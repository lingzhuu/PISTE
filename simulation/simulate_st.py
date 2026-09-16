"""Simulate spatial transcriptomics data with pathway-level ground truth.

Design:
- 10,000 genes: 50 test pathways x 100 genes + 5,000 background genes
- 5,000 background genes: 20% SVGs (1000) with independently generated patterns, 80% non-SVGs.
- Categories: 25 coherent (regional/gradient/boundary/multi-hotspot x5 each + 5 mixed), 
  15 heterogeneous (all regional; direction ratios 50:50 / 60:40 / 75:25, five each), 10 null.
- Visium-like hexagonal lattice, 100 um center-to-center (~70x70 = 4,900 spots).
- Raw counts:  Y_ig ~ NB(m_g (1 + d_g b_g H_ig), theta_g),
  m_g/theta_g sampled from empirical reference (calibrate_reference.py).
- Gaussian:    X_ig ~ N(a_g + d_g delta_g H_ig, sigma_g^2),
  delta_g = r_g sigma_g, a_g/sigma2_g from the same reference.
"""


import json
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.spatial import cKDTree

# ---------------- configuration  ----------------
CFG = dict(
    seed=20260915,
    n_genes=10000,
    n_pathways=50,
    pathway_size=100,
    frac_responsive=0.30,          # fraction of responsive genes per pathway
    bg_frac_svg=0.20,              # fraction of SVGs among background genes
    hex_rows=70, hex_cols=70, spacing=100.0,   # um, Visium-like
    coherent_patterns=(["regional"] * 5 + ["gradient"] * 5 + ["boundary"] * 5
                       + ["multi_hotspot"] * 5 + ["mixed"] * 5),
    n_hotspots=3,                  # for multi-hotspot pattern
    het_ratios=[0.50] * 5 + [0.60] * 5 + [0.75] * 5,   # all regional
    n_null=10,
    beta_ab=(2.0, 2.0),            # b_g ~ Beta(2,2): median 0.5, ~50% per-gene detection
    beta_ar=(2.0, 5.0),            # r_g ~ Beta(2,5): median 0.26
    ref_params="data/reference_params.npz",
    outdir="output",
    chunk=1000,                    # genes per sampling chunk
)

work_dir = Path("../data")
output_dir = work_dir / CFG["outdir"]
output_dir.mkdir(exist_ok=True)
rng = np.random.default_rng(CFG["seed"])

# Check that the numbers of coherent, heterogeneous, and null pathways add up to the total number of pathways
assert len(CFG["coherent_patterns"]) + len(CFG["het_ratios"]) + CFG["n_null"] == CFG["n_pathways"]
N_PW_GENES = CFG["n_pathways"] * CFG["pathway_size"]
N_BG = CFG["n_genes"] - N_PW_GENES
N_BG_SVG = int(round(N_BG * CFG["bg_frac_svg"]))

# ---------------- reference parameters ----------------
ref = np.load(work_dir / CFG["ref_params"])
M_EMP, TH_EMP = ref["m"], ref["theta"]
A_EMP, S2_EMP = ref["a"], ref["sigma2"]

# ---------------- spatial coordinates: hexagonal lattice ----------------
rows, cols = np.meshgrid(
    np.arange(CFG["hex_rows"]),
    np.arange(CFG["hex_cols"]),
    indexing="ij"
)
rows, cols = rows.ravel(), cols.ravel()
# Convert row/column indices into 2D hexagonal-grid coordinates
coords = np.c_[
    cols * CFG["spacing"] + (rows % 2) * CFG["spacing"] / 2,
    rows * CFG["spacing"] * np.sqrt(3) / 2
]
N_SPOTS = len(coords)
EXTENT = np.ptp(coords, axis=0)
print(f"spots={N_SPOTS}, extent={EXTENT[0]:.0f} x {EXTENT[1]:.0f} um")

# ---------------- spatial pattern generators -> A in [0,1]^n_spots ----------------
# Randomly choose a center point within the spatial domain
def _rand_center(c):
    lo, hi = c.min(axis=0), c.max(axis=0)
    return lo + rng.random(2) * (hi - lo)

def pat_regional(c):
    ctr, r = _rand_center(c), rng.uniform(0.15, 0.30) * EXTENT.max()
    d = np.linalg.norm(c - ctr, axis=1) # L-2 norm
    # Convert distance from the regional center into a smooth activity score between 0 and 1
    return 1.0 / (1.0 + np.exp((d - r) / 150.0))

def pat_gradient(c):
    u = rng.normal(size=2)
    u /= np.linalg.norm(u)
    t = c @ u
    return (t - t.min()) / (t.max() - t.min())

def pat_boundary(c):
    u = rng.normal(size=2)
    u /= np.linalg.norm(u)
    mid = (c.max(axis=0) + c.min(axis=0)) / 2
    d = np.abs((c - mid) @ u - rng.uniform(-500, 500))
    s = rng.uniform(200, 400)
    return np.exp(-d**2 / (2 * s**2))

def pat_multi_hotspot(c, n_hotspots=3):
    A = np.zeros(len(c))
    for _ in range(n_hotspots):
        ctr, r = _rand_center(c), rng.uniform(0.06, 0.12) * EXTENT.max()
        d = np.linalg.norm(c - ctr, axis=1)
        A = np.maximum(A, 1.0 / (1.0 + np.exp((d - r) / 100.0)))
    return A

# Generate a mixed pattern combining regional and gradient effects
def pat_mixed(c, regional_weight=0.5, gradient_weight=0.5, boundary_weight=0, hotspot_weight=0):
    return np.clip(
        regional_weight * pat_regional(c) +
        gradient_weight * pat_gradient(c) +
        boundary_weight * pat_boundary(c) +
        hotspot_weight * pat_multi_hotspot(c),
        0,
        1
    )

PATTERNS = {"regional": pat_regional, "gradient": pat_gradient,
            "boundary": pat_boundary, "multi_hotspot": pat_multi_hotspot,
            "mixed": pat_mixed}

# ---------------- gene-level arrays ----------------
G = CFG["n_genes"]
m_g = np.empty(G)                 # NB baseline mean
theta_g = np.empty(G)             # NB dispersion
a_g = np.empty(G)                 # Gaussian baseline
s2_g = np.empty(G)                # Gaussian variance
d_g = np.zeros(G)                 # response direction in {-1, +1}
b_g = np.zeros(G)                 # NB spatial effect magnitude
r_g = np.zeros(G)                 # Gaussian standardized effect
H = np.zeros((N_SPOTS, G), dtype=np.float32)
gt = []                           # ground-truth rows

A_B, B_B = CFG["beta_ab"]
A_R, B_R = CFG["beta_ar"]

def beta_b(k):
    return rng.beta(A_B, B_B, size=k)

def beta_r(k):
    return rng.beta(A_R, B_R, size=k)

def fill_basic(idx):
    m_g[idx] = rng.choice(M_EMP, size=len(idx), replace=True)
    theta_g[idx] = rng.choice(TH_EMP, size=len(idx), replace=True)
    a_g[idx] = rng.choice(A_EMP, size=len(idx), replace=True)
    s2_g[idx] = rng.choice(S2_EMP, size=len(idx), replace=True)

K_RESP = int(round(CFG["pathway_size"] * CFG["frac_responsive"]))

# ---------------- test pathways ----------------
for p in range(CFG["n_pathways"]):
    idx = np.arange(p * CFG["pathway_size"], (p + 1) * CFG["pathway_size"])
    fill_basic(idx)
    pw = f"PW{p + 1:02d}"

    if p < len(CFG["coherent_patterns"]):
        cat, pattern = "coherent", CFG["coherent_patterns"][p]
        ratio = ""
    elif p < len(CFG["coherent_patterns"]) + len(CFG["het_ratios"]):
        cat, pattern = "heterogeneous", "regional"
        ratio = CFG["het_ratios"][p - len(CFG["coherent_patterns"])]
    else:
        cat, pattern, ratio = "null", "", ""

    if cat == "null":
        for g in idx:
            gt.append((f"Gene_{g + 1:05d}", pw, cat, "", "", False, 0, 0.0, 0.0))
        continue

    A = PATTERNS[pattern](coords)
    resp = rng.choice(CFG["pathway_size"], size=K_RESP, replace=False)
    resp_g = idx[resp]

    if cat == "coherent":
        b_g[resp_g] = beta_b(K_RESP)
        r_g[resp_g] = beta_r(K_RESP)
        d_g[resp_g] = 1.0
        H[:, resp_g] = A[:, None]
    else:  # heterogeneous: matched pairs share (m, b)
        n_pos = int(round(K_RESP * ratio))
        n_neg = K_RESP - n_pos
        n_pair = min(n_pos, n_neg)
        pos_g, neg_g = resp_g[:n_pos], resp_g[n_pos:]
        # paired genes: identical m and b -> contributions cancel
        for k in range(n_pair):
            gm, gb = rng.choice(M_EMP, size=1)[0], beta_b(1)[0]
            m_g[pos_g[k]] = m_g[neg_g[k]] = gm
            b_g[pos_g[k]] = b_g[neg_g[k]] = gb
            r_g[pos_g[k]] = r_g[neg_g[k]] = beta_r(1)[0]
        # excess majority genes: independent draws
        extra = pos_g[n_pair:]
        if len(extra):
            b_g[extra] = beta_b(len(extra))
            r_g[extra] = beta_r(len(extra))
        d_g[pos_g], d_g[neg_g] = 1.0, -1.0
        H[:, resp_g] = A[:, None]

    resp_set = set(resp_g.tolist())
    for g in idx:
        gt.append((f"Gene_{g + 1:05d}", pw, cat, pattern,
                   ratio if cat == "heterogeneous" else "",
                   g in resp_set, int(d_g[g]) if g in resp_set else 0,
                   float(b_g[g]), float(r_g[g])))

# ---------------- background genes ----------------
bg_idx = np.arange(N_PW_GENES, G)
fill_basic(bg_idx)
svg_g = rng.choice(bg_idx, size=N_BG_SVG, replace=False)
for g in svg_g:
    pat = PATTERNS[rng.choice(["regional", "gradient", "boundary", "multi_hotspot"])]
    H[:, g] = pat(coords)
    b_g[g] = beta_b(1)[0]
    r_g[g] = beta_r(1)[0]
    d_g[g] = 1.0
bg_set = set(svg_g.tolist())
for g in bg_idx:
    gt.append((f"Gene_{g + 1:05d}", "background",
               "background_svg" if g in bg_set else "background_nonsvg",
               "", "", g in bg_set, int(d_g[g]), float(b_g[g]), float(r_g[g])))

# ---------------- raw-count simulation (Gamma-Poisson mixture = NB) ----------------
blocks = []
for s in range(0, G, CFG["chunk"]):
    e = min(s + CFG["chunk"], G)
    MU = m_g[s:e] * (1.0 + d_g[s:e] * b_g[s:e] * H[:, s:e])
    th = theta_g[s:e]
    lam = rng.gamma(shape=th[None, :], scale=np.maximum(MU, 1e-12) / th[None, :])
    blocks.append(sp.csr_matrix(rng.poisson(lam).astype(np.int32)))
Y = sp.hstack(blocks, format="csr")
print(f"Y: {Y.shape}, median depth/spot={np.median(Y.sum(axis=1).A.ravel()):.0f}")

# ---------------- Gaussian simulation (shared H, independent parameters) ----------------
delta = r_g * np.sqrt(s2_g)
X = (a_g[None, :] + d_g[None, :] * delta[None, :] * H
     + rng.normal(size=(N_SPOTS, G)) * np.sqrt(s2_g)[None, :]).astype(np.float32)

# ---------------- save ----------------
pd.DataFrame(coords, columns=["x_um", "y_um"]).to_csv(output_dir / "coords.csv", index=False)
pd.DataFrame(gt, columns=["gene", "set", "category", "pattern", "ratio",
                          "responsive", "direction", "b", "r"]).to_csv(
    output_dir / "ground_truth.csv", index=False)
sp.save_npz(output_dir / "counts_raw.npz", Y)
sp.save_npz(output_dir / "H_maps.npz", sp.csr_matrix(H))
np.savez_compressed(output_dir / "gaussian.npz", X=X)
json.dump({k: (list(v) if isinstance(v, tuple) else v) for k, v in CFG.items()},
          open(output_dir / "config.json", "w"), indent=2)

# ---------------- sanity checks ----------------
def morans_I(v, W):
    z = v - v.mean()
    return float(z @ (W @ z) / max(z @ z, 1e-12))


tree = cKDTree(coords)
_, nn = tree.query(coords, k=7)
W = sp.csr_matrix((np.ones(N_SPOTS * 6),
                   (np.repeat(np.arange(N_SPOTS), 6), nn[:, 1:].ravel())),
                  shape=(N_SPOTS, N_SPOTS))
W = sp.diags(1.0 / W.sum(axis=1).A.ravel()) @ (W + W.T > 0).astype(float)

depth = np.asarray(Y.sum(axis=1)).ravel()
L = np.log1p(Y.multiply(1e4 / depth[:, None]).toarray())

pw_idx = lambda p: np.arange(p * 100, (p + 1) * 100)
mu_path = lambda p: (m_g[pw_idx(p)] * (1 + d_g[pw_idx(p)] * b_g[pw_idx(p)]
                     * H[:, pw_idx(p)])).mean(axis=1)
mean_path = lambda p: L[:, pw_idx(p)].mean(axis=1)

p_coh, p_h50, p_h75, p_null = 0, 25, 35, 45
flat = mu_path(p_h50)

# ---------------- sanity figure ----------------
import matplotlib
matplotlib.use("Agg")
import sys
sys.path.insert(0, str(Path(sys.executable).parent.parent.parent))
from daimon_runtime import setup_plot
import matplotlib.pyplot as plt
setup_plot()

panels = [("pattern: regional", PATTERNS["regional"](coords)),
          ("pattern: gradient", PATTERNS["gradient"](coords)),
          ("pattern: boundary", PATTERNS["boundary"](coords)),
          ("pattern: multi-hotspot", PATTERNS["multi_hotspot"](coords)),
          ("pattern: mixed", PATTERNS["mixed"](coords)),
          ("coherent pathway mean", mean_path(p_coh)),
          ("hetero 50:50 mean (cancels)", mean_path(p_h50)),
          ("hetero 75:25 mean", mean_path(p_h75)),
          ("null pathway mean", mean_path(p_null))]
fig, axes = plt.subplots(3, 3, figsize=(13, 12))
for ax, (title, v) in zip(axes.ravel(), panels):
    sc = ax.scatter(coords[:, 0], coords[:, 1], c=v, s=3, cmap="viridis")
    ax.set_title(title, fontsize=11)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.colorbar(sc, ax=ax, fraction=0.046)
fig.suptitle("Simulation sanity check (10k genes, 50 pathways, 4,900 spots)", fontsize=13)
fig.savefig(output_dir / "sanity_check.png", bbox_inches="tight", dpi=150)
print(f"saved outputs to {output_dir}")
