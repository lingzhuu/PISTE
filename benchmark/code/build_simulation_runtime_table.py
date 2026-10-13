"""Supplementary table of run time and compute for the simulation benchmark.

Reads only saved outputs and logs of a finished run (nothing is rerun):
  SCT              $seconds in the SCT log (preprocess_manifest)
  GSVA             GSVA_SCT/runtime.csv + permutation_runtime.csv
  GSVA bootstrap   benchmark log: last RUN line to log modification time
  STpathway        benchmark log span
  GSDensity        'seconds:' line of the benchmark log
  GSDensity target benchmark log span
  SpaPath          fit_seconds of the latest successful fit per pathway (status.csv)
Usage: python build_simulation_runtime_table.py SIM_ROOT OUTPUT.xlsx
"""
import csv
import datetime as dt
import os
import re
import statistics as st
import sys
from pathlib import Path

from openpyxl import Workbook
from openpyxl.styles import Alignment, Border, Font, Side

root, out = Path(sys.argv[1]), Path(sys.argv[2])
logs = root / "benchmark_logs"
seeds = list(range(2020, 2027))


def text(path):
    return path.read_bytes().decode("utf-8", "replace") if path.is_file() else ""


def read(path):
    with path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def log_span(label):
    path = logs / f"{label}.log"
    runs = re.findall(r"^RUN (\S+) attempt", text(path), re.M)
    start = dt.datetime.fromisoformat(runs[-1])
    return (dt.datetime.fromtimestamp(os.path.getmtime(path)) - start).total_seconds()


def quantile(v, p):
    v = sorted(v)
    k = (len(v) - 1) * p
    i = int(k)
    return v[i] + (v[min(i + 1, len(v) - 1)] - v[i]) * (k - i)


seconds = {}  # (method, seed) -> seconds
fits = []     # SpaPath per-pathway fits
for s in seeds:
    base = root / f"seed{s}"
    seconds["SCT", s] = float(re.search(r"\$seconds\s*\n\[1\]\s*([\d.]+)", text(logs / f"{s}_SCT.log")).group(1))
    seconds["GSVA", s] = sum(float(r["seconds"]) for f in ("runtime.csv", "permutation_runtime.csv")
                             for r in read(base / "GSVA_SCT" / f))
    seconds["GSVA_bootstrap", s] = log_span(f"{s}_GSVA_bootstrap")
    seconds["STpathway", s] = log_span(f"{s}_STpathway")
    seconds["GSDensity", s] = float(re.findall(r"^seconds:\s*([\d.]+)", text(logs / f"{s}_GSDensity.log"), re.M)[-1])
    seconds["GSDensity_target", s] = log_span(f"{s}_GSDensity_target")
    truth = {r["pathway"]: r for r in read(base / "pathway_truth.csv")}
    latest = {r["pathway"]: float(r["fit_seconds"])
              for r in read(base / "SpaPath_dev_0.6.3" / "Gaussian" / "status.csv") if r["status"] == "ok"}
    assert set(latest) == set(truth), f"seed {s}: SpaPath fits incomplete ({len(latest)}/50)"
    seconds["SpaPath", s] = sum(latest.values())
    fits += [dict(truth[p], seed=s, seconds=v) for p, v in latest.items()]

METHODS = [
    ("SCT", "SCTransform v2 (shared preprocessing)", "Raw counts",
     "1 data set (10,000 genes x 3,639 spots)", "1"),
    ("GSVA", "GSVA, target vs rest", "SCT Pearson residuals",
     "50 gene sets; 1,999 spot-label permutations", "1"),
    ("GSVA_bootstrap", "GSVA spatial block bootstrap", "Saved GSVA scores",
     "50 gene sets x 3 block widths x 1,999 draws", "1"),
    ("STpathway", "STpathway predefined-area test", "SCT Pearson residuals",
     "50 gene sets; dcov.test with 199 permutations", "1"),
    ("GSDensity", "GSDensity (MCA, coherence, RWR, spatial KLD)", "SCT corrected counts",
     "50 gene sets; 100 grids x 100 resamples", "1"),
    ("GSDensity_target", "GSDensity regional specificity", "GSDensity PAL",
     "Gene sets passing coherence q < 0.05", "1"),
    ("SpaPath", "SpaPath 0.6.3 (PathwayLGM, INLA)", "SCT Pearson residuals",
     "50 gene sets, one joint fit of 100 genes each", "2"),
]
PEAK = {"SpaPath": "about 4.3"}  # Rscript peak working set 2.74 GB (sampled every 30 s) + inla.exe 1.6 GB (spot check)

TITLE = Font(name="Calibri", size=11, bold=True)
HEAD = Font(name="Calibri", size=10, bold=True)
BODY = Font(name="Calibri", size=10)
NOTE = Font(name="Calibri", size=9, italic=False)
LINE = Border(bottom=Side(style="thin"))
TOP = Border(top=Side(style="thin"))
WRAP = Alignment(wrap_text=True, vertical="center")
CENTER = Alignment(horizontal="center", vertical="center", wrap_text=True)


def sheet(ws, title, header, rows, notes, widths, center_from=2):
    n = len(header)
    ws.append([title])
    ws.merge_cells(start_row=1, start_column=1, end_row=1, end_column=n)
    ws.cell(1, 1).font = TITLE
    ws.append(header)
    for c in range(1, n + 1):
        cell = ws.cell(2, c)
        cell.font, cell.border, cell.alignment = HEAD, LINE, CENTER
    for row in rows:
        ws.append(row)
        r = ws.max_row
        for c in range(1, n + 1):
            ws.cell(r, c).font = BODY
            ws.cell(r, c).alignment = CENTER if c >= center_from else WRAP
    for c in range(1, n + 1):
        ws.cell(ws.max_row, c).border = LINE
    ws.append([])
    for note in notes:
        ws.append([note])
        r = ws.max_row
        ws.merge_cells(start_row=r, start_column=1, end_row=r, end_column=n)
        ws.cell(r, 1).font = NOTE
        ws.cell(r, 1).alignment = Alignment(wrap_text=True, vertical="top")
        ws.row_dimensions[r].height = 15 * max(1, len(note) // (sum(widths) * 1.1) + 1)
    for i, w in enumerate(widths, 1):
        ws.column_dimensions[chr(64 + i)].width = w
    ws.freeze_panes = "A3"
    ws.sheet_view.showGridLines = False


def fmt_min(v):
    return f"{v / 60:.2f}"


COMMON_NOTES = [
    "Hardware: Intel Core i7-13700K (16 cores, 24 threads), 32 GB RAM, Windows 11 x64. Software: R 4.6.1, "
    "Seurat 5.5.1, sctransform 0.4.3, GSVA 2.6.6, energy 1.7-12 (STpathway), gsdensity 0.1.0 with CelliD 1.19.0, "
    "PathwayLGM 0.6.3 (SpaPath) with INLA 26.08.07; Python 3.12 with NumPy 2.5 for the bootstrap.",
    "All jobs ran one at a time, each restricted to two CPU cores. SpaPath used two threads "
    "(OMP_NUM_THREADS = 2, INLA threads '2:1'); all other methods were single-threaded (BLAS threads = 1).",
    "Each simulated data set has 10,000 genes, 3,639 spots and 50 gene sets of 100 genes; 7 data sets (seeds 2020-2026).",
]

wb = Workbook()
ws = wb.active
ws.title = "Runtime summary"
rows = []
for key, name, inp, unit, threads in METHODS:
    v = [seconds[key, s] for s in seeds]
    rows.append([name, inp, unit, threads, f"{st.median(v) / 60:.2f} ({min(v) / 60:.2f}-{max(v) / 60:.2f})",
                 round(sum(v) / 3600, 2) if sum(v) >= 36 else "<0.01", PEAK.get(key, "NA")])
total = sum(seconds.values())
sp = sum(seconds["SpaPath", s] for s in seeds)
sheet(ws, "Supplementary Table X. Run time and computing resources of each method on the simulated data sets.",
      ["Method", "Input", "Work per data set", "Threads", "Wall time per data set, min: median (range)",
       "Total wall time, 7 data sets (h)", "Peak memory (GB)"],
      rows,
      [f"Wall time is elapsed time per data set; SpaPath is the sum of its 50 per-gene-set fits. "
       f"SpaPath accounted for {100 * sp / total:.1f}% of the total {total / 3600:.1f} h.",
       "Timing sources: SCT, GSVA, GSDensity and SpaPath from the methods' own timers; GSVA bootstrap, STpathway "
       "and GSDensity regional specificity from the start and end times of their logs (includes R start-up and data loading).",
       "Peak memory: largest working set of the R process (sampled every 30 s over the run, 2.7 GB) plus the INLA child process (1.6 GB in a spot check); "
       "NA, not recorded for the other methods."] + COMMON_NOTES,
      [36, 20, 34, 9, 22, 16, 12], center_from=4)

ws = wb.create_sheet("Per data set")
order = [m for m, *_ in METHODS]
names = {m: n for m, n, *_ in METHODS}
rows = [[s] + [round(seconds[m, s] / 60, 2) for m in order] for s in seeds]
rows.append(["Total"] + [round(sum(seconds[m, s] for s in seeds) / 60, 2) for m in order])
sheet(ws, "Supplementary Table X (continued). Wall time (min) of each method for each simulated data set.",
      ["Seed"] + [names[m] for m in order], rows, COMMON_NOTES, [9] + [17] * len(order))
for c in range(1, len(order) + 2):
    ws.cell(2 + len(rows), c).font = HEAD

ws = wb.create_sheet("SpaPath per fit")
GROUPS = [("K0", "0 (null)"), ("K1_10", "1-10"), ("K11_30", "11-30"), ("K31_60", "31-60"), ("K61_100", "61-100")]
PROFILES = [("broad_weak", "broad, weak-effect"), ("narrow_strong", "narrow, strong-effect")]
rows = []
for pkey, pname in PROFILES:
    for gkey, gname in GROUPS:
        v = [f["seconds"] for f in fits if f["effect_profile"] == pkey and f["response_group"] == gkey]
        rows.append([pname, gname, len(v), f"{st.median(v):.0f} ({quantile(v, .25):.0f}-{quantile(v, .75):.0f})",
                     f"{min(v):.0f}-{max(v):.0f}"])
v = [f["seconds"] for f in fits]
rows.append(["All", "All", len(v), f"{st.median(v):.0f} ({quantile(v, .25):.0f}-{quantile(v, .75):.0f})",
             f"{min(v):.0f}-{max(v):.0f}"])
sheet(ws, "Supplementary Table X (continued). SpaPath fitting time per gene set.",
      ["Effect distribution", "Responding genes K", "Gene sets", "Seconds per fit: median (IQR)", "Range (s)"],
      rows,
      ["One fit is the joint INLA model of the 100 genes of one gene set (Gaussian observed mode, oracle "
       "target-versus-rest regions, SPDE mesh with max edge 0.25 and cutoff 0.08, CCD integration, two threads)."]
      + COMMON_NOTES, [22, 18, 11, 26, 13], center_from=2)

wb.save(out)
print("saved", out, "total h", round(total / 3600, 2), "SpaPath share", round(100 * sp / total, 1))
