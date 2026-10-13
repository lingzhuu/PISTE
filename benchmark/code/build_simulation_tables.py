"""Fill Table 1 (simulation_summary.xlsx) and Supplementary Table 1
(simulation_per_pathway_result.xlsx) with the results of a simulation run.

The two workbooks are used as templates: titles, headers, notes and formatting
stay as they are and only the result cells are rewritten.
Input: tables_input/per_pathway_all.csv written by tables_input_v3.R.
Usage: python build_simulation_tables.py per_pathway_all.csv RESULT_DIR
"""
import csv
import statistics as st
import sys
from pathlib import Path

from openpyxl import load_workbook

rows = list(csv.DictReader(open(sys.argv[1], encoding="utf-8")))
out = Path(sys.argv[2])
assert len(rows) == 350
seeds = sorted({int(r["seed"]) for r in rows})
GROUPS = ["K0", "K1_10", "K11_30", "K31_60", "K61_100"]
PROFILES = ["broad_weak", "narrow_strong"]
PROFILE_LABEL = {"broad_weak": "broad, weak-effect", "narrow_strong": "narrow, strong-effect"}


def num(x):
    return None if x in ("", "NA", None) else float(x)


def small(x, digits):
    """x.xxe-n for tiny values (as in the templates), else fixed decimals."""
    m, e = f"{x:.2e}".split("e")
    return f"{m}e-{abs(int(e))}"


def prob(x):
    """Probability with 4 decimals; tails beyond 4 decimals as 1-x.xxe-n or x.xxe-n."""
    if 1 - x < 5e-5:
        return "1.0000" if 1 - x <= 1e-15 else "1−" + small(1 - x, 2)
    if x < 5e-5:
        return small(x, 2)
    return f"{x:.4f}"


def pval(x):
    return f"{x:.2e}" if x < 1e-3 else f"{x:.4f}"


def tested(r):
    return num(r["gsd_chisq_p"]) is not None


# ---------------------------------------------------------- Supplementary Table 1
wb = load_workbook(out / "simulation_per_pathway_result.xlsx")
assert wb.sheetnames == [f"Seed {s}" for s in seeds], wb.sheetnames
for s in seeds:
    ws = wb[f"Seed {s}"]
    d = sorted((r for r in rows if int(r["seed"]) == s), key=lambda r: r["pathway"])
    assert len(d) == 50
    for i, r in enumerate(d):
        row = 3 + i
        assert ws.cell(row, 1).value == r["pathway"]
        passed = num(r["gsd_deviation_q"]) < .05
        vals = [r["pathway"], r["spatial_class"], PROFILE_LABEL[r["effect_profile"]],
                f'{r["K"]} ({r["n_up"]}/{r["n_down"]})',
                f'{float(r["spapath_pi"]):.4f} / {prob(float(r["spapath_evidence"]))}',
                float(r["stpathway_dcor"]), float(r["gsva_delta"]),
                f'{float(r["gsva_boot_p"]):.4f} / {float(r["gsva_boot_q"]):.4f}',
                float(r["gsd_deviation_q"]),
                float(r["gsd_spatial_log_kld"]) if passed and num(r["gsd_spatial_log_kld"]) is not None else "NA",
                f'{int(float(r["gsd_high_target"]))}/{int(float(r["gsd_n_target"]))} / '
                f'{int(float(r["gsd_high_rest"]))}/{int(float(r["gsd_n_rest"]))}' if tested(r) else "NA",
                pval(float(r["gsd_chisq_p"])) if tested(r) else "NA",
                pval(float(r["gsd_chisq_q"])) if tested(r) else "NA"]
        for c, v in enumerate(vals, 1):
            ws.cell(row, c).value = v
wb.save(out / "simulation_per_pathway_result.xlsx")

# ---------------------------------------------------------------------- Table 1
wb = load_workbook(out / "simulation_summary.xlsx")
ws, wd = wb["Median summary"], wb["Detection counts"]
med = lambda v: st.median(v) if v else None
for i, (p, g) in enumerate((p, g) for p in PROFILES for g in GROUPS):
    d = [r for r in rows if r["effect_profile"] == p and r["response_group"] == g]
    assert len(d) == 35
    ps = [r for r in d if float(r["gsd_deviation_q"]) < .05]
    tt = [r for r in ps if tested(r)]
    kld = [float(r["gsd_spatial_log_kld"]) for r in ps if num(r["gsd_spatial_log_kld"]) is not None]
    row = 3 + i
    ws.cell(row, 3).value = (f'{med([float(r["spapath_pi"]) for r in d]):.4f} / '
                             f'{prob(med([float(r["spapath_evidence"]) for r in d]))}')
    ws.cell(row, 4).value = med([float(r["stpathway_dcor"]) for r in d])
    ws.cell(row, 5).value = med([float(r["gsva_delta"]) for r in d])
    ws.cell(row, 6).value = (f'{med([float(r["gsva_boot_p"]) for r in d]):.4f} / '
                             f'{med([float(r["gsva_boot_q"]) for r in d]):.4f}')
    ws.cell(row, 7).value = f"{len(ps)} / 35"
    ws.cell(row, 8).value = med(kld) if kld else "NA"
    if tt:
        ft = med([float(r["gsd_high_target"]) / float(r["gsd_n_target"]) for r in tt])
        fr = med([float(r["gsd_high_rest"]) / float(r["gsd_n_rest"]) for r in tt])
        ws.cell(row, 9).value = f"{ft:.4f} / {fr:.4f}"
        ws.cell(row, 10).value = (f'{pval(med([float(r["gsd_chisq_p"]) for r in tt]))} / '
                                  f'{pval(med([float(r["gsd_chisq_q"]) for r in tt]))}')
    else:
        ws.cell(row, 9).value = ws.cell(row, 10).value = "NA"
    wd.cell(row, 3).value = f'{sum(float(r["gsva_boot_q"]) < .05 for r in d)} / 35'
    wd.cell(row, 4).value = f'{sum(tested(r) and float(r["gsd_chisq_q"]) < .05 for r in d)} / 35'
wb.save(out / "simulation_summary.xlsx")
print("filled", out)
