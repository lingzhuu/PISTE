"""Supplementary table of SpaPath run time on the simulated data sets.

Reads only the saved fits of a finished run (nothing is rerun): fit_seconds of
the latest successful fit of every gene set in
seed*/SpaPath_dev_0.6.3/Gaussian/status.csv, and pathway_truth.csv for the groups.
Usage: python build_simulation_runtime_table.py SIM_ROOT OUTPUT.xlsx
"""
import csv
import statistics as st
import sys
from pathlib import Path

from openpyxl import Workbook
from openpyxl.styles import Alignment, Border, Font, Side

root, out = Path(sys.argv[1]), Path(sys.argv[2])
seeds = list(range(2020, 2027))


def read(path):
    with path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def quantile(v, p):
    v = sorted(v)
    k = (len(v) - 1) * p
    i = int(k)
    return v[i] + (v[min(i + 1, len(v) - 1)] - v[i]) * (k - i)


def mid(v):
    return f"{st.median(v):.0f} ({quantile(v, .25):.0f}-{quantile(v, .75):.0f})"


fits = []
for s in seeds:
    base = root / f"seed{s}"
    truth = {r["pathway"]: r for r in read(base / "pathway_truth.csv")}
    latest = {r["pathway"]: float(r["fit_seconds"])
              for r in read(base / "SpaPath_dev_0.6.3" / "Gaussian" / "status.csv") if r["status"] == "ok"}
    assert set(latest) == set(truth), f"seed {s}: SpaPath fits incomplete ({len(latest)}/50)"
    fits += [dict(truth[p], seed=s, seconds=v) for p, v in latest.items()]

TITLE = Font(name="Calibri", size=11, bold=True)
HEAD = Font(name="Calibri", size=10, bold=True)
BODY = Font(name="Calibri", size=10)
NOTE = Font(name="Calibri", size=9)
LINE = Border(bottom=Side(style="thin"))
CENTER = Alignment(horizontal="center", vertical="center", wrap_text=True)
LEFT = Alignment(vertical="center", wrap_text=True)


def sheet(ws, title, header, rows, notes, widths, bold_last=False):
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
            ws.cell(r, c).alignment = CENTER if c > 1 else LEFT
    last = ws.max_row
    for c in range(1, n + 1):
        ws.cell(last, c).border = LINE
        if bold_last:
            ws.cell(last, c).font = HEAD
    ws.append([])
    for note in notes:
        ws.append([note])
        r = ws.max_row
        ws.merge_cells(start_row=r, start_column=1, end_row=r, end_column=n)
        ws.cell(r, 1).font = NOTE
        ws.cell(r, 1).alignment = Alignment(wrap_text=True, vertical="top")
        ws.row_dimensions[r].height = 13 * (len(note) // int(sum(widths) * 1.1) + 1)
    for i, w in enumerate(widths, 1):
        ws.column_dimensions[chr(64 + i)].width = w
    ws.freeze_panes = "A3"
    ws.sheet_view.showGridLines = False


NOTES = [
    "One fit is the joint model of the 100 genes of one gene set (Gaussian observed mode, SCT Pearson residuals, "
    "oracle target-versus-rest regions, SPDE mesh with maximum edge 0.25 and cutoff 0.08, CCD integration). "
    "Time is the elapsed time of the fit and of computing pi and evidence; reading the data is not included.",
    "Each simulated data set has 10,000 genes, 3,639 spots and 50 gene sets of 100 genes; 7 data sets (seeds 2020-2026). "
    "Gene sets were fitted one after another.",
    "Hardware and software: Intel Core i7-13700K, 32 GB RAM, Windows 11 x64; R 4.6.1, PathwayLGM 0.6.3, INLA 26.08.07. "
    "Peak memory of one fit was about 4.3 GB (R process 2.7 GB, sampled every 30 s, plus the INLA process, about 1.6 GB in a spot check).",
]

wb = Workbook()
ws = wb.active
ws.title = "Per data set"
rows = []
for s in seeds:
    v = [f["seconds"] for f in fits if f["seed"] == s]
    rows.append([s, len(v), round(sum(v) / 60, 1), mid(v)])
v = [f["seconds"] for f in fits]
rows.append(["All", len(v), round(sum(v) / 60, 1), mid(v)])
sheet(ws, "Supplementary Table X. SpaPath run time on the simulated data sets.",
      ["Seed", "Gene sets", "Total time (min)", "Seconds per gene set: median (IQR)"],
      rows, NOTES, [10, 11, 16, 30], bold_last=True)

ws = wb.create_sheet("By group")
GROUPS = [("K0", "0 (null)"), ("K1_10", "1-10"), ("K11_30", "11-30"), ("K31_60", "31-60"), ("K61_100", "61-100")]
PROFILES = [("broad_weak", "broad, weak-effect"), ("narrow_strong", "narrow, strong-effect")]
rows = []
for pkey, pname in PROFILES:
    for gkey, gname in GROUPS:
        v = [f["seconds"] for f in fits if f["effect_profile"] == pkey and f["response_group"] == gkey]
        rows.append([pname, gname, len(v), mid(v), f"{min(v):.0f}-{max(v):.0f}"])
sheet(ws, "Supplementary Table X (continued). SpaPath run time per gene set by effect distribution and number of responding genes K.",
      ["Effect distribution", "Responding genes K", "Gene sets", "Seconds per gene set: median (IQR)", "Range (s)"],
      rows, NOTES, [22, 18, 11, 30, 13])

wb.save(out)
print("saved", out, "fits", len(fits), "total h", round(sum(f["seconds"] for f in fits) / 3600, 2))
