"""
Builds (1) a formatted p-value / BH workbook and (2) publication figures
from stats_results.xlsx (produced by stats_tests.py).

Outputs
  Statistical_tests_summary.xlsx
  figures/Fig1_variant_mismatch_heatmap.(png|pdf)
  figures/Fig2_mismatch_over_time.(png|pdf)
  figures/Fig3_pairwise_variant_BH.(png|pdf)
  figures/Fig4_effect_size_vs_significance.(png|pdf)
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm, BoundaryNorm, ListedColormap
from openpyxl import load_workbook
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.formatting.rule import CellIsRule, FormulaRule
from openpyxl.utils import get_column_letter

from stats_tests import NAME_MAP  # assay names (edit there)

SRC, OUT_XLSX, FIG = "stats_results.xlsx", "Statistical_tests_summary.xlsx", "figures"
os.makedirs(FIG, exist_ok=True)
ALPHA = 0.05
P_FLOOR = 1e-300            # numerical underflow limit
MC_FLOOR = 1 / 200_001      # Monte-Carlo resolution of the exact test

rd = lambda s: pd.read_excel(SRC, s)
v_omni, v_pair = rd("Variant_omnibus"), rd("Variant_pairwise_Fisher")
q_omni, q_cons, q_trend = rd("Quarter_omnibus"), rd("Quarter_consecutive_Fisher"), rd("Quarter_trend_CA")
mm_v = pd.read_excel(SRC, "Mismatch_counts_variant", index_col=0)
tot_v = pd.read_excel(SRC, "Totals_variant", index_col=0)
mm_q = pd.read_excel(SRC, "Mismatch_counts_quarter", index_col=0)
tot_q = pd.read_excel(SRC, "Totals_quarter", index_col=0)
COLS = list(tot_v.columns)
NAMES = [NAME_MAP[c] for c in COLS]
VARIANTS = ["Other", "Alpha", "Delta", "Omicron"]


def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    r = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    out = np.empty(n); out[o] = np.minimum(r, 1); return out


def primary(df):
    """Primary omnibus p: chi-square when expected counts are adequate, otherwise exact (FFH)."""
    use_exact = df["pct_cells_exp_lt5"] > 20
    p = np.where(use_exact, df["p_fisher_exact"], df["p_chi2"])
    test = np.where(use_exact, "Fisher-Freeman-Halton exact (MC)", "Pearson chi-square")
    p = np.maximum(p, P_FLOOR)
    return test, p


def stars(p):
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"


def fmt_p(p, exact_mc=False):
    if exact_mc and p <= MC_FLOOR * 1.01:
        return "<5.0E-06"
    if p <= P_FLOOR * 1e3:  # BH of an underflowed p
        return "<1.0E-300"
    return f"{p:.2E}" if p < 0.001 else f"{p:.3f}"


# =================================================================== tables
def omnibus_table(df, label):
    test, p = primary(df)
    padj = bh(p)
    out = pd.DataFrame({
        "Assay": [NAME_MAP[c] for c in df["column"]],
        "Comparison": label,
        "Test used": test,
        "Chi-square (χ²)": df["chi2"].round(2),
        "df": df["df"].astype(int),
        "Cramér's V": df["cramers_v"].round(3),
        "% cells expected <5": df["pct_cells_exp_lt5"].round(1),
        "p (chi-square)": np.maximum(df["p_chi2"], P_FLOOR),
        "p (exact, MC)": df["p_fisher_exact"],
        "p (primary)": p,
        "p BH-adjusted": padj,
        "p Bonferroni": np.minimum(p * len(p), 1),
        "Significant (BH<0.05)": np.where(padj < ALPHA, "Yes", "No"),
        "Significance": [stars(x) for x in padj],
    })
    return out


T_var = omnibus_table(v_omni, "Variants (Other/Alpha/Delta/Omicron)")
T_qtr = omnibus_table(q_omni, "Quarters (12 sampled quarters)")
for c in VARIANTS:
    T_var[f"Mismatch {c}"] = v_omni[f"mm_{c}"]

T_trend = pd.DataFrame({
    "Assay": [NAME_MAP[c] for c in q_trend["column"]], "Test": "Cochran-Armitage trend",
    "Z": q_trend["z"].round(2), "p": np.maximum(q_trend["p_trend"], P_FLOOR),
    "p BH-adjusted": np.maximum(q_trend["p_trend_BH"], P_FLOOR),
    "p Bonferroni": np.maximum(q_trend["p_trend_Bonf"], P_FLOOR),
    "Direction": np.where(q_trend["z"] > 0, "Increasing", "Decreasing"),
})
T_trend["Significant (BH<0.05)"] = np.where(T_trend["p BH-adjusted"] < ALPHA, "Yes", "No")
T_trend["Significance"] = [stars(x) for x in T_trend["p BH-adjusted"]]


def pair_table(df, g1, g2, pct1, pct2, extra=None):
    out = pd.DataFrame({
        "Assay": [NAME_MAP[c] for c in df["column"]],
        "Group 1": df[g1], "Group 2": df[g2],
        "Mismatch % group 1": df[pct1].round(3), "Mismatch % group 2": df[pct2].round(3),
        "Odds ratio (1 vs 2)": df["odds_ratio"], "OR 95% CI low": df["OR_95CI_low"],
        "OR 95% CI high": df["OR_95CI_high"],
        "p (Fisher exact)": np.maximum(df["p_fisher"], P_FLOOR),
        "p BH-adjusted": np.maximum(df["p_fisher_BH"], P_FLOOR),
        "p Bonferroni": np.maximum(df["p_fisher_Bonf"], P_FLOOR),
    })
    if extra is not None:
        out.insert(3, "Gap (quarters)", extra)
    out["Significant (BH<0.05)"] = np.where(out["p BH-adjusted"] < ALPHA, "Yes", "No")
    out["Significance"] = [stars(x) for x in out["p BH-adjusted"]]
    return out


T_vpair = pair_table(v_pair, "group1", "group2", "mm_pct_1", "mm_pct_2")
T_qcons = pair_table(q_cons, "from_q", "to_q", "mm_pct_from", "mm_pct_to", q_cons["gap_quarters"])
T_qcons = T_qcons.rename(columns={"Group 1": "From quarter", "Group 2": "To quarter"})

pct_v = (100 * mm_v.loc[VARIANTS] / tot_v.loc[VARIANTS]).T
pct_q = (100 * mm_q / tot_q).T
pct_v.index, pct_q.index = NAMES, NAMES

# ---------------------------------------------------------------- write xlsx
methods = [
    ["Statistical analysis of primer/probe mismatch frequencies"],
    [""],
    ["Data", "result.xlsx: Lineage / Lineage_Match and Date_Count / Date_match_Count. Mismatch = total − matched."],
    ["Variant test", "r×2 table (Other, Alpha, Delta, Omicron × mismatch/match) per assay. Pearson χ² when ≤20% of cells have expected count <5; otherwise Fisher–Freeman–Halton exact test (Monte-Carlo, 200,000 tables). Effect size: Cramér's V."],
    ["Quarter test", "12×2 table (2019.Q4–2023.Q1; 2020.Q2 and 2021.Q4 not sampled). Same rule as above. Plus Cochran–Armitage test for linear trend using calendar-spaced quarter scores."],
    ["Pairwise tests", "Two-sided Fisher's exact test (2×2). Odds ratio with Haldane–Anscombe correction (+0.5) and Woolf 95% CI."],
    ["Multiple testing", "Benjamini–Hochberg FDR (primary) and Bonferroni (sensitivity), applied separately within each family: variant omnibus (12), quarter omnibus (12), trend (12), pairwise variants (72), consecutive quarters (132)."],
    ["Significance codes", "*** BH p<0.001; ** <0.01; * <0.05; ns ≥0.05. Green cells = BH p<0.05."],
    ["Display limits", "p<1.0E-300 = below floating-point limit; exact Monte-Carlo p has resolution ≈5.0E-06."],
    ["Caution", "Sequences are a non-random GISAID sample with very unequal group sizes; with n=25,658, tiny differences reach significance — interpret together with Cramér's V / odds ratios."],
]
sheets = [
    ("README", None),
    ("1_Variant_omnibus", T_var),
    ("2_Variant_pairwise", T_vpair),
    ("3_Quarter_omnibus", T_qtr),
    ("4_Quarter_trend", T_trend),
    ("5_Quarter_consecutive", T_qcons),
    ("6_Mismatch_pct_variant", pct_v.round(3).reset_index(names="Assay")),
    ("7_Mismatch_pct_quarter", pct_q.round(3).reset_index(names="Assay")),
]
with pd.ExcelWriter(OUT_XLSX, engine="openpyxl") as w:
    pd.DataFrame(methods).to_excel(w, sheet_name="README", index=False, header=False)
    for name, df in sheets[1:]:
        df.to_excel(w, sheet_name=name, index=False)

wb = load_workbook(OUT_XLSX)
hdr_fill = PatternFill("solid", fgColor="1C5CAB")
sig_fill = PatternFill("solid", fgColor="D5F0DC")
ns_fill = PatternFill("solid", fgColor="F2F2F2")
thin = Side(style="thin", color="D0D0D0")
ws = wb["README"]
ws["A1"].font = Font(bold=True, size=14)
ws.column_dimensions["A"].width = 20; ws.column_dimensions["B"].width = 120
for row in ws.iter_rows(min_row=3):
    row[0].font = Font(bold=True)
    for c in row:
        c.alignment = Alignment(wrap_text=True, vertical="top")

for name, df in sheets[1:]:
    ws = wb[name]
    ws.freeze_panes = "B2"
    for j, col in enumerate(df.columns, 1):
        cell = ws.cell(1, j)
        cell.font = Font(bold=True, color="FFFFFF"); cell.fill = hdr_fill
        cell.alignment = Alignment(wrap_text=True, horizontal="center", vertical="center")
        width = max(len(str(col)), *(len(str(v)) for v in df[col].head(200))) + 2
        ws.column_dimensions[get_column_letter(j)].width = min(max(width, 9), 45)
        is_p = str(col).startswith("p ") or str(col) == "p"
        for i in range(2, len(df) + 2):
            c = ws.cell(i, j)
            c.border = Border(top=thin, bottom=thin, left=thin, right=thin)
            if is_p:
                c.number_format = "0.00E+00"
            elif isinstance(c.value, float):
                c.number_format = "0.000"
        if is_p or col == "Significant (BH<0.05)":
            pass
    ws.row_dimensions[1].height = 45
    # highlight significance on the BH column and the Yes/No column
    last = len(df) + 1
    if "p BH-adjusted" in df.columns:
        L = get_column_letter(list(df.columns).index("p BH-adjusted") + 1)
        rng = f"{L}2:{L}{last}"
        ws.conditional_formatting.add(rng, CellIsRule(operator="lessThan", formula=[str(ALPHA)], fill=sig_fill,
                                                      font=Font(bold=True, color="0B6E2E")))
        ws.conditional_formatting.add(rng, CellIsRule(operator="greaterThanOrEqual", formula=[str(ALPHA)], fill=ns_fill))
    if "Significant (BH<0.05)" in df.columns:
        L = get_column_letter(list(df.columns).index("Significant (BH<0.05)") + 1)
        ws.conditional_formatting.add(f"{L}2:{L}{last}", FormulaRule(formula=[f'{L}2="Yes"'], fill=sig_fill))
    if name.startswith(("6_", "7_")):   # sequential shading for mismatch %
        from openpyxl.formatting.rule import ColorScaleRule
        ws.conditional_formatting.add(
            f"B2:{get_column_letter(len(df.columns))}{last}",
            ColorScaleRule(start_type="num", start_value=0, start_color="FFFFFF",
                           mid_type="num", mid_value=5, mid_color="9EC5F4",
                           end_type="num", end_value=100, end_color="184F95"))
    ws.auto_filter.ref = ws.dimensions
wb.save(OUT_XLSX)

# =================================================================== figures
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9, "axes.spines.top": False,
    "axes.spines.right": False, "axes.edgecolor": "#8a8986", "axes.labelcolor": "#0b0b0b",
    "xtick.color": "#52514e", "ytick.color": "#52514e", "axes.titleweight": "bold",
    "axes.titlesize": 10, "savefig.dpi": 300, "savefig.bbox": "tight",
})
BLUE_SEQ = ["#ffffff", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
seq_cmap = LinearSegmentedColormap.from_list("seq", BLUE_SEQ)
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#6250d6"]  # validated categorical slots
INK, INK2 = "#0b0b0b", "#52514e"


def save(fig, stem):
    for ext in ("png", "pdf"):
        fig.savefig(f"{FIG}/{stem}.{ext}")
    plt.close(fig)


# ---- Fig 1: mismatch % heatmap (assay × variant) + BH-adjusted omnibus p
fig, (ax, axp) = plt.subplots(1, 2, figsize=(7.2, 5.2), gridspec_kw=dict(width_ratios=[4, 1.35], wspace=0.05))
data = pct_v.to_numpy()
norm = LogNorm(vmin=0.01, vmax=100)
im = ax.imshow(np.clip(data, 0.01, 100), cmap=seq_cmap, norm=norm, aspect="auto")
for i in range(data.shape[0]):
    for j in range(data.shape[1]):
        v = data[i, j]
        ax.text(j, i, f"{v:.2f}" if v < 10 else f"{v:.1f}", ha="center", va="center", fontsize=7.5,
                color="white" if v >= 10 else INK)
n_lbl = [f"{v}\n(n={tot_v.loc[v].iloc[0]:,})" for v in VARIANTS]
ax.set_xticks(range(4), n_lbl); ax.set_yticks(range(len(NAMES)), NAMES)
ax.xaxis.tick_top(); ax.tick_params(length=0)
for s in ax.spines.values(): s.set_visible(False)
ax.set_ylabel("Primer–probe configuration")
axp.set_xlim(0, 1); axp.set_ylim(len(NAMES) - 0.5, -0.5); axp.axis("off")
axp.text(0.5, -0.85, "Omnibus test\nBH-adjusted p", ha="center", va="bottom", fontsize=8, color=INK2)
for i, r in T_var.iterrows():
    txt = fmt_p(r["p BH-adjusted"]) + f"  {r['Significance']}"
    axp.text(0.5, i, txt, ha="center", va="center", fontsize=7.5,
             color=INK if r["Significance"] != "ns" else "#8a8986",
             fontweight="bold" if r["Significance"] != "ns" else "normal")
cb = fig.colorbar(im, ax=[ax, axp], orientation="horizontal", fraction=0.04, pad=0.03, aspect=40)
cb.set_label("Sequences with mismatch (%, log scale)")
cb.set_ticks([0.01, 0.1, 1, 10, 100], labels=["≤0.01", "0.1", "1", "10", "100"])
fig.suptitle("Mismatch frequency by SARS-CoV-2 variant", x=0.42, y=1.02, fontweight="bold")
save(fig, "Fig1_variant_mismatch_heatmap")

# ---- Fig 2: mismatch % over time (small multiples, all assays) with gaps + significant jumps
qlab = list(tot_q.index)
qnum = np.array([int(s[:4]) + (int(s[-1]) - 1) / 4 for s in qlab])
fig, axes = plt.subplots(3, 4, figsize=(10, 6.6), sharex=True)
for k, (col, ax) in enumerate(zip(COLS, axes.flat)):
    y = pct_q.loc[NAME_MAP[col]].to_numpy()
    for miss in (2020.25, 2021.75):  # 2020.Q2, 2021.Q4 not sampled
        ax.axvspan(miss - 0.125, miss + 0.125, color="#ececea", lw=0, zorder=0)
    big = y.max() >= 5
    ax.plot(qnum, y, color=SERIES[0], lw=2, marker="o", ms=4, zorder=3)
    sub = T_qcons[(T_qcons["Assay"] == NAME_MAP[col]) & (T_qcons["Significant (BH<0.05)"] == "Yes")]
    for _, r in sub.iterrows():
        x1 = int(r["To quarter"][:4]) + (int(r["To quarter"][-1]) - 1) / 4
        ax.plot(x1, r["Mismatch % group 2"], "o", ms=8, mfc="none", mec=SERIES[1], mew=1.6, zorder=4)
    tr = T_trend[T_trend["Assay"] == NAME_MAP[col]].iloc[0]
    ax.set_title(f"{NAME_MAP[col]}", loc="left")
    pt = fmt_p(tr['p BH-adjusted']); pt = pt if pt.startswith("<") else "= " + pt
    ax.set_title(f"trend {'↑' if tr['Z'] > 0 else '↓'}  BH p {pt} {tr['Significance']}",
                 loc="right", fontsize=6.5, color=INK2, fontweight="normal")
    ax.set_ylim(-3 if big else -0.05 * max(y.max(), 0.5), 105 if big else max(y.max(), 0.5) * 1.35)
    ax.grid(axis="y", color="#ececea", lw=0.6)
    if k % 4 == 0: ax.set_ylabel("Mismatch (%)")
for ax in axes[-1]:
    ax.set_xticks([2020, 2021, 2022, 2023]); ax.set_xlim(2019.6, 2023.25)
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
fig.legend(handles=[Line2D([], [], color=SERIES[0], lw=2, marker="o", ms=4, label="Mismatch % per quarter"),
                    Line2D([], [], ls="none", marker="o", ms=8, mfc="none", mec=SERIES[1], mew=1.6,
                           label="Significant change vs previous sampled quarter (Fisher, BH<0.05)"),
                    Patch(color="#ececea", label="Quarter not sampled (2020.Q2, 2021.Q4)")],
           loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, -0.04), fontsize=8)
fig.suptitle("Mismatch frequency over time by primer–probe configuration (note: independent y-scales)",
             fontweight="bold", y=1.0)
fig.tight_layout(rect=(0, 0.03, 1, 0.98))
save(fig, "Fig2_mismatch_over_time")

# ---- Fig 3: pairwise variant comparisons, -log10(BH p), annotated with OR
pairs = [f"{a} vs {b}" for a, b in zip(T_vpair["Group 1"], T_vpair["Group 2"])]
T_vpair["Pair"] = pairs
order_pairs = list(dict.fromkeys(pairs))
M = T_vpair.pivot(index="Assay", columns="Pair", values="p BH-adjusted").loc[NAMES, order_pairs]
OR = T_vpair.pivot(index="Assay", columns="Pair", values="Odds ratio (1 vs 2)").loc[NAMES, order_pairs]
bounds = [0, 1.30103, 2, 3, 10, 50, 301]  # -log10 p: ns | <0.05 | <0.01 | <0.001 | <1e-10 | <1e-50 | ...
cols = ["#f2f2f2", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#1c5cab"]
cmap3, norm3 = ListedColormap(cols), BoundaryNorm(bounds, len(cols))
L = -np.log10(M.to_numpy())
fig, ax = plt.subplots(figsize=(7.4, 5.2))
im = ax.imshow(L, cmap=cmap3, norm=norm3, aspect="auto")
for i in range(L.shape[0]):
    for j in range(L.shape[1]):
        p = M.iat[i, j]; o = OR.iat[i, j]
        lab = "ns" if p >= ALPHA else (f"OR {o:.2g}" if 0.01 <= o < 1000 else f"OR {o:.0e}")
        ax.text(j, i, lab, ha="center", va="center", fontsize=6.8,
                color="white" if L[i, j] >= 10 else (INK if p < ALPHA else "#8a8986"))
ax.set_xticks(range(len(order_pairs)), order_pairs, rotation=30, ha="left")
ax.set_yticks(range(len(NAMES)), NAMES); ax.xaxis.tick_top(); ax.tick_params(length=0)
ax.set_xticks(np.arange(-.5, len(order_pairs)), minor=True); ax.set_yticks(np.arange(-.5, len(NAMES)), minor=True)
ax.grid(which="minor", color="white", lw=2); ax.tick_params(which="minor", length=0)
for s in ax.spines.values(): s.set_visible(False)
cb = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02, ticks=[(a + b) / 2 for a, b in zip(bounds[:-1], bounds[1:])])
cb.ax.set_yticklabels(["ns (≥0.05)", "<0.05", "<0.01", "<0.001", "<1e-10", "<1e-50"])
cb.set_label("BH-adjusted p (pairwise Fisher's exact test)")
ax.set_ylabel("Primer–probe configuration")
ax.set_title("Pairwise variant comparisons of mismatch frequency\n(cell text = odds ratio of mismatch, group 1 vs group 2)",
             pad=12, fontsize=9.5)
save(fig, "Fig3_pairwise_variant_BH")

# ---- Fig 4: effect size vs significance (separates statistical from practical significance)
fig, axes = plt.subplots(1, 2, figsize=(9, 3.9), sharey=True)
for ax, T, ttl in [(axes[0], T_var, "Across variants"), (axes[1], T_qtr, "Across quarters")]:
    T = T.set_index("Assay").loc[NAMES[::-1]]
    y = np.arange(len(T))
    sig = T["p BH-adjusted"] < ALPHA
    ax.hlines(y, 0, T["Cramér's V"], color="#cfcfcb", lw=2)
    ax.plot(T["Cramér's V"][sig], y[sig], "o", ms=8, color=SERIES[0], mec="white", mew=1.5, label="BH p < 0.05")
    ax.plot(T["Cramér's V"][~sig], y[~sig], "o", ms=8, mfc="white", mec="#8a8986", mew=1.5, label="BH p ≥ 0.05 (ns)")
    for yy, (v, s) in zip(y, zip(T["Cramér's V"], T["Significance"])):
        ax.text(v + 0.02, yy, f"{v:.3f} {s}", va="center", fontsize=7, color=INK2)
    for x, lab in [(0.1, "small"), (0.3, "medium"), (0.5, "large")]:
        ax.axvline(x, color="#ececea", lw=1, zorder=0)
        ax.text(x, len(T) - 0.3, lab, ha="center", fontsize=6.5, color="#8a8986")
    ax.set_yticks(y, T.index); ax.set_xlim(0, 1.05); ax.set_title(ttl, loc="left", pad=16)
    ax.set_xlabel("Effect size (Cramér's V)")
axes[0].legend(loc="upper right", bbox_to_anchor=(1.0, 0.9), frameon=False, fontsize=7.5)
fig.suptitle("Statistical vs practical significance of omnibus tests", fontweight="bold")
fig.tight_layout()
save(fig, "Fig4_effect_size_vs_significance")

print("wrote", OUT_XLSX, "and figures/:", sorted(os.listdir(FIG)))
