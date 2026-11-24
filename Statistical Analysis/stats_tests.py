"""
Inferential tests for primer/probe mismatch frequencies (reviewer comments 1-3).

Input : result.xlsx  (sheets Lineage / Lineage_Match, Date_Count / Date_match_Count)
Output: stats_results.xlsx

For every primer-probe configuration (columns A..R):
  * Variants  : omnibus Pearson chi-square + Fisher-Freeman-Halton exact test
                (Monte-Carlo, 200k tables) + Cramer's V;
                pairwise 2x2 Fisher exact tests with odds ratios.
  * Quarters  : omnibus chi-square + Fisher-Freeman-Halton exact,
                Cochran-Armitage test for linear trend,
                2x2 Fisher exact tests between consecutive sampled quarters.
  * Multiple testing: Benjamini-Hochberg FDR (primary) and Bonferroni
                (sensitivity), applied within each family of tests.
"""
import numpy as np
import pandas as pd
from scipy import stats

XLSX = "result.xlsx"
OUT = "stats_results.xlsx"
N_SIM = 200_000
rng = np.random.default_rng(20261120)

# Fill in the assay name for each column letter (e.g. "B": "China CDC-N ...").
NAME_MAP = {c: c for c in list("ABCDEFGHIPQR")}

MAIN_VARIANTS = ["Other", "Alpha", "Delta", "Omicron"]
MINOR_VARIANTS = ["Beta", "Gamma", "Epsilon", "Zeta", "Lota", "Kappa", "Eta", "Mu"]


def read_block(sheet):
    df = pd.read_excel(XLSX, sheet, header=0, index_col=1)
    df = df.drop(columns=df.columns[0]).dropna(how="all")
    df.index = df.index.astype(str).str.strip()
    return df.astype(int)


# ---------------------------------------------------------------- statistics
def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(ranked, 1)
    return out


def bonf(p):
    return np.minimum(np.asarray(p, float) * len(p), 1)


def log_table_prob(mm, totals, M):
    """log P(mismatch counts mm | row totals, total mismatches M) – multivariate hypergeometric."""
    from scipy.special import gammaln
    lc = lambda n, k: gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1)
    return lc(totals, mm).sum(axis=-1) - lc(totals.sum(), M)


def ffh_exact_mc(mm, totals):
    """Fisher-Freeman-Halton exact test for an r x 2 table, Monte-Carlo p-value."""
    M = int(mm.sum())
    if M == 0 or M == totals.sum():
        return np.nan
    obs = log_table_prob(mm, totals, M)
    sims = stats.multivariate_hypergeom.rvs(m=totals, n=M, size=N_SIM, random_state=rng)
    lp = log_table_prob(sims, totals, M)
    k = np.sum(lp <= obs + 1e-7)
    return (k + 1) / (N_SIM + 1)


def omnibus(groups, mm, totals):
    table = np.column_stack([mm, totals - mm])
    keep = totals > 0
    table, mm_k, tot_k = table[keep], mm[keep], totals[keep]
    if mm_k.sum() == 0:
        return dict(chi2=np.nan, df=len(tot_k) - 1, p_chi2=np.nan, cramers_v=np.nan,
                    min_expected=np.nan, pct_cells_exp_lt5=np.nan, p_fisher_exact=np.nan,
                    note="no mismatches in any group – test not applicable")
    chi2, p, dof, exp = stats.chi2_contingency(table, correction=False)
    v = np.sqrt(chi2 / (table.sum() * (min(table.shape) - 1)))
    p_ex = ffh_exact_mc(mm_k, tot_k)
    lt5 = (exp < 5).mean() * 100
    return dict(chi2=chi2, df=dof, p_chi2=p, cramers_v=v, min_expected=exp.min(),
                pct_cells_exp_lt5=lt5, p_fisher_exact=p_ex,
                note="use exact p (expected <5 in >20% of cells)" if lt5 > 20 else "")


def cochran_armitage(mm, totals, scores):
    x, n, s = map(lambda a: np.asarray(a, float), (mm, totals, scores))
    N, X = n.sum(), x.sum()
    pbar = X / N
    if X == 0 or X == N:
        return np.nan, np.nan
    T = np.sum(s * (x - n * pbar))
    var = pbar * (1 - pbar) * (np.sum(n * s ** 2) - np.sum(n * s) ** 2 / N)
    z = T / np.sqrt(var)
    return z, 2 * stats.norm.sf(abs(z))


def fisher2x2(m1, n1, m2, n2):
    t = [[m1, n1 - m1], [m2, n2 - m2]]
    or_, p = stats.fisher_exact(t)
    # Haldane-Anscombe corrected OR + Woolf 95% CI (robust to zero cells)
    a, b, c, d = m1 + .5, n1 - m1 + .5, m2 + .5, n2 - m2 + .5
    lor = np.log(a * d / (b * c)); se = np.sqrt(1/a + 1/b + 1/c + 1/d)
    return p, np.exp(lor), np.exp(lor - 1.96 * se), np.exp(lor + 1.96 * se)


# ---------------------------------------------------------------- analyses
def group_analysis(total_df, match_df, groups, label):
    omni, pair = [], []
    for col in total_df.columns:
        tot = total_df.loc[groups, col].to_numpy()
        mm = tot - match_df.loc[groups, col].to_numpy()
        r = omnibus(groups, mm, tot)
        r.update(assay=NAME_MAP[col], column=col, comparison=label,
                 **{f"mm_{g}": f"{m}/{t} ({100*m/t:.2f}%)" for g, m, t in zip(groups, mm, tot)})
        omni.append(r)
        for i in range(len(groups)):
            for j in range(i + 1, len(groups)):
                p, OR, lo, hi = fisher2x2(mm[i], tot[i], mm[j], tot[j])
                pair.append(dict(assay=NAME_MAP[col], column=col, group1=groups[i], group2=groups[j],
                                 mm_pct_1=100 * mm[i] / tot[i], mm_pct_2=100 * mm[j] / tot[j],
                                 odds_ratio=OR, OR_95CI_low=lo, OR_95CI_high=hi, p_fisher=p))
    omni, pair = pd.DataFrame(omni), pd.DataFrame(pair)
    for df, pcol in [(omni, "p_fisher_exact"), (omni, "p_chi2"), (pair, "p_fisher")]:
        ok = df[pcol].notna()
        df.loc[ok, pcol + "_BH"] = bh(df.loc[ok, pcol])
        df.loc[ok, pcol + "_Bonf"] = bonf(df.loc[ok, pcol])
    lead = ["assay", "column", "comparison"]
    omni = omni[lead + [c for c in omni.columns if c not in lead]]
    return omni, pair


def main():
    lin_t, lin_m = read_block("Lineage"), read_block("Lineage_Match")
    dat_t, dat_m = read_block("Date_Count"), read_block("Date_match_Count")
    assert (lin_t.sum() == dat_t.sum()).all()

    # pooled minor VOC/VOI group for the sensitivity analysis
    for d in (lin_t, lin_m):
        d.loc["Other VOC/VOI"] = d.loc[MINOR_VARIANTS].sum()

    minor_counts = lin_t.loc[MINOR_VARIANTS].sum()

    print("Minor lineage counts per assay:")
    print(minor_counts)

    print("\nTotal minor lineage sequences:")
    print(minor_counts.iloc[0])

    v_omni, v_pair = group_analysis(lin_t, lin_m, MAIN_VARIANTS, "Variants (Other, Alpha, Delta, Omicron)")
    v5_omni, _ = group_analysis(lin_t, lin_m, MAIN_VARIANTS + ["Other VOC/VOI"],
                                "Sensitivity: + pooled minor VOC/VOI")

    # chronological ordering of quarters
    q = sorted(dat_t.index, key=lambda s: (int(s[:4]), int(s[-1])))
    dat_t, dat_m = dat_t.loc[q], dat_m.loc[q]
    t_omni, t_pair_all = group_analysis(dat_t, dat_m, q, "Quarters")

    # consecutive sampled quarters + Cochran-Armitage trend
    qi = np.array([int(s[:4]) * 4 + int(s[-1]) for s in q]); qi -= qi.min()  # gaps respected
    cons, trend = [], []
    for col in dat_t.columns:
        tot = dat_t[col].to_numpy(); mm = tot - dat_m[col].to_numpy()
        z, p = cochran_armitage(mm, tot, qi)
        trend.append(dict(assay=NAME_MAP[col], column=col, z=z, p_trend=p))
        for i in range(len(q) - 1):
            p2, OR, lo, hi = fisher2x2(mm[i], tot[i], mm[i + 1], tot[i + 1])
            cons.append(dict(assay=NAME_MAP[col], column=col, from_q=q[i], to_q=q[i + 1],
                             gap_quarters=qi[i + 1] - qi[i],
                             mm_pct_from=100 * mm[i] / tot[i], mm_pct_to=100 * mm[i + 1] / tot[i + 1],
                             odds_ratio=OR, OR_95CI_low=lo, OR_95CI_high=hi, p_fisher=p2))
    cons, trend = pd.DataFrame(cons), pd.DataFrame(trend)
    for df, pcol in [(cons, "p_fisher"), (trend, "p_trend")]:
        ok = df[pcol].notna()
        df.loc[ok, pcol + "_BH"] = bh(df.loc[ok, pcol]); df.loc[ok, pcol + "_Bonf"] = bonf(df.loc[ok, pcol])

    # contingency tables for transparency
    mm_lin = (lin_t - lin_m).loc[MAIN_VARIANTS + MINOR_VARIANTS + ["Other VOC/VOI"]]
    mm_dat = dat_t - dat_m

    with pd.ExcelWriter(OUT) as w:
        v_omni.to_excel(w, sheet_name="Variant_omnibus", index=False)
        v5_omni.to_excel(w, sheet_name="Variant_omnibus_sens", index=False)
        v_pair.to_excel(w, sheet_name="Variant_pairwise_Fisher", index=False)
        t_omni.to_excel(w, sheet_name="Quarter_omnibus", index=False)
        trend.to_excel(w, sheet_name="Quarter_trend_CA", index=False)
        cons.to_excel(w, sheet_name="Quarter_consecutive_Fisher", index=False)
        mm_lin.to_excel(w, sheet_name="Mismatch_counts_variant")
        lin_t.to_excel(w, sheet_name="Totals_variant")
        mm_dat.to_excel(w, sheet_name="Mismatch_counts_quarter")
        dat_t.to_excel(w, sheet_name="Totals_quarter")

    pd.set_option("display.width", 250, "display.max_columns", 30, "display.max_colwidth", 22)
    print(v_omni.drop(columns=["comparison"]).to_string(float_format=lambda x: f"{x:.3g}"))
    print(t_omni[["column", "chi2", "df", "p_chi2", "cramers_v", "pct_cells_exp_lt5",
                  "p_fisher_exact", "p_fisher_exact_BH"]].to_string(float_format=lambda x: f"{x:.3g}"))
    print(trend.to_string(float_format=lambda x: f"{x:.3g}"))
    print(v_pair[v_pair.column.isin(["B", "F", "Q"])].to_string(float_format=lambda x: f"{x:.3g}"))
    print(cons[cons.column.isin(["B", "F", "Q"])].to_string(float_format=lambda x: f"{x:.3g}"))


if __name__ == "__main__":
    main()
