"""Result 1: do the two labs agree better after dividing by PMMoV?

Same-day pairs at site 1179 where both labs detected the target.
gap = log10(WastewaterSCAN) - log10(CDC_Verily), raw and divided by each lab's own PMMoV.
The labs report in different units, so the mean gap is a unit conversion; only the spread is interpreted.
"""
import json
import numpy as np, pandas as pd
from scipy import stats
from common import load_ww, LABS, SEED, N_BOOT, OUT, stationary_boot_idx

MEAN_BLOCK_DAYS = 20   # mean block = 20 paired sampling days (about 10 weeks)


def paired(pathogen, normalized_field="computed"):
    d = load_ww(pathogen)
    d = d[d["detected"] & (d["pcr_target_avg_conc_lin"] > 0) & (d["hum_frac_mic_conc"] > 0)]
    if normalized_field == "cdc_rounded":          # sensitivity: CDC's rounded normalized field
        d = d[d["pcr_target_mic_lin"] > 0]
    d = d.assign(lraw=np.log10(d["pcr_target_avg_conc_lin"]), lpm=np.log10(d["hum_frac_mic_conc"]))
    p = d.pivot_table(index="date", columns="source", values=["lraw", "lpm"], aggfunc="mean").dropna()
    w, v = LABS
    out = pd.DataFrame({"raw": p[("lraw", w)] - p[("lraw", v)],
                        "g": p[("lpm", w)] - p[("lpm", v)]}).sort_index()
    out["divided"] = out["raw"] - out["g"]           # dividing becomes subtracting on the log scale
    return out


def summarize(p, rng):
    raw, g, dv = p["raw"].to_numpy(), p["g"].to_numpy(), p["divided"].to_numpy()
    sd_raw, sd_div = raw.std(ddof=1), dv.std(ddof=1)
    var_raw, var_g, cov = raw.var(ddof=1), g.var(ddof=1), np.cov(raw, g, ddof=1)[0, 1]
    r, pval = stats.pearsonr(raw, g)
    boots = []
    for _ in range(N_BOOT):
        i = stationary_boot_idx(len(raw), MEAN_BLOCK_DAYS, rng)
        boots.append(dv[i].std(ddof=1) / raw[i].std(ddof=1))
    lo, hi = np.percentile(boots, [2.5, 97.5])
    centred = raw - np.median(raw)
    return {
        "paired_days": len(p),
        "sd_raw": sd_raw, "sd_divided": sd_div,
        "sd_ratio": sd_div / sd_raw, "sd_ratio_ci": [lo, hi],
        "loa_fold_raw": 10 ** (2 * 1.96 * sd_raw), "loa_fold_divided": 10 ** (2 * 1.96 * sd_div),
        "median_fold_disagreement": 10 ** np.median(np.abs(centred)),
        "pearson_r_gap_vs_pmmov_gap": r, "pearson_p": pval,
        "shared_part_cancelled": cov / var_raw,             # Cov(raw, g) / Var(raw)
        "pmmov_noise_added": (var_g - cov) / var_raw,       # [Var(g) - Cov(raw, g)] / Var(raw)
        "best_possible_removal_r2": r ** 2,
    }


if __name__ == "__main__":
    res = {}
    for pathogen in ["SARS-CoV-2", "Influenza A"]:
        rng = np.random.default_rng(SEED)
        res[pathogen] = summarize(paired(pathogen), rng)
        rng = np.random.default_rng(SEED)
        res[pathogen + " (CDC rounded field)"] = summarize(paired(pathogen, "cdc_rounded"), rng)
    rsv = load_ww("RSV")
    both = rsv.pivot_table(index="date", columns="source", values="detected", aggfunc="max").dropna()
    res["RSV shared days / both detected"] = [len(both), int(both.all(axis=1).sum())]
    for pathogen in ["SARS-CoV-2", "Influenza A"]:
        d = load_ww(pathogen)
        m = d.pivot_table(index="date", columns="source", values="detected", aggfunc="max").dropna()
        res[pathogen + " shared days / both detected"] = [len(m), int(m.all(axis=1).sum())]
    json.dump(res, open(OUT / "result1.json", "w"), indent=2, default=float)
    for k, v in res.items():
        if isinstance(v, dict):
            print(f"\n{k}")
            for kk, vv in v.items():
                print(f"  {kk:28s} {np.round(vv, 4) if not isinstance(vv, int) else vv}")
        else:
            print(f"\n{k}: {v}")
