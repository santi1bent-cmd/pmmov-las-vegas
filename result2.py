"""Result 2: does dividing by PMMoV change how well wastewater tracks clinical data?

dR = Spearman(divided, clinical) - Spearman(raw, clinical), zero lag, on weekly series
smoothed with a centred five-week moving average.

Wastewater weekly values: mean of log10 concentrations within each MMWR week (a geometric
mean), then the five-week average of those weekly log values. The raw and divided series
are averaged the same way. Non-detects are set to the detection limit / sqrt(2).
Clinical series: NHSN weekly COVID-19 admissions (smoothed, then weeks with < 75% of the usual
47 Nevada hospitals reporting are removed) and NSSP weekly percent of ED visits, both statewide.
"""
import json
import numpy as np, pandas as pd
from scipy.stats import spearmanr
from common import (RAW, OUT, LABS, SEED, N_BOOT, HOSPITALS_USUAL, MIN_REPORTING, NHSN_COLS,
                    NSSP_NAMES, load_ww, mmwr_week_end, stationary_boot_idx)


def smooth5(s):
    """Centred five-week moving average over the weeks that have data (edges use what is available)."""
    s = s.dropna().sort_index()
    return s.asfreq("7D").rolling(5, center=True, min_periods=1).mean().reindex(s.index)

MEAN_BLOCK_WEEKS = 10


def clinical(pathogen, filter_reporting=True):
    n = pd.DataFrame(json.load(open(RAW / "nhsn.json")))
    n["wk"] = pd.to_datetime(n["weekendingdate"])
    n = n.set_index("wk").sort_index()
    adm_col, rep_col = NHSN_COLS[pathogen]
    adm = pd.to_numeric(n[adm_col], errors="coerce")
    ok = pd.to_numeric(n[rep_col], errors="coerce") >= MIN_REPORTING * HOSPITALS_USUAL
    if not filter_reporting:
        ok = adm.notna()
    e = pd.DataFrame(json.load(open(RAW / "nssp.json")))
    e = e[e["pathogen"] == NSSP_NAMES[pathogen]].copy()
    e["wk"] = pd.to_datetime(e["week_end"])
    ed = pd.to_numeric(e.set_index("wk").sort_index()["percent_visits"], errors="coerce")
    # Admissions are smoothed on the full series, then weeks with incomplete reporting are removed.
    A = smooth5(adm)
    return A[ok.reindex(A.index).fillna(False).astype(bool)], smooth5(ed), ok


def wastewater(pathogen, lab, nondetect="lod_sqrt2"):
    d = load_ww(pathogen)
    d = d[d["source"] == lab].copy()
    c = d["pcr_target_avg_conc_lin"].copy()
    nd = ~d["detected"]
    c[nd] = d.loc[nd, "lod_sewage"] / (np.sqrt(2) if nondetect == "lod_sqrt2" else 1.0)
    d["raw"] = c
    d["divided"] = c / d["hum_frac_mic_conc"]
    d = d[(d["raw"] > 0) & (d["divided"] > 0)]
    d["wk"] = mmwr_week_end(d["date"])
    return pd.DataFrame({k: smooth5(d.groupby("wk")[k].apply(lambda v: np.mean(np.log10(v))))
                         for k in ["raw", "divided"]})


def delta_r(j, col):
    return spearmanr(j["divided"], j[col])[0] - spearmanr(j["raw"], j[col])[0]


def with_ci(j, col, rng=None):
    rng = np.random.default_rng(SEED)        # fresh, fixed seed for every interval
    boots = [delta_r(j.iloc[stationary_boot_idx(len(j), MEAN_BLOCK_WEEKS, rng)], col) for _ in range(N_BOOT)]
    lo, hi = np.nanpercentile(boots, [2.5, 97.5])   # a resample can be constant; skip those
    return {"weeks": len(j), "first_week": str(j.index.min().date()), "last_week": str(j.index.max().date()),
            "raw_R": spearmanr(j["raw"], j[col])[0], "dR": delta_r(j, col), "ci": [lo, hi]}


def windows(pathogen, lab, rng, nondetect="lod_sqrt2", filter_reporting=True):
    A, E, ok = clinical(pathogen, filter_reporting)
    W = wastewater(pathogen, lab, nondetect)
    keep = lambda j: j[ok.reindex(j.index).fillna(False).astype(bool)]
    full = keep(pd.concat([W, A.rename("adm")], axis=1, sort=True).dropna())
    same = keep(pd.concat([W, A.rename("adm"), E.rename("ed")], axis=1, sort=True).dropna())
    ed_all = pd.concat([W, E.rename("ed")], axis=1, sort=True).dropna()
    out = {}
    for name, j, col in [("admissions, full window", full, "adm"), ("admissions, same weeks as ED", same, "adm"),
                         ("ED share, same weeks", same, "ed"), ("ED share, every ED week", ed_all, "ed")]:
        if len(j) >= 10:
            out[name] = with_ci(j, col, rng)
    return out


if __name__ == "__main__":
    res = {}
    for pathogen in ["SARS-CoV-2", "Influenza A", "RSV"]:
        for lab in LABS:
            res[f"{pathogen} | {lab}"] = windows(pathogen, lab, np.random.default_rng(SEED))
    json.dump(res, open(OUT / "result2.json", "w"), indent=2, default=float)
    for k, v in res.items():
        print(f"\n{k}")
        for w, s in v.items():
            print(f"  {w:30s} {s['weeks']:4d} wk  raw R {s['raw_R']:.3f}  dR {s['dR']:+.3f}  "
                  f"[{s['ci'][0]:+.3f}, {s['ci'][1]:+.3f}]  {s['first_week']} to {s['last_week']}")
