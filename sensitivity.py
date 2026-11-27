"""Sensitivity checks reported in Sections 4.3 and 4.4 (SARS-CoV-2 unless stated).

1. Justen et al.'s smoothing kernel (log1p -> weekly mean -> centred 5-week mean -> expm1)
   and CDC's rounded normalized field (pcr_target_mic_lin), alone and together.
2. Leaving in the NHSN weeks with incomplete hospital reporting.
3. Non-detects at the detection limit instead of the limit / sqrt(2).
"""
import json
import numpy as np, pandas as pd
from scipy.stats import spearmanr
from common import OUT, LABS, load_ww, mmwr_week_end
from result2 import clinical, wastewater, delta_r, smooth5


def justen_kernel(per_sample):
    """Their smoothing.py kernel: weekly mean of log1p(x), centred 5-week mean on the week grid, expm1."""
    w = per_sample.groupby(level=0).apply(lambda v: np.mean(np.log1p(v))).sort_index()
    return np.expm1(w.asfreq("7D").rolling(5, center=True, min_periods=1).mean().reindex(w.index))


def ww_variant(lab, kernel, field):
    d = load_ww("SARS-CoV-2")
    d = d[d["source"] == lab].copy()
    c = d["pcr_target_avg_conc_lin"].copy()
    nd = ~d["detected"]
    c[nd] = d.loc[nd, "lod_sewage"] / np.sqrt(2)
    d["raw"] = c
    d["divided"] = d["pcr_target_mic_lin"] if field == "cdc_rounded" else c / d["hum_frac_mic_conc"]
    d["wk"] = mmwr_week_end(d["date"])
    d = d.set_index("wk")
    if kernel == "justen":
        d = d[(d["raw"] >= 0) & (d["divided"] >= 0)]
        return pd.DataFrame({k: justen_kernel(d[k]) for k in ["raw", "divided"]})
    out = {}
    for k in ["raw", "divided"]:
        s = d.loc[d[k] > 0, k]                       # log scale: zeros from rounding cannot be logged
        out[k] = smooth5(np.log10(s).groupby(level=0).mean())
    return pd.DataFrame(out)


def run():
    res = {}
    A, E, ok = clinical("SARS-CoV-2")
    keep = lambda j: j[ok.reindex(j.index).fillna(False).astype(bool)]
    for lab in LABS:
        r = {}
        for name, kernel, field in [("this analysis", "log", "computed"),
                                    ("CDC rounded field only", "log", "cdc_rounded"),
                                    ("Justen smoothing only", "justen", "computed"),
                                    ("Justen smoothing + CDC rounded field", "justen", "cdc_rounded")]:
            W = ww_variant(lab, kernel, field)
            full = keep(pd.concat([W, A.rename("adm")], axis=1, sort=True).dropna())
            same = keep(pd.concat([W, A.rename("adm"), E.rename("ed")], axis=1, sort=True).dropna())
            r[name] = {"dR admissions, full window": delta_r(full, "adm"),
                       "dR ED share, same weeks": delta_r(same, "ed")}
        # NHSN weeks with incomplete reporting left in
        A_all, _, _ = clinical("SARS-CoV-2", filter_reporting=False)
        W = wastewater("SARS-CoV-2", lab)
        jf = keep(pd.concat([W, A.rename("adm")], axis=1, sort=True).dropna())
        ja = pd.concat([W, A_all.rename("adm")], axis=1, sort=True).dropna()
        r["raw R vs admissions, full window: filtered / gap weeks left in"] = [
            spearmanr(jf["raw"], jf["adm"])[0], spearmanr(ja["raw"], ja["adm"])[0]]
        # non-detects at the detection limit
        W1 = wastewater("SARS-CoV-2", lab, nondetect="lod")
        j1 = keep(pd.concat([W1, A.rename("adm"), E.rename("ed")], axis=1, sort=True).dropna())
        j1f = keep(pd.concat([W1, A.rename("adm")], axis=1, sort=True).dropna())
        r["non-detects at LOD: dR adm full / adm same weeks / ED same weeks"] = [
            delta_r(j1f, "adm"), delta_r(j1, "adm"), delta_r(j1, "ed")]
        res[lab] = r
    return res


if __name__ == "__main__":
    res = run()
    json.dump(res, open(OUT / "sensitivity.json", "w"), indent=2, default=float)
    for lab, r in res.items():
        print(f"\n{lab}")
        for k, v in r.items():
            if isinstance(v, dict):
                print(f"  {k:40s} " + "  ".join(f"{kk}: {vv:+.3f}" for kk, vv in v.items()))
            else:
                print(f"  {k}: " + ", ".join(f"{x:+.3f}" for x in v))
