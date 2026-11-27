"""Shared loading and helper functions."""
import json, pathlib
import numpy as np, pandas as pd

ROOT = pathlib.Path(__file__).parent
RAW = ROOT / "data" / "raw"
OUT = ROOT / "output"
OUT.mkdir(exist_ok=True)

LABS = ["WastewaterSCAN", "CDC_Verily"]
SEED = 20261126
N_BOOT = 2000
HOSPITALS_USUAL = 47          # Nevada hospitals normally reporting to NHSN
MIN_REPORTING = 0.75          # keep NHSN weeks with >= 75% of the usual hospitals

WW_FILES = {"SARS-CoV-2": "sars", "Influenza A": "flua", "RSV": "rsv"}
NHSN_COLS = {"SARS-CoV-2": ("totalconfc19newadm", "totalconfc19newadmhosprep"),
             "Influenza A": ("totalconfflunewadm", "totalconfflunewadmhosprep"),
             "RSV": ("totalconfrsvnewadm", "totalconfrsvnewadmhosprep")}
NSSP_NAMES = {"SARS-CoV-2": "COVID-19", "Influenza A": "Influenza", "RSV": "RSV"}


def load_ww(pathogen):
    """One row per lab per sample day, WastewaterSCAN and CDC_Verily only (CDC_Biobot excluded)."""
    d = pd.DataFrame(json.load(open(RAW / f"{WW_FILES[pathogen]}.json")))
    d = d[d["source"].isin(LABS)].copy()
    d["date"] = pd.to_datetime(d["sample_collect_date"])
    for c in ["pcr_target_avg_conc_lin", "hum_frac_mic_conc", "pcr_target_mic_lin", "lod_sewage"]:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    d["detected"] = d["pcr_target_detect"].str.lower().eq("yes")
    return d


def mmwr_week_end(dates):
    """Saturday ending the MMWR week that contains each date."""
    return dates + pd.to_timedelta((5 - dates.dt.dayofweek) % 7, unit="D")


def smooth5(s):
    """Centred five-week moving average on a regular weekly index."""
    return s.asfreq("7D").rolling(5, center=True, min_periods=3).mean()


def stationary_boot_idx(n, mean_block, rng):
    """Politis-Romano stationary bootstrap: blocks of geometric length, wrapping around."""
    p = 1.0 / mean_block
    idx = np.empty(n, dtype=int)
    i = rng.integers(n)
    for k in range(n):
        idx[k] = i
        i = rng.integers(n) if rng.random() < p else (i + 1) % n
    return idx
