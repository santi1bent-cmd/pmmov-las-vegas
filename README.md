# Does dividing by PMMoV fix wastewater data? A check in Las Vegas

Code and data for the write-up *Does dividing by PMMoV fix wastewater data? A check in my hometown, Las Vegas* (Santi Bent, November 2026).

CDC wastewater site 1179, the Flamingo Water Resource Center in Las Vegas (about 1.25 million people served), was measured on the same days from November 2023 to April 2026 by two separate pipelines: WastewaterSCAN (solids, copies per gram of dry sludge) and CDC/Verily (liquid, copies per litre). This repository asks two questions:

1. **Do the two labs agree better after dividing by PMMoV?** For SARS-CoV-2, no: the spread of the between-lab gap changes by a ratio of 1.02 (95% interval 0.95 to 1.12). Dividing cancels a shared part worth 32% of the disagreement but adds PMMoV's own noise worth 36%.
2. **Does dividing change how well wastewater tracks clinical data?** Against Nevada hospital admissions, no (ΔR −0.001, interval −0.035 to +0.026, 164 weeks). The comparator "sign flip" that Justen et al. (2026) report for four Chicago sites does not appear here.

Scope: one sewershed, PCR rather than sequencing, PMMoV rather than the tobamovirus marker Justen et al. found best.

## Run it

```bash
pip install -r requirements.txt
python fetch.py        # optional: re-download the CDC data into data/raw/
python result1.py      # Test 1  -> output/result1.json
python result2.py      # Test 2  -> output/result2.json
python sensitivity.py  # checks in Sections 4.3-4.4 -> output/sensitivity.json
python figures.py      # Figures 3-7 -> output/figures/
```

`data/raw/` holds the snapshot the write-up uses (fetched 27 November 2026; see `data/raw/FETCH_DATE`). Run from that snapshot, the scripts reproduce every Las Vegas number and Figures 3–7 in the write-up. The Chicago and WastewaterSCAN-site values come from the public tables in Justen et al.'s repository, and the Illinois and Massachusetts reporting counts from the same NHSN dataset for those states. Running `fetch.py` again downloads the current data, which can differ slightly if CDC revises it.

## Files

| File | What it does |
| --- | --- |
| `fetch.py` | Downloads five public CDC datasets from data.cdc.gov: j9g8-acpt (SARS-CoV-2), ymmh-divb (influenza A) and 45cq-cw4i (RSV) for site 1179; ua7e-t2fy (NHSN admissions, Nevada); 7mra-9cq9 (NSSP emergency department visits, Nevada). |
| `common.py` | Loading, MMWR weeks, the stationary block bootstrap, shared settings (seed 20261126, 2,000 resamples). |
| `result1.py` | Test 1. Same-day pairs where both labs detected the target. SD of the log10 between-lab gap before and after dividing, its ratio with a block-bootstrap interval (mean block 20 paired days), limits of agreement, and the split into "shared part cancelled" and "PMMoV's own noise added". Also the check using CDC's rounded normalized field. |
| `result2.py` | Test 2. ΔR = Spearman(divided, clinical) − Spearman(raw, clinical) on weekly series with a centred five-week average, with block-bootstrap intervals (mean block 10 weeks). Windows: admissions over the full window, both clinical series on the same weeks, and ED on every ED week. |
| `sensitivity.py` | Justen et al.'s smoothing kernel and CDC's rounded field (alone and together), leaving in the weeks with incomplete hospital reporting, and non-detects at the detection limit. |
| `figures.py` | Figures 3 to 7 of the write-up. |

## Method notes

- **Normalization.** Raw concentration (`pcr_target_avg_conc_lin`) divided by the same sample's PMMoV (`hum_frac_mic_conc`), computed from unrounded fields. CDC's own normalized field (`pcr_target_mic_lin`) is rounded and turns some low values into zero, so it is used only as a sensitivity check.
- **Labs.** WastewaterSCAN and CDC_Verily only. The earlier CDC_Biobot series (May 2022 to September 2023) used a different method (qPCR) and sample type (raw wastewater) and is excluded.
- **Weekly wastewater values** are the mean of the log10 readings in each MMWR week, followed by a centred five-week average of those weekly values. The raw and divided series are averaged the same way. Non-detects are set to the detection limit divided by √2 (Test 2 only; Test 1 uses days both labs detected the target).
- **Hospital reporting gap.** Hospital reporting of COVID-19 admissions was voluntary from 1 May to 1 November 2024 and Nevada reporting collapsed. Admissions are smoothed on the full series, then weeks with fewer than 75% of the usual 47 Nevada hospitals reporting are removed.
- **Intervals** are percentile intervals from a stationary block bootstrap, each drawn with a fresh generator seeded 20261126. They are not adjusted for multiple comparisons.

## Data and licence

CDC NWSS, NHSN and NSSP data are public and come from data.cdc.gov. Code is under the MIT licence (`LICENSE`).
