"""Figures 3-7 of the write-up, from the raw data and the outputs of result1.py and result2.py.
(Figure 1 is a photograph and Figure 2 a design diagram, so neither is generated here.)"""
import json
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from common import OUT, RAW, LABS, load_ww, HOSPITALS_USUAL
from result1 import paired

FIG = OUT / "figures"; FIG.mkdir(exist_ok=True)
BLUE, ORANGE, VIOLET, GRAY = "#257de0", "#f47230", "#6a4ce2", "#777777"
plt.rcParams.update({"font.family": "Liberation Serif", "font.size": 11, "axes.edgecolor": "#999",
                     "axes.labelcolor": "#555", "xtick.color": "#555", "ytick.color": "#555",
                     "axes.spines.top": False, "axes.spines.right": False})
R1 = json.load(open(OUT / "result1.json"))
R2 = json.load(open(OUT / "result2.json"))
GAP = (pd.Timestamp("2024-05-01"), pd.Timestamp("2024-11-01"))


def fig3():
    d = load_ww("SARS-CoV-2")
    days = paired("SARS-CoV-2").index                 # the paired days used in Result 1
    d = d[d["detected"] & d["date"].isin(days)]
    start, end = days.min(), days.max()
    fig, ax = plt.subplots(figsize=(5.855, 2.44), dpi=200)
    for lab, col, name in [("WastewaterSCAN", BLUE, "WastewaterSCAN"), ("CDC_Verily", ORANGE, "CDC / Verily")]:
        x = d[(d.source == lab) & d.date.between(start, end)].sort_values("date")
        y = x.pcr_target_avg_conc_lin / x.pcr_target_avg_conc_lin.median()
        ax.plot(x.date, y, "-o", color=col, lw=1.2, ms=2.5)
        ax.text(0.2, 0.93 if lab == LABS[0] else 0.81, name, transform=ax.transAxes, color=col, fontweight="bold", fontsize=10)
    ax.set_yscale("log"); ax.set_ylim(0.05, 100); ax.set_ylabel("SARS-CoV-2 level\n(× each lab’s own median)", fontsize=10)
    ax.xaxis.set_major_locator(mdates.MonthLocator((1, 7))); ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %Y"))
    ax.grid(axis="y", color="#e5e5e5"); fig.tight_layout(); fig.savefig(FIG / "fig3_two_labs.png"); plt.close(fig)


def fig4():
    p = paired("SARS-CoV-2")
    raw, g, dv = (p[c] - p[c].median() for c in ["raw", "g", "divided"])
    s = R1["SARS-CoV-2"]
    fig, (a, b) = plt.subplots(1, 2, figsize=(6.375, 2.775), dpi=200)
    a.scatter(g, raw, s=9, color=VIOLET, alpha=0.6, lw=0)
    k, c0 = np.polyfit(g, raw, 1); xx = np.linspace(g.min(), g.max(), 2); a.plot(xx, k * xx + c0, color="#222", lw=1.2)
    a.axhline(0, color="#ccc", lw=0.8); a.axvline(0, color="#ccc", lw=0.8)
    a.text(0.05, 0.95, f"r = +{s['pearson_r_gap_vs_pmmov_gap']:.2f}\n{s['paired_days']} paired days", transform=a.transAxes, va="top", fontsize=9)
    a.set_xlabel("Gap between labs in PMMoV reading", fontsize=9.5); a.set_ylabel("Gap between labs in\nSARS-CoV-2 reading", fontsize=9.5)
    a.set_title("a  When the labs disagree on PMMoV,\n   they tend to disagree on the virus too", loc="left", fontsize=9.5, fontweight="bold")
    bins = np.arange(-1.25, 1.3, 0.1)
    b.hist(raw, bins=bins, color="#d0d0d0", edgecolor=GRAY, histtype="stepfilled", lw=1.2)
    b.hist(dv, bins=bins, color=VIOLET, histtype="step", lw=1.6)
    b.text(0.05, 0.95, f"Raw: SD {s['sd_raw']:.3f}", transform=b.transAxes, color=GRAY, fontweight="bold", fontsize=9, va="top")
    b.text(0.05, 0.84, f"Divided by PMMoV: SD {s['sd_divided']:.3f}", transform=b.transAxes, color=VIOLET, fontweight="bold", fontsize=9, va="top")
    b.set_xlabel("Gap between labs (log10 units, centred)", fontsize=9.5); b.set_ylabel("Number of days", fontsize=9.5)
    b.set_title("b  Dividing by PMMoV leaves the spread\n   almost unchanged", loc="left", fontsize=9.5, fontweight="bold")
    b.set_ylim(0, 52); b.grid(axis="y", color="#e5e5e5"); b.set_axisbelow(True)
    fig.tight_layout(w_pad=2); fig.savefig(FIG / "fig4_gap.png"); plt.close(fig)


def fig5():
    fig, (a, b) = plt.subplots(1, 2, figsize=(6.375, 2.425), dpi=200)
    rows = [(f"SARS-CoV-2\n({R1['SARS-CoV-2']['paired_days']} days)", R1["SARS-CoV-2"]),
            (f"Influenza A\n({R1['Influenza A']['paired_days']} days)", R1["Influenza A"])]
    for i, (lab, s) in enumerate(rows):
        y, m, (lo, hi) = 1 - i, s["sd_ratio"], s["sd_ratio_ci"]
        a.plot([lo, hi], [y, y], color=VIOLET, lw=3, solid_capstyle="round")
        a.plot(m, y, "o", color=VIOLET, ms=7, mec="white", mew=1)
        a.text(hi + 0.02, y, f"{m:.2f}  [{lo:.2f}, {hi:.2f}]", va="center", fontsize=9.5)
    a.axvline(1, ls="--", color=GRAY, lw=1); a.text(1, -0.55, "no change", ha="center", fontsize=9, color="#666")
    a.set_yticks([1, 0]); a.set_yticklabels([r[0] for r in rows], fontsize=10.5)
    a.set_ylim(-0.7, 1.6); a.set_xlim(0.85, 1.75); a.set_xticks([1.0, 1.2, 1.4, 1.6])
    a.set_xlabel("Spread after ÷ PMMoV, relative to raw\n(1.0 = no change)", fontsize=10)
    a.set_title("a  Spread of the between-lab gap", loc="left", fontsize=10.5, fontweight="bold")
    a.grid(axis="x", color="#e5e5e5"); a.set_axisbelow(True)
    shared = [round(100 * R1[k]["shared_part_cancelled"]) for k in ["SARS-CoV-2", "Influenza A"]]
    own = [round(100 * R1[k]["pmmov_noise_added"]) for k in ["SARS-CoV-2", "Influenza A"]]
    w, x = 0.36, np.array([0, 1])
    b1 = b.bar(x - w / 2, shared, w, color=VIOLET, label="Shared part cancelled")
    b2 = b.bar(x + w / 2, own, w, color="white", edgecolor=VIOLET, hatch="///", lw=1.3, label="PMMoV’s own noise added")
    for bars, vals in [(b1, shared), (b2, own)]:
        for r, v in zip(bars, vals):
            b.text(r.get_x() + r.get_width() / 2, v + 1.5, f"{v}%", ha="center", fontsize=9.5)
    b.set_xticks(x); b.set_xticklabels(["SARS-CoV-2", "Influenza A"], fontsize=10.5)
    b.set_ylim(0, 100); b.set_yticks([0, 25, 50, 75, 100]); b.set_ylabel("% of the labs’ disagreement", fontsize=10)
    b.set_title("b  What dividing by PMMoV does", loc="left", fontsize=10.5, fontweight="bold")
    b.legend(frameon=False, fontsize=9, loc="upper left"); b.grid(axis="y", color="#e5e5e5"); b.set_axisbelow(True)
    fig.tight_layout(w_pad=2); fig.savefig(FIG / "fig5_why_nil.png"); plt.close(fig)


# Justen et al. medians for SARS-CoV-2, from tables/pcr_normalization_summary_wide.tsv in their repository
JUSTEN = {"adm": [("Justen et al., 4 Chicago sites", 0.197), ("Justen et al., 5 WastewaterSCAN sites", -0.049)],
          "ed": [("Justen et al., 4 Chicago sites", -0.071), ("Justen et al., 5 WastewaterSCAN sites", 0.037)]}


def fig6():
    W, V = R2["SARS-CoV-2 | WastewaterSCAN"], R2["SARS-CoV-2 | CDC_Verily"]
    panels = {
        "adm": ("a  Against hospital admissions (NHSN)", [
            (f"Las Vegas, WastewaterSCAN ({W['admissions, same weeks as ED']['weeks']} weeks)", W["admissions, same weeks as ED"], BLUE),
            (f"Las Vegas, CDC/Verily ({V['admissions, same weeks as ED']['weeks']} weeks)", V["admissions, same weeks as ED"], ORANGE),
            (f"Las Vegas, WastewaterSCAN (all {W['admissions, full window']['weeks']} weeks)", W["admissions, full window"], BLUE),
            (f"Las Vegas, CDC/Verily (all {V['admissions, full window']['weeks']} weeks)", V["admissions, full window"], ORANGE)]),
        "ed": ("b  Against ED visits (Las Vegas) or Epic Cosmos clinical encounters (Justen et al.)", [
            (f"Las Vegas, WastewaterSCAN ({W['ED share, same weeks']['weeks']} weeks)", W["ED share, same weeks"], BLUE),
            (f"Las Vegas, CDC/Verily ({V['ED share, same weeks']['weeks']} weeks)", V["ED share, same weeks"], ORANGE)])}
    fig, axes = plt.subplots(2, 1, figsize=(7.205, 4.145), dpi=200, gridspec_kw={"height_ratios": [6, 4]})
    for ax, key in zip(axes, ["adm", "ed"]):
        title, rows = panels[key]
        rows = rows + [(n, v, None) for n, v in JUSTEN[key]]
        for i, (name, s, col) in enumerate(rows):
            y = len(rows) - 1 - i
            if col:
                ax.plot(s["ci"], [y, y], color=col, lw=2.2); ax.plot(s["dR"], y, "o", color=col, ms=7, mec="white", mew=1)
                ax.text(0.27, y, f"{s['dR']:+.3f}  [{s['ci'][0]:+.3f}, {s['ci'][1]:+.3f}]", va="center", fontsize=8)
            else:
                ax.plot(s, y, "D", mfc="white", mec="#555", ms=8, mew=1.5)
                ax.text(0.27, y, f"{s:+.3f}  (point estimate only)", va="center", fontsize=8, color="#555")
        ax.set_yticks(range(len(rows))); ax.set_yticklabels([r[0] for r in rows][::-1], fontsize=8)
        ax.axvline(0, ls="--", color=GRAY, lw=1); ax.set_xlim(-0.16, 0.5); ax.set_xticks([-0.1, 0, 0.1, 0.2])
        ax.set_ylim(-0.6, len(rows) - 0.4); ax.grid(axis="x", color="#e5e5e5"); ax.set_axisbelow(True)
        ax.spines["left"].set_color("#bbb")
        ax.set_title(title, loc="left", fontsize=9, fontweight="bold")
    axes[1].set_xlabel("ΔR: change in correlation with clinical data after dividing by PMMoV  (right of dashed line = normalizing helps)", fontsize=8)
    fig.tight_layout(h_pad=1.5); fig.savefig(FIG / "fig6_delta_r.png"); plt.close(fig)


def fig7():
    n = pd.DataFrame(json.load(open(RAW / "nhsn.json")))
    n["wk"] = pd.to_datetime(n["weekendingdate"]); n = n.set_index("wk").sort_index()
    n = n.loc["2023-10-01":"2025-06-30"]
    rep = pd.to_numeric(n["totalconfc19newadmhosprep"]); adm = pd.to_numeric(n["totalconfc19newadm"])
    fig, (a, b) = plt.subplots(2, 1, figsize=(5.83, 3.375), dpi=200, sharex=True)
    for ax in (a, b):
        ax.axvspan(*GAP, color="#ebebeb", lw=0)
    a.step(rep.index, rep, where="mid", color=BLUE, lw=1.8)
    a.text(pd.Timestamp("2024-01-15"), HOSPITALS_USUAL + 3, f"{HOSPITALS_USUAL} hospitals", ha="center", fontsize=9)
    a.text(pd.Timestamp("2024-05-28"), 10, "3–7 hospitals", fontsize=9)
    a.text(pd.Timestamp("2024-08-01"), 32, "voluntary reporting\n1 May – 1 Nov 2024", ha="center", fontsize=8.5, color="#777")
    a.set_ylim(0, 58); a.set_ylabel("Nevada hospitals\nreporting", fontsize=9.5); a.grid(axis="y", color="#e5e5e5")
    a.set_title("Weekly COVID-19 admissions reporting in Nevada (NHSN)", loc="left", fontsize=9.5, fontweight="bold")
    b.plot(adm.index, adm, "-o", color=BLUE, ms=2.5, lw=1.2)
    b.set_ylabel("New admissions\nreported per week", fontsize=9.5); b.set_ylim(0, None); b.grid(axis="y", color="#e5e5e5")
    b.xaxis.set_major_locator(mdates.MonthLocator((1, 4, 7, 10))); b.xaxis.set_major_formatter(mdates.DateFormatter("%b %Y"))
    fig.tight_layout(); fig.savefig(FIG / "fig7_reporting_gap.png"); plt.close(fig)


if __name__ == "__main__":
    for f in (fig3, fig4, fig5, fig6, fig7):
        f(); print("made", f.__name__)
