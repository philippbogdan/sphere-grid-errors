"""Render relaxation, resolution and restart/twin charts from the archived records."""
from __future__ import annotations

import json
import math
import sys

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.ticker import FixedLocator, NullFormatter, NullLocator  # noqa: E402

sys.path.insert(0, str(__import__("pathlib").Path(__file__).parent))
from common import (BLUE, GREY, GUIDE, INK, MAIN_C, MAIN_SCHEME, OUT, QA, VERM,  # noqa: E402
                    audit_text, rel, save_chart, select_one, style_rc)

# Marker geometry (points). Unfilled, nested sizes: a diamond fits inside the square and a circle inside the
# diamond, so coincident or near-coincident samples stay separately visible without moving any centre.
MS = {"s": 9.0, "D": 6.3, "o": 5.3, "^": 7.6}
LW = 1.6          # data lines (pt): 2.4 CSS px at 850 px display width
MEW = 1.45        # marker edge width (pt)
THIN = 1.1        # guides (pt): 1.6 CSS px
LABEL_PT = 11.0
TICK_PT = 10.5

L2 = r"$\mathit{l}_2$"     # typeset form of the whitelisted symbol "l2" (plain sans l reads as capital I)
TYPESET = {
    "Normalised height l2 vs T85 (1)": f"Normalised height {L2} vs T85 (1)",
    "Day-14 normalised height l2 (1)": f"Day-14 normalised height {L2} (1)",
    "Normalised height l2 vs T213 (1)": f"Normalised height {L2} vs T213 (1)",
    "Tetrahedral l = 3": r"Tetrahedral $\mathit{l}$ = 3",
}


def T(s):
    return TYPESET.get(s, s)


def log_ticks(ax, axis, values):
    """Fixed log ticks with explicit 'a x 10^b' labels; minor decade subdivisions unlabelled."""
    labels = []
    for v in values:
        e = math.floor(math.log10(v) + 1e-12)
        m = v / 10 ** e
        m_s = f"{m:.0f}" if abs(m - round(m)) < 1e-9 else f"{m:g}"
        labels.append(f"$10^{{{e}}}$" if m_s == "1" else f"${m_s}\\times10^{{{e}}}$")
    a = ax.yaxis if axis == "y" else ax.xaxis
    a.set_major_locator(FixedLocator(values))
    a.set_major_formatter(mpl.ticker.FixedFormatter(labels))
    a.set_minor_formatter(NullFormatter())


def style_axes(ax, grid_axis="y"):
    ax.tick_params(labelsize=TICK_PT, pad=4)
    ax.tick_params(which="minor", length=3)
    if grid_axis:
        ax.grid(True, axis=grid_axis, which="major", color="#E6E7E8", lw=1.05, zorder=0)
    ax.set_axisbelow(True)


def series(ax, x, y, color, marker, ls, label=None, z=3):
    return ax.plot(x, y, color=color, lw=LW, ls=ls, marker=marker, ms=MS[marker], mew=MEW,
                   mfc="none", mec=color, label=label, zorder=z, clip_on=False,
                   dash_capstyle="round", solid_capstyle="round")[0]


def fit_ylabels(fig, axes, margin_px=6):
    """Keep each rotated y label inside the canvas: slide it along its axis if it overruns an edge."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    H = fig.canvas.get_width_height()[1]
    for ax in axes:
        lab = ax.yaxis.label
        bb = lab.get_window_extent(r)
        ax_bb = ax.get_window_extent(r)
        shift = 0.0
        if bb.y1 > H - margin_px:
            shift = (H - margin_px) - bb.y1
        elif bb.y0 < margin_px:
            shift = margin_px - bb.y0
        if shift:
            x, _ = ax.yaxis.get_label_position(), None
            cx = ax.transAxes.inverted().transform((bb.x0 + bb.width / 2, 0))[0]
            cy = ax.transAxes.inverted().transform((0, (bb.y0 + bb.y1) / 2 + shift))[1]
            ax.yaxis.set_label_coords(cx, cy)


LS_DASH = (0, (5.0, 2.6))
LS_DOT = (0, (1.0, 2.2))


# ---------------------------------------------------------------- Visual 3: relaxation evidence
REL_ROWS = [  # label, region, par
    ("Gaussian 20 deg", "bp", 20),
    ("Disc 20 deg", "disc", 20),
    ("Annulus 10 to 20 deg", "annulus", 10),
    ("Disc 10 deg", "disc", 10),
    ("Gaussian 10 deg", "bp", 10),
    ("Disc 5 deg", "disc", 5),
    ("Polar caps 60 deg", "polar", 60),
]
DAYS = ["1", "3", "5", "7", "10"]


def relax_records():
    from common import EVID, load_jsonl
    rows = load_jsonl(EVID / "tc6_week2_relax.jsonl")
    assert len(rows) == 12, len(rows)
    out = {}
    for _, r in rows:
        assert r["M"] == 128 and r["scheme"] == MAIN_SCHEME and r["mode"] == "relax", r
    for _, r in rows:
        key = (r["region"], r["par"], r["control_rotation_deg"])
        rec, prov = select_one("tc6_week2_relax.jsonl",
                               {"M": 128, "region": key[0], "par": key[1], "control_rotation_deg": key[2],
                                "scheme": MAIN_SCHEME})
        out[key] = (rec, prov)
    assert len(out) == 12
    return out


def visual3():
    recs = relax_records()
    none = recs[("none", 0, 0)][0]
    base10 = none["l2_vs_T85_by_day"]["10"]
    red = {k: base10 / v[0]["l2_vs_T85_by_day"]["10"] for k, v in recs.items()}
    prov = {"inputs": [dict(p, used=["weighted_area", "l2_vs_T85_by_day"]) for (_, p) in recs.values()]}

    whitelist = ["Time (days)", "Normalised height l2 vs T85 (1)", "Day-10 error reduction (x)", "Main",
                 "Cones", "Control", "Polar caps", "None"] + [r[0] for r in REL_ROWS]
    with mpl.rc_context(style_rc(TICK_PT)):
        fig = plt.figure(figsize=(8.0, 6.0), dpi=300)
        x0, x1 = 0.215, 0.972
        axt = fig.add_axes([x0, 0.585, x1 - x0, 0.39])
        axb = fig.add_axes([x0, 0.085, x1 - x0, 0.365])

        # Top: recorded days only, straight guides between samples, no day-0 point.
        d = [int(k) for k in DAYS]
        spec = [(("none", 0, 0), "Main", MAIN_C, "s", "-"),
                (("bp", 20, 0), "Cones", BLUE, "D", "-"),
                (("bp", 20, 60), "Control", VERM, "o", LS_DASH)]
        handles = []
        for key, lab, c, m, ls in spec:
            y = [recs[key][0]["l2_vs_T85_by_day"][k] for k in DAYS]
            handles.append(series(axt, d, y, c, m, ls, lab, z=4 if m == "o" else 3))
        axt.set_yscale("log")
        axt.set_ylim(2e-5, 4e-3)
        axt.set_xlim(0, 10)
        axt.xaxis.set_major_locator(FixedLocator([0, 1, 3, 5, 7, 10]))
        axt.xaxis.set_minor_locator(NullLocator())
        log_ticks(axt, "y", [1e-4, 1e-3])
        axt.yaxis.set_minor_locator(mpl.ticker.LogLocator(base=10, subs=range(2, 10)))
        axt.set_xlabel(T("Time (days)"), fontsize=LABEL_PT, labelpad=5)
        axt.set_ylabel(T("Normalised height l2 vs T85 (1)"), fontsize=LABEL_PT, labelpad=8)
        style_axes(axt, "y")
        axt.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.0, 1.02), fontsize=TICK_PT,
                   handlelength=3.0, handletextpad=0.7, borderaxespad=0.2, labelspacing=0.45)

        # Bottom: all 12 configurations as day-10 reductions relative to none.
        for i, (lab, reg, par) in enumerate(REL_ROWS):
            axb.axhline(i, color="#EEEFF0", lw=1.05, zorder=0)
            cone = red.get((reg, par, 0))
            ctrl = red.get((reg, par, 60))
            if reg == "polar":
                axb.plot([cone], [i], ls="none", marker="^", ms=MS["^"], mew=MEW, mfc="none", mec=GREY,
                         zorder=4, clip_on=False)
                continue
            if ctrl is not None:
                axb.plot([ctrl, cone], [i, i], color=GUIDE, lw=THIN, zorder=2, solid_capstyle="butt")
                axb.plot([ctrl], [i], ls="none", marker="o", ms=MS["o"], mew=MEW, mfc="none", mec=VERM,
                         zorder=4, clip_on=False)
            axb.plot([cone], [i], ls="none", marker="D", ms=MS["D"], mew=MEW, mfc="none", mec=BLUE,
                     zorder=3, clip_on=False)
        axb.axvline(1.0, color=INK, lw=1.3, ls=LS_DOT, zorder=1)
        axb.set_xlim(0, 6)
        axb.set_ylim(len(REL_ROWS) - 0.45, -0.55)
        axb.yaxis.set_major_locator(FixedLocator(range(len(REL_ROWS))))
        axb.yaxis.set_major_formatter(mpl.ticker.FixedFormatter([T(r[0]) for r in REL_ROWS]))
        axb.tick_params(axis="y", length=0, pad=8)
        axb.spines["left"].set_visible(False)
        axb.xaxis.set_major_locator(FixedLocator(range(0, 7)))
        axb.set_xlabel(T("Day-10 error reduction (x)"), fontsize=LABEL_PT, labelpad=5)
        style_axes(axb, None)
        axb.tick_params(axis="y", length=0)
        lh = [Line2D([], [], ls="none", marker="D", ms=MS["D"], mew=MEW, mfc="none", mec=BLUE),
              Line2D([], [], ls="none", marker="o", ms=MS["o"], mew=MEW, mfc="none", mec=VERM),
              Line2D([], [], ls="none", marker="^", ms=MS["^"], mew=MEW, mfc="none", mec=GREY),
              Line2D([], [], color=INK, lw=1.3, ls=LS_DOT)]
        axb.legend(lh, ["Cones", "Control", "Polar caps", "None"], loc="lower right", fontsize=TICK_PT,
                   handlelength=2.0, handletextpad=0.6, borderaxespad=0.3, labelspacing=0.45)

        audit = audit_text(fig, whitelist, TYPESET)
        png, svg = save_chart(fig, "relaxation", (2400, 1800))
        plt.close(fig)

    plotted = {
        "top": {lab: {"days": d, "l2_vs_T85": [recs[key][0]["l2_vs_T85_by_day"][k] for k in DAYS]}
                for key, lab, *_ in spec},
        "bottom": {lab: {"cones_rot0": red.get((reg, par, 0)), "control_rot60": red.get((reg, par, 60))}
                   for lab, reg, par in REL_ROWS},
        "weighted_area": {f"{k[0]},{k[1]},{k[2]}": v[0]["weighted_area"] for k, v in recs.items()},
    }
    checks = {
        "cones_bp20_day10_reduction": {"value": red[("bp", 20, 0)], "expected": 5.33265185658},
        "control_bp20_rot60_day10_reduction": {"value": red[("bp", 20, 60)], "expected": 2.81015510331},
        "polar60_day10_reduction": {"value": red[("polar", 60, 0)], "expected_display": "1.40"},
        "pairs_favouring_cones": sum(red[(r, p, 0)] > red[(r, p, 60)] for _, r, p in REL_ROWS[:4]),
        "weighted_area_pct_cones_bp20": 100 * recs[("bp", 20, 0)][0]["weighted_area"],
        "weighted_area_pct_control_bp20": 100 * recs[("bp", 20, 60)][0]["weighted_area"],
        "configurations_plotted": 1 + 8 + 2 + 1,
    }
    for c in checks.values():
        if isinstance(c, dict) and "expected" in c:
            c["abs_diff"] = abs(c["value"] - c["expected"])
    return {"outputs": [png, svg], "audit": audit, "provenance": prov, "plotted": plotted, "checks": checks,
            "limits": {"top_x": [0, 10], "top_y_log": [2e-5, 4e-3], "bottom_x": [0, 6]}}


# ---------------------------------------------------------------- Visual 4: resolution
RES = [(128, 42, 900.0), (256, 85, 450.0), (512, 170, 225.0)]


def visual4():
    rows, prov = [], []
    for M, trunc, dt in RES:
        main, pm = select_one("accuracy_main_ref.jsonl",
                              {"case": "tc6", "days": 14.0, "grid": ["torus", M],
                               "spec": f"torus:{M}:tc6:{MAIN_SCHEME}"},
                              extra=lambda r: r["scheme"]["courant"] == 1.252)
        trsk, pt = select_one("mpas_ref.jsonl", {"solver": "trsk", "case": "tc6", "days": 14.0,
                                                 "grid": ["torus", M], "mesh": "torus"})
        spec_, ps = select_one("stswm_errors.jsonl", {"case": "tc6", "against": "T213dt22",
                                                      "truncation": trunc, "dt": dt})
        assert main["nunique"] == trsk["nunique"] == spec_["npts"] == M * M // 2, (main["nunique"], trsk["nunique"], spec_["npts"])
        rows.append({"M": M, "n": M * M // 2, "main": main["l2"], "trsk": trsk["l2"], "spectral": spec_["l2"],
                     "main_over_spectral": main["l2"] / spec_["l2"], "trsk_over_spectral": trsk["l2"] / spec_["l2"]})
        prov += [dict(pm, used=["l2", "nunique"]), dict(pt, used=["l2", "nunique"]), dict(ps, used=["l2", "npts"])]

    n = [r["n"] for r in rows]
    whitelist = ["Unique points (count)", "Day-14 normalised height l2 (1)", "Error / spectral error (1)", "Main",
                 "TRSK on torus points", "Spectral"]
    xlim = (n[0] / 1.45, n[-1] * 1.45)
    with mpl.rc_context(style_rc(TICK_PT)):
        fig = plt.figure(figsize=(8.0, 1600 / 300), dpi=300)
        x0, x1 = 0.155, 0.972
        axt = fig.add_axes([x0, 0.465, x1 - x0, 0.51])
        axb = fig.add_axes([x0, 0.115, x1 - x0, 0.29], sharex=axt)
        h = [series(axt, n, [r["main"] for r in rows], MAIN_C, "s", "-", "Main", z=3),
             series(axt, n, [r["trsk"] for r in rows], BLUE, "D", LS_DASH, "TRSK on torus points", z=4),
             series(axt, n, [r["spectral"] for r in rows], VERM, "o", LS_DOT, "Spectral", z=5)]
        axt.set_xscale("log")
        axt.set_yscale("log")
        axt.set_ylim(5e-4, 2e-2)
        log_ticks(axt, "y", [5e-4, 1e-3, 2e-3, 5e-3, 1e-2, 2e-2])
        axt.yaxis.set_minor_locator(mpl.ticker.LogLocator(base=10, subs=range(2, 10)))
        axt.set_ylabel(T("Day-14 normalised height l2 (1)"), fontsize=LABEL_PT, labelpad=8)
        style_axes(axt, "y")
        axt.tick_params(axis="x", labelbottom=False)
        axt.legend(handles=h, loc="upper right", fontsize=TICK_PT, handlelength=3.0, handletextpad=0.7,
                   borderaxespad=0.2, labelspacing=0.45)

        axb.axhline(1.0, color=VERM, lw=THIN, ls=LS_DOT, zorder=1)
        series(axb, n, [r["main_over_spectral"] for r in rows], MAIN_C, "s", "-", z=3)
        series(axb, n, [r["trsk_over_spectral"] for r in rows], BLUE, "D", LS_DASH, z=4)
        axb.set_ylim(0, 2.2)
        axb.yaxis.set_major_locator(FixedLocator([0, 0.5, 1.0, 1.5, 2.0]))
        axb.yaxis.set_major_formatter(mpl.ticker.FixedFormatter(["0", "0.5", "1", "1.5", "2"]))
        axb.set_xlim(*xlim)
        axb.xaxis.set_major_locator(FixedLocator(n))
        axb.xaxis.set_major_formatter(mpl.ticker.FixedFormatter([f"{v:,}" for v in n]))
        axb.xaxis.set_minor_locator(NullLocator())
        axt.xaxis.set_minor_locator(NullLocator())
        axb.set_xlabel(T("Unique points (count)"), fontsize=LABEL_PT, labelpad=5)
        axb.set_ylabel(T("Error / spectral error (1)"), fontsize=LABEL_PT, labelpad=8)
        style_axes(axb, "y")
        # align y labels of the two panels
        fig.align_ylabels([axt, axb])
        audit = audit_text(fig, whitelist, TYPESET)
        png, svg = save_chart(fig, "resolution", (2400, 1600))
        plt.close(fig)
    checks = {f"M{r['M']}": {"main_over_spectral": r["main_over_spectral"],
                             "trsk_over_spectral": r["trsk_over_spectral"]} for r in rows}
    checks["expected_display"] = {"main_over_spectral": ["1.86", "1.94", "1.92"],
                                  "trsk_over_spectral": ["1.91", "1.54", "1.046"],
                                  "fine_trsk_ratio_exact": 1.04582187688}
    checks["fine_trsk_ratio_abs_diff"] = abs(rows[-1]["trsk_over_spectral"] - 1.04582187688)
    return {"outputs": [png, svg], "audit": audit, "provenance": {"inputs": prov}, "plotted": rows,
            "checks": checks, "limits": {"x_log": list(xlim), "top_y_log": [5e-4, 2e-2], "bottom_y": [0, 2.2]}}


# ---------------------------------------------------------------- Visual 5: restart and twins
TWINS = [("random", "Smooth random", MAIN_C, "s", "-"),
         ("zonal_P2", "Zonal P2", BLUE, "o", LS_DASH),
         ("tetrahedral_l3", "Tetrahedral l = 3", VERM, "D", LS_DOT)]


def visual5():
    s0, p0 = select_one("tc6_week2_restart.jsonl", {"M": 128, "label": "main", "start_day": 0},
                        extra=lambda r: r["scheme"] == MAIN_SCHEME)
    s7, p7 = select_one("tc6_week2_restart.jsonl", {"M": 128, "label": "main", "start_day": 7},
                        extra=lambda r: r["scheme"] == MAIN_SCHEME)
    a_days = sorted(int(k) for k in s0["l2_by_day"])
    b_days = sorted(int(k) for k in s7["l2_by_day"])
    assert a_days == list(range(1, 8)) and b_days == list(range(8, 15)), (a_days, b_days)
    a = [(d, s0["l2_by_day"][str(d)]) for d in a_days]
    b = [(d - 7, s7["l2_by_day"][str(d)]) for d in b_days]
    tw = {}
    for shape, *_ in TWINS:
        rec, p = select_one("tc6_week2_twins7.jsonl", {"M": 128, "shape": shape})
        days = sorted(int(k) for k in rec["rms_difference_m_by_day"])
        assert days == list(range(8, 15)), days
        tw[shape] = (rec, p, [7] + days, [rec["initial_rms_m"]] + [rec["rms_difference_m_by_day"][str(d)] for d in days])

    whitelist = ["Elapsed time (days)", "Normalised height l2 vs T213 (1)", "Start at day 0", "Restart at day 7",
                 "Time (days)", "Twin height difference RMS (m)", "Smooth random", "Zonal P2", "Tetrahedral l = 3"]
    with mpl.rc_context(style_rc(TICK_PT)):
        fig = plt.figure(figsize=(8.0, 1600 / 300), dpi=300)
        x0, x1 = 0.155, 0.972
        axt = fig.add_axes([x0, 0.565, x1 - x0, 0.42])
        axb = fig.add_axes([x0, 0.105, x1 - x0, 0.345])
        h = [series(axt, [p[0] for p in a], [p[1] for p in a], MAIN_C, "s", "-", "Start at day 0"),
             series(axt, [p[0] for p in b], [p[1] for p in b], BLUE, "D", LS_DASH, "Restart at day 7", z=4)]
        axt.set_yscale("log")
        axt.set_ylim(2e-5, 1e-2)
        axt.set_xlim(0, 7)
        axt.xaxis.set_major_locator(FixedLocator(range(0, 8)))
        log_ticks(axt, "y", [1e-4, 1e-3, 1e-2])
        axt.yaxis.set_minor_locator(mpl.ticker.LogLocator(base=10, subs=range(2, 10)))
        axt.set_xlabel(T("Elapsed time (days)"), fontsize=LABEL_PT, labelpad=5)
        axt.set_ylabel(T("Normalised height l2 vs T213 (1)"), fontsize=LABEL_PT, labelpad=8)
        style_axes(axt, "y")
        axt.legend(handles=h, loc="lower right", fontsize=TICK_PT, handlelength=3.0, handletextpad=0.7,
                   borderaxespad=0.3, labelspacing=0.45)

        hb = []
        for shape, lab, c, m, ls in TWINS:
            _, _, x, y = tw[shape]
            hb.append(series(axb, x, y, c, m, ls, T(lab), z={"s": 3, "D": 4, "o": 5}[m]))
        axb.set_xlim(7, 14)
        axb.set_ylim(0, 3)
        axb.xaxis.set_major_locator(FixedLocator(range(7, 15)))
        axb.yaxis.set_major_locator(FixedLocator([0, 1, 2, 3]))
        axb.set_xlabel(T("Time (days)"), fontsize=LABEL_PT, labelpad=5)
        axb.set_ylabel(T("Twin height difference RMS (m)"), fontsize=LABEL_PT, labelpad=8)
        style_axes(axb, "y")
        axb.legend(handles=hb, loc="upper left", fontsize=TICK_PT, handlelength=3.0, handletextpad=0.7,
                   borderaxespad=0.2, labelspacing=0.45)
        fig.align_ylabels([axt, axb])
        fit_ylabels(fig, [axt, axb])
        audit = audit_text(fig, whitelist, TYPESET)
        png, svg = save_chart(fig, "restart-twins", (2400, 1600))
        plt.close(fig)
    ratio = s7["l2_by_day"]["10"] / s0["l2_by_day"]["3"]
    checks = {"restart_ratio_elapsed3": {"value": ratio, "expected": 12.7711668561,
                                         "abs_diff": abs(ratio - 12.7711668561)},
              "twin_day14_rms_m": {lab: tw[s][0]["rms_difference_m_by_day"]["14"] for s, lab, *_ in TWINS},
              "twin_initial_rms_m": {lab: tw[s][0]["initial_rms_m"] for s, lab, *_ in TWINS}}
    prov = {"inputs": [dict(p0, used=["l2_by_day"]), dict(p7, used=["l2_by_day"])] +
            [dict(tw[s][1], used=["initial_rms_m", "rms_difference_m_by_day"]) for s, *_ in TWINS]}
    plotted = {"top": {"Start at day 0": a, "Restart at day 7": b},
               "bottom": {lab: list(zip(tw[s][2], tw[s][3])) for s, lab, *_ in TWINS}}
    return {"outputs": [png, svg], "audit": audit, "provenance": prov, "plotted": plotted, "checks": checks,
            "limits": {"top_x": [0, 7], "top_y_log": [2e-5, 1e-2], "bottom_x": [7, 14], "bottom_y": [0, 3]}}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    QA.mkdir(parents=True, exist_ok=True)
    which = sys.argv[1:] or ["3", "4", "5"]
    res = {}
    for k in which:
        r = {"3": visual3, "4": visual4, "5": visual5}[k]()
        r["outputs"] = [rel(p) for p in r["outputs"]]
        res[k] = r
        print(k, "audit ok" if r["audit"]["ok"] else f"AUDIT FAIL {r['audit']['violations']}")
        print(json.dumps(r["checks"], indent=1, default=str))
    prev = QA / "charts-provenance.json"
    old = json.loads(prev.read_text()) if prev.exists() else {}
    old.update(res)
    prev.write_text(json.dumps(old, indent=1, default=str))


if __name__ == "__main__":
    main()
