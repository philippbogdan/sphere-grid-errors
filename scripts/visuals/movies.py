"""Render the flow and correction movies defined in docs/FIELDS.md.

Run movies.py all to check fields, render the day-10 hero poster and encode both movies.
Use --data, --out and --qa to select local input and output directories.
"""
from __future__ import annotations

import argparse
import glob
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.patheffects as pe  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402
from matplotlib.legend_handler import HandlerTuple  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Ellipse, Patch  # noqa: E402

sys.path.insert(0, str(Path(__file__).parent))
import fields as F  # noqa: E402
from common import (DATA, FONT, INK, OUT, QA, SCRATCH, audit_text, rel, save_rgb_png, select_one,  # noqa: E402
                    sha256)

FFMPEG = os.environ.get("FFMPEG", "ffmpeg")
FS_TICK = 18.0      # pt at 100 dpi = 25 px in the movie, 50 px in the poster
FS_LABEL = 19.0
EDGE = "#80868B"
GRAT = (0.5, 0.5, 0.5, 0.42)
TRACK = "#DADCE0"
CB_EDGE = "#9AA0A6"
LATC = -19.4712206
RUN_NAMES = ("main", "cones20", "controls20")
REFS = ("T85", "T213")
ARCHIVE_KEY = {"main": ("none", 0, 0), "cones20": ("bp", 20, 0), "controls20": ("bp", 20, 60)}


def rc():
    return {"font.family": FONT, "mathtext.fontset": "dejavusans", "mathtext.default": "regular",
            "axes.unicode_minus": True,
            "text.color": INK, "axes.labelcolor": INK, "xtick.color": "#5F6368", "xtick.labelcolor": INK,
            "ytick.color": "#5F6368", "figure.facecolor": "white", "path.simplify": False}


# ====================================================================== data
def load_inputs(data: Path):
    runs, side, files = {}, {}, {}
    for r in RUN_NAMES:
        p = data / f"{r}.npz"
        with np.load(p) as z:
            runs[r] = {k: z[k] for k in z.files}
        files[f"{r}.npz"] = {"path": rel(p), "sha256": sha256(p), "bytes": p.stat().st_size,
                             "arrays": {k: [list(v.shape), str(v.dtype)] for k, v in runs[r].items()}}
        js = data / f"{r}.json"
        if js.exists():
            side[r] = json.loads(js.read_text())
            files[f"{r}.json"] = {"path": rel(js), "sha256": sha256(js)}
    chunks = sorted(glob.glob(str(data / "t213-vorticity-*.npz")))
    if not chunks:
        raise FileNotFoundError("no t213-vorticity-*.npz chunks")
    T, Z, lat, lon = [], [], None, None
    for c in chunks:
        with np.load(c) as z:
            T.append(np.asarray(z["t_hours"], np.float64))
            Z.append(np.asarray(z["zeta_s_inv"]))
            if lat is None:
                lat, lon = np.asarray(z["lat_rad"], np.float64), np.asarray(z["lon_rad"], np.float64)
            else:
                assert np.array_equal(lat, z["lat_rad"]) and np.array_equal(lon, z["lon_rad"])
        files[Path(c).name] = {"path": rel(c), "sha256": sha256(c), "bytes": Path(c).stat().st_size,
                               "frames": int(T[-1].size)}
    vort = {"t_hours": np.concatenate(T), "zeta": np.concatenate(Z), "lat_rad": lat, "lon_rad": lon}
    return runs, side, vort, files


def check_inputs(runs, side, vort):
    """Renderer-side checks before any field is drawn. Returns a dict; raises on schema violations."""
    out = {"runs": {}}
    M = 128
    t0 = runs["main"]["t_hours"]
    for r, d in runs.items():
        nt = d["t_hours"].size
        for k, shp in (("X", (3, M, M)), ("weight", (M, M)), ("relax_weight", (M, M)), ("centres", (4, 3)),
                       ("err_vs_T85", (nt, M, M)), ("err_vs_T213", (nt, M, M))):
            assert d[k].shape == shp, (r, k, d[k].shape, shp)
        pr = F.Pairs(d["X"], d["weight"])
        rec = {"nt": int(nt), "t_first_last_h": [float(d["t_hours"][0]), float(d["t_hours"][-1])],
               "times_equal_main": bool(np.array_equal(d["t_hours"], t0)),
               "max_abs_time_diff_vs_main_h": float(np.abs(d["t_hours"] - t0).max()) if nt == t0.size else None,
               "pairs": pr.stats(),
               "nan_counts": {k: int(np.isnan(d[k]).sum()) for k in ("err_vs_T85", "err_vs_T213", "relax_weight")}}
        # actual W against the nearest-centre Gaussian of the stated centres
        P = np.asarray(d["X"], np.float64).reshape(3, -1).T
        dist = np.degrees(np.arccos(np.clip(P @ np.asarray(d["centres"], np.float64).T, -1, 1))).min(1)
        Wan = np.exp(-(dist / 20.0) ** 2) if r != "main" else np.zeros_like(dist)
        rec["relax_weight_max_abs_diff_vs_formula"] = float(np.abs(d["relax_weight"].ravel() - Wan).max())
        w = d["weight"].astype(np.float64).ravel()
        rec["weighted_area_from_fields"] = float((w * d["relax_weight"].ravel()).sum() / w.sum())
        if r in side:
            rec["weighted_area_sidecar"] = side[r].get("weighted_area")
        # l2 from the stored fields at exact day endpoints, all M^2 points, original weights
        rec["l2"] = {}
        if "h" in d:
            for ref, key in (("T85", "l2_vs_T85_by_day"), ("T213", "l2_vs_T213dt22_by_day")):
                e_all = d[f"err_vs_{ref}"]
                for day in (1, 3, 5, 7, 10, 14):
                    hit = np.nonzero(np.abs(d["t_hours"] - 24.0 * day) < 1e-6)[0]
                    if hit.size != 1:
                        continue
                    e = e_all[hit[0]].astype(np.float64).ravel()
                    href = d["h"][hit[0]].astype(np.float64).ravel() - e
                    l2 = math.sqrt((w * e * e).sum() / (w * href * href).sum())
                    item = {"from_float32_fields": l2}
                    sv = side.get(r, {}).get(key, {}).get(str(day))
                    if sv is not None:
                        item["sidecar"] = sv
                        item["rel_diff_vs_sidecar"] = abs(l2 - sv) / sv
                    if ref == "T85" and day <= 10:
                        reg, par, rot = ARCHIVE_KEY[r]
                        arc, prov = select_one("tc6_week2_relax.jsonl", {"M": 128, "region": reg, "par": par,
                                                                          "control_rotation_deg": rot})
                        av = arc["l2_vs_T85_by_day"][str(day)]
                        item["archive"] = av
                        item["archive_line"] = prov["line"]
                        item["rel_diff_vs_archive"] = abs(l2 - av) / av
                        if sv is not None:
                            item["sidecar_rel_diff_vs_archive"] = abs(sv - av) / av
                    rec["l2"][f"{ref}_day{day}"] = item
        out["runs"][r] = rec
    # day-10 reductions from the sidecars (headline arithmetic is against T85)
    try:
        out["sidecar_day10_reduction_vs_T85"] = {
            r: side["main"]["l2_vs_T85_by_day"]["10"] / side[r]["l2_vs_T85_by_day"]["10"] for r in ("cones20", "controls20")}
        out["sidecar_day10_reduction_vs_T213"] = {
            r: side["main"]["l2_vs_T213dt22_by_day"]["10"] / side[r]["l2_vs_T213dt22_by_day"]["10"]
            for r in ("cones20", "controls20")}
    except (KeyError, ZeroDivisionError, TypeError):
        pass
    tv = vort["t_hours"]
    z = vort["zeta"]
    assert z.shape[1:] == (vort["lat_rad"].size, vort["lon_rad"].size), z.shape
    out["vorticity"] = {"records": int(tv.size), "grid": [int(z.shape[2]), int(z.shape[1])],
                        "t_first_last_h": [float(tv[0]), float(tv[-1])],
                        "t_spacing_h": sorted(set(np.round(np.diff(tv), 9).tolist())),
                        "lat_ascending": bool(np.all(np.diff(vort["lat_rad"]) > 0)),
                        "lat_extremes_deg": [float(np.degrees(vort["lat_rad"][0])), float(np.degrees(vort["lat_rad"][-1]))],
                        "nan_count": int(np.isnan(z).sum()), "dtype": str(z.dtype)}
    return out


def aligned(series, t_src, t_dst):
    """series: (nt, ...) at t_src; returns a function k -> value at t_dst[k] (linear, no extrapolation)."""
    j, a = F.bracket(t_src, t_dst)

    def get(k):
        if a[k] == 0.0:
            return np.asarray(series[j[k]], np.float64)
        return (1 - a[k]) * np.asarray(series[j[k]], np.float64) + a[k] * np.asarray(series[j[k] + 1], np.float64)
    return get, j, a


class Model:
    """Everything shared by both movies: pairs, triangulation, aligned E_pair series, frozen scales, schedule."""

    def __init__(self, data: Path):
        self.runs, self.side, self.vort, self.files = load_inputs(data)
        self.checks = check_inputs(self.runs, self.side, self.vort)
        m = self.runs["main"]
        self.t = np.asarray(m["t_hours"], np.float64)
        self.pairs = F.Pairs(m["X"], m["weight"])
        for r in RUN_NAMES[1:]:
            d = self.runs[r]
            assert np.allclose(np.asarray(d["X"], np.float64), np.asarray(m["X"], np.float64), atol=1e-6)
        self.tri = F.SphereTri(self.pairs.U)
        self.pairs_by_run = {r: F.Pairs(self.runs[r]["X"], self.runs[r]["weight"]) for r in RUN_NAMES}
        # E_max over all three runs, both references, all saved times (before any alignment)
        emax = 0.0
        self.emax_by = {}
        for r in RUN_NAMES:
            for ref in REFS:
                E = self.pairs_by_run[r].epair(self.runs[r][f"err_vs_{ref}"])
                v = float(np.nanmax(E))
                self.emax_by[f"{r}/{ref}"] = v
                emax = max(emax, v)
        self.E_largest = emax
        self.E_max = F.nice_ceiling(emax)
        tv = self.vort["t_hours"]
        jz, az = F.bracket(tv, self.t)
        used = sorted(set(jz.tolist()) | set((jz[az > 0] + 1).tolist()))
        zmax = float(np.nanmax(np.abs(self.vort["zeta"][used])))
        self.zeta_largest, self.zeta_lim = zmax, F.nice_ceiling(zmax)
        self.vort_get, self.vort_j, self.vort_a = aligned(self.vort["zeta"], tv, self.t)
        self.err_get, self.align = {}, {}
        for r in RUN_NAMES:
            for ref in REFS:
                g, j, a = aligned(self.runs[r][f"err_vs_{ref}"], self.runs[r]["t_hours"], self.t)
                self.err_get[(r, ref)] = g
                self.align[r] = {"j": j.tolist(), "a": a.tolist(), "identity": bool(np.all(a == 0) and
                                                                                     np.array_equal(j, np.arange(self.t.size)))}
        self.tau, self.kf = F.display_schedule(self.t)
        self.k_day = {d: int(np.nonzero(np.abs(self.t - 24.0 * d) < 1e-6)[0][0]) for d in (0, 7, 10, 14)}
        self.cones = np.asarray(self.runs["cones20"]["centres"], np.float64)
        self.controls = np.asarray(self.runs["controls20"]["centres"], np.float64)

    def epair(self, run, ref, k):
        return self.pairs_by_run[run].epair(self.err_get[(run, ref)](k))

    def summary(self):
        f_day = {d: int(np.nonzero(self.kf == self.k_day[d])[0][0]) for d in (7, 10)}
        return {"E_largest_m": self.E_largest, "E_max_m": self.E_max, "E_largest_by_run_ref": self.emax_by,
                "zeta_largest_abs_s_inv": self.zeta_largest, "zeta_limit_s_inv": self.zeta_lim,
                "triangulation": self.tri.check, "pairs_main": self.pairs.stats(),
                "t_hours_main": self.t.tolist(), "snapshot_index_by_day": self.k_day,
                "first_frame_showing_day": f_day,
                "frame_schedule": {"tau_hours": self.tau.tolist(), "snapshot_index": self.kf.tolist()},
                "vorticity_alignment_to_main_times": {"j": self.vort_j.tolist(), "a": self.vort_a.tolist(),
                                                       "reference_t_hours": self.vort["t_hours"].tolist()},
                "run_alignment_to_main_times": {r: {"identity": v["identity"]} | ({} if v["identity"] else v)
                                                for r, v in self.align.items()}}


# ====================================================================== drawing helpers
def halo_marker(ax, x, y, marker, ms, z=20, transform=None):
    kw = dict(ls="none", clip_on=False, zorder=z, mfc="none")
    if transform is not None:
        kw["transform"] = transform
    ax.plot(x, y, marker=marker, ms=ms, mew=4.6, mec="white", **kw)
    ax.plot(x, y, marker=marker, ms=ms, mew=1.9, mec=INK, **kw)


def marker_handle(marker, ms):
    return (Line2D([], [], ls="none", marker=marker, ms=ms, mew=4.6, mec="white", mfc="none"),
            Line2D([], [], ls="none", marker=marker, ms=ms, mew=1.9, mec=INK, mfc="none"))


def draw_contour(ax, XX, YY, Z, level, ls, lw, clip):
    c = ax.contour(XX, YY, Z, levels=[level], colors=INK, linewidths=lw, zorder=7)
    c.set_linestyle(ls)
    c.set_path_effects([pe.Stroke(linewidth=lw + 2.2, foreground="white"), pe.Normal()])
    c.set_clip_path(clip)
    return c


def colourbar(fig, box, scale, ticks, labels, label, W, H, s):
    x, y, w, h = box
    ax = fig.add_axes([x / W, 1 - (y + h) / H, w / W, h / H])
    ax.imshow(scale.lut[None, :, :], extent=(0, 1, 0, 1), aspect="auto", interpolation="nearest")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.set_xticks([float(scale.u(t)) for t in ticks])
    ax.set_xticklabels(labels, fontsize=FS_TICK)
    ax.tick_params(axis="x", length=6, width=1.0, pad=3, color="#5F6368")
    for sp in ax.spines.values():
        sp.set_color(CB_EDGE)
        sp.set_linewidth(0.8)
    ax.set_xlabel(label, fontsize=FS_LABEL, labelpad=5)
    return ax


def asinh_ticks(E_max):
    vals = [v for v in (0, 1, 3, 10, 30, 100, 300, 1000) if v <= E_max]
    return vals, [str(v) for v in vals]


def vort_ticks(L):
    e = math.floor(math.log10(L))
    m = L / 10 ** e
    step = {1: 0.5, 2: 1, 5: 2.5}[int(round(m))]
    if step == 2.5:
        step = 2
    vals = [i * step * 10 ** e for i in range(-int(m // step), int(m // step) + 1)]

    def lab(v):
        if abs(v) < 1e-30:
            return "0"
        mm = v / 10 ** e
        ms = f"{mm:g}".replace("-", "−")
        return f"${ms}\\times10^{{{e}}}$"
    return vals, [lab(v) for v in vals]


class TimeAxis:
    def __init__(self, fig, box, W, H, label_right_px, legend=None):
        x, y, w, h = box
        self.ax = ax = fig.add_axes([x / W, 1 - (y + h) / H, w / W, h / H])
        ax.set_xlim(0, 14)
        ax.set_ylim(-1, 1)
        ax.axis("off")
        ax.plot([0, 14], [0, 0], color=TRACK, lw=5.0, solid_capstyle="round", zorder=1)
        self.done, = ax.plot([0, 0], [0, 0], color=INK, lw=5.0, solid_capstyle="round", zorder=2)
        for d in (0, 7, 10, 14):
            ax.plot([d, d], [-0.55, -1.0], color="#5F6368", lw=1.4, zorder=1, clip_on=False)
            ax.text(d, -1.25, str(d), ha="center", va="top", fontsize=FS_TICK, clip_on=False)
        self.cur, = ax.plot([0], [0], ls="none", marker="o", ms=13, mfc=INK, mec="white", mew=2.2, zorder=3,
                            clip_on=False)
        fig.text(label_right_px / W, 1 - (y + h / 2) / H, "Time (days)", ha="right", va="center", fontsize=FS_LABEL)

    def set(self, t_days):
        self.done.set_data([0, t_days], [0, 0])
        self.cur.set_data([t_days], [0])


def graticule_ortho(ax, lon0, lat0, step=30):
    segs = []
    s = np.radians(np.arange(-180, 180.01, 0.5))
    for la in range(-90 + step, 90, step):
        P = F.latlon_to_vec(np.full_like(s, math.radians(la)), s)
        segs += _visible_segments(P, lon0, lat0)
    p = np.radians(np.arange(-80, 80.01, 0.5))          # meridians stop short of the poles
    for lo in range(0, 360, step):
        P = F.latlon_to_vec(p, np.full_like(p, math.radians(lo)))
        segs += _visible_segments(P, lon0, lat0)
    ax.add_collection(LineCollection(segs, colors=[GRAT], linewidths=0.85, zorder=5, capstyle="round"))


def _visible_segments(P, lon0, lat0):
    x, y, d = F.ortho_forward(P, lon0, lat0)
    vis = d > 0
    out, cur = [], []
    for i in range(len(x)):
        if vis[i]:
            cur.append((x[i], y[i]))
        elif cur:
            if len(cur) > 1:
                out.append(np.array(cur))
            cur = []
    if len(cur) > 1:
        out.append(np.array(cur))
    return out


def graticule_moll(ax, lon0, step=30):
    segs = []
    t80 = F.moll_theta(np.array([math.radians(80.0)]))[0]
    th = np.linspace(-t80, t80, 721)                     # meridians stop short of the poles
    for k in range(-180 + step, 180, step):
        lam = math.radians(k)
        segs.append(np.stack([2 * F.SQ2 / np.pi * lam * np.cos(th), F.SQ2 * np.sin(th)], 1))
    for la in range(-90 + step, 90, step):
        t = F.moll_theta(np.array([math.radians(la)]))[0]
        xe = 2 * F.SQ2 * math.cos(t)
        segs.append(np.array([[-xe, F.SQ2 * math.sin(t)], [xe, F.SQ2 * math.sin(t)]]))
    ax.add_collection(LineCollection(segs, colors=[GRAT], linewidths=0.85, zorder=5))


def moll_points(P, lon0):
    """Projected marker positions; points exactly on the seam are drawn on both edges."""
    lat, lon = F.vec_to_latlon(P)
    x, y, lam = F.moll_forward(lat, lon, lon0)
    xs, ys = list(x), list(y)
    for i in range(len(x)):
        if abs(abs(lam[i]) - np.pi) < 1e-6 and abs(abs(lat[i]) - np.pi / 2) > 1e-6:
            xs.append(-x[i])
            ys.append(y[i])
    return np.array(xs), np.array(ys)


# ====================================================================== hero
class Hero:
    W = H = 1440
    D, GX = 570, 28
    VIEWS = ((0.0, 35.2643897), (180.0, -35.2643897))

    def __init__(self, model: Model, s: int, ss=2):
        self.m, self.s = model, s
        W, H, D = self.W, self.H, self.D
        x0 = (W - (2 * D + self.GX)) // 2
        self.globe_xy = {}
        tops = (14, 696)
        for row in (0, 1):
            for col in (0, 1):
                self.globe_xy[(row, col)] = (x0 + col * (D + self.GX), tops[row])
        self.x0, self.span = x0, 2 * D + self.GX
        t0 = time.time()
        self.grids = [F.DiscGrid(D * s, ss, *v) for v in self.VIEWS]
        self.vs = [F.GaussSampler(model.vort["lat_rad"], model.vort["lon_rad"], g.Q) for g in self.grids]
        self.ts = [F.TriSampler(model.tri, g.Q) for g in self.grids]
        for g in self.grids:
            del g.Q
        self.prep_s = time.time() - t0
        self.vscale = F.Scale("RdBu_r", "linear", model.zeta_lim)
        self.escale = F.Scale("magma", "asinh", model.E_max)
        self.fig = None

    def build(self):
        W, H, s, D = self.W, self.H, self.s, self.D
        with mpl.rc_context(rc()):
            fig = plt.figure(figsize=(W / 100, H / 100), dpi=100 * s)
        self.fig = fig
        self.ims = {}
        blank = np.full((D * s, D * s, 3), 255, np.uint8)
        with mpl.rc_context(rc()):
            for (row, col), (gx, gy) in self.globe_xy.items():
                lon0, lat0 = self.VIEWS[col]
                self.ims[(row, col)] = fig.figimage(blank, xo=gx * s, yo=(H - gy - D) * s, origin="upper", zorder=0)
                ax = fig.add_axes([gx / W, 1 - (gy + D) / H, D / W, D / H], zorder=2)
                ax.set_xlim(-1, 1)
                ax.set_ylim(-1, 1)
                ax.set_aspect("equal")
                ax.axis("off")
                graticule_ortho(ax, lon0, lat0)
                ax.add_patch(Ellipse((0, 0), 2, 2, fill=False, ec=EDGE, lw=1.0, zorder=6))
                x, y, d = F.ortho_forward(self.m.cones, lon0, lat0)
                vis = d > 0
                halo_marker(ax, x[vis], y[vis], "D", 12.5)
            vt, vl = vort_ticks(self.m.zeta_lim)
            colourbar(fig, (self.x0, 598, self.span, 16), self.vscale, vt, vl, "Relative vorticity (s$^{-1}$)", W, H, s)
            et, el = asinh_ticks(self.m.E_max)
            colourbar(fig, (self.x0, 1280, self.span, 16), self.escale, et, el, "Height error (m, asinh scale)", W, H, s)
            self.time = TimeAxis(fig, (self.x0 + 196, 1378, self.span - 196 - 170, 24), W, H, self.x0 + 172)
            h = marker_handle("D", 12.5)
            fig.legend([h], ["Cones"], handler_map={h: HandlerTuple(ndivide=1)}, loc="center left",
                       bbox_to_anchor=((self.x0 + self.span - 150) / W, 1 - 1390 / H), frameon=False,
                       fontsize=FS_LABEL, handletextpad=0.5, borderaxespad=0, borderpad=0)
        self.missing = False

    def render(self, k):
        m, s = self.m, self.s
        z = m.vort_get(k)
        E = m.epair("main", "T213", k)
        for col in (0, 1):
            g = self.grids[col]
            img, _ = F.compose(g, self.vs[col](z), self.vscale, (g.D, g.D))
            self.ims[(0, col)].set_data(img)
            img, miss = F.compose(g, self.ts[col](E), self.escale, (g.D, g.D))
            self.missing |= miss
            self.ims[(1, col)].set_data(img)
        self.time.set(m.t[k] / 24.0)
        self.fig.canvas.draw()
        return np.asarray(self.fig.canvas.buffer_rgba())[..., :3].copy()

    def regions(self):
        D, R = self.D, self.D / 2 - 3
        c = np.arange(D) + 0.5 - D / 2
        m = np.hypot(*np.meshgrid(c, c)) < R
        return {f"{['vorticity', 'error'][r]}-{['left', 'right'][c_]}": (x, y, D, D, m)
                for (r, c_), (x, y) in self.globe_xy.items()}

    whitelist = ["Relative vorticity (s^-1)", "Height error (m, asinh scale)", "Time (days)", "Cones", "Missing"]
    typeset = {"Relative vorticity (s^-1)": "Relative vorticity (s$^{-1}$)"}


# ====================================================================== relaxation fields
class FieldsMovie:
    W, H = 1600, 1400
    WM, HM, GX, GY = 724, 362, 16, 22
    X0, TOP = 124, 64
    LON0 = 60.0
    CBW = 760

    def __init__(self, model: Model, s: int, ss=2):
        self.m, self.s = model, s
        t0 = time.time()
        self.grid = F.MollGrid(self.WM * s, self.HM * s, ss, self.LON0)
        self.ts = F.TriSampler(model.tri, self.grid.Q)
        del self.grid.Q
        # W contour grids from the actual relax_weight (piecewise linear on the same triangulation)
        n = 2 * self.WM
        xs = np.linspace(-2 * F.SQ2 * 1.01, 2 * F.SQ2 * 1.01, n)
        ys = np.linspace(-F.SQ2 * 1.01, F.SQ2 * 1.01, n // 2)
        XX, YY = np.meshgrid(xs, ys)
        lat, lon = F.moll_inverse(XX.ravel(), YY.ravel(), self.LON0)
        sam = F.TriSampler(model.tri, F.latlon_to_vec(lat, lon))
        self.wgrid = {}
        for r in ("cones20", "controls20"):
            W_ = model.pairs_by_run[r].epair(model.runs[r]["relax_weight"])   # equal copies: returns W itself
            self.wgrid[r] = (XX, YY, sam(W_).reshape(XX.shape))
        self.prep_s = time.time() - t0
        self.escale = F.Scale("magma", "asinh", model.E_max)

    def box(self, row, col):
        return (self.X0 + col * (self.WM + self.GX), self.TOP + row * (self.HM + self.GY))

    def build(self):
        W, H, s = self.W, self.H, self.s
        with mpl.rc_context(rc()):
            fig = plt.figure(figsize=(W / 100, H / 100), dpi=100 * s)
        self.fig = fig
        self.ims = {}
        blank = np.full((self.HM * s, self.WM * s, 3), 255, np.uint8)
        a, b = 2 * F.SQ2, F.SQ2
        with mpl.rc_context(rc()):
            for row, run in enumerate(RUN_NAMES):
                for col, ref in enumerate(REFS):
                    bx, by = self.box(row, col)
                    self.ims[(run, ref)] = fig.figimage(blank, xo=bx * s, yo=(H - by - self.HM) * s, origin="upper",
                                                        zorder=0)
                    ax = fig.add_axes([bx / W, 1 - (by + self.HM) / H, self.WM / W, self.HM / H], zorder=2)
                    ax.set_xlim(-a, a)
                    ax.set_ylim(-b, b)
                    ax.axis("off")
                    graticule_moll(ax, self.LON0)
                    edge = Ellipse((0, 0), 2 * a, 2 * b, fill=False, ec=EDGE, lw=1.0, zorder=6)
                    ax.add_patch(edge)
                    clip = Ellipse((0, 0), 2 * a, 2 * b, transform=ax.transData)
                    if run != "main":
                        XX, YY, WG = self.wgrid[run]
                        ls = "solid" if run == "cones20" else (0, (5.0, 3.0))
                        draw_contour(ax, XX, YY, WG, math.exp(-1), ls, 1.7, clip)
                        draw_contour(ax, XX, YY, WG, 0.001, (0, (1.0, 2.4)), 1.1, clip)
                    x, y = moll_points(self.m.cones, self.LON0)
                    halo_marker(ax, x, y, "D", 12.5, z=20)
                    if run == "controls20":
                        x, y = moll_points(self.m.controls, self.LON0)
                        halo_marker(ax, x, y, "o", 11.5, z=21)
                lab_y = 1 - (self.box(row, 0)[1] + self.HM / 2) / H
                fig.text((self.X0 - 14) / W, lab_y, {"main": "Main", "cones20": "Cones", "controls20": "Control"}[run],
                         ha="right", va="center", fontsize=FS_LABEL)
            for col, ref in enumerate(REFS):
                bx, _ = self.box(0, col)
                fig.text((bx + self.WM / 2) / W, 1 - (self.TOP - 12) / H, ref, ha="center", va="bottom",
                         fontsize=FS_LABEL)
            et, el = asinh_ticks(self.m.E_max)
            cb_y = self.TOP + 3 * self.HM + 2 * self.GY + 22
            colourbar(fig, (self.X0, cb_y, self.CBW, 16), self.escale, et, el, "Height error (m, asinh scale)", W, H, s)
            halo = lambda w: [pe.Stroke(linewidth=w + 2.2, foreground="white"), pe.Normal()]  # noqa: E731
            hs = [marker_handle("D", 12.5), marker_handle("o", 11.5),
                  (Line2D([], [], color=INK, lw=1.7, ls="solid", path_effects=halo(1.7)),
                   Line2D([], [], color=INK, lw=1.7, ls=(0, (5.0, 3.0)), path_effects=halo(1.7))),
                  Line2D([], [], color=INK, lw=1.1, ls=(0, (1.0, 2.4)), path_effects=halo(1.1))]
            labs = ["Cone points", "Control points", "W = 0.368", "W = 0.001"]
            hm = {hs[0]: HandlerTuple(ndivide=1), hs[1]: HandlerTuple(ndivide=1), hs[2]: HandlerTuple(ndivide=None, pad=0.5)}
            fig.legend(hs, labs, handler_map=hm, ncol=2,
                       loc="upper right", bbox_to_anchor=((W - 14) / W, 1 - (cb_y - 8) / H),
                       frameon=False, fontsize=FS_TICK, handlelength=3.8, handletextpad=0.6, columnspacing=1.3,
                       labelspacing=0.7, borderaxespad=0, borderpad=0)
            ty = cb_y + 112
            self.time = TimeAxis(fig, (self.X0 + 196, ty, 1600 - 12 - (self.X0 + 196) - 24, 24), W, H, self.X0 + 172)
        self.missing = False

    def render(self, k):
        m = self.m
        g = self.grid
        for run in RUN_NAMES:
            for ref in REFS:
                img, miss = F.compose(g, self.ts(m.epair(run, ref, k)), self.escale, (g.Hm, g.Wm))
                self.missing |= miss
                self.ims[(run, ref)].set_data(img)
        self.time.set(m.t[k] / 24.0)
        self.fig.canvas.draw()
        return np.asarray(self.fig.canvas.buffer_rgba())[..., :3].copy()

    def regions(self):
        Wm, Hm = self.WM, self.HM
        X, Y = np.meshgrid((np.arange(Wm) + 0.5 - Wm / 2) / (Wm / 2 - 3), (np.arange(Hm) + 0.5 - Hm / 2) / (Hm / 2 - 3))
        m = X ** 2 + Y ** 2 < 1
        return {f"{run}-{ref}": (*self.box(r, c), Wm, Hm, m) for r, run in enumerate(RUN_NAMES)
                for c, ref in enumerate(REFS)}

    whitelist = ["Main", "Cones", "Control", "T85", "T213", "Cone points", "Control points", "W = 0.368",
                 "W = 0.001", "Height error (m, asinh scale)", "Time (days)", "Missing"]
    typeset = {}


VIS = {"hero": (Hero, "hero"), "fields": (FieldsMovie, "relaxation-fields")}


# ====================================================================== commands
def ffmpeg_exe():
    executable = shutil.which(FFMPEG)
    if executable is None:
        raise FileNotFoundError("Install FFmpeg with libx264, or set FFMPEG to its executable")
    return executable


def encode(seq_dir, out_mp4, bitrate, log_dir):
    exe = ffmpeg_exe()
    vf = "scale=out_color_matrix=bt709:out_range=tv:flags=accurate_rnd+full_chroma_int+lanczos,format=yuv420p"
    args = ["-framerate", "24", "-i", "%04d.png", "-an", "-vf", vf, "-c:v", "libx264", "-preset", "veryslow",
            "-profile:v", "high", "-b:v", bitrate, "-g", "96", "-color_primaries", "bt709", "-color_trc", "bt709",
            "-colorspace", "bt709", "-color_range", "tv"]
    run = [exe, "-v", "error", "-y"] + [str(seq_dir / a) if a == "%04d.png" else a for a in args] + [
        "-passlogfile", str(log_dir / "x264pass")]
    t0 = time.time()
    subprocess.run(run + ["-pass", "1", "-f", "mp4", os.devnull], check=True)
    tmp = out_mp4.with_suffix(".tmp.mp4")
    subprocess.run(run + ["-pass", "2", "-movflags", "+faststart", str(tmp)], check=True)
    os.replace(tmp, out_mp4)
    ver = subprocess.run([exe, "-version"], capture_output=True, text=True).stdout.split("\n")[0]
    return {"seconds": time.time() - t0, "ffmpeg": ver, "args": args + ["two-pass (-pass 1, -pass 2)",
                                                                       "-movflags", "+faststart"]}


def _lap(a):
    return (a[1:-1, :-2] + a[1:-1, 2:] + a[:-2, 1:-1] + a[2:, 1:-1] - 4 * a[1:-1, 1:-1])


def compression_check(mp4, fdir, kf, regions, frames, qa, stem, keep=()):
    """Decode chosen frames and compare each field region with the source render: PSNR (RGB), mean and 99.9th
    percentile absolute difference, and the ratio of decoded to source Laplacian energy of luminance (fine-structure
    retention; 1 means the encode kept the small-scale texture)."""
    from PIL import Image
    exe = ffmpeg_exe()
    res = {}
    for f in frames:
        dec_p = fdir / f"dec-{f:04d}.png"
        subprocess.run([exe, "-v", "error", "-y", "-i", str(mp4), "-vf",
                        f"select=eq(n\\,{f}),scale=in_color_matrix=bt709:in_range=tv:flags=accurate_rnd+full_chroma_int,"
                        "format=rgb24", "-fps_mode", "passthrough", "-frames:v", "1", str(dec_p)], check=True)
        src = np.asarray(Image.open(fdir / "seq" / f"{f:04d}.png").convert("RGB"), np.float64)
        dec = np.asarray(Image.open(dec_p).convert("RGB"), np.float64)
        row = {"snapshot": int(kf[f])}
        for name, (x, y, w, h, m) in regions.items():
            a, b = src[y:y + h, x:x + w][m], dec[y:y + h, x:x + w][m]
            mse = ((a - b) ** 2).mean()
            la = (src[y:y + h, x:x + w] @ [0.2126, 0.7152, 0.0722])
            lb = (dec[y:y + h, x:x + w] @ [0.2126, 0.7152, 0.0722])
            mm = m[1:-1, 1:-1]
            ea, eb = (_lap(la)[mm] ** 2).sum(), (_lap(lb)[mm] ** 2).sum()
            d = np.abs(a - b)
            row[name] = {"psnr_db": float(10 * np.log10(255 ** 2 / max(mse, 1e-12))), "mean_abs": float(d.mean()),
                         "p999_abs": float(np.percentile(d, 99.9)), "hf_energy_ratio": float(eb / max(ea, 1e-12))}
        res[str(f)] = row
        if f in keep:
            shutil.copy(fdir / "seq" / f"{f:04d}.png", qa / f"{stem}-source-f{f:04d}.png")
            shutil.copy(dec_p, qa / f"{stem}-decoded-f{f:04d}.png")
    return res


def do_check(model, qa):
    (qa / "movie-checks.json").write_text(json.dumps({"inputs": model.files, "checks": model.checks,
                                                      "summary": model.summary()}, indent=1, default=float))


def do_poster(model, which, out, qa, k=None, s=2, frame=False):
    cls, stem = VIS[which]
    k = model.k_day[10] if k is None else k
    if not frame:
        assert abs(model.t[k] - 240.0) < 1e-6, model.t[k]
        if which == "hero":
            assert model.vort_a[k] == 0.0 and abs(model.vort["t_hours"][model.vort_j[k]] - 240.0) < 1e-6
    v = cls(model, s)
    v.build()
    t1 = time.time()
    rgb = v.render(k)
    audit = audit_text(v.fig, v.whitelist, v.typeset)
    dest = (out / f"{stem}-poster.png") if not frame else (qa / f"{stem}-frame-k{k:03d}-s{s}.png")
    save_rgb_png(rgb, dest, dpi=300)
    info = {"file": dest.name, "snapshot_index": k, "t_hours": float(model.t[k]), "px": list(rgb.shape[:2][::-1]),
            "prep_s": v.prep_s, "render_s": time.time() - t1, "audit": audit, "missing_present": v.missing}
    (qa / f"{stem}-{'frame' if frame else 'poster'}.json").write_text(json.dumps(info, indent=1, default=float))
    plt.close(v.fig)
    print(json.dumps({k_: v_ for k_, v_ in info.items() if k_ != "audit"}, default=float),
          "audit ok" if audit["ok"] else audit["violations"], flush=True)
    return info


def do_movie(model, which, out, qa, bitrate):
    from PIL import Image
    cls, stem = VIS[which]
    v = cls(model, 1)
    v.build()
    audit = audit_text(v.fig, v.whitelist, v.typeset)
    if not audit["ok"]:
        raise SystemExit(f"text audit failed: {audit['violations']}")
    root = SCRATCH / "frames"
    root.mkdir(parents=True, exist_ok=True)
    fdir = Path(tempfile.mkdtemp(prefix=f"{stem}-", dir=root))
    (fdir / "u").mkdir()
    (fdir / "seq").mkdir()
    try:
        t1 = time.time()
        uniq = sorted(set(model.kf.tolist()))
        for i, k in enumerate(uniq):
            Image.fromarray(v.render(k)).save(fdir / "u" / f"{k:03d}.png", compress_level=1)
            if i % 20 == 0:
                print(f"{stem} frame {i + 1}/{len(uniq)} k={k} {time.time() - t1:.0f}s", flush=True)
        for f, k in enumerate(model.kf):
            os.link(fdir / "u" / f"{k:03d}.png", fdir / "seq" / f"{f:04d}.png")
        render_s = time.time() - t1
        mp4 = out / f"{stem}.mp4"
        enc = encode(fdir / "seq", mp4, bitrate, fdir)
        f7 = int(np.nonzero(model.kf == model.k_day[7])[0][0])
        f10 = int(np.nonzero(model.kf == model.k_day[10])[0][0])
        comp = compression_check(mp4, fdir, model.kf, v.regions(), [0, 150, f7, f10, 450, 575], qa, stem,
                                 keep=[f10, 575])
    finally:
        shutil.rmtree(fdir, ignore_errors=True)
    info = {"file": mp4.name, "bytes": mp4.stat().st_size, "sha256": sha256(mp4), "unique_frames": len(uniq),
            "render_s": render_s, "prep_s": v.prep_s, "encode": enc, "compression": comp, "audit": audit,
            "missing_present": v.missing, "day7_frame": int(np.nonzero(model.kf == model.k_day[7])[0][0]),
            "day10_frame": int(np.nonzero(model.kf == model.k_day[10])[0][0]),
            "frame_snapshot": model.kf.tolist()}
    (qa / f"{stem}-frames.json").write_text(json.dumps(info, indent=1, default=float))
    plt.close(v.fig)
    print(json.dumps({k_: v_ for k_, v_ in info.items() if k_ not in ("audit", "frame_snapshot")}, default=float),
          flush=True)
    return info


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["check", "poster", "posters", "movie", "frame", "all"])
    ap.add_argument("which", nargs="?", default="hero", choices=list(VIS))
    ap.add_argument("--data", default=str(DATA))
    ap.add_argument("--synthetic", action="store_true", help="generate stand-in data in a temporary directory")
    ap.add_argument("--out", default=str(OUT), help="deliverables (posters, MP4s)")
    ap.add_argument("--qa", default=str(QA), help="audit records")
    ap.add_argument("--bitrate", default="2900k")
    ap.add_argument("--k", type=int, default=None)
    ap.add_argument("--scale", type=int, default=None)
    a = ap.parse_args()
    out, qa = Path(a.out), Path(a.qa)
    out.mkdir(parents=True, exist_ok=True)
    qa.mkdir(parents=True, exist_ok=True)
    data = Path(a.data)
    if a.synthetic:
        import synth
        data = Path(tempfile.mkdtemp(prefix="synthetic-"))
        t0 = time.time()
        synth.main(data)
        print(f"synthetic data {time.time() - t0:.0f}s", flush=True)
    t0 = time.time()
    model = Model(data)
    print(f"model {time.time() - t0:.1f}s  E_max={model.E_max} (largest {model.E_largest:.4g})  "
          f"zeta_lim={model.zeta_lim} (largest {model.zeta_largest:.4g})", flush=True)
    do_check(model, qa)
    if a.cmd == "check":
        print(json.dumps(model.checks, indent=1, default=float)[:5000])
    elif a.cmd in ("poster", "frame"):
        do_poster(model, a.which, out, qa, a.k, a.scale or (2 if a.cmd == "poster" else 1), frame=a.cmd == "frame")
    elif a.cmd == "posters":
        do_poster(model, "hero", out, qa)
    elif a.cmd == "movie":
        do_movie(model, a.which, out, qa, a.bitrate)
    else:
        do_poster(model, "hero", out, qa)
        for w in ("hero", "fields"):
            do_movie(model, w, out, qa, a.bitrate)
    print(f"done {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
