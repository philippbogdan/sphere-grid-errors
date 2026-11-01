"""Run one archived relaxation configuration with two-hourly display snapshots."""

import argparse
import json
import math
from pathlib import Path
import time
from evidence import CONFIGURATIONS, REVISION, sha256

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("--region", choices=["none", "bp", "disc", "annulus", "polar"], default="none",
               help="Correction weight: none, Gaussian, smooth disc, smooth annulus or polar caps.")
p.add_argument("--par", type=float, default=0, help="Angular scale in degrees; use an archived configuration.")
p.add_argument("--rot", type=float, default=0, help="Zero for cones, 60 for the control centres.")
p.add_argument("--M", type=int, default=128, help="Even torus width; the archived intervention uses 128.")
p.add_argument("--t85", required=True, help="T85 reference binary: 45 s steps and two-hourly output.")
p.add_argument("--t213", default="", help="Optional T213 binary: 22.5 s steps and two-hourly output.")
p.add_argument("--days", default="1,3,5,7,10,14", choices=["1,3,5,7,10", "1,3,5,7,10,14"],
               help="Preserve archived stepping boundaries, with an optional day-14 extension.")
p.add_argument("--out", required=True, help="Output stem for the NPZ fields and JSON scores.")
a = p.parse_args()
if (a.region, a.par, a.rot) not in CONFIGURATIONS:
    p.error("Choose one of the 12 combinations listed in docs/METHOD.md")
if a.M < 16 or a.M % 2:
    p.error("M must be even and at least 16")
Path(a.out).parent.mkdir(parents=True, exist_ok=True)
if any(Path(a.out + suffix).exists() for suffix in (".json", ".npz")):
    p.error("The output stem already exists; choose a new one")
from runtime import configure, reference_read
configure()
import numpy as np
from polarswe import grid, testcases as tc, reference as ref, lwgg, configs, experiments as ex

VERTS = np.array([[0, 0, 1.0]] + [[math.cos(math.radians(-19.4712206)) * math.cos(math.radians(l)),
                                    math.cos(math.radians(-19.4712206)) * math.sin(math.radians(l)),
                                    math.sin(math.radians(-19.4712206))] for l in (0, 120, -120)])

nlon, nlat, g, recs = reference_read(a.t85)
times = np.array([r[0] for r in recs]) * 3600.0
R213 = None
if a.t213:
    _, _, g2, recs2 = reference_read(a.t213); R213 = {round(r[0], 6): r for r in recs2}
    if recs2[-1][0] < 24 * int(a.days.split(",")[-1]):
        raise ValueError("The T213 reference does not cover the requested duration")
if times[-1] < 24 * int(a.days.split(",")[-1]) * 3600:
    raise ValueError("The T85 reference does not cover the requested duration")
M = a.M; X = grid.torus_X(M); lat, _ = tc.latlon(X); latd = np.degrees(lat); e, n = tc.east_north(X)
X_, _ = ex.make_grid(("torus", M), ex.scheme_of(configs.MAIN).get("orient"))
F, S0, h, cmax, dt = ex.build(X_, tc.TC6(), configs.MAIN)
w = S0.weight
V = VERTS
if a.rot:
    c_, s_ = math.cos(math.radians(a.rot)), math.sin(math.radians(a.rot))
    Rz = np.array([[c_, -s_, 0], [s_, c_, 0], [0, 0, 1.0]])
    V = np.array([[0, 0, -1.0]] + [Rz @ v_ for v_ in VERTS[1:]])
dist = np.degrees(np.arccos(np.clip(np.tensordot(V, X, axes=1), -1, 1))).min(0)
par = a.par
W = {"none": lambda: np.zeros_like(dist), "bp": lambda: np.exp(-(dist / max(par, 1e-9)) ** 2),
     "disc": lambda: 0.5 * (1 - np.tanh((dist - par) / 2.0)),
     "annulus": lambda: 0.25 * (1 + np.tanh((dist - par) / 2.0)) * (1 - np.tanh((dist - 2 * par) / 2.0)),
     "polar": lambda: 0.5 * (1 + np.tanh((np.abs(latd) - par) / 3.0))}[a.region]()
sel = W > 1e-3; Ws = W[sel]; tau = 3600.0
cache = {}
def truth_at(t):
    """Interpolate the target at the active sites, using the original two-hour indexing."""
    j = min(int(t // 7200.0), len(times) - 2); al = (t - times[j]) / (times[j + 1] - times[j])
    for i_ in (j, j + 1):
        if i_ not in cache:
            _, u, vv, hh = recs[i_]
            cache[i_] = (ref.SphereInterp(g, hh)(X[:, sel]) - S0.hs[sel],
                         ref.SphereInterp(g, u, vector=True)(X[:, sel]) * e[:, sel] + ref.SphereInterp(g, vv, vector=True)(X[:, sel]) * n[:, sel])
    for i_ in [i_ for i_ in cache if i_ < j]:
        del cache[i_]
    return (1 - al) * cache[j][0] + al * cache[j + 1][0], (1 - al) * cache[j][1] + al * cache[j + 1][1]

def href_t(rs, gg, t_h):
    """Interpolate reference height in space and between the original saved times."""

    j = int(t_h // 2.0); al = (t_h - 2.0 * j) / 2.0
    h0 = ref.SphereInterp(gg, rs[j][3])(X)
    if al < 1e-9 or j + 1 >= len(rs):
        return h0
    return (1 - al) * h0 + al * ref.SphereInterp(gg, rs[j + 1][3])(X)

recs2_list = sorted(R213.values(), key=lambda r: r[0]) if R213 else None
D = np.ascontiguousarray(tc.TC6().height(X) - S0.hs); v = np.ascontiguousarray(tc.TC6().velocity(X))
t = 0.0; res = {}; res213 = {}
snaps_t, snaps85, snaps213, snapsh = [], [], [], []
next_snap = 0.0
def snap(t):
    """Save actual time and float32 fields; endpoint scores remain float64."""
    th = t / 3600.0
    snaps_t.append(th)
    snapsh.append((D + S0.hs).astype(np.float32))
    snaps85.append((D + S0.hs - href_t(recs, g, th)).astype(np.float32))
    if recs2_list:
        snaps213.append((D + S0.hs - href_t(recs2_list, g2, th)).astype(np.float32))
snap(0.0); next_snap = 2 * 3600.0
segments = []
t0 = time.time()
for day in [int(x) for x in a.days.split(",")]:
    T = day * tc.DAY; nn = int(np.ceil((T - t) / dt)); hh = (T - t) / nn
    segments.append(dict(start_s=t, end_s=T, nsteps=nn, step_s=hh))
    for i_ in range(nn):
        D, v = F.step(D, v, hh); t += hh
        if sel.any():
            Dt, vt = truth_at(t); c = (hh / tau) * Ws
            D[sel] += c * (Dt - D[sel]); v[:, sel] += c * (vt - v[:, sel])
        if t >= next_snap - 1e-6 and t < T - 1e-6:
            snap(t); next_snap += 2 * 3600.0
    t = T
    hr = ref.SphereInterp(g, [r for r in recs if abs(r[0] - 24 * day) < 1e-6][0][3])(X)
    res[str(day)] = float(lwgg.norms(D + S0.hs, hr, w)[1])
    if R213:
        res213[str(day)] = float(lwgg.norms(D + S0.hs, ref.SphereInterp(g2, R213[round(24.0 * day, 6)][3])(X), w)[1])
    snap(t); next_snap = t + 2 * 3600.0
    print(day, res[str(day)], res213.get(str(day)), f"{time.time() - t0:.0f}s", flush=True)
rec = dict(mode="relax", M=M, region=a.region, par=a.par, control_rotation_deg=a.rot,
           weighted_area=float((w * W).sum() / w.sum()), l2_vs_T85_by_day=res, l2_vs_T213dt22_by_day=res213,
           dt_s=dt, wall_s=time.time() - t0, scheme=configs.MAIN, public_solver_revision=REVISION,
           references={"T85dt45": {"sha256": sha256(a.t85)},
                       "T213dt22": {"sha256": sha256(a.t213)} if a.t213 else None},
           stepping_segments=segments, snapshot_rule="First step crossing each two-hour threshold; exact day endpoints")
json.dump(rec, open(a.out + ".json", "w"), indent=1)
np.savez_compressed(a.out + ".npz", t_hours=np.array(snaps_t), h=np.array(snapsh), err_vs_T85=np.array(snaps85),
                    err_vs_T213=np.array(snaps213) if snaps213 else np.zeros(0), X=X.astype(np.float32),
                    weight=w.astype(np.float32), relax_weight=W.astype(np.float32), centres=V)
print(json.dumps(rec))
