"""Prepare display fields, regional scores and vorticity from local numerical outputs."""

import argparse
import csv
import json
import math
import os
from pathlib import Path
from evidence import EVIDENCE, ROOT, REVISION, sha256

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("--runs", type=Path, default=ROOT / "build/runs",
               help="Directory containing none, cones and controls NPZ/JSON pairs.")
p.add_argument("--t85", type=Path, default=ROOT / "build/stswm/runs/T85dt45/grid_out.bin",
               help="T85 reference with 45 s steps and two-hourly output.")
p.add_argument("--t213", type=Path, default=ROOT / "build/stswm/runs/T213dt22/grid_out.bin",
               help="T213 reference with 22.5 s steps and two-hourly output.")
p.add_argument("--out", type=Path, default=ROOT / "build/fields", help="Destination for the local render pack.")
a = p.parse_args()
from runtime import configure, reference_read
configure()
import numpy as np
from polarswe import reference as ref, testcases as tc
import vort
RUNS, T85, T213, OUT = a.runs, a.t85, a.t213, a.out
OUT.mkdir(parents=True, exist_ok=False)

NAMES = {"main": "none", "cones20": "cones", "controls20": "controls"}
sha = sha256
val = {"inputs": {}, "runs": {}}

def ref_series(fn):
    """Keep the reference archive on disk and expose its original record values."""
    nlon, nlat, g, recs = reference_read(fn)
    return g, recs

def href_at(g, recs, th, X, cache):
    """Interpolate reference height in space and between adjacent two-hour records."""
    if th < recs[0][0] - 1e-9 or th > recs[-1][0] + 1e-9:
        raise ValueError("Field time is outside the reference interval")
    j = int(math.floor(th / 2.0 + 1e-9)); al = th / 2.0 - j
    def H(i):
        if i not in cache:
            cache.clear() if len(cache) > 4 else None
            cache[i] = ref.SphereInterp(g, recs[i][3])(X)
        return cache[i]
    assert abs(recs[j][0] - 2.0 * j) < 1e-6
    if al < 1e-9: return H(j), (j, j, 0.0)
    return (1 - al) * H(j) + al * H(j + 1), (j, j + 1, float(al))

g213, r213 = ref_series(T213); g85, r85 = ref_series(T85)
if len(r213) != 169 or len(r85) != 169:
    raise ValueError("Expected 169 reference samples from day zero to day 14")
if r213[0][3].shape != (320, 640) or r85[0][3].shape != (128, 256):
    raise ValueError("The reference grids differ from T213 and T85")
val["inputs"]["T213dt22_tc6_2h"] = dict(sha256=sha(T213), records=len(r213), t_first=r213[0][0], t_last=r213[-1][0], shape=list(r213[0][3].shape))
val["inputs"]["T85dt45_tc6_2h"] = dict(sha256=sha(T85), records=len(r85), t_first=r85[0][0], t_last=r85[-1][0], shape=list(r85[0][3].shape))

rows = []
dist_cache = {}
for name, tag in NAMES.items():
    z = np.load(os.path.join(RUNS, f"{tag}.npz")); meta = json.load(open(os.path.join(RUNS, f"{tag}.json")))
    X = z["X"].astype(np.float64); w = z["weight"].astype(np.float64); t = z["t_hours"]; h = z["h"].astype(np.float64)
    M = X.shape[1]
    if M != 128 or t[0] != 0 or abs(t[-1] - 336) > 1e-6 or not np.all(np.diff(t) > 0):
        raise ValueError("The visual pack requires M=128 and increasing times spanning days zero to 14")
    val["inputs"][name] = {"field_sha256": sha(RUNS / f"{tag}.npz"),
                           "sidecar_sha256": sha(RUNS / f"{tag}.json")}

    # Reconstruct errors from float32 saved heights exactly as in the display pack.
    # Full-precision endpoint scores remain in the JSON sidecar.
    e85 = np.empty(h.shape, np.float32); e213 = np.empty(h.shape, np.float32); al213 = []
    c85, c213 = {}, {}
    hr85 = []; hr213 = []
    for k, th in enumerate(t):
        a, _ = href_at(g85, r85, th, X, c85); b, s = href_at(g213, r213, th, X, c213)
        e85[k] = h[k] - a; e213[k] = h[k] - b; al213.append(s)
        hr85.append(a.astype(np.float32)); hr213.append(b.astype(np.float32))
    np.savez_compressed(os.path.join(OUT, f"{name}.npz"), t_hours=t, X=z["X"], weight=z["weight"], h=z["h"],
                        err_vs_T85=e85, err_vs_T213=e213, relax_weight=z["relax_weight"], centres=z["centres"])
    side = dict(meta); side["public_solver_revision"] = REVISION; side.update(run=name, reference_T85="STSWM T85, 45 s step, HDC 3.1106E+14, 2-hourly output",
                                   reference_T213="STSWM T213, 22.5 s step, HDC 8.0E+12, 2-hourly output",
                                   snapshot_rule="First solver step crossing each two-hour threshold; actual hours retained; day endpoints exact",
                                   reference_time_alignment="Linear between the two bracketing reference records; no extrapolation")
    json.dump(side, open(os.path.join(OUT, f"{name}.json"), "w"), indent=1)

    I, J = np.meshgrid(np.arange(M), np.arange(M), indexing="ij"); Ip, Jp = (-1 - I) % M, (-1 - J) % M
    rep = (I * M + J) < (Ip * M + Jp)
    pos = np.linalg.norm(X[:, rep] - X[:, Ip[rep], Jp[rep]], axis=0).max()
    wdiff = np.abs(w[rep] - w[Ip[rep], Jp[rep]]).max()
    vr = dict(position_max_diff=float(pos), weight_max_abs_diff=float(wdiff), n_pairs=int(rep.sum()),
              t_hours_first=float(t[0]), t_hours_last=float(t[-1]), n_snapshots=int(len(t)),
              max_offset_from_2h_grid_h=float(np.abs(t - 2 * np.round(t / 2)).max()))

    for refname, E, HR in (("T85dt45", e85, hr85), ("T213dt22", e213, hr213)):
        pd = []
        for k, th in enumerate(t):
            e = E[k].astype(np.float64); hr = HR[k].astype(np.float64)
            ei, ep = e[rep], e[Ip[rep], Jp[rep]]; wi, wp = w[rep], w[Ip[rep], Jp[rep]]
            prms = math.sqrt(((wi + wp) * (ei - ep) ** 2).sum() / (wi + wp).sum()); pmax = float(np.abs(ei - ep).max())
            pd.append(pmax)
            for region in ("global", "cones20", "controls20"):
                if region == "global": q = np.ones_like(w)
                else:
                    src = "cones" if region == "cones20" else "controls"
                    if src not in dist_cache:
                        zz = np.load(os.path.join(RUNS, f"{src}.npz")); C = zz["centres"]
                        d = np.degrees(np.arccos(np.clip(np.tensordot(C, X, axes=1), -1, 1))).min(0)
                        dist_cache[src] = np.exp(-(d / 20.0) ** 2)
                    q = dist_cache[src]
                sw = w.sum(); area = (w * q).sum() / sw; wes = (w * q * e * e).sum() / sw; rhs = (w * q * hr * hr).sum() / sw
                tot = (w * e * e).sum() / sw
                share = wes / tot if tot > 0 else float("nan")
                rows.append(dict(run=name, reference=refname, t_hours=float(th), region=region, area_share=area,
                                 weighted_error_sq_m2=wes, reference_height_sq_m2=rhs, l2=math.sqrt(wes / rhs) if rhs > 0 else float("nan"),
                                 rms_m=math.sqrt(wes / area) if area > 0 else float("nan"), error_energy_share=share,
                                 enrichment=share / area if area > 0 else float("nan"),
                                 pair_rms_difference_m=prms, pair_max_difference_m=pmax))
        vr[f"pair_max_abs_error_diff_m_{refname}"] = float(max(pd))

    vr["reference_alignment_T213"] = al213
    vr["new_l2_vs_T85_by_day"] = meta["l2_vs_T85_by_day"]
    arch = [json.loads(l) for l in open(EVIDENCE / "tc6_week2_relax.jsonl")]
    key = {"none": ("none", 0, 0), "cones": ("bp", 20, 0), "controls": ("bp", 20, 60)}[tag]
    a = [r for r in arch if (r["region"], r["par"], r["control_rotation_deg"]) == key]; assert len(a) == 1
    vr["archived_l2_vs_T85_by_day"] = a[0]["l2_vs_T85_by_day"]
    vr["relative_difference_new_vs_archived"] = {d: meta["l2_vs_T85_by_day"][d] / a[0]["l2_vs_T85_by_day"][d] - 1 for d in a[0]["l2_vs_T85_by_day"]}
    val["runs"][name] = vr
    print(name, "done", flush=True)

import csv
cols = "run,reference,t_hours,region,area_share,weighted_error_sq_m2,reference_height_sq_m2,l2,rms_m,error_energy_share,enrichment,pair_rms_difference_m,pair_max_difference_m".split(",")
with open(os.path.join(OUT, "field-metrics.csv"), "w", newline="") as f:
    wr = csv.DictWriter(f, fieldnames=cols); wr.writeheader()
    for r in rows: wr.writerow({k: (f"{r[k]:.10g}" if isinstance(r[k], float) else r[k]) for k in cols})

o = np.argsort(g213); phi = g213[o]; nlon = r213[0][1].shape[1]
lam = 2 * np.pi * np.arange(nlon) / nlon
n = len(r213); chunks = [(0, 56), (57, 112), (113, n - 1)]
zmax = 0.0; chk = {}
for a_, b_ in chunks:
    Z = []; T = []
    for i in range(a_, b_ + 1):
        t_, u, v, _ = r213[i]
        z = vort.zeta(u[o], v[o], phi); Z.append(z.astype(np.float32)); T.append(t_)
        zmax = max(zmax, float(np.abs(z).max()))
        if abs(t_ - 24 * round(t_ / 24)) < 1e-6 and round(t_ / 24) in (0, 7, 10, 14):
            z2 = vort.zeta_fd(u[o], v[o], phi); dd = z2 - z
            chk[f"day{round(t_ / 24)}"] = dict(rms_diff_vs_fd2=float(np.sqrt((dd ** 2).mean())), max_diff_vs_fd2=float(np.abs(dd).max()),
                                              polar_max_diff_vs_fd2=float(np.abs(dd[np.abs(phi) > np.radians(80)]).max()), max_abs_zeta=float(np.abs(z).max()))
    fn = os.path.join(OUT, f"t213-vorticity-{a_:03d}-{b_:03d}.npz")
    np.savez_compressed(fn, t_hours=np.array(T), lat_rad=phi, lon_rad=lam, zeta_s_inv=np.array(Z))
    if os.path.getsize(fn) >= 50_000_000:
        raise ValueError("Vorticity chunk exceeds 50 MB")
    print(fn, os.path.getsize(fn) / 1e6, "MB", flush=True)

P, L = np.meshgrid(phi, lam, indexing="ij"); A = vort.A; w_ = K = 7.848e-6; R = 4
u0 = A * w_ * np.cos(P) + A * K * np.cos(P) ** (R - 1) * (R * np.sin(P) ** 2 - np.cos(P) ** 2) * np.cos(R * L)
v0 = -A * K * R * np.cos(P) ** (R - 1) * np.sin(P) * np.sin(R * L)
zt = 2 * w_ * np.sin(P) - K * np.sin(P) * np.cos(P) ** R * (R ** 2 + 3 * R + 2) * np.cos(R * L)
d = vort.zeta(u0, v0, phi) - zt
z0 = vort.zeta(r213[0][1][o], r213[0][2][o], phi) - zt
val["vorticity"] = dict(method="Longitude: periodic Fourier derivative. Latitude: cubic spline of u*cos(phi), with eight rows continued across each pole and shifted by half a longitude turn. Radius from the pinned solver.",
                        analytic_tc6_rms_err=float(np.sqrt((d ** 2).mean())), analytic_tc6_max_err=float(np.abs(d).max()),
                        analytic_tc6_polar_max_err=float(np.abs(d[np.abs(phi) > np.radians(80)]).max()),
                        stswm_day0_vs_analytic_rms=float(np.sqrt((z0 ** 2).mean())), stswm_day0_vs_analytic_max=float(np.abs(z0).max()),
                        second_derivative_checks=chk, max_abs_zeta_all_frames=zmax, chunks=[list(c) for c in chunks])
json.dump(val, open(os.path.join(OUT, "validation.json"), "w"), indent=1)
print("pack done")

assets = {}
for path in sorted(OUT.glob("*")):
    if path.is_file():
        assets[path.name] = {"bytes": path.stat().st_size, "sha256": sha256(path)}
(OUT / "assets.json").write_text(json.dumps({"schema": 1, "public_solver_revision": REVISION,
    "assets": assets}, indent=2) + "\n")
