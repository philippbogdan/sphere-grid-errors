"""Generate synthetic fixtures for renderer tests; these are not study results."""
from __future__ import annotations

import json
import math
import sys
import tempfile
from pathlib import Path

import numpy as np
from scipy.spatial import SphericalVoronoi

sys.path.insert(0, str(Path(__file__).parent))
from common import SCRATCH  # noqa: E402

from polarswe import grid  # noqa: E402

grid.CACHE = str(SCRATCH / "solver-cache") if SCRATCH.exists() else tempfile.mkdtemp(prefix="solver-cache-")
# Keep generated caches separate from the evidence.

LATC = -19.4712206


def centres(control=False):
    cl = math.cos(math.radians(LATC))
    sl = math.sin(math.radians(LATC))
    if not control:
        return np.array([[0, 0, 1.0]] + [[cl * math.cos(math.radians(l)), cl * math.sin(math.radians(l)), sl]
                                         for l in (0, 120, -120)])
    return np.array([[0, 0, -1.0]] + [[cl * math.cos(math.radians(l)), cl * math.sin(math.radians(l)), sl]
                                      for l in (60, 180, -60)])


def times(rng):
    t = [0.0]
    ends = {24, 72, 120, 168, 240, 336}
    for k in range(1, 169):
        nom = 2.0 * k
        t.append(nom if nom in ends else nom + rng.uniform(0.01, 0.2))
    return np.array(t)


def smooth_field(P, rng, kmax, n=160, amp_decay=1.0):
    """Random superposition of plane waves restricted to the sphere: smooth, signed, O(1)."""
    k = rng.normal(size=(n, 3))
    k /= np.linalg.norm(k, axis=1)[:, None]
    mag = rng.uniform(2, kmax, size=n)
    ph = rng.uniform(0, 2 * np.pi, size=n)
    a = rng.normal(size=n) * mag ** (-amp_decay)
    return k * mag[:, None], ph, a


def eval_waves(P, waves, shift=0.0):
    k, ph, a = waves
    return (a[None] * np.cos(P @ k.T + ph[None] + shift)).sum(1) / np.sqrt((a ** 2).sum())


def main(out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(20261101)
    M = 128
    X = grid.torus_X(M)
    P = X.reshape(3, -1).T
    Mi, Mj = np.meshgrid(np.arange(M), np.arange(M), indexing="ij")
    flat = (Mi * M + Mj).ravel()
    pflat = (((-1 - Mi) % M) * M + ((-1 - Mj) % M)).ravel()
    rep = flat < pflat
    sv = SphericalVoronoi(P[flat[rep]], radius=1.0)
    area = sv.calculate_areas()
    w = np.zeros(M * M)
    w[flat[rep]] = area / 2
    w[pflat[rep]] = area / 2
    w = w.reshape(M, M)
    t = times(rng)
    lat = np.arcsin(P[:, 2])
    h0 = 8000 + 2500 * np.cos(lat) ** 2
    waves85 = smooth_field(P, rng, 40)
    waves213 = smooth_field(P, rng, 40)
    asym = smooth_field(P, rng, 20)
    C = centres()
    Cc = centres(True)
    dC = np.degrees(np.arccos(np.clip(P @ C.T, -1, 1)))
    dK = np.degrees(np.arccos(np.clip(P @ Cc.T, -1, 1)))
    bump = np.exp(-(dC.min(1) / 9.0) ** 2)
    for name, cen, d in (("main", C, None), ("cones20", C, dC), ("controls20", Cc, dK)):
        W = np.zeros(M * M) if d is None else np.exp(-(d.min(1) / 20.0) ** 2)
        damp = 1 - 0.8 * W
        E85, E213, H = [], [], []
        for th in t:
            day = th / 24
            A = 0.02 * math.exp(day / 1.75)
            s = 0.25 * day
            base = (0.35 * eval_waves(P, waves85, s) + 3.0 * bump * np.cos(25 * P[:, 0] + s)) * A * damp
            base213 = base + 0.15 * A * eval_waves(P, waves213, s)
            pair_noise = 1e-3 * A * eval_waves(P, asym, s)
            e85 = base.copy()
            e213 = base213.copy()
            e85[pflat[rep]] = base[flat[rep]] + pair_noise[flat[rep]]
            e213[pflat[rep]] = base213[flat[rep]] + pair_noise[flat[rep]]
            E85.append(e85.reshape(M, M))
            E213.append(e213.reshape(M, M))
            H.append((h0 + e85).reshape(M, M))
        np.savez_compressed(out / f"{name}.npz", t_hours=t, X=X.astype(np.float32), weight=w.astype(np.float32),
                            err_vs_T85=np.array(E85, np.float32), err_vs_T213=np.array(E213, np.float32),
                            h=np.array(H, np.float32), relax_weight=W.reshape(M, M).astype(np.float32),
                            centres=cen)
        side = {"mode": "relax", "M": M, "region": "none" if d is None else "bp", "par": 0 if d is None else 20,
                "control_rotation_deg": 60 if name == "controls20" else 0,
                "weighted_area": float((w.ravel() * W).sum() / w.sum()),
                "l2_vs_T85_by_day": {}, "l2_vs_T213dt22_by_day": {}, "dt_s": 740.0, "SYNTHETIC": True}
        (out / f"{name}.json").write_text(json.dumps(side))
        print(name, flush=True)

    # T213-shaped vorticity: analytic RH wave advected at its phase speed plus growing fine structure
    nlat, nlon = 320, 640
    x, _ = np.polynomial.legendre.leggauss(nlat)
    phi = np.arcsin(np.sort(x))
    lam = 2 * np.pi * np.arange(nlon) / nlon
    PH, LA = np.meshgrid(phi, lam, indexing="ij")
    G = np.stack([np.cos(PH) * np.cos(LA), np.cos(PH) * np.sin(LA), np.sin(PH)], -1).reshape(-1, 3)
    w_ = K = 7.848e-6
    R = 4
    nu = (R * (3 + R) * w_ - 2 * 7.292e-5) / ((1 + R) * (2 + R))
    fine = smooth_field(G, rng, 160, n=200, amp_decay=0.6)
    tref = np.arange(169) * 2.0
    Z = []
    for th in tref:
        ts = th * 3600
        z = 2 * w_ * np.sin(PH) - K * np.sin(PH) * np.cos(PH) ** R * (R * R + 3 * R + 2) * np.cos(R * (LA - nu * ts))
        g = min(1.0, (th / 24 / 9.0) ** 2)
        z = z + g * 6e-5 * eval_waves(G, fine, th / 24 * 0.9).reshape(nlat, nlon) * np.cos(PH)
        Z.append(z.astype(np.float32))
    Z = np.array(Z)
    for a, b in ((0, 56), (57, 112), (113, 168)):
        np.savez_compressed(out / f"t213-vorticity-{a:03d}-{b:03d}.npz", t_hours=tref[a:b + 1], lat_rad=phi,
                            lon_rad=lam, zeta_s_inv=Z[a:b + 1])
    print("vorticity", np.abs(Z).max())


if __name__ == "__main__":
    main(sys.argv[1])
