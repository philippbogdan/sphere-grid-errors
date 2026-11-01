"""Fast checks of the field geometry. python test_fields.py"""
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import fields as F  # noqa: E402

rng = np.random.default_rng(1)
res = {}

# Mollweide round trip, including the seam at lon0 +- 180
lat = rng.uniform(-1.5, 1.5, 20000)
lon = rng.uniform(-np.pi, np.pi, 20000)
x, y, lam = F.moll_forward(lat, lon, 60.0)
la2, lo2 = F.moll_inverse(x, y, 60.0)
d = np.linalg.norm(F.latlon_to_vec(lat, lon) - F.latlon_to_vec(la2, lo2), axis=1)
res["moll_roundtrip_max_rad"] = float(d.max())
assert d.max() < 1e-12
# a point exactly on the seam (lon -120 when centred on 60E) projects to the boundary
xs, ys, ls = F.moll_forward(np.array([math.radians(-19.4712206)]), np.array([math.radians(-120.0)]), 60.0)
g = (xs / (2 * F.SQ2)) ** 2 + (ys / F.SQ2) ** 2
res["seam_point_on_boundary"] = float(g[0])
assert abs(g[0] - 1) < 1e-12 and abs(abs(ls[0]) - np.pi) < 1e-12

# orthographic: inverse then forward
xx = rng.uniform(-1, 1, 20000)
yy = rng.uniform(-1, 1, 20000)
ok = xx * xx + yy * yy < 1
for lon0, lat0 in ((0, 35.2643897), (180, -35.2643897)):
    P = F.ortho_inverse(xx[ok], yy[ok], lon0, lat0)
    x2, y2, dep = F.ortho_forward(P, lon0, lat0)
    assert np.abs(x2 - xx[ok]).max() < 1e-12 and np.abs(y2 - yy[ok]).max() < 1e-12 and (dep >= 0).all()
    # geographic north is up: a point slightly north of the centre projects above it
    c = F.latlon_to_vec(np.radians([lat0 + 1.0]), np.radians([lon0]))
    xn, yn, _ = F.ortho_forward(c, lon0, lat0)
    assert yn[0] > 0 and abs(xn[0]) < 1e-12
    # east is right
    c = F.latlon_to_vec(np.radians([lat0]), np.radians([lon0 + 1.0]))
    xe, ye, _ = F.ortho_forward(c, lon0, lat0)
    assert xe[0] > 0
res["ortho_roundtrip"] = "ok"

# cone visibility: each view shows two cones, at 54.7 degrees from the view centre
LATC = -19.4712206
cones = F.latlon_to_vec(np.radians([90, LATC, LATC, LATC]), np.radians([0, 0, 120, -120]))
for lon0, lat0 in ((0, 35.2643897), (180, -35.2643897)):
    _, _, dep = F.ortho_forward(cones, lon0, lat0)
    res[f"cone_depths_view_{lon0}"] = dep.round(6).tolist()
    assert (dep > 0.57).sum() == 2

# Gaussian-grid sampler: smooth fields, pole consistency
nlat, nlon = 320, 640
xg, _ = np.polynomial.legendre.leggauss(nlat)
phi = np.arcsin(np.sort(xg))
lamg = 2 * np.pi * np.arange(nlon) / nlon
PH, LA = np.meshgrid(phi, lamg, indexing="ij")
G = np.stack([np.cos(PH) * np.cos(LA), np.cos(PH) * np.sin(LA), np.sin(PH)], -1)
a = np.array([0.3, -0.7, 0.5])
f_grid = (G @ a) + 0.2 * G[..., 0] * G[..., 2]
Q = rng.normal(size=(200000, 3))
Q /= np.linalg.norm(Q, axis=1)[:, None]
# add points very near both poles, all longitudes
pl = np.radians(89.9)
lonr = rng.uniform(0, 2 * np.pi, 2000)
Q = np.concatenate([Q, F.latlon_to_vec(np.full(2000, pl), lonr), F.latlon_to_vec(np.full(2000, -pl), lonr),
                    np.array([[0, 0, 1.0], [0, 0, -1.0]])])
gs = F.GaussSampler(phi, lamg, Q)
f_true = Q @ a + 0.2 * Q[:, 0] * Q[:, 2]
err = np.abs(gs(f_grid) - f_true)
res["gauss_sampler_max_err"] = float(err.max())
res["gauss_sampler_cap_counts"] = gs.counts
assert err.max() < 2e-4, err.max()
# value at the pole equals the ring mean regardless of longitude
NP = gs(f_grid)[-2]
res["pole_value_minus_ring_mean"] = float(NP - f_grid[-1].mean())
assert abs(NP - f_grid[-1].mean()) < 1e-6

# triangulation: linear field reproduced to second order; every point inside a facet
n = 8192
i = np.arange(n) + 0.5
U = np.stack([np.cos(np.pi * (1 + 5 ** 0.5) * i) * np.sqrt(1 - (1 - 2 * i / n) ** 2),
              np.sin(np.pi * (1 + 5 ** 0.5) * i) * np.sqrt(1 - (1 - 2 * i / n) ** 2), 1 - 2 * i / n], 1)
tri = F.SphereTri(U)
res["tri_check"] = tri.check
ts = F.TriSampler(tri, Q)
e = np.abs(ts(U @ a) - Q @ a)
res["tri_linear_max_err"] = float(e.max())
assert e.max() < 2e-3
# values never exceed the node range (no invented extrema)
v = rng.uniform(0, 1, n)
out = ts(v)
assert out.min() >= v.min() - 1e-6 and out.max() <= v.max() + 1e-6

# schedule
t = [0.0] + [2.0 * k + (0 if 2 * k in (24, 72, 120, 168, 240, 336) else 0.1) for k in range(1, 169)]
tau, k = F.display_schedule(np.array(t))
res["schedule"] = {"frames": int(k.size), "unique": int(np.unique(k).size), "first": int(k[0]), "last": int(k[-1]),
                   "frames_at_day0": int((k == 0).sum()), "frames_at_day14": int((k == 168).sum())}
assert k.size == 576 and k[:24].max() == 0 and (k[528:] == 168).all() and np.all(np.diff(k) >= 0)
assert np.all(np.array(t)[k] <= tau + 1e-9)
j, al = F.bracket(np.arange(5.0), np.array([0.0, 0.5, 4.0]))
assert j.tolist() == [0, 0, 4] and al.tolist() == [0.0, 0.5, 0.0]
assert F.nice_ceiling(0.00031) == 5e-4 and F.nice_ceiling(870) == 1000 and F.nice_ceiling(1000) == 1000
assert F.nice_ceiling(1.7e-4) == 2e-4
print(res)
print("all field geometry tests passed")
