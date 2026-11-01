"""Field mathematics for the two movies: double-cover pairing, spherical triangulation, Gaussian-grid sampling,
orthographic and Mollweide projections, fixed colour scales and time alignment (figure-pack/brief.md)."""
from __future__ import annotations

import math

import numpy as np
from scipy.spatial import ConvexHull, cKDTree

DEG = math.pi / 180.0


# ------------------------------------------------------------------ double cover
class Pairs:
    """Sites (i,j) and ((-1-i) mod M, (-1-j) mod M) share a position. One representative per pair."""

    def __init__(self, X, weight):
        X = np.asarray(X, np.float64)
        M = X.shape[1]
        assert X.shape == (3, M, M) and weight.shape == (M, M), (X.shape, weight.shape)
        I, J = np.meshgrid(np.arange(M), np.arange(M), indexing="ij")
        flat = (I * M + J).ravel()
        pflat = (((-1 - I) % M) * M + ((-1 - J) % M)).ravel()
        rep = flat < pflat
        assert rep.sum() * 2 == M * M, "every site must have a distinct partner"
        self.M = M
        self.i = flat[rep]
        self.p = pflat[rep]
        Xf = X.reshape(3, -1).T
        self.pos_diff = np.linalg.norm(Xf[self.i] - Xf[self.p], axis=1)
        U = Xf[self.i]
        self.norm_dev = np.abs(np.linalg.norm(U, axis=1) - 1.0)
        self.U = U / np.linalg.norm(U, axis=1)[:, None]
        w = np.asarray(weight, np.float64).ravel()
        self.wi, self.wp = w[self.i], w[self.p]
        self.w_rel_diff = np.abs(self.wi - self.wp) / np.maximum(np.abs(self.wi), 1e-300)

    def epair(self, e):
        """Weighted RMS of the two copies' signed errors, per unique point. e: (..., M, M)."""
        e = np.asarray(e, np.float64)
        ef = e.reshape(e.shape[:-2] + (-1,))
        ei, ep = ef[..., self.i], ef[..., self.p]
        return np.sqrt((self.wi * ei ** 2 + self.wp * ep ** 2) / (self.wi + self.wp))

    def stats(self):
        return {"pairs": int(self.i.size), "max_position_discrepancy": float(self.pos_diff.max()),
                "max_norm_deviation": float(self.norm_dev.max()),
                "max_weight_rel_difference": float(self.w_rel_diff.max())}


# ------------------------------------------------------------------ triangulation
class SphereTri:
    """Convex-hull triangulation of unique unit vectors; piecewise-linear interpolation on its facets."""

    def __init__(self, U):
        U = np.asarray(U, np.float64)
        hull = ConvexHull(U)
        tri = hull.simplices.copy()
        A, B, C = U[tri[:, 0]], U[tri[:, 1]], U[tri[:, 2]]
        nrm = np.cross(B - A, C - A)
        flip = np.einsum("ij,ij->i", nrm, A + B + C) < 0
        tri[flip] = tri[flip][:, [0, 2, 1]]
        self.U, self.tri = U, tri
        A, B, C = U[tri[:, 0]], U[tri[:, 1]], U[tri[:, 2]]
        Mt = np.stack([A, B, C], axis=-1)            # columns a, b, c
        self.Minv = np.linalg.inv(Mt)
        cen = A + B + C
        self.cen = cen / np.linalg.norm(cen, axis=1)[:, None]
        self.tree = cKDTree(self.cen)
        N, F = len(U), len(tri)
        edges = np.concatenate([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [2, 0]]])
        ang = np.degrees(np.arccos(np.clip(np.einsum("ij,ij->i", U[edges[:, 0]], U[edges[:, 1]]), -1, 1)))
        d_nn, _ = cKDTree(U).query(U, k=2)
        nn = np.degrees(2 * np.arcsin(d_nn[:, 1] / 2))
        self.check = {"points": N, "faces": F, "euler_faces_expected": 2 * N - 4,
                      "all_points_are_vertices": bool(np.unique(tri).size == N),
                      "edge_deg_max": float(ang.max()), "edge_deg_p99": float(np.percentile(ang, 99)),
                      "edge_deg_median": float(np.median(ang)),
                      "nn_deg_max": float(nn.max()), "nn_deg_median": float(np.median(nn)),
                      "max_edge_over_max_nn": float(ang.max() / nn.max())}
        assert F == 2 * N - 4 and self.check["all_points_are_vertices"], self.check

    def locate(self, Q, k=8, chunk=250_000):
        """Facet index and barycentric weights (radial projection onto the facet) for unit vectors Q (n,3)."""
        Q = np.asarray(Q, np.float64)
        n = len(Q)
        fid = np.full(n, -1, np.int64)
        bary = np.zeros((n, 3), np.float64)
        for s in range(0, n, chunk):
            q = Q[s:s + chunk]
            f, b = self._locate_chunk(q, k)
            fid[s:s + chunk], bary[s:s + chunk] = f, b
        assert (fid >= 0).all()
        return fid, bary

    def _locate_chunk(self, q, k):
        m = len(q)
        fid = np.full(m, -1, np.int64)
        bary = np.zeros((m, 3))
        todo = np.arange(m)
        kk = k
        while todo.size:
            _, cand = self.tree.query(q[todo], k=kk)
            cand = np.atleast_2d(cand)
            found = np.zeros(todo.size, bool)
            for j in range(cand.shape[1]):
                t = cand[:, j]
                x = np.einsum("nij,nj->ni", self.Minv[t], q[todo])
                ok = (x >= -1e-12).all(1) & ~found
                if ok.any():
                    sel = todo[ok]
                    fid[sel] = t[ok]
                    xs = np.clip(x[ok], 0, None)
                    bary[sel] = xs / xs.sum(1)[:, None]
                    found |= ok
            todo = todo[~found]
            kk = min(kk * 4, len(self.tri))
            if kk == len(self.tri) and todo.size:
                raise RuntimeError(f"{todo.size} points not located")
        return fid, bary


class TriSampler:
    def __init__(self, tri: SphereTri, Q):
        f, b = tri.locate(Q)
        self.idx = tri.tri[f].astype(np.int32)
        self.w = b.astype(np.float32)

    def __call__(self, vals):
        return (np.asarray(vals, np.float32)[self.idx] * self.w).sum(1)


# ------------------------------------------------------------------ Gaussian grid sampling (vorticity)
class GaussSampler:
    """Bilinear in (latitude, longitude) on the native Gaussian grid, periodic in longitude. Poleward of the
    outermost Gaussian latitude the value blends linearly in polar distance from the ring value to a single pole
    value (the ring mean), so the field is single-valued at each pole and continuous across it."""

    def __init__(self, lat_rad, lon_rad, Q):
        lat_rad = np.asarray(lat_rad, np.float64)
        nlat, nlon = lat_rad.size, np.asarray(lon_rad).size
        assert np.all(np.diff(lat_rad) > 0), "latitudes must ascend"
        assert np.allclose(lon_rad, 2 * np.pi * np.arange(nlon) / nlon, atol=1e-12, rtol=0)
        self.nlat, self.nlon = nlat, nlon
        Q = np.asarray(Q, np.float64)
        lat = np.arcsin(np.clip(Q[:, 2], -1, 1))
        lon = np.arctan2(Q[:, 1], Q[:, 0]) % (2 * np.pi)
        fl = lon / (2 * np.pi / nlon)
        i0f = np.floor(fl)
        fx = fl - i0f
        i0 = i0f.astype(np.int64) % nlon
        i1 = (i0 + 1) % nlon
        n = len(Q)
        idx = np.zeros((n, 4), np.int64)
        w = np.zeros((n, 4), np.float64)
        NP, SP = nlat * nlon, nlat * nlon + 1
        j = np.searchsorted(lat_rad, lat) - 1
        mid = (j >= 0) & (j < nlat - 1)
        jm = j[mid]
        fy = (lat[mid] - lat_rad[jm]) / (lat_rad[jm + 1] - lat_rad[jm])
        a, b = i0[mid], i1[mid]
        idx[mid] = np.stack([jm * nlon + a, jm * nlon + b, (jm + 1) * nlon + a, (jm + 1) * nlon + b], 1)
        fxm = fx[mid]
        w[mid] = np.stack([(1 - fxm) * (1 - fy), fxm * (1 - fy), (1 - fxm) * fy, fxm * fy], 1)
        north = lat >= lat_rad[-1]
        north &= ~mid
        s = (np.pi / 2 - lat[north]) / (np.pi / 2 - lat_rad[-1])
        r = nlat - 1
        idx[north] = np.stack([r * nlon + i0[north], r * nlon + i1[north], np.full(s.size, NP), np.full(s.size, NP)], 1)
        w[north] = np.stack([s * (1 - fx[north]), s * fx[north], 1 - s, np.zeros_like(s)], 1)
        south = (lat <= lat_rad[0]) & ~mid & ~north
        s = (lat[south] + np.pi / 2) / (lat_rad[0] + np.pi / 2)
        idx[south] = np.stack([i0[south], i1[south], np.full(s.size, SP), np.full(s.size, SP)], 1)
        w[south] = np.stack([s * (1 - fx[south]), s * fx[south], 1 - s, np.zeros_like(s)], 1)
        assert (mid | north | south).all()
        self.idx, self.w = idx.astype(np.int32), w.astype(np.float32)
        self.counts = {"interior": int(mid.sum()), "north_cap": int(north.sum()), "south_cap": int(south.sum())}

    def __call__(self, z):
        z = np.asarray(z, np.float64)
        assert z.shape == (self.nlat, self.nlon)
        ext = np.concatenate([z.ravel(), [z[-1].mean(), z[0].mean()]]).astype(np.float32)
        return (ext[self.idx] * self.w).sum(1)


# ------------------------------------------------------------------ projections
def frame_vectors(lon0_deg, lat0_deg):
    l, p = lon0_deg * DEG, lat0_deg * DEG
    C = np.array([math.cos(p) * math.cos(l), math.cos(p) * math.sin(l), math.sin(p)])
    E = np.array([-math.sin(l), math.cos(l), 0.0])
    N = np.array([-math.sin(p) * math.cos(l), -math.sin(p) * math.sin(l), math.cos(p)])
    return C, E, N


def ortho_inverse(x, y, lon0, lat0):
    """Screen (x right, y up) in the unit disc -> unit vectors. Points outside the disc are put on the limb."""
    C, E, N = frame_vectors(lon0, lat0)
    r = np.hypot(x, y)
    sc = np.where(r > 1, 1 / np.maximum(r, 1e-300), 1.0)
    x, y = x * sc, y * sc
    z = np.sqrt(np.clip(1 - x * x - y * y, 0, None))
    return x[:, None] * E + y[:, None] * N + z[:, None] * C


def ortho_forward(P, lon0, lat0):
    C, E, N = frame_vectors(lon0, lat0)
    P = np.atleast_2d(P)
    return P @ E, P @ N, P @ C           # x, y, depth (visible where > 0)


SQ2 = math.sqrt(2.0)


def moll_inverse(x, y, lon0):
    """Mollweide (x in [-2 sqrt2, 2 sqrt2], y in [-sqrt2, sqrt2]) -> lat, lon (radians). Outside points are put on
    the boundary ellipse."""
    g = np.sqrt((x / (2 * SQ2)) ** 2 + (y / SQ2) ** 2)
    sc = np.where(g > 1, 1 / np.maximum(g, 1e-300), 1.0)
    x, y = x * sc, y * sc
    th = np.arcsin(np.clip(y / SQ2, -1, 1))
    lat = np.arcsin(np.clip((2 * th + np.sin(2 * th)) / np.pi, -1, 1))
    c = np.cos(th)
    lam = np.where(c > 1e-12, np.pi * x / (2 * SQ2 * np.maximum(c, 1e-12)), 0.0)
    lam = np.clip(lam, -np.pi, np.pi)
    return lat, lon0 * DEG + lam


def moll_theta(lat):
    lat = np.asarray(lat, np.float64)
    th = lat.copy()
    target = np.pi * np.sin(lat)
    for _ in range(60):
        f = 2 * th + np.sin(2 * th) - target
        d = 2 + 2 * np.cos(2 * th)
        step = np.where(d > 1e-12, f / np.maximum(d, 1e-12), 0.0)
        th = th - step
    th = np.where(np.abs(np.abs(lat) - np.pi / 2) < 1e-9, np.sign(lat) * np.pi / 2, th)
    return th


def moll_forward(lat, lon, lon0):
    """lat, lon in radians. Returns x, y and the relative longitude in (-pi, pi]."""
    lam = (np.asarray(lon) - lon0 * DEG + np.pi) % (2 * np.pi) - np.pi
    th = moll_theta(lat)
    return 2 * SQ2 / np.pi * lam * np.cos(th), SQ2 * np.sin(th), lam


def latlon_to_vec(lat, lon):
    lat, lon = np.asarray(lat, np.float64), np.asarray(lon, np.float64)
    return np.stack([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)], -1)


def vec_to_latlon(P):
    P = np.atleast_2d(P)
    return np.arcsin(np.clip(P[:, 2], -1, 1)), np.arctan2(P[:, 1], P[:, 0])


# ------------------------------------------------------------------ pixel samplers with analytic limb coverage
class DiscGrid:
    """Active pixels of a D x D globe image, ss x ss subsamples each, and analytic disc coverage (alpha)."""

    def __init__(self, D, ss, lon0, lat0):
        R = D / 2.0
        c = np.arange(D) + 0.5 - R
        X, Y = np.meshgrid(c, -c)                      # row 0 at the top
        r = np.hypot(X, Y)
        act = r <= R + 1.0
        self.rows, self.cols = np.nonzero(act)
        self.alpha = np.clip(R - r[act] + 0.5, 0, 1).astype(np.float32)
        o = (np.arange(ss) + 0.5) / ss - 0.5
        ox, oy = np.meshgrid(o, -o)
        sx = (X[act][:, None] + ox.ravel()[None]) / R
        sy = (Y[act][:, None] + oy.ravel()[None]) / R
        self.D, self.ss, self.n = D, ss, act.sum()
        self.Q = ortho_inverse(sx.ravel(), sy.ravel(), lon0, lat0)


class MollGrid:
    """Active pixels of a Wm x Hm (2:1) Mollweide map, subsamples, analytic ellipse coverage."""

    def __init__(self, Wm, Hm, ss, lon0):
        assert Wm == 2 * Hm
        a_px, b_px = Wm / 2.0, Hm / 2.0
        cx = np.arange(Wm) + 0.5 - a_px
        cy = -(np.arange(Hm) + 0.5 - b_px)
        X, Y = np.meshgrid(cx, cy)
        g = np.sqrt((X / a_px) ** 2 + (Y / b_px) ** 2)
        grad = np.sqrt((X / a_px ** 2) ** 2 + (Y / b_px ** 2) ** 2) / np.maximum(g, 1e-12)
        dist = (1 - g) / np.maximum(grad, 1e-12)
        act = dist >= -1.0
        self.rows, self.cols = np.nonzero(act)
        self.alpha = np.clip(dist[act] + 0.5, 0, 1).astype(np.float32)
        o = (np.arange(ss) + 0.5) / ss - 0.5
        ox, oy = np.meshgrid(o, -o)
        k = 2 * SQ2 / a_px                               # projection units per pixel
        sx = (X[act][:, None] + ox.ravel()[None]) * k
        sy = (Y[act][:, None] + oy.ravel()[None]) * k
        lat, lon = moll_inverse(sx.ravel(), sy.ravel(), lon0)
        self.Wm, self.Hm, self.ss, self.n = Wm, Hm, ss, act.sum()
        self.Q = latlon_to_vec(lat, lon)


# ------------------------------------------------------------------ colour scales
def nice_ceiling(v):
    """Next 1, 2 or 5 times a power of ten at or above v (v > 0)."""
    e = math.floor(math.log10(v))
    for m in (1, 2, 5, 10):
        c = m * 10 ** e
        if c >= v * (1 - 1e-12):
            return float(c)
    raise AssertionError


class Scale:
    def __init__(self, cmap_name, kind, vmax, n=2048):
        import matplotlib as mpl
        self.lut = mpl.colormaps[cmap_name](np.linspace(0, 1, n))[:, :3].astype(np.float32)
        self.kind, self.vmax, self.n = kind, float(vmax), n
        self.cmap_name = cmap_name

    def u(self, v):
        v = np.asarray(v, np.float64)
        if self.kind == "asinh":                     # E >= 0, 1 m linear width
            return np.arcsinh(v / 1.0) / np.arcsinh(self.vmax / 1.0)
        return (v + self.vmax) / (2 * self.vmax)      # symmetric linear

    def rgb(self, v):
        u = self.u(v)
        bad = ~np.isfinite(u)
        k = np.clip(np.round(np.nan_to_num(u, nan=0.0) * (self.n - 1)), 0, self.n - 1).astype(np.int32)
        out = self.lut[k]
        return out, bad


HATCH_A = np.array([0.62, 0.62, 0.62], np.float32)
HATCH_B = np.array([0.78, 0.78, 0.78], np.float32)


def compose(grid, vals, scale, shape, missing_ok=False):
    """Colour subsample values, average to pixels (anti-aliasing), composite over white. Returns uint8 RGB."""
    ss2 = grid.ss ** 2
    col, bad = scale.rgb(vals)
    if bad.any():
        if not missing_ok:
            raise RuntimeError(f"{int(bad.sum())} missing samples: stop and diagnose")
        r = np.repeat(grid.rows, ss2)
        c = np.repeat(grid.cols, ss2)
        stripe = ((r + c) // 6) % 2 == 0
        col[bad] = np.where(stripe[bad, None], HATCH_A, HATCH_B)
    pix = col.reshape(grid.n, ss2, 3).mean(1)
    a = grid.alpha[:, None]
    img = np.ones(shape + (3,), np.float32)
    img[grid.rows, grid.cols] = a * pix + (1 - a)
    return (np.clip(img, 0, 1) * 255 + 0.5).astype(np.uint8), bool(bad.any())


# ------------------------------------------------------------------ time
N_FRAMES, FPS, HOLD0, SEG, HOLD1 = 576, 24, 24, 504, 48
T_END_H = 336.0


def display_schedule(t_hours):
    """Frame -> displayed snapshot index. 1 s at day 0, 21 s mapping day 0 to day 14 linearly (frames 24 to 527
    inclusive), 2 s at day 14. Each frame shows the most recent saved snapshot at or before its mapped time."""
    t = np.asarray(t_hours, np.float64)
    assert np.all(np.diff(t) > 0) and abs(t[0]) < 1e-9 and abs(t[-1] - T_END_H) < 1e-6, (t[0], t[-1])
    tau = np.empty(N_FRAMES)
    f = np.arange(N_FRAMES)
    tau[f < HOLD0] = 0.0
    seg = (f >= HOLD0) & (f < HOLD0 + SEG)
    tau[seg] = T_END_H * (f[seg] - HOLD0) / (SEG - 1)
    tau[f >= HOLD0 + SEG] = T_END_H
    k = np.searchsorted(t, tau + 1e-9, side="right") - 1
    return tau, k


def bracket(t_src, t_dst):
    """Linear alignment of a saved series to target times: indices j and weights a with
    value(t_dst) = (1-a) v[j] + a v[j+1]. No extrapolation."""
    t_src = np.asarray(t_src, np.float64)
    t_dst = np.asarray(t_dst, np.float64)
    if t_dst.min() < t_src[0] - 1e-9 or t_dst.max() > t_src[-1] + 1e-9:
        raise RuntimeError("alignment would extrapolate")
    j = np.clip(np.searchsorted(t_src, t_dst, side="right") - 1, 0, len(t_src) - 2)
    a = (t_dst - t_src[j]) / (t_src[j + 1] - t_src[j])
    a = np.where(np.abs(a) < 1e-12, 0.0, a)
    # an exact hit on the upper node
    hit = np.abs(a - 1) < 1e-12
    j = np.where(hit, j + 1, j)
    a = np.where(hit, 0.0, a)
    j = np.minimum(j, len(t_src) - 1)
    return j, np.clip(a, 0, 1)
