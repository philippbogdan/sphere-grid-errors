"""Relative vorticity from Gaussian-grid winds, with a pole-consistent spline derivative."""

import numpy as np
from scipy.interpolate import CubicSpline
from polarswe.testcases import A

def zeta(u, v, phi, ext=8):
    """Return curl in s^-1 for winds (nlat, nlon) and ascending latitude in radians."""

    nlat, nlon = u.shape
    k = np.fft.rfftfreq(nlon, 1.0 / nlon)
    dv = np.fft.irfft(1j * k * np.fft.rfft(v, axis=1), n=nlon, axis=1)
    g = u * np.cos(phi)[:, None]
    # The half-turn and vector sign make u*cos(phi) even across each pole.
    gs = np.roll(g, nlon // 2, axis=1)
    lat_e = np.concatenate([-np.pi - phi[:ext][::-1], phi, np.pi - phi[-ext:][::-1]])
    g_e = np.concatenate([gs[:ext][::-1], g, gs[-ext:][::-1]], axis=0)
    dg = CubicSpline(lat_e, g_e, axis=0)(phi, 1)
    return (dv - dg) / (A * np.cos(phi)[:, None])

def zeta_fd(u, v, phi):
    """Cross-check the latitude derivative using second-order finite differences."""

    nlat, nlon = u.shape
    k = np.fft.rfftfreq(nlon, 1.0 / nlon)
    dv = np.fft.irfft(1j * k * np.fft.rfft(v, axis=1), n=nlon, axis=1)
    dg = np.gradient(u * np.cos(phi)[:, None], phi, axis=0, edge_order=2)
    return (dv - dg) / (A * np.cos(phi)[:, None])

if __name__ == "__main__":

    for nlat in (128, 320):
        nlon = 2 * nlat
        x, _ = np.polynomial.legendre.leggauss(nlat); phi = np.arcsin(np.sort(x))
        lam = 2 * np.pi * np.arange(nlon) / nlon
        P, L = np.meshgrid(phi, lam, indexing="ij")
        w = K = 7.848e-6; R = 4
        u = A * w * np.cos(P) + A * K * np.cos(P) ** (R - 1) * (R * np.sin(P) ** 2 - np.cos(P) ** 2) * np.cos(R * L)
        v = -A * K * R * np.cos(P) ** (R - 1) * np.sin(P) * np.sin(R * L)
        z_true = 2 * w * np.sin(P) - K * np.sin(P) * np.cos(P) ** R * (R ** 2 + 3 * R + 2) * np.cos(R * L)
        for name, f in (("spline", zeta), ("fd2", zeta_fd)):
            d = f(u, v, phi) - z_true
            pol = np.abs(np.degrees(phi)) > 80
            print(nlat, name, "rms", np.sqrt((d ** 2).mean()), "max", np.abs(d).max(), "polar max", np.abs(d[pol]).max(), "zmax", np.abs(z_true).max())
