"""Read the included records with explicit, unique selectors using the standard library."""

import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
EVIDENCE = DATA / "evidence"
REVISION = "376f20bda091acee109abd183bc4684080288fc9"
SOLVER_URL = "https://github.com/philippbogdan/sphere-on-a-torus"
ARCHIVE_URL = "https://www.netlib.org/chammp/stswm.tgz"
ARCHIVE_SHA256 = "385b001d0b116f86026e8df9f8c27f7749ed0533231c8b42541b1b8e2af2d9f1"
MAIN = ("order=8,integrator=ls10o4,courant=1.252,hv_mode=step,kappa=0.0001,"
        "cone=8_3_4_lines_cc_dw10,cone_filter=0,tracer=skew")
CONFIGURATIONS = [
    ("none", 0, 0), ("bp", 10, 0), ("bp", 20, 0),
    ("disc", 5, 0), ("disc", 10, 0), ("disc", 20, 0),
    ("annulus", 10, 0), ("polar", 60, 0), ("bp", 20, 60),
    ("disc", 20, 60), ("annulus", 10, 60), ("disc", 10, 60),
]
METRIC_COLUMNS = ("run,reference,t_hours,region,area_share,weighted_error_sq_m2,"
                  "reference_height_sq_m2,l2,rms_m,error_energy_share,enrichment,"
                  "pair_rms_difference_m,pair_max_difference_m").split(",")


def sha256(path):
    """Hash bytes incrementally, including large locally generated fields."""
    h = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def records(name):
    """Return source line numbers and decoded records without altering the source."""
    return [(i, json.loads(line)) for i, line in
            enumerate((EVIDENCE / name).read_text().splitlines(), 1) if line.strip()]


def one(name, selector):
    """Require exactly one record whose named fields equal the selector."""
    hits = [(i, r) for i, r in records(name)
            if all(r.get(k) == v for k, v in selector.items())]
    if len(hits) != 1:
        raise ValueError(f"{name}: expected one record for {selector}; found {len(hits)}")
    line, row = hits[0]
    return row, {"file": f"data/evidence/{name}", "line": line,
                 "selector": selector, "sha256": sha256(EVIDENCE / name)}


def relaxation(region, par, rotation):
    return one("tc6_week2_relax.jsonl", {
        "mode": "relax", "M": 128, "region": region, "par": par,
        "control_rotation_deg": rotation, "scheme": MAIN})


def resolution(M, truncation, dt):
    """Select main, TRSK and spectral endpoints at equal unique-point counts."""
    rows = [
        one("accuracy_main_ref.jsonl", {"case": "tc6", "days": 14.0,
            "grid": ["torus", M], "spec": f"torus:{M}:tc6:{MAIN}"}),
        one("mpas_ref.jsonl", {"solver": "trsk", "case": "tc6", "days": 14.0,
            "grid": ["torus", M], "mesh": "torus"}),
        one("stswm_errors.jsonl", {"case": "tc6", "against": "T213dt22",
            "truncation": truncation, "dt": dt}),
    ]
    main, trsk, spectral = [r for r, _ in rows]
    if not (main["nunique"] == trsk["nunique"] == spectral["npts"] == M * M // 2):
        raise ValueError("The selected resolutions have different unique-point counts")
    if main["scheme"]["courant"] != 1.252:
        raise ValueError("The main endpoint uses a different Courant number")
    return rows
