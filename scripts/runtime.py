"""Keep solver caches local and require the recorded public solver revision."""

import importlib.metadata
import json
import os
from pathlib import Path

from evidence import MAIN, REVISION, ROOT


def configure():
    """Use one numerical thread and verify the pip VCS installation before a run."""
    for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
                "NUMBA_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        os.environ[key] = "1"
    os.environ["POLARSWE_DEVICE"] = "cpu"
    os.environ["POLARSWE_LATTICE"] = "hex"
    cache = ROOT / "build/cache"
    cache.mkdir(parents=True, exist_ok=True)
    os.environ["NUMBA_CACHE_DIR"] = str(cache / "numba")
    dist = importlib.metadata.distribution("sphere-on-a-torus")
    provenance = json.loads(dist.read_text("direct_url.json") or "{}")
    if provenance.get("vcs_info", {}).get("commit_id") != REVISION:
        raise RuntimeError("Install requirements.txt: the solver must come from the pinned git revision")
    from polarswe import configs, grid
    if configs.MAIN != MAIN:
        raise RuntimeError("The installed main scheme differs from the recorded configuration")
    grid.CACHE = str(cache / "solver")


def reference_read(path):
    """Expose the original binary values as memory-mapped record views.

    This keeps the solver reader's order and precision while avoiding a second
    resident copy of every reference wind and height field.
    """
    import numpy as np
    path = Path(path)
    with path.open("rb") as stream:
        header = np.fromfile(stream, dtype=np.int32, count=2)
        if header.size != 2:
            raise ValueError("Reference file has no complete grid header")
        nlon, nlat = map(int, header)
        if nlon <= 0 or nlat <= 0:
            raise ValueError("Reference dimensions must be positive")
        glat = np.fromfile(stream, dtype=np.float64, count=nlat)
    offset = 8 + 8 * nlat
    width = 1 + 3 * nlat * nlon
    remaining = path.stat().st_size - offset
    if glat.size != nlat or remaining <= 0 or remaining % (8 * width):
        raise ValueError("Reference file ends inside a record")
    raw = np.memmap(path, mode="r", dtype=np.float64, offset=offset,
                    shape=(remaining // (8 * width), width))
    recs = []
    for row in raw:
        u, v, h = row[1:].reshape(3, nlat, nlon)
        recs.append((float(row[0]), u, v, h))
    times = np.array([r[0] for r in recs])
    if not np.allclose(times, 2 * np.arange(len(times)), rtol=0, atol=1e-6):
        raise ValueError("Expected reference records every two hours, starting at zero")
    return nlon, nlat, glat, recs
