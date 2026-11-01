"""Record figure hashes, numerical transformations and visual checks."""
from __future__ import annotations

import datetime
import json
import platform
import os
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import FONT, WS, OUT, QA, SRC, rel, sha256  # noqa: E402
from gif import recipe  # noqa: E402


def load(p):
    p = Path(p)
    return json.loads(p.read_text()) if p.exists() else None


def out_entry(name):
    p = OUT / name
    if not p.exists():
        return {"file": name, "present": False}
    e = {"file": name, "bytes": p.stat().st_size, "sha256": sha256(p), "below_50MB": p.stat().st_size < 50_000_000}
    if p.suffix == ".png":
        from PIL import Image
        with Image.open(p) as im:
            e |= {"px": list(im.size), "mode": im.mode, "dpi": [round(x) for x in im.info.get("dpi", (0, 0))],
                  "srgb_icc": bool(im.info.get("icc_profile"))}
    elif p.suffix == ".gif":
        from PIL import Image
        with Image.open(p) as im:
            e |= {"px": list(im.size), "frames": im.n_frames,
                  "palette_colours": len(im.getpalette()) // 3,
                  "loop": "forever" if im.info.get("loop") == 0 else im.info.get("loop")}
            delays = []
            for index in range(im.n_frames):
                im.seek(index)
                delays.append(im.info["duration"] // 10)
            e |= {"frame_delays_cs": delays, "duration_s": sum(delays) / 100,
                  "below_10MB": p.stat().st_size < 10_000_000}
    return e


def versions():
    import matplotlib
    import numpy
    import PIL
    import scipy
    exe = shutil.which(os.environ.get("FFMPEG", "ffmpeg"))
    ff = subprocess.run([exe, "-version"], capture_output=True, text=True).stdout.split("\n")[0] if exe else None
    return {"python": platform.python_version(), "numpy": numpy.__version__, "scipy": scipy.__version__,
            "matplotlib": matplotlib.__version__, "pillow": PIL.__version__, "ffmpeg": ff}


def main():
    charts = load(QA / "charts-provenance.json") or {}
    mc = load(QA / "movie-checks.json")
    man = {
        "schema": "sphere-grid-errors figure pack manifest v1",
        "generated_utc": datetime.datetime.now(datetime.UTC).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "specification": {"brief": rel(WS / "docs/FIELDS.md"), "sha256": sha256(WS / "docs/FIELDS.md")},
        "renderer": {"source": {rel(p): sha256(p) for p in sorted(SRC.glob("*.py"))}, "versions": versions(),
                     "font": FONT + " (text drawn as glyph outlines in SVG)",
                     "png": "opaque 8-bit sRGB, embedded sRGB ICC profile, 300 dpi stored as 11811 px/m"},
        "visuals": {},
    }
    names = {"3": ("relaxation", ["relaxation.png", "relaxation.svg"]),
             "4": ("resolution", ["resolution.png", "resolution.svg"]),
             "5": ("restart-twins", ["restart-twins.png", "restart-twins.svg"])}
    for k, (stem, files) in names.items():
        c = charts.get(k)
        if not c:
            continue
        man["visuals"][stem] = {
            "outputs": [out_entry(f) for f in files],
            "inputs": c["provenance"]["inputs"],
            "plotted_values": c["plotted"],
            "limits": c["limits"],
            "checks": c["checks"],
            "text": sorted({t["text"] for t in c["audit"]["texts"]}),
            "text_audit_ok": c["audit"]["ok"],
        }
    if mc:
        s = mc["summary"]
        shared = {
            "inputs": mc["inputs"],
            "adapter_mapping": {f"build/fields/{n}": v["path"] for n, v in mc["inputs"].items()},
            "renderer_input_checks": mc["checks"],
            "error_colour_scale": {"quantity": "E_pair = sqrt((w_i e_i^2 + w_p e_p^2)/(w_i + w_p)) per unique point",
                                   "transform": "asinh(E / 1 m) / asinh(E_max / 1 m)", "cmap": "magma",
                                   "E_largest_m": s["E_largest_m"], "E_max_m": s["E_max_m"],
                                   "E_largest_by_run_ref": s["E_largest_by_run_ref"]},
            "vorticity_colour_scale": {"cmap": "RdBu_r", "linear_symmetric_limit_s_inv": s["zeta_limit_s_inv"],
                                       "largest_abs_s_inv": s["zeta_largest_abs_s_inv"]},
            "interpolation": {"error": "piecewise linear on convex-hull facets of the unique points (radial projection),"
                                       " 2x2 supersampled per output pixel, analytic limb coverage",
                              "vorticity": "bilinear in latitude and longitude on the native Gaussian grid, periodic"
                                           " longitude; poleward of the outermost Gaussian latitude a linear blend in"
                                           " polar distance to the ring-mean pole value",
                              "triangulation": s["triangulation"], "pairs": s["pairs_main"]},
            "time": {"t_hours_main": s["t_hours_main"], "snapshot_index_by_day": s["snapshot_index_by_day"],
                     "frame_schedule": s["frame_schedule"],
                     "timing": "576 frames at 24 fps: frames 0-23 hold day 0; frames 24-527 map day 0 to day 14"
                               " linearly (inclusive); frames 528-575 hold day 14; each frame shows the most recent"
                               " saved snapshot, cursor at that snapshot's stored time",
                     "vorticity_alignment": s["vorticity_alignment_to_main_times"],
                     "run_alignment": s["run_alignment_to_main_times"]},
        }
        for stem, extra in (("hero", {"views_lon_lat_deg": [[0.0, 35.2643897], [180.0, -35.2643897]],
                                      "projection": "orthographic, geographic north up, fixed cameras"}),
                            ("relaxation-fields", {"projection": "Mollweide, central longitude 60E, seam at -120",
                                                   "contours": "actual relax_weight interpolated on the triangulation;"
                                                               " W = exp(-1) solid (cones) / dashed (controls);"
                                                               " W = 0.001 dotted"})):
            fr = load(QA / f"{stem}-frames.json")
            po = load(QA / f"{stem}-poster.json") if stem == "hero" else None
            if not fr and not po:
                continue
            files = [f"{stem}.mp4"] + (["hero-poster.png"] if stem == "hero" else [])
            v = {"outputs": [out_entry(name) for name in files]} | extra
            if fr:
                v["movie"] = {k: fr[k] for k in ("unique_frames", "encode", "day7_frame", "day10_frame",
                                                 "missing_present") if k in fr}
                v["text"] = sorted({t["text"] for t in fr["audit"]["texts"]})
                v["text_audit_ok"] = fr["audit"]["ok"]
            if po:
                v["poster"] = {k: po[k] for k in ("snapshot_index", "t_hours", "px", "missing_present")}
                v["poster_text_audit_ok"] = po["audit"]["ok"]
            qa = load(QA / "qa-report.json") or {}
            if stem in qa.get("movies", {}):
                v["probe"] = qa["movies"][stem]["probe"]
            man["visuals"][stem] = v
        man["fields"] = shared
    if (OUT / "relaxation-fields.gif").exists():
        v = man["visuals"].setdefault("relaxation-fields", {"outputs": []})
        v["outputs"].append(out_entry("relaxation-fields.gif"))
        source = OUT / "relaxation-fields.mp4"
        if not source.exists():
            source = WS / "figures/relaxation-fields.mp4"
        v["gif"] = load(QA / "relaxation-fields-gif.json") or recipe(source)
    (OUT / "manifest.json").write_text(json.dumps(man, indent=1, default=float))
    print("manifest", (OUT / "manifest.json").stat().st_size, "bytes;", list(man["visuals"]))


if __name__ == "__main__":
    main()
