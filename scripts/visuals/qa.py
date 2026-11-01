"""Audit images and movies, writing metadata and visual previews to build/qa/.

python qa.py charts | movies | all
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
from PIL import Image

sys.path.insert(0, str(Path(__file__).parent))
from common import OUT, QA, sha256  # noqa: E402

FFPROBE = os.environ.get("FFPROBE", "ffprobe")
FFMPEG = os.environ.get("FFMPEG", "ffmpeg")

# Machado, Oliveira and Fernandes (2009), severity 1.0, applied in linear sRGB.
CVD = {
    "protanopia": np.array([[0.152286, 1.052583, -0.204868], [0.114503, 0.786281, 0.099216],
                            [-0.003882, -0.048116, 1.051998]]),
    "deuteranopia": np.array([[0.367322, 0.860646, -0.227968], [0.280085, 0.672501, 0.047413],
                              [-0.011820, 0.042940, 0.968881]]),
    "tritanopia": np.array([[1.255528, -0.076749, -0.178779], [-0.078411, 0.930809, 0.147602],
                            [0.004733, 0.691367, 0.303900]]),
}


def to_lin(a):
    a = a / 255.0
    return np.where(a <= 0.04045, a / 12.92, ((a + 0.055) / 1.055) ** 2.4)


def to_srgb(a):
    a = np.clip(a, 0, 1)
    return np.round(255 * np.where(a <= 0.0031308, 12.92 * a, 1.055 * a ** (1 / 2.4) - 0.055)).astype(np.uint8)


def variants(img: Image.Image, width: int):
    h = round(img.height * width / img.width)
    small = img.convert("RGB").resize((width, h), Image.LANCZOS)
    a = np.asarray(small).astype(np.float64)
    lin = to_lin(a)
    out = {"colour": small}
    y = lin @ np.array([0.2126, 0.7152, 0.0722])
    out["greyscale"] = Image.fromarray(np.repeat(to_srgb(y)[..., None], 3, axis=2))
    for k, m in CVD.items():
        out[k] = Image.fromarray(to_srgb(lin @ m.T))
    return out


def sheet(images, cols=2, pad=12):
    w = max(i.width for i in images)
    h = max(i.height for i in images)
    rows = (len(images) + cols - 1) // cols
    S = Image.new("RGB", (cols * w + (cols + 1) * pad, rows * h + (rows + 1) * pad), (200, 200, 200))
    for k, im in enumerate(images):
        r, c = divmod(k, cols)
        S.paste(im, (pad + c * (w + pad), pad + r * (h + pad)))
    return S


def image_report(png: Path):
    name = png.stem
    d = QA / "previews"
    d.mkdir(parents=True, exist_ok=True)
    img = Image.open(png)
    v850 = variants(img, 850)
    v850["colour"].save(d / f"{name}-850.png")
    img.convert("RGB").resize((390, round(img.height * 390 / img.width)), Image.LANCZOS).save(d / f"{name}-390.png")
    sheet([v850[k] for k in ("colour", "greyscale", "protanopia", "deuteranopia", "tritanopia")]).save(
        d / f"{name}-cvd-sheet.png")
    meta = QA / "metadata" / f"{name}.png.json"
    meta.parent.mkdir(parents=True, exist_ok=True)
    W, H = img.size
    checks = {"png": img.format == "PNG", "rgb": img.mode == "RGB",
              "opaque": "transparency" not in img.info,
              "dpi": min(img.info.get("dpi", (0, 0))) >= 299.99,
              "below_50MB": png.stat().st_size < 50_000_000,
              "srgb_profile": bool(img.info.get("icc_profile"))}
    meta.write_text(json.dumps(checks, indent=2) + "\n")
    if not all(checks.values()):
        raise ValueError(f"PNG metadata check failed: {png.name}: {checks}")
    return {"file": png.name, "px": [W, H], "dpi": img.info.get("dpi"), "mode": img.mode,
            "icc": bool(img.info.get("icc_profile")), "bytes": png.stat().st_size, "sha256": sha256(png),
            "metadata_screen_rc": 0, "metadata_report": str(meta.relative_to(QA.parent))}


def svg_report(svg: Path):
    s = svg.read_text()
    meta = QA / "metadata" / f"{svg.stem}.svg.json"
    root = ET.fromstring(s)
    checks = {"svg": root.tag == "{http://www.w3.org/2000/svg}svg",
              "below_50MB": svg.stat().st_size < 50_000_000,
              "outlined_text": "<text" not in s, "vector_only": "<image" not in s}
    meta.parent.mkdir(parents=True, exist_ok=True)
    meta.write_text(json.dumps(checks, indent=2) + "\n")
    if not all(checks.values()):
        raise ValueError(f"SVG metadata check failed: {svg.name}: {checks}")
    head = s[:600]
    return {"file": svg.name, "bytes": svg.stat().st_size, "sha256": sha256(svg),
            "text_elements": s.count("<text"), "image_elements": s.count("<image"),
            "root": head[head.find("<svg"):head.find(">", head.find("<svg")) + 1],
            "metadata_screen_rc": 0}


def palette_reports():
    d = QA / "palette"
    d.mkdir(parents=True, exist_ok=True)
    palettes = {
        "charts-graphical-on-white": ["333333", "0072B2", "D55E00", "6E6E6E", "8A8F94", "202124", "5F6368"],
        "chart-text-on-white": ["202124"],
    }
    res = {}
    for name, cols in palettes.items():
        role = "normal-text" if "text" in name else "graphical"
        threshold = 4.5 if role == "normal-text" else 3.0
        contrasts = {}
        for c in cols:
            rgb = np.array([int(c[i:i + 2], 16) for i in (0, 2, 4)], dtype=float)
            luminance = float(to_lin(rgb) @ np.array([0.2126, 0.7152, 0.0722]))
            contrasts[c] = 1.05 / (luminance + 0.05)
        report = {"background": "FFFFFF", "role": role, "threshold": threshold,
                  "contrast": contrasts, "ok": all(v >= threshold for v in contrasts.values())}
        (d / f"{name}.json").write_text(json.dumps(report, indent=2) + "\n")
        res[name] = {"rc": 0 if report["ok"] else 1,
                     "report": str((d / f"{name}.json").relative_to(QA.parent)), "stderr": ""}
    return res


def probe(mp4: Path):
    r = subprocess.run([FFPROBE, "-v", "error", "-show_streams", "-show_format", "-count_frames", "-of", "json",
                        str(mp4)], capture_output=True, text=True, check=True)
    j = json.loads(r.stdout)
    v = [s for s in j["streams"] if s["codec_type"] == "video"]
    a = [s for s in j["streams"] if s["codec_type"] == "audio"]
    raw = mp4.read_bytes()
    moov, mdat = raw.find(b"moov"), raw.find(b"mdat")
    s = v[0]
    if not (len(v) == 1 and not a and s["codec_name"] == "h264" and s["pix_fmt"] == "yuv420p"
            and s["r_frame_rate"] == "24/1" and int(s.get("nb_read_frames", -1)) == 576
            and abs(float(j["format"]["duration"]) - 24) < 1e-3
            and mp4.stat().st_size < 10_000_000 and 0 <= moov < mdat):
        raise ValueError(f"Movie encoding check failed: {mp4.name}")
    return {"file": mp4.name, "bytes": mp4.stat().st_size, "below_10MB": mp4.stat().st_size < 10_000_000,
            "sha256": sha256(mp4), "codec": s["codec_name"], "profile": s.get("profile"),
            "level": s.get("level"), "pix_fmt": s["pix_fmt"], "width": s["width"], "height": s["height"],
            "r_frame_rate": s["r_frame_rate"], "avg_frame_rate": s["avg_frame_rate"],
            "nb_read_frames": int(s.get("nb_read_frames", -1)), "duration_s": float(j["format"]["duration"]),
            "bit_rate": int(j["format"]["bit_rate"]), "audio_streams": len(a),
            "colour": {k: s.get(k) for k in ("color_range", "color_space", "color_transfer", "color_primaries")},
            "faststart_moov_before_mdat": 0 <= moov < mdat}


def movie_frames(mp4: Path, frames):
    """Decode the given frame indices exactly (one pass, BT.709 limited range to sRGB)."""
    d = QA / "movie-frames"
    d.mkdir(parents=True, exist_ok=True)
    for old in d.glob(f"{mp4.stem}-*.png"):
        old.unlink()
    sel = "+".join(f"eq(n\\,{f})" for f in frames)
    tmp = d / f"{mp4.stem}-tmp-%03d.png"
    subprocess.run([FFMPEG, "-v", "error", "-y", "-i", str(mp4), "-vf",
                    f"select={sel},scale=in_color_matrix=bt709:in_range=tv:flags=accurate_rnd+full_chroma_int,"
                    "format=rgb24", "-fps_mode", "passthrough", str(tmp)], check=True)
    outs = []
    for i, f in enumerate(frames, 1):
        src = d / f"{mp4.stem}-tmp-{i:03d}.png"
        dst = d / f"{mp4.stem}-f{f:04d}.png"
        src.rename(dst)
        outs.append(dst)
    return outs


def main():
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    rep_path = QA / "qa-report.json"
    rep = json.loads(rep_path.read_text()) if rep_path.exists() else {}
    if what in ("charts", "all"):
        rep["charts"] = {}
        for stem in ("relaxation", "resolution", "restart-twins"):
            rep["charts"][stem] = {"png": image_report(OUT / f"{stem}.png"), "svg": svg_report(OUT / f"{stem}.svg")}
        rep["palette"] = palette_reports()
    if what in ("movies", "all"):
        rep["movies"] = {}
        for stem in ("hero", "relaxation-fields"):
            mp4 = OUT / f"{stem}.mp4"
            if not mp4.exists():
                continue
            pr = probe(mp4)
            # frame f is displayed during [f/24, (f+1)/24): sample mid-frame
            fr = json.loads((QA / f"{stem}-frames.json").read_text())
            pick = sorted(set([0, fr["day7_frame"], fr["day10_frame"], 575] + list(range(60, 560, 90))))
            frames = movie_frames(mp4, pick)
            ims = [Image.open(p).convert("RGB") for p in frames]
            ims = [im.resize((480, round(im.height * 480 / im.width)), Image.LANCZOS) for im in ims]
            sheet(ims, cols=4).save(QA / "movie-frames" / f"{stem}-sheet.png")
            rep["movies"][stem] = {"probe": pr, "inspected_frames": [str(p.relative_to(QA.parent)) for p in frames],
                                   "inspected_frame_indices": pick}
            poster = OUT / f"{stem}-poster.png"
            if poster.exists():
                rep["movies"][stem]["poster"] = image_report(poster)
    rep_path.write_text(json.dumps(rep, indent=1, default=str))
    print(json.dumps(rep, indent=1, default=str)[:6000])


if __name__ == "__main__":
    main()
