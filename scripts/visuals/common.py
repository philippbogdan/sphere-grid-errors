"""Shared paths, evidence selection, style, export and audit helpers for the figure pack."""
from __future__ import annotations

import hashlib
import io
import json
import os
import re
from pathlib import Path

import numpy as np

WS = Path(__file__).resolve().parents[2]
FP = WS / "build"
SRC = Path(__file__).resolve().parent
OUT = Path(os.environ.get("SPHERE_FIGURES_OUT", str(WS / "build/figures")))
QA = Path(os.environ.get("SPHERE_FIGURES_QA", str(WS / "build/qa")))
DATA = WS / "build/fields"
EVID = WS / "data/evidence"
SCRATCH = WS / "build/render"
os.environ.setdefault("MPLCONFIGDIR", str(WS / "build/cache/matplotlib"))

MAIN_SCHEME = ("order=8,integrator=ls10o4,courant=1.252,hv_mode=step,kappa=0.0001,"
               "cone=8_3_4_lines_cc_dw10,cone_filter=0,tracer=skew")

# Palette (brief): charcoal labels, Okabe-Ito blue and vermilion, neutral greys.
INK = "#202124"
MAIN_C = "#333333"
BLUE = "#0072B2"
VERM = "#D55E00"
GREY = "#6E6E6E"
GUIDE = "#8A8F94"
SPINE = "#5F6368"

FONT = "DejaVu Sans"   # bundled with matplotlib: identical rendering on every machine


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def rel(path) -> str:
    p = Path(path).resolve()
    try:
        return str(p.relative_to(WS))
    except ValueError:
        return p.name


def load_jsonl(path):
    rows = []
    with open(path) as f:
        for i, line in enumerate(f, 1):
            line = line.strip()
            if line:
                rows.append((i, json.loads(line)))
    return rows


def select_one(fname, selector: dict, extra=None):
    """Exactly one record of data/evidence/<fname> whose fields equal `selector` (and satisfy `extra`).

    Returns (record, provenance). Raises on zero or several matches."""
    path = EVID / fname
    hits = []
    for ln, rec in load_jsonl(path):
        ok = all(rec.get(k) == v for k, v in selector.items())
        if ok and extra is not None:
            ok = extra(rec)
        if ok:
            hits.append((ln, rec))
    if len(hits) != 1:
        raise RuntimeError(f"{fname}: selector {selector} matched {len(hits)} rows "
                           f"(lines {[h[0] for h in hits]})")
    ln, rec = hits[0]
    prov = {"path": rel(path), "line": ln, "selector": selector, "sha256": sha256(path)}
    return rec, prov


def style_rc(base_pt=11.0):
    """Matplotlib rc for the static charts (sizes in points at the 300 dpi export)."""
    return {
        "font.family": FONT,
        "font.size": base_pt,
        "mathtext.fontset": "dejavusans",
        "mathtext.default": "regular",
        "axes.unicode_minus": True,
        "text.color": INK,
        "axes.labelcolor": INK,
        "axes.edgecolor": SPINE,
        "axes.linewidth": 1.05,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "xtick.color": SPINE,
        "ytick.color": SPINE,
        "xtick.labelcolor": INK,
        "ytick.labelcolor": INK,
        "xtick.major.width": 1.05,
        "ytick.major.width": 1.05,
        "xtick.minor.width": 1.05,
        "ytick.minor.width": 1.05,
        "xtick.major.size": 5.0,
        "ytick.major.size": 5.0,
        "xtick.minor.size": 3.0,
        "ytick.minor.size": 3.0,
        "legend.frameon": False,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
        "svg.fonttype": "path",
        "svg.hashsalt": "sphere-grid-errors",
        "path.simplify": False,
    }


def srgb_icc() -> bytes:
    from PIL import ImageCms
    return ImageCms.ImageCmsProfile(ImageCms.createProfile("sRGB")).tobytes()


def save_rgb_png(rgb: np.ndarray, path, dpi=300):
    """Write an opaque 8-bit sRGB PNG with dpi metadata and an embedded sRGB profile."""
    from PIL import Image
    assert rgb.dtype == np.uint8 and rgb.ndim == 3 and rgb.shape[2] == 3, rgb.shape
    img = Image.fromarray(rgb, "RGB")
    path = Path(path)
    tmp = path.with_suffix(".tmp.png")
    img.save(tmp, format="PNG", dpi=(dpi, dpi), icc_profile=srgb_icc(), optimize=True)
    os.replace(tmp, path)


def fig_to_rgb(fig, dpi) -> np.ndarray:
    """Render a figure at `dpi` without bbox changes and return opaque RGB (white composite)."""
    buf = io.BytesIO()
    fig.savefig(buf, format="rgba", dpi=dpi, facecolor="white")
    w, h = fig.get_size_inches()
    W, H = int(round(w * dpi)), int(round(h * dpi))
    a = np.frombuffer(buf.getvalue(), np.uint8).reshape(H, W, 4).astype(np.float32)
    al = a[..., 3:4] / 255.0
    rgb = a[..., :3] * al + 255.0 * (1 - al)
    return np.clip(np.round(rgb), 0, 255).astype(np.uint8)


def save_chart(fig, stem, size_px):
    """PNG (2x design size, 300 dpi metadata) and SVG with the same intrinsic pixel dimensions."""
    W, H = size_px
    w, h = fig.get_size_inches()
    dpi = W / w
    assert abs(dpi - 300) < 1e-6 and abs(h * dpi - H) < 1e-6, (w, h, dpi)
    rgb = fig_to_rgb(fig, dpi)
    assert rgb.shape == (H, W, 3), rgb.shape
    png = OUT / f"{stem}.png"
    save_rgb_png(rgb, png, dpi=300)
    buf = io.StringIO()
    fig.savefig(buf, format="svg", metadata={"Date": None, "Creator": None})
    svg = buf.getvalue()
    # Same intrinsic size as the PNG so both scale to the full content width identically.
    svg, n = re.subn(r'<svg([^>]*?) width="[^"]+" height="[^"]+"',
                     lambda m: f'<svg{m.group(1)} width="{W}" height="{H}"', svg, count=1)
    assert n == 1
    svgp = OUT / f"{stem}.svg"
    svgp.write_text(svg)
    return png, svgp


NUMERIC_TICK = re.compile(
    r"^(\$?[−\-]?(\d+([.,]\d+)*)?(\\times)?\s*(10\^\{[−\-]?\d+\})?\$?|"
    r"[−\-]?\d{1,3}(,\d{3})*(\.\d+)?)$")


def is_numeric_tick(s: str) -> bool:
    s = s.strip()
    if not s:
        return True
    core = s.replace("$", "").replace("\\mathdefault", "").replace("{", "").replace("}", "")
    core = core.replace("\\times", "x").replace("×", "x").replace("−", "-").replace(" ", "")
    return re.fullmatch(r"-?(\d+(\.\d+)?|\d{1,3}(,\d{3})+)?(x?10\^-?\d+)?", core) is not None and any(
        ch.isdigit() for ch in core)


DASHES = tuple(map(chr, (0x2013, 0x2014, 0x2012, 0x2015)))


def audit_text(fig, whitelist, typeset=None):
    """Every visible text artist must be a whitelisted string (or its declared typeset form) or a numeric tick."""
    import matplotlib.text as mtext
    fig.canvas.draw()
    typeset = typeset or {}
    allowed = set(whitelist) | set(typeset.values())
    seen, bad, boxes = [], [], []
    r = fig.canvas.get_renderer()
    W, H = fig.canvas.get_width_height()
    for t in fig.findobj(mtext.Text):
        if not t.get_visible():
            continue
        s = t.get_text()
        if not s.strip():
            continue
        bb = t.get_window_extent(r)
        kind = "label" if s in allowed else ("tick" if is_numeric_tick(s) else "VIOLATION")
        if any(d in s for d in DASHES):
            kind = "VIOLATION(dash)"
        inside = bb.x0 >= -0.5 and bb.y0 >= -0.5 and bb.x1 <= W + 0.5 and bb.y1 <= H + 0.5
        seen.append({"text": s, "kind": kind, "bbox_px": [round(bb.x0, 1), round(bb.y0, 1), round(bb.x1, 1),
                                                         round(bb.y1, 1)], "inside_canvas": bool(inside)})
        if kind.startswith("VIOLATION"):
            bad.append(s)
        if not inside:
            bad.append(f"CLIPPED: {s}")
        boxes.append((s, bb))
    overlaps = []
    for i in range(len(boxes)):
        for j in range(i + 1, len(boxes)):
            a, b = boxes[i][1], boxes[j][1]
            if min(a.x1, b.x1) - max(a.x0, b.x0) > 1 and min(a.y1, b.y1) - max(a.y0, b.y0) > 1:
                overlaps.append([boxes[i][0], boxes[j][0]])
    bad += [f"OVERLAP: {o}" for o in overlaps]
    return {"texts": seen, "violations": bad, "overlaps": overlaps, "ok": not bad}
