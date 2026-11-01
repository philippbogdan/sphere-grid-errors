"""Build the inline correction GIF from an MP4 without numerical Python packages."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]
HOURS = tuple(range(0, 337, 8))
SCALE = "scale=1000:-1:flags=lanczos"
PALETTE = ("split[a][b];[a]palettegen=max_colors=128:stats_mode=full[p];"
           "[b][p]paletteuse=dither=none:diff_mode=rectangle")
COMPRESSION = ["-O3", "--lossy=140", "--loopcount=forever"]
DELAYS = ["-d", "25", "#0-41", "-d", "200", "#42"]


def recipe(source):
    """Record the source and the exact sampling, palette and timing settings."""
    source = Path(source).resolve()
    try:
        name = str(source.relative_to(ROOT))
    except ValueError:
        name = str(source)
    return {
        "source": {"file": name, "bytes": source.stat().st_size,
                   "sha256": hashlib.sha256(source.read_bytes()).hexdigest()},
        "renderer": {"file": "scripts/visuals/gif.py",
                     "sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
        "sample_hours": list(HOURS),
        "seek_formula_seconds": "1 + 21 * h / 336 + 0.06",
        "seek_seconds": [round(1 + 21 * h / 336 + 0.06, 2) for h in HOURS],
        "scale_filter": SCALE,
        "input_frame_rate": 4,
        "palette_filter": PALETTE,
        "gifsicle_options": COMPRESSION,
        "gifsicle_frame_selection": DELAYS,
        "frame_delays_cs": [25] * 42 + [200],
        "loop": "forever",
        "note": "Lossy GIF compression changes pixel values slightly; the MP4 is the source.",
    }


def executable(variable, default):
    value = os.environ.get(variable, default)
    found = shutil.which(value)
    if not found:
        raise FileNotFoundError(f"Install {default}, or set {variable} to its executable")
    return found


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=ROOT / "figures/relaxation-fields.mp4")
    parser.add_argument("--out", type=Path, default=ROOT / "build/figures")
    parser.add_argument("--qa", type=Path, default=Path(os.environ.get(
        "SPHERE_FIGURES_QA", str(ROOT / "build/qa"))))
    args = parser.parse_args()
    source = args.input.resolve()
    provenance = recipe(source)
    ffmpeg = executable("FFMPEG", "ffmpeg")
    gifsicle = executable("GIFSICLE", "gifsicle")
    args.out.mkdir(parents=True, exist_ok=True)
    args.qa.mkdir(parents=True, exist_ok=True)
    scratch = ROOT / "build/render"
    scratch.mkdir(parents=True, exist_ok=True)
    dest = args.out / "relaxation-fields.gif"
    with tempfile.TemporaryDirectory(prefix="correction-gif-", dir=scratch) as temporary:
        frames = Path(temporary)
        for index, seconds in enumerate(provenance["seek_seconds"]):
            subprocess.run([ffmpeg, "-v", "error", "-y", "-ss", f"{seconds:.2f}",
                            "-i", str(source), "-frames:v", "1", "-vf", SCALE,
                            str(frames / f"{index:03d}.png")], check=True)
        raw = frames / "raw.gif"
        subprocess.run([ffmpeg, "-v", "error", "-y", "-framerate", "4", "-i",
                        str(frames / "%03d.png"), "-vf", PALETTE, str(raw)], check=True)
        encoded = frames / "relaxation-fields.gif"
        subprocess.run([gifsicle, *COMPRESSION, str(raw), *DELAYS, "-o", str(encoded)], check=True)
        shutil.copyfile(encoded, dest)
    provenance["versions"] = {
        "ffmpeg": subprocess.check_output([ffmpeg, "-version"], text=True).splitlines()[0],
        "gifsicle": subprocess.check_output([gifsicle, "--version"], text=True).splitlines()[0],
    }
    (args.qa / "relaxation-fields-gif.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(f"GIF: {dest} ({dest.stat().st_size:,} bytes; 43 frames, 12.5 s, loops forever)")


if __name__ == "__main__":
    main()
