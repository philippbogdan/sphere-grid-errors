"""Render charts, movies and the inline GIF from saved records and local fields."""

import argparse
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
VISUALS = Path(__file__).resolve().parent / "visuals"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=["charts", "check", "posters", "movies", "gif", "qa", "manifest", "all"])
    parser.add_argument("--data", type=Path, default=ROOT / "build/fields", help="Local field pack.")
    parser.add_argument("--out", type=Path, default=ROOT / "build/figures", help="Destination for rendered files.")
    parser.add_argument("--input", type=Path, help="Source MP4 for gif mode; defaults to the bundled movie.")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, SPHERE_FIGURES_OUT=str(args.out.resolve()),
               MPLCONFIGDIR=str(ROOT / "build/cache/matplotlib"))

    def run(script, *argv):
        python = [sys.executable] + (["-S"] if sys.flags.no_site else [])
        subprocess.run([*python, str(VISUALS / script), *map(str, argv)], env=env, check=True)

    if args.mode in ("charts", "all"):
        run("charts.py")
    if args.mode in ("check", "posters", "movies", "all"):
        run("movies.py", "check", "--data", args.data)
    if args.mode in ("posters", "movies", "all"):
        run("movies.py", "poster", "hero", "--data", args.data, "--out", args.out)
    if args.mode in ("movies", "all"):
        for which in ("hero", "fields"):
            run("movies.py", "movie", which, "--data", args.data, "--out", args.out)
    if args.mode in ("gif", "all"):
        source = args.out / "relaxation-fields.mp4" if args.mode == "all" else (
            args.input or ROOT / "figures/relaxation-fields.mp4")
        run("gif.py", "--input", source, "--out", args.out)
    if args.mode in ("qa", "all"):
        run("qa.py", "all")
    if args.mode in ("charts", "posters", "movies", "qa", "manifest", "all"):
        run("make_manifest.py")


if __name__ == "__main__":
    main()
