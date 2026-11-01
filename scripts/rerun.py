"""Build local references or run the archived relaxation configurations."""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time

from evidence import CONFIGURATIONS, ROOT, sha256

SCRIPTS = Path(__file__).resolve().parent


def namelist(dt, hdc, output_hours=2.0):
    """Use the same TC6 namelist and decimal diffusion coefficient as the visual runs."""
    return f""" $PHYVAR
 ALPHA  =  0.0,
 HDC    =  {hdc:.4E},
 $END
 $EXPDEF
 DT     =   {dt},
 EGYFRQ =   24.0,
 ERRFRQ =   999.0,
 SPCFRQ =   999.0,
 TAUE   =   336.0,
 TAUO   =   {output_hours},
 GPHFRQ =   999.0,
 SITS   =   .TRUE.,
 FORCED =   .FALSE.,
 MOMENT =   .FALSE.,
 ICOND  =   6,
 $END
 $PLTDEF
 LGPHS  = .FALSE.,
 $END
 $FNAMES
 $END
"""


def references(args):
    root = args.build.resolve() / "stswm"
    subprocess.run([sys.executable, str(SCRIPTS / "fetch_reference.py"),
                    "--out", str(root / "stswm.tgz")], check=True)
    configurations = {"T85": (85, 128, 256, 45.0, 8e12*(213*214/(85*86))**2, "T85dt45"),
                      "T213": (213, 320, 640, 22.5, 8e12, "T213dt22")}
    for name in (["T85", "T213"] if args.which == "both" else [args.which]):
        trunc, nlat, nlon, dt, hdc, tag = configurations[name]
        executable = root / f"build_T{trunc}/stswm"
        if not executable.exists():
            subprocess.run([sys.executable, str(SCRIPTS / "build_reference.py"), str(root),
                            str(trunc), str(nlat), str(nlon)], check=True)
        dest = root / "runs" / tag
        dest.mkdir(parents=True, exist_ok=False)
        config = dest / "exp"
        output_hours = 2.0
        config.write_text(namelist(dt, hdc, output_hours))
        start = time.perf_counter()
        with config.open("rb") as inp, (dest / "stdout").open("wb") as out, (dest / "stderr").open("wb") as err:
            subprocess.run([str(executable)], cwd=dest, stdin=inp, stdout=out, stderr=err, check=True)
        wall = time.perf_counter() - start
        field = dest / "grid_out.bin"
        if not field.exists():
            raise RuntimeError("The reference run produced no grid_out.bin; inspect its local logs")
        metadata = {"model": "NCAR STSWM 2.0 derivative", "case": "tc6", "truncation": trunc,
                    "nlat": nlat, "nlon": nlon, "dt_s": dt, "hdc": float(f"{hdc:.4E}"),
                    "duration_hours": 336, "output_interval_hours": output_hours, "wall_s": wall,
                    "namelist_sha256": sha256(config), "executable_sha256": sha256(executable),
                    "field_sha256": sha256(field), "field_bytes": field.stat().st_size}
        (dest / "reference.json").write_text(json.dumps(metadata, indent=2) + "\n")
        print(f"{tag}: {wall:.1f} s; {field}", flush=True)


def relax(args):
    if args.configuration == "visual":
        chosen = [(("none", 0, 0), "none"), (("bp", 20, 0), "cones"), (("bp", 20, 60), "controls")]
    elif args.configuration == "all":
        chosen = [(cfg, f"config-{i:02d}-{cfg[0]}-{cfg[1]}-{cfg[2]}") for i, cfg in enumerate(CONFIGURATIONS)]
    else:
        index = int(args.configuration)
        cfg = CONFIGURATIONS[index]
        chosen = [(cfg, f"config-{index:02d}-{cfg[0]}-{cfg[1]}-{cfg[2]}")]
    t85 = args.build / "stswm/runs/T85dt45/grid_out.bin"
    t213 = args.build / "stswm/runs/T213dt22/grid_out.bin"
    if not t85.is_file():
        raise FileNotFoundError("Generate the T85 reference before running relaxation")
    extra = ["--t213", str(t213)] if t213.is_file() else []
    for (region, par, rotation), tag in chosen:
        subprocess.run([sys.executable, str(SCRIPTS / "run_fields.py"),
                        "--region", region, "--par", str(par), "--rot", str(rotation),
                        "--t85", str(t85), *extra,
                        "--days", "1,3,5,7,10,14" if args.end_day == 14 else "1,3,5,7,10",
                        "--out", str(args.build / "runs" / tag)], check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build", type=Path, default=ROOT / "build", help="Untracked output root.")
    sub = parser.add_subparsers(dest="mode", required=True)
    ref = sub.add_parser("references", help="Fetch, build and run TC6 at T85 and T213 through day 14.")
    ref.add_argument("--which", choices=["both", "T85", "T213"], default="both")
    runs = sub.add_parser("relax", help="Run the visual trio, all 12 archived configurations or one index.")
    runs.add_argument("--configuration", default="visual", choices=["visual", "all"] + [str(i) for i in range(12)],
                      help="Configuration indices are listed in docs/METHOD.md.")
    runs.add_argument("--end-day", type=int, choices=[10, 14], default=14)
    args = parser.parse_args()
    args.build = args.build.resolve()
    for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
        os.environ[key] = "1"
    {"references": references, "relax": relax}[args.mode](args)


if __name__ == "__main__":
    main()
