"""Apply the TC6 output and precision changes to a locally fetched NCAR archive."""

import argparse
import os
from pathlib import Path
import shutil
import subprocess
import tarfile

from evidence import ARCHIVE_SHA256, sha256


def replace_once(path, before, after):
    text = path.read_text()
    if text.count(before) != 1:
        raise ValueError(f"Expected exactly one patch location in {path.name}")
    path.write_text(text.replace(before, after, 1))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("truncation", type=int)
    parser.add_argument("nlat", type=int)
    parser.add_argument("nlon", type=int)
    a = parser.parse_args()
    root = a.directory.resolve()
    archive = root / "stswm.tgz"
    if sha256(archive) != ARCHIVE_SHA256:
        raise SystemExit("NCAR archive hash differs from the recorded input")
    dest = root / f"build_T{a.truncation}"
    dest.mkdir(parents=True, exist_ok=False)
    # Extract fresh sources from the verified archive, never an edited source tree.
    unpacked = dest / "upstream"
    with tarfile.open(archive) as stream:
        stream.extractall(unpacked, filter="data")
    for pattern in ("src/*.f", "src/*.i", "src/netcdf.inc", "lib/adquad.f",
                    "lib/fft99f.f", "lib/netcdf.f", "lib/ncarg.f"):
        for p in unpacked.glob(pattern):
            shutil.copyfile(p, dest / p.name)
    shutil.copyfile(Path(__file__).with_name("grdout.f"), dest / "grdout.f")
    p = dest / "netcdf.inc"
    p.write_text(p.read_text().replace("1.797693134862315900e+308", "1.0e+30")
                 .replace("3.402823466e+38", "1.0e+30"))
    p = dest / "params.i"
    replace_once(p, "      PARAMETER (MM=42, NN=42, KK=42, NLAT=64, NLON=128)",
                 f"      PARAMETER (MM={a.truncation}, NN={a.truncation}, KK={a.truncation}, NLAT={a.nlat}, NLON={a.nlon})")
    replace_once(p, "      PARAMETER (MAXH=42)", f"      PARAMETER (MAXH={a.truncation})")
    replace_once(p, "      PARAMETER (RTRUNC=42)", f"      PARAMETER (RTRUNC={a.truncation})")
    p = dest / "stswm.f"
    replace_once(p, "         CALL OUTPTP (ZETASC,DIVSC,PHISC,PHIBAR,LDIAG,.FALSE.,OCTR)",
                 "         CALL GRDOUT(WS1,WS2,WS3,TAU)")
    replace_once(p, "         CALL OUTPTP (ZETASC,DIVSC,PHISC,PHIBAR,LDIAG,.TRUE.,OCTR)",
                 "         CONTINUE")
    # At day zero WS1 and WS2 still contain winds, before conversion to divergence and curl.
    key = "      CALL ZD(WS1,WS2,DIVSC,ZETASC)"
    text = p.read_text()
    index = text.index(key)
    p.write_text(text[:index] + "      IF (TAUO .LE. TAUE) CALL GRDOUT(WS1,WS2,WS3,0.0)\n" + text[index:])
    flags = ["-O2", "-std=legacy", "-fdefault-real-8", "-fdefault-double-8", "-w",
             "-fallow-argument-mismatch", "-ffixed-line-length-none"]
    subprocess.run([os.environ.get("FC", "gfortran"), *flags, "-o", "stswm",
                    *[p.name for p in sorted(dest.glob("*.f"))]], cwd=dest, check=True)
    print(f"Built NCAR STSWM derivative at T{a.truncation}: {dest / 'stswm'}")


if __name__ == "__main__":
    main()
