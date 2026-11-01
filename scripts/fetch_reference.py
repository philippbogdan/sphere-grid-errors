"""Download the checked NCAR source archive for local research or education."""

import argparse
from pathlib import Path
import urllib.request

from evidence import ARCHIVE_SHA256, ARCHIVE_URL, ROOT, sha256


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=ROOT / "build/stswm/stswm.tgz",
                        help="Destination for the source archive, outside the tracked file set.")
    args = parser.parse_args()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    if args.out.exists():
        if sha256(args.out) != ARCHIVE_SHA256:
            raise SystemExit("Existing archive has a different hash; it has not been replaced.")
    else:
        temporary = args.out.with_suffix(".download")
        with urllib.request.urlopen(ARCHIVE_URL, timeout=60) as source, temporary.open("wb") as dest:
            while block := source.read(1 << 20):
                dest.write(block)
        if sha256(temporary) != ARCHIVE_SHA256:
            raise SystemExit("Downloaded archive has a different hash; do not build it.")
        temporary.replace(args.out)
    print(f"Verified {args.out.name}: {ARCHIVE_SHA256}")
    print("NCAR STSWM retains its research and education terms. Read THIRD_PARTY.md.")


if __name__ == "__main__":
    main()
