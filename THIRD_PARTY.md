# Reuse and reference software

Our code, result records, tables, text and figures are dedicated under [CC0 1.0](LICENSE). Third-party software retains its own terms. No NCAR source, archive, executable or reference-field archive is included; the solver is installed separately at its pinned public revision.

## NCAR STSWM 2.0

We acknowledge NCAR STSWM as the source of the spectral model used for the reference fields and comparisons. These simulations use a derivative of NCAR STSWM 2.0, described by Hack and Jakob in *Description of a Global Shallow Water Model Based on the Spectral Transform Method*, NCAR Technical Note TN-343+STR (1992).

The [Netlib distribution](https://www.netlib.org/chammp/) supplies the source archive. Its SHA-256 is `385b001d0b116f86026e8df9f8c27f7749ed0533231c8b42541b1b8e2af2d9f1`. The archive's README gives the following copyright and distribution conditions, reproduced unchanged:

```text
SPECTRAL TRANSFORM SHALLOW WATER MODEL (Version 2.0)
Copyright (C) 1992
University Corporation for Atmospheric Research
All Rights Reserved

This spectral transform shallow water model (STSWM) is made available
``as is'' without warranty, expressed or implied, as to its suitability
for general purpose applications. The source code has been made available 
for scientific research and educational purposes only, and should not
be construed to be in the public domain. UCAR does not indemnify any 
infringement of copyright, patent, or trademark through use or 
modification of this software. NCAR/CGD Division support for the software 
is limited to copies of the source code and the distribution of 
a technical note (i.e. no consulting support will be provided).

NCAR should be acknowledged as the source of the software used in any
resulting research or publications. If users make substantial
modifications to the model source code, any references in publications
of simulation results should clearly indicate that it is a derivative
of the NCAR STSWM. Users are requested to send reprints of their work 
with this model when available to the NCAR CGD Division Office.
```

The local build changes truncation, uses 64-bit real arithmetic and substitutes our grid writer for spectral netCDF output. It writes the initial winds before the model reuses their work arrays, and uses the archive's graphics and netCDF stubs. The TC6 equations, initial condition and time integration remain those of STSWM. The build utility and grid writer are our contributions; they do not extend CC0 to upstream code.

Large reference-derived fields are regenerated locally and are not offered for download. The distribution notice does not establish general redistribution terms for those fields, so this repository makes no such grant. Keep the upstream README beside the fetched archive and consult it before using the model.

## Solver and dependencies

The [sphere-on-a-torus solver](https://github.com/philippbogdan/sphere-on-a-torus/tree/376f20bda091acee109abd183bc4684080288fc9) uses CC0 1.0. Its TRSK implementation follows published mathematics; this repository includes neither MPAS source nor an MPAS executable. The historical filename `mpas_ref.jsonl` identifies records from that TRSK implementation.

NumPy, SciPy and Numba use BSD-style licences. Matplotlib uses its PSF-based licence, and Pillow uses the HPND licence. FFmpeg and the selected H.264 encoder retain the terms of the installed build. [gifsicle](https://www.lcdf.org/gifsicle/) uses GNU GPL version 2 and compresses the inline GIF. FFmpeg, the encoder and gifsicle are external executables, not bundled software. Dependencies and fonts are outside our dedication. SVG charts contain glyph outlines; rebuilding uses DejaVu Sans, bundled with Matplotlib.
