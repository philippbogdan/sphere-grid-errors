# Figures

The correction GIF plays inline and links to its 24-second source movie, covering days 0 to 14. The flow movie has a day-10 poster. All visuals use DejaVu Sans; PNGs are opaque sRGB with 300 dpi metadata, and SVG text is stored as outlines.

| Visual | What it shows | Numerical source |
|---|---|---|
| [Correction GIF](relaxation-fields.gif), [source movie](relaxation-fields.mp4) | Whole-sphere height error for main, cone correction and control correction, against T85 and T213 | [Visual runs](../data/visual-runs.json), [field metrics](../data/field-metrics.csv) |
| [Flow and error movie](hero.mp4), [day-10 poster](hero-poster.png) | T213 relative vorticity above main-scheme height error against T213, with two hemispheres per row | [Field assets](../data/field-assets.json), [validation](../data/validation.json) |
| Relaxation: [PNG](relaxation.png), [SVG](relaxation.svg) | Three error histories and all 12 day-10 configurations | [relaxation.csv](../data/relaxation.csv), `tc6_week2_relax.jsonl` |
| Resolution: [PNG](resolution.png), [SVG](resolution.svg) | Main, TRSK on torus points and spectral at three equal unique-point counts | [resolution.csv](../data/resolution.csv), `accuracy_main_ref.jsonl`, `mpas_ref.jsonl`, `stswm_errors.jsonl` |
| Restart and twins: [PNG](restart-twins.png), [SVG](restart-twins.svg) | Equal elapsed-time restart comparison and three initially 1 m perturbations | [restart.csv](../data/restart.csv), [twins.csv](../data/twins.csv), corresponding original records |

The correction movie shows an early penalty followed by lower late error, with cones better than controls. The flow-and-error movie shows broad error patches, not four bright cone-centred spots.

Height error uses one fixed asinh colour scale from 0 to 500 m in both movies. Relative vorticity uses a fixed linear scale from -1e-4 to +1e-4 s^-1; [FIELDS.md](../docs/FIELDS.md) defines pairing, interpolation, actual snapshot times, contours and the visible lattice ripple.

| File | Pixels | Bytes |
|---|---:|---:|
| `hero.mp4` | 1440 by 1440 | 8,775,486 |
| `hero-poster.png` | 2880 by 2880 | 2,694,683 |
| `relaxation-fields.mp4` | 1600 by 1400 | 8,565,646 |
| `relaxation-fields.gif` | 1000 by 875 | 9,354,994 |
| `relaxation.png` | 2400 by 1800 | 194,117 |
| `resolution.png` | 2400 by 1600 | 180,428 |
| `restart-twins.png` | 2400 by 1600 | 197,557 |

Both movies contain 576 H.264 frames at 24 frames/s, use `yuv420p` with BT.709 colour metadata and have no audio. The GIF samples 43 frames at eight-hour intervals: 42 frames last 0.25 s each, and the last lasts 2 s, for a 12.5-second loop that repeats forever. It uses 128 colours, no dithering and lossy compression. [manifest.json](manifest.json) records exact input/output hashes, plotted values, fixed scales, time alignment, encoding and checks; the CSVs provide numerical alternatives to the images.

[REPRODUCE.md](../REPRODUCE.md) gives regeneration and visual-check commands. The fields are not bundled or hosted; their measured sizes and hashes are in [field-assets.json](../data/field-assets.json).
