# Field data and visual transformations

The field workflow turns locally regenerated TC6 solutions into two movies, while retaining the numerical scores and the transformations used to draw them. [REPRODUCE.md](../REPRODUCE.md) gives the commands. Large fields, references, caches and temporary frames stay under the untracked `build/` directory.

## Inputs and schemas

The reference writer stores two native-endian int32 values, `nlon,nlat`, followed by `nlat` float64 Gaussian latitudes in radians. Each record then contains a float64 time in hours followed by float64 eastward wind u, northward wind v and height h, in that order. Each field has shape `(nlat,nlon)`, with longitude varying fastest. Winds are in m/s and height is in metres. The T85 and T213 configurations are specified in the reproduction guide.

The raw solver runs contain `none.npz`, `cones.npz` and `controls.npz`, plus JSON sidecars. `prepare_fields.py` maps these to `main.npz`, `cones20.npz` and `controls20.npz` in `build/fields/`:

| Array | Shape | Meaning |
|---|---|---|
| `t_hours` | `(nt,)` | Actual saved time in hours, float64 |
| `X` | `(3,128,128)` | Cartesian unit-sphere positions, float32 |
| `weight` | `(128,128)` | Numerical area weights, float32 |
| `h` | `(nt,128,128)` | Saved height in metres, float32 |
| `err_vs_T85` | `(nt,128,128)` | Signed height difference from T85dt45 in metres, float32 |
| `err_vs_T213` | `(nt,128,128)` | Signed height difference from T213dt22 in metres, float32 |
| `relax_weight` | `(128,128)` | Actual correction W, float32; zero for main |
| `centres` | `(4,3)` | Cone or control unit vectors, float64 |

The bundled [run sidecars](../data/field-runs/main.json) retain full-precision T85 endpoint scores, timestep and correction area. Their T213 score objects are empty; the checked T213 scores in [visual-runs.json](../data/visual-runs.json) come from the prepared field arrays. Fresh runs also retain the scheme, public solver revision, reference hashes and stepping segments. Full-precision scores are distinct from quantities reconstructed from float32 fields. The pack builder recomputes errors using the saved float32 heights and positions, then stores those errors as float32. It keeps the original height for consistency checks. No image is used to calculate a norm.

Native T213 vorticity is stored in `t213-vorticity-000-056.npz`, `t213-vorticity-057-112.npz` and `t213-vorticity-113-168.npz`. Each contains `t_hours`, ascending `lat_rad`, `lon_rad = 2*pi*arange(640)/640` and `zeta_s_inv` with shape `(nframes,320,640)`. The derivative is calculated in float64 and stored in float32. All 169 reference samples from 0 to 336 hours are retained; each chunk is below 50 MB.

`assets.json` records hashes and byte sizes of the local products. `validation.json` records reference hashes, input schemas, pair consistency, discrepancies from archived scores, actual times and vorticity checks. Render checks add triangulation coverage, scales and time alignment under `build/qa/`.

[field-assets.json](../data/field-assets.json) records computed SHA-256 hashes and byte sizes for all six NPZ files, about 114 MB in total. They are not committed or hosted: regenerate them with [REPRODUCE.md](../REPRODUCE.md). The bundled [validation.json](../data/validation.json) and [figure manifest](../figures/manifest.json) retain derivative, pair, geometry and time-alignment checks.

## Regional scores

[field-metrics.csv](../data/field-metrics.csv) contains 3,042 rows: three runs, two references, 169 actual times and three regions. These are derived measurements from the visual extension, separate from the archived 12-configuration study.

Rows identify `run`, `reference`, actual `t_hours` and `region`. Runs are `main`, `cones20` and `controls20`. References are `T85dt45` and `T213dt22`. Regions are global, the cone Gaussian 20 degree weight and the control Gaussian 20 degree weight. Let e be signed height error, w the original numerical weights and q the regional weight, with q = 1 globally:

```text
area_share               = sum(w*q) / sum(w)
weighted_error_sq_m2     = sum(w*q*e^2) / sum(w)
reference_height_sq_m2   = sum(w*q*h_ref^2) / sum(w)
l2                       = sqrt(weighted_error_sq_m2 / reference_height_sq_m2)
rms_m                    = sqrt(sum(w*q*e^2) / sum(w*q))
error_energy_share       = sum(w*q*e^2) / sum(w*e^2)
enrichment               = error_energy_share / area_share
```

The regional Gaussians include their unthresholded tails and overlap. Their energy shares cannot be added as a partition of the sphere. Enrichment above one indicates more squared error per weighted area than the global average; it does not locate an error source. Zero denominators produce missing values, encoded as `nan` in generated CSV files. The pack builder retains its ten-significant-digit CSV formatting; endpoint JSON scores carry greater precision.

For the unrelaxed main run against T213:

| Day | Cone enrichment | Control enrichment |
|---:|---:|---:|
| 3 | 0.78 | 1.04 |
| 7 | 0.67 | 0.83 |
| 10 | 1.10 | 1.15 |
| 14 | 1.37 | 1.34 |

The squared-error share is broadly comparable to the area share in both sets of regions. The field does not show error concentrated at four cones; the intervention identifies where correction pays off most, not where error is largest.

Against T213, the cone/control reductions are 3.94x/2.52x at day 10 and 4.38x/2.11x at day 14. [visual-runs.json](../data/visual-runs.json) records the scores, ratios and their precision; the archived 5.3x and 2.8x numbers belong only to T85 scoring. The renderer's reconstruction from float32 fields differs from the full-precision T85 sidecar scores by at most `5.125e-6` relative across the checked endpoints. The separate main day-14 endpoint record agrees to about `2.4e-7` relative.

## Double-cover error

The partner of `(i,j)` is `((-1-i) mod M,(-1-j) mod M)`. The positions must coincide within stored precision. One representative per pair supplies the geometry, but both signed errors enter the displayed magnitude:

```text
E_pair = sqrt((w_i*e_i^2 + w_p*e_p^2) / (w_i+w_p))
```

This is a nonnegative error magnitude in metres. Averaging the signed errors first could cancel a disagreement between the two computed copies. Pair differences remain recoverable from the retained signed arrays. The metrics table reports

```text
pair_rms_difference_m = sqrt(sum((w_i+w_p)*(e_i-e_p)^2) / sum(w_i+w_p))
pair_max_difference_m = max(abs(e_i-e_p))
```

In all three saved runs, paired heights and signed errors agree exactly at every snapshot; positions and weights also agree exactly. Each sum counts a pair once. These are global consistency diagnostics, repeated unchanged across the three region rows for a given run, reference and time. Position and weight discrepancies are reported separately. The display geometry normalises the float32 position vectors to unit length before triangulation.

Every unique sphere point enters an outward convex-hull triangulation. Error values are interpolated linearly on its planar facets using radial barycentric coordinates. The triangulation must retain every point and satisfy the spherical face count `F = 2*N - 4`. Edge lengths and nearest-neighbour distances are recorded to expose long connections across a fold. No spatial smoothing is applied. Two by two subpixel samples provide raster antialiasing, not additional flow resolution.

## Relative vorticity

For longitude lambda, latitude phi and radius a:

```text
zeta = [d(v)/d(lambda) - d(u*cos(phi))/d(phi)] / [a*cos(phi)]
```

The radius comes from the pinned solver: a = 6.37122e6 m. `vort.py` uses a periodic Fourier derivative in longitude and a cubic-spline latitude derivative of `u*cos(phi)`. Eight Gaussian rows are extended across each pole with a half-turn in longitude. The meridional continuation changes the signs of both u and cos(phi), so their product has even parity. This is relative vorticity; no planetary rotation term or zonal-mean subtraction is added.

On analytic TC6 initial winds at the T213 grid, the checked RMS error is `1.3e-10 s^-1`, with maximum `2.4e-9 s^-1`, against a field maximum of `7.5e-5 s^-1`. These are derivative-validation errors on that initial field, not uncertainty bounds for the later simulation.

The analytic curl used for this check is

```text
zeta_exact = 2*omega*sin(phi)
             - K*sin(phi)*cos(phi)^R*(R^2 + 3*R + 2)*cos(R*lambda)
omega = K = 7.848e-6 s^-1, R = 4
```

The preparation script compares the spline and finite-difference results at days 0, 7, 10 and 14, including absolute latitudes above 80 degrees. It also distinguishes the analytic initial wind from the stored STSWM initial wind. At day 7 the maximum difference between the two derivative versions is `2.9e-6 s^-1`, compared with a field maximum of `1.0e-4 s^-1`. The full checks, including the polar maxima, are in [validation.json](../data/validation.json).

The renderer samples vorticity on its native Gaussian grid with periodic bilinear interpolation in latitude and longitude. Poleward of the outermost Gaussian ring, it blends linearly in polar distance to that ring's mean at the pole. This defines a continuous, single-valued scalar pole without inserting a zero cap. The reference is never downsampled to the M = 128 solver grid.

## Time and scales

Each run has 169 snapshots. They use the first solver step crossing each two-hour threshold, up to 0.2013 h later; saved day endpoints are exact. Display times come from the main run. Other runs and the reference vorticity are linearly aligned between their bracketing samples, with indices and weights retained. No extrapolation is allowed. The hero poster uses the exact 240-hour state.

Each movie has 576 frames at 24 frames/s. Frames 0 to 23 hold day 0. Frames 24 to 527 map linearly from day 0 to day 14, inclusive. Frames 528 to 575 hold day 14. During the middle segment, each frame holds the latest aligned saved snapshot, and the cursor marks that snapshot's actual time. There is no flow morphing, optical interpolation or fade back to day 0.

All height-error panels share one `magma` scale:

```text
colour = asinh(E_pair / 1 m) / asinh(E_max / 1 m)
```

E_max is the next 1, 2 or 5 times a power of ten at or above the largest saved paired error across all three runs, both references and all times. The delivered error limit is 500 m, above the measured maximum of 416.6555 m in main against T213. Zero remains visible. The 1 m linear width and fixed maximum apply to both movies, the GIF and the hero poster; there is no percentile clipping or per-frame brightening. Vorticity uses `RdBu_r` with fixed symmetric linear limits of +/-1e-4 s^-1, rounded upwards by the same rule from the measured maximum magnitude of `9.9877e-5 s^-1`.

Unexpected non-finite samples stop rendering. The low-level colour routine can depict explicitly permitted missing samples with a distinct grey hatch; the production movie path does not silently substitute zeros or discard a run.

## Views and labels

The hero uses four orthographic globes: T213 relative vorticity above main-scheme height error against T213. Each row looks towards `(longitude,latitude) = (0,35.2643897)` and its antipode `(180,-35.2643897)`, with geographic north upwards. Together they show the sphere, with two cones inside each visible disc. Cameras are fixed; graticules are spaced every 30 degrees. The movie is 1440 by 1440 pixels and its poster 2880 by 2880.

The spatial intervention uses a three-row, two-column array of complete-sphere Mollweide maps. Rows are main, cones and control; columns score against T85 and T213. Every panel has central longitude 60 degrees east and seam longitude -120 degrees. Seam markers appear at both boundaries where needed, without duplicating field contributions. The movie is 1600 by 1400 pixels; its inline GIF is 1000 by 875 and samples every eight simulated hours using the [GIF recipe](../REPRODUCE.md#regenerate-only-the-inline-gif).

Hollow diamonds mark cones; hollow circles mark controls. Corrected rows show the actual W field's `exp(-1)` contour, labelled `W = 0.368`, with solid lines for cones and dashed lines for controls. A thinner dotted contour marks the active-region boundary `W = 0.001`. The uncorrected row has no correction contour. All maps use a white background without surface lighting.

Text inside the movies is restricted to colourbar units, time ticks and minimal marker, row, column and contour labels. Static charts likewise contain axes, ticks and legends only. They use distinct markers and line styles as well as colour, and italic l in the l2 norm. Captions and interpretation belong outside the images. The CSV tables supply accessible numerical alternatives to the plots.

All visuals use DejaVu Sans. The delivered movies are H.264, `yuv420p`, 24 seconds at 24 frames/s, with BT.709 colour metadata and no audio. `hero.mp4` is 8,775,486 bytes; `relaxation-fields.mp4` is 8,565,646 bytes. [The figure guide](../figures/README.md) gives dimensions and formats, and [the manifest](../figures/manifest.json) records hashes, scales, encoding and actual frame times.

## Reading the visible structure

The hero's height error forms broad patches reaching 100 to 300 m late in the run. Thin dark lines follow crossings through zero signed error. They are not seams or missing data. Faint diagonal stripes about two grid cells wide in low-error regions are lattice-scale ripple, most visible beside the longitude-0 cone. The checked ripple estimate is about 1.5 m RMS, about 7% of error within 10 degrees of a cone versus 5% farther away; this small component is not the broad late error. Pale polar regions in vorticity reflect the wave's cos^4(latitude) amplitude, rather than inserted blank caps.

In the correction movie the corrected rows are initially brighter than main, then cones become darkest and controls intermediate. Both columns and all times share the same height-error scale. The inspected frames include the beginning, days 7 and 10, the endpoint and intermediate states; the [visual-check commands](../REPRODUCE.md#regenerate-the-fields-and-visuals) regenerate sample sheets in `build/qa/movie-frames/`. The projection check includes all 8,192 unique points and 16,380 triangles, with no non-finite field samples or unequal paired positions.
