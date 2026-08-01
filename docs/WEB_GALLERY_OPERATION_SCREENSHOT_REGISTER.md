# Web Gallery operation screenshot register

Last updated: 2026-08-01

This register records task-specific decisions for Gallery operation media.
Capture metadata remains the executable source of truth in each tutorial JSON.

## Coordinate-scale visibility (#311 and #315)

| Tutorial | Operation media | Decision | Required capture state | Status |
| --- | --- | --- | --- | --- |
| `lambda_basic_linear` | `manual-03-03-scale-style-ruler.webp` | Recapture | Exact Lambda session; Linear mode; **Show Coordinate Scale** selected; **Scale Style** set to **Ruler (Ticks)**; both controls inside the crop | Recaptured and visually verified (990×768, DSF 3, quality 94) |
| `hepatoplasmataceae_orthogroup` | `manual-05-01-layout-gc-skew-ruler.webp` | Recapture | Exact orthogroup session; Linear mode; **Show Coordinate Scale** selected; **Scale Style** set to **Ruler (Ticks)**; both controls inside the crop | Recaptured and visually verified (990×768, DSF 3, quality 94) |
| `hepatoplasmataceae_collinear` | `manual-05-01-layout-gc-skew-ruler.webp` | Recapture | Exact collinear session; Linear mode; **Show Coordinate Scale** selected; **Scale Style** set to **Ruler (Ticks)**; both controls inside the crop | Recaptured and visually verified (990×768, DSF 3, quality 94) |
| `BGC0000708-BGC0000713` | `manual-05-02-scale-style-ruler.webp` | Recapture | Exact BGC session; Linear mode; **Show Coordinate Scale** selected; **Scale Style** set to **Ruler (Ticks)**; both controls inside the crop | Recaptured and visually verified (990×768, DSF 3, quality 94) |
| `HmmtDNA_ATskew` | `manual-07-01-tick-track-context.webp` | Keep | Exact mitochondrial session; custom stack enabled; final slot is enabled `ticks` with `label_in_tick_out` | Current image remains truthful; deterministic metadata added |
| `vibrio-harveyi-group-collinear` | `manual-06-01-linear-layout.webp` | Keep | Layout card only; the settings table retains the visible ruler configuration | No Axis & Scale controls in the crop |
| `WSSV_genome_comparison` | `manual-03-01-basic-circular-layout.webp` | Keep | Layout card only; the settings table retains the 20,000 bp interval and the default visible scale | No Axis & Scale controls in the crop |

The four recaptures replace stale crops from the pre-checkbox Axis & Scale
card. Final previews and thumbnails remain unchanged because coordinate scales
remain visible by default.

## Explicit Linear comparison plan (#316)

| Tutorial | Operation media | Decision | Required capture state | Status |
| --- | --- | --- | --- | --- |
| `lambda_basic_linear` | `manual-02-03-no-comparison.webp` | Add | Exact Lambda session; Linear; **No comparison** selected; focused global comparison choices | Captured and visually verified (1008×300, DSF 3, quality 94) |
| `BGC0000708-BGC0000713` | `manual-02-01-record-text-row.webp` | Recapture | Exact BGC session; row 1 GenBank, organism, and subtitle controls; omit the retired per-record LOSAT filename area | Recaptured and visually verified (1008×1263, DSF 3, quality 94) |
| `BGC0000708-BGC0000713` | `manual-03-01-open-pairwise.webp` | Recapture | Exact BGC session; Linear; **Run LOSAT**; LOSATP; Similarity groups; thread and threshold values from the tutorial table | Recaptured and visually verified (1068×1686, DSF 3, quality 94) |
| `hepatoplasmataceae_orthogroup` | `manual-02-01-upload-row-context.webp` | Recapture | Exact orthogroup session; row 1 GenBank and neighboring record controls; omit the retired per-record LOSAT filename area | Recaptured and visually verified (1008×1263, DSF 3, quality 94) |
| `hepatoplasmataceae_orthogroup` | `manual-03-01-browser-losat.webp` | Recapture | Exact orthogroup session; Linear; **Run LOSAT**; LOSATP; comparison selector and LOSAT settings both visible | Recaptured and visually verified (996×1686, DSF 3, quality 94) |
| `majanivirus_orthogroup` | `manual-02-01-upload-row-label-context.webp` | Recapture | Exact Majanivirus session; row 1 GenBank and organism label; omit the retired per-record LOSAT filename area | Recaptured and visually verified (1008×1263, DSF 3, quality 94) |
| `majanivirus_orthogroup` | `manual-03-01-open-pairwise.webp` | Recapture | Exact Majanivirus session; `adjacent + losat`; LOSATP Similarity groups; requested 32 threads per run; tutorial thresholds | Recaptured and visually verified (996×1686, DSF 3, quality 94) |
| `hepatoplasmataceae_collinear` | `manual-02-01-upload-row-context.webp` | Recapture | Exact collinear session; row 1 GenBank and neighboring record controls; omit the retired per-record LOSAT filename area | Recaptured and visually verified (1008×1263, DSF 3, quality 94) |
| `hepatoplasmataceae_collinear` | `manual-03-01-open-pairwise.webp` | Replace | Convert the stale manual crop to a declarative exact-session capture with `adjacent + losat`, LOSATP, and Collinear blocks | Recaptured and visually verified (996×2058, DSF 3, quality 94) |
| `vibrio-harveyi-group-collinear` | `manual-04-03-adjacent-pairs.webp` | Inspect | Exact Vibrio session; LOSATP Collinear blocks; Evidence scope **Adjacent pairs**; retain only if the new global selector stays outside the crop | Inspected; keep existing crop (840×483), which does not show the global selector |

The comparison-panel captures declare the exact saved session, comparison-plan
state, LOSAT program, and visible selected controls. Record-upload crops no
longer describe the retired per-record raw-result filename field.
