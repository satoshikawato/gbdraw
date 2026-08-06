# Web Gallery operation screenshot register

Last updated: 2026-08-05

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

## Worked-example final previews

| Tutorial | Operation media | Decision | Required capture state | Status |
| --- | --- | --- | --- | --- |
| `HmmtDNA_basic_circular` | `manual-04-01-final-preview.webp` | Add | Exact current session; labeled features; GC content and skew; ticks; metadata; right legend; tight SVG crop | Captured and visually verified (3072×1944, DSF 3, quality 94); replaces the card-thumbnail reference only in the tutorial |
| `lambda_basic_linear` | `manual-05-01-final-preview.webp` | Add | Exact current session; both strand lanes; all labels; ruler; metadata; left legend; tight SVG crop | Captured and visually verified (4182×1440, DSF 3, quality 94); full ruler and rightmost labels remain in frame |
| `hepatoplasmataceae_collinear` | `manual-06-01-collinear-overview.webp` | Replace | Exact current session; five records; both GC tracks; orientation-and-identity blocks; ruler; right legend without overlap | Captured and visually verified (4182×1278, DSF 3, quality 94); same-size comparison confirms the legend no longer overlaps the GC tracks |
| `majanivirus_orthogroup` | `manual-07-01-orthogroup-preview.webp` | Recrop | Exact current session; nine labeled records; product colors; similarity-group ribbons; right legend; SVG-aspect crop | Captured and visually verified (4182×1020, DSF 3, quality 94); same-size comparison confirms removal of the letterboxed app chrome |
| `tobacco-chloroplast` | `manual-08-01-chloroplast-preview.webp` | Replace | Exact current session; LSC/IRb/SSC/IRa brackets; radial labels; GC track; metadata; one entry per legend category | Captured and visually verified (3072×2250, DSF 3, quality 94); same-size comparison confirms the duplicate legend category is gone |

Before these captures, the Gallery refresh path synchronizes the legacy
`config.form.legend` control with the canonical render-request legend. The four
affected sessions (`HmmtDNA_basic_circular`, both Hepatoplasmataceae examples,
and `majanivirus_orthogroup`) must restore the same legend position that their
saved result renders.

The refresh path also copies canonical Circular track slots into the restored
Web draft. This repairs the chloroplast session's region-annotation track state,
so the documented session opens successfully before the final preview is
captured.

## Diagram-layout overhaul recapture (WP7)

The pre-capture contract is
`docs/WEB_GALLERY_DIAGRAM_LAYOUT_RECAPTURE_PLAN.md`. Documentation images are
full-viewport PNGs owned by `docs/capture/run_all.py`; Gallery media are compact
operation/result WebP crops owned by tutorial capture metadata. The initial
documentation owner pass replaced 18 result PNGs and kept seven input/settings
PNGs byte-stable. Six accepted outputs cover `T-GUI-01` and `H-GUI-03`; the
other 12 need the corrected recapture described below. All six Gallery result
WebPs remain pending.

| Scenario | Documentation outputs | Gallery relation | Required state | Decision | Status |
| --- | --- | --- | --- | --- | --- |
| `T-GUI-01` | Six files under `docs/images/t-gui-01/` | Exact: `HmmtDNA_basic_circular/manual-04-01-final-preview.webp` | Fresh HmmtDNA flow for docs; exact `HmmtDNA_basic_circular` session for Gallery; Circular Middle, Out labels, right legend | Re-run all owner outputs; expect input-only PNG unchanged; recapture generated-result/full-preview views | Documentation pass accepted: five results replaced, input view unchanged; Gallery WebP pending |
| `T-GUI-02` | Four files under `docs/images/t-gui-02/` | Exact: `lambda_basic_linear/manual-05-01-final-preview.webp` | Fresh Lambda flow for docs; exact `lambda_basic_linear` session for Gallery; Linear, no comparison, all labels, ruler, left legend | Re-run all owner outputs; expect input-only PNG unchanged; recapture generated-result/full-preview views | Three result PNGs replaced provisionally; corrected recapture and Gallery WebP pending |
| `T-GUI-05` | Five files under `docs/images/t-gui-05/` | Exact: `tobacco-chloroplast/manual-08-01-chloroplast-preview.webp` | Fresh NC_001879.2 flow for docs; exact `tobacco-chloroplast` session for Gallery; three-slot stack, four region labels, upper-left legend | Re-run all owner outputs and recapture the exact-session Gallery result | Four result PNGs replaced provisionally; corrected recapture and Gallery WebP pending |
| `H-GUI-02` | `grid-settings.png`, `grid-result.png` | Contextual only: `Vnig_TUMSAT-TG-2018/manual-06-01-multirecord-preview.webp` | Docs use four complete mitochondrial records in a 2×2 equal-size grid; Gallery uses six Vibrio replicons with Auto sizing, left legend, bottom title | Keep distinction explicit; add exact-session metadata for the Gallery result before recapture | Settings unchanged; result PNG replaced provisionally; corrected recapture and Gallery WebP pending |
| `H-GUI-03` | `record-layout.png`, `orientation-result.png` | Contextual only: `vibrio-harveyi-group-collinear/manual-08-01-collinear-overview.webp` | Docs use two comparison-free phage rows at 24 px; Gallery uses 11 comparison records in five rows at 48 px with bottom legend | Recapture docs result; strengthen Gallery app-state/visible-text assertions before its result recapture | Documentation pass accepted: result replaced, settings unchanged; Gallery WebP pending |
| `H-GUI-09` | `track-settings.png`, `track-result.png` | Contextual only: `HmmtDNA_ATskew/manual-09-01-atskew-preview.webp` | Docs use AP027133 depth + GC content/skew; Gallery uses HmmtDNA GC/AT skew without depth | Recapture docs result; add exact-session metadata for the manual-only Gallery result | Settings unchanged; result PNG replaced provisionally; corrected recapture and Gallery WebP pending |
| `H-GUI-10` | `slot-settings.png`, `annotation-result.png` | Contextual only: tobacco chloroplast final preview | Docs use alternating annotation lanes, an outside annotation slot, AT skew, top title, right legend; Gallery uses one inside region lane, no AT skew, no title, upper-left legend | Recapture docs result; reuse the T-GUI-05 Gallery recapture only as contextual coverage | Settings unchanged; result PNG replaced provisionally; corrected recapture pending |
| `H-GUI-11` | `style-settings.png`, `style-result.png` | Contextual only: HmmtDNA basic final preview | Docs use soft_pastels, whitelist, selected labels, top title, and exact legend order; Gallery basic example has a different style/title state | Recapture both docs views; reuse the T-GUI-01 Gallery recapture only as contextual right-legend coverage | Both PNGs replaced provisionally; corrected recapture pending |

The documentation environment is pinned to 1440×900 at DSF 1 with the light
theme. Planned Gallery result crops use their declared viewports, DSF 3, WebP
quality 94, and the exact selectors/actions/captions/alt text in the recapture
plan. Existing compact control-operation WebPs remain `Keep` unless strict
capture validation proves the visible control or selected value is stale.
