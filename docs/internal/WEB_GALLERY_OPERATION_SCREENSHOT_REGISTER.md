# Web Gallery operation screenshot register

Last updated: 2026-08-11

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

## Linear comparison workflow (Work package B)

This section supersedes the #316 comparison-panel and row-boundary entries. The
saved sessions remain the operation data source. Capture actions may activate a
comparison plan for the screenshot; they do not rewrite the session file.

| Tutorial | Operation media | Decision | Required capture state | Status |
| --- | --- | --- | --- | --- |
| `lambda_basic_linear` | `manual-02-03-no-comparison.webp` | Recapture | Exact Lambda session; **No comparison** command; **Current: No comparison** status; closed **Settings** and **Selected pairs (0)** | Captured at DSF 3; strict validation passed; visually accepted |
| `BGC0000708-BGC0000713` | `manual-02-01-record-text-row.webp`, `manual-04-02-reverse-bgc0000713.webp` | Recapture | Exact BGC session; uploader before an open **Record options** disclosure; row-specific organism, subtitle, region, and reverse-complement controls | Captured at DSF 3; strict validation passed; visually accepted |
| `BGC0000708-BGC0000713` | `manual-03-01-open-pairwise.webp`, `manual-03-02-select-losatp-orthogroups.webp` | Recapture | Exact BGC session; **Run LOSAT** command and current status; open **Settings**; all three **LOSAT Mode** buttons with **LOSATP** pressed; open **LOSATP mode** menu in UI order with **Similarity groups** selected; tutorial filters | Button UI captured at DSF 3; strict validation passed; visually accepted |
| `BGC0000708-BGC0000713` | `manual-03-02-runtime-reproducibility.webp`, `manual-03-04-first-raw-result.webp` | Add | Open **Advanced comparison and layout**; capture **Runtime and reproducibility** separately from the first pair's raw filename, status, and download | Captured at DSF 3; strict validation passed; visually accepted |
| `BGC0000708-BGC0000713` | `manual-03-03-first-comparison-boundary.webp` | Recapture | Open **Selected pairs (4)**; first `record-1->record-2` pair; LOSAT source and endpoint identities; no raw-result owner in this crop | Captured at DSF 3; strict validation passed; visually accepted |
| `hepatoplasmataceae_orthogroup` | `manual-02-01-upload-row-context.webp` | Recapture | Exact orthogroup session; first uploader and open **Record options** disclosure | Captured at DSF 3; strict validation passed; visually accepted |
| `hepatoplasmataceae_orthogroup` | `manual-03-01-browser-losat.webp`, `manual-04-01-orthogroups-mode.webp`, `manual-04-02-orthogroup-settings.webp` | Recapture | **Run LOSAT** command and status; open **Settings**; all three **LOSAT Mode** buttons with **LOSATP** pressed; open **LOSATP mode** menu in UI order with **Similarity groups** selected; tutorial filters | Button UI captured at DSF 3; strict validation passed; visually accepted |
| `hepatoplasmataceae_orthogroup` | `manual-03-02-runtime-reproducibility.webp` | Add | Open **Advanced comparison and layout**; Auto execution, Safe total threads, and automatic run allocation | Captured at DSF 3; strict validation passed; visually accepted |
| `hepatoplasmataceae_collinear` | `manual-02-01-upload-row-context.webp` | Recapture | Exact collinear session; first uploader and open **Record options** disclosure | Captured at DSF 3; strict validation passed; visually accepted |
| `hepatoplasmataceae_collinear` | `manual-03-01-open-pairwise.webp`, `manual-04-01-collinear-reduction.webp`, `manual-04-02-orientation-identity.webp` | Recapture | **Run LOSAT** command and status; open **Settings**; all three **LOSAT Mode** buttons with **LOSATP** pressed; open **LOSATP mode** menu in UI order with **Collinear blocks** selected; **Color mode: Orientation + identity** | Button UI captured at DSF 3; strict validation passed; visually accepted |
| `hepatoplasmataceae_collinear` | `manual-03-02-runtime-reproducibility.webp`, `manual-04-03-advanced-collinear.webp` | Add | Open **Advanced comparison and layout**; capture runtime and **Advanced collinear search** as separate operations | Captured at DSF 3; strict validation passed; visually accepted |
| `majanivirus_orthogroup` | `manual-02-01-upload-row-label-context.webp` | Recapture | Exact Majanivirus session; first uploader and open **Record options** with organism text | Captured at DSF 3; strict validation passed; visually accepted |
| `majanivirus_orthogroup` | `manual-03-01-open-pairwise.webp`, `manual-03-02-losatp-orthogroups.webp` | Recapture | **Run LOSAT** command and status; open **Settings**; all three **LOSAT Mode** buttons with **LOSATP** pressed; open **LOSATP mode** menu in UI order with **Similarity groups** selected; tutorial filters | Button UI captured at DSF 3; strict validation passed; visually accepted |
| `majanivirus_orthogroup` | `manual-03-03-runtime-reproducibility.webp` | Add | Open **Advanced comparison and layout** for the requested thread count | Captured at DSF 3; strict validation passed; visually accepted |
| `vibrio-harveyi-group-collinear` | `manual-02-01-record-row.webp`, `manual-03-01-record-layout.webp` | Recapture | Exact Vibrio session; open per-record **Record options** and top-level **Advanced comparison and layout > Record Layout** | Captured at DSF 3; strict validation passed; visually accepted |
| `vibrio-harveyi-group-collinear` | `manual-04-00-run-adjacent-losat.webp`, `manual-04-01-search-method-collinear.webp`, `manual-04-02-color-mode-orientation-identity.webp`, `manual-04-03-adjacent-pairs.webp` | Add/recapture | Capture action only: activate all-adjacent LOSAT; command/status; all three **LOSAT Mode** buttons with **LOSATP** pressed; open **LOSATP mode** menu in UI order with **Collinear blocks** selected; color and evidence scope | Button UI captured at DSF 3; strict validation passed; visually accepted |
| `vibrio-harveyi-group-collinear` | `manual-04-04-advanced-collinear.webp` | Add | Capture action only: activate all-adjacent LOSAT; open **Advanced comparison and layout > Advanced collinear search** | Captured at DSF 3; strict validation passed; visually accepted |

The Vibrio editable session stays CLI-only with
`linearComparisonPlan.mode = none`. Every data-dependent comparison operation
asserts `adjacent + losat` after its capture action. **LOSAT Mode** is a
three-button group: LOSATN, LOSATP, and TLOSATX. When LOSATP is active,
**LOSATP mode** lists Similarity groups, Collinear blocks, and Pairwise matches.
Similarity groups always runs
an all-vs-all protein search. The collinear tutorials retain their explicit
saved **Adjacent pairs** evidence scope.

## Worked-example final previews

| Tutorial | Operation media | Decision | Required capture state | Status |
| --- | --- | --- | --- | --- |
| `HmmtDNA_basic_circular` | `manual-04-01-final-preview.webp` | Add | Exact current session; labeled features; GC content and skew; ticks; metadata; right legend; tight SVG crop | Captured and visually verified (3072×2049, DSF 3, quality 94); replaces the card-thumbnail reference only in the tutorial |
| `lambda_basic_linear` | `manual-05-01-final-preview.webp` | Add | Exact current session; both strand lanes; all labels; ruler; metadata; left legend; tight SVG crop | Captured and visually verified (4182×1452, DSF 3, quality 94); full ruler and rightmost labels remain in frame |
| `hepatoplasmataceae_collinear` | `manual-06-01-collinear-overview.webp` | Replace | Exact current session; five records; both GC tracks; orientation-and-identity blocks; ruler; right legend without overlap | Captured and visually verified (4182×1278, DSF 3, quality 94); same-size comparison confirms the legend no longer overlaps the GC tracks |
| `majanivirus_orthogroup` | `manual-07-01-orthogroup-preview.webp` | Recrop | Exact current session; nine labeled records; product colors; similarity-group ribbons; right legend; SVG-aspect crop | Captured and visually verified (4182×1020, DSF 3, quality 94); same-size comparison confirms removal of the letterboxed app chrome |
| `tobacco-chloroplast` | `manual-08-01-chloroplast-preview.webp` | Replace | Exact current session; LSC/IRb/SSC/IRa brackets; radial labels; GC track; metadata; one entry per legend category | Captured and visually verified (3072×2187, DSF 3, quality 94); same-size comparison confirms the duplicate legend category is gone |

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
`docs/internal/WEB_GALLERY_DIAGRAM_LAYOUT_RECAPTURE_PLAN.md`. Documentation images are
full-viewport PNGs owned by `docs/capture/run_all.py`; Gallery media are compact
operation/result WebP crops owned by tutorial capture metadata. The initial
documentation owner pass was followed by corrected recaptures and exact checks.
All six linked Gallery result WebPs were refreshed from their exact sessions.

| Scenario | Documentation outputs | Gallery relation | Required state | Decision | Status |
| --- | --- | --- | --- | --- | --- |
| `T-GUI-01` | Six files under `docs/images/t-gui-01/` | Exact: `HmmtDNA_basic_circular/manual-04-01-final-preview.webp` | Fresh HmmtDNA flow for docs; exact `HmmtDNA_basic_circular` session for Gallery; Circular Middle, Out labels, right legend | Re-run all owner outputs; expect input-only PNG unchanged; recapture generated-result/full-preview views | Documentation images accepted; exact-session Gallery WebP recaptured at 3072×2049 |
| `T-GUI-02` | Four files under `docs/images/t-gui-02/` | Exact: `lambda_basic_linear/manual-05-01-final-preview.webp` | Fresh Lambda flow for docs; exact `lambda_basic_linear` session for Gallery; Linear, no comparison, all labels, ruler, left legend | Re-run all owner outputs; expect input-only PNG unchanged; recapture generated-result/full-preview views | Corrected documentation images accepted; exact-session Gallery WebP recaptured at 4182×1452 |
| `T-GUI-05` | Five files under `docs/images/t-gui-05/` | Exact: `tobacco-chloroplast/manual-08-01-chloroplast-preview.webp` | Fresh NC_001879.2 flow for docs; exact `tobacco-chloroplast` session for Gallery; three-slot stack, four region labels, upper-left legend | Re-run all owner outputs and recapture the exact-session Gallery result | Corrected documentation images accepted; exact-session Gallery WebP recaptured at 3072×2187 |
| `H-GUI-02` | `grid-settings.png`, `grid-result.png` | Contextual only: `Vnig_TUMSAT-TG-2018/manual-06-01-multirecord-preview.webp` | Docs use four complete mitochondrial records in a 2×2 equal-size grid; Gallery uses six Vibrio replicons with Auto sizing, left legend, bottom title | Keep distinction explicit; add exact-session metadata for the Gallery result before recapture | Documentation images accepted; contextual Gallery WebP recaptured at 3072×2016 |
| `H-GUI-03` | `record-layout.png`, `orientation-result.png` | Contextual only: `vibrio-harveyi-group-collinear/manual-08-01-collinear-overview.webp` | Docs use two comparison-free phage rows at 24 px; Gallery uses 11 comparison records in five rows at 48 px with bottom legend | Recapture docs result; strengthen Gallery app-state/visible-text assertions before its result recapture | Documentation images accepted; contextual Gallery WebP recaptured at 3804×1200 |
| `H-GUI-09` | `track-settings.png`, `track-result.png` | Contextual only: `HmmtDNA_ATskew/manual-09-01-atskew-preview.webp` | Docs use AP027133 depth + GC content/skew; Gallery uses HmmtDNA GC/AT skew without depth | Recapture docs result; add exact-session metadata for the manual-only Gallery result | Documentation result accepted at 70%; contextual Gallery WebP recaptured at 3072×2304 |
| `H-GUI-10` | `slot-settings.png`, `annotation-result.png` | Contextual only: tobacco chloroplast final preview | Docs use alternating annotation lanes, an outside annotation slot, AT skew, top title, right legend; Gallery uses one inside region lane, no AT skew, no title, upper-left legend | Recapture docs result; reuse the T-GUI-05 Gallery recapture only as contextual coverage | Documentation result accepted at 70%; contextual Gallery coverage verified |
| `H-GUI-11` | `style-settings.png`, `style-result.png` | Contextual only: HmmtDNA basic final preview | Docs use soft_pastels, whitelist, selected labels, top title, and exact legend order; Gallery basic example has a different style/title state | Recapture both docs views; reuse the T-GUI-01 Gallery recapture only as contextual right-legend coverage | Documentation images accepted at the largest complete-fit scale; contextual Gallery coverage verified |

The documentation environment is pinned to 1440×900 at DSF 1 with the light
theme. Planned Gallery result crops use their declared viewports, DSF 3, WebP
quality 94, and the exact selectors/actions/captions/alt text in the recapture
plan. Existing compact control-operation WebPs remain `Keep` unless strict
capture validation proves the visible control or selected value is stale.

## Circular screenshot scale audit

- Use 70% whenever the complete title, plot, labels, and legend remain visible.
  This applies to `T-GUI-10`, `T-GUI-12`, `H-GUI-09`, and `H-GUI-10`, plus
  the intermediate `T-GUI-05` states.
- Use 60% only where 70% clips required content. This applies to `T-GUI-01`,
  `T-GUI-09`, `H-GUI-11`, `H-GUI-12`, `H-GUI-13`, `H-GUI-14`, and
  `H-GUI-15`.
- Use 50% only for the dense final tobacco chloroplast figure (`T-GUI-05`) and
  the three-comparison-ring figures (`T-GUI-06` and `H-GUI-06`), where 60%
  clips the title, labels, legend, or outer comparison ring.
- The HmmtDNA feature-highlight result uses `Middle`, strand separation off,
  70%, and gene labels for all 13 mitochondrial CDS features, including
  `COX1`.
