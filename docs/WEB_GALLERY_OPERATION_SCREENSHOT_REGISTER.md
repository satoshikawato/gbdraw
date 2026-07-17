# Web Gallery Tutorial Audit Register

Audit date: 2026-07-17

Scope: all ready tutorials under `gbdraw/web/gallery/tutorials/`, their referenced WebP media, current web-app controls, saved Gallery sessions, and desktop/mobile Gallery rendering.

## Decisions

| Tutorial | Writing and structure | Media decision |
|---|---|---|
| `BGC0000708-BGC0000713` | Move orthogroup alignment after the first generation, use exact italic markup and current UI labels, and document the orange palette, `#dddddd` CDS fallback, and definition-line styles required by the saved figure. | Keep the compact title crop, recapture the rule rows and selected `livE`/`og_1` context, and remove the redundant threshold crop plus generic export/session/editor media. |
| `HmmtDNA_ATskew` | Use literal species markup and current UI labels, correct the Middle template to four slots, and explicitly move ticks after the added AT-skew slot. | Add one compact Middle-template crop; retain the earlier compact AT-skew and ticks crops; recrop Qualifier Priority to `CDS → gene`; recapture RNR2 with the highlighted l-rRNA fully visible beside its popup; remove generic export/session media. |
| `Vnig_TUMSAT-TG-2018` | Preserve literal italic markup and use current Multi-Record Canvas field names. | Consolidate the overlapping multi-record screenshots, recapture a compact title-control crop and the selected-feature controls, and keep the existing clicked-feature context. |
| `WSSV_genome_comparison` | Make the bundled session the exact-reproduction path, provide available sequence accessions, use current UI labels, remove Linear-only thread instructions, and combine per-ring file/label/color data. | Reuse the compact comparison-series and VP180 Details crops, move zoom guidance before feature inspection, and remove redundant per-ring, browser-LOSAT, export, and session media. |
| `hepatoplasmataceae_collinear` | Add the missing collinear overview, document all collinear thresholds, and separate result inspection from popup inspection. | Replace the Generate-button-only crop with the real overview, recapture the exact block and feature targets, add a layout overview, and remove generic legend/export/session media. |
| `hepatoplasmataceae_orthogroup` | Correct Minimum Identity from 30 to 0, describe all-record orthogroup evidence, and use current layout labels. | Merge Generate with the overview, recapture the exact orthogroup and feature targets, add a layout overview, and remove generic legend/export/session media. |
| `majanivirus_orthogroup` | Correct the intended LOSATP, layout, legend, pairwise, and user-enterable hex color values; explain that requested threads are capped by the Safe budget. | Consolidate the overlapping color-rule crops, replace stale LOSAT/final/popup states, add focused Right-legend and Pairwise Match crops, and remove generic post-generation media. |
| `_common-color-rule-guide` | Use the actual `Legend Caption` label and anchored exact-match examples. | No media needed. |
| `tobacco-chloroplast` | Replace the single-step summary with a full Circular workflow: input and supporting TSVs, core layout and label settings, the four coordinate rows, the annotation track binding, and generated-result inspection. Keep exact values in tables and explain the bundled-session path. | Replace the thumbnail-only operation with focused real-UI crops for mode, upload, basic settings, layout, labels, supporting TSVs, the four-row Region Annotations set, the three-row Custom Track Slots stack, and the completed preview. |

## Cross-tutorial cleanup

- Keep mode-selector reuse; it is intentional and readable.
- Treat a good existing screenshot as the semantic-crop baseline. Compare replacements side by side at the same displayed size, and reject wider crops or popup crops that hide the selected feature/ribbon without adding necessary information.
- Remove generic Export and Save Session strips whose controls render too small in operation cards and add no example-specific state.
- Keep tables as the authoritative source for exact repeated values; use screenshots to identify the real control or visible result.
- The current tutorial set references 94 unique operation images: 91 files under `gallery/media/` plus three intentional Gallery thumbnails. Delete only media that has become unreferenced and rerun strict media validation plus desktop/mobile browser checks.

## Verification

- JSON parsing and strict capture/media validation pass for 100 operations and 100 operation media entries.
- The Gallery media inventory is exact: 91 referenced WebP files under `gallery/media/`, 91 files present, 0 missing, and 0 unreferenced. Three operations intentionally use Gallery thumbnails.
- The tobacco chloroplast contact sheet and the full-size region, track-stack, and final-preview images were inspected after capture.
- The focused tobacco chloroplast Chromium test passes, and a 390 px mobile check loads all 11 images with no horizontal document overflow or page errors.
