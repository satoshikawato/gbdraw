# Web Gallery Tutorial Audit Register

Audit date: 2026-07-10

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
| `tobacco-chloroplast` | Add a focused chloroplast region-annotation walkthrough with the four coordinate rows as a table and keep the saved-session path explicit. | Add the generated final-result thumbnail only. The input values remain copyable in the table and command; no reconstructed or redundant operation screenshots are needed. |

## Cross-tutorial cleanup

- Keep mode-selector reuse; it is intentional and readable.
- Treat a good existing screenshot as the semantic-crop baseline. Compare replacements side by side at the same displayed size, and reject wider crops or popup crops that hide the selected feature/ribbon without adding necessary information.
- Remove generic Export and Save Session strips whose controls render too small in operation cards and add no example-specific state.
- Keep tables as the authoritative source for exact repeated values; use screenshots to identify the real control or visible result.
- The final tutorial set references 80 unique operation images. Delete only media that has become unreferenced and rerun strict media validation plus desktop/mobile browser checks.

## Verification

- JSON parsing and strict capture/media validation pass for 80 operations and 80 operation images.
- The tutorial media inventory is exact: 80 referenced WebP files, 80 files present, 0 missing, and 0 unreferenced.
- The seven final contact sheets were inspected after recapture; AT skew and ticks retain their prior compact bitmaps, while RNR2 is now a shorter landscape crop with the selected l-rRNA highlighted beside the popup.
- `tests/web/gallery-tutorial.playwright.spec.js` passes all 18 Chromium desktop/mobile tests.
