# Web Gallery operation screenshot register

Last updated: 2026-10-10

This register records task-specific decisions for Gallery operation media.
Capture metadata remains the executable source of truth in each tutorial JSON.

## BGC record orientation

MIBiG stores BGC0000713 on the opposite strand from the other four clusters.
Since 2026-07-15 the BGC Gallery Session stores a reverse-complemented copy of
that file instead of the original file with **Reverse complement** on, so a
reader who starts from the MIBiG files still needs the step. The recipe reads
the records first, so the **Record** select shows its loaded options rather
than **Loading records...**.

| Tutorial | Operation media | Decision | Required capture state | Status |
| --- | --- | --- | --- | --- |
| `BGC0000708-BGC0000713` | `manual-04-02-reverse-bgc0000713.webp` | Recapture | Exact BGC Session; the BGC0000713 file card from its file name to **Region (optional)**, with **Record options** open and **Reverse complement** checked | Captured at DSF 3; accepted |
| `BGC0000708-BGC0000713` | `manual-04-04-pairwise-style-curve.webp` | Delete | None; no step referenced it | Deleted |

## BGC og_6 alignment and arrow shaft width (Gallery revision 2026-10-10)

The BGC Gallery Session now uses Shaft Width Ratio 0.6 and aligns on og_6
(neoU of BGC0000709, livU kept for BGC0000708). Preview recipes compute their
pan from the target's rectangle, so a layout change no longer leaves a stale
fixed offset. The anchors-dialog crops show the current **Select alignment
anchors** dialog; the dialog redesign recaptures `manual-08-02` and
`manual-08-03`.

| Tutorial | Operation media | Decision | Required capture state | Status |
| --- | --- | --- | --- | --- |
| `BGC0000708-BGC0000713` | `manual-04-04-arrow-shaft-width.webp` | Add | Exact BGC Session; **Features** open; **Arrow Geometry** box with Head Length Ratio at Auto and Shaft Width Ratio 0.6 highlighted | Captured at DSF 3; accepted |
| `BGC0000708-BGC0000713` | `manual-08-01-align-og1.webp` → `manual-08-01-align-og6.webp` | Replace | Clicked neoU of BGC0000709 (og_6) beside the popup with **Align…** and **Review alignment options…** | Captured at DSF 3; accepted |
| `BGC0000708-BGC0000713` | `manual-08-02-anchor-candidates.webp` | Add | **Review alignment options…** open; the BGC0000708 card with the transport protein CAG38692.1 and the selected, **Recommended** livU CAG38700.1 | Captured at DSF 3; accepted |
| `BGC0000708-BGC0000713` | `manual-08-03-alignment-direction.webp` | Add | Same dialog; **Alignment direction** with **Keep current directions** selected | Captured at DSF 3; accepted |
| `BGC0000708-BGC0000713` | `manual-08-01-bgc-preview.webp` → `manual-08-04-aligned-preview.webp` | Replace | The old image had no capture recipe; the new one is the restored Session's preview at fixed width, preview controls hidden | Captured at DSF 3; accepted |
| `BGC0000708-BGC0000713` | `manual-09-01-orthogroup-popup.webp`, `manual-10-01-feature-popup.webp` | Recapture; recipe corrected | Pans computed from the og_18 ribbon and the livE feature rectangles | Captured at DSF 3; accepted |

## LOSAT runtime controls under Comparison Settings (OV-366)

The runtime controls moved from **Advanced comparison and layout** into
**Comparison › Settings**; the recipes waited for a section that no longer
exists there. They also set the drawing's retired `losat.*` execution fields
instead of the app-level `losatExecution`.

| Tutorial | Operation media | Decision | Required capture state | Status |
| --- | --- | --- | --- | --- |
| `BGC0000708-BGC0000713`, `hepatoplasmataceae_collinear`, `hepatoplasmataceae_orthogroup` | `manual-03-02-runtime-reproducibility.webp` | Recapture; recipe and text corrected | Open **Settings**; **Runtime and reproducibility** section; Auto execution, Safe total threads, Auto threads per run | Captured at DSF 3; accepted |
| `majanivirus_orthogroup` | `manual-03-03-runtime-reproducibility.webp` | Recapture; recipe and text corrected | Same section with 32 threads per run | Captured at DSF 3; accepted |

## Comparison pressed state and Align help-tip (GUI remediation S07)

The removed **Current: …** status line and the always-on Align paragraph were
visible in these images. Each recapture was compared with the committed image
at the same display size.

| Tutorial | Operation media | Decision | Required capture state | Status |
| --- | --- | --- | --- | --- |
| `BGC0000708-BGC0000713` | `manual-08-01-align-og1.webp` | Recapture; recipe corrected | Exact BGC Session; the recipe pans with `app.canvasPan` (the old container transform no longer moves the canvas); clicked livE beside the popup; **Align…**, **Review alignment options…** and the help-tip in one row | DSF 3; livE, its highlight and label are no longer covered; accepted |
| `BGC0000708-BGC0000713`, `majanivirus_orthogroup` | `manual-03-01-open-pairwise.webp` | Recapture; taller viewport | **Run LOSAT** pressed; open **Settings** with the asserted filters; viewport 1600 × 1600 so the whole card fits without the sticky header | Taller because Settings now also shows Runtime and reproducibility and Comparison appearance; no status line; accepted |
| `hepatoplasmataceae_collinear`, `hepatoplasmataceae_orthogroup`, `vibrio-harveyi-group-collinear` | `manual-03-01-open-pairwise.webp`, `manual-03-01-browser-losat.webp`, `manual-04-00-run-adjacent-losat.webp` | Recapture; recipe corrected | **Run LOSAT** pressed; the recipe closes the Settings disclosure that the app opens for LOSAT, as the alt text states | Compact card; the old Vibrio image was cut off under the header; accepted |
| `lambda_basic_linear` | `manual-02-03-no-comparison.webp` | Recapture | **No comparison** pressed; closed **Settings** and **Selected pairs (0)** | No status line; accepted |

The other Gallery media are unchanged. The Vibrio linear-layout image never
showed the Lock explanation paragraph, and other popup images do not show the
Align area.

## Popup framing with `app.canvasPan` (GUI remediation follow-up)

These recipes panned with the preview container's `style.transform`. The
preview binds that transform to `canvasPan` and `zoom`, so the next zoom change
overwrote it and the step had no effect. Each recipe now sets `app.canvasPan`
and waits 500 ms, so the clicked target is inside the visible canvas before the
click. Each recapture was compared with the committed image at the same scale.

| Tutorial | Operation media | Decision | Required capture state | Status |
| --- | --- | --- | --- | --- |
| `BGC0000708-BGC0000713` | `manual-09-01-orthogroup-popup.webp` | Recapture; recipe corrected | Zoom 1.4; the og_18 ribbon is highlighted and the popup opens beside it; the crop runs to the first member row | The ribbon is no longer hidden behind the popup and the footer is gone; legend fragments remain at the lower left, as before; accepted |
| `BGC0000708-BGC0000713` | `manual-10-01-feature-popup.webp` | Recapture; recipe corrected | livE, its highlight and label left of the **Qualifiers** popup; drawer toggle hidden | No app-header fragments; accepted |
| `Vnig_TUMSAT-TG-2018` | `manual-08-01-feature-popup.webp` | Recapture; recipe corrected | The clicked dnaA CDS is highlighted beside the 720 px popup; the whole 4.0 Mbp tick label is visible | dnaA was not visible before; accepted |
| `hepatoplasmataceae_collinear` | `manual-07-01-collinear-block-popup.webp` | Recapture; recipe corrected | The whole highlighted block_0024 is to the right of the popup | Fragments of the AP027133.1 record label remain between popup and block because that label ends 4 px left of the block; accepted |
| `hepatoplasmataceae_collinear`, `hepatoplasmataceae_orthogroup` | `post-01-01-feature-popup.webp`, `post-02-01-feature-popup.webp` | Recapture; recipe corrected | The highlighted DnaA at the start of AP027078.1 is to the right of the popup; the popup is below the toolbar | Taller because the popup now shows the **Feature placement** row; no Run info or toolbar fragments; accepted |
| `majanivirus_orthogroup` | `manual-08-01-orthogroup-popup.webp` | Recapture; recipe corrected | The whole og_31 ribbon and its highlight are visible; the other og_31 members show the group outline | The ribbon was cut off at the left edge before; accepted |

The Vibrio feature-popup caption says the popup reports qualifiers and the
sequence, but the image shows the **Details** tab, where those appear only as
tab names. This predates the recapture and is left for a caption review.

## Circular record discovery and loaded-preview inspection (#597 S07)

| Tutorial | Operation media | Decision | Required capture state | Status |
| --- | --- | --- | --- | --- |
| `tobacco-chloroplast` | `manual-06-01-region-annotations.webp` | Keep; recipe corrected | Exact tobacco Session; **Inspect source records** replaces the removed Circular **Load record rotation controls** click before the four region rows are asserted | Corrected recipe captured and passed every declared assertion; the new crop truncates the label cells, so the existing bitmap is kept |
| `vibrio-harveyi-group-collinear` | `manual-02-01-record-row.webp` | Keep; recipe corrected | Exact Vibrio Session; the removed click is dropped because the loaded Linear record rows already exist | Corrected recipe captured the same row and values; only the region label gained its later **Applies on Generate** note, so the existing bitmap is kept |

Tutorial text and captions are unchanged. Public documentation images for the
changed Circular upload and Session steps come from `T-GUI-01` and `T-GUI-09`;
see `results/S07_RESULT.md` in the Issue #597 plan for the remaining stale
Circular Tutorial captures that this change did not regenerate.

## Circular Width/Radius numeric and unit controls (#619 S03)

| Tutorial | Operation media | Decision | Required capture state | Status |
| --- | --- | --- | --- | --- |
| `tobacco-chloroplast` | `manual-07-01-custom-track-slots.webp` | Replace | Exact tobacco Session; features Auto; plastome_regions 20 px / 0.65 ×R; GC content 0.08 ×R / 0.56 ×R; explicit numeric/unit controls | Recaptured at DSF 3 / quality 94; eight value/unit controls passed; equal-width visual comparison accepted |

Keep the three-row stack, annotation binding and all finished-preview media.
Update the existing table and caption to name numeric values and units. A taller viewport keeps the added unit/help controls fully visible without widening
the crop. The operation declares its own Session, app state and eight visible
controls. The existing GUI Tutorial's corresponding track-controls crop uses
its original `T-GUI-05` owner recipe; unrelated images remain unchanged.

## Alignment direction and Reset (#598 S04)

The existing BGC tutorial keeps its public reader route and media. The alignment
instruction now names exclusive Keep/right/left/Custom and both Reset scopes.
`manual-08-01-align-og1.webp` shows the feature-popup entry points, with no
retired Match control. The final preview uses automatic Keep alignment. Both
images remain truthful: **Keep**, with no recapture or new public smoke figure.
The five-source GUI recipe regenerates the finished BGC example and verifies
optional minority-reference reversal, both Reset scopes and Undo; see
`docs/capture/README.md`. Internal acceptance captures remain internal.

## Optional Similarity Group alignment review (#586)

| Tutorial | Operation media | Decision | Required capture state | Status |
| --- | --- | --- | --- | --- |
| `BGC0000708-BGC0000713` | `manual-08-01-align-og1.webp` | Recapture | Exact BGC Session; clicked livE feature in og_1; popup shows both **Align…** and **Review alignment options…** in one crop | Recaptured at DSF 3, quality 94; old/new reviewed at equal display size; current controls and clicked feature verified |

## Two-species Vibrio collinearity example

The Product owner selected a smaller public comparison containing only
*Vibrio parahaemolyticus* RIMD 2210633 and *Vibrio alginolyticus* NBRC 15630.
The stable Gallery ID remains `vibrio-harveyi-group-collinear`, while its owner
records table, generated session/SVG/thumbnail, tutorial copy, and every
data-dependent operation capture are regenerated for four chromosomes in two
rows (`1,1,2,2`). The LOSATP scope remains Adjacent pairs, now covering one
species boundary and four cross-record combinations.

| Tutorial | Operation media | Decision | Required capture state | Status |
| --- | --- | --- | --- | --- |
| `vibrio-harveyi-group-collinear` | all referenced media under its own media directory | Recapture | Exact regenerated two-species session; 4 Linear records; rows `1,1,2,2`; *V. parahaemolyticus* and *V. alginolyticus* identity; one adjacent collinear boundary | Recaptured at DSF 3 and quality 94; strict validation passed; final 3804×591 overview and operation crops visually accepted |

## Linear input workflow cleanup (#568 Phase 1)

The Lambda tutorial is the only Gallery example that shows the renamed Linear
GenBank uploader and the single-record comparison state. Update its capture
metadata and replace only the two affected operation crops; keep the remaining
control and final-preview media.

| Tutorial | Operation media | Decision | Required capture state | Status |
| --- | --- | --- | --- | --- |
| `lambda_basic_linear` | `manual-02-01-genbank-upload.webp` | Recapture | Exact Lambda session; one Linear row; **GenBank / DDBJ File** uploader containing `NC_001416.gb`; no record-reorder or source-card Remove action | Recaptured at DSF 3, quality 94; visually accepted |
| `lambda_basic_linear` | `manual-02-03-no-comparison.webp` | Recapture | Exact Lambda session; **No comparison** command; one **Current: No comparison** status; disabled **Run LOSAT** with its two-input requirement | Recaptured at DSF 3, quality 94; visually accepted |

## Titles, Record Labels, and Legend regrouping (#562)

The combined **Title & Legend** card was replaced by **Titles & Record Labels**
and a separate **Legend · position** card. Recapture only operations that showed
the replaced controls; retain input, comparison, and generated-preview media.

| Tutorial | Operation media | Decision | Required capture state | Status |
| --- | --- | --- | --- | --- |
| `HmmtDNA_ATskew` | `manual-03-03-legend-position-left.webp` | Replace | Separate Legend card; Position menu open with Left selected | Recaptured at DSF 3, quality 94; visually accepted |
| `lambda_basic_linear` | `manual-03-04-legend-position-left.webp` | Replace | Separate Legend card; Position menu open with Left selected | Recaptured at DSF 3, quality 94; visually accepted |
| `BGC0000708-BGC0000713` | `manual-05-03-legend-position-bottom.webp`, `manual-06-01-title-record-text.webp`, `manual-06-02-record-labels.webp` | Replace/add | Bottom Legend menu; focused Plot Title subsection; focused Record Labels rows with Show/Show visibility | Recaptured at DSF 3, quality 94; old/new same-size review accepted the smaller focused crops |
| `Vnig_TUMSAT-TG-2018` | `manual-03-02-legend-position-left.webp`, `manual-05-01-bottom-title.webp` | Replace | Separate Legend card; focused Plot Title subsection with Bottom selected | Recaptured at DSF 3, quality 94; old/new review accepted |
| `majanivirus_orthogroup` | `manual-04-02-legend-position-right.webp` | Replace | Separate Legend card; Position menu open with Right selected | Recaptured at DSF 3, quality 94; visually accepted |
| `vibrio-harveyi-group-collinear` | `manual-07-02-bottom-legend.webp` | Replace | Separate Legend card; Position menu open with Bottom selected | Recaptured at DSF 3, quality 94; visually accepted |

The tobacco tutorial changes table labels only; its referenced operations do
not show the replaced controls. No generated preview or session artifact is
recaptured for this UI-only organization change.

## Comparison commands in History (05A4-02)

Keep the existing comparison-operation media, captions, values, and crops in
the six Linear tutorials. Their capture scripts now await the comparison
History transaction before setting LOSAT modes or reading its state. This
changes execution ordering only; the declared screenshot state is unchanged.
Replay verified all 23 affected operations. Visual comparison also found
pre-existing image drift in additional comparison/collinear controls and the
Vibrio selected-pair count. Those images are retained under this session's
History-only scope; this entry does not accept them as a completed image refresh.

## Annotation TSV download reconciliation

The pre-freeze comparison uses source `e2f867087c9bab1f129484dff8ec1df463392c06`.
The existing chloroplast annotation screenshots omit the current **Download TSV**
button. Keep both tutorial pages and their capture owners. Recapture the
`T-GUI-05` annotation view and
`tobacco-chloroplast/manual-06-01-region-annotations.webp`, preserving their
existing scopes. Add the download button to the capture assertions and retain
the four region rows. Keep other Gallery media unless a separate visible
mismatch is demonstrated.

All 26 manifest-owned GUI scenarios were then executed against this source.
The review found stale Record options, comparison settings, feature editors,
layout controls, and last-successful-result UI in 84 of their 88 images. Those
images were regenerated by their existing owners; four matching images were
retained. Capture repairs use the current **Comparison match style** label and
**Undo Generate diagram** history entry, scroll the relevant input panels into
view, and wait for layout-edit zoom before dragging the legend. Finished views
keep their titles visible and the edited legend clear of gene labels.

The Gallery annotation crop was visually accepted with **Download TSV**, all
four region IDs, and `NC_001879.2` selected. Its owner first opens **Load record
rotation controls** so lazy record discovery has finished before checking the
target selector. The other 111 Gallery WebP files and all tracked Gallery SVGs
were retained. Strict capture metadata validation passed for 113 operations.

The full CLI/Python recipe run also found old comparison metadata in the
`T-CLI-10` and `T-PY-07` SVGs, plus a session 40/schema 5 download and an old
Interactive SVG script in `T-CLI-11`. Those four artifacts were regenerated
through the existing recipe owners and the session assertion now uses 41/schema 7.
The collinear SVGs retain their geometry, text, and colors; the handoff's static
SVG and replay remain byte-identical. No renderer or compatibility path changes
are needed.

## Coordinate-scale visibility (#311 and #315)

| Tutorial | Operation media | Decision | Required capture state | Status |
| --- | --- | --- | --- | --- |
| `lambda_basic_linear` | `manual-03-03-scale-style-ruler.webp` | Recapture | Exact Lambda session; Linear mode; **Show Coordinate Scale** selected; **Scale Style** set to **Ruler (Ticks)**; both controls inside the crop | Recaptured and visually verified (990×768, DSF 3, quality 94) |
| `hepatoplasmataceae_orthogroup` | `manual-05-01-layout-gc-skew-ruler.webp` | Recapture | Exact orthogroup session; Linear mode; **Show Coordinate Scale** selected; **Scale Style** set to **Ruler (Ticks)**; both controls inside the crop | Recaptured and visually verified (990×768, DSF 3, quality 94) |
| `hepatoplasmataceae_collinear` | `manual-05-01-layout-gc-skew-ruler.webp` | Recapture | Exact collinear session; Linear mode; **Show Coordinate Scale** selected; **Scale Style** set to **Ruler (Ticks)**; both controls inside the crop | Recaptured and visually verified (990×768, DSF 3, quality 94) |
| `BGC0000708-BGC0000713` | `manual-05-02-scale-style-ruler.webp` | Recapture | Exact BGC session; Linear mode; **Show Coordinate Scale** selected; **Scale Style** set to **Ruler (Ticks)**; both controls inside the crop | Recaptured and visually verified (990×768, DSF 3, quality 94) |
| `HmmtDNA_ATskew` | `manual-07-01-tick-track-context.webp` | Keep | Exact mitochondrial session; custom stack enabled; final slot is enabled `ticks` with `label_in_tick_out` | Current image remains truthful; deterministic metadata added |
| `vibrio-harveyi-group-collinear` | `manual-06-01-linear-layout.webp` | Keep | Layout card only; the settings table retains the visible ruler configuration | No Axis & Scale controls in the crop |

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
| `H-GUI-03` | `record-layout.png`, `orientation-result.png` | Contextual only: `vibrio-harveyi-group-collinear/manual-08-01-collinear-overview.webp` | Docs use two comparison-free phage rows at 24 px; Gallery uses four comparison records in two rows at 48 px with bottom legend | Recapture docs result; strengthen Gallery app-state/visible-text assertions before its result recapture | Documentation images accepted; contextual Gallery WebP recaptured at 3804×591 |
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
- Use 60% only where 70% clips required content. This applies to `H-GUI-11`,
  `H-GUI-12`, `H-GUI-13`, `H-GUI-14`, and `H-GUI-15`.
- Use 50% only for the dense final tobacco chloroplast figure (`T-GUI-05`) and
  the three-comparison-ring figures (`T-GUI-06` and `H-GUI-06`), where 60%
  clips the title, labels, legend, or outer comparison ring. `T-GUI-01` and
  `T-GUI-09` also use 50% since #597 S07: the Generation status row above the
  Preview leaves too little canvas height at 60%, and bottom labels such as
  `tRNA-Asp` fall behind the preview toolbar.
- The HmmtDNA feature-highlight result uses `Middle`, strand separation off,
  70%, and gene labels for all 13 mitochondrial CDS features, including
  `COX1`.
