# Web Gallery operation screenshot register

Last updated: 2026-10-24

This register records task-specific decisions for Gallery operation media.
Capture metadata remains the executable source of truth in each tutorial JSON.

## Progressive hepatoplasmataceae_collinear tutorial and fresh-upload Session (2026-10-24)

The tutorial now starts from the Web defaults and Generates four times: the
defaults, LOSATP collinear blocks (the only LOSAT run; its settings come before
it), Orientation + identity, and the layout. The Gallery Session was saved from
a fresh upload of the five pinned files (`tutorial-data/`) followed by these
steps, then published with `tools/refresh_gallery_sessions.py`, so the Session,
figure, captures, and a reader's run show the inferred organism name lines
(`<i>Ca.</i> Tyloplasma litorale Fukuoka2020`, OV-1698, OV-1720). Publication
now keeps the File card rows that Save writes (`webFiles.bindings.linearSeqs`),
so a loaded Session holds each name as inferred, not typed, and **Reset
Settings** keeps it (OV-1701). The steps therefore replay from the Session
alone, with no file chosen again. The four replayed previews match the
fresh-upload run element for element.

| Tutorial | Operation media | Decision | Required capture state | Status |
| --- | --- | --- | --- | --- |
| `hepatoplasmataceae_collinear` | `manual-04-01-default-preview.webp`, `manual-05-04-collinear-preview.webp`, `manual-06-02-identity-preview.webp`, `manual-07-03-collinear-overview.webp` | Add (the last replaces `manual-06-01-collinear-overview.webp`) | Replay through steps 4, 5, 6, 7; SVG at 1370 px in a 1900x1000 viewport; name lines visible | Captured at DSF 3; accepted |
| `hepatoplasmataceae_collinear` | `manual-03-01-add-sequence.webp` | Add | Replay through step 3 (Session alone); the File 5 card (NZ_CP006932.gb) and **Add sequence**, clipped together | Captured at DSF 3; accepted |
| `hepatoplasmataceae_collinear` | `manual-05-01-open-pairwise.webp`, `manual-05-02-collinear-reduction.webp`, `manual-05-03-runtime-reproducibility.webp`, `manual-06-01-orientation-identity.webp`, `manual-07-01-layout-gc-skew-ruler.webp`, `manual-07-02-layout-overview.webp`, `manual-08-01-collinear-block-popup.webp` | Recapture; renamed from `manual-03-01`, `04-01`, `03-02`, `04-02`, `05-01`, `05-02`, `07-01` (the step numbers moved) | Session alone, as before | Captured at DSF 3; accepted |
| `hepatoplasmataceae_collinear` | `manual-02-01-upload-row-context.webp` | Recapture | Replay through step 2 (Session alone); File 1's card found by its record uid; empty name fields with the inferred name as placeholder | Captured at DSF 3; accepted |
| `hepatoplasmataceae_collinear` | `manual-01-01-linear-mode.webp`, `post-01-01-feature-popup.webp` | Recapture | Session alone, as before | Captured at DSF 3; accepted |
| `hepatoplasmataceae_collinear` | `manual-04-03-advanced-collinear.webp` | Delete | None: it showed Web defaults (Diagonal drift 0, Merge conflicts 1), now named in the step 7 note | Deleted |

## Progressive HmmtDNA_ATskew, tobacco-chloroplast, and Vnig_TUMSAT-TG-2018 tutorials (2026-10-23)

The three Circular tutorials now start from the Web defaults and Generate after
each change, as the BGC and Vibrio tutorials do. Every manual step carries
`replay`, and each preview is captured by replay through its own step
(Gallery Session, **Reset Settings**, steps 1..N). The tobacco steps choose
their tables with the capture tool's `upload` action (the files under
`gallery/files/`, the reader's downloads). Each tutorial now lists every
setting its Gallery figure needs: Multi-Record Canvas off and Separate Strands
off for HmmtDNA_ATskew (OV-1695); Multi-Record Canvas off, the outline colors
and widths, `rep_origin`, and Circular Label Spacing for tobacco (OV-1699).
Settings that equal the Web defaults (Track Preset Tuckin, Legend position
Left on Circular, Separate Strands on for Vnig) are not taught. Sessions,
source SVGs, thumbnails, and `examples.json` are unchanged.

| Tutorial | Operation media | Decision | Required capture state | Status |
| --- | --- | --- | --- | --- |
| `Vnig_TUMSAT-TG-2018` | `manual-03-01-default-preview.webp`, `manual-04-02-plasmid-row-preview.webp`, `manual-05-02-orchid-preview.webp`, `manual-06-02-multirecord-preview.webp` | Add (the last replaces `manual-06-01-multirecord-preview.webp`) | Replay through steps 3, 4, 5, 6; Record Order rows and palette asserted per stage | Captured at DSF 3; accepted |
| `Vnig_TUMSAT-TG-2018` | `manual-04-01-record-order.webp` | Replace `manual-04-00-multirecord-toggle.webp` | Replay through step 4; the row select of record #3 (found by its `Row for #3 (` name) highlighted with value 2 | Captured at DSF 3; accepted |
| `Vnig_TUMSAT-TG-2018` | `manual-05-01-palette-orchid.webp`, `manual-06-01-bottom-title.webp`, `manual-07-01-feature-popup.webp` | Recapture; renamed from `manual-03-01-…`, `manual-05-01-…`, `manual-08-01-…` | Session alone, as before | Captured at DSF 3; accepted |
| `Vnig_TUMSAT-TG-2018` | `manual-03-00-circular-layout-settings.webp`, `manual-03-02-legend-position-left.webp`, `manual-03-03-selected-features.webp` | Delete | None: they showed Web defaults | Deleted |
| `HmmtDNA_ATskew` | `manual-03-01-default-preview.webp`, `manual-04-03-labels-preview.webp`, `manual-05-01-restyle-preview.webp`, `manual-06-03-layout-preview.webp`, `manual-07-02-middle-stack-preview.webp`, `manual-09-02-atskew-preview.webp` | Add (the last replaces `manual-09-01-atskew-preview.webp`) | Replay through steps 3, 4, 5, 6, 7, 9; the final asserts the five-slot stack and the a_skew_2 parameters | Captured at DSF 3; accepted |
| `HmmtDNA_ATskew` | `manual-06-01-layout-options.webp` | Add | Session; Separate Strands and Multi-Record Canvas cleared (`visibleControls`) | Captured at DSF 3; accepted |
| `HmmtDNA_ATskew` | `manual-08-01-at-skew-slot-context.webp` | Recapture; renamed from `manual-06-01-…`; recipe corrected | Replay through step 8, so ticks is still second and a_skew_2 last (slot 4); disclosure open from the replay | Captured at DSF 3; accepted |
| `HmmtDNA_ATskew` | `manual-04-01-label-mode-out.webp`, `manual-04-02-cds-label-priority.webp`, `manual-06-02-window-step.webp`, `manual-07-01-middle-template.webp`, `manual-09-01-tick-track-context.webp` | Recapture or rename (renumbered from `manual-03-02`, `08-01`, `04-01`, `05-01`, `07-01`) | Session alone, as before; `window-step` has no recipe and keeps its bytes | Accepted |
| `HmmtDNA_ATskew` | `manual-03-01-track-preset-middle.webp`, `manual-03-03-legend-position-left.webp` | Delete | None: the Middle template sets the track preset, and Left is the Circular default | Deleted |
| `tobacco-chloroplast` | `manual-03-01-default-preview.webp`, `manual-04-04-labels-preview.webp`, `manual-05-02-colors-preview.webp`, `manual-06-03-regions-preview.webp`, `manual-07-02-tint-preview.webp`, `manual-08-02-chloroplast-preview.webp` | Add (the last replaces `manual-08-01-chloroplast-preview.webp`) | Replay through steps 3..8 at 1400x1400 with the footer hidden; the final also waits for the eight annotation ids | Captured at DSF 3; accepted |
| `tobacco-chloroplast` | `manual-04-03-qualifier-priority.webp`, `manual-05-01-specific-color-table.webp` | Recapture; recipe replaced | Replay through step 4 or 5: the file is chosen with `upload`, so the crop shows the reader's file name without renaming the restored File | Captured at DSF 3; accepted |
| `tobacco-chloroplast` | `manual-06-02-region-ring-slots.webp` | Add | Replay through step 6: features, plastome_regions, gc_content with their Width/Radius values | Captured at DSF 3; accepted |
| `tobacco-chloroplast` | `manual-08-01-layout-options.webp` | Replace `manual-04-01-layout-options.webp` | Replay through step 8; Separate Strands on, Multi-Record Canvas off | Captured at DSF 3; accepted |
| `tobacco-chloroplast` | `manual-02-02-basic-settings.webp`, `manual-04-01-label-mode-both.webp`, `manual-04-02-label-placement-radial.webp` | Recapture; renamed from `manual-03-01`, `04-02`, `04-03` | Session alone, as before | Captured at DSF 3; accepted |

## Edit-surface redesign: feature popup, search bar, Legend tab (2026-10-12)

Operations that show the redesigned feature popup, the Feature search bar, or the Editor drawer's Legend tab.
Recaptured at the tool default scale; an image is kept only when it visibly differs from the committed one.

| Tutorial | Operation media | Decision | Required capture state | Status |
| --- | --- | --- | --- | --- |
| `vibrio-harveyi-group-collinear` | `manual-06-01-search-dnaa.webp`, `manual-07-01-search-rctb.webp` | Recapture; recipe corrected | Search bar with **Regular expression** and the status `1 of 2 matches`; the replay waits on that status; the open feature popup is closed first | Accepted; no popup fragment |
| `vibrio-harveyi-group-collinear` | `manual-06-02-rotate-start.webp` | Recapture; text corrected | **Edit** tab: header gene, product, record and location; **Similarity group** above **Layout**; **Apply when I Generate** and **Apply now** | Accepted |
| `BGC0000708-BGC0000713` | `manual-07-01-align-og6.webp` | Recapture; recipe corrected | neoU popup with gene and product in the header, **Similarity group**, **Align…**, **Review alignment options…**; the recipe sets the popup's own width to 500 px (the size state it wrote before is not on `window.__GBDRAW_APP__`, so the popup stayed 720 px wide and left the window), crops from the whole popup header, and stops when the popup does not fit the window | Accepted at 1926x858: close button, tabs, and **Align…** inside the crop |
| `BGC0000708-BGC0000713` | `manual-12-01-feature-popup.webp` | Recapture; recipe corrected | **Qualifiers** tab as a name and value list (no row borders): the recipe marks each name and the value beside it and scrolls the name to the top | Accepted |
| `BGC0000708-BGC0000713` | `manual-11-01-orthogroup-popup.webp` | Keep | Similarity-group popup unchanged apart from sub-percent noise | Committed bytes kept |
| `HmmtDNA_ATskew` | `manual-10-01-feature-popup.webp` | Recapture; recipe corrected | RNR2 **Details** popup left of the highlighted l-rRNA; the recipe pans the feature to a fixed place in the preview before the click, where it left the place to the scroll that the click does on a feature outside the visible preview (three different outcomes seen, one capture run timed out) | Accepted at 3153x903; two runs byte-identical |
| `HmmtDNA_ATskew` | `post-02-01-legend-editor.webp` | Recapture; recipe added | Exact Session; **Editor** drawer open on **Legend**; list scrolled to the AT skew rows; crop from the drawer header to the last row | Accepted; the AT skew swatches show the GC skew colors (Product bug, see the handoff) |

The capture tool's `resetSettings` action confirms the in-app **Reset all settings?** dialog; the tool no longer expects a browser dialog.

## BGC progressive tutorial (Gallery revision 2026-10-11)

The BGC Gallery Session is rebuilt from the five MIBiG files in
`gbdraw/web/tutorial-data/aminoglycoside-bgc-five` through the Web app (D-08):
BGC0000713 is the original file with record-level **Reverse complement** on, so
**Reset Settings** returns to the reader's starting state. The tutorial starts
from the Web defaults and shows the preview after each Generate (D-01): every
manual step carries `replay`, and each preview capture uses `replayThrough` of
its own step. The route that builds the Session runs the same `replay` actions
(task folder `scripts/b_route.py`). Before each upload it gives the empty slot
the key `record-<n>`, the key CLI-written Sessions use, instead of the random
`linear-seq-<uuid>` key; nothing refers to the key yet, so the LOSAT cache and
protein identities are built from it. Renaming keys after the fact breaks
those hashes. Earlier BGC media were renamed to their new step numbers or
deleted.

| Tutorial | Operation media | Decision | Required capture state | Status |
| --- | --- | --- | --- | --- |
| `BGC0000708-BGC0000713` | `manual-02-01-file1-defaults.webp`, `manual-03-01-file5-defaults.webp` | Add (replace `manual-02-01-record-text-row.webp`) | Replay through step 2 or 3; File 1 (or File 5 with **Add sequence**) from its header to **File defaults** | Captured at DSF 3; accepted |
| `BGC0000708-BGC0000713` | `manual-04-01-preview-defaults.webp`, `manual-05-03-preview-losatp.webp`, `manual-06-02-preview-reversed.webp`, `manual-07-04-preview-aligned.webp`, `manual-08-02-preview-colors.webp`, `manual-09-02-preview-labels.webp`, `manual-10-04-preview-final.webp` | Add (replace `manual-08-04-aligned-preview.webp`) | Replay through the step; the whole SVG at 1100 px, preview controls hidden; state asserted per stage | Captured at DSF 3; accepted |
| `BGC0000708-BGC0000713` | `manual-05-01-comparison-settings.webp` | Replace `manual-03-01-open-pairwise.webp` | Final Session; **Settings** open; LOSATP, Similarity groups, Minimum identity 30 | Captured at DSF 3; accepted |
| `BGC0000708-BGC0000713` | `manual-07-01-align-og6.webp`, `manual-07-02-anchor-candidates.webp`, `manual-07-03-alignment-direction.webp` | Recapture (renamed from `manual-08-0x`) | Replay through step 6 (reversed, not yet aligned, default colors) | Captured at DSF 3; accepted. 2026-10-12: `-02` and `-03` recaptured for the docked Review alignment panel (After Apply strip, no candidate numbers); `-01` (popup opener) is byte-identical |
| `BGC0000708-BGC0000713` | `manual-08-01-specific-rules-all.webp`, `manual-09-01-show-labels-first.webp`, `manual-10-01-arrow-shaft-width.webp`, `manual-10-02-scale-style-ruler.webp`, `manual-10-03-plot-title.webp`, `manual-11-01-orthogroup-popup.webp`, `manual-12-01-feature-popup.webp` | Recapture (renumbered) | Final Session. The og_18 popup now shows the computed group name: the old Session carried a rename of og_18 to livA that no step taught | Captured at DSF 3; accepted |
| `BGC0000708-BGC0000713` | `manual-03-02-runtime-reproducibility.webp`, `manual-03-03-first-comparison-boundary.webp`, `manual-03-04-first-raw-result.webp`, `manual-04-03-track-layout-middle.webp`, `manual-06-02-record-labels.webp` | Delete | None: the progressive tutorial shows only controls the reader changes (Track Layout Middle and the record-label visibility are defaults) | Deleted |

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

### Progressive Vibrio tutorial (2026-10-11)

The tutorial now teaches in stages: add File 1 and fill its File defaults, add
File 2 the same way, Generate with the Web defaults, run LOSATP Collinear blocks
once, rotate the four chromosomes to *dnaA* and *rctB*, Generate again from the
cached LOSAT results, then finish the figure. Each Generate step shows its own
preview, captured by replay (`replayThrough` = that step): the tool loads the
Gallery Session, clicks **Reset Settings**, and performs the steps' `replay`
actions. Every preview recipe fails if a record label shows an accession,
replicon, or length (D-07). The final replayed figure matches the Gallery SVG
in layout, rotations, and all 116 comparison ribbons; only the ruler tick and
label color differs (black instead of dimgray, OV-440).

| Tutorial | Operation media | Decision | Required capture state | Status |
| --- | --- | --- | --- | --- |
| `vibrio-harveyi-group-collinear` | `manual-02-01-file1-defaults.webp`, `manual-03-01-file2-defaults.webp` | Add (replace `manual-02-01-record-row.webp`) | Replay through step 2 or 3; File card with the GBFF name, **File defaults** filled, **Number of records: 2**; File 2 crop includes **Add sequence** | Captured at DSF 3; accepted |
| `vibrio-harveyi-group-collinear` | `manual-04-01-default-preview.webp` | Add | Replay through step 4; Web defaults, rows `1,1,2,2`, no comparison, species and strain only | Captured at DSF 3; accepted |
| `vibrio-harveyi-group-collinear` | `manual-05-05-unrotated-blocks.webp` | Add | Replay through step 5; LOSATP Collinear blocks on unrotated records; diagonal ribbons | Captured at DSF 3; accepted |
| `vibrio-harveyi-group-collinear` | `manual-08-01-rotated-blocks.webp` | Add | Replay through step 8; rotated records, search cleared, vertical blocks | Captured at DSF 3; accepted |
| `vibrio-harveyi-group-collinear` | `manual-08-01-collinear-overview.webp` → `manual-10-01-collinear-overview.webp` | Replace | Replay through step 10; the finished figure, legend gap check | Captured at DSF 3, 3804×591; accepted |
| `vibrio-harveyi-group-collinear` | `manual-03-10-search-dnaa.webp`, `manual-03-11-rotate-start.webp`, `manual-03-13-search-rctb.webp` → `manual-06-01-…`, `manual-06-02-…`, `manual-07-01-…` | Recapture; renamed | Replay through step 5 (dnaA) or 6 (RctB), so the search and rotation act on unrotated records | Captured at DSF 3; accepted |
| `vibrio-harveyi-group-collinear` | `manual-04-00`, `manual-04-01`, `manual-04-02`, `manual-04-04` → `manual-05-01` … `manual-05-04`; `manual-05-01-rectangle-features.webp`, `manual-06-01-linear-layout.webp` → `manual-09-01-…`, `manual-09-02-…` | Recapture; renamed | Session alone; the values equal the final Session | Captured at DSF 3; accepted |
| `vibrio-harveyi-group-collinear` | `manual-03-01-record-layout.webp`, `manual-04-03-adjacent-pairs.webp`, `manual-07-02-bottom-legend.webp` | Delete | None: Arrange in rows, Adjacent pairs, and Legend Bottom are Web defaults | Deleted |

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
| `H-GUI-03` | `record-layout.png`, `orientation-result.png` | Contextual only: `vibrio-harveyi-group-collinear/manual-10-01-collinear-overview.webp` | Docs use two comparison-free phage rows at 24 px; Gallery uses four comparison records in two rows at 48 px with bottom legend | Recapture docs result; strengthen Gallery app-state/visible-text assertions before its result recapture | Documentation images accepted; contextual Gallery WebP recaptured at 3804×591 |
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

## Compact sidebar input cards (Gallery revision 2026-10-12)

The Web sidebar's File cards, Record options, Colors, and slot rows were
compacted. The File card's defaults block now reads **Names for all records in
this file**, and the tutorial text, `visibleText`, alt text, and captions follow
that label. Media were recaptured only where the surface visibly changed;
captures that differ from the committed bytes by noise only keep the committed
bytes.

| Tutorial | Operation media | Decision | Required capture state | Status |
| --- | --- | --- | --- | --- |
| `BGC0000708-BGC0000713`, `vibrio-harveyi-group-collinear` | `manual-02-01-file1-defaults.webp`, `manual-03-01-file5-defaults.webp` (Vibrio: `manual-03-01-file2-defaults.webp`) | Recapture | Replay through the step; File card with **Names for all records in this file** and the rendered taxon-name preview line | Recaptured at DSF 3; compared by eye |
| `BGC0000708-BGC0000713` | `manual-06-01-reverse-bgc0000713.webp` | Recapture | File 5 card from its header through the compact **Record options** (record header, **Using file default** names, **Reverse complement**) | Recaptured at DSF 3; compared by eye |
| `BGC0000708-BGC0000713` | `manual-08-01-specific-rules-all.webp` | Recapture; scroll corrected | **SPECIFIC RULES** heading scrolled to `block: 'start'` with a 120 px scroll margin so the app header does not cover the **TSV** / **Clear All** row; the four rule rows are opened | No app-header overlap; visually accepted |
| `BGC0000708-BGC0000713` | `manual-10-03-plot-title.webp` | Recapture | **Titles & Record Labels** card; Plot Title with its markup preview line | Recaptured at DSF 3; compared by eye |
| `BGC0000708-BGC0000713` | `manual-05-01-comparison-settings.webp`, `manual-09-01-show-labels-first.webp`, `manual-10-01-arrow-shaft-width.webp`, `manual-10-02-scale-style-ruler.webp`; Vibrio `manual-09-01-rectangle-features.webp`, `manual-09-02-linear-layout.webp` | Keep | Surfaces without a visible change | Recapture differed by noise only |
| `BGC0000708-BGC0000713` | `manual-06-01-reverse-bgc0000713.webp` (record wait) | Recipe corrected | The wait for the record-5 select proves the records were read: the select `Record selector for sequence 5` is enabled and lists `BGC0000713` | The old wait on the absence of "Loading records" proved nothing |
| `HmmtDNA_ATskew` | `manual-06-01-at-skew-slot-context.webp` | Recapture; recipe corrected | Exact mitochondrial session; custom stack enabled; the capture opens the slot's `[data-track-slot-more='a_skew_2']` disclosure; `assertAppState` (slot 3 id, renderer, nt, legend label, both colors) and six `visibleControls` (id, renderer, dinucleotide AT, legend label AT skew, both colors) | Data-dependent; recaptured 810x2319; disclosure open, both swatches inside the highlight |

## Short/long size boundary: BGC previews and step 10 (Gallery revision 2026-10-12)

Linear Auto sizes use the short styles for a longest record under 100 kb (it
was 50 kb), so the previews that replay the reader's steps from the Web
defaults now draw the five clusters (longest 50,466 bp) with the short sizes.
The published figure does not change: its Session fixes Head Length Ratio
0.64, and step 10 types that value. The step 10 replay gives the published
SVG byte for byte, so the finished-figure and popup captures keep their bytes.

| Tutorial | Operation media | Decision | Required capture state | Status |
| --- | --- | --- | --- | --- |
| `BGC0000708-BGC0000713` | `manual-04-01-preview-defaults.webp`, `manual-05-03-preview-losatp.webp`, `manual-06-02-preview-reversed.webp`, `manual-07-04-preview-aligned.webp`, `manual-08-02-preview-colors.webp`, `manual-09-02-preview-labels.webp` | Recapture | Replay through the step, recipes unchanged; short Auto sizes (80 px features on two strand lanes, 24 px definition lines) | Recaptured at DSF 3; each compared with the previous image side by side |
| `BGC0000708-BGC0000713` | `manual-07-01-align-og6.webp` | Recapture | Replay through step 6, recipe unchanged; the highlighted neoU arrow is drawn at the short feature height beside the same popup | Recaptured at DSF 3, same 2502x792 crop |
| `BGC0000708-BGC0000713` | `manual-10-01-arrow-shaft-width.webp` | Recapture; contract extended | Final Session; **Arrow Geometry** box with Head Length Ratio 0.64 and Shaft Width Ratio 0.6, both highlighted; `assertAppState` and `visibleControls` name both ratios | Recaptured at DSF 3, same 1422x315 crop |
| `BGC0000708-BGC0000713` | `manual-09-02-preview-labels.webp` caption; step 10 text, replay, and both tables | Text | The caption no longer says that the shafts hide the heads of short genes; step 10 has the Head Length Ratio row and says what each group of rows does against the Auto values the reader sees | Tables checked against the Auto notes read after the step 9 replay |
| `BGC0000708-BGC0000713` | `manual-10-04-preview-final.webp`, `manual-11-01-orthogroup-popup.webp`, `manual-12-01-feature-popup.webp` | Keep | The pinned Session and the step 10 replay draw the published figure | Not recaptured |

## Record name fields on one-record File cards (Gallery revision 2026-10-12)

A Linear File with one record and no File-level name shows that record's own
**Organism / strain (optional)** and **Subtitle / title (optional)** on its
card, in place of **Names for all records in this file**; **Record options** no
longer repeat them. The BGC tutorial therefore types the names into the record
fields, and its Gallery Session was saved again through the Web app with the
five names on the records instead of the Files. The published figure is
unchanged: the source SVG, the rendered example, the thumbnail and
`examples.json` keep their bytes. This section replaces the capture state of
the same media in the rows of the earlier sections; those rows stay as history.

| Tutorial | Operation media | Decision | Required capture state | Status |
| --- | --- | --- | --- | --- |
| `BGC0000708-BGC0000713` | `manual-02-01-file1-defaults.webp`, `manual-03-01-file5-defaults.webp` | Recapture; recipe changed | Replay through the step (step 2 `replayThrough` 1 -> 2); clip `:scope > [data-linear-record-names]` of the File card; File card from its header to the record's **Organism / strain (optional)** and **Subtitle / title (optional)** with the markup preview line; no **Names for all records in this file**; File 5 with **Add sequence** | Recaptured at DSF 3, 1164x768 and 1203x1233; compared by eye |
| `BGC0000708-BGC0000713` | `manual-06-01-reverse-bgc0000713.webp` | Recapture | File 5 card with the name fields on the card; open **Record options** hold the record header, Region and **Reverse complement** on; no **Using file default** | Recaptured at DSF 3, 1704x1839; compared by eye |
| `hepatoplasmataceae_collinear`, `hepatoplasmataceae_orthogroup` | `manual-02-01-upload-row-context.webp` (the two files are byte-identical) | Recapture; alt and caption changed | Exact Session; AP027078.gb File card with empty name fields on the card; open **Record options** with the region controls | Recaptured, 1422x1689; compared by eye |
| `majanivirus_orthogroup` | `manual-02-01-upload-row-label-context.webp` | Recapture; recipe changed | The capture finds the field by `[data-linear-record-names="record-1"]`; MjeNMV.gb File card with its Organism / strain value on the card; **Record options** closed | Recaptured, 1782x1002; compared by eye |
| `lambda_basic_linear` | `manual-02-01-genbank-upload.webp` | Keep | The clip is the card header and the uploader; the names area is outside it | Recaptured and reverted |
| `BGC0000708-BGC0000713` | `manual-04-01-preview-defaults.webp`, `manual-10-04-preview-final.webp` | Keep | The replay with the names in the record fields draws the same previews (step 4 SVG byte-identical, step 10 within raster noise) | Recaptured as proof and reverted |

## File card fold line and long upload names (Gallery revision 2026-10-13)

Each Linear File card has a header line with a caret that folds the card, and
the upload zone cuts a long file name at its end with an ellipsis. Two upload
crops were held back from the first recapture of the File cards and are
recorded here; the rows of the earlier sections for the same media stay as
history.

The capture tool now lets running CSS transitions end before it measures and
takes a crop (`wait_for_css_transitions`). A caret turns over 200 ms after its
card or section opens, and the Lambda upload crop, taken right after the tool
opened the loaded Session's folded card, showed the caret mid-turn.

| Tutorial | Operation media | Decision | Required capture state | Status |
| --- | --- | --- | --- | --- |
| `lambda_basic_linear` | `manual-02-01-genbank-upload.webp` | Recapture; recipe unchanged | Exact Lambda Session; open File 1 card from its header line, with the caret at rest, to the **GenBank / DDBJ File** uploader containing `NC_001416.gb` | Recaptured at DSF 3, 864x336 (was 930x324; the same two elements are 22 px narrower and 4 px taller in the current pane); compared side by side |
| `Vnig_TUMSAT-TG-2018` | `manual-02-01-gbff-upload.webp` | Recapture; recipe added | Exact Vnig Session; `dataDependent`; clip from the **Input Genomes** card header to the **GenBank/DDBJ File** uploader; `assertAppState` and `visibleText` name `GCF_015097735.1_ASM1509773v1_genomic.gbff`; the zone shows the name cut at its end with an ellipsis. The rest of the card is outside the clip: the restored Session shows the record hints of the finished figure (Multi-Record Canvas on), not the state at the upload step | Recaptured at DSF 3, 942x483 (was a manual 656x542 crop); compared side by side |

## The qualifier step before the colors (Gallery revision 2026-10-26)

Issue #1385. In the `BGC0000708-BGC0000713` tutorial the qualifier inspection
(old step 12) is now step 8, before the colors; the old steps 8-11 are steps
9-12. The colors, the palette and the color tables stay as they were. The
Gallery Session was saved again from a fresh upload, so it keeps the File card
rows (`webFiles.bindings`) and its record keys are the uids `linear-seq-…`, not
`record-N`; the recipes of steps 6 and 7 find the File 5 card and the
alignment review cards by those uids. This section replaces the capture state
of the same media in the rows of the earlier sections; those rows stay as
history.

| Tutorial | Operation media | Decision | Required capture state | Status |
| --- | --- | --- | --- | --- |
| `BGC0000708-BGC0000713` | `manual-08-01-feature-popup.webp` (was `manual-12-01-feature-popup.webp`) | Recapture; recipe changed | Replay through step 7 (every CDS still blue); `livE` typed in **Search features** and Enter; the livE CDS of BGC0000708 clicked; popup on **Qualifiers** with `gene` `livE` and `gene_kind` `biosynthetic` | Recaptured at DSF 3, 2511x945; compared by eye |
| `BGC0000708-BGC0000713` | `manual-09-02-preview-colors.webp`, `manual-10-02-preview-labels.webp`, `manual-11-04-preview-final.webp` (were `manual-08-02`, `manual-09-02`, `manual-10-04`) | Recapture | The same previews one step later; the current app draws the legend 3 px apart from the committed crops | Recaptured at DSF 3, 3372x1332, 3372x1521 and 3372x1128; compared by eye |
| `BGC0000708-BGC0000713` | `manual-09-01-specific-rules-all.webp`, `manual-10-01-show-labels-first.webp`, `manual-11-01-arrow-shaft-width.webp`, `manual-11-02-scale-style-ruler.webp`, `manual-11-03-plot-title.webp`, `manual-12-01-orthogroup-popup.webp` (were `manual-08-01`, `manual-09-01`, `manual-10-01`, `manual-10-02`, `manual-10-03`, `manual-11-01`) | Rename; keep bytes | Recaptures pixel-identical or within raster noise | Renamed with the committed bytes |
| `BGC0000708-BGC0000713` | `manual-04-01-preview-defaults.webp`, `manual-05-03-preview-losatp.webp`, `manual-06-02-preview-reversed.webp`, `manual-07-04-preview-aligned.webp` | Recapture | Previews of steps 4-7; the legend 3 px apart, as above | Recaptured at DSF 3; compared by eye |
| `BGC0000708-BGC0000713` | `manual-05-01-comparison-settings.webp`, `manual-07-01-align-og6.webp` | Recapture | Current labels: **Max hits per protein for groups**; current tab widths in the alignment panel | Recaptured at DSF 3, 1422x3204 and 1926x858; compared by eye |
| `BGC0000708-BGC0000713` | `manual-07-02-anchor-candidates.webp` and the other media of steps 1-7 | Keep | Recaptures within raster noise (`07-02` maximum channel difference 8) | Recaptured and reverted |
