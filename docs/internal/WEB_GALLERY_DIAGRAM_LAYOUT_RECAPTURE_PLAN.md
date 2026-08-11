# Web and Gallery diagram-layout recapture plan

Status: Completed; owner outputs and linked Gallery results verified
Prepared: 2026-08-04
Implementation owner: WP7 diagram-layout artifact refresh

## Purpose and scope

This is the pre-capture contract for the layout-sensitive GUI documentation
identified by `IMPLEMENTATION_PLAN_DIAGRAM_LAYOUT_OVERHAUL.md`: `T-GUI-01`,
`T-GUI-02`, `T-GUI-05`, `H-GUI-02`, `H-GUI-03`, `H-GUI-09`, `H-GUI-10`, and
`H-GUI-11`. It also identifies the Gallery result crops that give exact or
contextual coverage for those journeys.

The plan does not authorize hand-editing a screenshot, Gallery session, source
SVG, example SVG, thumbnail, or `examples.json`. Documentation PNGs belong to
`docs/capture/run_all.py`. Gallery WebP files belong to tutorial capture
metadata and `tools/capture_gallery_tutorial_screenshots.py`. Gallery sessions
and rendered assets must be refreshed by `tools/refresh_gallery_sessions.py`
before any result crop is replaced.

The Circular record-local GC/skew-radius defect is resolved. Regression and
visual checks verified the corrected renderer before capture. Screenshots must
use that renderer and may not hide or normalize renderer defects.

The final owner run reviewed every in-scope PNG at 1440×900, accepted the
corrected documentation captures, refreshed all 11 Gallery sessions, and
recaptured the six linked Gallery result WebPs from their exact sessions. The
strict Gallery audit now checks 117 operations and media entries, and the
desktop/mobile tutorial and session-regeneration browser suites pass.

## Capture classes

| Class | Meaning | Required framing |
| --- | --- | --- |
| `D-OP` | Documentation operation/state image | Pinned 1440×900 full viewport. Keep the operated control and the neighboring preview/context produced by the real flow. Do not crop. |
| `D-RESULT` | Documentation generated-result image | Pinned 1440×900 full viewport after the flow's real zoom/pan helper. The complete title, primary diagram, labels, and legend must remain readable. |
| `D-EXPORT` | Documentation export action | Pinned 1440×900 full viewport with the real export control visible or hovered. |
| `G-RESULT` | Gallery generated-result media | Exact saved session; tight selector crop around the real SVG; WebP at DSF 3 and quality 94; no app-shell padding except the declared crop padding. |
| `G-OP` | Gallery input/control media | Keep the current compact operation crop unless the operated control or its truthful state changed. A renderer-only geometry change is not a reason to widen or replace it. |

All 25 in-scope documentation PNGs are currently 1440×900 RGB images. The
documentation owner always rewrites every manifest-owned image for a selected
scenario; `Keep expected` below means the regenerated pixels are expected to
remain unchanged and must be investigated if they do not.

## Pinned environments

### Documentation PNGs

- one fresh browser context per scenario; no saved session or browser-storage
  seeding;
- loopback-only app server, service workers blocked, and no external requests;
- Python Playwright 1.61.0 and Chromium 149.0.7827.55 (revision v1228);
- viewport 1440×900, device scale factor 1, RGB PNG;
- `en-US`, UTC, light color scheme, reduced motion, and vendored fonts;
- every screenshot is a full viewport captured by
  `flows.web_capture.capture_screenshot`.

### Gallery WebP

- load the exact `capture.session` listed below after the Gallery session
  refresh; never synthesize a different example's data;
- default light browser UI, device scale factor 3, WebP quality 94;
- use the operation's declared viewport, selector/clip selectors, crop padding,
  state assertions, visible controls, and visible text;
- hide floating preview controls only with capture-time actions already allowed
  by the screenshot-maintenance contract;
- preserve the old bitmap under `/tmp` and compare old/new at identical display
  size before accepting the replacement.

## Documentation capture inventory

The `Caption intent` text is the manifest `reason`. The `Alt text` is the exact
manifest target. The consuming Markdown now uses the same value. Eight
references had drifted from that target:

| Scenario / output | Previous Markdown alt text | Target manifest alt text |
| --- | --- | --- |
| T-GUI-05 / `01-input-ready.png` | Complete tobacco plastome selected as the Circular input | Tobacco plastome GenBank file ready in Circular mode |
| T-GUI-05 / `02-first-diagram.png` | First complete tobacco plastome before Gallery presentation settings | First circular tobacco plastome diagram |
| T-GUI-05 / `03-annotation-table.png` | The four Gallery plastome regions in one annotation set | Annotation table containing LSC, SSC, IRa, and IRb regions |
| T-GUI-05 / `04-track-settings.png` | Gallery circular stack with one inner region lane and GC content | Circular custom-track controls for features, one plastome-region lane, and GC content |
| T-GUI-05 / `05-finished-diagram.png` | Gallery-quality tobacco plastome map | Gallery-style tobacco plastome with functional colors, radial labels, structural regions, GC content, and upper-left legend |
| H-GUI-03 / `record-layout.png` | Linear record rows with complete coordinate ranges | Linear record rows with selected coordinate regions |
| H-GUI-03 / `orientation-result.png` | Linear diagram with complete DE3 reverse complemented | Linear diagram with one reversed record and a coordinate ruler |
| H-GUI-11 / `style-settings.png` | CDS gene-label whitelist and qualifier-priority controls | Feature color, label, title, and legend controls |

### T-GUI-01 — first Circular diagram

Session/state source: fresh context; upload
`gbdraw/web/tutorial-data/human-mitochondrion/HmmtDNA.gbk` and
`gbdraw/web/tutorial-data/shared/cds_gene_qualifier_priority.tsv`. Final state:
Circular single record `NC_012920.1`; output prefix `human_mitochondrion`;
species `<i>Homo sapiens</i>`; Middle; separate strands; GC content/skew shown;
labels Out; right legend.

| Output | Class / decision | Selector or action state | Caption intent | Alt text |
| --- | --- | --- | --- | --- |
| `docs/images/t-gui-01/01-input-ready.png` | `D-OP`; Keep expected | Select `Circular`, `GenBank`, then upload through `GenBank/DDBJ File`; file-selection group contains `HmmtDNA.gbk`. | Confirm the selected tutorial file before generation. | Circular GenBank input showing HmmtDNA.gbk selected |
| `docs/images/t-gui-01/02-first-diagram.png` | `D-RESULT`; Recapture | Click `Generate Diagram`; assert the first `NC_012920.1` SVG before publication settings. | Show the first successful diagram in Step 2. | First circular human mitochondrial genome diagram |
| `docs/images/t-gui-01/03-publication-label.png` | `D-RESULT`; Recapture | Fill `Output Prefix` and `Species`, regenerate, and assert rendered `Homo sapiens` markup. | Verify the species label in the preview. | Circular preview labeled Homo sapiens |
| `docs/images/t-gui-01/04-layout-settings.png` | `D-OP`; Recapture | Select `Track Preset=middle`, check `Separate Strands`, uncheck both Hide GC controls, set `Legend Position=right`, `Label Mode=out`, upload the priority TSV, and scroll `Track Preset` into view. | Record the small set of layout choices used by the final figure. | Circular layout controls set to Middle, Labels Out, and Legend Right |
| `docs/images/t-gui-01/04-finished-diagram.png` | `D-RESULT`; Recapture | Regenerate; assert final semantics; use four real `Zoom out` clicks to 60%; pan the result left so the right legend remains in frame. | Provide the finished result and tutorial hero image. | Finished circular human mitochondrial genome diagram with external labels and a right legend |
| `docs/images/t-gui-01/05-export-svg.png` | `D-EXPORT`; Recapture | Hover the real `SVG` button after the finished preview; the subsequent click must download `human_mitochondrion.svg`. | Locate the standard SVG export action. | SVG download button below the finished result preview |

Gallery relation: exact-parity coverage is
`HmmtDNA_basic_circular/manual-04-01-final-preview.webp`.

### T-GUI-02 — first Linear diagram

Session/state source: fresh context; upload
`gbdraw/web/tutorial-data/lambda/NC_001416.gb` and the shared CDS gene-priority
TSV. Final state: Linear one-record layout; `No comparison`; prefix
`lambda_linear`; Middle; separate strands; all labels; ruler; left legend.

| Output | Class / decision | Selector or action state | Caption intent | Alt text |
| --- | --- | --- | --- | --- |
| `docs/images/t-gui-02/01-input-ready.png` | `D-OP`; Keep expected | Select `Linear`, `GenBank`, and `No comparison`; upload with `GenBank File`; selection contains `NC_001416.gb`. | Confirm the Lambda input and Linear mode. | Linear GenBank input showing NC_001416.gb selected |
| `docs/images/t-gui-02/02-first-diagram.png` | `D-RESULT`; Recapture | Generate the first diagram and call `fit_complete_linear_preview`. | Show the first Linear result in Step 2. | First linear Lambda genome diagram |
| `docs/images/t-gui-02/03-layout-settings.png` | `D-OP`; Recapture | Set Middle, separate strands, all labels, priority TSV, left legend, `Show Coordinate Scale (Linear)=on`, and `Linear scale style=ruler`; keep `Axis & Scale` open and scroll the scale selector into view. | Record the label and ruler settings. | Linear layout controls configured for labels and a ruler |
| `docs/images/t-gui-02/04-finished-diagram.png` | `D-RESULT`; Recapture | Regenerate, assert the labeled/ruler semantics, and fit the complete Linear preview; SVG download must be `lambda_linear.svg`. | Verify the labeled final map before export. | Finished linear Lambda genome diagram with labels and ruler |

Gallery relation: exact-parity coverage is
`lambda_basic_linear/manual-05-01-final-preview.webp`.

### T-GUI-05 — Gallery-quality tobacco plastome

Session/state source: fresh context; upload `NC_001879.gbk`, the four-region TSV,
the chloroplast specific-color table, and qualifier-priority table. Final state:
Circular `NC_001879.2`; Tuckin; separate strands; functional colors; both radial
label bands; features / `plastome_regions` / GC-content slots; no skew or ticks;
upper-left legend; no plot title.

| Output | Class / decision | Selector or action state | Caption intent | Alt text |
| --- | --- | --- | --- | --- |
| `docs/images/t-gui-05/01-input-ready.png` | `D-OP`; Keep expected | Resize the settings pane, select Circular/GenBank, upload `NC_001879.gbk`, and fill prefix `annotated_chloroplast_map`. | Confirm the complete tobacco plastome input. | Tobacco plastome GenBank file ready in Circular mode |
| `docs/images/t-gui-05/02-first-diagram.png` | `D-RESULT`; Recapture | Generate the plain complete plastome, fit to 70%, and hide the floating feature-search control. | Show the plain plastome before custom tracks. | First circular tobacco plastome diagram |
| `docs/images/t-gui-05/03-annotation-table.png` | `D-OP`; Recapture | Apply Gallery presentation settings; open `Region Annotations`; import the TSV; assert set id `plastome_regions` and four lane controls; scroll the section into view. | Record the four imported structural regions. | Annotation table containing LSC, SSC, IRa, and IRb regions |
| `docs/images/t-gui-05/04-track-settings.png` | `D-OP`; Recapture | Enable custom stack; remove ticks/skew; bind `plastome_regions` to an inside annotations slot at radius 0.65 and width 20 px; set GC content to radius 0.56 and width 0.08; scroll group `Circular track slot plastome_regions` into view. | Record the Gallery three-slot stack. | Circular custom-track controls for features, one plastome-region lane, and GC content |
| `docs/images/t-gui-05/05-finished-diagram.png` | `D-RESULT`; Recapture | Generate and assert 147 logical/197 rendered features, labels LSC/IRb/SSC/IRa, required legend entries, and no skew; fit to 50% and close any feature dialog. | Show the Gallery-equivalent finished figure. | Gallery-style tobacco plastome with functional colors, radial labels, structural regions, GC content, and upper-left legend |

Gallery relation: exact-parity coverage is
`tobacco-chloroplast/manual-08-01-chloroplast-preview.webp`.

### H-GUI-02 — multi-record Circular grid

Session/state source: fresh context; concatenate the four checksum-pinned complete
mitochondrial records into `complete_metazoan_mitochondria.gb`. Final state:
Circular grid; rows `#1@1,#2@1,#3@2,#4@2`; equal size; minimum radius 0.75;
column gap 0.40; row gap 0.08; top title; full definitions; Out gene labels.

| Output | Class / decision | Selector or action state | Caption intent | Alt text |
| --- | --- | --- | --- | --- |
| `docs/images/h-gui-02/grid-settings.png` | `D-OP`; Keep expected | Enable `Multi-Record Canvas`; set the four `Row for #N (record)` selectors and sizing/gap controls; set title/definitions/labels; scroll `Record Size Mode` into view and require row 4 visible. | Record the grid placement choices. | Circular multi-record grid settings |
| `docs/images/h-gui-02/grid-result.png` | `D-RESULT`; Recapture | Generate one four-record SVG; assert IDs, 147 features, labels on every record, and non-overlap; zoom to 30%, pan the complete 2×2 grid into frame, and clear text selection. | Verify the labeled, non-overlapping four-record canvas. | Four complete mitochondrial genomes with gene labels arranged in a grid |

Gallery relation: `Vnig_TUMSAT-TG-2018` is contextual multi-record coverage only;
it uses six Vibrio replicons, Auto sizing, different gaps, a left legend, and a
bottom title, so it is not parity evidence for H-GUI-02.

### H-GUI-03 — Linear rows, ranges, and orientation

Session/state source: fresh context; upload complete `NC_001416.gb` and
`NC_042057.1.gb`. Final state: two Linear rows; ranges 1–48,502 and 1–42,925;
DE3 reverse complemented; 24 px record gap; Above layout; ruler on axis; no
comparison; prefix `linear_regions_orientation`.

| Output | Class / decision | Selector or action state | Caption intent | Alt text |
| --- | --- | --- | --- | --- |
| `docs/images/h-gui-03/record-layout.png` | `D-OP`; Keep expected | Add two sequences; fill both definitions/ranges; check reverse complement only for sequence 2; enable row layout; enter rows 1 and 2 and gap 24; set Above and the ruler; scroll `Arrange linear records in rows` into view. | Show record order, rows, and regions together. | Linear record rows with selected coordinate regions |
| `docs/images/h-gui-03/orientation-result.png` | `D-RESULT`; Recapture | Generate; assert record order, range text `1-48502` / `42925-1`, two non-overlapping rows, no matches, and ruler ticks; fit to 30%. | Verify reverse-complement orientation and the ruler. | Linear diagram with one reversed record and a coordinate ruler |

Gallery relation: `vibrio-harveyi-group-collinear` is contextual multi-row
Linear coverage only; its 11 records, comparison blocks, 48 px gap, bottom
legend, and title-free state are not parity evidence for H-GUI-03.

### H-GUI-09 — depth, GC content, and skew

Session/state source: fresh context; upload complete `AP027133.gb` and the 607-row
`AP027133.DRR394922.depth-1kb.tsv`. Final state: custom slot order features,
ticks, depth, GC content, GC skew; depth range 0–80 with 20/10 ticks; GC percent
axis/ticks; top title; right legend.

| Output | Class / decision | Selector or action state | Caption intent | Alt text |
| --- | --- | --- | --- | --- |
| `docs/images/h-gui-09/track-settings.png` | `D-OP`; Keep expected | Resize the sidebar; open `Depth TSV tracks`; upload depth TSV; set window/step/min/max, axis/ticks/title and tick intervals; set dinucleotide mode Percent with axis/ticks; enable custom stack; scroll the depth section into view. | Show the selected quantitative series and scale. | Depth, GC content, and skew track settings |
| `docs/images/h-gui-09/track-result.png` | `D-RESULT`; Recapture | Generate; assert slot order, depth and percent ticks, `AP027133.1`, title, legends, and 576 features; fit Circular preview to 70%. | Verify series order, axes, and legend. | Genome diagram with depth, GC content, and skew tracks |

Gallery relation: `HmmtDNA_ATskew` is contextual quantitative-track coverage
only; it has no depth series and adds AT skew, so it is not parity evidence for
H-GUI-09.

### H-GUI-10 — region annotations and custom slots

Session/state source: fresh context; upload complete `NC_001879.gbk` and
`nicotiana-tabacum-regions.tsv`. Final state: Middle; annotation lanes 0,1,0,1;
annotations slot outside the axis with labels; AT-skew slot; top title; right
legend.

| Output | Class / decision | Selector or action state | Caption intent | Alt text |
| --- | --- | --- | --- | --- |
| `docs/images/h-gui-10/slot-settings.png` | `D-OP`; Keep expected | Import `plastome_regions`; set legend label and lanes; enable custom stack; change `gc_skew` to AT/`AT skew`; add annotations, bind the set, move it outside, enable labels, and scroll group `Circular track slot annotations` into view. | Show annotation and numeric slot placement. | Custom track-slot settings with annotation and numeric tracks |
| `docs/images/h-gui-10/annotation-result.png` | `D-RESULT`; Recapture | Generate; assert slot order, outside annotation placement, visible LSC/IRb/SSC/IRa and AT-skew legend, and 145 logical/195 rendered features; fit to 70%. | Verify annotation lanes and visible labels. | Genome diagram with region annotations in custom slots |

Gallery relation: `tobacco-chloroplast` is contextual coverage on the same
record and annotation set; its single inside lane, no AT skew, upper-left
legend, and title-free state are not parity evidence for H-GUI-10.

### H-GUI-11 — styles, labels, title, and legend

Session/state source: fresh context; upload `HmmtDNA.gbk`,
`modified_default_colors.tsv`, and `cds_gene_qualifier_priority.tsv`. Final
state: Circular Middle; separate strands; soft_pastels; Out labels; CDS/gene
whitelist pattern `^(COX[123]|ATP[68]|CYTB|ND4L|ND4)$`; top title
`Human mitochondrial genome: selected gene labels`; right legend in the exact
order tRNA, rRNA, CDS, GC content, GC skew (+), GC skew (-).

| Output | Class / decision | Selector or action state | Caption intent | Alt text |
| --- | --- | --- | --- | --- |
| `docs/images/h-gui-11/style-settings.png` | `D-OP`; Recapture | Generate and assert final state first; park feature search; zoom once to 90%; open `Labels`; scroll the priority-file selection and `Whitelist` control into view. | Show the focused style and label choices. | Feature color, label, title, and legend controls |
| `docs/images/h-gui-11/style-result.png` | `D-RESULT`; Recapture | Close `Labels`; reach 60%; pan left to reveal the right legend; assert title, eight selected gene labels, all fills, and exact legend order. | Verify colors, labels, title, and legend together. | Styled genome diagram with custom labels, title, and legend |

Gallery relation: `HmmtDNA_basic_circular` is contextual right-legend coverage
on the same record. Its palette, filter, title, and label subset differ, so it is
not parity evidence for H-GUI-11.

## Linked Gallery result crops

Only these generated-result media are in the eight-scenario crosswalk. Other
Gallery examples remain part of the broader Gallery refresh and visual audit,
but they do not become parity evidence for these documentation scenarios.

| Scenario coverage | Relation | Media / exact session | Required state and capture | Caption / alt | Decision |
| --- | --- | --- | --- | --- | --- |
| T-GUI-01; H-GUI-11 | Exact for T-GUI-01; contextual for H-GUI-11 | `media/HmmtDNA_basic_circular/manual-04-01-final-preview.webp`; `sessions/HmmtDNA_basic_circular.gbdraw-session.json` | Circular single; right legend; Middle; Out; `Homo sapiens`, accession, GC content/skew visible. Viewport 1400×1000; hide preview chrome; reset transform; render SVG at 1000 px width; selector `svg[data-capture-hmmt-basic-final='true']`; padding 12. | Caption: “The completed figure retains labeled features, GC content, GC skew, coordinate ticks, record metadata, and the right legend.” Alt: “Completed circular human mitochondrial genome diagram with labeled features, GC content, GC skew, coordinate ticks, record metadata, and a right legend.” | Recaptured from refreshed exact session and visually verified at 3072×2049. |
| T-GUI-02 | Exact | `media/lambda_basic_linear/manual-05-01-final-preview.webp`; `sessions/lambda_basic_linear.gbdraw-session.json` | Linear one record; comparison mode none; left legend; all labels; ruler; accession/length/CDS visible. Viewport 1900×900; hide preview chrome; reset transform; render SVG at 1370 px width; selector `svg[data-capture-lambda-final='true']`; padding 12. | Caption: “The completed figure retains both strand lanes, all feature labels, the ruler, record metadata, and the left legend.” Alt: “Completed linear lambda phage diagram with two strand lanes, feature labels, a coordinate ruler, record metadata, and a left CDS legend.” | Recaptured from refreshed exact session and visually verified at 4182×1452. |
| T-GUI-05; H-GUI-10 | Exact for T-GUI-05; contextual for H-GUI-10 | `media/tobacco-chloroplast/manual-08-01-chloroplast-preview.webp`; `sessions/tobacco-chloroplast.gbdraw-session.json` | Circular; upper-left legend; set `plastome_regions`; custom stack enabled; LSC/IRb/SSC/IRa and GC content visible. Viewport 1400×1000; hide preview chrome; reset transform; render SVG at 1000 px width; selector `svg[data-capture-chloroplast-final='true']`; padding 12. | Caption: “The completed map shows all four labeled region brackets, chloroplast feature colors, radial gene labels, the inner GC-content track, and one entry per legend category.” Alt: “Completed Nicotiana tabacum chloroplast map with four region brackets, chloroplast feature colors, radial gene labels, an inner GC-content track, record metadata, and one entry per legend category.” | Recaptured from refreshed exact session and visually verified at 3072×2187. |
| H-GUI-02 | Contextual only | `media/Vnig_TUMSAT-TG-2018/manual-06-01-multirecord-preview.webp`; `sessions/Vnig_TUMSAT-TG-2018.gbdraw-session.json.gz` | Require Circular grid, six records, the documented Auto sizing and row positions, left legend, bottom title, and chromosome/plasmid identity text. Use viewport 1600×1100; hide preview chrome; reset pan/zoom; selector `svg[data-capture-vnig-multirecord-final='true']`; padding 12. | Caption: “The generated preview uses automatic radius scaling for two chromosomes and four plasmids.” Alt: “Generated Vibrio multi-record preview with two chromosomes and four plasmids.” | Exact-session metadata added; recaptured and visually verified at 3072×2016. |
| H-GUI-03 | Contextual only | `media/vibrio-harveyi-group-collinear/manual-08-01-collinear-overview.webp`; `sessions/vibrio-harveyi-group-collinear.gbdraw-session.json.gz` | Linear 11 records in five rows; 48 px gap; bottom legend; no plot title; species and Collinear/Inverted text visible. Viewport 1800×1100; selector `svg[data-capture-vibrio-overview='true']`; padding 12. | Caption: “The completed title-free diagram retains all 11 replicons, centers each 18 px bold species name and 16 px strain name on its feature row, preserves the Axis-sized gap around each collinear block band, and places the legend below the final row.” Alt: “Generated five-row Vibrio Harveyi-group diagram with bold species names, 11 replicons, and orientation-and-identity collinear blocks.” | Exact-session assertions strengthened; recaptured and visually verified at 3804×1200. |
| H-GUI-09 | Contextual only | `media/HmmtDNA_ATskew/manual-09-01-atskew-preview.webp`; `sessions/HmmtDNA_ATskew.gbdraw-session.json` | Require Circular single; left legend; enabled slot order `features,gc_content,gc_skew,a_skew_2,ticks`; AT-skew colors/label; `label_in_tick_out`; Homo sapiens and GC/AT skew legend text. Use viewport 1400×1000; hide preview chrome; reset transform; selector `svg[data-capture-hmmt-atskew-final='true']`; padding 12. | Caption: “The generated preview includes the custom track stack and left legend.” Alt: “Generated human mitochondrial AT skew preview with the custom track stack and left legend.” | Exact-session metadata added; recaptured and visually verified at 3072×2304. |

All adjacent `G-OP` controls remain `Keep` unless the strict metadata check or
side-by-side review proves that the operated control/state is stale. Do not
replace a compact control crop merely because its saved result SVG changed.

## Execution order

1. Resolve the renderer blocker and complete WP1–WP6 gates.
2. Run each documentation `--check` command to inventory pixel changes without
   overwriting the committed PNGs.
3. Preserve the 25 old PNGs in `/tmp`, run each scenario through its owner, and
   compare old/new at the same 1440×900 display size.
4. Refresh Gallery sessions and rendered assets with
   `python tools/refresh_gallery_sessions.py`; do not use
   `--skip-session-refresh`.
5. Add or strengthen the three Gallery declarative contracts identified above
   before replacing their WebP files.
6. Preserve each old WebP under `/tmp`, run strict checks, capture only the
   planned result operation at DSF 3 / quality 94, and compare at the Gallery's
   rendered thumbnail size and at natural size.
7. Run the focused documentation, Gallery, session-regeneration, packaging, and
   browser checks. Record actual dimensions and review status in
   `WEB_GALLERY_OPERATION_SCREENSHOT_REGISTER.md`.

Focused documentation commands:

```bash
python docs/capture/run_all.py --scenario T-GUI-01 --tier core --check
python docs/capture/run_all.py --scenario T-GUI-02 --tier core --check
python docs/capture/run_all.py --scenario T-GUI-05 --tier extended --check
python docs/capture/run_all.py --scenario H-GUI-02 --tier extended --check
python docs/capture/run_all.py --scenario H-GUI-03 --tier extended --check
python docs/capture/run_all.py --scenario H-GUI-09 --tier extended --check
python docs/capture/run_all.py --scenario H-GUI-10 --tier extended --check
python docs/capture/run_all.py --scenario H-GUI-11 --tier extended --check
```

After the review-only pass, repeat the same commands without `--check` to let
the owner replace its files, then repeat with `--check` to prove
reproducibility.

Focused Gallery commands after session refresh and metadata updates:

```bash
python tools/capture_gallery_tutorial_screenshots.py --example HmmtDNA_basic_circular --check --strict
python tools/capture_gallery_tutorial_screenshots.py --example lambda_basic_linear --check --strict
python tools/capture_gallery_tutorial_screenshots.py --example tobacco-chloroplast --check --strict
python tools/capture_gallery_tutorial_screenshots.py --example Vnig_TUMSAT-TG-2018 --check --strict
python tools/capture_gallery_tutorial_screenshots.py --example vibrio-harveyi-group-collinear --check --strict
python tools/capture_gallery_tutorial_screenshots.py --example HmmtDNA_ATskew --check --strict
```

Capture a planned result with the same example and an `--operation` filter,
plus `--quality 94 --device-scale-factor 3`; do not run an unreviewed broad
media rewrite.

## Acceptance checks

- All documentation PNGs remain 1440×900 RGB and are below the pinned size
  limit; input-only screenshots remain pixel-stable unless an explained UI
  change is present.
- No painted SVG bounds cross the viewBox. Titles remain centered on primary
  content, docked legends retain the shared spacing policy, and labels or
  floating preview controls do not collide with title, legend, or records.
- The T/H scenario semantic assertions still prove record identities, complete
  lengths, feature counts, slot/row order, orientation, axes, and export safety.
- Exact Gallery parity is claimed only for T-GUI-01, T-GUI-02, and T-GUI-05.
  H-GUI-02/03/09/10/11 mappings remain explicitly contextual.
- Every data-dependent Gallery result uses its own session, state assertions,
  and visible identity text. No cross-example data crop is reused.
- Old/new result images are reviewed side by side at identical displayed size;
  tighter canvas geometry must not reduce label readability or biological
  context.
- Captions and alt text above remain accurate after replacement. If the image
  truth changes, update tutorial text and capture metadata in the same later
  change rather than accepting a misleading image.
- `examples/gbdraw_social_preview.png` and all other manually managed assets
  remain untouched.

Final verification after capture:

```bash
python -m pytest tests/test_documentation_capture_contracts.py \
  tests/test_documentation_scenario_contracts.py \
  tests/test_gui_presentation_capture_contracts.py \
  tests/test_gui_tracks_capture_contracts.py \
  tests/test_gui_nucleotide_comparison_capture_contracts.py -v
node --check gbdraw/web/gallery/gallery.js
npx playwright test tests/web/gallery-tutorial.playwright.spec.js --project=chromium
npx playwright test tests/web/gallery-session-regeneration.playwright.spec.js --project=chromium
```

## Completion evidence

- All manifest-owned documentation screenshots reproduce through their capture
  scenarios; all GUI PNGs are 1440×900 RGB.
- The strict Gallery audit validates 117 media entries, 117 operations, and 117
  operation-media entries after refreshing all 11 sessions.
- Documentation capture and scenario contract tests pass, including bounded
  comparison for the two Chromium complex-SVG raster cases.
- Gallery desktop/mobile tutorial tests and session-regeneration tests pass.
- Circular radial-layout tests explicitly require the plus strand to be outside
  the minus strand for separated `Tuckin`, `Middle`, and `Spreadout` layouts.
- Renderer output-comparison and circular feature-width regression suites pass.
