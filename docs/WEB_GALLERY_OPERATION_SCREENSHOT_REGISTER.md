# Web Gallery Operation Screenshot Register

This register controls operation-level screenshot coverage for the seven ready Gallery tutorials. It is the source of truth for capture order, expected filenames, and whether an existing crop can be reused during JSON migration.

Conventions:

- Source values: Gallery page, Web app session, generated preview, popup, Files tab.
- `reuse` means an existing screenshot may be attached to an operation if it still meets the crop and accessibility rules.
- `new` means capture a new WebP under `gbdraw/web/gallery/media/<example-id>/`.
- Highlight notes are selector-oriented where possible; exact selectors may be added to tutorial JSON `capture` metadata during migration.

## hepatoplasmataceae_collinear

| Section | Step | Op | Target UI area | Source | Expected media | Highlight / note | Status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| quick | 1 | 1 | Gallery action toolbar (`.viewer-actions`) with Session visible | Gallery page | `quick-01-01-download-session.webp` | actual `.viewer-actions` toolbar crop, not an isolated button | reuse `01-gallery-session.webp` if crop remains valid |
| quick | 2 | 1 | Web app header/session toolbar (`.app-toolbar`) | Web app session | `quick-02-01-open-session-controls.webp` | actual `.app-toolbar` crop with Load Session visible | new |
| quick | 3 | 1 | Restored generated diagram | generated preview | `quick-03-01-restored-diagram.webp` | preview viewport | new |
| manual | 1 | 1 | Mode segmented control | Web app session | `manual-01-01-linear-mode.webp` | Linear selected | new |
| manual | 2 | 1 | Linear GenBank uploader with file order | Web app session | `manual-02-01-five-file-upload-order.webp` | file chips in order | new |
| manual | 3 | 1 | Pairwise Comparisons panel | Web app session | `manual-03-01-open-pairwise.webp` | expanded Pairwise Comparisons | reuse `02-web-linear-losat-settings.webp` if crop remains valid |
| manual | 3 | 2 | LOSAT run settings | Web app session | `manual-03-02-browser-losat.webp` | Run LOSAT, LOSATP | reuse `02-web-linear-losat-settings.webp` if crop remains valid |
| manual | 4 | 1 | LOSATP collinear controls | Web app session | `manual-04-01-collinear-reduction.webp` | Collinear blocks selected | new |
| manual | 4 | 2 | Collinear color mode | Web app session | `manual-04-02-orientation-identity.webp` | Orientation + identity | new |
| manual | 5 | 1 | Linear layout controls | Web app session | `manual-05-01-layout-gc-skew-ruler.webp` | align center, GC/skew, ruler | new |
| manual | 6 | 1 | Generate button and processing context | Web app session | `manual-06-01-generate.webp` | Generate Diagram button | new |
| manual | 7 | 1 | Collinear block popup | popup | `manual-07-01-collinear-block-popup.webp` | ribbon + popup | reuse `03-collinearity-popup.webp` |
| post | 1 | 1 | Feature popup | popup | `post-01-01-feature-popup.webp` | feature + popup | new |
| post | 2 | 1 | Legend editor | Web app session | `post-02-01-legend-editor.webp` | legend editor drawer | new |
| post | 3 | 1 | Preview export toolbar (`.preview-toolbar`) | Web app session | `post-03-01-export-controls.webp` | actual `.preview-toolbar` crop with DPI, SVG, Interactive SVG, PNG, and PDF controls | new |
| post | 4 | 1 | Web app header/session toolbar (`.app-toolbar`) | Web app session | `post-04-01-save-session.webp` | actual `.app-toolbar` crop with Save Session visible | new |

## hepatoplasmataceae_orthogroup

| Section | Step | Op | Target UI area | Source | Expected media | Highlight / note | Status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| quick | 1 | 1 | Gallery action toolbar (`.viewer-actions`) with Session visible | Gallery page | `quick-01-01-download-session.webp` | actual `.viewer-actions` toolbar crop, not an isolated button | reuse `01-gallery-session.webp` if crop remains valid |
| quick | 2 | 1 | Web app header/session toolbar (`.app-toolbar`) | Web app session | `quick-02-01-open-session-controls.webp` | actual `.app-toolbar` crop with Load Session visible | new |
| quick | 3 | 1 | Contrast with collinear example | Gallery page | `quick-03-01-compare-collinear.webp` | sample list entries | new |
| manual | 1 | 1 | Mode segmented control | Web app session | `manual-01-01-linear-mode.webp` | Linear selected | new |
| manual | 2 | 1 | Linear GenBank uploader with file order | Web app session | `manual-02-01-five-file-upload-order.webp` | file chips in order | new |
| manual | 3 | 1 | Pairwise Comparisons panel | Web app session | `manual-03-01-browser-losat.webp` | Run LOSAT, LOSATP | new |
| manual | 4 | 1 | Orthogroups mode controls | Web app session | `manual-04-01-orthogroups-mode.webp` | blastp mode Orthogroups | new |
| manual | 4 | 2 | Orthogroup membership and hit settings | Web app session | `manual-04-02-orthogroup-settings.webp` | membership, max hits | new |
| manual | 5 | 1 | Linear layout controls | Web app session | `manual-05-01-layout-gc-skew-ruler.webp` | align center, GC/skew, ruler | new |
| manual | 6 | 1 | Generate button and processing context | Web app session | `manual-06-01-generate.webp` | Generate Diagram button | new |
| manual | 7 | 1 | Orthogroup overview in preview | generated preview | `manual-07-01-orthogroup-overview.webp` | ribbon-dense region | reuse `02-orthogroup-preview.webp` |
| manual | 8 | 1 | Orthogroup match popup | popup | `manual-08-01-orthogroup-popup.webp` | ribbon + popup | reuse `03-orthogroup-popup.webp` |
| post | 1 | 1 | Pairwise ribbon popup | popup | `post-01-01-match-popup.webp` | ribbon + popup | reuse `03-orthogroup-popup.webp` if crop remains valid |
| post | 2 | 1 | CDS feature popup | popup | `post-02-01-feature-popup.webp` | feature + popup | new |
| post | 3 | 1 | Legend editor | Web app session | `post-03-01-legend-editor.webp` | legend editor drawer | new |
| post | 4 | 1 | Preview export toolbar (`.preview-toolbar`) | Web app session | `post-04-01-export-controls.webp` | actual `.preview-toolbar` crop with DPI, SVG, Interactive SVG, PNG, and PDF controls | new |
| post | 5 | 1 | Web app header/session toolbar (`.app-toolbar`) | Web app session | `post-05-01-save-session.webp` | actual `.app-toolbar` crop with Save Session visible | new |

## BGC0000708-BGC0000713

| Section | Step | Op | Target UI area | Source | Expected media | Highlight / note | Status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| quick | 1 | 1 | Gallery action toolbar (`.viewer-actions`) with Session visible | Gallery page | `quick-01-01-download-session.webp` | actual `.viewer-actions` toolbar crop, not an isolated button | reuse `01-gallery-session.webp` if crop remains valid |
| quick | 2 | 1 | Web app header/session toolbar (`.app-toolbar`) | Web app session | `quick-02-01-open-session-controls.webp` | actual `.app-toolbar` crop with Load Session visible | new |
| quick | 3 | 1 | Restored BGC settings summary | Web app session | `quick-03-01-restored-settings.webp` | Pairwise and color settings | new |
| manual | 1 | 1 | Mode segmented control | Web app session | `manual-01-01-linear-mode.webp` | Linear selected | new |
| manual | 2 | 1 | Linear GenBank uploader with five BGC files | Web app session | `manual-02-01-five-file-upload-order.webp` | file chips in order | new |
| manual | 3 | 1 | Pairwise Comparisons panel | Web app session | `manual-03-01-open-pairwise.webp` | expanded panel | new |
| manual | 3 | 2 | LOSATP orthogroups mode | Web app session | `manual-03-02-select-losatp-orthogroups.webp` | Orthogroups selected | new |
| manual | 3 | 3 | Threshold controls | Web app session | `manual-03-03-set-thresholds.webp` | identity, bitscore, e-value | new |
| manual | 4 | 1 | Orthogroup alignment field | Web app session | `manual-04-01-align-og1.webp` | `og_1` value | new |
| manual | 4 | 2 | Reverse-complement selector | Web app session | `manual-04-02-reverse-bgc0000713.webp` | BGC0000713 selected | new |
| manual | 5 | 1 | Label and ruler layout controls | Web app session | `manual-05-01-label-layout.webp` | first-record labels, ruler, legend bottom | new |
| manual | 6 | 1 | Title and record text controls | Web app session | `manual-06-01-title-record-text.webp` | plot title and labels | new |
| manual | 7 | 1 | Rule 1: core biosynthetic genes | Web app session | `manual-07-01-color-rule-core.webp` | gene_kind biosynthetic$ | new |
| manual | 7 | 2 | Rule 2: additional biosynthetic genes | Web app session | `manual-07-02-color-rule-additional.webp` | gene_kind biosynthetic-additional | new |
| manual | 7 | 3 | Rule 3: transport genes | Web app session | `manual-07-03-color-rule-transport.webp` | gene_kind transport | new |
| manual | 7 | 4 | Rule 4: regulatory genes | Web app session | `manual-07-04-color-rule-regulatory.webp` | gene_kind regulatory | new |
| manual | 8 | 1 | Generated BGC comparison preview | generated preview | `manual-08-01-bgc-preview.webp` | five rows + bottom legend | reuse `02-bgc-preview.webp` |
| manual | 9 | 1 | Orthogroup ribbon popup | popup | `manual-09-01-orthogroup-popup.webp` | ribbon + popup | reuse `03-orthogroup-popup.webp` |
| manual | 10 | 1 | Colored feature popup | popup | `manual-10-01-feature-popup.webp` | feature + popup | reuse `04-feature-popup.webp` |
| post | 1 | 1 | Colored feature popup metadata | popup | `manual-10-01-feature-popup.webp` | gene_kind and orthogroup fields | reuse manual popup crop |
| post | 2 | 1 | Legend editor | Web app session | `post-02-01-legend-editor.webp` | legend editor drawer | new |
| post | 3 | 1 | Preview export toolbar (`.preview-toolbar`) | Web app session | `post-03-01-export-controls.webp` | actual `.preview-toolbar` crop with DPI, SVG, Interactive SVG, PNG, and PDF controls | new |
| post | 4 | 1 | Web app header/session toolbar (`.app-toolbar`) | Web app session | `post-04-01-save-session.webp` | actual `.app-toolbar` crop with Save Session visible | new |

## majanivirus_orthogroup

| Section | Step | Op | Target UI area | Source | Expected media | Highlight / note | Status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| quick | 1 | 1 | Gallery action toolbar (`.viewer-actions`) with Session visible | Gallery page | `quick-01-01-download-session.webp` | actual `.viewer-actions` toolbar crop, not an isolated button | reuse `01-gallery-session.webp` if crop remains valid |
| quick | 2 | 1 | Web app header/session toolbar (`.app-toolbar`) | Web app session | `quick-02-01-open-session-controls.webp` | actual `.app-toolbar` crop with Load Session visible | new |
| quick | 3 | 1 | Restored large-session preview | generated preview | `quick-03-01-restored-large-session.webp` | preview viewport | new |
| manual | 1 | 1 | Mode segmented control | Web app session | `manual-01-01-linear-mode.webp` | Linear selected | new |
| manual | 2 | 1 | Linear GenBank uploader with nine viral files | Web app session | `manual-02-01-nine-file-upload-order.webp` | file chips in order | new |
| manual | 3 | 1 | Pairwise Comparisons panel | Web app session | `manual-03-01-open-pairwise.webp` | expanded panel | new |
| manual | 3 | 2 | LOSATP orthogroups mode | Web app session | `manual-03-02-losatp-orthogroups.webp` | Orthogroups selected | new |
| manual | 3 | 3 | Large-session thread and threshold settings | Web app session | `manual-03-03-thread-threshold-settings.webp` | threads, bitscore, e-value, identity | new |
| manual | 4 | 1 | Linear layout controls | Web app session | `manual-04-01-linear-layout.webp` | align center and compact layout | new |
| manual | 5 | 1 | Record label editor | Web app session | `manual-05-01-record-labels.webp` | nine descriptive labels | new |
| manual | 6 | 1 | Rule 1: WSSV-like proteins | Web app session | `manual-06-01-color-rule-wssv-like.webp` | product regex | new |
| manual | 6 | 2 | Rule 2: BIRP | Web app session | `manual-06-02-color-rule-birp.webp` | product regex | new |
| manual | 6 | 3 | Rule 3: tyrosine recombinase | Web app session | `manual-06-03-color-rule-recombinase.webp` | product regex | new |
| manual | 6 | 4 | Default CDS color | Web app session | `manual-06-04-default-cds-color.webp` | default CDS swatch | new |
| manual | 7 | 1 | Generated orthogroup preview | generated preview | `manual-07-01-orthogroup-preview.webp` | nine rows + right legend | reuse `02-orthogroup-preview.webp` |
| manual | 8 | 1 | Orthogroup ribbon popup | popup | `manual-08-01-orthogroup-popup.webp` | ribbon + popup | reuse `03-orthogroup-popup.webp` |
| manual | 9 | 1 | Gallery tab bar plus Files panel (`.tab-bar` + `#files-panel`) | Files tab | `manual-09-01-files-tab.webp` | actual `.tab-bar` crop with Files selected plus input names and artifacts | reuse `04-files.webp` |
| post | 1 | 1 | Match ribbon popup | popup | `post-01-01-match-popup.webp` | ribbon + popup | reuse `03-orthogroup-popup.webp` if crop remains valid |
| post | 2 | 1 | CDS feature popup | popup | `post-02-01-feature-popup.webp` | feature + popup | new |
| post | 3 | 1 | Identity legend in preview | generated preview | `post-03-01-identity-legend.webp` | identity legend | new |
| post | 4 | 1 | Preview export toolbar (`.preview-toolbar`) | Web app session | `post-04-01-export-controls.webp` | actual `.preview-toolbar` crop with DPI, SVG, Interactive SVG, PNG, and PDF controls | new |
| post | 5 | 1 | Web app header/session toolbar (`.app-toolbar`) | Web app session | `post-05-01-save-session.webp` | actual `.app-toolbar` crop with Save Session visible | new |

## WSSV_genome_comparison

| Section | Step | Op | Target UI area | Source | Expected media | Highlight / note | Status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| quick | 1 | 1 | Gallery action toolbar (`.viewer-actions`) with Session visible | Gallery page | `quick-01-01-download-session.webp` | actual `.viewer-actions` toolbar crop, not an isolated button | reuse `01-gallery-session.webp` if crop remains valid |
| quick | 2 | 1 | Web app header/session toolbar (`.app-toolbar`) | Web app session | `quick-02-01-open-session-controls.webp` | actual `.app-toolbar` crop with Load Session visible | new |
| quick | 3 | 1 | Restored conservation stack | generated preview | `quick-03-01-restored-conservation-stack.webp` | preview viewport | new |
| manual | 1 | 1 | Mode segmented control | Web app session | `manual-01-01-circular-mode.webp` | Circular selected | new |
| manual | 2 | 1 | Circular GenBank uploader | Web app session | `manual-02-01-reference-genbank-upload.webp` | AP027280.gb file chip | new |
| manual | 3 | 1 | Basic circular layout controls | Web app session | `manual-03-01-basic-circular-layout.webp` | spreadout, legend left, suppress GC/skew | new |
| manual | 4 | 1 | Circular Conservation panel | Web app session | `manual-04-01-open-conservation-panel.webp` | panel expanded | new |
| manual | 4 | 2 | Conservation FASTA uploader | Web app session | `manual-04-02-upload-fasta-comparisons.webp` | comparison FASTA list | new |
| manual | 5 | 1 | Ring label controls | Web app session | `manual-05-01-ring-labels.webp` | labels textarea/list | new |
| manual | 6 | 1 | Ring color controls | Web app session | `manual-06-01-ring-colors.webp` | color swatches | new |
| manual | 7 | 1 | Conservation thresholds | Web app session | `manual-07-01-thresholds.webp` | bitscore, e-value, identity, length | new |
| manual | 8 | 1 | Browser LOSAT execution controls | Web app session | `manual-08-01-browser-losat-run.webp` | Auto, Safe, blastn megablast | new |
| manual | 9 | 1 | Generated conservation rings | generated preview | `manual-09-01-conservation-rings.webp` | 20 rings + legend | reuse `03-conservation-rings.webp` |
| manual | 10 | 1 | Gallery tab bar plus Files panel (`.tab-bar` + `#files-panel`) | Files tab | `manual-10-01-files-tab.webp` | actual `.tab-bar` crop with Files selected plus input names and artifacts | reuse `02-input-files.webp` if crop is reduced |
| manual | 11 | 1 | CDS feature popup | popup | `manual-11-01-feature-popup.webp` | outer CDS + popup | reuse `04-feature-popup.webp` |
| post | 1 | 1 | Floating preview toolbar (`.preview-controls`) | Web app session | `post-01-01-zoom-controls.webp` | actual floating toolbar crop with zoom, reset layout, layout edit, and canvas padding controls | new |
| post | 2 | 1 | CDS feature popup metadata | popup | `manual-11-01-feature-popup.webp` | note, product, protein_id fields | reuse manual popup crop |
| post | 3 | 1 | Conservation legend | generated preview | `post-03-01-conservation-legend.webp` | identity endpoints visible | new |
| post | 4 | 1 | Preview export toolbar (`.preview-toolbar`) | Web app session | `post-04-01-export-controls.webp` | actual `.preview-toolbar` crop with DPI, SVG, Interactive SVG, PNG, and PDF controls | new |
| post | 5 | 1 | Web app header/session toolbar (`.app-toolbar`) | Web app session | `post-05-01-save-session.webp` | actual `.app-toolbar` crop with Save Session visible | new |

## HmmtDNA_ATskew

| Section | Step | Op | Target UI area | Source | Expected media | Highlight / note | Status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| quick | 1 | 1 | Gallery action toolbar (`.viewer-actions`) with Session visible | Gallery page | `quick-01-01-download-session.webp` | actual `.viewer-actions` toolbar crop, not an isolated button | reuse `01-gallery-session.webp` if crop remains valid |
| quick | 2 | 1 | Web app header/session toolbar (`.app-toolbar`) | Web app session | `quick-02-01-open-session-controls.webp` | actual `.app-toolbar` crop with Load Session visible | new |
| quick | 3 | 1 | Restored custom track stack | Web app session | `quick-03-01-restored-custom-stack.webp` | custom track slots | new |
| manual | 1 | 1 | Mode segmented control | Web app session | `manual-01-01-circular-mode.webp` | Circular selected | new |
| manual | 2 | 1 | Circular GenBank uploader | Web app session | `manual-02-01-genbank-upload.webp` | HmmtDNA.gbk file chip | new |
| manual | 3 | 1 | Basic circular layout controls | Web app session | `manual-03-01-basic-circular-layout.webp` | middle, left legend, labels out | new |
| manual | 4 | 1 | Dinucleotide window controls | Web app session | `manual-04-01-window-step.webp` | GC, 500, 50 | new |
| manual | 5 | 1 | Custom Track Slots panel | Web app session | `manual-05-01-open-custom-track-slots.webp` | panel expanded | new |
| manual | 5 | 2 | Middle template slot stack | Web app session | `manual-05-02-middle-template-stack.webp` | features, GC, GC skew slots | new |
| manual | 6 | 1 | Add AT skew slot | Web app session | `manual-06-01-at-skew-slot-values.webp` | AT, label, width, colors | reuse `02-atskew-preview.webp` only for result, not control |
| manual | 7 | 1 | Tick track slot | Web app session | `manual-07-01-tick-track-slot.webp` | ticks renderer and layout | new |
| manual | 8 | 1 | CDS label priority controls | Web app session | `manual-08-01-cds-label-priority.webp` | CDS, gene | new |
| manual | 9 | 1 | Generated AT skew preview | generated preview | `manual-09-01-atskew-preview.webp` | track stack + left legend | reuse `02-atskew-preview.webp` |
| manual | 10 | 1 | Feature popup | popup | `manual-10-01-feature-popup.webp` | feature + popup | reuse `03-feature-popup.webp` |
| post | 1 | 1 | Feature popup metadata | popup | `manual-10-01-feature-popup.webp` | gene, product, db_xref fields | reuse manual popup crop |
| post | 2 | 1 | Legend editor | Web app session | `post-02-01-legend-editor.webp` | AT skew entries visible | new |
| post | 3 | 1 | Preview export toolbar (`.preview-toolbar`) | Web app session | `post-03-01-export-controls.webp` | actual `.preview-toolbar` crop with DPI, SVG, Interactive SVG, PNG, and PDF controls | new |
| post | 4 | 1 | Web app header/session toolbar (`.app-toolbar`) | Web app session | `post-04-01-save-session.webp` | actual `.app-toolbar` crop with Save Session visible | new |

## Vnig_TUMSAT-TG-2018

| Section | Step | Op | Target UI area | Source | Expected media | Highlight / note | Status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| quick | 1 | 1 | Gallery action toolbar (`.viewer-actions`) with Session visible | Gallery page | `quick-01-01-download-session.webp` | actual `.viewer-actions` toolbar crop, not an isolated button | reuse `01-gallery-session.webp` if crop remains valid |
| quick | 2 | 1 | Web app header/session toolbar (`.app-toolbar`) | Web app session | `quick-02-01-open-session-controls.webp` | actual `.app-toolbar` crop with Load Session visible | new |
| quick | 3 | 1 | Restored multi-record layout | generated preview | `quick-03-01-restored-multirecord-layout.webp` | preview viewport | new |
| manual | 1 | 1 | Mode segmented control | Web app session | `manual-01-01-circular-mode.webp` | Circular selected | new |
| manual | 2 | 1 | Multi-record GBFF uploader | Web app session | `manual-02-01-gbff-upload.webp` | GBFF file chip | new |
| manual | 2 | 2 | Detected record list | Web app session | `manual-02-02-record-list.webp` | chromosomes and plasmids | new |
| manual | 3 | 1 | Basic circular layout controls | Web app session | `manual-03-01-basic-circular-layout.webp` | tuckin, orchid, strands, legend | new |
| manual | 4 | 1 | Multi-record canvas controls | Web app session | `manual-04-01-multirecord-canvas.webp` | enabled, auto size | new |
| manual | 4 | 2 | Record position controls | Web app session | `manual-04-02-record-positions.webp` | `#1@1` through `#6@2` | new |
| manual | 5 | 1 | Bottom title controls | Web app session | `manual-05-01-bottom-title.webp` | title text and bottom position | new |
| manual | 6 | 1 | Generated multi-record preview | generated preview | `manual-06-01-multirecord-preview.webp` | two chromosomes + plasmids | reuse `02-multirecord-preview.webp` |
| manual | 7 | 1 | Gallery tab bar plus Files panel (`.tab-bar` + `#files-panel`) | Files tab | `manual-07-01-files-tab.webp` | actual `.tab-bar` crop with Files selected plus GBFF and artifacts | reuse `03-files.webp` |
| manual | 8 | 1 | CDS feature popup | popup | `manual-08-01-feature-popup.webp` | feature + popup | reuse `04-feature-popup.webp` |
| post | 1 | 1 | Feature popup with record ownership | popup | `post-01-01-record-feature-popup.webp` | record ID visible | reuse `04-feature-popup.webp` if crop remains valid |
| post | 2 | 1 | Compact left legend | generated preview | `post-02-01-left-legend-compactness.webp` | legend + chromosome 1 edge | new |
| post | 3 | 1 | Preview export toolbar (`.preview-toolbar`) | Web app session | `post-03-01-export-controls.webp` | actual `.preview-toolbar` crop with DPI, SVG, Interactive SVG, PNG, and PDF controls | new |
| post | 4 | 1 | Web app header/session toolbar (`.app-toolbar`) | Web app session | `post-04-01-save-session.webp` | actual `.app-toolbar` crop with Save Session visible | new |
