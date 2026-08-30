[Documentation home](./DOCS.md) | [Tutorials](./TUTORIALS/README.md) | [Technical documentation](./REFERENCE/README.md) | **FAQ** | [Gallery](./GALLERY.md)

# Frequently asked questions

## Choosing a workflow

### Should I use the hosted app, local GUI, CLI, or Python?

Use the hosted app for interactive work without installing Python. Use
`gbdraw gui` for the same browser interface from a local installation, the CLI
for repeatable commands and batches, and Python when a pipeline already owns
Biopython records or needs output bytes. Integrations that need explicit
planning and session conversion can use the typed request API.

See the exact boundaries for the [web
app](./REFERENCE/web-app.md#execution-privacy-and-offline-use), [command
line](./REFERENCE/command-line.md), [package-root Python
API](./REFERENCE/python-api.md), and [typed
requests](./REFERENCE/typed-requests.md). The [Tutorials](./TUTORIALS/README.md)
provide a first project for each surface.

### Should I use Circular or Linear?

Use Circular for a complete circular replicon and radial tracks. Circular
multi-record layout places complete circular records on one canvas. Use Linear
when record order, cropped regions, reverse-complement display, rows, or
pairwise links are central to the figure.

Changing modes does not translate mode-specific placement. Circular grid
positions do not become Linear rows, and Linear row assignments, crops, and
reverse-complement states do not become Circular settings.

The [diagram-layout contract](./REFERENCE/palettes-feature-rules-labels-shapes-and-tracks.md#diagram-layout)
defines the two layouts. Start with the [first Circular
Tutorial](./TUTORIALS/GUI/first-circular-genome-diagram.md) or [first Linear
Tutorial](./TUTORIALS/GUI/first-linear-genome-diagram.md).

### Which comparison method should I use?

Use an uploaded BLAST table when the search ran elsewhere or the evidence must
stay fixed. LOSATN compares nucleotide sequence; TLOSATX is useful when coding
similarity remains after nucleotide divergence. In Linear **Settings**,
the three **LOSAT Mode** buttons choose LOSATN, LOSATP, or TLOSATX. When LOSATP is selected,
**LOSATP mode** chooses Similarity groups, Collinear blocks, or Pairwise
matches. Pairwise matches show individual protein matches, Similarity groups
show search-derived membership, and Collinear blocks emphasize compatible
ordered anchors. Similarity groups always searches all loaded record pairs.
Fresh Collinear settings and **Reset Settings** default to **Adjacent pairs**,
while a saved **All records** scope remains explicit. Circular rings place
retained evidence around one reference. Linear comparisons connect selected
query and subject record endpoints.

The [comparison capability
matrix](./REFERENCE/comparison-programs-thresholds-and-results.md#capability-matrix)
lists interface availability, filters, direction, and scientific limits.

### Is my data uploaded, and can I work offline?

The hosted page loads from `gbdraw.app`, but genome parsing, LOSAT searches,
rendering, session handling, and export run in the browser without a gbdraw
application-server upload. The hosted site may collect aggregate page-use
analytics; genome files and generated diagrams are not analytics payloads.
`gbdraw gui` can run the normal workflow offline after installation.

Sessions embed input resources and should be protected like the source data.
See [execution, privacy, and offline
use](./REFERENCE/web-app.md#execution-privacy-and-offline-use) for the full
boundary and browser-performance constraints.

### How should I prepare a figure for publication?

Keep an SVG master when the submission workflow accepts vector artwork, proof
the final derivative at its placement size, and retain the inputs, software
versions, comparison evidence, render instructions, and any manual-edit notes.

The [publication and reproducible-handoff
contract](./REFERENCE/output-formats-and-export.md#publication-and-reproducible-handoff)
lists the format choices and archive contents.

## Inputs, layout, and presentation

### Why do my CLI and browser renders differ slightly?

Small differences in label placement and legend sizing are expected. The CLI
uses kerning-aware font metrics, while the web app uses browser text metrics.

### How do I hide the coordinate scale without hiding the genome axis?

Use `--hide_scale` on the command line or clear **Show Coordinate Scale** in
the web app. Explicit custom track slots can own a separate `ticks` renderer,
and quantitative tracks own their own axes. See [diagram
layout](./REFERENCE/palettes-feature-rules-labels-shapes-and-tracks.md#diagram-layout)
and the [CLI layout
options](./REFERENCE/command-line.md#record-selection-and-layout).

### Can I use a GFF3 file by itself?

No. Pair each GFF3 input with the matching FASTA sequence. The [sequence and
annotation contract](./REFERENCE/input-formats-and-tsv-schemas.md#sequence-and-annotation-files)
defines ID matching, coordinates, phase, and translation behavior.

### My labels overlap. What should I do?

Reduce the label size, filter the label set, or change the Circular track
placement. The [feature-presentation
contract](./REFERENCE/palettes-feature-rules-labels-shapes-and-tracks.md#feature-presentation)
defines qualifier priority, whitelist and blacklist behavior, overrides, and
overlap controls.

### How do I change the color of one specific gene?

Use a specific-color rule that matches the intended feature qualifier. See the
[styling-table schemas](./REFERENCE/input-formats-and-tsv-schemas.md#styling-tables)
and [feature-rule precedence](./REFERENCE/palettes-feature-rules-labels-shapes-and-tracks.md#feature-presentation).

### Why does a web color edit create a qualifier rule for some labels and `hash` rules for others?

The editor uses one exact qualifier rule only when that rule represents the
selected scope without ambiguity. Otherwise it keeps feature-specific rules;
identical duplicate biological identities cannot be reconstructed as one
instance after regeneration. See [preview, search, and
editor](./REFERENCE/web-app.md#preview-search-and-editor) for the exact
conditions and scope controls.

### How do I mark a coordinate range or a group of features?

Use an annotation table and bind its `set_id` to an `annotations` track
slot. See [annotation table
fields](./REFERENCE/input-formats-and-tsv-schemas.md#annotation-table-fields)
and [track placement](./REFERENCE/palettes-feature-rules-labels-shapes-and-tracks.md#tracks-axes-and-annotations).

### Can I use gene names instead of product descriptions for labels?

Yes. Set qualifier priority so that `gene` precedes `product`. The
[feature-presentation
contract](./REFERENCE/palettes-feature-rules-labels-shapes-and-tracks.md#feature-presentation)
and [qualifier-priority
schema](./REFERENCE/input-formats-and-tsv-schemas.md#styling-tables) define the
behavior.

## Comparisons and sessions

### My comparative diagram has no links. What should I check?

Check the table format, record direction and identifiers, display thresholds,
selected edges, and protein translations. Circular rings also require the
intended reference side. The [result meanings and
limits](./REFERENCE/comparison-programs-thresholds-and-results.md#result-meanings-and-limits)
section lists the empty-result cases.

### How do I draw several Linear records without comparing them in the web app?

Fresh Linear pages and **Reset Settings** start with **No comparison**. If a
comparison is active, select **No comparison** in the **Comparison** command
group. The three buttons there apply one choice to every adjacent pair; the
separate **Current:** status reports the effective plan. Open **Selected pairs
(N)** to inspect or edit a custom pair plan. Choosing **No comparison** keeps
retained pair files and raw-result names inactive for later reuse. See the [web
comparison controls](./REFERENCE/web-app.md#comparison-surfaces) and
[selected-edge contract](./REFERENCE/comparison-programs-thresholds-and-results.md#selected-linear-edges).

### Why did gbdraw rerun LOSATP after I loaded a session?

A cached search is reused only when its biological inputs, direction, program,
and meaningful search settings still match. See [saved comparison results and
cache reuse](./REFERENCE/session-and-request-compatibility.md#saved-comparison-results-and-cache-reuse).

### Why does Save Raw LOSAT TSV not contain the internal `h_` IDs?

Generated protein results export stable readable aliases instead of
session-only handles. Uploaded tables are not rewritten. See [raw results and
cache identity](./REFERENCE/comparison-programs-thresholds-and-results.md#raw-results-and-cache-identity).

### Can pairwise comparison links be curved?

Yes. `curve` bends the same mapped spans that `ribbon` draws as filled
links; it does not change the evidence. See [result meanings and
limits](./REFERENCE/comparison-programs-thresholds-and-results.md#result-meanings-and-limits).

## Depth, content, and skew

### What if one record has no Depth TSV for a sample?

Keep an explicit missing entry for that record. gbdraw does not substitute
another file or draw zero coverage. See [tracks, axes, and
annotations](./REFERENCE/palettes-feature-rules-labels-shapes-and-tracks.md#tracks-axes-and-annotations)
and the [CLI depth input](./REFERENCE/command-line.md#tracks-annotations-and-presentation).

### How do I make the GC content track smoother?

Use a larger window and step. Larger values smooth the trace and reduce local
detail; the [numeric-track
contract](./REFERENCE/palettes-feature-rules-labels-shapes-and-tracks.md#tracks-axes-and-annotations)
defines the relationship.

### Can I plot AT instead of GC?

Yes. Select the `AT` dinucleotide. The [tracks, axes, and annotations
contract](./REFERENCE/palettes-feature-rules-labels-shapes-and-tracks.md#tracks-axes-and-annotations)
defines content and skew behavior for reversed pairs.

## Output and limitations

### Why does SVG export work but PNG/PDF/EPS/PS export fail?

Command-line and Python conversion to those formats requires CairoSVG and its
runtime libraries. See [output names, dependencies, and
overwrite](./REFERENCE/output-formats-and-export.md#names-dependencies-and-overwrite).

### Are there known visualization or conversion limitations?

Yes. See the [feature-presentation
limits](./REFERENCE/palettes-feature-rules-labels-shapes-and-tracks.md#feature-presentation)
and [static, interactive, and raster
output](./REFERENCE/output-formats-and-export.md#static-interactive-and-raster-output).

[Documentation home](./DOCS.md) | [Tutorials](./TUTORIALS/README.md) | [Technical documentation](./REFERENCE/README.md) | **FAQ** | [Gallery](./GALLERY.md)
