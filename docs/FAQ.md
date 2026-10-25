[Documentation home](./DOCS.md) | [Tutorials](./TUTORIALS/README.md) | [Technical documentation](./REFERENCE/README.md) | **FAQ** | [Gallery](./GALLERY.md)

# Frequently asked questions

## Choosing a workflow

### Should I use the hosted app, local GUI, CLI, or Python?

Use the hosted app for interactive work without installing Python. Use
`gbdraw gui` for the same browser interface from a local installation, the CLI
for repeatable commands and batches, and Python when your pipeline already
holds Biopython records or needs the output bytes. If you integrate gbdraw and
need explicit planning and session conversion, use the typed request API.

The exact differences are in the [web
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

Switching modes does not carry placement over. Circular grid positions do not
become Linear rows, and Linear row assignments, crops, and reverse-complement
states do not become Circular settings.

[Diagram layout](./REFERENCE/palettes-feature-rules-labels-shapes-and-tracks.md#diagram-layout)
describes the two layouts. Start with the [first Circular
Tutorial](./TUTORIALS/first-circular-genome-diagram.md) or [first Linear
Tutorial](./TUTORIALS/first-linear-genome-diagram.md).

### Which comparison method should I use?

Use an uploaded BLAST table when the search ran elsewhere or the hits must not
change. LOSATN compares nucleotide sequence. TLOSATX is useful when coding
similarity remains after nucleotide divergence.

In Linear **Settings**, the three **LOSAT Mode** buttons choose LOSATN, LOSATP,
or TLOSATX. When LOSATP is selected, **LOSATP mode** chooses Similarity groups,
Collinear blocks, or Pairwise matches:

- Pairwise matches show individual protein matches.
- Similarity groups show membership derived from the search. They always
  search all loaded record pairs.
- Collinear blocks emphasize compatible ordered anchors. Fresh Collinear
  settings and **Reset Settings** default to **Adjacent pairs**, while a saved
  **All records** scope stays as saved.

Circular rings place the hits that pass the filters around one reference.
Linear comparisons connect selected query and subject record endpoints.

On the command line, `gbdraw linear --losat losatn`, `--losat tlosatx`, or
`--losat losatp` runs the same Linear searches. `--comparisons_table` rows with
`source` `losat` or `table` mix searched and uploaded edges. `gbdraw circular
--losat losatn` or `--losat tlosatx` with `--conservation_sequence` (FASTA,
GenBank, or DDBJ) runs the ring searches against the displayed reference.

The [comparison capability
matrix](./REFERENCE/comparison-programs-thresholds-and-results.md#capability-matrix)
lists interface availability, filters, direction, and scientific limits.

### Is my data uploaded, and can I work offline?

Your genome files are not uploaded. The hosted page loads from `gbdraw.app`,
but genome parsing, LOSAT searches, rendering, Session save and load, and
export all run in your browser. The hosted site may collect aggregate page-use
analytics; genome files and generated diagrams are never sent with them. To
work offline, install gbdraw and run `gbdraw gui`.

A saved Session embeds your input files, so protect it like the source data.
[Execution, privacy, and offline
use](./REFERENCE/web-app.md#execution-privacy-and-offline-use) gives the full
details and the browser-performance limits.

### How should I prepare a figure for publication?

Keep an SVG master when the submission workflow accepts vector artwork. Proof
the final file at its placement size. Keep the inputs, software versions,
comparison results, render instructions, and any manual-edit notes.

[Publication and reproducible
handoff](./REFERENCE/output-formats-and-export.md#publication-and-reproducible-handoff)
lists the format choices and what to archive.

## Inputs, layout, and presentation

### Why do my CLI and browser renders differ slightly?

The web app draws with the same Python engine as the CLI, so a difference has
an input cause: a font that exists on one machine but not the other (custom or
installed fonts; see [Fonts](./REFERENCE/output-formats-and-export.md#fonts)),
different gbdraw versions, or settings that differ because the web app starts
from its own defaults (for example legend side and **Arrange in rows**; see
[Browser defaults](./REFERENCE/web-app.md#browser-defaults)).

### Can I use my own font, and will others see it?

Yes. Choose it in **Font** in the web app (a font file, or an installed font in
Chrome and Edge), or use `--font_family` with `--font_file` on the command line.
Exported SVG embeds the characters it uses of the font, so others see it
without installing it, unless the font's license forbids embedding; then it
shows the font only where it is installed. PNG always draws the font. PDF
embeds it too, but draws a font whose license forbids embedding in the nearest
Liberation font. A Session keeps only the font's
name, so whoever opens it needs the same font or gbdraw draws the nearest
Liberation font. See [Fonts](./REFERENCE/output-formats-and-export.md#fonts).

### How do I hide the coordinate scale without hiding the genome axis?

Use `--hide_scale` on the command line or clear **Show Coordinate Scale** in
the web app. Explicit custom track slots can have a separate `ticks` renderer,
and quantitative tracks have their own axes. See [diagram
layout](./REFERENCE/palettes-feature-rules-labels-shapes-and-tracks.md#diagram-layout)
and the [CLI layout
options](./REFERENCE/command-line.md#record-selection-and-layout).

### Can I use a GFF3 file by itself?

No. Pair each GFF3 input with the matching FASTA sequence. [Sequence and
annotation files](./REFERENCE/input-formats-and-tsv-schemas.md#sequence-and-annotation-files)
explains ID matching, coordinates, phase, and translation.

### My labels overlap. What should I do?

Reduce the label size, filter the label set, or change the Circular track
placement. [Feature
presentation](./REFERENCE/palettes-feature-rules-labels-shapes-and-tracks.md#feature-presentation)
covers qualifier priority, whitelist and blacklist behavior, overrides, and
overlap controls.

### Why did gbdraw draw nothing and say my labels do not fit?

Some labels would overlap, or radial labels needed more room than the circular
tracks leave, so gbdraw stopped instead of drawing them. Run the command again
with one of the settings the message lists, for example
`--outer_label_y_radius_offset 1.15` to move outer labels outward,
`--label_font_size 12` to make the labels smaller, or `--labels out` to move
inner labels outside. If only leader lines cross, labels cover the definition
text, or the Circular definition text overlaps a track, the figure is drawn and
the same settings are shown as a warning (for the definition text, for example
`--definition_font_size 10.5` or `--center_reserved_radius 231.99`). A long
organism name in the centre of a Circular map is never wrapped onto more
lines; if it overlaps a track, type a shorter name in **Species**
(`--species`), such as the abbreviation `MjeNMV` for "Marsupenaeus japonicus
endogenous nimavirus". See
[Labels that do not
fit](./REFERENCE/palettes-feature-rules-labels-shapes-and-tracks.md#labels-that-do-not-fit).

### How do I change the color of one specific gene?

Use a specific-color rule that matches the intended feature qualifier. See the
[styling tables](./REFERENCE/input-formats-and-tsv-schemas.md#styling-tables)
and [which feature rule wins](./REFERENCE/palettes-feature-rules-labels-shapes-and-tracks.md#feature-presentation).

### Why does a web color edit create a qualifier rule for some labels and `hash` rules for others?

The editor writes one exact qualifier rule only when that rule covers the
selected scope without ambiguity. Otherwise it keeps feature-specific rules,
because features with identical biological identities cannot be told apart as
single instances after you generate the diagram again. See [preview, search, and
editor](./REFERENCE/web-app.md#preview-search-and-editor) for the exact
conditions and scope controls.

### How do I mark a coordinate range or a group of features?

Use an annotation table and bind its `set_id` to an `annotations` track
slot. See [annotation table
fields](./REFERENCE/input-formats-and-tsv-schemas.md#annotation-table-fields)
and [track placement](./REFERENCE/palettes-feature-rules-labels-shapes-and-tracks.md#tracks-axes-and-annotations).

### Can I use gene names instead of product descriptions for labels?

Yes. Set qualifier priority so that `gene` precedes `product`. See [feature
presentation](./REFERENCE/palettes-feature-rules-labels-shapes-and-tracks.md#feature-presentation)
and the [qualifier-priority
table](./REFERENCE/input-formats-and-tsv-schemas.md#styling-tables).

## Comparisons and sessions

### My comparative diagram has no links. What should I check?

Check the table format, record direction and identifiers, display thresholds,
selected edges, and protein translations. Circular rings also need the
intended reference side. [Result meanings and
limits](./REFERENCE/comparison-programs-thresholds-and-results.md#result-meanings-and-limits)
lists the cases that give an empty result.

### How do I draw several Linear records without comparing them in the web app?

Fresh Linear pages and **Reset Settings** start with **No comparison**. If a
comparison is active, select **No comparison** among the three buttons at the top of **Comparison**.
They apply one choice to every adjacent pair, and the
button that matches the effective plan is pressed. Open **Selected pairs
(N)** to inspect or edit a custom pair plan. Choosing **No comparison** keeps
the pair files and raw-result names you already loaded, inactive, so you can
reuse them later. See the [web
comparison controls](./REFERENCE/web-app.md#comparison-surfaces) and
[selected Linear edges](./REFERENCE/comparison-programs-thresholds-and-results.md#selected-linear-edges).

### Why did gbdraw rerun LOSATP after I loaded a Session?

gbdraw reuses a saved search only when its sequences, direction, program, and
search settings that affect the result still match; otherwise it searches again. See [saved comparison results and
cache reuse](./REFERENCE/session-and-request-compatibility.md#saved-comparison-results-and-cache-reuse).

### Why does a loaded Circular Session say Records not inspected?

**Load Session** shows the saved Result without reading the embedded source
again. Select **Inspect source records**, or **Generate Diagram**, to list the
records and show their rotation rows. A newly uploaded file is inspected at
once. See [Save and Load Sessions](./REFERENCE/web-app.md#save-and-load-sessions)
and [Circular source records](./REFERENCE/web-app.md#circular-source-records-and-one-record-settings).

### Why are controls unavailable while a Session saves or loads?

Save and Load need the Session to stay unchanged, so edits, **Generate Diagram**,
and the other Session button wait until **Saving session…** or **Loading
session…** disappears. Scrolling, preview pan and zoom, and feature search stay
available. If Load reports **Session not loaded**, check that a plain Session
file is at most 200 MiB and that a gzip Session expands to at most 512 MiB. See
[Save and Load Sessions](./REFERENCE/web-app.md#save-and-load-sessions).

### Why does Save Raw LOSAT TSV not contain the internal `h_` IDs?

Generated protein results export stable, readable aliases instead of
handles that exist only inside a Session. Uploaded tables are not rewritten. See [raw results and
cache identity](./REFERENCE/comparison-programs-thresholds-and-results.md#raw-results-and-cache-identity).

### Can pairwise comparison links be curved?

Yes. `curve` bends the same mapped spans that `ribbon` draws as filled
links; it does not change which regions are matched. See [result meanings and
limits](./REFERENCE/comparison-programs-thresholds-and-results.md#result-meanings-and-limits).

## Depth, content, and skew

### What if one record has no Depth TSV for a sample?

Keep an explicit missing entry for that record. gbdraw does not use another
file in its place and does not draw zero coverage. See [tracks, axes, and
annotations](./REFERENCE/palettes-feature-rules-labels-shapes-and-tracks.md#tracks-axes-and-annotations)
and the [CLI depth input](./REFERENCE/command-line.md#tracks-annotations-and-presentation).

### How do I make the GC content track smoother?

Use a larger window and step. Larger values smooth the trace and reduce local
detail. [Tracks, axes, and
annotations](./REFERENCE/palettes-feature-rules-labels-shapes-and-tracks.md#tracks-axes-and-annotations)
gives the relationship.

### Can I plot AT instead of GC?

Yes. Select the `AT` dinucleotide. [Tracks, axes, and
annotations](./REFERENCE/palettes-feature-rules-labels-shapes-and-tracks.md#tracks-axes-and-annotations)
describes content and skew for reversed pairs.

## Output and limitations

### Why does SVG export work but PNG/PDF/EPS/PS export fail?

Command-line and Python conversion to those formats requires CairoSVG and its
runtime libraries. In a conda environment, PDF, EPS, or PS export can also abort
with `free(): invalid pointer` when pandas loads conda-forge `pyarrow` (its `orc`
2.1.2 ships its own zlib). Remove `pyarrow` from that environment
or export SVG; a pip install is not affected. In Python, `import cairocffi` before
`import pandas` or `import gbdraw` also avoids the abort. See [output names, dependencies, and
overwrite](./REFERENCE/output-formats-and-export.md#names-dependencies-and-overwrite).

### Are there known visualization or conversion limitations?

Yes. See the [feature-presentation
limits](./REFERENCE/palettes-feature-rules-labels-shapes-and-tracks.md#feature-presentation)
and [static, interactive, and raster
output](./REFERENCE/output-formats-and-export.md#static-interactive-and-raster-output).

[Documentation home](./DOCS.md) | [Tutorials](./TUTORIALS/README.md) | [Technical documentation](./REFERENCE/README.md) | **FAQ** | [Gallery](./GALLERY.md)
