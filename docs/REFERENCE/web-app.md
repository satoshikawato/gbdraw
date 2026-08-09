[Documentation home](../DOCS.md) | [Web app Tutorials](../TUTORIALS/GUI/README.md) | [Technical documentation](README.md) | [FAQ](../FAQ.md) | [Gallery](../GALLERY.md)

# Web app

The hosted application at [gbdraw.app](https://gbdraw.app/) and local
`gbdraw gui` command expose the same single-page interface.

## Execution, privacy, and offline use

The hosted site downloads the application and its packaged runtimes from
`gbdraw.app`. After those assets load, the browser parses genome files, runs
LOSAT searches, renders diagrams, manages sessions, and exports files. These
operations do not require an application-server upload. The hosted site uses
Google Analytics 4 for aggregate page-usage metrics, but uploaded genome files
and generated diagrams are not analytics payloads.

`gbdraw gui` serves the application and runtimes from the installed Python
package. Its normal drawing workflow works without an internet connection once
gbdraw is installed. Browser security policy still applies. Threaded LOSAT
also requires a cross-origin-isolated page. Select **Serial** under
**Execution** when that capability is unavailable. Where threaded execution is
available, **Safe** is the conservative **Total threads** choice.

A saved session embeds input resources and should be protected like those
inputs. Browser extensions, proxies, the operating system, and downloaded
files remain outside gbdraw's privacy boundary. Runtime cost depends on record
length, feature density, record and comparison counts, search mode and
thresholds, memory, and CPU. There is no universal maximum genome size.
For a reproducible comparison, begin with **Serial** and one worker, then raise
concurrency only after measuring the same inputs. Large assemblies, dense
all-pairs searches, and large depth tables are better suited to a controlled
local environment or prepared command-line evidence.

## Modes and inputs

| Mode | Primary browser input | Result |
|---|---|---|
| Circular | One GenBank/GBFF container or one matched GFF3 + FASTA pair | One Circular result, separate results, or one multi-record canvas |
| Linear | Ordered GenBank rows or matched GFF3 + FASTA rows | One Linear result with an independent comparison plan |

A Circular GenBank upload uses **GenBank/DDBJ File**. Each Linear record card
uses **GenBank File**; **Add Seq** has the accessible name **Add sequence**.
A GenBank file may contain several biological records. GFF3 input requires the
matching FASTA sequence and exact sequence-ID agreement. See [Input formats and
TSV schemas](input-formats-and-tsv-schemas.md) for the file contract.

## Main workflow

1. Select **Circular** or **Linear**.
2. Add the biological inputs and resolve any validation message.
3. Set layout, tracks, comparisons, labels, and output values.
4. Select **Generate Diagram**.
5. Inspect **Result Preview**, then export a file or select **Save Session**.

Settings are committed to a result only after generation succeeds. The form
may contain newer draft values than the visible result, so a session records
the committed request together with supported editable state.

## Circular multi-record canvas

Turn on **Multi-Record Canvas** to place every selected record from one
multi-record Circular input on a shared canvas. The browser's Circular uploader
accepts one input container; use the command line or Python API when the records
must remain in separate source files. Use Circular layout only for complete
records whose biological topology is circular.

To combine separate files for browser upload, concatenate the complete GenBank
flat files in the intended record order without editing their contents. Every
record must retain its terminating `//`. If an expected accession is absent
from **Record Order**, stop and repair the container before generating.

```bash
cat record-a.gb record-b.gb > combined.gb
```

**Record Size Mode** accepts **Auto (Default)**, **Linear**, and **Equal**.
**Equal** gives every record the same radius, **Linear** scales radius with
record length, and **Auto (Default)** applies bounded automatic scaling.
**Min Radius Ratio** accepts 0.01 through 1. **Column Gap Ratio** defaults to
0.10 and **Row Gap Ratio** defaults to 0.05; both add clearance between visible
record bounds.

**Record Order** controls row assignment and left-to-right order. Every loaded
record appears exactly once. A top or bottom plot title can use **Keep Full
Definition with Plot Title** to keep each record definition inside its own
circle.

## Record selection and layout

Each Linear input card owns its record selector, inclusive **Start** and
**End** coordinates, **Reverse complement** state, definition, and row
placement. A region changes the displayed interval, not the source file.
Reverse complementation changes displayed coordinates, feature orientation,
and comparison endpoint mapping without rewriting the input.

Turn on **Arrange in rows** to assign records to rows. Record-card order is the
left-to-right order within a row, and **Record gap (px)** separates records in
that row. Records that share a row use one bp-per-pixel scale. Row placement is
independent of the comparison plan: **No comparison** draws the records without
links.

**Show Coordinate Scale** controls coordinate ticks and labels while retaining
the record axes. **Ruler on Axis** uses each record axis as its ruler only when
the scale is visible, **Ruler (Ticks)** is selected, and **Track Layout** is
**Above** or **Below**. Titles, record definitions, and legends do not change
biological coordinates.

## Browser defaults

| Setting | Circular | Linear |
|---|---|---|
| Separate strands | On | On |
| Legend | Left | Bottom |
| Feature placement | Tuckin preset | Features on axis |
| Comparison | No rings until configured | No comparison until a plan is selected |

Loading a session restores its saved values instead of applying these
defaults. Turning off **Use custom stack** preserves its draft slots. **Reset
Settings** rebuilds settings from the browser defaults while retaining uploaded
files and the current generated result.

The custom-stack editor cannot add an `annotations` slot until an annotation
table has been imported and an annotation set is available.

## Comparison surfaces

The browser runs LOSATN, TLOSATX, and LOSATP. LOSATP has **Pairwise**,
**Similarity groups**, and **Collinear blocks** modes. **Upload BLAST TSV**
uses prepared BLAST-compatible rows for Linear links or Circular rings.

For Linear input, **Apply to all adjacent gaps** sets **No comparison**, **Run
LOSAT**, or **Upload BLAST TSV** across positional neighbors. Each comparison
boundary can override that bulk choice. **Advanced pair setup** adds selected
record pairs, including non-adjacent edges. An uploaded edge participates only
when it has an active file; an omitted edge draws no link and starts no search.

**LOSAT Mode** selects **LOSATN**, **LOSATP**, or **TLOSATX**. **Execution**
selects **Auto (...)**, **Serial**, or **Threaded**. Threaded runs allocate
**Total threads** through **Safe (...)**, **Available (...)**, or an explicit
number, then expose **Parallel runs** and **Threads per run**. LOSATN **Task**
is `megablast`, `blastn`, or `dc-megablast`.

TLOSATX translates each sequence with its selected genetic code. In Linear
mode, each card's **Gencode (this entry)** control, with accessible name
**TLOSATX gencode for sequence N**, supplies the code for that endpoint. In
Circular mode, **Reference gencode** applies to the displayed subject. Each
comparison-FASTA row's visible **Subject gencode** control applies to that
comparison sequence, even though the search passes the sequence as its query.

LOSATP **blastp mode** selects **Pairwise**, **Similarity groups**, or
**Collinear blocks**. Pairwise presentation uses **Pairwise Match Style**
(**Ribbon** or **Curve**) and **Pairwise Match Height**. Collinear settings
include **Max unit gap**, **Min block genes**, **Diagonal drift**, **Merge
conflicts**, **Evidence scope** (**Adjacent pairs** or **All records**), and
**Color mode** (**Average identity**, **Orientation**, or **Orientation +
identity**). The common filters are **Bitscore**, **E-value**, **Minimum
Identity**, and **Minimum Length**; scientific meanings and limits belong to
the comparison reference below.

Circular **Pairwise Comparisons** selects **Run LOSAT** or **Upload BLAST**.
Uploaded evidence uses **BLAST outfmt 6/7 files** and **Reference side**
(**Auto (...)**, **Query**, or **Subject**). A browser-generated Circular
comparison uses the displayed Circular record as the search subject and each
**Comparison FASTA** as a query. **Ring Width** and **Ring Gap** control the
ordered evidence tracks. **Save Raw LOSAT TSV** exports generated search rows.

See [Comparison programs, thresholds, and result
semantics](comparison-programs-thresholds-and-results.md) for search boundaries,
filters, direction, and interpretation.

## Preview, search, and editor

**Result Preview** exposes records, feature and quantitative tracks,
comparisons, definitions, ticks, legends, and annotations as semantic SVG
objects. Feature search can target **All**, **Label**, **Feature type**,
**Record ID**, **Location**, **Strand**, or **Similarity group**. When rich
feature popups are enabled, the **Field** menu also includes **Qualifier key**,
**Qualifier value**, **Nucleotide**, and **Amino acid**. Search may use a
literal value or **Regex**, and the previous and next controls move through
rendered matches.

A normal feature click opens its identity, location, strand, qualifiers, and
available sequence actions. Match popups report mapped endpoints and evidence;
Similarity-group and Collinear popups add member or anchor context. Sequence
downloads are available only when the required source sequence and metadata
are present.

Ctrl-click selects features for bulk color, legend-caption, visibility, and
stroke edits. Use the visible **Apply** action to make an edit part of the
editor state. **Apply to all label** and **Apply to all source label** become
one anchored qualifier rule only when the selected features share one feature
type, qualifier, and value and that rule matches exactly the intended loaded
features. Otherwise the editor keeps one exact `hash` rule per biological
feature. Identical duplicate records can share the same hash, so a regenerated
diagram cannot preserve a one-instance-only rule for indistinguishable
duplicates.

**Layout edit** moves supported legends and other layout objects without
changing record geometry. **Undo** and **Redo** traverse supported form and
editor changes. **Reset Settings** is broader than undo and requires
confirmation. Regenerate after rule-backed edits when the exported figure must
match the current form state.

The export actions and session handoff rules are documented in [Output formats
and export](output-formats-and-export.md) and [Session and request
compatibility](session-and-request-compatibility.md).

## Accessibility

Primary controls have stable accessible names: **Circular**, **Linear**,
**GenBank/DDBJ File**, **GenBank File**, **Add sequence**, **Output Prefix**,
**Species**, **Track Preset**, **Separate Strands**, **Hide GC Content**,
**Hide GC Skew**, **Label Mode**, **Legend Position**, **Generate Diagram**,
**Result Preview**, and **SVG**. The visible **Show Coordinate Scale** control
has the mode-qualified accessible name **Show Coordinate Scale (Circular)** or
**Show Coordinate Scale (Linear)**. The visible annotation **Labels** toggle
uses **Show annotation labels**. Mode buttons expose pressed state, file
controls are labelled, and the preview is a named region.
