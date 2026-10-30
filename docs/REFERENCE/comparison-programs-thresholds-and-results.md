[Documentation home](../DOCS.md) | [Tutorials](../TUTORIALS/README.md) | [Technical documentation](README.md) | [Comparison FAQ](../FAQ.md#which-comparison-method-should-i-use)

# Comparison programs, thresholds, and result semantics

## Capability matrix

| Comparison input or mode | Web app | Command line and Python | Result shown |
|---|---|---|---|
| Uploaded BLAST outfmt 6/7 | **Upload BLAST TSV** | Prepared table input | One retained row per Linear link or Circular span |
| LOSATN | Browser nucleotide search | Linear and Circular rings: `--losat losatn` (Python `losat="losatn"`) | Local nucleotide-alignment spans |
| TLOSATX | Browser six-frame translated-nucleotide search | Linear and Circular rings: `--losat tlosatx` (Python `losat="tlosatx"`) | Translated query and subject spans |
| LOSATP Pairwise | Browser protein search | `--losat losatp --losatp_mode pairwise` (Python `losat="losatp", losatp_mode="pairwise"`) | Retained individual protein matches |
| LOSATP **Similarity groups** | Browser protein search | `--losat losatp --losatp_mode similarity_groups` (Python `losatp_mode="similarity_groups"`) | Search-derived group membership and links |
| LOSATP **Collinear blocks** | Browser protein search | `--losat losatp --losatp_mode collinear` (Python `losatp_mode="collinear"`) | Ordered blocks built from compatible anchors |
| Selected mixed edges | Per-edge browser plan | `--comparisons_table` rows with `source` `losat` or `table` | Only the selected record pairs |
| Circular similarity rings | Uploaded or browser-generated results | `--losat losatn\|tlosatx` with `--conservation_sequence`, or prepared `--conservation_blast` | Ordered rings around one reference |

LOSATN compares nucleotide sequence directly. TLOSATX translates both sides
and is useful when coding similarity remains after nucleotide divergence.
LOSATP searches CDS protein sequences, taken from `/translation` or translated
as described in [Input formats](input-formats-and-tsv-schemas.md#sequence-and-annotation-files).
A CDS whose protein is unusable, such as one with an internal stop codon,
cannot contribute protein matches.

The command line and Python API run LOSATN, TLOSATX, and LOSATP for Linear
diagrams, with native LOSAT or the matching NCBI BLAST+ program. They search
the same source-file databases as the web app and write the same raw-cache
keys, so a Session saved from either reuses its LOSAT results in the other.
For Circular similarity rings they run LOSATN or TLOSATX: each comparison
genome (`--conservation_sequence`, one file per ring) is the query and the
subject database is the whole reference file, as in the web app: every record
of the file that has a sequence, uncropped and on the file's strand, also
records that are not drawn (for GFF3 input, every entry of the FASTA file,
also one the GFF3 file does not name). A record without a sequence (a CON
record) is left out. Turning records off, drawing one record of the file, or
a crop does not change ring E-values, and the subject gets the web app's
raw-cache key. With several reference files, each file is searched whole, in
input order; records passed to the Python API in memory are searched as
given. The search reports positions on the whole reference file; gbdraw puts
each row in the search frame of the record it draws (the record after its
crop, on its source strand, 1-based), the frame of every ring table, and leaves
out rows of records that are not drawn or wholly outside a crop. The web app
converts the rows with the same Python code before Generate renders them. The
ring tables of `--losat_output_dir`, `conservation.tsv`, and a Session hold
these rows; the raw LOSAT cache and **Save Raw LOSAT TSV** keep the whole-file
search output. A comparison
FASTA with ID-only headers and 60-column uppercase lines, like the tutorial
files, gets the web app's raw-cache key; other FASTA layouts and GenBank or
DDBJ files give the same rows under a different key. Precomputed rings read
exported rows (`--conservation_blast`).
The web app can run all three LOSAT program families. The web app offers the
LOSATN task `dc-megablast`; a native runtime that does not list that task (the
released LOSAT 0.1.0) stops before searching. Threaded browser search requires cross-origin
isolation.

For Linear browser searches, the three **LOSAT Mode** buttons select
**LOSATN**, **LOSATP**, or **TLOSATX**. When LOSATP is active, the **LOSATP
mode** menu selects **Similarity groups**, **Collinear blocks**, or **Pairwise
matches**. LOSATN, TLOSATX, and LOSATP Pairwise matches can run on selected
LOSAT edges. Similarity groups
and Collinear blocks require an all-adjacent LOSAT plan because they expand
the search across the whole record set. Changing either mode control preserves inactive
program, presentation, and appearance drafts. It does not set **Match style**
as a side effect. **Comparison appearance** remains reachable for every active
Linear comparison, so set **Curve** directly after selecting the intended
LOSATP mode.

Similarity groups always uses all-vs-all search results across the loaded
records and does not read the Collinear **Search between** field. For Collinear
blocks, fresh state and **Reset Settings** use **Adjacent pairs**. A saved
**All records** value is restored from a session. Search scope and display
topology are separate, so an all-record search can still produce links
only between adjacent display rows.

## Filters and direction

Display filters are e-value, bit score, identity, and alignment length. Values
must be finite and non-negative; identity is limited to 0 through 100, and
alignment length is an integer. Protein modes can also limit candidates or
hits per query. Collinear mode adds anchor, unit-gap, diagonal-drift, scope,
and conflict rules. Its unit mode is `auto`, `cds`, or `locus`
(`--collinear_unit_mode`); anchor mode is `rbh` or `all`
(`--collinear_anchor_mode`); and merge orientation is `strand` or `order`
(`--collinear_merge_orientation`). `strand` requires the anchor strands to
follow the block orientation as well as the gene order; `order` reads the gene
order only. A saved `either`, an earlier name that behaved as `strand`, reads
as `strand`. A saved `one_to_one` anchor mode, which was removed, reads as
`rbh`, and Web Load lists the replacement. For each record pair and orientation, gbdraw takes the longest chain
of compatible anchors (most anchors, then the highest summed bit score), removes
it, and repeats; each anchor left over is a singleton block. Two neighbouring
anchors of a chain may have at most `max_unit_gap` units between them, on both
records; at `0` they must be adjacent units.
`max_conflicts` defaults to `1` and limits how many other anchors may lie
strictly inside both order intervals between two neighbouring anchors of a
chain, so it has no effect when `max_unit_gap` is `0`. Those other anchors stay
in the result in their own blocks. `max_diagonal_drift` limits the
difference between the two unit gaps of neighbouring anchors, so at its
default `0` a gap is bridged only when both records skip the same number of
units; the drift can never exceed `max_unit_gap`, so it also has no effect when
`max_unit_gap` is `0`.

Filters apply to display or derived results after raw search rows exist. A raw
TSV can therefore contain more rows than the finished diagram. Record the
program and version, inputs and order, task, thresholds, search scope,
scheduling, and thread count when an exact result must be reproduced.

LOSATP Pairwise search keeps at most one HSP for each query-subject protein
combination in its raw result. The display-stage `max_hits` setting
(`--losatp_max_hits` on the command line) is a separate limit: it
defaults to `5` and keeps the strongest distinct subject proteins for each
query protein after the display thresholds are applied. Similarity-group and
Collinear construction use **Max hits per protein for groups**
(`--losatp_member_max_hits`), a separate limit on threshold-qualified hits: in
each directed record-pair search it keeps, for each query protein, the
strongest distinct subject proteins other than the protein itself. Changing a
derived limit does not change raw search jobs or raw cache keys.

In Web **Settings**, **Max target seqs** sets the raw LOSATP
`--max-target-seqs` limit (`--losatp_max_target_seqs` on the command line).
A blank value leaves the search unbounded. In a within-record search, the
protein's hit to itself counts toward this limit, so a limit of N leaves at
most N - 1 other proteins of the same record for each query protein. Both
**Max target seqs** and **Max hits per protein for groups** are unbounded by
default in every LOSATP mode, on the Web, the command line, and in Python.
Each LOSATP mode remembers its own values when switching modes, including after
a Session round trip; a Session keeps the values it saved.

Collinear blocks builds no similarity groups by default: **Build similarity
groups with within-record searches** on the Web, `--collinear_infer_orthogroups
off` on the command line, `collinear_infer_orthogroups=False` in typed requests,
and `collinearity_infer_orthogroups=False` in `LinearComparisonOptions`. OFF
searches between records only and builds blocks from the selected anchors. ON
adds within-record searches and groups proteins, paralogs included, into
similarity groups before blocks are built. Similarity groups always builds
these groups. A saved Collinear Session or request without the setting, and a
recorded Collinear command line without the flag, keep ON, the behavior they
were made with. **Search between** (`--collinear_search_scope`) and **Max hits
per protein for groups** (`--losatp_member_max_hits`) change the group evidence:
with groups OFF, **Search between** changes nothing, and the member limit
changes the anchors only with the `all` anchor mode, so the Web shows each only
when it can change the result.

The Web app, the CLI, and the Python API treat one source file as one genome
when they run LOSATN, TLOSATX, or LOSATP. Records from two different files are
searched as one job per directed file pair, and the E-value database is the
subject file. A comparison between two records of the same file searches that
file without the query record. A requested within-record search (Similarity
groups, or Collinear with inference ON) searches the record alone. A record
therefore never searches itself unless that search was requested. Records of
one multi-record file therefore get the same E-values and raw cache keys in all
three. Run Info states this scope. Settings shows the job count from the same
plan that Generate runs.

In the Web app, a Similarity group name or description stays with the exact set
of member proteins. When regrouping changes the members, the name is kept as a
saved name that waits for those members, returns when the same group forms
again, and is listed with **Clear** in the Similarity group panel. Session files
store saved names in `webEdits.orthogroupDormantOverrides`; older Sessions have
none.

Prepared rows keep `qseqid`, `sseqid`, `pident`, `length`, `mismatch`,
`gapopen`, `qstart`, `qend`, `sstart`, `send`, `evalue`, and `bitscore` from
the first 12 columns; extra columns are ignored. Query and subject are
directional. In Linear mode, a row that names the other endpoint or another
displayed record is rejected. A version suffix difference is accepted. IDs
that match no displayed record keep the positional pair with a warning (logged
by the CLI, shown beside the Result in the web app). See
[Comparison and numeric tables](input-formats-and-tsv-schemas.md#comparison-and-numeric-tables).
A start greater than its end marks reverse orientation; it does not by itself
mean that a hit crosses a Circular origin.

For a Circular ring, the selected reference side is `query`, `subject`, or
`auto`. The renderer paints retained spans from that side on the displayed
record. It does not infer binned coverage or a wraparound hit from one BLAST
row.

## Selected Linear edges

The Linear **Comparison** command group has **No comparison**, **Run LOSAT**,
and **Upload BLAST TSV** buttons. Each button applies one choice to all
positional adjacent pairs. The button that matches the resolved plan is pressed;
a selected or mixed plan presses none of them and shows a **Custom** badge.

Open **Selected pairs (N)** to change the source or uploaded file for one pair,
omit a pair, or select **Add** and define an explicit non-adjacent pair. The
pair sections are grouped there instead of appearing between record cards.
Each included uploaded edge needs its own active file. An omitted edge
contributes no comparison results and does not reserve a search job. Retained inactive
files and raw-result names do not reactivate an edge.

A selected or mixed plan can use LOSATN, TLOSATX, or LOSATP Pairwise matches.
When **LOSAT Mode** is LOSATP, its **Similarity groups** and **Collinear
blocks** presentation options are disabled and marked **(all adjacent pairs only)**.
Select **Run LOSAT** at the top of the Comparison card, or **Run LOSAT for all
pairs** in its issue list, to change to the compatible all-adjacent topology; opening or closing
**Settings** does not change the plan.

Search scope and display topology are separate. A Similarity group can include
members from several records even when links are drawn only between adjacent
rows. Collinear search can use adjacent pairs or all records while the
finished diagram limits ribbon endpoints for readability.

## Result meanings and limits

- A Pairwise link is one retained search row. It does not establish that two
  complete genes or genomes are equivalent.
- A Similarity group is a connected group derived from retained protein-search
  relationships. It is not a phylogenetic orthogroup. Generated `og_*` IDs are
  stable only for the same inputs, input order, program and runtime version,
  and settings; they are not universal biological identifiers.
- A Collinear block is an ordered run of compatible retained anchors under the
  selected gap, drift, and scope rules. It is not a whole-genome synteny call.
- A Circular ring is a view of retained spans against one reference. It does
  not infer evolutionary conservation for every base inside a span.

Pairwise links use filled `ribbon` geometry by default. `curve` bends the same
mapped spans; it does not change the underlying matches. Collinear color modes can encode
orientation, average identity, or both. Changing an anchor, gap, scope, or
color calculation setting creates a different derived result.

An empty comparison can mean that the search produced no rows, every row
failed a filter, endpoints did not map to displayed records or features, the
wrong Circular reference side was selected, or CDS translations were unusable.
Inspect warnings and the raw search rows before relaxing thresholds.

## Raw results and cache identity

**Save Raw LOSAT TSV** and the Linear `--losat_output_dir` preserve raw search
rows. For Circular rings, the raw rows are the search output on the whole
reference file, and the ring tables are in the search frame of the drawn
records.
Generated protein results replace session-only runtime handles with stable
percent-encoded protein or feature aliases. Export fails if a handle cannot be
resolved safely. An uploaded comparison table is not rewritten.

A reusable cache entry binds sequence content, selected proteins, record and
feature identities, direction, program, and meaningful search arguments.
Changing display labels or filenames does not invalidate the biological
search; changing input content, selected features, direction, program, or
search settings does. Reuse is an optimization. Keep the exported raw results
as the record of what was searched and found.
