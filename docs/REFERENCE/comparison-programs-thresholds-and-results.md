[Documentation home](../DOCS.md) | [Tutorials](../TUTORIALS/README.md) | [Technical documentation](README.md) | [Comparison FAQ](../FAQ.md#which-comparison-method-should-i-use)

# Comparison programs, thresholds, and result semantics

## Capability matrix

| Evidence or mode | Web app | Command line and Python | Result shown |
|---|---|---|---|
| Uploaded BLAST outfmt 6/7 | **Upload BLAST TSV** | Prepared table input | One retained row per Linear link or Circular span |
| LOSATN | Browser nucleotide search | Read its exported table as prepared input | Local nucleotide-alignment spans |
| TLOSATX | Browser six-frame translated-nucleotide search | Read its exported table as prepared input | Translated query and subject spans |
| LOSATP Pairwise | Browser protein search | `pairwise` protein-search mode | Retained individual protein matches |
| LOSATP **Similarity groups** | Browser protein search | `orthogroup` compatibility token | Search-derived group membership and links |
| LOSATP **Collinear blocks** | Browser protein search | `collinear` protein-search mode | Ordered blocks built from compatible anchors |
| Selected mixed edges | Per-edge browser plan | Explicit prepared-table endpoints | Only the selected record pairs |
| Circular similarity rings | Uploaded or browser-generated evidence | Prepared conservation inputs | Ordered evidence tracks around one reference |

LOSATN compares nucleotide sequence directly. TLOSATX translates both sides
and is useful when coding similarity remains after nucleotide divergence.
LOSATP searches annotated CDS translations. Unusable or missing CDS
translations cannot contribute protein matches.

The command line can run LOSATP or a compatible BLASTP runtime. It does not run
LOSATN or TLOSATX; use `--blast`, `--comparisons_table`, or
`--conservation_blast` for their exported rows. The web app can run all three
LOSAT program families. Threaded browser search requires cross-origin
isolation.

For Linear browser searches, the three **LOSAT Mode** buttons select
**LOSATN**, **LOSATP**, or **TLOSATX**. When LOSATP is active, the **LOSATP
mode** menu selects **Similarity groups**, **Collinear blocks**, or **Pairwise
matches**. LOSATN, TLOSATX, and LOSATP Pairwise matches can run on selected
LOSAT edges. Similarity groups
and Collinear blocks require an all-adjacent LOSAT plan because they expand
evidence across the record set. Changing either mode control preserves inactive
program, presentation, and appearance drafts. It does not set **Match style**
as a side effect. When a Similarity-group or Collinear recipe requires
**Curve** from fresh state, select LOSATP Pairwise matches, set that style,
then select the final LOSATP mode.

Similarity groups always uses all-vs-all search evidence across the loaded
records and does not read the Collinear evidence-scope field. For Collinear
blocks, fresh state and **Reset Settings** use **Adjacent pairs**. A saved
**All records** value is restored from a session. Search evidence and display
topology are separate, so all-record search evidence can still produce links
only between adjacent display rows.

## Filters and direction

Display filters are e-value, bit score, identity, and alignment length. Values
must be finite and non-negative; identity is limited to 0 through 100, and
alignment length is an integer. Protein modes can also limit candidates or
hits per query. Collinear mode adds anchor, unit-gap, diagonal-drift, scope,
and conflict rules. Its unit mode is `auto`, `cds`, or `locus`; anchor mode is
`all`, `one_to_one`, or `rbh`; and merge orientation is `strand`, `order`, or
`either`. `max_conflicts` defaults to `1` and limits how many retained singleton
anchors may lie inside both order intervals when two compatible clusters
merge. Those singleton anchors remain in the result.

Filters apply to display or derived results after raw search rows exist. A raw
TSV can therefore contain more rows than the finished diagram. Record the
program and version, inputs and order, task, thresholds, search scope,
scheduling, and thread count when an exact result must be reproduced.

LOSATP Pairwise search keeps at most one HSP for each query-subject protein
combination in its raw result. The display-stage `max_hits` setting
(`--protein_blastp_max_hits` on the command line) is a separate limit: it
defaults to `5` and retains the strongest distinct subject proteins for each
query protein after the display thresholds are applied. Similarity-group and
Collinear construction use a separate member-hits limit, also `5` by default,
for the directional rows considered during grouping. Changing either derived
limit does not change raw search jobs or raw cache keys.

Prepared rows retain `qseqid`, `sseqid`, `pident`, `length`, `mismatch`,
`gapopen`, `qstart`, `qend`, `sstart`, `send`, `evalue`, and `bitscore`.
Query and subject are directional. They must map to the intended displayed
records, including version suffixes. A start greater than its end marks reverse
orientation; it does not by itself mean that a hit crosses a Circular origin.

For a Circular ring, the selected reference side is `query`, `subject`, or
`auto`. The renderer paints retained spans from that side on the displayed
record. It does not infer binned coverage or a wraparound hit from one BLAST
row.

## Selected Linear edges

The Linear **Comparison** command group has **No comparison**, **Run LOSAT**,
and **Upload BLAST TSV** buttons. Each button applies one choice to all
positional adjacent pairs. These are bulk commands rather than a complete list
of possible current states. The separate **Current:** status reports the
resolved state; a selected or mixed plan appears as **Current: Selected pairs
(N; ...)** with a **Custom** badge.

Open **Selected pairs (N)** to change the source or uploaded file for one pair,
omit a pair, or select **Add** and define an explicit non-adjacent pair. The
pair sections are grouped there instead of appearing between record cards.
Each included uploaded edge needs its own active file. An omitted edge
contributes no evidence and does not reserve a search job. Retained inactive
files and raw-result names do not reactivate an edge.

A selected or mixed plan can use LOSATN, TLOSATX, or LOSATP Pairwise matches.
When **LOSAT Mode** is LOSATP, its **Similarity groups** and **Collinear
blocks** presentation options are disabled. Select **Use all adjacent LOSAT**
to change to the compatible all-adjacent topology; opening or closing
**Settings** does not change the plan.

Search scope and display topology are separate. A Similarity group can include
members from several records even when links are drawn only between adjacent
rows. Collinear search can use adjacent or all-record evidence while the
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
mapped spans; it does not change the evidence. Collinear color modes can encode
orientation, average identity, or both. Changing an anchor, gap, scope, or
color calculation setting creates a different derived result.

An empty comparison can mean that the search produced no rows, every row
failed a filter, endpoints did not map to displayed records or features, the
wrong Circular reference side was selected, or CDS translations were unusable.
Inspect warnings and raw evidence before relaxing thresholds.

## Raw results and cache identity

**Save Raw LOSAT TSV** and `--protein_blastp_output` preserve raw search rows.
Generated protein results replace session-only runtime handles with stable
percent-encoded protein or feature aliases. Export fails if a handle cannot be
resolved safely. An uploaded comparison table is not rewritten.

A reusable cache entry binds sequence content, selected proteins, record and
feature identities, direction, program, and meaningful search arguments.
Changing display labels or filenames does not invalidate the biological
search; changing input content, selected features, direction, program, or
search settings does. Reuse is an optimization. Keep exported raw results as
the evidence record.
