[Documentation home](../DOCS.md) | [Tutorials](../TUTORIALS/README.md) | [Input reference](input-formats-and-tsv-schemas.md)

# Tutorial fixture provenance

The machine-readable authority is [`gbdraw/web/tutorial-data/manifest.json`](../../gbdraw/web/tutorial-data/manifest.json). Its schema version 2 records each file's fixture ID, package tier, byte size, SHA-256 checksum, biological identity, source, derivation, intended scenarios, and expected semantics.

`retrievedOn` is the date on which the source was acquired. Some pre-policy fixtures have no recoverable historical acquisition date; they retain `retrievedOn: null` and explicitly declare `retrievalDateStatus: unknown-legacy`. This marker does not estimate or silently substitute a date. Their `repositoryAddedOn` value, repository commit when known, current byte size, and verified SHA-256 preserve the auditable facts that are available. Files added after 2026-08-04 must record an ISO retrieval date and cannot use the legacy exception.

Core fixtures are small public biological records shared by the hosted GUI, local GUI, CLI recipes, Python recipes, and documentation tests. A fixture has one canonical copy under `gbdraw/web/tutorial-data/`; public procedures do not link to `tests/test_inputs`.

Derived files must preserve biological identity and document their source checksum, tool version, arguments, and output statistics. A naturally single sequence is never divided into artificial contigs for a Tutorial. GFF3 + FASTA examples use a complete biological sequence, and multi-record examples use genuine biological records or replicons.

Circular multi-record examples use four complete mitochondrial RefSeq records from
*Homo sapiens*, *Danio rerio*, *Drosophila melanogaster*, and *Caenorhabditis elegans*.
Circular conservation examples use the same complete records. Three frozen
TLOSATX tables compare zebrafish, fruit fly, and nematode mtDNA as queries with
human mtDNA as the subject reference. The reproducible builder pins the native
LOSAT 0.1.0 executable, mitochondrial genetic codes, one-thread scheduling,
two byte-identical runs, and the displayed thresholds.
Track examples use the complete `NC_001879.2` tobacco plastome and complete `AP027133.1`
genome. Only the AP027133 depth measurements are reduced: the committed builder records
1 kbp arithmetic means without altering or cropping the matching genome record.

The frozen nucleotide-comparison fixture uses the complete Lambda
[`NC_001416.1`](https://www.ncbi.nlm.nih.gov/nuccore/NC_001416.1) and DE3
[`NC_042057.1`](https://www.ncbi.nlm.nih.gov/nuccore/NC_042057.1) records. Its
LOSATN and TLOSATX tables are reproducible derivatives of those two whole
records. The fixture manifest binds both inputs and both result tables by
SHA-256 checksum and records the pinned program, scheduling, arguments, row
counts, coordinate coverage, and orientations.
Lambda and DE3 are naturally linear genomes, so public diagrams use this
evidence only in Linear layouts. Circular comparison figures use the complete
metazoan mitochondrial records described above.

The protein-comparison fixture keeps the five native BGC database records
`BGC0000708`, `BGC0000709`, `BGC0000711`, `BGC0000712`, and `BGC0000713`
unchanged. They are linear BGC-region records, not complete chromosomes, and
the documentation never renders them as circular genomes. LOSATP connects only
adjacent records in their displayed order; Similarity groups and Collinear
blocks are documented as separate result models.

Package tiers are:

- `core`: included in the beginner package, with a total uncompressed budget of 1 MiB;
- `extended`: deterministic material with up to 5 MiB additional budget;
- `gallery-nightly`: large data kept outside the beginner bundle and exercised separately.

Tests recalculate checksums, parse biological records, verify accession/version, sequence length, feature/CDS counts and strands, validate derived-result schemas, and reject duplicate canonical copies. The manifest is the source for exact values; this page does not duplicate numbers that would drift.
