[Home](./DOCS.md) | [Tutorials](./TUTORIALS/TUTORIALS.md) | [Workflow guide](./WORKFLOW_GUIDE.md) | **GFF3 + FASTA** | [CLI Reference](./CLI_Reference.md)

# Draw an annotated GFF3 + FASTA assembly

Pair a GFF3 annotation with its FASTA sequence while preserving record identity, strand, phase, and CDS translation behavior. The examples use a small two-contig fixture extracted from the RefSeq lambda phage record `NC_001416.1`. The fixture retains gene/CDS annotations on both strands and is stored under [`examples/gff3_lambda`](../examples/gff3_lambda/). NCBI sequence records are public-domain data; the extraction is included for documentation and tests.

## 1. Check record IDs

The first GFF3 column and each FASTA header must use the same record ID:

```text
lambda_left  RefSeq-extract  CDS  ...
>lambda_left source=NC_001416.1:1-18000
```

This fixture contains `lambda_left` and `lambda_right` in both files. A missing or misspelled FASTA ID leaves the annotation without its sequence.

## 2. Generate a multi-contig diagram

```bash
gbdraw linear \
  --gff examples/gff3_lambda/lambda_two_contigs.gff3 \
  --fasta examples/gff3_lambda/lambda_two_contigs.fna \
  --features CDS,gene \
  --show_labels first \
  --separate_strands \
  -o lambda_gff3_two_contigs \
  -f svg
```

The two GFF3/FASTA records become two displayed linear records. For several independent GFF3 + FASTA pairs, pass matching path lists in the same order or use the [GFF3 columns in `--records_table`](./TUTORIALS/5_Table_Driven_Inputs.md#3-linear---records_table-for-gff3--fasta-rows).

## 3. Preserve CDS semantics

- Column 7 is `+` or `-`; it controls feature direction and translation strand.
- Column 8 is CDS phase `0`, `1`, or `2`, not a generic frame label.
- `ID` values should be unique within the file, and `Parent` should connect a CDS to its gene or transcript when the annotation model uses that relationship.
- A `translation` attribute makes the protein sequence explicit. Without it, gbdraw may translate a valid CDS from the FASTA sequence, strand, phase/codon start, and genetic code, but incomplete or biologically exceptional CDS features need careful validation.

Protein-search modes require usable CDS translations. Check warnings and feature counts before interpreting missing similarity links.

## 4. Common validation failures

| Symptom | Check | Fix |
|---|---|---|
| Record has no sequence | GFF3 column 1 differs from FASTA ID | Rename one side so IDs match exactly. |
| CDS points the wrong way | Strand is `.` or incorrect | Export strand-aware CDS annotations from the annotation tool. |
| Translation is shifted | CDS phase is missing or wrong | Preserve the caller's phase and confirm it against the encoded protein. |
| Features disappear | Requested `--features` omits their type | Include the exact GFF3 feature type or adjust the visibility table. |
| Protein comparison is empty | CDS has no valid translation | Supply translation attributes or correct sequence/phase/genetic-code metadata. |
| Contigs are merged unexpectedly | IDs are reused | Give each contig a unique, stable ID in both files. |

For production data, keep the annotation tool/version, source assembly accession, GFF3 checksum, FASTA checksum, and any normalization commands beside the figure recipe.

[Home](./DOCS.md) | [Tutorials](./TUTORIALS/TUTORIALS.md) | [Workflow guide](./WORKFLOW_GUIDE.md) | **GFF3 + FASTA**
