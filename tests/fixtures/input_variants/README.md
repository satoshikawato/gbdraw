# Input variants

Small synthetic inputs that reproduce quirks of files other tools write, taken from issues of
competing tools. No attachment from an issue tracker is copied. Each file is a few kilobytes,
generated from a random sequence, and carries only the quirk named below.
`tests/test_input_variants.py` reads every file listed here through `gbdraw.read_genbank` or
`gbdraw.read_gff`, then draws it with `gbdraw circular` and `gbdraw linear`. A file must read
as the test expects, or fail with a named `ERROR:` message. It must never read as empty or
crash with a bare exception.

A `.gff3` row below lists its `.fasta` partner under the same name stem, except where noted.
The test checks that every file in this folder appears in the table.

## GenBank and EMBL

| File | Origin | Exercises |
|---|---|---|
| `prokka_locus_tag_only.gb` | Prokka output, many | CDS with `/locus_tag` and `/inference` only, no `/gene`; tRNA and rRNA |
| `rast_seed_dbxref.gb` | RAST output | `/db_xref="SEED:fig\|..."`, `misc_feature` |
| `geneious_locus_spaces.gb` | Geneious export | LOCUS name containing spaces; a `/note` that spans three lines |
| `artemis_colour.gb` | Artemis-saved file | `/colour="2 125 0"` qualifier |
| `mibig_gene_kind.gb` | MIBiG entries | `/gene_kind`, `/gene_functions`, a `protocluster` feature |
| `pgap_partial_locations.gb` | NCBI PGAP | partial locations `<1..300`, `2501..>3000`, `complement(<801..1200)` |
| `pgap_contig_no_sequence.gb` | NCBI PGAP, WGS master | `CONTIG` line and no sequence |
| `geseq_plastome.gb` | GeSeq (OGDRAW family) | `repeat_region` with `/rpt_type=inverted` (IRa, IRb); trans-spliced `rps12` as `join(complement(..),..)` with `/trans_splicing`; intron-containing genes as `join()` |
| `ddbj_flavour.gb` | DDBJ | DDBJ header (REFERENCE, COMMENT), `/protein_id`, `taxon:` db_xref |
| `easyfig72_embl.embl` | Easyfig#72 | EMBL flat file instead of GenBank |
| `bakta_cog_dbxref.gb` | Bakta | repeated `/db_xref` with `COG:COG1136`, `COG:V`, `RefSeq:`, `UniRef:`, `SO:`, `KEGG:`; `ncRNA` |
| `pharokka_phrog.gb` | Pharokka | `/function` and `/phrog` qualifiers |
| `easyfig59_tabs.gb` | Easyfig#59 | tab characters inside a qualifier value and as qualifier-line indent |
| `crlf_line_endings.gb` | Windows editors, many | CRLF line endings |
| `cgviewjs3_one_bp_feature.gb` | CGView.js#3 | 1 bp features on both strands, and a 1 bp `variation` |
| `pseudo_and_pseudogene.gb` | NCBI annotation | `/pseudo` and `/pseudogene="unitary"` on gene and CDS |
| `dfv85_complement_join.gb` | DnaFeaturesViewer#85 | `complement(join(..))` and `join(complement(..),complement(..))` in both part orders |
| `dfv29_origin_spanning_join.gb` | DnaFeaturesViewer#29 | joins across the origin on the plus strand, the minus strand, and as `join(complement(..),complement(..))` |
| `cgview7_multi_record.gbff` | CGView#7 | one `.gbff` with 3 records of 2000, 3500 and 800 bp |
| `zero_features_record.gb` | many | a record with no `FEATURES` table |
| `long_feature_over_10kb.gb` | many | a 10.8 kb CDS in a 12 kb record |
| `lowercase_sequence_n_runs.gb` | many | lower-case sequence with runs of `N` and `n` |

## GFF3

| File | Origin | Exercises |
|---|---|---|
| `ggn35_embedded_fasta.gff3` | gggenomes#35 | `##FASTA` section in the same file (passed as its own FASTA partner; no `.fasta` file) |
| `ggn51_feature_without_id.gff3` | gggenomes#51 | features without `ID=` |
| `ggn67_hash_in_attribute.gff3` | gggenomes#67 | `#` inside attribute values |
| `ggn238_lowercase_cds_type.gff3` | gggenomes#238 | lower-case `cds` type |
| `dfv70_sequence_region_longer.gff3` | DnaFeaturesViewer#70 | `##sequence-region` far longer than the last feature |
| `ggn117_crlf.gff3` | gggenomes#117 | CRLF line endings in GFF3 and FASTA |
| `ggn227_fasta_header_tab.gff3` | gggenomes#227 | FASTA header holding a tab and `key=value` words |
| `ncbi_gene_mrna_cds.gff3` | NCBI RefSeq GFF | gene, mRNA, exon, and a 3-exon CDS sharing one `ID` with phase 0, 1, 2 |
| `prokka_embedded_fasta.gff3` | Prokka GFF | `##FASTA` section; `inference`, `locus_tag`; a score on rRNA (passed as its own FASTA partner) |
| `ggn65_augustus_exon_hierarchy.gff3` | gggenomes#65 | Augustus gene, transcript, exon, intron, start and stop codon hierarchy |
| `ggn93_exonerate_style.gff3` | gggenomes#93 | exonerate `gene`, `cds`, `exon`, `splice5`, `intron`, `splice3`, `similarity`; attributes separated by ` ; ` |
| `score_and_phase_dot.gff3` | many | score `.` and phase `.` on CDS; strand `.` on a `misc_feature` |
| `pgv54_no_strand.gff3` | pyGenomeViz#54 | strand `.` and `?` |
| `clinker67_single_line_fasta.gff3` | clinker#67 | the FASTA sequence on one line, no trailing newline |
| `gff3_fasta_missing.gff3` | many | a GFF3 whose `.fasta` partner does not exist; the CLI must say so |
