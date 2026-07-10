[Home](../DOCS.md) | [Installation](../INSTALL.md) | [Quickstart](../QUICKSTART.md) | [Tutorials](./TUTORIALS.md) | [Recipes](../RECIPES.md) | [CLI Reference](../CLI_Reference.md) | [Gallery](../GALLERY.md) | [FAQ](../FAQ.md) | [About](../ABOUT.md)

[< Back to the Tutorials Index](./TUTORIALS.md)
[< Back to Tutorial 4](./4_Protein_Comparisons.md) | [Go to Tutorial 6 >](./6_Depth_Quantitative_Tracks.md)

# Tutorial 5: Table-Driven CLI Inputs

**Goal:** use TSV manifests when row-coupled inputs are clearer than long order-sensitive option lists.

Relative paths in table files resolve against the table file, not against the shell's current directory.

## 1. Prepare Example GenBank Files

```bash
wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=LC738868.1&rettype=gbwithparts&retmode=text" -O MjeNMV.gb
wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=LC738874.1&rettype=gbwithparts&retmode=text" -O MelaMJNV.gb
wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=LC738870.1&rettype=gbwithparts&retmode=text" -O PemoMJNVA.gb
```

## 2. Linear `--records_table` for GenBank Rows

Create `linear_records.tsv`:

```tsv
gbk	record_label	record_subtitle	record_id	region	reverse_complement	order
MjeNMV.gb	Marsupenaeus japonicus endogenous nimavirus	Ginoza2017	LC738868.1		0	1
MelaMJNV.gb	Melicertus latisulcatus majanivirus	Okinawa2016	LC738874.1	1-160000	0	2
PemoMJNVA.gb	Penaeus monodon majanivirus A	Mikawa2016	LC738870.1	1-160000	1	3
```

Then run:

```bash
gbdraw linear \
  --records_table linear_records.tsv \
  --protein_blastp_mode pairwise \
  --pairwise_match_style curve \
  -o majani_records_table \
  -f svg
```

`--records_table` replaces `--gbk`, `--gff`, and `--fasta`. In linear mode, put per-record labels, subtitles, selectors, crops, and orientation in the table instead of combining `--records_table` with `--record_label`, `--record_subtitle`, `--record_id`, `--region`, or `--reverse_complement`.

## 3. Linear `--records_table` for GFF3 + FASTA Rows

Create a tiny pair of GFF3 + FASTA inputs:

```bash
cat > toy_a.fna <<'EOF'
>toyA
ATGAAACCCGGGTTTAAATTTCCCGGGAAATTTCCCGGGAAATTTCCCGGGAAATTTCCCGGGAAATTTCCCGGGAAATTTCCCTAA
EOF

cat > toy_a.gff3 <<'EOF'
##gff-version 3
toyA	gbdraw	CDS	1	90	.	+	0	ID=toyA_cds1;product=toy protein A
EOF

cat > toy_b.fna <<'EOF'
>toyB
ATGAAACCCGGGTTTAAATTTCCCGGGAAATTTCCCGGGAAATTTCCCGGGAAATTTCCCGGGAAATTTCCCGGGAAATTTCCCTAA
EOF

cat > toy_b.gff3 <<'EOF'
##gff-version 3
toyB	gbdraw	CDS	1	90	.	+	0	ID=toyB_cds1;product=toy protein B
EOF
```

Create `gff_records.tsv`:

```tsv
gff	fasta	record_label	record_id	region	reverse_complement	order
toy_a.gff3	toy_a.fna	Toy A	toyA		0	1
toy_b.gff3	toy_b.fna	Toy B	toyB	1-90	0	2
```

```bash
gbdraw linear \
  --records_table gff_records.tsv \
  -o toy_gff_records_table \
  -f svg
```

A table must use either all GenBank rows or all GFF3/FASTA rows. Do not mix `gbk` with `gff`/`fasta` in the same table.

## 4. Circular Multi-Record Placement

Create `circular_records.tsv`:

```tsv
gbk	record_id	order	row	column
MjeNMV.gb	LC738868.1	1	1	1
MelaMJNV.gb	LC738874.1	2	1	2
PemoMJNVA.gb	LC738870.1	3	2	1
```

```bash
gbdraw circular \
  --records_table circular_records.tsv \
  --multi_record_canvas \
  --multi_record_size_mode auto \
  -o majani_circular_grid \
  -f svg
```

When a records table provides `row` and `column`, do not also use `--multi_record_position`.

## 5. Circular Conservation Rings with `--conservation_table`

Create a minimal BLAST outfmt 6 file and `conservation.tsv`:

```bash
cat > MjeNMV.MelaMJNV.tblastx.out <<'EOF'
LC738868.1	LC738874.1	91.2	1000	88	0	15000	16000	14800	15800	1e-80	300
EOF
```

Create `conservation.tsv`:

```tsv
blast	label	color
MjeNMV.MelaMJNV.tblastx.out	MelaMJNV	#4E79A7
```

Use the table like this:

```bash
gbdraw circular \
  --gbk MjeNMV.gb \
  --conservation_table conservation.tsv \
  --conservation_reference query \
  --identity 70 \
  --alignment_length 100 \
  -o MjeNMV_conservation_table \
  -f svg
```

`--conservation_table` cannot be combined with `--conservation_blast`, `--conservation_labels`, or `--conservation_colors`.

## 6. Circular Track Slots with `--circular_track_table`

Create `circular_tracks.tsv`:

```tsv
id	renderer	side	r	w	params
features	features	axis
gc_content	dinucleotide_content	inside		0.1	nt=GC
gc_skew	dinucleotide_skew	inside		0.1	nt=GC
at_skew	dinucleotide_skew	inside		0.1	nt=AT,positive_color=#deaf6e,negative_color=#7294e3,legend_label=AT skew
ticks	ticks	inside			tick_label_layout=label_in_tick_out
```

```bash
gbdraw circular \
  --gbk MjeNMV.gb \
  --circular_track_table circular_tracks.tsv \
  --track_type middle \
  --window 500 \
  --step 50 \
  -o MjeNMV_track_table \
  -f svg
```

`--circular_track_table` cannot be combined with inline circular track slot options such as `--circular_track_order`, `--circular_track_slot`, or `--circular_track_axis_index`.

[< Back to the Tutorials Index](./TUTORIALS.md)
[< Back to Tutorial 4](./4_Protein_Comparisons.md) | [Go to Tutorial 6 >](./6_Depth_Quantitative_Tracks.md)

[Home](../DOCS.md) | [Installation](../INSTALL.md) | [Quickstart](../QUICKSTART.md) | [Tutorials](./TUTORIALS.md) | [Recipes](../RECIPES.md) | [CLI Reference](../CLI_Reference.md) | [Gallery](../GALLERY.md) | [FAQ](../FAQ.md) | [About](../ABOUT.md)
