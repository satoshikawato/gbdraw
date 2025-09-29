[Home](./DOCS.md) | [Installation](./INSTALL.md) | [Quickstart](./QUICKSTART.md) | [Tutorials](./TUTORIALS/TUTORIALS.md) | **Gallery** | [FAQ](./FAQ.md) | [ABOUT](./ABOUT.md)
# Gallery
This gallery showcases a variety of plots created with `gbdraw` and the commands used to generate them.

### *Escherichia coli* K-12

A standard circular plot with separate strands and custom title formatting.
```bash
gbdraw circular \
  --gbk NC_000913.gb \
  --species "<i>Escherichia coli</i>" \
  --strain "K-12" \
  -f svg \
  --separate_strands \
  -o ecoli_k12
```

###  *Vibrio cholerae* O395
`gbdraw` automatically processes GenBank files with multiple records (e.g., chromosomes) and saves each as a separate file.

```bash
# Download and decompress data
# wget https://.../GCF_000016245.1_ASM1624v1_genomic.gbff.gz
# gunzip GCF_000016245.1_ASM1624v1_genomic.gbff.gz

gbdraw circular \
  --gbk GCF_000016245.1_ASM1624v1_genomic.gbff \
  --species "<i>Vibrio cholerae</i>" \
  --strain "O395" \
  -f svg \
  --separate_strands \
  --track_type middle
```

**Chromosome I**

**Chromosome II**


### Human Herpesvirus 6 (HHV-6) Comparison
A linear plot comparing two virus genomes (HHV-6A and HHV-6B) using BLASTN results.

```bash
# First, run BLASTN with outfmt 7
# blastn -query NC_001664.fasta -subject NC_000898.fasta -outfmt 7 -out hhv6_comparison.blast

gbdraw linear \
  --gbk NC_001664.gb NC_000898.gb \
  -b hhv6_comparison.blast \
  --align_center \
  --separate_strands \
  -o HHV-6 \
  -f svg
```

### Majaniviruses Multi-Genome Comparison

A more complex linear comparison showing sequence similarity across ten different virus genomes. This demonstrates how `gbdraw` can visualize relationships within a viral family.

```bash
# This requires running pairwise TBLASTX between all adjacent genomes
# (e.g., genome1 vs genome2, genome2 vs genome3, etc.)

# Pairwise TBLASTX search 
tblastx -query MjeNMV.fasta -subject MelaMJNV.fasta -outfmt 7 -out MjeNMV.MelaMJNV.tblastx.out 
tblastx -query MelaMJNV.fasta -subject PemoMJNVA.fasta -outfmt 7 -out MelaMJNV.PemoMJNVA.tblastx.out 
tblastx -query PemoMJNVA.fasta -subject PeseMJNV.fasta -outfmt 7 -out PemoMJNVA.PeseMJNV.tblastx.out 
tblastx -query PeseMJNV.fasta -subject PemoMJNVB.fasta -outfmt 7 -out PeseMJNV.PemoMJNVB.tblastx.out 
tblastx -query PemoMJNVB.fasta -subject LvMJNV.fasta -outfmt 7 -out PemoMJNVB.LvMJNV.tblastx.out 
tblastx -query LvMJNV.fasta -subject TrcuMJNV.fasta -outfmt 7 -out LvMJNV.TrcuMJNV.tblastx.out 
tblastx -query TrcuMJNV.fasta -subject MellatMJNV.fasta -outfmt 7 -out TrcuMJNV.MellatMJNV.tblastx.out 
tblastx -query MellatMJNV.fasta -subject MeenMJNV.fasta -outfmt 7 -out MellatMJNV.MeenMJNV.tblastx.out 
tblastx -query MeenMJNV.fasta -subject MejoMJNV.fasta -outfmt 7 -out MeenMJNV.MejoMJNV.tblastx.out 

# gbdraw
gbdraw linear \
--gbk \
./in_gbk/MjeNMV.gb \
./in_gbk/MelaMJNV.gb \
./in_gbk/PemoMJNVA.gb \
./in_gbk/PeseMJNV.gb \
./in_gbk/PemoMJNVB.gb \
./in_gbk/LvMJNV.gb \
./in_gbk/TrcuMJNV.gb \
./in_gbk/MetlamMJNV.gb \
./in_gbk/MeenMJNV.gb \
./in_gbk/MejoMJNV.gb \
-b \
./in_fna/MjeNMV.MelaMJNV.tblastx.out \
./in_fna/MelaMJNV.PemoMJNVA.tblastx.out \
./in_fna/PemoMJNVA.PeseMJNV.tblastx.out \
./in_fna/PeseMJNV.PemoMJNVB.tblastx.out \
./in_fna/PemoMJNVB.LvMJNV.tblastx.out \
./in_fna/LvMJNV.TrcuMJNV.tblastx.out \
./in_fna/TrcuMJNV.MetlamMJNV.tblastx.out \
./in_fna/MetlamMJNV.MeenMJNV.tblastx.out \
./in_fna/MeenMJNV.MejoMJNV.tblastx.out \
-t majani_custom_color_table.tsv \
-d modified_default_colors.tsv \
--block_stroke_width 1 \
--block_stroke_color gray \
--align_center \
--separate_strands \
-o majani -f svg
```
![majaniviruses](https://github.com/satoshikawato/gbdraw/blob/main/examples/majani.svg)

#### <i>Sorangium cellulosum</i> So ce56 (label whitelist)

```bash
gbdraw circular \
--gbk NC_010162.gb \
-f svg \
--palette edelweiss \
--show_labels \
--separate_strands \
-t NC_010162.feature-specific_table.tsv \
--label_whitelist NC_010162.whitelist.tsv \
-t NC_010162.feature-specific_table.tsv \
--species "<i>Sorangium cellulosum</i>" \
--strain "So ce56"
```
##### Lable whitelist example
| feature type | target qualifier | qualifier value regex (Python) |
| ------ | ------- | ------- |
| CDS | old_locus_tag | sce4138 |
| CDS | old_locus_tag | sce4137 |

[NC_010162.whitelist.tsv](https://github.com/satoshikawato/gbdraw/blob/main/examples/NC_010162.whitelist.tsv) (part)
```
CDS	old_locus_tag	sce4138
CDS	old_locus_tag	sce4137
CDS	old_locus_tag	sce4136
CDS	old_locus_tag	sce4135
CDS	old_locus_tag	sce4134
CDS	old_locus_tag	sce4133
CDS	old_locus_tag	sce4132
...
```
[NC_010162.feature-specific_table.tsv](https://github.com/satoshikawato/gbdraw/blob/main/examples/NC_010162.feature-specific_table.tsv) (part)
```
CDS	old_locus_tag	sce4138	#a4d8a7	Chivosazol biosynthesis
CDS	old_locus_tag	sce4137	#a4d8a7	Chivosazol biosynthesis
CDS	old_locus_tag	sce4136	#a4d8a7	Chivosazol biosynthesis
CDS	old_locus_tag	sce4135	#a4d8a7	Chivosazol biosynthesis
CDS	old_locus_tag	sce4134	#a4d8a7	Chivosazol biosynthesis
CDS	old_locus_tag	sce4133	#a4d8a7	Chivosazol biosynthesis
CDS	old_locus_tag	sce4132	#a4d8a7	Chivosazol biosynthesis
...
```
![NC_010162_edelweiss](https://github.com/satoshikawato/gbdraw/blob/main/examples/NC_010162_edelweiss.svg)


#### Human mitochondrial genome (feature qualifier priority)
The following `qualifier_priority.tsv` designated by `--qualifier_priority` opton specifies which qualifier should be used for the label text of a given feature type. Other features remain the same as default:

```modified_default_colors.tsv
CDS gene
```
| feature type | qualifier |
| ------ | ------- |
| CDS | gene |

```bash
gbdraw circular \
--gbk NC_012920.gb \
-f svg --track_type middle \
--species "<i>Homo sapiens</i>" \
--block_stroke_width 2 \
--axis_stroke_width 5 \
--allow_inner_labels \
--show_labels \
--qualifier_priority qualifier_priority.tsv \
-o NC_012920_middle_qualifier_priority_inner_axis5_def28_italic \
--definition_font_size 28
```
![NC_012920_middle_qualifier_priority_inner_axis5_def28_italic](https://github.com/satoshikawato/gbdraw/blob/main/examples/NC_012920_middle_qualifier_priority_inner_axis5_def28_italic.svg)

#### <i>Nicotiana tabacum</i> chloroplast genome

```bash
gbdraw circular \
--gbk NC_001879.gbk \
--separate_strands \
-f svg \
-o NC_001879_color \
-k CDS,rRNA,tRNA,tmRNA,ncRNA,misc_RNA,rep_origin \
-t 2025-09-19_chloroplast.tsv \
--block_stroke_width 1 \
--block_stroke_color black \
--axis_stroke_width 3 \
--line_stroke_width 2 \
--suppress_gc \
--suppress_skew \
-p default \
--track tuckin \
--show_labels \
--allow_inner_labels \
--qualifier_priority qualifier_priority.tsv \
--outer_label_x_radius_offset 0.90 \
--outer_label_y_radius_offset 0.90 \
--inner_label_x_radius_offset 0.975 \
--inner_label_y_radius_offset 0.975 \
--species "<i>Nicotiana tabacum</i>" \
--definition_font_size 28 \
--legend upper_left
```

![NC_001879_color.svg](https://github.com/satoshikawato/gbdraw/blob/main/examples/NC_001879_color.svg)

#### Lambda phage

```bash
gbdraw linear \
--gbk NC_001416.gb \
-o NC_001416 \
-f svg \
--show_labels \
--separate_strands \
--legend left \
-d cds_white.tsv \
-t lambda_specific_table.tsv \
--block_stroke_width 2 \
--axis_stroke_width 5 \
--definition_font_size 24 
```
![NC_001416.svg](https://github.com/satoshikawato/gbdraw/blob/main/examples/NC_001416.svg)


#### <i>Ca.</i> Sukunaarchaeum mirabile M16-5
```bash
gbdraw circular \
--gbk M16-5.gb \
--separate_strands \
-f svg \
-o M16-5_fugaku \
--block_stroke_width 1 \
--axis_stroke_width 1 \
-p fugaku \
--track middle \
--species "<i>Ca.</i> Sukunaarchaeum mirabile" \
--definition_font_size 22 \
--legend upper_right
```
![M16-5_fugaku.svg](https://github.com/satoshikawato/gbdraw/blob/main/examples/M16-5_fugaku.svg)


#### <i>Pandoravirus_salinus</i>
```bash
gbdraw circular \
--gbk Pandoravirus_salinus.gb \
--separate_strands \
-f svg \
-o Pandoravirus_salinus_forest \
-p forest \
--track tuckin \
--species "<i>Pandoravirus salinus</i>" \
--definition_font_size 22 \
--legend upper_right
```
![Pandoravirus_salinus_forest.svg](https://github.com/satoshikawato/gbdraw/blob/main/examples/Pandoravirus_salinus_forest.svg)

#### Candidatus Pelagibacter ubique HTCC1062
```bash
gbdraw circular \
--gbk NC_007205.gb \
-f svg \
--separate_strands \
--species "<i>Ca. </i> Pelagibacter ubique" \
--strain "HTCC1062" \
--legend none \
--palette oceanic_voyage \
-o NC_007205_oceanic_voyage
```
![NC_007205_oceanic_voyage.svg](https://github.com/satoshikawato/gbdraw/blob/main/examples/NC_007205_oceanic_voyage.svg)

#### <i>Prochlorococcus marinus</i> CCMP1375
```bash
gbdraw circular \
--gbk NC_005042.gb \
-f svg \
--separate_strands \
--species "<i>Prochlorococcus marinus</i>" \
--strain "CCMP1375" \
--legend none \
--palette pine_reflection \
-o NC_005042_pine_reflection
```
![NC_005042_pine_reflection.svg](https://github.com/satoshikawato/gbdraw/blob/main/examples/NC_005042_pine_reflection.svg)

#### <i>Flavobacterium columnare</i> ATCC 49512
```bash
gbdraw circular \
--gbk NC_016510.gb \
-f svg \
--separate_strands \
--species "<i>Flavobacterium columnare</i>" \
--strain "ATCC 49512" \
--legend none \
--palette mint \
-o NC_016510_mint
```
![NC_016510_mint.svg](https://github.com/satoshikawato/gbdraw/blob/main/examples/NC_016510_mint.svg)

#### <i>Thermus aquaticus</i> Y51MC23
```bash
gbdraw circular \
--gbk NZ_CP010822.gb \
-f svg \
--separate_strands \
--species "<i>Thermus aquaticus</i>" \
--strain "Y51MC23" \
--legend none \
--palette orange \
-o NZ_CP010822_orange
```
![NZ_CP010822_orange.svg](https://github.com/satoshikawato/gbdraw/blob/main/examples/NZ_CP010822_orange.svg)

#### <i>Helicobacter pylori</i> J99
```bash
gbdraw circular \
--gbk NC_000921.gb \
-f svg \
--separate_strands \
--species "<i>Helicobacter pylori</i>" \
--strain "J99" \
--legend none \
--palette spring \
-o NC_000921_spring
```
![NC_000921_spring.svg](https://github.com/satoshikawato/gbdraw/blob/main/examples/NC_000921_spring.svg)


#### <i>Mycobacterium tuberculosis</i> H37Rv
```bash
gbdraw circular \
--gbk NC_000962.gb \
-f svg \
--separate_strands \
--species "<i>Mycobacterium tuberculosis</i>" \
--strain "H37Rv" \
--legend none \
--palette psyche \
-o NC_000962_psyche
```
![NC_000962_psyche.svg](https://github.com/satoshikawato/gbdraw/blob/main/examples/NC_000962_psyche.svg)


[Home](./DOCS.md) | [Installation](./INSTALL.md) | [Quickstart](./QUICKSTART.md) | [Tutorials](./TUTORIALS/TUTORIALS.md) | **Gallery** | [FAQ](./FAQ.md) | [ABOUT](./ABOUT.md)