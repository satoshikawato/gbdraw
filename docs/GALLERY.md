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

#### Lable whitelist

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


#### Feature qualifier priority
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

[Home](./DOCS.md) | [Installation](./INSTALL.md) | [Quickstart](./QUICKSTART.md) | [Tutorials](./TUTORIALS/TUTORIALS.md) | **Gallery** | [FAQ](./FAQ.md) | [ABOUT](./ABOUT.md)