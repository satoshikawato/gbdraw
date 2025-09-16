[Home](./README.md) | [Installation](./INSTALL.md) | [Quickstart](./QUICKSTART.md) | [Tutorials](./TUTORIALS/) | **Gallery**
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

gbdraw linear \
  --gbk ./in_gbk/MjeNMV.gb ./in_gbk/MelaMJNV.gb ... \
  -b ./in_fna/MjeNMV.MelaMJNV.tblastx.out ./in_fna/MelaMJNV.PemoMJNVA.tblastx.out ... \
  -t majani_custom_color_table.tsv \
  -d modified_default_colors.tsv \
  --align_center \
  --separate_strands \
  -o majani \
  -f svg
```