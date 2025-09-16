[Home](../README.md) | [Installation](./INSTALL.md) | [Usage](./USAGE.md) | [Gallery](./GALLERY.md) | [RECIPES](./RECIPES.md) | [Advanced](./ADVANCED.md) | [FAQ](./FAQ.md)
#### <i>Haemophilus influenzae</i>
`gbdraw` automatically identifies and displays the organism and strain name from the sequence record. However, these names are not italicized by default. For example:
```bash
wget -c https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/931/575/GCF_000931575.1_ASM93157v1/GCF_000931575.1_ASM93157v1_genomic.gbff.gz # download genome
gunzip GCF_000931575.1_ASM93157v1_genomic.gbff.gz # decompress .gz file
gbdraw circular --gbk GCF_000931575.1_ASM93157v1_genomic.gbff -o Haemophilus_influenzae -f svg
```
![hinfluenzae](https://github.com/satoshikawato/gbdraw/blob/main/examples/Haemophilus_influenzae.svg)
#### <i>Escherichia coli</i> K-12
To italicize a portion of the organism name, you can use the <i></i> tags in the --species and --strain parameters. This will format the specified text in italics. The following command will render the species name "_Escherichia coli_" in italics, while keeping "K-12" in standard text (the organim name will be overridden):
```bash
gbdraw circular --gbk NC_000913.gb --species "<i>Escherichia coli</i>" --strain "K-12" -f svg --separate_strands
```
![ecoli](https://github.com/satoshikawato/gbdraw/blob/main/examples/NC_000913.svg)

#### <i>Vibrio cholerae</i> Strain O395 (GCF_000016245.1)
For GenBank files containing multiple entries, `gbdraw` saves each entry as a separate file. Here's how to do it for the _Vibrio cholerae_ strain O395 genome, which has two chromosomes:
```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/016/245/GCF_000016245.1_ASM1624v1/GCF_000016245.1_ASM1624v1_genomic.gbff.gz # download genome
gunzip GCF_000016245.1_ASM1624v1_genomic.gbff.gz # decomperss .gz file
gbdraw circular --gbk GCF_000016245.1_ASM1624v1_genomic.gbff --species "<i>Vibrio cholerae</i>" --strain "O395" -f svg --separate_strands --track_type middle # Draw genome; results in "NC_009457.svg" for Chromosome I and "NC_009456.svg" for Chromosome II

```
<i>Vibrio cholerae</i> Chromosome I
![Vibrio cholerae chromosome I](https://github.com/satoshikawato/gbdraw/blob/main/examples/NC_009457.svg)
<i>Vibrio cholerae</i> Chromosome II
![Vibrio cholerae chromosome II](https://github.com/satoshikawato/gbdraw/blob/main/examples/NC_009456.svg)

### Human herpesvirus 6 (HHV-6)
`gbdraw` can draw linear genomes and pairwise matches depicting similar genomic regions.
Perform a pariwise BLASTN analysis between the genomes of HHV-6A and HHV-6B:
```bash
blastn -query NC_000898.fasta -subject NC_001664.fasta -outfmt 7 -out NC_000898_NC_001664.blastn.out
```
Use `gbdraw linear` to visualize the comparison, aligning genomes to the center and separating forward and reverse strands:
```bash
umamba activate blast-2.16.0
blastn -query NC_000898.fasta -subject NC_001664.fasta -outfmt 7 -out NC_000898_NC_001664.blastn.out
umamba deactivate
umamba activate gbdraw-0.1.0
gbdraw linear --gbk NC_000898.gb NC_001664.gb -b NC_000898_NC_001664.blastn.out --resolve_overlaps --align_center --separate_strands -o HHV-6 -f svg --block_stroke_width 0.5

umamba deactivate
```
![HHV-6](https://github.com/satoshikawato/gbdraw/blob/main/examples/HHV-6.svg)

### <i>Candidatus</i> Hepatoplasmataceae (mollicutes)
```bash
tblastx -query Fukuoka2020.fasta -subject Av-JP.fasta  -outfmt 7 -out Fukuoka2020_Av-JP.tblastx.out
tblastx -query Av-JP.fasta -subject Ps-JP.fasta  -outfmt 7 -out Av-JP_Ps-JP.tblastx.out
tblastx -query Ps-JP.fasta -subject Tokyo2021.fasta  -outfmt 7 -out Ps-JP_Tokyo2021.tblastx.out
tblastx -query Tokyo2021.fasta -subject Av.fasta  -outfmt 7 -out Tokyo2021_Av.tblastx.out
gbdraw linear --gbk Fukuoka2020.gb Av-JP.gb Ps-JP.gb Tokyo2021.gb Av.gb -b Fukuoka2020_Av-JP.tblastx.out  Av-JP_Ps-JP.tblastx.out  Ps-JP_Tokyo2021.tblastx.out  Tokyo2021_Av.tblastx.out -o hepatoplasmataceae --align_center --bitscore 50 --evalue 1e-3 --separate_strands -f svg
```
![hepatoplasmataceae](https://github.com/satoshikawato/gbdraw/blob/main/examples/hepatoplasmataceae.svg)

[Home](../README.md) | [Installation](./INSTALL.md) | [Usage](./USAGE.md) | [Gallery](./GALLERY.md) | [RECIPES](./RECIPES.md) | [Advanced](./ADVANCED.md) | [FAQ](./FAQ.md)