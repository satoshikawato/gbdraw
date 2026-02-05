[Home](../README.md) | [Installation](../INSTALL.md) | [Quickstart](../QUICKSTART.md) | [**Tutorials**](../TUTORIALS/TUTORIALS.md) | [Gallery](../GALLERY.md) | [FAQ](../FAQ.md) | [ABOUT](../ABOUT.md)


[< Back to the Index of Tutorials](./TUTORIALS.md)
[< Back to Tutorial 1: Customizing Your Plot](./1_Customizing_Plots.md)　　　　　　[Go to Tutorial 3: Advanced Customization >](./3_Advanced_Customization.md) 


# Tutorial 2: Comparative Genomics with BLAST
**Goal**: Learn how to use `gbdraw linear` to visualize sequence similarity between two or more genomes using BLAST results.

---

### 1. Required Inputs

To create a comparative plot, you need two types of files:

1.  **Genome Files**: Two or more genome annotation files in GenBank or GFF3+FASTA format.
2.  **BLAST Result File**: A tab- or comma-separated file (`-outfmt 7` or `-outfmt 6`) containing the results of a BLAST comparison between the genomes.

---

### 2. BLAST Output Format

Let's prepare an example using *Shigella*, the causative agent of bacterial dysentery and is genetically nested within the species *E.coli*.

1.  **Download the sample data (GenBank and FASTA).**
    ```bash
    # Escherichia coli str. K-12 substr. MG1655 (NC_000913.3)
    wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_000913.3&rettype=gbwithparts&retmode=text" -O Escherichia_coli.gbk
    wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_000913.3&rettype=fasta&retmode=text" -O Escherichia_coli.fasta

    # Shigella dysenteriae SWHEFF_49 (NZ_CP055055.1)
    wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NZ_CP055055.1&rettype=gbwithparts&retmode=text" -O Shigella_dysenteriae.gbk
    wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NZ_CP055055.1&rettype=fasta&retmode=text" -O Shigella_dysenteriae.fasta
    ```

2.  **Run BLAST.**
    ```bash
    # Run blastn and save the results in outfmt 7 or 6
    blastn -query Escherichia_coli.fasta -subject Shigella_dysenteriae.fasta -outfmt 7 -out Escherichia_coli-Shigella_dysenteriae.blastn.out

    # blastn -query Escherichia_coli.fasta -subject Shigella_dysenteriae.fasta -outfmt 6 -out Escherichia_coli-Shigella_dysenteriae.blastn.out

    ```
This creates the `Escherichia_coli-Shigella_dysenteriae.blastn.out` file we need.

---

### 3. Generate the Comparative Plot

Now, provide both the GenBank files and the BLAST result file to `gbdraw linear`.

* `--gbk`: List all genome files.
* `-b` or `--blast`: List the BLAST result file that connects them.

```bash
gbdraw linear \
  --gbk Escherichia_coli.gbk Shigella_dysenteriae.gbk \
  -b Escherichia_coli-Shigella_dysenteriae.blastn.out \
  --align_center \
  --separate_strands \
  -o Escherichia_Shigella_pairwise \
  -f svg
```

The resulting file, `Escherichia_Shigella_pair.svg`, will show the two genomes as tracks, with ribbons connecting the regions of sequence similarity identified by BLAST.
![Escherichia_Shigella_pair.svg](../../examples/Escherichia_Shigella_pair.svg)

### 4. Selecting Records and Regions (Advanced Input Control)

When a GenBank or GFF3+FASTA file contains multiple records, you can target
specific inputs for linear mode:

- `--record_id`: Select a record by ID or `#index` per input file (repeatable).
- `--reverse_complement`: Reverse-complement a record per input file.
- `--region`: Crop a region using `record_id:start-end[:rc]`.

Example (single file, select one record and crop a region):
```bash
gbdraw linear \
  --gbk Escherichia_coli.gbk \
  --record_id NC_000913.3 \
  --region NC_000913.3:100000-250000 \
  -o Escherichia_coli_region \
  -f svg
```

Example (multi-file, reverse-complement the second record):
```bash
gbdraw linear \
  --gbk Genome1.gbk Genome2.gbk \
  --reverse_complement false true \
  -o Genome1_Genome2_rc \
  -f svg
```

If you use `--region` or `--reverse_complement` with BLAST, ensure the BLAST
coordinates match the cropped/reversed records.

### 5. Comparing Multiple Genomes
You can compare more than two genomes by providing them in sequence. Ensure that you provide a BLAST file for each adjacent pair.

**Example**: To compare `Genome1 -> Genome2 -> Genome3`, you need:

- BLAST result for Genome1 vs. Genome2

- BLAST result for Genome2 vs. Genome3


```bash
# Shigella flexneri 2a str. 301 (NC_004337.2)
wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_004337.2&rettype=gbwithparts&retmode=text" -O Shigella_flexneri.gbk
wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_004337.2&rettype=fasta&retmode=text" -O Shigella_flexneri.fasta

# Shigella sonnei ATCC 29930 (NZ_CP026802.1)
wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NZ_CP026802.1&rettype=gbwithparts&retmode=text" -O Shigella_sonnei.gbk
wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NZ_CP026802.1&rettype=fasta&retmode=text" -O Shigella_sonnei.fasta
```


```bash
blastn -query Shigella_dysenteriae.fasta -subject Shigella_flexneri.fasta -outfmt 7 -out Shigella_dysenteriae-Shigella_flexneri.blastn.out
blastn -query Shigella_flexneri.fasta -subject Shigella_sonnei.fasta -outfmt 7 -out Shigella_flexneri-Shigella_sonnei.blastn.out
```
> [!IMPORTANT]
> The order of the `--gbk` and `-b` files must correspond correctly.

> [!TIP]
> You can filter BLAST ribbons by setting thresholds for E-value (`--evalue`), bitscore (`--bitscore`), or identity (`--identity`).

```bash
gbdraw linear \
  --gbk Escherichia_coli.gbk Shigella_dysenteriae.gbk Shigella_flexneri.gbk Shigella_sonnei.gbk \
  -b Escherichia_coli-Shigella_dysenteriae.blastn.out Shigella_dysenteriae-Shigella_flexneri.blastn.out Shigella_flexneri-Shigella_sonnei.blastn.out \
  --align_center \
  --separate_strands \
  --evalue 1e-99 --bitscore 5000 \
  -o Escherichia_Shigella_multi \
  -f svg
```

![Escherichia_Shigella_multi.svg](../../examples/Escherichia_Shigella_multi.svg)

[< Back to the Index of Tutorials](./TUTORIALS.md)
[< Back to Tutorial 1: Customizing Your Plot](./1_Customizing_Plots.md)　　　　　　[Go to Tutorial 3: Advanced Customization >](./3_Advanced_Customization.md) 

[Home](../README.md) | [Installation](../INSTALL.md) | [Quickstart](../QUICKSTART.md) | [**Tutorials**](../TUTORIALS/TUTORIALS.md) | [Gallery](../GALLERY.md) | [FAQ](../FAQ.md) | [ABOUT](../ABOUT.md)
