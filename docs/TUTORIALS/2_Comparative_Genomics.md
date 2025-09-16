# Tutorial 2: Comparative Genomics with BLAST

[< Back to Tutorials](./)

**Goal**: Learn how to use `gbdraw linear` to visualize sequence similarity between two or more genomes using BLAST results.

---

### 1. Required Inputs

To create a comparative plot, you need two types of files:

1.  **Genome Files**: Two or more genome annotation files in GenBank (`.gbk`) or GFF3+FASTA format.
2.  **BLAST Result File**: A file containing the results of a BLAST comparison between the genomes.

---

### 2. BLAST Output Format

This is a critical step. `gbdraw` requires the BLAST output to be in the tabular format with comment lines, specified as **`outfmt 7`**. These comment lines contain the names of the query and subject sequences, which `gbdraw` uses to map the results to the correct genomes.

Let's prepare an example using two Human herpesvirus genomes.

1.  **Download the sample data (GenBank and FASTA).**
    ```bash
    # HHV-6A
    wget "[https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?db=nuccore&dopt=gbwithparts&id=NC_001664.4](https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?db=nuccore&dopt=gbwithparts&id=NC_001664.4)" -O NC_001664.gbk
    wget "[https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?db=nuccore&dopt=fasta&id=NC_001664.4](https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?db=nuccore&dopt=fasta&id=NC_001664.4)" -O NC_001664.fasta
    # HHV-6B
    wget "[https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?db=nuccore&dopt=gbwithparts&id=NC_000898.1](https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?db=nuccore&dopt=gbwithparts&id=NC_000898.1)" -O NC_000898.gbk
    wget "[https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?db=nuccore&dopt=fasta&id=NC_000898.1](https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?db=nuccore&dopt=fasta&id=NC_000898.1)" -O NC_000898.fasta
    ```

2.  **Run BLAST.**
    ```bash
    # Create a BLAST database for the subject genome
    makeblastdb -in NC_000898.fasta -dbtype nucl

    # Run blastn and save the results in outfmt 7
    blastn -query NC_001664.fasta -db NC_000898.fasta -outfmt 7 -out hhv6a_vs_hhv6b.blast
    ```
This creates the `hhv6a_vs_hhv6b.blast` file we need.

---

### 3. Generate the Comparative Plot

Now, provide both the GenBank files and the BLAST result file to `gbdraw linear`.

* `--gbk`: List all genome files.
* `-b` or `--blast`: List the BLAST result file that connects them.

```bash
gbdraw linear \
  --gbk NC_001664.gbk NC_000898.gbk \
  -b hhv6a_vs_hhv6b.blast \
  --align_center \
  -o HHV6_comparison \
  -f svg
```

The resulting file, `HHV6_comparison.svg`, will show the two genomes as tracks, with ribbons connecting the regions of sequence similarity identified by BLAST.

### 4. Comparing Multiple Genomes
You can compare more than two genomes by providing them in sequence. Ensure that you provide a BLAST file for each adjacent pair.

**Example**: To compare `Genome1 -> Genome2 -> Genome3`, you need:

- BLAST result for Genome1 vs. Genome2

- BLAST result for Genome2 vs. Genome3

The order of the `--gbk` and `-b` files must correspond correctly.
```bash
gbdraw linear \
  --gbk genome1.gbk genome2.gbk genome3.gbk \
  -b g1_vs_g2.blast g2_vs_g3.blast \
  -o multi_comparison
```

