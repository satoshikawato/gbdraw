# Quickstart: Your First Plot in 5 Minutes

[Home](./README.md) | [Installation](./INSTALL.md) | **Quickstart** | [Tutorials](./TUTORIALS/) | [Gallery](./GALLERY.md)

This tutorial will guide you from a fresh installation to generating your first circular genome plot.

### 1. Prerequisites

Ensure you have `gbdraw` installed locally via one of the methods described on the **[Installation](./INSTALL.md)** page. Make sure your conda environment is activated.

### 2. Get Sample Data

For this tutorial, we will use the GenBank file for *Escherichia coli* K-12. Download and decompress it with the following commands:

```bash
wget [https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gbff.gz](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gbff.gz)
gunzip GCF_000005845.2_ASM584v2_genomic.gbff.gz
```

This will give you the file GCF_000005845.2_ASM584v2_genomic.gbff.


### 3. Generate the Plot
In the same directory as the downloaded file, run the following command:
```bash
gbdraw circular --gbk GCF_000005845.2_ASM584v2_genomic.gbff -o ecoli_k12_plot -f svg
```

This command tells gbdraw to:

- `circular`: Create a circular diagram.

- `--gbk ...`: Use the specified GenBank file as input.

- `-o ecoli_k12_plot`: Set the prefix for the output filename.

- `-f svg`: Set the output format to SVG (a scalable vector format).

### 4. Check Your Output

A new file named `ecoli_k12_plot.svg` will appear in your directory. Open it in a web browser or vector graphics editor (like Inkscape or Illustrator). You should see a complete genome map of *E. coli*!

### 5. Next Steps

Congratulations on creating your first plot!

To learn how to change colors, add titles, and show labels, continue to [Tutorial 1: Customizing Your Plot](./TUTORIALS.md).

To see more examples of what `gbdraw` can do, check out the [Gallery](./GALLERY.md).