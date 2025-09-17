[Home](./README.md) | [Installation](./INSTALL.md) | [Quickstart](./QUICKSTART.md) | [Tutorials](./TUTORIALS/) | **FAQ**


# Frequently Asked Questions (FAQ)

Here are answers to some common questions about `gbdraw`.

---

### Q: Can I generate a plot if I only have a GFF3 file?
**A:** Yes. You can use a GFF3 file by providing the annotation with `--gff` and the corresponding sequence with `--fasta`.
```bash
gbdraw circular --gff my_genome.gff --fasta my_genome.fasta -o my_plot
```
### Q: How can I make the GC content graph smoother or more detailed?
**A:** You can control the graph's resolution with the --window and --step options. A larger window size will produce a smoother plot, while a smaller size will show more detail.

### Q: My feature labels are overlapping and unreadable. How can I fix this?

**A:** You have several options to deal with label clutter:

1.  **Reduce Font Size**: The simplest fix is to make the labels smaller using `--label_font_size`.
    ```bash
    gbdraw circular --show_labels --label_font_size 6 ...
    ```
2.  **Blacklist Common Labels**: Use `--label_blacklist` to hide non-informative labels like "hypothetical protein". This often cleans up the plot significantly.
    ```bash
    gbdraw circular --show_labels --label_blacklist "hypothetical protein" ...
    ```
3.  **Whitelist Key Labels**: If you only care about a few specific genes, use `--label_whitelist` to show *only* those labels.

For more details, see the **[Advanced Customization Tutorial](./TUTORIALS/3_Advanced_Customization.md)**.

---

### Q: How can I change the color of just one specific gene?

**A:** Use a feature-specific color table with the `-t` option. This lets you create rules to color individual features based on their annotations (like product name or locus tag).

For a complete guide, see the "Feature-Specific Colors" section in the **[Advanced Customization Tutorial](./TUTORIALS/3_Advanced_Customization.md)**.

---

### Q: My comparative plot (using `-b`) isn't showing any similarity ribbons. What's wrong?

**A:** This is almost always due to one of two issues:

1.  **Incorrect BLAST Format**: `gbdraw` strictly requires BLAST tabular output with headers, which is **`outfmt 7`**. Standard `outfmt 6` will not work. Re-run your BLAST search with the `-outfmt 7` flag.
2.  **Incorrect File Order**: The order of genome files (`--gbk`) and BLAST files (`-b`) must correspond. For `gbdraw linear --gbk A B C -b A_vs_B B_vs_C`, the first BLAST file must be the comparison between A and B, and the second between B and C.

See the **[Comparative Genomics Tutorial](./TUTORIALS/2_Comparative_Genomics.md)** for a working example.

---

### Q: Can I use the gene name (e.g., *dnaA*) for labels instead of the full product description?

**A:** Yes. By default, `gbdraw` prefers the `product` qualifier. You can change this behavior with the `--qualifier_priority` option.

Create a simple TSV file (e.g., `priority.tsv`) that tells `gbdraw` to look for the `gene` qualifier first for CDS features:
```tsv
# priority.tsv
CDS	gene
```
Then, use it in your command:
```bash
gbdraw circular --show_labels --qualifier_priority priority.tsv ...
```
