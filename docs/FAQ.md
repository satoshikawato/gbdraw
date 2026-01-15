[Home](./DOCS.md) | [Installation](./INSTALL.md) | [Quickstart](./QUICKSTART.md) | [Tutorials](./TUTORIALS/TUTORIALS.md) | [Gallery](./GALLERY.md) | **FAQ** | [ABOUT](./ABOUT.md)


# Frequently Asked Questions (FAQ)

Here are answers to some common questions about `gbdraw`.

---

## Known issues
- **Trans-introns** are not currently visualized.
- **Mixed-format text** (e.g., combining italic and block elements like `<i>Ca.</i> Tyloplasma litorale`) cannot be reliably converted from SVG to PDF/PNG/EPS/PS.  
  → As a workaround, export to **SVG format** and convert to other formats using external tools like [**Inkscape**](https://inkscape.org/).

---

### Q: Is there a web GUI? Do I need Streamlit?
**A:** Use [https://gbdraw.app/](https://gbdraw.app/) for the hosted web app, or run `gbdraw gui` locally after installing `gbdraw`. The `gbdraw gui` command starts a local HTTP server on a free port and opens the web UI in your browser. Streamlit is not required.

---

### Q: My labels or legend spacing changed after upgrading. Is that expected?

**A:** Yes. The CLI now uses kerning-aware font metrics for bounding box calculations, which can slightly change label placement and legend sizing. The web UI uses browser text metrics, so small differences between CLI and browser renders are expected.

---

### Q: Can I generate a plot if I only have a GFF3 file containing the genome sequence?
**A:** No. The genome sequence must be provided in a separate file. gbdraw cannot read the embedded sequence from the `##FASTA` section of a GFF3 file.
To generate a diagram, you must provide the GFF3 annotation with `--gff` and the corresponding sequence in a FASTA file with `--fasta`.
```bash
gbdraw circular --gff my_genome.gff --fasta my_genome.fasta -o my_plot
```

---

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

1.  **Incorrect BLAST Format**: `gbdraw` requires BLAST tabular output in **`outfmt 6`** or **`outfmt 7`** format. Both formats work correctly. Re-run your BLAST search with the `-outfmt 6` or `-outfmt 7` flag.
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

---

### Q: How can I make the GC content graph smoother/finer?
**A:** You can control the graph's resolution with the `--window (-w)` and `--step (-s)` options. A larger window size will produce a smoother plot, while a smaller size will show more detail.
```bash
gbdraw circular --window 10000 --step 1000 ...
```
![window_step_comparison.png](../examples/window_step_comparison.png)

---

### Q: Can I plot the AT content instead of GC content?
**A:** Yes. Use the --nt option to switch from the default GC to AT or any other dinucleotide.
```bash
gbdraw circular --nt AT ...
```
![skew_comparison.png](../examples/skew_comparison.png)

---

[Home](./DOCS.md) | [Installation](./INSTALL.md) | [Quickstart](./QUICKSTART.md) |[Tutorials](./TUTORIALS/TUTORIALS.md) | [Gallery](./GALLERY.md) | **FAQ** | [ABOUT](./ABOUT.md)
