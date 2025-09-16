# Tutorial 1: Customizing Your Plot

[< Back to Tutorials](./)

**Goal**: Learn the basics of customizing plot appearance, including colors, titles, and labels, using the plot from the Quickstart.

---

### 1. Changing the Color Scheme

`gbdraw` comes with [55 built-in color palettes](../examples/color_palette_examples.md). You can specify one using the `-p` or `--palette` option. Let's try the `ajisai` palette.

```bash
gbdraw circular \
  --gbk GCF_000005845.2_ASM584v2_genomic.gbff \
  -o ecoli_ajisai \
  -f svg \
  -p ajisai
```

The new file `ecoli_ajisai.svg` will be generated with a different color scheme. You can see examples of all available palettes [here](../examples/color_palette_examples.md).


### 2. Adding a Title
Use the `--species` and `--strain` options to add a title to the center of your plot. You can use HTML `<i>` tags for *italics*.

```bash
gbdraw circular \
  --gbk GCF_000005845.2_ASM584v2_genomic.gbff \
  -o ecoli_with_title \
  -f svg \
  -p ajisai \
  --species "<i>Escherichia coli</i>" \
  --strain "K-12"
```
The plot will now have a formatted title in the center.

### 3. Showing Feature Labels
By default, gene labels are hidden. Use the `--show_labels` option to display them.
```bash
gbdraw circular \
  --gbk GCF_000005845.2_ASM584v2_genomic.gbff \
  -o ecoli_with_labels \
  -f svg \
  -p ajisai \
  --species "<i>Escherichia coli</i>" \
  --strain "K-12" \
  --show_labels
```
For large genomes, labels may overlap. You can adjust their size with `--label_font_size` or use advanced filtering techniques covered in a later tutorial.

### Next Steps
**Tutorial 2**: Comparative Genomics: Learn how to compare multiple genomes.

**Tutorial 3**: Advanced Customization: Learn how to control colors and labels with configuration files.