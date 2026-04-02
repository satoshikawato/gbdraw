[Home](./DOCS.md) | [Installation](./INSTALL.md) | [Quickstart](./QUICKSTART.md) | [Tutorials](./TUTORIALS/TUTORIALS.md) | [Recipes](./RECIPES.md) | [CLI Reference](./CLI_Reference.md) | [Gallery](./GALLERY.md) | **FAQ** | [About](./ABOUT.md)

# Frequently Asked Questions

## Is there a web GUI? Do I need Streamlit?

Use [https://gbdraw.app/](https://gbdraw.app/) for the hosted app, or run `gbdraw gui` locally after installation. Streamlit is not required.

## Why do my CLI and browser renders differ slightly?

Small differences in label placement and legend sizing are expected. The CLI uses kerning-aware font metrics, while the web UI uses browser text metrics.

## Can I use a GFF3 file by itself?

No. `gbdraw` requires both annotation and sequence data. When using GFF3 input, provide the matching FASTA file with `--fasta`.

```bash
gbdraw circular --gff my_genome.gff --fasta my_genome.fasta -o my_plot
```

## My labels overlap. What should I do?

Common fixes:

1. Reduce `--label_font_size`
2. Hide noisy labels with `--label_blacklist`
3. Keep only important labels with `--label_whitelist`
4. Switch to `--track_type middle` for circular plots or reduce the number of displayed labels

See [Tutorial 3](./TUTORIALS/3_Advanced_Customization.md) for examples.

## How do I change the color of one specific gene?

Use a feature-specific color table with `-t`. This matches selected features by qualifier values and assigns a color and legend label.

See [Tutorial 3](./TUTORIALS/3_Advanced_Customization.md) and [Recipes](./RECIPES.md).

## My comparative plot has no ribbons. What is usually wrong?

The most common causes are:

1. The BLAST file is not in outfmt 6 or 7
2. The BLAST file order does not match the genome input order
3. Filtering thresholds such as `--evalue`, `--bitscore`, `--identity`, or `--alignment_length` are too strict

See [Tutorial 2](./TUTORIALS/2_Comparative_Genomics.md) for a working example.

## Can I use gene names instead of product descriptions for labels?

Yes. Use `--qualifier_priority` to prefer `gene`, `locus_tag`, or other qualifiers.

```tsv
CDS	gene
```

```bash
gbdraw circular --gbk genome.gb --labels --qualifier_priority priority.tsv -o output -f svg
```

## How do I make the GC curve smoother?

Increase the window and step sizes:

```bash
gbdraw circular --gbk genome.gb --window 10000 --step 1000 -o output -f svg
```

![window_step_comparison.png](../examples/window_step_comparison.png)

## Can I plot AT instead of GC?

Yes. Use `--nt AT`.

```bash
gbdraw circular --gbk genome.gb --nt AT -o output -f svg
```

![skew_comparison.png](../examples/skew_comparison.png)

## Why does SVG export work but PNG/PDF/EPS/PS export fail?

Non-SVG export requires CairoSVG. Install the optional export dependency and, if needed on your platform, the system Cairo/Pango libraries.

## Are there known visualization limitations?

- Trans-introns are not currently visualized.
- Mixed-format text such as `<i>Ca.</i> Tyloplasma litorale` does not reliably survive SVG-to-PNG/PDF/EPS/PS conversion. Use SVG if you need exact mixed formatting.

[Home](./DOCS.md) | [Installation](./INSTALL.md) | [Quickstart](./QUICKSTART.md) | [Tutorials](./TUTORIALS/TUTORIALS.md) | [Recipes](./RECIPES.md) | [CLI Reference](./CLI_Reference.md) | [Gallery](./GALLERY.md) | **FAQ** | [About](./ABOUT.md)
