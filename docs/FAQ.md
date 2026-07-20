[Home](./DOCS.md) | [Installation](./INSTALL.md) | [Quickstart](./QUICKSTART.md) | [Tutorials](./TUTORIALS/TUTORIALS.md) | [Recipes](./RECIPES.md) | [CLI Reference](./CLI_Reference.md) | [Gallery](./GALLERY.md) | **FAQ** | [About](./ABOUT.md)

# Frequently asked questions

## Is there a web GUI? Do I need Streamlit?

Use [https://gbdraw.app/](https://gbdraw.app/) for the hosted app, or run `gbdraw gui` locally after installation. Streamlit is not required. Local GUI analysis runs on your machine; the interactive gallery examples are hosted at [https://gbdraw.app/gallery/](https://gbdraw.app/gallery/).

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
3. Keep only important labels with `--label_whitelist` regex patterns
4. Use the `--track_type middle` circular preset or reduce the number of displayed labels

See [Set feature colors and labels](./TUTORIALS/3_Advanced_Customization.md) for examples.

## How do I change the color of one specific gene?

Use a feature-specific color table with `-t`. This matches selected features by qualifier values and assigns a color and legend label.

See [Set feature colors and labels](./TUTORIALS/3_Advanced_Customization.md) and [Recipes](./RECIPES.md).

## How do I mark a coordinate range or a group of features?

Use `--annotation_table` and bind each `set_id` to an `annotations` custom track slot. Coordinate targets are 1-based and inclusive. Feature targets use the qualifiers already loaded from GenBank or GFF3, such as `locus_tag=ABC_001`.

The same table works in Circular and Linear mode. See [Region annotation tables](./TUTORIALS/5_Table_Driven_Inputs.md#7-region-annotation-tables).

## My comparative diagram has no ribbons. What is usually wrong?

The most common causes are:

1. The BLAST file is not in outfmt 6 or 7
2. The BLAST file order does not match the genome input order
3. Filtering thresholds such as `--evalue`, `--bitscore`, `--identity`, or `--alignment_length` are too strict

See [Draw genome comparison links from precomputed BLAST results](./TUTORIALS/2_Comparative_Genomics.md) for a working example.

## What if one record has no Depth TSV for a sample?

Keep the record's position in the repeated `--depth_track` group with a quoted empty argument:

```bash
gbdraw linear \
  --gbk record-a.gb record-b.gb \
  --depth_track record-a.depth.tsv '' \
  -o depth-partial \
  -f svg
```

Use `--depth_track '' record-b.depth.tsv` when only the second record has data. The empty argument means that the logical series is missing for that record. gbdraw does not substitute another file or draw zero coverage. In Linear mode the missing cell also reserves no vertical space, while the series identity remains stable for records that do contain data. Each group must contain at least one real file. See [Plot read depth and other numeric tracks](./TUTORIALS/6_Depth_Quantitative_Tracks.md#3-compare-depth-across-records) for a runnable example.

## Why is my circular BLAST similarity ring empty?

Check that the BLAST file is outfmt 6 or 7, the displayed circular record ID appears on the side selected by `--conservation_reference`, and the thresholds are not too strict. When BLAST was generated as `blastn -query comparison.fasta -subject reference.fasta`, use `--conservation_reference subject`.

These rings draw raw HSP spans; they do not infer evolutionary conservation. A BLAST row where the selected reference start is greater than the selected reference end is treated as reverse orientation, not as a hit crossing the circular origin. The current implementation does not infer binned or wraparound hits.

## Can pairwise comparison links be curved?

Yes. In linear mode, `--pairwise_match_style ribbon` draws straight filled ribbons by default. Use `--pairwise_match_style curve` to bend the same match spans; curved links can be easier to distinguish in a dense comparison diagram.

## Can I use gene names instead of product descriptions for labels?

Yes. Use `--qualifier_priority` to prefer `gene`, `locus_tag`, or other qualifiers.

```tsv
CDS	gene
```

```bash
gbdraw circular --gbk genome.gb --labels --qualifier_priority priority.tsv -o output -f svg
```

## How do I make the GC content track smoother?

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
