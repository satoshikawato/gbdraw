[Home](../DOCS.md) | [Installation](../INSTALL.md) | [Quickstart](../QUICKSTART.md) | [Tutorials](./TUTORIALS.md) | [Recipes](../RECIPES.md) | [CLI Reference](../CLI_Reference.md) | [Gallery](../GALLERY.md) | [FAQ](../FAQ.md) | [About](../ABOUT.md)

[< Back to Quickstart](../QUICKSTART.md) | [Go to Tutorial 2 >](./2_Comparative_Genomics.md)

[< Back to the Tutorials Index](./TUTORIALS.md)

# Tutorial 1: Customize a circular plot

Start here for the main circular-plot styling controls: palettes, track presets, centered definition text, titles, and labels.

## 1. Change the color scheme

`gbdraw` ships with many built-in palettes. Use `-p` or `--palette` to select one.

If you did not complete the [Quickstart](../QUICKSTART.md), download the *Escherichia coli* K-12 record first:

```bash
wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_000913.3&rettype=gbwithparts&retmode=text" -O NC_000913.gbk
```

```bash
gbdraw circular \
  --gbk NC_000913.gbk \
  -o ecoli_orchid \
  -f svg \
  --separate_strands \
  -p orchid
```

The SVG uses the orchid palette while keeping forward- and reverse-strand features on separate tracks.

![Circular E. coli K-12 diagram in the orchid palette with separate strand tracks](../../examples/ecoli_orchid.svg)

See the [built-in palette examples](../../examples/color_palette_examples.md) for the full palette list.

## 2. Choose a circular preset

These options control circular feature-track layout:

- `--track_type`: preset name, `tuckin`, `middle`, or `spreadout`
- `--separate_strands`: split forward and reverse features into different tracks

`tuckin` is the compact default. `middle` is often easier to label, while `spreadout` provides the widest track separation.

Custom track slots can reorder tracks without replacing the preset geometry. When a built-in slot omits its radius, width, physical gaps, placement, or standard renderer parameters, it inherits them from `--track_type`. Values set on the slot override those defaults. Use `inner_gap_px` and `outer_gap_px` to set radial gaps.

The circular axis radius is fixed at `canvas.circular.radius`; it cannot be moved or hidden with a circular track slot.

The montage compares `tuckin`, `middle`, and `spreadout` from left to right. All three plots use `--separate_strands`.

![Circular WSSV diagrams comparing tuckin, middle, and spreadout presets with separate strand tracks](../../examples/track_layout_separate_strands.png)

## 3. Add centered organism text or a plot title

Use `--species` and `--strain` to control the centered definition text:

```bash
gbdraw circular \
  --gbk NC_000913.gbk \
  -o ecoli_with_title \
  -f svg \
  --separate_strands \
  --species "<i>Escherichia coli</i>" \
  --strain "K-12"
```

The species name is italicized and shown with the strain in the centered definition text.

![Circular E. coli K-12 diagram with an italic species name and centered strain text](../../examples/ecoli_with_title.svg)

Use `--plot_title` with `--plot_title_position top` or `bottom` to add a title above or below the plot.

> [!CAUTION]
> Mixed-format text such as `<i>Ca.</i> Tyloplasma litorale` does not reliably survive SVG-to-PNG/PDF/EPS/PS conversion. Use SVG when exact formatting matters.

## 4. Show feature labels

Circular labels are hidden by default. Use `--labels` to show outer labels, or `--labels both` to show outer and inner labels.

This example uses the white spot syndrome virus genome so that individual labels remain legible:

```bash
wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=AP027280.1&rettype=gbwithparts&retmode=text" -O AP027280.gb
```

```bash
gbdraw circular \
  --gbk AP027280.gb \
  -o WSSV_with_labels \
  -f svg \
  --block_stroke_width 1 \
  --track_type middle \
  --labels
```

With the default `--label_rendering auto` mode, labels that fit are embedded in feature bodies and the remaining labels are routed outside the circle.

![Circular white spot syndrome virus genome with feature labels placed within and around the map](../../examples/WSSV_with_labels.svg)

> [!WARNING]
> Avoid `--labels` or `--labels both` on feature-dense genomes unless you also filter labels with `--label_blacklist` or `--label_whitelist`.

## 5. Simplify the plot

Suppress the GC tracks and legend to leave more space for annotated features:

```bash
gbdraw circular \
  --gbk AP027280.gb \
  -o WSSV_filtered \
  -f svg \
  --block_stroke_width 1 \
  --suppress_gc \
  --suppress_skew \
  --separate_strands \
  --labels \
  --legend none
```

The output remains a labeled, strand-separated feature map.

![Simplified circular white spot syndrome virus map without GC tracks, GC skew, or a legend](../../examples/WSSV_filtered.svg)

## 6. When you move to linear mode

Linear mode has its own input selectors:

- `--record_id`
- `--reverse_complement`
- `--region`

See [Tutorial 2](./2_Comparative_Genomics.md) and the [CLI Reference](../CLI_Reference.md) for details.

[< Back to Quickstart](../QUICKSTART.md) | [Go to Tutorial 2 >](./2_Comparative_Genomics.md)

[< Back to the Tutorials Index](./TUTORIALS.md)

[Home](../DOCS.md) | [Installation](../INSTALL.md) | [Quickstart](../QUICKSTART.md) | [Tutorials](./TUTORIALS.md) | [Recipes](../RECIPES.md) | [CLI Reference](../CLI_Reference.md) | [Gallery](../GALLERY.md) | [FAQ](../FAQ.md) | [About](../ABOUT.md)
