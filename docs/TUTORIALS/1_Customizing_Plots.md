[Home](../DOCS.md) | [Installation](../INSTALL.md) | [Quickstart](../QUICKSTART.md) | [Tutorials](./TUTORIALS.md) | [Recipes](../RECIPES.md) | [CLI Reference](../CLI_Reference.md) | [Gallery](../GALLERY.md) | [FAQ](../FAQ.md) | [About](../ABOUT.md)

[< Back to Quickstart](../QUICKSTART.md) | [Next: Draw precomputed BLAST links >](./2_Comparative_Genomics.md)

[< Back to the guide index](./TUTORIALS.md)

# Style a circular genome diagram

Use palettes, track presets, center labels, plot titles, and feature labels to change a circular genome diagram.

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

The montage compares `tuckin`, `middle`, and `spreadout` from left to right. All three diagrams use `--separate_strands`.

![Circular WSSV diagrams comparing tuckin, middle, and spreadout presets with separate strand tracks](../../examples/track_layout_separate_strands.png)

## 3. Add a center label or plot title

Use `--species` and `--strain` to label the center of the diagram:

```bash
gbdraw circular \
  --gbk NC_000913.gbk \
  -o ecoli_with_title \
  -f svg \
  --separate_strands \
  --species "<i>Escherichia coli</i>" \
  --strain "K-12"
```

The species name is italicized and shown with the strain at the center of the diagram.

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

External circular labels are horizontal by default. Use radial placement when each external label should follow the circle radius:

```bash
gbdraw circular \
  --gbk tests/test_inputs/NC_001879.gbk \
  --labels both \
  --label_placement radial \
  --qualifier_priority examples/qualifier_priority.tsv \
  -o circular_radial_labels
```

Radial placement works with `auto`, `embedded_only`, and `external_only`. Unlike horizontal `--labels both`, radial placement keeps requested GC and GC-skew tracks and expands the radial/canvas geometry when inner or long labels need more room.

![Circular white spot syndrome virus genome with feature labels placed within and around the map](../../examples/WSSV_with_labels.svg)

> [!WARNING]
> Avoid `--labels` or `--labels both` on feature-dense genomes unless you also filter labels with `--label_blacklist` or `--label_whitelist`.

## 5. Simplify the diagram

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

## 6. Mark chloroplast genome regions

Region annotations are independent of GenBank features. This example adds LSC, SSC, IRa, and IRb brackets to the *Nicotiana tabacum* chloroplast map while retaining its feature colors and inner GC-content track.

The repository includes the NC_001879.2 GenBank record, chloroplast color rules, label priority, and region table. The four intervals follow the JLB, JSB, JSA, and JLA boundaries in the record:

| Region | Start | End | Length |
| --- | ---: | ---: | ---: |
| LSC | 1 | 86,686 | 86,686 bp |
| IRb | 86,687 | 112,029 | 25,343 bp |
| SSC | 112,030 | 130,600 | 18,571 bp |
| IRa | 130,601 | 155,943 | 25,343 bp |

Run the example from the repository root:

```bash
gbdraw circular \
  --gbk tests/test_inputs/NC_001879.gbk \
  --annotation_table examples/nicotiana-tabacum-regions.tsv \
  -t examples/chloroplast_specific_table.tsv \
  --qualifier_priority examples/qualifier_priority.tsv \
  -k CDS,rRNA,tRNA,tmRNA,ncRNA,misc_RNA,rep_origin \
  --separate_strands \
  --block_stroke_width 1 \
  --block_stroke_color black \
  --axis_stroke_width 3 \
  --line_stroke_width 2 \
  --track_type tuckin \
  --labels both \
  --label_placement radial \
  --outer_label_x_radius_offset 0.90 \
  --outer_label_y_radius_offset 0.90 \
  --inner_label_x_radius_offset 0.975 \
  --inner_label_y_radius_offset 0.975 \
  --suppress_skew \
  --species "<i>Nicotiana tabacum</i>" \
  --definition_font_size 28 \
  --legend upper_left \
  --circular_track_slot 'features:features@side=overlay,lane_direction=split' \
  --circular_track_slot 'plastome_regions:annotations@set_id=plastome_regions,side=inside,r=0.65,w=20px,show_labels=true,padding_px=1,overflow=compress,inner_gap_px=1,outer_gap_px=1' \
  --circular_track_slot 'gc_content:dinucleotide_content@side=inside,r=0.56,w=0.08' \
  -o NC_001879_regions \
  -f svg
```

![Nicotiana tabacum chloroplast map with bracket annotations for LSC, SSC, IRa, and IRb](../../examples/NC_001879_regions.svg)

The brackets share one compact annotation lane between the feature map and GC-content track. Radial feature labels allow the GC track to remain visible with `--labels both`. Open the [interactive Gallery version](https://gbdraw.app/gallery/#tobacco-chloroplast) to inspect the same plot as an interactive SVG or restore its saved session.

Coordinates in annotation tables are 1-based and inclusive. Use `wraps_origin=true` with `start > end` for an origin-spanning circular annotation. See [Use TSV manifests for CLI inputs](./5_Table_Driven_Inputs.md#7-region-annotation-tables) for feature targets, styles, and the complete table schema.

In the web app, open **Region Annotations**, import `nicotiana-tabacum-regions.tsv`, then choose **Annotations** in **Custom Track Slots** and bind the track to `plastome_regions`.

## 7. When you move to linear mode

Linear mode has its own input selectors:

- `--record_id`
- `--reverse_complement`
- `--region`

See [Arrange linear tracks, record labels, and rulers](./7_Linear_Layout.md#7-select-records-regions-and-orientation) and the [CLI Reference](../CLI_Reference.md) for details.

[< Back to Quickstart](../QUICKSTART.md) | [Next: Draw precomputed BLAST links >](./2_Comparative_Genomics.md)

[< Back to the guide index](./TUTORIALS.md)

[Home](../DOCS.md) | [Installation](../INSTALL.md) | [Quickstart](../QUICKSTART.md) | [Tutorials](./TUTORIALS.md) | [Recipes](../RECIPES.md) | [CLI Reference](../CLI_Reference.md) | [Gallery](../GALLERY.md) | [FAQ](../FAQ.md) | [About](../ABOUT.md)
