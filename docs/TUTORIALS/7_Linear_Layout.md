[Home](../DOCS.md) | [Installation](../INSTALL.md) | [Quickstart](../QUICKSTART.md) | [Tutorials](./TUTORIALS.md) | [Recipes](../RECIPES.md) | [CLI Reference](../CLI_Reference.md) | [Gallery](../GALLERY.md) | [FAQ](../FAQ.md) | [About](../ABOUT.md)

[< Back to the Tutorials Index](./TUTORIALS.md)
[< Back to Tutorial 6](./6_Depth_Quantitative_Tracks.md) | [Go to Tutorial 8 >](./8_Interactive_SVG_Sessions.md)

# Tutorial 7: Linear Layout, Definitions, and Rulers

**Goal:** customize linear plots with track placement, rulers, record definitions, titles, and custom track slots.

## 1. Prepare Inputs

```bash
wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=LC738868.1&rettype=gbwithparts&retmode=text" -O MjeNMV.gb
wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=LC738874.1&rettype=gbwithparts&retmode=text" -O MelaMJNV.gb
```

## 2. Place Tracks Above, Middle, or Below

`--track_layout above`, `middle`, and `below` control where the feature track sits relative to the record axis.

```bash
gbdraw linear \
  --gbk MjeNMV.gb MelaMJNV.gb \
  --track_layout below \
  --track_axis_gap auto \
  --show_gc \
  --show_skew \
  -o majani_tracks_below \
  -f svg
```

Use `--track_axis_gap 12` when you want an explicit pixel gap instead of automatic spacing.

## 3. Use a Ruler on the Axis

`--ruler_on_axis` is effective when `--scale_style ruler` is used with `--track_layout above` or `below`.

```bash
gbdraw linear \
  --gbk MjeNMV.gb MelaMJNV.gb \
  --track_layout below \
  --scale_style ruler \
  --ruler_on_axis \
  --scale_interval 50000 \
  -o majani_axis_ruler \
  -f svg
```

## 4. Set Labels, Subtitles, and Plot Title

`--record_label` and `--record_subtitle` are repeatable and order-sensitive. Their order follows the input records unless you use `--records_table`, where labels and subtitles belong in table columns.

```bash
gbdraw linear \
  --gbk MjeNMV.gb MelaMJNV.gb \
  --record_label "Marsupenaeus japonicus endogenous nimavirus" \
  --record_label "Melicertus latisulcatus majanivirus" \
  --record_subtitle "Ginoza2017" \
  --record_subtitle "Okinawa2016" \
  --plot_title "Majanivirus comparison" \
  --plot_title_position top \
  -o majani_titles \
  -f svg
```

## 5. Tune Definition Lines

Use definition controls when you need compact or publication-specific record labels.

```bash
gbdraw linear \
  --gbk MjeNMV.gb MelaMJNV.gb \
  --record_label "MjeNMV" \
  --record_label "MelaMJNV" \
  --show_replicon \
  --hide_accession \
  --hide_length \
  --keep_definition_left_aligned \
  --definition_line_style name:weight=bold,size=18 \
  --definition_line_style 'subtitle:size=14,color=#555555' \
  -o majani_definition_lines \
  -f svg
```

`--show_replicon` adds inferred replicon labels. `--hide_accession` and `--hide_length` remove default metadata lines.

## 6. Customize Linear Track Slots

For simple ordering, use `--linear_track_order`:

```bash
gbdraw linear \
  --gbk MjeNMV.gb \
  --show_gc \
  --show_skew \
  --linear_track_order features,gc_content,gc_skew \
  -o MjeNMV_linear_track_order \
  -f svg
```

Use `--linear_track_slot` when a track needs explicit height, spacing, side, or renderer parameters:

```bash
gbdraw linear \
  --gbk MjeNMV.gb \
  --show_gc \
  --show_skew \
  --linear_track_slot features:features@side=overlay,h=60px \
  --linear_track_slot gc_content:gc_content@side=below,h=24px,spacing=8px \
  --linear_track_slot gc_skew:gc_skew@side=below,h=24px,spacing=8px \
  --linear_track_axis_index 0 \
  -o MjeNMV_linear_slots \
  -f svg
```

[< Back to the Tutorials Index](./TUTORIALS.md)
[< Back to Tutorial 6](./6_Depth_Quantitative_Tracks.md) | [Go to Tutorial 8 >](./8_Interactive_SVG_Sessions.md)

[Home](../DOCS.md) | [Installation](../INSTALL.md) | [Quickstart](../QUICKSTART.md) | [Tutorials](./TUTORIALS.md) | [Recipes](../RECIPES.md) | [CLI Reference](../CLI_Reference.md) | [Gallery](../GALLERY.md) | [FAQ](../FAQ.md) | [About](../ABOUT.md)
