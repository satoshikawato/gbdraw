[Home](../DOCS.md) | [Installation](../INSTALL.md) | [Quickstart](../QUICKSTART.md) | [Tutorials](./TUTORIALS.md) | [Recipes](../RECIPES.md) | [CLI Reference](../CLI_Reference.md) | [Gallery](../GALLERY.md) | [FAQ](../FAQ.md) | [About](../ABOUT.md)

[< Back to the Tutorials Index](./TUTORIALS.md)
[< Back to Tutorial 2](./2_Comparative_Genomics.md)

# Tutorial 3: Advanced Customization

**Goal:** use tables and styling options for fine-grained control over colors, labels, and plot appearance.

## Part 1: Advanced Color Control

`gbdraw` supports two complementary table formats:

| Method | Option | Use case |
| --- | --- | --- |
| Default override table | `-d` | Replace the default color for an entire feature type |
| Feature-specific table | `-t` | Highlight selected features based on qualifier matches |

### 1. Override Default Colors

Create `modified_default_colors.tsv`:

```tsv
CDS	#d3d3d3
```

```bash
wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=LC738868.1&rettype=gbwithparts&retmode=text" -O MjeNMV.gbk
```

```bash
gbdraw circular \
  --gbk MjeNMV.gbk \
  --separate_strands \
  --track_type middle \
  -f svg \
  --block_stroke_width 1 \
  -d modified_default_colors.tsv \
  -o MjeNMV_modified_default_colors
```

![MjeNMV_modified_default_colors.svg](../../examples/MjeNMV_modified_default_colors.svg)

### 2. Highlight Specific Features

Create `feature_specific_colors.tsv`:

```tsv
CDS	product	wsv.*-like protein	#47b8f8	WSSV-like proteins
CDS	product	baculoviral IAP repeat-containing protein	yellow	BIRP
CDS	product	tyrosine recombinase	red	tyrosine recombinase
```

Then combine `-d` and `-t`:

```bash
gbdraw circular \
  --gbk MjeNMV.gbk \
  --separate_strands \
  --track_type middle \
  -f svg \
  --block_stroke_width 1 \
  -d modified_default_colors.tsv \
  -t feature_specific_colors.tsv \
  --labels \
  -o MjeNMV_feature_specific_colors_with_labels
```

![MjeNMV_feature_specifc_colors_with_labels.svg](../../examples/MjeNMV_feature_specifc_colors_with_labels.svg)

## Part 2: Advanced Label Control

### 1. Blacklist Uninformative Labels

```bash
gbdraw circular \
  --gbk MjeNMV.gbk \
  --separate_strands \
  --track_type middle \
  -f svg \
  --block_stroke_width 1 \
  -d modified_default_colors.tsv \
  -t feature_specific_colors.tsv \
  --labels \
  --label_blacklist "hypothetical" \
  -o MjeNMV_feature_specific_colors_with_labels_blacklist
```

![MjeNMV_feature_specifc_colors_with_labels_blacklist.svg](../../examples/MjeNMV_feature_specifc_colors_with_labels_blacklist.svg)

You can also pass a file to `--label_blacklist`, with one term per line.

### 2. Whitelist Only the Labels You Need

Create `stx_whitelist.tsv`:

```tsv
CDS	gene	^stx1A$
CDS	gene	^stx1B$
CDS	gene	^stx2A$
CDS	gene	^stx2B$
```

Whitelist patterns use the same case-insensitive Python regex `search(...)` matching as `-t`.
Use `^...$` for exact matches, or broader patterns such as `wsv.*-like protein` when you want to keep a label family.

```bash
gbdraw circular \
  --gbk O157_H7.gbk \
  -o O157_H7_stx_whitelist \
  --labels \
  --separate_strands \
  --species "<i>Escherichia coli</i> O157:H7" \
  --strain "Sakai" \
  --label_whitelist stx_whitelist.tsv \
  --label_font_size 16 \
  -f svg
```

![O157_H7_stx_whitelist.svg](../../examples/O157_H7_stx_whitelist.svg)

### 3. Change Which Qualifier Supplies the Label

Create `qualifier_priority.tsv`:

```tsv
CDS	gene
```

```bash
gbdraw circular \
  --gbk HmmtDNA.gbk \
  -f svg \
  --track_type middle \
  --species "<i>Homo sapiens</i>" \
  --block_stroke_width 2 \
  --axis_stroke_width 5 \
  --labels both \
  --palette soft_pastels \
  --definition_font_size 28 \
  --label_font_size 18 \
  --qualifier_priority qualifier_priority.tsv \
  -o HmmtDNA_qualifier_priority_soft_pastels
```

![HmmtDNA_qualifier_priority_soft_pastels.svg](../../examples/HmmtDNA_qualifier_priority_soft_pastels.svg)

### 4. Override Label Text After Filtering

`--label_table` runs after whitelist, blacklist, and qualifier-priority processing.

Example table:

```tsv
NC_010162.1	CDS	label	^ATP synthase subunit alpha$	ATP synthase alpha
*	*	label	^hypothetical protein$	HP
```

Circular mode:

```bash
gbdraw circular \
  --gbk MjeNMV.gbk \
  --labels \
  --label_table label_override.tsv \
  -o MjeNMV_label_override
```

Linear mode:

```bash
gbdraw linear \
  --gbk MjeNMV.gbk \
  --show_labels all \
  --label_table label_override.tsv \
  -o MjeNMV_linear_label_override
```

## Part 3: Fine-Tune Plot Appearance

Use the current stroke and text options for publication-oriented styling:

- `--block_stroke_color` and `--block_stroke_width`
- `--axis_stroke_color` and `--axis_stroke_width`
- `--line_stroke_color` and `--line_stroke_width`
- `--definition_font_size` and `--label_font_size`
- `--outer_label_x_radius_offset`, `--outer_label_y_radius_offset`
- `--inner_label_x_radius_offset`, `--inner_label_y_radius_offset`

Example:

```bash
gbdraw circular \
  --gbk MjeNMV.gbk \
  --separate_strands \
  --track_type middle \
  --labels both \
  --block_stroke_color black \
  --block_stroke_width 1 \
  --axis_stroke_width 5 \
  --line_stroke_color gray \
  --line_stroke_width 1 \
  --definition_font_size 24 \
  --label_font_size 12 \
  --outer_label_x_radius_offset 1.05 \
  --outer_label_y_radius_offset 1.05 \
  -o MjeNMV_tuned \
  -f svg
```

Reference images:

- ![definition_font_size_comparison.png](../../examples/definition_font_size_comparison.png)
- ![label_font_size_comparison.png](../../examples/label_font_size_comparison.png)
- ![outer_label_offset_comparison.png](../../examples/outer_label_offset_comparison.png)

## Part 4: Linear Mode Input Selectors

For linear plots, the main selectors are:

- `--record_id`
- `--reverse_complement`
- `--region`

See [Tutorial 2](./2_Comparative_Genomics.md) and the [CLI Reference](../CLI_Reference.md) for the full syntax.

[< Back to the Tutorials Index](./TUTORIALS.md)
[< Back to Tutorial 2](./2_Comparative_Genomics.md)

[Home](../DOCS.md) | [Installation](../INSTALL.md) | [Quickstart](../QUICKSTART.md) | [Tutorials](./TUTORIALS.md) | [Recipes](../RECIPES.md) | [CLI Reference](../CLI_Reference.md) | [Gallery](../GALLERY.md) | [FAQ](../FAQ.md) | [About](../ABOUT.md)
