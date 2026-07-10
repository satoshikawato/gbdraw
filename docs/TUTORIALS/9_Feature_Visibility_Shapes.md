[Home](../DOCS.md) | [Installation](../INSTALL.md) | [Quickstart](../QUICKSTART.md) | [Tutorials](./TUTORIALS.md) | [Recipes](../RECIPES.md) | [CLI Reference](../CLI_Reference.md) | [Gallery](../GALLERY.md) | [FAQ](../FAQ.md) | [About](../ABOUT.md)

[< Back to the Tutorials Index](./TUTORIALS.md)
[< Back to Tutorial 8](./8_Interactive_SVG_Sessions.md)

# Tutorial 9: Feature Visibility and Shapes

**Goal:** override which features are drawn or analyzed, and change feature shapes without changing the input annotation.

## 1. Prepare an Input

```bash
wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=LC738868.1&rettype=gbwithparts&retmode=text" -O MjeNMV.gb
```

## 2. Override Feature Shapes

`--feature_shape TYPE=SHAPE` is repeatable. Supported shapes are `arrow` and `rectangle`.

```bash
gbdraw circular \
  --gbk MjeNMV.gb \
  --feature_shape CDS=rectangle \
  --feature_shape repeat_region=arrow \
  --track_type middle \
  -o MjeNMV_feature_shapes \
  -f svg
```

Shape overrides apply by feature type. They do not change colors, labels, or feature selection by themselves.

## 3. Create a Feature Visibility Table

Create `feature_visibility.tsv`:

```tsv
record_id	feature_type	qualifier	value	action
LC738868.1	CDS	product	^hypothetical protein$	off
*	CDS	product	wsv.*-like protein	show
*	CDS	product	baculoviral IAP repeat-containing protein	exclude_matching
```

Columns:

- `record_id`: exact record ID, or `*` for any record
- `feature_type`: feature type such as `CDS`, or `*`
- `qualifier`: qualifier key, or special selectors `hash`, `location`, or `record_location`
- `value`: case-insensitive Python regular expression
- `action`: `show`, `off`, or `exclude_matching`

## 4. Apply Visibility Rules

```bash
gbdraw circular \
  --gbk MjeNMV.gb \
  --feature_visibility_table feature_visibility.tsv \
  --labels \
  --track_type middle \
  -o MjeNMV_visibility \
  -f svg
```

Actions:

- `show` draws matching features even when the baseline `-k/--features` list would not include them.
- `off` hides matching features and removes them from downstream analysis inputs.
- `exclude_matching` keeps the normal drawing decision, but removes matching features from downstream analysis inputs such as protein comparison grouping.

Use only the current action names listed above.

## 5. Combine with Feature Type, Color, and Label Controls

`-k/--features` sets the baseline feature types. Color tables and label tables still act on the features that remain visible.

```bash
gbdraw linear \
  --gbk MjeNMV.gb \
  -k CDS,rRNA,tRNA,repeat_region \
  --feature_visibility_table feature_visibility.tsv \
  --feature_shape CDS=rectangle \
  --show_labels all \
  --label_blacklist hypothetical \
  -o MjeNMV_visibility_linear \
  -f svg
```

For precise targeting, prefer stable selectors such as `record_id`, `feature_type`, `protein_id`, `locus_tag`, `location`, or `record_location`. Use broad product regexes only when the annotation text is consistent across records.

[< Back to the Tutorials Index](./TUTORIALS.md)
[< Back to Tutorial 8](./8_Interactive_SVG_Sessions.md)

[Home](../DOCS.md) | [Installation](../INSTALL.md) | [Quickstart](../QUICKSTART.md) | [Tutorials](./TUTORIALS.md) | [Recipes](../RECIPES.md) | [CLI Reference](../CLI_Reference.md) | [Gallery](../GALLERY.md) | [FAQ](../FAQ.md) | [About](../ABOUT.md)
