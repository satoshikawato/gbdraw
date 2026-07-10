[Home](../DOCS.md) | [Installation](../INSTALL.md) | [Quickstart](../QUICKSTART.md) | [Tutorials](./TUTORIALS.md) | [Recipes](../RECIPES.md) | [CLI Reference](../CLI_Reference.md) | [Gallery](../GALLERY.md) | [FAQ](../FAQ.md) | [About](../ABOUT.md)

[< Back to the Tutorials Index](./TUTORIALS.md)
[< Back to Tutorial 5](./5_Table_Driven_Inputs.md) | [Go to Tutorial 7 >](./7_Linear_Layout.md)

# Tutorial 6: Coverage Depth and Quantitative Tracks

**Goal:** add coverage/depth tracks, absolute GC percent tracks, and quantitative axes to circular and linear diagrams.

## 1. Prepare Inputs

```bash
wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=LC738868.1&rettype=gbwithparts&retmode=text" -O MjeNMV.gb
wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=LC738874.1&rettype=gbwithparts&retmode=text" -O MelaMJNV.gb
```

Depth TSV files use three tab-separated columns: `reference`, `position`, and `depth`. A header is optional; the parser normalizes the first three columns to `reference_name`, `position`, and `depth`.

Create two small depth files:

```bash
cat > MjeNMV.depth.tsv <<'EOF'
reference	position	depth
LC738868.1	1	18
LC738868.1	5000	24
LC738868.1	10000	42
LC738868.1	20000	8
LC738868.1	40000	60
EOF

cat > MelaMJNV.depth.tsv <<'EOF'
reference	position	depth
LC738874.1	1	12
LC738874.1	5000	20
LC738874.1	10000	35
LC738874.1	20000	16
LC738874.1	40000	55
EOF
```

## 2. Add One Logical Depth Track

Use `--depth` when you need one logical depth track.

```bash
gbdraw circular \
  --gbk MjeNMV.gb \
  --depth MjeNMV.depth.tsv \
  --depth_width 45 \
  --show_depth_axis \
  --show_depth_ticks \
  --depth_large_tick_interval 20 \
  --depth_small_tick_interval 10 \
  -o MjeNMV_depth_circular \
  -f svg
```

Circular mode uses `--depth_width` for the radial thickness of the depth track.

## 3. Add Multiple Logical Depth Tracks

Use repeatable `--depth_track` when you need multiple logical depth tracks. Each `--depth_track` group accepts one shared file or one file per displayed record.

```bash
gbdraw linear \
  --gbk MjeNMV.gb MelaMJNV.gb \
  --depth_track MjeNMV.depth.tsv MelaMJNV.depth.tsv \
  --depth_track MjeNMV.depth.tsv MelaMJNV.depth.tsv \
  --depth_track_label "Run A" "Run B" \
  --depth_track_color "#4E79A7" "#E15759" \
  --depth_height 36 \
  --show_depth_axis \
  --show_depth_ticks \
  --share_depth_axis \
  -o majani_depth_tracks \
  -f svg
```

Linear mode uses `--depth_height` for the vertical height of depth tracks. `--depth` and `--depth_track` are alternatives and cannot be used in the same command.

Track-specific overrides are available when each logical track needs its own axis styling:

- `--depth_track_large_tick_interval`
- `--depth_track_small_tick_interval`
- `--depth_track_tick_font_size`

## 4. Control Scaling

Log scaling is useful when a few high-depth bins would otherwise flatten the rest of the track.

```bash
gbdraw linear \
  --gbk MjeNMV.gb \
  --depth MjeNMV.depth.tsv \
  --depth_height 40 \
  --depth_log_scale \
  --depth_min 1 \
  --depth_max 100 \
  --show_depth_axis \
  -o MjeNMV_depth_log \
  -f svg
```

Use `--no_depth_log_scale` to force linear scaling when a config file or saved session enables log scaling.

## 5. Use GC Percent Mode

The default GC content track is mean-centered deviation. Use `--gc_content_mode percent` when the y-axis should show absolute GC percent.

```bash
gbdraw linear \
  --gbk MjeNMV.gb \
  --show_gc \
  --gc_content_mode percent \
  --gc_content_min_percent 25 \
  --gc_content_max_percent 75 \
  --gc_content_large_tick_interval 10 \
  --gc_content_small_tick_interval 5 \
  --show_gc_content_axis \
  --show_gc_content_ticks \
  -o MjeNMV_gc_percent \
  -f svg
```

The same percent-mode options are available in circular mode, together with circular track geometry such as `--gc_content_width` and `--gc_content_radius`.

## 6. Add Another Skew Track

Custom track slots can add a second skew track with a different dinucleotide. This example keeps the standard GC skew and adds AT skew below it:

```bash
gbdraw linear \
  --gbk MjeNMV.gb \
  --show_skew \
  --linear_track_slot features:features@side=overlay \
  --linear_track_slot gc_skew:gc_skew@side=below,h=24px,spacing=8px \
  --linear_track_slot at_skew:dinucleotide_skew@side=below,h=24px,spacing=8px,nt=AT,positive_color=#deaf6e,negative_color=#7294e3 \
  --linear_track_axis_index 0 \
  -o MjeNMV_two_skew_tracks \
  -f svg
```

Use the same `nt=AT` pattern in circular track tables or circular track slots when you need an additional circular skew ring.

[< Back to the Tutorials Index](./TUTORIALS.md)
[< Back to Tutorial 5](./5_Table_Driven_Inputs.md) | [Go to Tutorial 7 >](./7_Linear_Layout.md)

[Home](../DOCS.md) | [Installation](../INSTALL.md) | [Quickstart](../QUICKSTART.md) | [Tutorials](./TUTORIALS.md) | [Recipes](../RECIPES.md) | [CLI Reference](../CLI_Reference.md) | [Gallery](../GALLERY.md) | [FAQ](../FAQ.md) | [About](../ABOUT.md)
