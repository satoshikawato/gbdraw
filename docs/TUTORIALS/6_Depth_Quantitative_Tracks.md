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

In a source checkout, both files are also available under `examples/`.

Depth TSV files use the first three columns of samtools depth output: `reference`, 1-based `position`, and non-negative `depth`. A header is optional; gbdraw normalizes these columns to `reference_name`, `position`, and `depth`.

Create two compact synthetic depth files. They contain one sample every 25,000 bp, so the commands below use `--depth_window 1 --depth_step 25000`. For ordinary per-base samtools depth output, omit those two options or choose aggregation values that fit your data.

```bash
cat > MjeNMV.depth.tsv <<'EOF'
reference	position	depth
LC738868.1	1	18
LC738868.1	25001	24
LC738868.1	50001	42
LC738868.1	75001	8
LC738868.1	100001	60
LC738868.1	125001	35
LC738868.1	150001	22
LC738868.1	175001	48
LC738868.1	200001	30
LC738868.1	225001	70
LC738868.1	250001	40
LC738868.1	275001	28
LC738868.1	300001	55
EOF

cat > MelaMJNV.depth.tsv <<'EOF'
reference	position	depth
LC738874.1	1	12
LC738874.1	25001	20
LC738874.1	50001	35
LC738874.1	75001	16
LC738874.1	100001	55
LC738874.1	125001	31
LC738874.1	150001	18
LC738874.1	175001	44
LC738874.1	200001	27
LC738874.1	225001	63
LC738874.1	250001	38
LC738874.1	275001	50
EOF
```

## 2. Add One Logical Depth Track

Use `--depth` when you need one logical depth track.

```bash
gbdraw circular \
  --gbk MjeNMV.gb \
  --depth MjeNMV.depth.tsv \
  --depth_width 45 \
  --depth_window 1 \
  --depth_step 25000 \
  --show_depth_axis \
  --show_depth_ticks \
  --depth_large_tick_interval 20 \
  --depth_small_tick_interval 10 \
  -o MjeNMV_depth_circular \
  -f svg
```

This writes `MjeNMV_depth_circular.svg`. Circular mode uses `--depth_width` for the radial thickness of the depth track. Because this record is annotated as linear, gbdraw also prints the expected topology warning for this circular-only example.

![Circular MjeNMV diagram with a blue depth ring and quantitative depth ticks](../../examples/tutorial-6-depth-circular.svg)

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
  --depth_window 1 \
  --depth_step 25000 \
  --show_depth_axis \
  --show_depth_ticks \
  --share_depth_axis \
  -o majani_depth_tracks \
  -f svg
```

This writes `majani_depth_tracks.svg`. The repeated input groups keep the example setup short; replace the second `--depth_track` group with another sample's files in a real comparison.

![Two majanivirus records with blue Run A and red Run B depth tracks on shared axes](../../examples/tutorial-6-depth-tracks.svg)

Linear mode uses `--depth_height` for the vertical height of depth tracks. `--share_depth_axis` uses the same automatically determined y-axis maximum across records for each logical track. `--depth` and `--depth_track` are alternatives and cannot be used in the same command.

Use `--depth_track_height` when each logical track needs its own height. Track-specific axis overrides are also available:

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
  --depth_window 1 \
  --depth_step 25000 \
  --depth_log_scale \
  --depth_min 1 \
  --depth_max 100 \
  --show_depth_axis \
  -o tutorial-depth-log-axis \
  -f svg
```

This writes `tutorial-depth-log-axis.svg`, with a log-scaled depth axis spanning the requested 1x to 100x range.

![Linear majanivirus diagram with a blue log-scaled depth track and quantitative axis](../../examples/tutorial-depth-log-axis.svg)

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

This writes `MjeNMV_gc_percent.svg`.

![Linear MjeNMV diagram with an absolute GC percent track and quantitative axis](../../examples/tutorial-6-gc-percent.svg)

The same percent-mode options are available in circular mode, together with circular track geometry such as `--gc_content_width` and `--gc_content_radius`.

## 6. Add Another Skew Track

Custom track slots can add a second skew track with a different dinucleotide. This example keeps the standard GC skew and adds AT skew below it:

```bash
gbdraw linear \
  --gbk MjeNMV.gb \
  --show_skew \
  --linear_track_slot 'features:features@side=overlay' \
  --linear_track_slot 'gc_skew:gc_skew@side=below,h=24px,spacing=8px' \
  --linear_track_slot 'at_skew:dinucleotide_skew@side=below,h=24px,spacing=8px,nt=AT,positive_color=#deaf6e,negative_color=#7294e3' \
  --linear_track_axis_index 0 \
  -o MjeNMV_two_skew_tracks \
  -f svg
```

This writes `MjeNMV_two_skew_tracks.svg`.

![Linear MjeNMV diagram with separate GC skew and AT skew tracks](../../examples/tutorial-6-two-skew-tracks.svg)

Use the same `nt=AT` pattern in circular track tables or circular track slots when you need an additional circular skew ring.

[< Back to the Tutorials Index](./TUTORIALS.md)
[< Back to Tutorial 5](./5_Table_Driven_Inputs.md) | [Go to Tutorial 7 >](./7_Linear_Layout.md)

[Home](../DOCS.md) | [Installation](../INSTALL.md) | [Quickstart](../QUICKSTART.md) | [Tutorials](./TUTORIALS.md) | [Recipes](../RECIPES.md) | [CLI Reference](../CLI_Reference.md) | [Gallery](../GALLERY.md) | [FAQ](../FAQ.md) | [About](../ABOUT.md)
