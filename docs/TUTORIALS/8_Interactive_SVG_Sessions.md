[Home](../DOCS.md) | [Installation](../INSTALL.md) | [Quickstart](../QUICKSTART.md) | [Tutorials](./TUTORIALS.md) | [Recipes](../RECIPES.md) | [CLI Reference](../CLI_Reference.md) | [Gallery](../GALLERY.md) | [FAQ](../FAQ.md) | [About](../ABOUT.md)

[< Back to the Tutorials Index](./TUTORIALS.md)
[< Back to Tutorial 7](./7_Linear_Layout.md) | [Go to Tutorial 9 >](./9_Feature_Visibility_Shapes.md)

# Tutorial 8: Interactive SVG and Session Round Trips

**Goal:** move between CLI output, standalone interactive SVG, saved GUI sessions, and regenerated diagrams.

## 1. Prepare Inputs

```bash
wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=LC738868.1&rettype=gbwithparts&retmode=text" -O MjeNMV.gb
wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=LC738874.1&rettype=gbwithparts&retmode=text" -O MelaMJNV.gb
```

## 2. Export Standalone Interactive SVG

Use canonical `interactive-svg` spelling:

```bash
gbdraw circular \
  --gbk MjeNMV.gb \
  --labels \
  -o MjeNMV_interactive \
  -f interactive-svg
```

Expected outputs:

- `MjeNMV_interactive.svg`
- `MjeNMV_interactive.interactive.svg`

Interactive SVG output does not require CairoSVG. Open the `.interactive.svg` file in a browser; some desktop SVG viewers block embedded scripts.

## 3. What the Interactive SVG Contains

Interactive SVG embeds feature metadata for rendered features. Linear comparison plots can also include pairwise match metadata, and generated orthogroup workflows can include orthogroup metadata.

```bash
gbdraw linear \
  --gbk MjeNMV.gb MelaMJNV.gb \
  --protein_blastp_mode orthogroup \
  --show_labels orthogroup_top \
  --pairwise_match_style curve \
  -o MjeNMV_MelaMJNV_interactive \
  -f interactive-svg
```

## 4. Save a GUI Session Sidecar

`--save_session` writes a GUI-loadable `.gbdraw-session.json` next to the diagram.

```bash
gbdraw circular \
  --gbk MjeNMV.gb \
  --labels \
  --save_session \
  -o MjeNMV_session_demo \
  -f interactive-svg
```

Expected session output:

- `MjeNMV_session_demo.gbdraw-session.json`

Use `--session_output` when you want an explicit session path:

```bash
gbdraw circular \
  --gbk MjeNMV.gb \
  --labels \
  --session_output roundtrip.gbdraw-session.json \
  -o MjeNMV_session_named \
  -f svg
```

## 5. Regenerate from a Session

Use `--session` to regenerate a diagram from a saved GUI session JSON. Only output and format overrides are supported in this mode.

```bash
gbdraw circular \
  --session roundtrip.gbdraw-session.json \
  -o MjeNMV_regenerated \
  -f svg
```

## 6. Open the Local Web UI

Launch the browser UI from an installed environment:

```bash
gbdraw gui
```

Load the `.gbdraw-session.json` file in the GUI to inspect or adjust the same run.

[< Back to the Tutorials Index](./TUTORIALS.md)
[< Back to Tutorial 7](./7_Linear_Layout.md) | [Go to Tutorial 9 >](./9_Feature_Visibility_Shapes.md)

[Home](../DOCS.md) | [Installation](../INSTALL.md) | [Quickstart](../QUICKSTART.md) | [Tutorials](./TUTORIALS.md) | [Recipes](../RECIPES.md) | [CLI Reference](../CLI_Reference.md) | [Gallery](../GALLERY.md) | [FAQ](../FAQ.md) | [About](../ABOUT.md)
