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

If you are working from a source checkout, the same files are available under `examples/`.

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

Hover over a feature for a compact summary, or click it to open the **Details**, **Qualifiers**, and **Sequence** tabs. The **+**, **-**, **Original**, and **Search** controls are another quick way to distinguish the interactive file from the static SVG.

![Standalone interactive SVG with zoom and search controls and a feature details popup open](./images/tutorial-8-interactive-feature-popup.png)

## 3. What the Interactive SVG Contains

Interactive SVG embeds feature metadata for rendered features. Linear comparison plots can also include pairwise match metadata. Create a minimal BLAST outfmt 6 file for one comparison ribbon:

```bash
cat > MjeNMV.MelaMJNV.tblastx.out <<'EOF'
LC738868.1	LC738874.1	91.2	1000	88	0	15000	16000	14800	15800	1e-80	300
EOF
```

```bash
gbdraw linear \
  --gbk MjeNMV.gb MelaMJNV.gb \
  --blast MjeNMV.MelaMJNV.tblastx.out \
  --pairwise_match_style curve \
  -o MjeNMV_MelaMJNV_interactive \
  -f interactive-svg
```

Click the comparison ribbon in `MjeNMV_MelaMJNV_interactive.interactive.svg` to inspect its match metadata. Generated orthogroup workflows also include orthogroup metadata; complete the runtime setup in [Tutorial 4](./4_Protein_Comparisons.md) before using `--protein_blastp_mode orthogroup`.

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

As in Step 2, this command also writes `MjeNMV_session_demo.svg` and `MjeNMV_session_demo.interactive.svg` because the requested format is `interactive-svg`.

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

Use `--session` with the same diagram subcommand that created the session. Alongside `--session`, the CLI accepts output and format overrides plus optional `--save_session` or `--session_output`; other diagram options are rejected.

```bash
gbdraw circular \
  --session roundtrip.gbdraw-session.json \
  -o MjeNMV_regenerated \
  -f svg
```

This writes `MjeNMV_regenerated.svg` without requiring the original GenBank file next to the session.

## 6. Open the Local Web UI

Launch the browser UI from an installed environment:

```bash
gbdraw gui
```

Click **Load Session**, then choose the `.gbdraw-session.json` file to restore its embedded inputs, settings, and saved result. Load the JSON sidecar, not the `.interactive.svg` file.

[< Back to the Tutorials Index](./TUTORIALS.md)
[< Back to Tutorial 7](./7_Linear_Layout.md) | [Go to Tutorial 9 >](./9_Feature_Visibility_Shapes.md)

[Home](../DOCS.md) | [Installation](../INSTALL.md) | [Quickstart](../QUICKSTART.md) | [Tutorials](./TUTORIALS.md) | [Recipes](../RECIPES.md) | [CLI Reference](../CLI_Reference.md) | [Gallery](../GALLERY.md) | [FAQ](../FAQ.md) | [About](../ABOUT.md)
