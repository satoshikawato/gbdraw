[Home](../DOCS.md) | [Installation](../INSTALL.md) | [Quickstart](../QUICKSTART.md) | [Tutorials](./TUTORIALS.md) | [Recipes](../RECIPES.md) | [CLI Reference](../CLI_Reference.md) | [Gallery](../GALLERY.md) | [FAQ](../FAQ.md) | [About](../ABOUT.md)

[< Back to the guide index](./TUTORIALS.md)
[< Previous: Arrange linear tracks and labels](./7_Linear_Layout.md) | [Next: Control feature visibility and shapes >](./9_Feature_Visibility_Shapes.md)

# Create interactive SVGs and restore saved sessions

Write interactive SVG files, save GUI sessions, and regenerate diagrams from saved sessions.

## 1. Prepare inputs

```bash
wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=LC738868.1&rettype=gbwithparts&retmode=text" -O MjeNMV.gb
wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=LC738874.1&rettype=gbwithparts&retmode=text" -O MelaMJNV.gb
```

If you are working from a source checkout, the same files are available under `examples/`.

## 2. Export standalone interactive SVG

The format name is `interactive_svg`:

```bash
gbdraw circular \
  --gbk MjeNMV.gb \
  --labels \
  -o MjeNMV_interactive \
  -f interactive_svg
```

Expected outputs:

- `MjeNMV_interactive.svg`
- `MjeNMV_interactive.interactive.svg`

Interactive SVG output does not require CairoSVG. Open the `.interactive.svg` file in a browser; some desktop SVG viewers block embedded scripts.

Hover over a feature for a compact summary, or click it to open the **Details**, **Qualifiers**, and **Sequence** tabs. The **+**, **-**, **Original**, and **Search** controls are another quick way to distinguish the interactive file from the static SVG.

![Standalone interactive SVG with zoom and search controls and a feature details popup open](./images/tutorial-8-interactive-feature-popup.png)

## 3. What the interactive SVG contains

Interactive SVG embeds feature metadata for rendered features. Linear comparison plots can also include pairwise match metadata. Reuse the precomputed BLAST outfmt 7 table maintained in `examples/`:

New exports use compact metadata schema v2. The browser reconstructs FASTA text and match-popup rows when needed, so qualifier display, sequence copy/download, search, and match popups remain available without storing duplicate pre-rendered values. Existing files that use schema v1 remain compatible with the embedded runtime and do not need conversion.

Click a Linear pairwise ribbon or collinear block, or a Circular Homology-ring HSP, to open the shared match popup. The sequence section can copy or download the query/reference span, the subject/comparison span, or both spans as nucleotide FASTA. Collinear export uses the complete block envelopes, which may include intergenic bases and genes that are not anchors. Feature and orthogroup-member actions remain separate.

These exports are ungapped genomic spans, not reconstructed alignments. Coordinates stay 1-based and inclusive in the FASTA header. A reversed coordinate pair is sliced and reverse-complemented. Circular comparison spans require a companion FASTA supplied through `--conservation_fasta`, the `comparison_fasta` table column, or the web app's **Comparison FASTA (optional)** control. Web sessions retain each optional companion file with its Homology source.

```bash
cp examples/MjeNMV.MelaMJNV.tblastx.out .
```

You can also download [the same existing BLAST table](../../examples/MjeNMV.MelaMJNV.tblastx.out) directly.

```bash
gbdraw linear \
  --gbk MjeNMV.gb MelaMJNV.gb \
  --blast MjeNMV.MelaMJNV.tblastx.out \
  --identity 98 \
  --alignment_length 1400 \
  --pairwise_match_style curve \
  -o MjeNMV_MelaMJNV_interactive \
  -f interactive_svg
```

The thresholds leave one retained ribbon. Click it in `MjeNMV_MelaMJNV_interactive.interactive.svg` to inspect its match metadata:

![Standalone interactive SVG with a precomputed TBLASTX ribbon and its pairwise match metadata popup](./images/tutorial-8-interactive-match-popup.png)

Comparisons generated with gbdraw's `orthogroup` mode also include gbdraw similarity-group metadata. Complete the runtime setup in [Draw protein matches from annotated CDS features](./4_Protein_Comparisons.md) before using `--protein_blastp_mode orthogroup`.

## 4. Save a GUI session sidecar

`--save_session` writes a GUI-loadable `.gbdraw-session.json` next to the diagram.

```bash
gbdraw circular \
  --gbk MjeNMV.gb \
  --labels \
  --save_session \
  -o MjeNMV_session_demo \
  -f interactive_svg
```

Expected session output:

- `MjeNMV_session_demo.gbdraw-session.json`

As in Step 2, this command also writes `MjeNMV_session_demo.svg` and `MjeNMV_session_demo.interactive.svg` because the requested format is `interactive_svg`.

Use `--session_output` when you want an explicit session path:

```bash
gbdraw circular \
  --gbk MjeNMV.gb \
  --labels \
  --session_output roundtrip.gbdraw-session.json \
  -o MjeNMV_session_named \
  -f svg
```

Add a `.gz` suffix, such as `roundtrip.gbdraw-session.json.gz`, to write the same Session JSON with lossless gzip compression. `--save_session` keeps the uncompressed `.gbdraw-session.json` default for compatibility.

## 5. Regenerate from a session

Use `--session` with the same diagram subcommand that created the session. It accepts both `.gbdraw-session.json` and `.gbdraw-session.json.gz`. Alongside `--session`, the CLI accepts output and format overrides plus optional `--save_session` or `--session_output`; other diagram options are rejected.

```bash
gbdraw circular \
  --session roundtrip.gbdraw-session.json \
  -o MjeNMV_regenerated \
  -f svg
```

This writes `MjeNMV_regenerated.svg` without requiring the original GenBank file next to the session.

## 6. Open the local web UI

Launch the browser UI from an installed environment:

```bash
gbdraw gui
```

Click **Load Session**, then choose the `.gbdraw-session.json` file to restore its embedded inputs, settings, and saved result. Load the JSON sidecar, not the `.interactive.svg` file.

The web app's **Save Session** action downloads a lossless gzip-compressed `.gbdraw-session.json.gz` file. **Load Session** accepts both this compressed form and the uncompressed `.gbdraw-session.json` files written by the CLI.

Current Python and Web writers use session version 35 with canonical `renderRequest` schema 3. Readers accept versions 27 through 35; public typed conversion starts at version 31, while versions 27 through 30 remain CLI replay inputs.

For Linear protein comparisons, version 35 stores stable protein identities in a schema-1 manifest, current protein raw cache entries as schema 3, and derived comparison payloads as schema 2. Nucleotide raw cache entries remain schema 2. Older protein cache entries are isolated as legacy candidates instead of being treated as current hits. Saving immediately after loading an older session preserves those candidates, even before **Generate Diagram**; generation promotes only candidates that can be verified against the restored proteins and search settings.

![Local gbdraw web app after loading a session, with the embedded GenBank input, circular settings, and saved result restored](./images/tutorial-8-loaded-session.png)

[< Back to the guide index](./TUTORIALS.md)
[< Previous: Arrange linear tracks and labels](./7_Linear_Layout.md) | [Next: Control feature visibility and shapes >](./9_Feature_Visibility_Shapes.md)

[Home](../DOCS.md) | [Installation](../INSTALL.md) | [Quickstart](../QUICKSTART.md) | [Tutorials](./TUTORIALS.md) | [Recipes](../RECIPES.md) | [CLI Reference](../CLI_Reference.md) | [Gallery](../GALLERY.md) | [FAQ](../FAQ.md) | [About](../ABOUT.md)
