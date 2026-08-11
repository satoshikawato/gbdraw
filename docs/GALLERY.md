[Documentation home](./DOCS.md) | [Tutorials](./TUTORIALS/README.md) | [Technical documentation](./REFERENCE/README.md) | [FAQ](./FAQ.md) | **Gallery** | [Installation](./INSTALL.md) | [About](./ABOUT.md)

# Gallery

This is a curated showcase of figures made with gbdraw. Click any figure to open the full-size output. A tracked image is not, by itself, a reproducible recipe. Entries are described as reproducible only when they link to a version-pinned fixture and an executable Tutorial.

For zooming, feature popups, match inspection, and downloadable sessions, open the [interactive gallery](https://gbdraw.app/gallery/).

## Circular genome maps

<table>
  <tr>
    <td width="50%" valign="top">
      <a href="../examples/NC_001879_regions.svg"><img src="../examples/NC_001879_regions.svg" alt="Circular Nicotiana tabacum chloroplast map with LSC, SSC, IRa, and IRb region brackets" width="100%"></a><br>
      <strong><em>Nicotiana tabacum</em> chloroplast</strong><br>
      Feature labels, a GC-content track, and LSC, SSC, IRa, and IRb region brackets. Follow the <a href="./TUTORIALS/GUI/build-an-annotated-chloroplast-map.md">annotated chloroplast Tutorial</a> or open the <a href="https://gbdraw.app/gallery/#tobacco-chloroplast">interactive example</a>.
    </td>
    <td width="50%" valign="top">
      <a href="../examples/HmmtDNA_qualifier_priority_soft_pastels.svg"><img src="../examples/HmmtDNA_qualifier_priority_soft_pastels.svg" alt="Circular human mitochondrial genome with labels placed inside and outside the feature ring" width="100%"></a><br>
      <strong>Human mitochondrial genome</strong><br>
      Qualifier-based labels placed inside and outside the ring with the <code>soft_pastels</code> palette. See the <a href="./REFERENCE/palettes-feature-rules-labels-shapes-and-tracks.md#feature-presentation">feature-presentation technical documentation</a>.
    </td>
  </tr>
  <tr>
    <td width="50%" valign="top">
      <a href="../examples/M16-5_fugaku.svg"><img src="../examples/M16-5_fugaku.svg" alt="Compact circular genome map of Candidatus Sukunaarchaeum mirabile" width="100%"></a><br>
      <strong><em>Ca.</em> Sukunaarchaeum mirabile</strong><br>
      A compact archaeal genome using separated strands, a centered feature track, and the <code>fugaku</code> palette. This input is not bundled with gbdraw; browse the <a href="https://gbdraw.app/gallery/palettes/">Circular Palette Explorer</a> to try <code>fugaku</code> on your own genome.
    </td>
    <td width="50%" valign="top">
      <a href="../examples/Pandoravirus_salinus_forest.svg"><img src="../examples/Pandoravirus_salinus_forest.svg" alt="Circular Pandoravirus salinus genome map with dense feature tracks" width="100%"></a><br>
      <strong><em>Pandoravirus salinus</em></strong><br>
      A large viral genome with dense forward- and reverse-strand annotations using the <code>forest</code> palette. This input is not bundled with gbdraw; browse the <a href="https://gbdraw.app/gallery/palettes/">Circular Palette Explorer</a> to try <code>forest</code> on your own genome.
    </td>
  </tr>
</table>

## Comparative genomics

<table>
  <tr>
    <td width="50%" valign="top">
      <a href="../examples/Escherichia_Shigella_pair.svg"><img src="../examples/Escherichia_Shigella_pair.svg" alt="Linear comparison of Escherichia coli and Shigella dysenteriae with nucleotide match ribbons" width="100%"></a><br>
      <strong><em>Escherichia coli</em> and <em>Shigella dysenteriae</em></strong><br>
      A two-record Linear diagram with nucleotide-match ribbons and separated feature strands. See the <a href="./REFERENCE/web-app.md#comparison-surfaces">web app comparison documentation</a>.
    </td>
    <td width="50%" valign="top">
      <a href="../examples/Escherichia_Shigella_multi.svg"><img src="../examples/Escherichia_Shigella_multi.svg" alt="Four-record linear comparison of Escherichia and Shigella genomes" width="100%"></a><br>
      <strong>Four-record <em>Escherichia</em>–<em>Shigella</em> comparison</strong><br>
      Adjacent nucleotide comparisons across four bacterial records on one canvas.
    </td>
  </tr>
  <tr>
    <td colspan="2" valign="top">
      <a href="https://gbdraw.app/gallery/#BGC0000708-BGC0000713">Open the interactive Gallery example</a><br>
      <strong>Five aminoglycoside biosynthetic gene clusters</strong><br>
      The complete BGC0000708, BGC0000709, BGC0000711, BGC0000712, and BGC0000713 records are compared with LOSATP <em>Similarity groups</em>. The ribbons show retained protein-search relationships used for visualization; they are not phylogenetic orthology calls. Open the <a href="https://gbdraw.app/gallery/#BGC0000708-BGC0000713">Gallery tutorial</a> to inspect the five-record result and its session.
    </td>
  </tr>
  <tr>
    <td colspan="2" valign="top">
      <a href="https://gbdraw.app/gallery/#vibrio-harveyi-group-collinear">Open the interactive Gallery example</a><br>
      <strong><em>Vibrio</em> Harveyi group multi-record collinearity</strong><br>
      Five RefSeq assemblies are arranged as one species per row, retaining all 11 chromosomes and plasmids. LOSATP searches 18 cross-record combinations between adjacent rows.<br><br>
      The blocks are colored by orientation and identity. All features are rectangular; species and strain appear as a two-line left definition. The legend sits below the diagram, which has no plot title.<br><br>
      The Gallery provides an interactive SVG with internally gzip-compressed metadata and a separate gzip-compressed Session JSON. Follow the <a href="https://gbdraw.app/gallery/#vibrio-harveyi-group-collinear">Gallery tutorial</a> to reproduce the workflow.
    </td>
  </tr>
  <tr>
    <td colspan="2" valign="top">
      <a href="../examples/majani.svg">Open the full-size SVG</a><br>
      <strong>Majanivirus comparison</strong><br>
      Ten viral records connected by translated-nucleotide matches, with product-based feature colors. For a protein-search version, open the <a href="https://gbdraw.app/gallery/#majanivirus_orthogroup">interactive case study</a>.
    </td>
  </tr>
</table>

## Labels, tracks, and visual styles

<table>
  <tr>
    <td width="50%" valign="top">
      <a href="../examples/O157_H7_stx_whitelist.svg"><img src="../examples/O157_H7_stx_whitelist.svg" alt="Circular Escherichia coli O157 H7 genome with selected virulence-feature labels" width="100%"></a><br>
      <strong>Selected virulence-feature labels</strong><br>
      A label whitelist keeps attention on selected <em>E. coli</em> O157:H7 features. See the <a href="./REFERENCE/palettes-feature-rules-labels-shapes-and-tracks.md#feature-presentation">feature-presentation technical documentation</a>.
    </td>
    <td width="50%" valign="top">
      <a href="../examples/tutorial-6-depth-circular.svg"><img src="../examples/tutorial-6-depth-circular.svg" alt="Circular bacterial genome with a blue read-depth track and quantitative tick labels" width="100%"></a><br>
      <strong>Read-depth track</strong><br>
      A circular depth profile with a quantitative axis and evenly spaced tick labels. Follow the <a href="./TUTORIALS/GUI/build-a-quantitative-genome-map.md">quantitative genome map Tutorial</a>.
    </td>
  </tr>
  <tr>
    <td width="50%" valign="top">
      <a href="../examples/tutorial-9-feature-shapes.svg"><img src="../examples/tutorial-9-feature-shapes.svg" alt="Human mitochondrial genome with rectangular CDS, rRNA, and tRNA features" width="100%"></a><br>
      <strong>Feature-shape overrides</strong><br>
      CDS, rRNA, and tRNA features rendered as rectangles instead of directional arrows. Follow the <a href="./TUTORIALS/GUI/highlight-mitochondrial-features.md">mitochondrial feature Tutorial</a>.
    </td>
    <td width="50%" valign="top">
      <a href="https://gbdraw.app/gallery/palettes/"><img src="../examples/AP027078_tuckin_separate_strands_default.svg" alt="Circular genome map in the default gbdraw color palette" width="100%"></a><br>
      <strong>Built-in color palettes</strong><br>
      Choose any built-in palette and recolor one full-size Circular SVG immediately in the <a href="https://gbdraw.app/gallery/palettes/">Circular Palette Explorer</a>. The <a href="../examples/color_palette_examples.md">palette reference</a> lists the underlying colors.
    </td>
  </tr>
</table>

## Reproduce or adapt a figure

- Use the [web app](https://gbdraw.app/) to create a diagram in the browser.
- Follow the [Quickstart](./QUICKSTART.md) to choose a checked first workflow.
- Browse the [Tutorials](./TUTORIALS/README.md) for complete workflows, [Technical documentation](./REFERENCE/README.md) for exact behavior, or [Recipes](./RECIPES.md) for command templates.

[Documentation home](./DOCS.md) | [Tutorials](./TUTORIALS/README.md) | [Technical documentation](./REFERENCE/README.md) | [FAQ](./FAQ.md) | **Gallery** | [Installation](./INSTALL.md) | [About](./ABOUT.md)
