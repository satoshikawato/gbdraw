[Home](./DOCS.md) | [Installation](./INSTALL.md) | [Quickstart](./QUICKSTART.md) | [Tutorials](./TUTORIALS/TUTORIALS.md) | [Recipes](./RECIPES.md) | [CLI Reference](./CLI_Reference.md) | **Gallery** | [FAQ](./FAQ.md) | [About](./ABOUT.md)

# Gallery

Use this page to choose a visual capability, then use the reproducible recipes below when you need the exact inputs and settings. For diagrams with feature popups, match inspection, zoom controls, downloadable sessions, and guided browser tutorials, open the [interactive gallery](https://gbdraw.app/gallery/).

## Capabilities

### Circular genome overview

Circular mode can separate strands, place labels inside and outside the ring, and add GC content or dinucleotide-skew tracks.

![Circular human mitochondrial genome with qualifier-based labels](../examples/HmmtDNA_qualifier_priority_soft_pastels.svg)

### Linear genome comparison

Linear mode can align several records and connect nucleotide or protein matches. Use precomputed BLAST-family tables when the search is performed outside gbdraw, or use the built-in protein-search modes described in the [comparison tutorials](./TUTORIALS/2_Comparative_Genomics.md).

![Linear Escherichia and Shigella comparison](../examples/Escherichia_Shigella_pair.svg)

### Feature filtering and color rules

Whitelist and blacklist tables control visible features. Default and qualifier-specific color tables keep the visual encoding explicit and reusable.

![Filtered E. coli virulence-feature diagram](../examples/O157_H7_stx_whitelist.svg)

### Quantitative tracks

Depth, GC content, and dinucleotide-skew tracks can use quantitative axes and shared scales. See [Plot read depth and other numeric tracks](./TUTORIALS/6_Depth_Quantitative_Tracks.md).

![Circular comparison of window and step settings](../examples/window_step_comparison.png)

### Layout and typography

Track placement, strand separation, label offsets, legend position, and font sizes can be adjusted independently.

![Comparison of circular definition font sizes](../examples/definition_font_size_comparison.png)

### Palettes

gbdraw includes palettes for different visual moods and contrast needs. The [palette reference](../examples/color_palette_examples.md) shows every palette in circular and linear contexts and lists its underlying colors.

![Circular palette contact sheet](../examples/palettes_combined_image_1.png)

## Reproducible recipes

The repository contains a machine-readable figure manifest in [`tools/reproduce_examples_manifest.py`](../tools/reproduce_examples_manifest.py). The reproducer resolves inputs from `tests/test_inputs/` and `examples/`, creates small support TSV files, and reports any unavailable source instead of silently substituting it.

List recipe readiness without rendering:

```bash
python tools/reproduce_examples.py --list
```

Render one recipe into `_reproduced/examples/`:

```bash
python tools/reproduce_examples.py --figure HmmtDNA_qualifier_priority_soft_pastels
```

Render all documentation figures:

```bash
python tools/reproduce_examples.py --group docs
```

To deliberately refresh tracked figures after reviewing the output, add `--output-root .`. Do not use that option merely to preview a change.

### Human mitochondrial genome labels

Question: how can CDS labels use the `gene` qualifier while other feature types keep their defaults?

```bash
python tools/reproduce_examples.py --figure HmmtDNA_qualifier_priority_soft_pastels
```

The recipe generates the qualifier-priority TSV and renders the circular figure. The browser equivalent is the [human mtDNA AT-skew tutorial](https://gbdraw.app/gallery/#HmmtDNA_ATskew).

### Escherichia and Shigella nucleotide matches

Question: how do adjacent-genome BLAST HSPs appear in a two-record or multi-record layout?

```bash
python tools/reproduce_examples.py \
  --figure Escherichia_Shigella_pair \
  --figure Escherichia_Shigella_multi
```

The manifest records query and subject order for each search preparation. For interpretation details, see [Compare genomes with precomputed BLAST-family results](./TUTORIALS/2_Comparative_Genomics.md).

### Majanivirus TBLASTX comparison

Question: how can ten viral records be compared with precomputed translated nucleotide matches and product-based colors?

```bash
python tools/reproduce_examples.py --figure majani
```

The static recipe uses adjacent-pair TBLASTX output. The [interactive nine-genome case study](https://gbdraw.app/gallery/#majanivirus_orthogroup) instead runs a protein search and assigns CDS-derived proteins to gbdraw similarity groups; those are distinct workflows.

### Additional single-genome styles

These recipes demonstrate chloroplast annotation colors, compact archaeal circular layout, and large viral genomes:

```bash
python tools/reproduce_examples.py \
  --figure NC_001879_color \
  --figure M16-5_fugaku \
  --figure Pandoravirus_salinus_forest
```

## Advanced interactive case studies

The hosted gallery keeps the larger session-based workflows separate from the static recipes:

- [Vibrio multi-record circular layout](https://gbdraw.app/gallery/#Vnig_TUMSAT-TG-2018) — six records on one canvas.
- [Aminoglycoside biosynthetic gene clusters](https://gbdraw.app/gallery/#BGC0000708-BGC0000713) — protein search, similarity groups, qualifier colors, and legend editing.
- [Hepatoplasmataceae similarity groups](https://gbdraw.app/gallery/#hepatoplasmataceae_orthogroup) and [collinear blocks](https://gbdraw.app/gallery/#hepatoplasmataceae_collinear) — the same records rendered with two different match representations.
- [WSSV nucleotide-similarity rings](https://gbdraw.app/gallery/#WSSV_genome_comparison) — an advanced session-first case study. Its displayed command documents provenance but is not labeled as directly runnable because the complete exact public input bundle is unavailable.

[Home](./DOCS.md) | [Installation](./INSTALL.md) | [Quickstart](./QUICKSTART.md) | [Tutorials](./TUTORIALS/TUTORIALS.md) | [Recipes](./RECIPES.md) | [CLI Reference](./CLI_Reference.md) | **Gallery** | [FAQ](./FAQ.md) | [About](./ABOUT.md)
