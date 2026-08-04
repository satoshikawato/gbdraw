[Documentation home](../DOCS.md) | [How-to guides](../HOW_TO/README.md) | [Explanation](../EXPLANATION/README.md) | [Reference](../REFERENCE/README.md)

# gbdraw Tutorials

Tutorials take you through one complete result. Choose a learning goal; each
project links to focused How-to and Reference material for optional variations.

## One figure, three ways to build it

A Tutorial project is defined by its finished biological figure, not by its
interface. Each project should have GUI, CLI, and Python API variants that use
the same fixtures, record order, comparisons, feature selection, labels,
colors, track geometry, and finished SVG semantics. The variants may differ
only in the steps and syntax used to produce that figure. Surface-specific SVG
wrappers, such as browser-added outer canvas padding and serialization
metadata, are allowed when they do not change the finished diagram.

Every complete project starts with a **Choose how to build this figure** table.
GitHub Markdown cannot provide portable interactive tabs, so these tables are
the source-compatible interface switch. A future documentation site can render
the same project metadata as tabs without changing the tutorial content.

Do not fill a missing interface with a different or biologically weaker
example. If a surface cannot express the shared figure, record the missing
capability explicitly. Existing single-surface chapters are migration backlog;
new Tutorial projects must start with all three variants.

| Finished figure | GUI | CLI | Python API | Parity |
| --- | --- | --- | --- | --- |
| Human mitochondrial first circular map | [GUI](GUI/first-circular-genome-diagram.md) | [CLI](CLI/first-circular-genome-diagram.md) | [Python](PYTHON/first-genome-diagram.md) | Verified |
| Lambda first linear map | [GUI](GUI/first-linear-genome-diagram.md) | [CLI](CLI/first-linear-genome-diagram.md) | Not yet migrated | Migration |
| Lambda–DE3 LOSATN comparison | [GUI](GUI/compare-genomes-losatn.md) | Not yet migrated | Not yet migrated | Migration |
| Aminoglycoside BGC Similarity groups | [GUI](GUI/compare-proteins-losatp.md) | Not yet migrated | Not yet migrated | Migration |
| Gallery tobacco chloroplast map | [GUI](GUI/build-an-annotated-chloroplast-map.md) | [CLI](CLI/build-an-annotated-chloroplast-map.md) | [Python](PYTHON/build-an-annotated-chloroplast-map.md) | Verified |
| Human mitochondrial comparison rings | [GUI](GUI/add-precomputed-circular-comparison-rings.md) | Not yet migrated | Not yet migrated | Migration |
| Hepatoplasmataceae Collinear map | [GUI](GUI/compare-proteins-losatp-collinear.md) | Not yet migrated | Not yet migrated | Migration |
| Interactive session handoff | [GUI](GUI/create-and-resume-an-interactive-figure.md) | Not yet migrated | Not yet migrated | Migration |
| Human mitochondrial feature presentation | Not yet migrated | [CLI](CLI/highlight-mitochondrial-features.md) | Not yet migrated | Migration |
| Majanivirus table-driven comparison | Not yet migrated | [CLI](CLI/build-a-table-driven-genome-comparison.md) | Not yet migrated | Migration |
| Quantitative Hepatoplasma map | Not yet migrated | [CLI](CLI/build-a-quantitative-genome-map.md) | Not yet migrated | Migration |

## Start with one genome

- [Create and export your first circular genome diagram](GUI/first-circular-genome-diagram.md) — web app
- [Create and export your first linear genome diagram](GUI/first-linear-genome-diagram.md) — web app
- [Create a reproducible circular diagram from the command line](CLI/first-circular-genome-diagram.md)
- [Create a reproducible linear diagram from the command line](CLI/first-linear-genome-diagram.md)
- [Draw and save your first genome diagram from Python](PYTHON/first-genome-diagram.md)

## Customize a complete biological map

- [Recreate the Interactive SVG Gallery chloroplast map](GUI/build-an-annotated-chloroplast-map.md) — web app
- [Recreate the Gallery chloroplast map from the command line](CLI/build-an-annotated-chloroplast-map.md) — CLI
- [Recreate the Gallery chloroplast map from Python](PYTHON/build-an-annotated-chloroplast-map.md) — Python API
- [Highlight mitochondrial features without editing the GenBank file](CLI/highlight-mitochondrial-features.md)
- [Build a quantitative genome map with depth, GC content, and skew](CLI/build-a-quantitative-genome-map.md)

## Compare genomes and proteins

- [Compare two genomes with LOSATN in the browser](GUI/compare-genomes-losatn.md)
- [Create protein Similarity groups with LOSATP in the browser](GUI/compare-proteins-losatp.md)
- [Add circular comparison rings from precomputed results](GUI/add-precomputed-circular-comparison-rings.md)
- [Find conserved gene order from the Gallery LOSATP Collinear project](GUI/compare-proteins-losatp-collinear.md)
- [Build a biological-pair genome comparison from TSV manifests](CLI/build-a-table-driven-genome-comparison.md)

## Hand off reproducible work

- [Create an interactive figure and reproduce it from a saved session](GUI/create-and-resume-an-interactive-figure.md)

Browse the same projects by interface in the [web app](GUI/README.md),
[command-line](CLI/README.md), and [Python](PYTHON/README.md) indexes.

## Compatibility routes for the retired numbered guides

The migration to focused Tutorial, How-to, Reference, and Explanation pages is complete. The former numbered URLs remain as short compatibility routes for one renovation release:

1. [Style a circular genome diagram](1_Customizing_Plots.md)
2. [Draw genome comparison links from precomputed BLAST results](2_Comparative_Genomics.md)
3. [Set feature colors and labels](3_Advanced_Customization.md)
4. [Draw protein matches from annotated CDS features](4_Protein_Comparisons.md)
5. [Use TSV manifests for CLI inputs](5_Table_Driven_Inputs.md)
6. [Plot read depth and other numeric tracks](6_Depth_Quantitative_Tracks.md)
7. [Arrange linear tracks, record labels, and rulers](7_Linear_Layout.md)
8. [Create interactive SVGs and restore saved sessions](8_Interactive_SVG_Sessions.md)
9. [Control feature visibility and shapes](9_Feature_Visibility_Shapes.md)

Each compatibility page points to the current canonical destinations and contains no duplicate runnable procedure. Complete option inventories, schemas, defaults, and compatibility history live in [Reference](../REFERENCE/README.md).

[Documentation home](../DOCS.md) | [How-to guides](../HOW_TO/README.md) | [Explanation](../EXPLANATION/README.md) | [Reference](../REFERENCE/README.md)
