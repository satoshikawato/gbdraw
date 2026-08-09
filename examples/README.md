# About this directory

This directory holds pre-rendered output files (SVG/PNG/PDF), input fixtures
(GenBank, FASTA, TSV), and BLAST/LOSAT result files used as image sources and
downloadable inputs by [Gallery](../docs/GALLERY.md), the
[palette reference](./color_palette_examples.md), and a small number of
Tutorial and technical documentation pages. It is a media and fixture store,
not a curated "getting started" example gallery.

## Where to find runnable examples

For copyable, currently tested example code and commands, use:

- [Tutorials](../docs/TUTORIALS/README.md): complete projects for the web app, CLI, and Python.
- [Technical documentation](../docs/REFERENCE/README.md): exact controls, options, schemas, and API contracts.
- [Recipes](../docs/RECIPES.md): short copyable CLI command templates.

The documentation test suite checks Tutorial programs and command recipes
against the current CLI, GUI, and Python API. Files in this directory are not
covered by that check, and some predate
the current CLI option names. Do not copy commands from files in this
directory; copy them from a Tutorial, Technical documentation page, or Recipe
instead.

## Provenance

Most files here are named after the source accession or fixture they were
generated from (for example, `AP027078_tuckin_separate_strands_default.svg`).
A file referenced from a current documentation page is regenerated when that
page's recipe changes. Other files remain only when a manifest or an asset
inventory still requires them; removed outputs remain available in Git history.
