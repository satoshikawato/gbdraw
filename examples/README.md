# About this directory

This directory holds pre-rendered output files (SVG/PNG/PDF), input fixtures
(GenBank, FASTA, TSV), and BLAST/LOSAT result files used as image sources and
downloadable inputs by [Gallery](../docs/GALLERY.md), the
[palette reference](./color_palette_examples.md), and a small number of
How-to and Tutorial pages. It is a media and fixture store, not a curated
"getting started" example gallery.

## Where to find runnable examples

For copyable, currently tested example code and commands, use:

- [Tutorials](../docs/TUTORIALS/README.md) — complete first-diagram walkthroughs for the web app, CLI, and Python.
- [How-to guides](../docs/HOW_TO/README.md) — one focused task per page, with the exact working CLI/GUI/Python steps.
- [Recipes](../docs/RECIPES.md) — short copyable CLI command templates.

Every command and code block in those pages is checked by the test suite
(`pytest tests/ -k documentation`) against the current CLI, GUI, and Python
API. Files in this directory are not covered by that check, and some predate
the current CLI option names. Do not copy commands from files in this
directory; copy them from a Tutorial, How-to guide, or Recipe instead.

## Provenance

Most files here are named after the source accession or fixture they were
generated from (for example, `AP027078_tuckin_separate_strands_default.svg`).
A file referenced from a current documentation page is regenerated when that
page's recipe changes; a file not referenced from any current page is
historical output kept for the commit history and may not reflect the
current renderer.
