[Home](./DOCS.md) | [Installation](./INSTALL.md) | [Quickstart](./QUICKSTART.md) | [Tutorials](./TUTORIALS/TUTORIALS.md) | [Recipes](./RECIPES.md) | [CLI Reference](./CLI_Reference.md) | [Gallery](./GALLERY.md) | [FAQ](./FAQ.md) | [About](./ABOUT.md)

# gbdraw documentation

`gbdraw` is available as a hosted web app, a local browser-based GUI, and a command-line tool for circular and linear genome diagrams.

Use [https://gbdraw.app/](https://gbdraw.app/) to create a plot without installing `gbdraw`, or install it locally for pipeline integration, non-SVG export, or offline use.

## Start here

- [Installation](./INSTALL.md): choose between the hosted web app, Bioconda, or a local development install.
- [Quickstart](./QUICKSTART.md): create a first circular plot in a few minutes.
- [Recipes](./RECIPES.md): copy-paste command patterns for common tasks.

## Tutorials

- [Tutorial 1: Customize a circular plot](./TUTORIALS/1_Customizing_Plots.md)
- [Tutorial 2: Compare genomes with BLAST and protein searches](./TUTORIALS/2_Comparative_Genomics.md)
- [Tutorial 3: Control colors, labels, and appearance](./TUTORIALS/3_Advanced_Customization.md)
- [More workflows (Tutorials 4-9)](./TUTORIALS/TUTORIALS.md#more-workflows): protein comparisons, TSV manifests, coverage/depth tracks, linear layout, interactive SVG sessions, and feature visibility.

## Reference

- [CLI Reference](./CLI_Reference.md): current command help for `gbdraw`, `gbdraw circular`, and `gbdraw linear`.
- [Gallery](./GALLERY.md): example plots and the commands used to generate them.
- [Interactive SVG gallery](https://gbdraw.app/gallery/): hosted JavaScript-enabled SVG examples embedded in sandboxed iframes, with tutorial tabs for selected Hepatoplasmataceae, aminoglycoside BGC, WSSV conservation, majanivirus, Vibrio multi-record, and human mitochondrial entries.
- [FAQ](./FAQ.md): common questions, limitations, and workarounds.
- [About](./ABOUT.md): citation information and project background.

## Interface guide

- Hosted web app: [https://gbdraw.app/](https://gbdraw.app/)
- Hosted interactive SVG gallery and web tutorials: [https://gbdraw.app/gallery/](https://gbdraw.app/gallery/)
- Local GUI: `gbdraw gui`
- CLI entry point: `gbdraw circular ...` and `gbdraw linear ...`
- Protein BLASTP comparisons: `--protein_blastp_mode pairwise` draws adjacent CDS protein matches, `orthogroup` groups CDS proteins across the input records before drawing links, and `collinear` draws gene-order-preserving blocks from compatible orthogroup anchors. gbdraw uses LOSAT when available and can fall back to NCBI BLAST+ `blastp`; see [Tutorial 2](./TUTORIALS/2_Comparative_Genomics.md#5-interpreting-orthogroups-and-collinear-blocks) for the biological meaning of orthogroups and collinear blocks.

[Home](./DOCS.md) | [Installation](./INSTALL.md) | [Quickstart](./QUICKSTART.md) | [Tutorials](./TUTORIALS/TUTORIALS.md) | [Recipes](./RECIPES.md) | [CLI Reference](./CLI_Reference.md) | [Gallery](./GALLERY.md) | [FAQ](./FAQ.md) | [About](./ABOUT.md)
