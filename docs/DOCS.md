[Home](./DOCS.md) | [Installation](./INSTALL.md) | [Quickstart](./QUICKSTART.md) | [Tutorials](./TUTORIALS/TUTORIALS.md) | [Recipes](./RECIPES.md) | [CLI Reference](./CLI_Reference.md) | [Gallery](./GALLERY.md) | [FAQ](./FAQ.md) | [About](./ABOUT.md)

# gbdraw Documentation

`gbdraw` is available as a hosted web app, a local browser-based GUI, and a command-line tool for circular and linear genome diagrams.

Use [https://gbdraw.app/](https://gbdraw.app/) if you want the fastest path to a plot, or install `gbdraw` locally if you need pipeline integration, non-SVG export, or offline use.

## Start Here

- [Installation](./INSTALL.md): choose between the hosted web app, Bioconda, or a local development install.
- [Quickstart](./QUICKSTART.md): create a first circular plot in a few minutes.
- [Recipes](./RECIPES.md): copy-paste command patterns for common tasks.

## Tutorials

- [Tutorial 1: Customizing Your Plot](./TUTORIALS/1_Customizing_Plots.md)
- [Tutorial 2: Comparative Genomics](./TUTORIALS/2_Comparative_Genomics.md)
- [Tutorial 3: Advanced Customization](./TUTORIALS/3_Advanced_Customization.md)

## Reference

- [CLI Reference](./CLI_Reference.md): current command help for `gbdraw`, `gbdraw circular`, and `gbdraw linear`.
- [Gallery](./GALLERY.md): example plots and the commands used to generate them.
- [Interactive SVG Gallery](https://gbdraw.app/gallery/): JavaScript-enabled SVG examples embedded in sandboxed iframes.
- [FAQ](./FAQ.md): common questions, limitations, and workarounds.
- [About](./ABOUT.md): citation information and project background.

## Interface Guide

- Hosted web app: [https://gbdraw.app/](https://gbdraw.app/)
- Interactive SVG gallery: [https://gbdraw.app/gallery/](https://gbdraw.app/gallery/)
- Local GUI: `gbdraw gui`
- CLI entry point: `gbdraw circular ...` and `gbdraw linear ...`
- Collinear linear comparisons: use `--protein_blastp_mode collinear`; `--collinear_min_anchors` sets the minimum anchors/genes required for a rendered Collinear block. The default 1 includes singleton links. Collinear blocks use RBH anchors. In the web UI, Evidence scope controls which record pairs provide Collinear evidence: Adjacent pairs uses neighboring records and reports local collinear gene groups, while All records uses every record pair for collinearity-backed global metadata; rendered ribbons remain adjacent display pairs. Orthogroups mode remains all-vs-all global orthogroup inference. `--collinear_color_mode orientation_identity` uses separate forward and inverted identity gradients.

[Home](./DOCS.md) | [Installation](./INSTALL.md) | [Quickstart](./QUICKSTART.md) | [Tutorials](./TUTORIALS/TUTORIALS.md) | [Recipes](./RECIPES.md) | [CLI Reference](./CLI_Reference.md) | [Gallery](./GALLERY.md) | [FAQ](./FAQ.md) | [About](./ABOUT.md)
