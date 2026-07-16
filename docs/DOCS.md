[Home](./DOCS.md) | [Installation](./INSTALL.md) | [Quickstart](./QUICKSTART.md) | [Tutorials](./TUTORIALS/TUTORIALS.md) | [Recipes](./RECIPES.md) | [CLI Reference](./CLI_Reference.md) | [Gallery](./GALLERY.md) | [FAQ](./FAQ.md) | [About](./ABOUT.md)

# gbdraw documentation

`gbdraw` is available as a hosted web app, a local browser-based GUI, and a command-line tool for circular and linear genome diagrams.

Use [https://gbdraw.app/](https://gbdraw.app/) to create a diagram without installing `gbdraw`, or install it locally for pipeline integration, non-SVG export, or offline use.

## Start here

- [Installation](./INSTALL.md): choose between the hosted web app, Bioconda, or a local development install.
- [Quickstart](./QUICKSTART.md): draw a circular genome diagram from a GenBank record.
- [Workflow guide](./WORKFLOW_GUIDE.md): choose the shortest path for one genome, comparisons, quantitative tracks, automation, or the web app.
- [Recipes](./RECIPES.md): copy-paste command patterns for common tasks.

## Tutorials

- [Style a circular genome diagram](./TUTORIALS/1_Customizing_Plots.md)
- [Draw genome comparison links from precomputed BLAST results](./TUTORIALS/2_Comparative_Genomics.md)
- [Set feature colors and labels](./TUTORIALS/3_Advanced_Customization.md)
- [More command-line guides](./TUTORIALS/TUTORIALS.md): CDS protein matches, TSV manifests, read-depth tracks, linear layout, interactive SVG sessions, and feature visibility.

## Reference

- [CLI Reference](./CLI_Reference.md): current command help for `gbdraw`, `gbdraw circular`, and `gbdraw linear`.
- [Python API](./PYTHON_API.md): render circular and linear diagrams from Python.
- [0.14.0b0 release notes](./RELEASE_NOTES_0.14.0b0.md): Python API additions, behavior corrections, and migration guidance.
- [GFF3 + FASTA input](./GFF3_FASTA.md): prepare, validate, and draw annotation files without GenBank conversion.
- [Export formats](./EXPORT.md): choose SVG, interactive SVG, PNG, PDF, EPS, or PS and understand paired interactive output files.
- [Gallery](./GALLERY.md): example diagrams and the commands used to generate them.
- [Interactive SVG gallery](https://gbdraw.app/gallery/): hosted genome diagrams with feature and match popups. Selected examples include web-app tutorials.
- [FAQ](./FAQ.md): common questions, limitations, and workarounds.
- [About](./ABOUT.md): citation information and project background.

## Interface guide

- Hosted web app: [https://gbdraw.app/](https://gbdraw.app/)
- Hosted interactive SVG gallery and web tutorials: [https://gbdraw.app/gallery/](https://gbdraw.app/gallery/)
- Local GUI: `gbdraw gui`
- CLI entry point: `gbdraw circular ...` and `gbdraw linear ...`
- Protein-search modes: `--protein_blastp_mode pairwise` draws adjacent CDS protein matches, `orthogroup` assigns CDS-derived proteins to gbdraw similarity groups, and `collinear` combines compatible runs of protein-match anchors into collinear blocks. gbdraw uses LOSAT when available and can fall back to NCBI BLAST+ `blastp`. The similarity groups are not phylogeny-based orthogroups; see [the protein-comparison guide](./TUTORIALS/4_Protein_Comparisons.md#4-gbdraw-similarity-group-ribbons-orthogroup-mode).

[Home](./DOCS.md) | [Installation](./INSTALL.md) | [Quickstart](./QUICKSTART.md) | [Tutorials](./TUTORIALS/TUTORIALS.md) | [Recipes](./RECIPES.md) | [CLI Reference](./CLI_Reference.md) | [Gallery](./GALLERY.md) | [FAQ](./FAQ.md) | [About](./ABOUT.md)
