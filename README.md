[![Static Badge](https://img.shields.io/badge/gbdraw%20webapp-8A2BE2)](https://gbdraw.app/)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/gbdraw/badges/version.svg)](https://anaconda.org/bioconda/gbdraw)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/gbdraw/badges/platforms.svg)](https://anaconda.org/bioconda/gbdraw)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/gbdraw/README.html)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/gbdraw/badges/latest_release_date.svg)](https://anaconda.org/bioconda/gbdraw)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/gbdraw/badges/downloads.svg)](https://anaconda.org/bioconda/gbdraw)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/gbdraw/badges/license.svg)](https://anaconda.org/bioconda/gbdraw)
[![Ask DeepWiki](https://deepwiki.com/badge.svg)](https://deepwiki.com/satoshikawato/gbdraw)

# gbdraw

![gbdraw](https://github.com/satoshikawato/gbdraw/blob/main/examples/gbdraw_social_preview.png)

`gbdraw` is a command-line and browser-based tool for publication-quality genome diagrams. It accepts GenBank/DDBJ flatfiles or GFF3 + FASTA pairs and produces circular or linear plots in SVG, PNG, PDF, EPS, or PS.

## Features

- Circular and linear genome diagrams from annotated sequence files
- Comparative genomics tracks from BLAST outfmt 6 or 7
- Browser-based GUI at [gbdraw.app](https://gbdraw.app/) and local GUI via `gbdraw gui`
- Fine-grained control over colors, labels, legends, and track layout
- SVG-first workflow with optional CairoSVG-based export to PNG, PDF, EPS, and PS

## Documentation

| Guide | Description |
| --- | --- |
| [Full documentation](./docs/DOCS.md) | Main entry point for the published docs set. |
| [Installation](./docs/INSTALL.md) | Hosted app, Bioconda, and local development installation. |
| [Quickstart](./docs/QUICKSTART.md) | First circular plot in a few minutes. |
| [Tutorials](./docs/TUTORIALS/TUTORIALS.md) | Step-by-step guides for styling, comparisons, and advanced customization. |
| [Recipes](./docs/RECIPES.md) | Copy-paste command patterns for common tasks. |
| [CLI Reference](./docs/CLI_Reference.md) | Current command help for the CLI. |
| [Gallery](./docs/GALLERY.md) | Example plots and commands. |
| [FAQ](./docs/FAQ.md) | Common questions and known limitations. |
| [About](./docs/ABOUT.md) | Citation information and project background. |

## Use Without Local Installation

The hosted web app runs entirely in your browser:

[https://gbdraw.app/](https://gbdraw.app/)

Your data stays on your machine, and no build step or Streamlit install is required.

## Local Installation

For regular use, Bioconda is the main installation path:

```bash
mamba create -n gbdraw -c conda-forge -c bioconda gbdraw
conda activate gbdraw
gbdraw -h
gbdraw gui
```

For development from source:

```bash
git clone https://github.com/satoshikawato/gbdraw.git
cd gbdraw
python -m pip install -e ".[dev]"
```

If you need PNG/PDF/EPS/PS export from a source install, add the optional export dependency:

```bash
python -m pip install -e ".[dev,export]"
```

See [Installation](./docs/INSTALL.md) for details and platform notes.

## Bug Reports and Suggestions

Please open an issue if you find a bug or want to propose an improvement:

https://github.com/satoshikawato/gbdraw/issues
