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
`gbdraw` is a command-line/GUI tool designed for creating detailed diagrams of microbial genomes. 
`gbdraw` accepts GenBank/DDBJ-format annotated genomes or GFF3+FASTA file pairs as input and outputs a visual representation of the genomes in SVG/PNG/PDF/EPS/PS formats.

**Try gbdraw Web App!** [![Static Badge](https://img.shields.io/badge/gbdraw%20webapp-8A2BE2)](https://gbdraw.app/)

## Features
- **Circular and linear diagrams:** Generates both circular and linear representations of genome structures.
- **Multiple Input Formats:** Supports standard GenBank/DDBJ files as well as GFF3 + FASTA file pairs.
- **Dual Interface:** Available as a powerful command-line tool and an interactive local/web GUI.
- **Various output formats:** Vector and raster graphics suitable for publication and further editing.
- **Flexible Label Control:** Provides advanced control over feature labels, including priority rules, blacklists, and whitelists.
- **Comparative genomics:** Visualizes sequence similarity between genomes using BLAST results.
- **Accurate text metrics:** CLI uses kerning-aware font measurements for improved label placement and legend sizing.

## 📖 Documentation & Tutorials

For detailed instructions on how to use `gbdraw`, please visit our main documentation page. You'll find step-by-step guides covering everything from installation to advanced customization.

**[➡️ Go to Full Documentation](./docs/DOCS.md)**

| Section | Description |
| :--- | :--- |
| **[Installation](./docs/INSTALL.md)** | Get started instantly with the web app or a full local installation. |
| **[Quickstart](./docs/QUICKSTART.md)** | Create your first plot in under 5 minutes. |
| **[Tutorials](./docs/TUTORIALS/TUTORIALS.md)** | Learn everything from basic customization to comparative genomics. |
| **[Gallery](./docs/GALLERY.md)** | See example plots and the commands used to create them. |
| **[FAQ](./docs/FAQ.md)** | Find answers to common questions. |

## Use without local installation
### Web App
A GUI web app of `gbdraw` (latest commit on `main` branch) is available without any local installation:

[![Static Badge](https://img.shields.io/badge/gbdraw%20webapp-8A2BE2)](https://gbdraw.app/)

The web app runs entirely in your browser.

### Local Web UI
If you have `gbdraw` installed and prefer to run the GUI locally, use:

```bash
gbdraw gui
```

This starts a local HTTP server on a free port and opens the app in your browser.

## Bug reports and suggestions
Please feel free to submit a new issue if you find a bug or have a suggestion:
https://github.com/satoshikawato/gbdraw/issues







