[![Anaconda-Server Badge](https://anaconda.org/bioconda/gbdraw/badges/version.svg)](https://anaconda.org/bioconda/gbdraw)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/gbdraw/badges/platforms.svg)](https://anaconda.org/bioconda/gbdraw)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/gbdraw/README.html)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/gbdraw/badges/latest_release_date.svg)](https://anaconda.org/bioconda/gbdraw)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/gbdraw/badges/license.svg)](https://anaconda.org/bioconda/gbdraw)
[![Ask DeepWiki](https://deepwiki.com/badge.svg)](https://deepwiki.com/satoshikawato/gbdraw)

# gbdraw
![gbdraw](https://github.com/satoshikawato/gbdraw/blob/main/examples/gbdraw_social_preview.png)
`gbdraw` is a command-line/GUI tool designed for creating detailed diagrams of microbial genomes. 
`gbdraw` accepts GenBank/DDBJ-format annotated genomes or GFF3+FASTA file pairs as input and outputs a visual representation of the genomes in SVG/PNG/PDF/EPS/PS formats.

**Try gbdraw Web App!** [https://gbdraw.streamlit.app/](https://gbdraw.streamlit.app/)

**Try gbdraw on Colab Notebook!** [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/satoshikawato/gbdraw/blob/main/gbdraw_colab.ipynb)

## Features
- Circular and linear diagrams: Generates both circular and linear representations of genome structures.
- **Multiple Input Formats:** Supports standard **GenBank/DDBJ** files as well as **GFF3 + FASTA** file pairs.
- **Dual Interface:** Available as a powerful command-line tool and an interactive **local/web GUI**.
- Various output formats: Vector and raster graphics suitable for publication and further editing.
- **Flexible Label Control:** Provides advanced control over feature labels, including priority rules, blacklists, and **whitelists**.
- Comparative genomics: Visualizes sequence similarity between genomes using BLAST results.


### Optional dependencies (GUI)
- [Streamlit](https://streamlit.io/) (for the `gbdraw gui` command)
## Use without local installation
### Streamlit Web App
A GUI web app of `gbdraw` (latest commit on `main` branch) is available on Streamlit without any local installation:
[https://gbdraw.streamlit.app/](https://gbdraw.streamlit.app/)

### Colab Notebook (Google account required)
You can try `gbdraw` (latest release) on Google Colaboratory without any local installation:

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/satoshikawato/gbdraw/blob/main/gbdraw_colab.ipynb)


## Bug reports and suggestions
Please feel free to submit a new issue if you find a bug or have a suggestion:
https://github.com/satoshikawato/gbdraw/issues





