[Home](../README.md) | [Installation](./INSTALL.md) | [Usage](./USAGE.md) | [Gallery](./GALLERY.md) | [RECIPES](./RECIPES.md) | [Advanced](./ADVANCED.md) | [FAQ](./FAQ.md)
# Install `gbdraw`

**Prerequisite:** Make sure you have a [conda](https://docs.conda.io/en/latest/)-compatible package manager—[mamba](https://github.com/mamba-org/mamba) ,[micromamba](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html), [miniforge](https://github.com/conda-forge/miniforge) or plain conda—already installed and on your `$PATH`. All steps below assume you run the commands in such an environment.

### Which Installation Method is Right for You?
`gbdraw` offers several installation methods to suit different needs. Use this table to choose the best option for your use case.
| Method | Ease of Use | Access to Latest Version | Performance | Customization | Ideal User / Use Case |
| ------ | ------- | ------- | ------- | ------- | ------- |
| Streamlit Web App | ★★★★★ | Stable | Medium | Low | Researchers wanting quick visualization with zero installation; educational use. | 
| Google Colab | ★★★★☆ | Stable | Medium | Medium | Users with a Google account who want to try gbdraw in a cloud environment. | 
| Bioconda | ★★★☆☆ | Stable | High | High | Bioinformaticians needing a reproducible environment; pipeline integration. | 
 | Local Build | ★☆☆☆☆ | Development | High | High | Developers wanting to test the latest features or contribute to the project. | 
 

## Installation Instructions

### gbdraw Web App
The gbdraw GUI is available as a web app, ready to use without any local installation.

**Simply visit**: [https://gbdraw.streamlit.app/](https://gbdraw.streamlit.app/)

### Colab Notebook
`gbdraw` also provides a Google Colaboratory notebook. A Google account is required.

**Open in Colab**: [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/satoshikawato/gbdraw/blob/main/gbdraw_colab.ipynb)

### Local installation
For regular use and integration into analysis pipelines, we recommend installing `gbdraw` into a dedicated conda environment.

#### Bioconda (recommended)
`gbdraw` is available on the Bioconda channel.
```bash
mamba create -n gbdraw-0.5.1 -y -c conda-forge -c bioconda gbdraw=0.5.1
```
If you also want to use the graphical interface (`gbdraw gui`), install streamlit into the same environment.
```bash
mamba create -n gbdraw-0.5.1 -c conda-forge -c bioconda gbdraw=0.5.1 streamlit # Install streamlit if you want to use GUI mode. Streamlit can also be installed later
```

#### Local build (development version)
To use the latest development version locally, clone the repository yourself using `git` and build the package locally with [conda-build](https://anaconda.org/anaconda/conda-build).
```bash
# 1. Clone the source
git clone https://github.com/satoshikawato/gbdraw.git
cd gbdraw/

# 2. Make sure conda-build is installed
mamba install -y conda-build             # or: conda install conda-build

# 3. Build the package locally
conda-build .

# 4. Create an isolated environment from the locally built package
mamba create -n gbdraw -y  -c conda-forge -c bioconda -c local gbdraw

# 5. Activate the environment
mamba activate gbdraw
```
[Home](../README.md) | [Installation](./INSTALL.md) | [Usage](./USAGE.md) | [Gallery](./GALLERY.md) | [RECIPES](./RECIPES.md) | [Advanced](./ADVANCED.md) | [FAQ](./FAQ.md)