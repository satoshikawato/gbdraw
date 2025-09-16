[Home](../README.md) | [Installation](./INSTALL.md) | [Usage](./USAGE.md) | [Gallery](./GALLERY.md) | [RECIPES](./RECIPES.md) | [Advanced](./ADVANCED.md) | [FAQ](./FAQ.md)
# Install gbdraw
## Dependencies
- [Python](https://www.python.org/) >=3.10
- [Biopython](https://biopython.org/)
- [bcbio-gff](https://github.com/chapmanb/bcbb/tree/master/gff)
- [pandas](https://pandas.pydata.org/)
- [svgwrite](https://github.com/mozman/svgwrite)
- [CairoSVG](https://cairosvg.org/)
- [Liberation Fonts](https://github.com/liberationfonts/liberation-fonts) (bundled; SIL Open Font_License 1.1)

## Local Installation
**Prerequisite:** Make sure you have a [conda](https://docs.conda.io/en/latest/)-compatible package manager—[mamba](https://github.com/mamba-org/mamba) ,[micromamba](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html), [miniforge](https://github.com/conda-forge/miniforge) or plain conda—already installed and on your `$PATH`. All steps below assume you run the commands in such an environment.
### Bioconda (recommended)
`gbdraw` is available on the Bioconda channel.
```bash
mamba create -n gbdraw-0.5.1 -y -c conda-forge -c bioconda gbdraw=0.5.1
mamba create -n gbdraw-0.5.1 -c conda-forge -c bioconda gbdraw=0.5.1 streamlit # Install streamlit if you want to use GUI mode. Streamlit can also be installed later

mamba activate gbdraw-0.5.1
```
### Local build (development version)
To use the latest development version, clone the repository yourself using `git` and build the package locally with [conda-build](https://anaconda.org/anaconda/conda-build).
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