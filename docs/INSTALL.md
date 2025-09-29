[Home](./README.md) | **Installation** | [Quickstart](./QUICKSTART.md) | [Tutorials](./TUTORIALS/TUTORIALS.md)) | [Gallery](./GALLERY.md) | [FAQ](./FAQ.md)

[< Back to gbdraw Documentation](./DOCS.md)　　　　　　[Go to Quickstart >](./QUICKSTART.md)

# `gbdraw` Installation

`gbdraw` offers several installation methods to suit different needs. Use this table to choose the best option for your use case.
| Method | Ease of Use | Access to Latest Version | Performance | Customization | Ideal User / Use Case |
| ------ | ------- | ------- | ------- | ------- | ------- |
| Streamlit Web App | ★★★★★ | Development | Medium | High | Researchers wanting quick visualization with zero installation; educational use. | 
| Google Colab | ★★★★★ | Stable | Medium | Medium | Users with a Google account who want to try gbdraw in a cloud environment. | 
| Bioconda | ★★★★☆ | Stable | High | High | Bioinformaticians needing a reproducible environment; pipeline integration. | 
 | Local Build | ★☆☆☆☆ | Development | High | High | Developers wanting to test the latest features or contribute to the project. | 
 
---

## 1. Use Without Installation

### Streamlit Web App

The `gbdraw` GUI is available as a web app, ready to use without any local installation.

**Simply visit**: [https://gbdraw.streamlit.app/](https://gbdraw.streamlit.app/)

### Google Colab

`gbdraw` also provides a Google Colaboratory notebook, although the functionality is limited. A Google account is required.
**Open in Colab**: [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/satoshikawato/gbdraw/blob/main/gbdraw_colab.ipynb)

---

## 2. Local Installation

For regular use and integration into analysis pipelines, installing `gbdraw` into a dedicated conda environment is recommended.

### Prerequisite

You must have a [conda](https://docs.conda.io/en/latest/)-compatible package manager installed, such as [mamba](https://github.com/mamba-org/mamba), miniforge, or micromamba. We recommend `mamba` or `micromamba` for its speed.

### Bioconda (Recommended)

`gbdraw` is available from the Bioconda channel. Creating a dedicated environment prevents dependency conflicts. If you also want to use the graphical interface (`gbdraw gui`), install streamlit into the same environment.

1.  **Create and activate a new conda environment for `gbdraw`.**
    * **For both CLI and GUI:**
        ```bash
        mamba create -n gbdraw -c conda-forge -c bioconda gbdraw streamlit
        conda activate gbdraw
        ```
    * **For CLI only:**
        ```bash
        mamba create -n gbdraw -c conda-forge -c bioconda gbdraw
        conda activate gbdraw
        ```
2.  **Verify the installation.**
    ```bash
    gbdraw -h
    ```

Installation is complete!

### Build from Source (For Developers)

To use the latest development version locally, clone the repository yourself using `git` and build the package locally with [conda-build](https://anaconda.org/anaconda/conda-build).

```bash
# 1. Clone the source code
git clone https://github.com/satoshikawato/gbdraw.git

cd gbdraw/

# 2. Install conda-build
mamba install -y conda-build

# 3. Build the package locally
conda build .

# 4. Create an environment from the locally built package
mamba create -n gbdraw-dev -y -c conda-forge -c bioconda -c local gbdraw

# 5. Activate the environment
conda activate gbdraw-dev
```

[< Back to gbdraw Documentation](./DOCS.md)　　　　　　[Go to Quickstart >](./QUICKSTART.md)


[Home](./README.md) | **Installation** | [Quickstart](./QUICKSTART.md) | [Tutorials](./TUTORIALS/TUTORIALS.md)) | [Gallery](./GALLERY.md) | [FAQ](./FAQ.md)