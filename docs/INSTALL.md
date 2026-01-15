[Home](./DOCS.md) | **Installation** | [Quickstart](./QUICKSTART.md) | [Tutorials](./TUTORIALS/TUTORIALS.md) | [Gallery](./GALLERY.md) | [FAQ](./FAQ.md) | [ABOUT](./ABOUT.md)

[< Back to gbdraw Documentation](./DOCS.md)　　　　　　[Go to Quickstart >](./QUICKSTART.md)

# `gbdraw` Installation

`gbdraw` offers several installation methods to suit different needs. Use this table to choose the best option for your use case.
| Method | Ease of Use | Access to Latest Version | Performance | Customization | Ideal User / Use Case |
| ------ | ------- | ------- | ------- | ------- | ------- |
| Web App (gbdraw.app) | ★★★★★ | Development | Medium | High | Researchers wanting quick visualization with zero installation; educational use. | 
| Bioconda | ★★★★☆ | Stable | High | High | Bioinformaticians needing a reproducible environment; pipeline integration. | 
| Local Build | ★☆☆☆☆ | Development | High | High | Developers wanting to test the latest features or contribute to the project. | 
 
---

## 1. Use Without Installation

### Web App (gbdraw.app)

The `gbdraw` GUI is available as a browser-based web app, ready to use without any local installation.

**Simply visit**: [https://gbdraw.app/](https://gbdraw.app/)

The legacy Streamlit deployment now redirects to the new web app. If you have a local install, you can also run the same web UI with `gbdraw gui`.

---

## 2. Local Installation

For regular use and integration into analysis pipelines, installing `gbdraw` into a dedicated conda environment is recommended.

### Prerequisite

You must have a [conda](https://docs.conda.io/en/latest/)-compatible package manager installed, such as [mamba](https://github.com/mamba-org/mamba), miniforge, or micromamba. We recommend `mamba` or `micromamba` for its speed.

### Bioconda (Recommended)

`gbdraw` is available from the Bioconda channel. Creating a dedicated environment prevents dependency conflicts. The graphical interface (`gbdraw gui`) runs as a local web app and does not require Streamlit.

1.  **Create and activate a new conda environment for `gbdraw`.**
    ```bash
    mamba create -n gbdraw -c conda-forge -c bioconda gbdraw
    conda activate gbdraw
    ```
2.  **(Optional) Launch the local web UI.**
    ```bash
    gbdraw gui
    ```
3.  **Verify the installation.**
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


[Home](./DOCS.md) | **Installation** | [Quickstart](./QUICKSTART.md) | [Tutorials](./TUTORIALS/TUTORIALS.md) | [Gallery](./GALLERY.md) | [FAQ](./FAQ.md) | [ABOUT](./ABOUT.md)
