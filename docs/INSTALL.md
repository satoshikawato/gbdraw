[Home](./DOCS.md) | **Installation** | [Quickstart](./QUICKSTART.md) | [Tutorials](./TUTORIALS/TUTORIALS.md) | [Recipes](./RECIPES.md) | [CLI Reference](./CLI_Reference.md) | [Gallery](./GALLERY.md) | [FAQ](./FAQ.md) | [About](./ABOUT.md)

[< Back to Home](./DOCS.md) | [Go to Quickstart >](./QUICKSTART.md)

# Installation

`gbdraw` supports three common ways of working:

| Method | Best for | Notes |
| --- | --- | --- |
| Hosted web app | Quick plotting with zero install | Runs at [gbdraw.app](https://gbdraw.app/) in your browser. |
| Bioconda | Daily local use and reproducible environments | Recommended for most users. |
| Source install | Development and testing current repository state | Uses `pip install -e ".[dev]"`. |

## 1. Hosted Web App

If you want the fastest path to a figure, open:

[https://gbdraw.app/](https://gbdraw.app/)

The hosted app is served as a static site on GitHub Pages behind the `gbdraw.app`
custom domain. The same `gbdraw/web/` directory is also used by `gbdraw gui` for
local browser launches.

Pyodide and the main browser-side assets used by the hosted app are vendored and
self-hosted from the repository, so the web UI does not need to fetch those runtime
dependencies from third-party CDNs.

The same UI can also be launched locally after installation with:

```bash
gbdraw gui
```

## 2. Bioconda Installation

Bioconda is the recommended local installation path for normal use.

```bash
mamba create -n gbdraw -c conda-forge -c bioconda gbdraw
conda activate gbdraw
gbdraw -h
```

Optional: launch the local GUI:

```bash
gbdraw gui
```

## 3. Source Installation for Development

Use a source install when you want the current repository state, need to run tests, or plan to contribute.

```bash
git clone https://github.com/satoshikawato/gbdraw.git
cd gbdraw
python -m pip install -U pip
python -m pip install -e ".[dev]"
```

Verify the install:

```bash
python -m gbdraw.cli -h
pytest tests/ -v -m "not slow"
```

## Optional: Non-SVG Export Support

SVG export works with the base install. PNG, PDF, EPS, and PS export require CairoSVG:

```bash
python -m pip install -e ".[dev,export]"
```

Depending on your platform, CairoSVG may also require system Cairo/Pango libraries.

## Related Commands

```bash
gbdraw circular --gbk genome.gb -o output -f svg
gbdraw linear --gbk genome1.gb genome2.gb -b comparison.out -o output -f svg
gbdraw gui
```

[< Back to Home](./DOCS.md) | [Go to Quickstart >](./QUICKSTART.md)

[Home](./DOCS.md) | **Installation** | [Quickstart](./QUICKSTART.md) | [Tutorials](./TUTORIALS/TUTORIALS.md) | [Recipes](./RECIPES.md) | [CLI Reference](./CLI_Reference.md) | [Gallery](./GALLERY.md) | [FAQ](./FAQ.md) | [About](./ABOUT.md)
