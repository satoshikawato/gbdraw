[Documentation home](./DOCS.md) | [Tutorials](./TUTORIALS/README.md) | [Technical documentation](./REFERENCE/README.md) | [FAQ](./FAQ.md) | [Gallery](./GALLERY.md) | **Installation** | [About](./ABOUT.md)

[< Back to Home](./DOCS.md) | [Go to Quickstart >](./QUICKSTART.md)

# Installation

Choose an installation route according to the version and interface you need:

| Method | Best for | Notes |
| --- | --- | --- |
| Hosted web app | Making a diagram without a local installation | Runs at [gbdraw.app](https://gbdraw.app/) in your browser. |
| Bioconda | Routine command-line use and reproducible environments | Recommended for most users. |
| PyPI (after publication) | Installing into an existing Python environment | Future release route: `python -m pip install gbdraw`. |
| Source install | Developing or testing the current checkout | Uses `pip install -e ".[dev]"`. |

## 1. Hosted web app

To create a figure without installing `gbdraw`, open:

[https://gbdraw.app/](https://gbdraw.app/)

The hosted app is served by Cloudflare Workers with static assets at the `gbdraw.app`
custom domain. GitHub Pages is not used for deployment.

Pyodide and the main browser-side assets used by the hosted app are vendored and
self-hosted from the repository, so the web UI does not need to fetch those runtime
dependencies from third-party CDNs.

Uploaded genomic data is processed locally in the browser. The hosted deployment uses
Google Analytics 4 for aggregate page-usage metrics; gbdraw does not send uploaded
genome files or generated diagrams to Google Analytics.

The same UI can also be launched locally after installation with:

```bash
gbdraw gui
```

Local `gbdraw gui` analysis runs on your machine. Its packaged Web assets and
browser wheel have no hosted Google Analytics injection. Prepared local installs
include the browser runtime assets and GUI palette data needed for offline
analysis; obtaining the package and its dependencies initially requires network
access or a local package source. The interactive Gallery examples are hosted
separately at [gbdraw.app/gallery](https://gbdraw.app/gallery/) and are not bundled
into local installs.

## 2. Bioconda installation

Bioconda is the recommended local installation path for routine command-line use.
The command installs the version available on that channel; it does not select
the unreleased 0.14.0 checkout.

```bash
mamba create -n gbdraw -c conda-forge -c bioconda gbdraw
conda activate gbdraw
gbdraw -h
```

Optional: launch the local GUI:

```bash
gbdraw gui
```

## 3. PyPI installation

**Not yet published:** the package version remains `0.14.0b0`, and 0.14.0 has
not been published to PyPI. The Trusted Publishing workflow is prepared;
publisher setup and the release transaction must complete before this route is
available. To test the current implementation now, use a source install below.

After publication, install the released package in an activated environment
using a supported Python version (3.10, 3.11, or 3.12):

```bash
python -m pip install gbdraw
gbdraw -h
```

Use an isolated virtual environment rather than modifying the system Python installation.
An unpinned command selects the version available on PyPI, not an unpublished
candidate. For the release's migration notes, see
[0.14.0 release notes (unreleased)](./RELEASE_NOTES_0.14.0.md).

## 4. Source installation for development

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

## Optional: non-SVG export support

SVG export works with the base install. PNG, PDF, EPS, and PS export require
CairoSVG. After PyPI publication:

```bash
python -m pip install "gbdraw[export]"
```

For an editable source checkout, use `python -m pip install -e ".[dev,export]"` instead.

Depending on your platform, CairoSVG may also require system Cairo/Pango
libraries. On Ubuntu, the tested system packages are `libcairo2-dev` and
`libpango1.0-dev`.

## Supported and verified environments

The supported Python versions are 3.10, 3.11, and 3.12. The 0.14.0 development
package has passed isolated wheel/sdist installation, CLI, Python API, session
replay, and export checks on Linux across those versions. The installed local
GUI has also been checked with Chromium. These checks do not establish
Windows/macOS installation, Edge, or later-Python validation.

## Related commands

```bash
gbdraw circular --gbk genome.gb -o output -f svg
gbdraw linear --gbk genome1.gb genome2.gb -b comparison.out -o output -f svg
gbdraw gui
```

[< Back to Home](./DOCS.md) | [Go to Quickstart >](./QUICKSTART.md)

[Documentation home](./DOCS.md) | [Tutorials](./TUTORIALS/README.md) | [Technical documentation](./REFERENCE/README.md) | [FAQ](./FAQ.md) | [Gallery](./GALLERY.md) | **Installation** | [About](./ABOUT.md)
