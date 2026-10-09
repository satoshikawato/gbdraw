[Tutorials](./TUTORIALS/README.md) | [Technical documentation](./REFERENCE/README.md) | [FAQ](./FAQ.md) | [Gallery](./GALLERY.md)

# gbdraw documentation

gbdraw creates publication-quality Circular and Linear genome diagrams in a hosted web app, a local browser app, the command line, or Python.

## Choose a route

| I want to... | Go to | What you will find |
| --- | --- | --- |
| Make my first genome diagram | [Draw a labeled circular map](./TUTORIALS/first-circular-genome-diagram.md) or [a labeled linear map](./TUTORIALS/first-linear-genome-diagram.md) | One figure, built step by step in the web app, on the command line, or in Python |
| Build a complete figure for a specific task | [Tutorials](./TUTORIALS/README.md) | Ten complete figures, from genome comparisons to quantitative tracks and saved sessions |
| Look up exactly what a control, option, or API does | [Technical documentation](./REFERENCE/README.md) | Controls, options, schemas, APIs, compatibility rules, and output formats |
| Choose an approach or fix a problem | [FAQ](./FAQ.md) | Short answers about layouts, interfaces, comparison methods, privacy, publication, and troubleshooting |
| See what gbdraw figures can look like | [Gallery](./GALLERY.md) | Finished static and interactive figures, with links to the Tutorial or settings behind each one |

## Tutorials

Each Tutorial starts from named inputs and ends with a figure you can check.
One page covers the web app, the command line, and Python. Begin with a first
[Circular](./TUTORIALS/first-circular-genome-diagram.md) or
[Linear](./TUTORIALS/first-linear-genome-diagram.md) figure, then choose
another project from [all Tutorials](./TUTORIALS/README.md).

## Technical documentation

Use [Technical documentation](./REFERENCE/README.md) for exact controls, CLI
options, schemas, APIs, compatibility rules, formats, SVG hooks, and source
provenance. Common entry points include:

- [Web app](./REFERENCE/web-app.md)
- [Command line](./REFERENCE/command-line.md)
- [Generated CLI option inventory](./CLI_Reference.md)
- [Python API](./REFERENCE/python-api.md)
- [Typed requests](./REFERENCE/typed-requests.md)
- [Input formats and TSV schemas](./REFERENCE/input-formats-and-tsv-schemas.md)
- [Comparison programs, thresholds, and results](./REFERENCE/comparison-programs-thresholds-and-results.md)
- [Palettes, feature rules, labels, shapes, and tracks](./REFERENCE/palettes-feature-rules-labels-shapes-and-tracks.md)
- [Session and request compatibility](./REFERENCE/session-and-request-compatibility.md)
- [Output formats and export](./REFERENCE/output-formats-and-export.md)
- [Recipes](./RECIPES.md)
- [SVG semantic hooks](./SVG_SEMANTIC_HOOKS.md)

For Linear Similarity Group alignment, see the [Web controls](./REFERENCE/web-app.md#similarity-group-alignment-in-linear-view),
[strict CLI behavior](./REFERENCE/command-line.md#strict-similarity-group-alignment),
[typed Python example](./REFERENCE/python-api.md#typed-linear-similarity-group-alignment),
and [Session compatibility](./REFERENCE/session-and-request-compatibility.md#similarity-alignment-request-ownership).

## FAQ

The [FAQ](./FAQ.md) answers layout, interface, comparison method, privacy,
publication, and troubleshooting questions. Each answer links to the full
technical documentation when needed.

## Gallery

Use the [Gallery](./GALLERY.md) to compare finished figures. Its entries link
back to reproducible Tutorials or the relevant technical documentation.

## Supporting documents

- [Installation](./INSTALL.md)
- [Get the tutorial inputs](./GETTING_TUTORIAL_DATA.md)
- [Palette Explorer](https://gbdraw.app/gallery/palettes/)
- [About and citation](./ABOUT.md)
- [0.14.0 release notes](./RELEASE_NOTES_0.14.0.md)
- [0.14.0b0 beta history](./RELEASE_NOTES_0.14.0b0.md)

## Entry points

- Hosted web app: [gbdraw.app](https://gbdraw.app/)
- Local web app: `gbdraw gui`
- Command line: `gbdraw circular` and `gbdraw linear`
- Python: package-root drawing API or typed `gbdraw.api` requests

[Tutorials](./TUTORIALS/README.md) | [Technical documentation](./REFERENCE/README.md) | [FAQ](./FAQ.md) | [Gallery](./GALLERY.md)
