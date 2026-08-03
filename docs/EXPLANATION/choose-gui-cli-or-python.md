[Documentation home](../DOCS.md) | [Tutorials](../TUTORIALS/README.md) | [How-to guides](../HOW_TO/README.md) | [Reference](../REFERENCE/README.md)

# Choose the GUI, CLI, or Python interface

All three interfaces use the same drawing engine. Choose the surface that fits how the figure will be created, reviewed, and reproduced.

| Interface | Good fit | Reproducible handoff |
|---|---|---|
| Hosted web app | Interactive exploration without installing Python | Download the input-backed session and final output. |
| Local `gbdraw gui` | The browser workflow with local package assets and no hosted-site dependency | Save a session; keep the package version with it. |
| CLI | Repeatable commands, batch work, and scripted exports | Keep the command, inputs, version, and generated files. |
| Simple Python API | Pipelines that already use Biopython records or need output bytes | Keep the script and pinned environment. |
| Typed request API | Applications that need explicit validation, resource planning, and session codecs | Persist the typed request or session with its schema version. |

The hosted app runs genome parsing, search, rendering, and export in the browser after its static assets load. Uploaded data is not sent to a gbdraw application server. The local GUI serves the same application from the installed package. See [Browser privacy, offline execution, and performance](browser-privacy-offline-execution-and-performance.md) for the boundary.

The CLI exposes built-in protein search through LOSATP or BLASTP. Browser LOSATN and TLOSATX are browser workflows; the CLI consumes their frozen tabular results as prepared comparison input. The [comparison-method guide](choose-a-genome-comparison-method.md) lists these surface differences.

For a first result, use one of the surface-specific Tutorials rather than translating a command by hand: [GUI](../TUTORIALS/GUI/first-circular-genome-diagram.md), [CLI](../TUTORIALS/CLI/first-circular-genome-diagram.md), or [Python](../TUTORIALS/PYTHON/first-genome-diagram.md).
