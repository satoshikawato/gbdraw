[Home](../../DOCS.md) | [Tutorials](../README.md) | [Get tutorial files](../../GETTING_TUTORIAL_DATA.md) | [Typed requests](../../REFERENCE/typed-requests.md)

# Create an interactive figure and reproduce it from Python

## Choose how to build this figure

| GUI | CLI | Python API |
| --- | --- | --- |
| [Web-app workflow](../GUI/create-and-resume-an-interactive-figure.md) | [CLI workflow](../CLI/create-and-resume-an-interactive-figure.md) | **This page** |

Build the human mitochondrial figure as a typed request, export its static and
Interactive SVG forms, save a compressed session, and replay that session with
a new output prefix.

## What you'll need

| Kind | Filename | What to do |
| --- | --- | --- |
| Download | `HmmtDNA.gbk` ([NCBI accession NC_012920.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_012920.1)) | Download the complete GenBank record and save it as `HmmtDNA.gbk`. |
| Download | [`cds_gene_qualifier_priority.tsv`](../../../gbdraw/web/tutorial-data/shared/cds_gene_qualifier_priority.tsv) | Save the supplied label rule as `cds_gene_qualifier_priority.tsv`. |
| Create | `interactive_handoff.py` | Save the complete Python program from Step 2 with this filename. |
| Generated | `interactive_human_mitochondrion.svg` | The first render writes this static master. |
| Generated | `interactive_human_mitochondrion.interactive.svg` | The first render writes this offline Interactive SVG. |
| Generated | `interactive_handoff.gbdraw-session.json.gz` | The program saves this compressed session and then reloads it. |
| Generated | `restored_interactive_figure.svg` | The replay writes this static SVG. |
| Reference result | [`restored_interactive_figure.svg`](../../images/t-py-08/restored_interactive_figure.svg) | Compare the restored Generated SVG with this versioned result. |

Install gbdraw in the active Python environment before starting.

## Step 1: Prepare the working directory

Create and enter an empty directory:

```bash
mkdir gbdraw-python-interactive-session
cd gbdraw-python-interactive-session
```

Acquire the sequence from its authoritative accession link and download the
supplied label rule. Save both with the exact filenames shown. See [Get the
tutorial files](../../GETTING_TUTORIAL_DATA.md) for browser, `curl`, and
PowerShell instructions for authoritative sequence acquisition.

## Step 2: Save and run the Python program

Save the following complete program as `interactive_handoff.py`:

<!-- executable:T-PY-08:start -->
```python
from datetime import datetime, timezone
from pathlib import Path

from gbdraw.api import (
    CircularDiagramOptions,
    CircularDiagramRequest,
    CircularOutputOptions,
    GenBankInputSource,
    RecordInput,
    RenderOutputRequest,
    RequestRenderResult,
    load_session_document,
    materialize_session,
    render_request,
    save_session_document,
    session_to_request,
    with_request_output,
)


request = CircularDiagramRequest(
    records=(
        RecordInput(
            source=GenBankInputSource("HmmtDNA.gbk"),
            record_key="human-mitochondrion",
        ),
    ),
    options=CircularDiagramOptions(
        qualifier_priority_file="cds_gene_qualifier_priority.tsv",
        species="<i>Homo sapiens</i>",
        output=CircularOutputOptions(legend="right"),
        config_overrides={
            "canvas.strandedness": True,
            "canvas.circular.track_type": "middle",
            "canvas.show_gc": True,
            "canvas.show_skew": True,
            "labels.circular.scope": "outer",
            "labels.circular.placement": "horizontal",
        },
    ),
    output=RenderOutputRequest(
        output_prefix="interactive_human_mitochondrion",
        formats=("svg", "interactive_svg"),
    ),
)

result = render_request(request)
assert isinstance(result, RequestRenderResult)

session_path = Path("interactive_handoff.gbdraw-session.json.gz")
session_document = save_session_document(
    session_path,
    request,
    title="Interactive human mitochondrial handoff",
    created_at=datetime(2026, 8, 4, tzinfo=timezone.utc),
)
loaded_document = load_session_document(session_path)

with materialize_session(loaded_document, output_directory=Path(".")) as materialized:
    replay_request = session_to_request(materialized)
    restored_request = with_request_output(
        replay_request,
        output_prefix="restored_interactive_figure",
        formats=("svg",),
    )
    restored_result = render_request(restored_request)

assert isinstance(restored_result, RequestRenderResult)
print("Exported the interactive figure and session")
print("Restored restored_interactive_figure.svg")
```
<!-- executable:T-PY-08:end -->

Before running it, your working directory should contain:

```text
gbdraw-python-interactive-session/
├── HmmtDNA.gbk
├── cds_gene_qualifier_priority.tsv
└── interactive_handoff.py
```

### Export and replay the session

```bash
python interactive_handoff.py
```

Expected output: the program prints
`Exported the interactive figure and session`, followed by
`Restored restored_interactive_figure.svg`. It writes all four Generated files
listed above.

## Step 3: Inspect the restored figure

Open the Generated `restored_interactive_figure.svg`. It should reproduce the
original static export exactly, using the embedded resources and feature
metadata carried by the saved session.

![Human mitochondrial map restored from Python](../../images/t-py-08/restored_interactive_figure.svg)

The image above is the Reference result. Compare the record identity, labels,
tracks, and legend with your restored SVG. The two static Generated SVGs should
also be byte-identical.

## Next steps

- [Build typed requests and round-trip sessions](../../HOW_TO/PYTHON/build-typed-requests-and-round-trip-sessions.md)
- [Review typed-request fields](../../REFERENCE/typed-requests.md)
- [Review session compatibility](../../REFERENCE/session-and-request-compatibility.md)
