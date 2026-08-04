[Home](../../DOCS.md) | [Tutorials](../README.md) | [Typed requests](../../REFERENCE/typed-requests.md)

# Create an interactive figure and reproduce it from Python

## Choose how to build this figure

| GUI | CLI | Python API |
| --- | --- | --- |
| [Web-app workflow](../GUI/create-and-resume-an-interactive-figure.md) | [CLI workflow](../CLI/create-and-resume-an-interactive-figure.md) | **This page** |

Build the human mitochondrial figure as a typed request, export its static and
Interactive SVG forms, save a compressed session, and replay that session with
a new output prefix.

## What you'll need

Place `HmmtDNA.gbk` and `cds_gene_qualifier_priority.tsv` in an empty working
directory. Both files are bundled under `gbdraw/web/tutorial-data/`.

## Step 1: Save and run the program

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

## Step 2: Verify the handoff

```bash
python docs/recipes/run_python_scenarios.py --scenario T-PY-08 --check
```

The runner executes the block in a clean directory and checks the typed return
values, current session and request schemas, embedded resources, interactive
controls and feature metadata, and exact static-SVG reproduction.

![Human mitochondrial map restored from Python](../../images/t-py-08/restored_interactive_figure.svg)

## What you built

The request owns the figure settings, the Interactive SVG owns browser-side
inspection, and the compressed session owns the portable replay state. The
session is materialized only while its embedded resources are in use.

## Next steps

- [Build typed requests and round-trip sessions](../../HOW_TO/PYTHON/build-typed-requests-and-round-trip-sessions.md)
- [Review typed-request fields](../../REFERENCE/typed-requests.md)
- [Review session compatibility](../../REFERENCE/session-and-request-compatibility.md)
