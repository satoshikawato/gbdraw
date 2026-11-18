"""The frozen browser LOSAT evidence of the Lambda and DE3 tutorial data matches its manifest entries and the shipped
serial LOSAT file, without a browser (``tools/qualify_lambda_de3_comparison_fixture.py``)."""

from __future__ import annotations

import copy
import dataclasses
import importlib.util
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
_SPEC = importlib.util.spec_from_file_location(
    "qualify_lambda_de3_comparison_fixture", REPO_ROOT / "tools" / "qualify_lambda_de3_comparison_fixture.py"
)
assert _SPEC is not None and _SPEC.loader is not None
qualifier = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = qualifier
_SPEC.loader.exec_module(qualifier)


def test_the_frozen_evidence_was_qualified_with_the_shipped_serial_losat() -> None:
    # The swap to LOSAT 0.2.0 (#1540) left the evidence qualified with 0.1.0, and no check ran.
    summaries = qualifier._check()
    assert sorted(summaries) == ["lambda-de3-losatn", "lambda-de3-tlosatx"]


def test_writing_the_evidence_writes_every_frozen_value_of_its_entry() -> None:
    # --write records a TSV's summary into its manifest entry; recording the frozen TSV again changes nothing, so the
    # entry holds no value that --write leaves stale.
    manifest = json.loads(qualifier.FIXTURE_MANIFEST_PATH.read_text(encoding="utf-8"))
    build = json.loads(qualifier.BUILD_MANIFEST_PATH.read_text(encoding="utf-8"))
    for item in qualifier._evidence(manifest):
        written = dataclasses.replace(item, entry=copy.deepcopy(item.entry))
        summary = qualifier._summarize(item.output_path.read_bytes(), item)
        qualifier._record(written, summary, build, item.entry["derivation"]["arguments"]["runs"])
        assert written.entry == item.entry
        assert set(qualifier._frozen(item.entry)) <= set(summary)
