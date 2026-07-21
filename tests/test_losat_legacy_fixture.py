"""Characterize the real schema-v2 LOSATP migration fixture."""

from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path


FIXTURE_DIR = Path(__file__).parent / "fixtures" / "sessions"
SESSION_PATH = FIXTURE_DIR / "BGC0000708-BGC0000713.schema-v2.gbdraw-session.json.gz"
EXPECTED_PATH = FIXTURE_DIR / "BGC0000708-BGC0000713.schema-v2.expected.json"


def test_real_legacy_losat_fixture_matches_its_historical_oracle() -> None:
    expected = json.loads(EXPECTED_PATH.read_text(encoding="utf-8"))
    source_bytes = gzip.decompress(SESSION_PATH.read_bytes())
    session = json.loads(source_bytes)

    assert hashlib.sha256(source_bytes).hexdigest() == expected["sourceSha256"]
    assert expected["sourceCommit"] == "c64ff8c3b42fc315975dc36cf1942f471c93e847"
    assert session["version"] == expected["sessionVersion"] == 33
    assert session["renderRequest"]["schema"] == expected["renderRequestSchema"] == 2

    entries = session["losatCache"]["entries"]
    assert len(entries) == expected["storedRawEntries"] == 34
    assert all(entry["schema"] == 2 for entry in entries)
    assert all(entry["kind"] == "raw-losat" for entry in entries)
    assert all(entry["program"] == "blastp" for entry in entries)

    record_count = len(session["renderRequest"]["records"])
    assert session["config"]["losat"]["blastp"]["mode"] == "orthogroup"
    assert record_count * record_count == expected["totalPairs"] == 25
    assert expected["cacheHits"] == 25
    assert expected["cacheMisses"] == expected["uniqueJobs"] == expected["workerCalls"] == 0
