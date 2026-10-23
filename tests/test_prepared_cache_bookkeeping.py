"""The Worker's prepared-input cache keeps its bookkeeping without walking entries.

OV-1583: a transaction does not read the entries it does not use to report
the retained bytes or to find entries of a removed resource.
OV-1588: the retained bytes count what an entry holds (a comparison frame,
the nucleotide edge text beside it), not the size of the file it came from.
OV-1581: a cached record names its source where the Worker keeps that
resource, not in the workspace of the Generate that parsed it.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from gbdraw.api.prepared import PreparedBiologicalInputCache, PreparedResourceIdentity
from gbdraw.exceptions import ValidationError
from gbdraw.comparisons.losat_jobs import losat_source_ids
from gbdraw.session_request_codec import encode_canonical_request
from tests.test_web_helper_prepared_lookup import (
    _encoded_resources,
    _linear_comparison_request,
    _Stage,
)


class _CountingKey:
    """A cache key that counts how often the cache hashes it."""

    hashes = 0

    def __init__(self, index: int) -> None:
        self.index = index

    def __hash__(self) -> int:
        type(self).hashes += 1
        return hash(("counting-key", self.index))

    def __eq__(self, other: object) -> bool:
        return isinstance(other, _CountingKey) and other.index == self.index


def _resource(tmp_path: Path) -> tuple[Path, PreparedResourceIdentity]:
    path = tmp_path / "source.bin"
    path.write_bytes(b"x" * 64)
    return path, PreparedResourceIdentity("resource-1", "render-resource-1", 64)


def _all_entries_bytes(cache: PreparedBiologicalInputCache) -> int:
    return sum(entry.retained_bytes for entries in cache._layers.values() for entry in entries.values())


def _assert_bookkeeping_matches_a_recount(cache: PreparedBiologicalInputCache) -> None:
    """The running total equals a full recount; the held set covers every entry."""

    assert cache._retained_total == _all_entries_bytes(cache)
    named = {identity for entries in cache._layers.values() for entry in entries.values()
             for identity in entry.resource_identities}
    assert named <= cache._held_identities


class _Failure(Exception):
    pass


def test_a_write_that_fails_after_an_eviction_keeps_the_running_total(tmp_path: Path) -> None:
    """A failed write transaction does not commit, so its evictions keep the total without a recount."""

    first, second = tmp_path / "first.bin", tmp_path / "second.bin"
    first.write_bytes(b"x" * 64)
    second.write_bytes(b"y" * 32)
    kept = PreparedResourceIdentity("resource-1", "render-resource-1", 64)
    replaced = PreparedResourceIdentity("resource-2", "render-resource-2", 32)
    cache = PreparedBiologicalInputCache()
    with cache.transaction(resource_paths={first: kept, second: replaced}, diagnostics=None):
        cache.get_or_build("decoded", ("kept",), frozenset({kept}), object)
        cache.get_or_build("decoded", ("replaced",), frozenset({replaced}), object)
        cache.get_or_build("decoded", ("both",), frozenset({kept, replaced}), object)
    _assert_bookkeeping_matches_a_recount(cache)

    # The second resource gets a new token; the transaction then fails.
    renewed = PreparedResourceIdentity("resource-2", "render-resource-3", 32)
    diagnostics: dict[str, dict[str, int]] = {"metrics": {}}
    with pytest.raises(_Failure):
        with cache.transaction(resource_paths={first: kept, second: renewed}, diagnostics=diagnostics):
            raise _Failure()
    # The entries that name the old token were evicted when the transaction began.
    assert list(cache._layers["decoded"]) == [("kept",)]
    assert diagnostics["metrics"]["preparedInputCacheEvictionCount"] == 2
    _assert_bookkeeping_matches_a_recount(cache)
    assert diagnostics["metrics"]["preparedInputCacheRetainedBytes"] == _all_entries_bytes(cache) > 0


def test_a_lookup_that_evicts_a_mutated_entry_keeps_the_running_total(tmp_path: Path) -> None:
    """A read-only helper lookup never commits; its owner eviction keeps the total."""

    path, identity = _resource(tmp_path)
    cache = PreparedBiologicalInputCache()

    def records() -> tuple[SeqRecord, ...]:
        return (SeqRecord(Seq("ACGT" * 10), id="record", features=[
            SeqFeature(FeatureLocation(0, 9), type="CDS", qualifiers={"product": ["p"]}),
        ]),)

    with cache.transaction(resource_paths={path: identity}, diagnostics=None):
        cached = cache.get_or_build("parsed", ("mutated",), frozenset({identity}), records)
        cache.get_or_build("decoded", ("kept",), frozenset({identity}), object)
    _assert_bookkeeping_matches_a_recount(cache)

    cached[0].features.append(SeqFeature(FeatureLocation(0, 1), type="misc_feature"))
    diagnostics: dict[str, dict[str, int]] = {"metrics": {}}
    with pytest.raises(ValidationError, match="ownership was mutated"):
        with cache.transaction(resource_paths={path: identity}, diagnostics=diagnostics, read_only=True):
            cache.get_or_build("parsed", ("mutated",), frozenset({identity}), records)
    assert not cache._layers["parsed"]
    assert diagnostics["metrics"]["preparedInputCacheEvictionCount"] == 1
    _assert_bookkeeping_matches_a_recount(cache)
    assert diagnostics["metrics"]["preparedInputCacheRetainedBytes"] == _all_entries_bytes(cache) > 0


def test_a_transaction_does_not_walk_the_entries_it_does_not_use(tmp_path: Path) -> None:
    path, identity = _resource(tmp_path)
    cache = PreparedBiologicalInputCache()
    with cache.transaction(resource_paths={path: identity}, diagnostics=None):
        for index in range(50):
            cache.get_or_build("parsed", _CountingKey(index), frozenset({identity}), object)
    assert len(cache._layers["parsed"]) == 50

    _CountingKey.hashes = 0
    diagnostics: dict[str, object] = {"metrics": {}}
    with cache.transaction(resource_paths={path: identity}, diagnostics=diagnostics, read_only=True):
        pass
    assert _CountingKey.hashes == 0
    with cache.transaction(resource_paths={path: identity}, diagnostics=diagnostics):
        pass
    assert _CountingKey.hashes == 0
    # The write transaction used none of them, so its commit released them.
    assert not cache._layers["parsed"]
    assert cache._retained_total == 0


def test_retained_bytes_count_the_comparison_frame_and_the_edge_text(tmp_path: Path) -> None:
    path, identity = _resource(tmp_path)
    frame = pd.DataFrame({"query": ["q"] * 100, "subject": ["s"] * 100, "bitscore": [1.0] * 100})
    text = "q\ts\t1.0\n" * 100
    cache = PreparedBiologicalInputCache()
    diagnostics: dict[str, dict[str, int]] = {"metrics": {}}
    with cache.transaction(resource_paths={path: identity}, diagnostics=diagnostics):
        cache.get_or_build("decoded", ("frame",), frozenset({identity}), lambda: frame)
        cache.get_or_build("decoded", ("edge",), frozenset({identity}), lambda: (text, frame))
    held = diagnostics["metrics"]["preparedInputCacheRetainedBytes"]
    shallow = int(frame.memory_usage(index=True, deep=False).sum())
    # Each object cell holds a string beside its pointer.
    assert cache._layers["decoded"][("frame",)].retained_bytes > shallow
    assert cache._layers["decoded"][("edge",)].retained_bytes >= (
        cache._layers["decoded"][("frame",)].retained_bytes + len(text)
    )
    assert held == _all_entries_bytes(cache) == cache._retained_total


def test_the_running_total_follows_evictions_and_commits(tmp_path: Path) -> None:
    encoded = encode_canonical_request(_linear_comparison_request(tmp_path))
    stage = _Stage(tmp_path / "worker", _encoded_resources(encoded))
    cache = PreparedBiologicalInputCache()
    first = stage.render(cache, encoded.payload)
    assert first["metrics"]["preparedInputCacheRetainedBytes"] == _all_entries_bytes(cache) > 0
    # A new token for one table evicts the entries that name the old one.
    table_id = next(item.resource_id for item in encoded.resources if item.source_path and item.source_path.name == "pair.tsv")
    stage.identities[table_id] = {**stage.identities[table_id], "cacheToken": f"{99:064x}"}
    second = stage.render(cache, encoded.payload)
    assert second["metrics"]["preparedInputCacheEvictionCount"] > 0
    assert second["metrics"]["preparedInputCacheRetainedBytes"] == _all_entries_bytes(cache) == cache._retained_total


def test_a_cached_record_names_the_resource_where_the_worker_keeps_it(tmp_path: Path) -> None:
    encoded = encode_canonical_request(_linear_comparison_request(tmp_path))
    stage = _Stage(tmp_path / "worker", _encoded_resources(encoded))
    cache = PreparedBiologicalInputCache()
    stage.render(cache, encoded.payload)
    second = stage.render(cache, encoded.payload)
    assert second["metrics"]["resolvedRecordCacheHitCount"] == 1
    # Both workspaces are gone; the Worker's resource files stay while their tokens do.
    kept = {str(path.resolve()) for path in stage.cache_paths.values()}
    (collection,) = (entry.value for entry in cache._layers["resolved"].values())
    parsed = [record for entry in cache._layers["parsed"].values() for record in entry.value.records]
    named = {path for item in collection.provenance for path in item.source_paths}
    named |= {path for paths in losat_source_ids(collection.records) for path in paths}
    named |= {record.annotations["gbdraw_source_file"] for record in (*collection.records, *parsed)}
    assert named and named <= kept
    assert all(Path(path).is_file() for path in named)
