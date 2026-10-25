"""A Cancel stops a Worker request only at the prepared-cache checkpoints.

The page writes SIGINT into the Worker's interrupt buffer; the Worker's
``WorkerCancelRequest`` records it, and ``gbdraw.api.prepared.cancel_checkpoint``
raises ``KeyboardInterrupt`` at the next point between builds. A render
stopped at any checkpoint leaves a cache whose next render equals a render
from a fresh cache (an interrupted build is never kept).
"""

from __future__ import annotations

import signal
import sys
from pathlib import Path

import pytest

from gbdraw.api.prepared import (
    PreparedBiologicalInputCache,
    PreparedResourceIdentity,
    cancel_checkpoint,
    get_or_build_parsed_source,
)
from gbdraw.web_support.worker_cancel import WorkerCancelRequest
from tests.test_web_render_preparation_cache import (
    _document,
    _metrics,
    _record,
    _render,
    _replace_first_feature,
    _write_genbank,
)


class _StopAt:
    """A cancel check that stops at its ``stop_at``-th call; ``None`` never stops."""

    def __init__(self, stop_at: int | None = None) -> None:
        self.calls = 0
        self.stop_at = stop_at

    def __call__(self) -> None:
        self.calls += 1
        if self.stop_at is not None and self.calls == self.stop_at:
            raise KeyboardInterrupt


def _comparable(result: dict[str, object]) -> dict[str, object]:
    """The reply without the catalog revision, which names one cache's build."""

    return {key: value for key, value in result.items() if key != "featureCatalogRevision"}


def _sweep(tmp_path: Path, documents: list[dict[str, object]], canceled: dict[str, object]) -> int:
    """Stop ``canceled`` at each checkpoint after ``documents`` committed, then render it again.

    Every second render must equal the render of a fresh cache. Returns the
    number of checkpoints swept.
    """

    request_ids = iter(range(1, 100_000))
    expected, _ = _render(tmp_path, PreparedBiologicalInputCache(), canceled, request_id=next(request_ids))
    counting = _StopAt()
    warm = PreparedBiologicalInputCache(cancel_check=counting)
    for document in documents:
        _render(tmp_path, warm, document, request_id=next(request_ids))
    counting.calls = 0
    _render(tmp_path, warm, canceled, request_id=next(request_ids))
    checkpoints = counting.calls
    assert checkpoints > 0
    for stop_at in range(1, checkpoints + 1):
        stop = _StopAt()
        cache = PreparedBiologicalInputCache(cancel_check=stop)
        for document in documents:
            _render(tmp_path, cache, document, request_id=next(request_ids))
        stop.calls = 0
        stop.stop_at = stop_at
        with pytest.raises(KeyboardInterrupt):
            _render(tmp_path, cache, canceled, request_id=next(request_ids))
        stop.stop_at = None
        again, _ = _render(tmp_path, cache, canceled, request_id=next(request_ids))
        assert _comparable(again) == _comparable(expected), f"checkpoint {stop_at} of {checkpoints}"
    return checkpoints


def test_a_render_stopped_at_any_checkpoint_renders_again_like_a_fresh_cache(tmp_path: Path) -> None:
    source = tmp_path / "source.gbk"
    _write_genbank(source, (_record("cancel"),))
    assert _sweep(tmp_path, [], _document(source, title="First")) > 5
    # A re-render after a committed one: stages and inputs are published, and
    # the canceled render reuses some and builds others.
    assert _sweep(tmp_path, [_document(source, title="First")], _document(source, title="Second")) > 0
    assert _sweep(
        tmp_path, [_document(source, title="First")], _document(source, title="First", reverse=True)
    ) > 5


def test_a_checkpoint_outside_a_transaction_does_nothing() -> None:
    # The CLI and the Python API never open a transaction: no check runs.
    cancel_checkpoint()


def test_the_cancel_request_records_sigint_and_raises_only_at_a_check() -> None:
    cancel = WorkerCancelRequest()
    previous = signal.signal(signal.SIGINT, cancel.on_signal)
    try:
        signal.raise_signal(signal.SIGINT)
        with pytest.raises(KeyboardInterrupt):
            cancel.check()
        # It stays requested until the drain clears it.
        with pytest.raises(KeyboardInterrupt):
            cancel.check()
        cancel.clear()
        cancel.check()
    finally:
        signal.signal(signal.SIGINT, previous)


def test_a_stop_silently_evicts_a_cached_input_the_canceled_build_changed(tmp_path: Path) -> None:
    path = tmp_path / "resource.bin"
    path.write_bytes(b"resource")
    identity = PreparedResourceIdentity("record-1", "render-resource-1", 8)
    identities = frozenset({identity})
    stop = _StopAt()
    cache = PreparedBiologicalInputCache(cancel_check=stop)
    cached = (_record("cached"),)
    with cache.transaction(resource_paths={path: identity}, diagnostics=None):
        get_or_build_parsed_source(("source", identity), identities, lambda: cached)

    diagnostics: dict[str, object] = {"metrics": {}}
    with pytest.raises(KeyboardInterrupt):
        with cache.transaction(resource_paths={path: identity}, diagnostics=diagnostics):
            assert get_or_build_parsed_source(
                ("source", identity), identities, lambda: pytest.fail("the committed source is reused")
            ) is cached
            # A builder changed the record it read, then the Cancel stopped it.
            _replace_first_feature(cached)
            stop.stop_at = stop.calls + 1
            get_or_build_parsed_source(("other", identity), identities, lambda: pytest.fail("stopped first"))
    assert ("source", identity) not in cache._layers["parsed"]
    assert _metrics(diagnostics)["preparedInputCacheEvictionCount"] == 1

    # The next request rebuilds it instead of failing on "ownership was mutated".
    stop.stop_at = None
    rebuilt = (_record("rebuilt"),)
    with cache.transaction(resource_paths={path: identity}, diagnostics=None):
        assert get_or_build_parsed_source(("source", identity), identities, lambda: rebuilt) is rebuilt


class _CallSites:
    """A cancel check that records which function called ``cancel_checkpoint``."""

    def __init__(self) -> None:
        self.sites: set[str] = set()

    def __call__(self) -> None:
        # check_cancel <- cancel_checkpoint <- the build that asked.
        checkpoint = sys._getframe(2)
        if checkpoint.f_code is cancel_checkpoint.__code__:
            caller = checkpoint.f_back
            assert caller is not None
            self.sites.add(caller.f_code.co_name)


def test_a_cold_render_checks_inside_the_long_builds(tmp_path: Path) -> None:
    # Tier 2: the builds that took longest between checkpoints check per feature
    # or per phase, so a Cancel does not wait for a whole build.
    source = tmp_path / "source.gbk"
    _write_genbank(source, (_record("cold"),))
    sites = _CallSites()
    _render(tmp_path, PreparedBiologicalInputCache(cancel_check=sites), _document(source), request_id=1)
    assert {
        "_build_feature_layers",
        "extract_features_from_records_payload",
        "build_interactive_svg_context",
        "_normalized_biological_features",
        "_feature_catalog_item",
    } <= sites.sites
