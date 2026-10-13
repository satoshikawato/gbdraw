"""Unit tests of the render-stage mechanism (gbdraw.render_stages) and its cache layer."""

from __future__ import annotations

import dataclasses
from pathlib import Path
from typing import Iterator

import pandas as pd
import pytest

import gbdraw.render_stages as stages_module
from gbdraw.api.config import apply_config_overrides
from gbdraw.api.prepared import (
    PreparedBiologicalInputCache,
    PreparedResourceIdentity,
    register_prepared_resource_value,
)
from gbdraw.config.models.root import GbdrawConfig
from gbdraw.render_stages import render_stage

CALLS: list[str] = []


@pytest.fixture(autouse=True)
def _reset_calls() -> Iterator[None]:
    CALLS.clear()
    yield
    CALLS.clear()


@pytest.fixture
def declared(monkeypatch: pytest.MonkeyPatch):
    """Set the declared sections of named stages for one test."""

    def apply(**sections: frozenset[str]) -> None:
        table = dict(stages_module.RENDER_STAGE_SETTINGS)
        table.update(sections)
        monkeypatch.setattr(stages_module, "RENDER_STAGE_SETTINGS", table)

    return apply


def _config(**overrides: object) -> GbdrawConfig:
    return apply_config_overrides(None, dict(overrides) or None)


def _transaction(cache: PreparedBiologicalInputCache, metrics: dict[str, object] | None = None):
    return cache.transaction(resource_paths={}, diagnostics=metrics)


def _metrics(diagnostics: dict[str, object]) -> dict[str, int]:
    metrics = diagnostics["metrics"]
    assert isinstance(metrics, dict)
    return metrics


def double(value: int) -> tuple[int, ...]:
    CALLS.append("double")
    return (value, value)


def total(values: tuple[int, ...]) -> int:
    CALLS.append("total")
    return sum(values)


def count(values: tuple[object, ...]) -> int:
    CALLS.append("count")
    return len(values)


def outer_total(value: int) -> int:
    CALLS.append("outer")
    return render_stage("feature_dom_index", total, render_stage("gc_skew_data", double, value))


def font_size(config: GbdrawConfig) -> float:
    CALLS.append("font_size")
    return config.labels.font_size.long


def font_size_and_gc_window(config: GbdrawConfig) -> tuple[float, int]:
    CALLS.append("font_size_and_gc_window")
    return config.labels.font_size.long, config.objects.sliding_window.default[0]


def outer_font_size(config: GbdrawConfig) -> float:
    CALLS.append("outer_font_size")
    return render_stage("linear_labels", font_size, config)


@dataclasses.dataclass(frozen=True)
class Profile:
    config: GbdrawConfig
    scale: float
    derived: float = dataclasses.field(init=False, default=0.0)

    def __post_init__(self) -> None:
        CALLS.append("profile_post_init")
        object.__setattr__(self, "derived", self.scale * 2)


def profile_font_size(profile: Profile) -> float:
    CALLS.append("profile_font_size")
    return profile.config.labels.font_size.long * profile.derived


def table_rows(table: pd.DataFrame) -> int:
    CALLS.append("table_rows")
    return len(table)


def test_without_a_transaction_every_call_runs() -> None:
    assert render_stage("gc_skew_data", double, 2) == (2, 2)
    assert render_stage("gc_skew_data", double, 2) == (2, 2)
    assert CALLS == ["double", "double"]


def test_a_later_transaction_reuses_a_stage_with_the_same_arguments() -> None:
    cache = PreparedBiologicalInputCache()
    cold: dict[str, object] = {}
    with _transaction(cache, cold):
        assert render_stage("gc_skew_data", double, 3) == (3, 3)
    warm: dict[str, object] = {}
    with _transaction(cache, warm):
        assert render_stage("gc_skew_data", double, 3) == (3, 3)
    changed: dict[str, object] = {}
    with _transaction(cache, changed):
        assert render_stage("gc_skew_data", double, 4) == (4, 4)
    assert CALLS == ["double", "double"]
    assert _metrics(cold)["renderStageBuild:gc_skew_data"] == 1
    assert _metrics(cold)["renderStageHitCount"] == 0
    assert _metrics(warm)["renderStageHit:gc_skew_data"] == 1
    assert _metrics(warm)["renderStageBuildCount"] == 0
    assert _metrics(changed)["renderStageBuildCount"] == 1


def test_an_int_and_a_float_of_the_same_value_have_different_keys() -> None:
    cache = PreparedBiologicalInputCache()
    with _transaction(cache):
        render_stage("gc_skew_data", double, 1)
    with _transaction(cache):
        assert render_stage("gc_skew_data", double, 1.0) == (1.0, 1.0)
    assert CALLS == ["double", "double"]


def test_an_entry_the_transaction_did_not_use_is_dropped_at_commit() -> None:
    cache = PreparedBiologicalInputCache()
    with _transaction(cache):
        render_stage("gc_skew_data", double, 1)
    with _transaction(cache):
        render_stage("gc_skew_data", double, 2)
    with _transaction(cache):
        render_stage("gc_skew_data", double, 1)
    assert CALLS == ["double", "double", "double"]


def test_a_failed_transaction_publishes_nothing() -> None:
    cache = PreparedBiologicalInputCache()
    with pytest.raises(RuntimeError):
        with _transaction(cache):
            render_stage("gc_skew_data", double, 5)
            raise RuntimeError("render failed")
    with _transaction(cache):
        render_stage("gc_skew_data", double, 5)
    assert CALLS == ["double", "double"]


def test_a_stage_result_is_keyed_by_its_stage_and_an_outer_hit_keeps_inner_stages() -> None:
    cache = PreparedBiologicalInputCache()
    with _transaction(cache):
        assert render_stage("legend_table", outer_total, 2) == 4
    assert CALLS == ["outer", "double", "total"]
    CALLS.clear()
    warm: dict[str, object] = {}
    with _transaction(cache, warm):
        assert render_stage("legend_table", outer_total, 2) == 4
    assert CALLS == []
    # The outer hit marked both inner stages used, so they survive the commit.
    inner: dict[str, object] = {}
    with _transaction(cache, inner):
        pair = render_stage("gc_skew_data", double, 2)
        assert render_stage("feature_dom_index", total, pair) == 4
    assert CALLS == []
    assert _metrics(inner)["renderStageHitCount"] == 2


def test_a_result_passed_on_is_keyed_by_identity_not_by_its_value() -> None:
    cache = PreparedBiologicalInputCache()
    with _transaction(cache):
        pair = render_stage("gc_skew_data", double, 2)
        render_stage("feature_dom_index", total, pair)
    CALLS.clear()
    with _transaction(cache):
        # An equal tuple built outside a stage has a value key, not the stage key.
        render_stage("feature_dom_index", total, (2, 2))
    assert CALLS == ["total"]


def test_an_unkeyable_argument_makes_the_stage_and_its_consumers_uncacheable() -> None:
    cache = PreparedBiologicalInputCache()
    marker = object()
    for _ in range(2):
        metrics: dict[str, object] = {}
        with _transaction(cache, metrics):
            pair = render_stage("gc_skew_data", double, marker)  # type: ignore[arg-type]
            render_stage("feature_dom_index", count, pair)
        assert _metrics(metrics)["renderStageUncacheable:gc_skew_data"] == 1
        assert _metrics(metrics)["renderStageUncacheable:feature_dom_index"] == 1
        assert _metrics(metrics)["renderStageUncacheableCount"] == 2


class Doubler:
    def double(self, value: int) -> tuple[int, ...]:
        return (value, value)


def test_only_module_level_functions_and_classes_are_cached() -> None:
    import functools

    cache = PreparedBiologicalInputCache()

    def closure(value: int) -> tuple[int, ...]:
        return (value, value)

    metrics: dict[str, object] = {}
    with _transaction(cache, metrics):
        render_stage("gc_skew_data", closure, 1)
        render_stage("gc_skew_data", (lambda value: (value,)), 1)
        render_stage("gc_skew_data", Doubler().double, 1)
        render_stage("gc_skew_data", functools.partial(double), 1)
        render_stage("gc_skew_data", double, 1)
    assert _metrics(metrics)["renderStageUncacheableCount"] == 4
    assert _metrics(metrics)["renderStageBuild:gc_skew_data"] == 1


def test_a_declared_section_is_in_the_key_and_others_are_not(declared) -> None:
    declared(linear_labels=frozenset({"labels.font_size"}))
    cache = PreparedBiologicalInputCache()
    base = _config()
    with _transaction(cache):
        assert render_stage("linear_labels", font_size, base) == base.labels.font_size.long
    other_section = _config(**{"objects.sliding_window.default": [77, 7]})
    warm: dict[str, object] = {}
    with _transaction(cache, warm):
        render_stage("linear_labels", font_size, other_section)
    assert _metrics(warm)["renderStageHit:linear_labels"] == 1
    bigger = _config(**{"labels.font_size.long": base.labels.font_size.long + 3})
    with _transaction(cache):
        assert render_stage("linear_labels", font_size, bigger) == bigger.labels.font_size.long
    assert CALLS == ["font_size", "font_size"]


def test_an_undeclared_read_is_returned_counted_and_never_stored(declared) -> None:
    declared(linear_labels=frozenset({"labels.font_size"}))
    cache = PreparedBiologicalInputCache()
    config = _config()
    for _ in range(2):
        metrics: dict[str, object] = {}
        with _transaction(cache, metrics):
            value = render_stage("linear_labels", font_size_and_gc_window, config)
            render_stage("feature_dom_index", count, value)
        assert value == (config.labels.font_size.long, config.objects.sliding_window.default[0])
        assert _metrics(metrics)["renderStageUndeclaredRead:linear_labels"] == 1
        assert _metrics(metrics)["renderStageUndeclaredReadCount"] == 1
        # The unstored result has no identity, so its consumer is not cached either.
        assert _metrics(metrics)["renderStageUncacheable:feature_dom_index"] == 1
    assert CALLS.count("font_size_and_gc_window") == 2


def test_an_inner_stage_that_declares_more_marks_the_outer_stage(declared) -> None:
    declared(
        legend_table=frozenset({"objects.sliding_window"}),
        linear_labels=frozenset({"labels.font_size"}),
    )
    cache = PreparedBiologicalInputCache()
    metrics: dict[str, object] = {}
    with _transaction(cache, metrics):
        render_stage("legend_table", outer_font_size, _config())
    assert _metrics(metrics)["renderStageUndeclaredRead:legend_table"] == 1
    # The inner stage declared what it read and was stored.
    assert "renderStageUndeclaredRead:linear_labels" not in _metrics(metrics)
    assert _metrics(metrics)["renderStageUndeclaredReadCount"] == 1


def test_an_undeclared_read_in_an_inner_stage_marks_every_enclosing_stage(declared) -> None:
    declared(
        legend_table=frozenset({"labels.font_size"}),
        linear_labels=frozenset({"labels.font_size"}),
    )
    cache = PreparedBiologicalInputCache()

    metrics: dict[str, object] = {}
    with _transaction(cache, metrics):
        render_stage("legend_table", outer_font_size, _config())
    assert _metrics(metrics)["renderStageUndeclaredReadCount"] == 0
    assert _metrics(metrics)["renderStageBuildCount"] == 2


def outer_window_through_inner(config: GbdrawConfig) -> tuple[float, int]:
    return render_stage("linear_labels", font_size_and_gc_window, config)


def test_an_inner_read_undeclared_by_both_marks_both(declared) -> None:
    declared(
        legend_table=frozenset({"labels.font_size"}),
        linear_labels=frozenset({"labels.font_size"}),
    )
    cache = PreparedBiologicalInputCache()
    metrics: dict[str, object] = {}
    with _transaction(cache, metrics):
        render_stage("legend_table", outer_window_through_inner, _config())
    assert _metrics(metrics)["renderStageUndeclaredRead:legend_table"] == 1
    assert _metrics(metrics)["renderStageUndeclaredRead:linear_labels"] == 1


def test_a_config_in_a_frozen_dataclass_gets_a_view_without_rerunning_post_init(declared) -> None:
    declared(linear_labels=frozenset({"labels.font_size"}))
    cache = PreparedBiologicalInputCache()
    profile = Profile(_config(), 1.5)
    CALLS.clear()
    metrics: dict[str, object] = {}
    with _transaction(cache, metrics):
        result = render_stage("linear_labels", profile_font_size, profile)
    assert result == profile.config.labels.font_size.long * 3.0
    assert CALLS == ["profile_font_size"]
    assert _metrics(metrics)["renderStageUndeclaredReadCount"] == 0
    with _transaction(cache):
        render_stage("linear_labels", profile_font_size, Profile(profile.config, 1.5))
    assert CALLS == ["profile_font_size", "profile_post_init"]


def test_a_config_inside_a_tuple_is_not_keyable(declared) -> None:
    declared(linear_labels=frozenset({"labels.font_size"}))
    cache = PreparedBiologicalInputCache()
    metrics: dict[str, object] = {}
    with _transaction(cache, metrics):
        render_stage("linear_labels", count, (_config(),))
    assert _metrics(metrics)["renderStageUncacheable:linear_labels"] == 1


def test_a_table_is_keyed_by_its_resource_identity(tmp_path: Path) -> None:
    path = tmp_path / "colors.tsv"
    path.write_text("x\n", encoding="utf-8")
    identity = PreparedResourceIdentity("colors", "token-1", 2)
    cache = PreparedBiologicalInputCache()
    for _ in range(2):
        metrics: dict[str, object] = {}
        with cache.transaction(resource_paths={path: identity}, diagnostics=metrics):
            table = register_prepared_resource_value(pd.DataFrame({"a": [1, 2]}), path)
            assert render_stage("color_rules", table_rows, table) == 2
            render_stage("color_rules", table_rows, pd.DataFrame({"a": [1]}))
        assert _metrics(metrics)["renderStageUncacheable:color_rules"] == 1
    assert CALLS == ["table_rows", "table_rows", "table_rows"]
    changed: dict[str, object] = {}
    other = PreparedResourceIdentity("colors", "token-2", 2)
    with cache.transaction(resource_paths={path: other}, diagnostics=changed):
        table = register_prepared_resource_value(pd.DataFrame({"a": [1, 2]}), path)
        render_stage("color_rules", table_rows, table)
    assert _metrics(changed)["renderStageBuild:color_rules"] == 1


def test_unknown_stage_names_are_rejected() -> None:
    with pytest.raises(ValueError):
        render_stage("not_a_stage", double, 1)  # type: ignore[arg-type]


def test_label_filter_maps_leave_the_config_mapping_unchanged() -> None:
    from gbdraw.labels.filtering import DERIVED_LABEL_FILTERING_KEYS, label_filter_maps

    config = _config()
    raw = config.labels.filtering.as_dict()
    before = dict(raw)
    maps = label_filter_maps(raw)
    assert raw == before
    assert DERIVED_LABEL_FILTERING_KEYS <= set(maps)
    cache = PreparedBiologicalInputCache()
    for _ in range(2):
        metrics: dict[str, object] = {}
        with _transaction(cache, metrics):
            label_filter_maps(raw)
            label_filter_maps(dict(raw))
        assert _metrics(metrics)["renderStageHit:label_filter_maps"] >= 1
    assert raw == before
