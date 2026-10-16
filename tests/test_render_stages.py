"""Unit tests of the render-stage mechanism (gbdraw.render_stages) and its cache layer."""

from __future__ import annotations

import dataclasses
from pathlib import Path
from typing import Iterator

import pandas as pd
import pytest
from fontTools.ttLib import TTFont

import gbdraw.render_stages as stages_module
from gbdraw.api.config import apply_config_overrides
from gbdraw.api.prepared import (
    PreparedBiologicalInputCache,
    PreparedResourceIdentity,
    register_prepared_resource_value,
)
from gbdraw.config.models.root import GbdrawConfig
from gbdraw.core import fonts
from gbdraw.core.text import measurement_font_file
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
    return render_stage("legend_measurement", font_size, config)


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
    declared(legend_measurement=frozenset({"labels.font_size"}))
    cache = PreparedBiologicalInputCache()
    base = _config()
    with _transaction(cache):
        assert render_stage("legend_measurement", font_size, base) == base.labels.font_size.long
    other_section = _config(**{"objects.sliding_window.default": [77, 7]})
    warm: dict[str, object] = {}
    with _transaction(cache, warm):
        render_stage("legend_measurement", font_size, other_section)
    assert _metrics(warm)["renderStageHit:legend_measurement"] == 1
    bigger = _config(**{"labels.font_size.long": base.labels.font_size.long + 3})
    with _transaction(cache):
        assert render_stage("legend_measurement", font_size, bigger) == bigger.labels.font_size.long
    assert CALLS == ["font_size", "font_size"]


def test_an_undeclared_read_is_returned_counted_and_never_stored(declared) -> None:
    declared(legend_measurement=frozenset({"labels.font_size"}))
    cache = PreparedBiologicalInputCache()
    config = _config()
    for _ in range(2):
        metrics: dict[str, object] = {}
        with _transaction(cache, metrics):
            value = render_stage("legend_measurement", font_size_and_gc_window, config)
            render_stage("feature_dom_index", count, value)
        assert value == (config.labels.font_size.long, config.objects.sliding_window.default[0])
        assert _metrics(metrics)["renderStageUndeclaredRead:legend_measurement"] == 1
        assert _metrics(metrics)["renderStageUndeclaredReadCount"] == 1
        # The unstored result has no identity, so its consumer is not cached either.
        assert _metrics(metrics)["renderStageUncacheable:feature_dom_index"] == 1
    assert CALLS.count("font_size_and_gc_window") == 2


def test_an_inner_stage_that_declares_more_marks_the_outer_stage(declared) -> None:
    declared(
        legend_table=frozenset({"objects.sliding_window"}),
        legend_measurement=frozenset({"labels.font_size"}),
    )
    cache = PreparedBiologicalInputCache()
    metrics: dict[str, object] = {}
    with _transaction(cache, metrics):
        render_stage("legend_table", outer_font_size, _config())
    # The inner stage declared what it read and was stored: only the outer stage is marked.
    assert {
        name: count
        for name, count in _metrics(metrics).items()
        if name.startswith("renderStageUndeclaredRead:")
    } == {"renderStageUndeclaredRead:legend_table": 1}
    assert _metrics(metrics)["renderStageBuild:legend_measurement"] == 1
    assert _metrics(metrics)["renderStageUncacheableCount"] == 0
    assert _metrics(metrics)["renderStageUndeclaredReadCount"] == 1


def test_an_inner_read_declared_by_both_marks_neither(declared) -> None:
    declared(
        legend_table=frozenset({"labels.font_size"}),
        legend_measurement=frozenset({"labels.font_size"}),
    )
    cache = PreparedBiologicalInputCache()

    metrics: dict[str, object] = {}
    with _transaction(cache, metrics):
        render_stage("legend_table", outer_font_size, _config())
    assert _metrics(metrics)["renderStageUndeclaredReadCount"] == 0
    assert _metrics(metrics)["renderStageBuildCount"] == 2


def outer_window_through_inner(config: GbdrawConfig) -> tuple[float, int]:
    return render_stage("legend_measurement", font_size_and_gc_window, config)


def test_an_inner_read_undeclared_by_both_marks_both(declared) -> None:
    declared(
        legend_table=frozenset({"labels.font_size"}),
        legend_measurement=frozenset({"labels.font_size"}),
    )
    cache = PreparedBiologicalInputCache()
    metrics: dict[str, object] = {}
    with _transaction(cache, metrics):
        render_stage("legend_table", outer_window_through_inner, _config())
    assert _metrics(metrics)["renderStageUndeclaredRead:legend_table"] == 1
    assert _metrics(metrics)["renderStageUndeclaredRead:legend_measurement"] == 1


def test_a_config_in_a_frozen_dataclass_gets_a_view_without_rerunning_post_init(declared) -> None:
    declared(legend_measurement=frozenset({"labels.font_size"}))
    cache = PreparedBiologicalInputCache()
    profile = Profile(_config(), 1.5)
    CALLS.clear()
    metrics: dict[str, object] = {}
    with _transaction(cache, metrics):
        result = render_stage("legend_measurement", profile_font_size, profile)
    assert result == profile.config.labels.font_size.long * 3.0
    assert CALLS == ["profile_font_size"]
    assert _metrics(metrics)["renderStageUndeclaredReadCount"] == 0
    with _transaction(cache):
        render_stage("legend_measurement", profile_font_size, Profile(profile.config, 1.5))
    assert CALLS == ["profile_font_size", "profile_post_init"]


def test_a_config_inside_a_tuple_is_not_keyable(declared) -> None:
    declared(legend_measurement=frozenset({"labels.font_size"}))
    cache = PreparedBiologicalInputCache()
    metrics: dict[str, object] = {}
    with _transaction(cache, metrics):
        render_stage("legend_measurement", count, (_config(),))
    assert _metrics(metrics)["renderStageUncacheable:legend_measurement"] == 1


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


def test_a_label_fit_override_keeps_the_label_filtering_objects() -> None:
    # The label fit verification draws with overrides; the label tables in
    # labels.filtering keep their prepared identity there (no deep copy).
    from gbdraw.config.modify import with_config_overrides

    config = _config()
    fitted = with_config_overrides(config, {"labels.spacing.circular": 0.5})
    assert fitted.labels.spacing.circular == 0.5
    assert fitted.labels.filtering is config.labels.filtering
    filtered = with_config_overrides(config, {"labels.filtering.blacklist_keywords": ["x"]})
    assert filtered.labels.filtering is not config.labels.filtering
    assert filtered.labels.filtering.blacklist_keywords == ["x"]


# --- review fixes (REVIEW-PR2-MECHANISM.md) ----------------------------------


def section_long(section: object) -> float:
    CALLS.append("section_long")
    return section.long  # type: ignore[attr-defined]


def outer_section_long(config: GbdrawConfig) -> float:
    return render_stage("legend_measurement", section_long, config.labels.font_size)


def test_a_section_passed_to_an_inner_stage_marks_the_outer_stage() -> None:
    cache = PreparedBiologicalInputCache()
    base = _config()
    for _ in range(2):
        metrics: dict[str, object] = {}
        with _transaction(cache, metrics):
            assert render_stage("legend_table", outer_section_long, base) == base.labels.font_size.long
        # The outer stage did not declare labels.font_size, on a build and on an inner hit.
        assert _metrics(metrics)["renderStageUndeclaredRead:legend_table"] == 1
    bigger = _config(**{"labels.font_size.long": base.labels.font_size.long + 5})
    with _transaction(cache):
        assert render_stage("legend_table", outer_section_long, bigger) == bigger.labels.font_size.long


def test_a_section_passed_by_value_is_keyed_like_the_section(declared) -> None:
    declared(legend_table=frozenset({"labels.font_size"}))
    cache = PreparedBiologicalInputCache()
    config = _config()
    for _ in range(2):
        with _transaction(cache):
            render_stage("legend_table", outer_section_long, config)
    # The inner stage is keyed by the section's value, not by the view's class.
    with _transaction(cache):
        assert render_stage("legend_measurement", section_long, config.labels.font_size) == config.labels.font_size.long
    assert CALLS == ["section_long"]


@dataclasses.dataclass(frozen=True)
class Holder:
    config: GbdrawConfig


def hold(config: GbdrawConfig) -> Holder:
    return Holder(config)


def read_held_show_gc(holder: Holder) -> bool:
    return holder.config.canvas.show_gc


def test_a_result_that_holds_a_config_is_never_stored() -> None:
    cache = PreparedBiologicalInputCache()
    for show_gc in (True, False):
        metrics: dict[str, object] = {}
        with _transaction(cache, metrics):
            holder = render_stage("legend_table", hold, _config(**{"canvas.show_gc": show_gc}))
            assert render_stage("feature_dom_index", read_held_show_gc, holder) is show_gc
        assert _metrics(metrics)["renderStageUndeclaredRead:legend_table"] == 1


def compare_font_size(config: GbdrawConfig, other: object) -> bool:
    return config.labels.font_size == other


def deep_copied_font_size(config: GbdrawConfig) -> float:
    import copy

    return copy.deepcopy(config).labels.font_size.long


def test_a_view_compares_like_its_section_and_the_comparison_is_a_read() -> None:
    cache = PreparedBiologicalInputCache()
    config = _config()
    metrics: dict[str, object] = {}
    with _transaction(cache, metrics):
        assert render_stage("legend_table", compare_font_size, config, config.labels.font_size) is True
    assert _metrics(metrics)["renderStageUndeclaredRead:legend_table"] == 1


def test_a_deep_copy_of_a_view_still_flags_reads() -> None:
    cache = PreparedBiologicalInputCache()
    metrics: dict[str, object] = {}
    with _transaction(cache, metrics):
        render_stage("legend_table", deep_copied_font_size, _config())
    assert _metrics(metrics)["renderStageUndeclaredRead:legend_table"] == 1


def wrap_values(raw: dict[str, object]) -> dict[str, object]:
    CALLS.append("wrap_values")
    return {**raw, "compiled": True}


def test_a_hit_keeps_the_keys_of_values_passed_through(tmp_path: Path) -> None:
    path = tmp_path / "table.tsv"
    path.write_text("x\n", encoding="utf-8")
    identity = PreparedResourceIdentity("table", "token-1", 2)
    cache = PreparedBiologicalInputCache()
    for _ in range(2):
        metrics: dict[str, object] = {}
        with cache.transaction(resource_paths={path: identity}, diagnostics=metrics):
            table = register_prepared_resource_value(pd.DataFrame({"a": [1, 2]}), path)
            wrapped = render_stage("label_filter_maps", wrap_values, {"values": [1, 2], "table": table})
            render_stage("feature_dom_index", count, wrapped["values"])
            render_stage("legend_table", table_rows, wrapped["table"])
    assert CALLS == ["wrap_values", "count", "table_rows"]
    assert _metrics(metrics)["renderStageHitCount"] == 3


def test_a_class_defined_in_a_function_is_not_cached() -> None:
    @dataclasses.dataclass(frozen=True)
    class Local:
        value: int

    cache = PreparedBiologicalInputCache()
    metrics: dict[str, object] = {}
    with _transaction(cache, metrics):
        render_stage("gc_skew_data", Local, 1)
    assert _metrics(metrics)["renderStageUncacheable:gc_skew_data"] == 1


def test_container_subclasses_and_signed_zeros_have_their_own_keys() -> None:
    import collections

    Pair = collections.namedtuple("Pair", "a b")
    cache = PreparedBiologicalInputCache()
    with _transaction(cache):
        render_stage("feature_dom_index", count, (1, 2))
        render_stage("feature_dom_index", count, Pair(1, 2))
        render_stage("gc_skew_data", double, 0.0)
        render_stage("gc_skew_data", double, -0.0)
    assert CALLS == ["count", "count", "double", "double"]


def test_the_settings_table_names_only_existing_sections() -> None:
    stages_module._validate_stage_settings(stages_module.RENDER_STAGE_SETTINGS)
    with pytest.raises(ValueError):
        stages_module._validate_stage_settings({"legend_table": frozenset({"labels.no_such_section"})})


# --- record identities (design 2.1, M1) ---------------------------------------


def _genbank(path: Path, *record_ids: str) -> Path:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    records = []
    for record_id in record_ids:
        record = SeqRecord(Seq("ATGC" * 150), id=record_id, name=record_id, description=record_id)
        record.annotations["molecule_type"] = "DNA"
        records.append(record)
    SeqIO.write(records, path, "genbank")
    return path


def _resolve(paths: tuple[Path, ...], *regions: str):
    from gbdraw.api import CircularDiagramRequest, GenBankInputSource, RecordInput
    from gbdraw.api.prepared import prepared_value_identity
    from gbdraw.api.request_render import normalize_request_records
    from gbdraw.api.requests import RecordCardinality, RecordCollectionOptions
    from gbdraw.io.regions import parse_region_spec

    request = CircularDiagramRequest(
        records=tuple(
            RecordInput(
                source=GenBankInputSource(path),
                cardinality=RecordCardinality.ALL,
                record_key=f"record-{index}",
            )
            for index, path in enumerate(paths)
        ),
        record_options=RecordCollectionOptions(
            regions=tuple(parse_region_spec(region) for region in regions)
        ),
    )
    identities = {
        path: PreparedResourceIdentity(path.name, path.name, path.stat().st_size) for path in paths
    }
    with PreparedBiologicalInputCache().transaction(resource_paths=identities, diagnostics=None):
        records = normalize_request_records(request)
        return [(record.id, len(record.seq), prepared_value_identity(record)) for record in records]


def test_a_file_scoped_region_crops_the_record_of_the_named_file(tmp_path: Path) -> None:
    """OV-700: the region key keeps its file selector."""

    a = _genbank(tmp_path / "a.gb", "a1")
    b = _genbank(tmp_path / "b.gb", "b1")
    cache = PreparedBiologicalInputCache()
    identities = {
        path: PreparedResourceIdentity(path.name, path.name, path.stat().st_size) for path in (a, b)
    }
    from gbdraw.api import CircularDiagramRequest, GenBankInputSource, RecordInput
    from gbdraw.api.request_render import normalize_request_records
    from gbdraw.api.requests import RecordCollectionOptions
    from gbdraw.io.regions import parse_region_spec

    def lengths(region: str) -> dict[str, int]:
        request = CircularDiagramRequest(
            records=(
                RecordInput(source=GenBankInputSource(a), record_key="a"),
                RecordInput(source=GenBankInputSource(b), record_key="b"),
            ),
            record_options=RecordCollectionOptions(regions=(parse_region_spec(region),)),
        )
        with cache.transaction(resource_paths=identities, diagnostics=None):
            return {record.id: len(record.seq) for record in normalize_request_records(request)}

    assert lengths("a.gb:#1:1-100") == {"a1": 100, "b1": 600}
    assert lengths("b.gb:#1:1-100") == {"a1": 600, "b1": 100}


def test_a_collection_transform_gives_records_the_identity_of_their_place(tmp_path: Path) -> None:
    """With collection options a record's identity depends on every input (review #4)."""

    one = _genbank(tmp_path / "one.gb", "a1")
    two = _genbank(tmp_path / "two.gb", "a1", "a2")
    b = _genbank(tmp_path / "b.gb", "b1")
    first = _resolve((one, b), "#2:1-100")
    second = _resolve((two, b), "#2:1-100")
    assert first[-1][:2] == ("b1", 100)
    assert second[-1][:2] == ("b1", 600)
    assert first[-1][2] != second[-1][2]


def test_the_default_colors_key_names_the_palette_of_an_overrides_file(tmp_path: Path) -> None:
    """The merged default colors depend on the palette, not only the file (review #6)."""

    from gbdraw.api import CircularDiagramRequest, GenBankInputSource, RecordInput
    from gbdraw.api.options import CircularDiagramOptions, ColorOptions
    from gbdraw.api.prepared import prepared_value_identity
    from gbdraw.api.request_render import plan_request

    source = _genbank(tmp_path / "a.gb", "a1")
    overrides = tmp_path / "default_colors.tsv"
    overrides.write_text("tRNA\t#123456\n", encoding="utf-8")
    identities = {
        path: PreparedResourceIdentity(path.name, path.name, path.stat().st_size)
        for path in (source, overrides)
    }
    cache = PreparedBiologicalInputCache()

    def resolve(palette: str) -> tuple[dict[str, str], object]:
        request = CircularDiagramRequest(
            records=(RecordInput(source=GenBankInputSource(source), record_key="a"),),
            options=CircularDiagramOptions(
                colors=ColorOptions(default_colors_file=str(overrides), default_colors_palette=palette)
            ),
        )
        with cache.transaction(resource_paths=identities, diagnostics=None):
            plan = plan_request(request)
            assert plan.inputs is not None
            features = plan.inputs.features
            return dict(features.default_color_map), prepared_value_identity(features.default_colors)

    default, default_identity = resolve("default")
    mint, mint_identity = resolve("mint")
    assert (default["CDS"], mint["CDS"], mint["tRNA"]) == ("#54bcf8", "#9cd9cf", "#123456")
    # A stage that takes the table keys it by this identity.
    assert default_identity is not None and mint_identity is not None
    assert default_identity != mint_identity


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_a_render_without_a_cache_compiles_the_label_maps_once_per_request(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    mode: str,
) -> None:
    """The label filter maps compile a fixed number of times, not once per record (review #13)."""

    import gbdraw.circular
    import gbdraw.labels.filtering as filtering
    import gbdraw.linear
    from Bio import SeqIO

    compile_maps = filtering.preprocess_label_filtering
    calls: list[int] = []

    def counting(label_filtering):
        calls.append(1)
        return compile_maps(label_filtering)

    monkeypatch.setattr(filtering, "preprocess_label_filtering", counting)

    def compiles(record_count: int) -> int:
        from Bio.Seq import Seq
        from Bio.SeqFeature import FeatureLocation, SeqFeature
        from Bio.SeqRecord import SeqRecord

        paths = []
        for index in range(record_count):
            record = SeqRecord(Seq("ATGC" * 300), id=f"r{index}", name=f"r{index}", description="r")
            record.annotations["molecule_type"] = "DNA"
            record.features = [
                SeqFeature(FeatureLocation(30, 300, strand=1), type="CDS", qualifiers={"product": ["p"]})
            ]
            path = tmp_path / f"{record_count}-{index}.gb"
            SeqIO.write([record], path, "genbank")
            paths.append(str(path))
        calls.clear()
        out = tmp_path / f"out-{len(list(tmp_path.iterdir()))}"
        if mode == "linear":
            gbdraw.linear.linear_main(["--gbk", *paths, "-o", str(out), "--show_labels", "all", "-f", "svg"])
        else:
            gbdraw.circular.circular_main(
                ["--gbk", paths[0], "-o", str(out), "--labels", "out", "-f", "svg"]
            )
        return len(calls)

    if mode == "circular":
        assert compiles(1) == 1
    else:
        assert compiles(3) == compiles(1)


def test_feature_inputs_report_the_color_table_error_before_the_visibility_error() -> None:
    """The error order of 3d748b78: color rules before visibility rules (review #8)."""

    import re

    from gbdraw.api.diagram import resolve_feature_inputs
    from gbdraw.exceptions import ValidationError

    color_table = pd.DataFrame(
        {"feature_type": ["CDS"], "qualifier_key": ["product"], "value": ["("], "color": ["red"], "caption": [""]}
    )
    visibility_table = pd.DataFrame(
        {"record_id": [""], "feature_type": ["CDS"], "qualifier": ["product"], "value": ["x"], "action": ["show"]}
    )
    with pytest.raises(re.error):
        resolve_feature_inputs(
            color_table=color_table,
            default_colors=None,
            feature_visibility_table=visibility_table,
        )
    with pytest.raises(ValidationError, match="record_id"):
        resolve_feature_inputs(
            color_table=None,
            default_colors=None,
            feature_visibility_table=visibility_table,
        )


def test_a_linear_render_compiles_the_color_table_only_in_its_stage(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The linear feature passes use the configurator's compiled maps (design 2.6)."""

    import gbdraw.features.colors as feature_colors
    import gbdraw.linear
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqFeature import FeatureLocation, SeqFeature
    from Bio.SeqRecord import SeqRecord

    recompiles: list[int] = []

    def counting(*args, **kwargs):
        recompiles.append(1)
        return compile_rules(*args, **kwargs)

    compile_rules = feature_colors.build_specific_color_rules
    monkeypatch.setattr(feature_colors, "build_specific_color_rules", counting)
    record = SeqRecord(Seq("ATGC" * 300), id="r1", name="r1", description="r")
    record.annotations["molecule_type"] = "DNA"
    record.features = [SeqFeature(FeatureLocation(30, 300, strand=1), type="CDS", qualifiers={"product": ["p"]})]
    source = tmp_path / "r1.gb"
    SeqIO.write([record], source, "genbank")
    table = tmp_path / "colors.tsv"
    table.write_text("CDS\tproduct\tp\tred\tP\n", encoding="utf-8")
    gbdraw.linear.linear_main(
        ["--gbk", str(source), "-t", str(table), "-o", str(tmp_path / "out"), "--show_labels", "all", "-f", "svg"]
    )
    assert recompiles == []


def measured_face(family: str) -> str | None:
    CALLS.append("measured_face")
    return measurement_font_file(family)


def _font_named(family: str, target: Path, source: str = "LiberationSans-Regular") -> Path:
    font = TTFont(str(Path(fonts.__file__).parents[1] / "data" / f"{source}.ttf"))
    name_table = font["name"]
    name_table.names = [record for record in name_table.names if record.nameID not in {1, 2, 4, 6, 16, 17}]
    for name_id, value in ((1, family), (2, "Regular"), (4, family), (6, family.replace(" ", "")), (16, family)):
        name_table.setName(value, name_id, 3, 1, 0x409)
    font.save(str(target))
    return target


def test_a_font_registered_between_renders_rebuilds_every_stage(tmp_path: Path) -> None:
    """D-11: the registered fonts are part of every stage key, so a cached measurement is not reused."""

    cache = PreparedBiologicalInputCache()
    family = "Render Stage Token Sans"
    try:
        with _transaction(cache):
            before = render_stage("legend_measurement", measured_face, family)
        fonts.register_font_file(_font_named(family, tmp_path / "token.ttf"))
        with _transaction(cache):
            after = render_stage("legend_measurement", measured_face, family)
        fonts.register_font_file(_font_named(family, tmp_path / "token.ttf", "LiberationSerif-Regular"))
        with _transaction(cache):
            replaced = render_stage("legend_measurement", measured_face, family)
    finally:
        fonts.clear_registered_fonts()
    assert CALLS == ["measured_face"] * 3
    assert len({before, after, replaced}) == 3
    assert after is not None and after.startswith("user-")


def catalog_part(value: int, name: str) -> dict[str, object]:
    CALLS.append(f"catalog_part:{value}")
    return {"value": value, "name": name}


def test_a_stage_revision_names_the_builds_of_one_cache_in_order() -> None:
    """Step 3a: a hit keeps its build number; a rebuild, even of an evicted equal entry, gets a new one."""

    cache = PreparedBiologicalInputCache()
    with _transaction(cache):
        first = [render_stage("catalog_item", catalog_part, index, "r") for index in (1, 2)]
        cold = stages_module.render_stage_revision("catalog_item", first)
        assert stages_module.render_stage_revision("catalog_item", list(reversed(first))) != cold
    assert cold is not None
    epoch, builds = cold.split(":")
    assert len(epoch) == 16 and int(epoch, 16) >= 0
    assert builds == "1.2"
    with _transaction(cache):
        warm = stages_module.render_stage_revision(
            "catalog_item", [render_stage("catalog_item", catalog_part, index, "r") for index in (1, 2)]
        )
    assert warm == cold
    with _transaction(cache):
        changed = stages_module.render_stage_revision(
            "catalog_item", [render_stage("catalog_item", catalog_part, index, "r") for index in (1, 3)]
        )
    assert changed == f"{epoch}:1.3"
    # The build of 2 was not used by the last commit, so it was dropped: the
    # same inputs build again under a new number.
    with _transaction(cache):
        again = stages_module.render_stage_revision(
            "catalog_item", [render_stage("catalog_item", catalog_part, index, "r") for index in (1, 2)]
        )
    assert again == f"{epoch}:1.4"
    assert CALLS == ["catalog_part:1", "catalog_part:2", "catalog_part:3", "catalog_part:2"]


def test_two_caches_name_the_same_inputs_differently() -> None:
    revisions = []
    for _cache in range(2):
        cache = PreparedBiologicalInputCache()
        with _transaction(cache):
            revisions.append(
                stages_module.render_stage_revision(
                    "catalog_item", [render_stage("catalog_item", catalog_part, 1, "r")]
                )
            )
    assert None not in revisions
    assert revisions[0] != revisions[1]


def test_a_stage_revision_is_none_unless_every_result_is_a_stored_build_of_that_stage() -> None:
    class Local:
        pass

    assert stages_module.render_stage_revision("catalog_item", [{"value": 1}]) is None
    cache = PreparedBiologicalInputCache()
    with _transaction(cache):
        stored = render_stage("catalog_item", catalog_part, 1, "r")
        unstored = render_stage("catalog_item", catalog_part, 2, Local())
        other_stage = render_stage("gc_skew_data", double, 4)
        assert stages_module.render_stage_revision("catalog_item", []) is None
        assert stages_module.render_stage_revision("catalog_item", [stored, unstored]) is None
        assert stages_module.render_stage_revision("catalog_item", [other_stage]) is None
        assert stages_module.render_stage_revision("catalog_item", [{"value": 1, "name": "r"}]) is None
        assert stages_module.render_stage_revision("catalog_item", [stored]) is not None
    assert stages_module.render_stage_revision("catalog_item", [stored]) is None


@pytest.mark.parametrize("error", [NotImplementedError, OSError])
def test_a_cache_without_a_random_source_names_no_build(
    monkeypatch: pytest.MonkeyPatch, error: type[Exception]
) -> None:
    """No epoch from a time or an id (D-38 OQ-5): the catalog is then always sent."""

    import gbdraw.api.prepared as prepared_module

    def no_random_source(_size: int) -> str:
        raise error("no random source")

    monkeypatch.setattr(prepared_module.secrets, "token_hex", no_random_source)
    cache = PreparedBiologicalInputCache()
    with _transaction(cache):
        result = render_stage("catalog_item", catalog_part, 1, "r")
        assert stages_module.render_stage_revision("catalog_item", [result]) is None
    with _transaction(cache):
        assert render_stage("catalog_item", catalog_part, 1, "r") is result


def _bookkeeping_transaction(cache: PreparedBiologicalInputCache, stage_args: tuple[int, ...], decoded: tuple[str, ...]):
    metrics: dict[str, object] = {}
    with _transaction(cache, metrics):
        for value in stage_args:
            render_stage("gc_skew_data", double, value)
        for key in decoded:
            cache.get_or_build("decoded", ("decoded", key), frozenset(), lambda key=key: {"value": key})
    return metrics


def test_render_stage_entries_are_not_owner_checked(monkeypatch: pytest.MonkeyPatch) -> None:
    """A render-stage entry has one constant stamp, so its owner check compared it with itself."""

    checked: list[str] = []
    real = PreparedBiologicalInputCache._verify_owner

    def spy(self, layer, *args, **kwargs):
        checked.append(layer)
        return real(self, layer, *args, **kwargs)

    monkeypatch.setattr(PreparedBiologicalInputCache, "_verify_owner", spy)
    cache = PreparedBiologicalInputCache()
    _bookkeeping_transaction(cache, (1, 2), ("x",))
    _bookkeeping_transaction(cache, (1, 2), ("x",))

    assert "render_stage" not in checked
    assert "decoded" in checked


def test_commit_counts_evictions_of_held_entries_only() -> None:
    import gbdraw.api.prepared as prepared_module

    cache = PreparedBiologicalInputCache()
    _bookkeeping_transaction(cache, (1, 2, 3), ("x", "y"))
    # One stage entry and one decoded entry are kept, one of each is new
    # (pending, not an eviction), two stage entries and one decoded entry go.
    metrics = _bookkeeping_transaction(cache, (1, 4), ("x", "z"))

    assert _metrics(metrics)["preparedInputCacheEvictionCount"] == 3
    assert not hasattr(prepared_module, "_layer_stamp")
