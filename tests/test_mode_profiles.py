from __future__ import annotations

from dataclasses import asdict
import math

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from svgwrite import Drawing

import gbdraw.api.diagram as api_diagram
import gbdraw.circular as circular_cli
import gbdraw.interface as interface
import gbdraw.linear as linear_cli
from gbdraw.api.options import (
    CircularDiagramOptions,
    CircularMultiRecordOptions,
    DiagramOptions,
    LinearDiagramOptions,
    LinearMultiRecordOptions,
    TrackOptions,
    resolve_diagram_options_for_mode,
)
from gbdraw.api.requests import (
    CircularDiagramRequest,
    InMemoryRecordSource,
    LinearDiagramRequest,
    RecordInput,
)
from gbdraw.exceptions import ValidationError
from gbdraw.mode_profiles import (
    CIRCULAR_MODE_PROFILE,
    DEFAULT_FEATURE_TYPES,
    LINEAR_MODE_PROFILE,
    ComparisonThresholds,
)
from gbdraw.tracks import CircularTrackSlot, LinearTrackSlot


EXPECTED_THRESHOLDS = {
    "circular": {
        "evalue": 1e-5,
        "bitscore": 50.0,
        "identity": 70.0,
        "alignment_length": 0,
    },
    "linear": {
        "evalue": 1e-2,
        "bitscore": 50.0,
        "identity": 0.0,
        "alignment_length": 0,
    },
}


def _record() -> SeqRecord:
    return SeqRecord(Seq("ATGC" * 25), id="record")


def _record_input() -> RecordInput:
    return RecordInput(source=InMemoryRecordSource(_record()))


def _threshold_dict(value: object) -> dict[str, float | int]:
    return {
        "evalue": getattr(value, "evalue"),
        "bitscore": getattr(value, "bitscore"),
        "identity": getattr(value, "identity"),
        "alignment_length": getattr(value, "alignment_length"),
    }


def test_mode_profiles_define_the_approved_defaults() -> None:
    assert asdict(CIRCULAR_MODE_PROFILE.comparison) == EXPECTED_THRESHOLDS["circular"]
    assert asdict(LINEAR_MODE_PROFILE.comparison) == EXPECTED_THRESHOLDS["linear"]
    assert (CIRCULAR_MODE_PROFILE.show_gc, CIRCULAR_MODE_PROFILE.show_skew) == (
        True,
        True,
    )
    assert (LINEAR_MODE_PROFILE.show_gc, LINEAR_MODE_PROFILE.show_skew) == (
        False,
        False,
    )
    assert CIRCULAR_MODE_PROFILE.feature_types == DEFAULT_FEATURE_TYPES
    assert LINEAR_MODE_PROFILE.feature_types == DEFAULT_FEATURE_TYPES
    assert "misc_RNA" in DEFAULT_FEATURE_TYPES
    assert LINEAR_MODE_PROFILE.linear_axis_color == "lightgray"
    assert LINEAR_MODE_PROFILE.linear_ruler_axis_color == "dimgray"


@pytest.mark.parametrize(
    ("mode", "root_options", "request_type", "cli_module"),
    [
        ("circular", interface.CircularOptions(), CircularDiagramRequest, circular_cli),
        ("linear", interface.LinearOptions(), LinearDiagramRequest, linear_cli),
    ],
)
def test_fresh_python_request_and_cli_threshold_defaults_match(
    mode: str,
    root_options: interface.CircularOptions | interface.LinearOptions,
    request_type: type[CircularDiagramRequest] | type[LinearDiagramRequest],
    cli_module: object,
) -> None:
    expected = EXPECTED_THRESHOLDS[mode]
    request = request_type(records=(_record_input(),))
    cli_args = cli_module._get_args(["--gbk", "record.gbk"])

    assert _threshold_dict(root_options.thresholds) == expected
    assert _threshold_dict(request.options) == expected
    assert _threshold_dict(cli_args) == expected
    assert tuple(root_options.features.types or ()) == DEFAULT_FEATURE_TYPES
    assert tuple(request.options.selected_features_set or ()) == DEFAULT_FEATURE_TYPES
    assert tuple(str(cli_args.features).split(",")) == DEFAULT_FEATURE_TYPES


def test_mode_resolution_preserves_explicit_values_and_applies_track_defaults() -> None:
    explicit = DiagramOptions(
        evalue=0.25,
        config_overrides={"show_gc": True, "linear_ruler_on_axis": True},
    )

    resolved = resolve_diagram_options_for_mode(explicit, mode="linear")

    assert resolved.evalue == 0.25
    assert resolved.bitscore == 50.0
    assert resolved.identity == 0.0
    assert resolved.config_overrides == {
        "show_gc": True,
        "show_skew": False,
        "linear_axis_stroke_color": "dimgray",
        "linear_ruler_on_axis": True,
    }


@pytest.mark.parametrize(
    ("options_type", "expected"),
    [
        (interface.CircularOptions, EXPECTED_THRESHOLDS["circular"]),
        (interface.LinearOptions, EXPECTED_THRESHOLDS["linear"]),
    ],
)
def test_partial_root_thresholds_resolve_omitted_values_for_the_parent_mode(
    options_type: type[interface.CircularOptions] | type[interface.LinearOptions],
    expected: dict[str, float | int],
) -> None:
    options = options_type(thresholds=interface.Thresholds(identity=88))

    assert _threshold_dict(options.thresholds) == {
        **expected,
        "identity": 88.0,
    }


def test_none_config_overrides_do_not_mask_linear_profile_defaults() -> None:
    resolved = resolve_diagram_options_for_mode(
        DiagramOptions(
            config_overrides={
                "show_gc": None,
                "linear_ruler_on_axis": True,
                "linear_axis_stroke_color": None,
            }
        ),
        mode="linear",
    )

    assert resolved.config_overrides == {
        "show_gc": False,
        "show_skew": False,
        "linear_axis_stroke_color": "dimgray",
        "linear_ruler_on_axis": True,
    }
    svg = api_diagram.build_linear_diagram(
        [_record()],
        options=DiagramOptions(config_overrides={"show_gc": None}),
    ).tostring()
    assert 'id="gc_content"' not in svg


def test_builders_forward_resolved_mode_profiles(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    captured: dict[str, dict[str, object]] = {}

    def fake_circular(_record: SeqRecord, **kwargs: object) -> Drawing:
        captured["circular"] = kwargs
        return Drawing("circular.svg")

    def fake_linear(_records: object, **kwargs: object) -> Drawing:
        captured["linear"] = kwargs
        return Drawing("linear.svg")

    monkeypatch.setattr(
        api_diagram,
        "assemble_circular_diagram_from_record",
        fake_circular,
    )
    monkeypatch.setattr(
        api_diagram,
        "assemble_linear_diagram_from_records",
        fake_linear,
    )

    api_diagram.build_circular_diagram(_record())
    api_diagram.build_linear_diagram([_record()])

    for mode in ("circular", "linear"):
        kwargs = captured[mode]
        assert _threshold_dict(type("Forwarded", (), kwargs)()) == EXPECTED_THRESHOLDS[
            mode
        ]
        assert kwargs["selected_features_set"] == DEFAULT_FEATURE_TYPES
    assert captured["circular"]["config_overrides"] == {
        "show_gc": True,
        "show_skew": True,
    }
    assert captured["linear"]["config_overrides"] == {
        "show_gc": False,
        "show_skew": False,
        "linear_axis_stroke_color": "lightgray",
    }


def test_paired_default_renders_apply_mode_track_and_axis_profiles() -> None:
    circular_svg = api_diagram.build_circular_diagram(_record()).tostring()
    linear_svg = api_diagram.build_linear_diagram([_record()]).tostring()

    assert 'id="gc_content"' in circular_svg
    assert 'id="skew"' in circular_svg
    assert 'id="gc_content"' not in linear_svg
    assert 'id="skew"' not in linear_svg
    assert 'stroke="lightgray"' in linear_svg


@pytest.mark.parametrize(
    "kwargs",
    [
        {"evalue": -1},
        {"evalue": math.nan},
        {"evalue": math.inf},
        {"evalue": "1e-5"},
        {"evalue": True},
        {"bitscore": -1},
        {"bitscore": math.nan},
        {"bitscore": math.inf},
        {"bitscore": "50"},
        {"identity": -1},
        {"identity": 101},
        {"identity": math.nan},
        {"identity": "70"},
        {"alignment_length": -1},
        {"alignment_length": 1.5},
        {"alignment_length": True},
    ],
)
def test_threshold_models_reject_invalid_public_domains(
    kwargs: dict[str, object],
) -> None:
    values: dict[str, object] = {
        "evalue": 1e-5,
        "bitscore": 50,
        "identity": 70,
        "alignment_length": 0,
    }
    values.update(kwargs)

    with pytest.raises(ValidationError):
        ComparisonThresholds(**values)
    with pytest.raises(ValidationError):
        interface.Thresholds(**values)
    with pytest.raises(ValidationError):
        DiagramOptions(**kwargs)


def test_threshold_overflow_is_normalized_to_validation_error() -> None:
    with pytest.raises(ValidationError):
        ComparisonThresholds(
            evalue=10**10000,
            bitscore=50,
            identity=70,
            alignment_length=0,
        )


@pytest.mark.parametrize(
    "factory",
    [
        lambda: DiagramOptions(colors="bad"),  # type: ignore[arg-type]
        lambda: DiagramOptions(tracks="bad"),  # type: ignore[arg-type]
        lambda: DiagramOptions(output="bad"),  # type: ignore[arg-type]
        lambda: DiagramOptions(annotations="bad"),  # type: ignore[arg-type]
        lambda: DiagramOptions(selected_features_set="CDS"),
        lambda: DiagramOptions(selected_features_set=("CDS", "")),
        lambda: DiagramOptions(config_overrides=[]),  # type: ignore[arg-type]
        lambda: DiagramOptions(feature_shapes=[]),  # type: ignore[arg-type]
    ],
)
def test_diagram_options_reject_invalid_nested_structure(factory: object) -> None:
    with pytest.raises(ValidationError):
        factory()


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_identity_100_is_valid_at_every_fresh_entry(mode: str) -> None:
    thresholds = ComparisonThresholds(
        evalue=1e-5,
        bitscore=50,
        identity=100,
        alignment_length=0,
    )
    request_type = CircularDiagramRequest if mode == "circular" else LinearDiagramRequest
    options_type = (
        CircularDiagramOptions if mode == "circular" else LinearDiagramOptions
    )
    cli_module = circular_cli if mode == "circular" else linear_cli

    request = request_type(
        records=(_record_input(),),
        options=options_type(identity=100),
    )
    cli_args = cli_module._get_args(["--gbk", "record.gbk", "--identity", "100"])

    assert thresholds.identity == 100
    assert request.options.identity == 100
    assert cli_args.identity == 100


@pytest.mark.parametrize(
    ("option", "value"),
    [
        ("evalue", "-1"),
        ("evalue", "nan"),
        ("bitscore", "-1"),
        ("bitscore", "inf"),
        ("identity", "-1"),
        ("identity", "101"),
        ("alignment_length", "-1"),
    ],
)
@pytest.mark.parametrize("cli_module", [circular_cli, linear_cli])
def test_cli_threshold_domain_errors_are_argparse_errors(
    cli_module: object,
    option: str,
    value: str,
) -> None:
    with pytest.raises(SystemExit, match="2"):
        cli_module._get_args(["--gbk", "record.gbk", f"--{option}", value])


def test_wrong_mode_nested_values_fail_during_construction() -> None:
    with pytest.raises(ValidationError, match="CircularTrackOptions"):
        interface.CircularOptions(tracks=interface.LinearTrackOptions())  # type: ignore[arg-type]
    with pytest.raises(ValidationError, match="LinearTrackOptions"):
        interface.LinearOptions(tracks=interface.CircularTrackOptions())  # type: ignore[arg-type]
    with pytest.raises(ValidationError, match="Circular title position"):
        interface.CircularOptions(title=interface.TitleOptions(position="center"))
    with pytest.raises(ValidationError, match="Linear title position"):
        interface.LinearOptions(title=interface.TitleOptions(position="none"))
    with pytest.raises(ValidationError, match="linear diagrams"):
        interface.CircularOptions(
            depth_tracks=(interface.DepthTrackOptions("depth.tsv", height=20),)
        )
    with pytest.raises(ValidationError, match="CircularDiagramOptions"):
        CircularDiagramRequest(
            records=(_record_input(),),
            options=LinearDiagramOptions(),  # type: ignore[arg-type]
        )
    with pytest.raises(ValidationError, match="CircularTrackSlot"):
        interface.CircularTrackOptions(
            slots=(LinearTrackSlot(id="wrong", renderer="features"),)
        )
    with pytest.raises(ValidationError, match="LinearTrackSlot"):
        interface.LinearTrackOptions(
            slots=(CircularTrackSlot(id="wrong", renderer="features"),)
        )
    with pytest.raises(ValidationError, match="CircularTrackSlot"):
        TrackOptions(
            circular_track_slots=(
                LinearTrackSlot(id="wrong", renderer="features"),
            )
        )


def test_root_facade_rejects_wrong_mode_layout_objects() -> None:
    records = (_record(), _record())

    with pytest.raises(ValidationError, match="CircularLayout"):
        interface.draw_circular(
            records,
            layout=interface.LinearLayout(),  # type: ignore[arg-type]
        )
    with pytest.raises(ValidationError, match="LinearLayout"):
        interface.draw_linear(
            records,
            layout=interface.CircularLayout(),  # type: ignore[arg-type]
        )
    with pytest.raises(ValidationError, match="CircularMultiRecordOptions"):
        api_diagram.build_circular_multi_diagram(
            records,
            layout=LinearMultiRecordOptions(),  # type: ignore[arg-type]
        )
    with pytest.raises(ValidationError, match="LinearMultiRecordOptions"):
        api_diagram.build_linear_diagram(
            records,
            layout=CircularMultiRecordOptions(),  # type: ignore[arg-type]
        )


def test_cli_track_defaults_are_sourced_from_mode_profiles() -> None:
    circular_args = circular_cli._get_args(["--gbk", "record.gbk"])
    linear_args = linear_cli._get_args(["--gbk", "record.gbk"])

    assert circular_args.show_gc is CIRCULAR_MODE_PROFILE.show_gc
    assert circular_args.show_skew is CIRCULAR_MODE_PROFILE.show_skew
    assert linear_args.show_gc is LINEAR_MODE_PROFILE.show_gc
    assert linear_args.show_skew is LINEAR_MODE_PROFILE.show_skew
    for args in (circular_args, linear_args):
        assert not hasattr(args, "depth")
        assert not hasattr(args, "show_depth")
        assert not hasattr(args, "suppress_gc")
        assert not hasattr(args, "suppress_skew")


@pytest.mark.parametrize(
    ("cli_module", "profile"),
    [
        (circular_cli, CIRCULAR_MODE_PROFILE),
        (linear_cli, LINEAR_MODE_PROFILE),
    ],
)
def test_cli_gc_and_skew_use_symmetric_boolean_pairs(cli_module, profile) -> None:
    base_args = ["--gbk", "record.gbk"]

    assert cli_module._get_args([*base_args, "--gc"]).show_gc is True
    assert cli_module._get_args([*base_args, "--no-gc"]).show_gc is False
    assert cli_module._get_args([*base_args, "--skew"]).show_skew is True
    assert cli_module._get_args([*base_args, "--no-skew"]).show_skew is False
    assert cli_module._get_args(base_args).show_gc is profile.show_gc
    assert cli_module._get_args(base_args).show_skew is profile.show_skew

    with pytest.raises(SystemExit):
        cli_module._get_args([*base_args, "--gc", "--no-gc"])
    with pytest.raises(SystemExit):
        cli_module._get_args([*base_args, "--skew", "--no-skew"])


@pytest.mark.parametrize(
    ("cli_module", "legacy_args"),
    [
        (circular_cli, ["--suppress_gc"]),
        (circular_cli, ["--suppress_skew"]),
        (circular_cli, ["--depth", "depth.tsv"]),
        (circular_cli, ["--show_depth"]),
        (linear_cli, ["--show_gc"]),
        (linear_cli, ["--show_skew"]),
        (linear_cli, ["--depth", "depth.tsv"]),
        (linear_cli, ["--show_depth"]),
    ],
)
def test_fresh_cli_rejects_removed_track_flags(cli_module, legacy_args) -> None:
    with pytest.raises(SystemExit):
        cli_module._get_args(["--gbk", "record.gbk", *legacy_args])


@pytest.mark.parametrize(
    "layout",
    [
        lambda: interface.CircularLayout(min_radius_ratio=math.nan),
        lambda: interface.CircularLayout(column_gap_ratio=-1),
        lambda: interface.CircularLayout(row_gap_ratio=True),
        lambda: interface.CircularLayout(positions=123),  # type: ignore[arg-type]
        lambda: interface.LinearLayout(record_gap=math.inf),
        lambda: interface.LinearLayout(record_gap=-1),
        lambda: interface.LinearLayout(record_gap="24"),
        lambda: interface.LinearLayout(positions=123),  # type: ignore[arg-type]
        lambda: CircularMultiRecordOptions(
            multi_record_positions=123,  # type: ignore[arg-type]
        ),
        lambda: LinearMultiRecordOptions(
            multi_record_positions=123,  # type: ignore[arg-type]
        ),
    ],
)
def test_root_layouts_validate_eagerly(layout: object) -> None:
    with pytest.raises(ValidationError):
        layout()
