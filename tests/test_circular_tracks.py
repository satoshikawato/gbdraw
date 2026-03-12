from __future__ import annotations

import json
from pathlib import Path

import pytest
from Bio import SeqIO
from svgwrite import Drawing

import gbdraw.circular as circular_cli_module
import gbdraw.diagrams.circular.assemble as circular_assemble_module
from gbdraw.api.diagram import assemble_circular_diagram_from_record
from gbdraw.config.models import GbdrawConfig
from gbdraw.config.modify import modify_config_dict
from gbdraw.config.toml import load_config_toml
from gbdraw.exceptions import ValidationError
from gbdraw.features.objects import FeatureLocationPart, FeatureObject
from gbdraw.tracks import (
    CircularTrackPlacement,
    ScalarSpec,
    TrackSpec,
    iter_compiled_custom_rules,
    load_circular_track_specs,
    normalize_circular_track_specs,
)
from gbdraw.diagrams.circular.track_manager import (
    assign_features_to_feature_tracks,
    build_feature_track_entries,
    build_ordered_annular_track_specs,
)


SELECTED_FEATURES = ["CDS", "rRNA", "tRNA", "tmRNA", "ncRNA", "misc_RNA", "repeat_region"]


def _load_record():
    input_path = Path(__file__).parent / "test_inputs" / "HmmtDNA.gbk"
    return SeqIO.read(str(input_path), "genbank")


def _make_config_dict(*, show_labels: bool) -> dict:
    config_dict = load_config_toml("gbdraw.data", "config.toml")
    return modify_config_dict(
        config_dict,
        show_labels=show_labels,
        allow_inner_labels=False,
        resolve_overlaps=False,
        strandedness=True,
        track_type="tuckin",
    )


def _make_feature(
    feature_id: str,
    feature_type: str,
    strand: str,
    start: int,
    end: int,
    *,
    qualifiers: dict[str, object] | None = None,
) -> FeatureObject:
    location = [FeatureLocationPart("block", "001", strand, start, end, True)]
    return FeatureObject(
        feature_id=feature_id,
        location=location,
        is_directional=True,
        color="#54bcf8",
        note="",
        label_text="",
        coordinates=location,
        type=feature_type,
        qualifiers=qualifiers or {},
    )


def _write_track_file(tmp_path: Path, payload: list[dict]) -> Path:
    path = tmp_path / "tracks.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path


def _capture_annular_annuli(track_specs: list[TrackSpec]) -> dict[str, tuple[float, float]]:
    record = _load_record()
    config_dict = _make_config_dict(show_labels=False)
    cfg_model = GbdrawConfig.from_dict(config_dict)
    captured: dict[str, tuple[float, float]] = {}

    def fake_add_record_group_on_canvas(
        canvas,
        gb_record,
        canvas_config,
        feature_config,
        config_dict,
        *,
        cfg=None,
        precomputed_feature_dict=None,
        precalculated_labels=None,
        feature_track_ratio_factor_override=None,
        group_id=None,
        radius_override=None,
    ):
        if precomputed_feature_dict:
            feature_ratio_factor = (
                float(feature_track_ratio_factor_override)
                if feature_track_ratio_factor_override is not None
                else float(cfg_model.canvas.circular.track_ratio_factors[str(canvas_config.length_param)][0])
            )
            bounds = circular_assemble_module._compute_feature_band_bounds_px(
                precomputed_feature_dict,
                len(gb_record.seq),
                base_radius_px=float(canvas_config.radius),
                track_ratio=float(canvas_config.track_ratio),
                length_param=str(canvas_config.length_param),
                track_ratio_factor=feature_ratio_factor,
                cfg=cfg or canvas_config._cfg,
            )
            if bounds is not None:
                captured["features"] = bounds
        return canvas

    def fake_add_gc_content_group_on_canvas(
        canvas,
        gb_record,
        gc_df,
        canvas_config,
        gc_config,
        config_dict,
        *,
        track_width_override=None,
        norm_factor_override=None,
        group_id=None,
        track_id_override=None,
        extra_attribs=None,
        cfg=None,
    ):
        width = (
            float(track_width_override)
            if track_width_override is not None
            else (
                float(canvas_config.radius)
                * float(canvas_config.track_ratio)
                * float((cfg or canvas_config._cfg).canvas.circular.track_ratio_factors[str(canvas_config.length_param)][1])
            )
        )
        center = (
            float(norm_factor_override) * float(canvas_config.radius)
            if norm_factor_override is not None
            else float(
                circular_assemble_module._default_track_center_radius_px(
                    "gc_content",
                    canvas_config=canvas_config,
                    cfg=cfg or canvas_config._cfg,
                )
            )
        )
        track_key = (
            str((extra_attribs or {}).get("data-track-id") or "")
            or (str(group_id)[6:] if str(group_id).startswith("track_") else "")
            or "gc_content"
        )
        captured[track_key] = circular_assemble_module._annulus_from_center_and_width(center, width)
        return canvas

    def fake_add_gc_skew_group_on_canvas(
        canvas,
        gb_record,
        gc_df,
        canvas_config,
        skew_config,
        config_dict,
        *,
        track_width_override=None,
        norm_factor_override=None,
        group_id=None,
        track_id_override=None,
        extra_attribs=None,
        cfg=None,
    ):
        width = (
            float(track_width_override)
            if track_width_override is not None
            else (
                float(canvas_config.radius)
                * float(canvas_config.track_ratio)
                * float((cfg or canvas_config._cfg).canvas.circular.track_ratio_factors[str(canvas_config.length_param)][2])
            )
        )
        center = (
            float(norm_factor_override) * float(canvas_config.radius)
            if norm_factor_override is not None
            else float(
                circular_assemble_module._default_track_center_radius_px(
                    "gc_skew",
                    canvas_config=canvas_config,
                    cfg=cfg or canvas_config._cfg,
                )
            )
        )
        track_key = (
            str((extra_attribs or {}).get("data-track-id") or "")
            or (str(group_id)[6:] if str(group_id).startswith("track_") else "")
            or "gc_skew"
        )
        captured[track_key] = circular_assemble_module._annulus_from_center_and_width(center, width)
        return canvas

    original_record = circular_assemble_module.add_record_group_on_canvas
    original_gc = circular_assemble_module.add_gc_content_group_on_canvas
    original_skew = circular_assemble_module.add_gc_skew_group_on_canvas
    original_axis = circular_assemble_module.add_axis_group_on_canvas
    original_labels = circular_assemble_module.add_labels_group_on_canvas
    original_definition = circular_assemble_module.add_record_definition_group_on_canvas
    original_legend = circular_assemble_module.add_legend_group_on_canvas
    try:
        circular_assemble_module.add_record_group_on_canvas = fake_add_record_group_on_canvas
        circular_assemble_module.add_gc_content_group_on_canvas = fake_add_gc_content_group_on_canvas
        circular_assemble_module.add_gc_skew_group_on_canvas = fake_add_gc_skew_group_on_canvas
        circular_assemble_module.add_axis_group_on_canvas = lambda canvas, *args, **kwargs: canvas
        circular_assemble_module.add_labels_group_on_canvas = lambda canvas, *args, **kwargs: canvas
        circular_assemble_module.add_record_definition_group_on_canvas = (
            lambda canvas, gb_record, canvas_config, species, strain, config_dict, **kwargs: canvas
        )
        circular_assemble_module.add_legend_group_on_canvas = (
            lambda canvas, canvas_config, legend_config, legend_table: canvas
        )
        assemble_circular_diagram_from_record(
            record,
            config_dict=config_dict,
            selected_features_set=SELECTED_FEATURES,
            legend="none",
            track_specs=track_specs,
        )
    finally:
        circular_assemble_module.add_record_group_on_canvas = original_record
        circular_assemble_module.add_gc_content_group_on_canvas = original_gc
        circular_assemble_module.add_gc_skew_group_on_canvas = original_skew
        circular_assemble_module.add_axis_group_on_canvas = original_axis
        circular_assemble_module.add_labels_group_on_canvas = original_labels
        circular_assemble_module.add_record_definition_group_on_canvas = original_definition
        circular_assemble_module.add_legend_group_on_canvas = original_legend

    return captured


def _order_annuli_outer_to_inner(annuli: dict[str, tuple[float, float]]) -> list[str]:
    return [
        name
        for name, _bounds in sorted(
            annuli.items(),
            key=lambda item: float(item[1][1]),
            reverse=True,
        )
    ]


def test_load_circular_track_specs_accepts_valid_config(tmp_path: Path) -> None:
    track_file = _write_track_file(
        tmp_path,
        [
            {
                "id": "features",
                "kind": "features",
                "show": True,
                "placement": {"width": "48px"},
                "params": {"feature_types": ["CDS", "tRNA"]},
            },
            {
                "id": "membrane",
                "kind": "custom",
                "show": True,
                "placement": {"radius": "0.91", "width": "16px"},
                "params": {
                    "caption": "Membrane",
                    "feature_types": ["CDS"],
                    "strand_mode": "positive",
                    "rules": [{"qualifier": "product", "pattern": "membrane"}],
                },
            },
            {
                "id": "at_content",
                "kind": "analysis",
                "show": True,
                "placement": {"radius": "0.70", "width": "18px"},
                "params": {"caption": "AT content", "metric": "content", "dinucleotide": "AT"},
            },
            {
                "id": "gc_content",
                "kind": "gc_content",
                "show": False,
                "placement": {"radius": "0.74", "width": "22px"},
            },
        ],
    )

    specs = load_circular_track_specs(track_file)

    assert [spec.id for spec in specs] == ["features", "membrane", "at_content", "gc_content"]
    assert specs[0].params == {"feature_types": ["CDS", "tRNA"]}
    assert specs[0].placement is not None
    assert specs[0].placement.width is not None
    assert specs[0].placement.width.resolve(400.0) == pytest.approx(48.0)
    assert specs[1].params == {
        "caption": "Membrane",
        "feature_types": ["CDS"],
        "strand_mode": "positive",
        "rules": [{"qualifier": "product", "pattern": "membrane"}],
        "match_all": True,
    }
    assert specs[1].placement is not None
    assert specs[1].placement.width is not None
    assert specs[1].placement.width.resolve(400.0) == pytest.approx(16.0)
    compiled_rules = iter_compiled_custom_rules(specs[1])
    assert len(compiled_rules) == 1
    assert compiled_rules[0].regex.search("Inner membrane protein")
    assert specs[2].params == {
        "caption": "AT content",
        "metric": "content",
        "dinucleotide": "AT",
    }
    assert specs[2].placement is not None
    assert specs[2].placement.radius is not None
    assert specs[2].placement.radius.resolve(400.0) == pytest.approx(280.0)


def test_load_circular_track_specs_preserves_blank_analysis_caption(tmp_path: Path) -> None:
    track_file = _write_track_file(
        tmp_path,
        [
            {
                "id": "at_content",
                "kind": "analysis",
                "show": True,
                "params": {"caption": "", "metric": "content", "dinucleotide": "AT"},
            }
        ],
    )

    specs = load_circular_track_specs(track_file)

    assert len(specs) == 1
    assert specs[0].params == {
        "caption": "",
        "metric": "content",
        "dinucleotide": "AT",
    }


@pytest.mark.parametrize(
    ("payload", "message_fragment"),
    [
        (
            [
                {"id": "features", "kind": "features", "show": True},
                {"id": "features", "kind": "custom", "show": True, "params": {"caption": "dup", "feature_types": ["CDS"], "strand_mode": "all", "rules": []}},
            ],
            "Duplicate circular track id",
        ),
        (
            [
                {"id": "features_a", "kind": "features", "show": True},
                {"id": "features_b", "kind": "features", "show": True},
            ],
            "Duplicate circular built-in track kind",
        ),
        (
            [{"id": "bad", "kind": "axis", "show": True}],
            "unsupported kind",
        ),
        (
            [{"id": "bad", "kind": "custom", "show": True, "params": {"caption": "Bad", "feature_types": ["CDS"], "strand_mode": "sideways", "rules": []}}],
            "strand_mode",
        ),
        (
            [{"id": "bad", "kind": "custom", "show": True, "params": {"caption": "Bad", "feature_types": [], "strand_mode": "all", "rules": []}}],
            "feature_types",
        ),
        (
            [{"id": "bad", "kind": "analysis", "show": True, "params": {"caption": "Bad", "metric": "noise", "dinucleotide": "AT"}}],
            "params.metric",
        ),
        (
            [{"id": "bad", "kind": "analysis", "show": True, "params": {"caption": "Bad", "metric": "content", "dinucleotide": "A"}}],
            "params.dinucleotide",
        ),
        (
            [{"id": "bad", "kind": "custom", "show": True, "params": {"caption": "Bad", "feature_types": ["CDS"], "strand_mode": "all", "rules": [{"qualifier": "product", "pattern": "["}]}}],
            "Invalid regex",
        ),
    ],
)
def test_load_circular_track_specs_rejects_invalid_payloads(
    tmp_path: Path,
    payload: list[dict],
    message_fragment: str,
) -> None:
    track_file = _write_track_file(tmp_path, payload)

    with pytest.raises(ValidationError, match=message_fragment):
        load_circular_track_specs(track_file)


def test_normalize_circular_track_specs_accepts_direct_custom_trackspec() -> None:
    specs = normalize_circular_track_specs(
        [
            TrackSpec(
                id="membrane",
                kind="custom",
                mode="circular",
                show=True,
                params={
                    "caption": "Membrane",
                    "feature_types": ["CDS"],
                    "strand_mode": "all",
                    "rules": [{"qualifier": "product", "pattern": "membrane"}],
                },
            )
        ]
    )

    assert len(specs) == 1
    assert specs[0].params is not None
    assert specs[0].params["match_all"] is True
    assert iter_compiled_custom_rules(specs[0])[0].regex.search("membrane protein")


def test_assign_features_to_feature_tracks_first_match_wins_and_fallback() -> None:
    base_feature_dict = {
        "f1": _make_feature("f1", "CDS", "positive", 10, 80, qualifiers={"product": ["membrane protein"]}),
        "f2": _make_feature("f2", "CDS", "positive", 100, 160, qualifiers={"product": ["enzyme"]}),
        "f3": _make_feature("f3", "tRNA", "positive", 180, 220),
    }
    annular_specs = build_ordered_annular_track_specs(
        [
            TrackSpec(
                id="membrane",
                kind="custom",
                mode="circular",
                show=True,
                params={
                    "caption": "Membrane",
                    "feature_types": ["CDS"],
                    "strand_mode": "all",
                    "rules": [{"qualifier": "product", "pattern": "membrane"}],
                },
            ),
            TrackSpec(
                id="cds_rest",
                kind="custom",
                mode="circular",
                show=True,
                params={
                    "caption": "Other CDS",
                    "feature_types": ["CDS"],
                    "strand_mode": "all",
                    "rules": [],
                },
            ),
            TrackSpec(
                id="features",
                kind="features",
                mode="circular",
                show=True,
                params={"feature_types": ["CDS", "tRNA"]},
            ),
        ],
        show_gc=True,
        show_skew=True,
    )
    feature_track_entries = build_feature_track_entries(
        annular_specs,
        fallback_feature_types=SELECTED_FEATURES,
    )

    rendered_track_dicts = assign_features_to_feature_tracks(
        base_feature_dict,
        feature_track_entries,
        separate_strands=True,
        resolve_overlaps=False,
        split_overlaps_by_strand=False,
        genome_length=500,
    )

    assert set(rendered_track_dicts["membrane"]) == {"f1"}
    assert set(rendered_track_dicts["cds_rest"]) == {"f2"}
    assert set(rendered_track_dicts["features"]) == {"f3"}
    assert base_feature_dict["f1"].circular_track_id == "membrane"
    assert base_feature_dict["f2"].circular_track_id == "cds_rest"
    assert base_feature_dict["f3"].circular_track_id == "features"


def test_assign_features_to_feature_tracks_hidden_builtin_hides_unclaimed() -> None:
    base_feature_dict = {
        "f1": _make_feature("f1", "CDS", "positive", 10, 80, qualifiers={"product": ["membrane protein"]}),
        "f2": _make_feature("f2", "tRNA", "positive", 100, 160),
    }
    annular_specs = build_ordered_annular_track_specs(
        [
            TrackSpec(
                id="membrane",
                kind="custom",
                mode="circular",
                show=True,
                params={
                    "caption": "Membrane",
                    "feature_types": ["CDS"],
                    "strand_mode": "all",
                    "rules": [{"qualifier": "product", "pattern": "membrane"}],
                },
            ),
            TrackSpec(
                id="features",
                kind="features",
                mode="circular",
                show=False,
                params={"feature_types": ["CDS", "tRNA"]},
            ),
        ],
        show_gc=True,
        show_skew=True,
    )
    feature_track_entries = build_feature_track_entries(
        annular_specs,
        fallback_feature_types=SELECTED_FEATURES,
    )

    rendered_track_dicts = assign_features_to_feature_tracks(
        base_feature_dict,
        feature_track_entries,
        separate_strands=True,
        resolve_overlaps=False,
        split_overlaps_by_strand=False,
        genome_length=500,
    )

    assert set(rendered_track_dicts["membrane"]) == {"f1"}
    assert "features" not in rendered_track_dicts or rendered_track_dicts["features"] == {}
    assert base_feature_dict["f2"].circular_track_id is None


def test_cli_track_file_precedence_ignores_conflicting_builtin_shorthands(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
) -> None:
    record = _load_record()
    track_file = _write_track_file(
        tmp_path,
        [
            {
                "id": "features",
                "kind": "features",
                "show": True,
                "placement": {"width": "55px"},
            },
            {
                "id": "gc_content",
                "kind": "gc_content",
                "show": True,
                "placement": {"radius": "0.72", "width": "22px"},
            },
        ],
    )
    captured: dict[str, object] = {}

    monkeypatch.setattr(circular_cli_module, "load_gbks", lambda paths, mode: [record])
    monkeypatch.setattr(circular_cli_module, "read_color_table", lambda _path: None)
    monkeypatch.setattr(circular_cli_module, "load_default_colors", lambda _path, _palette: None)
    monkeypatch.setattr(circular_cli_module, "save_figure", lambda canvas, formats: None)

    def fake_assemble(*args, **kwargs):
        captured["track_specs"] = kwargs.get("track_specs")
        return Drawing(filename=str(tmp_path / "dummy.svg"))

    monkeypatch.setattr(circular_cli_module, "assemble_circular_diagram_from_record", fake_assemble)

    with caplog.at_level("WARNING"):
        circular_cli_module.circular_main(
            [
                "--gbk",
                "dummy.gb",
                "--track_file",
                str(track_file),
                "--feature_width",
                "42",
                "--gc_content_width",
                "18",
                "--gc_content_radius",
                "0.74",
                "--format",
                "svg",
                "-o",
                str(tmp_path / "out"),
            ]
        )

    track_specs = captured["track_specs"]
    assert isinstance(track_specs, list)
    assert len(track_specs) == 2
    assert [track_spec.id for track_spec in track_specs] == ["features", "gc_content"]
    assert any("Ignoring --feature_width" in message for message in caplog.messages)
    assert any("Ignoring --gc_content_width/--gc_content_radius" in message for message in caplog.messages)


def test_assemble_circular_diagram_renders_multiple_feature_track_groups() -> None:
    record = _load_record()
    config_dict = _make_config_dict(show_labels=False)

    custom_tracks = [
        TrackSpec(
            id="protein_track",
            kind="custom",
            mode="circular",
            show=True,
            params={
                "caption": "Proteins",
                "feature_types": ["CDS"],
                "strand_mode": "all",
                "rules": [],
            },
        ),
        TrackSpec(
            id="trna_track",
            kind="custom",
            mode="circular",
            show=True,
            params={
                "caption": "tRNAs",
                "feature_types": ["tRNA"],
                "strand_mode": "all",
                "rules": [],
            },
        ),
        TrackSpec(
            id="features",
            kind="features",
            mode="circular",
            show=True,
            params={"feature_types": ["rRNA", "tmRNA", "ncRNA", "misc_RNA", "repeat_region"]},
        ),
    ]

    canvas = assemble_circular_diagram_from_record(
        record,
        config_dict=config_dict,
        selected_features_set=SELECTED_FEATURES,
        legend="none",
        track_specs=custom_tracks,
    )
    svg_text = canvas.tostring()

    assert 'id="track_protein_track"' in svg_text
    assert 'id="track_trna_track"' in svg_text
    assert 'id="track_features"' in svg_text
    assert svg_text.index('id="track_protein_track"') < svg_text.index('id="track_trna_track"')
    assert svg_text.index('id="track_trna_track"') < svg_text.index('id="track_features"')


def test_assemble_circular_diagram_renders_extra_analysis_track_groups() -> None:
    record = _load_record()
    config_dict = _make_config_dict(show_labels=False)

    track_specs = [
        TrackSpec(id="features", kind="features", mode="circular", show=True),
        TrackSpec(
            id="at_content",
            kind="analysis",
            mode="circular",
            show=True,
            placement=CircularTrackPlacement(
                radius=ScalarSpec(0.74, "factor"),
                width=ScalarSpec(18.0, "px"),
            ),
            params={"caption": "AT content", "metric": "content", "dinucleotide": "AT"},
        ),
        TrackSpec(
            id="at_skew",
            kind="analysis",
            mode="circular",
            show=True,
            placement=CircularTrackPlacement(
                radius=ScalarSpec(0.66, "factor"),
                width=ScalarSpec(18.0, "px"),
            ),
            params={"caption": "AT skew", "metric": "skew", "dinucleotide": "AT"},
        ),
    ]

    canvas = assemble_circular_diagram_from_record(
        record,
        config_dict=config_dict,
        selected_features_set=SELECTED_FEATURES,
        legend="right",
        track_specs=track_specs,
    )
    svg_text = canvas.tostring()

    assert 'id="track_at_content"' in svg_text
    assert 'id="track_at_skew"' in svg_text
    assert 'data-track-id="at_content"' in svg_text
    assert 'data-track-id="at_skew"' in svg_text
    assert 'data-track-metric="content"' in svg_text
    assert 'data-track-metric="skew"' in svg_text
    assert 'AT content' in svg_text
    assert 'AT skew (+)' in svg_text
    assert 'AT skew (-)' in svg_text


def test_assemble_circular_diagram_preserves_blank_analysis_caption_without_legend_entry() -> None:
    record = _load_record()
    config_dict = _make_config_dict(show_labels=False)

    track_specs = [
        TrackSpec(id="features", kind="features", mode="circular", show=True),
        TrackSpec(
            id="at_content",
            kind="analysis",
            mode="circular",
            show=True,
            placement=CircularTrackPlacement(
                radius=ScalarSpec(0.74, "factor"),
                width=ScalarSpec(18.0, "px"),
            ),
            params={"caption": "", "metric": "content", "dinucleotide": "AT"},
        ),
        TrackSpec(id="gc_content", kind="gc_content", mode="circular", show=False),
        TrackSpec(id="gc_skew", kind="gc_skew", mode="circular", show=False),
    ]

    canvas = assemble_circular_diagram_from_record(
        record,
        config_dict=config_dict,
        selected_features_set=SELECTED_FEATURES,
        legend="right",
        track_specs=track_specs,
    )
    svg_text = canvas.tostring()

    assert 'id="track_at_content"' in svg_text
    assert 'data-track-id="at_content"' in svg_text
    assert 'data-track-metric="content"' in svg_text
    assert 'AT content' not in svg_text
    assert 'AT content (+)' not in svg_text
    assert 'AT content (-)' not in svg_text


def test_ordered_annular_layout_respects_builtin_track_order() -> None:
    default_annuli = _capture_annular_annuli(
        [
            TrackSpec(id="features", kind="features", mode="circular", show=True),
            TrackSpec(id="gc_content", kind="gc_content", mode="circular", show=True),
            TrackSpec(id="gc_skew", kind="gc_skew", mode="circular", show=True),
        ]
    )
    reordered_annuli = _capture_annular_annuli(
        [
            TrackSpec(id="gc_skew", kind="gc_skew", mode="circular", show=True),
            TrackSpec(id="gc_content", kind="gc_content", mode="circular", show=True),
            TrackSpec(id="features", kind="features", mode="circular", show=True),
        ]
    )

    assert _order_annuli_outer_to_inner(default_annuli) == ["features", "gc_content", "gc_skew"]
    assert _order_annuli_outer_to_inner(reordered_annuli) == ["gc_skew", "gc_content", "features"]


def test_ordered_annular_layout_keeps_explicit_gc_radius_when_reordered() -> None:
    annuli = _capture_annular_annuli(
        [
            TrackSpec(
                id="gc_content",
                kind="gc_content",
                mode="circular",
                show=True,
                placement=CircularTrackPlacement(
                    radius=ScalarSpec(0.92, "factor"),
                    width=ScalarSpec(22.0, "px"),
                ),
            ),
            TrackSpec(id="gc_skew", kind="gc_skew", mode="circular", show=True),
            TrackSpec(id="features", kind="features", mode="circular", show=True),
        ]
    )

    cfg = GbdrawConfig.from_dict(_make_config_dict(show_labels=False))
    gc_inner, gc_outer = annuli["gc_content"]
    gc_center = 0.5 * (float(gc_inner) + float(gc_outer))

    assert gc_center == pytest.approx(0.92 * float(cfg.canvas.circular.radius))
    assert _order_annuli_outer_to_inner(annuli) == ["gc_content", "gc_skew", "features"]


def test_ordered_annular_layout_skips_hidden_builtin_rows() -> None:
    annuli = _capture_annular_annuli(
        [
            TrackSpec(id="gc_skew", kind="gc_skew", mode="circular", show=False),
            TrackSpec(id="gc_content", kind="gc_content", mode="circular", show=True),
            TrackSpec(id="features", kind="features", mode="circular", show=True),
        ]
    )

    assert set(annuli) == {"features", "gc_content"}
    assert _order_annuli_outer_to_inner(annuli) == ["gc_content", "features"]


def test_ordered_annular_layout_moves_gc_inside_multiple_feature_annuli(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    record = _load_record()
    config_dict = load_config_toml("gbdraw.data", "config.toml")
    config_dict = modify_config_dict(
        config_dict,
        show_labels=False,
        allow_inner_labels=False,
        resolve_overlaps=True,
        strandedness=False,
        track_type="tuckin",
    )
    captured_gc_norms: list[float | None] = []

    def fake_add_gc_content_group_on_canvas(
        canvas,
        gb_record,
        gc_df,
        canvas_config,
        gc_config,
        config_dict,
        *,
        track_width_override=None,
        norm_factor_override=None,
        radius_override=None,
        group_id=None,
        track_id_override=None,
        extra_attribs=None,
        cfg=None,
    ):
        captured_gc_norms.append(norm_factor_override)
        return canvas

    monkeypatch.setattr(
        circular_assemble_module,
        "add_gc_content_group_on_canvas",
        fake_add_gc_content_group_on_canvas,
    )

    assemble_circular_diagram_from_record(
        record,
        config_dict=config_dict,
        selected_features_set=SELECTED_FEATURES,
        legend="none",
        track_specs=["features@w=96px"],
    )
    assemble_circular_diagram_from_record(
        record,
        config_dict=config_dict,
        selected_features_set=SELECTED_FEATURES,
        legend="none",
        track_specs=[
            TrackSpec(
                id="protein_track",
                kind="custom",
                mode="circular",
                show=True,
                params={
                    "caption": "Proteins",
                    "feature_types": ["CDS"],
                    "strand_mode": "all",
                    "rules": [],
                },
            ),
            TrackSpec(
                id="features",
                kind="features",
                mode="circular",
                show=True,
                placement=CircularTrackPlacement(width=ScalarSpec(96.0, "px")),
                params={"feature_types": ["rRNA", "tRNA", "tmRNA", "ncRNA", "misc_RNA", "repeat_region"]},
            ),
        ],
    )

    assert len(captured_gc_norms) == 2
    assert captured_gc_norms[0] is not None
    assert captured_gc_norms[1] is not None
    assert float(captured_gc_norms[1]) < float(captured_gc_norms[0])
