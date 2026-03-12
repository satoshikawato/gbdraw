from __future__ import annotations

import json
import math
from pathlib import Path
import xml.etree.ElementTree as ET

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
        extra_attribs=None,
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


_SVG_NS = {"svg": "http://www.w3.org/2000/svg"}


def _get_group_attr_float(group: ET.Element | None, attr_name: str) -> float | None:
    if group is None:
        return None
    raw = str(group.attrib.get(attr_name, "")).strip()
    if not raw:
        return None
    return float(raw)


def _find_track_group(root: ET.Element, track_id: str) -> ET.Element | None:
    return (
        root.find(f".//svg:g[@data-track-id='{track_id}']", _SVG_NS)
        or root.find(f".//svg:g[@id='track_{track_id}']", _SVG_NS)
    )


def _make_definition_clearance_track_specs(width_px: float | None = None) -> list[TrackSpec]:
    placement_kwargs: dict[str, ScalarSpec] = {
        "radius": ScalarSpec(0.345, "factor"),
    }
    if width_px is not None:
        placement_kwargs["width"] = ScalarSpec(float(width_px), "px")

    return [
        TrackSpec(
            id="inner_custom",
            kind="custom",
            mode="circular",
            show=True,
            placement=CircularTrackPlacement(**placement_kwargs),
            params={
                "caption": "Positive CDS",
                "feature_types": ["CDS"],
                "strand_mode": "positive",
                "rules": [],
            },
        ),
        TrackSpec(
            id="features",
            kind="features",
            mode="circular",
            show=False,
            params={"feature_types": SELECTED_FEATURES},
        ),
        TrackSpec(id="gc_content", kind="gc_content", mode="circular", show=False),
        TrackSpec(id="gc_skew", kind="gc_skew", mode="circular", show=False),
    ]


def _make_center_definition_autofit_track_specs(
    *,
    include_analysis: bool = False,
    analysis_width_px: float | None = None,
) -> list[TrackSpec]:
    track_specs = [
        TrackSpec(
            id="features",
            kind="features",
            mode="circular",
            show=True,
            params={"feature_types": SELECTED_FEATURES},
        ),
        TrackSpec(id="gc_content", kind="gc_content", mode="circular", show=True),
        TrackSpec(id="gc_skew", kind="gc_skew", mode="circular", show=True),
    ]
    if include_analysis:
        analysis_placement = (
            CircularTrackPlacement(width=ScalarSpec(float(analysis_width_px), "px"))
            if analysis_width_px is not None
            else None
        )
        track_specs.append(
            TrackSpec(
                id="at_skew",
                kind="analysis",
                mode="circular",
                show=True,
                placement=analysis_placement,
                params={"caption": "AT skew", "metric": "skew", "dinucleotide": "AT"},
            )
        )
    return track_specs


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

    canvas = assemble_circular_diagram_from_record(
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
    root = ET.fromstring(canvas.tostring())

    protein_group = _find_track_group(root, "protein_track")
    feature_group = _find_track_group(root, "features")
    gc_group = _find_track_group(root, "gc_content")
    skew_group = _find_track_group(root, "gc_skew")

    assert protein_group is not None
    assert feature_group is not None
    assert gc_group is not None
    assert skew_group is not None

    protein_inner = _get_group_attr_float(protein_group, "data-track-inner-radius")
    protein_outer = _get_group_attr_float(protein_group, "data-track-outer-radius")
    feature_inner = _get_group_attr_float(feature_group, "data-track-inner-radius")
    feature_outer = _get_group_attr_float(feature_group, "data-track-outer-radius")
    gc_inner = _get_group_attr_float(gc_group, "data-track-inner-radius")
    gc_outer = _get_group_attr_float(gc_group, "data-track-outer-radius")
    skew_outer = _get_group_attr_float(skew_group, "data-track-outer-radius")

    assert protein_inner is not None
    assert protein_outer is not None
    assert feature_inner is not None
    assert feature_outer is not None
    assert gc_inner is not None
    assert gc_outer is not None
    assert skew_outer is not None
    assert protein_outer > feature_outer > gc_outer > skew_outer
    assert protein_inner > feature_outer
    assert feature_inner > gc_outer
    assert gc_inner > skew_outer

    canvas = assemble_circular_diagram_from_record(
        record,
        config_dict=config_dict,
        selected_features_set=SELECTED_FEATURES,
        legend="none",
        track_specs=["features@w=96px"],
    )
    root = ET.fromstring(canvas.tostring())
    single_feature_group = _find_track_group(root, "features")
    single_gc_group = _find_track_group(root, "gc_content")

    assert single_feature_group is not None
    assert single_gc_group is not None

    single_feature_inner = _get_group_attr_float(single_feature_group, "data-track-inner-radius")
    single_gc_outer = _get_group_attr_float(single_gc_group, "data-track-outer-radius")

    assert single_feature_inner is not None
    assert single_gc_outer is not None
    assert single_feature_inner > single_gc_outer


def test_auto_width_custom_track_shrinks_to_clear_center_definition() -> None:
    record = _load_record()
    config_dict = _make_config_dict(show_labels=False)
    cfg = GbdrawConfig.from_dict(config_dict)
    base_radius = float(cfg.canvas.circular.radius)
    default_width = (
        base_radius
        * float(cfg.canvas.circular.track_ratio)
        * float(cfg.canvas.circular.track_ratio_factors["short"][0])
    )
    track_gap_px = max(4.0, min(16.0, default_width * 0.15))

    canvas = assemble_circular_diagram_from_record(
        record,
        config_dict=config_dict,
        selected_features_set=SELECTED_FEATURES,
        legend="none",
        track_specs=_make_definition_clearance_track_specs(),
    )
    root = ET.fromstring(canvas.tostring())

    track_group = root.find(".//svg:g[@id='track_inner_custom']", _SVG_NS)
    definition_group = root.find(f".//svg:g[@id='{record.id}_definition']", _SVG_NS)

    assert track_group is not None
    assert definition_group is not None

    center_radius = _get_group_attr_float(track_group, "data-track-center-radius")
    inner_radius = _get_group_attr_float(track_group, "data-track-inner-radius")
    width_px = _get_group_attr_float(track_group, "data-track-width")
    definition_max_radius = _get_group_attr_float(definition_group, "data-definition-max-radius")

    assert center_radius is not None
    assert inner_radius is not None
    assert width_px is not None
    assert definition_max_radius is not None
    assert math.isclose(center_radius, 0.345 * base_radius, rel_tol=1e-6, abs_tol=1e-6)
    assert width_px < default_width
    assert inner_radius >= (definition_max_radius + track_gap_px - 1e-3)


def test_auto_width_builtin_tracks_shrink_to_clear_center_definition() -> None:
    record = _load_record()
    config_dict = _make_config_dict(show_labels=False)
    cfg = GbdrawConfig.from_dict(config_dict)
    feature_default_bounds = _capture_annular_annuli(
        [
            TrackSpec(
                id="features",
                kind="features",
                mode="circular",
                show=True,
                params={"feature_types": SELECTED_FEATURES},
            )
        ]
    )["features"]
    default_feature_width = float(feature_default_bounds[1]) - float(feature_default_bounds[0])
    base_radius = float(cfg.canvas.circular.radius)
    default_gc_width = (
        base_radius
        * float(cfg.canvas.circular.track_ratio)
        * float(cfg.canvas.circular.track_ratio_factors["short"][1])
    )
    default_skew_width = (
        base_radius
        * float(cfg.canvas.circular.track_ratio)
        * float(cfg.canvas.circular.track_ratio_factors["short"][2])
    )
    track_gap_px = max(4.0, min(16.0, max(default_feature_width, default_gc_width, default_skew_width) * 0.15))

    canvas = assemble_circular_diagram_from_record(
        record,
        config_dict=config_dict,
        selected_features_set=SELECTED_FEATURES,
        legend="none",
        track_specs=_make_center_definition_autofit_track_specs(),
    )
    root = ET.fromstring(canvas.tostring())

    feature_group = _find_track_group(root, "features")
    gc_group = _find_track_group(root, "gc_content")
    skew_group = _find_track_group(root, "gc_skew")
    definition_group = root.find(f".//svg:g[@id='{record.id}_definition']", _SVG_NS)

    assert feature_group is not None
    assert gc_group is not None
    assert skew_group is not None
    assert definition_group is not None

    feature_width = _get_group_attr_float(feature_group, "data-track-width")
    gc_width = _get_group_attr_float(gc_group, "data-track-width")
    skew_width = _get_group_attr_float(skew_group, "data-track-width")
    skew_inner = _get_group_attr_float(skew_group, "data-track-inner-radius")
    definition_max_radius = _get_group_attr_float(definition_group, "data-definition-max-radius")

    assert feature_width is not None
    assert gc_width is not None
    assert skew_width is not None
    assert skew_inner is not None
    assert definition_max_radius is not None
    assert feature_width < default_feature_width
    assert gc_width < default_gc_width
    assert skew_width < default_skew_width
    assert skew_inner >= (definition_max_radius + track_gap_px - 1e-3)


def test_auto_width_analysis_track_shares_width_budget_and_clears_definition() -> None:
    record = _load_record()
    config_dict = _make_config_dict(show_labels=False)
    cfg = GbdrawConfig.from_dict(config_dict)
    feature_default_bounds = _capture_annular_annuli(
        [
            TrackSpec(
                id="features",
                kind="features",
                mode="circular",
                show=True,
                params={"feature_types": SELECTED_FEATURES},
            )
        ]
    )["features"]
    default_feature_width = float(feature_default_bounds[1]) - float(feature_default_bounds[0])
    base_radius = float(cfg.canvas.circular.radius)
    default_gc_width = (
        base_radius
        * float(cfg.canvas.circular.track_ratio)
        * float(cfg.canvas.circular.track_ratio_factors["short"][1])
    )
    default_skew_width = (
        base_radius
        * float(cfg.canvas.circular.track_ratio)
        * float(cfg.canvas.circular.track_ratio_factors["short"][2])
    )
    default_analysis_width = default_skew_width
    track_gap_px = max(
        4.0,
        min(16.0, max(default_feature_width, default_gc_width, default_skew_width, default_analysis_width) * 0.15),
    )

    canvas = assemble_circular_diagram_from_record(
        record,
        config_dict=config_dict,
        selected_features_set=SELECTED_FEATURES,
        legend="none",
        track_specs=_make_center_definition_autofit_track_specs(include_analysis=True),
    )
    root = ET.fromstring(canvas.tostring())

    feature_group = _find_track_group(root, "features")
    gc_group = _find_track_group(root, "gc_content")
    skew_group = _find_track_group(root, "gc_skew")
    analysis_group = _find_track_group(root, "at_skew")
    definition_group = root.find(f".//svg:g[@id='{record.id}_definition']", _SVG_NS)

    assert feature_group is not None
    assert gc_group is not None
    assert skew_group is not None
    assert analysis_group is not None
    assert definition_group is not None

    feature_width = _get_group_attr_float(feature_group, "data-track-width")
    gc_width = _get_group_attr_float(gc_group, "data-track-width")
    skew_width = _get_group_attr_float(skew_group, "data-track-width")
    analysis_width = _get_group_attr_float(analysis_group, "data-track-width")
    analysis_inner = _get_group_attr_float(analysis_group, "data-track-inner-radius")
    definition_max_radius = _get_group_attr_float(definition_group, "data-definition-max-radius")

    assert feature_width is not None
    assert gc_width is not None
    assert skew_width is not None
    assert analysis_width is not None
    assert analysis_inner is not None
    assert definition_max_radius is not None
    assert feature_width < default_feature_width
    assert gc_width < default_gc_width
    assert skew_width < default_skew_width
    assert analysis_width < default_analysis_width
    assert analysis_inner >= (definition_max_radius + track_gap_px - 1e-3)
    assert math.isclose(analysis_width, skew_width, rel_tol=1e-6, abs_tol=1e-3)


def test_explicit_analysis_width_is_preserved_while_auto_tracks_shrink() -> None:
    record = _load_record()
    config_dict = _make_config_dict(show_labels=False)
    cfg = GbdrawConfig.from_dict(config_dict)
    feature_default_bounds = _capture_annular_annuli(
        [
            TrackSpec(
                id="features",
                kind="features",
                mode="circular",
                show=True,
                params={"feature_types": SELECTED_FEATURES},
            )
        ]
    )["features"]
    default_feature_width = float(feature_default_bounds[1]) - float(feature_default_bounds[0])
    base_radius = float(cfg.canvas.circular.radius)
    default_gc_width = (
        base_radius
        * float(cfg.canvas.circular.track_ratio)
        * float(cfg.canvas.circular.track_ratio_factors["short"][1])
    )
    default_skew_width = (
        base_radius
        * float(cfg.canvas.circular.track_ratio)
        * float(cfg.canvas.circular.track_ratio_factors["short"][2])
    )
    explicit_width = 60.0
    track_gap_px = max(
        4.0,
        min(16.0, max(default_feature_width, default_gc_width, default_skew_width, explicit_width) * 0.15),
    )

    canvas = assemble_circular_diagram_from_record(
        record,
        config_dict=config_dict,
        selected_features_set=SELECTED_FEATURES,
        legend="none",
        track_specs=_make_center_definition_autofit_track_specs(
            include_analysis=True,
            analysis_width_px=explicit_width,
        ),
    )
    root = ET.fromstring(canvas.tostring())

    feature_group = _find_track_group(root, "features")
    gc_group = _find_track_group(root, "gc_content")
    skew_group = _find_track_group(root, "gc_skew")
    analysis_group = _find_track_group(root, "at_skew")
    definition_group = root.find(f".//svg:g[@id='{record.id}_definition']", _SVG_NS)

    assert feature_group is not None
    assert gc_group is not None
    assert skew_group is not None
    assert analysis_group is not None
    assert definition_group is not None

    feature_width = _get_group_attr_float(feature_group, "data-track-width")
    gc_width = _get_group_attr_float(gc_group, "data-track-width")
    skew_width = _get_group_attr_float(skew_group, "data-track-width")
    analysis_width = _get_group_attr_float(analysis_group, "data-track-width")
    analysis_inner = _get_group_attr_float(analysis_group, "data-track-inner-radius")
    definition_max_radius = _get_group_attr_float(definition_group, "data-definition-max-radius")

    assert feature_width is not None
    assert gc_width is not None
    assert skew_width is not None
    assert analysis_width is not None
    assert analysis_inner is not None
    assert definition_max_radius is not None
    assert feature_width < default_feature_width
    assert gc_width < default_gc_width
    assert skew_width < default_skew_width
    assert math.isclose(analysis_width, explicit_width, rel_tol=1e-6, abs_tol=1e-3)
    assert analysis_inner >= (definition_max_radius + track_gap_px - 1e-3)


def test_explicit_custom_track_width_is_not_shrunk_for_definition_clearance() -> None:
    record = _load_record()
    config_dict = _make_config_dict(show_labels=False)
    cfg = GbdrawConfig.from_dict(config_dict)
    base_radius = float(cfg.canvas.circular.radius)
    explicit_width = (
        base_radius
        * float(cfg.canvas.circular.track_ratio)
        * float(cfg.canvas.circular.track_ratio_factors["short"][0])
    )
    track_gap_px = max(4.0, min(16.0, explicit_width * 0.15))

    canvas = assemble_circular_diagram_from_record(
        record,
        config_dict=config_dict,
        selected_features_set=SELECTED_FEATURES,
        legend="none",
        track_specs=_make_definition_clearance_track_specs(width_px=explicit_width),
    )
    root = ET.fromstring(canvas.tostring())

    track_group = root.find(".//svg:g[@id='track_inner_custom']", _SVG_NS)
    definition_group = root.find(f".//svg:g[@id='{record.id}_definition']", _SVG_NS)

    assert track_group is not None
    assert definition_group is not None

    width_px = _get_group_attr_float(track_group, "data-track-width")
    inner_radius = _get_group_attr_float(track_group, "data-track-inner-radius")
    definition_max_radius = _get_group_attr_float(definition_group, "data-definition-max-radius")

    assert width_px is not None
    assert inner_radius is not None
    assert definition_max_radius is not None
    assert math.isclose(width_px, explicit_width, rel_tol=1e-6, abs_tol=1e-6)
    assert inner_radius < (definition_max_radius + track_gap_px)


def test_track_and_definition_groups_emit_geometry_metadata() -> None:
    record = _load_record()
    config_dict = _make_config_dict(show_labels=False)

    canvas = assemble_circular_diagram_from_record(
        record,
        config_dict=config_dict,
        selected_features_set=SELECTED_FEATURES,
        legend="right",
        track_specs=[
            TrackSpec(id="features", kind="features", mode="circular", show=True),
            TrackSpec(
                id="at_content",
                kind="analysis",
                mode="circular",
                show=True,
                placement=CircularTrackPlacement(
                    radius=ScalarSpec(0.72, "factor"),
                    width=ScalarSpec(18.0, "px"),
                ),
                params={"caption": "AT content", "metric": "content", "dinucleotide": "AT"},
            ),
        ],
    )
    root = ET.fromstring(canvas.tostring())

    feature_group = root.find(".//svg:g[@data-track-id='features']", _SVG_NS)
    analysis_group = root.find(".//svg:g[@id='track_at_content']", _SVG_NS)
    definition_group = root.find(f".//svg:g[@id='{record.id}_definition']", _SVG_NS)

    assert feature_group is not None
    assert analysis_group is not None
    assert definition_group is not None
    assert feature_group.attrib["data-track-id"] == "features"
    assert feature_group.attrib["data-track-kind"] == "features"
    assert "data-track-center-radius" in feature_group.attrib
    assert "data-track-inner-radius" in feature_group.attrib
    assert "data-track-outer-radius" in feature_group.attrib
    assert "data-track-width" in feature_group.attrib
    assert analysis_group.attrib["data-track-id"] == "at_content"
    assert analysis_group.attrib["data-track-kind"] == "analysis"
    assert analysis_group.attrib["data-track-metric"] == "content"
    assert "data-track-center-factor" in analysis_group.attrib
    assert "data-definition-max-radius" in definition_group.attrib
