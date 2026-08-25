from __future__ import annotations

import base64
from collections.abc import Callable
from functools import cache
import math
from pathlib import Path
import re
from tempfile import TemporaryDirectory
from typing import Any
from xml.etree import ElementTree as ET

import pytest

from gbdraw.api import (
    load_session_document,
    materialize_session,
    session_to_request,
)
from gbdraw.api.session_compat import build_session_compatible_request_diagram
from gbdraw.session_request_codec import CANONICAL_REQUEST_SCHEMA
from tools.prepare_interactive_gallery_assets import EXAMPLES, GallerySessionExample


pytestmark = pytest.mark.gallery


_TAG_RE = re.compile(r"<[^>]+>")
_SVG_NUMBER = r"-?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"
_PATH_POINT_RE = re.compile(
    rf"[ML]\s*({_SVG_NUMBER})\s*,?\s*({_SVG_NUMBER})"
)
_TRANSLATE_COMPONENT_RE = re.compile(
    rf"translate\(\s*({_SVG_NUMBER})\s*[, ]\s*({_SVG_NUMBER})\s*\)"
)


def _translate_chain_y(transform: str) -> float:
    matches = list(_TRANSLATE_COMPONENT_RE.finditer(transform))
    assert matches, transform
    cursor = 0
    for match in matches:
        assert not transform[cursor : match.start()].strip(), transform
        cursor = match.end()
    assert not transform[cursor:].strip(), transform
    return sum(float(match.group(2)) for match in matches)


@cache
def _example(example_id: str) -> GallerySessionExample:
    return next(example for example in EXAMPLES if example.id == example_id)


@pytest.fixture(scope="module")
def gallery_sessions(
    load_cached_gallery_session: Callable[[Path], dict[str, object]],
) -> dict[str, tuple[GallerySessionExample, dict[str, object]]]:
    return {
        example.id: (example, load_cached_gallery_session(example.session_path))
        for example in EXAMPLES
    }


def _request(session: dict[str, object]) -> dict[str, object]:
    request = session["renderRequest"]
    assert isinstance(request, dict)
    assert request["schema"] in {5, CANONICAL_REQUEST_SCHEMA}
    return request


@pytest.mark.parametrize(
    "example",
    EXAMPLES,
    ids=lambda example: example.id,
)
def test_gallery_sessions_do_not_persist_overwrite_permission(
    example: GallerySessionExample,
    gallery_sessions: dict[str, tuple[GallerySessionExample, dict[str, object]]],
) -> None:
    _, session = gallery_sessions[example.id]
    output = _request(session)["output"]
    outputs = output if isinstance(output, list) else [output]

    assert outputs
    assert all(
        isinstance(item, dict) and item.get("overwrite") is False
        for item in outputs
    )


def test_lambda_basic_linear_gallery_session_opts_out_of_comparisons(
    gallery_sessions: dict[str, tuple[GallerySessionExample, dict[str, object]]],
) -> None:
    _, session = gallery_sessions["lambda_basic_linear"]

    plan = session["config"]["linearComparisonPlan"]
    assert plan["mode"] == "none"
    assert plan["edges"] == []
    assert _request(session)["comparisons"] == []


def _resource_bytes(
    session: dict[str, object], ref: dict[str, object]
) -> bytes:
    resource_id = ref["resourceId"]
    assert isinstance(resource_id, str)
    resources = session["resources"]
    assert isinstance(resources, dict)
    assert resource_id in resources
    resource = resources[resource_id]
    assert isinstance(resource, dict)
    data = resource.get("data", "")
    assert isinstance(data, str) and data
    if resource.get("encoding") == "base64":
        return base64.b64decode(data)
    return data.encode("utf-8")


def _resource_text(
    session: dict[str, object], ref: dict[str, object]
) -> str:
    return _resource_bytes(session, ref).decode("utf-8")


def _option_resource_ref(
    options: dict[str, object], *names: str
) -> dict[str, object]:
    for name in names:
        value = options.get(name)
        if isinstance(value, dict) and value.get("resourceId"):
            return value
    raise AssertionError(f"missing resource-backed option: {' or '.join(names)}")


@pytest.fixture(scope="module")
def gallery_visual_roots(
    gallery_sessions: dict[str, tuple[GallerySessionExample, dict[str, object]]],
    load_cached_svg_root: Callable[[Path], ET.Element],
) -> Callable[[str, bool], tuple[tuple[str, ET.Element], ...]]:
    @cache
    def load(
        example_id: str,
        include_gallery: bool = False,
    ) -> tuple[tuple[str, ET.Element], ...]:
        example, session = gallery_sessions[example_id]
        results = session["results"]
        assert isinstance(results, list) and results
        result = results[0]
        assert isinstance(result, dict)
        content = result.get("content")
        assert isinstance(content, str) and content
        roots = (
            ("session result", ET.fromstring(content)),
            ("gallery source", load_cached_svg_root(example.source_svg_path)),
        )
        if include_gallery:
            roots += (
                ("gallery example", load_cached_svg_root(example.gallery_svg_path)),
            )
        return roots

    return load


def _texts(root: ET.Element) -> list[str]:
    return [
        text
        for element in root.iter()
        if element.tag.endswith("text")
        if (text := "".join(element.itertext()).strip())
    ]


def _plain_label(value: object) -> str:
    return _TAG_RE.sub("", str(value or "")).strip()


def _group_translate_y(root: ET.Element, group_id: str) -> float:
    group = next(
        (element for element in root.iter() if element.get("id") == group_id),
        None,
    )
    assert group is not None, group_id
    transform = str(group.get("transform") or "")
    return _translate_chain_y(transform)


def _record_group_translate_y(
    root: ET.Element,
    *,
    record_id: str,
    record_index: int,
) -> float:
    group = next(
        (
            element
            for element in root.iter()
            if element.tag.endswith("g")
            and str(element.get("id") or "").startswith("record_group_")
            and element.get("data-gbdraw-record-id") == record_id
            and element.get("data-gbdraw-record-index") == str(record_index)
        ),
        None,
    )
    if group is None:
        legacy_group_id = f"{record_id}_record_{record_index + 1}"
        group = next(
            (
                element
                for element in root.iter()
                if element.get("id") == legacy_group_id
            ),
            None,
        )
    assert group is not None, (record_id, record_index)
    transform = str(group.get("transform") or "")
    return _translate_chain_y(transform)


def _circular_axis_and_feature_band(
    root: ET.Element,
) -> tuple[float, float, float]:
    axis = next(
        (element for element in root.iter() if element.get("id") == "Axis"),
        None,
    )
    assert axis is not None
    axis_circle = next(
        (element for element in axis.iter() if element.tag.endswith("circle")),
        None,
    )
    assert axis_circle is not None

    radii = [
        math.hypot(float(x), float(y))
        for element in root.iter()
        if element.get("data-gbdraw-feature-id")
        for path_data in [element.get("d")]
        if path_data
        for x, y in _PATH_POINT_RE.findall(path_data)
    ]
    assert radii
    return float(axis_circle.get("r", "nan")), min(radii), max(radii)


@cache
def _materialized_track_geometry(example_id: str) -> dict[str, Any]:
    example = _example(example_id)
    document = load_session_document(example.session_path)
    with TemporaryDirectory(prefix=f"gbdraw-gallery-{example_id}-") as output_dir:
        with materialize_session(
            document,
            output_directory=output_dir,
        ) as materialized:
            prepared = build_session_compatible_request_diagram(
                session_to_request(materialized),
                document.to_dict(),
            )
    geometry = getattr(prepared.drawing, "_gbdraw_track_slot_geometry", None)
    assert isinstance(geometry, dict)
    assert geometry.get("mode") == "linear"
    return geometry


def _local_track_geometry_signature(
    geometry: dict[str, Any],
) -> tuple[tuple[object, ...], ...]:
    records = geometry.get("records")
    assert isinstance(records, list)

    def band_signature(value: object) -> tuple[float, float] | None:
        if value is None:
            return None
        assert isinstance(value, dict)
        return float(value["topPx"]), float(value["bottomPx"])

    return tuple(
        (
            int(record["recordIndex"]),
            str(record["recordId"]),
            tuple(
                (
                    str(slot["slotId"]),
                    str(slot["renderer"]),
                    str(slot["side"]),
                    float(slot["resolvedOriginPx"]),
                    float(slot["heightPx"]),
                    float(slot["spacingAfterPx"]),
                    bool(slot["dataAvailable"]),
                    band_signature(slot["paintBand"]),
                    band_signature(slot["reserveBand"]),
                )
                for slot in record["slots"]
            ),
        )
        for record in records
    )


def test_hmmt_at_skew_session_keeps_tracks_palette_and_gene_labels(
    gallery_sessions: dict[str, tuple[GallerySessionExample, dict[str, object]]],
    gallery_visual_roots: Callable[[str, bool], tuple[tuple[str, ET.Element], ...]],
) -> None:
    _example, session = gallery_sessions["HmmtDNA_ATskew"]
    request = _request(session)
    options = request["diagramOptions"]
    assert isinstance(options, dict)

    colors = options["colors"]
    assert isinstance(colors, dict)
    assert colors["defaultColorsPalette"] == "ajisai"

    tracks = options["tracks"]
    assert isinstance(tracks, dict)
    slots = tracks["circularTrackSlots"]
    assert isinstance(slots, list)
    assert [slot["id"] for slot in slots] == [
        "features",
        "gc_content",
        "gc_skew",
        "a_skew_2",
        "ticks",
    ]
    at_skew = next(slot for slot in slots if slot["id"] == "a_skew_2")
    assert at_skew["renderer"] == "dinucleotide_skew"
    assert at_skew["params"] == {
        "nt": "AT",
        "positive_color": "#deaf6e",
        "negative_color": "#7294e3",
        "legend_label": "AT skew",
    }

    priority_ref = _option_resource_ref(
        options, "qualifierPriorityFile", "qualifierPriorityTable"
    )
    assert "CDS\tgene" in _resource_text(session, priority_ref)

    for location, root in gallery_visual_roots("HmmtDNA_ATskew", True):
        slot_ids = {
            element.get("data-gbdraw-slot-id")
            for element in root.iter()
            if element.tag.endswith("g")
        }
        texts = set(_texts(root))
        fills = {element.get("fill") for element in root.iter()}
        assert "a_skew_2" in slot_ids, location
        assert {"AT skew (+)", "AT skew (-)", "ND1"} <= texts, location
        assert {
            "#84b9ec",
            "#7cecd5",
            "#ddce76",
            "#deaf6e",
            "#7294e3",
        } <= fills, location


def test_hmmt_at_skew_session_and_visuals_keep_middle_feature_placement(
    gallery_sessions: dict[str, tuple[GallerySessionExample, dict[str, object]]],
    gallery_visual_roots: Callable[[str, bool], tuple[tuple[str, ET.Element], ...]],
) -> None:
    _example, session = gallery_sessions["HmmtDNA_ATskew"]
    request = _request(session)
    options = request["diagramOptions"]
    assert isinstance(options, dict)
    tracks = options["tracks"]
    assert isinstance(tracks, dict)
    slots = tracks["circularTrackSlots"]
    assert isinstance(slots, list)

    feature_slot = next(slot for slot in slots if slot["id"] == "features")
    assert feature_slot["renderer"] == "features"
    assert feature_slot["params"]["lane_direction"] == "split"

    for location, root in gallery_visual_roots("HmmtDNA_ATskew", True):
        axis_radius, feature_inner, feature_outer = (
            _circular_axis_and_feature_band(root)
        )
        assert feature_inner < axis_radius < feature_outer, location


def test_bgc_gallery_session_keeps_curated_presentation_and_styles(
    gallery_sessions: dict[str, tuple[GallerySessionExample, dict[str, object]]],
    gallery_visual_roots: Callable[[str, bool], tuple[tuple[str, ET.Element], ...]],
) -> None:
    expected_labels = (
        "Streptomyces lividus CBS 844.73",
        "Streptomyces fradiae ATCC 10745",
        "Streptomyces fradiae MCIMB 8233",
        "Streptomyces rimosus subsp. paromomycinus NRRL 2455",
        "Streptomyces ribosidificus ATCC 21294",
    )
    expected_subtitles = (
        "Lividomycin biosynthetic gene cluster",
        "Neomycin biosynthetic gene cluster",
        "Neomycin biosynthetic gene cluster",
        "Paromomycin biosynthetic gene cluster",
        "Ribostamycin biosynthetic gene",
    )
    title = "Aminoglycoside biosynthetic gene clusters from Streptomyces spp."
    rule_colors = {"#d03535", "#f787a9", "#577edb", "#57b767"}
    rule_captions = {
        "Core biosynthetic genes",
        "Additional biosynthetic genes",
        "Transport-related genes",
        "Regulatory genes",
    }

    _example, session = gallery_sessions["BGC0000708-BGC0000713"]
    request = _request(session)
    records = request["records"]
    assert isinstance(records, list)
    presentations = [record["presentation"] for record in records]
    assert tuple(_plain_label(item["label"]) for item in presentations) == expected_labels
    assert tuple(item["subtitle"] for item in presentations) == expected_subtitles

    options = request["diagramOptions"]
    assert isinstance(options, dict)
    assert _plain_label(options["plotTitle"]) == title
    colors = options["colors"]
    assert isinstance(colors, dict)
    assert colors["defaultColorsPalette"] == "orange"
    rules_ref = _option_resource_ref(colors, "colorTableFile", "colorTable")
    rules = _resource_text(session, rules_ref)
    assert rule_colors <= set(re.findall(r"#[0-9a-fA-F]{6}", rules))
    assert all(caption in rules for caption in rule_captions)
    priority_ref = _option_resource_ref(
        options, "qualifierPriorityFile", "qualifierPriorityTable"
    )
    assert "CDS\tgene" in _resource_text(session, priority_ref)

    for location, root in gallery_visual_roots("BGC0000708-BGC0000713", False):
        texts = set(_texts(root))
        fills = {element.get("fill") for element in root.iter()}
        group_ids = {element.get("id") for element in root.iter()}
        definitions = [
            element
            for element in root.iter()
            if "_definition" in (element.get("id") or "")
        ]

        assert "plot_title" in group_ids, location
        assert title in texts, location
        assert set(expected_labels) <= texts, location
        assert set(expected_subtitles) <= texts, location
        assert rule_captions <= texts, location
        assert {"#d3d3d3", *rule_colors} <= fills, location
        assert "livZ" in texts, location
        assert len(definitions) == 5, location
        for definition in definitions:
            lines = [
                element
                for element in definition.iter()
                if element.tag.endswith("text")
            ]
            assert len(lines) == 4, (location, definition.get("id"))
            assert lines[0].get("font-size") == "20.0"
            assert lines[0].get("font-weight") == "bold"
            assert lines[1].get("font-size") == "20.0"
            assert all(line.get("text-anchor") == "start" for line in lines)
            assert all(line.get("font-size") == "20.0" for line in lines[2:])
            assert all(line.get("fill") == "#7b7c7d" for line in lines[2:])


def test_majanivirus_gallery_session_keeps_record_labels_and_color_rules(
    gallery_sessions: dict[str, tuple[GallerySessionExample, dict[str, object]]],
    gallery_visual_roots: Callable[[str, bool], tuple[tuple[str, ET.Element], ...]],
) -> None:
    expected_labels = (
        "Marsupenaeus japonicus endogenous nimavirus",
        "Melicertus latisulcatus majanivirus",
        "Penaeus monodon majanivirus A",
        "Penaeus semisulcatus majanivirus",
        "Penaeus monodon majanivirus B",
        "Litopenaeus vannamei majanivirus",
        "Trachysalambria curvirostris majanivirus",
        "Metapenaeus ensis majanivirus",
        "Metapenaeus joyneri majanivirus",
    )
    captions = {"WSSV-like proteins", "BIRP", "tyrosine recombinase"}
    expected_colors = {"#89d1fa", "yellow", "red"}

    _example, session = gallery_sessions["majanivirus_orthogroup"]
    request = _request(session)
    records = request["records"]
    assert isinstance(records, list)
    assert tuple(
        _plain_label(record["presentation"]["label"]) for record in records
    ) == expected_labels

    options = request["diagramOptions"]
    assert isinstance(options, dict)
    colors = options["colors"]
    assert isinstance(colors, dict)
    rules_ref = _option_resource_ref(colors, "colorTableFile", "colorTable")
    rules = _resource_text(session, rules_ref)
    assert all(caption in rules for caption in captions)
    assert all(color in rules for color in expected_colors)

    for location, root in gallery_visual_roots("majanivirus_orthogroup", False):
        texts = set(_texts(root))
        fills = {element.get("fill") for element in root.iter()}
        definitions = [
            element
            for element in root.iter()
            if "_definition" in (element.get("id") or "")
        ]
        assert set(expected_labels) <= texts, location
        assert captions <= texts, location
        assert expected_colors <= fills, location
        assert len(definitions) == 9, location


def test_wssv_gallery_session_keeps_all_twenty_conservation_rings(
    gallery_sessions: dict[str, tuple[GallerySessionExample, dict[str, object]]],
    gallery_visual_roots: Callable[[str, bool], tuple[tuple[str, ET.Element], ...]],
) -> None:
    labels = (
        "CN01",
        "WSSV-TW",
        "WSSV-CN",
        "WSSV-TH",
        "JP01A",
        "JP01B",
        "Pc2020",
        "E1",
        "0722-1",
        "CN03",
        "CN04",
        "WSSV-AU",
        "EU129",
        "GCF7",
        "MES-753",
        "Shantou2019",
        "POMZ1",
        "POMZ4",
        "MG18PR-0187-N40S",
        "Angostura2013",
    )
    colors = (
        "#6e91b7",
        "#f4a251",
        "#77b26f",
        "#e67577",
        "#8fc4c0",
        "#f0d369",
        "#be92b2",
        "#ffafb7",
        "#ae8e7c",
        "#c6bebb",
        "#6e91b7",
        "#f4a251",
        "#e67577",
        "#8fc4c0",
        "#bcb4ca",
        "#f0d369",
        "#be92b2",
        "#ffafb7",
        "#ae8e7c",
        "#c6bebb",
    )

    _example, session = gallery_sessions["WSSV_genome_comparison"]
    request = _request(session)
    options = request["diagramOptions"]
    assert isinstance(options, dict)
    conservation = options["conservationBlastFiles"]
    assert isinstance(conservation, list) and len(conservation) == 20
    assert all(_resource_bytes(session, ref) for ref in conservation)
    assert tuple(options["conservationLabels"]) == labels
    assert tuple(options["conservationColors"]) == colors
    # The saved gallery state uses automatic BLAST-side detection.  The
    # canonical encoder may omit this default rather than serializing "auto".
    assert options.get("conservationReference") in {None, "auto"}
    assert options["conservationRingWidth"] == 5
    assert options["conservationRingGap"] == 2
    palette = options["colors"]
    assert isinstance(palette, dict)
    assert palette["defaultColorsPalette"] == "royal_gala"

    for location, root in gallery_visual_roots("WSSV_genome_comparison", False):
        group_labels = {
            element.get("data-track-label")
            for element in root.iter()
            if element.tag.endswith("g")
            and element.get("data-track-label")
        }
        texts = set(_texts(root))
        fills = {element.get("fill") for element in root.iter()}
        assert group_labels == set(labels), location
        assert set(labels) <= texts, location
        assert "#dc7078" in fills, location


@pytest.mark.parametrize(
    ("example_id", "record_count", "precomputed_count", "orthogroup_count"),
    [
        ("hepatoplasmataceae_collinear", 5, 4, 554),
        ("hepatoplasmataceae_orthogroup", 5, 4, 577),
        ("majanivirus_orthogroup", 9, 8, 152),
    ],
)
def test_gallery_sessions_keep_precomputed_comparisons_and_orthogroups(
    example_id: str,
    record_count: int,
    precomputed_count: int,
    orthogroup_count: int,
    gallery_sessions: dict[str, tuple[GallerySessionExample, dict[str, object]]],
) -> None:
    _, session = gallery_sessions[example_id]
    request = _request(session)
    records = request["records"]
    comparisons = request["comparisons"]
    assert isinstance(records, list) and len(records) == record_count
    assert isinstance(comparisons, list)

    precomputed = [
        item
        for item in comparisons
        if item["kind"] == "precomputedProteinComparison"
    ]
    orthogroups = [
        item for item in comparisons if item["kind"] == "orthogroupResult"
    ]
    assert len(precomputed) == precomputed_count
    assert len(orthogroups) == 1
    assert all(_resource_bytes(session, item) for item in precomputed)
    assert _resource_bytes(session, orthogroups[0])

    state = session["orthogroupState"]
    assert isinstance(state, dict)
    assert "groups" not in state

    editor_state = session["editorState"]
    assert isinstance(editor_state, dict)
    catalog = editor_state["featureCatalog"]
    assert isinstance(catalog, dict)
    items = catalog["items"]
    assert isinstance(items, list) and len(items) == 1
    groups = items[0]["orthogroups"]
    assert isinstance(groups, list)
    assert len(groups) == orthogroup_count
    assert all(group.get("members") for group in groups)


@pytest.mark.parametrize(
    "example_id",
    (
        "hepatoplasmataceae_collinear",
        "hepatoplasmataceae_orthogroup",
    ),
)
def test_hepatoplasmataceae_gallery_keeps_shared_track_spacing(
    example_id: str,
    gallery_sessions: dict[str, tuple[GallerySessionExample, dict[str, object]]],
    gallery_visual_roots: Callable[[str, bool], tuple[tuple[str, ET.Element], ...]],
) -> None:
    _example, session = gallery_sessions[example_id]
    request = _request(session)
    options = request["diagramOptions"]
    assert isinstance(options, dict)
    config_overrides = options["configOverrides"]
    assert isinstance(config_overrides, dict)
    declared_spacing = float(
        config_overrides.get("canvas.linear.track_spacing", 0.0)
    )
    assert declared_spacing == pytest.approx(0.0)

    geometry = _materialized_track_geometry(example_id)
    records = geometry["records"]
    assert isinstance(records, list)
    assert len(records) == 5

    for record in records:
        assert isinstance(record, dict)
        slots = {
            str(slot["slotId"]): slot
            for slot in record["slots"]
        }
        features = slots["features"]
        gc_content = slots["gc_content"]
        gc_skew = slots["gc_skew"]
        feature_reserve = features["reserveBand"]
        gc_reserve = gc_content["reserveBand"]
        skew_reserve = gc_skew["reserveBand"]
        comparison_exclusion = record["comparisonExclusionBand"]
        assert isinstance(feature_reserve, dict)
        assert isinstance(gc_reserve, dict)
        assert isinstance(skew_reserve, dict)
        assert isinstance(comparison_exclusion, dict)
        assert float(features["spacingAfterPx"]) == pytest.approx(declared_spacing)
        assert float(gc_content["spacingAfterPx"]) == pytest.approx(declared_spacing)
        assert (
            float(gc_reserve["topPx"])
            - float(feature_reserve["bottomPx"])
        ) == pytest.approx(declared_spacing)
        assert (
            float(skew_reserve["topPx"])
            - float(gc_reserve["bottomPx"])
        ) == pytest.approx(declared_spacing)
        assert (
            float(comparison_exclusion["bottomPx"])
            - float(skew_reserve["bottomPx"])
        ) == pytest.approx(declared_spacing)

    for location, root in gallery_visual_roots(example_id, True):
        for record in records:
            record_index = int(record["recordIndex"])
            record_number = record_index + 1
            record_id = str(record["recordId"])
            slots = {
                str(slot["slotId"]): slot
                for slot in record["slots"]
            }
            assert _record_group_translate_y(
                root,
                record_id=record_id,
                record_index=record_index,
            ) == pytest.approx(float(record["axisYpx"])), (
                location,
                record_id,
            )
            expected_group_y = {
                f"gc_content_record_{record_number}": float(
                    slots["gc_content"]["finalYOffsetPx"]
                ),
                f"gc_skew_record_{record_number}": float(
                    slots["gc_skew"]["finalYOffsetPx"]
                ),
            }
            for group_id, expected_y in expected_group_y.items():
                assert _group_translate_y(root, group_id) == pytest.approx(
                    expected_y
                ), (location, group_id)


def test_hepatoplasmataceae_comparison_modes_share_record_local_track_geometry() -> None:
    collinear = _materialized_track_geometry(
        "hepatoplasmataceae_collinear"
    )
    orthogroup = _materialized_track_geometry(
        "hepatoplasmataceae_orthogroup"
    )

    assert _local_track_geometry_signature(
        collinear
    ) == _local_track_geometry_signature(orthogroup)
