from __future__ import annotations

import base64
import re
from xml.etree import ElementTree as ET

import pytest

from gbdraw.session_io import load_session
from tools.prepare_interactive_gallery_assets import EXAMPLES, GallerySessionExample


_TAG_RE = re.compile(r"<[^>]+>")


def _example(example_id: str) -> GallerySessionExample:
    return next(example for example in EXAMPLES if example.id == example_id)


def _session(example_id: str) -> tuple[GallerySessionExample, dict[str, object]]:
    example = _example(example_id)
    return example, load_session(example.session_path)


def _request(session: dict[str, object]) -> dict[str, object]:
    request = session["renderRequest"]
    assert isinstance(request, dict)
    assert request["schema"] == 3
    return request


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


def _visual_roots(
    example: GallerySessionExample, session: dict[str, object]
) -> tuple[tuple[str, ET.Element], ...]:
    results = session["results"]
    assert isinstance(results, list) and results
    result = results[0]
    assert isinstance(result, dict)
    content = result.get("content")
    assert isinstance(content, str) and content
    return (
        ("session result", ET.fromstring(content)),
        ("gallery source", ET.parse(example.source_svg_path).getroot()),
    )


def _texts(root: ET.Element) -> list[str]:
    return [
        text
        for element in root.iter()
        if element.tag.endswith("text")
        if (text := "".join(element.itertext()).strip())
    ]


def _plain_label(value: object) -> str:
    return _TAG_RE.sub("", str(value or "")).strip()


def test_hmmt_at_skew_session_keeps_tracks_palette_and_gene_labels() -> None:
    example, session = _session("HmmtDNA_ATskew")
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
    assert [str(slot).split(":", 1)[0] for slot in slots] == [
        "features",
        "gc_content",
        "gc_skew",
        "a_skew_2",
        "ticks",
    ]
    at_skew = next(str(slot) for slot in slots if str(slot).startswith("a_skew_2:"))
    assert "a_skew_2:dinucleotide_skew" in at_skew
    assert "nt=AT" in at_skew
    assert "legend_label=AT skew" in at_skew
    assert "positive_color=#deaf6e" in at_skew
    assert "negative_color=#7294e3" in at_skew

    priority_ref = _option_resource_ref(
        options, "qualifierPriorityFile", "qualifierPriorityTable"
    )
    assert "CDS\tgene" in _resource_text(session, priority_ref)

    for location, root in _visual_roots(example, session):
        group_ids = {element.get("id") for element in root.iter()}
        texts = set(_texts(root))
        fills = {element.get("fill") for element in root.iter()}
        assert "a_skew_2" in group_ids, location
        assert {"AT skew (+)", "AT skew (-)", "ND1"} <= texts, location
        assert {
            "#84b9ec",
            "#7cecd5",
            "#ddce76",
            "#deaf6e",
            "#7294e3",
        } <= fills, location


def test_bgc_gallery_session_keeps_curated_presentation_and_styles() -> None:
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

    example, session = _session("BGC0000708-BGC0000713")
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

    for location, root in _visual_roots(example, session):
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
        assert {"#dddddd", *rule_colors} <= fills, location
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


def test_majanivirus_gallery_session_keeps_record_labels_and_color_rules() -> None:
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

    example, session = _session("majanivirus_orthogroup")
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

    for location, root in _visual_roots(example, session):
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


def test_wssv_gallery_session_keeps_all_twenty_conservation_rings() -> None:
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

    example, session = _session("WSSV_genome_comparison")
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

    for location, root in _visual_roots(example, session):
        group_ids = {
            element.get("id")
            for element in root.iter()
            if (element.get("id") or "").startswith("conservation_")
            and element.get("id") != "conservation_identity_legend"
        }
        texts = set(_texts(root))
        fills = {element.get("fill") for element in root.iter()}
        assert group_ids == {f"conservation_{label}" for label in labels}, location
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
) -> None:
    _, session = _session(example_id)
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
    groups = state["groups"]
    assert isinstance(groups, list)
    assert len(groups) == orthogroup_count
    assert all(group.get("members") for group in groups)
