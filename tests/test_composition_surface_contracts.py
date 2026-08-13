from __future__ import annotations

import json
from pathlib import Path
import re
from types import SimpleNamespace
from typing import Literal
import xml.etree.ElementTree as ET

from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
import pytest

import gbdraw.render.export as export_module
from gbdraw.api import (
    CircularDiagramOptions,
    CircularDiagramRequest,
    CircularOutputOptions,
    InMemoryRecordSource,
    LinearDiagramOptions,
    LinearDiagramRequest,
    LinearOutputOptions,
    PreparedDiagramRequest,
    RecordInput,
    build_request_diagram,
    render_to_bytes,
    save_figure_to,
)
from gbdraw.api.diagram import build_circular_diagram, build_linear_diagram
from gbdraw.exceptions import ValidationError
from gbdraw.render.composition import (
    COMPOSITION_METADATA_ATTRIBUTE,
    COMPOSITION_ROLE_ATTRIBUTE,
    COMPOSITION_SCHEMA_ATTRIBUTE,
)
from gbdraw.session_io import CURRENT_SESSION_VERSION
from gbdraw.session_request_codec import (
    CANONICAL_REQUEST_SCHEMA,
    EncodedCanonicalRequest,
    decode_canonical_request,
    encode_canonical_request,
)


CIRCULAR_LEGEND_VALUES = (
    "left",
    "right",
    "top",
    "bottom",
    "upper_left",
    "upper_right",
    "lower_left",
    "lower_right",
    "none",
)
LINEAR_LEGEND_VALUES = ("left", "right", "top", "bottom", "none")
CIRCULAR_TITLE_VALUES = ("none", "top", "bottom")
LINEAR_TITLE_VALUES = ("center", "top", "bottom")


def _record() -> SeqRecord:
    record = SeqRecord(
        Seq("ATG" * 600),
        id="contract",
        name="contract",
        description="composition surface contract",
        annotations={"molecule_type": "DNA", "organism": "Contract example"},
    )
    record.features = [
        SeqFeature(
            FeatureLocation(60, 900, strand=1),
            type="CDS",
            qualifiers={"product": ["composition protein"]},
        ),
        SeqFeature(
            FeatureLocation(1080, 1260, strand=-1),
            type="tRNA",
            qualifiers={"product": ["tRNA-Lys"]},
        ),
    ]
    return record


def _request(
    mode: Literal["circular", "linear"],
    *,
    legend: str,
    title_position: str,
    plot_title: str | None = None,
) -> CircularDiagramRequest | LinearDiagramRequest:
    record = RecordInput(source=InMemoryRecordSource(_record()))
    if mode == "circular":
        output = CircularOutputOptions(
            legend=legend,
            plot_title_position=title_position,
        )
        return CircularDiagramRequest(
            records=(record,),
            options=CircularDiagramOptions(
                plot_title=plot_title,
                output=output,
            ),
        )
    output = LinearOutputOptions(
        legend=legend,
        plot_title_position=title_position,
    )
    return LinearDiagramRequest(
        records=(record,),
        options=LinearDiagramOptions(
            plot_title=plot_title,
            output=output,
        ),
    )


def _materialize_resources(
    encoded: EncodedCanonicalRequest,
    directory: Path,
) -> dict[str, Path]:
    directory.mkdir(parents=True)
    paths: dict[str, Path] = {}
    for resource in encoded.resources:
        path = directory / resource.name
        content = (
            resource.content
            if resource.content is not None
            else resource.source_path.read_bytes()
        )
        path.write_bytes(content)
        paths[resource.resource_id] = path
    return paths


def _round_trip_request(
    request: CircularDiagramRequest | LinearDiagramRequest,
    tmp_path: Path,
) -> tuple[EncodedCanonicalRequest, CircularDiagramRequest | LinearDiagramRequest]:
    encoded = encode_canonical_request(request)
    replay = decode_canonical_request(
        encoded.payload,
        resource_paths=_materialize_resources(encoded, tmp_path / "resources"),
        output_directory=tmp_path / "output",
    )
    assert isinstance(replay, (CircularDiagramRequest, LinearDiagramRequest))
    return encoded, replay


def _build_root(
    request: CircularDiagramRequest | LinearDiagramRequest,
) -> ET.Element:
    prepared = build_request_diagram(request)
    assert isinstance(prepared, PreparedDiagramRequest)
    return ET.fromstring(prepared.drawing.tostring())


def _composition(root: ET.Element) -> dict[str, object]:
    assert root.attrib[COMPOSITION_SCHEMA_ATTRIBUTE] == "1"
    return json.loads(root.attrib[COMPOSITION_METADATA_ATTRIBUTE])


def _normalized_target_transforms(
    root: ET.Element,
) -> tuple[tuple[str, tuple[tuple[str, tuple[float, ...]], ...]], ...]:
    targets: list[tuple[str, tuple[tuple[str, tuple[float, ...]], ...]]] = []
    for element in root.iter():
        role = element.attrib.get(COMPOSITION_ROLE_ATTRIBUTE)
        if role is None:
            continue
        operations = tuple(
            (
                name,
                tuple(
                    round(float(number), 9)
                    for number in re.findall(
                        r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?",
                        arguments,
                    )
                ),
            )
            for name, arguments in re.findall(
                r"([A-Za-z]+)\(([^)]*)\)",
                element.attrib.get("transform", ""),
            )
        )
        targets.append((role, operations))
    return tuple(sorted(targets))


@pytest.mark.parametrize("legend", CIRCULAR_LEGEND_VALUES)
def test_circular_legend_values_cross_current_schema_and_reach_composition(
    legend: str,
    tmp_path: Path,
) -> None:
    request = _request(
        "circular",
        legend=legend,
        title_position="none",
    )
    encoded, replay = _round_trip_request(request, tmp_path)

    assert isinstance(request.options.output, CircularOutputOptions)
    assert encoded.payload["schema"] == CANONICAL_REQUEST_SCHEMA
    assert encoded.payload["diagramOptions"]["output"]["legend"] == legend
    assert isinstance(replay.options.output, CircularOutputOptions)
    assert replay.options.output.legend == legend
    metadata = _composition(_build_root(replay))
    assert metadata["legendSide"] == legend
    assert (metadata["legend"] is None) is (legend == "none")


@pytest.mark.parametrize("legend", LINEAR_LEGEND_VALUES)
def test_linear_legend_values_cross_current_schema_and_reach_composition(
    legend: str,
    tmp_path: Path,
) -> None:
    request = _request(
        "linear",
        legend=legend,
        title_position="bottom",
    )
    encoded, replay = _round_trip_request(request, tmp_path)

    assert isinstance(request.options.output, LinearOutputOptions)
    assert encoded.payload["schema"] == CANONICAL_REQUEST_SCHEMA
    assert encoded.payload["diagramOptions"]["output"]["legend"] == legend
    assert isinstance(replay.options.output, LinearOutputOptions)
    assert replay.options.output.legend == legend
    metadata = _composition(_build_root(replay))
    assert metadata["legendSide"] == legend
    assert (metadata["legend"] is None) is (legend == "none")


@pytest.mark.parametrize("title_position", CIRCULAR_TITLE_VALUES)
def test_circular_title_values_cross_schema_5_and_reach_composition(
    title_position: str,
    tmp_path: Path,
) -> None:
    request = _request(
        "circular",
        legend="none",
        title_position=title_position,
        plot_title="Circular composition contract",
    )
    encoded, replay = _round_trip_request(request, tmp_path)

    assert encoded.payload["diagramOptions"]["output"]["plotTitlePosition"] == (
        title_position
    )
    assert isinstance(replay.options.output, CircularOutputOptions)
    assert replay.options.output.plot_title_position == title_position
    metadata = _composition(_build_root(replay))
    assert metadata["titleSide"] == title_position
    assert (metadata["title"] is None) is (title_position == "none")


@pytest.mark.parametrize("title_position", LINEAR_TITLE_VALUES)
def test_linear_title_values_cross_schema_5_and_reach_composition(
    title_position: str,
    tmp_path: Path,
) -> None:
    request = _request(
        "linear",
        legend="none",
        title_position=title_position,
        plot_title="Linear composition contract",
    )
    encoded, replay = _round_trip_request(request, tmp_path)

    assert encoded.payload["diagramOptions"]["output"]["plotTitlePosition"] == (
        title_position
    )
    assert isinstance(replay.options.output, LinearOutputOptions)
    assert replay.options.output.plot_title_position == title_position
    metadata = _composition(_build_root(replay))
    assert metadata["titleSide"] == title_position
    assert metadata["title"] is not None


@pytest.mark.parametrize(
    ("mode", "expected_title_side", "has_title"),
    (
        ("circular", "none", False),
        ("linear", "bottom", True),
    ),
)
def test_output_option_defaults_remain_mode_specific_after_current_schema_replay(
    mode: Literal["circular", "linear"],
    expected_title_side: str,
    has_title: bool,
    tmp_path: Path,
) -> None:
    record = RecordInput(source=InMemoryRecordSource(_record()))
    request = (
        CircularDiagramRequest(
            records=(record,),
            options=CircularDiagramOptions(
                plot_title="Default placement",
                output=CircularOutputOptions(),
            ),
        )
        if mode == "circular"
        else LinearDiagramRequest(
            records=(record,),
            options=LinearDiagramOptions(
                plot_title="Default placement",
                output=LinearOutputOptions(),
            ),
        )
    )
    encoded, replay = _round_trip_request(request, tmp_path)

    assert encoded.payload["diagramOptions"]["output"] == {
        "legend": "right",
        "plotTitlePosition": None,
    }
    metadata = _composition(_build_root(replay))
    assert metadata["legendSide"] == "right"
    assert metadata["titleSide"] == expected_title_side
    assert (metadata["title"] is not None) is has_title


@pytest.mark.parametrize(
    ("mode", "legend", "title_position"),
    (
        ("circular", "lower_right", "bottom"),
        ("linear", "top", "center"),
    ),
)
def test_fresh_and_current_schema_replay_have_identical_automatic_composition(
    mode: Literal["circular", "linear"],
    legend: str,
    title_position: str,
    tmp_path: Path,
) -> None:
    request = _request(
        mode,
        legend=legend,
        title_position=title_position,
        plot_title="Replay equivalence",
    )
    _, replay = _round_trip_request(request, tmp_path)

    fresh_root = _build_root(request)
    replay_root = _build_root(replay)

    assert replay_root.attrib["viewBox"] == fresh_root.attrib["viewBox"]
    assert replay_root.attrib[COMPOSITION_METADATA_ATTRIBUTE] == fresh_root.attrib[
        COMPOSITION_METADATA_ATTRIBUTE
    ]
    assert _normalized_target_transforms(
        replay_root
    ) == _normalized_target_transforms(fresh_root)


@pytest.mark.parametrize(
    ("options_type", "invalid_legend"),
    (
        (CircularOutputOptions, "center"),
        (CircularOutputOptions, " upper_left"),
        (CircularOutputOptions, "UPPER_LEFT"),
        (CircularOutputOptions, None),
        (LinearOutputOptions, "upper_left"),
        (LinearOutputOptions, "lower_right"),
        (LinearOutputOptions, "center"),
        (LinearOutputOptions, ""),
    ),
)
def test_output_options_reject_legend_values_outside_the_mode_contract(
    options_type: type[CircularOutputOptions] | type[LinearOutputOptions],
    invalid_legend: object,
) -> None:
    with pytest.raises(ValidationError, match="legend must be one of"):
        options_type(legend=invalid_legend)  # type: ignore[arg-type]


def test_linear_cli_rejects_circular_corner_legend_before_rendering(
    gbdraw_runner,
    test_inputs_dir: Path,
    tmp_path: Path,
) -> None:
    return_code, output, svg_path = gbdraw_runner.run(
        "linear",
        [test_inputs_dir / "BGC0000708.gbk"],
        "invalid-corner-legend",
        tmp_path,
        extra_args=["--legend", "upper_left"],
    )

    assert return_code == 1
    assert "Linear legend must be one of" in output
    assert not svg_path.exists()


def test_current_request_and_session_schema_versions_are_explicit() -> None:
    assert CANONICAL_REQUEST_SCHEMA == 6
    assert CURRENT_SESSION_VERSION == 41


@pytest.mark.parametrize("mode", ("circular", "linear"))
def test_real_exports_use_the_composed_svg_viewbox(
    mode: Literal["circular", "linear"],
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    if mode == "circular":
        drawing = build_circular_diagram(
            _record(),
            options=CircularDiagramOptions(
                plot_title="Export contract",
                output=CircularOutputOptions(
                    legend="left",
                    plot_title_position="top",
                ),
            ),
        )
    else:
        drawing = build_linear_diagram(
            [_record()],
            options=LinearDiagramOptions(
                plot_title="Export contract",
                output=LinearOutputOptions(
                    legend="top",
                    plot_title_position="bottom",
                ),
            ),
        )

    original_root = ET.fromstring(drawing.tostring())
    expected_viewbox = original_root.attrib["viewBox"]
    real_cairosvg = export_module.get_cairosvg()
    converted_sources: list[bytes] = []

    def capture(method_name: str, *, bytestring: bytes, write_to) -> object:
        converted_sources.append(bytestring)
        return getattr(real_cairosvg, method_name)(
            bytestring=bytestring,
            write_to=write_to,
        )

    proxy = SimpleNamespace(
        svg2png=lambda **kwargs: capture("svg2png", **kwargs),
        svg2pdf=lambda **kwargs: capture("svg2pdf", **kwargs),
        svg2ps=lambda **kwargs: capture("svg2ps", **kwargs),
    )
    monkeypatch.setattr(export_module, "get_cairosvg", lambda: proxy)

    paths = save_figure_to(
        drawing,
        ("svg", "png", "pdf", "eps", "ps"),
        output_dir=str(tmp_path),
        output_prefix=mode,
    )
    by_suffix = {Path(path).suffix: Path(path) for path in paths}

    assert set(by_suffix) == {".svg", ".png", ".pdf", ".eps", ".ps"}
    assert all(path.is_file() and path.stat().st_size > 0 for path in by_suffix.values())
    assert by_suffix[".svg"].read_bytes().lstrip().startswith(b"<?xml")
    assert by_suffix[".png"].read_bytes().startswith(b"\x89PNG\r\n\x1a\n")
    assert by_suffix[".pdf"].read_bytes().startswith(b"%PDF-")
    assert by_suffix[".eps"].read_bytes().startswith(b"%!PS-Adobe")
    assert by_suffix[".ps"].read_bytes().startswith(b"%!PS-Adobe")

    saved_root = ET.parse(by_suffix[".svg"]).getroot()
    assert saved_root.attrib["viewBox"] == expected_viewbox
    assert saved_root.attrib[COMPOSITION_SCHEMA_ATTRIBUTE] == "1"
    assert len(converted_sources) == 4
    for source in converted_sources:
        converted_root = ET.fromstring(source)
        assert converted_root.attrib["viewBox"] == expected_viewbox
        assert converted_root.attrib[COMPOSITION_SCHEMA_ATTRIBUTE] == "1"
        assert (
            converted_root.attrib[COMPOSITION_METADATA_ATTRIBUTE]
            == original_root.attrib[COMPOSITION_METADATA_ATTRIBUTE]
        )


@pytest.mark.parametrize("mode", ("circular", "linear"))
def test_interactive_svg_preserves_original_viewbox_and_composition_schema(
    mode: Literal["circular", "linear"],
) -> None:
    drawing = (
        build_circular_diagram(
            _record(),
            options=CircularDiagramOptions(
                output=CircularOutputOptions(legend="bottom")
            ),
        )
        if mode == "circular"
        else build_linear_diagram(
            [_record()],
            options=LinearDiagramOptions(
                output=LinearOutputOptions(legend="bottom")
            ),
        )
    )
    source_root = ET.fromstring(drawing.tostring())
    interactive_root = ET.fromstring(render_to_bytes(drawing, "interactive_svg"))

    assert interactive_root.attrib["viewBox"] == source_root.attrib["viewBox"]
    assert interactive_root.attrib[COMPOSITION_SCHEMA_ATTRIBUTE] == "1"
    assert interactive_root.attrib[COMPOSITION_METADATA_ATTRIBUTE] == source_root.attrib[
        COMPOSITION_METADATA_ATTRIBUTE
    ]
