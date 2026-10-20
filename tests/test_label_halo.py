"""How a feature label with a leader is drawn (LP-D98, LP-D100, LP-D103, LP-D106): its halo, its box, its
Linear underline (trimmed at its ends, on a white casing) and gap-stop, a Circular underline's joint, and the PDF
text of a figure with halos."""

from __future__ import annotations

import math
import re
import xml.etree.ElementTree as ET
import zlib
from collections.abc import Iterator
from pathlib import Path

import pytest

import gbdraw.circular as circular_cli
import gbdraw.linear as linear_cli
from gbdraw.core.text import dominant_baseline_shift_em, get_font_vertical_metrics, text_ink_offsets
from gbdraw.labels.label_style import (
    BOX_OPACITY,
    GAP_STOP_EM,
    HALO_COLOR,
    HALO_WIDTH_EM,
    UNDERLINE_GAP_EM,
    UNDERLINE_TRIM_EM,
    leader_width_px,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
INPUTS = REPO_ROOT / "tests" / "test_inputs"
LAMBDA = INPUTS / "NC_001416.gb"
SVG = "{http://www.w3.org/2000/svg}"
XLINK_HREF = "{http://www.w3.org/1999/xlink}href"
FEATURE = "data-label-feature-id"


def _svg(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, mode: str, *argv: str, fmt: str = "svg") -> Path:
    monkeypatch.chdir(tmp_path)
    name = f"{mode}-{abs(hash(argv))}"
    main = linear_cli.linear_main if mode == "linear" else circular_cli.circular_main
    main([*argv, "-o", name, "-f", fmt])
    return tmp_path / f"{name}.{fmt}"


def _walk(element: ET.Element, parents: tuple[ET.Element, ...] = ()) -> Iterator[tuple[ET.Element, tuple[ET.Element, ...]]]:
    yield element, parents
    for child in element:
        yield from _walk(child, (*parents, element))


def _label_texts(root: ET.Element) -> list[ET.Element]:
    return [element for element, _ in _walk(root) if element.tag == f"{SVG}text" and element.get(FEATURE)]


def _halos(root: ET.Element) -> list[ET.Element]:
    return [element for element, _ in _walk(root) if element.tag == f"{SVG}use"]


def _record_of(root: ET.Element) -> dict[int, ET.Element | None]:
    """Each element's record group (the nearest ``data-gbdraw-record-index`` or ``record_*`` ancestor)."""

    out: dict[int, ET.Element | None] = {}
    for element, parents in _walk(root):
        record = next(
            (parent for parent in reversed(parents) if parent.get("data-gbdraw-record-index") is not None),
            None,
        )
        out[id(element)] = record
    return out


def _check_halos(root: ET.Element) -> None:
    """One halo per drawn label with a leader: white, unfilled, 0.08 em wide, before its letters and after
    every leader, pointing to a text of its own record subtree; ids are unique."""

    ids = [element.get("id") for element, _ in _walk(root) if element.get("id")]
    assert len(ids) == len(set(ids))
    texts = {text.get("id"): text for text in _label_texts(root) if text.get("id")}
    halos = _halos(root)
    assert halos and len(halos) == len(texts)
    records = _record_of(root)
    for halo in halos:
        text = texts[halo.get(XLINK_HREF, "")[1:]]
        assert halo.get("fill") == "none" and halo.get("stroke") == "#ffffff"
        assert halo.get("stroke-linejoin") == "round" and halo.get(FEATURE) is None
        assert float(halo.get("stroke-width", "")) == pytest.approx(HALO_WIDTH_EM * float(text.get("font-size", "")))
        assert records[id(halo)] is records[id(text)]
    order = {id(element): index for index, (element, _) in enumerate(_walk(root))}
    for record in {id(records[id(halo)]): records[id(halo)] for halo in halos}.values():
        mine = [element for element, _ in _walk(root) if records[id(element)] is record and element.get(FEATURE)]
        letters = [order[id(text)] for text in texts.values() if records[id(text)] is record]
        own = [order[id(halo)] for halo in halos if records[id(halo)] is record]
        lines = [line for line in mine if line.tag == f"{SVG}line"]
        leaders = [order[id(line)] for line in lines if _is_leader(line)]
        casings = [order[id(line)] for line in lines if _is_casing(line)]
        unders = [order[id(line)] for line in lines if _is_underline(line)]
        # Leaders, then the halos, the casings, the leaders' tails, the underlines and the letters.
        tails = [index for index in leaders if index > min(own)]
        assert max(own) < min(letters) and min(own) > max(index for index in leaders if index < min(own))
        assert all(max(own) < index for index in casings + tails)
        assert max(casings, default=-1) < min(tails, default=len(order)) and max(tails, default=-1) < min(unders, default=len(order))
        assert max(unders, default=-1) < min(letters)


def _is_casing(line: ET.Element) -> bool:
    return line.get("stroke") == HALO_COLOR


def _is_underline(line: ET.Element) -> bool:
    return not _is_casing(line) and line.get("y1") == line.get("y2") and line.get("x1") != line.get("x2")


def _is_leader(line: ET.Element) -> bool:
    return not _is_casing(line) and not _is_underline(line)


@pytest.mark.parametrize(
    ("mode", "argv"),
    [
        ("linear", ("--gbk", str(LAMBDA), "--show_labels", "--label_protection", "halo")),
        (
            "linear",
            ("--gbk", str(INPUTS / "NC_001416.gb"), str(LAMBDA), "--show_labels", "all", "--separate_strands", "--label_protection", "halo"),
        ),
        ("circular", ("--gbk", str(LAMBDA), "--labels", "both", "--label_protection", "halo")),
        ("circular", ("--gbk", str(LAMBDA), str(LAMBDA), "--labels", "--multi_record_canvas", "--label_protection", "halo")),
    ],
    ids=["linear", "linear-two-records", "circular", "circular-multi-canvas"],
)
def test_each_label_with_a_leader_has_one_halo_behind_its_letters(
    mode: str, argv: tuple[str, ...], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    root = ET.parse(_svg(tmp_path, monkeypatch, mode, *argv)).getroot()
    _check_halos(root)
    # Embedded labels (drawn on their features) have no leader and no halo.
    for text in _label_texts(root):
        assert bool(text.get("id")) != _embedded(text, root)
    if mode == "circular":
        assert any(_embedded(text, root) for text in _label_texts(root)), "the case should embed some labels"


@pytest.mark.parametrize(
    ("mode", "argv"),
    [
        ("linear", ("--gbk", str(LAMBDA), "--show_labels")),
        ("circular", ("--gbk", str(LAMBDA), "--labels", "both", "--label_protection", "box")),
    ],
    ids=["linear-default", "circular"],
)
def test_protection_none_draws_no_halo_and_box_draws_a_box_per_label(
    mode: str, argv: tuple[str, ...], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Linear draws a box by default (LP-D106): one per label with a leader, white at ``BOX_OPACITY``, bound
    to its label so that hiding the label hides it, and, like a halo or a casing, taking no pointer events."""

    plain = ET.parse(_svg(tmp_path, monkeypatch, mode, *argv, "--label_protection", "none")).getroot()
    assert not _halos(plain)
    boxed = ET.parse(_svg(tmp_path, monkeypatch, mode, *argv)).getroot()
    boxes = [element for element, _ in _walk(boxed) if element.tag == f"{SVG}rect" and element.get("fill-opacity") == str(BOX_OPACITY)]
    with_leader = [text for text in _label_texts(boxed) if not _embedded(text, boxed)]
    assert not _halos(boxed) and len(boxes) == len(with_leader) > 0
    assert sorted(box.get(FEATURE, "") for box in boxes) == sorted(text.get(FEATURE, "") for text in with_leader)
    assert all(
        (box.get("fill"), box.get("pointer-events"), box.get("aria-hidden")) == ("#ffffff", "none", "true") for box in boxes
    )


def _embedded(text: ET.Element, root: ET.Element) -> bool:
    feature = text.get(FEATURE)
    return not any(
        line.get(FEATURE) == feature for line, _ in _walk(root) if line.tag == f"{SVG}line"
    )


def _lines_of(root: ET.Element, feature: str | None) -> list[ET.Element]:
    return [line for line, _ in _walk(root) if line.tag == f"{SVG}line" and line.get(FEATURE) == feature]


def _ink(text: ET.Element) -> tuple[float, float, float, float]:
    size = float(text.get("font-size", ""))
    offsets = text_ink_offsets(
        text.text or "", text.get("font-family", ""), size,
        dominant_baseline=text.get("dominant-baseline", ""), text_anchor=text.get("text-anchor", ""),
    )
    assert offsets is not None
    x, y = float(text.get("x", "")), float(text.get("y", ""))
    return x + offsets[0], x + offsets[1], y + offsets[2], y + offsets[3]


@pytest.mark.parametrize("font_size", ["5", "24"])
def test_a_linear_label_above_its_feature_is_underlined_and_one_below_stops_its_leader_above_it(
    font_size: str, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    argv = ("--gbk", str(LAMBDA), "--show_labels", "--label_font_size", font_size)
    above_svg = ET.parse(_svg(tmp_path, monkeypatch, "linear", *argv)).getroot()
    below_svg = ET.parse(_svg(tmp_path, monkeypatch, "linear", *argv, "--track_layout", "below")).getroot()
    size = float(font_size)
    width = leader_width_px(size)
    metrics = get_font_vertical_metrics(_label_texts(above_svg)[0].get("font-family", ""))
    assert metrics is not None
    shift = dominant_baseline_shift_em("central", metrics)
    assert shift is not None
    above = below = 0
    for root in (above_svg, below_svg):
        for text in _label_texts(root):
            lines = _lines_of(root, text.get(FEATURE))
            if not lines:
                continue
            left, right, top, bottom = _ink(text)
            y = float(text.get("y", ""))
            parts = [line for line in lines if _is_leader(line)]
            unders = [line for line in lines if _is_underline(line)]
            casings = [line for line in lines if _is_casing(line)]
            assert all(float(part.get("stroke-width", "")) == pytest.approx(width) for part in parts)
            if float(parts[0].get("y1", "")) > y:
                # Above its feature: the leader meets the line under the letters at a T.
                above += 1
                (under,), (casing,) = unders, casings
                assert under.get("stroke") == parts[0].get("stroke")
                assert float(under.get("stroke-width", "")) == pytest.approx(width)
                line_y = float(under.get("y1", ""))
                descent_y = y + (shift + metrics.descent) * size
                assert line_y - width / 2 - descent_y == pytest.approx(UNDERLINE_GAP_EM * size)
                # Each end stops 0.12 em inside the letters (at most a quarter of the line).
                trim = min(UNDERLINE_TRIM_EM * size, 0.25 * (right - left))
                ends = sorted((float(under.get("x1", "")), float(under.get("x2", ""))))
                assert ends == pytest.approx([left + trim, right - trim])
                # The leader's parts (up to its casing, then its tail) join, and it ends on the line.
                points = [((float(p.get("x1", "")), float(p.get("y1", ""))), (float(p.get("x2", "")), float(p.get("y2", "")))) for p in parts]
                assert all(first[1] == pytest.approx(second[0]) for first, second in zip(points, points[1:]))
                assert points[-1][1][1] == pytest.approx(line_y) and ends[0] - 1e-6 <= points[-1][1][0] <= ends[1] + 1e-6
                # The casing: the line's width plus 0.08 em each side, each end 0.08 em and one line width out.
                margin = HALO_WIDTH_EM * size
                assert casing.get("stroke-linecap") == "butt" and casing.get("y1") == casing.get("y2") == under.get("y1")
                assert float(casing.get("stroke-width", "")) == pytest.approx(width + 2 * margin)
                assert sorted((float(casing.get("x1", "")), float(casing.get("x2", "")))) == pytest.approx(
                    [ends[0] - margin - width, ends[1] + margin + width]
                )
            else:
                # Below its feature: the leader would cross the letters, so it stops above them; no line.
                below += 1
                assert not unders and not casings and len(parts) == 1
                assert float(parts[0].get("y2", "")) == pytest.approx(top - GAP_STOP_EM * size)
    assert above and below


def test_the_linear_underline_takes_no_room_in_the_label_rows(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """The underline is not counted in a Linear label's row box: every label and the canvas are where
    ``--label_line none`` puts them, and each line still clears the letters of the next row."""

    base = ("--gbk", str(LAMBDA), "--show_labels")
    roots = [ET.parse(_svg(tmp_path, monkeypatch, "linear", *base, "--label_line", line)).getroot() for line in ("none", "under")]
    places = [sorted((text.text, text.get("x"), text.get("y")) for text in _label_texts(root)) for root in roots]
    without, root = roots
    assert len(places[0]) > 2 and places[0] == places[1]
    assert (without.get("width"), without.get("height")) == (root.get("width"), root.get("height"))
    rows = sorted({float(text.get("y", "")) for text in _label_texts(root) if not _embedded(text, root)})
    for text in _label_texts(root):
        unders = [line for line in _lines_of(root, text.get(FEATURE)) if _is_underline(line)]
        row = rows.index(float(text.get("y", "")))
        if unders and row + 1 < len(rows):
            # The line clears the letters of the next row toward the track.
            nearer = [other for other in _label_texts(root) if float(other.get("y", "")) == rows[row + 1]]
            line_bottom = float(unders[0].get("y1", "")) + float(unders[0].get("stroke-width", "")) / 2
            assert line_bottom < min(_ink(other)[2] for other in nearer)


@pytest.mark.parametrize(
    ("mode", "argv"),
    [
        ("linear", ("--gbk", str(LAMBDA), "--show_labels", "--label_placement", "above_feature", "--label_rotation", "45")),
        (
            "linear",
            (
                "--gbk", str(LAMBDA), "--show_labels", "all", "--separate_strands",
                "--label_placement", "above_feature", "--label_rotation", "45",
            ),
        ),
        ("circular", ("--gbk", str(LAMBDA), "--labels", "--label_placement", "radial")),
    ],
    ids=["linear", "linear-separate-strands", "circular-radial"],
)
def test_a_turned_label_draws_no_underline_and_its_leader_ends_as_without_one(
    mode: str, argv: tuple[str, ...], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """LP-D112: a feature label drawn turned (not horizontal) gets no underline and no casing, and its leader
    and its box are where ``--label_line none`` draws them, so the whole figure is the same."""

    under, none = (
        _svg(tmp_path, monkeypatch, mode, *argv, "--label_line", line) for line in ("under", "none")
    )
    root = ET.parse(under).getroot()
    turned = [text for text in _label_texts(root) if "rotate(" in text.get("transform", "") and _lines_of(root, text.get(FEATURE))]
    assert len(turned) > 2
    assert under.read_text(encoding="utf-8") == none.read_text(encoding="utf-8")


@pytest.mark.parametrize(
    "argv",
    [
        ("--show_labels", "--label_rotation", "45"),
        ("--show_labels", "all", "--label_rotation", "30", "--separate_strands"),
        ("--show_labels", "all", "--label_rotation", "-60", "--label_font_size", "24"),
    ],
    ids=["45", "all-30-separate-strands", "all-minus-60-large"],
)
def test_every_leader_of_a_turned_linear_label_is_vertical_in_the_default_style(
    argv: tuple[str, ...], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """LP-D112: the leader of a turned Linear label runs straight up or down from its feature to the letters, as
    a horizontal label's does; no part of it is slanted or bent along the text."""

    root = ET.parse(
        _svg(tmp_path, monkeypatch, "linear", "--gbk", str(LAMBDA), "--label_placement", "above_feature", *argv)
    ).getroot()
    turned = [text for text in _label_texts(root) if "rotate(" in text.get("transform", "")]
    checked = 0
    for text in turned:
        for line in _lines_of(root, text.get(FEATURE)):
            assert not _is_casing(line), text.text
            x1, _y1, x2, _y2 = _line_points(line)
            assert x1 == pytest.approx(x2, abs=1e-9), text.text
            checked += 1
    assert checked > 2


def _line_points(line: ET.Element) -> tuple[float, float, float, float]:
    return tuple(float(line.get(key, "")) for key in ("x1", "y1", "x2", "y2"))  # type: ignore[return-value]


def _uncut_crossings(root: ET.Element) -> list[tuple[str, str]]:
    """Each (underlined label, other label) where a part of the other label's leader crosses the underline
    and is not cut there by the underline's casing: drawn after the casing, or past its ends. A leader within
    half a pixel of the underlined label's own leader lies on it and is not counted."""

    order = {id(element): index for index, (element, _) in enumerate(_walk(root))}
    parent = {id(child): element for element, _ in _walk(root) for child in element}
    lines = [element for element, _ in _walk(root) if element.tag == f"{SVG}line" and element.get(FEATURE)]
    found = []
    for under in (line for line in lines if _is_underline(line)):
        owner = under.get(FEATURE)
        (casing,) = [line for line in lines if _is_casing(line) and line.get(FEATURE) == owner]
        own_x = [_line_points(line)[0] for line in lines if _is_leader(line) and line.get(FEATURE) == owner]
        ux1, line_y, ux2, _ = _line_points(under)
        cx1, _, cx2, _ = _line_points(casing)
        for part in lines:
            if not _is_leader(part) or part.get(FEATURE) == owner or parent[id(part)] is not parent[id(under)]:
                continue
            x1, y1, x2, y2 = _line_points(part)
            if not min(y1, y2) < line_y < max(y1, y2):
                continue
            x = x1 + (x2 - x1) * (line_y - y1) / (y2 - y1)
            half = float(part.get("stroke-width", "")) / 2
            if not min(ux1, ux2) - half < x < max(ux1, ux2) + half or any(abs(x - own) < 0.5 for own in own_x):
                continue
            cut = order[id(part)] < order[id(casing)] and min(cx1, cx2) <= x - half and x + half <= max(cx1, cx2)
            if not cut:
                found.append((str(owner), str(part.get(FEATURE))))
    return found


@pytest.mark.parametrize("font_size", ["8", "24"])
def test_a_casing_cuts_other_leaders_at_an_underline(font_size: str, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Where another label's leader crosses an underline, the underline's casing cuts it, so it never reads as
    a joint (LP-D103); without casings the crossings are joints."""

    root = ET.parse(_svg(tmp_path, monkeypatch, "linear", "--gbk", str(LAMBDA), "--show_labels", "--label_font_size", font_size)).getroot()
    assert _uncut_crossings(root) == []
    for casing in [line for line, _ in _walk(root) if line.tag == f"{SVG}line" and _is_casing(line)]:
        casing.set("x2", casing.get("x1", ""))
    assert font_size == "8" or len(_uncut_crossings(root)) > 50


def _raster(root: ET.Element, drop) -> "object":
    import io

    cairosvg = pytest.importorskip("cairosvg")
    from PIL import Image

    copy = ET.fromstring(ET.tostring(root))
    for element, parents in list(_walk(copy)):
        if parents and drop(element):
            parents[-1].remove(element)
    return Image.open(io.BytesIO(cairosvg.svg2png(bytestring=ET.tostring(copy), scale=2.0, background_color="white"))).convert("L")


@pytest.mark.parametrize(
    "argv",
    [
        ("--gbk", str(LAMBDA), "--show_labels", "--label_font_size", "24"),
        ("--gbk", str(INPUTS / "HmmtDNA.gbk"), "--show_labels", "all", "--track_layout", "above", "--separate_strands"),
    ],
    ids=["lambda-24pt", "hmmtdna"],
)
def test_a_casing_hides_nothing_but_leaders(argv: tuple[str, ...], tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """With every leader taken out, the figure with casings and the one without are the same picture: a casing
    hides no letters, underlines, feature shapes or other marks."""

    from PIL import ImageChops

    root = ET.parse(_svg(tmp_path, monkeypatch, "linear", *argv)).getroot()
    assert [line for line, _ in _walk(root) if line.tag == f"{SVG}line" and _is_casing(line)]

    def leader(element: ET.Element) -> bool:
        return element.tag == f"{SVG}line" and bool(element.get(FEATURE)) and _is_leader(element)

    with_casings = _raster(root, leader)
    without = _raster(root, lambda element: leader(element) or (element.tag == f"{SVG}line" and _is_casing(element)))
    # A pixel the casings make lighter is hidden ink.
    assert ImageChops.subtract(with_casings, without).getbbox() is None
    # With the leaders kept, the casings do cut them.
    casing = lambda element: element.tag == f"{SVG}line" and _is_casing(element)  # noqa: E731
    assert ImageChops.subtract(_raster(root, lambda element: False), _raster(root, casing)).getbbox() is not None


def test_a_circular_underline_continues_its_leader_without_a_new_bend(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    base = ("--gbk", str(LAMBDA), "--labels", "both")
    default = ET.parse(_svg(tmp_path, monkeypatch, "circular", *base)).getroot()
    under = ET.parse(_svg(tmp_path, monkeypatch, "circular", *base, "--label_line", "under")).getroot()

    def groups(root: ET.Element) -> tuple[ET.Element, ET.Element]:
        found = {element.get("id"): element for element, _ in _walk(root) if element.get("id") in {"label_leaders", "label_text"}}
        return found["label_leaders"], found["label_text"]

    leaders, text_layer = groups(default)
    assert not [child for child in text_layer if child.tag == f"{SVG}line"]
    under_leaders, under_text = groups(under)
    # Two segments per label, as before: the last one is re-aimed at the line's end.
    assert len(list(under_leaders)) == len(list(leaders))
    ends = {}
    for line in under_text:
        if line.tag == f"{SVG}line":
            ends.setdefault(line.get(FEATURE), []).append(
                {(float(line.get("x1", "")), float(line.get("y1", ""))), (float(line.get("x2", "")), float(line.get("y2", "")))}
            )
    assert ends
    for leader in under_leaders:
        if leader.get(FEATURE) in ends:
            end = (float(leader.get("x2", "")), float(leader.get("y2", "")))
            if any(end in points for points in ends[leader.get(FEATURE)]):
                break
    else:
        pytest.fail("no Circular leader ends on its underline")


def _pdf_text_runs(path: Path) -> int:
    data = path.read_bytes()
    runs = 0
    for match in re.finditer(rb"<<(.*?)>>\s*stream\r?\n", data, re.S):
        start = match.end()
        stream = data[start:data.find(b"endstream", start)]
        if b"FlateDecode" in match.group(1):
            try:
                stream = zlib.decompress(stream)
            except zlib.error:
                continue
        runs += len(re.findall(rb"\)\s*Tj|\]\s*TJ|>\s*Tj", stream))
    return runs


def test_a_cairosvg_pdf_carries_each_haloed_label_once(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    pytest.importorskip("cairosvg")
    base = ("--gbk", str(LAMBDA), "--show_labels")
    with_halos = _svg(tmp_path, monkeypatch, "linear", *base, "--label_protection", "halo", fmt="pdf")
    without = _svg(tmp_path, monkeypatch, "linear", *base, "--label_protection", "none", fmt="pdf")
    assert _pdf_text_runs(with_halos) == _pdf_text_runs(without) > 0


def _feature_y_extent(root: ET.Element, feature: str | None) -> tuple[float, float]:
    """The top and bottom y of a feature's drawn shapes: their path points and half their outline."""

    tops, bottoms = [], []
    for element, _ in _walk(root):
        if element.get("data-gbdraw-feature-id") == feature and element.tag == f"{SVG}path":
            numbers = [float(value) for value in re.findall(r"-?\d+(?:\.\d+)?(?:e-?\d+)?", element.get("d", ""))]
            half = 0.5 * float(element.get("stroke-width", "0") or 0)
            tops.append(min(numbers[1::2]) - half)
            bottoms.append(max(numbers[1::2]) + half)
    assert tops, feature
    return min(tops), max(bottoms)


def _box_corners(box: ET.Element) -> list[tuple[float, float]]:
    """The four corners of a ``rect``, turned by its ``rotate(angle,cx,cy)`` transform."""

    x, y, width, height = (float(box.get(key, "")) for key in ("x", "y", "width", "height"))
    corners = [(x, y), (x + width, y), (x + width, y + height), (x, y + height)]
    match = re.fullmatch(r"rotate\(([-\d.e]+),([-\d.e]+),([-\d.e]+)\)", box.get("transform", "") or "rotate(0,0,0)")
    assert match is not None, box.get("transform")
    angle, cx, cy = (float(value) for value in match.groups())
    cos, sin = math.cos(math.radians(angle)), math.sin(math.radians(angle))
    return [(cx + (px - cx) * cos - (py - cy) * sin, cy + (px - cx) * sin + (py - cy) * cos) for px, py in corners]


@pytest.mark.parametrize(
    "argv",
    [
        ("--show_labels", "--label_rotation", "45"),
        ("--show_labels", "all", "--label_rotation", "30", "--separate_strands"),
        ("--show_labels", "all", "--label_rotation", "-60", "--label_font_size", "24"),
    ],
    ids=["45", "all-30-separate-strands", "all-minus-60-large"],
)
def test_a_turned_label_on_its_feature_gets_its_box_which_stays_off_the_feature(
    argv: tuple[str, ...], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """OV-1361 (LP-D119): a turned Linear label drawn on its feature (no leader) also gets its box, so another
    label's leader that passes through its letters is faded there. The box stays on the label's side of the
    feature's outline, and every box lies under every label's letters."""

    root = ET.parse(
        _svg(tmp_path, monkeypatch, "linear", "--gbk", str(LAMBDA), "--label_placement", "above_feature", *argv)
    ).getroot()
    order = {id(element): index for index, (element, _) in enumerate(_walk(root))}
    boxes = {
        element.get(FEATURE): element
        for element, _ in _walk(root)
        if element.tag == f"{SVG}rect" and element.get("fill-opacity") == str(BOX_OPACITY)
    }
    texts = _label_texts(root)
    turned = [text for text in texts if "rotate(" in text.get("transform", "")]
    on_feature = [text for text in turned if not _lines_of(root, text.get(FEATURE))]
    assert len(on_feature) > 2
    assert sorted(boxes) == sorted(text.get(FEATURE) for text in turned)
    for text in on_feature:
        feature = text.get(FEATURE)
        top, bottom = _feature_y_extent(root, feature)
        ys = [y for _x, y in _box_corners(boxes[feature])]
        if float(text.get("y", "")) < top:
            assert max(ys) <= top + 1e-6, text.text
        else:
            assert min(ys) >= bottom - 1e-6, text.text
    assert max(order[id(box)] for box in boxes.values()) < min(order[id(text)] for text in texts)


@pytest.mark.parametrize(
    "argv",
    [("--show_labels",), ("--show_labels", "all", "--separate_strands")],
    ids=["one-side", "all-separate-strands"],
)
def test_a_horizontal_above_feature_label_has_its_leader_underline_and_box_as_in_auto(
    argv: tuple[str, ...], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """LP-D120 (OV-1396): every horizontal ``above_feature`` label, the lowest row too, is drawn as an auto
    label in its rows: a leader from its feature, the underline it ends on (a label below its feature: the
    leader stops above the letters, no line), and its white box. Before, the lowest row had none of the
    three (lambda: "terminase small subunit", "head decoration")."""

    root = ET.parse(
        _svg(tmp_path, monkeypatch, "linear", "--gbk", str(LAMBDA), "--label_placement", "above_feature", *argv)
    ).getroot()
    boxes = {
        element.get(FEATURE)
        for element, _ in _walk(root)
        if element.tag == f"{SVG}rect" and element.get("fill-opacity") == str(BOX_OPACITY)
    }
    texts = _label_texts(root)
    assert {"terminase small subunit", "head decoration"} <= {text.text for text in texts}
    above = below = 0
    for text in texts:
        feature = text.get(FEATURE)
        lines = _lines_of(root, feature)
        leaders = [line for line in lines if _is_leader(line)]
        unders = [line for line in lines if _is_underline(line)]
        top, bottom = _feature_y_extent(root, feature)
        assert leaders and feature in boxes, text.text
        start_y = float(leaders[0].get("y1", ""))
        if float(text.get("y", "")) < top:
            above += 1
            assert start_y == pytest.approx(top, abs=1.0) and len(unders) == 1, text.text
        else:
            below += 1
            assert start_y == pytest.approx(bottom, abs=1.0) and not unders, text.text
    assert above > 5 and (below > 5) is ("--separate_strands" in argv)
