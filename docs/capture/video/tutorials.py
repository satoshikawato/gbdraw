"""Record the Gallery tutorials as human-paced Web walkthroughs.

Each journey follows one tutorial in ``gbdraw/web/gallery/tutorials`` from an
empty app to the exported SVG. Captions use the tutorial's step numbers and
its UI labels, so a viewer can follow the written steps beside the video.
Generate is pressed at a few extra checkpoints to show each change as it
lands; the final settings are the tutorial's own.
"""

from __future__ import annotations

import re
from dataclasses import replace
from pathlib import Path

from playwright.sync_api import Locator, expect

from config import ACTION_TIMEOUT_MS, GENERATION_TIMEOUT_MS, REPO_ROOT, tutorial_file_identity
from flows.web_capture import assert_fixture_identity
from video.walkthrough import (PORTRAIT_FIGURE_WIDTH, _PREVIEW, Journey, Recorder, _center, export_interactive_svg,
                               open_in_browser)
from video.walkthrough_render import VERTICAL


def _fixture(file_id: str) -> Path:
    return tutorial_file_identity(file_id)[0]


def _verify(*file_ids: str) -> None:
    """Record only the frozen, checksum-verified tutorial inputs."""

    for file_id in file_ids:
        path, size, digest = tutorial_file_identity(file_id)
        assert_fixture_identity(path, expected_size=size, expected_sha256=digest)


def _labels(svg: Path) -> set[str]:
    return set(re.findall(r">([^<>]{1,80})</(?:text|textPath|tspan)>", svg.read_text(encoding="utf-8")))


def _fills(svg: Path) -> set[str]:
    return {fill.lower() for fill in re.findall(r'fill="(#[0-9a-fA-F]{6})"', svg.read_text(encoding="utf-8"))}


def _require(report: dict, ok: bool, svg: Path) -> dict:
    if not ok:
        raise AssertionError(f"{svg.name} lacks the tutorial's edits: {report}")
    return report


_HIT = """
(selector) => {
  // In an opened Interactive SVG the drawing is the whole page.
  const region = document.querySelector('[role="region"][aria-label="Result Preview"]') || document.documentElement;
  const view = region.getBoundingClientRect();
  for (const node of region.querySelectorAll(selector)) {
    const b = node.getBoundingClientRect();
    if (!b.width && !b.height) continue;
    const points = [];
    for (let i = 1; i < 10; i++) for (let j = 1; j < 10; j++) points.push([b.left + b.width * j / 10, b.top + b.height * i / 10]);
    const cx = b.left + b.width / 2, cy = b.top + b.height / 2;
    points.sort((p, q) => Math.hypot(p[0] - cx, p[1] - cy) - Math.hypot(q[0] - cx, q[1] - cy));
    for (const [x, y] of points) {
      if (x < view.left + 8 || x > view.right - 8 || y < view.top + 8 || y > view.bottom - 8) continue;
      const hit = document.elementFromPoint(x, y);
      if (hit === node || node.contains(hit)) return {x, y, box: {x: b.left, y: b.top, width: b.width, height: b.height}};
    }
  }
  return null;
}
"""


# -- shared steps -------------------------------------------------------------

def _select_mode(rec: Recorder, name: str) -> None:
    button = rec.page.get_by_role("button", name=name, exact=True)
    rec.focus(button, zoom=2.2)
    rec.hold(0.4)
    rec.click(button, pause=0.6)
    expect(button).to_have_attribute("aria-pressed", "true")


def _fit_preview(rec: Recorder, *, margin: float = 24) -> None:
    """Zoom out with the preview's own control until the figure fits above its buttons, then center it.

    A portrait take also narrows the figure until the tall window can show all
    of it. It frames the canvas meanwhile and skims the repeated clicks under
    the fast-forward badge, so the new figure is not seen cut off.
    """

    page = rec.page
    zoom_out = page.get_by_role("button", name="Zoom out", exact=True)
    skimming = False
    for _ in range(16):
        geometry = page.evaluate(_PREVIEW)
        canvas, svg = geometry["canvas"], geometry["svg"]
        width = min(canvas["width"] - margin, PORTRAIT_FIGURE_WIDTH) if rec.portrait else canvas["width"] - margin
        if svg["width"] <= width and svg["height"] <= canvas["height"] - margin:
            break
        if rec.portrait and not skimming:
            rec.camera(canvas, duration=0.5)
            rec.fast_forward(factor=3.0)
            skimming = True
        rec.click(zoom_out, pause=0.35)
    else:
        raise AssertionError(f"The preview did not fit the figure: {geometry}")
    page.wait_for_timeout(300)
    geometry = page.evaluate(_PREVIEW)
    canvas, svg = geometry["canvas"], geometry["svg"]
    dx = canvas["x"] + canvas["width"] / 2 - (svg["x"] + svg["width"] / 2)
    dy = canvas["y"] + canvas["height"] / 2 - (svg["y"] + svg["height"] / 2)
    if abs(dx) > 24 or abs(dy) > 24:
        # Leaving the canvas ends a pan, so the drag must end inside it too.
        inside = [(x, y) for x, y in geometry["blank"]
                  if canvas["x"] + 16 <= x + dx <= canvas["x"] + canvas["width"] - 16
                  and canvas["y"] + 16 <= y + dy <= canvas["y"] + canvas["height"] - 16]
        if not inside:
            raise AssertionError(f"No blank preview area to pan from: {geometry}")
        start = inside[0]
        rec.drag(start, (start[0] + dx, start[1] + dy))
        page.wait_for_timeout(300)
        geometry = page.evaluate(_PREVIEW)
        canvas, svg = geometry["canvas"], geometry["svg"]
        off = canvas["x"] + canvas["width"] / 2 - (svg["x"] + svg["width"] / 2)
        if abs(off) > 40:
            raise AssertionError(f"The preview pan did not center the figure: {geometry}")
    if skimming:
        rec.normal_speed()


def _generate(rec: Recorder, *, target: float = 1.4, look: float = 2.0,
              timeout_ms: int = GENERATION_TIMEOUT_MS) -> None:
    button = rec.page.get_by_role("button", name="Generate Diagram", exact=True)
    rec.focus(button, zoom=1.9)
    rec.hold(0.4)
    rec.generate(target=target, timeout_ms=timeout_ms)
    _fit_preview(rec)
    rec.camera(duration=1.0)
    rec.hold(look)


def _field(rec: Recorder, field: Locator, text: str, *, zoom: float = 2.2) -> None:
    rec.reveal(field)
    rec.focus(field, zoom=zoom, duration=0.6)
    rec.type(field, text)


def _value(rec: Recorder, field: Locator, text: str, *, zoom: float = 2.2) -> None:
    """Type ``text``, or only point at the field when it already holds it."""

    if field.input_value() != text:
        _field(rec, field, text, zoom=zoom)
        return
    _point(rec, field, zoom=zoom)


def _point(rec: Recorder, control: Locator, *, zoom: float = 2.2) -> None:
    """Point at a control that already holds the tutorial's value."""

    rec.reveal(control)
    rec.focus(control, zoom=zoom, duration=0.5)
    rec.move_to(*_center(control.bounding_box()))
    rec.hold(0.5)


def _pick(rec: Recorder, select: Locator, value: str, *, zoom: float = 2.2) -> None:
    rec.reveal(select)
    rec.focus(select, zoom=zoom, duration=0.6)
    rec.choose(select, value)


def _press(rec: Recorder, button: Locator, *, zoom: float = 2.2, pause: float = 0.5) -> None:
    rec.reveal(button)
    rec.focus(button, zoom=zoom, duration=0.6)
    rec.click(button, pause=pause)


def _look(rec: Recorder, box: dict, *, zoom: float, seconds: float, pad: float = 0) -> None:
    if pad:
        box = {"x": box["x"] - pad, "y": box["y"] - pad,
               "width": box["width"] + pad * 2, "height": box["height"] + pad * 2}
    rec.camera(box, zoom=zoom, duration=1.1)
    rec.hold(seconds)


def _priority_rule(rec: Recorder, feature: str, order: str) -> None:
    """Labels › Qualifier priority: one feature type, one qualifier order."""

    field = rec.page.locator('input[placeholder="product,gene,locus_tag"]')
    rec.reveal(field)
    rec.focus(field, zoom=2.2)
    rec.choose(field.locator("xpath=preceding-sibling::select[1]"), feature)
    rec.type(field, order)
    rec.click(field.locator("xpath=following-sibling::button[1]"), pause=0.6)


def _color_rules(rec: Recorder, rules: tuple[tuple[str, str, str, str, str], ...]) -> None:
    """Colors › Specific rules; the first is shown at pace, the rest at 3×."""

    regex = rec.page.get_by_label("New color rule pattern", exact=True)
    rule_box = regex.locator('xpath=ancestor::div[contains(@class,"bg-slate-50")][1]')
    for index, (feature, qualifier_name, pattern, color, legend) in enumerate(rules):
        rec.reveal(rule_box, anchor=0.45)
        rec.mark(f"rule-{index + 1}")
        rec.focus(rule_box, zoom=2.2, duration=0.6)
        if index == 1:
            rec.fast_forward(factor=3.0)
        rec.choose(rule_box.get_by_label("New color rule feature type", exact=True), feature)
        qualifier = rule_box.get_by_label("Qualifier name", exact=True)
        if qualifier.input_value() != qualifier_name:
            rec.type(qualifier, qualifier_name)
        rec.type(regex, pattern)
        rec.pick_color(rule_box.get_by_label("New color rule color", exact=True), color)
        rec.type(rule_box.get_by_label("New color rule legend caption", exact=True), legend)
        rec.click(rule_box.get_by_role("button", name="Add color rule", exact=True), pause=0.45)
        # Each rule is validated before it is committed; a quicker next rule
        # would be validated without this one.
        rec.page.wait_for_function(
            "(count) => (window.__GBDRAW_APP__?.manualSpecificRules || []).length === count", arg=index + 1)
    rec.normal_speed()
    rec.hold(0.5)


def _feature_types(rec: Recorder, *, remove: tuple[str, ...] = (), add: tuple[str, ...] = ()) -> None:
    """Features › Feature Types: drop and add the drawn feature types."""

    page = rec.page
    summary = page.locator('summary[aria-label="Features"]')
    rec.open_section(summary)
    card = summary.locator("xpath=..")
    for name in remove:
        chip = card.locator("span.inline-flex").filter(
            has=page.locator("span.font-medium", has_text=re.compile(rf"^{re.escape(name)}$")))
        _press(rec, chip.locator("button"), zoom=2.2, pause=0.5)
        expect(chip).to_have_count(0)
    for name in add:
        _pick(rec, card.locator(f'select:has(option[value="{name}"])'), name)
        # Icon glyphs precede the button text in the accessible name.
        rec.click(card.get_by_role("button", name=re.compile(r"\bAdd$")), pause=0.5)
    features = page.evaluate("() => [...(window.__GBDRAW_APP__?.adv?.features || [])]")
    if set(remove) & set(features) or set(add) - set(features):
        raise AssertionError(f"Feature types were not updated: {features}")


def _color_value(rec: Recorder, label: str, color: str, *, zoom: float = 2.2) -> None:
    """Switch an Auto / None / Color control to Color and pick ``color``."""

    page = rec.page
    _pick(rec, page.get_by_label(f"{label} mode", exact=True), "color", zoom=zoom)
    rec.pick_color(page.get_by_label(label, exact=True), color)


def _hit(rec: Recorder, selector: str) -> dict:
    """A visible point that really lands on a drawn element matching ``selector``."""

    point = rec.page.evaluate(_HIT, selector)
    if point is None:
        raise AssertionError(f"No clickable preview element matches {selector}")
    return point


def _svg_id(rec: Recorder, feature_filter: str, arg: object) -> str:
    """The drawn id of the feature that ``feature_filter`` picks; exports keep it."""

    svg_id = rec.page.evaluate(
        f"(arg) => ((window.__GBDRAW_APP__?.extractedFeatures || []).find({feature_filter}) || {{}}).svg_id", arg)
    if not svg_id:
        raise AssertionError(f"No feature matches {feature_filter} ({arg!r})")
    return svg_id


def _feature_point(rec: Recorder, feature_filter: str, arg: object) -> dict:
    """A clickable point on the drawn feature that ``feature_filter`` picks."""

    return _hit(rec, f'svg [data-gbdraw-feature-id="{_svg_id(rec, feature_filter, arg)}"]')


def _feature_tab(rec: Recorder, popup: Locator, tab: str, *, pause: float) -> None:
    # Icon glyphs precede the tab text in the accessible name.
    rec.click(popup.get_by_role("button", name=re.compile(rf"\b{tab}$")), pause=pause)


def _open_feature(rec: Recorder, point: dict, tab: str) -> Locator:
    """Click a drawn feature, switch its popup to ``tab``, and frame both."""

    page = rec.page
    rec.click_point(point["x"], point["y"], pause=0.4)
    popup = page.get_by_role("dialog", name=re.compile(r"^Feature details"))
    expect(popup).to_be_visible()
    _feature_tab(rec, popup, tab, pause=0.5)
    box = popup.bounding_box()
    left = min(box["x"], point["x"]) - 40
    right = max(box["x"] + box["width"], point["x"]) + 40
    # The tall window cannot hold both, so it reads the popup alone.
    rec.camera({"x": left, "y": box["y"] - 40, "width": right - left, "height": box["height"] + 80},
               zoom=1.55, duration=1.0, portrait={"x": box["x"] - 30, "y": box["y"] - 30,
                                                  "width": box["width"] + 60, "height": box["height"] + 60})
    return popup


def _scroll_popup_to(rec: Recorder, popup: Locator, text: str) -> None:
    """Wheel-scroll the popup until the ``text`` row sits in its lower part."""

    target = popup.get_by_text(text, exact=True).first
    box = popup.bounding_box()
    rec.move_to(box["x"] + box["width"] * 0.6, box["y"] + box["height"] * 0.6)
    for _ in range(60):
        row = target.bounding_box()
        if row["y"] + row["height"] <= box["y"] + box["height"] * 0.8:
            break
        rec.page.mouse.wheel(0, 40)
        rec.page.wait_for_timeout(40)
    else:
        raise AssertionError(f"Could not scroll {text!r} into the feature popup")
    rec.hold(0.3)


def _close_feature(rec: Recorder) -> None:
    rec.click(rec.page.get_by_label("Close feature popup", exact=True), pause=0.3)
    rec.camera(duration=0.9)
    rec.hold(0.5)


def _open_custom_slots(rec: Recorder) -> None:
    page = rec.page
    rec.open_section(page.get_by_label("Layout", exact=True))
    _press(rec, page.get_by_role("button", name="Custom Track Slots", exact=True), zoom=2.0)
    rec.set_checked(page.get_by_role("checkbox", name="Use custom stack", exact=True), True)


def _slot(rec: Recorder, slot_id: str) -> Locator:
    return rec.page.get_by_role("group", name=f"Circular track slot {slot_id}", exact=True)


def _run_losat(rec: Recorder) -> None:
    """Comparison › Run LOSAT for all adjacent pairs; the pressed button confirms the choice."""

    button = rec.page.get_by_role("button", name="Run LOSAT for all adjacent pairs", exact=True)
    _press(rec, button, zoom=2.0)
    expect(button).to_have_attribute("aria-pressed", "true")


def _slot_control(slot: Locator, name: str) -> Locator:
    """A control in a track-slot row; its accessible name starts with the row's slot id."""

    return slot.get_by_label(re.compile(rf" {re.escape(name)}$"))


def _measure(rec: Recorder, slot: Locator, name: str, value: str, unit: str = "factor") -> None:
    """Type a Width or Radius value in ``unit`` (px or factor, shown as ×R).

    The unit is picked first: on an Auto field it applies to the next input.
    """

    select = _slot_control(slot, f"{name} unit")
    if select.input_value() != unit:
        _pick(rec, select, unit, zoom=2.0)
    _field(rec, _slot_control(slot, f"{name} value"), value, zoom=2.0)


def _slot_ids(rec: Recorder) -> list[str]:
    return rec.page.evaluate("() => (window.__GBDRAW_APP__?.adv?.circular_track_slots || []).map((s) => s.id)")


# -- the exported Interactive SVG ---------------------------------------------

_SEARCH_STATUS = """
() => document.querySelector('foreignObject#gbdraw-feature-search-controls')?.textContent || ''
"""


def _open_export(rec: Recorder, svg: Path, step: int) -> list[str]:
    rec.caption(step, "Open it in any web browser", "The Interactive SVG needs no server and no install")
    errors = open_in_browser(rec, svg)
    rec.mark("interactive")
    rec.camera(duration=0.6)
    rec.move_to(1500, 760)
    rec.hold(2.2)
    return errors


def _frame_popup(rec: Recorder, popup: Locator) -> None:
    """Frame the opened file's popup alone, large enough to read its rows."""

    box = popup.bounding_box()
    rec.camera({"x": box["x"] - 30, "y": box["y"] - 30, "width": box["width"] + 60, "height": box["height"] + 60},
               zoom=1.7, duration=1.0)


def _svg_popup(rec: Recorder, point: dict, *, overview: float = 0) -> Locator:
    """Click ``point`` in the opened file and frame the popup it opens.

    ``overview`` first holds on the whole drawing, where the clicked feature's
    similarity-group members light up.
    """

    rec.click_point(point["x"], point["y"], pause=0.5)
    popup = rec.page.locator("foreignObject#gbdraw-feature-popup")
    expect(popup).to_be_visible()
    if overview:
        rec.camera(duration=0.8)
        rec.hold(overview)
    _frame_popup(rec, popup)
    return popup


def _svg_tab(rec: Recorder, popup: Locator, name: str, *, seconds: float) -> None:
    rec.click(popup.get_by_role("button", name=name, exact=True), pause=0.4)
    rec.hold(seconds)


def _svg_close(rec: Recorder, popup: Locator) -> None:
    rec.click(popup.get_by_role("button", name="x", exact=True), pause=0.3)
    expect(popup).to_be_hidden()
    rec.camera(duration=0.9)
    rec.hold(0.4)


def _svg_feature(rec: Recorder, svg_id: str, step: int, *,
                 first: str = "Details", overview: float = 0, shows: str = "") -> Locator:
    """Show one feature's popup in the opened file, then its embedded sequences."""

    popup = _svg_popup(rec, _hit(rec, f'[data-gbdraw-feature-id="{svg_id}"]'), overview=overview)
    _svg_tab(rec, popup, first, seconds=2.6)
    if shows:
        expect(popup).to_contain_text(shows)
    _svg_sequences(rec, popup, step)
    return popup


def _svg_sequences(rec: Recorder, popup: Locator, step: int) -> None:
    rec.caption(step, "Its sequences are in the file", "Nucleotide and amino acid FASTA, ready to copy or download")
    _svg_tab(rec, popup, "Sequence", seconds=3.2)
    if not re.search(r"NUCLEOTIDE|Nucleotide", popup.evaluate("node => node.textContent")):
        raise AssertionError("The Sequence tab shows no nucleotide sequence")


def _svg_search(rec: Recorder, query: str, step: int, *, sequences: bool = False,
                qualifier: str | None = None) -> None:
    """Search › Label in the opened file, open the first match, then return to the whole drawing.

    ``sequences`` also shows the opened match's Sequence tab. ``qualifier``
    searches Qualifier value in that key instead of the label.
    """

    page = rec.page
    expand = page.get_by_role("button", name="Expand feature search", exact=True)
    rec.focus(expand, zoom=2.4, duration=0.8)
    rec.click(expand, pause=0.4)
    controls = page.locator("foreignObject#gbdraw-feature-search-controls")
    rec.camera(controls.bounding_box(), zoom=2.6, duration=0.6)
    rec.choose(page.locator('select[aria-label="Search field"]'), "qualifier-value" if qualifier else "label",
               in_svg=True)
    if qualifier:
        rec.type(page.locator('input[aria-label="Qualifier key for qualifier value search"]'), qualifier, in_svg=True)
    rec.type(page.locator('input[aria-label="Search features"]'), query, in_svg=True)
    rec.click(page.get_by_role("button", name="Search features", exact=True), pause=0.5)
    status = re.search(r"(\d+) / (\d+) features", page.evaluate(_SEARCH_STATUS))
    if not status or int(status.group(1)) < 1:
        raise AssertionError(f"Searching for {query!r} matched nothing: {page.evaluate(_SEARCH_STATUS)!r}")
    rec.hold(0.8)
    rec.camera(duration=1.0)
    rec.hold(2.2)
    rec.caption(step, "Open a match", "Open zooms to the feature; Original returns to the whole drawing")
    rec.camera(controls.bounding_box(), zoom=2.6, duration=0.8)
    rec.click(page.get_by_role("button", name="Open active feature", exact=True), pause=0.8)
    popup = page.locator("foreignObject#gbdraw-feature-popup")
    expect(popup).to_be_visible()
    rec.camera(duration=1.0)
    rec.hold(1.2)
    _frame_popup(rec, popup)
    rec.hold(2.4)
    if sequences:
        _svg_sequences(rec, popup, step)
        rec.caption(step, "Back to the whole drawing", "Clear the search, then Original")
    rec.click(popup.get_by_role("button", name="x", exact=True), pause=0.3)
    rec.click(page.get_by_role("button", name="Clear search", exact=True), pause=0.3)
    original = page.get_by_role("button", name="Return to original view", exact=True)
    rec.focus(original, zoom=2.4, duration=0.8)
    rec.click(original, pause=0.6)
    rec.camera(duration=1.0)
    rec.hold(1.6)


def _svg_done(rec: Recorder, errors: list[str]) -> None:
    if errors:
        raise AssertionError(f"The Interactive SVG logged errors: {errors}")
    rec.hold(0.8)


# -- Human mitochondrial genome (AT skew) ------------------------------------

HMMT = _fixture("human-mitochondrion-genbank")
AT_SKEW = ("#deaf6e", "#7294e3")
MT_GENES = {"ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6", "COX1", "COX2", "COX3", "ATP6", "ATP8", "CYTB"}


def _hmmtdna_atskew(rec: Recorder, downloads: Path) -> Path:
    _verify("human-mitochondrion-genbank")
    page = rec.page
    preview = page.get_by_role("region", name="Result Preview", exact=True)

    rec.camera(duration=0)
    rec.caption(1, "Switch to Circular mode", "Open gbdraw.app and choose the circular layout")
    rec.hold(1.6)
    _select_mode(rec, "Circular")

    rec.caption(2, "Upload the GenBank file", "HmmtDNA.gbk, the human mitochondrial genome NC_012920.1")
    upload = page.get_by_role("button", name="Choose GenBank/DDBJ File", exact=True)
    rec.focus(upload, zoom=2.1)
    rec.hold(0.4)
    rec.upload(upload, HMMT)
    expect(page.get_by_role("group", name="GenBank/DDBJ File selection", exact=True)).to_contain_text(HMMT.name)
    rec.hold(0.6)
    rec.caption(2, "Name the output and the species", "Type the <i> tags to italicize the center label")
    _field(rec, page.get_by_label("Output Prefix", exact=True), "HmmtDNA_ATskew")
    _field(rec, page.get_by_label("Species", exact=True), "<i>Homo sapiens</i>")
    rec.caption(2, "Generate a first map", "The first run starts Python inside the browser")
    _generate(rec, target=1.8)
    rec.camera(preview.bounding_box(), zoom=1.12, duration=1.4)
    rec.hold(1.6)

    rec.caption(3, "Set the basic circular layout", "Layout › Track Preset: Middle")
    rec.open_section(page.get_by_label("Layout", exact=True))
    _pick(rec, page.get_by_label("Track Preset", exact=True), "middle")
    rec.caption(3, "Put both strands on one track", "Layout › Separate Strands: off, as in the Gallery session")
    strands = page.get_by_label("Separate Strands", exact=True)
    rec.reveal(strands)
    rec.focus(strands, zoom=2.2, duration=0.6)
    rec.set_checked(strands, False)
    # The multi-record canvas drops custom-slot legend entries such as AT skew.
    rec.caption(3, "Draw the record on its own canvas", "Layout › Multi-Record Canvas: off, as in the Gallery session")
    single = page.get_by_label("Multi-Record Canvas", exact=True)
    rec.reveal(single)
    rec.focus(single, zoom=2.2, duration=0.6)
    rec.set_checked(single, False)
    rec.caption(3, "Show labels outside the ring", "Labels › Label Mode: Out · Label Font Size: 20")
    rec.open_section(page.get_by_label("Labels", exact=True))
    _pick(rec, page.locator("#circular-label-mode"), "out")
    _field(rec, page.get_by_label("Label Font Size", exact=True), "20")
    rec.caption(3, "Place the legend and size the center text", "Legend position: Left · Default font size: 28")
    rec.open_section(page.get_by_label("Legend settings", exact=True))
    _pick(rec, page.locator("#legend-position"), "left")
    rec.open_section(page.get_by_label("Titles and Record Labels", exact=True))
    _field(rec, page.get_by_label("Default font size", exact=True), "28")
    rec.caption(3, "Choose the palette", "Colors › Palette: ajisai")
    rec.open_section(page.get_by_label("Colors", exact=True))
    _pick(rec, page.get_by_label("Palette", exact=True), "ajisai")
    rec.caption(3, "Draw the Gallery feature types", "Features › Feature Types: remove misc_RNA")
    _feature_types(rec, remove=("misc_RNA",))
    rec.caption(3, "Generate to check the layout", "")
    _generate(rec)

    rec.caption(4, "Set the dinucleotide window", "Dinucleotide content/skew › Window 500 · Step 50")
    dinucleotide = page.get_by_label("Dinucleotide content/skew", exact=True)
    rec.open_section(dinucleotide)
    windows = dinucleotide.locator("xpath=..").locator("div.grid-cols-3").first.locator("input")
    _field(rec, windows.nth(0), "500")
    _field(rec, windows.nth(1), "50")
    expect(windows.nth(2)).to_have_value("GC")

    rec.caption(5, "Open Custom Track Slots", "Use custom stack, then Reset to preset › Middle")
    _open_custom_slots(rec)
    _press(rec, page.get_by_title("Replace the current custom stack with the Middle template", exact=True))
    page.wait_for_function(
        "() => (window.__GBDRAW_APP__?.adv?.circular_track_slots || []).map((s) => s.id).join()"
        " === 'features,ticks,gc_content,gc_skew'")
    rec.caption(5, "Set the GC track widths", "gc_content and gc_skew › Width: 0.1")
    for slot_id in ("gc_content", "gc_skew"):
        _measure(rec, _slot(rec, slot_id), "Width", "0.1")

    rec.caption(6, "Add the AT skew track", "New track: Dinucleotide skew › Add track")
    _pick(rec, page.get_by_label("New circular track renderer", exact=True), "dinucleotide_skew")
    rec.click(page.get_by_role("button", name=re.compile(r"Add track$")), pause=0.6)
    # Typing the slot id renames the row, so address it by position.
    added = page.locator('[role="group"][aria-label^="Circular track slot "]').nth(4)
    expect(added).to_have_attribute("aria-label", "Circular track slot gc_skew_2")
    rec.caption(6, "Use it for AT skew", "Slot id a_skew_2 · Dinucleotide AT · Legend label AT skew")
    _field(rec, added.get_by_placeholder("slot_id"), "a_skew_2", zoom=2.0)
    _field(rec, _slot_control(added, "Track dinucleotide"), "AT", zoom=2.0)
    _field(rec, _slot_control(added, "Track legend label"), "AT skew", zoom=2.0)
    _measure(rec, added, "Width", "0.1")
    rec.caption(6, "Pick the skew colors", "Positive #deaf6e · Negative #7294e3")
    for title, color in zip(("Positive skew color", "Negative skew color"), AT_SKEW):
        swatch = added.get_by_title(title, exact=True)
        rec.reveal(swatch)
        rec.focus(swatch, zoom=2.2, duration=0.6)
        rec.pick_color(swatch, color)

    rec.caption(7, "Move the ticks to the inside", "Ticks › Move inward until it follows a_skew_2")
    ticks = _slot(rec, "ticks")
    while _slot_ids(rec)[-1] != "ticks":
        _press(rec, ticks.get_by_title("Move inward", exact=True), zoom=2.0, pause=0.6)
    rec.caption(7, "Labels inside, ticks outside", "Ticks › Tick label layout")
    _pick(rec, ticks.locator('select:has(option[value="label_in_tick_out"])'), "label_in_tick_out", zoom=2.0)
    if _slot_ids(rec) != ["features", "gc_content", "gc_skew", "a_skew_2", "ticks"]:
        raise AssertionError(f"Unexpected AT skew slot stack: {_slot_ids(rec)}")

    rec.caption(8, "Prefer gene labels for CDS", "Labels › Qualifier priority: CDS › gene")
    rec.open_section(page.get_by_label("Labels", exact=True))
    _priority_rule(rec, "CDS", "gene")

    rec.caption(9, "Generate the diagram", "GC content, GC skew, and AT skew rings inside the axis")
    _generate(rec, look=1.2)
    rec.mark("atskew")
    _look(rec, rec.text_bounds(r"^AT skew \([+-]\)$"), zoom=2.0, seconds=1.8, pad=60)
    rec.camera(preview.bounding_box(), zoom=1.3, duration=1.2)
    rec.hold(2.0)
    rec.camera(duration=1.0)
    rec.hold(0.6)

    rec.caption(10, "Inspect a feature", "Click RNR2 (l-rRNA), then open Details")
    point = _feature_point(rec, "(item) => item.gene === arg", "RNR2")
    _open_feature(rec, point, "Details")
    rec.hold(3.0)
    _close_feature(rec)

    cytb = _svg_id(rec, "(item) => item.gene === arg", "CYTB")
    exported = export_interactive_svg(rec, downloads, 11)
    errors = _open_export(rec, exported, 11)
    rec.caption(11, "Click a feature for its details", "CYTB: record, type, location, and protein ID")
    _svg_close(rec, _svg_feature(rec, cytb, 11))
    rec.caption(11, "Search the features", "Search › Label: COX finds COX1, COX2, and COX3")
    _svg_search(rec, "COX", 11)
    _svg_done(rec, errors)
    return exported


def _inspect_hmmtdna_atskew(svg: Path) -> dict:
    labels, fills = _labels(svg), _fills(svg)
    report = {
        "gene_symbols": sorted(MT_GENES & labels),
        "product_labels_present": any(label.startswith("NADH dehydrogenase subunit") for label in labels),
        "skew_legend": sorted(label for label in labels if re.fullmatch(r"(AT|GC) skew \([+-]\)", label)),
        "at_skew_fills": sorted(set(AT_SKEW) & fills),
        "species": "Homo sapiens" in labels,
    }
    return _require(report, set(report["gene_symbols"]) == MT_GENES and not report["product_labels_present"]
                    and {"AT skew (+)", "AT skew (-)", "GC skew (+)", "GC skew (-)"} <= set(report["skew_legend"])
                    and len(report["at_skew_fills"]) == 2 and report["species"], svg)


# -- Nicotiana tabacum chloroplast genome regions ----------------------------

PLASTOME = _fixture("tobacco-plastome-genbank")
PLASTOME_COLORS = _fixture("tobacco-plastome-specific-colors")
PLASTOME_PRIORITY = _fixture("tobacco-plastome-qualifier-priority")
PLASTOME_REGIONS = _fixture("tobacco-plastome-regions-table")
REGIONS = ("LSC", "IRb", "SSC", "IRa")


def _tobacco_chloroplast(rec: Recorder, downloads: Path) -> Path:
    _verify("tobacco-plastome-genbank", "tobacco-plastome-specific-colors",
            "tobacco-plastome-qualifier-priority", "tobacco-plastome-regions-table")
    page = rec.page
    preview = page.get_by_role("region", name="Result Preview", exact=True)

    rec.camera(duration=0)
    rec.caption(1, "Switch to Circular mode", "Open gbdraw.app and choose the circular layout")
    rec.hold(1.6)
    _select_mode(rec, "Circular")

    rec.caption(2, "Upload the chloroplast GenBank record", "NC_001879.2, whose coordinates the region table uses")
    upload = page.get_by_role("button", name="Choose GenBank/DDBJ File", exact=True)
    rec.focus(upload, zoom=2.1)
    rec.hold(0.4)
    rec.upload(upload, PLASTOME)
    expect(page.get_by_role("group", name="GenBank/DDBJ File selection", exact=True)).to_contain_text(PLASTOME.name)
    rec.hold(0.6)

    rec.caption(3, "Set the output name and center label", "Type the <i> tags to italicize the species")
    _field(rec, page.get_by_label("Output Prefix", exact=True), "tobacco-chloroplast")
    _field(rec, page.get_by_label("Species", exact=True), "<i>Nicotiana tabacum</i>")
    rec.set_checked(page.get_by_label("Rich Feature Popup", exact=True), True)
    rec.caption(3, "Generate a first map", "The first run starts Python inside the browser")
    _generate(rec, target=1.8)
    rec.camera(preview.bounding_box(), zoom=1.12, duration=1.4)
    rec.hold(1.4)

    rec.caption(4, "Configure the circular layout", "Tuckin preset, separate strands, GC content without GC skew")
    rec.open_section(page.get_by_label("Layout", exact=True))
    _pick(rec, page.get_by_label("Track Preset", exact=True), "tuckin")
    for name, checked in (("Separate Strands", True), ("Hide GC Content", False), ("Hide GC Skew", True)):
        box = page.get_by_label(name, exact=True)
        rec.focus(box, zoom=2.2, duration=0.5)
        rec.set_checked(box, checked, pause=0.4)
    rec.caption(4, "Draw the record on its own canvas", "Layout › Multi-Record Canvas: off, as in the Gallery session")
    single = page.get_by_label("Multi-Record Canvas", exact=True)
    rec.reveal(single)
    rec.focus(single, zoom=2.2, duration=0.6)
    rec.set_checked(single, False)
    rec.caption(4, "Draw the Gallery feature types", "Features › Feature Types: rep_origin instead of repeat_region")
    _feature_types(rec, remove=("repeat_region",), add=("rep_origin",))
    rec.caption(4, "Outline the features", "Block stroke: black, width 1 · Line stroke width: 2")
    _color_value(rec, "Block stroke color", "#000000")
    _field(rec, page.get_by_label("Block Stroke Width", exact=True), "1")
    _field(rec, page.get_by_label("Line Stroke Width", exact=True), "2")
    rec.caption(4, "Thicken the axis", "Axis & Scale › Axis Stroke Width: 3")
    rec.open_section(page.locator("summary").filter(has_text=re.compile(r"^\s*Axis & Scale")))
    _field(rec, page.get_by_label("Axis Stroke Width", exact=True), "3")
    rec.caption(4, "Label both sides of the ring", "Labels › Label Mode: Both (Out + Inner)")
    rec.open_section(page.get_by_label("Labels", exact=True))
    _pick(rec, page.locator("#circular-label-mode"), "both")
    # The offsets are edited before Radial placement, which hides them.
    rec.caption(4, "Set the label offsets", "Outer X / Y 0.9 · Inner X / Y 0.975")
    for label, value in (("Outer X Offset", "0.9"), ("Outer Y Offset", "0.9"),
                         ("Inner X Offset", "0.975"), ("Inner Y Offset", "0.975")):
        _field(rec, page.get_by_text(label, exact=True).locator("xpath=..").locator('input[type="number"]'),
               value)
    rec.caption(4, "Place the labels radially", "Circular Label Placement: Radial")
    _pick(rec, page.get_by_text("Circular Label Placement", exact=True).locator("xpath=..").locator("select"),
          "radial")
    rec.caption(4, "Place the legend and size the text", "Legend position: Upper Left · Default font size: 28")
    rec.open_section(page.get_by_label("Legend settings", exact=True))
    _pick(rec, page.locator("#legend-position"), "upper_left")
    rec.open_section(page.get_by_label("Titles and Record Labels", exact=True))
    _field(rec, page.get_by_label("Default font size", exact=True), "28")

    rec.caption(5, "Load the color table", "Colors › Specific Table (-t): chloroplast_specific_table.tsv")
    rec.open_section(page.get_by_label("Colors", exact=True))
    table = page.get_by_role("button", name="Choose Specific Table (-t)", exact=True)
    rec.reveal(table)
    rec.focus(table, zoom=2.2)
    rec.upload(table, PLASTOME_COLORS)
    rec.hold(0.8)
    rec.caption(5, "Load the label priority table", "Labels › Priority File (TSV): qualifier_priority.tsv")
    rec.open_section(page.get_by_label("Labels", exact=True))
    priority = page.get_by_role("button", name="Choose Priority File (TSV)", exact=True)
    rec.reveal(priority)
    rec.focus(priority, zoom=2.2)
    rec.upload(priority, PLASTOME_PRIORITY)
    rec.hold(0.8)

    rec.caption(6, "Import the four plastome regions", "Region Annotations › Import TSV: nicotiana-tabacum-regions.tsv")
    rec.open_section(page.get_by_label("Region Annotations", exact=True))
    importer = page.get_by_role("button", name="Import TSV", exact=True)
    rec.reveal(importer)
    rec.focus(importer, zoom=2.2)
    rec.upload(importer, PLASTOME_REGIONS)
    expect(page.get_by_label("Annotation set id", exact=True)).to_have_value("plastome_regions")
    expect(page.get_by_label("Annotation lane", exact=True)).to_have_count(4)
    rec.caption(6, "One bracket per region", "LSC, IRb, SSC, and IRa, 1-based and inclusive")
    first = page.get_by_placeholder("Start (1-based)").first
    rec.reveal(first, anchor=0.3)
    rec.focus(page.get_by_label("Annotation set id", exact=True), zoom=1.6)
    rec.move_to(*_center(first.bounding_box()))
    rec.hold(2.4)

    rec.caption(7, "Bind the regions to a track", "Layout › Custom Track Slots › Use custom stack")
    _open_custom_slots(rec)
    rec.caption(7, "Keep only the tracks the map needs", "Remove ticks and gc_skew")
    for slot_id in ("ticks", "gc_skew"):
        if slot_id in _slot_ids(rec):
            _press(rec, _slot(rec, slot_id).get_by_title("Remove", exact=True), zoom=2.0)
    features = _slot(rec, "features")
    rec.caption(7, "Put the features on the axis", "features › Feature on axis")
    _pick(rec, features.locator('select:has(option[value="split"])'), "split", zoom=2.0)
    rec.caption(7, "Add an annotation track", "New track: Annotations › Add track")
    _pick(rec, page.get_by_label("New circular track renderer", exact=True), "annotations")
    rec.click(page.get_by_role("button", name=re.compile(r"Add track$")), pause=0.6)
    annotations = page.locator('[role="group"][aria-label^="Circular track slot "]').last
    expect(annotations).to_have_attribute("aria-label", "Circular track slot annotations")
    rec.caption(7, "Show plastome_regions inside the features", "Annotation set · Placement: inside")
    _pick(rec, _slot_control(annotations, "Annotation set"), "plastome_regions", zoom=2.0)
    _pick(rec, _slot_control(annotations, "Annotation placement"), "inside", zoom=2.0)
    _field(rec, annotations.get_by_placeholder("slot_id"), "plastome_regions", zoom=2.0)
    regions = _slot(rec, "plastome_regions")
    while _slot_ids(rec).index("plastome_regions") > 1:
        _press(rec, regions.get_by_title("Move outward", exact=True), zoom=2.0, pause=0.6)
    rec.caption(7, "Size the bracket lane", "Radius 0.65 ×R · Width 20 px · Labels on · Overflow: Compress")
    _measure(rec, regions, "Radius", "0.65")
    _measure(rec, regions, "Width", "20", "px")
    rec.fast_forward(factor=2.0)
    _field(rec, regions.get_by_title("Inner gap", exact=True), "1", zoom=2.0)
    _field(rec, regions.get_by_title("Outer gap", exact=True), "1", zoom=2.0)
    rec.normal_speed()
    show = _slot_control(regions, "Show annotation labels")
    rec.reveal(show)
    rec.set_checked(show, True)
    _pick(rec, regions.locator('select:has(option[value="compress"])'), "compress", zoom=2.0)
    rec.fast_forward(factor=2.0)
    _field(rec, _slot_control(regions, "Padding"), "1", zoom=2.0)
    rec.normal_speed()
    rec.caption(7, "Move GC content farther inward", "gc_content › Radius 0.56 ×R · Width 0.08 ×R")
    gc = _slot(rec, "gc_content")
    _measure(rec, gc, "Radius", "0.56")
    _measure(rec, gc, "Width", "0.08")
    if _slot_ids(rec) != ["features", "plastome_regions", "gc_content"]:
        raise AssertionError(f"Unexpected chloroplast slot stack: {_slot_ids(rec)}")

    rec.caption(8, "Generate and inspect the chloroplast map", "One labeled bracket lane between the genes and GC content")
    _generate(rec, look=1.2)
    rec.mark("regions")
    _look(rec, rec.text_bounds(r"^(LSC|IRb|SSC|IRa)$"), zoom=1.5, seconds=2.4, pad=40)
    rec.camera(preview.bounding_box(), zoom=1.3, duration=1.2)
    rec.hold(1.8)
    rec.camera(duration=1.0)
    rec.hold(0.6)

    rbcl = _svg_id(rec, "(item) => item.gene === arg", "rbcL")
    exported = export_interactive_svg(rec, downloads, 9)
    errors = _open_export(rec, exported, 9)
    rec.caption(9, "Click a feature for its details", "rbcL: every qualifier from the GenBank record")
    _svg_close(rec, _svg_feature(rec, rbcl, 9, first="Qualifiers"))
    rec.caption(9, "Search the features", "Search › Label: psb highlights the photosystem II genes")
    _svg_search(rec, "psb", 9)
    _svg_done(rec, errors)
    return exported


def _inspect_tobacco_chloroplast(svg: Path) -> dict:
    labels = _labels(svg)
    text = svg.read_text(encoding="utf-8")
    report = {
        "regions": sorted(set(REGIONS) & labels),
        "annotations": len(set(re.findall(r'data-gbdraw-annotation-id="([^"]+)"', text))),
        "genes": sorted({"matK", "psaA", "rbcL"} & labels),
        "species": "Nicotiana tabacum" in labels,
        "gc_content": "GC content" in labels,
        "skew_legend": sorted(label for label in labels if "skew" in label),
    }
    return _require(report, len(report["regions"]) == 4 and report["annotations"] == 4
                    and len(report["genes"]) == 3 and report["species"] and report["gc_content"]
                    and not report["skew_legend"], svg)


# -- Aminoglycoside biosynthetic gene clusters from Streptomyces spp. ---------

BGC_FILES = tuple(f"bgc-{number}-genbank" for number in ("0000708", "0000709", "0000711", "0000712", "0000713"))
BGC_SOURCES = tuple(_fixture(file_id) for file_id in BGC_FILES)
BGC_RECORDS = (
    ("<i>Streptomyces lividus</i> CBS 844.73", "Lividomycin biosynthetic gene cluster"),
    ("<i>Streptomyces fradiae</i> ATCC 10745", "Neomycin biosynthetic gene cluster"),
    ("<i>Streptomyces fradiae</i> MCIMB 8233", "Neomycin biosynthetic gene cluster"),
    ("<i>Streptomyces rimosus</i> subsp. <i>paromomycinus</i> NRRL 2455", "Paromomycin biosynthetic gene cluster"),
    ("<i>Streptomyces ribosidificus</i> ATCC 21294", "Ribostamycin biosynthetic gene"),
)
BGC_TITLE = "Aminoglycoside biosynthetic gene clusters from <i>Streptomyces</i> spp."
BGC_RULES = (
    ("CDS", "gene_kind", "biosynthetic$", "#d03535", "Core biosynthetic genes"),
    ("CDS", "gene_kind", "biosynthetic-additional", "#f787a9", "Additional biosynthetic genes"),
    ("CDS", "gene_kind", "transport", "#577edb", "Transport-related genes"),
    ("CDS", "gene_kind", "regulatory", "#57b767", "Regulatory genes"),
)
BGC_DEFINITION_STYLES = (
    ("Name / Species", True, None),
    ("Subtitle", False, None),
    ("Accession", False, "#7b7c7d"),
    ("Length / Coordinates", False, "#7b7c7d"),
)


def _record_options(rec: Recorder, index: int) -> Locator:
    button = rec.page.get_by_role("button", name=f"Record options for sequence {index}", exact=True)
    definition = rec.page.get_by_label(f"Definition for sequence {index}", exact=True)
    if not definition.is_visible():
        _press(rec, button, zoom=2.0, pause=0.4)
    return definition


def _align_to_group(rec: Recorder, popup: Locator, group_id: str) -> None:
    """Align on the popup's similarity group through Review alignment options…."""

    page = rec.page
    review = popup.get_by_role("button", name=re.compile(r"^Review alignment options"))
    _scroll_popup_to(rec, popup, "Review alignment options…")
    section = popup.locator('[data-feature-popup-group="similarity"]')
    box = section.bounding_box()
    rec.camera(box, zoom=1.9, duration=0.9, portrait={"x": box["x"] - 20, "y": box["y"] - 20,
                                                      "width": box["width"] + 40, "height": box["height"] + 40})
    for name in ("Show in editor", "Align…"):
        rec.focus(popup.get_by_role("button", name=name, exact=True), zoom=2.2, duration=0.6)
        rec.hold(0.5)
    rec.caption(8, "Review the alignment options", "Review alignment options… lists every candidate per record")
    rec.focus(review, zoom=2.4, duration=0.6)
    rec.click(review, pause=0.4)
    dialog = page.get_by_role("dialog", name="Review alignment", exact=True)
    expect(dialog).to_be_visible(timeout=ACTION_TIMEOUT_MS)
    page.wait_for_timeout(600)
    rec.camera(dialog.bounding_box(), zoom=1.4, duration=0.8)
    rec.hold(1.6)
    # BGC0000708 has two og_6 members: the recommended livU and the coortholog CAG38692.1.
    card = dialog.locator("[data-alignment-record-key]").filter(has_text="CAG38692.1")
    expect(card).to_have_count(1)
    candidates = card.get_by_role("radio", name=re.compile(r"^Select "))
    if candidates.count() != 2:
        raise AssertionError(f"Expected two og_6 candidates for BGC0000708, found {candidates.count()}")
    recommended = card.locator("label").filter(has_text="Recommended")
    expect(recommended).to_have_count(1)
    expect(recommended).to_contain_text("livU")
    expect(card).to_contain_text("CAG38692.1")
    rec.camera(card.bounding_box(), zoom=1.8, duration=0.8)
    other = card.locator("label").filter(has_not_text="Recommended")
    spot = other.first.bounding_box()
    rec.move_to(spot["x"] + spot["width"] / 2, spot["y"] + spot["height"] / 2)
    rec.hold(1.4)
    rec.click(recommended.get_by_role("radio"), pause=0.5)
    rec.hold(0.8)
    rec.caption(8, "Keep the arrow direction", "Arrow direction: Keep as is · All right · All left")
    direction = dialog.locator("[data-similarity-alignment-direction-mode]")
    rec.camera(direction.bounding_box(), zoom=2.0, duration=0.8)
    expect(direction.get_by_role("radio", name="Keep as is")).to_be_checked()
    rec.hold(2.0)
    rec.camera(dialog.bounding_box(), zoom=1.4, duration=0.8)
    rec.click(dialog.get_by_role("button", name="Apply", exact=True), pause=0.2)
    expect(dialog).to_be_hidden(timeout=ACTION_TIMEOUT_MS)
    rec.camera(duration=0.8)
    rec.fast_forward(target=1.4)
    # Apply regenerates the preview through the alignment owner, which leaves the run marker as it was.
    page.wait_for_function(
        """(group) => { const app = window.__GBDRAW_APP__;
          return !app?.processing && app?.similarityAlignmentStatus === 'idle' && !app?.errorLog
            && app?.similarityAlignmentPlanInspector?.groupId === group; }""",
        arg=group_id, timeout=GENERATION_TIMEOUT_MS)
    rec.hold(0.3)
    rec.normal_speed()


def _bgc_aminoglycosides(rec: Recorder, downloads: Path) -> Path:
    _verify(*BGC_FILES)
    page = rec.page

    rec.camera(duration=0)
    rec.caption(1, "Switch to Linear mode", "Five clusters, one row each")
    rec.hold(1.6)
    _select_mode(rec, "Linear")

    rec.caption(2, "Upload the five GenBank files", "BGC0000708, 709, 711, 712, and 713 from MIBiG, in this order")
    add = page.get_by_role("button", name="Add sequence", exact=True)
    rec.focus(add, zoom=2.0)
    for _ in range(len(BGC_SOURCES) - 1):
        _press(rec, add, zoom=2.0, pause=0.35)
    for index, source in enumerate(BGC_SOURCES, start=1):
        if index == 2:
            rec.fast_forward(factor=3.0)
        row = page.get_by_role("region", name=f"Linear input file {index}", exact=True)
        chooser = row.get_by_role("button", name=re.compile(r"^Choose GenBank"))
        rec.reveal(chooser)
        rec.focus(chooser, zoom=2.0, duration=0.6)
        rec.upload(chooser, source)
        expect(row).to_contain_text(source.name)
    rec.normal_speed()
    rec.caption(2, "Fill in each record's text", "Record options › Organism / strain and Subtitle / title")
    for index, (definition, subtitle) in enumerate(BGC_RECORDS, start=1):
        if index == 2:
            rec.fast_forward(factor=3.0)
        field = _record_options(rec, index)
        _field(rec, field, definition, zoom=2.0)
        _field(rec, page.get_by_label(f"Subtitle / title for sequence {index}", exact=True), subtitle, zoom=2.0)
    rec.normal_speed()

    rec.caption(3, "Run LOSAT in the browser", "Comparison › Run LOSAT for all adjacent pairs")
    _run_losat(rec)
    settings = page.locator('details[data-linear-comparison-disclosure="settings"]')
    summary = settings.get_by_role("button", name="Comparison Settings", exact=True)
    if settings.get_attribute("open") is None:
        _press(rec, summary, zoom=2.0)
    rec.caption(3, "Group proteins by similarity", "LOSAT Mode: LOSATP · LOSATP mode: Similarity groups")
    losatp = settings.get_by_role("group", name="LOSAT Mode", exact=True).get_by_role("button", name="LOSATP",
                                                                                      exact=True)
    _press(rec, losatp, zoom=2.2)
    expect(losatp).to_have_attribute("aria-pressed", "true")
    _pick(rec, settings.get_by_role("combobox", name="LOSATP mode", exact=True), "orthogroup")
    rec.caption(3, "Set the match filters", "Bitscore 50 · E-value 0.01 · Identity 30 · Length 0")
    for label, value in (("minimum bitscore", "50"), ("maximum e-value", "0.01"),
                         ("minimum identity", "30"), ("minimum alignment length", "0")):
        _value(rec, settings.get_by_label(f"Linear comparison {label}", exact=True), value)
    rec.caption(3, "Let the browser choose the runtime", "Execution: Auto · Total threads: Safe")
    for name, value in (("LOSAT execution", "auto"), ("LOSAT total threads", "safe")):
        _pick(rec, settings.get_by_role("combobox", name=name, exact=True), value, zoom=2.0)

    rec.caption(4, "Reverse-complement BGC0000713", "Record options for sequence 5 › Reverse complement")
    _record_options(rec, 5)
    reverse = page.get_by_label("Reverse complement for sequence 5", exact=True)
    rec.reveal(reverse)
    rec.focus(reverse, zoom=2.2, duration=0.6)
    rec.set_checked(reverse, True)
    rec.caption(4, "Keep features in the middle of each axis",
                "Layout › Track Layout: Middle · Lock Definition Column")
    rec.open_section(page.get_by_label("Layout", exact=True))
    _pick(rec, page.get_by_label("Track Layout", exact=True), "middle")
    lock = page.get_by_label("Lock Definition Column", exact=True)
    rec.reveal(lock)
    rec.focus(lock, zoom=2.2, duration=0.6)
    rec.set_checked(lock, True)
    rec.caption(4, "Put both strands on one track", "Layout › Separate Strands: off, as in the Gallery session")
    strands = page.get_by_label("Separate Strands", exact=True)
    rec.reveal(strands)
    rec.focus(strands, zoom=2.2, duration=0.6)
    rec.set_checked(strands, False)
    rec.caption(4, "Narrow the arrow shafts", "Features › Arrow Geometry › Shaft Width Ratio: 0.6")
    rec.open_section(page.locator('summary[aria-label="Features"]'))
    _field(rec, page.get_by_label("Arrow shaft width ratio", exact=True), "0.6", zoom=2.2)

    rec.caption(5, "Label the first record", "Labels › Show Labels: First Record Only · CDS › gene")
    rec.open_section(page.get_by_label("Labels", exact=True))
    _pick(rec, page.get_by_label("Show Labels", exact=True), "first")
    _priority_rule(rec, "CDS", "gene")
    rec.caption(5, "Rotate the labels", "Label Placement: Above Feature · Label Rotation: 45 · Font size: 18")
    _pick(rec, page.get_by_label("Label Placement", exact=True), "above_feature")
    _field(rec, page.get_by_label("Label Rotation", exact=True), "45")
    _field(rec, page.get_by_label("Label Font Size", exact=True), "18")
    rec.caption(5, "Gray out the other CDS features", "Colors › CDS: #d3d3d3 on the default palette")
    rec.open_section(page.get_by_label("Colors", exact=True))
    cds = page.get_by_label("CDS feature color", exact=True)
    rec.reveal(cds)
    rec.focus(cds, zoom=2.2, duration=0.6)
    rec.pick_color(cds, "#d3d3d3")
    rec.caption(5, "Draw the Gallery feature types", "Features › Feature Types: remove misc_RNA")
    _feature_types(rec, remove=("misc_RNA",))
    rec.caption(5, "Size and outline the features", "Height 75 · Block stroke #262626, width 2 · Line stroke width 2")
    _field(rec, page.get_by_label("Feature Height", exact=True), "75")
    _color_value(rec, "Block stroke color", "#262626")
    _field(rec, page.get_by_label("Block Stroke Width", exact=True), "2")
    _field(rec, page.get_by_label("Line Stroke Width", exact=True), "2")
    rec.caption(5, "Add a ruler on a gray axis", "Axis & Scale › Axis gray, width 5 · Scale Style: Ruler · Font size 20")
    rec.open_section(page.get_by_label("Axis & Scale", exact=True))
    show_scale = page.get_by_label("Show Coordinate Scale (Linear)", exact=True)
    rec.reveal(show_scale)
    rec.set_checked(show_scale, True)
    _color_value(rec, "Axis stroke color", "#808080")
    _field(rec, page.get_by_label("Axis Stroke Width", exact=True), "5")
    _pick(rec, page.get_by_label("Linear scale style", exact=True), "ruler")
    _field(rec, page.get_by_label("Linear scale font size", exact=True), "20")
    ruler = page.get_by_label("Linear ruler label font size", exact=True)
    if ruler.input_value() != "20":
        _field(rec, ruler, "20")
    rec.caption(5, "Put the legend at the bottom", "Legend position: Bottom · Swatch size 20 · Font size 20")
    rec.open_section(page.get_by_label("Legend settings", exact=True))
    _pick(rec, page.locator("#legend-position"), "bottom")
    _field(rec, page.get_by_label("Legend swatch size", exact=True), "20")
    _field(rec, page.get_by_label("Legend font size", exact=True), "20")

    rec.caption(6, "Add the plot title", "Titles & Record Labels › Plot Title, Position: Bottom (Default)")
    rec.open_section(page.get_by_label("Titles and Record Labels", exact=True))
    _field(rec, page.get_by_role("textbox", name="Plot Title", exact=True), BGC_TITLE, zoom=2.0)
    _pick(rec, page.get_by_label("Plot Title Position", exact=True), "bottom")
    rec.caption(6, "Show every definition line", "Accession and Length / Coordinates visibility: Show")
    for label in ("Accession", "Length / Coordinates"):
        _pick(rec, page.get_by_label(f"{label} visibility", exact=True), "show", zoom=2.0)
    rec.caption(6, "Style the definition lines", "Size 20 · bold names · gray accession and length")
    rec.fast_forward(factor=3.0)
    for label, bold, color in BGC_DEFINITION_STYLES:
        _press(rec, page.get_by_label(f"{label} style", exact=True), zoom=2.0, pause=0.3)
        _field(rec, page.get_by_label(f"{label} style size", exact=True), "20", zoom=2.0)
        weight = page.get_by_role("button", name=f"{label} {'bold' if bold else 'normal'} weight", exact=True)
        _press(rec, weight, zoom=2.0, pause=0.3)
        if color:
            _field(rec, page.get_by_label(f"{label} definition line fill value", exact=True), color, zoom=2.0)
    rec.normal_speed()

    rec.caption(7, "Color genes by antiSMASH category", "Colors › Specific rules on the gene_kind qualifier")
    rec.open_section(page.get_by_label("Colors", exact=True))
    _color_rules(rec, BGC_RULES)

    rec.caption(8, "Generate the comparison", "LOSATP searches every pair of records in the browser")
    _generate(rec, target=2.4, look=2.4)
    groups = page.evaluate("() => (window.__GBDRAW_APP__?.orthogroups || []).length")
    if groups != 23:
        raise AssertionError(f"Expected 23 similarity groups, found {groups}")
    rec.caption(8, "Align the rows on og_6", "Click the neoU CDS of BGC0000709, then Align…")
    neou = _feature_point(rec, "(item) => item.gene === arg && item.orthogroupId === 'og_6' && /^BGC0000709/.test(item.record_id || '')",
                          "neoU")
    rec.mark("align")
    rec.click_point(neou["x"], neou["y"], pause=0.4)
    popup = page.get_by_role("dialog", name=re.compile(r"^Feature details"))
    expect(popup).to_be_visible()
    expect(popup).to_contain_text("og_6")
    box = popup.bounding_box()
    rec.camera({"x": min(box["x"], neou["x"]) - 40, "y": box["y"] - 40,
                "width": max(box["x"] + box["width"], neou["x"]) - min(box["x"], neou["x"]) + 80,
                "height": box["height"] + 80}, zoom=1.5, duration=0.9,
               portrait={"x": box["x"] - 30, "y": box["y"] - 30, "width": box["width"] + 60,
                         "height": box["height"] + 60})
    rec.hold(1.2)
    _align_to_group(rec, popup, "og_6")
    _fit_preview(rec)
    rec.camera(duration=1.0)
    rec.hold(2.4)

    rec.caption(9, "Inspect a similarity-group ribbon", "og_18: group ID, 6 genes in 5 records, and the member table")
    ribbon = _hit(rec, '[data-match-kind="orthogroup"][data-orthogroup-id="og_18"]')
    rec.click_point(ribbon["x"], ribbon["y"], pause=0.4)
    details = page.get_by_role("dialog", name="Pairwise match details", exact=True)
    expect(details).to_be_visible()
    expect(details).to_contain_text("og_18")
    expect(details).to_contain_text("6 genes in 5 records")
    rec.camera(details.bounding_box(), zoom=1.5, duration=1.0)
    rec.hold(3.2)
    rec.click(details.get_by_label("Close match popup", exact=True), pause=0.3)
    rec.camera(duration=0.9)
    rec.hold(0.5)

    rec.caption(10, "Inspect the feature qualifiers", "livE: gene supplies the label, gene_kind the color")
    live = _feature_point(rec, "(item) => item.gene === arg && item.orthogroupId === 'og_1'", "livE")
    popup = _open_feature(rec, live, "Qualifiers")
    rec.hold(1.0)
    _scroll_popup_to(rec, popup, "gene_kind")
    rec.hold(3.0)
    _close_feature(rec)

    live_id = _svg_id(rec, "(item) => item.gene === arg && item.orthogroupId === 'og_1'", "livE")
    exported = export_interactive_svg(rec, downloads, 11)
    errors = _open_export(rec, exported, 11)
    rec.caption(11, "Click a feature for its similarity group", "livE: the other og_1 members light up in every cluster")
    _svg_close(rec, _svg_feature(rec, live_id, 11, overview=2.6, shows="og_1"))
    rec.caption(11, "Click a similarity-group link", "og_1: every member with its record, coordinates, and product")
    popup = _svg_popup(rec, _hit(rec, '[data-match-kind="orthogroup"][data-orthogroup-id="og_1"]'))
    expect(popup).to_contain_text("og_1")
    rec.hold(3.0)
    rec.caption(11, "Download the group's sequences", "DL aa saves every member protein in one FASTA file")
    save = popup.get_by_role("button", name=re.compile(r"^DL aa \((\d+)\)$"))
    members = int(re.search(r"\((\d+)\)", save.inner_text()).group(1))
    rec.focus(save, zoom=2.4, duration=0.8)
    with rec.page.expect_download(timeout=ACTION_TIMEOUT_MS) as info:
        rec.click(save, pause=0.2)
    fasta = downloads / info.value.suggested_filename
    info.value.save_as(fasta)
    if fasta.read_text(encoding="utf-8").count(">") != members:
        raise AssertionError(f"{fasta.name} does not hold the {members} og_1 proteins")
    rec.toast(f"{fasta.name} downloaded ({max(1, round(fasta.stat().st_size / 1024))} KB)", "download")
    rec.hold(2.4)
    _svg_close(rec, popup)
    _svg_done(rec, errors)
    return exported


def _inspect_bgc_aminoglycosides(svg: Path) -> dict:
    labels, fills = _labels(svg), _fills(svg)
    text = svg.read_text(encoding="utf-8")
    report = {
        "rule_fills": sorted({color for *_, color, _ in BGC_RULES} & fills),
        "legend": sorted({legend for *_, legend in BGC_RULES} & labels),
        "similarity_links": text.count('data-match-kind="orthogroup"'),
        "livE": "livE" in labels,
        "records": sorted({"BGC0000708", "BGC0000709", "BGC0000711", "BGC0000712", "BGC0000713"}
                          & {label.split(".")[0] for label in labels}),
    }
    return _require(report, len(report["rule_fills"]) == 4 and len(report["legend"]) == 4
                    and report["similarity_links"] > 0 and report["livE"] and len(report["records"]) == 5, svg)


# -- Vibrio nigripulchritudo chromosomes and plasmids on one canvas -----------

# The GBFF is not in the tutorial data package. The Gallery session was built
# from this RefSeq file, so its size and digest are pinned here instead.
VNIG = REPO_ROOT / "tests/test_inputs/GCF_015097735.1_ASM1509773v1_genomic.gbff"
VNIG_IDENTITY = (17494876, "a2f696d584451866032339eb17de4eadad55372a36bc87b9eb85583a482ad027")
VNIG_TITLE = "<i>Vibrio nigripulchritudo</i> TUMSAT-TG-2018, complete genome"
VNIG_RECORDS = (("Chromosome 1", "NZ_AP024087.1"), ("Chromosome 2", "NZ_AP024088.1"), ("pVNTG1", "NZ_AP024089.1"),
                ("pVNTG2", "NZ_AP024090.1"), ("pVNTG3", "NZ_AP024091.1"), ("pVNTG4", "NZ_AP024092.1"))
VNIG_ROWS = ("1", "1", "2", "2", "2", "2")
_RECORD_ORDER = """
() => (window.__GBDRAW_APP__?.adv?.multi_record_positions || []).map((entry) => `${entry.selector}@${entry.row}`).join()
"""


def _vnig_multirecord(rec: Recorder, downloads: Path) -> Path:
    assert_fixture_identity(VNIG, expected_size=VNIG_IDENTITY[0], expected_sha256=VNIG_IDENTITY[1])
    page = rec.page
    preview = page.get_by_role("region", name="Result Preview", exact=True)

    rec.camera(duration=0)
    rec.caption(1, "Switch to Circular mode", "Open gbdraw.app and choose the circular layout")
    rec.hold(1.6)
    _select_mode(rec, "Circular")

    rec.caption(2, "Upload the GenBank file", "One RefSeq GBFF holds two chromosomes and four plasmids")
    upload = page.get_by_role("button", name="Choose GenBank/DDBJ File", exact=True)
    rec.focus(upload, zoom=2.1)
    rec.hold(0.4)
    rec.upload(upload, VNIG)
    expect(page.get_by_role("group", name="GenBank/DDBJ File selection", exact=True)).to_contain_text(VNIG.name)
    rec.hold(0.6)
    rec.caption(2, "Generate a first map", "The first run starts Python inside the browser")
    _generate(rec, target=1.8)
    rec.camera(preview.bounding_box(), zoom=1.12, duration=1.4)
    rec.hold(1.6)

    rec.caption(3, "Set the circular layout", "Output Prefix · Track Preset: Tuckin · Separate Strands: on")
    _field(rec, page.get_by_label("Output Prefix", exact=True), "Vnig_TUMSAT-TG-2018")
    rec.open_section(page.get_by_label("Layout", exact=True))
    _pick(rec, page.get_by_label("Track Preset", exact=True), "tuckin")
    strands = page.get_by_label("Separate Strands", exact=True)
    rec.reveal(strands)
    rec.focus(strands, zoom=2.2, duration=0.6)
    rec.set_checked(strands, True)
    rec.caption(3, "Choose the palette", "Colors › Palette: orchid")
    rec.open_section(page.get_by_label("Colors", exact=True))
    _pick(rec, page.get_by_label("Palette", exact=True), "orchid")
    rec.caption(3, "Put the legend on the left", "Legend position: Left")
    rec.open_section(page.get_by_label("Legend settings", exact=True))
    _pick(rec, page.locator("#legend-position"), "left")
    rec.caption(3, "Draw the Gallery feature types", "Features › Feature Types: remove misc_RNA")
    _feature_types(rec, remove=("misc_RNA",))

    rec.caption(4, "Enable Multi-Record Canvas", "Layout › Multi-Record Canvas: on · Record Size Mode: Auto")
    rec.open_section(page.get_by_label("Layout", exact=True))
    canvas = page.get_by_label("Multi-Record Canvas", exact=True)
    rec.reveal(canvas)
    rec.focus(canvas, zoom=2.2, duration=0.6)
    rec.set_checked(canvas, True)
    _pick(rec, page.get_by_label("Record Size Mode", exact=True), "auto")
    rec.caption(4, "Set the record sizes and gaps", "Min Radius Ratio 0.55 · Column Gap Ratio 0.1 · Row Gap Ratio 0.05")
    for label, value in (("Min Radius Ratio", "0.55"), ("Column Gap Ratio", "0.1"), ("Row Gap Ratio", "0.05")):
        _value(rec, page.get_by_label(label, exact=True), value)
    rec.caption(4, "Chromosomes on top, plasmids below", "Record Order: #1 and #2 in row 1, #3 to #6 in row 2")
    rows = page.locator('[data-capture="multi-record-order"] select')
    expect(rows).to_have_count(len(VNIG_ROWS))
    changed = 0
    for index, row in enumerate(VNIG_ROWS):
        select = rows.nth(index)
        if select.input_value() == row:
            continue
        changed += 1
        if changed == 2:
            rec.fast_forward(factor=3.0)
        _pick(rec, select, row, zoom=2.0)
    rec.normal_speed()
    expected = ",".join(f"#{index}@{row}" for index, row in enumerate(VNIG_ROWS, start=1))
    if page.evaluate(_RECORD_ORDER) != expected:
        raise AssertionError(f"Unexpected record order: {page.evaluate(_RECORD_ORDER)}")

    rec.caption(5, "Add the plot title", "Titles & Record Labels › Plot Title, Position: Bottom")
    rec.open_section(page.get_by_label("Titles and Record Labels", exact=True))
    _field(rec, page.get_by_role("textbox", name="Plot Title", exact=True), VNIG_TITLE, zoom=2.0)
    _pick(rec, page.get_by_label("Plot Title Position", exact=True), "bottom")

    rec.caption(6, "Generate the diagram", "Two chromosomes above four plasmids, each with GC content and GC skew")
    _generate(rec, look=1.2)
    rec.mark("canvas")
    _look(rec, rec.text_bounds(r"^pVNTG[1-4]$"), zoom=1.5, seconds=2.2, pad=150)
    _look(rec, rec.text_bounds(r"TUMSAT-TG-2018, complete genome$"), zoom=2.0, seconds=1.8, pad=30)
    rec.camera(duration=1.0)
    rec.hold(0.6)

    rec.caption(7, "Inspect a CDS feature", "Click topB on pVNTG3, then open Details")
    topb = _feature_point(rec, "(item) => item.gene === arg && item.record_id === 'NZ_AP024091.1'", "topB")
    popup = _open_feature(rec, topb, "Details")
    rec.hold(3.0)
    rec.caption(7, "Its qualifiers and sequence", "The Qualifiers and Sequence tabs of the same popup")
    for tab in ("Qualifiers", "Sequence"):
        _feature_tab(rec, popup, tab, pause=0.4)
        rec.hold(2.6)
    _close_feature(rec)

    exported = export_interactive_svg(rec, downloads, 8)
    errors = _open_export(rec, exported, 8)
    rec.caption(8, "Find a gene on a chromosome", "Search › Label: dnaA, the replication initiator on chromosome 1")
    _svg_search(rec, "dnaA", 8, sequences=True)
    _svg_done(rec, errors)
    return exported


VNIG_GALLERY = REPO_ROOT / "gbdraw/web/gallery/sources/Vnig_TUMSAT-TG-2018.svg"


def _record_offsets(text: str) -> dict[str, tuple[float, float]]:
    """Where each record sits on the multi-record canvas, keyed by record ID."""

    offsets = {}
    for tag in re.findall(r'<g [^>]*\bid="record_\d+"[^>]*>', text):
        record_id = re.search(r'data-gbdraw-record-id="([^"]+)"', tag).group(1)
        pairs = re.findall(r"translate\(([-\d.]+),([-\d.]+)\)", re.search(r'transform="([^"]+)"', tag).group(1))
        offsets[record_id] = (round(sum(float(x) for x, _ in pairs), 1), round(sum(float(y) for _, y in pairs), 1))
    return offsets


def _inspect_vnig_multirecord(svg: Path) -> dict:
    labels, fills = _labels(svg), _fills(svg)
    text = svg.read_text(encoding="utf-8")
    offsets = _record_offsets(text)
    tops = {name: offsets[accession][1] for name, accession in VNIG_RECORDS if accession in offsets}
    report = {
        "records": sorted(name for name, accession in VNIG_RECORDS if {name, accession} <= labels),
        # Each row centers its records vertically, so their tops differ.
        "chromosomes_above_plasmids": len(tops) == 6 and max(tops[name] for name, _ in VNIG_RECORDS[:2])
                                      < min(tops[name] for name, _ in VNIG_RECORDS[2:]),
        "gallery_placement": offsets == _record_offsets(VNIG_GALLERY.read_text(encoding="utf-8")),
        "title": {"Vibrio nigripulchritudo", " TUMSAT-TG-2018, complete genome"} <= labels,
        "gc_tracks": {"GC content", "GC skew (+)", "GC skew (-)"} <= labels,
        "orchid_cds": "#c25792" in fills,
        "legend_left": "&quot;legendSide&quot;:&quot;left&quot;" in text[:20000],
    }
    return _require(report, len(report["records"]) == 6
                    and all(value for key, value in report.items() if key != "records"), svg)


# -- Hepatoplasmataceae collinear protein-match blocks ------------------------

HEP_FILES = ("hepatoplasmataceae-ap027078-genbank", "hepatoplasmataceae-ap027131-genbank", "ap027133-genbank",
             "hepatoplasmataceae-ap027132-genbank", "hepatoplasmataceae-nz-cp006932-genbank")
HEP_SOURCES = tuple(_fixture(file_id) for file_id in HEP_FILES)
HEP_RECORDS = ("AP027078.1", "AP027131.1", "AP027133.1", "AP027132.1", "NZ_CP006932.1")
HEP_GALLERY = REPO_ROOT / "gbdraw/web/gallery/sources/hepatoplasmataceae_collinear.svg"
# The tutorial's example block: AP027078.1 to AP027131.1, plus orientation.
HEP_BLOCK = "block_0024"


def _blocks(text: str) -> set[str]:
    return set(re.findall(r'data-collinearity-block-id="([^"]+)"', text))


def _hepatoplasmataceae_collinear(rec: Recorder, downloads: Path) -> Path:
    _verify(*HEP_FILES)
    page = rec.page

    rec.camera(duration=0)
    rec.caption(1, "Switch to Linear mode", "Five genomes, one row each")
    rec.hold(1.6)
    _select_mode(rec, "Linear")

    rec.caption(2, "Upload the five GenBank files", "AP027078, AP027131, AP027133, AP027132, and NZ_CP006932, in this order")
    add = page.get_by_role("button", name="Add sequence", exact=True)
    rec.focus(add, zoom=2.0)
    for _ in range(len(HEP_SOURCES) - 1):
        _press(rec, add, zoom=2.0, pause=0.35)
    for index, source in enumerate(HEP_SOURCES, start=1):
        if index == 2:
            rec.fast_forward(factor=3.0)
        row = page.get_by_role("region", name=f"Linear input file {index}", exact=True)
        chooser = row.get_by_role("button", name=re.compile(r"^Choose GenBank"))
        rec.reveal(chooser)
        rec.focus(chooser, zoom=2.0, duration=0.6)
        rec.upload(chooser, source)
        expect(row).to_contain_text(source.name)
    rec.normal_speed()

    rec.caption(3, "Run LOSAT in the browser", "Comparison › Run LOSAT for all adjacent pairs")
    _run_losat(rec)
    settings = page.locator('details[data-linear-comparison-disclosure="settings"]')
    if settings.get_attribute("open") is None:
        _press(rec, settings.get_by_role("button", name="Comparison Settings", exact=True), zoom=2.0)
    losatp = settings.get_by_role("group", name="LOSAT Mode", exact=True).get_by_role("button", name="LOSATP",
                                                                                      exact=True)
    rec.caption(3, "Let the browser choose the runtime",
                "Comparison Settings › Runtime and reproducibility: Execution Auto · Total threads Safe")
    for name, value in (("LOSAT execution", "auto"), ("LOSAT total threads", "safe")):
        _pick(rec, page.get_by_role("combobox", name=name, exact=True), value, zoom=2.0)
    # Both already read Auto, and Threads per run stays locked until LOSATP.
    if page.evaluate("() => { const losat = window.__GBDRAW_APP__?.losatExecution;"
                     " return losat?.parallelWorkers == null && losat?.threadsPerJob === 'auto'; }") is not True:
        raise AssertionError("Parallel runs and Threads per run are not Auto")
    for name in ("LOSAT parallel runs", "LOSAT threads per run"):
        _point(rec, page.get_by_role("combobox", name=name, exact=True), zoom=2.0)

    rec.caption(4, "Build collinear blocks", "LOSAT Mode: LOSATP · LOSATP mode: Collinear blocks")
    _press(rec, losatp, zoom=2.2)
    expect(losatp).to_have_attribute("aria-pressed", "true")
    _pick(rec, settings.get_by_role("combobox", name="LOSATP mode", exact=True), "collinear")
    # The Web app leaves this off; the Gallery command's CLI default is on,
    # and the tutorial's group IDs and group counts come from it.
    rec.caption(4, "Group the proteins, as in the Gallery session", "Infer orthogroups with self-comparisons: on")
    infer = settings.get_by_label("Infer orthogroups with self-comparisons", exact=True)
    rec.reveal(infer)
    rec.focus(infer, zoom=2.2, duration=0.6)
    rec.set_checked(infer, True)
    rec.caption(4, "Set the block rules", "Max unit gap 0 · Min block genes 1 · Evidence scope: Adjacent pairs")
    _value(rec, settings.get_by_label("Collinear max unit gap", exact=True), "0")
    _value(rec, settings.get_by_label("Collinear minimum block genes", exact=True), "1")
    _pick(rec, settings.get_by_label("Collinear evidence scope", exact=True), "adjacent")
    rec.caption(4, "Color blocks by orientation and identity", "Color mode: Orientation + identity")
    _pick(rec, settings.get_by_label("Collinear color mode", exact=True), "orientation_identity")
    rec.caption(4, "Set the match filters", "Bitscore 50 · E-value 1e-2 · Identity 0 · Length 0")
    for label, value in (("minimum bitscore", "50"), ("maximum e-value", "1e-2"),
                         ("minimum identity", "0"), ("minimum alignment length", "0")):
        _value(rec, settings.get_by_label(f"Linear comparison {label}", exact=True), value)
    rec.caption(4, "Limit block merging", "Advanced comparison and layout › Diagonal drift 0 · Merge conflicts 1")
    advanced = page.locator('details[data-linear-comparison-disclosure="advanced"]')
    if advanced.get_attribute("open") is None:
        _press(rec, advanced.get_by_role("button", name="Advanced comparison and layout", exact=True), zoom=2.0)
    _value(rec, advanced.get_by_label("Collinear diagonal drift", exact=True), "0")
    _value(rec, advanced.get_by_label("Collinear merge conflicts", exact=True), "1")

    rec.caption(5, "Set the diagram layout", "Layout › Track Layout: Features on axis · Align to Center · Separate Strands")
    rec.open_section(page.get_by_label("Layout", exact=True))
    _pick(rec, page.get_by_label("Track Layout", exact=True), "middle")
    for name in ("Align to Center", "Separate Strands"):
        box = page.get_by_label(name, exact=True)
        rec.reveal(box)
        rec.focus(box, zoom=2.2, duration=0.5)
        rec.set_checked(box, True, pause=0.4)
    rec.caption(5, "Show both GC tracks", "Layout › Show GC Content and Show GC Skew: on")
    for name in ("Show GC Content", "Show GC Skew"):
        box = page.get_by_label(name, exact=True)
        rec.reveal(box)
        rec.focus(box, zoom=2.2, duration=0.5)
        rec.set_checked(box, True, pause=0.4)
    rec.caption(5, "Add a ruler", "Axis & Scale › Show Coordinate Scale · Scale Style: Ruler (Ticks)")
    rec.open_section(page.get_by_label("Axis & Scale", exact=True))
    show_scale = page.get_by_label("Show Coordinate Scale (Linear)", exact=True)
    rec.reveal(show_scale)
    rec.focus(show_scale, zoom=2.2, duration=0.5)
    rec.set_checked(show_scale, True)
    _pick(rec, page.get_by_label("Linear scale style", exact=True), "ruler")
    rec.caption(5, "Choose the palette", "Colors › Palette: ajisai")
    rec.open_section(page.get_by_label("Colors", exact=True))
    _pick(rec, page.get_by_label("Palette", exact=True), "ajisai")

    rec.caption(6, "Generate and inspect the diagram", "LOSATP compares the proteins in the browser; the wait is fast-forwarded")
    # Four 3-4 Mbp proteomes: slower than the five small Hepatoplasmataceae genomes.
    _generate(rec, target=2.4, look=2.4, timeout_ms=900_000)
    rec.mark("blocks")
    _look(rec, rec.text_bounds(r"^(Collinear|Inverted)$"), zoom=2.2, seconds=1.8, pad=40)
    rec.camera(duration=1.0)
    rec.hold(0.6)

    rec.caption(7, "Inspect collinear blocks", "Click a block ribbon: its records, intervals, orientation, identity, and gene pairs")
    ribbon = _hit(rec, f'[data-collinearity-block-id="{HEP_BLOCK}"]')
    rec.click_point(ribbon["x"], ribbon["y"], pause=0.4)
    details = page.get_by_role("dialog", name="Pairwise match details", exact=True)
    expect(details).to_be_visible()
    expect(details).to_contain_text(HEP_BLOCK)
    rec.camera(details.bounding_box(), zoom=1.5, duration=1.0)
    rec.hold(3.6)
    rec.click(details.get_by_label("Close match popup", exact=True), pause=0.3)
    rec.camera(duration=0.9)
    rec.hold(0.5)

    rec.caption(8, "Inspect a CDS feature", "Click DnaA, then open Details: its record, location, and group ID")
    dnaa = "(item) => item.gene === arg && item.record_id === 'AP027078.1'"
    popup = _open_feature(rec, _feature_point(rec, dnaa, "dnaA"), "Details")
    expect(popup).to_contain_text("og_1")
    rec.hold(3.2)
    _close_feature(rec)

    dnaa_id = _svg_id(rec, dnaa, "dnaA")
    exported = export_interactive_svg(rec, downloads, 9)
    errors = _open_export(rec, exported, 9)
    rec.caption(9, "Click a feature for its similarity group", "DnaA: og_1 has one member in each of the five genomes")
    _svg_close(rec, _svg_feature(rec, dnaa_id, 9, overview=2.6, shows="og_1"))
    # The opened file's block popup reports "Match feature endpoint identity is
    # invalid." for both spans, so the video inspects blocks in the app only.
    rec.caption(9, "Search the features", "Search › Label: ribosomal protein, in all five genomes")
    _svg_search(rec, "ribosomal protein", 9)
    _svg_done(rec, errors)
    return exported


def _inspect_hepatoplasmataceae_collinear(svg: Path) -> dict:
    labels, fills = _labels(svg), _fills(svg)
    text = svg.read_text(encoding="utf-8")
    gallery = _blocks(HEP_GALLERY.read_text(encoding="utf-8"))
    report = {
        "records": sorted(set(HEP_RECORDS) & labels),
        "blocks": len(_blocks(text)),
        "gallery_blocks": len(gallery),
        "example_block": HEP_BLOCK in _blocks(text),
        "legend": sorted({"Collinear", "Inverted"} & labels),
        "gc_tracks": {"GC content", "GC skew (+)", "GC skew (-)"} <= labels,
        "ruler": "600 kbp" in labels,
        "ajisai_cds": "#84b9ec" in fills,
    }
    return _require(report, len(report["records"]) == 5 and _blocks(text) == gallery and report["example_block"]
                    and len(report["legend"]) == 2 and report["gc_tracks"] and report["ruler"]
                    and report["ajisai_cds"], svg)


# -- Vibrio harveyi-group collinear blocks ------------------------------------

VIBRIO_FILES = (
    (REPO_ROOT / "tests/test_inputs/GCF_000196095.1_ASM19609v1_genomic.gbff", 13513689,
     "c785f2feed7dbef0a7ad881aef4a43860f65318e1a44db66b032fc65c9cc34dc"),
    (REPO_ROOT / "tests/test_inputs/GCF_000354175.2_ASM35417v2_genomic.gbff", 13030087,
     "4357443d248ab3d5c08b1854dcbcca1f59adbf064640fd99d5504f6ee4d59096"),
)
VIBRIO_SOURCES = tuple(path for path, _, _ in VIBRIO_FILES)
VIBRIO_RECORDS = ("NC_004603.1", "NC_004605.1", "NC_022349.1", "NC_022359.1")
VIBRIO_ROWS = ("1", "1", "2", "2")
VIBRIO_SPECIES = ((1, "<i>Vibrio parahaemolyticus</i>", "RIMD 2210633"),
                  (3, "<i>Vibrio alginolyticus</i>", "NBRC 15630"))
# (qualifier key, query, the record of each match, in the order the search lists them)
VIBRIO_INITIATORS = (("gene", "dnaA", ("NC_004603.1", "NC_022349.1")),
                     ("product", "RctB", ("NC_004605.1", "NC_022359.1")))
VIBRIO_GALLERY = REPO_ROOT / "gbdraw/web/gallery/sources/vibrio-harveyi-group-collinear.svg"
_LEFTMOST = """
([recordId, svgId]) => {
  const svg = document.querySelector('[role="region"][aria-label="Result Preview"] svg');
  const left = (node) => node.getBoundingClientRect().left;
  const own = [...svg.querySelectorAll('path[data-gbdraw-feature-id]')]
    .filter((node) => node.getAttribute('data-gbdraw-record-id') === recordId);
  const target = svg.querySelector(`[data-gbdraw-feature-id="${svgId}"]`);
  if (!own.length || !target) return null;
  const record = Math.min(...own.map(left));
  return {offset: left(target) - record, width: svg.getBoundingClientRect().width};
}
"""


def _verify_vibrio() -> None:
    for path, size, digest in VIBRIO_FILES:
        assert_fixture_identity(path, expected_size=size, expected_sha256=digest)


def _search_palette(rec: Recorder, shown: bool) -> None:
    """Show or hide the floating feature search, which the recorder hides by default."""

    rec.page.evaluate("(shown) => { if (shown) document.documentElement.dataset.recSearch = '';"
                      " else delete document.documentElement.dataset.recSearch; }", shown)


def _search_status(rec: Recorder) -> str:
    return rec.page.get_by_role("status", name="Feature search status", exact=True).inner_text().strip()


def _frame_box(rec: Recorder, box: dict, *, zoom: float = 1.6, pad: float = 30) -> None:
    area = {"x": box["x"] - pad, "y": box["y"] - pad, "width": box["width"] + 2 * pad, "height": box["height"] + 2 * pad}
    rec.camera(area, zoom=zoom, duration=0.9, portrait=area)


def _rotate_record(rec: Recorder, record_id: str) -> None:
    """In the open feature popup: Rotate record using this feature › Start of the record › Apply on Generate."""

    page = rec.page
    popup = page.get_by_role("dialog", name=re.compile(r"^Feature details"))
    expect(popup).to_be_visible()
    if page.evaluate("() => window.__GBDRAW_APP__?.clickedFeature?.feat?.record_id") != record_id:
        raise AssertionError(f"The popup is not for {record_id}")
    toggle = popup.locator("#feature-record-actions-toggle")
    toggle.evaluate("node => node.scrollIntoView({block: 'center'})")
    page.wait_for_timeout(300)
    _frame_box(rec, popup.bounding_box())
    rec.click(toggle, pause=0.4)
    start = popup.locator("input[name='feature-record-rotation-position'][value='start']")
    start.evaluate("node => node.scrollIntoView({block: 'center'})")
    page.wait_for_timeout(300)
    rec.set_checked(start, True, pause=0.6)
    preview = popup.locator("#feature-record-rotation-preview")
    expect(preview).to_contain_text(f"{record_id} will start at")
    apply = popup.get_by_role("button", name="Apply on Generate", exact=True)
    apply.evaluate("node => node.scrollIntoView({block: 'end'})")
    page.wait_for_timeout(300)
    _frame_box(rec, popup.bounding_box())
    rec.hold(1.6)
    rec.click(apply, pause=0.8)
    rec.hold(1.0)
    close = page.get_by_label("Close feature popup", exact=True)
    if close.count() and close.is_visible():
        rec.click(close, pause=0.3)
    rec.camera(duration=0.9)
    rec.hold(0.4)


def _search_initiator(rec: Recorder, key: str, query: str) -> None:
    """Search › Qualifier value, ``key`` and ``query``: the initiator of both species."""

    page = rec.page
    _search_palette(rec, True)
    palette = page.locator(".preview-feature-search")
    rec.camera(palette.bounding_box(), zoom=2.4, duration=0.9, portrait=palette.bounding_box())
    rec.choose(palette.get_by_label("Search field", exact=True), "qualifier-value")
    rec.type(palette.get_by_label("Qualifier key", exact=True), key)
    rec.type(palette.locator("input[aria-label='Search features']"), query)
    rec.click(palette.locator("button.gfs-button--search"), pause=0.6)
    expect(page.get_by_role("status", name="Feature search status", exact=True)).to_have_text(
        "1 / 2 features", timeout=ACTION_TIMEOUT_MS)
    rec.hold(1.4)


def _rotate_match(rec: Recorder, record_id: str, *, second: bool) -> None:
    """Open the active match (the second one after Next match) and start its record at the feature."""

    page = rec.page
    palette = page.locator(".preview-feature-search")
    if second:
        rec.click(palette.get_by_role("button", name="Next match", exact=True), pause=0.5)
        expect(page.get_by_role("status", name="Feature search status", exact=True)).to_have_text(
            "2 / 2 features", timeout=ACTION_TIMEOUT_MS)
        rec.hold(1.0)
    rec.click(palette.get_by_role("button", name="Open active feature", exact=True), pause=0.6)
    _search_palette(rec, False)
    _rotate_record(rec, record_id)
    _search_palette(rec, True)


def _vibrio_harveyi_group_collinear(rec: Recorder, downloads: Path) -> Path:
    _verify_vibrio()
    page = rec.page

    rec.camera(duration=0)
    rec.caption(1, "Switch to Linear mode", "Two Vibrio genomes, two chromosomes each")
    rec.hold(1.6)
    _select_mode(rec, "Linear")

    rec.caption(2, "Add the two assembly files", "GCF_000196095.1 as File 1 and GCF_000354175.2 as File 2")
    add = page.get_by_role("button", name="Add sequence", exact=True)
    _press(rec, add, zoom=2.0, pause=0.35)
    for index, source in enumerate(VIBRIO_SOURCES, start=1):
        row = page.get_by_role("region", name=f"Linear input file {index}", exact=True)
        chooser = row.get_by_role("button", name=re.compile(r"^Choose GenBank"))
        rec.reveal(chooser)
        rec.focus(chooser, zoom=2.0, duration=0.6)
        rec.upload(chooser, source)
        expect(row).to_contain_text(source.name)
    page.wait_for_function("() => window.__GBDRAW_APP__?.linearSeqs?.length === 4", timeout=GENERATION_TIMEOUT_MS)
    rec.caption(2, "Name the species", "Record options › Organism / strain and Subtitle / title of each File's first record")
    for number, species, strain in VIBRIO_SPECIES:
        file_number = 1 if number == 1 else 2
        records = page.get_by_role("button", name=f"Records for file {file_number}", exact=True)
        rec.open_section(records)
        options = page.get_by_role("button", name=f"Record options for sequence {number}", exact=True)
        rec.open_section(options)
        _field(rec, page.get_by_label(f"Definition for sequence {number}", exact=True), species, zoom=2.0)
        _field(rec, page.get_by_label(f"Subtitle / title for sequence {number}", exact=True), strain, zoom=2.0)

    rec.caption(3, "Arrange one species per row", "Advanced comparison and layout › Arrange in rows · Record gap 48 px")
    advanced = page.locator('details[data-linear-comparison-disclosure="advanced"]')
    rec.open_section(advanced.get_by_role("button", name="Advanced comparison and layout", exact=True))
    rec.set_checked(page.get_by_label("Arrange linear records in rows", exact=True), True)
    for number, row in enumerate(VIBRIO_ROWS, start=1):
        field = page.get_by_label(f"Linear record row for sequence {number}", exact=True)
        if number == 3:
            rec.fast_forward(factor=3.0)
        _value(rec, field, row, zoom=2.0)
    rec.normal_speed()
    _value(rec, page.get_by_label("Linear record gap", exact=True), "48", zoom=2.0)

    rec.caption(4, "Generate a first preview", "Without a comparison it is ready in seconds; the first run starts Python in the browser")
    _generate(rec, target=2.4, look=1.6)

    rec.caption(4, "Search for dnaA", "Search field: Qualifier value · Qualifier key gene · dnaA")
    _search_initiator(rec, "gene", "dnaA")
    rec.caption(5, "Start chromosome I at dnaA",
                "Open › Rotate record using this feature › Start of the record › Apply on Generate")
    _rotate_match(rec, VIBRIO_INITIATORS[0][2][0], second=False)
    rec.caption(5, "Do the same for the second match", "Next match › Open › Start of the record › Apply on Generate")
    _rotate_match(rec, VIBRIO_INITIATORS[0][2][1], second=True)
    rec.caption(6, "Search for rctB", "Search field: Qualifier value · Qualifier key product · RctB")
    _search_initiator(rec, "product", "RctB")
    rec.caption(6, "Start chromosome II at rctB",
                "Open › Rotate record using this feature › Start of the record › Apply on Generate")
    _rotate_match(rec, VIBRIO_INITIATORS[1][2][0], second=False)
    rec.caption(6, "Do the same for the second match", "Next match › Open › Start of the record › Apply on Generate")
    _rotate_match(rec, VIBRIO_INITIATORS[1][2][1], second=True)
    _search_palette(rec, False)

    rec.caption(6, "Generate to apply the new starts", "A checkpoint before the comparison: every record begins at its initiator")
    _generate(rec, target=2.4, look=2.0)
    initiators = (("(item) => item.gene === arg && item.record_id === 'NC_004603.1'", "dnaA", "NC_004603.1"),
                  ("(item) => item.gene === arg && item.record_id === 'NC_022349.1'", "dnaA", "NC_022349.1"),
                  ("(item) => /RctB/i.test(item.product || '') && item.record_id === 'NC_004605.1'", "", "NC_004605.1"),
                  ("(item) => /RctB/i.test(item.product || '') && item.record_id === 'NC_022359.1'", "", "NC_022359.1"))
    for feature_filter, arg, record_id in initiators:
        found = page.evaluate(_LEFTMOST, [record_id, _svg_id(rec, feature_filter, arg)])
        if found is None or found["offset"] > 3:
            raise AssertionError(f"{record_id} does not start at its initiator: {found}")

    rec.caption(7, "Run LOSAT in the browser", "Comparison › Run LOSAT for all adjacent pairs")
    _run_losat(rec)
    settings = page.locator('details[data-linear-comparison-disclosure="settings"]')
    if settings.get_attribute("open") is None:
        _press(rec, settings.get_by_role("button", name="Comparison Settings", exact=True), zoom=2.0)
    losatp = settings.get_by_role("group", name="LOSAT Mode", exact=True).get_by_role("button", name="LOSATP",
                                                                                      exact=True)
    rec.caption(7, "Build collinear blocks", "LOSAT Mode: LOSATP · LOSATP mode: Collinear blocks")
    _press(rec, losatp, zoom=2.2)
    expect(losatp).to_have_attribute("aria-pressed", "true")
    _pick(rec, settings.get_by_role("combobox", name="LOSATP mode", exact=True), "collinear")
    # The Web app leaves this off; the Gallery session has it on.
    rec.caption(7, "Group the proteins, as in the Gallery session", "Infer orthogroups with self-comparisons: on")
    infer = settings.get_by_label("Infer orthogroups with self-comparisons", exact=True)
    rec.reveal(infer)
    rec.focus(infer, zoom=2.2, duration=0.6)
    rec.set_checked(infer, True)
    rec.caption(7, "Set the block rules", "Max unit gap 2 · Min block genes 3 · Evidence scope: Adjacent pairs")
    _value(rec, settings.get_by_label("Collinear max unit gap", exact=True), "2")
    _value(rec, settings.get_by_label("Collinear minimum block genes", exact=True), "3")
    _pick(rec, settings.get_by_label("Collinear evidence scope", exact=True), "adjacent")
    rec.caption(7, "Color blocks by orientation and identity", "Color mode: Orientation + identity")
    _pick(rec, settings.get_by_label("Collinear color mode", exact=True), "orientation_identity")
    rec.caption(7, "Set the match filters", "Bitscore 50 · E-value 1e-2 · Identity 0 · Length 0")
    for label, value in (("minimum bitscore", "50"), ("maximum e-value", "1e-2"),
                         ("minimum identity", "0"), ("minimum alignment length", "0")):
        _value(rec, settings.get_by_label(f"Linear comparison {label}", exact=True), value)
    rec.caption(7, "Limit block merging", "Advanced comparison and layout › Diagonal drift 2 · Merge conflicts 1")
    if advanced.get_attribute("open") is None:
        _press(rec, advanced.get_by_role("button", name="Advanced comparison and layout", exact=True), zoom=2.0)
    _value(rec, advanced.get_by_label("Collinear diagonal drift", exact=True), "2")
    _value(rec, advanced.get_by_label("Collinear merge conflicts", exact=True), "1")

    rec.caption(8, "Render every feature as a rectangle", "Features › Rendering of all seven feature types: Rectangle")
    rec.open_section(page.locator('summary[aria-label="Features"]'))
    shapes = page.locator('select[aria-label^="Rendering for "]')
    expect(shapes).to_have_count(7)
    for index in range(7):
        if index == 1:
            rec.fast_forward(factor=3.0)
        _pick(rec, shapes.nth(index), "rectangle", zoom=2.0)
    rec.normal_speed()

    rec.caption(9, "Set the layout", "Layout › Track Layout: Above · Separate Strands · Lock Definition Column")
    rec.open_section(page.get_by_label("Layout", exact=True))
    _pick(rec, page.get_by_label("Track Layout", exact=True), "above")
    for name in ("Separate Strands", "Lock Definition Column"):
        box = page.get_by_label(name, exact=True)
        rec.reveal(box)
        rec.focus(box, zoom=2.2, duration=0.5)
        rec.set_checked(box, True, pause=0.4)
    rec.caption(9, "Draw a 750,000 bp ruler on each axis", "Axis & Scale › Scale Style: Ruler (Ticks) · Ruler on Axis · Scale Interval 750000")
    rec.open_section(page.get_by_label("Axis & Scale", exact=True))
    _pick(rec, page.get_by_label("Linear scale style", exact=True), "ruler")
    on_axis = page.get_by_label("Ruler on Axis", exact=True)
    rec.reveal(on_axis)
    rec.focus(on_axis, zoom=2.2, duration=0.5)
    rec.set_checked(on_axis, True)
    _value(rec, page.get_by_label("Linear ruler interval", exact=True), "750000")
    rec.caption(9, "Show each definition as species and strain", "Titles & Record Labels › Default font size 16 · Replicon off · Accession and Length hidden")
    rec.open_section(page.get_by_label("Titles and Record Labels", exact=True))
    _value(rec, page.get_by_label("Default font size", exact=True), "16")
    for summary, size, weight in (("Name / Species", "18", "bold"), ("Subtitle", "16", "normal")):
        rec.open_section(page.locator(f'summary[aria-label="{summary} style"]'))
        _value(rec, page.get_by_label(f"{summary} style size", exact=True), size, zoom=2.0)
        _press(rec, page.get_by_label(f"{summary} {weight} weight", exact=True), zoom=2.0, pause=0.4)
    rec.set_checked(page.get_by_label("Replicon visibility", exact=True), False)
    for name in ("Accession visibility", "Length / Coordinates visibility"):
        _pick(rec, page.get_by_label(name, exact=True), "hide")

    rec.caption(10, "Place the legend at the bottom", "Plot Title stays empty · Legend position: Bottom")
    _point(rec, page.get_by_role("textbox", name="Plot Title", exact=True), zoom=2.0)
    rec.open_section(page.get_by_label("Legend settings", exact=True))
    _pick(rec, page.locator("#legend-position"), "bottom")

    rec.caption(11, "Generate once and inspect", "LOSATP compares the proteins in the browser; the wait is fast-forwarded")
    # Four 3-4 Mbp proteomes: slower than the five small Hepatoplasmataceae genomes.
    _generate(rec, target=2.4, look=2.4, timeout_ms=900_000)
    rec.mark("blocks")
    _look(rec, rec.text_bounds(r"^(Collinear|Inverted)$"), zoom=2.2, seconds=1.8, pad=40)
    rec.camera(duration=1.0)
    rec.hold(0.6)

    exported = export_interactive_svg(rec, downloads, 12)
    errors = _open_export(rec, exported, 12)
    rec.caption(12, "Find an initiator in the file", "Search › Qualifier value: gene dnaA, in both species")
    _svg_search(rec, "dnaA", 12, qualifier="gene")
    _svg_done(rec, errors)
    return exported


def _inspect_vibrio_harveyi_group_collinear(svg: Path) -> dict:
    labels = _labels(svg)
    text = svg.read_text(encoding="utf-8")
    gallery = _blocks(VIBRIO_GALLERY.read_text(encoding="utf-8"))
    records = set(re.findall(r'data-gbdraw-record-id="([^"]+)"', text))
    report = {
        "records": sorted(set(VIBRIO_RECORDS) & records),
        "blocks": len(_blocks(text)),
        "gallery_blocks": len(gallery),
        "species": sorted({"Vibrio parahaemolyticus", "Vibrio alginolyticus", "RIMD 2210633", "NBRC 15630"} & labels),
        "legend": sorted({"Collinear", "Inverted"} & labels),
        "ruler": {"0.75 Mbp", "1.5 Mbp"} <= labels,
        "same_blocks": _blocks(text) == gallery,
    }
    return _require(report, len(report["records"]) == 4 and len(report["species"]) == 4
                    and len(report["legend"]) == 2 and report["ruler"] and report["blocks"] > 0
                    and report["same_blocks"], svg)


# Each tutorial ends in its exported Interactive SVG, opened in the browser,
# rather than in the vector zoom.
_LANDSCAPE_TAKES = (
    Journey("HmmtDNA_ATskew", (HMMT,), _hmmtdna_atskew, _inspect_hmmtdna_atskew, zoom_text="8 kbp",
            vector_zoom=False),
    Journey("tobacco-chloroplast", (PLASTOME, PLASTOME_COLORS, PLASTOME_PRIORITY, PLASTOME_REGIONS),
            _tobacco_chloroplast, _inspect_tobacco_chloroplast, zoom_text="IRb", vector_zoom=False),
    Journey("BGC0000708-BGC0000713", BGC_SOURCES, _bgc_aminoglycosides, _inspect_bgc_aminoglycosides,
            zoom_text="livE", vector_zoom=False),
    Journey("Vnig_TUMSAT-TG-2018", (VNIG,), _vnig_multirecord, _inspect_vnig_multirecord, zoom_text="pVNTG3",
            vector_zoom=False),
    Journey("hepatoplasmataceae_collinear", HEP_SOURCES, _hepatoplasmataceae_collinear,
            _inspect_hepatoplasmataceae_collinear, zoom_text="AP027078.1", vector_zoom=False),
    Journey("vibrio-harveyi-group-collinear", VIBRIO_SOURCES, _vibrio_harveyi_group_collinear,
            _inspect_vibrio_harveyi_group_collinear, zoom_text="NC_004603.1", vector_zoom=False),
)
# The tall video is recorded again, so the preview figure and the opened file
# can take a tall window's shape; the steps and checks are the same.
TUTORIALS = _LANDSCAPE_TAKES + tuple(replace(journey, name=journey.name + VERTICAL, portrait=True)
                                     for journey in _LANDSCAPE_TAKES)
