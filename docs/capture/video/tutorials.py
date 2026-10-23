"""Record the Gallery tutorials as human-paced Web walkthroughs.

Each journey follows one tutorial in ``gbdraw/web/gallery/tutorials`` from an
empty app to the exported SVG. Captions use the tutorial's step numbers and
its UI labels, so a viewer can follow the written steps beside the video.
Generate is pressed at a few extra checkpoints to show each change as it
lands; the final settings are the tutorial's own.
"""

from __future__ import annotations

import re
from collections.abc import Callable
from dataclasses import replace
from pathlib import Path

from playwright.sync_api import Locator, expect

from config import ACTION_TIMEOUT_MS, GENERATION_TIMEOUT_MS, REPO_ROOT, tutorial_file_identity
from flows.web_capture import assert_fixture_identity
from video.walkthrough import (PORTRAIT_FIGURE_WIDTH, PREVIEW_REGION, VIEWPORT, _PREVIEW, Journey, Recorder, _center,
                               export_interactive_svg, open_in_browser)
from video.walkthrough_render import PORTRAIT, TUTORIAL_LANGUAGES, take_name


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


# The union of the preview drawing's elements that match a selector, clipped
# to the visible canvas, in page CSS px.
_DETAIL = """
(selector) => {
  const region = document.querySelector('$REGION');
  const view = region.getBoundingClientRect();
  let left = Infinity, top = Infinity, right = -Infinity, bottom = -Infinity;
  for (const node of region.querySelectorAll(selector)) {
    const b = node.getBoundingClientRect();
    if (!node.closest('svg') || !b.width || !b.height) continue;
    left = Math.min(left, b.left); top = Math.min(top, b.top);
    right = Math.max(right, b.right); bottom = Math.max(bottom, b.bottom);
  }
  left = Math.max(left, view.left); top = Math.max(top, view.top);
  right = Math.min(right, view.right); bottom = Math.min(bottom, view.bottom);
  return right > left && bottom > top ? {x: left, y: top, width: right - left, height: bottom - top} : null;
}
""".replace("$REGION", PREVIEW_REGION)

_HIT = """
(selector) => {
  // In an opened Interactive SVG the drawing is the whole page.
  const region = document.querySelector('$REGION') || document.documentElement;
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
""".replace("$REGION", PREVIEW_REGION)


# -- shared steps -------------------------------------------------------------

def _css(rec: Recorder, attribute: str, text: str, params: dict | None = None) -> str:
    """A CSS attribute test for the UI text ``text``, in the app's current language."""

    value = rec.t(text, params).replace("\\", "\\\\").replace('"', '\\"')
    return f'[{attribute}="{value}"]'


def _summary(rec: Recorder, name: str, params: dict | None = None) -> Locator:
    """A settings disclosure by its (translated) name."""

    return rec.page.locator(f"summary{_css(rec, 'aria-label', name, params)}")


def _label(rec: Recorder, name: str, params: dict | None = None, *, within: Locator | None = None) -> Locator:
    return (within or rec.page).get_by_label(rec.t(name, params), exact=True)


def _button(rec: Recorder, name: str, params: dict | None = None, *, within: Locator | None = None) -> Locator:
    return (within or rec.page).get_by_role("button", name=rec.t(name, params), exact=True)


def _select_mode(rec: Recorder, name: str) -> None:
    button = _button(rec, name)
    rec.focus(button, zoom=2.2)
    rec.hold(0.4)
    rec.click(button, pause=0.6)
    expect(button).to_have_attribute("aria-pressed", "true")


def _fit_preview(rec: Recorder, *, margin: float = 24, skim: float | None = None) -> None:
    """Zoom out with the preview's own control until the figure fits above its buttons, then center it.

    A portrait take also narrows the figure until the tall window can show all
    of it. It frames the canvas meanwhile and skims the repeated clicks under
    the fast-forward badge, so the new figure is not seen cut off: at 3x, or
    in about ``skim`` seconds when given.
    """

    page = rec.page
    zoom_out = _button(rec, "Zoom out")
    skimming = False
    for _ in range(16):
        geometry = page.evaluate(_PREVIEW)
        canvas, svg = geometry["canvas"], geometry["svg"]
        width = min(canvas["width"] - margin, PORTRAIT_FIGURE_WIDTH) if rec.portrait else canvas["width"] - margin
        if svg["width"] <= width and svg["height"] <= canvas["height"] - margin:
            break
        if rec.portrait and not skimming:
            rec.camera(canvas, duration=0.5)
            rec.fast_forward(factor=None if skim else 3.0, target=skim)
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
        _pan_preview(rec, dx, dy)
        geometry = page.evaluate(_PREVIEW)
        canvas, svg = geometry["canvas"], geometry["svg"]
        off = canvas["x"] + canvas["width"] / 2 - (svg["x"] + svg["width"] / 2)
        if abs(off) > 40:
            raise AssertionError(f"The preview pan did not center the figure: {geometry}")
    if skimming:
        rec.normal_speed()


def _pan_preview(rec: Recorder, dx: float, dy: float, *, duration: float = 0.9) -> None:
    """Drag the preview's blank background by ``(dx, dy)`` CSS px, as a person pans it.

    Leaving the canvas ends a pan, so each drag ends inside it; a pan longer
    than one drag can reach takes another from a fresh blank spot.
    """

    for _ in range(4):
        if abs(dx) < 2 and abs(dy) < 2:
            return
        geometry = rec.page.evaluate(_PREVIEW)
        canvas = geometry["canvas"]
        left, top = canvas["x"] + 16, canvas["y"] + 16
        right, bottom = canvas["x"] + canvas["width"] - 16, canvas["y"] + canvas["height"] - 16

        def reach(x: float, y: float) -> float:
            """The share of the pan a drag from ``(x, y)`` can make without leaving the canvas."""

            share = 1.0
            for start, delta, low, high in ((x, dx, left, right), (y, dy, top, bottom)):
                if delta:
                    share = min(share, max(0.0, ((high if delta > 0 else low) - start) / delta))
            return share

        share, start = max(((reach(*point), point) for point in geometry["blank"]), key=lambda row: row[0],
                           default=(0.0, None))
        if share < 0.05:
            raise AssertionError(f"No blank preview area to pan from: {geometry}")
        rec.drag(start, (start[0] + dx * share, start[1] + dy * share), duration=duration * max(share, 0.5))
        rec.page.wait_for_timeout(300)
        dx, dy = dx * (1 - share), dy * (1 - share)
    raise AssertionError(f"The preview pan did not finish: {dx:.0f}, {dy:.0f} px left")


def _generate(rec: Recorder, *, target: float | None = 1.4, look: float = 2.0,
              timeout_ms: int = GENERATION_TIMEOUT_MS, stages: bool = False) -> None:
    """Generate, fit the figure, and show the whole page for ``look`` seconds.

    ``stages`` frames the progress overlay's stage line meanwhile; ``target``
    None plays the wait at normal speed.
    """

    button = _button(rec, "Generate Diagram")
    rec.focus(button, zoom=1.9)
    rec.hold(0.4)
    rec.generate(target=target, timeout_ms=timeout_ms, stages=stages)
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

    field = _label(rec, "New priority rule qualifier order")
    rec.reveal(field)
    rec.focus(field, zoom=2.2)
    rec.choose(_label(rec, "New priority rule feature type"), feature)
    rec.type(field, order)
    rec.click(_button(rec, "Add priority rule"), pause=0.6)


def _color_rules(rec: Recorder, rules: tuple[tuple[str, str, str, str, str], ...]) -> None:
    """Colors › Specific rules; the first is shown at pace, the rest at 3×."""

    regex = rec.page.locator("input[data-new-specific-rule-pattern]")
    rule_box = rec.page.locator("div:has(> div > input[data-new-specific-rule-pattern])")
    for index, (feature, qualifier_name, pattern, color, legend) in enumerate(rules):
        rec.reveal(rule_box, anchor=0.45)
        rec.mark(f"rule-{index + 1}")
        rec.focus(rule_box, zoom=2.2, duration=0.6)
        if index == 1:
            rec.fast_forward(factor=3.0)
        rec.choose(_label(rec, "New color rule feature type"), feature)
        qualifier = _label(rec, "Qualifier name", within=rule_box)
        if qualifier.input_value() != qualifier_name:
            rec.type(qualifier, qualifier_name)
        rec.type(regex, pattern)
        rec.pick_color(_label(rec, "New color rule color"), color)
        rec.type(_label(rec, "New color rule legend name"), legend)
        rec.click(_button(rec, "Add color rule"), pause=0.45)
        # Each rule is validated before it is committed; a quicker next rule
        # would be validated without this one.
        rec.page.wait_for_function(
            "(count) => (window.__GBDRAW_APP__?.manualSpecificRules || []).length === count", arg=index + 1)
    rec.normal_speed()
    rec.hold(0.5)


def _feature_types(rec: Recorder, *, remove: tuple[str, ...] = (), add: tuple[str, ...] = ()) -> None:
    """Features › Feature Types: drop and add the drawn feature types."""

    page = rec.page
    summary = _summary(rec, "Features")
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

    _pick(rec, _label(rec, "{name} mode", {"name": rec.t(label)}), "color", zoom=zoom)
    rec.pick_color(_label(rec, label), color)


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
    rec.click(popup.get_by_role("button", name=re.compile(rf"(^|\s){re.escape(rec.t(tab))}$")), pause=pause)


def _open_feature(rec: Recorder, point: dict, tab: str) -> Locator:
    """Click a drawn feature, switch its popup to ``tab``, and frame both."""

    page = rec.page
    rec.click_point(point["x"], point["y"], pause=0.4)
    popup = page.locator(".feature-popup")
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
    rec.click(_label(rec, "Close feature popup"), pause=0.3)
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

    button = _button(rec, "Run LOSAT for all adjacent pairs")
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


def _open_export(rec: Recorder, svg: Path, step: int | None) -> list[str]:
    rec.caption(step, "Open it in any web browser", "No server, no install")
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
    # Back to the whole drawing, and long enough there that the take never
    # ends while the camera is still easing out of the popup's zoom. The
    # close button sits over the drawing; the pointer leaves it so no hover
    # tooltip from an unrelated ribbon or feature shows meanwhile.
    rec.camera(duration=0.9)
    rec.move_to(60, rec.page.viewport_size["height"] - 40)
    rec.hold(1.0)


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
    preview = page.locator(PREVIEW_REGION)

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
    preview = page.locator(PREVIEW_REGION)

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


# -- shared steps of the progressive Linear tutorials ---------------------------
#
# These journeys follow their tutorial JSON step by step: Generate with the
# defaults first, look at each preview, name what it lacks, change it, and
# Generate again. Every UI name in a locator goes through ``rec.t``, so the same
# steps record the English and the translated app.

def _frame_figure(rec: Recorder) -> dict:
    """Frame the whole preview figure; returns its stage bounds."""

    figure = rec.figure()
    if figure is None:
        raise AssertionError("No preview figure to show")
    zoom = min(VIEWPORT[0] / figure["width"], VIEWPORT[1] / figure["height"]) / 1.08
    rec.camera(figure, zoom=max(1.0, min(zoom, 2.2)), duration=1.0)
    return figure


def _show_figure(rec: Recorder, step: int | None, title: str, subtitle: str = "", *, seconds: float = 3.4,
                 detail: Callable[[Recorder], dict] | None = None) -> None:
    """Name what the new preview shows while the camera holds the whole figure.

    ``detail`` returns the part of a figure too small to read whole (a page
    CSS px rect, such as one column of records): the wide window then holds
    it about 3 s at a camera zoom where it reads, and the tall window reads
    it closely in the app (``_read_closely``). Then the whole figure again.
    """

    rec.caption(step, title, subtitle)
    figure = _frame_figure(rec)
    bottom = figure["y"] + figure["height"]
    if figure["x"] <= rec.x <= figure["x"] + figure["width"] and figure["y"] <= rec.y <= bottom:
        # A pointer resting on the legend or a label opens the preview's
        # "Turn on Layout edit" hint over the figure being shown.
        rec.move_to(rec.x - rec.origin[0], bottom - rec.origin[1] + 10)
    rec.hold(seconds)
    if detail:
        part = detail(rec)
        if rec.portrait:
            _read_closely(rec, part)
        else:
            zoom = min(VIEWPORT[0] / part["width"], VIEWPORT[1] / part["height"]) / 1.08
            rec.camera(part, zoom=max(1.0, min(zoom, 3.0)), duration=1.0)
            # The move eases for 1 s, so the part then stays still about 3 s.
            rec.hold(4.0)
        _frame_figure(rec)
        rec.hold(1.6)


# A portrait detail hold zooms the preview until the record labels' capitals
# are this tall in the output, as they read in the landscape Vibrio detail.
DETAIL_CAP_PX = 11
_CAP_HEIGHT = 0.716   # Liberation Sans capitals, in em
# The smallest record label whose center lies in a page rect: its font size in page CSS px.
_RECORD_LABEL = """
(box) => {
  const region = document.querySelector('$REGION');
  let size = Infinity;
  for (const node of region.querySelectorAll('svg [data-gbdraw-role^="record-definition"] text')) {
    const b = node.getBoundingClientRect(), user = node.getBBox().height;
    const x = b.left + b.width / 2, y = b.top + b.height / 2;
    if (!b.height || !user || x < box.x || x > box.x + box.width || y < box.y || y > box.y + box.height) continue;
    size = Math.min(size, parseFloat(getComputedStyle(node).fontSize) * b.height / user);
  }
  return Number.isFinite(size) ? size : null;
}
""".replace("$REGION", PREVIEW_REGION)


def _detail_view(part: dict, label: float) -> tuple[float, float]:
    """The page height a tall-window view of ``part`` shows, and the width of its stops.

    The view shows the part's full height and is tall enough that ``label``
    (a font size in page CSS px) reads at ``DETAIL_CAP_PX``; each stop is the
    widest target that view shows whole.
    """

    readable = label * _CAP_HEIGHT * PORTRAIT.card[3] / DETAIL_CAP_PX
    rows = {"x": part["x"], "y": part["y"], "width": 0, "height": part["height"]}
    height = max(PORTRAIT.portrait_view(VIEWPORT, rows, "bottom")[2], readable)
    return height, PORTRAIT.widest_target(height)


def _read_closely(rec: Recorder, part: dict, *, margin: float = 16) -> None:
    """Read ``part`` in the tall window: zoom in with the preview's own control, pan across, fit again.

    The fitted figure is too small for the tall window to read, and the
    frames hold only 2 device px per CSS px, so the app draws the part
    larger: "Zoom in" (skimmed under the fast-forward badge) to the largest
    zoom whose view still fits inside the canvas, a pan that puts the part's
    left end in the view, about 3 s still, a slow pan to its right end and
    another 3 s, then ``_fit_preview``. The camera holds one place meanwhile.
    """

    page = rec.page
    geometry = page.evaluate(_PREVIEW)
    canvas, svg, zoom = geometry["canvas"], geometry["svg"], geometry["zoom"]
    label = page.evaluate(_RECORD_LABEL, part)
    if not label:
        raise AssertionError(f"No record label in the detail part {part}")
    # The part as shares of the drawing, which zooming and panning keep.
    share = {"x": (part["x"] - svg["x"]) / svg["width"], "y": (part["y"] - svg["y"]) / svg["height"],
             "width": part["width"] / svg["width"], "height": part["height"] / svg["height"]}
    clicks = 0
    for more in range(1, 49):
        if zoom + more / 10 > 5.0 + 1e-6:
            break
        scale = (zoom + more / 10) / zoom
        height, _ = _detail_view({"x": 0, "y": 0, "width": part["width"] * scale,
                                  "height": part["height"] * scale}, label * scale)
        if height > canvas["height"] - 2 * margin or height * PORTRAIT.aspect > canvas["width"] - 2 * margin:
            break
        clicks = more
    rec.camera(canvas, duration=0.5, portrait=canvas)
    rec.fast_forward(factor=3.0)
    zoom_in = _button(rec, "Zoom in")
    for _ in range(clicks):
        rec.click(zoom_in, pause=0.35)
    page.wait_for_timeout(300)

    def located() -> dict:
        svg = page.evaluate(_PREVIEW)["svg"]
        return {"x": svg["x"] + share["x"] * svg["width"], "y": svg["y"] + share["y"] * svg["height"],
                "width": share["width"] * svg["width"], "height": share["height"] * svg["height"]}

    now = located()
    label = page.evaluate(_RECORD_LABEL, now)
    height, width = _detail_view(now, label)
    # The stops' left edges: the part's two ends, or its middle when one view holds it.
    stops = ([now["x"], now["x"] + now["width"] - width] if now["width"] > width
             else [now["x"] + (now["width"] - width) / 2])

    def stop(part: dict, x: float) -> dict:
        return {"x": x, "y": part["y"], "width": width, "height": part["height"]}

    cx, cy, _ = PORTRAIT.portrait_view(VIEWPORT, stop(now, stops[0]), "bottom")
    _pan_preview(rec, canvas["x"] + canvas["width"] / 2 - cx, canvas["y"] + canvas["height"] / 2 - cy)
    # The pointer waits on the zoom buttons, below the view and off the drawing.
    rec.move_to(*_center(zoom_in.bounding_box()))
    rec.normal_speed()
    moved = located()
    first = stop(moved, moved["x"] + stops[0] - now["x"])
    rec.event("detail", zoom=round(zoom + clicks / 10, 1), clicks=clicks, stops=len(stops),
              cap_px=round(label * _CAP_HEIGHT * PORTRAIT.card[3] / PORTRAIT.portrait_view(VIEWPORT, first, "bottom")[2], 1))
    rec.camera(first, duration=1.0, portrait=first)
    # The move eases for 1 s, so the part then stays still about 3 s (the
    # pointer's way to the next pan's start, beside the view, adds to it).
    rec.hold(4.0 if len(stops) == 1 else 3.0)
    if len(stops) > 1:
        _pan_preview(rec, stops[0] - stops[1], 0, duration=1.2)
        rec.move_to(*_center(zoom_in.bounding_box()))
        rec.hold(3.0)
    _fit_preview(rec)


def _preview_part(selector: str) -> Callable[[Recorder], dict]:
    """A ``detail`` of the preview drawing's elements that match ``selector``."""

    def part(rec: Recorder) -> dict:
        box = rec.page.evaluate(_DETAIL, selector)
        if not box:
            raise AssertionError(f"No preview element matches {selector!r}")
        return box
    return part


def _linear_file(rec: Recorder, index: int, source: Path) -> Locator:
    """Upload ``source`` into File ``index`` (one-based) and return its card."""

    page = rec.page
    card = page.get_by_role("region", name=rec.t("Linear input file {n}", {"n": index}), exact=True)
    chooser = card.locator(".upload-zone").first
    rec.reveal(chooser)
    rec.focus(chooser, zoom=2.0, duration=0.6)
    rec.hold(0.3)
    rec.upload(chooser, source)
    expect(card).to_contain_text(source.name, timeout=GENERATION_TIMEOUT_MS)
    rec.hold(0.6)
    return card


def _add_sequence(rec: Recorder) -> None:
    _press(rec, _button(rec, "Add sequence"), zoom=2.0, pause=0.5)


def _losat_settings(rec: Recorder) -> Locator:
    """Comparison › Run LOSAT, open its Settings, and press LOSATP."""

    page = rec.page
    _run_losat(rec)
    settings = page.locator('details[data-linear-comparison-disclosure="settings"]')
    rec.open_section(_summary(rec, "Comparison Settings"))
    losatp = settings.locator("[data-linear-comparison-losat-mode-option='blastp']")
    _press(rec, losatp, zoom=2.2)
    expect(losatp).to_have_attribute("aria-pressed", "true")
    return settings


def _check(rec: Recorder, box: Locator, checked: bool = True) -> None:
    rec.reveal(box)
    rec.focus(box, zoom=2.2, duration=0.5)
    rec.set_checked(box, checked, pause=0.5)


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
# Titles & Record Labels › definition-line styles: (line, font size, bold, fill)
BGC_DEFINITION_STYLES = (
    ("Name / Species", "20", True, None),
    ("Subtitle", "20", False, None),
    ("Accession", "20", False, "#7b7c7d"),
    ("Length / Coordinates", "20", False, "#7b7c7d"),
)
_ALIGNMENT_READY = """
() => { const app = window.__GBDRAW_APP__;
  const row = app?.similarityAlignmentDraft?.rows?.find((item) => item.recordDetail === 'BGC0000708');
  const chosen = row?.candidates?.find((candidate) => candidate.key === row?.choice?.candidateKey);
  return Boolean(row && !app.similarityAlignmentBusy && row.choice?.kind === 'select'
    && row.choice.candidateKey === row.recommendedKey && String(chosen?.name || '') === 'livU'
    && document.querySelector('[data-similarity-alignment-dialog]')); }
"""
_ALIGNMENT_APPLIED = """
(before) => { const app = window.__GBDRAW_APP__; const content = String(app?.results?.[0]?.content || '');
  return Boolean(app && !app.processing && !app.similarityAlignmentBusy
    && !document.querySelector('[data-similarity-alignment-dialog]') && content && content !== before); }
"""


def _bgc_record(rec: Recorder, card: Locator, index: int) -> None:
    """Organism / strain and Subtitle / title on the File card of a one-record file."""

    definition, subtitle = BGC_RECORDS[index - 1]
    _field(rec, _label(rec, "Definition for sequence {n}", {"n": index}, within=card), definition, zoom=2.0)
    _field(rec, _label(rec, "Subtitle / title for sequence {n}", {"n": index}, within=card), subtitle, zoom=2.0)


def _bgc_align(rec: Recorder) -> None:
    """Step 7: Review alignment options… on neoU (og_6), keep livU, Apply."""

    page = rec.page
    rec.caption(7, "Align the rows on og_6", "Click neoU of BGC0000709, in row 2")
    neou = _feature_point(rec, "(item) => item.gene === arg && item.orthogroupId === 'og_6'"
                               " && /^BGC0000709/.test(item.record_id || '')", "neoU")
    rec.mark("align")
    rec.click_point(neou["x"], neou["y"], pause=0.4)
    popup = page.locator(".feature-popup")
    expect(popup).to_be_visible()
    expect(popup).to_contain_text("og_6")
    review = popup.locator("[data-similarity-alignment-popup-review]")
    review.evaluate("node => node.scrollIntoView({block: 'center'})")
    page.wait_for_timeout(300)
    box = popup.bounding_box()
    rec.camera({"x": min(box["x"], neou["x"]) - 40, "y": box["y"] - 40,
                "width": max(box["x"] + box["width"], neou["x"]) - min(box["x"], neou["x"]) + 80,
                "height": box["height"] + 80}, zoom=1.5, duration=0.9,
               portrait={"x": box["x"] - 30, "y": box["y"] - 30, "width": box["width"] + 60,
                         "height": box["height"] + 60})
    rec.hold(1.6)
    rec.caption(7, "Review alignment options…", "One card per record, with its og_6 candidates")
    rec.focus(review, zoom=2.4, duration=0.6)
    before = page.evaluate("() => String(window.__GBDRAW_APP__.results?.[0]?.content || '')")
    rec.click(review, pause=0.4)
    dialog = page.locator("[data-similarity-alignment-dialog]")
    expect(dialog).to_be_visible(timeout=ACTION_TIMEOUT_MS)
    page.wait_for_function(_ALIGNMENT_READY, timeout=GENERATION_TIMEOUT_MS)
    page.wait_for_timeout(600)
    rec.camera(dialog.bounding_box(), zoom=1.4, duration=0.8)
    rec.hold(2.4)
    # BGC0000708 has two og_6 members: the Recommended livU and CAG38692.1.
    card = dialog.locator("[data-alignment-record-key]").filter(has_text="CAG38692.1")
    expect(card).to_have_count(1)
    expect(card).to_contain_text("livU")
    rec.caption(7, "Keep livU (Recommended)", "CAG38692.1 is the other og_6 member")
    rec.camera(card.bounding_box(), zoom=1.8, duration=0.8)
    rec.hold(3.0)
    direction = dialog.locator("[data-similarity-alignment-direction-mode]")
    rec.caption(7, "Arrow direction: Keep as is", "Apply regenerates the preview")
    rec.camera(direction.bounding_box(), zoom=2.0, duration=0.8)
    expect(direction.get_by_role("radio", name=rec.t("Keep as is"), exact=True)).to_be_checked()
    rec.hold(2.0)
    apply = dialog.locator('button[aria-describedby="similarity-alignment-apply-reason"]')
    rec.focus(apply, zoom=2.0, duration=0.6)
    rec.click(apply, pause=0.2)
    rec.camera(duration=0.8)
    rec.fast_forward(target=1.4)
    page.wait_for_function(_ALIGNMENT_APPLIED, arg=before, timeout=GENERATION_TIMEOUT_MS)
    rec.hold(0.3)
    rec.normal_speed()
    # The neoU popup, hidden during the review, comes back after Apply.
    expect(popup).to_be_visible()
    rec.caption(7, "Close the neoU popup", "")
    _close_feature(rec)
    expect(popup).to_be_hidden()
    _fit_preview(rec)
    _show_figure(rec, 7, "The og_6 genes line up", "One column across the five clusters", detail=_bgc_detail(rec))


def _bgc_finish(rec: Recorder) -> None:
    """Step 10: arrow geometry, sizes, strokes, ruler, legend, title, and definition lines."""

    page = rec.page
    rec.caption(10, "Finish the figure", "Arrow Geometry: Head 0.64 · Shaft 0.6")
    rec.open_section(_summary(rec, "Features"))
    _field(rec, _label(rec, "Arrow head length ratio"), "0.64")
    _field(rec, _label(rec, "Arrow shaft width ratio"), "0.6")
    rec.caption(10, "Fix the feature sizes", "Feature Height 75 · strokes #262626, 2 and 2")
    rec.fast_forward(factor=3.0)
    _field(rec, _label(rec, "Feature Height"), "75")
    _color_value(rec, "Block stroke color", "#262626")
    _field(rec, _label(rec, "Block Stroke Width"), "2")
    _field(rec, _label(rec, "Line Stroke Width"), "2")
    rec.normal_speed()
    rec.caption(10, "One row of arrows per cluster", "Layout › Separate Strands: Off")
    rec.open_section(_summary(rec, "Layout"))
    _check(rec, _label(rec, "Separate Strands"), False)
    rec.caption(10, "A gray axis with a ruler", "Axis & Scale › #808080, 5 · Ruler (Ticks) · 20")
    rec.open_section(_summary(rec, "Axis & Scale"))
    rec.fast_forward(factor=3.0)
    _color_value(rec, "Axis stroke color", "#808080")
    _field(rec, _label(rec, "Axis Stroke Width"), "5")
    _pick(rec, _label(rec, "Linear scale style"), "ruler")
    _field(rec, _label(rec, "Linear scale font size"), "20")
    rec.normal_speed()
    rec.caption(10, "Fix the text sizes", "Label Font Size 18 · Legend swatch and font 20")
    rec.open_section(_summary(rec, "Labels"))
    _field(rec, _label(rec, "Label Font Size"), "18")
    rec.open_section(_summary(rec, "Legend settings"))
    rec.fast_forward(factor=3.0)
    _field(rec, _label(rec, "Legend swatch size"), "20")
    _field(rec, _label(rec, "Legend font size"), "20")
    rec.normal_speed()
    rec.caption(10, "Add the plot title", "Titles & Record Labels › Plot Title")
    rec.open_section(_summary(rec, "Titles and Record Labels"))
    _field(rec, page.get_by_role("textbox", name=rec.t("Plot Title"), exact=True), BGC_TITLE, zoom=2.0)
    rec.caption(10, "Style the definition lines", "Size 20 · bold names · gray accession and length")
    rec.fast_forward(factor=3.0)
    for line, size, bold, fill in BGC_DEFINITION_STYLES:
        rec.open_section(_summary(rec, "{label} style", {"label": rec.t(line)}))
        _field(rec, _label(rec, "{label} font size", {"label": rec.t(line)}), size, zoom=2.0)
        if bold:
            _press(rec, _button(rec, "{label} bold weight", {"label": rec.t(line)}), zoom=2.0, pause=0.3)
        if fill:
            _field(rec, _label(rec, "{label} definition line fill value", {"label": rec.t(line)}), fill, zoom=2.0)
    rec.normal_speed()


# The five BGC rows: record labels, genes, gene labels, and the ribbons between them (no legend).
BGC_ROWS = _preview_part('[data-gbdraw-role="record-definition"], [data-gbdraw-feature-id], '
                         '[data-label-feature-id], [data-match-kind]')


def _bgc_detail(rec: Recorder) -> Callable[[Recorder], dict] | None:
    """The tall window shows the five fitted rows too small to read, so it reads them closely."""

    return BGC_ROWS if rec.portrait else None


def _bgc_aminoglycosides(rec: Recorder, downloads: Path) -> Path:
    _verify(*BGC_FILES)
    page = rec.page

    rec.camera(duration=0)
    rec.caption(1, "Switch to Linear mode", "One row per gene cluster")
    rec.hold(1.6)
    _select_mode(rec, "Linear")

    rec.caption(2, "Add the first cluster", "File 1: upload BGC0000708.gbk")
    card = _linear_file(rec, 1, BGC_SOURCES[0])
    rec.caption(2, "Name the cluster", "Organism / strain · Subtitle / title")
    _bgc_record(rec, card, 1)

    rec.caption(3, "Add the other four clusters", "Add sequence, then upload and name each file")
    for index, source in enumerate(BGC_SOURCES[1:], start=2):
        if index == 3:
            rec.caption(3, "Add the other four clusters", "BGC0000711, 712, and 713, in this order")
            rec.fast_forward(factor=3.0)
        _add_sequence(rec)
        _bgc_record(rec, _linear_file(rec, index, source), index)
    rec.normal_speed()

    rec.caption(4, "Generate with the defaults", "The first run starts Python in the browser")
    _generate(rec, target=1.8, look=0.6)
    _show_figure(rec, 4, "Five clusters, one row each", "Nothing shows yet which genes they share",
                 detail=_bgc_detail(rec))

    rec.caption(5, "Link the clusters with LOSATP", "Comparison › Run LOSAT, then Settings")
    settings = _losat_settings(rec)
    rec.caption(5, "Group similar proteins", "LOSATP mode: Similarity groups")
    _pick(rec, _label(rec, "LOSATP mode", within=settings), "orthogroup")
    rec.caption(5, "Drop the weak matches", "Min identity (%): 30")
    _field(rec, _label(rec, "Linear comparison minimum identity", within=settings), "30")
    rec.caption(5, "Run LOSAT once", "Generate shows each stage while it runs")
    _generate(rec, target=3.0, look=0.6, stages=True)
    if page.evaluate("() => (window.__GBDRAW_APP__?.orthogroups || []).length") != 23:
        raise AssertionError("Expected 23 similarity groups")
    _show_figure(rec, 5, "Ribbons link similar proteins", "The X: BGC0000713 runs the other way",
                 detail=_bgc_detail(rec))

    rec.caption(6, "Reverse-complement File 5", "BGC0000713 › Record options › Reverse complement")
    rec.open_section(_summary(rec, "Record options for sequence {n}", {"n": 5}))
    _check(rec, _label(rec, "Reverse complement for sequence {n}", {"n": 5}))
    rec.caption(6, "LOSAT results are reused", "Generate redraws without a new search")
    _generate(rec, target=None, look=0.6, stages=True)
    _show_figure(rec, 6, "The X is gone", "All five clusters now run the same way", detail=_bgc_detail(rec))

    _bgc_align(rec)

    rec.caption(8, "Color by antiSMASH gene kind", "Colors › Palette: orange · CDS: #dddddd")
    rec.open_section(_summary(rec, "Colors"))
    _pick(rec, _label(rec, "Palette"), "orange")
    cds = _label(rec, "{type} feature color", {"type": "CDS"})
    rec.reveal(cds)
    rec.focus(cds, zoom=2.2, duration=0.6)
    rec.pick_color(cds, "#dddddd")
    rec.caption(8, "Add four SPECIFIC RULES", "gene_kind: core, additional, transport, regulatory")
    _color_rules(rec, BGC_RULES)
    _generate(rec, look=0.6)
    _show_figure(rec, 8, "Red core genes stand out", "against the blue ribbons", detail=_bgc_detail(rec))

    rec.caption(9, "Label the first cluster", "Labels › Show Labels: First Record Only")
    rec.open_section(_summary(rec, "Labels"))
    _pick(rec, page.locator("#linear-show-labels"), "first")
    rec.caption(9, "Use the short gene names", "Qualifier priority: CDS › gene")
    _priority_rule(rec, "CDS", "gene")
    rec.caption(9, "Above the features, at 45°", "Label Placement · Label Rotation")
    _pick(rec, _label(rec, "Label Placement"), "above_feature")
    _field(rec, _label(rec, "Label Rotation"), "45")
    _generate(rec, look=0.6)
    _show_figure(rec, 9, "Gene names on the first row", "The other rows stay uncrowded", detail=_bgc_detail(rec))

    _bgc_finish(rec)
    rec.caption(10, "Generate the finished figure", "It keeps og_6 and reuses the LOSAT results")
    _generate(rec, target=None, look=0.6, stages=True)
    _show_figure(rec, 10, "The finished figure", "", seconds=4.0, detail=_bgc_detail(rec))

    rec.caption(11, "Inspect a similarity ribbon", "og_18: 6 genes in 5 records, and its members")
    ribbon = _hit(rec, '[data-match-kind="orthogroup"][data-orthogroup-id="og_18"]')
    rec.click_point(ribbon["x"], ribbon["y"], pause=0.4)
    details = page.get_by_role("dialog", name=rec.t("Pairwise match details"), exact=True)
    expect(details).to_be_visible()
    expect(details).to_contain_text("og_18")
    rec.camera(details.bounding_box(), zoom=1.5, duration=1.0)
    rec.hold(4.0)
    rec.click(_label(rec, "Close match popup", within=details), pause=0.3)
    rec.camera(duration=0.9)
    rec.hold(0.5)

    rec.caption(12, "Inspect the feature qualifiers", "livE: gene gives the label, gene_kind the color")
    live = "(item) => item.gene === arg && item.orthogroupId === 'og_1'"
    popup = _open_feature(rec, _feature_point(rec, live, "livE"), "Qualifiers")
    rec.hold(1.0)
    _scroll_popup_to(rec, popup, "gene_kind")
    rec.hold(3.2)
    _close_feature(rec)

    live_id = _svg_id(rec, live, "livE")
    exported = export_interactive_svg(rec, downloads, None)
    errors = _open_export(rec, exported, None)
    rec.caption(None, "Click a feature in the file", "livE: its og_1 members light up")
    popup = _svg_popup(rec, _hit(rec, f'[data-gbdraw-feature-id="{live_id}"]'), overview=2.6)
    expect(popup).to_contain_text("og_1")
    rec.hold(3.0)
    _svg_close(rec, popup)
    _svg_done(rec, errors)
    return exported


def _inspect_bgc_aminoglycosides(svg: Path) -> dict:
    labels, fills = _labels(svg), _fills(svg)
    text = svg.read_text(encoding="utf-8")
    report = {
        "rule_fills": sorted({color for *_, color, _ in BGC_RULES} & fills),
        "legend": sorted({legend for *_, legend in BGC_RULES} & labels),
        "default_cds": "#dddddd" in fills,
        "similarity_links": text.count('data-match-kind="orthogroup"'),
        "livE": "livE" in labels,
        "title": "Aminoglycoside biosynthetic gene clusters from " in labels,
        "records": sorted({"BGC0000708", "BGC0000709", "BGC0000711", "BGC0000712", "BGC0000713"}
                          & {label.split(".")[0] for label in labels}),
    }
    return _require(report, len(report["rule_fills"]) == 4 and len(report["legend"]) == 4 and report["default_cds"]
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
    preview = page.locator(PREVIEW_REGION)

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
# File defaults: (Organism / strain, Subtitle / title) for all records of each file.
VIBRIO_SPECIES = (("<i>Vibrio parahaemolyticus</i>", "RIMD 2210633"), ("<i>Vibrio alginolyticus</i>", "NBRC 15630"))
# (qualifier key, query, the record of each match in the order the search lists them)
VIBRIO_INITIATORS = (("gene", "dnaA", ("NC_004603.1", "NC_022349.1")),
                     ("product", "RctB", ("NC_004605.1", "NC_022359.1")))
# Chromosome I of both rows with their species labels: the part of the wide
# figure the detail holds frame, so its labels and blocks read.
VIBRIO_CHROMOSOME_I = _preview_part('[data-record-column="0"], [data-gbdraw-role="record-definition-row"]')
VIBRIO_SHAPES = ("CDS", "rRNA", "tRNA", "tmRNA", "ncRNA", "misc_RNA", "repeat_region")
VIBRIO_GALLERY = REPO_ROOT / "gbdraw/web/gallery/sources/vibrio-harveyi-group-collinear.svg"
_LEFTMOST = """
([recordId, svgId]) => {
  const svg = document.querySelector('$REGION svg');
  const left = (node) => node.getBoundingClientRect().left;
  const own = [...svg.querySelectorAll('path[data-gbdraw-feature-id]')]
    .filter((node) => node.getAttribute('data-gbdraw-record-id') === recordId);
  const target = svg.querySelector(`[data-gbdraw-feature-id="${svgId}"]`);
  if (!own.length || !target) return null;
  const record = Math.min(...own.map(left));
  return {offset: left(target) - record, width: svg.getBoundingClientRect().width};
}
""".replace("$REGION", PREVIEW_REGION)


def _verify_vibrio() -> None:
    for path, size, digest in VIBRIO_FILES:
        assert_fixture_identity(path, expected_size=size, expected_sha256=digest)


def _frame_box(rec: Recorder, box: dict, *, zoom: float = 1.6, pad: float = 30) -> None:
    area = {"x": box["x"] - pad, "y": box["y"] - pad, "width": box["width"] + 2 * pad, "height": box["height"] + 2 * pad}
    rec.camera(area, zoom=zoom, duration=0.9, portrait=area)


def _vibrio_file(rec: Recorder, index: int) -> None:
    """Upload File ``index`` and fill its Names for all records in this file."""

    page = rec.page
    card = _linear_file(rec, index, VIBRIO_SOURCES[index - 1])
    page.wait_for_function("(count) => window.__GBDRAW_APP__?.linearSeqs?.length === count",
                           arg=2 * index, timeout=GENERATION_TIMEOUT_MS)
    species, strain = VIBRIO_SPECIES[index - 1]
    rec.caption(index + 1, "Name both of its records", "Names for all records in this file")
    _field(rec, _label(rec, "Default definition for file {n}", {"n": index}, within=card), species, zoom=2.0)
    _field(rec, _label(rec, "Default subtitle for file {n}", {"n": index}, within=card), strain, zoom=2.0)


def _search_matches(rec: Recorder, current: int, count: int = 2) -> str:
    return rec.tn(count, "{current} of {count} match", "{current} of {count} matches", {"current": current})


def _rotate_open_feature(rec: Recorder, step: int, record_id: str) -> None:
    """In the open feature popup: Rotate record using this feature › Start of the record › Apply when I Generate."""

    page = rec.page
    popup = page.locator(".feature-popup")
    expect(popup).to_be_visible()
    page.wait_for_function("(id) => window.__GBDRAW_APP__?.clickedFeature?.feat?.record_id === id", arg=record_id)
    toggle = popup.locator("#feature-record-actions-toggle")
    toggle.evaluate("node => node.scrollIntoView({block: 'center'})")
    page.wait_for_timeout(300)
    _frame_box(rec, popup.bounding_box())
    rec.hold(0.8)
    rec.click(toggle, pause=0.4)
    rec.caption(step, "Start of the record", "then Apply when I Generate")
    start = popup.locator("input[name='feature-record-rotation-position'][value='start']")
    start.evaluate("node => node.scrollIntoView({block: 'center'})")
    page.wait_for_timeout(300)
    rec.set_checked(start, True, pause=0.6)
    expect(popup.locator("#feature-record-rotation-preview")).to_contain_text(record_id)
    apply = _button(rec, "Apply when I Generate", within=popup)
    apply.evaluate("node => node.scrollIntoView({block: 'end'})")
    page.wait_for_timeout(300)
    _frame_box(rec, popup.bounding_box())
    rec.hold(1.6)
    rec.click(apply, pause=0.8)
    rec.hold(0.6)


def _rotate_initiators(rec: Recorder, step: int, key: str, query: str, records: tuple[str, str]) -> None:
    """Search the preview for ``key`` = ``query``, then start both matches' records at them."""

    page = rec.page
    search = page.locator(".preview-feature-search")
    rec.camera(search.bounding_box(), zoom=2.4, duration=0.9, portrait=search.bounding_box())
    rec.choose(_label(rec, "Field to search", within=search), "qualifier-value")
    rec.type(_label(rec, "Qualifier key", within=search), key)
    rec.type(search.locator("input.gfs-query"), query)
    rec.click(search.locator("button.gfs-button--search"), pause=0.6)
    status = search.locator(".gfs-count")
    expect(status).to_have_text(_search_matches(rec, 1), timeout=ACTION_TIMEOUT_MS)
    rec.hold(1.6)
    rec.caption(step, "Open the first match", "Rotate record using this feature")
    rec.click(search.locator("button.gfs-button--open"), pause=0.6)
    _rotate_open_feature(rec, step, records[0])
    rec.caption(step, "Now the second match", "> moves the open popup to it")
    rec.camera(search.bounding_box(), zoom=2.4, duration=0.8, portrait=search.bounding_box())
    rec.click(_button(rec, "Next match", within=search), pause=0.5)
    expect(status).to_have_text(_search_matches(rec, 2), timeout=ACTION_TIMEOUT_MS)
    rec.hold(0.8)
    _rotate_open_feature(rec, step, records[1])
    _close_feature(rec)


def _vibrio_harveyi_group_collinear(rec: Recorder, downloads: Path) -> Path:
    _verify_vibrio()
    page = rec.page

    rec.camera(duration=0)
    rec.caption(1, "Switch to Linear mode", "Two Vibrio genomes, two chromosomes each")
    rec.hold(1.6)
    _select_mode(rec, "Linear")

    rec.caption(2, "Add the first assembly", "File 1: GCF_000196095.1, two chromosomes")
    _vibrio_file(rec, 1)
    rec.caption(3, "Add the second assembly", "Add sequence, then GCF_000354175.2")
    _add_sequence(rec)
    _vibrio_file(rec, 2)

    rec.caption(4, "Generate with the defaults", "The first run starts Python in the browser")
    _generate(rec, target=2.4, look=0.6)
    _show_figure(rec, 4, "One species per row", "Each row shows only its species and strain",
                 detail=VIBRIO_CHROMOSOME_I)

    rec.caption(5, "Compare with Collinear blocks", "Comparison › Run LOSAT, then Settings")
    settings = _losat_settings(rec)
    rec.caption(5, "LOSATP mode: Collinear blocks", "Color mode: Orientation + identity")
    _pick(rec, _label(rec, "LOSATP mode", within=settings), "collinear")
    _pick(rec, _label(rec, "Collinear color mode", within=settings), "orientation_identity")
    rec.caption(5, "Set the block rules", "Max unit gap 2 · Min block genes 3")
    _field(rec, _label(rec, "Collinear max unit gap", within=settings), "2")
    _field(rec, _label(rec, "Collinear minimum block genes", within=settings), "3")
    rec.caption(5, "Build similarity groups", "with within-record searches: On")
    _check(rec, _label(rec, "Build similarity groups with within-record searches", within=settings))
    rec.caption(5, "Diagonal drift: 2", "Advanced comparison and layout")
    rec.open_section(_summary(rec, "Advanced comparison and layout"))
    _field(rec, _label(rec, "Collinear diagonal drift"), "2")
    rec.caption(5, "Run LOSAT once", "Generate shows each stage while it runs")
    # Four 3-5 Mbp proteomes: the search takes minutes.
    _generate(rec, target=3.0, look=0.6, stages=True, timeout_ms=900_000)
    _show_figure(rec, 5, "The ribbons run diagonally", "NCBI starts each chromosome elsewhere",
                 detail=VIBRIO_CHROMOSOME_I)

    rec.caption(6, "Start chromosome I at dnaA", "Search: Qualifier value · gene · dnaA")
    _rotate_initiators(rec, 6, *VIBRIO_INITIATORS[0])
    rec.caption(7, "Start chromosome II at rctB", "Search: Qualifier value · product · RctB")
    _rotate_initiators(rec, 7, *VIBRIO_INITIATORS[1])
    rec.caption(7, "Clear the search", "The other features stop dimming")
    search = page.locator(".preview-feature-search")
    rec.camera(search.bounding_box(), zoom=2.4, duration=0.8, portrait=search.bounding_box())
    rec.click(search.locator("button.gfs-button--clear"), pause=0.8)

    rec.caption(8, "Generate rotated chromosomes", "Generate reuses the LOSAT results")
    # Redrawing four chromosomes still takes about 20 s, so it is fast-forwarded too.
    _generate(rec, target=4.0, look=0.6, stages=True, timeout_ms=900_000)
    initiators = (("(item) => item.gene === arg && item.record_id === 'NC_004603.1'", "dnaA", "NC_004603.1"),
                  ("(item) => item.gene === arg && item.record_id === 'NC_022349.1'", "dnaA", "NC_022349.1"),
                  ("(item) => /RctB/i.test(item.product || '') && item.record_id === 'NC_004605.1'", "", "NC_004605.1"),
                  ("(item) => /RctB/i.test(item.product || '') && item.record_id === 'NC_022359.1'", "", "NC_022359.1"))
    for feature_filter, arg, record_id in initiators:
        found = page.evaluate(_LEFTMOST, [record_id, _svg_id(rec, feature_filter, arg)])
        if found is None or found["offset"] > 3:
            raise AssertionError(f"{record_id} does not start at its initiator: {found}")
    _show_figure(rec, 8, "Blocks run straight down", "Inversions stand out in red", detail=VIBRIO_CHROMOSOME_I)

    rec.caption(9, "Finish the figure", "Features › Rendering: Rectangle for all seven")
    rec.open_section(_summary(rec, "Features"))
    for index, feature in enumerate(VIBRIO_SHAPES):
        if index == 1:
            rec.fast_forward(factor=3.0)
        _pick(rec, _label(rec, "Rendering for {feature}", {"feature": feature}), "rectangle", zoom=2.0)
    rec.normal_speed()
    rec.caption(9, "Thin the outlines", "Block Stroke Width 0 · Line Stroke Width 1")
    _field(rec, _label(rec, "Block Stroke Width"), "0")
    _field(rec, _label(rec, "Line Stroke Width"), "1")
    rec.caption(9, "Layout › Track Layout: Above", "Separate Strands and Lock Definition Column stay on")
    rec.open_section(_summary(rec, "Layout"))
    _pick(rec, page.locator("#linear-track-layout"), "above")
    for name in ("Separate Strands", "Lock Definition Column"):
        box = _label(rec, name)
        expect(box).to_be_checked()
        _point(rec, box)
    rec.caption(9, "A ruler on each axis", "Axis & Scale › Ruler (Ticks) every 750,000 bp")
    rec.open_section(_summary(rec, "Axis & Scale"))
    _field(rec, _label(rec, "Axis Stroke Width"), "2")
    _pick(rec, _label(rec, "Linear scale style"), "ruler")
    _check(rec, page.locator("input[aria-describedby='help-linear-axis-ruler-on-axis']"))
    _field(rec, _label(rec, "Linear ruler interval"), "750000")
    rec.caption(9, "Set the record-label fonts", "Default 16 · Name / Species 18, bold · Subtitle 16")
    rec.open_section(_summary(rec, "Titles and Record Labels"))
    _field(rec, _label(rec, "Default font size"), "16")
    rec.fast_forward(factor=2.0)
    for line, size, bold in (("Name / Species", "18", True), ("Subtitle", "16", False)):
        label = {"label": rec.t(line)}
        rec.open_section(_summary(rec, "{label} style", label))
        _field(rec, _label(rec, "{label} font size", label), size, zoom=2.0)
        if bold:
            _press(rec, _button(rec, "{label} bold weight", label), zoom=2.0, pause=0.3)
    rec.normal_speed()

    rec.caption(10, "Generate and check the result", "Generate again reuses the LOSAT results")
    _generate(rec, target=4.0, look=0.6, stages=True, timeout_ms=900_000)
    rec.mark("blocks")
    _show_figure(rec, 10, "Rows start at dnaA and rctB", "Each chromosome begins at its initiator",
                 detail=VIBRIO_CHROMOSOME_I)
    _look(rec, rec.text_bounds(r"^(Collinear|Inverted)$"), zoom=2.2, seconds=2.4, pad=40)
    _show_figure(rec, 10, "Red blocks: inversions", "within chromosome I or within chromosome II")
    _show_figure(rec, 10, "Few blocks cross over", "from one species' chromosome I to the other's II")

    dnaa = _svg_id(rec, "(item) => item.gene === arg && item.record_id === 'NC_004603.1'", "dnaA")
    exported = export_interactive_svg(rec, downloads, None)
    errors = _open_export(rec, exported, None)
    rec.caption(None, "Click a feature in the file", "dnaA of Vibrio parahaemolyticus")
    popup = _svg_popup(rec, _hit(rec, f'[data-gbdraw-feature-id="{dnaa}"]'))
    expect(popup).to_contain_text("dnaA")
    rec.hold(3.0)
    _svg_close(rec, popup)
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
        # Auto hides Accession and Length / Coordinates for two records per row; Replicon is off.
        "accession_or_length_drawn": sorted(label for label in labels
                                            if label in VIBRIO_RECORDS or re.fullmatch(r"[\d,]+ bp", label)
                                            or re.fullmatch(r"(chromosome|Chromosome) (1|2|I|II)", label)),
    }
    return _require(report, len(report["records"]) == 4 and len(report["species"]) == 4
                    and len(report["legend"]) == 2 and report["ruler"] and report["blocks"] > 0
                    and report["same_blocks"] and not report["accession_or_length_drawn"], svg)


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
# Each take is one journey: the tall video is recorded again, so the preview
# figure and the opened file can take a tall window's shape, and a take in
# another language records the app in that language; the steps and checks are
# the same.
TUTORIALS = tuple(replace(journey, name=take_name(journey.name, vertical=vertical, language=language),
                          portrait=vertical)
                  for journey in _LANDSCAPE_TAKES for vertical in (False, True)
                  for language in TUTORIAL_LANGUAGES[journey.name])
