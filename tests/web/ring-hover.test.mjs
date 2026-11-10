// The similarity-ring hover (gbdraw/web/js/app/feature-editor/ring-hover.js, PD-OI-110 scenario revision 3): the rules
// that dim the other rings to 35% and bold the hovered series' Legend row, the rings of a record outermost first, every
// ring's identity at a drawn position from the spans Python wrote, the position in source coordinates (R5) with the
// radial line's ends, and what a pointer is over. The axis comes from the vectors of tests/test_region_axis_geometry.py.
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import test from 'node:test';

const VECTORS = JSON.parse(readFileSync(new URL('../fixtures/region_axis_geometry_vectors.json', import.meta.url), 'utf8')).cases;
const axisOf = (name) => VECTORS.find((entry) => entry.name === name).records[0].axis;

class FakeElement {
  constructor(attributes = {}, { children = [], box = null, tag = 'g', text = '' } = {}) {
    this.attributes = attributes;
    this.children = children;
    this.box = box;
    this.tag = tag;
    this.textContent = text;
    this.parent = null;
    children.forEach((child) => { child.parent = this; });
  }

  getAttribute(name) { return name in this.attributes ? String(this.attributes[name]) : null; }

  // The selectors the module asks, by what they name.
  matches(selector) {
    if (selector === 'text') return this.tag === 'text';
    if (selector.includes('sequence_conservation')) return this.attributes['data-gbdraw-slot-renderer'] === 'sequence_conservation';
    if (selector.includes('conservation_identity_legend')) return this.attributes.legendRow === true;
    if (selector.startsWith('path[data-display-start]')) return 'data-display-start' in this.attributes;
    if (selector === '[data-gbdraw-record-index]') return 'data-gbdraw-record-index' in this.attributes;
    if (selector.startsWith('g[data-gbdraw-record-index="')) return selector.includes(`"${this.attributes['data-gbdraw-record-index']}"`);
    return false;
  }

  descendants() { return this.children.flatMap((child) => [child, ...child.descendants()]); }

  querySelectorAll(selector) { return this.descendants().filter((element) => element.matches(selector)); }

  querySelector(selector) { return this.querySelectorAll(selector)[0] ?? null; }

  closest(selector) {
    for (let element = this; element; element = element.parent) if (element.matches(selector)) return element;
    return null;
  }

  contains(other) { return other === this || this.descendants().includes(other); }

  getBoundingClientRect() { return this.box; }

  getScreenCTM() { return { a: 1, b: 0, c: 0, d: 1, e: 0, f: 0, inverse() { return this; } }; }
}
Object.assign(globalThis, {
  Element: FakeElement,
  CSS: { escape: (value) => String(value).replace(/["\\]/g, '\\$&') },
  DOMPoint: class {
    constructor(x, y) { Object.assign(this, { x, y }); }

    matrixTransform(m) { return { x: m.a * this.x + m.c * this.y + m.e, y: m.b * this.x + m.d * this.y + m.f }; }
  }
});

const {
  DIMMED_RING_OPACITY, identityText, ringHitAt, ringHoverStyle, ringHoverTarget, ringPointAt, ringPopupShows,
  ringPopupSource, ringsOfRecord
} = await import('../../gbdraw/web/js/app/feature-editor/ring-hover.js');

/** Coordinates to a millionth of a pixel. @param {Record<string, number>} line */
const rounded = (line) => Object.fromEntries(Object.entries(line).map(([key, value]) => [key, Math.round(value * 1e6) / 1e6 || 0]));

const hit = (start, end, identity, fill = '#8b9cc1', opacity = null) => new FakeElement({
  'data-match-kind': 'homology', 'data-identity': identity, 'data-display-start': start, 'data-display-end': end, fill,
  ...(opacity === null ? {} : { 'fill-opacity': opacity })
}, { tag: 'path' });
const ring = (sourceIndex, label, radii, hits, trackIndex = Number(sourceIndex) + 1) => new FakeElement({
  'data-gbdraw-slot-renderer': 'sequence_conservation', 'data-source-index': sourceIndex, 'data-track-index': trackIndex,
  'data-track-label': label, ...(radii ? { 'data-ring-inner-radius': radii[0], 'data-ring-outer-radius': radii[1] } : {})
}, { children: hits });
/** A ring's Legend row: its drawn name (a `<text>` with that box) after its identity ramp. */
const legendRow = (sourceIndex, drawnName, box = null) => new FakeElement({ legendRow: true, 'data-source-index': sourceIndex }, {
  children: [new FakeElement({}, { tag: 'path' }), new FakeElement({}, { tag: 'text', text: drawnName, box })]
});
const legend = (...rows) => new FakeElement({ id: 'conservation_identity_legend' }, { children: rows });
/** A Result: its record groups and, unless null, its ring Legend. */
const result = (records, rows = null) => new FakeElement({}, { children: [...records, ...(rows ? [legend(...rows)] : [])] });
const recordGroup = (index, rings) => new FakeElement({ 'data-gbdraw-record-index': index }, { children: rings });
/** The client point of a display position on a ring of radius `radius` (frame and client coordinates are one here). */
const pointOn = (axis, position, radius) => {
  const angle = (axis.startAngleDeg + axis.degreesPerBp * position) * Math.PI / 180;
  return { clientX: radius * Math.cos(angle), clientY: radius * Math.sin(angle) };
};

test('a hover dims every other ring of the series to 35% and bolds its Legend row, on screen in the preview only', () => {
  assert.equal(DIMMED_RING_OPACITY, 0.35);
  const css = ringHoverStyle('2');
  assert.match(css, /^\.gbdraw-preview-surface g\[data-gbdraw-slot-renderer="sequence_conservation"\]:not\(\[data-source-index="2"\]\) \{ opacity: 0\.35; \}$/m);
  assert.match(css, /^\.gbdraw-preview-surface #conservation_identity_legend g\[data-source-index\]\[data-source-index="2"\] text \{ font-weight: 700; \}$/m);
  // Every rule names the preview surface, so no other figure on the page (a Gallery card) changes, and none prints.
  const lines = css.trim().split('\n');
  assert.equal(lines[0], '@media screen {');
  assert.equal(lines.at(-1), '}');
  assert.ok(lines.slice(1, -1).every((rule) => rule.startsWith('.gbdraw-preview-surface ')));
});

test('the rings of a record are listed outermost first by the radii Python wrote, named as the Legend draws them', () => {
  const rings = [ring('0', 'WSSV', [100, 110], []), ring('1', 'WSSV', [130, 140], []), ring('2', 'orig', [115, 125], [])];
  const record = recordGroup(0, rings);
  // Two rings of one label are numbered by the Legend, and a renamed row shows its new name (<i> runs kept).
  const root = result([record], [legendRow('1', 'WSSV (1)'), legendRow('2', 'Renamed <i>E. coli</i>'), legendRow('0', 'WSSV (2)')]);
  root.descendants().filter((element) => element.textContent === 'Renamed <i>E. coli</i>').forEach((label) => {
    label.textContent = '';
    label.children = [new FakeElement({}, { tag: 'tspan', text: 'Renamed ' }), new FakeElement({ 'font-style': 'italic' }, { tag: 'tspan', text: 'E. coli' })];
  });
  assert.deepEqual(ringsOfRecord(root, record).map((entry) => [entry.name, entry.sourceIndex, entry.inner, entry.outer]), [
    ['WSSV (1)', '1', 130, 140], ['Renamed <i>E. coli</i>', '2', 115, 125], ['WSSV (2)', '0', 100, 110]
  ]);
  // Without a Legend for the rings, the Legend's naming rule over their labels, in the Legend's order (outer first).
  const outerFirst = [ring('1', 'WSSV', [130, 140], []), ring('0', 'WSSV', [100, 110], [])];
  const bare = recordGroup(0, outerFirst);
  assert.deepEqual(ringsOfRecord(result([bare]), bare).map((entry) => entry.name), ['WSSV (1)', 'WSSV (2)']);
  // A Result drawn before Python wrote the radii has no rings to list.
  const old = recordGroup(0, [ring('0', 'a', null, []), ring('1', 'b', null, [])]);
  assert.equal(ringsOfRecord(result([old]), old), null);
});

test('a ring reads its identity and drawn colour at a position from its hits\' drawn spans, the last drawn on top', () => {
  const drawn = ring('0', 'A', [100, 110], [hit(10, 60, 95, '#a0b0c0'), hit(40, 80, 99.5, '#8b9cc1'), hit(100, 120, 80, '#e0e4ee')]);
  assert.deepEqual(ringHitAt(drawn, 10), { identity: 95, color: '#a0b0c0' });
  assert.deepEqual(ringHitAt(drawn, 50.5), { identity: 99.5, color: '#8b9cc1' });
  assert.equal(ringHitAt(drawn, 79.99).identity, 99.5);
  // The end is the base after the hit's last one.
  assert.equal(ringHitAt(drawn, 80), null);
  assert.deepEqual(ringHitAt(drawn, 119), { identity: 80, color: '#e0e4ee' });
  // Python paints the later path over the earlier: the later one is shown whatever its identity.
  const painted = ring('0', 'A', [100, 110], [hit(0, 100, 99, '#000000'), hit(20, 30, 88, '#808080')]);
  assert.equal(ringHitAt(painted, 25).identity, 88);
  // The swatch shows the fill at the opacity the hit is drawn with.
  assert.equal(ringHitAt(ring('0', 'A', [1, 2], [hit(0, 10, 90, '#ff8000', 0.5)]), 5).color, 'rgba(255, 128, 0, 0.5)');
  assert.equal(ringHitAt(ring('0', 'A', [1, 2], [hit(0, 10, 90, '#ff8000', 1)]), 5).color, '#ff8000');
  // A hit without the spans (a Result drawn before Python wrote them) gives no identity.
  const old = ring('0', 'A', null, [new FakeElement({ 'data-match-kind': 'homology', 'data-identity': 99 })]);
  assert.equal(ringHitAt(old, 10), null);
});

test('an identity shows one decimal, and a value below 100 never shows as 100.0', () => {
  assert.equal(identityText(99.5), '99.5%');
  assert.equal(identityText(87.04), '87.0%');
  assert.equal(identityText(87.06), '87.1%');
  assert.equal(identityText(99.96), '99.9%');
  assert.equal(identityText(100), '100.0%');
  assert.equal(identityText(null), 'No hit');
});

test('the position is the source base under the pointer (R5), and the line runs across the rings and the features', () => {
  const axis = axisOf('circular-display-start');
  const rings = [ring('1', 'outer', [300, 310], []), ring('0', 'inner', [280, 290], [])];
  const root = result([recordGroup(0, rings)]);
  const host = { left: 10, top: 20 };
  const source = ringPopupSource(root, rings[0], [{ recordIndex: 0, axis }]);
  const at = (clientX, clientY) => ringPointAt({ source, clientX, clientY, host });
  // Straight up is display index 0, which shows source base 16000 (the drawing starts there).
  const top = at(0, -305);
  assert.equal(top.base, 16000);
  assert.equal(top.position, 0);
  assert.deepEqual(rounded(top.line), { x1: -10, y1: -300, x2: -10, y2: -406.1 });
  // A quarter turn clockwise: display index 4142, source base 3573 across the origin.
  const right = at(305, 0);
  assert.equal(Math.floor(right.position), 4142);
  assert.equal(right.base, 3573);
  assert.deepEqual(rounded(right.line), { x1: 270, y1: -20, x2: 376.1, y2: -20 });
});

test('on a record drawn from another start, the popup names the exact source base and every ring\'s identity there', () => {
  const axis = axisOf('circular-display-start');
  // Python's fragments of a hit across the drawn origin (display 16560..16569 and 0..10), and hits further on.
  const outer = ring('1', 'outer', [300, 310], [hit(16560, 16569, 92), hit(0, 10, 92), hit(50, 200, 97.25), hit(90, 120, 88)]);
  const inner = ring('0', 'inner', [280, 290], [hit(16560, 16569, 99.04), hit(0, 10, 99.04)]);
  const root = result([recordGroup(0, [inner, outer])], [legendRow('1', 'outer'), legendRow('0', 'inner')]);
  const source = ringPopupSource(root, inner.children[0], [{ recordIndex: 0, axis }]);
  const popupAt = (position) => {
    const point = ringPointAt({ source, ...pointOn(axis, position, 285), host: { left: 0, top: 0 } });
    return { base: point.base, rings: source.rings.map((entry) => [entry.name, identityText(ringHitAt(entry.element, point.position)?.identity ?? null)]) };
  };
  // Display 100 is source base 16100; display 16568 is 15999 and display 5 is 16005, either side of the origin.
  assert.deepEqual(popupAt(100.5), { base: 16100, rings: [['outer', '88.0%'], ['inner', 'No hit']] });
  assert.deepEqual(popupAt(60.5), { base: 16060, rings: [['outer', '97.3%'], ['inner', 'No hit']] });
  assert.deepEqual(popupAt(16568.5), { base: 15999, rings: [['outer', '92.0%'], ['inner', '99.0%']] });
  assert.deepEqual(popupAt(5.5), { base: 16005, rings: [['outer', '92.0%'], ['inner', '99.0%']] });
});

test('the ring popup shows, and replaces the match hover summary, only where the Result has the axis and the radii', () => {
  const axis = axisOf('circular-display-start');
  const drawn = ring('0', 'A', [300, 310], [hit(0, 10, 90)]);
  const root = result([recordGroup(0, [drawn])]);
  assert.equal(ringPopupShows(root, drawn.children[0], [{ recordIndex: 0, axis }]), true);
  // A Result loaded from a Session the gbdraw command saved has no axis: the hit keeps its match hover summary.
  assert.equal(ringPopupShows(root, drawn.children[0], []), false);
  // A Result drawn by gbdraw 0.14.5 or earlier has no radii: the same.
  const old = ring('0', 'A', null, [new FakeElement({ 'data-match-kind': 'homology', 'data-identity': 90 })]);
  assert.equal(ringPopupShows(result([recordGroup(0, [old])]), old.children[0], [{ recordIndex: 0, axis }]), false);
  // An element of no ring, or of a ring of another figure, shows no ring popup.
  assert.equal(ringPopupShows(root, root, [{ recordIndex: 0, axis }]), false);
  assert.equal(ringPopupShows(root, ring('0', 'B', [1, 2], []), [{ recordIndex: 0, axis }]), false);
});

test('each record of the Multi-Record Canvas reads its own rings and axis; the series is the same in every record', () => {
  const axis = axisOf('circular-multi-record-canvas');
  const first = [ring('0', 'A', [150, 160], [hit(0, 100, 91)])];
  const second = [ring('0', 'A', [150, 160], [hit(0, 100, 97)])];
  const records = [first, second].map((rings, index) => recordGroup(index, rings));
  const root = result(records);
  const source = ringPopupSource(root, second[0], [{ recordIndex: 0, axis: null }, { recordIndex: 1, axis }]);
  const point = ringPointAt({ source, clientX: 0, clientY: -155, host: { left: 0, top: 0 } });
  assert.equal(point.base, 1);
  assert.equal(ringHitAt(second[0], point.position).identity, 97);
  assert.equal(ringHoverTarget(root, second[0].children[0], 0, 0).sourceIndex, first[0].getAttribute('data-source-index'));
});

test('the pointer is over a ring by its target, and over a Legend name by the name\'s box, not its identity ramp', () => {
  const drawnRing = ring('3', 'A', [100, 110], [hit(0, 10, 90)]);
  const root = result([drawnRing], [legendRow('3', 'A', { left: 540, right: 600, top: 40, bottom: 60 })]);
  assert.deepEqual(ringHoverTarget(root, drawnRing.children[0], 0, 0), { kind: 'ring', sourceIndex: '3', ring: drawnRing });
  // Between the glyphs of a name the target is the root; the name's box still holds the point.
  assert.deepEqual(ringHoverTarget(root, root, 550, 50), { kind: 'legend', sourceIndex: '3' });
  // The row's ramp, left of the name, is not the name.
  assert.equal(ringHoverTarget(root, root, 520, 50), null);
  assert.equal(ringHoverTarget(root, root, 450, 50), null);
  // A ring of another figure is not this Result's.
  assert.equal(ringHoverTarget(root, ring('3', 'B', [1, 2], []), 0, 0), null);
});
