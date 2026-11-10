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
  constructor(attributes = {}, { children = [], box = null } = {}) {
    this.attributes = attributes;
    this.children = children;
    this.box = box;
    this.parent = null;
    children.forEach((child) => { child.parent = this; });
  }

  getAttribute(name) { return name in this.attributes ? String(this.attributes[name]) : null; }

  // The selectors the module asks, by what they name.
  matches(selector) {
    if (selector.includes('sequence_conservation')) return this.attributes['data-gbdraw-slot-renderer'] === 'sequence_conservation';
    if (selector.includes('conservation_identity_legend')) return this.attributes.legendRow === true;
    if (selector.startsWith('path[data-display-start]')) return 'data-display-start' in this.attributes;
    if (selector === '[data-gbdraw-record-index]') return 'data-gbdraw-record-index' in this.attributes;
    if (selector.startsWith('g[data-gbdraw-record-index="')) return selector.includes(`"${this.attributes['data-gbdraw-record-index']}"`);
    return false;
  }

  descendants() { return this.children.flatMap((child) => [child, ...child.descendants()]); }

  querySelectorAll(selector) { return this.descendants().filter((element) => element.matches(selector)); }

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
  DIMMED_RING_OPACITY, identityText, ringHitAt, ringHoverStyle, ringHoverTarget, ringPointAt, ringsOfRecord
} = await import('../../gbdraw/web/js/app/feature-editor/ring-hover.js');

/** Coordinates to a millionth of a pixel. @param {Record<string, number>} line */
const rounded = (line) => Object.fromEntries(Object.entries(line).map(([key, value]) => [key, Math.round(value * 1e6) / 1e6 || 0]));

const hit = (start, end, identity, fill = '#8b9cc1') => new FakeElement({
  'data-match-kind': 'homology', 'data-identity': identity, 'data-display-start': start, 'data-display-end': end, fill
});
const ring = (sourceIndex, label, radii, hits) => new FakeElement({
  'data-gbdraw-slot-renderer': 'sequence_conservation', 'data-source-index': sourceIndex, 'data-track-label': label,
  ...(radii ? { 'data-ring-inner-radius': radii[0], 'data-ring-outer-radius': radii[1] } : {})
}, { children: hits });

test('a hover dims every other ring of the series to 35% and bolds its Legend row, in the preview only', () => {
  assert.equal(DIMMED_RING_OPACITY, 0.35);
  const css = ringHoverStyle('2');
  assert.match(css, /^\.gbdraw-preview-surface g\[data-gbdraw-slot-renderer="sequence_conservation"\]:not\(\[data-source-index="2"\]\) \{ opacity: 0\.35; \}$/m);
  assert.match(css, /^\.gbdraw-preview-surface #conservation_identity_legend g\[data-source-index\]\[data-source-index="2"\] text \{ font-weight: 700; \}$/m);
  // Every rule names the preview surface, so no other figure on the page (a Gallery card) changes.
  assert.ok(css.trim().split('\n').every((rule) => rule.startsWith('.gbdraw-preview-surface ')));
});

test('the rings of a record are listed outermost first by the radii Python wrote', () => {
  const record = new FakeElement({ 'data-gbdraw-record-index': 0 }, { children: [
    ring('0', 'inner', [100, 110], []), ring('1', 'outer', [130, 140], []), ring('2', 'middle', [115, 125], [])
  ] });
  assert.deepEqual(ringsOfRecord(record).map((entry) => [entry.name, entry.sourceIndex, entry.inner, entry.outer]), [
    ['outer', '1', 130, 140], ['middle', '2', 115, 125], ['inner', '0', 100, 110]
  ]);
  // A Result drawn before Python wrote the radii keeps the document order and no radius.
  const old = new FakeElement({}, { children: [ring('0', 'a', null, []), ring('1', 'b', null, [])] });
  assert.deepEqual(ringsOfRecord(old).map((entry) => [entry.name, entry.outer]), [['a', null], ['b', null]]);
});

test('a ring reads its identity and drawn colour at a position from its hits\' drawn spans, the highest one on top', () => {
  const drawn = ring('0', 'A', [100, 110], [hit(10, 60, 95, '#a0b0c0'), hit(40, 80, 99.5, '#8b9cc1'), hit(100, 120, 80, '#e0e4ee')]);
  assert.deepEqual(ringHitAt(drawn, 10), { identity: 95, color: '#a0b0c0' });
  assert.deepEqual(ringHitAt(drawn, 50.5), { identity: 99.5, color: '#8b9cc1' });
  assert.equal(ringHitAt(drawn, 79.99).identity, 99.5);
  // The end is the base after the hit's last one.
  assert.equal(ringHitAt(drawn, 80), null);
  assert.deepEqual(ringHitAt(drawn, 119), { identity: 80, color: '#e0e4ee' });
  assert.equal(identityText(99.5), '99.5%');
  assert.equal(identityText(100), '100.0%');
  assert.equal(identityText(null), 'No hit');
  // A hit without the spans (a Result drawn before Python wrote them) gives no identity.
  const old = ring('0', 'A', null, [new FakeElement({ 'data-match-kind': 'homology', 'data-identity': 99 })]);
  assert.equal(ringHitAt(old, 10), null);
});

test('the position is the source base under the pointer (R5), and the line runs across the rings and the features', () => {
  const axis = axisOf('circular-display-start');
  const rings = [ring('1', 'outer', [300, 310], []), ring('0', 'inner', [280, 290], [])];
  const frame = new FakeElement({ 'data-gbdraw-record-index': 0 }, { children: rings });
  const root = new FakeElement({}, { children: [frame] });
  const host = { left: 10, top: 20 };
  const at = (clientX, clientY) => ringPointAt({
    root, ring: rings[0], rings: ringsOfRecord(frame), axisRecords: [{ recordIndex: 0, axis }], clientX, clientY, host
  });
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
  // No axis for the record, or rings without radii: no position, so no line and no popup.
  assert.equal(ringPointAt({ root, ring: rings[0], rings: ringsOfRecord(frame), axisRecords: [], clientX: 0, clientY: -305, host }), null);
  const old = [ring('0', 'a', null, [])];
  const oldRoot = new FakeElement({}, { children: [new FakeElement({ 'data-gbdraw-record-index': 0 }, { children: old })] });
  assert.equal(ringPointAt({ root: oldRoot, ring: old[0], rings: ringsOfRecord(oldRoot), axisRecords: [{ recordIndex: 0, axis }], clientX: 0, clientY: -305, host }), null);
});

test('each record of the Multi-Record Canvas reads its own rings and axis; the series is the same in every record', () => {
  const axis = axisOf('circular-multi-record-canvas');
  const first = [ring('0', 'A', [150, 160], [hit(0, 100, 91)])];
  const second = [ring('0', 'A', [150, 160], [hit(0, 100, 97)])];
  const records = [first, second].map((rings, index) => new FakeElement({ 'data-gbdraw-record-index': index }, { children: rings }));
  const root = new FakeElement({}, { children: records });
  const point = ringPointAt({
    root, ring: second[0], rings: ringsOfRecord(records[1]), axisRecords: [{ recordIndex: 0, axis: null }, { recordIndex: 1, axis }],
    clientX: 0, clientY: -155, host: { left: 0, top: 0 }
  });
  assert.equal(point.base, 1);
  assert.equal(ringHitAt(second[0], point.position).identity, 97);
  assert.equal(ringHoverTarget(root, second[0].children[0], 0, 0).sourceIndex, first[0].getAttribute('data-source-index'));
});

test('the pointer is over a ring by its target, and over a Legend name by the row\'s box', () => {
  const drawnRing = ring('3', 'A', [100, 110], [hit(0, 10, 90)]);
  const row = new FakeElement({ legendRow: true, 'data-source-index': '3' }, { box: { left: 500, right: 600, top: 40, bottom: 60 } });
  const root = new FakeElement({}, { children: [drawnRing, new FakeElement({ id: 'conservation_identity_legend' }, { children: [row] })] });
  assert.deepEqual(ringHoverTarget(root, drawnRing.children[0], 0, 0), { kind: 'ring', sourceIndex: '3', ring: drawnRing });
  // Between the glyphs of a name the target is the root; the row's box still holds the point.
  assert.deepEqual(ringHoverTarget(root, root, 550, 50), { kind: 'legend', sourceIndex: '3' });
  assert.equal(ringHoverTarget(root, root, 450, 50), null);
  // A ring of another figure is not this Result's.
  assert.equal(ringHoverTarget(root, ring('3', 'B', [1, 2], []), 0, 0), null);
});
