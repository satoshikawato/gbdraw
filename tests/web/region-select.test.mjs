// Select region's owner (gbdraw/web/js/app/region-select.js, PD-OI-115): the toggle and its exclusivity with Layout
// edit, the reason it is unavailable, and the drag gesture (pointer capture, cancel, a record that is not circular).
// The axis comes from the same vectors as tests/test_region_axis_geometry.py.
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import test from 'node:test';

const VECTORS = JSON.parse(readFileSync(new URL('../fixtures/region_axis_geometry_vectors.json', import.meta.url), 'utf8')).cases;
const axisOf = (name) => VECTORS.find((entry) => entry.name === name).records[0].axis;

// The page the owner reads: an identity screen matrix for the record's frame, frames at once.
class FakeElement {
  constructor() {
    this.captured = [];
    this.released = [];
  }

  setPointerCapture(pointerId) { this.captured.push(pointerId); }

  releasePointerCapture(pointerId) { this.released.push(pointerId); }
}
const identity = { a: 1, b: 0, c: 0, d: 1, e: 0, f: 0, inverse() { return this; } };
Object.assign(globalThis, {
  Element: FakeElement,
  DOMPoint: class {
    constructor(x, y) { Object.assign(this, { x, y }); }

    matrixTransform(m) { return { x: m.a * this.x + m.c * this.y + m.e, y: m.b * this.x + m.d * this.y + m.f }; }
  },
  window: { innerHeight: 1000, requestAnimationFrame: (callback) => { callback(); return 0; } },
  document: { querySelector: () => null, querySelectorAll: () => [], addEventListener: () => {}, removeEventListener: () => {} }
});

const { createRegionSelect } = await import('../../gbdraw/web/js/app/region-select.js');

const frame = { contains: () => false, getScreenCTM: () => identity };
const root = { querySelectorAll: (selector) => (selector.startsWith('g[data-gbdraw-record-index') ? [frame] : []) };

/** An owner over one record of `axis`; `records` empties when the shown Result has no axis. */
const setup = ({ axis = axisOf('circular-display-start'), withAxis = true } = {}) => {
  const layoutEdit = { value: false };
  const opened = [];
  let open = null;
  const owner = createRegionSelect({
    readAxisRecords: () => (withAxis ? [{ resultIndex: 0, recordIndex: 0, recordKey: 'rec', recordId: 'REC1', axis }] : []),
    shownRoot: () => /** @type {any} */ (root),
    layoutEdit,
    openSelection: () => open,
    openRegionPopup: (selection, _element, anchor) => {
      opened.push({ selection, anchor });
      open = selection;
      return true;
    },
    closeRegionPopup: () => { open = null; },
    drawnFeatures: () => [],
    selectFeatures: () => {},
    addRegionAnnotation: () => null,
    newDrawingFromRegion: async () => null,
    readRegionFasta: async () => ''
  });
  owner.bindBand({});
  return { owner, layoutEdit, opened, canvas: new FakeElement() };
};

/** The point on the middle of the feature band of display position `position`. */
const pointAt = (axis, position) => {
  const radius = (axis.bandInnerPx + axis.bandOuterPx) / 2;
  const angle = (axis.startAngleDeg + axis.degreesPerBp * position) * Math.PI / 180;
  return { clientX: radius * Math.cos(angle), clientY: radius * Math.sin(angle) };
};
const press = (canvas, point, extra = {}) => ({
  pointerId: 7, button: 0, isPrimary: true, currentTarget: canvas, preventDefault() {}, stopPropagation() {}, ...point, ...extra
});

test('Select region and Layout edit turn each other off; without an axis the toggle is off with its reason', () => {
  const { owner, layoutEdit } = setup();
  assert.equal(owner.on.value, false, 'off at first');
  assert.equal(owner.unavailableReason.value, '');
  layoutEdit.value = true;
  owner.toggle();
  assert.deepEqual([owner.active.value, layoutEdit.value], [true, false]);
  owner.toggleLayoutEdit();
  assert.deepEqual([owner.active.value, layoutEdit.value], [false, true]);
  owner.toggleLayoutEdit();
  assert.deepEqual([owner.active.value, layoutEdit.value], [false, false], 'turning Layout edit off leaves Select region off');

  // A Result saved before Python wrote the axis (Q-A1).
  const old = setup({ withAxis: false });
  assert.equal(old.owner.available.value, false);
  assert.equal(old.owner.unavailableReason.value, 'Generate to enable region selection');
  old.owner.toggle();
  assert.equal(old.owner.active.value, false);
});

test('a drag captures the pointer once it moves and releases it on release; a press that stays put captures nothing', () => {
  const axis = axisOf('circular-display-start');
  const { owner, opened, canvas } = setup({ axis });
  owner.toggle();
  owner.pointerDown(press(canvas, pointAt(axis, 4_000.5)));
  owner.pointerUp(press(canvas, pointAt(axis, 4_000.5)));
  assert.deepEqual([canvas.captured, opened.length], [[], 0], 'a click is left to the click');

  owner.pointerDown(press(canvas, pointAt(axis, 4_000.5)));
  owner.pointerMove(press(canvas, pointAt(axis, 5_000.5)));
  assert.deepEqual(canvas.captured, [7]);
  owner.pointerMove(press(canvas, pointAt(axis, 6_000.5)));
  owner.pointerUp(press(canvas, pointAt(axis, 6_000.5)));
  assert.deepEqual(canvas.released, [7]);
  assert.equal(opened.length, 1);
  const [{ selection, anchor }] = opened;
  assert.deepEqual(selection.span, { start: 4_000, length: 2_001 });
  // The popup's limits come with the anchor (no toolbar on this page).
  assert.deepEqual(anchor.obstacles, []);
  assert.equal(anchor.within, null);
});

test('a lost capture or a cancelled pointer ends the drag with no selection', () => {
  const axis = axisOf('circular-display-start');
  for (const end of ['lostpointercapture', 'pointercancel']) {
    const { owner, opened, canvas } = setup({ axis });
    owner.toggle();
    owner.pointerDown(press(canvas, pointAt(axis, 100.5)));
    owner.pointerMove(press(canvas, pointAt(axis, 900.5)));
    assert.ok(owner.overlay.value, `${end}: the band shows while dragging`);
    // A touch press is captured by the element under it; its capture moving to the canvas is no loss (OV-2483).
    owner.pointerCancel(press(canvas, {}, { type: 'lostpointercapture', target: new FakeElement() }));
    assert.ok(owner.overlay.value, `${end}: the implicit touch capture moving to the canvas keeps the drag`);
    owner.pointerCancel(press(canvas, {}, { type: end, target: canvas }));
    assert.equal(owner.overlay.value, null, `${end}: the band goes`);
    // The release that follows opens nothing.
    owner.pointerUp(press(canvas, pointAt(axis, 900.5)));
    assert.equal(opened.length, 0, end);
    assert.deepEqual(canvas.released, [7], `${end}: the capture is released`);
  }
});

test('a drag on a record that is not circular stops at its last base instead of crossing the origin', () => {
  const axis = axisOf('circular-linear-topology');
  const length = axis.source.drawnLength;
  const { owner, opened, canvas } = setup({ axis });
  owner.toggle();
  // From about 11 o'clock clockwise through 12 o'clock to about 1 o'clock.
  owner.pointerDown(press(canvas, pointAt(axis, length - 3_000.5)));
  for (const position of [length - 1_500.5, length - 10.5, 1_000.5, 3_000.5]) owner.pointerMove(press(canvas, pointAt(axis, position)));
  owner.pointerUp(press(canvas, pointAt(axis, 3_000.5)));
  const [{ selection }] = opened;
  assert.deepEqual(selection.source, { start: length - 3_000, end: length, length: 3_001, wrapsOrigin: false });
  assert.equal(owner.popup.value.newDrawingReason, '');
  assert.equal(owner.popup.value.featureText, '0 features');
});

test('a span across the origin of a circular record says so and disables New drawing from region with its reason', () => {
  const axis = axisOf('circular-display-start');
  const { owner, opened, canvas } = setup({ axis });
  owner.toggle();
  // Display start 16,000: display index 569 shows base 16,569, 570 base 1.
  owner.pointerDown(press(canvas, pointAt(axis, 500.5)));
  owner.pointerMove(press(canvas, pointAt(axis, 560.5)));
  owner.pointerMove(press(canvas, pointAt(axis, 700.5)));
  owner.pointerUp(press(canvas, pointAt(axis, 700.5)));
  assert.equal(opened[0].selection.source.wrapsOrigin, true);
  assert.equal(owner.popup.value.featureText, '0 features · crosses the origin');
  assert.equal(owner.popup.value.newDrawingReason, 'A region that crosses the origin cannot be drawn on its own.');
});
