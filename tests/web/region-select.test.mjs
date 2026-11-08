// Select region's owner (gbdraw/web/js/app/region-select.js, PD-OI-115): the toggle and its exclusivity with Layout
// edit, the reason it is unavailable, and the drag gesture (pointer capture, cancel, a record that is not circular).
// The axis comes from the same vectors as tests/test_region_axis_geometry.py. The actions that wait on Python (one at
// a time per selection, and a failure's text) are here too.
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import test from 'node:test';

const VECTORS = JSON.parse(readFileSync(new URL('../fixtures/region_axis_geometry_vectors.json', import.meta.url), 'utf8')).cases;
const axisOf = (name) => VECTORS.find((entry) => entry.name === name).records[0].axis;

// The page the owner reads: an identity screen matrix for the record's frame, frames at once, and the window's
// listeners (Esc and blur during a drag).
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
  window: {
    innerHeight: 1000,
    requestAnimationFrame: (callback) => { callback(); return 0; },
    listeners: new Map(),
    addEventListener(type, listener) { this.listeners.set(type, [...(this.listeners.get(type) ?? []), listener]); },
    removeEventListener(type, listener) { this.listeners.set(type, (this.listeners.get(type) ?? []).filter((l) => l !== listener)); },
    dispatch(type, event = {}) { for (const listener of this.listeners.get(type) ?? []) listener(event); }
  },
  document: { querySelector: () => null, querySelectorAll: () => [], addEventListener: () => {}, removeEventListener: () => {} }
});

const { createRegionSelect } = await import('../../gbdraw/web/js/app/region-select.js');

const frame = { contains: () => false, getScreenCTM: () => identity, isConnected: true };
const root = { querySelectorAll: (selector) => (selector.startsWith('g[data-gbdraw-record-index') ? [frame] : []) };

/** An owner over one record of `axis`; `records` empties when the shown Result has no axis. */
const setup = ({ axis = axisOf('circular-display-start'), withAxis = true, ports = {} } = {}) => {
  const layoutEdit = { value: false };
  const opened = [];
  const pans = [];
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
    readRegionFasta: async () => '',
    panStart: (event) => pans.push(['start', event.pointerId]),
    panMove: (event) => pans.push(['move', event.pointerId]),
    panEnd: (event) => pans.push(['end', event.pointerId]),
    ...ports
  });
  owner.bindBand({});
  return { owner, layoutEdit, opened, pans, canvas: new FakeElement() };
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
  assert.equal('newDrawingReason' in owner.popup.value, false);
  assert.equal(owner.popup.value.featureText, '0 features');
});

test('a span across the origin of a circular record says so and leaves New drawing from region to draw it (D-34)', () => {
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
  assert.equal('newDrawingReason' in owner.popup.value, false, 'no reason: Python draws it as one window across the origin');
});

test('Esc, leaving the window, turning the tool off, and another Result shown cancel a drag: no popup, capture released', () => {
  const axis = axisOf('circular-display-start');
  const cancels = {
    escape: () => {
      let stopped = false;
      window.dispatch('keydown', { key: 'Escape', preventDefault() {}, stopPropagation() { stopped = true; } });
      assert.ok(stopped, 'the Esc of a drag reaches no other handler');
    },
    blur: () => window.dispatch('blur'),
    'tool off': (owner) => owner.toggle(),
    'Layout edit on': (owner) => owner.toggleLayoutEdit(),
    // A drawing switch (or a Result switch, Generate, Load) mounts another Result: the record's frame leaves the page.
    'Result replaced': () => { frame.isConnected = false; }
  };
  for (const [name, cancel] of Object.entries(cancels)) {
    const { owner, opened, canvas } = setup({ axis });
    owner.toggle();
    owner.pointerDown(press(canvas, pointAt(axis, 100.5)));
    owner.pointerMove(press(canvas, pointAt(axis, 900.5)));
    assert.ok(owner.overlay.value, `${name}: the band shows while dragging`);
    cancel(owner);
    owner.pointerMove(press(canvas, pointAt(axis, 1_200.5)));
    owner.pointerUp(press(canvas, pointAt(axis, 1_200.5)));
    frame.isConnected = true;
    assert.equal(opened.length, 0, `${name}: no popup`);
    assert.equal(owner.overlay.value, null, `${name}: the band goes`);
    assert.deepEqual(canvas.released, [7], `${name}: the capture is released`);
    assert.deepEqual([window.listeners.get('keydown'), window.listeners.get('blur')], [[], []], `${name}: no listener stays`);
  }
});

test('with the tool on, a touch outside the tracks pans through the pan owner; a mouse press there and a touch on them do not', () => {
  const axis = axisOf('circular-display-start');
  const { owner, opened, pans, canvas } = setup({ axis });
  owner.toggle();
  const outside = { clientX: 0, clientY: 0 };
  const touch = { pointerType: 'touch', target: canvas };
  owner.pointerDown(press(canvas, outside, touch));
  owner.pointerMove(press(canvas, { clientX: 40, clientY: 30 }, touch));
  owner.pointerUp(press(canvas, { clientX: 40, clientY: 30 }, touch));
  assert.deepEqual(pans, [['start', 7], ['move', 7], ['end', 7]]);
  assert.equal(opened.length, 0);

  pans.length = 0;
  owner.pointerDown(press(canvas, outside, { pointerType: 'mouse', target: canvas }));
  owner.pointerMove(press(canvas, { clientX: 40, clientY: 30 }, { pointerType: 'mouse', target: canvas }));
  owner.pointerUp(press(canvas, { clientX: 40, clientY: 30 }, { pointerType: 'mouse', target: canvas }));
  assert.deepEqual(pans, [], 'a mouse press pans by its own mouse events');

  owner.pointerDown(press(canvas, pointAt(axis, 100.5), touch));
  owner.pointerMove(press(canvas, pointAt(axis, 900.5), touch));
  owner.pointerUp(press(canvas, pointAt(axis, 900.5), touch));
  assert.deepEqual(pans, [], 'a touch on the tracks selects');
  assert.equal(opened.length, 1);
});

/** A promise and its settle functions. */
const deferred = () => {
  let resolve;
  let reject;
  const promise = new Promise((done, fail) => { resolve = done; reject = fail; });
  return { promise, resolve, reject };
};
/** An owner whose popup is open on a span of `axis`; `ports` replaces the action ports. */
const openPopup = (ports) => {
  const axis = axisOf('circular-display-start');
  const made = setup({ axis, ports });
  made.owner.toggle();
  made.owner.pointerDown(press(made.canvas, pointAt(axis, 100.5)));
  made.owner.pointerMove(press(made.canvas, pointAt(axis, 900.5)));
  made.owner.pointerUp(press(made.canvas, pointAt(axis, 900.5)));
  assert.ok(made.owner.popup.value, 'the Region popup is open');
  return made;
};

test('while an action of the selection runs, every action of its popup is disabled and refused (OV-2561)', async () => {
  const add = deferred();
  const calls = [];
  const { owner } = openPopup({
    addRegionAnnotation: (region) => { calls.push(['add', region]); return add.promise; },
    readRegionFasta: async () => { calls.push(['copy']); return '>r\nACGT\n'; },
    newDrawingFromRegion: async () => { calls.push(['new']); return null; },
    selectFeatures: () => calls.push(['select'])
  });
  assert.equal(owner.popup.value.busy, false);
  const adding = owner.addToAnnotationSet();
  // The buttons read this one state; the second click is refused instead of being dropped behind an enabled button.
  assert.deepEqual([owner.popup.value.busy, owner.popup.value.add.busy, owner.popup.value.copy.busy], [true, true, false]);
  assert.equal(await owner.copyFasta(), false);
  assert.equal(await owner.newDrawing(), null);
  owner.selectFeaturesInside();
  assert.equal(await owner.addToAnnotationSet(), false);
  assert.deepEqual(calls.map(([name]) => name), ['add']);
  add.resolve(null);
  assert.equal(await adding, true);
  assert.equal(owner.popup.value.busy, false);
});

test('Add to annotation set that fails shows the normalized error under its button, tagged with its operation (OV-2562)', async () => {
  const { normalizeUserFacingError } = await import('../../gbdraw/web/js/utils/error-normalization.js');
  const { errorText } = await import('../../gbdraw/web/js/utils/error-text.js');
  const failure = normalizeUserFacingError({ code: 'INPUT_UNREADABLE' }, { operation: 'listSequenceRecords', stage: 'helper' });
  const { owner } = openPopup({ addRegionAnnotation: async () => { throw failure; } });
  assert.equal(await owner.addToAnnotationSet(), false);
  assert.deepEqual(owner.popup.value.add, { busy: false, status: errorText(failure.summaryParts), failed: true });
  assert.equal(owner.popup.value.busy, false);
  // The operation the action names is one the normalizer knows, not 'unknown'.
  assert.equal(normalizeUserFacingError(failure, { operation: 'addRegionAnnotation', stage: 'helper' }).operation, 'addRegionAnnotation');
});

// OV-2717 parity: Python's catalog of a record drawn in two windows, admitted by the Web, binds each drawn feature to
// its record (the key a Web feature edit names, which Python binds: tests/test_record_windows.py) and to its window.
test('a span of a window counts the features that window draws, by the record and the drawn unit', async () => {
  const { evaluatePython } = await import('./helpers/python-evaluator.mjs');
  const { admitFeatureCatalog } = await import('../../gbdraw/web/js/services/feature-catalog.js');
  const { readRegionAxisRecords } = await import('../../gbdraw/web/js/services/region-axis.js');
  const rendered = JSON.parse(await evaluatePython(`
import json, tempfile
from pathlib import Path
from gbdraw.api import GenBankInputSource, LinearDiagramRequest, RecordInput, resolve_request
from gbdraw.io.regions import RecordWindow
from gbdraw.session_request_codec import encode_canonical_request
from gbdraw.web_support.request_render import render_canonical_web_request
directory = Path(tempfile.mkdtemp())
request = LinearDiagramRequest(records=(RecordInput(GenBankInputSource(Path('tests/test_inputs/HmmtDNA.gbk')),
    record_key='mito', windows=(RecordWindow(1000, 6500), RecordWindow(7000, 9000))),))
encoded = encode_canonical_request(resolve_request(request))
paths = {}
for index, resource in enumerate(encoded.resources):
    path = directory / f'r{index}'
    path.write_bytes(resource.content if resource.content is not None else Path(str(resource.source_path)).read_bytes())
    paths[resource.resource_id] = str(path)
reply = render_canonical_web_request(encoded.payload, resource_paths=paths, output_directory=directory / 'out')
result = json.dumps({'catalog': reply['featureCatalog'], 'geometry': reply['metadata']['trackSlotGeometry'],
    'name': reply['results'][0]['name']})
`));
  const { extractedFeatures } = admitFeatureCatalog(rendered.catalog, [{ name: rendered.name }], { mode: 'linear' });
  const geometry = { ...rendered.geometry, records: rendered.geometry.records.map((record) => ({ ...record, resultIndex: 0 })) };
  const axisRecords = readRegionAxisRecords(geometry, 0);
  assert.deepEqual(axisRecords.map(({ recordIndex, recordKey, windowIndex }) => ({ recordIndex, recordKey, windowIndex })),
    [{ recordIndex: 0, recordKey: 'mito', windowIndex: 0 }, { recordIndex: 1, recordKey: 'mito', windowIndex: 1 }]);
  const rows = rendered.catalog.items[0].features;
  for (const axisRecord of axisRecords) {
    let open = null;
    const { owner } = setup({ ports: {
      readAxisRecords: () => axisRecords, openSelection: () => open, drawnFeatures: () => extractedFeatures
    } });
    // The whole window.
    const span = { start: 0, length: axisRecord.axis.source.drawnLength };
    open = { resultIndex: 0, recordIndex: axisRecord.recordIndex, recordKey: 'mito', windowIndex: axisRecord.windowIndex, span,
      source: { start: 0, end: 0, length: span.length, wrapsOrigin: false } };
    const drawnHere = rows.filter((row) => row.recordKey === 'mito' && row.recordIndex === axisRecord.recordIndex);
    // Every feature the window draws whole; one cut by the window's edge has a part outside it.
    const whole = drawnHere.filter((row) => !extractedFeatures.find((feature) => feature.svg_id === row.svgId)
      .location_parts?.some((part) => part.start < axisRecord.axis.source.anchorBase - 1));
    assert.ok(whole.length > 0, `window ${axisRecord.windowIndex + 1} draws features`);
    assert.equal(owner.popup.value.featureCount, whole.length, `the features of window ${axisRecord.windowIndex + 1}`);
  }
});
