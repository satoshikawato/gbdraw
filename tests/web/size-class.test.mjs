import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import test from 'node:test';

import { createAutoValueDisplay } from '../../gbdraw/web/js/app/auto-value-display.js';
import { circularDiagramRecordLengths } from '../../gbdraw/web/js/services/record-options.js';
import { sizeClassForRecordLengths } from '../../gbdraw/web/js/services/track-slot-display.js';
import { withDrawings as withFixtureDrawings } from './helpers/drawing-state.mjs';

// Every drawing holds the settings its Session preserved (state.js); a case
// that keeps a boundary passes its own.
const withDrawings = (/** @type {Record<string, unknown>} */ fixture) => (
  withFixtureDrawings({ unmanagedConfigOverrides: {}, ...fixture })
);

test('the Web size class uses the packaged short/long boundary of each mode', () => {
  // The boundaries reach the Web through mode-profiles.generated.js (labels.length_threshold).
  for (const [mode, boundary] of /** @type {const} */ ([['circular', 75_000], ['linear', 100_000]])) {
    const drawing = { unmanagedConfigOverrides: {} };
    assert.equal(sizeClassForRecordLengths(mode, [3000, boundary - 1], drawing), 'short', mode);
    assert.equal(sizeClassForRecordLengths(mode, [3000, boundary], drawing), 'long', mode);
    assert.equal(sizeClassForRecordLengths(mode, [], drawing), null, mode);
  }
});

test('the Circular label font note follows the shared Auto size vectors, mixed size classes included', () => {
  const vectors = JSON.parse(readFileSync(new URL('../fixtures/auto_size_vectors.json', import.meta.url), 'utf8'));
  const cases = vectors.autoValues.filter((entry) => entry.mode === 'circular' && 'labels.font_size' in entry.expected);
  assert.ok(cases.some((entry) => entry.recordLengths.length > 1), 'a vector with records of both size classes');
  for (const entry of cases) {
    // The records of a vector share one canvas: one diagram.
    const display = createAutoValueDisplay(withDrawings({
      mode: { value: 'circular' },
      form: {},
      adv: {}
    }), { circularDiagramLengths: () => [entry.recordLengths] });
    assert.equal(display.autoValueText('labelFontSize'), `${entry.expected['labels.font_size']} px (auto)`, entry.name);
  }
});

// A Session that stores its whole configuration keeps the boundary it was saved
// with (labels.length_threshold.<mode> among the drawing's preserved settings),
// and Generate draws with it. The Auto notes and the track preview read the same
// boundary: 60,000 bp is short by both packaged boundaries and long by a kept 50 kb.
test('a drawing that keeps a stored boundary shows the Auto values Generate draws', () => {
  const lengths = [60_000];
  const drawing = (/** @type {Record<string, unknown>} */ unmanagedConfigOverrides) => ({
    form: {}, adv: {}, unmanagedConfigOverrides
  });
  const display = (/** @type {'circular' | 'linear'} */ mode, /** @type {Record<string, unknown>} */ overrides) => (
    createAutoValueDisplay(withDrawings({
      mode: { value: mode },
      ...drawing(overrides)
    }), { circularDiagramLengths: () => [lengths], linearRecordLengths: () => lengths })
  );

  for (const mode of /** @type {const} */ (['circular', 'linear'])) {
    const kept = { [`labels.length_threshold.${mode}`]: 50_000 };
    assert.equal(sizeClassForRecordLengths(mode, lengths, drawing({})), 'short', mode);
    assert.equal(sizeClassForRecordLengths(mode, lengths, drawing(kept)), 'long', mode);
    // The other mode's boundary is not this drawing's.
    const other = mode === 'circular' ? 'linear' : 'circular';
    assert.equal(
      sizeClassForRecordLengths(mode, lengths, drawing({ [`labels.length_threshold.${other}`]: 50_000 })),
      'short',
      mode
    );
  }

  const circular = { packaged: display('circular', {}), kept: display('circular', { 'labels.length_threshold.circular': 50_000 }) };
  assert.equal(circular.packaged.autoValueText('labelFontSize'), '14 px (auto)');
  assert.equal(circular.kept.autoValueText('labelFontSize'), '8 px (auto)');
  assert.equal(circular.kept.autoValueText('axisStrokeWidth'), '1 px (auto)');
  assert.equal(circular.kept.autoValueText('circularFeatureWidth'), '18.5 px (auto)');
  assert.equal(circular.packaged.autoValueText('circularFeatureWidth'), '37 px (auto)');

  const linear = { packaged: display('linear', {}), kept: display('linear', { 'labels.length_threshold.linear': 50_000 }) };
  assert.equal(linear.packaged.autoValueText('labelFontSize'), '24 px (auto)');
  assert.equal(linear.kept.autoValueText('labelFontSize'), '5 px (auto)');
  assert.equal(linear.kept.autoValueText('linearFeatureHeight'), '20 px (auto)');
  assert.equal(linear.kept.autoValueText('legendFontSize'), '16 px (auto)');
});

// OV-494: a file with a 30,000 bp and a 300,000 bp record. The request draws
// the ON records of a Multi-Record Canvas as one diagram, separate diagrams one
// record each, and a single choice its record; the notes read those records,
// not every record discovered in the file.
const SHORT_AND_LONG = Object.freeze([
  { selector: '#1', record_id: 'short', record_length: 30_000 },
  { selector: '#2', record_id: 'long', record_length: 300_000 }
]);
const circularDrawing = (/** @type {Record<string, unknown>} */ form = {}, /** @type {string[]} */ recordsOff = []) => ({
  form: {
    circular_record_selector: '', multi_record_canvas: false,
    circular_region_start: null, circular_region_end: null, ...form
  },
  adv: { circular_grouping_intent: 'auto' },
  recordsOff,
  unmanagedConfigOverrides: {}
});
const circularNotes = (/** @type {ReturnType<typeof circularDrawing>} */ drawing) => createAutoValueDisplay(
  withDrawings({ mode: { value: 'circular' }, ...drawing }),
  { circularDiagramLengths: () => circularDiagramRecordLengths(SHORT_AND_LONG, drawing) }
);

test('OV-494: the drawn Circular records are those of the request, one list per diagram', () => {
  const lengths = (/** @type {Record<string, unknown>} */ form, /** @type {string[]} */ recordsOff = []) => (
    circularDiagramRecordLengths(SHORT_AND_LONG, circularDrawing(form, recordsOff))
  );
  assert.deepEqual(lengths({ multi_record_canvas: true }), [[30_000, 300_000]]);
  assert.deepEqual(lengths({ multi_record_canvas: true }, ['#2']), [[30_000]]);
  assert.deepEqual(lengths({}), [[30_000], [300_000]], 'separate diagrams');
  assert.deepEqual(lengths({}, ['#2']), [[30_000]]);
  assert.deepEqual(lengths({ circular_record_selector: 'short' }), [[30_000]], 'a single choice');
  assert.deepEqual(lengths({ circular_record_selector: 'long' }, ['#2']), [[300_000]], 'a single choice draws its record ON or OFF');
  // Not known: a region crop of the one record, a selector that matches no record, every record OFF, nothing read.
  assert.deepEqual(lengths({ circular_record_selector: 'long', circular_region_start: 1, circular_region_end: 30_000 }), []);
  assert.deepEqual(lengths({ circular_record_selector: 'missing' }), []);
  assert.deepEqual(lengths({ multi_record_canvas: true }, ['#1', '#2']), []);
  assert.deepEqual(circularDiagramRecordLengths([], circularDrawing()), []);
  assert.deepEqual(circularDiagramRecordLengths([{ selector: '#1' }], circularDrawing()), []);
});

test('OV-494: the Circular Auto notes show the sizes of the records Generate draws', () => {
  const canvas = circularNotes(circularDrawing({ multi_record_canvas: true }));
  assert.equal(canvas.autoValueText('labelFontSize'), '8 px (auto)', 'a canvas of both size classes is long');
  assert.equal(canvas.autoValueText('circularFeatureWidth'), '18.5 px (auto)');
  // Each record of a canvas keeps the scale interval and the window of its own length.
  assert.equal(canvas.autoValueText('scaleInterval'), 'by record length (auto)');
  assert.equal(canvas.autoValueText('dinucleotideWindow'), '1 kbp (auto)');

  for (const longOff of [circularNotes(circularDrawing({ multi_record_canvas: true }, ['#2'])), circularNotes(circularDrawing({}, ['#2']))]) {
    assert.equal(longOff.autoValueText('labelFontSize'), '14 px (auto)');
    assert.equal(longOff.autoValueText('axisStrokeWidth'), '3 px (auto)');
    assert.equal(longOff.autoValueText('circularFeatureWidth'), '37 px (auto)');
    assert.equal(longOff.autoValueText('scaleInterval'), '1 kbp (auto)');
  }

  // Separate diagrams: the 30 kb record is drawn short and the 300 kb record long.
  const separate = circularNotes(circularDrawing());
  assert.equal(separate.autoValueText('labelFontSize'), 'Auto (by genome size)');
  assert.equal(separate.autoValueText('circularFeatureWidth'), 'Auto (by genome size)');
  assert.equal(separate.autoValueText('scaleInterval'), 'by record length (auto)');

  assert.equal(circularNotes(circularDrawing({ circular_record_selector: 'short' })).autoValueText('labelFontSize'), '14 px (auto)');
  assert.equal(circularNotes(circularDrawing({ circular_record_selector: 'long' })).autoValueText('labelFontSize'), '8 px (auto)');
});

// OV-493: Linear, File 1 = 40,579 bp, File 2 = 306,008 bp. Generate sizes the
// labels by the records in the label scope (Show Labels) and every other size
// by the longest record; which records a scope names comes from Python
// (mode-profiles.generated.js).
test('OV-493: the Linear Label Font Size note follows the records in the label scope', () => {
  const notes = (/** @type {string} */ show_labels_linear, /** @type {number[]} */ lengths = [40_579, 306_008]) => (
    createAutoValueDisplay(
      withDrawings({ mode: { value: 'linear' }, form: { show_labels_linear }, adv: {} }),
      { linearRecordLengths: () => lengths }
    )
  );
  assert.equal(notes('first').autoValueText('labelFontSize'), '24 px (auto)');
  assert.equal(notes('first').autoValueText('linearFeatureHeight'), '20 px (auto)', 'the other sizes stay long');
  assert.equal(notes('first', [306_008, 40_579]).autoValueText('labelFontSize'), '5 px (auto)');
  for (const scope of ['all', 'orthogroup_top', 'none']) {
    assert.equal(notes(scope).autoValueText('labelFontSize'), '5 px (auto)', scope);
  }
  assert.equal(notes('first', []).autoValueText('labelFontSize'), 'Auto (by genome size)');
});

