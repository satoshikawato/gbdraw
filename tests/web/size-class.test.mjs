import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import test from 'node:test';

import { createAutoValueDisplay } from '../../gbdraw/web/js/app/auto-value-display.js';
import { createCircularTrackSlotEditor } from '../../gbdraw/web/js/app/circular-track-slots.js';
import {
  circularDrawingRegionNote, circularEstimateLengths, circularFormRegion, circularRequestDiagrams, linearCardRegion,
  regionClampNote
} from '../../gbdraw/web/js/services/record-options.js';
import { circularSingleRecordKey } from '../../gbdraw/web/js/services/record-draw-selection.js';
import { createDefaultAdv, createDefaultForm } from '../../gbdraw/web/js/services/session-active-config-contract.js';
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
/** @typedef {Parameters<typeof circularRequestDiagrams>[0]} DiscoveredRecords */
/** @param {DiscoveredRecords} records @param {ReturnType<typeof circularDrawing>} drawing */
const circularDiagramRecordLengths = (records, drawing) => circularRequestDiagrams(records, drawing).map(({ lengths }) => lengths);
const circularNotes = (/** @type {ReturnType<typeof circularDrawing>} */ drawing, records = SHORT_AND_LONG) => createAutoValueDisplay(
  withDrawings({ mode: { value: 'circular' }, ...drawing }),
  { circularDiagramLengths: () => circularDiagramRecordLengths(records, drawing) }
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
  // A selector that matches no record, or every record OFF: the records the drawing offers, as before a choice.
  assert.deepEqual(lengths({ circular_record_selector: 'missing' }), [[30_000], [300_000]]);
  assert.deepEqual(lengths({ multi_record_canvas: true }, ['#1', '#2']), [[30_000, 300_000]]);
  // Not known: nothing read, or no length read.
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

// OV-518: Generate crops the one record of a single choice to its region
// before it sizes the diagram, and the notes and the track preview estimate
// the diagram at the lengths they know. One 16,569 bp record, Region 1..5000.
const MITOGENOME = Object.freeze([{ selector: '#1', record_id: 'mito', record_length: 16_569 }]);
/** @param {DiscoveredRecords} records @param {ReturnType<typeof circularDrawing>} drawing */
const featureWidthEstimate = (records, drawing, shownRequest = null, resultIndex = 0) => {
  const state = withDrawings({
    mode: { value: 'circular' }, ...drawing,
    form: { ...createDefaultForm(), ...drawing.form }, adv: { ...createDefaultAdv('circular'), ...drawing.adv }
  });
  const editor = createCircularTrackSlotEditor({
    state,
    estimateRecordLengths: () => circularEstimateLengths(circularRequestDiagrams(records, drawing), shownRequest, resultIndex)
  });
  editor.normalizeCircularTrackSlots();
  const features = state.adv.circular_track_slots.findIndex((/** @type {{ renderer: string }} */ slot) => slot.renderer === 'features');
  return editor.circularTrackSlotGeometryAutoText(state.adv.circular_track_slots[features], features, 'width');
};

test('OV-518: a record cropped to a region is drawn at the length of the region', () => {
  const cropped = circularDrawing({ circular_region_start: 1, circular_region_end: 5000 });
  assert.deepEqual(circularDiagramRecordLengths(MITOGENOME, cropped), [[5000]]);
  assert.equal(featureWidthEstimate(MITOGENOME, cropped), '≈ 37 px (estimate)');
  const notes = circularNotes(cropped, MITOGENOME);
  assert.equal(notes.autoValueText('labelFontSize'), '14 px (auto)');
  assert.equal(notes.autoValueText('scaleInterval'), '1 kbp (auto)');
  assert.equal(notes.autoValueText('dinucleotideWindow'), '1 kbp (auto)');

  // A long record cropped to a short region is drawn short.
  const longCropped = circularDrawing({ circular_record_selector: 'long', circular_region_start: 1, circular_region_end: 30_000 });
  assert.deepEqual(circularDiagramRecordLengths(SHORT_AND_LONG, longCropped), [[30_000]]);
  assert.equal(circularNotes(longCropped).autoValueText('labelFontSize'), '14 px (auto)');
  // A region Generate cannot draw yet (End not typed) leaves the record whole.
  assert.deepEqual(circularDiagramRecordLengths(MITOGENOME, circularDrawing({ circular_region_start: 1 })), [[16_569]]);
  // Only a record offered alone is cropped: separate diagrams draw whole records.
  assert.deepEqual(
    circularDiagramRecordLengths(SHORT_AND_LONG, circularDrawing({ circular_region_start: 1, circular_region_end: 5000 })),
    [[30_000], [300_000]]
  );
});

// OV-1097 (Owner decision 2026-10-18): Python draws a region whose End is past
// the record end to the record end; the Web predicts it once, and the lengths,
// the note beside the Region fields, and the request read that one prediction.
test('OV-1097: an End past the record end is drawn to the record end, and the note says so', () => {
  const past = circularDrawing({ circular_region_start: 10_000, circular_region_end: 20_000 });
  assert.deepEqual(circularDiagramRecordLengths(MITOGENOME, past), [[6_570]]);
  assert.equal(
    circularDrawingRegionNote(MITOGENOME, past),
    'Region reaches 20,000, past the record end (16,569 bp); Generate draws to 16,569.'
  );
  assert.deepEqual(circularFormRegion(past.form, MITOGENOME[0]), { start: 10_000, end: 20_000, drawnEnd: 16_569, recordLength: 16_569 });
  // A region within the record, or a record of unknown length, has no note.
  assert.equal(circularDrawingRegionNote(MITOGENOME, circularDrawing({ circular_region_start: 1, circular_region_end: 16_569 })), '');
  assert.equal(circularDrawingRegionNote([{ selector: '#1', record_id: 'unread' }], past), '');
  // Separate diagrams draw whole records, so no region and no note.
  assert.equal(circularDrawingRegionNote(SHORT_AND_LONG, circularDrawing({ circular_region_start: 1, circular_region_end: 400_000 })), '');
  // A Start past the record end stays an error, as Python refuses it.
  const late = circularFormRegion(circularDrawing({ circular_region_start: 16_570, circular_region_end: 20_000 }).form, MITOGENOME[0]);
  assert.ok(late && 'error' in late);
  assert.equal(/** @type {{ code?: string }} */ (late.error).code, 'REGION_INVALID');
  assert.deepEqual(/** @type {{ context?: unknown }} */ (late.error).context, { field: 'region', reason: 'RECORD_BOUNDS' });
  assert.deepEqual(circularDiagramRecordLengths(MITOGENOME, circularDrawing({ circular_region_start: 16_570, circular_region_end: 20_000 })), [[16_569]]);
});

test('a Region Generate cannot draw is a typed REGION_INVALID diagnostic, named by its reason', () => {
  const reasonOf = (start, end) => {
    const region = circularFormRegion(circularDrawing({ circular_region_start: start, circular_region_end: end }).form, MITOGENOME[0]);
    assert.ok(region && 'error' in region);
    const { code, context } = /** @type {{ code?: string, context?: unknown }} */ (region.error);
    return [code, context];
  };
  assert.deepEqual(reasonOf(5, null), ['REGION_INVALID', { field: 'region', reason: 'BOTH_ENDPOINTS' }]);
  assert.deepEqual(reasonOf(0, 9), ['REGION_INVALID', { field: 'region', reason: 'POSITIVE_INTEGER' }]);
  assert.deepEqual(reasonOf(9, 5), ['REGION_INVALID', { field: 'region', reason: 'ORDER' }]);
});

// The Linear Region fields have one reader: Generate's request, the sequences a
// match reads, and the note beside the fields all take its result.
test('OV-1097: a Linear card region is read once, with the span Generate draws', () => {
  assert.equal(linearCardRegion({ region_start: '', region_end: null }, 2, 48_502), null);
  const reasonOf = (start, end) => {
    const region = linearCardRegion({ region_start: start, region_end: end }, 2, 48_502);
    assert.ok(region && 'error' in region);
    return /** @type {{ context?: unknown }} */ (region.error).context;
  };
  assert.deepEqual(reasonOf(5, ''), { inputOrdinal: 2, reason: 'BOTH_ENDPOINTS' });
  assert.deepEqual(reasonOf('x', 9), { inputOrdinal: 2, reason: 'FINITE' });
  assert.deepEqual(reasonOf(1.5, 9), { inputOrdinal: 2, reason: 'INTEGER' });
  assert.deepEqual(reasonOf(0, 9), { inputOrdinal: 2, reason: 'POSITIVE_INTEGER' });
  const past = linearCardRegion({ region_start: 60_000, region_end: 40_000 }, 2, 48_502);
  assert.deepEqual(past, { start: 60_000, end: 40_000, drawn: { start: 40_000, end: 60_000, drawnEnd: 48_502, recordLength: 48_502 } });
  assert.equal(regionClampNote(past && 'drawn' in past ? past.drawn : null),
    'Region reaches 60,000, past the record end (48,502 bp); Generate draws to 48,502.');
  assert.deepEqual(linearCardRegion({ region_start: 50_000, region_end: 60_000 }, 2, 48_502), { start: 50_000, end: 60_000, drawn: null });
  // Without a known length the region is drawn as typed.
  assert.deepEqual(linearCardRegion({ region_start: '3', region_end: '9' }, 1), { start: 3, end: 9, drawn: { start: 3, end: 9, drawnEnd: 9 } });
});

test('OV-518: a record of unknown length leaves the lengths that are known', () => {
  const partlyRead = [{ selector: '#1', record_id: 'short', record_length: 30_000 }, { selector: '#2', record_id: 'unread' }];
  assert.deepEqual(circularDiagramRecordLengths(partlyRead, circularDrawing()), [[30_000]]);
  assert.deepEqual(circularDiagramRecordLengths(partlyRead, circularDrawing({ multi_record_canvas: true })), [[30_000]]);
  assert.equal(featureWidthEstimate(partlyRead, circularDrawing()), '≈ 37 px (estimate)');
  assert.equal(circularNotes(circularDrawing(), partlyRead).autoValueText('labelFontSize'), '14 px (auto)');
  // Long only while no length is known at all.
  assert.equal(featureWidthEstimate([{ selector: '#1', record_id: 'unread' }], circularDrawing()), '≈ 18.5 px (estimate)');
});

// OV-519: three separate records of 30 kb, 300 kb and 30 kb, drawn as a batch;
// Result 2 is shown and record 1 is turned off. The next Generate draws records
// 2 and 3, so the estimate is that of record 2 (the shown Result's record), not
// that of the second diagram of the next Generate (record 3).
test('OV-519: the track preview estimates the shown Result by its record, not its position', () => {
  const three = [
    { selector: '#1', record_id: 'a', record_length: 30_000 },
    { selector: '#2', record_id: 'b', record_length: 300_000 },
    { selector: '#3', record_id: 'c', record_length: 30_000 }
  ];
  const batch = { mode: 'circular', grouping: 'batch', records: [{ recordKey: 'record-1' }, { recordKey: 'record-2' }, { recordKey: 'record-3' }] };
  const firstOff = circularDrawing({}, ['#1']);
  assert.equal(featureWidthEstimate(three, firstOff, batch, 1), '≈ 18.5 px (estimate)');
  assert.equal(featureWidthEstimate(three, firstOff, batch, 2), '≈ 37 px (estimate)');
  // The shown Result's record is OFF now: the first diagram of the next Generate.
  assert.equal(featureWidthEstimate(three, circularDrawing({}, ['#2']), batch, 1), '≈ 37 px (estimate)');
  // A Result of one record chosen alone has the key of a record drawn alone;
  // the next Generate draws the same record among the others of a batch.
  const single = { mode: 'circular', grouping: 'single', records: [{ recordKey: circularSingleRecordKey({ recordId: 'b', selector: '#2' }) }] };
  assert.equal(featureWidthEstimate(three, circularDrawing(), single, 0), '≈ 18.5 px (estimate)');
  // A Linear Result names no Circular record.
  assert.equal(featureWidthEstimate(three, firstOff, { ...batch, mode: 'linear' }, 1), '≈ 18.5 px (estimate)');
});

/** @param {string} text "500 kbp (auto)" -> 500000 */
const shownBp = (text) => {
  const match = /^([\d.]+) (Mbp|kbp|bp)/.exec(text);
  assert.ok(match, text);
  return Number(match[1]) * { Mbp: 1_000_000, kbp: 1000, bp: 1 }[/** @type {'Mbp' | 'kbp' | 'bp'} */ (match[2])];
};

// Python's tiers reach the Web hints through mode-profiles.generated.js; every
// tier edge (the bound and its neighbours) is pinned to Python's value.
test('the scale interval and sliding window hints follow Python at every length-tier edge', () => {
  const vectors = JSON.parse(readFileSync(new URL('../fixtures/auto_size_vectors.json', import.meta.url), 'utf8'));
  assert.ok(vectors.autoHints.length > 40, 'edge vectors exist');
  for (const entry of vectors.autoHints) {
    const display = createAutoValueDisplay(withDrawings({
      mode: { value: entry.mode },
      form: {},
      adv: {}
    }), { circularDiagramLengths: () => [[entry.length]], linearRecordLengths: () => [entry.length] });
    const name = `${entry.mode} ${entry.length}`;
    assert.equal(shownBp(display.autoValueText('scaleInterval')), entry.scaleInterval, name);
    assert.equal(shownBp(display.autoValueText('dinucleotideWindow')), entry.window, name);
    assert.equal(shownBp(display.autoValueText('dinucleotideStep')), entry.step, name);
  }
});
