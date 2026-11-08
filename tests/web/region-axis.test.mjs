// Select region's converter (gbdraw/web/js/services/region-axis.js) reads the
// axis Python writes and the same vectors as tests/test_region_axis_geometry.py:
// each probe is a point in a record's frame with the display index and source
// base Python says it shows.
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import test from 'node:test';

import {
  axisPosition, circularDragSpan, displayIndexOfRecordOffset, featureParts, linearDragSpan, linearRowPosition, onTheWay,
  partsInsideSpan, partsMidpoint, readRegionAxisRecords, regionBandPath, sourceBaseAt, sourceSpanOf
} from '../../gbdraw/web/js/services/region-axis.js';

const VECTORS = JSON.parse(readFileSync(new URL('../fixtures/region_axis_geometry_vectors.json', import.meta.url), 'utf8')).cases;
const caseRecords = (name) => VECTORS.find((entry) => entry.name === name).records;
const geometryOf = (name, resultIndex = 0) => ({
  schema: 1, mode: caseRecords(name)[0].axis.kind, source: 'resolved',
  records: caseRecords(name).map((record) => ({ ...record, resultIndex, recordId: `id-${record.recordKey}`, slots: [] }))
});

test('a point Python names shows the display index and source base Python names', () => {
  for (const { name, records } of VECTORS) {
    for (const record of records) {
      for (const probe of record.probes) {
        const [x, y] = probe.point;
        const { position, inBand } = axisPosition(record.axis, { x, y });
        assert.equal(Math.floor(position), probe.displayIndex, `${name} ${record.recordKey} display ${probe.displayIndex}`);
        assert.equal(sourceBaseAt(record.axis.source, probe.displayIndex), probe.sourceBase, `${name} source base`);
        assert.equal(displayIndexOfRecordOffset(record.axis.source, probe.recordOffset), probe.displayIndex, `${name} record offset`);
        if (record.axis.kind === 'circular') assert.ok(inBand, `${name}: a probe lies on the feature tracks`);
      }
    }
  }
});

test('the reader takes the records of one Result and drops a record without an axis', () => {
  const geometry = geometryOf('circular-multi-record-canvas');
  assert.deepEqual(readRegionAxisRecords(geometry, 0).map((record) => record.recordKey), ['mito', 'lambda']);
  assert.deepEqual(readRegionAxisRecords(geometry, 1), []);
  // A Result generated before Python wrote the axis has none: Select region is unavailable.
  const old = { ...geometry, records: geometry.records.map(({ axis: _axis, ...rest }) => rest) };
  assert.deepEqual(readRegionAxisRecords(old, 0), []);
  assert.deepEqual(readRegionAxisRecords(null, 0), []);
  const broken = { ...geometry, records: [{ ...geometry.records[0], axis: { ...geometry.records[0].axis, source: { step: 2 } } }] };
  assert.deepEqual(readRegionAxisRecords(broken, 0), []);
  // The record's circularity is required: without it a drag could not tell whether it may cross the origin.
  const { isCircular: _circular, ...withoutCircularity } = geometry.records[0].axis.source;
  const unknownTopology = { ...geometry, records: [{ ...geometry.records[0], axis: { ...geometry.records[0].axis, source: withoutCircularity } }] };
  assert.deepEqual(readRegionAxisRecords(unknownTopology, 0), []);
});

test('a Circular drag across the display origin gives a span that wraps the source origin', () => {
  const [mito] = caseRecords('circular-display-start');
  const { axis } = mito;
  // Display start 16,000: display index 569 shows base 16,569 and 570 base 1. A drag from 560 to 580.
  const span = circularDragSpan(axis, 560.5, 20);
  assert.deepEqual(span, { start: 560, length: 21 });
  assert.deepEqual(sourceSpanOf(axis.source, span), { start: 16_560, end: 11, length: 21, wrapsOrigin: true });
  // The same arc dragged the other way gives the same span.
  assert.deepEqual(circularDragSpan(axis, 580.5, -20), span);
  // Across display index 0 (12 o'clock) on a whole record.
  assert.deepEqual(circularDragSpan(axis, 16_560.5, 20), { start: 16_560, length: 21 });
  assert.deepEqual(sourceSpanOf(axis.source, { start: 16_560, length: 21 }), { start: 15_991, end: 16_011, length: 21, wrapsOrigin: false });
});

// The origin-crossing twin (tools/twins.json): the Web says a span crosses the origin from wrapsOrigin; Python draws
// the same source span as one window across the origin (tests/test_region_drawing.py reads the same vectors).
const ORIGIN_CROSSING = JSON.parse(readFileSync(new URL('../fixtures/region_origin_crossing_vectors.json', import.meta.url), 'utf8')).cases;

test('a span crosses the origin on the Web exactly when Python draws it as one window across the origin', () => {
  assert.ok(ORIGIN_CROSSING.some(({ pythonOutcome }) => pythonOutcome === 'acrossOrigin'));
  assert.ok(ORIGIN_CROSSING.some(({ pythonOutcome }) => pythonOutcome === 'withinRecord'));
  for (const { name, axisSource, displaySpan, sourceSpan, pythonOutcome } of ORIGIN_CROSSING) {
    const source = sourceSpanOf(axisSource, displaySpan);
    assert.deepEqual(source, sourceSpan, name);
    assert.equal(source.wrapsOrigin, pythonOutcome === 'acrossOrigin', name);
  }
});

test('a whole turn selects the whole record, 1 to its length, not a span across the origin', () => {
  const [mito] = caseRecords('circular-display-start');
  const { axis } = mito;
  const whole = { start: 0, length: 16_569 };
  for (const [start, sweep] of [[10, -20_000], [4_142.5, 16_569], [100.5, -16_570]]) {
    assert.deepEqual(circularDragSpan(axis, start, sweep), whole, `from ${start}, sweep ${sweep}`);
  }
  // Display start 16,000 draws base 16,000 first; the whole record is still 1-16,569.
  assert.deepEqual(sourceSpanOf(axis.source, whole), { start: 1, end: 16_569, length: 16_569, wrapsOrigin: false });
  // Every feature lies inside it, also one that crosses the display origin or the source origin.
  assert.ok(partsInsideSpan(axis.source, whole, [{ start: 15_990, end: 16_010 }]));
  assert.ok(partsInsideSpan(axis.source, whole, [{ start: 16_560, end: 16_569 }, { start: 0, end: 9 }]));
});

test('a Circular drag on a record that is not circular, or on a crop of every base of one, stops at its ends', () => {
  for (const name of ['circular-linear-topology', 'circular-fasta-only', 'circular-whole-crop']) {
    const [record] = caseRecords(name);
    const { axis } = record;
    const length = axis.source.drawnLength;
    assert.equal(axis.source.isCircular, false, name);
    assert.equal(length, axis.source.sourceLength, `${name} draws the whole record`);
    // Clockwise through 12 o'clock: the span ends at the last base.
    const late = circularDragSpan(axis, length - 999.5, 2_000);
    assert.deepEqual(late, { start: length - 1_000, length: 1_000 }, name);
    assert.deepEqual(sourceSpanOf(axis.source, late), { start: length - 999, end: length, length: 1_000, wrapsOrigin: false }, name);
    // Counter-clockwise through 12 o'clock: the span starts at the first base.
    assert.deepEqual(circularDragSpan(axis, 500.5, -1_000), { start: 0, length: 501 }, name);
    // A whole turn is the whole record.
    assert.deepEqual(sourceSpanOf(axis.source, circularDragSpan(axis, 0.5, length + 10)), { start: 1, end: length, length, wrapsOrigin: false }, name);
  }
});

test('a Linear drag stays on its row, and a reversed crop gives its span on the source strand', () => {
  const [wrapped] = caseRecords('linear-wrap-reverse');
  // The second row shows display 20,000..39,999; a drag past its end stops at its last base.
  const start = linearRowPosition(wrapped.axis, 1, 100).position;
  const end = linearRowPosition(wrapped.axis, 1, 99_999).position;
  const span = linearDragSpan(start, end);
  assert.equal(span.start + span.length, 40_000);
  const source = sourceSpanOf(wrapped.axis.source, span);
  assert.equal(source.end, 48_502 - Math.floor(start));
  assert.equal(source.start, 48_502 - 39_999);
  assert.equal(source.wrapsOrigin, false);

  const [crop] = caseRecords('linear-two-records-crop-reverse');
  // The crop 1001-9000:rc draws 9,000 first: its first 100 drawn bases are 8,901..9,000.
  assert.deepEqual(sourceSpanOf(crop.axis.source, { start: 0, length: 100 }), { start: 8_901, end: 9_000, length: 100, wrapsOrigin: false });
});

test('a feature counts only when every part lies inside the span', () => {
  const [mito] = caseRecords('circular-display-start');
  const { source } = mito.axis;
  const span = { start: 560, length: 21 }; // source 16,560..11
  // Record bases are 0-based, end excluded: base 16,561..16,569 is [16_560, 16_569).
  assert.ok(partsInsideSpan(source, span, [{ start: 16_560, end: 16_569 }]));
  assert.ok(partsInsideSpan(source, span, [{ start: 16_565, end: 16_569 }, { start: 0, end: 5 }]));
  assert.ok(!partsInsideSpan(source, span, [{ start: 16_565, end: 16_569 }, { start: 0, end: 20 }]));
  assert.ok(!partsInsideSpan(source, span, []));
});

// OV-2718: the catalog gives a feature's parts in the input file's coordinates; on a cropped or reverse-complemented
// record they are not the drawn record's own bases.
test('a feature of a cropped, reverse-complemented record counts where its source bases are drawn', () => {
  const [crop] = caseRecords('linear-two-records-crop-reverse');
  const { source } = crop.axis;
  // The crop 1001-9000:rc draws 9,000 first: bases 8,901..9,000, [8_900, 9_000), are display 0..99.
  assert.ok(partsInsideSpan(source, { start: 0, length: 100 }, [{ start: 8_900, end: 9_000 }]));
  assert.ok(!partsInsideSpan(source, { start: 0, length: 100 }, [{ start: 8_899, end: 9_000 }]));
  assert.ok(partsInsideSpan(source, { start: 7_900, length: 100 }, [{ start: 1_000, end: 1_100 }]));
  assert.equal(partsMidpoint(source, [{ start: 8_900, end: 9_000 }]), 50);
});

test('the band covers the span on the feature tracks', () => {
  const [lambda] = caseRecords('linear-single');
  const row = lambda.axis.rows[0];
  const path = regionBandPath(lambda.axis, { start: 0, length: 48_502 });
  assert.equal(path, `M ${row.xStartPx} ${row.bandTopPx} h ${row.xEndPx - row.xStartPx} V ${row.bandBottomPx} H ${row.xStartPx} Z`);
  const [wrapped] = caseRecords('linear-wrap-reverse');
  // A span on the second row draws only there.
  assert.equal(regionBandPath(wrapped.axis, { start: 20_100, length: 50 }).split('M').length - 1, 1);
  const [mito] = caseRecords('circular-display-start');
  assert.match(regionBandPath(mito.axis, { start: 0, length: 100 }), /^M .* A .* L .* A .* Z$/);
});

test('a Shift+click range runs clockwise on a Circular record and between the two features on a Linear one (OV-2481)', () => {
  const [mito] = caseRecords('circular-display-start');
  const { axis } = mito;
  const middle = (feature) => partsMidpoint(axis.source, featureParts(feature));
  // Four features on one arc, the second on the other strand: the range from the first to the third holds the second.
  const [f1, f2, f3, f4] = [[100, 300], [400, 600], [700, 900], [1_000, 1_200]].map(([start, end]) => middle({ start, end }));
  assert.deepEqual([f1, f2, f3, f4].map((position) => onTheWay(axis, f1, f3, position)), [true, true, true, false]);
  // Clockwise past the display origin, and the way back the other way round.
  const late = middle({ start: 16_400, end: 16_500 });
  assert.equal(onTheWay(axis, late, f2, f1), true);
  assert.equal(onTheWay(axis, late, f2, f4), false);
  assert.equal(onTheWay(axis, f2, late, f4), true);
  assert.equal(onTheWay(axis, f2, late, f1), false);
  // A joined feature's middle lies between its first and its last part along the axis.
  assert.equal(middle({ location_parts: [{ start: 16_560, end: 16_569 }, { start: 0, end: 9 }] }), mod(displayIndexOfRecordOffset(axis.source, 16_560) + 9, 16_569));
  const [lambda] = caseRecords('linear-single');
  assert.deepEqual([10, 50, 90].map((position) => onTheWay(lambda.axis, 80, 20, position)), [false, true, false]);
});

/** @param {number} value @param {number} modulus */
const mod = (value, modulus) => ((value % modulus) + modulus) % modulus;
