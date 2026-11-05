// Select region's converter (gbdraw/web/js/services/region-axis.js) reads the
// axis Python writes and the same vectors as tests/test_region_axis_geometry.py:
// each probe is a point in a record's frame with the display index and source
// base Python says it shows.
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import test from 'node:test';

import {
  axisPosition, circularDragSpan, displayIndexOfRecordOffset, linearDragSpan, linearRowPosition, partsInsideSpan,
  readRegionAxisRecords, regionBandPath, sourceBaseAt, sourceSpanOf
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
  // A sweep of a whole turn or more is the whole record.
  assert.equal(circularDragSpan(axis, 10, -20_000).length, axis.source.drawnLength);
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
