import assert from 'node:assert/strict';
import test from 'node:test';
import {
  buildRecordDisplayRows, parseRecordDisplayStart, reconcileRecordDisplayDrafts,
  recordDisplaySurface, selectedFeatureDisplayStart
} from '../../gbdraw/web/js/app/record-display-options.js';
import { parseSequenceRecordText } from '../../gbdraw/web/js/app/record-discovery.js';

const source = {};
const records = [
  { selector: '#1', recordId: 'same', recordLength: 100, detectedTopology: 'circular' },
  { selector: '#2', recordId: 'same', recordLength: 100, detectedTopology: 'linear' }
];
const rows = buildRecordDisplayRows({ scope: 'linear', sourceUid: 'card-1', source, records });

test('discovery observes topology without inferring it from the diagram or record name', () => {
  const parsed = parseSequenceRecordText('LOCUS       same 100 bp DNA circular\n//\nLOCUS       same 100 bp DNA linear\n//\nLOCUS       circular 100 bp DNA\n//\n', 'genbank');
  assert.deepEqual(parsed.map(r => r.detectedTopology), ['circular', 'linear', 'unknown']);
  assert.equal(parseSequenceRecordText('>circular\nACGT\n', 'fasta')[0].detectedTopology, 'unknown');
});

test('ALL rows and exact-one use instance selectors, with no duplicate-ID fallback', () => {
  assert.equal(rows.length, 2);
  assert.notEqual(rows[0].key, rows[1].key);
  assert.deepEqual(buildRecordDisplayRows({ scope: 'linear', sourceUid: 'card-1', source, records, selector: '#2' }), [rows[1]]);
  for (const selector of ['same', 'missing', '#3']) {
    assert.throws(() => buildRecordDisplayRows({ scope: 'linear', sourceUid: 'card-1', records, selector }), /selector/);
  }
});

test('source replacement purges drafts; inactive selectors survive rediscovery and reorder', () => {
  const drafts = rows.map(r => ({ scope: r.scope, sourceUid: r.sourceUid, selector: r.selector, topologyOverride: null, startCoordinate: 25 }));
  assert.deepEqual(reconcileRecordDisplayDrafts(drafts, [...rows].reverse()), drafts);
  assert.deepEqual(reconcileRecordDisplayDrafts(drafts, rows, ['card-1']), []);
  assert.deepEqual(reconcileRecordDisplayDrafts(drafts, rows.map(r => ({ ...r, sourceUid: 'replacement' }))), []);
});

test('blank remains null, explicit 1 survives, malformed input is never clamped', () => {
  for (const value of [null, '', '  ']) assert.equal(parseRecordDisplayStart(value), null);
  for (const value of [1, 100, '25']) assert.equal(parseRecordDisplayStart(value), Number(value));
  for (const value of [true, false, undefined, 0, -1, 1.5, 'oops', Infinity, [], {}]) {
    assert.throws(() => parseRecordDisplayStart(value), /integer/);
  }
});

test('detected topology, nullable reset, crop lock, and RC default stay separate from raw draft', () => {
  const draft = { topologyOverride: false, startCoordinate: 25 };
  assert.equal(recordDisplaySurface(rows[0], draft).startEnabled, false);
  assert.equal(draft.startCoordinate, 25);
  assert.equal(recordDisplaySurface(rows[0], { ...draft, topologyOverride: null }).startEnabled, true);
  assert.equal(recordDisplaySurface(rows[1]).startEnabled, false);
  assert.equal(recordDisplaySurface(rows[1], { topologyOverride: true }).startEnabled, true);
  assert.equal(recordDisplaySurface(rows[0], {}, { cropped: true }).startEnabled, false);
  assert.equal(recordDisplaySurface(rows[0], {}, { reverse: true }).currentStart, 100);
});

for (const strand of ['+', '-']) {
  for (const parts of [[[0, 1]], [[99, 100]], [[10, 14]], [[90, 100], [0, 5]], [[10, 13], [40, 43]]]) {
    test(`covered-base shortcut oracle ${strand} ${JSON.stringify(parts)}`, () => {
      const ordered = strand === '-' ? [...parts].reverse() : parts;
      const bases = ordered.flatMap(([start, end]) => {
        const sequence = Array.from({ length: end - start }, (_, i) => start + i + 1);
        return strand === '-' ? sequence.reverse() : sequence;
      });
      const feature = { record_key: 'instance', strand, location_parts: ordered.map(([start, end]) => ({ start, end, strand })) };
      const args = { row: rows[0], committedRow: { ...rows[0], recordKey: 'instance' }, selectedFeatures: [feature] };
      assert.equal(selectedFeatureDisplayStart({ ...args, shortcut: 'five-prime' }), bases[0]);
      assert.equal(selectedFeatureDisplayStart({ ...args, shortcut: 'midpoint' }), bases[Math.floor((bases.length - 1) / 2)]);
    });
  }
}

test('shortcuts reject stale, unbound, cross-record, empty, and unknown/mixed-strand selections', () => {
  const feature = { record_key: 'instance', strand: '+', location_parts: [{ start: 5, end: 10, strand: '+' }] };
  const args = { row: rows[0], committedRow: { ...rows[0], recordKey: 'instance' }, selectedFeatures: [feature], shortcut: 'midpoint' };
  for (const changed of [
    { committedRow: null }, { committedRow: { ...args.committedRow, source: {} } },
    { committedRow: { ...args.committedRow, selector: '#2' } },
    { selectedFeatures: [] }, { selectedFeatures: [feature, feature] },
    { selectedFeatures: [{ ...feature, record_key: 'other' }] },
    { selectedFeatures: [{ ...feature, strand: 'unknown' }] },
    { selectedFeatures: [{ ...feature, location_parts: [] }] },
    { selectedFeatures: [{ ...feature, location_parts: [{ start: 5, end: 10, strand: '-' }] }] }
  ]) assert.throws(() => selectedFeatureDisplayStart({ ...args, ...changed }));
});
