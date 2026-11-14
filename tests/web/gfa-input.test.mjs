// A GFA 1 file as a drawing's input type (D-86, D-126): its card slot, its record source, and its readers. Only the
// Python loader (gbdraw/io/gfa.py through load_gfa) reads a GFA, so every browser reader declines it to the Worker.
import assert from 'node:assert/strict';
import test from 'node:test';

const {
  INPUT_TYPES, inputSequenceSource, inputSourceSlotFiles, inputTypeOfSourceKind, normalizeInputType, recordSourceOf,
  recordSourceSlotIds
} = await import('../../gbdraw/web/js/services/input-types.js');
const {
  discoverInputRecords, discoverWorkerInputRecords, readAllSequenceRecords, readParsedSequenceRecords
} = await import('../../gbdraw/web/js/services/record-discovery.js');
const { countSequenceRecords } = await import('../../gbdraw/web/js/services/genbank-header.js');
const { RECORD_SOURCE_FILES } = await import('../../gbdraw/web/js/contracts/schema.generated.js');

const GRAPH = 'H\tVN:Z:1.0\nS\tcircle\tACGTACGTAC\nS\tleft\tGGGG\nL\tcircle\t+\tcircle\t+\t2M\nP\tp1\tcircle+\t*\n';
/** @param {string} text @param {string} [name] */
const gfaFile = (text, name = 'graph.gfa') => ({
  name, size: text.length, lastModified: 0, arrayBuffer: async () => new TextEncoder().encode(text).buffer
});

test('a GFA input holds its file in the sequence slot and projects a gfa record source', () => {
  assert.equal(normalizeInputType('gfa'), 'gfa');
  assert.deepEqual(INPUT_TYPES.gfa, { primary: 'fasta', paired: null, sourceKind: 'gfa' });
  assert.equal(inputTypeOfSourceKind('gfa'), 'gfa');
  assert.deepEqual(inputSourceSlotFiles('gfa').map(({ resource, field, slot }) => ({ resource, field, slot })),
    [{ ...RECORD_SOURCE_FILES.gfa[0], slot: 'fasta' }]);
  assert.deepEqual(recordSourceOf('gfa', () => 'graph-id'), { kind: 'gfa', gfaResourceId: 'graph-id' });
  assert.deepEqual(recordSourceSlotIds({ kind: 'gfa', gfaResourceId: 'graph-id' }), { gb: '', gff: '', fasta: 'graph-id' });
  assert.equal(inputSequenceSource('gfa').format, 'gfa');
  assert.equal(inputSequenceSource('fasta').format, 'fasta');
  assert.equal(inputSequenceSource('gff').format, 'fasta');
});

test('the Worker lists a GFA input\'s records through load_gfa; the browser readers decline it', async () => {
  /** @type {{ operation: string, format: string, roles: string[] }[]} */
  const calls = [];
  const records = await discoverWorkerInputRecords({
    file: gfaFile(GRAPH), format: 'gfa',
    runHelperOperation: async (operation, payload) => {
      calls.push({ operation, format: payload.format, roles: payload.files.map(({ role }) => role) });
      return { result: { records: [
        { selector: '#1', record_id: 'circle', record_length: 8, topology: 'circular' },
        { selector: '#2', record_id: 'left', record_length: 4, topology: 'linear' }
      ] } };
    }
  });
  assert.deepEqual(records.map(({ recordId, detectedTopology }) => [recordId, detectedTopology]),
    [['circle', 'circular'], ['left', 'linear']]);
  assert.deepEqual(calls, [{ operation: 'listSequenceRecords', format: 'gfa', roles: ['source'] }]);
  assert.equal(typeof discoverInputRecords, 'function');

  const file = gfaFile(GRAPH, 'browser.gfa');
  assert.deepEqual(await readParsedSequenceRecords(file, 'gfa'), { records: null, parsed: false });
  /** @type {string[]} */
  const formats = [];
  const all = await readAllSequenceRecords(file, 'gfa', async (_operation, payload) => {
    formats.push(payload.format);
    return { result: { records: [{ id: 'circle', aliases: ['circle'], sequence: 'ACGTACGT' }] } };
  });
  assert.deepEqual(all.map(({ id, sequence }) => [id, sequence]), [['circle', 'ACGTACGT']]);
  assert.deepEqual(formats, ['gfa']);
});

test('a GFA file counts one record per segment (S line)', () => {
  assert.equal(countSequenceRecords(GRAPH, 'gfa'), 2);
  assert.equal(countSequenceRecords('>a\nAC\n>b\nGT\n', 'fasta'), 2);
  assert.equal(countSequenceRecords('LOCUS       a\n//\n', 'genbank'), 1);
});
