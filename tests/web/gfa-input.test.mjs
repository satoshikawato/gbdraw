// A GFA 1 file as a drawing's input type (D-86, D-126): its card slot, its record source, and its readers. Only the
// Python loader (gbdraw/io/gfa.py through load_gfa) reads a GFA, so every browser reader declines it to the Worker.
import assert from 'node:assert/strict';
import test from 'node:test';

const {
  INPUT_TYPES, inputSequenceSource, inputSourceSlotFiles, inputTypeOfSourceKind, normalizeInputType, recordSourceOf,
  recordSourceSlotIds
} = await import('../../gbdraw/web/js/services/input-types.js');
const {
  countWorkerInputRecords, discoverInputRecords, discoverWorkerInputRecords, readAllSequenceRecords, readParsedSequenceRecords
} = await import('../../gbdraw/web/js/services/record-discovery.js');
const { readCanonicalResourceRecordCount } = await import('../../gbdraw/web/js/services/session-request.js');
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

// A Source recipe names a record by #index only when its file holds more than one. The browser counts GenBank and
// FASTA records; a GFA's count is the Worker's listing (load_gfa), so no JavaScript copy of the segment rule exists.
test('a GFA resource counts the records the Worker lists; GenBank and FASTA count in the browser', async () => {
  /** @param {string} name @param {string} text */
  const descriptor = (name, text) => ({
    kind: 'source', name, type: 'text/plain', size: text.length, lastModified: 0, encoding: 'base64',
    data: Buffer.from(text).toString('base64')
  });
  const resources = { graph: descriptor('graph.gfa', GRAPH), reads: descriptor('reads.fasta', '>a\nAC\n>b\nGT\n') };
  /** @type {{ operation: string, format: string, name: string, text: string }[]} */
  const calls = [];
  /** @param {string} operation @param {{ format: string, files: { name: string, blob: Blob }[] }} payload */
  const runHelperOperation = async (operation, payload) => {
    calls.push({ operation, format: payload.format, name: payload.files[0].name, text: await payload.files[0].blob.text() });
    return { result: { records: [
      { selector: '#1', record_id: 'circle', record_length: 8, topology: 'circular' },
      { selector: '#2', record_id: 'left', record_length: 4, topology: 'linear' },
      { selector: '#3', record_id: 'right', record_length: 4, topology: 'linear' }
    ] } };
  };
  /** @param {object} file @param {'gfa'} format */
  const countWorkerRecords = (file, format) => countWorkerInputRecords(file, format, runHelperOperation);
  assert.equal(await readCanonicalResourceRecordCount(resources, 'graph', 'gfa', countWorkerRecords), 3);
  assert.deepEqual(calls, [{ operation: 'listSequenceRecords', format: 'gfa', name: 'graph.gfa', text: GRAPH }]);
  assert.equal(await readCanonicalResourceRecordCount(resources, 'reads', 'fasta', countWorkerRecords), 2);
  assert.equal(calls.length, 1);
  await assert.rejects(readCanonicalResourceRecordCount(resources, 'graph', 'gfa'), /Only the Worker counts/);
  assert.equal(countSequenceRecords('LOCUS       a\n//\n', 'genbank'), 1);
});
