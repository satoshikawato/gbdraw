import assert from 'node:assert/strict';
import { cp, mkdtemp, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';

const repoRoot = process.cwd();
const sourceRoot = join(repoRoot, 'gbdraw', 'web', 'js', 'app');
const tempRoot = await mkdtemp(join(tmpdir(), 'gbdraw-record-selector-'));
await writeFile(join(tempRoot, 'package.json'), '{"type":"module"}', 'utf8');
await cp(join(sourceRoot, 'linear-record-selector.js'), join(tempRoot, 'linear-record-selector.js'));
await cp(join(sourceRoot, 'record-options.js'), join(tempRoot, 'record-options.js'));

const {
  AUTOMATIC_RECORD_OPTION_LABEL,
  buildRecordOptions,
  createLinearRecordSelector,
  formatRecordLength
} = await import(pathToFileURL(join(tempRoot, 'linear-record-selector.js')));
const {
  discoverGffFastaRecords,
  discoverSequenceRecords,
  normalizeSequenceRecords,
  parseSequenceRecordText
} = await import(pathToFileURL(join(sourceRoot, 'record-discovery.js')));

assert.equal(formatRecordLength(4641652), '4,641,652 bp');
assert.equal(formatRecordLength(null), 'length unavailable');
assert.equal(formatRecordLength(0), 'length unavailable');

const records = [
  { selector: '#1', recordId: 'RecA', recordLength: 1234 },
  { selector: '#2', recordId: 'RecB', recordLength: 567 }
];
assert.deepEqual(buildRecordOptions(records), [
  { value: '', label: AUTOMATIC_RECORD_OPTION_LABEL, synthetic: false },
  { value: 'RecA', label: 'RecA (1,234 bp)', synthetic: false },
  { value: 'RecB', label: 'RecB (567 bp)', synthetic: false }
]);
assert.deepEqual(
  buildRecordOptions([records[1], records[0]]).map((option) => option.value),
  ['', 'RecB', 'RecA']
);

const duplicateOptions = buildRecordOptions([
  { selector: '#1', recordId: 'RecA', recordLength: 1234 },
  { selector: '#2', recordId: 'RecA', recordLength: 1100 }
]);
assert.deepEqual(duplicateOptions.slice(1), [
  { value: '#1', label: 'RecA (1,234 bp) [#1]', synthetic: false },
  { value: '#2', label: 'RecA (1,100 bp) [#2]', synthetic: false }
]);
assert.deepEqual(
  buildRecordOptions([
    { selector: '#1', recordId: 'null', recordLength: 10 },
    { selector: '#2', recordId: '#named', recordLength: 20 }
  ]).slice(1).map((option) => option.value),
  ['#1', '#2']
);

const unmatchedOptions = buildRecordOptions(records, 'LegacyRec');
assert.deepEqual(unmatchedOptions[1], {
  value: 'LegacyRec',
  label: 'LegacyRec (not found in current file)',
  synthetic: true
});
assert.equal(buildRecordOptions(records, 'RecB').some((option) => option.synthetic), false);
assert.equal(buildRecordOptions([{ ...records[0], recordLength: undefined }])[1].label, 'RecA (length unavailable)');

assert.deepEqual(normalizeSequenceRecords({
  records: [
    { selector: '#1', record_id: 'RecA', record_length: 4 },
    { selector: '#2', record_id: '', record_length: 0 }
  ]
}), [
  { selector: '#1', recordId: 'RecA', recordLength: 4 },
  { selector: '#2', recordId: 'Record_2', recordLength: null }
]);
assert.throws(() => normalizeSequenceRecords({ records: [] }), /No records found/);
assert.throws(() => normalizeSequenceRecords({ error: 'Unsupported format: embl' }), /Unsupported format/);

assert.deepEqual(parseSequenceRecordText(`LOCUS       FIRST 12 bp DNA\nACCESSION   A1\nVERSION     A1.2\n//\nLOCUS       SECOND 7 bp DNA\nACCESSION   B1\n//\n`, 'genbank'), [
  { selector: '#1', recordId: 'A1.2', recordLength: 12 },
  { selector: '#2', recordId: 'B1', recordLength: 7 }
]);
assert.deepEqual(parseSequenceRecordText('>alpha note\nACGT\nAA\n>beta\nTTT\n', 'fasta'), [
  { selector: '#1', recordId: 'alpha', recordLength: 6 },
  { selector: '#2', recordId: 'beta', recordLength: 3 }
]);
let fastPathStagingCalls = 0;
const fastDiscovered = await discoverSequenceRecords({
  file: { text: async () => 'LOCUS       FAST 9 bp DNA\nVERSION     FAST.1\n//\n' },
  format: 'genbank',
  runHelperOperation: async () => { fastPathStagingCalls += 1; }
});
assert.equal(fastDiscovered[0].recordId, 'FAST.1');
assert.equal(fastPathStagingCalls, 0);

let sequenceHelperPayload = null;
const discovered = await discoverSequenceRecords({
  file: {
    name: 'records.fa',
    arrayBuffer: async () => new TextEncoder().encode('>FastaRec\nACGT\n').buffer
  },
  format: 'fasta',
  runHelperOperation: async (_operation, payload) => {
    sequenceHelperPayload = payload;
    return {
      result: {
        records: [{ selector: '#1', record_id: 'FastaRec', record_length: 10 }]
      }
    };
  }
});
assert.equal(discovered[0].recordId, 'FastaRec');
assert.equal(sequenceHelperPayload.files[0].role, 'source');
assert.equal(sequenceHelperPayload.files[0].bytes instanceof ArrayBuffer, true);

let pairedHelperPayload = null;
const paired = await discoverGffFastaRecords({
  gffFile: {
    name: 'records.gff',
    arrayBuffer: async () => new TextEncoder().encode('##gff-version 3\n').buffer
  },
  fastaFile: {
    name: 'records.fa',
    arrayBuffer: async () => new TextEncoder().encode('>GffOwned\nACGT\n').buffer
  },
  runHelperOperation: async (_operation, payload) => {
    pairedHelperPayload = payload;
    return {
      result: {
        records: [{ selector: '#1', record_id: 'GffOwned', record_length: 25 }]
      }
    };
  }
});
assert.equal(paired[0].recordId, 'GffOwned');
assert.deepEqual(pairedHelperPayload.files.map(({ role }) => role), ['gff', 'fasta']);

await assert.rejects(
  discoverSequenceRecords({
    file: {
      name: 'bad.gb',
      arrayBuffer: async () => new TextEncoder().encode('bad').buffer
    },
    format: 'genbank',
    runHelperOperation: async () => ({ result: { error: 'parse failed' } })
  }),
  /parse failed/
);

const ref = (value) => ({ value });
const fileA = { name: 'a.gb' };
const fileB = { name: 'b.gb' };
const state = {
  mode: ref('linear'),
  lInputType: ref('gb'),
  linearSeqs: [{ uid: 'row-a', gb: fileA, fasta: null, region_record_id: '' }]
};
const pending = new Map();
const recordReader = ({ primaryFile }) => new Promise((resolve) => pending.set(primaryFile, resolve));
const controller = createLinearRecordSelector({
  state,
  reactive: (value) => value,
  recordReader,
  logger: { warn: () => {} }
});

const loadingRefresh = controller.refresh();
assert.deepEqual(controller.optionsFor(state.linearSeqs[0]), [
  { value: '', label: 'Loading records...', synthetic: false }
]);
state.linearSeqs[0].region_record_id = 'LegacyRec';
assert.deepEqual(controller.optionsFor(state.linearSeqs[0]), [
  { value: 'LegacyRec', label: 'Loading records...', synthetic: true }
]);
state.linearSeqs[0].region_record_id = '';
pending.get(fileA)([{ selector: '#1', recordId: 'Initial', recordLength: 9 }]);
await loadingRefresh;

const firstRefresh = controller.refresh();
state.linearSeqs[0].gb = fileB;
const secondRefresh = controller.refresh();
pending.get(fileB)([{ selector: '#1', recordId: 'Current', recordLength: 20 }]);
await secondRefresh;
pending.get(fileA)([{ selector: '#1', recordId: 'Stale', recordLength: 10 }]);
await firstRefresh;
assert.equal(controller.optionsFor(state.linearSeqs[0])[1].value, 'Current');

state.linearSeqs.push({ uid: 'row-b', gb: fileA, fasta: null, region_record_id: 'Missing' });
const readyRefresh = controller.refresh();
assert.equal(controller.optionsFor(state.linearSeqs[1])[0].label, 'Loading records...');
pending.get(fileB)([{ selector: '#1', recordId: 'Current', recordLength: 20 }]);
await Promise.resolve();
pending.get(fileA)([{ selector: '#1', recordId: 'Other', recordLength: 30 }]);
await readyRefresh;
assert.match(controller.warningFor(state.linearSeqs[1]), /not found/);
assert.equal(controller.optionsFor(state.linearSeqs[1])[1].synthetic, true);

state.linearSeqs.reverse();
assert.equal(controller.optionsFor(state.linearSeqs[0]).at(-1).value, 'Other');
state.linearSeqs.splice(0, 1);
controller.purgeInactiveState();
assert.equal(Object.hasOwn(controller.selectorStateByUid, 'row-b'), false);
assert.equal(Object.hasOwn(controller.selectorStateByUid, 'row-a'), true);

const errorState = {
  mode: ref('linear'),
  lInputType: ref('gb'),
  linearSeqs: [{ uid: 'error-row', gb: null, fasta: null, region_record_id: 'LegacyRec' }]
};
const loggedErrors = [];
const errorController = createLinearRecordSelector({
  state: errorState,
  reactive: (value) => value,
  recordReader: async () => { throw new Error('parse details'); },
  logger: { warn: (...args) => loggedErrors.push(args) }
});
await errorController.refresh();
assert.deepEqual(errorController.optionsFor(errorState.linearSeqs[0]), [
  { value: 'LegacyRec', label: 'Upload a sequence file first', synthetic: true }
]);
errorState.linearSeqs[0].gb = { name: 'bad.gb' };
await errorController.refresh();
assert.equal(errorController.errorFor(errorState.linearSeqs[0]), 'Could not read records from this file.');
assert.deepEqual(errorController.optionsFor(errorState.linearSeqs[0]), [
  { value: 'LegacyRec', label: 'Records could not be loaded', synthetic: true }
]);
assert.equal(loggedErrors.length, 1);

const lazyState = {
  mode: ref('linear'),
  lInputType: ref('gb'),
  linearSeqs: [{ uid: 'lazy-row', gb: fileA, fasta: null, region_record_id: '' }]
};
let lazyReaderCalls = 0;
const lazyController = createLinearRecordSelector({
  state: lazyState,
  reactive: (value) => value,
  recordReader: async () => {
    lazyReaderCalls += 1;
    return [{ selector: '#1', recordId: 'Lazy', recordLength: 42 }];
  }
});
const lazyRefresh = lazyController.refresh();
const duplicateLazyRefresh = lazyController.refresh();
await Promise.all([lazyRefresh, duplicateLazyRefresh]);
assert.equal(lazyReaderCalls, 1);
assert.equal(lazyController.optionsFor(lazyState.linearSeqs[0])[1].value, 'Lazy');

const rollbackState = {
  mode: ref('linear'),
  lInputType: ref('gb'),
  semanticFileWatchersSuppressed: ref(false),
  sessionImportRollbackInProgress: ref(false),
  linearSeqs: [{
    uid: 'rollback-row',
    gb: fileA,
    fasta: null,
    region_record_id: 'Preserved'
  }]
};
let rejectStaleSelectorRead;
const rollbackController = createLinearRecordSelector({
  state: rollbackState,
  reactive: (value) => value,
  recordReader: () => new Promise((_resolve, reject) => {
    rejectStaleSelectorRead = reject;
  }),
  logger: { warn: () => {} }
});
const staleSelectorRefresh = rollbackController.refresh();
await Promise.resolve();
const selectorStateBeforeRollback = structuredClone(
  rollbackController.selectorStateByUid['rollback-row']
);
rollbackState.semanticFileWatchersSuppressed.value = true;
await rollbackController.refresh({ suppress: true });
rejectStaleSelectorRead(new Error('injected restored-file selector failure'));
await staleSelectorRefresh;
assert.deepEqual(
  rollbackController.selectorStateByUid['rollback-row'],
  selectorStateBeforeRollback
);
rollbackState.semanticFileWatchersSuppressed.value = false;

console.log('record selector tests passed');
