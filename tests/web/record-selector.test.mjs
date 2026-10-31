import assert from 'node:assert/strict';
import { cp, mkdtemp, readFile, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';

const repoRoot = process.cwd();
const sourceRoot = join(repoRoot, 'gbdraw', 'web', 'js', 'app');
const tempRoot = await mkdtemp(join(tmpdir(), 'gbdraw-record-selector-'));
await writeFile(join(tempRoot, 'package.json'), '{"type":"module"}', 'utf8');
await cp(join(sourceRoot, 'linear-record-selector.js'), join(tempRoot, 'linear-record-selector.js'));
await writeFile(join(tempRoot, 'linear-record-selector.js'), (await readFile(join(sourceRoot, 'linear-record-selector.js'), 'utf8'))
  .replace('../utils/i18n.js', pathToFileURL(join(sourceRoot, '../utils/i18n.js')).href)
  .replace('../utils/error-normalization.js', pathToFileURL(join(sourceRoot, '../utils/error-normalization.js')).href)
  .replace('../utils/error-text.js', pathToFileURL(join(sourceRoot, '../utils/error-text.js')).href)
  .replace('../services/record-options.js', pathToFileURL(join(sourceRoot, '../services/record-options.js')).href)
  .replace('../services/record-discovery.js', pathToFileURL(join(sourceRoot, '../services/record-discovery.js')).href), 'utf8');

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
  parseSequenceRecordText,
  readAllSequenceRecords
} = await import(pathToFileURL(join(sourceRoot, '..', 'services', 'record-discovery.js')));

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
  { selector: '#1', recordId: 'RecA', recordLength: 4, detectedTopology: 'unknown' },
  { selector: '#2', recordId: 'Record_2', recordLength: null, detectedTopology: 'unknown' }
]);
assert.throws(() => normalizeSequenceRecords({ records: [] }), /No records found/);
assert.throws(() => normalizeSequenceRecords({ error: 'Unsupported format: embl' }), /Unsupported format/);

assert.deepEqual(parseSequenceRecordText(`LOCUS       FIRST 12 bp DNA\nACCESSION   A1\nVERSION     A1.2\n//\nLOCUS       SECOND 7 bp DNA\nACCESSION   B1\n//\n`, 'genbank'), [
  { selector: '#1', recordId: 'A1.2', recordLength: 12, detectedTopology: 'unknown' },
  { selector: '#2', recordId: 'B1', recordLength: 7, detectedTopology: 'unknown' }
]);
assert.deepEqual(parseSequenceRecordText('>alpha note\nACGT\nAA\n>beta\nTTT\n', 'fasta'), [
  { selector: '#1', recordId: 'alpha', recordLength: 6, detectedTopology: 'unknown' },
  { selector: '#2', recordId: 'beta', recordLength: 3, detectedTopology: 'unknown' }
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

// A LOCUS line the browser reader cannot read exactly (a Geneious name with spaces) and an EMBL file
// decline to the Worker, whose reader (read_flat_file_records) reads both (R4).
for (const [name, text] of [
  ['geneious.gb', 'LOCUS       Contig 1 edited by me       4 bp    DNA     linear   UNK 01-JAN-2020\nORIGIN\n        1 acgt\n//\n'],
  ['records.embl', 'ID   EMBL1; SV 1; linear; DNA; STD; PRO; 4 BP.\nXX\nSQ   Sequence 4 BP;\n     acgt 4\n//\n']
]) {
  const helperCalls = [];
  const all = await readAllSequenceRecords(
    { name, arrayBuffer: async () => new TextEncoder().encode(text).buffer },
    'genbank',
    async (operation) => {
      helperCalls.push(operation);
      return { result: { records: [{ id: 'FROM_WORKER', aliases: ['FROM_WORKER'], sequence: 'ACGT' }] } };
    }
  );
  assert.deepEqual(all.map(({ id }) => id), ['FROM_WORKER'], name);
  assert.deepEqual(helperCalls, ['readAllSequenceRecords'], name);
}

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
  recordReader: async () => { throw {code:'INPUT_UNREADABLE',operation:'listSequenceRecords',stage:'helper',message:'PRIVATE_PARSE_SENTINEL'}; },
  logger: { warn: (...args) => loggedErrors.push(args) }
});
await errorController.refresh();
assert.deepEqual(errorController.optionsFor(errorState.linearSeqs[0]), [
  { value: 'LegacyRec', label: 'Upload a sequence file first', synthetic: true }
]);
errorState.linearSeqs[0].gb = { name: 'bad.gb' };
await errorController.refresh();
assert.match(errorController.errorFor(errorState.linearSeqs[0]), /input could not be read/i);
assert.equal(errorController.errorModelFor(errorState.linearSeqs[0]).code, 'INPUT_UNREADABLE');
assert.doesNotMatch(JSON.stringify(loggedErrors), /PRIVATE_PARSE_SENTINEL/);
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

// B4: a settled read of the same source files is reused by every later
// refresh (upload, Generate preparation, catalog check). A new File instance
// reads again, and so does a failure that is not a property of the bytes.
const reuseFiles = { gff: { name: 'm23.gff3' }, fasta: { name: 'mismatch.fasta' } };
const reuseState = {
  mode: ref('linear'),
  lInputType: ref('gff'),
  linearSeqs: [{ uid: 'reuse-a', gff: reuseFiles.gff, fasta: reuseFiles.fasta, region_record_id: '' }]
};
const reuseReads = [];
let reuseOutcome = () => {
  throw { code: 'FASTA_REQUIRED', operation: 'listGffFastaRecords', stage: 'helper', context: { reason: 'GFF_FASTA_MATCH' } };
};
const reuseController = createLinearRecordSelector({
  state: reuseState,
  reactive: (value) => value,
  recordReader: async ({ primaryFile, pairedFile }) => {
    reuseReads.push([primaryFile, pairedFile]);
    return reuseOutcome();
  },
  logger: { warn: () => {} }
});
for (let refresh = 0; refresh < 3; refresh += 1) await reuseController.refresh();
assert.equal(reuseReads.length, 1);
assert.deepEqual(
  (({ code, operation, stage, context }) => ({ code, operation, stage, context }))(
    reuseController.errorModelFor(reuseState.linearSeqs[0])
  ),
  { code: 'FASTA_REQUIRED', operation: 'listGffFastaRecords', stage: 'helper', context: { reason: 'GFF_FASTA_MATCH' } }
);
assert.match(reuseController.errorFor(reuseState.linearSeqs[0]), /matching FASTA/);

reuseState.linearSeqs[0].fasta = { name: 'mismatch.fasta' };
await reuseController.refresh();
assert.equal(reuseReads.length, 2);

reuseOutcome = () => { throw { code: 'WORKER_INIT', operation: 'listGffFastaRecords', stage: 'initialization' }; };
reuseState.linearSeqs[0].fasta = { name: 'retry.fasta' };
await reuseController.refresh();
await reuseController.refresh();
assert.equal(reuseReads.length, 4);
reuseOutcome = () => { throw new Error('Diagram generation worker failed.'); };
await reuseController.refresh();
await reuseController.refresh();
assert.equal(reuseReads.length, 6);

const discoveredRows = [];
const readyController = createLinearRecordSelector({
  state: reuseState,
  reactive: (value) => value,
  recordReader: async ({ primaryFile, pairedFile }) => {
    reuseReads.push([primaryFile, pairedFile]);
    return [{ selector: '#1', recordId: 'NC_001416.1', recordLength: 48502 }];
  },
  onRecordsDiscovered: ({ uid, records }) => { discoveredRows.push([uid, records[0].recordId]); }
});
reuseState.linearSeqs.push({ uid: 'reuse-b', gff: reuseState.linearSeqs[0].gff, fasta: reuseState.linearSeqs[0].fasta, region_record_id: '' });
await readyController.refresh();
await readyController.refresh();
assert.equal(reuseReads.length, 7);
assert.deepEqual(discoveredRows, [
  ['reuse-a', 'NC_001416.1'], ['reuse-b', 'NC_001416.1'],
  ['reuse-a', 'NC_001416.1'], ['reuse-b', 'NC_001416.1']
]);
assert.equal(readyController.optionsFor(reuseState.linearSeqs[1])[1].value, 'NC_001416.1');

// OV-1558: every reader of a source's records shares one discovery per File
// (or GFF3 and FASTA pair). A Linear -> Circular -> Linear trip and a draft
// Save read nothing again, also after a Generate transfer released the bytes.
{
  const { releaseFileContentCache } = await import(pathToFileURL(join(sourceRoot, '..', 'services', 'file-content-cache.js')));
  const { File } = await import('node:buffer');
  let fileReads = 0;
  const genbank = (id) => {
    const file = new File([`LOCUS       ${id}  10 bp    DNA     linear   BCT 01-JAN-2000\nORIGIN\n        1 acgtacgtac\n//\n`], `${id}.gb`);
    const read = file.arrayBuffer.bind(file);
    file.arrayBuffer = async () => { fileReads += 1; return read(); };
    return file;
  };
  const fileA = genbank('RecA');
  const fileB = genbank('RecB');
  const tripState = {
    lInputType: { value: 'gb' }, mode: { value: 'linear' },
    linearSeqs: [{ uid: 'trip-a', gb: fileA }, { uid: 'trip-b', gb: fileB }]
  };
  const trip = createLinearRecordSelector({
    state: tripState,
    reactive: (value) => value,
    recordReader: ({ primaryFile }) => discoverSequenceRecords({ file: primaryFile, format: 'genbank' })
  });
  const generateTransfer = () => [fileA, fileB].forEach((file) => releaseFileContentCache(file));
  await trip.refresh();
  assert.equal(fileReads, 2);
  generateTransfer();
  tripState.mode.value = 'circular';
  await trip.refresh();
  tripState.mode.value = 'linear';
  await trip.refresh();
  assert.equal(fileReads, 2, 'a mode round trip reads no card again');
  assert.equal(trip.optionsFor(tripState.linearSeqs[1])[1].value, 'RecB');
  // The draft Save's private discovery is a lookup.
  const draft = await discoverSequenceRecords({ file: fileA, format: 'genbank' });
  assert.equal(fileReads, 2, 'a draft Save reads no card again');
  // Each caller owns its copy of the records.
  draft[0].recordId = 'changed';
  assert.equal((await discoverSequenceRecords({ file: fileA, format: 'genbank' }))[0].recordId, 'RecA');

  // A failed read is not kept: the next caller reads again, and a success is kept.
  const bom = new File([new Uint8Array([0xef, 0xbb, 0xbf]), '>x\nACGT\n'], 'bom.fa');
  let helperCalls = 0;
  const runHelperOperation = async () => {
    helperCalls += 1;
    if (helperCalls === 1) throw Object.assign(new Error('Worker start failed.'), { code: 'UNKNOWN' });
    return { result: { records: [{ selector: '#1', record_id: 'x', record_length: 4 }] } };
  };
  await assert.rejects(discoverSequenceRecords({ file: bom, format: 'fasta', runHelperOperation }));
  assert.equal((await discoverSequenceRecords({ file: bom, format: 'fasta', runHelperOperation }))[0].recordId, 'x');
  await discoverSequenceRecords({ file: bom, format: 'fasta', runHelperOperation });
  assert.equal(helperCalls, 2);

  // A GFF3 and FASTA pair is one key.
  let gffReads = 0;
  const gff = new File(['##gff-version 3\nchr1\tsrc\tgene\t1\t4\t.\t+\t.\tID=g1\n'], 'a.gff');
  const fasta = new File(['>chr1\nACGT\n'], 'a.fa');
  [gff, fasta].forEach((file) => {
    const read = file.arrayBuffer.bind(file);
    file.arrayBuffer = async () => { gffReads += 1; return read(); };
  });
  await discoverGffFastaRecords({ gffFile: gff, fastaFile: fasta });
  [gff, fasta].forEach((file) => releaseFileContentCache(file));
  assert.equal((await discoverGffFastaRecords({ gffFile: gff, fastaFile: fasta }))[0].recordId, 'chr1');
  assert.equal(gffReads, 2);
}

// OV-1599: a Session resource view reads its bytes with the rule every browser reader applies. One that starts
// with a UTF-8 byte order mark goes to the Python reader, which reads the mark as content and skips the first
// record; the browser's text decoder would drop the mark and list that record.
{
  const { adoptCurrentSessionResources, createSessionResourceFileView } = await import(
    pathToFileURL(join(sourceRoot, '..', 'services', 'session-resource-backing.js'))
  );
  const bytes = new Uint8Array([0xef, 0xbb, 0xbf, ...new TextEncoder().encode('>F1\nAAAA\n>F2\nCCCC\n')]);
  const table = adoptCurrentSessionResources({
    'bom-fasta': {
      kind: 'fasta', name: 'bom.fa', type: 'text/plain', size: bytes.byteLength, lastModified: 7,
      encoding: 'base64', data: Buffer.from(bytes).toString('base64')
    }
  });
  const view = createSessionResourceFileView(table, 'bom-fasta');
  const formats = [];
  const records = await discoverSequenceRecords({
    file: view,
    format: 'fasta',
    runHelperOperation: async (_operation, payload) => {
      formats.push(payload.format);
      return { result: { records: [{ selector: '#1', record_id: 'F2', record_length: 4 }] } };
    }
  });
  assert.deepEqual(formats, ['fasta'], 'the Session resource goes to the Python reader');
  assert.deepEqual(records.map(({ recordId }) => recordId), ['F2']);
}

console.log('record selector tests passed');
