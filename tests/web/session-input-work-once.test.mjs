// Save and Load do each piece of input work once (OV-1559, OV-1567): a File
// Generate described is bound by its descriptor with no hash or decode and
// encoded once per Save (OV-2612),
// and the Worker check of a slice's unmanaged overrides runs once per distinct
// request across Load and every Save.
import { savedDrawingView } from './helpers/saved-drawing-view.mjs';
import assert from 'node:assert/strict';
import { test } from 'node:test';
import { webcrypto } from 'node:crypto';
import { File } from 'node:buffer';
import { installSessionImportWorker } from './helpers/session-import-node.mjs';
import { installFakeSvgDom } from './fake-svg-dom.mjs';
import { bindingDrawings, boundState } from './helpers/input-files.mjs';
import { modeDrawing } from './helpers/drawing-state.mjs';
import { callerSession } from './helpers/session-caller.mjs';
import { sessionRollbackPort } from './helpers/session-rollback-port.mjs';

if (!globalThis.crypto) globalThis.crypto = webcrypto;
globalThis.window = {
  Vue: {
    ref: (value) => ({ value }), reactive: (value) => value,
    computed: (getter) => ({ get value() { return getter(); } }),
    nextTick: async () => {}
  },
  DOMPurify: { sanitize: (value) => value }
};
globalThis.document = { createElement: () => ({ click: () => {} }) };
installFakeSvgDom();
globalThis.File = File;
globalThis.alert = () => {};
installSessionImportWorker();

// One Circular File card holding a GenBank file.
const card = (gb) => ({ uid: 'record-1', gb, gff: null, fasta: null });
const metrics = [];
globalThis.__GBDRAW_TEST_HOOKS__ = { onStructuralMetric: (metric) => metrics.push(metric) };
const count = (name) => metrics.filter((metric) => metric.name === name).length;

const sessionConfig = await import('../../gbdraw/web/js/services/config.js');
const { buildSessionResources } = await import('../../gbdraw/web/js/services/session-resources.js');
const { gallerySessionFile, writtenResourcesWith } = await import('./helpers/written-resources.mjs');
const { readSessionContainer } = await import('../../gbdraw/web/js/services/session-container.js');
// What Save writes for each resource, read as base64 descriptors.
const encodeSessionResources = writtenResourcesWith(await import('../../gbdraw/web/js/services/session-resource-backing.js'));
/** Save's resource table: allocated, then encoded by the Session file writer. */
const savedSession = async (...args) => {
  const built = await buildSessionResources(...args);
  return { ...built, resources: await encodeSessionResources(built.resources) };
};
const { adoptRuntimeCanonicalSession } = await import('../../gbdraw/web/js/services/session-authority.js');
const { contentTokenOfBytes, holdReadBytes, inputContentBlob, readFileText, releaseFileContentCache, residentBytesOf, useInputCopyStore } = await import('../../gbdraw/web/js/services/file-content-cache.js');
const { createDiagramResourceTransport } = await import('../../gbdraw/web/js/services/diagram-resource-staging.js');
const { discoverSequenceRecords, readOncePerFile, readParsedSequenceRecords } = await import('../../gbdraw/web/js/services/record-discovery.js');
const { CANONICAL_REQUEST_SCHEMA } = await import('../../gbdraw/web/js/contracts/schema.generated.js');
const { transformLegacyResultSvg } = await import('../../gbdraw/web/js/app/app-setup.js');
const { state } = await import('../../gbdraw/web/js/state.js');
const { exportSession, importSession, serializeActiveRenderFiles, setUnmanagedConfigOverrideValidator } = callerSession(sessionConfig, state, await sessionRollbackPort(state));

test('Save binds a File Generate described with no hash or decode and reads it once (OV-1559, OV-2612)', async () => {
  const text = `>x\n${'ACGT'.repeat(250_000)}\n`;
  const file = new File([text], 'g.fa', { type: '' });
  let reads = 0;
  let streams = 0;
  const read = file.arrayBuffer.bind(file);
  file.arrayBuffer = async () => { reads += 1; return read(); };
  const stream = file.stream.bind(file);
  file.stream = () => { streams += 1; return stream(); };
  // Generate: the request describes the File and stages its bytes (the transfer releases the byte cache).
  const sourceState = boundState({ files: {}, circularSources: [card(file)], linearSeqs: [] });
  const releaseReads = holdReadBytes();
  const described = (await serializeActiveRenderFiles(sourceState, { ...modeDrawing(state, 'circular'), files: sourceState.files, circularSources: sourceState.circularSources })).circularSources[0].gb;
  await inputContentBlob(file);
  releaseFileContentCache(file);
  releaseReads();
  const committed = adoptRuntimeCanonicalSession({
    renderRequest: {
      schema: CANONICAL_REQUEST_SCHEMA, mode: 'circular',
      records: [{ recordKey: 'r1', source: { kind: 'genbank', resourceId: 'src' }, presentation: { label: 'x' } }],
      diagramOptions: { trackSlots: [], annotations: { sets: [] } }, layout: {}, comparisons: [],
      output: { prefix: 'x', formats: ['svg'], overwrite: false }
    },
    resources: { src: described },
    webFiles: {}
  });
  assert.equal(reads, 1);
  assert.equal(count('base64EncodeCount'), 0);
  // Each Save streams the File once and never reads it whole: the page keeps no input bytes between Saves.
  for (const save of [1, 2]) {
    metrics.length = 0;
    const saved = await savedSession(sourceState, committed, bindingDrawings(sourceState));
    assert.deepEqual({
      reads,
      streams,
      hashes: count('resourceIdentityHashCount'),
      decodes: count('base64DecodeCount'),
      byteReads: count('resourceByteReadCount'),
      encodes: count('base64EncodeCount')
    }, { reads: 1, streams: save, hashes: 0, decodes: 0, byteReads: 0, encodes: 0 }, `Save ${save}`);
    assert.deepEqual(Object.keys(saved.resources), ['src']);
    assert.equal(saved.resources.src.data, Buffer.from(text).toString('base64'));
    assert.equal(saved.webFiles.bindings.circularSources[0].gb.resourceId, 'src');
  }

  // A File the committed request does not name is written by the Save, as bytes.
  metrics.length = 0;
  const other = new File(['other bytes'], 'other.tsv', { type: '' });
  const otherState = boundState({ files: {}, circularSources: [card(other)], linearSeqs: [] }, state.inputFiles);
  await serializeActiveRenderFiles(otherState, { ...modeDrawing(state, 'circular'), files: otherState.files, circularSources: otherState.circularSources });
  metrics.length = 0;
  const saved = await savedSession({ inputFiles: state.inputFiles }, committed, bindingDrawings({
    files: { d_color: state.inputFiles.registerFile(other) },
    circularSources: [card(state.inputFiles.registerFile(file))]
  }));
  const otherId = saved.webFiles.bindings.d_color.resourceId;
  assert.notEqual(otherId, 'src');
  assert.deepEqual(saved.resources[otherId], {
    kind: 'web-file', name: `${otherId}-other.tsv`, type: 'application/octet-stream',
    size: 11, lastModified: other.lastModified, encoding: 'base64', data: Buffer.from('other bytes').toString('base64')
  });
  assert.equal(count('base64EncodeCount'), 0);
});

test('Load and two Saves check a slice\'s unmanaged overrides once (OV-1567)', async () => {
  /** @type {unknown[]} */
  const checks = [];
  setUnmanagedConfigOverrideValidator(async (payload) => {
    checks.push(payload);
    return { result: { overrides: payload.configOverrides } };
  });
  // A Session file container (receipt session.file-container): Load reads each input from its entry.
  const sessionFile = await gallerySessionFile('lambda_basic_linear', (session) => {
    session.drawings.find((entry) => entry.id === session.ui.activeDrawingId).config.unmanagedConfigOverrides = { 'objects.gc_content.percent_background_opacity': 0.42 };
  });
  const input = { target: { files: [new File([sessionFile], 'session.gbdraw')], value: 'selected' } };
  const loaded = await importSession(input, { transformLegacyResultSvg, confirmLoadRepairs: async () => true });
  assert.equal(loaded.status, 'ok', JSON.stringify(loaded.error));
  const afterLoad = checks.length;
  for (const save of [1, 2]) {
    const result = await exportSession(`OV-1567 ${save}`, { askChoice: async () => 'download' });
    assert.equal(result?.status, 'saved', JSON.stringify(result?.error ?? result));
    const document = savedDrawingView(/** @type {Record<string, any>} */ ((await readSessionContainer(result.blob)).document), 'linear');
    assert.deepEqual(document.modes.linear.config.unmanagedConfigOverrides,
      { 'objects.gc_content.percent_background_opacity': 0.42 });
  }
  assert.deepEqual({ afterLoad, afterTwoSaves: checks.length }, { afterLoad: 1, afterTwoSaves: 1 });

  // A refusal is not kept: the same request asks the port again.
  checks.length = 0;
  setUnmanagedConfigOverrideValidator(async (payload) => {
    checks.push(payload);
    return checks.length === 1
      ? { result: { error: 'refused' } }
      : { result: { overrides: payload.configOverrides } };
  });
  await assert.rejects(exportSession('refused', { askChoice: async () => 'download' }).then((result) => {
    if (result?.status !== 'saved') throw new Error(String(result?.status));
  }));
  const retried = await exportSession('retried', { askChoice: async () => 'download' });
  assert.equal(retried?.status, 'saved');
  assert.equal(checks.length, 2);
});

// OV-2753, OV-2754: a Session file's input rests in the copy its Load wrote (PD-OI-121), as a chosen File's does.
// Generate stages the copy file itself, under the token the Load hashed while it wrote the copy, so the page reads
// and hashes nothing; the page keeps no bytes once a read of the input ends, so History counts none.
/**
 * Loads the Linear lambda Gallery Session file with a copy store that counts the reads of each copy it wrote.
 * @returns {Promise<{ view: object, copy: { token: string, blob: Blob, reads: number } }>}
 */
const loadWithCountedCopies = async () => {
  /** @type {{ token: string, blob: Blob, reads: number }[]} */
  const copies = [];
  useInputCopyStore({
    openWriter: async () => {
      const pieces = [];
      return {
        write: async (buffer) => { pieces.push(new Uint8Array(buffer).slice()); return buffer; },
        finish: async (token) => {
          const copy = { token, blob: new Blob(pieces), reads: 0 };
          copies.push(copy);
          return { read: async () => { copy.reads += 1; return new Uint8Array(await copy.blob.arrayBuffer()); },
            blob: async () => copy.blob, release: () => {} };
        },
        abort: () => {}
      };
    }
  });
  const sessionFile = await gallerySessionFile('lambda_basic_linear');
  const input = { target: { files: [new File([sessionFile], 'session.gbdraw')], value: 'selected' } };
  const loaded = await importSession(input, { transformLegacyResultSvg, confirmLoadRepairs: async () => true });
  assert.equal(loaded.status, 'ok', JSON.stringify(loaded.error));
  const view = state.inputFiles.fileOf(modeDrawing(state, 'linear').linearSeqs[0].gb);
  const copy = copies.find((entry) => entry.blob.size === view.size);
  assert.ok(copy, 'the Load wrote the GenBank input to a copy');
  return { view, copy };
};

/** Generate's staging of the active Linear drawing's GenBank input, as the diagram Worker receives it. */
const stageLinearGenBank = async () => {
  const linear = modeDrawing(state, 'linear');
  const described = (await serializeActiveRenderFiles({ inputFiles: state.inputFiles, files: linear.files,
    linearSeqs: linear.linearSeqs, circularSources: [] }, linear)).linearSeqs[0].gb;
  return createDiagramResourceTransport().prepare({ request: {
    schema: CANONICAL_REQUEST_SCHEMA, mode: 'linear',
    records: [{ recordKey: 'r1', source: { kind: 'genbank', resourceId: 'src' }, presentation: { label: 'x' } }],
    diagramOptions: { trackSlots: [], annotations: { sets: [] } }, layout: {}, comparisons: [],
    output: { prefix: 'x', formats: ['svg'], overwrite: false }
  }, resources: { src: described } });
};

test('Generate after a Load stages each input as its copy file, reading and hashing nothing; History counts none (OV-2753, OV-2754)', async () => {
  try {
    const { view, copy } = await loadWithCountedCopies();
    assert.equal(residentBytesOf(view), 0, 'History keeps no bytes alive for an input that rests in its copy');
    const text = await readFileText(view);
    const readsBefore = copy.reads;
    assert.equal(await readFileText(view), text);
    assert.equal(copy.reads, readsBefore + 1, 'a read that ended keeps nothing: the next one reads the copy again');

    const readsBeforeGenerate = copy.reads;
    metrics.length = 0;
    const prepared = await stageLinearGenBank();
    assert.equal(prepared.stagedResources[0].blob, copy.blob, 'the staged content is the copy file itself');
    assert.equal(prepared.resourceManifest[0].cacheToken, copy.token, 'under the token the Load hashed');
    assert.deepEqual({ copyReads: copy.reads - readsBeforeGenerate, byteReads: count('resourceByteReadCount') },
      { copyReads: 0, byteReads: 0 }, 'Generate reads no input bytes in the page');
  } finally {
    useInputCopyStore(null);
  }
});

// OV-2756: the readers of one reading (a once-per-File read) share the read of a loaded input's copy, as they share
// a chosen File's, and the page keeps it no longer than the reading.
test('One reading of a loaded Session input reads its copy once and keeps nothing after it (OV-2756)', async () => {
  try {
    const { view, copy } = await loadWithCountedCopies();
    const readsBefore = copy.reads;
    const records = await readOncePerFile([view], 'ov2756', async () => {
      const { records: parsed } = await readParsedSequenceRecords(view, 'genbank');
      return { parsed, listed: await discoverSequenceRecords({ file: view, format: 'genbank' }) };
    });
    assert.ok(records.parsed.length > 0 && records.listed.length > 0);
    assert.equal(copy.reads - readsBefore, 1, 'the parsed records and the record list take one read of the copy');
    const readsAfter = copy.reads;
    await readFileText(view);
    assert.equal(copy.reads, readsAfter + 1, 'the read ended with its reading: the next one reads the copy again');
  } finally {
    useInputCopyStore(null);
  }
});

// A Session loaded where the browser lets no copy be written (no private file system): the input rests in the
// Session file, whose entry may be deflated, so Generate stages a Blob streamed from the entry under the token the
// Load hashed during its checked read of the entry; the page holds no whole read, and History counts none.
test('Generate after a Load with no input copies streams each input from the Session file, reading and hashing nothing whole', async () => {
  useInputCopyStore(null);
  const sessionFile = await gallerySessionFile('lambda_basic_linear');
  const input = { target: { files: [new File([sessionFile], 'session.gbdraw')], value: 'selected' } };
  const loaded = await importSession(input, { transformLegacyResultSvg, confirmLoadRepairs: async () => true });
  assert.equal(loaded.status, 'ok', JSON.stringify(loaded.error));
  const view = state.inputFiles.fileOf(modeDrawing(state, 'linear').linearSeqs[0].gb);
  assert.equal(residentBytesOf(view), 0, 'History keeps no bytes alive for an input that rests in the Session file');
  const text = await readFileText(view);

  metrics.length = 0;
  const prepared = await stageLinearGenBank();
  assert.equal(count('resourceByteReadCount'), 0, 'Generate reads no input bytes whole in the page');
  const staged = new Uint8Array(await prepared.stagedResources[0].blob.arrayBuffer());
  assert.equal(new TextDecoder().decode(staged), text, 'the staged content is the stored input');
  assert.equal(prepared.resourceManifest[0].cacheToken, await contentTokenOfBytes(staged), 'under the token of its bytes');
});
