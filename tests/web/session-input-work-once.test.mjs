// Save and Load do each piece of input work once (OV-1559, OV-1567): a File
// Generate described is bound by its descriptor with no hash or decode and
// encoded once per Save (OV-2612),
// and the Worker check of a slice's unmanaged overrides runs once per distinct
// request across Load and every Save.
import { savedDrawingView } from './helpers/saved-drawing-view.mjs';
import assert from 'node:assert/strict';
import { test } from 'node:test';
import { readFileSync } from 'node:fs';
import { gunzipSync } from 'node:zlib';
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
const { encodeSessionResources } = await import('../../gbdraw/web/js/services/session-file.js');
/** Save's resource table: allocated, then encoded by the Session file writer. */
const savedSession = async (...args) => {
  const built = await buildSessionResources(...args);
  return { ...built, resources: await encodeSessionResources(built.resources) };
};
const { adoptRuntimeCanonicalSession } = await import('../../gbdraw/web/js/services/session-authority.js');
const { holdReadBytes, takeFileBytesForTransfer } = await import('../../gbdraw/web/js/services/file-content-cache.js');
const { CANONICAL_REQUEST_SCHEMA } = await import('../../gbdraw/web/js/contracts/schema.generated.js');
const { transformLegacyResultSvg } = await import('../../gbdraw/web/js/app/app-setup.js');
const { state } = await import('../../gbdraw/web/js/state.js');
const { exportSession, importSession, serializeActiveRenderFiles, setUnmanagedConfigOverrideValidator } = callerSession(sessionConfig, state, await sessionRollbackPort(state));

test('Save binds a File Generate described with no hash or decode and encodes it once (OV-1559, OV-2612)', async () => {
  const text = `>x\n${'ACGT'.repeat(250_000)}\n`;
  const file = new File([text], 'g.fa', { type: '' });
  let reads = 0;
  const read = file.arrayBuffer.bind(file);
  file.arrayBuffer = async () => { reads += 1; return read(); };
  // Generate: the request describes the File and stages its bytes (the transfer releases the byte cache).
  const sourceState = boundState({ files: {}, circularSources: [card(file)], linearSeqs: [] });
  const releaseReads = holdReadBytes();
  const described = (await serializeActiveRenderFiles(sourceState, { ...modeDrawing(state, 'circular'), files: sourceState.files, circularSources: sourceState.circularSources })).circularSources[0].gb;
  await takeFileBytesForTransfer(file);
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
  // Each Save reads the File once: the page keeps no input bytes between Saves.
  for (const save of [1, 2]) {
    metrics.length = 0;
    const saved = await savedSession(sourceState, committed, bindingDrawings(sourceState));
    assert.deepEqual({
      reads,
      hashes: count('resourceIdentityHashCount'),
      decodes: count('base64DecodeCount'),
      byteReads: count('resourceByteReadCount'),
      encodes: count('base64EncodeCount')
    }, { reads: 1 + save, hashes: 0, decodes: 0, byteReads: 0, encodes: 1 }, `Save ${save}`);
    assert.deepEqual(Object.keys(saved.resources), ['src']);
    assert.equal(saved.resources.src.data, Buffer.from(text).toString('base64'));
    assert.equal(saved.webFiles.bindings.circularSources[0].gb.resourceId, 'src');
  }

  // A File the committed request does not name is encoded by the Save.
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
  assert.equal(count('base64EncodeCount'), 2);
});

test('Load and two Saves check a slice\'s unmanaged overrides once (OV-1567)', async () => {
  /** @type {unknown[]} */
  const checks = [];
  setUnmanagedConfigOverrideValidator(async (payload) => {
    checks.push(payload);
    return { result: { overrides: payload.configOverrides } };
  });
  const session = JSON.parse(readFileSync(new URL('../../gbdraw/web/gallery/sessions/lambda_basic_linear.gbdraw-session.json', import.meta.url), 'utf8'));
  session.drawings.find((entry) => entry.id === session.ui.activeDrawingId).config.unmanagedConfigOverrides = { 'objects.gc_content.percent_background_opacity': 0.42 };
  const input = { target: { files: [new File([JSON.stringify(session)], 'session.json', { type: 'application/json' })], value: 'selected' } };
  const loaded = await importSession(input, { transformLegacyResultSvg, confirmLoadRepairs: async () => true });
  assert.equal(loaded.status, 'ok', JSON.stringify(loaded.error));
  const afterLoad = checks.length;
  for (const save of [1, 2]) {
    const result = await exportSession(`OV-1567 ${save}`, { askChoice: async () => 'download' });
    assert.equal(result?.status, 'saved', JSON.stringify(result?.error ?? result));
    const document = savedDrawingView(JSON.parse(gunzipSync(Buffer.from(await result.blob.arrayBuffer()))), 'linear');
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
