// Generate describes an input File without encoding it to base64 (OV-2612): the
// base64 text of a large input exceeds the V8 string limit, and every Generate
// held a 4/3-size copy of each input. A Session writer encodes it, once, and
// refuses a Session over its limit with SESSION_SIZE_LIMIT.
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

const metrics = [];
const hooks = { onStructuralMetric: (metric) => metrics.push(metric) };
globalThis.__GBDRAW_TEST_HOOKS__ = hooks;
const count = (name) => metrics.filter((metric) => metric.name === name).length;

const sessionConfig = await import('../../gbdraw/web/js/services/config.js');
const { buildSessionResources } = await import('../../gbdraw/web/js/services/session-resources.js');
const { encodeSessionResources } = await import('../../gbdraw/web/js/services/session-file.js');
const { adoptRuntimeCanonicalSession } = await import('../../gbdraw/web/js/services/session-authority.js');
const { takeFileBytesForTransfer } = await import('../../gbdraw/web/js/services/file-content-cache.js');
const { createDiagramResourceTransport } = await import('../../gbdraw/web/js/services/diagram-resource-staging.js');
const backing = await import('../../gbdraw/web/js/services/session-resource-backing.js');
const { CANONICAL_REQUEST_SCHEMA } = await import('../../gbdraw/web/js/contracts/schema.generated.js');
const { state } = await import('../../gbdraw/web/js/state.js');
const { serializeActiveRenderFiles } = callerSession(sessionConfig, state, await sessionRollbackPort(state));

const request = (resourceId) => ({
  schema: CANONICAL_REQUEST_SCHEMA, mode: 'circular',
  records: [{ recordKey: 'r1', source: { kind: 'genbank', resourceId }, presentation: { label: 'x' } }],
  diagramOptions: { trackSlots: [], annotations: { sets: [] } }, layout: {}, comparisons: [],
  output: { prefix: 'x', formats: ['svg'], overwrite: false }
});

/** Save's resource table: allocated, then encoded by the Session file writer. */
const savedSession = async (sourceState, committed, drawings) => {
  const built = await buildSessionResources(sourceState, committed, drawings);
  return { ...built, resources: await encodeSessionResources(built.resources) };
};

/** Generate's description of a File in a Circular drawing. */
const generateDescriptor = async (file) => {
  const sourceState = boundState({ files: { c_gb: file }, linearSeqs: [] }, state.inputFiles);
  const described = await serializeActiveRenderFiles(sourceState, { ...modeDrawing(state, 'circular'), files: sourceState.files });
  return { sourceState, descriptor: described.c_gb };
};

test('Generate describes an input File with no base64 text (OV-2612)', async () => {
  metrics.length = 0;
  const file = new File(['>x\nACGT\n'], 'g.fa', { type: '', lastModified: 5 });
  const { descriptor } = await generateDescriptor(file);
  assert.equal(Object.hasOwn(descriptor, 'data'), false);
  assert.equal(Object.hasOwn(descriptor, 'encoding'), false);
  assert.deepEqual({ ...descriptor }, { name: 'g.fa', type: '', size: 8, lastModified: 5 });
  assert.equal(count('base64EncodeCount'), 0);
});

test('The Worker receives the bytes of a file-backed descriptor', async () => {
  const file = new File(['>x\nACGTACGT\n'], 'g.fa', { type: '' });
  const { descriptor } = await generateDescriptor(file);
  const prepared = await createDiagramResourceTransport().prepare({ request: request('src'), resources: { src: descriptor } });
  assert.equal(prepared.stagedResources.length, 1);
  assert.equal(Buffer.from(prepared.stagedResources[0].bytes).toString('utf8'), '>x\nACGTACGT\n');
});

test('A second Generate describes a view of a file-backed resource under its shown name with its bytes', async () => {
  // After Generate 1, LOSAT conservation binds each FASTA as a view over the committed,
  // file-backed resource, shown under the name the user chose.
  const bytes = new TextEncoder().encode('>s\nACGT\n');
  const id = 'conservation-losat-fasta-files-1';
  const committed = backing.fileBackedDescriptor({ kind: 'conservation-fasta-file', name: `${id}-phage2.fasta`,
    type: '', size: bytes.byteLength, lastModified: 0 }, async () => bytes);
  const view = backing.createSessionResourceFileView(backing.adoptCurrentSessionResources({ [id]: committed }), id, { name: 'phage2.fasta' });
  const { descriptor } = await generateDescriptor(view);
  assert.equal(descriptor.name, 'phage2.fasta');
  assert.ok(backing.fileBackedReaderOf(descriptor), 'the renamed descriptor keeps its reader');
  const adopted = adoptRuntimeCanonicalSession({ renderRequest: request('src'), resources: { src: descriptor }, webFiles: {} });
  const prepared = await createDiagramResourceTransport().prepare({ request: request('src'), resources: adopted.resources });
  assert.equal(Buffer.from(prepared.stagedResources[0].bytes).toString('utf8'), '>s\nACGT\n');
});

test('A combined Session view is the source of the file-backed descriptor Generate made of it', async () => {
  // A composite GenBank card (two Session resources read as one file): after Generate, its
  // features are the displayed Result's only while the view matches the committed descriptor.
  const encoded = (text) => ({ kind: 'genbank', name: 'part.gb', type: 'text/plain', size: Buffer.byteLength(text),
    lastModified: 0, encoding: 'base64', data: Buffer.from(text).toString('base64') });
  const table = backing.adoptCurrentSessionResources({ a: encoded('LOCUS a\n//\n'), b: encoded('LOCUS b\n//\n') });
  const view = backing.createCombinedSessionResourceFileView(table, [
    { resourceId: 'a', name: 'a.gb', type: 'text/plain', lastModified: 0 },
    { resourceId: 'b', name: 'b.gb', type: 'text/plain', lastModified: 0 }
  ], { name: 'both.gb', type: 'text/plain', lastModified: 0 });
  const { descriptor } = await generateDescriptor(view);
  assert.equal(Object.hasOwn(descriptor, 'data'), false);
  assert.equal(backing.matchesSessionResourceDescriptor(view, descriptor), true);
  assert.equal(backing.matchesSessionResourceDescriptor(view, { ...encoded('LOCUS c\n//\n') }), false);
});

test('Save encodes a File Generate described, once, under the committed ID', async () => {
  const text = '>x\nACGTTGCA\n';
  const file = new File([text], 'g.fa', { type: '' });
  const { sourceState, descriptor } = await generateDescriptor(file);
  await takeFileBytesForTransfer(file);
  const committed = adoptRuntimeCanonicalSession({ renderRequest: request('src'), resources: { src: descriptor }, webFiles: {} });
  for (const save of [1, 2]) {
    metrics.length = 0;
    const saved = await savedSession(sourceState, committed, bindingDrawings(sourceState));
    assert.deepEqual(Object.keys(saved.resources), ['src'], `Save ${save}`);
    assert.equal(saved.resources.src.encoding, 'base64');
    assert.equal(Buffer.from(saved.resources.src.data, 'base64').toString('utf8'), text);
    assert.equal(saved.webFiles.bindings.c_gb.resourceId, 'src');
    assert.equal(count('resourceIdentityHashCount'), 0, `Save ${save}`);
    assert.equal(count('base64EncodeCount'), 1, `Save ${save}`);
  }
  // The committed request keeps its file-backed descriptor: Save writes a copy.
  assert.equal(Object.hasOwn(committed.resources.src, 'data'), false);
});

test('Two Files with equal metadata stay two resources in a Save', async () => {
  const first = new File(['AAAA'], 'same.tsv', { type: '', lastModified: 1 });
  const second = new File(['CCCC'], 'same.tsv', { type: '', lastModified: 1 });
  const one = await generateDescriptor(first);
  const two = await generateDescriptor(second);
  const saved = await encodeSessionResources({ a: one.descriptor, b: two.descriptor });
  assert.deepEqual(['a', 'b'].map((id) => Buffer.from(saved.a === saved[id] ? saved.a.data : saved[id].data, 'base64').toString()),
    ['AAAA', 'CCCC']);
  const committed = adoptRuntimeCanonicalSession({ renderRequest: { ...request('a'), records: [
    { recordKey: 'r1', source: { kind: 'genbank', resourceId: 'a' }, presentation: { label: 'x' } },
    { recordKey: 'r2', source: { kind: 'genbank', resourceId: 'b' }, presentation: { label: 'y' } }
  ] }, resources: { a: one.descriptor, b: two.descriptor }, webFiles: {} });
  const written = await savedSession({ inputFiles: state.inputFiles }, committed, bindingDrawings({ files: {} }));
  assert.deepEqual(Object.keys(written.resources).sort(), ['a', 'b']);
  assert.notEqual(written.resources.a.data, written.resources.b.data);
});

test('A Save over the Session limit fails with SESSION_SIZE_LIMIT before it encodes', async () => {
  const file = new File(['x'.repeat(3000)], 'big.tsv', { type: '' });
  const { sourceState, descriptor } = await generateDescriptor(file);
  const committed = adoptRuntimeCanonicalSession({ renderRequest: request('src'), resources: { src: descriptor }, webFiles: {} });
  hooks.maxExpandedSessionBytes = 2048;
  metrics.length = 0;
  try {
    await assert.rejects(savedSession(sourceState, committed, bindingDrawings(sourceState)),
      (error) => error.code === 'SESSION_SIZE_LIMIT' && error.operation === 'session-save'
        && error.context.sessionByteLimit === 2048 && error.context.largestInputs === 'big.tsv');
    assert.equal(count('base64EncodeCount'), 0);
  } finally {
    delete hooks.maxExpandedSessionBytes;
  }
});
