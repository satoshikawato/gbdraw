// Generate describes an input File without encoding it to base64 (OV-2612): the
// base64 text of a large input exceeds the V8 string limit, and every Generate
// held a 4/3-size copy of each input. A Session writer encodes it, once, and
// refuses a Session over its limit with SESSION_SIZE_LIMIT.
import assert from 'node:assert/strict';
import { test } from 'node:test';
import { webcrypto } from 'node:crypto';
import { File } from 'node:buffer';
import { gzipSync } from 'node:zlib';
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
const { writtenResourcesWith } = await import('./helpers/written-resources.mjs');
// What Save writes for each resource, read as base64 descriptors.
const encodeSessionResources = writtenResourcesWith(await import('../../gbdraw/web/js/services/session-resource-backing.js'));
const { adoptRuntimeCanonicalSession } = await import('../../gbdraw/web/js/services/session-authority.js');
const { inputContentBlob, readFileBytes, readFileText, releaseFileContentCache } = await import('../../gbdraw/web/js/services/file-content-cache.js');
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
  const sourceState = boundState({
    files: {}, linearSeqs: [], circularSources: [{ uid: 'record-1', gb: file, gff: null, fasta: null }]
  }, state.inputFiles);
  const described = await serializeActiveRenderFiles(sourceState, {
    ...modeDrawing(state, 'circular'), files: sourceState.files, circularSources: sourceState.circularSources
  });
  return { sourceState, descriptor: described.circularSources[0].gb };
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
  assert.equal(await prepared.stagedResources[0].blob.text(), '>x\nACGTACGT\n');
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
  assert.equal(await prepared.stagedResources[0].blob.text(), '>s\nACGT\n');
});

test('Save encodes a File Generate described, once, under the committed ID', async () => {
  const text = '>x\nACGTTGCA\n';
  const file = new File([text], 'g.fa', { type: '' });
  const { sourceState, descriptor } = await generateDescriptor(file);
  // Generate stages the File and lets go of its bytes (diagram-resource-staging.js).
  await inputContentBlob(file);
  releaseFileContentCache(file);
  const committed = adoptRuntimeCanonicalSession({ renderRequest: request('src'), resources: { src: descriptor }, webFiles: {} });
  for (const save of [1, 2]) {
    metrics.length = 0;
    const saved = await savedSession(sourceState, committed, bindingDrawings(sourceState));
    assert.deepEqual(Object.keys(saved.resources), ['src'], `Save ${save}`);
    assert.equal(saved.resources.src.encoding, 'base64');
    assert.equal(Buffer.from(saved.resources.src.data, 'base64').toString('utf8'), text);
    assert.equal(saved.webFiles.bindings.circularSources[0].gb.resourceId, 'src');
    assert.equal(count('resourceIdentityHashCount'), 0, `Save ${save}`);
    // Save streams the bytes into the Session file; it never makes a base64 text (receipt session.file-container).
    assert.equal(count('base64EncodeCount'), 0, `Save ${save}`);
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

// Receipt session.file-container: the inputs of a Session file have no cap; only the settings and results do.
test('A Save writes an input larger than the settings-and-results cap, streamed from its File', async () => {
  const text = 'x'.repeat(3000);
  const file = new File([text], 'big.tsv', { type: '' });
  const { descriptor } = await generateDescriptor(file);
  const { sessionContainerResources, writeSessionFile } = await import('../../gbdraw/web/js/services/session-file.js');
  const { openSessionContainerEntry, readSessionContainer } = await import('../../gbdraw/web/js/services/session-container.js');
  hooks.maxExpandedSessionBytes = 2048;
  metrics.length = 0;
  try {
    const document = { format: 'gbdraw-session', resources: sessionContainerResources({ src: descriptor }) };
    const { blob } = await writeSessionFile(document, { keep: null, readResource: async () => backing.sessionResourceStream(descriptor) });
    const container = await readSessionContainer(blob);
    assert.equal(await new Response(openSessionContainerEntry(blob, container.resources.src)).text(), text);
    assert.equal(count('base64EncodeCount'), 0);
  } finally {
    delete hooks.maxExpandedSessionBytes;
  }
});

// Receipt session.file-container: "gzip でない入力は圧縮して入れ、すでに gzip の入力はそのまま入れる".
test('A gzip File is described, sent and saved as its gzip bytes; after Load it reads expanded, nothing expanded kept, and Generate sends the stored bytes', async () => {
  const { openSessionContainerEntry, readSessionContainer } = await import('../../gbdraw/web/js/services/session-container.js');
  const text = `>x\n${'ACGT'.repeat(400)}\n`;
  const gzip = new Uint8Array(gzipSync(text));
  const file = new File([gzip], 'g.fa.gz', { type: '' });
  metrics.length = 0;
  const { descriptor } = await generateDescriptor(file);
  assert.equal(descriptor.size, gzip.byteLength, 'Generate counts the bytes as chosen');
  const sent = await createDiagramResourceTransport().prepare({ request: request('src'), resources: { src: descriptor } });
  assert.deepEqual(new Uint8Array(await sent.stagedResources[0].blob.arrayBuffer()), gzip, 'and sends them: Python expands a gzip input');
  assert.equal(count('inputExpansionCount'), 0, 'describing and sending a gzip File expand nothing');
  const { blob } = await sessionConfig.writeSessionDocumentFile({ format: 'gbdraw-session', resources: { src: descriptor } });
  const container = await readSessionContainer(blob);
  const entry = container.resources.src;
  assert.deepEqual({ method: entry.method, size: entry.size, declared: container.document.resources.src.size },
    { method: 0, size: gzip.byteLength, declared: gzip.byteLength });
  assert.deepEqual(new Uint8Array(await new Response(openSessionContainerEntry(blob, entry)).arrayBuffer()), gzip);

  // Load makes the entry file-backed as stored; the page's readers read it expanded (OV-2834).
  const loaded = backing.fileBackedDescriptor({ kind: 'genbank', name: 'g.fa.gz', type: '', size: gzip.byteLength, lastModified: 0 },
    async () => gzip.slice(), { size: gzip.byteLength, open: async () => openSessionContainerEntry(blob, entry), token: null });
  const view = backing.createSessionResourceFileView(backing.adoptCurrentSessionResources({ src: loaded }), 'src');
  metrics.length = 0;
  assert.equal(await readFileText(view), text);
  assert.equal(await readFileText(view), text);
  assert.equal(Buffer.from(await readFileBytes(view)).toString('utf8'), text);
  assert.equal(count('inputExpansionCount'), 3, 'each read expands again: a gzip view keeps no expanded bytes or text');
  // Generate after Load sends the view's stored bytes with the view's size; Python reads gzip by its signature.
  const { descriptor: viewDescriptor } = await generateDescriptor(view);
  assert.equal(viewDescriptor.size, gzip.byteLength);
  const prepared = await createDiagramResourceTransport().prepare({ request: request('src'), resources: { src: viewDescriptor } });
  assert.deepEqual(new Uint8Array(await prepared.stagedResources[0].blob.arrayBuffer()), gzip);
  // Saved again, the entry is the same gzip bytes.
  const again = await sessionConfig.writeSessionDocumentFile({ format: 'gbdraw-session', resources: { src: viewDescriptor } });
  const reread = await readSessionContainer(again.blob);
  assert.deepEqual({ method: reread.resources.src.method, size: reread.document.resources.src.size }, { method: 0, size: gzip.byteLength });
});

test('Save reads no file-backed input whole: a committed request\'s inputs and same-size candidates compare by stored token', async () => {
  const reads = [];
  const opened = [];
  const tokens = [];
  const loaded = (id, text) => {
    const bytes = new TextEncoder().encode(text);
    return backing.fileBackedDescriptor({ kind: 'genbank', name: `${id}.gb`, type: '', size: bytes.byteLength, lastModified: 0 },
      async () => { reads.push(id); return bytes.slice(); },
      { size: bytes.byteLength, open: async () => { opened.push(id); return new Blob([bytes]).stream(); },
        token: async () => { tokens.push(id); return `token-${text}`; } });
  };
  // Two inputs of one size (a Load's copies), in a request that is not adopted: a later drawing's request.
  const committed = { renderRequest: { ...request('a'), records: [
    { recordKey: 'r1', source: { kind: 'genbank', resourceId: 'a' }, presentation: { label: 'x' } },
    { recordKey: 'r2', source: { kind: 'genbank', resourceId: 'b' }, presentation: { label: 'y' } }
  ] }, resources: { a: loaded('a', 'LOCUS a\n'), b: loaded('b', 'LOCUS b\n') }, webFiles: {} };
  metrics.length = 0;
  const built = await buildSessionResources({ inputFiles: state.inputFiles }, committed, bindingDrawings({ files: {} }));
  assert.equal(Object.keys(built.resources).length, 2);
  assert.deepEqual({ reads, opened, hashes: count('resourceIdentityHashCount') }, { reads: [], opened: [], hashes: 0 });
  assert.deepEqual(tokens.sort(), ['a', 'b'], 'equal sizes compare by the tokens their copies were written under');
});
