import assert from 'node:assert/strict';

globalThis.location = { href: 'https://example.test/gbdraw/web/' };
delete globalThis.loadPyodide;

const { EXPECTED_WEB_RUNTIME_CAPABILITIES } = await import(
  '../../gbdraw/web/js/services/runtime-capabilities.js'
);

class FakeWorker {
  static instances = [];
  static deferHelpers = false;
  static failRuns = false;

  constructor() {
    this.listeners = new Map();
    this.messages = [];
    this.transferLists = [];
    this.terminated = false;
    FakeWorker.instances.push(this);
  }

  addEventListener(type, listener) {
    const listeners = this.listeners.get(type) || new Set();
    listeners.add(listener);
    this.listeners.set(type, listeners);
  }

  removeEventListener(type, listener) {
    this.listeners.get(type)?.delete(listener);
  }

  emit(data) {
    this.listeners.get('message')?.forEach((listener) => listener({ data }));
  }

  postMessage(message, transfer = []) {
    this.messages.push(message);
    this.transferLists.push(transfer);
    if (message.type === 'init') {
      queueMicrotask(() => this.emit({
        id: message.id,
        type: 'init',
        ok: true,
        capabilities: structuredClone(EXPECTED_WEB_RUNTIME_CAPABILITIES)
      }));
      return;
    }
    if (message.type === 'helper') {
      if (FakeWorker.deferHelpers) return;
      queueMicrotask(() => this.emit({
        requestId: message.requestId,
        type: 'helper',
        ok: true,
        result: {
          records: [{ selector: '#1', record_id: 'shared-worker', record_length: 4 }]
        }
      }));
      return;
    }
    if (message.type === 'run') {
      if (FakeWorker.failRuns) {
        queueMicrotask(() => this.emit({
          requestId: message.requestId,
          type: 'run',
          ok: false,
          error: {
            name: 'RenderStageError',
            message: 'diagram render stage failed',
            details: [{ label: 'Stage', text: 'typed render request' }]
          }
        }));
        return;
      }
      queueMicrotask(() => this.emit({
        requestId: message.requestId,
        type: 'run',
        ok: true,
        results: {
          results: [{ name: 'shared.svg', content: '<svg></svg>' }],
          metadata: {}
        }
      }));
    }
  }

  terminate() {
    this.terminated = true;
  }
}

globalThis.Worker = FakeWorker;

const {
  DIAGRAM_HELPER_OPERATIONS,
  disposeDiagramGenerationWorker,
  runDiagramGeneration,
  runDiagramHelperOperation
} = await import('../../gbdraw/web/js/services/diagram-generation.js');

await assert.rejects(
  runDiagramHelperOperation('runArbitraryPython', {}),
  /Unsupported diagram helper operation/
);
assert.equal(FakeWorker.instances.length, 0);
const sourceBuffer = new TextEncoder().encode('>shared-worker\nACGT\n').buffer;
const helper = await runDiagramHelperOperation(
  DIAGRAM_HELPER_OPERATIONS.LIST_SEQUENCE_RECORDS,
  {
    format: 'fasta',
    files: [{ role: 'source', bytes: sourceBuffer }]
  }
);
assert.equal(helper.result.records[0].record_id, 'shared-worker');
const helperWorker = FakeWorker.instances[0];
assert.equal(FakeWorker.instances.length, 1);
assert.equal(helperWorker.messages.filter(({ type }) => type === 'init').length, 1);
assert.equal(helperWorker.messages.filter(({ type }) => type === 'helper').length, 1);
assert.equal(helperWorker.messages.filter(({ type }) => type === 'run').length, 0);
const rendered = await runDiagramGeneration({ request: {}, resources: {} });
assert.equal(rendered.results[0].name, 'shared.svg');

const worker = FakeWorker.instances[0];
assert.equal(FakeWorker.instances.length, 1);
assert.equal(worker.messages.filter(({ type }) => type === 'init').length, 1);
assert.equal(worker.messages.filter(({ type }) => type === 'helper').length, 1);
assert.equal(worker.messages.filter(({ type }) => type === 'run').length, 1);
assert.equal(worker.transferLists[1][0], sourceBuffer);

await runDiagramHelperOperation(
  DIAGRAM_HELPER_OPERATIONS.HYDRATE_PROTEIN_LOSAT_TSV,
  { entry: {}, identityManifest: {} }
);
await runDiagramGeneration({ request: {}, resources: {} });
assert.equal(FakeWorker.instances.length, 1);
assert.equal(worker.messages.filter(({ type }) => type === 'init').length, 1);
FakeWorker.deferHelpers = true;
const pendingHelperA = runDiagramHelperOperation(
  DIAGRAM_HELPER_OPERATIONS.HYDRATE_PROTEIN_LOSAT_TSV,
  { entry: {}, identityManifest: {} }
);
const pendingHelperB = runDiagramHelperOperation(
  DIAGRAM_HELPER_OPERATIONS.MEASURE_LEGEND_TEXT,
  { caption: 'pending', fontFamily: 'Arial', fontSize: 14 }
);
const rejectedA = assert.rejects(pendingHelperA, /disposed/);
const rejectedB = assert.rejects(pendingHelperB, /disposed/);
await Promise.resolve();
await Promise.resolve();
disposeDiagramGenerationWorker();
await Promise.all([rejectedA, rejectedB]);
assert.equal(worker.terminated, true);

FakeWorker.deferHelpers = false;
const directRendered = await runDiagramGeneration({ request: {}, resources: {} });
assert.equal(directRendered.results[0].name, 'shared.svg');
const directWorker = FakeWorker.instances[1];
assert.equal(directWorker.messages.filter(({ type }) => type === 'init').length, 1);
assert.equal(directWorker.messages.filter(({ type }) => type === 'helper').length, 0);
assert.equal(directWorker.messages.filter(({ type }) => type === 'run').length, 1);
disposeDiagramGenerationWorker();

FakeWorker.failRuns = true;
await assert.rejects(
  runDiagramGeneration({ request: {}, resources: {} }),
  (error) => {
    assert.equal(error.name, 'RenderStageError');
    assert.match(error.message, /render stage failed/);
    assert.deepEqual(error.details, [{ label: 'Stage', text: 'typed render request' }]);
    return true;
  }
);
disposeDiagramGenerationWorker();
FakeWorker.failRuns = false;

const { cloneFileBytesForTransfer } = await import(
  '../../gbdraw/web/js/services/file-content-cache.js'
);
let fileReads = 0;
const cachedFile = {
  arrayBuffer: async () => {
    fileReads += 1;
    return new Uint8Array([1, 2, 3, 4]).buffer;
  }
};
const transferredClone = await cloneFileBytesForTransfer(cachedFile);
structuredClone(transferredClone, { transfer: [transferredClone] });
assert.equal(transferredClone.byteLength, 0);
assert.deepEqual(
  Array.from(new Uint8Array(await cloneFileBytesForTransfer(cachedFile))),
  [1, 2, 3, 4]
);
assert.equal(fileReads, 1);

const { createPaletteLoader } = await import('../../gbdraw/web/js/app/palettes.js');
const ref = (value) => ({ value });
const paletteState = {
  paletteDefinitions: ref({}),
  paletteNames: ref([]),
  selectedPalette: ref('default'),
  currentColors: ref({}),
  appliedPaletteName: ref('default'),
  appliedPaletteColors: ref({}),
  pendingPaletteName: ref(''),
  pendingPaletteColors: ref({}),
  normalizePaletteColors: (value) => value,
  normalizePaletteDefinitions: (value) => value
};
globalThis.fetch = async () => ({
  ok: true,
  status: 200,
  json: async () => ({
    palettes: {
      default: {
        CDS: '#54bcf8',
        repeat_region: '#d3d3d3'
      }
    }
  })
});
await createPaletteLoader({ state: paletteState }).loadPaletteAsset();
assert.equal(globalThis.loadPyodide, undefined);
assert.equal(paletteState.currentColors.value.repeat_region, '#d3d3d3');

console.log('single Worker Pyodide startup tests passed');
