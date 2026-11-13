import { modeDrawing } from './helpers/drawing-state.mjs';
import assert from 'node:assert/strict';
import { webcrypto } from 'node:crypto';
import { spawnSync } from 'node:child_process';
import test from 'node:test';
import { readFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { LOSAT_WASM_BUILD } from '../../gbdraw/web/js/services/losat-build.generated.js';
import { installFakeSvgDom } from './fake-svg-dom.mjs';
import { savedDrawingView } from './helpers/saved-drawing-view.mjs';
import { pythonLinearUnitPlan } from './helpers/linear-unit-plan.mjs';

globalThis.window = {
  location: { href: 'https://audit.invalid/' },
  Vue: {
    ref: (value) => ({ value }),
    reactive: (value) => value,
    computed: (getter) => ({ get value() { return getter(); } }),
    nextTick: async () => {}
  },
  DOMPurify: { sanitize: (value) => value }
};
let downloadedBlob = null;
let downloadedFilename = '';
// Every file name a download saved, in order (Save Raw LOSAT TSV saves one file per result).
const downloadedFilenames = [];
globalThis.document = {
  body: { appendChild: (node) => { node.parentNode = globalThis.document.body; } },
  createElement: () => ({
    addEventListener() {},
    click() { downloadedFilename = this.download; downloadedFilenames.push(this.download); },
    parentNode: null
  }),
  removeChild() {}
};
URL.createObjectURL = (blob) => {
  downloadedBlob = blob;
  return 'blob:audit-download';
};
URL.revokeObjectURL = () => {};
installFakeSvgDom();
globalThis.alert = () => {};
globalThis.confirm = () => true;
if (!globalThis.crypto) globalThis.crypto = webcrypto;

let activePrimaryReads = 0;
let activePrimaryStreams = 0;
let inactiveFileReads = 0;

/** A copy of a canonical Session whose resource descriptors stay shared, as `projectCommittedRecordTransform` copies one. */
const cloneSharingResources = ({ resources, ...rest }) => ({ ...structuredClone(rest), resources: { ...resources } });

class AuditFile extends Blob {
  constructor(parts, name, options = {}) {
    super(parts, options);
    this.name = String(name);
    this.lastModified = options.lastModified ?? 1;
  }

  async arrayBuffer() {
    activePrimaryReads += 1;
    return super.arrayBuffer();
  }

  stream() {
    activePrimaryStreams += 1;
    return super.stream();
  }
}

class AuditInactiveFile extends AuditFile {
  async arrayBuffer() {
    inactiveFileReads += 1;
    throw new Error(`inactive file was read: ${this.name}`);
  }
}

class AuditCacheMap extends Map {
  constructor(entries = []) {
    super(entries);
    this.probes = 0;
  }

  get(key) {
    this.probes += 1;
    return super.get(key);
  }

  has(key) {
    this.probes += 1;
    return super.has(key);
  }
}

globalThis.File = AuditFile;

const { EXPECTED_WEB_RUNTIME_CAPABILITIES } = await import(
  '../../gbdraw/web/js/services/runtime-capabilities.js'
);

const workerResponses = [];
const workerHelperResponses = [];
/** The last `planLinearUnits` reply, Python's. @type {{ units: { recordKey: string }[] } | null} */
let lastLinearUnitPlan = null;
const workerMessages = [];

/** The bytes of every resource the page staged, by content token: a held resource is not staged again. */
const stagedBytes = new Map();
const stagedBytesByToken = () => stagedBytes;
/** Reads the Blobs a message stages, as the Worker does before it answers (OV-2749). @param {any} message */
const readStagedBlobs = (message) => Promise.all((message.payload?.stagedResources || []).map(async ({ cacheToken, blob }) => {
  stagedBytes.set(cacheToken, new Uint8Array(await blob.arrayBuffer()));
}));
/** The latest ring table conversion message: its request, ring table IDs, and staged resources (base64 by ID). */
const ringTableConversionInput = () => {
  const { canonicalRequest, resourceIds, resourceManifest } = workerMessages
    .findLast(({ type, operation }) => type === 'helper' && operation === 'convertRingTablesToSearchFrame').payload;
  const bytes = stagedBytesByToken();
  const resources = Object.fromEntries(resourceManifest.map(({ resourceId, cacheToken }) => [
    resourceId, Buffer.from(bytes.get(cacheToken)).toString('base64')
  ]));
  return { canonicalRequest, resourceIds, resources };
};
/** The Worker's reply to a ring table conversion whose rows are already in the search frame. */
const ringTablesUnchanged = {
  then(resolve, reject) {
    try {
      const { resourceIds, resources } = ringTableConversionInput();
      resolve({ ok: true, result: { resources: Object.fromEntries(resourceIds.map((resourceId) => [
        resourceId, { data: resources[resourceId], size: Buffer.from(resources[resourceId], 'base64').length }
      ])) } });
    } catch (error) {
      reject(error);
    }
  }
};

const assertDownloadedFilesMatchCommands = (runInfo, entries) => {
  const expected = runInfo.helperFiles.map(({ name }) => name).sort();
  const downloaded = entries.map(({ name }) => name).sort();
  assert.deepEqual(downloaded, expected);
  const expectedSet = new Set(expected);
  const tableDependencies = entries.flatMap(({ text }) => {
    const rows = String(text || '').trim().split('\n').map((line) => line.split('\t'));
    const blastIndex = rows[0]?.indexOf('blast') ?? -1;
    return blastIndex < 0 ? [] : rows.slice(1).map((row) => row[blastIndex]).filter(Boolean);
  });
  const referenced = Array.from(new Set([
    ...(runInfo.sourceRecipe?.commandArgs || []),
    ...(runInfo.exactReplay?.commandArgs || []),
    ...tableDependencies
  ].filter((token) => expectedSet.has(token)))).sort();
  assert.deepEqual(referenced, expected);
};

// The Worker sends the feature catalog beside the metadata with its revision
// (step 3a); the fixtures below keep the catalog in `metadata` for the
// generationExecutor seam. A fixture without a catalog stands for one the page
// cannot admit, as a missing catalog was before. A payload that names its own
// `featureCatalogRevision` is sent as written (the "unchanged" cases).
let fakeCatalogBuilds = 0;
const fakeCatalogRevision = (revision) => ({ revision, sha256: '0'.repeat(64), bytes: 1 });
const asWorkerReply = (payload) => {
  if (!payload || typeof payload !== 'object' || Array.isArray(payload) || payload.error
    || Object.hasOwn(payload, 'featureCatalogRevision')) return payload;
  const { featureCatalog = {}, ...metadata } = payload.metadata || {};
  fakeCatalogBuilds += 1;
  return {
    ...payload,
    metadata,
    featureCatalog,
    featureCatalogRevision: fakeCatalogRevision(`fake:${fakeCatalogBuilds}`)
  };
};

class AuditSimplePathWorker {
  constructor() {
    this.listeners = new Map();
  }

  addEventListener(type, listener) {
    const listeners = this.listeners.get(type) || new Set();
    listeners.add(listener);
    this.listeners.set(type, listeners);
  }

  removeEventListener(type, listener) {
    this.listeners.get(type)?.delete(listener);
  }

  emit(type, data) {
    for (const listener of this.listeners.get(type) || []) {
      listener({ data });
    }
  }

  postMessage(message) {
    workerMessages.push(structuredClone(message));
    if (message.type === 'init') {
      queueMicrotask(() => this.emit('message', {
        type: 'init',
        id: message.id,
        ok: true,
        capabilities: structuredClone(EXPECTED_WEB_RUNTIME_CAPABILITIES)
      }));
      return;
    }
    if (message.type === 'run') {
      const response = workerResponses.shift();
      readStagedBlobs(message).then(() => response).then((payload) => queueMicrotask(() => this.emit('message', {
        type: 'run',
        requestId: message.requestId,
        ok: true,
        results: asWorkerReply(structuredClone(payload))
      })));
      return;
    }
    if (message.type === 'helper') {
      // The unit plan is Python's (`planLinearUnits`); a test queues no reply for it.
      if (message.operation === 'planLinearUnits') lastLinearUnitPlan = pythonLinearUnitPlan(message.payload);
      const response = message.operation === 'planLinearUnits'
        ? { ok: true, result: lastLinearUnitPlan } : workerHelperResponses.shift();
      if (!response) throw new Error(`missing helper worker response for ${message.operation}`);
      readStagedBlobs(message).then(() => response).then((payload) => this.emit('message', {
        type: 'helper',
        requestId: message.requestId,
        ...structuredClone(payload)
      }));
      return;
    }
    throw new Error(`unexpected worker message: ${String(message.type)}`);
  }

  terminate() {}
}

globalThis.Worker = AuditSimplePathWorker;

/**
 * The E-value scale port of these Generates: every unit of the last unit plan searched as a whole genome (no scale).
 * Python's scales for a cropped or windowed drawing are pinned by tests/test_genome_search_space.py and the
 * Playwright LOSAT cases; these tests read the args of whole-genome searches.
 */
const planWholeGenomeScales = async () => {
  if (!lastLinearUnitPlan) throw new Error('planGenomeScales ran before the unit plan');
  const recordKeys = lastLinearUnitPlan.units.map(({ recordKey }) => recordKey);
  return { recordKeys, scales: recordKeys.map(() => recordKeys.map(() => null)), warnings: [] };
};

const {
  afterFrame,
  afterPaint,
  afterTask,
  buildLosatCachePayload,
  createRunAnalysis,
  executeCanonicalRenderCandidate,
  requireLegacyProteinReferencesMigrate,
  showDiagramProgress
} = await import('../../gbdraw/web/js/app/run-analysis.js');
const {
  linearComparisonEdgeKey,
  resolveLinearComparisonPlan
} = await import('../../gbdraw/web/js/services/linear-comparisons.js');
const { cardRecordFacts } = await import('../../gbdraw/web/js/services/linear-card-windows.js');
const {
  applyEditorStateData,
  applyFeatureStateData,
  applyOrthogroupStateData,
  applyResultsData,
  applyRunStateData,
  applyUiStateData,
  buildEditorStateData,
  buildFeatureStateData,
  buildOrthogroupStateData,
  buildRunStateData,
  buildUiStateData,
  serializeActiveRenderFiles,
  serializeFile,
  SESSION_VERSION,
  installShownDrawing
} = await import('../../gbdraw/web/js/services/config.js');
const { CANONICAL_REQUEST_SCHEMA } = await import('../../gbdraw/web/js/services/session-request.js');
const { writeSessionDocumentFile } = await import('../../gbdraw/web/js/services/config.js');
const { readSessionContainer } = await import('../../gbdraw/web/js/services/session-container.js');
const { createInputFileStore } = await import(
  '../../gbdraw/web/js/services/history-files.js'
);
const { createHistoryManager } = await import('../../gbdraw/web/js/services/history.js');
const { createHistorySnapshotService } = await import(
  '../../gbdraw/web/js/services/history-snapshot.js'
);
const {
  adoptCurrentSessionResources,
  createSessionResourceFileView
} = await import('../../gbdraw/web/js/services/session-resource-backing.js');
const {
  featureStateFromCatalog
} = await import('../../gbdraw/web/js/services/feature-catalog.js');
const { createDrawingState, state } = await import('../../gbdraw/web/js/state.js');
await (await import('./helpers/drawing-state.mjs')).addStartLinearDrawing();
const { registerInputFileTree } = await import('../../gbdraw/web/js/services/drawing/inputs.js');
// A binding holds refs of the project's file store, as the app's writers bind them.
const bound = (value) => registerInputFileTree(state.inputFiles, value);
const BOUND_FIELDS = ['c_depth', 'c_conservation_blasts', 'c_conservation_fastas', 'c_conservation_sequence_sources'];
const { createCircularSource, replaceCircularSources } = await import('../../gbdraw/web/js/services/drawing/inputs.js');
/** The Circular drawing's File cards, one per `{ gb, gff, fasta }` given (each file bound as the app binds it). */
const setCircularCardsOf = (drawing, ...cards) => replaceCircularSources(drawing.slices.inputs, cards.map((card, index) => (
  createCircularSource({ uid: `record-${index + 1}`, gb: bound(card.gb ?? null), gff: bound(card.gff ?? null), fasta: bound(card.fasta ?? null) }))));
/** A record of card `record-1` as the record list holds it (`k` of `count` records in the card's file). */
const listedRecord = (recordId, recordLength = 10, k = 1, count = 1) => ({
  uid: 'record-1', k, identity: count > 1 ? `record-1:${k}` : 'record-1', selector: `#${k}`,
  record_id: recordId, record_length: recordLength, detectedTopology: 'circular'
});
const setCircularCards = (...cards) => setCircularCardsOf(modeDrawing(state, 'circular'), ...cards);
/** One card holding a GenBank file, or no card for `null`. */
const setCircularFile = (file, drawing = modeDrawing(state, 'circular')) => setCircularCardsOf(drawing, ...(file ? [{ gb: file }] : []));
const boundFiles = (fields) => Object.fromEntries(Object.entries(fields)
  .map(([key, value]) => [key, BOUND_FIELDS.includes(key) ? bound(value) : value]));
const {
  readComparisonSequence, readComparisonSequenceDigest, readOncePerFile
} = await import('../../gbdraw/web/js/services/record-discovery.js');

const artifactSnapshots = createHistorySnapshotService({
  state,
  fileStore: createInputFileStore(),
  nextTick: window.Vue.nextTick,
  buildUiStateData,
  applyUiStateData,
  buildFeatureStateData,
  applyFeatureStateData,
  buildEditorStateData,
  applyEditorStateData,
  buildOrthogroupStateData,
  applyOrthogroupStateData,
  currentResults: () => state.results.value,
  applyResultsData,
  buildRunStateData,
  applyRunStateData
});
const createImmediatePreviewRuntime = ({ onEvent = () => {} } = {}) => {
  let activeExpectation = null;
  let activeReceipt = null;
  let nextBindSequence = 1;
  return {
    registerReadinessExpectation(options) {
      const bindSequence = nextBindSequence++;
      const receipt = Object.freeze({
        artifactIdentity: options.artifactIdentity?.fingerprint || String(options.artifactIdentity || ''),
        resultIdentity: `test-result:${bindSequence}`,
        resultIndex: options.resultIndex,
        generationToken: String(options.generationToken),
        rootIdentity: bindSequence,
        rootGeneration: bindSequence,
        bindSequence,
        requiredBindingFlags: Object.freeze({ ready: true }),
        readyTimestamp: bindSequence,
        phase: options.phase
      });
      let rejectPromise;
      const promise = new Promise((resolve, reject) => {
        rejectPromise = reject;
        queueMicrotask(() => {
          if (activeExpectation?.receipt !== receipt) return;
          onEvent('test.preview-mounted');
          if (activeExpectation?.receipt !== receipt) return;
          onEvent('test.preview-bound');
          if (activeExpectation?.receipt !== receipt) return;
          activeReceipt = receipt;
          activeExpectation = null;
          onEvent('test.preview-ready');
          resolve(receipt);
        });
      });
      void promise.catch(() => {});
      activeExpectation = {
        generationToken: String(options.generationToken),
        receipt,
        reject: rejectPromise
      };
      return { promise };
    },
    invalidateReadinessExpectation(target, reason) {
      if (!activeExpectation || String(target) !== activeExpectation.generationToken) return false;
      const expectation = activeExpectation;
      activeExpectation = null;
      expectation.reject(reason || new Error('invalidated'));
      return true;
    },
    invalidateReadyReceipt(receipt) {
      if (activeReceipt === receipt) activeReceipt = null;
      return true;
    },
    getActiveRuntime() {
      return activeReceipt ? { readyReceipt: activeReceipt } : null;
    },
    async restorePreviousSelectedResult({ restore }) {
      onEvent('test.preview-restore-started');
      await restore();
      activeReceipt = Object.freeze({
        rootGeneration: nextBindSequence++,
        phase: 'rollback-restoration'
      });
      onEvent('test.preview-restore-ready');
      return activeReceipt;
    },
    selectResult(index) {
      if (state.selectedResultIndex.value === index) return false;
      state.selectedResultIndex.value = index;
      return true;
    }
  };
};
const immediatePreviewRuntime = createImmediatePreviewRuntime();
const generatedArtifactHandleOptions = {
  planGenomeScales: planWholeGenomeScales,
  serializeInputFile: serializeFile,
  captureGeneratedArtifactHandle: artifactSnapshots.captureGeneratedArtifactHandle,
  captureGeneratedArtifactOwnerSet: artifactSnapshots.captureGeneratedArtifactOwnerSet,
  installGeneratedArtifactOwnerSet: artifactSnapshots.installGeneratedArtifactOwnerSet,
  installRedrawnArtifact: artifactSnapshots.installRedrawnArtifact,
  displayedArtifactDrawing: artifactSnapshots.displayedArtifactDrawing,
  restoreGeneratedArtifactHandle: artifactSnapshots.restoreGeneratedArtifactHandle,
  setGeneratedArtifactIdentity: artifactSnapshots.setGeneratedArtifactIdentity,
  previewRuntime: immediatePreviewRuntime,
  // The label owner's ports (R13); these tests do not open the label dialog.
  closeLabelTextScopeDialog: () => {},
  clearLabelBuildNotices: () => {},
  moveDrawnLabelVisibility: () => {},
  nextTick: window.Vue.nextTick,
  waitForAfterPaint: async () => {}
};
const wireGeneratedArtifactRuntimeOwner = (runner) => {
  artifactSnapshots.setGeneratedArtifactRuntimeOwner({
    capture: runner.captureGeneratedArtifactRuntimeState,
    restore: runner.restoreGeneratedArtifactRuntimeState
  });
  return runner;
};

const result = (name, marker) => ({
  name,
  format: 'svg',
  content: `<svg data-marker="${marker}"></svg>`
});

const validCatalog = (name) => ({
  schema: 5,
  items: [{
    resultIndex: 0,
    resultName: name,
    recordKeys: ['record-1'],
    features: [{
      svgId: 'rendered-feature-1',
      recordKey: 'record-1',
      biologicalFeatureId: 'biological-feature-1',
      fillColor: '#abcdef',
      drawnSelector: { hash: 'rendered-feature-1', location: '0..9', recordLocation: 'record-1:0..9:+' }
    }],
    biologicalFeatures: [{
      recordKey: 'record-1',
      biologicalFeatureId: 'biological-feature-1',
      stableFeatureId: 'stable-feature-1',
      record_idx: 0,
      feature_index: 0,
      record_id: 'record-1',
      type: 'CDS',
      start: 0,
      end: 9,
      strand: 1,
      anchorProfile: {
        precision: 'exact', operator: 'single', partOrder: 'biological', strand: '+'
      },
      qualifiers: { product: ['audit protein'] }
    }],
    orthogroups: [],
    annotations: [],
    comparisonMatches: []
  }]
});

const response = (logicalResult, featureCatalog) => ({
  results: [logicalResult],
  metadata: featureCatalog === undefined ? {} : { featureCatalog, annotationWarnings: [{
    code: 'feature_selector_unmatched', setId: 's', annotationId: 'missing',
    recordId: 'record-1', recordIndex: 0, missingCount: 1,
    message: 'Skipped annotation: 1 feature selector(s) unmatched.',
    resultIndex: 0, resultName: logicalResult.name
  }] }
});

const storedZipEntries = async (blob) => {
  const bytes = new Uint8Array(await blob.arrayBuffer());
  const view = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
  const decoder = new TextDecoder();
  const entries = [];
  let offset = 0;
  while (offset + 30 <= bytes.byteLength && view.getUint32(offset, true) === 0x04034B50) {
    const size = view.getUint32(offset + 18, true);
    const nameLength = view.getUint16(offset + 26, true);
    const extraLength = view.getUint16(offset + 28, true);
    const nameStart = offset + 30;
    const dataStart = nameStart + nameLength + extraLength;
    entries.push({
      name: decoder.decode(bytes.subarray(nameStart, nameStart + nameLength)),
      text: decoder.decode(bytes.subarray(dataStart, dataStart + size)),
      bytes: bytes.slice(dataStart, dataStart + size)
    });
    offset = dataStart + size;
  }
  return entries;
};

const committedFeatureState = () => structuredClone({
  results: state.results.value,
  annotationWarnings: state.annotationWarnings.value,
  selectedResultIndex: state.selectedResultIndex.value,
  featureCatalog: state.featureCatalog.value,
  extractedFeatures: state.extractedFeatures.value,
  biologicalFeatures: state.biologicalFeatures.value,
  selectedFeatureIds: [...state.selectedFeatureIds.value],
  selectedFeatureAnchorId: state.selectedFeatureAnchorId.value,
  featureSelectionStatus: state.featureSelectionStatus.value,
  featureColorOverrides: state.activeDrawing().featureColorOverrides,
  featureStrokeOverrides: state.activeDrawing().featureStrokeOverrides,
  legendEntries: state.activeDrawing().legendEntries.value,
  featureEditorStatus: state.featureEditorStatus.value,
  featureExtractionPending: state.featureExtractionPending.value,
  featureRecordIds: state.featureRecordIds.value,
  selectedFeatureRecordIdx: state.selectedFeatureRecordIdx.value,
  orthogroups: state.orthogroups.value,
  resultGenerationKey: state.resultGenerationKey.value,
  lastRunInfo: state.lastRunInfo.value
});

const encodedResource = (kind, name, text) => {
  const bytes = new TextEncoder().encode(text);
  return {
    kind,
    name,
    type: 'text/plain',
    size: bytes.byteLength,
    lastModified: 0,
    encoding: 'base64',
    data: Buffer.from(bytes).toString('base64')
  };
};

const sha256Text = async (text) => {
  const digest = new Uint8Array(await crypto.subtle.digest(
    'SHA-256',
    new TextEncoder().encode(text)
  ));
  return [...digest].map((value) => value.toString(16).padStart(2, '0')).join('');
};

// Runs `generate` and returns its LOSAT timing line and how many SHA-256 digests of FASTA text it took, in
// all and by text.
const watchLosatPreparation = async (generate) => {
  const watched = { fastaDigests: 0, digestsByText: new Map(), timing: '' };
  const { subtle } = globalThis.crypto;
  const ownDigest = Object.hasOwn(subtle, 'digest');
  const digest = subtle.digest;
  const info = console.info;
  subtle.digest = function digestWatched(algorithm, data) {
    const text = new TextDecoder().decode(data);
    if (text.startsWith('>')) {
      watched.fastaDigests += 1;
      watched.digestsByText.set(text, (watched.digestsByText.get(text) || 0) + 1);
    }
    return digest.call(this, algorithm, data);
  };
  console.info = (...args) => {
    if (String(args[0]).startsWith('LOSAT timing')) watched.timing = String(args[0]);
    return info.apply(console, args);
  };
  try {
    await generate();
  } finally {
    if (ownDigest) subtle.digest = digest;
    else delete subtle.digest;
    console.info = info;
  }
  return watched;
};

// The pre-run check of a Generate that resolves an older Session's protein references through the current
// identity: it accepts or refuses as the rewrite it replaced did (OV-3352 follow-up; base answers recorded with
// that rewrite), and reads the values without changing them.
test('the legacy protein reference check answers as the rewrite it replaced', () => {
  const A = 'p_rec1_1_100_1_0123456789ab';
  const B = 'p_rec1_200_300_-1_0123456789ab';
  const duplicate = "Protein reference migration produced duplicate key 'h1'.";
  const remains = 'Legacy protein UI references remain after migration.';
  const cases = [
    ['every reference maps', { [A]: 'h1', [B]: 'h2' },
      [[{ id: 'og1', members: [A, B] }], A, [{ feature_id: 'f1', protein_id: B, note: `see ${A} here` }], [{ [A]: { label: B } }]], null],
    ['no reference', {}, [[], '', [{ feature_id: 'f1' }], null], null],
    ['a key maps onto a key the object has', { [A]: 'h1' }, [[], '', [{ [A]: 1, h1: 2 }], null], duplicate],
    ['a mapped handle reads as a legacy reference', { [A]: 'p_x_1_2_0_aaaaaaaaaaaa' }, [[A], '', [], null], remains],
    ['a key holds a reference inside other text', { [A]: 'h1' }, [[], '', [{ [`x${A}`]: 1 }], null], remains],
    ['a reference the map lacks', { [A]: 'h1' }, [[A, B], '', [], null], remains],
    ['the map is a list', [A], [[A], '', [], null], remains],
    ['a duplicate key deep in a list', { [A]: 'h1', [B]: 'h1' }, [[{ groups: [{ [A]: 1, [B]: 2 }] }], '', [], null], duplicate],
    ['a duplicate after a leftover reference', { [A]: 'h1' }, [[B], '', [{ [A]: 1, h1: 2 }], null], duplicate]
  ];
  for (const [name, map, values, refusal] of cases) {
    const before = structuredClone(values);
    const check = () => requireLegacyProteinReferencesMigrate(map, ...values);
    if (refusal) assert.throws(check, { message: refusal }, name);
    else assert.doesNotThrow(check, name);
    assert.deepEqual(values, before, `${name}: the check changes nothing`);
  }
});

test('afterPaint crosses one completed frame before resolving', async () => {
  const frames = [];
  let resolved = false;
  const painted = afterPaint({ requestFrame: (callback) => frames.push(callback) })
    .then(() => { resolved = true; });
  assert.equal(frames.length, 1);
  frames.shift()();
  assert.equal(resolved, false);
  assert.equal(frames.length, 1);
  frames.shift()();
  await painted;
  assert.equal(resolved, true);
});

// CAPACITY-WARN: the labels that may not fit stay under the status line until
// the run ends; the next Generate starts without them.
test('the processing overlay keeps each label notice Python sends while it renders', async () => {
  const { readFile } = await import('node:fs/promises');
  const overlay = { processingStatus: { value: '' }, labelCapacityNotices: { value: [] } };
  showDiagramProgress({ requestId: 1, stage: 'rendering' }, overlay);
  assert.deepEqual([overlay.processingStatus.value, overlay.labelCapacityNotices.value], ['Rendering diagram...', []]);
  const first = Object.freeze({ recordId: 'NC_1', labelKind: 'circular_outer', labelCount: 195, fontSize: 14, capacityRatio: 1.4 });
  const second = Object.freeze({ recordId: 'NC_2', labelKind: 'circular_inner', labelCount: 90, fontSize: 14, capacityRatio: 0.9 });
  showDiagramProgress({ requestId: 1, stage: 'rendering', notice: first }, overlay);
  showDiagramProgress({ requestId: 1, stage: 'rendering', notice: second }, overlay);
  showDiagramProgress({ requestId: 1, stage: 'finalizing' }, overlay);
  assert.deepEqual([overlay.processingStatus.value, overlay.labelCapacityNotices.value],
    ['Finalizing diagram results...', [first, second]]);
  // D-17: each pass of a figure-size fit sends its notice; a record's labels of one kind show once, the latest in place.
  const refitted = Object.freeze({ ...first, capacityRatio: 1.2 });
  showDiagramProgress({ requestId: 1, stage: 'rendering', notice: refitted }, overlay);
  assert.deepEqual(overlay.labelCapacityNotices.value, [refitted, second]);
  // Every run that shows the overlay (Generate, a target-record run, a Sequence search) starts it without the last
  // run's notices.
  const source = await readFile(new URL('../../gbdraw/web/js/app/run-analysis.js', import.meta.url), 'utf8');
  assert.equal(source.match(/processing\.value = true;\n\s*labelCapacityNotices\.value = \[\];/g)?.length, 3);
  assert.equal(source.match(/processing\.value = true;/g)?.length, 3);
  // The overlay shows the heading and each notice in the UI language (OV-2179), beside Cancel.
  const html = await readFile(new URL('../../gbdraw/web/index.html', import.meta.url), 'utf8');
  const start = html.indexOf('<div v-if="processing"');
  const markup = html.slice(start, html.indexOf('@click="cancelGeneration"', start));
  assert.match(markup, /<div v-if="labelCapacityNotices\.length"[^>]*data-generation-notice[^>]*>\s*<p[^>]*>\{\{ t\('Labels may not fit'\) \}\}<\/p>\s*<p v-for="\(notice, index\) in labelCapacityNotices" :key="index">\{\{ labelCapacityNoticeText\(notice\) \}\}<\/p>/);
});

test('afterFrame crosses one post-bind frame before resolving', async () => {
  const frames = [];
  let resolved = false;
  const framed = afterFrame({ requestFrame: (callback) => frames.push(callback) })
    .then(() => { resolved = true; });
  assert.equal(frames.length, 1);
  frames.shift()();
  await framed;
  assert.equal(resolved, true);
});

test('superseded canonical execution stops before catalog and SVG admission', async () => {
  let admissions = 0;
  let commits = 0;
  const outcome = await executeCanonicalRenderCandidate({
    canonical: { renderRequest: { schema: CANONICAL_REQUEST_SCHEMA }, resources: {} },
    mode: 'circular',
    kind: 'target-record-transform',
    generationExecutor: async () => ({ results: [], metadata: { annotationWarnings: ["stale-invalid-warning"] } }),
    shouldAdmit: () => false,
    catalogAdmission: () => { admissions += 1; },
    prepareCommit: () => { commits += 1; }
  });
  assert.equal(outcome.status, 'superseded');
  assert.equal(admissions, 0);
  assert.equal(commits, 0);
});

test('canonical execution carries the Worker comparison table ID warnings to the committed Result', async () => {
  const warning = {
    code: 'comparison_record_id_unmatched', queryRecordIndex: 0, subjectRecordIndex: 1,
    queryRecordId: 'R2', subjectRecordId: 'R3', rowCount: 1, exampleIds: ['contig_A'],
    message: 'Comparison between query record #1 \'R2\' and subject record #2 \'R3\': 1 row(s) use sequence IDs that match no displayed record (\'contig_A\'); these rows are drawn on the records assigned by position.',
    resultIndex: 0, resultName: 'out.svg'
  };
  const run = (metadata) => executeCanonicalRenderCandidate({
    canonical: { renderRequest: { schema: CANONICAL_REQUEST_SCHEMA }, resources: {} },
    mode: 'linear',
    generationExecutor: async () => ({ results: [{ name: 'out.svg', content: '<svg/>' }], metadata }),
    catalogAdmission: () => ({ catalog: null }),
    prepareCommit: () => ({})
  });
  const carried = await run({ comparisonWarnings: [warning] });
  assert.equal(carried.status, 'ok');
  assert.deepEqual(carried.comparisonWarnings, [warning]);
  assert.deepEqual((await run({})).comparisonWarnings, []);
  await assert.rejects(
    run({ comparisonWarnings: [{ ...warning, resultName: 'stale.svg' }] }),
    /Result metadata schema/
  );
});

test('a rerender admits the Results it draws again under their kept names (OV-136)', async () => {
  const admitted = [];
  const outcome = await executeCanonicalRenderCandidate({
    canonical: { renderRequest: { schema: CANONICAL_REQUEST_SCHEMA }, resources: {} },
    mode: 'circular',
    kind: 'reflow',
    resultNames: ['HmmtDNA_basic_circular'],
    generationExecutor: async () => ({
      results: [{ name: 'out.svg', content: '<svg/>' }],
      metadata: {
        featureCatalog: { items: [{ resultIndex: 0, resultName: 'out.svg' }] },
        annotationWarnings: [{
          code: 'empty_span', setId: 's', annotationId: 'a', recordId: 'NC_012920.1', recordIndex: 0,
          missingCount: 0, message: 'Empty span.', resultIndex: 0, resultName: 'out.svg'
        }]
      }
    }),
    catalogAdmission: (catalog, results) => {
      admitted.push({ items: catalog.items, names: results.map((result) => result.name) });
      return { catalog };
    },
    prepareCommit: ({ generationResponse }) => ({ results: generationResponse.results })
  });
  assert.equal(outcome.status, 'ok');
  assert.deepEqual(admitted, [{
    items: [{ resultIndex: 0, resultName: 'HmmtDNA_basic_circular' }],
    names: ['HmmtDNA_basic_circular']
  }]);
  assert.deepEqual(outcome.commit.results.map((result) => result.name), ['HmmtDNA_basic_circular']);
  assert.equal(outcome.annotationWarnings[0].resultName, 'HmmtDNA_basic_circular');
});

// One Circular GenBank input whose other slots hold files no run may read.
/** @param {AuditFile} primary @param {AuditInactiveFile} inactive */
const installSimpleCircularInput = (primary, inactive) => {
  installShownDrawing('circular');
  modeDrawing(state, 'circular').inputType.value = 'gb';
  setCircularFile(primary);
  modeDrawing(state, 'circular').files.c_depth = bound(inactive);
  modeDrawing(state, 'circular').files.c_conservation_blasts = bound([]);
  modeDrawing(state, 'circular').files.c_conservation_fastas = bound([]);
  modeDrawing(state, 'circular').files.c_conservation_sequence_sources = bound([]);
  modeDrawing(state, 'circular').files.d_color = inactive;
  modeDrawing(state, 'circular').files.t_color = inactive;
  modeDrawing(state, 'circular').files.blacklist = inactive;
  modeDrawing(state, 'circular').files.whitelist = inactive;
  modeDrawing(state, 'circular').files.qualifier_priority = inactive;
  state.activeDrawing().form.prefix = 'audit-simple';
  state.activeDrawing().form.show_depth = false;
  state.activeDrawing().form.multi_record_canvas = false;
  state.activeDrawing().adv.circular_track_slots_enabled = false;
  state.activeDrawing().circularConservation.enabled = false;
  state.activeDrawing().circularConservation.source = 'upload';
  state.activeDrawing().annotationSets.splice(0);
  modeDrawing(state, 'circular').circularRecordList.value = [listedRecord('audit')];
  Object.assign(modeDrawing(state, 'circular').circularRecordDiscovery, {
    status: 'ready',
    error: '',
    inputType: 'gb',
    sources: [{ uid: 'record-1', primaryFile: bound(primary), pairedFile: null, error: '' }]
  });
  state.drawingOf('linear').linearSeqs.splice(0, state.drawingOf('linear').linearSeqs.length, {
    uid: 'inactive-linear-row',
    gb: bound(inactive),
    gff: null,
    fasta: null,
    depth: bound(inactive),
    blast: null,
    losat_gencode: 1,
    losat_filename: '',
    definition: '',
    record_subtitle: '',
    region_record_id: '',
    region_start: null,
    region_end: null,
    region_reverse: false,
      windows: []
  });
};

test('audit-5 owner: direct simple createRunAnalysis path is worker-only and catalog-transactional', async () => {
  const structuralMetrics = {};
  const lifecycleEvents = [];
  globalThis.__GBDRAW_TEST_HOOKS__ = {
    onStructuralMetric(metric) {
      structuralMetrics[metric.name] = (
        Number(structuralMetrics[metric.name] || 0) + Number(metric.value || 0)
      );
    },
    onSessionLifecycleEvent({ name }) {
      lifecycleEvents.push(name);
    }
  };
  const primary = new AuditFile(['LOCUS       audit\nORIGIN\n//\n'], 'active.gb', {
    type: 'text/plain',
    lastModified: 7
  });
  const fallbackPrimary = new AuditFile(['record parser fallback'], 'fallback.gb', {
    type: 'application/octet-stream',
    lastModified: 8
  });
  // Discovery is shared per File (OV-1558): a scenario that needs its own
  // pending or failing read uses a File no earlier scenario read.
  const coalescedPrimary = new AuditFile(['record parser fallback'], 'coalesced.gb', {
    type: 'application/octet-stream', lastModified: 8
  });
  const stalePrimary = new AuditFile(['record parser fallback'], 'stale.gb', {
    type: 'application/octet-stream', lastModified: 8
  });
  const inactive = new AuditInactiveFile(['unused'], 'inactive.dat');

  installSimpleCircularInput(primary, inactive);

  let adoptedArtifacts = 0;
  let lastAdoptedCanonical = null;
  let failArtifactAdoption = false;
  let cancelDuringPreview = false;
  let failCandidateAdmission = false;
  let failDecorationCorrespondence = false;
  let failDecorationSnapshot = false;
  const activationOwnerSets = [];
  const captureForHistory = (...args) => {
    assert.equal(state.processing.value, true);
    assert.notEqual(state.processingStatus.value, '');
    lifecycleEvents.push('test.before-handle-capture');
    return artifactSnapshots.captureGeneratedArtifactHandle(args[0] ?? state.activeDrawing());
  };
  const generationHistory = createHistoryManager({
    buildIntent: async () => ({}),
    applyIntent: async () => {},
    buildCheckpoint: () => {
      throw new Error('Generate must not build a full artifact checkpoint.');
    },
    applyCheckpoint: async () => {
      throw new Error('Generate must not apply a full artifact checkpoint.');
    },
    captureGeneratedArtifactHandle: captureForHistory,
    restoreGeneratedArtifactHandle: artifactSnapshots.restoreGeneratedArtifactHandle,
  });
  await generationHistory.initializeIntentBaseline('Generated artifact baseline');
  const readinessRuntime = createImmediatePreviewRuntime({
    onEvent(name) {
      assert.equal(state.processing.value, true, name);
      if (name === 'test.preview-mounted' || name === 'test.preview-ready') {
        assert.equal(state.processingStatus.value, 'Preparing preview...', name);
      }
      lifecycleEvents.push(name);
      if (cancelDuringPreview && name === 'test.preview-mounted') {
        runner.cancelRunAnalysis();
      }
    }
  });
  // OV-36: Generate's start clears only the label build warning; a Generate
  // that commits its Results also clears the last live-edit failure.
  const labelNoticeClears = [];
  let runner;
  runner = wireGeneratedArtifactRuntimeOwner(createRunAnalysis({
    ...generatedArtifactHandleOptions,
    clearLabelBuildNotices: ({ rerender = false } = {}) => {
      labelNoticeClears.push(rerender ? 'rerender' : 'build');
    },
    captureGeneratedArtifactHandle: captureForHistory,
    installGeneratedArtifactOwnerSet: (ownerSet, options) => {
      activationOwnerSets.push(ownerSet);
      artifactSnapshots.installGeneratedArtifactOwnerSet(ownerSet, options);
      assert.equal(state.results.value, ownerSet.results);
      assert.equal(state.featureCatalog.value, ownerSet.featureCatalog);
      assert.equal(state.losatCache.value, ownerSet.losatCache);
      assert.equal(state.matchSequenceRegistry.captureTrustedOwner(), ownerSet.matchSequenceOwner);
    },
    waitForAfterPaint: async () => {
      lifecycleEvents.push('test.after-paint-returned');
    },
    waitForPostBindFrame: async () => {
      assert.equal(state.processing.value, true);
      lifecycleEvents.push('test.post-bind-frame-returned');
    },
    previewRuntime: readinessRuntime,
    runGeneratedArtifactReplacement: (...args) => (
      generationHistory.runUndoableArtifactReplacement(...args)
    ),
    state,
    captureDecorationContinuity: () => {
      if (failDecorationSnapshot) throw Object.assign(new Error('Invalid offset'), { code: 'DECORATION_CONTINUITY' });
      const mismatch = failDecorationCorrespondence;
      return mismatch ? () => {
        throw Object.assign(new Error('Cannot preserve title placement. Reset Layout.'), {
          code: 'DECORATION_CONTINUITY', stage: 'result-admission',
          context: { field: 'title', inputOrdinal: 1, reason: 'DECORATION_TARGET' }
        });
      } : null;
    },
    writeSessionDocumentFile,
    serializeCanonicalFiles: (_snapshot, _catalog, drawing) => {
      assert.equal(state.processingStatus.value, 'Preparing render inputs and session...');
      return serializeActiveRenderFiles(state, drawing);
    },
    adoptCanonicalRenderArtifacts: (canonical) => {
      if (failArtifactAdoption) {
        throw new Error('forced late canonical artifact adoption failure');
      }
      adoptedArtifacts += 1;
      // The committed Session shares its file descriptors (a File Generate described is file-backed).
      lastAdoptedCanonical = cloneSharingResources(canonical);
    },
    prepareCandidateCommit: ({
      generationResponse,
      catalogAdmission,
      featureColorOverrides
    }) => {
      if (failCandidateAdmission) {
        throw new Error('forced current Result admission failure');
      }
      const candidate = {
        results: structuredClone(generationResponse.results),
        featureState: featureStateFromCatalog(catalogAdmission.catalog),
        featureColorOverrides: structuredClone(featureColorOverrides)
      };
      return candidate;
    },
    resetPreviewViewport: ({ pan = null } = {}) => {
      state.canvasPan.x = Number(pan?.x) || 0;
      state.canvasPan.y = Number(pan?.y) || 0;
      if (!pan) state.zoom.value = 1;
    }
  }));

  const committedResult = result('audit-simple.svg', 'committed');
  const committedCatalog = validCatalog(committedResult.name);
  workerResponses.push(response(committedResult, committedCatalog));
  // OV-84: a stroke on a feature this Result does not draw (here, a record of
  // the other mode) stays in the draft through the Generate.
  const dormantStrokeKey = 'linear-seq-other\u0000feature-x';
  state.activeDrawing().featureStrokeOverrides[dormantStrokeKey] = { strokeColor: '#ff0000', strokeWidth: 3 };

  assert.deepEqual(
    await runner.runAnalysis(),
    { status: 'ok' },
    JSON.stringify(state.errorLog.value)
  );
  assert.deepEqual(state.activeDrawing().featureStrokeOverrides, {
    [dormantStrokeKey]: { strokeColor: '#ff0000', strokeWidth: 3 }
  });
  delete state.activeDrawing().featureStrokeOverrides[dormantStrokeKey];
  assert.deepEqual(labelNoticeClears.splice(0), ['build', 'rerender']);
  assert.equal(state.failedGeneratePreservedResult.value, false);
  assert.deepEqual(state.results.value, [committedResult]);
  assert.deepEqual(state.featureCatalog.value, committedCatalog);
  assert.equal(state.extractedFeatures.value.length, 1);
  assert.equal(state.biologicalFeatures.value.length, 1);
  assert.equal(state.lastRunInfo.value.sourceRecipe.available, true);
  assert.match(state.lastRunInfo.value.sourceRecipe.command, /--gbk active\.gb/);
  assert.doesNotMatch(state.lastRunInfo.value.sourceRecipe.command, /--session/);
  assert.match(state.lastRunInfo.value.exactReplay.command, /--session audit-simple\.gbdraw$/);
  assert.ok(
    lifecycleEvents.indexOf('generate.processing-published')
      < lifecycleEvents.indexOf('test.after-paint-returned')
  );
  assert.ok(
    lifecycleEvents.indexOf('generate.paint-opportunity-completed')
      < lifecycleEvents.indexOf('history.before-capture-started')
  );
  assert.ok(
    lifecycleEvents.indexOf('test.before-handle-capture')
      < lifecycleEvents.indexOf('history.before-capture-completed')
  );
  const orderedSuccessEvents = [
    'generate.processing-published',
    'generate.paint-opportunity-completed',
    'history.before-capture-completed',
    'result-admission-end',
    'artifact.activation-completed',
    'test.preview-mounted',
    'test.preview-bound',
    'test.preview-ready',
    'test.post-bind-frame-returned',
    'preview.post-bind-frame-completed',
    'history.finalization-completed',
    'artifact.finalization-completed',
    'generate.completed',
    'generate.processing-cleared'
  ];
  orderedSuccessEvents.slice(1).forEach((eventName, index) => {
    assert.ok(
      lifecycleEvents.indexOf(orderedSuccessEvents[index])
        < lifecycleEvents.indexOf(eventName),
      `${orderedSuccessEvents[index]} before ${eventName}`
    );
  });
  assert.equal(activationOwnerSets.length, 1);
  assert.equal(structuralMetrics.generatedArtifactCandidateBuildCount, 1);
  assert.equal(structuralMetrics.generatedArtifactActivationCount, 1);
  assert.equal(structuralMetrics.generatedArtifactFinalizeCount, 1);
  assert.equal(structuralMetrics.historyReplacementCount, 1);
  assert.equal(structuralMetrics.canonicalReplayFullSerializationCount || 0, 0);
  const readsBeforeReplay = activePrimaryReads;
  const streamsBeforeReplay = activePrimaryStreams;
  await runner.downloadCliHelperFiles();
  assert.equal(structuralMetrics.canonicalReplayFullSerializationCount, 1);
  // The replay Session streams the input into its entry once and never reads it whole.
  assert.deepEqual({ reads: activePrimaryReads - readsBeforeReplay, streams: activePrimaryStreams - streamsBeforeReplay },
    { reads: 0, streams: 1 });
  assert.equal(downloadedFilename, 'audit-simple-cli-files.zip');
  const helperEntries = await storedZipEntries(downloadedBlob);
  assertDownloadedFilesMatchCommands(state.lastRunInfo.value, helperEntries);
  // The exact replay is a Session file container, written by Save's writer (OV-2695).
  const replayEntry = helperEntries.find(({ name }) => name.endsWith('.gbdraw'));
  assert.ok(replayEntry);
  const replay = /** @type {Record<string, any>} */ ((await readSessionContainer(new Blob([replayEntry.bytes]))).document);
  assert.equal(replay.format, 'gbdraw-session');
  assert.equal(replay.version, SESSION_VERSION);
  // A replay Session over the Session limit is refused where the app shows errors, not lost.
  globalThis.__GBDRAW_TEST_HOOKS__.maxExpandedSessionBytes = 1;
  downloadedFilename = null;
  try {
    await runner.downloadCliHelperFiles();
  } finally {
    delete globalThis.__GBDRAW_TEST_HOOKS__.maxExpandedSessionBytes;
  }
  assert.equal(state.errorLog.value?.code, 'SESSION_SIZE_LIMIT');
  assert.equal(downloadedFilename, null);
  state.errorLog.value = null;
  // The exact replay is a Session 49 of the drawn drawing alone, by its ID and name (OV-1806).
  const drawn = state.activeDrawing();
  const drawnTab = state.drawingTabs.value.find((tab) => tab.id === drawn.id);
  assert.deepEqual(replay.drawings.map(({ id, name }) => ({ id, name })), [{ id: drawn.id, name: drawnTab.name }]);
  assert.equal(replay.ui.activeDrawingId, drawn.id);
  const replayView = savedDrawingView(replay, drawn.id);
  assert.equal(replayView.renderRequest.schema, CANONICAL_REQUEST_SCHEMA);
  assert.equal(Object.hasOwn(replay, 'renderRequest'), false);
  assert.ok(Object.keys(replay.resources).length > 0);

  const firstRunPayload = workerMessages.find(({ type }) => type === 'run').payload;
  // N7 (step 3a): the first Generate holds no catalog; later renders of the
  // shown drawing name the catalog revision it holds.
  assert.equal(firstRunPayload.heldFeatureCatalogRevision, null);
  const committedCatalogRevision = `fake:${fakeCatalogBuilds}`;
  assert.equal(Object.hasOwn(firstRunPayload, 'resources'), false);
  assert.ok(firstRunPayload.resourceManifest.length > 0);
  assert.ok(firstRunPayload.stagedResources.every(({ blob }) => blob instanceof Blob));

  let committedState = committedFeatureState();
  let committedExtractedFeatureIdentity = state.extractedFeatures.value;
  let committedBiologicalFeatureIdentity = state.biologicalFeatures.value;

  const targetCandidate = cloneSharingResources(lastAdoptedCanonical);
  targetCandidate.renderRequest.records[0].display.startCoordinate = 3;
  targetCandidate.renderRequest.records[0].presentation.reverseComplement = true;
  const committedOutputPrefix = targetCandidate.renderRequest.output.prefix;
  state.activeDrawing().form.prefix = 'unrelated-pending-prefix';
  let targetIntent = {
    startCoordinate: null,
    reverseComplementOverride: null,
    anchorIntent: null
  };
  const targetResult = result('target-only.svg', 'target-only');
  workerResponses.push(response(targetResult, validCatalog(targetResult.name)));
  const historyCountBeforeTarget = generationHistory.getUndoCount();
  let targetLosatExecutorJobs = 0;
  globalThis.__GBDRAW_LOSAT_EXECUTOR__ = async (jobs) => {
    targetLosatExecutorJobs += jobs.length;
    return [];
  };
  state.errorLog.value = { summary: 'Previous candidate failed; retry remains available.' };
  state.failedGeneratePreservedResult.value = true;
  assert.deepEqual(await runner.runCommittedCanonicalCandidate({
    canonical: targetCandidate,
    captureIntentCheckpoint: () => structuredClone(targetIntent),
    restoreIntentCheckpoint: (checkpoint) => { targetIntent = structuredClone(checkpoint); },
    commitIntent: () => {
      targetIntent = {
        startCoordinate: 3,
        reverseComplementOverride: true,
        anchorIntent: { schema: 1, recordKey: 'record-1' }
      };
    }
  }), { status: 'ok' });
  delete globalThis.__GBDRAW_LOSAT_EXECUTOR__;
  assert.equal(targetLosatExecutorJobs, 0);
  assert.equal(state.errorLog.value, null);
  assert.equal(state.failedGeneratePreservedResult.value, false);
  assert.equal(generationHistory.getUndoCount(), historyCountBeforeTarget + 1);
  assert.equal(state.activeDrawing().form.prefix, 'unrelated-pending-prefix');
  assert.equal(targetIntent.startCoordinate, 3);
  assert.deepEqual(state.results.value, [targetResult]);
  const targetRunPayload = workerMessages.filter(({ type }) => type === 'run').at(-1).payload;
  assert.equal(targetRunPayload.heldFeatureCatalogRevision?.revision, committedCatalogRevision);
  const targetCatalogRevision = `fake:${fakeCatalogBuilds}`;
  assert.equal(targetRunPayload.request.output.prefix, committedOutputPrefix);
  assert.notEqual(targetRunPayload.request.output.prefix, state.activeDrawing().form.prefix);
  assert.equal(structuralMetrics.canonicalCandidateExecutionCount, 2);
  await generationHistory.undo();
  assert.deepEqual(state.results.value, [committedResult]);
  assert.equal(targetIntent.startCoordinate, null);
  await generationHistory.redo();
  assert.deepEqual(state.results.value, [targetResult]);
  assert.equal(targetIntent.startCoordinate, 3);
  assert.equal(state.errorLog.value, null);
  assert.equal(state.failedGeneratePreservedResult.value, false);
  committedState = committedFeatureState();
  committedExtractedFeatureIdentity = state.extractedFeatures.value;
  committedBiologicalFeatureIdentity = state.biologicalFeatures.value;

  const targetHistoryCount = generationHistory.getUndoCount();
  const targetStateBeforeFailure = committedFeatureState();
  workerResponses.push(response(result('target-rejected.svg', 'rejected'), undefined));
  const rejectedTarget = await runner.runCommittedCanonicalCandidate({
    canonical: targetCandidate,
    captureIntentCheckpoint: () => structuredClone(targetIntent),
    restoreIntentCheckpoint: (checkpoint) => { targetIntent = structuredClone(checkpoint); },
    commitIntent: () => { throw new Error('must not commit rejected target intent'); }
  });
  assert.equal(rejectedTarget.status, 'error');
  assert.equal(rejectedTarget.error.code, 'FEATURE_METADATA');
  assert.equal(rejectedTarget.error.stage, 'result-admission');
  assert.equal(rejectedTarget.recovery, 'preserved');
  assert.deepEqual(committedFeatureState(), targetStateBeforeFailure);
  assert.equal(targetIntent.startCoordinate, 3);
  assert.equal(generationHistory.getUndoCount(), targetHistoryCount);

  const canceledTargetResult = result('target-canceled.svg', 'canceled');
  workerResponses.push(response(canceledTargetResult, validCatalog(canceledTargetResult.name)));
  cancelDuringPreview = true;
  assert.deepEqual(await runner.runCommittedCanonicalCandidate({
    canonical: targetCandidate,
    captureIntentCheckpoint: () => structuredClone(targetIntent),
    restoreIntentCheckpoint: (checkpoint) => { targetIntent = structuredClone(checkpoint); },
    commitIntent: () => { throw new Error('must not commit canceled target intent'); }
  }), { status: 'canceled' });
  cancelDuringPreview = false;
  assert.deepEqual(committedFeatureState(), targetStateBeforeFailure);
  assert.equal(targetIntent.startCoordinate, 3);
  assert.equal(generationHistory.getUndoCount(), targetHistoryCount);

  state.activeDrawing().form.prefix = 'audit-simple';

  workerResponses.push(response(result('missing.svg', 'missing'), undefined));
  const metricsBeforePreActivationFailure = { ...structuralMetrics };
  labelNoticeClears.length = 0;
  const failedRun = await runner.runAnalysis();
  assert.equal(
    workerMessages.filter(({ type }) => type === 'run').at(-1).payload.heldFeatureCatalogRevision?.revision,
    targetCatalogRevision
  );
  assert.equal(failedRun.status, 'error');
  assert.deepEqual(labelNoticeClears, ['build']);
  assert.equal(failedRun.error, state.errorLog.value);
  assert.equal(state.failedGeneratePreservedResult.value, true);
  assert.match(
    state.errorLog.value?.summary || '',
    /incompatible feature metadata/
  );
  assert.deepEqual(committedFeatureState(), committedState);
  assert.equal(
    structuralMetrics.generatedArtifactCandidateBuildCount,
    metricsBeforePreActivationFailure.generatedArtifactCandidateBuildCount
  );
  assert.equal(
    structuralMetrics.generatedArtifactActivationCount,
    metricsBeforePreActivationFailure.generatedArtifactActivationCount
  );
  assert.equal(
    structuralMetrics.generatedArtifactFinalizeCount,
    metricsBeforePreActivationFailure.generatedArtifactFinalizeCount
  );
  assert.equal(
    structuralMetrics.generatedArtifactRollbackCount || 0,
    metricsBeforePreActivationFailure.generatedArtifactRollbackCount || 0
  );

  const decorationFailureState = committedFeatureState();
  const decorationRequest = lastAdoptedCanonical;
  const decorationHistoryCount = generationHistory.getUndoCount();
  const decorationActivationCount = structuralMetrics.generatedArtifactActivationCount;
  failDecorationCorrespondence = true;
  const unmatchedResult = result('unmatched.svg', 'unmatched');
  workerResponses.push(response(unmatchedResult, validCatalog(unmatchedResult.name)));
  const unmatched = await runner.runAnalysis();
  failDecorationCorrespondence = false;
  assert.equal(unmatched.status, 'error');
  assert.equal(unmatched.error.code, 'DECORATION_CONTINUITY');
  assert.match(unmatched.error.summary, /Reset Layout.*Result 1.*title.*missing/);
  assert.deepEqual(committedFeatureState(), decorationFailureState);
  assert.equal(lastAdoptedCanonical, decorationRequest);
  assert.equal(generationHistory.getUndoCount(), decorationHistoryCount);
  assert.equal(structuralMetrics.generatedArtifactActivationCount, decorationActivationCount);

  failDecorationSnapshot = true;
  const snapshotFailure = await runner.runAnalysis();
  failDecorationSnapshot = false;
  assert.equal(snapshotFailure.status, 'error');
  assert.equal(snapshotFailure.recovery, 'preserved');
  assert.equal(state.failedGeneratePreservedResult.value, true);
  assert.deepEqual(committedFeatureState(), decorationFailureState);
  assert.equal(generationHistory.getUndoCount(), decorationHistoryCount);

  const admissionFailureState = committedFeatureState();
  const metricsBeforeAdmissionFailure = { ...structuralMetrics };
  failCandidateAdmission = true;
  const rejectedAdmissionResult = result('rejected-admission.svg', 'rejected');
  workerResponses.push(response(
    rejectedAdmissionResult,
    validCatalog(rejectedAdmissionResult.name)
  ));
  assert.equal((await runner.runAnalysis()).status, 'error');
  failCandidateAdmission = false;
  assert.equal(state.errorLog.value?.code, 'UNKNOWN');
  assert.doesNotMatch(JSON.stringify(state.errorLog.value), /forced current Result admission failure/);
  assert.deepEqual(committedFeatureState(), admissionFailureState);
  for (const metric of [
    'generatedArtifactCandidateBuildCount',
    'generatedArtifactActivationCount',
    'generatedArtifactRollbackCount',
    'generatedArtifactFinalizeCount',
    'historyReplacementCount'
  ]) {
    assert.equal(
      structuralMetrics[metric] || 0,
      metricsBeforeAdmissionFailure[metric] || 0,
      metric
    );
  }

  workerResponses.push(response(
    result('malformed.svg', 'malformed'),
    { schema: 2, items: [] }
  ));
  assert.equal((await runner.runAnalysis()).status, 'error');
  assert.equal(state.failedGeneratePreservedResult.value, true);
  assert.match(
    state.errorLog.value?.summary || '',
    /incompatible feature metadata/
  );
  assert.deepEqual(committedFeatureState(), committedState);

  state.zoom.value = 1.7;
  state.canvasPan.x = 31;
  state.canvasPan.y = -12;
  const bindFailureState = committedFeatureState();
  const bindHistoryCount = generationHistory.getUndoCount();
  const registerReadiness = readinessRuntime.registerReadinessExpectation;
  readinessRuntime.registerReadinessExpectation = options => {
    const expectation = registerReadiness(options);
    const promise = Promise.reject(new Error('forced candidate bind failure'));
    void promise.catch(() => {});
    return { ...expectation, promise };
  };
  const bindFailureResult = result('bind-failure.svg', 'bind-failure');
  workerResponses.push(response(bindFailureResult, validCatalog(bindFailureResult.name)));
  const bindFailure = await runner.runAnalysis();
  readinessRuntime.registerReadinessExpectation = registerReadiness;
  assert.equal(bindFailure.status, 'error');
  assert.equal(bindFailure.recovery, 'restored');
  assert.deepEqual(committedFeatureState(), bindFailureState);
  assert.equal(generationHistory.getUndoCount(), bindHistoryCount);

  const lateFailureState = committedFeatureState();
  const metricsBeforeLateFailure = { ...structuralMetrics };
  failArtifactAdoption = true;
  const lateFailureResult = result('late-failure.svg', 'late-failure');
  workerResponses.push(response(
    lateFailureResult,
    validCatalog(lateFailureResult.name)
  ));
  const rolledBack = await runner.runAnalysis();
  assert.equal(rolledBack.status, 'error');
  assert.equal(rolledBack.recovery, 'restored');
  failArtifactAdoption = false;
  assert.equal(state.errorLog.value?.code, 'UNKNOWN');
  assert.doesNotMatch(JSON.stringify(state.errorLog.value), /forced late canonical artifact adoption failure/);
  assert.deepEqual(committedFeatureState(), lateFailureState);
  assert.equal(state.extractedFeatures.value, committedExtractedFeatureIdentity);
  assert.equal(state.biologicalFeatures.value, committedBiologicalFeatureIdentity);
  assert.equal(state.zoom.value, 1.7);
  assert.deepEqual(state.canvasPan, { x: 31, y: -12 });
  assert.equal(
    structuralMetrics.generatedArtifactActivationCount,
    metricsBeforeLateFailure.generatedArtifactActivationCount + 1
  );
  assert.equal(
    structuralMetrics.generatedArtifactRollbackCount,
    (metricsBeforeLateFailure.generatedArtifactRollbackCount || 0) + 1
  );
  assert.equal(
    structuralMetrics.generatedArtifactFinalizeCount,
    metricsBeforeLateFailure.generatedArtifactFinalizeCount
  );
  const lateFailureEvents = lifecycleEvents.slice(
    lifecycleEvents.lastIndexOf('generate.processing-published')
  );
  assert.ok(
    lateFailureEvents.indexOf('artifact.rollback-started')
      < lateFailureEvents.indexOf('test.preview-restore-ready')
  );
  assert.equal(lateFailureEvents.includes('artifact.finalization-completed'), false);

  // A failed readiness restoration must retain the initiating cause, and must
  // never advertise completed rollback merely because the old Results exist.
  const originalRestore = readinessRuntime.restorePreviousSelectedResult;
  const checkpointBeforeRestoreFailure = await artifactSnapshots.captureGeneratedArtifactHandle(state.activeDrawing());
  readinessRuntime.restorePreviousSelectedResult = async () => { throw new Error('PRIVATE_ROLLBACK_SENTINEL'); };
  failArtifactAdoption = true;
  const restoreFailureResult = result('restore-failure.svg', 'restore-failure');
  workerResponses.push(response(restoreFailureResult, validCatalog(restoreFailureResult.name)));
  const failedRestore = await runner.runAnalysis();
  assert.equal(failedRestore.status, 'error');
  assert.equal(failedRestore.recovery, 'restore-failed');
  assert.equal(failedRestore.error.code, 'UNKNOWN');
  assert.equal(state.failedGeneratePreservedResult.value, false);
  assert.doesNotMatch(JSON.stringify(failedRestore.error), /PRIVATE_ROLLBACK_SENTINEL/);
  readinessRuntime.restorePreviousSelectedResult = originalRestore;
  failArtifactAdoption = false;
  await artifactSnapshots.restoreGeneratedArtifactHandle(checkpointBeforeRestoreFailure);

  const canceledState = committedFeatureState();
  const canceledResultIdentity = state.results.value;
  const metricsBeforePreviewCancel = { ...structuralMetrics };
  const canceledResult = result('canceled.svg', 'canceled');
  workerResponses.push(response(canceledResult, validCatalog(canceledResult.name)));
  cancelDuringPreview = true;
  assert.deepEqual(await runner.runAnalysis(), { status: 'canceled' });
  cancelDuringPreview = false;
  assert.equal(state.failedGeneratePreservedResult.value, true);
  assert.deepEqual(committedFeatureState(), canceledState);
  assert.equal(state.results.value, canceledResultIdentity);
  assert.deepEqual(state.results.value, canceledResultIdentity);
  assert.equal(state.results.value[0], canceledResultIdentity[0]);
  assert.equal(state.extractedFeatures.value, committedExtractedFeatureIdentity);
  assert.equal(state.biologicalFeatures.value, committedBiologicalFeatureIdentity);
  assert.equal(state.processingStatus.value, 'Canceled.');
  assert.equal(
    structuralMetrics.generatedArtifactActivationCount,
    metricsBeforePreviewCancel.generatedArtifactActivationCount + 1
  );
  assert.equal(
    structuralMetrics.generatedArtifactRollbackCount,
    metricsBeforePreviewCancel.generatedArtifactRollbackCount + 1
  );
  assert.equal(
    structuralMetrics.generatedArtifactFinalizeCount,
    metricsBeforePreviewCancel.generatedArtifactFinalizeCount
  );

  // GE-09: a Cancel with no Worker request keeps the warm Worker and its staged
  // primary source, so later runs do not read the file again. The reproducibility
  // download's replay Session streams the input and reads none of it (OV-2612).
  assert.equal(activePrimaryReads, 1);
  assert.equal(inactiveFileReads, 0);
  assert.equal(adoptedArtifacts, 2);
  assert.equal(workerMessages.filter(({ type }) => type === 'run').length, 12); // Includes the explicit failed-restoration attempt.
  assert.equal(workerMessages.filter(({ type }) => type === 'feature-extraction').length, 0);

  state.activeDrawing().form.multi_record_canvas = true;
  modeDrawing(state, 'circular').circularRecordList.value = [listedRecord('audit')];
  Object.assign(modeDrawing(state, 'circular').circularRecordDiscovery, {
    status: 'ready',
    error: '',
    inputType: 'gb',
    sources: [{ uid: 'record-1', primaryFile: bound(primary), pairedFile: null, error: '' }]
  });
  const readyResult = result('ready.svg', 'ready');
  workerResponses.push(response(readyResult, validCatalog(readyResult.name)));
  assert.deepEqual(await runner.runAnalysis(), { status: 'ok' });
  assert.equal(state.failedGeneratePreservedResult.value, false);
  assert.equal(activePrimaryReads, 1);

  // Cancel requested while preparation is pending ends before rendering.
  const runCount = () => workerMessages.filter(({ type }) => type === 'run').length;
  const runsBeforePreparedCancel = runCount();
  let releasePreparation;
  const canceledDuringPreparation = runner.runAnalysis(null, null, null, {
    prepareGenerate: () => new Promise((resolve) => { releasePreparation = resolve; })
  });
  while (!releasePreparation) await new Promise((resolve) => setImmediate(resolve));
  assert.equal(state.processing.value, true);
  runner.cancelRunAnalysis();
  releasePreparation({ status: 'ready', comparisonPlanSnapshot: null, comparisonExecution: null });
  assert.deepEqual(await canceledDuringPreparation, { status: 'canceled' });
  assert.equal(runCount(), runsBeforePreparedCancel);
  assert.equal(state.processing.value, false);
  assert.equal(state.failedGeneratePreservedResult.value, true);

  // A preparation that settles after a newer Generate is discarded as stale.
  let releaseStalePreparation;
  const stalePreparation = runner.runAnalysis(null, null, null, {
    prepareGenerate: () => new Promise((resolve) => { releaseStalePreparation = resolve; })
  });
  while (!releaseStalePreparation) await new Promise((resolve) => setImmediate(resolve));
  const newerResult = result('newer.svg', 'newer');
  workerResponses.push(response(newerResult, validCatalog(newerResult.name)));
  assert.deepEqual(await runner.runAnalysis(), { status: 'ok' });
  const runsAfterNewer = runCount();
  releaseStalePreparation({ status: 'ready', comparisonPlanSnapshot: null, comparisonExecution: null });
  assert.deepEqual(await stalePreparation, { status: 'stale' });
  assert.equal(runCount(), runsAfterNewer);
  assert.deepEqual(state.results.value, [newerResult]);
  assert.equal(state.processing.value, false);

  // A later operation owns the notification while this Worker response is held.
  const { normalizeUserFacingError } = await import('../../gbdraw/web/js/utils/error-normalization.js');
  let releaseOlderFailure;
  workerResponses.push(new Promise(resolve => { releaseOlderFailure = resolve; }));
  const runsBeforeOlderFailure = workerMessages.filter(message => message.type === 'run').length;
  const olderFailure = runner.runAnalysis();
  while (workerMessages.filter(message => message.type === 'run').length === runsBeforeOlderFailure) {
    await new Promise(resolve => setTimeout(resolve, 0));
  }
  const laterAlert = normalizeUserFacingError({ code: 'PDF_GLYPH', operation: 'export-pdf',
    stage: 'font-validation', context: { codepoint: 945 } });
  state.errorLog.value = laterAlert;
  releaseOlderFailure({ error: { code: 'REGEX_SYNTAX', operation: 'generate', stage: 'rule-validation',
    context: { reason: 'UNTERMINATED_SET' }, message: 'PRIVATE_OLDER_FAILURE' } });
  assert.deepEqual(await olderFailure, { status: 'stale' });
  assert.equal(state.errorLog.value, laterAlert);

  Object.assign(modeDrawing(state, 'circular').circularRecordDiscovery, {
    status: 'idle',
    error: '',
    inputType: '',
    sources: []
  });
  setCircularFile(fallbackPrimary);
  const workerRunCountBeforeDiscoveryFailure = workerMessages
    .filter(({ type }) => type === 'run')
    .length;
  const workerHelperCountBeforeDiscoveryFailure = workerMessages
    .filter(({ type }) => type === 'helper')
    .length;
  workerHelperResponses.push({
    ok: false,
    error: { code: 'INPUT_UNREADABLE', operation: 'listSequenceRecords', stage: 'helper', message: 'PRIVATE_DISCOVERY_SENTINEL' }
  });
  assert.equal((await runner.runAnalysis()).status, 'error');
  assert.equal(
    workerMessages.filter(({ type }) => type === 'helper').length,
    workerHelperCountBeforeDiscoveryFailure + 1
  );
  assert.equal(
    workerMessages.filter(({ type }) => type === 'run').length,
    workerRunCountBeforeDiscoveryFailure
  );
  assert.match(
    state.errorLog.value?.summary || '',
    /input could not be read/i
  );
  assert.doesNotMatch(
    JSON.stringify(state.errorLog.value),
    /Traceback|simple render must not initialize Pyodide/
  );

  state.activeDrawing().form.multi_record_canvas = false;
  state.activeDrawing().form.show_depth = false;
  modeDrawing(state, 'circular').files.c_depth = bound([[inactive]]);
  state.activeDrawing().adv.circular_track_slots_enabled = true;
  state.activeDrawing().adv.circular_track_slots = [{
    id: 'custom-depth',
    renderer: 'depth',
    enabled: true,
    side: 'inside',
    params: { track_index: 0 }
  }];
  Object.assign(modeDrawing(state, 'circular').circularRecordDiscovery, {
    status: 'idle',
    error: '',
    inputType: '',
    sources: []
  });
  workerHelperResponses.push({
    ok: false,
    error: { name: 'Error', message: 'injected depth record discovery failure' }
  });
  assert.equal((await runner.runAnalysis()).status, 'error');
  assert.equal(
    workerMessages.filter(({ type }) => type === 'run').length,
    workerRunCountBeforeDiscoveryFailure
  );

  let releaseCoalescedDiscovery;
  workerHelperResponses.push(new Promise((resolve) => {
    releaseCoalescedDiscovery = () => resolve({
      ok: true,
      result: {
        records: [{ selector: '#1', record_id: 'audit', record_length: 10 }]
      }
    });
  }));
  modeDrawing(state, 'circular').circularRecordList.value = [];
  setCircularFile(coalescedPrimary);
  Object.assign(modeDrawing(state, 'circular').circularRecordDiscovery, {
    status: 'loading',
    error: '',
    inputType: 'gb',
    sources: [{ uid: 'record-1', primaryFile: bound(coalescedPrimary), pairedFile: null, error: '' }]
  });
  const firstCoalescedDiscovery = runner.refreshCircularRecordOrder(undefined, state.activeDrawing());
  const secondCoalescedDiscovery = runner.refreshCircularRecordOrder(undefined, state.activeDrawing());
  assert.equal(firstCoalescedDiscovery, secondCoalescedDiscovery);
  await Promise.resolve();
  releaseCoalescedDiscovery();
  await Promise.all([firstCoalescedDiscovery, secondCoalescedDiscovery]);
  // The record's identity is its card's uid (a one-record file).
  assert.equal(modeDrawing(state, 'circular').circularRecordList.value[0].identity, 'record-1');

  // IN-05 (OIPC-C06): an inactive Circular source keeps its records and the
  // Multi-Record Canvas order through a mode round trip without rereading them.
  state.activeDrawing().adv.multi_record_positions.splice(
    0,
    state.activeDrawing().adv.multi_record_positions.length,
    { selector: '#1', row: 3 }
  );
  const circularStateBeforeRoundTrip = structuredClone({
    records: modeDrawing(state, 'circular').circularRecordList.value,
    positions: state.activeDrawing().adv.multi_record_positions,
    discovery: { ...modeDrawing(state, 'circular').circularRecordDiscovery }
  });
  const helperRequestsBeforeRoundTrip = workerMessages.filter(({ type }) => type === 'helper').length;
  installShownDrawing('linear');
  await runner.refreshCircularRecordOrder(undefined, state.activeDrawing());
  await runner.refreshCircularRecordOrder({ automatic: true }, state.activeDrawing());
  installShownDrawing('circular');
  await runner.refreshCircularRecordOrder({ automatic: true }, state.activeDrawing());
  assert.deepEqual(structuredClone({
    records: modeDrawing(state, 'circular').circularRecordList.value,
    positions: state.activeDrawing().adv.multi_record_positions,
    discovery: { ...modeDrawing(state, 'circular').circularRecordDiscovery }
  }), circularStateBeforeRoundTrip);
  assert.equal(
    workerMessages.filter(({ type }) => type === 'helper').length,
    helperRequestsBeforeRoundTrip
  );
  assert.equal(modeDrawing(state, 'circular').circularRecordList.value[0].identity, 'record-1');

  // A replaced or removed source is read again when Circular is shown, even if the change happened while it was inactive.
  installShownDrawing('linear');
  setCircularFile(new AuditFile(['LOCUS       audit 10 bp DNA circular\n//\n'], 'replacement.gb'));
  await runner.refreshCircularRecordOrder(undefined, state.activeDrawing());
  installShownDrawing('circular');
  await runner.refreshCircularRecordOrder({ automatic: true }, state.activeDrawing());
  assert.equal(modeDrawing(state, 'circular').circularRecordList.value[0].record_id, 'audit');
  installShownDrawing('linear');
  setCircularFile(null);
  await runner.refreshCircularRecordOrder(undefined, state.activeDrawing());
  assert.equal(modeDrawing(state, 'circular').circularRecordList.value[0].record_id, 'audit');
  installShownDrawing('circular');
  await runner.refreshCircularRecordOrder({ automatic: true }, state.activeDrawing());
  assert.deepEqual(modeDrawing(state, 'circular').circularRecordList.value, []);
  assert.deepEqual(state.activeDrawing().adv.multi_record_positions, []);
  setCircularFile(fallbackPrimary);

  state.activeDrawing().form.multi_record_canvas = true;
  modeDrawing(state, 'circular').files.c_depth = bound(null);
  state.activeDrawing().adv.circular_track_slots_enabled = false;
  state.activeDrawing().adv.circular_track_slots = [];
  modeDrawing(state, 'circular').circularRecordList.value = [];
  Object.assign(modeDrawing(state, 'circular').circularRecordDiscovery, {
    status: 'idle',
    error: '',
    inputType: '',
    sources: []
  });
  setCircularFile(fallbackPrimary);
  let releaseDiscovery;
  workerHelperResponses.push(new Promise((resolve) => {
    releaseDiscovery = () => resolve({
      ok: true,
      result: {
        records: [{ selector: '#1', record_id: 'audit', record_length: 10 }]
      }
    });
  }));
  const helperRequestsBeforeCancel = workerMessages.filter(({ type }) => type === 'helper').length;
  const canceledRun = runner.runAnalysis();
  for (
    let turn = 0;
    turn < 12 && workerMessages.filter(({ type }) => type === 'helper').length === helperRequestsBeforeCancel;
    turn += 1
  ) {
    await Promise.resolve();
  }
  await runner.cancelRunAnalysis();
  releaseDiscovery();
  assert.deepEqual(await canceledRun, { status: 'canceled' });
  assert.equal(state.processing.value, false);
  assert.equal(state.generationCancelRequested.value, false);
  assert.equal(
    workerMessages.filter(({ type }) => type === 'run').length,
    workerRunCountBeforeDiscoveryFailure
  );

  state.activeDrawing().form.multi_record_canvas = true;
  setCircularFile(stalePrimary);
  state.activeDrawing().adv.multi_record_positions.splice(
    0,
    state.activeDrawing().adv.multi_record_positions.length,
    { selector: '#1', row: 2 }
  );
  modeDrawing(state, 'circular').circularRecordList.value = [listedRecord('preserved-record', 123)];
  state.semanticFileWatchersSuppressed.value = false;
  state.sessionImportRollbackInProgress.value = false;
  let releaseStaleDiscovery;
  workerHelperResponses.push(new Promise((resolve) => {
    releaseStaleDiscovery = () => resolve({
      ok: false,
      error: { name: 'Error', message: 'injected restored-file reparse failure' }
    });
  }));
  const staleDiscovery = runner.refreshCircularRecordOrder(undefined, state.activeDrawing());
  await Promise.resolve();
  const discoveryStateBeforeRollback = {
    records: structuredClone(modeDrawing(state, 'circular').circularRecordList.value),
    positions: structuredClone(state.activeDrawing().adv.multi_record_positions),
    discovery: { ...modeDrawing(state, 'circular').circularRecordDiscovery }
  };
  state.semanticFileWatchersSuppressed.value = true;
  await runner.refreshCircularRecordOrder({ suppress: true }, state.activeDrawing());
  releaseStaleDiscovery();
  await staleDiscovery;
  assert.deepEqual(modeDrawing(state, 'circular').circularRecordList.value, discoveryStateBeforeRollback.records);
  assert.deepEqual(state.activeDrawing().adv.multi_record_positions, discoveryStateBeforeRollback.positions);
  assert.deepEqual(modeDrawing(state, 'circular').circularRecordDiscovery, discoveryStateBeforeRollback.discovery);
  state.semanticFileWatchersSuppressed.value = false;
});

test('a ring comparison file is read once, and a failed read is read again', async () => {
  const file = new File(['>ring\nacgt\n'], 'ring.fasta', { type: 'text/plain' });
  const operations = [];
  const replies = [
    { result: { error: new Error('injected ring read failure') } },
    { result: { fasta: '>ring\nACGT\n', recordLabel: 'Ring' } },
    { result: { fasta: '>ring\nACGT\n', recordLabel: 'Ring', recordIds: ['ring'], genomeFiles: [{ length: 4, recordLengths: [4] }] } }
  ];
  const runHelperOperation = async (operation) => {
    operations.push(operation);
    return replies.shift();
  };
  await assert.rejects(readComparisonSequence({ files: [file], runHelperOperation }), /injected ring read failure/);
  // A reply without its record IDs (or without each file's lengths) is a protocol error, and it is not kept either.
  await assert.rejects(readComparisonSequence({ files: [file], runHelperOperation }), { code: 'HELPER_PROTOCOL' });
  const read = await readComparisonSequence({ files: [file], runHelperOperation });
  assert.deepEqual(read, {
    fasta: '>ring\nACGT\n', recordLabel: 'Ring', recordIds: ['ring'], genomeFiles: [{ length: 4, recordLengths: [4] }]
  });
  assert.equal(await readComparisonSequence({ files: [file], runHelperOperation }), read);
  assert.deepEqual(operations, ['readComparisonSequence', 'readComparisonSequence', 'readComparisonSequence']);
});

// OV-1587: one per-File reading owner. Every reader gets the same frozen value,
// and the ring FASTA's SHA-256 is kept beside it, digested once per File.
test('a ring comparison read is frozen and its digest is kept with it', async () => {
  const file = new File(['>ring\nacgt\n'], 'ring-digest.fasta', { type: 'text/plain' });
  const operations = [];
  const runHelperOperation = async (operation) => {
    operations.push(operation);
    return { result: { fasta: '>ring\nACGT\n', recordLabel: null, recordIds: ['ring'], genomeFiles: [{ length: 4, recordLengths: [4] }] } };
  };
  const read = await readComparisonSequence({ files: [file], runHelperOperation });
  assert.equal(Object.isFrozen(read), true);
  assert.throws(() => { read.fasta = '>other\n'; }, TypeError);
  assert.equal((await readComparisonSequence({ files: [file], runHelperOperation })).fasta, '>ring\nACGT\n');
  const expected = await sha256Text('>ring\nACGT\n');
  const watched = await watchLosatPreparation(async () => {
    assert.equal(await readComparisonSequenceDigest({ files: [file], runHelperOperation }), expected);
    assert.equal(await readComparisonSequenceDigest({ files: [file], runHelperOperation }), expected);
  });
  assert.equal(watched.fastaDigests, 1);
  assert.deepEqual(operations, ['readComparisonSequence']);
});

// A ring's comparison genome may be several files: one read sends them all, in
// order, and is kept for that list (not for each file alone, nor another order).
test('a ring comparison genome of two files is one read of both, kept for those files in that order', async () => {
  const first = new File(['>chr1\nacgt\n'], 'chr1.fasta', { type: 'text/plain' });
  const second = new File(['>chr2\nttgg\n'], 'chr2.fasta', { type: 'text/plain' });
  const sent = [];
  const runHelperOperation = async (operation, payload) => {
    sent.push([operation, await Promise.all(payload.files.map(async ({ role, blob }) => [role, await blob.text()]))]);
    return { result: {
      fasta: `>read${sent.length}\nACGT\n`, recordLabel: null, recordIds: [`read${sent.length}`],
      genomeFiles: payload.files.map(() => ({ length: 4, recordLengths: [4] }))
    } };
  };
  const both = await readComparisonSequence({ files: [first, second], runHelperOperation });
  assert.deepEqual(sent, [['readComparisonSequence', [['source', '>chr1\nacgt\n'], ['source', '>chr2\nttgg\n']]]]);
  assert.equal(await readComparisonSequence({ files: [first, second], runHelperOperation }), both);
  assert.notEqual(await readComparisonSequence({ files: [first], runHelperOperation }), both);
  assert.notEqual(await readComparisonSequence({ files: [second, first], runHelperOperation }), both);
  assert.equal(sent.length, 3);
});

// #1777 review 9 and round 12 (OV-3203, OV-3204): Save Raw LOSAT TSV reads the shown drawing's results; both
// panels' "Clear N LOSAT results" count the searches Clear removes, the whole cache (PD-OI-022), a search only a
// deleted drawing names included (round 15, OV-3211).
test('Save Raw LOSAT TSV reads the shown drawing\'s results and both Clear buttons count the project\'s searches', () => {
  const html = readFileSync(new URL('../../gbdraw/web/index.html', import.meta.url), 'utf8');
  assert.equal(html.includes('losatCacheInfo'), false);
  assert.equal(html.match(/tn\(clearableLosatResultCount, 'Clear \{count\} LOSAT result'/g)?.length, 2);
  assert.equal(html.match(/@click="clearLosatCache"[^>]*:disabled="!semanticMutationAvailable \|\| \(clearableLosatResultCount === 0\)"/g)?.length, 2);
  assert.match(html, /@click="downloadLosatCache"[^>]*:disabled="losatResultRows\.length === 0" :aria-describedby="losatResultRows\.length === 0 \? 'circular-losat-cache-reason' : null"/);
  // Round 13 (OV-3206): each reason is true. The Circular hint gives Save Raw's reason: dev's hint while the project
  // has no LOSAT results (Clear is disabled too and points at it), else that the shown drawing has none (Clear is
  // enabled and points at nothing). The Linear hint gives Clear's reason.
  assert.ok(html.includes(`<p v-if="losatResultRows.length === 0" id="circular-losat-cache-reason" class="ui-hint">{{ clearableLosatResultCount === 0 ? t('No LOSAT results yet') : t('No LOSAT results in this drawing') }}</p>`));
  assert.ok(html.includes(`<p v-if="clearableLosatResultCount === 0" id="linear-losat-cache-reason" class="ui-hint">{{ t('No LOSAT results yet') }}</p>`));
  assert.equal(html.match(/:aria-describedby="clearableLosatResultCount === 0 \? '(circular|linear)-losat-cache-reason' : null"/g)?.length, 2);
});

test('neutral conservation replay delegates lazy resources to the shared reader', async () => {
  const resourceMetrics = [];
  const previousTestHooks = globalThis.__GBDRAW_TEST_HOOKS__;
  globalThis.__GBDRAW_TEST_HOOKS__ = {
    onStructuralMetric(metric) {
      resourceMetrics.push({ ...metric });
    }
  };
  const referenceText = [
    'LOCUS       reference-a 8 bp DNA',
    'ORIGIN',
    '        1 acgtacgt',
    '//',
    ''
  ].join('\n');
  const referenceFasta = '>reference-a\nACGTACGT\n';
  // The LOSAT query of each ring: the shared Python reader's canonical FASTA.
  const comparisonTexts = [
    '>comparison-b\nTTTTAAAA\n',
    '>comparison-c\nAAAATTTT\n',
    '>comparison-d\nCCCCGGGG\n'
  ];
  // Ring files as uploaded: a non-canonical FASTA, a GenBank flat file and a
  // canonical FASTA. Each gives the same raw key as its canonical FASTA (D12).
  const comparisonFileTexts = [
    '>comparison-b first ring\r\nttttaaaa\r\n',
    [
      'LOCUS       comparison-c 8 bp DNA',
      'DEFINITION  comparison c.',
      'ORIGIN',
      '        1 aaaatttt',
      '//',
      ''
    ].join('\n'),
    comparisonTexts[2]
  ];
  const resourceTable = adoptCurrentSessionResources({
    'reference-a': encodedResource('genbank', 'reference-a.gb', referenceText),
    'comparison-b': encodedResource(
      'conservation-fasta-file',
      'comparison-b.fasta',
      comparisonFileTexts[0]
    ),
    'comparison-c': encodedResource(
      'conservation-fasta-file',
      'comparison-c.gbk',
      comparisonFileTexts[1]
    ),
    'comparison-d': encodedResource(
      'conservation-fasta-file',
      'comparison-d.fasta',
      comparisonFileTexts[2]
    )
  });
  const referenceView = createSessionResourceFileView(resourceTable, 'reference-a');
  const comparisonViews = ['comparison-b', 'comparison-c', 'comparison-d'].map(
    (resourceId) => createSessionResourceFileView(resourceTable, resourceId)
  );
  const lazyFields = ['text', 'arrayBuffer', 'data', 'resourceId'];
  for (const file of [referenceView, ...comparisonViews]) {
    assert.equal(Object.isFrozen(file), true);
    lazyFields.forEach((field) => assert.equal(Object.hasOwn(file, field), false));
  }

  // The ring drawing: three LOSAT rings on one reference (also a copy of it, below).
  const setUpRingDrawing = (drawing) => {
    drawing.inputType.value = 'gb';
    drawing.form.prefix = 'lazy-circular';
    drawing.form.show_depth = false;
    drawing.form.multi_record_canvas = false;
    drawing.adv.circular_track_slots_enabled = false;
    drawing.adv.circular_track_slots = [];
    Object.assign(drawing.adv, {
      min_bitscore: 0,
      evalue: '1e-2',
      identity: 0,
      alignment_length: 0
    });
    Object.assign(drawing.circularConservation, {
      enabled: true,
      source: 'losat',
      losat_program: 'blastn',
      subject_gencode: 1,
      reference: 'auto',
      labels: '',
      ring_width: 5,
      ring_gap: 2
    });
    drawing.circularConservation.series.splice(
      0,
      drawing.circularConservation.series.length,
      ...comparisonViews.map((file, index) => ({
        sourceKey: `${file.name}|${file.size}|0|${index}`,
        fileName: file.name,
        sourceIndex: index,
        label: ['comparison-b', 'comparison-c', 'comparison-d'][index],
        color: ['#4e79a7', '#e15759', '#59a14f'][index],
        losat_gencode: 1
      }))
    );
    setCircularFile(referenceView, drawing);
    drawing.files.c_depth = bound(null);
    drawing.files.c_conservation_blasts = bound([]);
    drawing.files.c_conservation_fastas = bound(comparisonViews);
    drawing.files.c_conservation_sequence_sources = bound([]);
    drawing.annotationSets.splice(0);
    drawing.circularRecordList.value = [];
    Object.assign(drawing.circularRecordDiscovery, {
      status: 'idle', error: '', inputType: '', sources: []
    });
  };
  installShownDrawing('circular');
  setUpRingDrawing(modeDrawing(state, 'circular'));

  const subjectCanonicalHash = await sha256Text(referenceFasta);
  const cacheRows = [
    'comparison-b\treference-a\t100\t8\t0\t0\t1\t8\t1\t8\t1e-20\t40\n',
    'comparison-c\treference-a\t95\t8\t0\t0\t1\t8\t1\t8\t1e-18\t38\n',
    'comparison-d\treference-a\t90\t8\t0\t0\t1\t8\t1\t8\t1e-16\t36\n'
  ];
  const cacheEntries = [];
  for (let index = 0; index < comparisonTexts.length; index += 1) {
    const queryCanonicalHash = await sha256Text(comparisonTexts[index]);
    const cachePayload = {
      cacheSchema: 2,
      identityKind: 'nucleotide',
      program: 'blastn',
      outfmt: '6',
      args: ['--task', 'megablast'],
      queryCanonicalHash,
      subjectCanonicalHash,
      flow: 'circular-conservation'
    };
    const cacheKey = await sha256Text(JSON.stringify(cachePayload));
    cacheEntries.push([cacheKey, {
      schema: 2,
      kind: 'raw-losat',
      identityKind: 'nucleotide',
      text: cacheRows[index],
      program: 'blastn',
      flow: 'circular-conservation',
      outfmt: '6',
      args: ['--task', 'megablast'],
      queryCanonicalHash,
      subjectCanonicalHash,
      runtime: { kind: 'losat', source: 'wasm', version: LOSAT_WASM_BUILD.version, program: 'blastn' }
    }]);
  }
  state.losatCache.value = new Map(cacheEntries);
  state.losatDerivedCache.value = new Map();

  let failLateArtifactAdoption = false;
  let committedRenderRequest = null;
  let committedCanonicalSession = null;
  const recordDisplayRows = { value: [] };
  // The record selector's reads (OV-3227): a card with a selector reads the record it names, a card without one every
  // record of its file (`wholeFileRecords`); a card in `unreadCards` (a loaded Session's) is read by
  // `readLinearCardRecords`.
  const wholeFileRecords = new Map();
  const unreadCards = new Set();
  const runner = wireGeneratedArtifactRuntimeOwner(createRunAnalysis({
    ...generatedArtifactHandleOptions,
    state,
    recordDisplayRows,
    writeSessionDocumentFile,
    serializeCanonicalFiles: (_snapshot, _catalog, drawing) => serializeActiveRenderFiles(state, drawing),
    adoptCanonicalRenderArtifacts: (canonical) => {
      if (failLateArtifactAdoption) {
        throw new Error('injected LOSAT late artifact adoption failure');
      }
      committedRenderRequest = canonical.renderRequest;
      committedCanonicalSession = canonical;
    },
    getCommittedCanonicalSession: () => (committedCanonicalSession
      ? { ...committedCanonicalSession, renderRequest: committedRenderRequest }
      : null),
    prepareLinearRecordCatalog: async () => ({
      catalog: { mode: 'linear', status: 'ready', records: [] },
      error: ''
    }),
    linearCardRecord: (_drawing, card) => cardRecordFacts(unreadCards.has(card.uid) ? []
      : card.region_record_id ? [{ recordId: card.region_record_id, recordLength: 8 }] : wholeFileRecords.get(card.uid) ?? []),
    readLinearCardRecords: async () => { unreadCards.clear(); },
    prepareCandidateCommit: ({
      generationResponse,
      catalogAdmission,
      featureColorOverrides
    }) => ({
      results: structuredClone(generationResponse.results),
      featureState: featureStateFromCatalog(catalogAdmission.catalog),
      featureColorOverrides: structuredClone(featureColorOverrides)
    }),
    resetPreviewViewport: () => {}
  }));

  let losatCalls = 0;
  let capturedSequences = [];
  const previousLosatExecutor = globalThis.__GBDRAW_LOSAT_EXECUTOR__;
  globalThis.__GBDRAW_LOSAT_EXECUTOR__ = async (jobs, options) => {
    losatCalls += 1;
    capturedSequences = Array.from(options.sequences.values());
    return jobs.map((job) => ({
      cacheKey: job.cacheKey,
      text: 'MIDDLE\tTHIRD\t100\t8\t0\t0\t1\t8\t1\t8\t1e-10\t20\n'
    }));
  };

  try {
    const circularResult = result('lazy-circular.svg', 'lazy-circular');
    const helperMessageCount = workerMessages.filter(({ type }) => type === 'helper').length;
    const ringSequenceReply = (index) => ({
      ok: true,
      result: {
        fasta: comparisonTexts[index],
        format: index === 1 ? 'genbank' : 'fasta',
        recordLabel: null,
        recordIds: [comparisonTexts[index].slice(1).split('\n')[0]],
        genomeFiles: [{ length: 8, recordLengths: [8] }]
      }
    });
    workerHelperResponses.push(ringSequenceReply(0));
    // Add Seq reads the first ring file with the same reader before Generate,
    // which reuses that read (OV-1545).
    assert.equal((await readComparisonSequence({ files: [comparisonViews[0]] })).fasta, comparisonTexts[0]);
    // Generate reads the other rings for their keys, then plans the ring scales: one-record files keep LOSAT's
    // own scale.
    workerHelperResponses.push(ringSequenceReply(1), ringSequenceReply(2));
    workerHelperResponses.push({ ok: true, result: { scales: [null, null, null], warnings: [] } });
    workerHelperResponses.push(ringTablesUnchanged);
    workerResponses.push(response(circularResult, validCatalog(circularResult.name)));
    // A Circular Generate replaces only its own drawing's cache-info rows (#1777 review 8, OV-3201).
    const linearId = modeDrawing(state, 'linear').id;
    const linearRows = ['a->b', 'b->c'].map((edgeKey, index) => ({
      key: `linear-pair-${index}`, filename: `pair-${index}.tsv`, display: true, edgeKey, drawingId: linearId
    }));
    state.losatCacheInfo.value = [...linearRows, { key: 'stale-ring', filename: 'old.tsv', display: true, drawingId: state.activeDrawing().id }];
    // Each row names a cached search, as Generate and Load write them.
    state.losatCache.value = new Map([...state.losatCache.value, ...linearRows.map(({ key }) => [key, { kind: 'raw-losat', text: '' }])]);
    assert.deepEqual(
      await runner.runAnalysis(),
      { status: 'ok' },
      JSON.stringify(state.errorLog.value)
    );
    assert.equal(losatCalls, 0, 'verified circular cache entries must prevent LOSAT execution');
    assert.deepEqual(state.losatCacheInfo.value.map((row) => [row.key.startsWith('linear-pair') || row.key === 'stale-ring' ? row.key : 'ring', row.drawingId]),
      [['linear-pair-0', 'linear'], ['linear-pair-1', 'linear'], ...cacheEntries.map(() => ['ring', state.activeDrawing().id])]);

    // #1777 review 9 (OV-3203, OV-3204): Save Raw LOSAT TSV saves the LOSAT results of the drawing a panel shows,
    // one per cache key, through one reader (`losatResultRowsOf`); Clear counts the searches it removes, the whole
    // cache (`clearableLosatResultCount`, OV-3211), and clears them all, as on dev (PD-OI-022).
    const ringKeys = state.losatCacheInfo.value.filter(({ drawingId }) => drawingId === 'circular').map(({ key }) => key);
    assert.equal(ringKeys.length, 3);
    const panelOf = async (drawingId) => {
      installShownDrawing(drawingId);
      downloadedFilenames.splice(0);
      await runner.downloadLosatCache(state.activeDrawing());
      return { keys: runner.losatResultRowsOf(state.activeDrawing()).map(({ key }) => key), downloads: [...downloadedFilenames] };
    };
    const ringFiles = (await panelOf('circular')).downloads;
    assert.deepEqual(await panelOf('circular'), { keys: ringKeys, downloads: ringFiles },
      'the Circular panel says 3 and saves only its 3 rings, not the Linear pairs');
    assert.equal(new Set(ringFiles).size, 3, JSON.stringify(ringFiles));
    assert.deepEqual((await panelOf(linearId)).keys, linearRows.map(({ key }) => key), 'the Linear panel says 2');
    assert.equal(runner.clearableLosatResultCount(), 5, 'Clear removes the 2 pairs and the 3 ring searches');
    // Two rings of one genome in one drawing share a cache key: one result (as on dev, where Save writes one entry).
    const generatedRows = state.losatCacheInfo.value;
    state.losatCacheInfo.value = [...generatedRows, { ...generatedRows.at(-1), filename: 'same-genome-again.tsv' }];
    assert.deepEqual(await panelOf('circular'), { keys: ringKeys, downloads: ringFiles }, 'one row per cache key');
    state.losatCacheInfo.value = generatedRows;
    // A copy of the drawing, Generated: each drawing's panel says 3 and saves its 3 files once.
    const copy = createDrawingState('circular-copy', 'circular');
    state.addDrawing(copy, 'Circular copy', 2);
    try {
      setUpRingDrawing(copy);
      installShownDrawing(copy.id);
      // The copy's rings reuse the ring reads and tables of the first Generate: no helper runs.
      const copyResult = result('lazy-circular-copy.svg', 'lazy-circular-copy');
      workerResponses.push(response(copyResult, validCatalog(copyResult.name)));
      assert.deepEqual(await runner.runAnalysis(), { status: 'ok' }, JSON.stringify(state.errorLog.value));
      assert.equal(losatCalls, 0);
      assert.equal(state.losatCacheInfo.value.length, 8, 'each drawing keeps its own rows');
      for (const drawingId of [copy.id, 'circular']) {
        assert.deepEqual(await panelOf(drawingId), { keys: ringKeys, downloads: ringFiles }, drawingId);
      }
      // Two drawings hold the 3 ring searches: Clear still counts 5 searches, not 8 rows.
      assert.equal(runner.clearableLosatResultCount(), 5);
      // Clear from the copy's panel clears the project, as on dev: every drawing's rows and every search.
      const kept = { cache: state.losatCache.value, rows: state.losatCacheInfo.value, derived: state.losatDerivedCache.value };
      installShownDrawing(copy.id);
      runner.clearLosatCache();
      assert.deepEqual([state.losatCacheInfo.value, state.losatCache.value.size, runner.clearableLosatResultCount()], [[], 0, 0]);
      assert.deepEqual([...await Promise.all(['circular', linearId].map(async (drawingId) => (await panelOf(drawingId)).keys))], [[], []]);
      state.losatCache.value = kept.cache;
      state.losatCacheInfo.value = kept.rows;
      state.losatDerivedCache.value = kept.derived;
    } finally {
      installShownDrawing('circular');
      state.removeDrawing(copy.id, { takeSlot: () => null });
    }
    // Every ring file goes through the one Python reader (no JavaScript reader), once; Python plans the ring
    // scales once and then puts the ring rows in the search frame.
    assert.deepEqual(
      workerMessages.filter(({ type }) => type === 'helper').slice(helperMessageCount)
        .map(({ operation }) => operation),
      ['readComparisonSequence', 'readComparisonSequence', 'readComparisonSequence', 'planRingScales',
        'convertRingTablesToSearchFrame']
    );
    assert.equal(workerHelperResponses.length, 0);
    assert.deepEqual(
      resourceMetrics.filter(({ name }) => name === 'resourceTextReadCount')
        .map(({ resourceId }) => resourceId),
      ['reference-a']
    );
    assert.deepEqual(
      resourceMetrics.filter(({ name }) => name === 'resourceByteReadCount')
        .map(({ resourceId }) => resourceId)
        .filter((resourceId) => resourceId.startsWith('comparison-')),
      ['comparison-b', 'comparison-c', 'comparison-d']
    );
    assert.deepEqual(
      workerMessages.filter(({ type }) => type === 'run').at(-1)
        .payload.request.diagramOptions.conservationLabels,
      ['comparison-b', 'comparison-c', 'comparison-d']
    );
    assert.deepEqual(
      workerMessages.filter(({ type }) => type === 'run').at(-1)
        .payload.request.diagramOptions.conservationColors,
      ['#4e79a7', '#e15759', '#59a14f']
    );

    // A second Generate on the same files reads no ring file again (OV-1545):
    // neither the reference nor a comparison, and it sends the same LOSAT rows.
    const firstConservation = workerMessages.filter(({ type }) => type === 'run').at(-1)
      .payload.request.diagramOptions;
    const helperMessagesBeforeRepeat = workerMessages.filter(({ type }) => type === 'helper').length;
    resourceMetrics.splice(0);
    const repeatResult = result('lazy-circular-repeat.svg', 'lazy-circular-repeat');
    workerResponses.push(response(repeatResult, validCatalog(repeatResult.name)));
    const repeatWatch = await watchLosatPreparation(async () => assert.deepEqual(
      await runner.runAnalysis(),
      { status: 'ok' },
      JSON.stringify(state.errorLog.value)
    ));
    // OV-1565: the ring sequences' digests are kept with their FASTA, so no ring FASTA is digested again.
    assert.equal(repeatWatch.fastaDigests, 0);
    assert.equal(losatCalls, 0);
    assert.equal(
      workerMessages.filter(({ type }) => type === 'helper').length,
      helperMessagesBeforeRepeat
    );
    assert.deepEqual(
      resourceMetrics.filter(({ name, resourceId }) => name === 'resourceTextReadCount'
        || (name === 'resourceByteReadCount' && resourceId.startsWith('comparison-'))),
      []
    );
    assert.deepEqual(
      workerMessages.filter(({ type }) => type === 'run').at(-1)
        .payload.request.diagramOptions.conservationLabels,
      firstConservation.conservationLabels
    );

    // OV-1599: a reference that starts with a UTF-8 byte order mark. The browser reader declines it, and the
    // subject database of the rings is the records the Python loader reads: every record (OV-2354).
    const bomReference = new AuditFile([new Uint8Array([0xef, 0xbb, 0xbf]), ['BOMA', 'BOMB'].map((id, index) => [
      `LOCUS       ${id} 8 bp DNA`, 'ORIGIN', `        1 ${['acgtacgt', 'ccccgggg'][index]}`, '//'
    ].join('\n')).join('\n') + '\n'], 'bom-reference.gb');
    setCircularFile(bomReference);
    const helpersBeforeBomReference = workerMessages.filter(({ type }) => type === 'helper').length;
    // The record list of the new reference, the subject database, then the ring scales: the reference of two
    // records gives every ring the search space 16 x 8.
    workerHelperResponses.push(
      { ok: true, result: { records: [
        { selector: '#1', record_id: 'BOMA', record_length: 8 }, { selector: '#2', record_id: 'BOMB', record_length: 8 }
      ] } },
      { ok: true, result: { records: [
        { id: 'BOMA', aliases: ['BOMA'], sequence: 'ACGTACGT' }, { id: 'BOMB', aliases: ['BOMB'], sequence: 'CCCCGGGG' }
      ] } },
      { ok: true, result: { scales: ['128', '128', '128'], warnings: [] } },
      ringTablesUnchanged
    );
    const bomReferenceResult = result('bom-reference.svg', 'bom-reference');
    workerResponses.push(response(bomReferenceResult, validCatalog(bomReferenceResult.name)));
    assert.deepEqual(await runner.runAnalysis(), { status: 'ok' }, JSON.stringify(state.errorLog.value));
    assert.deepEqual(
      workerMessages.filter(({ type }) => type === 'helper').slice(helpersBeforeBomReference)
        .map(({ operation, payload }) => [operation, payload.format, payload.files?.map(({ role }) => role)]),
      [['listSequenceRecords', 'genbank', ['source']], ['readAllSequenceRecords', 'genbank', ['source']],
        ['planRingScales', undefined, undefined],
        ['convertRingTablesToSearchFrame', undefined, undefined]]
    );
    assert.equal(workerHelperResponses.length, 0);
    assert.ok(capturedSequences.some((fasta) => fasta.includes('>BOMA\nACGTACGT\n') && fasta.includes('>BOMB\nCCCCGGGG\n')),
      JSON.stringify(capturedSequences));
    // The new subject searched once; the Linear section below counts its own searches.
    assert.equal(losatCalls, 1);
    losatCalls = 0;
    setCircularFile(referenceView);

    // OV-3260 (D-200): no Generate and no Undo removes a cached search, so a later Generate reuses a search its last
    // Generate did not use. Each Generate below sends the rings it draws.
    const circularDrawing = modeDrawing(state, 'circular');
    const ringSeries = [...circularDrawing.circularConservation.series];
    // The copy above left the project without its rows (the drawing list's Delete holds them): no other row names the
    // ring searches.
    state.losatCacheInfo.value = state.losatCacheInfo.value.filter(({ drawingId }) => drawingId !== 'circular-copy');
    const shownRingKeys = () => state.losatCacheInfo.value.filter(({ drawingId }) => drawingId === circularDrawing.id).map(({ key }) => key);
    /**
     * A Generate of the reference-a drawing with the comparison rings `indexes`; `scales` is Python's reply for a
     * ring set whose scales were not planned yet (one-record files: LOSAT's own scale). A ring table conversion
     * held for the frame leaves its reply unread, so each Generate starts from no queued reply.
     * @param {number[]} indexes @param {(string | null)[] | null} [scales]
     */
    const generateRings = async (indexes, scales = null) => {
      circularDrawing.circularConservation.series.splice(0, circularDrawing.circularConservation.series.length,
        ...indexes.map((index) => ringSeries[index]));
      circularDrawing.files.c_conservation_fastas = bound(indexes.map((index) => comparisonViews[index]));
      workerHelperResponses.splice(0);
      if (scales) workerHelperResponses.push({ ok: true, result: { scales, warnings: [] } });
      workerHelperResponses.push(ringTablesUnchanged);
      const ringResult = result(`rings-${indexes.join('')}.svg`, `rings-${indexes.join('')}`);
      workerResponses.push(response(ringResult, validCatalog(ringResult.name)));
      assert.deepEqual(await runner.runAnalysis(), { status: 'ok' }, JSON.stringify(state.errorLog.value));
    };
    // Back on reference-a after the BOM reference's Generate: its ring searches are still cached.
    await generateRings([0, 1, 2]);
    assert.equal(losatCalls, 0, 'reference-a\'s ring searches outlive the BOM reference\'s Generate');
    // Ring d removed, then added again. The scales of the two rings are planned once.
    await generateRings([0, 1], [null, null]);
    const twoRings = shownRingKeys();
    await generateRings([0, 1, 2]);
    assert.equal(losatCalls, 0, 'ring d\'s search outlives the Generate without it');
    const ringDKey = shownRingKeys().find((key) => !twoRings.includes(key));
    // Through an Undo: ring d searched by a Generate whose Undo then runs, and Generated again.
    state.losatCache.value = new Map([...state.losatCache.value].filter(([key]) => key !== ringDKey));
    await generateRings([0, 1]);
    const beforeRingD = await artifactSnapshots.captureGeneratedArtifactHandle(state.activeDrawing());
    await generateRings([0, 1, 2]);
    assert.equal(losatCalls, 1, 'ring d searched');
    await artifactSnapshots.restoreGeneratedArtifactHandle(beforeRingD);
    await generateRings([0, 1, 2]);
    assert.equal(losatCalls, 1, 'the Undo of the Generate that searched ring d keeps its search');
    losatCalls = 0;

    // OV-2818: with Show Comparison Rings off, Generate reads no ring file, runs no LOSAT search (the cache is
    // empty) and sends no ring; the ring files and rows stay.
    const cachedRings = state.losatCache.value;
    state.losatCache.value = new Map();
    state.activeDrawing().circularConservation.enabled = false;
    const helpersBeforeRingsOff = workerMessages.filter(({ type }) => type === 'helper').length;
    const ringsOffResult = result('rings-off.svg', 'rings-off');
    workerResponses.push(response(ringsOffResult, validCatalog(ringsOffResult.name)));
    assert.deepEqual(await runner.runAnalysis(), { status: 'ok' }, JSON.stringify(state.errorLog.value));
    assert.equal(losatCalls, 0);
    assert.deepEqual(workerMessages.filter(({ type }) => type === 'helper').slice(helpersBeforeRingsOff)
      .map(({ operation }) => operation), []);
    const ringsOffOptions = workerMessages.filter(({ type }) => type === 'run').at(-1).payload.request.diagramOptions;
    assert.equal(ringsOffOptions.conservationLabels, undefined, JSON.stringify(ringsOffOptions));
    assert.equal(state.activeDrawing().circularConservation.enabled, false);
    assert.equal(state.activeDrawing().circularConservation.series.length, 3);
    assert.equal(modeDrawing(state, 'circular').files.c_conservation_fastas.length, 3);
    state.losatCache.value = cachedRings;
    state.activeDrawing().circularConservation.enabled = true;

    const invalidFile = Object.freeze({
      name: 'missing.gb', type: 'text/plain', size: 0, lastModified: 0
    });
    setCircularFile(invalidFile);
    modeDrawing(state, 'circular').circularRecordList.value = [
      listedRecord('missing', 1)
    ];
    Object.assign(modeDrawing(state, 'circular').circularRecordDiscovery, {
      status: 'ready',
      error: '',
      inputType: 'gb',
      sources: [{ uid: 'record-1', primaryFile: bound(invalidFile), pairedFile: null, error: '' }]
    });
    assert.equal((await runner.runAnalysis()).status, 'error');
    assert.match(
      state.errorLog.value?.summary || '',
      /Input resource preparation failed/
    );

    resourceMetrics.splice(0);
    const linearTable = adoptCurrentSessionResources({
      multi: encodedResource('genbank', 'alpha.gbk.tsv', [
        'LOCUS       FIRST 8 bp DNA',
        'ORIGIN',
        '        1 aaaaaaaa',
        '//',
        'LOCUS       SECOND 8 bp DNA',
        'ORIGIN',
        '        1 aaccggtt',
        '//',
        ''
      ].join('\n')),
      middle: encodedResource('genbank', 'middle.gbk', [
        'LOCUS       MIDDLE 8 bp DNA',
        'ORIGIN',
        '        1 ccccgggg',
        '//',
        ''
      ].join('\n')),
      third: encodedResource('genbank', 'third.gb', [
        'LOCUS       THIRD 8 bp DNA',
        'ORIGIN',
        '        1 ttttaaaa',
        '//',
        ''
      ].join('\n'))
    });
    const multiView = createSessionResourceFileView(linearTable, 'multi');
    const middleView = createSessionResourceFileView(linearTable, 'middle');
    const thirdView = createSessionResourceFileView(linearTable, 'third');
    installShownDrawing('linear');
    modeDrawing(state, 'linear').inputType.value = 'gb';
    state.activeDrawing().form.prefix = 'lazy-linear';
    state.activeDrawing().form.show_depth = false;
    state.activeDrawing().form.linear_track_layout = 'middle';
    state.activeDrawing().adv.linear_track_slots_enabled = false;
    state.activeDrawing().adv.linear_track_slots = [];
    state.activeDrawing().adv.comparison_height = null;
    state.activeDrawing().linearRecordLayoutEnabled.value = false;
    state.activeDrawing().linearRecordRows.splice(0);
    state.activeDrawing().linearSeqs.splice(0, state.activeDrawing().linearSeqs.length, {
      uid: 'multi',
      gb: bound(multiView),
      gff: null,
      fasta: null,
      depth: null,
      blast: null,
      losat_gencode: 1,
      losat_filename: '',
      definition: '',
      record_subtitle: '',
      region_record_id: 'SECOND',
      region_start: 2,
      region_end: 5,
      region_reverse: true,
      windows: []
    }, {
      uid: 'middle',
      gb: bound(middleView),
      gff: null,
      fasta: null,
      depth: null,
      blast: null,
      losat_gencode: 1,
      losat_filename: '',
      definition: '',
      record_subtitle: '',
      region_record_id: 'MIDDLE',
      region_start: null,
      region_end: null,
      region_reverse: false,
      windows: []
    }, {
      uid: 'third',
      gb: bound(thirdView),
      gff: null,
      fasta: null,
      depth: null,
      blast: null,
      losat_gencode: 1,
      losat_filename: '',
      definition: '',
      record_subtitle: '',
      region_record_id: 'THIRD',
      region_start: null,
      region_end: null,
      region_reverse: false,
      windows: []
    });
    Object.assign(state.activeDrawing().linearComparisonPlan, {
      mode: 'selected',
      defaultSource: 'losat'
    });
    state.activeDrawing().linearComparisonPlan.edges.splice(0, state.activeDrawing().linearComparisonPlan.edges.length, {
      id: 'lazy-linear-edge',
      queryUid: 'middle',
      subjectUid: 'third',
      included: true,
      source: 'losat',
      file: null,
      fileActive: false,
      losatFilename: 'alpha.gbk',
      losatFilenameActive: true
    });
    modeDrawing(state, 'linear').files.linearCanonicalComparisons = [];
    state.activeDrawing().losatProgram.value = 'blastn';
    state.activeDrawing().losat.blastn.task = 'megablast';
    state.losatExecution.executionMode = 'serial';
    state.losatCache.value = new Map();

    const comparisonPlanSnapshot = resolveLinearComparisonPlan({
      plan: state.activeDrawing().linearComparisonPlan,
      sequences: state.activeDrawing().linearSeqs,
      layout: [],
      losatProgram: state.activeDrawing().losatProgram.value,
      blastpMode: state.activeDrawing().losat.blastp.mode
    });
    state.activeDrawing().linearSeqs.find(({uid})=>uid==='middle').region_reverse=true;
    recordDisplayRows.value = [{ scope: 'linear', drawingId: state.activeDrawing().id, sourceUid: 'middle',
      key: JSON.stringify(['middle', '#1']), selector: '#1', recordId: 'MIDDLE', recordLength: 8,
      reverse: true, cropped: false, detectedTopology: 'linear' }];
    const linearResult = result('lazy-linear.svg', 'lazy-linear');
    workerResponses.push(response(linearResult, validCatalog(linearResult.name)));
    assert.deepEqual(
      await runner.runAnalysis(comparisonPlanSnapshot),
      { status: 'ok' },
      JSON.stringify({ error: state.errorLog.value, lastWorker: workerMessages.at(-1) })
    );
    assert.equal(losatCalls, 1);
    assert.ok(capturedSequences.includes('>MIDDLE\nCCCCGGGG\n'));
    assert.ok(capturedSequences.includes('>THIRD\nTTTTAAAA\n'));
    // D9/D10: each searched raw entry records the Web runtime identity, which
    // is not part of its key; Run Info lists the runtime of the displayed rows.
    const searchedEntries = Array.from(state.losatCache.value.entries());
    assert.ok(searchedEntries.length > 0);
    for (const [key, entry] of searchedEntries) {
      assert.deepEqual(entry.runtime, { kind: 'losat', source: 'wasm', version: LOSAT_WASM_BUILD.version, program: 'blastn' });
      const { runtime: _runtime, ...withoutRuntime } = entry;
      assert.equal(key, await sha256Text(JSON.stringify(buildLosatCachePayload(withoutRuntime))));
      assert.equal(key, await sha256Text(JSON.stringify(buildLosatCachePayload(entry))));
    }
    assert.deepEqual(state.lastRunInfo.value.losatRuntimes, [{
      program: 'blastn',
      kind: 'losat',
      version: LOSAT_WASM_BUILD.version,
      source: 'wasm',
      path: null,
      text: `blastn: LOSAT ${LOSAT_WASM_BUILD.version} (wasm)`,
      textParts: [{ text: `blastn: LOSAT ${LOSAT_WASM_BUILD.version} (wasm)`, verbatim: true }]
    }]);
    // Generated LOSATN rows stay in the search frame; the Python planner
    // projects the reversed record (PD-OI-073), so no Web conversion runs.
    assert.equal(workerMessages.some(({ type, operation }) => (
      type === 'helper' && /nucleotide/i.test(String(operation || ''))
    )), false);
    const projectedRun = workerMessages.findLast(({ type }) => type === 'run').payload.request;
    assert.equal(projectedRun.records.find(({ recordKey }) => recordKey === 'middle')
      .presentation.reverseComplement, true);
    assert.equal(state.activeDrawing().linearSeqs.find(({ uid }) => uid === 'middle').region_reverse, true);
    recordDisplayRows.value = [];
    assert.equal(
      resourceMetrics.filter(({ name }) => name === 'resourceByteReadCount').length,
      3
    );
    assert.equal(
      resourceMetrics.filter(({ name }) => name === 'resourceTextReadCount').length,
      0,
      'the explicit staged text fast path must not materialize file text again'
    );
    assert.equal(
      state.lastRunInfo.value.sourceRecipe.available,
      true,
      state.lastRunInfo.value.sourceRecipe.unavailableReason
    );
    assert.match(
      state.lastRunInfo.value.sourceRecipe.command,
      /--gbk alpha\.gbk\.tsv middle\.gbk third\.gb/
    );
    assert.match(
      state.lastRunInfo.value.sourceRecipe.command,
      /--comparisons_table comparisons\.tsv/
    );
    assert.doesNotMatch(state.lastRunInfo.value.sourceRecipe.command, /\/blast_0\.txt/);
    assert.equal(state.lastRunInfo.value.reproducibility.level, 'session-recommended');
    assert.match(
      state.lastRunInfo.value.reproducibility.notes.join('\n'),
      /Download reproducibility files/
    );
    assert.match(state.lastRunInfo.value.exactReplay.command, /--session lazy-linear\.gbdraw$/);
    await runner.downloadCliHelperFiles();
    assert.equal(downloadedFilename, 'lazy-linear-cli-files.zip');
    const linearHelperEntries = await storedZipEntries(downloadedBlob);
    assertDownloadedFilesMatchCommands(state.lastRunInfo.value, linearHelperEntries);
    assert.ok(linearHelperEntries.some(({ name }) => name.endsWith('.gbdraw')));
    const blastEntry = linearHelperEntries.find(({ name }) => name === 'alpha.gbk-2.tsv');
    const comparisonsEntry = linearHelperEntries.find(({ name }) => name === 'comparisons.tsv');
    assert.ok(blastEntry?.text.includes('MIDDLE\tTHIRD'));
    assert.ok(comparisonsEntry);
    const [comparisonHeader, comparisonRow] = comparisonsEntry.text.trim()
      .split('\n')
      .map((line) => line.split('\t'));
    assert.equal(comparisonRow[comparisonHeader.indexOf('blast')], 'alpha.gbk-2.tsv');
    assert.doesNotMatch(comparisonsEntry.text, /^alpha\.gbk\.tsv\t/m);
    assert.deepEqual(
      state.lastRunInfo.value.sourceRecipe.helperFiles.map(({ name }) => name).sort(),
      ['alpha.gbk-2.tsv', 'comparisons.tsv']
    );
    assert.doesNotMatch(
      [
        ...state.lastRunInfo.value.sourceRecipe.commandArgs,
        ...state.lastRunInfo.value.exactReplay.commandArgs,
        ...linearHelperEntries.map(({ name }) => name)
      ].join('\n'),
      /comparison-nucleotide-|\/blast_0\.txt/
    );

    const committedLosatCacheEntries = Array.from(state.losatCache.value.entries());
    const committedLosatCacheInfo = structuredClone(state.losatCacheInfo.value);
    const committedLosatTelemetry = structuredClone(
      globalThis.__GBDRAW_LAST_LOSAT_TELEMETRY__
    );
    state.activeDrawing().losat.blastn.task = 'blastn';
    failLateArtifactAdoption = true;
    resourceMetrics.splice(0);
    workerHelperResponses.push({
      ok: true,
      result: {
        tsv: 'MIDDLE\tTHIRD\t100\t8\t0\t0\t1\t8\t1\t8\t1e-10\t20\n'
      }
    });
    const failedLinearResult = result('lazy-linear-failed.svg', 'lazy-linear-failed');
    workerResponses.push(response(
      failedLinearResult,
      validCatalog(failedLinearResult.name)
    ));
    assert.equal((await runner.runAnalysis(comparisonPlanSnapshot)).status, 'error');
    // OV-1589: the LOSAT FASTA comes from the records parsed once per File and the
    // Worker holds the sources, so no input file is read again.
    assert.deepEqual(
      resourceMetrics.filter(({ name }) => name === 'resourceByteReadCount'),
      []
    );
    assert.equal(state.activeDrawing().linearSeqs.find(({ uid }) => uid === 'middle').region_reverse, true);
    failLateArtifactAdoption = false;
    state.activeDrawing().losat.blastn.task = 'megablast';
    assert.equal(state.errorLog.value?.code, 'UNKNOWN');
    assert.doesNotMatch(JSON.stringify(state.errorLog.value), /injected LOSAT late artifact adoption failure/);
    // D-200: the rollback puts back the committed rows and keeps every search, the failed run's completed one too
    // (PD-OI-022 retry): no install removes a cached search.
    assert.deepEqual(
      Array.from(state.losatCache.value.entries()).slice(0, committedLosatCacheEntries.length),
      committedLosatCacheEntries,
      'artifact rollback keeps the committed searches'
    );
    assert.equal(state.losatCache.value.size, committedLosatCacheEntries.length + 1, 'and the failed run\'s completed search');
    assert.deepEqual(state.losatCacheInfo.value, committedLosatCacheInfo);
    assert.deepEqual(
      globalThis.__GBDRAW_LAST_LOSAT_TELEMETRY__,
      committedLosatTelemetry,
      'failed candidate provenance must roll back with the committed artifact'
    );
    assert.equal(losatCalls, 2, 'the failed candidate must execute without becoming committed evidence');

    state.activeDrawing().losat.blastn.task = 'blastn';
    workerHelperResponses.push({
      ok: true,
      result: { tsv: 'MIDDLE\tTHIRD\t100\t8\t0\t0\t1\t8\t1\t8\t1e-10\t20\n' }
    });
    workerResponses.push(response(linearResult, validCatalog(linearResult.name)));
    assert.deepEqual(await runner.runAnalysis(comparisonPlanSnapshot), { status: 'ok' });
    assert.equal(losatCalls, 2, 'retry after artifact rollback must reuse completed raw searches');

    const warmCacheInfo = [{
      key: 'warm-cache', edgeKey: 'multi->third', queryUid: 'multi', subjectUid: 'third',
      queryIndex: 0, subjectIndex: 1, ordinal: 0, display: true
    }];
    state.losatCacheInfo.value = structuredClone(warmCacheInfo);
    state.activeDrawing().losatProgram.value = 'blastp';
    state.activeDrawing().losat.blastp.mode = 'pairwise';
    state.activeDrawing().losat.blastp.maxHits = 5;
    state.activeDrawing().losat.blastp.candidateLimit = null;
    const warmSettings = {
      proteinBlastpMaxHits: 5,
      proteinBlastpCandidateLimit: null
    };
    modeDrawing(state, 'linear').files.linearCanonicalComparisons = [{
      kind: 'precomputedProteinComparison',
      edgeKey: 'multi->third',
      ordinal: 0,
      queryRecordIndex: 0,
      subjectRecordIndex: 1,
      encoding: 'canonicalTsv',
      file: new AuditFile([
        'query\tsubject\tidentity\talignment_length\tmismatches\tgap_opens\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\n'
      ], 'warm.tsv')
    }, {
      kind: 'generatedProteinComparison', mode: 'none', pairs: [], settings: warmSettings
    }];
    committedRenderRequest = {
      records: structuredClone(committedRenderRequest.records),
      diagramOptions: {
        bitscore: Number(state.activeDrawing().adv.min_bitscore),
        evalue: Number(state.activeDrawing().adv.evalue),
        identity: Number(state.activeDrawing().adv.identity),
        alignmentLength: Number(state.activeDrawing().adv.alignment_length)
      },
      comparisons: [{
        kind: 'precomputedProteinComparison'
      }, {
        kind: 'generatedProteinComparison', mode: 'none', pairs: [], settings: warmSettings
      }]
    };
    const warmComparisonPlan = resolveLinearComparisonPlan({
      plan: state.activeDrawing().linearComparisonPlan,
      sequences: state.activeDrawing().linearSeqs,
      layout: [],
      losatProgram: state.activeDrawing().losatProgram.value,
      blastpMode: state.activeDrawing().losat.blastp.mode
    });
    const warmResult = result('lazy-linear-warm.svg', 'lazy-linear-warm');
    workerResponses.push(response(warmResult, validCatalog(warmResult.name)));
    assert.deepEqual(
      await runner.runAnalysis(warmComparisonPlan),
      { status: 'ok' },
      JSON.stringify(state.errorLog.value)
    );
    assert.equal(losatCalls, 2, 'resolved protein artifacts must bypass further LOSAT execution');
    assert.deepEqual(state.losatCacheInfo.value, warmCacheInfo);

    // OV-1128 B (D-114): a Similarity groups Generate draws, records and keys with the stored
    // paralog links per group. The committed result was drawn with 3, so it is reused and the
    // request records 3.
    const pairwiseCommittedRequest = committedRenderRequest;
    const pairwiseCanonicalComparisons = modeDrawing(state, 'linear').files.linearCanonicalComparisons;
    const groupSettings = {
      proteinBlastpCandidateLimit: null,
      orthogroupMembershipMode: state.activeDrawing().losat.blastp.orthogroupMembershipMode,
      orthogroupMemberMaxHits: null,
      collinearMaxParalogLinksPerOrthogroup: 3
    };
    // Similarity groups compares every adjacent pair.
    state.activeDrawing().linearComparisonPlan.mode = 'adjacent';
    state.activeDrawing().losat.blastp.mode = 'orthogroup';
    state.activeDrawing().losat.blastp.orthogroupMemberMaxHits = null;
    state.activeDrawing().losat.blastp.collinearMaxParalogLinksPerOrthogroup = 3;
    modeDrawing(state, 'linear').files.linearCanonicalComparisons = [{
      kind: 'orthogroupResult',
      encoding: 'canonicalJson',
      file: new AuditFile(['{"schema":3,"kind":"orthogroupResult","value":{}}\n'], 'groups.json')
    }, {
      kind: 'generatedProteinComparison', mode: 'none', pairs: [], settings: groupSettings
    }];
    committedRenderRequest = {
      ...pairwiseCommittedRequest,
      comparisons: [{ kind: 'orthogroupResult' }, {
        kind: 'generatedProteinComparison', mode: 'none', pairs: [], settings: groupSettings
      }]
    };
    const groupResult = result('lazy-linear-groups.svg', 'lazy-linear-groups');
    workerResponses.push(response(groupResult, validCatalog(groupResult.name)));
    assert.deepEqual(
      await runner.runAnalysis(resolveLinearComparisonPlan({
        plan: state.activeDrawing().linearComparisonPlan,
        sequences: state.activeDrawing().linearSeqs,
        layout: [],
        losatProgram: state.activeDrawing().losatProgram.value,
        blastpMode: state.activeDrawing().losat.blastp.mode
      })),
      { status: 'ok' },
      JSON.stringify(state.errorLog.value)
    );
    assert.equal(losatCalls, 2, 'the same paralog links reuse the committed groups');
    assert.equal(
      workerMessages.filter(({ type }) => type === 'run').at(-1).payload.request.comparisons
        .find(({ kind }) => kind === 'generatedProteinComparison').settings.collinearMaxParalogLinksPerOrthogroup,
      3,
      'a Similarity groups Generate records the stored paralog links per group'
    );
    state.activeDrawing().linearComparisonPlan.mode = 'selected';
    state.activeDrawing().losat.blastp.mode = 'pairwise';
    state.activeDrawing().losat.blastp.collinearMaxParalogLinksPerOrthogroup = 2;
    modeDrawing(state, 'linear').files.linearCanonicalComparisons = pairwiseCanonicalComparisons;
    committedRenderRequest = pairwiseCommittedRequest;
    committedRenderRequest = {
      ...committedRenderRequest,
      comparisons: [{ kind: 'precomputedProteinComparison' }, {
        kind: 'generatedProteinComparison', mode: 'none', pairs: [], settings: warmSettings
      }]
    };

    const alignmentAnchor = (recordKey, featureId, sourceFeatureIndex) => ({
      recordKey,
      biologicalFeatureId: featureId,
      sourceFeatureIndex,
      stableFeatureSvgId: `stable-${featureId}`
    });
    const alignmentReference = alignmentAnchor('multi', 'reference-feature', 0);
    const alignmentTarget = alignmentAnchor('middle', 'target-feature', 1);
    const alignmentPlan = {
      schema: 2,
      groupId: 'og-web-action',
      reference: alignmentReference,
      records: [
        {
          recordKey: 'multi', status: 'reference', rationale: 'reference',
          anchor: alignmentReference
        },
        {
          recordKey: 'middle', status: 'aligned', rationale: 'only_usable_candidate',
          anchor: alignmentTarget
        },
        {
          recordKey: 'third', status: 'skipped', rationale: 'skipped_no_candidate',
          anchor: null
        }
      ]
    };
    const alignmentTranslations = ['multi', 'middle', 'third'].map((recordKey) => ({
      recordKey, x: 0, y: 0
    }));
    const alignmentOrientations = ['multi', 'middle', 'third'].map((recordKey) => ({
      recordKey,
      reverseComplement: recordKey !== 'third'
    }));
    const helperRequestsBeforeAlignment = workerMessages.filter(({ type }) => type === 'helper')
      .length;
    const alignedResult = result('lazy-linear-aligned.svg', 'lazy-linear-aligned');
    workerResponses.push(response(alignedResult, validCatalog(alignedResult.name)));
    assert.deepEqual(
      await runner.runCommittedCanonicalCandidate({
        canonical:runner.projectCommittedSimilarityAlignment({committed:committedCanonicalSession,plan:alignmentPlan,translations:alignmentTranslations,orientations:alignmentOrientations,starts:[]}),
        label:'Align Similarity Group',
        commitIntent:()=>alignmentOrientations.forEach(({recordKey,reverseComplement})=>{state.activeDrawing().linearSeqs.find(s=>s.uid===recordKey).region_reverse=reverseComplement;})
      }),
      { status: 'ok' },
      JSON.stringify(state.errorLog.value)
    );
    assert.equal(losatCalls, 2, 'alignment Apply must reuse committed protein evidence');
    assert.deepEqual(state.similarityAlignmentPlan.value, alignmentPlan);
    assert.deepEqual(state.linearRecordTranslations.value, alignmentTranslations);
    const alignedRunRequest = workerMessages.filter(({ type }) => type === 'run').at(-1)
      .payload.request;
    assert.deepEqual(alignedRunRequest.layout.similarityAlignment, alignmentPlan);
    assert.deepEqual(alignedRunRequest.layout.recordTranslations, alignmentTranslations);
    const alignedOrientation = Object.fromEntries(alignedRunRequest.records.map((record) => [
      record.recordKey,
      Boolean(record.region?.reverseComplement ?? record.presentation?.reverseComplement)
    ]));
    assert.equal(alignedOrientation.middle, true);
    assert.equal(
      state.activeDrawing().linearSeqs.find(({ uid }) => uid === 'middle').region_reverse,
      true,
      'successful admission must install the materialized base orientation'
    );
    assert.equal(
      workerMessages.filter(({ type }) => type === 'helper').length,
      helperRequestsBeforeAlignment,
      'alignment Apply must not schedule protein extraction or group inference'
    );

    const committedAlignedResult = state.results.value;
    const committedAlignmentPlan = state.similarityAlignmentPlan.value;
    committedRenderRequest = {
      ...committedRenderRequest,
      comparisons: [{ kind: 'precomputedProteinComparison' }, {
        kind: 'generatedProteinComparison', mode: 'none', pairs: [], settings: warmSettings
      }]
    };
    failLateArtifactAdoption = true;
    const rejectedAlignedResult = result(
      'lazy-linear-aligned-rejected.svg',
      'lazy-linear-aligned-rejected'
    );
    workerResponses.push(response(
      rejectedAlignedResult,
      validCatalog(rejectedAlignedResult.name)
    ));
    assert.deepEqual(
      (await runner.runCommittedCanonicalCandidate({canonical:runner.projectCommittedSimilarityAlignment({committed:committedCanonicalSession,plan:{...alignmentPlan,groupId:'og-rejected-candidate'},translations:alignmentTranslations,orientations:alignmentOrientations,starts:[]})})).status,
      'error'
    );
    failLateArtifactAdoption = false;
    assert.equal(state.errorLog.value?.code, 'UNKNOWN');
    assert.doesNotMatch(JSON.stringify(state.errorLog.value), /injected LOSAT late artifact adoption failure/);
    assert.equal(state.results.value, committedAlignedResult);
    assert.equal(state.similarityAlignmentPlan.value, committedAlignmentPlan);
    assert.equal(losatCalls, 2, 'failed alignment admission must not schedule LOSATP');

    // OV-1556, OV-1565: one GenBank file of 13 records drawn as 13 cards, and a 14th card that draws the first
    // record again. Generate parses the file once and slices each card's record from it, and the next Generate
    // parses nothing and digests no FASTA again.
    const recordIds = Array.from({ length: 13 }, (_, index) => `REC${index + 1}`);
    const bases = ['a', 'c', 'g', 't'];
    const thirteenFile = new AuditFile([recordIds.map((id, index) => [
      `LOCUS       ${id} 8 bp DNA`,
      'ORIGIN',
      `        1 ${bases[index % 4].repeat(4)}${bases[(index + 1) % 4].repeat(4)}`,
      '//'
    ].join('\n')).join('\n') + '\n'], 'thirteen.gb');
    state.activeDrawing().linearSeqs.splice(0, state.activeDrawing().linearSeqs.length, ...[...recordIds, 'REC1'].map((id, index) => ({
      uid: `card-${index + 1}`, gb: bound(thirteenFile), gff: null, fasta: null, depth: null, blast: null,
      losat_gencode: 1, losat_filename: '', definition: '', record_subtitle: '', region_record_id: id,
      region_start: null, region_end: null, region_reverse: false,
      windows: []
    })));
    modeDrawing(state, 'linear').files.linearCanonicalComparisons = [];
    state.similarityAlignmentPlan.value = null;
    state.linearRecordTranslations.value = [];
    state.losatCacheInfo.value = [];
    state.losatCache.value = new Map();
    state.activeDrawing().losatProgram.value = 'blastn';
    state.activeDrawing().losat.blastn.task = 'megablast';
    Object.assign(state.activeDrawing().linearComparisonPlan, { mode: 'adjacent', defaultSource: 'losat' });
    state.activeDrawing().linearComparisonPlan.edges.splice(0);
    const thirteenPlan = resolveLinearComparisonPlan({
      plan: state.activeDrawing().linearComparisonPlan,
      sequences: state.activeDrawing().linearSeqs,
      layout: [],
      losatProgram: 'blastn',
      blastpMode: state.activeDrawing().losat.blastp.mode
    });
    globalThis.__GBDRAW_LOSAT_EXECUTOR__ = async (jobs) => {
      losatCalls += 1;
      return jobs.map((job) => ({ cacheKey: job.cacheKey, text: '' }));
    };
    const thirteenWatches = [];
    for (const name of ['thirteen-first', 'thirteen-second']) {
      const thirteenResult = result(`${name}.svg`, name);
      workerResponses.push(response(thirteenResult, validCatalog(thirteenResult.name)));
      thirteenWatches.push(await watchLosatPreparation(async () => assert.deepEqual(
        await runner.runAnalysis(thirteenPlan),
        { status: 'ok' },
        JSON.stringify(state.errorLog.value)
      )));
    }
    assert.equal(losatCalls, 3, 'the first Generate searched and the second reused its raw results');
    // OV-3227: the File name field's default is the name Generate gave each pair, the cards' record IDs.
    const thirteenNames = state.losatCacheInfo.value.map(({ filename }) => filename);
    assert.deepEqual(thirteenNames.slice(0, 2), ['REC1.REC2.losatn.tsv', 'REC2.REC3.losatn.tsv']);
    assert.deepEqual(state.losatCacheInfo.value.map(({ edgeKey }) => runner.getLosatPairDefaultName(state.activeDrawing(), edgeKey)),
      thirteenNames);
    assert.match(thirteenWatches[0].timing, /FASTA cache hits=1, .*JS FASTA=1, Worker FASTA=0/);
    assert.match(thirteenWatches[1].timing, /FASTA cache hits=14, .*JS FASTA=0, Worker FASTA=0/);
    assert.equal(thirteenWatches[1].fastaDigests, 0);
    // OV-1594: the held LOSAT FASTA of a record is frozen, and its SHA-256 is its own read beside it, so the
    // Generates read the entries without writing to them. A record's FASTA is digested once, whatever its card
    // count (the 14th card reads REC1's digest); its single-record query side reads that same digest (OV-1598).
    const heldLosatFasta = (id) => readOncePerFile([thirteenFile], `losat-fasta:${JSON.stringify({
      fmt: 'genbank', regionSpec: null, recordSelector: id, reverseFlag: '0'
    })}`, () => Promise.reject(new Error(`no held LOSAT FASTA for ${id}`)));
    for (const id of recordIds) {
      const entry = await heldLosatFasta(id);
      assert.ok(Object.isFrozen(entry), id);
      assert.deepEqual(Object.keys(entry).sort(), ['canonicalLength', 'fasta', 'recordId'], id);
      assert.equal(thirteenWatches[0].digestsByText.get(entry.fasta), 1, id);
    }

    // OV-1555, OV-1590: after a Generate, the blastp extraction of two of its records sends no
    // GenBank bytes. It stages the record source under the resource id and token the Worker
    // holds, so Python reads the records that Generate parsed; and it leaves the staging cache
    // as it was, so the next Generate sends nothing either.
    state.activeDrawing().linearSeqs.splice(2);
    const twoPlan = resolveLinearComparisonPlan({
      plan: state.activeDrawing().linearComparisonPlan, sequences: state.activeDrawing().linearSeqs, layout: [],
      losatProgram: 'blastn', blastpMode: state.activeDrawing().losat.blastp.mode
    });
    // The plan of the cards drawn now: later steps replace the cards, and a plan names the cards it was made for.
    const twoPlanNow = () => resolveLinearComparisonPlan({
      plan: state.activeDrawing().linearComparisonPlan, sequences: state.activeDrawing().linearSeqs, layout: [],
      losatProgram: 'blastn', blastpMode: state.activeDrawing().losat.blastp.mode
    });
    const twoResult = result('two-cards.svg', 'two-cards');
    workerResponses.push(response(twoResult, validCatalog(twoResult.name)));
    assert.deepEqual(await runner.runAnalysis(twoPlan), { status: 'ok' }, JSON.stringify(state.errorLog.value));
    const drawnRun = workerMessages.findLast(({ type }) => type === 'run').payload;
    const drawnToken = drawnRun.resourceManifest.find(({ resourceId }) => resourceId === 'record-1-genbank').cacheToken;
    state.activeDrawing().losatProgram.value = 'blastp';
    state.activeDrawing().losat.blastp.mode = 'pairwise';
    modeDrawing(state, 'linear').files.linearCanonicalComparisons = [];
    const helpersBeforeProteins = workerMessages.length;
    workerHelperResponses.splice(0);
    ['REC1', 'REC2'].forEach((recordId) => workerHelperResponses.push({ ok: true, result: {
      fasta: `>${recordId}_p1\nMK\n`, record_id: recordId, record_length: 8, protein_count: 1, protein_map: {},
      protein_set_hash: `set-${recordId}`, record_analysis_id: recordId, record_instance_key: recordId,
      runtime_binding_hash: `rt-${recordId}`, display_binding_hash: `dp-${recordId}`, identity_manifest: null
    } }));
    // The Generate stops after the extractions: these replies carry no protein identity manifest.
    assert.equal((await runner.runAnalysis(resolveLinearComparisonPlan({
      plan: state.activeDrawing().linearComparisonPlan, sequences: state.activeDrawing().linearSeqs, layout: [],
      losatProgram: 'blastp', blastpMode: 'pairwise'
    }))).status, 'error');
    const extractions = workerMessages.slice(helpersBeforeProteins)
      .filter(({ type, operation }) => type === 'helper' && operation === 'extractCdsProteinFasta');
    assert.equal(extractions.length, 2);
    for (const { payload } of extractions) {
      assert.deepEqual((payload.files || []).map(({ role }) => role).filter((role) => role !== 'visibility'), []);
      assert.deepEqual(payload.recordSource, { kind: 'genbank', resourceId: 'record-1-genbank' });
      assert.deepEqual(payload.resourceManifest.map(({ resourceId, cacheToken }) => [resourceId, cacheToken]),
        [['record-1-genbank', drawnToken]]);
      assert.deepEqual(payload.stagedResources, []);
    }
    workerHelperResponses.splice(0);
    state.activeDrawing().losatProgram.value = 'blastn';
    const afterProteins = result('two-cards-after.svg', 'two-cards-after');
    workerResponses.push(response(afterProteins, validCatalog(afterProteins.name)));
    assert.deepEqual(await runner.runAnalysis(twoPlan), { status: 'ok' }, JSON.stringify(state.errorLog.value));
    assert.deepEqual(workerMessages.findLast(({ type }) => type === 'run').payload.stagedResources, []);

    // OV-1595: the first blastp Generate after an upload sends the new file once. Its first extraction stages
    // it under the id the request names (record-1-genbank), so the other cards and the Generate that follows
    // send nothing, and Python parses it once for all of them.
    const uploadFile = new AuditFile([['UP1', 'UP2', 'UP3'].map((id) => [
      `LOCUS       ${id} 8 bp DNA`, 'ORIGIN', '        1 acgtacgt', '//'
    ].join('\n')).join('\n') + '\n'], 'upload.gb', { lastModified: 7 });
    state.activeDrawing().linearSeqs.splice(0, state.activeDrawing().linearSeqs.length, ...['UP1', 'UP2', 'UP3'].map((id, index) => ({
      uid: `upload-${index + 1}`, gb: bound(uploadFile), gff: null, fasta: null, depth: null, blast: null,
      losat_gencode: 1, losat_filename: '', definition: '', record_subtitle: '', region_record_id: id,
      region_start: null, region_end: null, region_reverse: false,
      windows: []
    })));
    const uploadPlan = (losatProgram) => resolveLinearComparisonPlan({
      plan: state.activeDrawing().linearComparisonPlan, sequences: state.activeDrawing().linearSeqs, layout: [],
      losatProgram, blastpMode: 'pairwise'
    });
    state.activeDrawing().losatProgram.value = 'blastp';
    const messagesBeforeUpload = workerMessages.length;
    ['UP1', 'UP2', 'UP3'].forEach((recordId) => workerHelperResponses.push({ ok: true, result: {
      fasta: `>${recordId}_p1\nMK\n`, record_id: recordId, record_length: 8, protein_count: 1, protein_map: {},
      protein_set_hash: `set-${recordId}`, record_analysis_id: recordId, record_instance_key: recordId,
      runtime_binding_hash: `rt-${recordId}`, display_binding_hash: `dp-${recordId}`, identity_manifest: null
    } }));
    assert.equal((await runner.runAnalysis(uploadPlan('blastp'))).status, 'error');
    const uploadExtractions = workerMessages.slice(messagesBeforeUpload)
      .filter(({ type, operation }) => type === 'helper' && operation === 'extractCdsProteinFasta');
    assert.equal(uploadExtractions.length, 3);
    const uploadToken = uploadExtractions[0].payload.resourceManifest?.[0]?.cacheToken;
    assert.notEqual(uploadToken, drawnToken);
    assert.deepEqual(uploadExtractions.map(({ payload }) => [
      payload.recordSource,
      (payload.files || []).map(({ role }) => role).filter((role) => role !== 'visibility'),
      (payload.resourceManifest ?? []).map(({ resourceId, cacheToken }) => [resourceId, cacheToken]),
      (payload.stagedResources ?? []).map(({ resourceId }) => resourceId)
    ]), [0, 1, 2].map((index) => [
      { kind: 'genbank', resourceId: 'record-1-genbank' },
      [],
      [['record-1-genbank', uploadToken]],
      index === 0 ? ['record-1-genbank'] : []
    ]));
    workerHelperResponses.splice(0);
    state.activeDrawing().losatProgram.value = 'blastn';
    const uploadResult = result('upload.svg', 'upload');
    workerResponses.push(response(uploadResult, validCatalog(uploadResult.name)));
    assert.deepEqual(await runner.runAnalysis(uploadPlan('blastn')), { status: 'ok' }, JSON.stringify(state.errorLog.value));
    const uploadRun = workerMessages.findLast(({ type }) => type === 'run').payload;
    assert.equal(
      uploadRun.resourceManifest.find(({ resourceId }) => resourceId === 'record-1-genbank').cacheToken,
      uploadToken
    );
    // Only the comparison table of the new third pair is new to the Worker.
    assert.deepEqual(uploadRun.stagedResources.map(({ resourceId }) => resourceId)
      .filter((resourceId) => resourceId.startsWith('record-')), []);
    const uploadBytesSent = workerMessages.slice(messagesBeforeUpload)
      .flatMap(({ payload }) => payload?.stagedResources ?? [])
      .filter(({ resourceId }) => resourceId.startsWith('record-'))
      .reduce((total, { blob }) => total + blob.size, 0);
    assert.equal(uploadBytesSent, uploadFile.size);

    // A new file inserted before the committed one: the next request names it record-1 and the committed file
    // record-2. The new file is staged once; the committed file's bytes stay held under their content token, so
    // its new id sends nothing, and the render sends neither again.
    const insertedFile = new AuditFile([
      'LOCUS       NEW1 8 bp DNA\nORIGIN\n        1 ggggcccc\n//\n'
    ], 'inserted.gb', { lastModified: 9 });
    state.activeDrawing().linearSeqs.unshift({
      uid: 'inserted-1', gb: bound(insertedFile), gff: null, fasta: null, depth: null, blast: null,
      losat_gencode: 1, losat_filename: '', definition: '', record_subtitle: '', region_record_id: 'NEW1',
      region_start: null, region_end: null, region_reverse: false,
      windows: []
    });
    state.activeDrawing().losatProgram.value = 'blastp';
    const messagesBeforeInsert = workerMessages.length;
    ['NEW1', 'UP1', 'UP2', 'UP3'].forEach((recordId) => workerHelperResponses.push({ ok: true, result: {
      fasta: `>${recordId}_p1\nMK\n`, record_id: recordId, record_length: 8, protein_count: 1, protein_map: {},
      protein_set_hash: `set-i-${recordId}`, record_analysis_id: `i-${recordId}`, record_instance_key: `i-${recordId}`,
      runtime_binding_hash: `rt-i-${recordId}`, display_binding_hash: `dp-i-${recordId}`, identity_manifest: null
    } }));
    assert.equal((await runner.runAnalysis(uploadPlan('blastp'))).status, 'error');
    assert.deepEqual(workerMessages.slice(messagesBeforeInsert)
      .filter(({ type, operation }) => type === 'helper' && operation === 'extractCdsProteinFasta')
      .map(({ payload }) => [payload.recordSource.resourceId, payload.stagedResources.map(({ resourceId }) => resourceId)]), [
      ['record-1-genbank', ['record-1-genbank']],
      ['record-2-genbank', []],
      ['record-2-genbank', []],
      ['record-2-genbank', []]
    ]);
    workerHelperResponses.splice(0);
    state.activeDrawing().losatProgram.value = 'blastn';
    const insertedResult = result('inserted.svg', 'inserted');
    workerResponses.push(response(insertedResult, validCatalog(insertedResult.name)));
    assert.deepEqual(await runner.runAnalysis(uploadPlan('blastn')), { status: 'ok' }, JSON.stringify(state.errorLog.value));
    const recordBytesSent = workerMessages.slice(messagesBeforeInsert)
      .flatMap(({ payload }) => payload?.stagedResources ?? [])
      .filter(({ resourceId }) => resourceId.startsWith('record-'))
      .reduce((total, { blob }) => total + blob.size, 0);
    assert.equal(recordBytesSent, insertedFile.size);
    state.activeDrawing().linearSeqs.splice(2);

    // A file that starts with a UTF-8 byte order mark: the record list leaves the file to the Python reader,
    // and the LOSAT FASTA does too.
    const bomFile = new AuditFile([new Uint8Array([0xef, 0xbb, 0xbf]), ['BOM1', 'BOM2'].map((id) => [
      `LOCUS       ${id} 8 bp DNA`, 'ORIGIN', '        1 acgtacgt', '//'
    ].join('\n')).join('\n') + '\n'], 'bom.gb');
    state.activeDrawing().linearSeqs.forEach((seq, index) => Object.assign(seq, { gb: bound(bomFile), region_record_id: `BOM${index + 1}` }));
    ['BOM1', 'BOM2'].forEach((id) => workerHelperResponses.push({ ok: true, result: {
      fasta: `>${id}\nACGTACGT\n`, record_id: id, record_length: 8
    } }));
    const helpersBeforeBom = workerMessages.length;
    const bomResult = result('bom.svg', 'bom');
    workerResponses.push(response(bomResult, validCatalog(bomResult.name)));
    const bomWatch = await watchLosatPreparation(async () => assert.deepEqual(
      await runner.runAnalysis(twoPlanNow()), { status: 'ok' }, JSON.stringify(state.errorLog.value)
    ));
    assert.match(bomWatch.timing, /JS FASTA=0, Worker FASTA=2/);
    assert.deepEqual(workerMessages.slice(helpersBeforeBom)
      .filter(({ type, operation }) => type === 'helper' && operation !== 'planLinearUnits')
      .map(({ operation, payload }) => [operation, payload.recordSelector]),
    [['extractFirstFasta', 'BOM1'], ['extractFirstFasta', 'BOM2']]);

    // The GFF3 and FASTA pair of a drawn record: after a Generate, the blastp extraction names the committed
    // GFF3 resource as the record's GFF3 file and the FASTA resource as its FASTA file, in that order.
    const gffFile = new AuditFile(['##gff-version 3\nG1\t.\tCDS\t1\t9\t.\t+\t0\tID=c1\n'], 'pair.gff');
    const pairedFasta = new AuditFile(['>G1\nATGAAATAA\n'], 'pair.fasta');
    modeDrawing(state, 'linear').inputType.value = 'gff';
    state.activeDrawing().losatProgram.value = 'blastn';
    state.activeDrawing().linearSeqs.forEach((seq) => Object.assign(seq, { gb: null, gff: bound(gffFile), fasta: bound(pairedFasta), region_record_id: 'G1' }));
    const pairResult = result('gff-pair.svg', 'gff-pair');
    workerResponses.push(response(pairResult, validCatalog(pairResult.name)));
    assert.deepEqual(await runner.runAnalysis(twoPlanNow()), { status: 'ok' }, JSON.stringify(state.errorLog.value));
    const pairRun = workerMessages.findLast(({ type }) => type === 'run').payload;
    const pairTokens = Object.fromEntries(pairRun.resourceManifest.map(({ resourceId, cacheToken }) => [resourceId, cacheToken]));
    state.activeDrawing().losatProgram.value = 'blastp';
    const helpersBeforePair = workerMessages.length;
    [0, 1].forEach((index) => workerHelperResponses.push({ ok: true, result: {
      fasta: `>G1_p${index}\nMK\n`, record_id: 'G1', record_length: 9, protein_count: 1, protein_map: {},
      protein_set_hash: `set-${index}`, record_analysis_id: `G1-${index}`, record_instance_key: `card-${index + 1}`,
      runtime_binding_hash: `rt-${index}`, display_binding_hash: `dp-${index}`, identity_manifest: null
    } }));
    assert.equal((await runner.runAnalysis(resolveLinearComparisonPlan({
      plan: state.activeDrawing().linearComparisonPlan, sequences: state.activeDrawing().linearSeqs, layout: [],
      losatProgram: 'blastp', blastpMode: 'pairwise'
    }))).status, 'error');
    const pairExtractions = workerMessages.slice(helpersBeforePair)
      .filter(({ type, operation }) => type === 'helper' && operation === 'extractCdsProteinFasta');
    assert.equal(pairExtractions.length, 2);
    for (const { payload } of pairExtractions) {
      assert.equal(payload.format, 'gff');
      assert.deepEqual((payload.files || []).map(({ role }) => role).filter((role) => role !== 'visibility'), []);
      assert.deepEqual(payload.recordSource,
        { kind: 'gffFasta', gffResourceId: 'record-1-gff3', fastaResourceId: 'record-1-fasta' });
      assert.deepEqual(payload.resourceManifest.map(({ resourceId, kind, cacheToken }) => [resourceId, kind, cacheToken]), [
        ['record-1-gff3', 'gff3', pairTokens['record-1-gff3']],
        ['record-1-fasta', 'fasta', pairTokens['record-1-fasta']]
      ]);
      assert.deepEqual(payload.stagedResources, []);
    }
    workerHelperResponses.splice(0);
    state.activeDrawing().losatProgram.value = 'blastn';

    // OV-1681: in GFF3 input `#N` counts the records the GFF3 names (Python's load_gff_fasta), not every FASTA
    // entry, so the LOSAT FASTA of `#1` and `#2` is the drawn G1 and G2, not the entry X the GFF3 does not name.
    const namedGff = new AuditFile([
      '##gff-version 3\nG1\t.\tCDS\t1\t9\t.\t+\t0\tID=c1\nG2\t.\tCDS\t1\t9\t.\t+\t0\tID=c2\n'
    ], 'named.gff');
    const extraFasta = new AuditFile(['>X\nCCCCCCCCC\n>G1\nATGAAATAA\n>G2\nATGGGGTAA\n'], 'extra.fasta');
    state.activeDrawing().linearSeqs.forEach((seq, index) => Object.assign(seq, {
      gff: bound(namedGff), fasta: bound(extraFasta), region_record_id: `#${index + 1}`
    }));
    let searched = [];
    globalThis.__GBDRAW_LOSAT_EXECUTOR__ = async (jobs, options) => {
      searched = [...options.sequences.values()];
      return jobs.map((job) => ({ cacheKey: job.cacheKey, text: '' }));
    };
    const namedResult = result('gff-named.svg', 'gff-named');
    workerResponses.push(response(namedResult, validCatalog(namedResult.name)));
    assert.deepEqual(await runner.runAnalysis(twoPlanNow()), { status: 'ok' }, JSON.stringify(state.errorLog.value));
    assert.deepEqual([...new Set(searched.map((fasta) => fasta.split('\n')[0]))].sort(), ['>G1', '>G2']);
    // The automatic record (an empty selector) is the first record the GFF3 names too: a FASTA file whose first
    // entry X has no GFF3 feature draws G1, so LOSAT searches G1 (OV-1681).
    const g1Gff = new AuditFile(['##gff-version 3\nG1\t.\tCDS\t1\t9\t.\t+\t0\tID=c1\n'], 'g1-only.gff');
    const xFirstFasta = new AuditFile(['>X\nCCCCCCCCC\n>G1\nATGAAATAA\n'], 'x-first.fasta');
    Object.assign(state.activeDrawing().linearSeqs[0], { gff: bound(g1Gff), fasta: bound(xFirstFasta), region_record_id: '' });
    const autoResult = result('gff-auto.svg', 'gff-auto');
    workerResponses.push(response(autoResult, validCatalog(autoResult.name)));
    assert.deepEqual(await runner.runAnalysis(twoPlanNow()), { status: 'ok' }, JSON.stringify(state.errorLog.value));
    assert.deepEqual([...new Set(searched.map((fasta) => fasta.split('\n')[0]))].sort(), ['>G1', '>G2']);
    modeDrawing(state, 'linear').inputType.value = 'gb';

    // OV-3227: a card with no record selector over a file of several records (a loaded Session's card) is named by
    // its number among the drawn cards, by Generate and the File name field alike, also when the card before it is
    // drawn in two windows and so searches as two units; the first record its search reads does not name it. A pair
    // row saved with no name downloads under that name once the cards' records are read.
    const pairsFile = new AuditFile([['P1', 'P2', 'P3'].map((id) => [
      `LOCUS       ${id} 8 bp DNA`, 'ORIGIN', '        1 acgtacgt', '//'
    ].join('\n')).join('\n') + '\n'], 'pairs.gb');
    const pairCard = (uid, recordId, windows = []) => ({
      uid, gb: bound(pairsFile), gff: null, fasta: null, depth: null, blast: null, losat_gencode: 1, losat_filename: '',
      definition: '', record_subtitle: '', region_record_id: recordId, region_start: null, region_end: null,
      region_reverse: false, windows
    });
    const span = (start, end) => ({ start, end, reverseComplement: false, column: null });
    state.activeDrawing().linearSeqs.splice(0, state.activeDrawing().linearSeqs.length,
      pairCard('pair-1', 'P1', [span(1, 4), span(5, 8)]), pairCard('pair-2', ''), pairCard('pair-3', 'P3'));
    wholeFileRecords.set('pair-2', ['P1', 'P2', 'P3'].map((recordId) => ({ recordId, recordLength: 8 })));
    state.activeDrawing().losatProgram.value = 'blastn';
    state.losatCacheInfo.value = [];
    state.losatCache.value = new Map();
    const pairsResult = result('whole-file-card.svg', 'whole-file-card');
    workerResponses.push(response(pairsResult, validCatalog(pairsResult.name)));
    assert.deepEqual(await runner.runAnalysis(twoPlanNow()), { status: 'ok' }, JSON.stringify(state.errorLog.value));
    const pairRows = state.losatCacheInfo.value;
    assert.deepEqual(pairRows.map(({ filename }) => filename).sort(),
      ['P1.seq_2.losatn.w1.tsv', 'P1.seq_2.losatn.w2.tsv', 'seq_2.P3.losatn.tsv']);
    const wholeFilePair = linearComparisonEdgeKey('pair-2', 'pair-3');
    assert.equal(runner.getLosatPairDefaultName(state.activeDrawing(), wholeFilePair), 'seq_2.P3.losatn.tsv');
    state.losatCacheInfo.value = pairRows.map((row) => (row.edgeKey === wholeFilePair ? { ...row, filename: '' } : row));
    unreadCards.add('pair-3');
    assert.equal(runner.getLosatPairDefaultName(state.activeDrawing(), wholeFilePair), '', 'no default before the read');
    downloadedFilename = '';
    await runner.downloadLosatPair(state.activeDrawing(), wholeFilePair, '');
    assert.equal(downloadedFilename, 'seq_2.P3.losatn.tsv');
    // #1911 review N1: a window pair whose card pair has a typed name downloads under that name without reading every
    // card, though neither the window pair nor its row has a name of its own.
    const windowRow = pairRows.find(({ filename }) => filename === 'P1.seq_2.losatn.w1.tsv');
    state.activeDrawing().linearComparisonPlan.edges.push({
      id: 'typed-card-pair', queryUid: 'pair-1', subjectUid: 'pair-2', included: true, source: 'losat', file: null,
      fileActive: false, losatFilename: 'typed.tsv', losatFilenameActive: true
    });
    state.losatCacheInfo.value = pairRows.map((row) => (row === windowRow ? { ...row, filename: '' } : row));
    // The drawn Result's units, as Python's catalog names them: the two windows of pair-1, then pair-2 and pair-3.
    const shownCatalog = state.featureCatalog.value;
    state.featureCatalog.value = { ...shownCatalog, items: [{ ...shownCatalog.items[0],
      recordKeys: ['pair-1~w1', 'pair-1~w2', 'pair-2', 'pair-3'], windowOf: { 'pair-1~w1': 'pair-1', 'pair-1~w2': 'pair-1' } }] };
    unreadCards.add('pair-3');
    downloadedFilename = '';
    await runner.downloadLosatPair(state.activeDrawing(), windowRow.edgeKey, '');
    assert.ok(unreadCards.has('pair-3'), 'a typed card pair name needs no read of the cards');
    assert.equal(downloadedFilename, runner.getLosatPairDefaultName(state.activeDrawing(), windowRow.edgeKey));
    assert.match(downloadedFilename, /^typed/);
    state.featureCatalog.value = shownCatalog;

    // OV-3326: a TLOSATX Gencode that is not an NCBI table LOSAT reads stops Generate before any search, and the
    // Generation Error panel names the card, its field and the tables.
    let searches = 0;
    globalThis.__GBDRAW_LOSAT_EXECUTOR__ = async (jobs) => {
      searches += jobs.length;
      return jobs.map((job) => ({ cacheKey: job.cacheKey, text: '' }));
    };
    state.activeDrawing().losatProgram.value = 'tblastx';
    state.losatCache.value = new Map();
    state.activeDrawing().linearSeqs[2].losat_gencode = 7;
    assert.equal((await runner.runAnalysis(resolveLinearComparisonPlan({
      plan: state.activeDrawing().linearComparisonPlan, sequences: state.activeDrawing().linearSeqs, layout: [],
      losatProgram: 'tblastx', blastpMode: state.activeDrawing().losat.blastp.mode
    }))).status, 'error');
    assert.equal(searches, 0);
    assert.equal(state.errorLog.value?.code, 'COMPARISON_INPUT');
    assert.match(state.errorLog.value?.summary, /Sequence 3\. Field: LOSAT Gencode\./);
    assert.match(state.errorLog.value?.summary, /Genetic codes LOSAT reads: 1–6, 9–16, 21–33\./);
    // The refusal leads to its cause only: no Retry, which refuses the same code again (OV-3317 rule).
    assert.deepEqual(state.errorLog.value?.actions, ['edit-comparison']);
    // Review B-2: LOSATN reads no genetic code, so the Gencode 7 left on the card (hidden under LOSATN) refuses nothing.
    state.activeDrawing().losatProgram.value = 'blastn';
    state.errorLog.value = null;
    await runner.runAnalysis(resolveLinearComparisonPlan({
      plan: state.activeDrawing().linearComparisonPlan, sequences: state.activeDrawing().linearSeqs, layout: [],
      losatProgram: 'blastn', blastpMode: state.activeDrawing().losat.blastp.mode
    }));
    assert.notEqual(state.errorLog.value?.code, 'COMPARISON_INPUT', state.errorLog.value?.summary);
    assert.ok(searches > 0, 'LOSATN searches');
    state.activeDrawing().linearSeqs[2].losat_gencode = 1;
    state.errorLog.value = null;
  } finally {
    state.similarityAlignmentPlan.value = null;
    state.linearRecordTranslations.value = [];
    if (previousTestHooks === undefined) {
      delete globalThis.__GBDRAW_TEST_HOOKS__;
    } else {
      globalThis.__GBDRAW_TEST_HOOKS__ = previousTestHooks;
    }
    if (previousLosatExecutor === undefined) {
      delete globalThis.__GBDRAW_LOSAT_EXECUTOR__;
    } else {
      globalThis.__GBDRAW_LOSAT_EXECUTOR__ = previousLosatExecutor;
    }
  }
});

test('the Circular ring subject database, its raw cache key and the ring rows equal the CLI ones (D-25, D-28)', async () => {
  // The Python oracle plans the CLI request that draws what each Web state below draws (records OFF, a single
  // presentation, a crop, GFF3 with an unnamed FASTA entry, a byte order mark) and returns the subject FASTA its
  // ring search reads, the raw cache key, the ring rows the CLI keeps for the same search output, and the Worker
  // helper's reply for a file the browser reader declines.
  const oracleOptions = {
    cwd: fileURLToPath(new URL('../../', import.meta.url)), encoding: /** @type {const} */ ('utf8'),
    maxBuffer: 16 * 1024 * 1024
  };
  const oracle = spawnSync(process.env.PYTHON || 'python', ['tests/web/circular-ring-subject-oracle.py'], oracleOptions);
  assert.equal(oracle.status, 0, `Python ring subject oracle failed:\n${oracle.stderr}`);
  const { query, queryReply, cases } = JSON.parse(oracle.stdout);
  const webStates = {
    'records not drawn': { selector: '', multi: true, off: ['#2'], region: [null, null] },
    'single presentation': { selector: '#2', multi: false, off: [], region: [null, null] },
    crop: { selector: '#2', multi: false, off: [], region: [11, 50] },
    'gff3 with an unnamed fasta entry': { selector: '', multi: false, off: [], region: [null, null] },
    'byte order mark': { selector: '#1', multi: false, off: [], region: [null, null] }
  };
  assert.deepEqual(cases.map(({ name }) => name), Object.keys(webStates));

  const queryFile = new AuditFile([query], 'ring-query.fasta');
  installShownDrawing('circular');
  const drawing = state.activeDrawing();
  drawing.form.prefix = 'ring-subject';
  drawing.form.show_depth = false;
  drawing.adv.circular_track_slots_enabled = false;
  drawing.adv.circular_track_slots = [];
  Object.assign(drawing.circularConservation, {
    enabled: true, source: 'losat', losat_program: 'blastn', subject_gencode: 1, reference: 'auto', labels: ''
  });
  drawing.circularConservation.series.splice(0, drawing.circularConservation.series.length, {
    sourceKey: `${queryFile.name}|${queryFile.size}|1|0`, fileName: queryFile.name, sourceIndex: 0,
    label: 'query', color: '#4e79a7', losat_gencode: 1
  });
  drawing.annotationSets.splice(0);
  Object.assign(modeDrawing(state, 'circular').files, boundFiles({
    c_depth: null, c_conservation_blasts: [], c_conservation_fastas: [queryFile],
    c_conservation_sequence_sources: []
  }));
  const runner = wireGeneratedArtifactRuntimeOwner(createRunAnalysis({
    ...generatedArtifactHandleOptions,
    state,
    recordDisplayRows: { value: [] },
    writeSessionDocumentFile,
    serializeCanonicalFiles: (_snapshot, _catalog, active) => serializeActiveRenderFiles(state, active),
    adoptCanonicalRenderArtifacts: () => {},
    getCommittedCanonicalSession: () => null,
    prepareLinearRecordCatalog: async () => ({ catalog: { mode: 'linear', status: 'ready', records: [] }, error: '' }),
    prepareCandidateCommit: ({ generationResponse, catalogAdmission, featureColorOverrides }) => ({
      results: structuredClone(generationResponse.results),
      featureState: featureStateFromCatalog(catalogAdmission.catalog),
      featureColorOverrides: structuredClone(featureColorOverrides)
    }),
    resetPreviewViewport: () => {}
  }));
  /** @type {{ subject: string, cacheKeys: string[] } | null} */
  let searched = null;
  let searchOutput = '';
  const previousLosatExecutor = globalThis.__GBDRAW_LOSAT_EXECUTOR__;
  globalThis.__GBDRAW_LOSAT_EXECUTOR__ = async (jobs, options) => {
    const subjectKey = [...options.sequences.keys()].find((key) => key.startsWith('circular-subject:'));
    searched = { subject: options.sequences.get(subjectKey), cacheKeys: jobs.map(({ cacheKey }) => cacheKey) };
    return jobs.map(({ cacheKey }) => ({ cacheKey, text: searchOutput }));
  };
  // The Worker's ring table conversion runs the Python helper on the request and the resources the page staged.
  const ringTableConversion = {
    then(resolve, reject) {
      try {
        const converted = spawnSync(process.env.PYTHON || 'python', ['tests/web/circular-ring-subject-oracle.py', '--convert'], {
          ...oracleOptions, input: JSON.stringify(ringTableConversionInput())
        });
        assert.equal(converted.status, 0, `Python ring table conversion failed:\n${converted.stderr}`);
        resolve({ ok: true, result: JSON.parse(converted.stdout) });
      } catch (error) {
        reject(error);
      }
    }
  };
  // The ring query is read once, by the Python reader, for every case.
  workerHelperResponses.push({ ok: true, result: queryReply });
  try {
    for (const [index, testCase] of cases.entries()) {
      const webState = webStates[testCase.name];
      const files = Object.fromEntries(Object.entries(testCase.files).map(([role, data]) => [
        role, new AuditFile([Buffer.from(data, 'base64')], `ring-subject-${index}-${role}`)
      ]));
      const gff = Boolean(files.fasta);
      modeDrawing(state, 'circular').inputType.value = gff ? 'gff' : 'gb';
      setCircularCards({ gb: gff ? null : files.source, gff: gff ? files.source : null, fasta: gff ? files.fasta : null });
      // One card holds the file: a record's identity is the card's uid, with its position in the file when it holds several.
      const identityOf = (selector) => (testCase.recordList.length > 1 && selector ? `record-1:${selector.slice(1)}` : selector && 'record-1');
      modeDrawing(state, 'circular').circularRecordList.value = testCase.recordList.map(({ id, length }, recordIndex) => ({
        uid: 'record-1', k: recordIndex + 1, identity: identityOf(`#${recordIndex + 1}`), selector: `#${recordIndex + 1}`,
        record_id: id, record_length: length, detectedTopology: 'circular'
      }));
      Object.assign(modeDrawing(state, 'circular').circularRecordDiscovery, {
        status: 'ready', error: '', inputType: gff ? 'gff' : 'gb',
        sources: [{ uid: 'record-1', primaryFile: bound(files.source), pairedFile: bound(gff ? files.fasta : null), error: '' }]
      });
      Object.assign(drawing.form, {
        circular_record_selector: identityOf(webState.selector), multi_record_canvas: webState.multi,
        circular_region_start: webState.region[0], circular_region_end: webState.region[1]
      });
      drawing.recordsOff.splice(0, drawing.recordsOff.length, ...webState.off.map(identityOf));
      state.losatCache.value = new Map();
      state.losatDerivedCache.value = new Map();
      // The page reads the file with its own reader, or, for a byte order mark, through the Worker helper.
      if (testCase.name === 'byte order mark') {
        workerHelperResponses.push({ ok: true, result: testCase.helperReply });
      }
      // The ring scales Python plans for the staged files (PD-OI-018 revision 6): the subject file of several
      // records searches with the CLI's -searchsp, so the raw cache key below is the CLI's with it.
      workerHelperResponses.push({ ok: true, result: testCase.ringScales });
      workerHelperResponses.push(ringTableConversion);
      searched = null;
      searchOutput = testCase.rawRows;
      const messagesBefore = workerMessages.length;
      const caseResult = result(`ring-subject-${index}.svg`, `ring-subject-${index}`);
      workerResponses.push(response(caseResult, validCatalog(caseResult.name)));
      assert.deepEqual(await runner.runAnalysis(), { status: 'ok' }, `${testCase.name}: ${JSON.stringify(state.errorLog.value)}`);
      assert.equal(workerHelperResponses.length, 0, testCase.name);
      assert.equal(searched?.subject, testCase.cliSubjectFasta, testCase.name);
      assert.deepEqual(searched?.cacheKeys, [testCase.cliCacheKey], testCase.name);
      assert.deepEqual([...state.losatCache.value.values()].map(({ args }) => args), [testCase.cliArgs], testCase.name);
      assert(testCase.cliArgs.includes('--searchsp'), testCase.name);
      // The ring sends the lengths the ring reader answered; the reference is staged as the request names it.
      const scaleRequest = workerMessages.slice(messagesBefore).find(({ operation }) => operation === 'planRingScales');
      assert.deepEqual(scaleRequest?.payload.rings, [{ names: [queryFile.name], files: queryReply.genomeFiles }], testCase.name);
      assert.deepEqual(scaleRequest?.payload.referenceSources, testCase.referenceSources, testCase.name);
      // One conversion per Generate; the request renders (and a Session saves) the CLI's ring rows, and the raw
      // LOSAT cache keeps the search output.
      const caseMessages = workerMessages.slice(messagesBefore);
      assert.deepEqual(caseMessages.filter(({ type, operation }) => type === 'helper'
        && operation === 'convertRingTablesToSearchFrame').map(({ payload }) => payload.resourceIds),
      [['conservation-resolved-blast-1']], testCase.name);
      const run = caseMessages.findLast(({ type }) => type === 'run').payload;
      const ring = run.resourceManifest.find(({ resourceId }) => resourceId === 'conservation-resolved-blast-1');
      assert.equal(Buffer.from(stagedBytesByToken().get(ring.cacheToken)).toString('utf8'), testCase.cliRingRows, testCase.name);
      assert.deepEqual([...state.losatCache.value.values()].map(({ text }) => text), [testCase.rawRows], testCase.name);
      if (testCase.name !== 'crop') continue;
      /** @param {string} name */
      const generateAgain = async (name) => {
        const before = workerMessages.length;
        const again = result(`${name}.svg`, name);
        workerResponses.push(response(again, validCatalog(again.name)));
        assert.deepEqual(await runner.runAnalysis(), { status: 'ok' }, `${name}: ${JSON.stringify(state.errorLog.value)}`);
        const messages = workerMessages.slice(before);
        const ringToken = messages.findLast(({ type }) => type === 'run').payload.resourceManifest
          .find(({ resourceId }) => resourceId === 'conservation-resolved-blast-1').cacheToken;
        return {
          conversions: messages.filter(({ type, operation }) => type === 'helper'
            && operation === 'convertRingTablesToSearchFrame').length,
          helpers: messages.filter(({ type }) => type === 'helper').length,
          rows: Buffer.from(stagedBytesByToken().get(ringToken)).toString('utf8')
        };
      };
      // An identical Generate sends only `run`: the page holds the converted rows of this search and frame.
      assert.deepEqual(await generateAgain('ring-subject-crop-again'), { conversions: 0, helpers: 0, rows: testCase.cliRingRows });
      // A new crop is a new frame: Python converts the held search output again, without a new search.
      drawing.form.circular_region_start = 21;
      workerHelperResponses.push(ringTableConversion);
      searched = null;
      const newCrop = await generateAgain('ring-subject-new-crop');
      assert.equal(newCrop.conversions, 1);
      assert.equal(newCrop.helpers, 1);
      assert.equal(searched, null);
      assert.match(newCrop.rows, /^QUERY1\tSUBB\.1\t99\.0\t20\t0\t0\t1\t20\t1\t20\t/);
      assert.equal(workerHelperResponses.length, 0);
    }
  } finally {
    workerHelperResponses.splice(0);
    drawing.recordsOff.splice(0);
    Object.assign(drawing.form, {
      circular_record_selector: '', multi_record_canvas: false, circular_region_start: null, circular_region_end: null
    });
    drawing.circularConservation.enabled = false;
    modeDrawing(state, 'circular').files.c_conservation_fastas = bound([]);
    if (previousLosatExecutor === undefined) {
      delete globalThis.__GBDRAW_LOSAT_EXECUTOR__;
    } else {
      globalThis.__GBDRAW_LOSAT_EXECUTOR__ = previousLosatExecutor;
    }
  }
});

test('Linear mode none ignores dormant comparison state while active depth and annotations render', async () => {
  const primaryReadsBefore = activePrimaryReads;
  const inactiveReadsBefore = inactiveFileReads;
  const workerRunCountBefore = workerMessages.filter(({ type }) => type === 'run').length;
  const first = new AuditFile(['LOCUS first\nORIGIN\n        1 acgtacgt\n//\n'], 'first.gb');
  const second = new AuditFile(['LOCUS second\nORIGIN\n        1 acgtacgt\n//\n'], 'second.gb');
  const depth = new AuditFile(['position\tdepth\n1\t5\n'], 'first-depth.tsv');
  const dormantUpload = new AuditInactiveFile(['unused upload'], 'dormant-upload.tsv');
  const staleCanonical = new AuditInactiveFile(['unused canonical'], 'stale-canonical.tsv');

  installShownDrawing('linear');
  modeDrawing(state, 'linear').inputType.value = 'gb';
  state.activeDrawing().form.prefix = 'linear-none';
  state.activeDrawing().form.show_depth = true;
  state.activeDrawing().form.multi_record_canvas = false;
  state.activeDrawing().adv.linear_track_slots_enabled = false;
  state.activeDrawing().adv.linear_track_slots = [];
  state.activeDrawing().adv.depth_tracks = [{
    label: 'Coverage',
    color: '#4A90E2',
    height: null,
    large_tick_interval: null,
    small_tick_interval: null,
    tick_font_size: null
  }];
  state.activeDrawing().linearRecordLayoutEnabled.value = false;
  state.activeDrawing().linearRecordRows.splice(0);
  state.activeDrawing().linearSeqs.splice(0, state.activeDrawing().linearSeqs.length,
    {
      uid: 'linear-none-first',
      gb: bound(first),
      gff: null,
      fasta: null,
      depth: bound(depth),
      blast: dormantUpload,
      losat_gencode: 1,
      definition: '',
      record_subtitle: '',
      region_record_id: '',
      region_start: null,
      region_end: null,
      region_reverse: false,
      windows: []
    },
    {
      uid: 'linear-none-second',
      gb: bound(second),
      gff: null,
      fasta: null,
      depth: null,
      blast: dormantUpload,
      losat_gencode: 1,
      definition: '',
      record_subtitle: '',
      region_record_id: '',
      region_start: null,
      region_end: null,
      region_reverse: false,
      windows: []
    }
  );
  state.activeDrawing().linearComparisonPlan.mode = 'none';
  state.activeDrawing().linearComparisonPlan.defaultSource = 'losat';
  state.activeDrawing().linearComparisonPlan.edges.splice(0, state.activeDrawing().linearComparisonPlan.edges.length, {
    id: 'dormant-upload-edge',
    queryUid: 'linear-none-first',
    subjectUid: 'linear-none-second',
    included: true,
    source: 'upload',
    file: dormantUpload,
    fileActive: true,
    losatFilename: '',
    losatFilenameActive: false
  });
  Object.assign(state.activeDrawing().adv, {
    comparison_height: 'dormant-invalid-height',
    min_bitscore: 'dormant-invalid-bitscore',
    evalue: 'dormant-invalid-evalue',
    identity: 'dormant-invalid-identity',
    alignment_length: 'dormant-invalid-alignment',
    pairwise_match_style: 'dormant-invalid-style'
  });
  state.activeDrawing().losatProgram.value = 'dormant-invalid-program';
  Object.assign(state.activeDrawing().losat, {
    totalThreadBudget: 'dormant-invalid-budget',
    threadsPerJob: 'dormant-invalid-threads',
    parallelWorkers: 'dormant-invalid-workers'
  });
  Object.assign(state.activeDrawing().losat.blastp, {
    mode: 'dormant-invalid-presentation',
    maxHits: 'dormant-invalid-max-hits',
    orthogroupMembershipMode: 'dormant-invalid-membership',
    orthogroupMemberMaxHits: 'dormant-invalid-member-hits',
    collinearMinAnchors: 'dormant-invalid-anchors',
    collinearMaxUnitGap: 'dormant-invalid-unit-gap',
    collinearMaxDiagonalDrift: 'dormant-invalid-drift',
    collinearMaxConflictsInMergeGap: 'dormant-invalid-conflicts',
    collinearMaxParalogLinksPerOrthogroup: 'dormant-invalid-paralogs',
    collinearColorMode: 'dormant-invalid-color',
    collinearAnchorMode: 'dormant-invalid-anchor-mode',
    collinearSearchScope: 'dormant-invalid-scope'
  });
  const dormantComparisonSettingsBefore = structuredClone({
    adv: {
      comparison_height: state.activeDrawing().adv.comparison_height,
      min_bitscore: state.activeDrawing().adv.min_bitscore,
      evalue: state.activeDrawing().adv.evalue,
      identity: state.activeDrawing().adv.identity,
      alignment_length: state.activeDrawing().adv.alignment_length,
      pairwise_match_style: state.activeDrawing().adv.pairwise_match_style
    },
    losatProgram: state.activeDrawing().losatProgram.value,
    losat: state.activeDrawing().losat
  });
  const dormantRawCache = new AuditCacheMap([['dormant-raw', { stale: true }]]);
  const dormantDerivedCache = new AuditCacheMap([['dormant-derived', { stale: true }]]);
  state.losatCache.value = dormantRawCache;
  state.losatDerivedCache.value = dormantDerivedCache;
  state.losatCacheInfo.value = [{ key: 'dormant-display-cache', drawingId: 'linear' }];
  modeDrawing(state, 'linear').files.linearCanonicalComparisons = [{
    kind: 'precomputed',
    edgeKey: 'linear-none-first->linear-none-second',
    queryIndex: 0,
    subjectIndex: 1,
    file: staleCanonical
  }, {
    kind: 'generated-protein',
    mode: 'pairwise',
    pairs: [{ queryIndex: 0, subjectIndex: 1 }]
  }];
  state.activeDrawing().annotationSets.splice(0, state.activeDrawing().annotationSets.length, {
    id: 'active-annotation',
    annotations: [{
      id: 'active-window',
      target: {
        kind: 'coordinateSpan',
        record: null,
        start: 1,
        end: 4,
        coordinateSpace: 'source',
        wrapsOrigin: false,
        outOfBounds: 'clip'
      },
      label: 'Active window',
      mark: 'band',
      lane: null,
      style: null,
      legendLabel: null,
      metadata: {}
    }],
    defaultStyle: { fill: null },
    legendLabel: null
  });

  const comparisonPlanSnapshot = resolveLinearComparisonPlan({
    plan: state.activeDrawing().linearComparisonPlan,
    sequences: state.activeDrawing().linearSeqs,
    layout: [],
    losatProgram: state.activeDrawing().losatProgram.value,
    blastpMode: state.activeDrawing().losat.blastp.mode
  });
  assert.equal(Object.isFrozen(comparisonPlanSnapshot), true);
  assert.equal(comparisonPlanSnapshot.mode, 'none');
  assert.deepEqual(comparisonPlanSnapshot.edges, []);

  let losatCalls = 0;
  let annotationValidationCalls = 0;
  let annotationValidationError = '';
  let serializeCalls = 0;
  let serializedSnapshot = null;
  let serializedRecordCatalog = null;
  const preparedRecordCatalog = { mode: 'linear', status: 'ready', records: [] };
  let prepareLinearRecordCatalogImpl = async () => ({
    catalog: preparedRecordCatalog,
    error: ''
  });
  const previousLosatExecutor = globalThis.__GBDRAW_LOSAT_EXECUTOR__;
  globalThis.__GBDRAW_LOSAT_EXECUTOR__ = async () => {
    losatCalls += 1;
    throw new Error('mode none must not execute LOSAT');
  };

  const runner = wireGeneratedArtifactRuntimeOwner(createRunAnalysis({
    ...generatedArtifactHandleOptions,
    state,
    writeSessionDocumentFile,
    serializeCanonicalFiles: (snapshot, recordCatalog, drawing) => {
      serializeCalls += 1;
      serializedSnapshot = snapshot;
      serializedRecordCatalog = recordCatalog;
      return serializeActiveRenderFiles(state, drawing, snapshot);
    },
    prepareLinearRecordCatalog: (...args) => prepareLinearRecordCatalogImpl(...args),
    adoptCanonicalRenderArtifacts: () => {},
    validateAnnotationTargets: () => {
      annotationValidationCalls += 1;
      return annotationValidationError;
    },
    prepareCandidateCommit: ({
      generationResponse,
      catalogAdmission,
      featureColorOverrides
    }) => ({
      results: structuredClone(generationResponse.results),
      featureState: featureStateFromCatalog(catalogAdmission.catalog),
      featureColorOverrides: structuredClone(featureColorOverrides)
    }),
    resetPreviewViewport: ({ pan = null } = {}) => {
      state.canvasPan.x = Number(pan?.x) || 0;
      state.canvasPan.y = Number(pan?.y) || 0;
      if (!pan) state.zoom.value = 1;
    }
  }));

  const linearResult = result('linear-none.svg', 'linear-none');
  workerResponses.push(response(linearResult, validCatalog(linearResult.name)));
  // OV-1589: nothing reads a staged copy of the annotation table; a Source
  // recipe writes its own from the request.
  const stagedAnnotationTables = [];
  const stagedColorTables = [];
  const NativeMap = globalThis.Map;
  globalThis.Map = class extends NativeMap {
    set(key, value) {
      if (key === '/web_annotations.tsv' && typeof value?.data === 'string') {
        stagedAnnotationTables.push(value.data);
      }
      if ((key === '/combined_d.tsv' || key === '/combined_t.tsv') && typeof value?.data === 'string') {
        stagedColorTables.push(key);
      }
      return super.set(key, value);
    }
  };
  try {
    assert.deepEqual(
      await runner.runAnalysis(comparisonPlanSnapshot),
      { status: 'ok' },
      JSON.stringify(state.errorLog.value)
    );
  } finally {
    globalThis.Map = NativeMap;
    if (previousLosatExecutor === undefined) {
      delete globalThis.__GBDRAW_LOSAT_EXECUTOR__;
    } else {
      globalThis.__GBDRAW_LOSAT_EXECUTOR__ = previousLosatExecutor;
    }
  }

  const workerRunMessages = workerMessages.filter(({ type }) => type === 'run');
  // OV-1573: once, and no color table copies: Run Info builds those from the request resources.
  assert.deepEqual(stagedAnnotationTables, [], 'Generate stages no annotation table copy');
  assert.deepEqual(stagedColorTables, []);
  // An inherited annotation fill has no lossless TSV form, so Run Info offers no Source recipe here.
  assert.match(state.lastRunInfo.value.sourceRecipe.unavailableReason, /canonical annotation semantics/);
  const payload = workerRunMessages.at(-1).payload;
  assert.equal(workerRunMessages.length, workerRunCountBefore + 1);
  assert.equal(serializedSnapshot, comparisonPlanSnapshot);
  assert.equal(serializedRecordCatalog, preparedRecordCatalog);
  assert.equal(losatCalls, 0);
  assert.equal(dormantRawCache.probes, 0);
  assert.equal(dormantDerivedCache.probes, 0);
  assert.deepEqual(state.losatCacheInfo.value, []);
  assert.equal(annotationValidationCalls, 1);
  assert.equal(activePrimaryReads, primaryReadsBefore + 3);
  assert.equal(inactiveFileReads, inactiveReadsBefore);
  assert.deepEqual(payload.request.comparisons, []);
  assert.deepEqual({
    adv: {
      comparison_height: state.activeDrawing().adv.comparison_height,
      min_bitscore: state.activeDrawing().adv.min_bitscore,
      evalue: state.activeDrawing().adv.evalue,
      identity: state.activeDrawing().adv.identity,
      alignment_length: state.activeDrawing().adv.alignment_length,
      pairwise_match_style: state.activeDrawing().adv.pairwise_match_style
    },
    losatProgram: state.activeDrawing().losatProgram.value,
    losat: state.activeDrawing().losat
  }, dormantComparisonSettingsBefore);
  for (const field of [
    'pairwiseMatchStyle',
    'evalue',
    'bitscore',
    'identity',
    'alignmentLength'
  ]) {
    assert.equal(Object.hasOwn(payload.request.diagramOptions, field), false);
  }
  assert.equal(
    Object.hasOwn(
      payload.request.diagramOptions.configOverrides,
      'canvas.linear.comparison_height'
    ),
    false
  );
  assert.equal(
    Object.hasOwn(
      payload.request.diagramOptions.configOverrides,
      'objects.blast_match.style'
    ),
    false
  );
  assert.equal(payload.request.diagramOptions.annotations.sets[0].id, 'active-annotation');
  assert.equal(payload.request.diagramOptions.depthTracks.length, 1);
  assert.equal(
    payload.resourceManifest.some((resource) => (
      String(resource?.kind || '').includes('comparison')
    )),
    false
  );

  annotationValidationError = 'injected annotation target failure';
  assert.equal((await runner.runAnalysis(comparisonPlanSnapshot)).status, 'error');
  assert.equal(state.errorLog.value?.code, 'UNKNOWN');
  assert.doesNotMatch(JSON.stringify(state.errorLog.value), /injected annotation target failure/);
  assert.equal(serializeCalls, 1, 'invalid annotations must fail before serialization');
  assert.equal(
    workerMessages.filter(({ type }) => type === 'run').length,
    workerRunCountBefore + 1
  );
  annotationValidationError = '';

  let releaseRecordCatalog;
  prepareLinearRecordCatalogImpl = () => new Promise((resolve) => {
    releaseRecordCatalog = resolve;
  });
  const committedResults = state.results.value;
  const workerRunsBeforeCancel = workerMessages.filter(({ type }) => type === 'run').length;
  const priorOrientations = state.activeDrawing().linearSeqs.map(({ uid, region_reverse }) => ({
    recordKey: uid, reverseComplement: Boolean(region_reverse)
  }));
  const canceledRun = runner.runAnalysis(comparisonPlanSnapshot);
  for (let turn = 0; turn < 4 && typeof releaseRecordCatalog !== 'function'; turn += 1) {
    await Promise.resolve();
  }
  await runner.cancelRunAnalysis();
  releaseRecordCatalog({ catalog: preparedRecordCatalog, error: '' });
  assert.deepEqual(await canceledRun, { status: 'canceled' });
  assert.deepEqual(state.activeDrawing().linearSeqs.map(({ uid, region_reverse }) => ({
    recordKey: uid, reverseComplement: Boolean(region_reverse)
  })), priorOrientations);
  assert.equal(state.results.value, committedResults);
  assert.deepEqual(state.results.value, committedResults);
  assert.equal(state.results.value[0], committedResults[0]);
  assert.equal(state.processing.value, false);
  assert.equal(state.generationCancelRequested.value, false);
  assert.equal(serializeCalls, 1);
  assert.equal(
    workerMessages.filter(({ type }) => type === 'run').length,
    workerRunsBeforeCancel
  );

  prepareLinearRecordCatalogImpl = async () => {
    throw new Error('injected record catalog failure');
  };
  assert.equal((await runner.runAnalysis(comparisonPlanSnapshot)).status, 'error');
  assert.equal(state.errorLog.value?.code, 'UNKNOWN');
  assert.doesNotMatch(JSON.stringify(state.errorLog.value), /injected record catalog failure/);
  assert.equal(state.processing.value, false);

  prepareLinearRecordCatalogImpl = async () => ({
    catalog: preparedRecordCatalog,
    error: ''
  });
  const retryResult = result('linear-none-retry.svg', 'linear-none-retry');
  workerResponses.push(response(retryResult, validCatalog(retryResult.name)));
  assert.deepEqual(await runner.runAnalysis(comparisonPlanSnapshot), { status: 'ok' });

  let releaseSupersededCatalog;
  prepareLinearRecordCatalogImpl = () => new Promise((resolve) => {
    releaseSupersededCatalog = resolve;
  });
  const supersededRun = runner.runAnalysis(comparisonPlanSnapshot);
  while (!releaseSupersededCatalog) {
    await new Promise((resolve) => setImmediate(resolve));
  }
  prepareLinearRecordCatalogImpl = async () => ({
    catalog: preparedRecordCatalog,
    error: ''
  });
  const newestResult = result('linear-none-newest.svg', 'linear-none-newest');
  workerResponses.push(response(newestResult, validCatalog(newestResult.name)));
  assert.deepEqual(await runner.runAnalysis(comparisonPlanSnapshot), { status: 'ok' });
  assert.deepEqual(state.results.value, [newestResult]);
  releaseSupersededCatalog({ catalog: preparedRecordCatalog, error: '' });
  assert.deepEqual(await supersededRun, { status: 'stale' });
  assert.deepEqual(state.activeDrawing().linearSeqs.map(({ uid, region_reverse }) => ({
    recordKey: uid, reverseComplement: Boolean(region_reverse)
  })), priorOrientations);
  assert.deepEqual(state.results.value, [newestResult]);
  assert.equal(state.failedGeneratePreservedResult.value, false);
  assert.equal(state.processing.value, false);
});


test('decoration candidate transform runs once after catalog validation for every candidate kind', async () => {
  for (const kind of ['generate', 'reflow', 'target-record-transform']) {
    const generatedResult = result('decoration.svg', 'candidate');
    const canonical = { renderRequest: { mode: 'linear' }, resources: {} };
    let transforms = 0;
    const candidate = await executeCanonicalRenderCandidate({ canonical, mode: 'linear', kind,
      generationExecutor: async () => response(generatedResult, validCatalog(generatedResult.name)),
      decorationContinuity: (actual, admission) => {
        assert.equal(actual, canonical);
        assert.equal(admission.catalog.items[0].resultName, generatedResult.name);
        transforms++;
        return [svg => { svg.decoration = [40, 20]; }];
      },
      prepareCommit: ({ resultTransforms }) => {
        const svg = { recordDelta: [9, 7] };
        resultTransforms?.[0](svg);
        return { svg };
      }
    });
    assert.equal(transforms, 1, kind);
    assert.deepEqual(candidate.commit.svg, { recordDelta: [9, 7], decoration: [40, 20] });
  }
});

test('a render that reuses the held catalog build keeps the catalog object the page admitted (step 3a)', async () => {
  const metrics = {};
  const events = [];
  const previousHooks = globalThis.__GBDRAW_TEST_HOOKS__;
  globalThis.__GBDRAW_TEST_HOOKS__ = {
    onStructuralMetric({ name, value }) { metrics[name] = (metrics[name] || 0) + Number(value || 0); },
    onSessionLifecycleEvent({ name }) { events.push(name); }
  };
  const canonical = { renderRequest: { schema: CANONICAL_REQUEST_SCHEMA }, resources: {} };
  const run = (options = {}) => executeCanonicalRenderCandidate({
    canonical,
    mode: 'circular',
    prepareCommit: ({ generationResponse }) => ({ generationResponse }),
    ...options
  });
  const runMessages = () => workerMessages.filter(({ type }) => type === 'run');
  const lastHeld = () => runMessages().at(-1).payload.heldFeatureCatalogRevision;
  const unchanged = (logicalResult, revision) => ({
    results: [logicalResult], metadata: {}, featureCatalogRevision: revision
  });
  const delta = (before, name) => (metrics[name] || 0) - (before[name] || 0);
  const name = 'reuse.svg';
  try {
    // N1: a sent catalog C with revision R; the next render names R, and an
    // "unchanged" reply puts C itself back.
    workerResponses.push(response(result(name, 'one'), validCatalog(name)));
    const first = await run();
    assert.equal(first.status, 'ok');
    assert.equal(lastHeld(), null);
    const catalog = first.catalog;
    assert.equal(first.generationResponse.metadata.featureCatalog, catalog);
    const heldRevision = fakeCatalogRevision(`fake:${fakeCatalogBuilds}`);
    workerResponses.push(unchanged(result(name, 'two'), heldRevision));
    const before = { ...metrics };
    events.length = 0;
    const reused = await run({ heldFeatureCatalog: catalog });
    assert.equal(reused.status, 'ok');
    assert.deepEqual(lastHeld(), heldRevision);
    assert.equal(reused.catalog, catalog);
    assert.equal(reused.generationResponse.metadata.featureCatalog, reused.catalogAdmission.catalog);
    assert.equal(delta(before, 'featureCatalogAdmissionCount'), 0);
    assert.equal(delta(before, 'featureCatalogReuseCount'), 1);
    assert.ok(events.includes('catalog.reused'));
    assert.equal(events.some((event) => event.startsWith('catalog.admission-')), false);

    // N2: an object the page did not receive from the Worker, or none, names nothing.
    for (const held of [structuredClone(catalog), validCatalog(name), null]) {
      workerResponses.push(response(result(name, 'copy'), validCatalog(name)));
      assert.equal((await run({ heldFeatureCatalog: held })).status, 'ok');
      assert.equal(lastHeld(), null);
    }

    // N3 (OV-136): items renamed after receipt; a rerender that keeps those
    // names names the revision, a Generate (engine names) or other names do not.
    workerResponses.push(response(result(name, 'renamed'), validCatalog(name)));
    const renamedSource = await run();
    const renamed = renamedSource.catalog;
    const renamedRevision = fakeCatalogRevision(`fake:${fakeCatalogBuilds}`);
    renamed.items[0].resultName = 'Saved name';
    workerResponses.push(unchanged(result(name, 'kept'), renamedRevision));
    const kept = await run({ heldFeatureCatalog: renamed, kind: 'reflow', resultNames: ['Saved name'] });
    assert.equal(kept.status, 'ok');
    assert.deepEqual(lastHeld(), renamedRevision);
    assert.equal(kept.catalog, renamed);
    for (const resultNames of [null, ['Other name']]) {
      workerResponses.push(response(result(name, 'resent'), validCatalog(name)));
      await run({ heldFeatureCatalog: renamed, resultNames }).catch(() => null);
      assert.equal(lastHeld(), null);
    }

    // N4: "unchanged" for another revision, with nothing held, or with no revision
    // fails the reply and ends the Worker.
    for (const [held, reply] of [
      [catalog, unchanged(result(name, 'other'), fakeCatalogRevision('fake:other'))],
      [null, unchanged(result(name, 'none-held'), heldRevision)],
      [catalog, { results: [result(name, 'neither')], metadata: {}, featureCatalogRevision: null }]
    ]) {
      const constructed = metrics.workerConstructionCount || 0;
      workerResponses.push(reply);
      await assert.rejects(run({ heldFeatureCatalog: held }), (error) => error.code === 'RESULT_INVALID');
      workerResponses.push(response(result(name, 'after-failure'), validCatalog(name)));
      assert.equal((await run()).status, 'ok');
      assert.equal((metrics.workerConstructionCount || 0), constructed + 1);
    }

    // N5: a sent catalog while one is held is admitted as a new object; C keeps its record.
    workerResponses.push(response(result(name, 'changed'), validCatalog(name)));
    const changed = await run({ heldFeatureCatalog: catalog });
    assert.deepEqual(lastHeld(), heldRevision);
    assert.notEqual(changed.catalog, catalog);
    workerResponses.push(unchanged(result(name, 'undo'), heldRevision));
    assert.equal((await run({ heldFeatureCatalog: catalog })).catalog, catalog);

    // N6: an engine error leaves the held catalog and its record.
    workerResponses.push({ error: { code: 'INPUT_INVALID', operation: 'generate', stage: 'render' } });
    assert.equal((await run({ heldFeatureCatalog: catalog })).status, 'engine-error');
    workerResponses.push(unchanged(result(name, 'after-engine-error'), heldRevision));
    await run({ heldFeatureCatalog: catalog });
    assert.deepEqual(lastHeld(), heldRevision);
  } finally {
    globalThis.__GBDRAW_TEST_HOOKS__ = previousHooks;
  }
});

// Step 3c B3 (D-53): after the first status a Generate yields one task, not two
// frames, at "Preparing render inputs and session..." and "Preparing preview...".
// The yield is the production one (`afterTask`, counted here), so a Cancel
// queued as a task while a status shows (a click during busy work) runs
// before the run continues (#1119 review N-2): GE-09 keeps
// the warm Worker before the post, and a Cancel during admission ends the run
// before it commits or mounts the candidate.
test('a Generate yields one frame wait and two task yields, and a Cancel queued at either later status ends the run there (3c B3)', async () => {
  const previousHooks = globalThis.__GBDRAW_TEST_HOOKS__;
  const events = [];
  globalThis.__GBDRAW_TEST_HOOKS__ = { onSessionLifecycleEvent({ name }) { events.push(name); } };
  installSimpleCircularInput(
    new AuditFile(['LOCUS       audit\nORIGIN\n//\n'], 'yield.gb', { type: 'text/plain', lastModified: 9 }),
    new AuditInactiveFile(['unused'], 'inactive.dat')
  );
  const waits = { afterPaint: 0, task: 0 };
  // The status a queued Cancel targets; the setter queues it as a task the
  // moment the run writes that status.
  let cancelAt = '';
  const cancelReturns = [];
  const statusRef = state.processingStatus;
  let statusValue = statusRef.value;
  Object.defineProperty(statusRef, 'value', {
    configurable: true,
    get: () => statusValue,
    set: (value) => {
      statusValue = value;
      if (cancelAt && value === cancelAt) {
        cancelAt = '';
        setTimeout(() => cancelReturns.push(runner.cancelRunAnalysis()), 0);
      }
    }
  });
  const runner = wireGeneratedArtifactRuntimeOwner(createRunAnalysis({
    ...generatedArtifactHandleOptions,
    state,
    previewRuntime: createImmediatePreviewRuntime({ onEvent: (name) => events.push(name) }),
    waitForAfterPaint: async () => { waits.afterPaint += 1; },
    yieldToTask: () => {
      waits.task += 1;
      return afterTask();
    },
    writeSessionDocumentFile,
    serializeCanonicalFiles: (_snapshot, _catalog, drawing) => serializeActiveRenderFiles(state, drawing),
    adoptCanonicalRenderArtifacts: () => {},
    prepareCandidateCommit: ({ generationResponse, catalogAdmission, featureColorOverrides }) => ({
      results: structuredClone(generationResponse.results),
      featureState: featureStateFromCatalog(catalogAdmission.catalog),
      featureColorOverrides: structuredClone(featureColorOverrides)
    }),
    resetPreviewViewport: () => {}
  }));
  const runPosts = () => workerMessages.filter(({ type }) => type === 'run').length;
  const statusOrder = (values) => values.filter((value) => [
    'Preparing input files...', 'Preparing render inputs and session...', 'Preparing preview...'
  ].includes(value)).filter((value, index, kept) => value !== kept[index - 1]);
  try {
    const committed = result('yield.svg', 'committed');
    workerResponses.push(response(committed, validCatalog(committed.name)));
    const statuses = [];
    const recordStatuses = Object.getOwnPropertyDescriptor(statusRef, 'value').set;
    Object.defineProperty(statusRef, 'value', {
      configurable: true,
      get: () => statusValue,
      set: (value) => { statuses.push(value); recordStatuses(value); }
    });
    assert.deepEqual(await runner.runAnalysis(), { status: 'ok' }, JSON.stringify(state.errorLog.value));
    assert.deepEqual(waits, { afterPaint: 1, task: 2 });
    assert.deepEqual(statusOrder(statuses), [
      'Preparing input files...', 'Preparing render inputs and session...', 'Preparing preview...'
    ]);
    assert.ok(events.includes('generate.paint-opportunity-completed'));
    const committedResults = state.results.value;

    // Cancel queued at "Preparing render inputs and session...": seen before
    // the post, with no Worker request to end (GE-09).
    const postsBefore = runPosts();
    events.length = 0;
    cancelAt = 'Preparing render inputs and session...';
    await runner.runAnalysis();
    assert.equal(runPosts(), postsBefore);
    assert.deepEqual(cancelReturns.splice(0), [false]);
    assert.equal(state.processingStatus.value, 'Canceled.');
    assert.equal(state.results.value, committedResults);

    // Cancel queued at "Preparing preview...": seen before the commit, so the
    // candidate is neither committed nor mounted.
    const candidate = result('yield.svg', 'candidate');
    workerResponses.push(response(candidate, validCatalog(candidate.name)));
    events.length = 0;
    cancelAt = 'Preparing preview...';
    await runner.runAnalysis();
    assert.equal(runPosts(), postsBefore + 1);
    assert.deepEqual(cancelReturns.splice(0), [false]);
    assert.equal(state.processingStatus.value, 'Canceled.');
    assert.equal(state.results.value, committedResults);
    assert.equal(events.includes('test.preview-mounted'), false);
    assert.deepEqual(waits, { afterPaint: 3, task: 5 });
  } finally {
    delete statusRef.value;
    statusRef.value = statusValue;
    globalThis.__GBDRAW_TEST_HOOKS__ = previousHooks;
  }
});

// Generate draws a Linear Depth row at its own height on every run while the
// row's series has no height of its own: the settings Generate gives a series
// leave its height Auto, so a second Generate from the same state draws the
// same figure (OV-1510, review 4 finding 2).
test('two Linear Generates from one state draw a Depth row at the same height', async () => {
  const genbank = new AuditFile(['LOCUS rowheight\nORIGIN\n        1 acgtacgt\n//\n'], 'row-height.gb');
  const depthFile = (n) => new AuditFile([`position\tdepth\n1\t${n}\n`], `row-height-depth-${n}.tsv`);
  installShownDrawing('linear');
  modeDrawing(state, 'linear').inputType.value = 'gb';
  const drawing = state.activeDrawing();
  drawing.form.prefix = 'row-height';
  drawing.form.show_depth = true;
  drawing.form.multi_record_canvas = false;
  drawing.linearRecordLayoutEnabled.value = false;
  drawing.linearRecordRows.splice(0);
  drawing.linearComparisonPlan.mode = 'none';
  drawing.linearComparisonPlan.edges.splice(0);
  drawing.annotationSets.splice(0);
  state.activeDrawing().linearSeqs.splice(0, state.activeDrawing().linearSeqs.length, {
    uid: 'row-height',
    gb: bound(genbank),
    gff: null,
    fasta: null,
    depth: bound([depthFile(1), depthFile(2)]),
    blast: null,
    losat_gencode: 1,
    definition: '',
    record_subtitle: '',
    region_record_id: '',
    region_start: null,
    region_end: null,
    region_reverse: false,
      windows: []
  });
  const depthRow = (id, trackIndex, height) => ({
    id, renderer: 'depth', enabled: true, side: 'below', height, spacing: '', z: 0, params: { track_index: trackIndex }
  });
  Object.assign(drawing.adv, {
    depth_height: 50,
    // Series 2 has no stored settings, as after a Session that saved one series.
    depth_tracks: [{
      label: '', color: '#4A90E2', height: null, large_tick_interval: null, small_tick_interval: null, tick_font_size: null
    }],
    linear_track_slots_enabled: true,
    linear_track_slots_axis_index: 0,
    linear_track_slots: [
      { id: 'features', renderer: 'features', enabled: true, side: 'overlay', height: '', spacing: '', z: 0, params: {} },
      depthRow('depth_1', 0, ''),
      depthRow('depth_2', 1, '30px')
    ]
  });
  const runner = wireGeneratedArtifactRuntimeOwner(createRunAnalysis({
    ...generatedArtifactHandleOptions,
    state,
    writeSessionDocumentFile,
    serializeCanonicalFiles: (snapshot, recordCatalog, activeDrawing) => (
      serializeActiveRenderFiles(state, activeDrawing, snapshot)
    ),
    prepareLinearRecordCatalog: async () => ({ catalog: { mode: 'linear', status: 'ready', records: [] }, error: '' }),
    adoptCanonicalRenderArtifacts: () => {},
    validateAnnotationTargets: () => '',
    prepareCandidateCommit: ({ generationResponse, catalogAdmission, featureColorOverrides }) => ({
      results: structuredClone(generationResponse.results),
      featureState: featureStateFromCatalog(catalogAdmission.catalog),
      featureColorOverrides: structuredClone(featureColorOverrides)
    }),
    resetPreviewViewport: () => {}
  }));
  const generate = async (name) => {
    const generated = result(`${name}.svg`, name);
    workerResponses.push(response(generated, validCatalog(generated.name)));
    assert.deepEqual(await runner.runAnalysis(), { status: 'ok' }, JSON.stringify(state.errorLog.value));
    return {
      row: drawing.adv.linear_track_slots.find((slot) => slot.id === 'depth_2')?.height,
      seriesHeights: drawing.adv.depth_tracks.map((track) => track.height)
    };
  };
  const first = await generate('row-height-1');
  const second = await generate('row-height-2');
  assert.deepEqual([first.row, second.row], ['30px', '30px']);
  assert.deepEqual([first.seriesHeights, second.seriesHeights], [[null, null], [null, null]]);
});

// OV-3010: the Run Info replay file of a Generate is built on demand from the
// request's resources. An input file's bytes stay in its File (the History file
// store counts them once), so the estimate of what an Undo step keeps for the
// replay file does not count them again.
test('the History estimate of the Run Info replay file leaves out the input file bytes (OV-3010)', async () => {
  const sequence = 'acgtacgtac'.repeat(100_000);
  const genbank = new AuditFile([`LOCUS largeinput\nORIGIN\n        1 ${sequence}\n//\n`], 'large-input.gb');
  installShownDrawing('linear');
  modeDrawing(state, 'linear').inputType.value = 'gb';
  const drawing = state.activeDrawing();
  drawing.form.prefix = 'large-input';
  drawing.form.show_depth = false;
  drawing.form.multi_record_canvas = false;
  drawing.linearRecordLayoutEnabled.value = false;
  drawing.linearRecordRows.splice(0);
  drawing.linearComparisonPlan.mode = 'none';
  drawing.linearComparisonPlan.edges.splice(0);
  drawing.annotationSets.splice(0);
  drawing.linearSeqs.splice(0, drawing.linearSeqs.length, {
    uid: 'large-input', gb: bound(genbank), gff: null, fasta: null, depth: null, blast: null, losat_gencode: 1,
    definition: '', record_subtitle: '', region_record_id: '', region_start: null, region_end: null,
    region_reverse: false, windows: []
  });
  Object.assign(drawing.adv, { linear_track_slots_enabled: false, linear_track_slots: [] });
  const runner = wireGeneratedArtifactRuntimeOwner(createRunAnalysis({
    ...generatedArtifactHandleOptions,
    state,
    writeSessionDocumentFile,
    serializeCanonicalFiles: (snapshot, recordCatalog, activeDrawing) => (
      serializeActiveRenderFiles(state, activeDrawing, snapshot)
    ),
    prepareLinearRecordCatalog: async () => ({ catalog: { mode: 'linear', status: 'ready', records: [] }, error: '' }),
    adoptCanonicalRenderArtifacts: () => {},
    validateAnnotationTargets: () => '',
    prepareCandidateCommit: ({ generationResponse, catalogAdmission, featureColorOverrides }) => ({
      results: structuredClone(generationResponse.results),
      featureState: featureStateFromCatalog(catalogAdmission.catalog),
      featureColorOverrides: structuredClone(featureColorOverrides)
    }),
    resetPreviewViewport: () => {}
  }));
  const generated = result('large-input.svg', 'large-input');
  workerResponses.push(response(generated, validCatalog(generated.name)));
  assert.deepEqual(await runner.runAnalysis(), { status: 'ok' }, JSON.stringify(state.errorLog.value));
  const runtime = runner.captureGeneratedArtifactRuntimeState();
  assert.deepEqual(runtime.latestCliHelperFiles.map((file) => file.name), ['large-input.gbdraw']);
  assert.ok(genbank.size > 1_000_000);
  assert.ok(runtime.retainedBytes < genbank.size, `${runtime.retainedBytes} counts the ${genbank.size} input bytes`);
});

// Circular Generate gives a Depth series without settings the color the
// series owner gives it, as the editor shows: series 2 the palette, not the
// Depth color of series 1 (OV-1510).
test('Circular Generate gives a Depth series without settings its default color', async () => {
  const depthFile = (n) => new AuditFile([`position\tdepth\n1\t${n}\n`], `circular-series-${n}.tsv`);
  installShownDrawing('circular');
  modeDrawing(state, 'circular').inputType.value = 'gb';
  const drawing = state.activeDrawing();
  drawing.form.prefix = 'circular-series';
  drawing.form.show_depth = true;
  drawing.form.multi_record_canvas = false;
  drawing.annotationSets.splice(0);
  drawing.circularConservation.enabled = false;
  drawing.circularConservation.series.splice(0);
  setCircularFile(new AuditFile(['LOCUS circularseries\nORIGIN\n        1 acgtacgt\n//\n'], 'circular-series.gb'));
  Object.assign(modeDrawing(state, 'circular').files, boundFiles({
    c_depth: [depthFile(1), depthFile(2)],
    c_conservation_blasts: [],
    c_conservation_fastas: [],
    c_conservation_sequence_sources: []
  }));
  modeDrawing(state, 'circular').circularRecordList.value = [];
  Object.assign(modeDrawing(state, 'circular').circularRecordDiscovery, {
    status: 'idle', error: '', inputType: '', sources: []
  });
  const depthRow = (id, trackIndex) => ({ id, renderer: 'depth', enabled: true, side: 'inside', params: { track_index: trackIndex } });
  Object.assign(drawing.adv, {
    depth_color: '#AA0000',
    // Series 2 has no stored settings, as after a Session that saved one series.
    depth_tracks: [{
      label: '', color: '#AA0000', height: null, large_tick_interval: null, small_tick_interval: null, tick_font_size: null
    }],
    circular_track_slots_enabled: true,
    circular_track_slots: [depthRow('depth_1', 0), depthRow('depth_2', 1)]
  });
  const runner = wireGeneratedArtifactRuntimeOwner(createRunAnalysis({
    ...generatedArtifactHandleOptions,
    state,
    recordDisplayRows: { value: [] },
    writeSessionDocumentFile,
    serializeCanonicalFiles: (_snapshot, _catalog, activeDrawing) => (
      serializeActiveRenderFiles(state, activeDrawing)
    ),
    adoptCanonicalRenderArtifacts: () => {},
    validateAnnotationTargets: () => '',
    prepareCandidateCommit: ({ generationResponse, catalogAdmission, featureColorOverrides }) => ({
      results: structuredClone(generationResponse.results),
      featureState: featureStateFromCatalog(catalogAdmission.catalog),
      featureColorOverrides: structuredClone(featureColorOverrides)
    }),
    resetPreviewViewport: () => {}
  }));
  // The record discovery of the new source; earlier tests may leave replies queued.
  workerHelperResponses.splice(0, workerHelperResponses.length, {
    ok: true,
    result: { records: [{ selector: '#1', record_id: 'circularseries', record_length: 8 }] }
  });
  const generated = result('circular-series.svg', 'circular-series');
  workerResponses.push(response(generated, validCatalog(generated.name)));
  assert.deepEqual(await runner.runAnalysis(), { status: 'ok' }, JSON.stringify(state.errorLog.value));
  const request = workerMessages.filter(({ type }) => type === 'run').at(-1).payload.request;
  assert.deepEqual(
    [drawing.adv.depth_tracks.map((track) => track.color), request.diagramOptions.depthTracks.map((track) => track.color)],
    [['#AA0000', '#E45756'], ['#AA0000', '#E45756']]
  );
});

// The Circular record refresh reads into the drawing its caller gives
// (Generate: the drawing it pinned), not the drawing named by the mode, and
// shares an in-flight read only with calls for the same drawing.
test('the Circular record refresh writes the drawing it reads and keys its read by the drawing', async () => {
  installShownDrawing('circular');
  const named = state.activeDrawing();
  const second = createDrawingState('circular-2', 'circular');
  state.addDrawing(second, 'Circular 2', 2);
  try {
    installShownDrawing(second.id);
    modeDrawing(state, 'circular').inputType.value = 'gb';
    setCircularFile(new AuditFile(['LOCUS refresh 10 bp DNA circular\n//\n'], 'refresh-owner.gb'), second);
    second.circularRecordList.value = [];
    Object.assign(second.circularRecordDiscovery, {
      status: 'idle', error: '', inputType: '', sources: []
    });
    named.adv.multi_record_positions.splice(0, named.adv.multi_record_positions.length);
    const runner = createRunAnalysis({ ...generatedArtifactHandleOptions, state, resetPreviewViewport: () => {} });
    const reply = { ok: true, result: { records: [{ selector: '#1', record_id: 'refresh', record_length: 10 }] } };
    workerHelperResponses.splice(0, workerHelperResponses.length, reply);
    await runner.refreshCircularRecordOrder({}, second);
    assert.deepEqual(
      [second.adv.multi_record_positions.map((position) => position.selector), named.adv.multi_record_positions],
      [['record-1'], []]
    );

    let release;
    workerHelperResponses.splice(0, workerHelperResponses.length, new Promise((resolve) => { release = () => resolve(reply); }), reply);
    second.circularRecordDiscovery.status = 'idle';
    const pinned = runner.refreshCircularRecordOrder({}, second);
    assert.equal(runner.refreshCircularRecordOrder({ automatic: true }, second), pinned, 'a call for the same drawing shares the read');
    assert.notEqual(runner.refreshCircularRecordOrder({}, named), pinned, 'another drawing of the mode reads on its own');
    release();
    await pinned;
  } finally {
    installShownDrawing('circular');
    state.removeDrawing(second.id, { takeSlot: () => null });
  }
});

// D-15-207: the read still running belongs to the drawing object too. A
// drawing added under the ID of a removed drawing whose read still runs (as a
// Load installs new drawings under reused IDs), with the same file bound,
// reads its own records: it does not join the old drawing's read, which
// writes only the old drawing.
test('a drawing re-added under its ID reads its own records while the old drawing\'s read runs', async () => {
  const id = 'circular-rejoined';
  const file = bound(new AuditFile(['LOCUS rejoined 8 bp DNA circular\nORIGIN\n        1 acgtacgt\n//\n'], 'rejoined.gb'));
  const reply = { ok: true, result: { records: [{ selector: '#1', record_id: 'rejoined', record_length: 8 }] } };
  const runner = createRunAnalysis({ ...generatedArtifactHandleOptions, state, resetPreviewViewport: () => {} });
  const idle = { status: 'idle', error: '', inputType: '', sources: [] };
  const old = createDrawingState(id, 'circular');
  state.addDrawing(old, 'Old', 2);
  old.inputType.value = 'gb';
  replaceCircularSources(old.slices.inputs, [createCircularSource({ uid: 'record-1', gb: file })]);
  Object.assign(old.circularRecordDiscovery, idle);
  let release = () => {};
  workerHelperResponses.splice(0, workerHelperResponses.length,
    new Promise((resolve) => { release = () => resolve(reply); }), reply);
  const oldRead = runner.refreshCircularRecordOrder({}, old);
  state.removeDrawing(id, { takeSlot: () => null });
  const fresh = createDrawingState(id, 'circular');
  state.addDrawing(fresh, 'Fresh', 2);
  try {
    fresh.inputType.value = 'gb';
    replaceCircularSources(fresh.slices.inputs, [createCircularSource({ uid: 'record-1', gb: file })]);
    Object.assign(fresh.circularRecordDiscovery, idle);
    const freshRead = runner.refreshCircularRecordOrder({}, fresh);
    assert.notEqual(freshRead, oldRead, 'the new drawing does not join the old drawing\'s read');
    release();
    await Promise.all([oldRead, freshRead]);
    assert.deepEqual(
      [fresh.circularRecordDiscovery.status, fresh.circularRecordList.value.map(({ selector }) => selector),
        fresh.adv.multi_record_positions.map(({ selector }) => selector)],
      ['ready', ['#1'], ['record-1']],
      'the new drawing shows its records and their Record Order'
    );
  } finally {
    state.removeDrawing(id, { takeSlot: () => null });
  }
});

// D-15-207: a drawing's record-read bookkeeping belongs to the drawing
// object. A drawing removed and re-added under its ID (as a Load installs new
// drawings under reused IDs) does not inherit the old drawing's read that the
// region dialog left unapplied: its Generate writes no Record Order from it.
test('a drawing re-added under its ID does not inherit the old drawing\'s unapplied record read', async () => {
  const id = 'circular-reused';
  const source = (name) => bound(new AuditFile([`LOCUS ${name} 8 bp DNA circular\nORIGIN\n        1 acgtacgt\n//\n`], `${name}.gb`));
  const reply = { ok: true, result: { records: [{ selector: '#1', record_id: 'reused', record_length: 8 }] } };
  const runner = wireGeneratedArtifactRuntimeOwner(createRunAnalysis({
    ...generatedArtifactHandleOptions,
    state,
    recordDisplayRows: { value: [] },
    readRecordIdentitySources: async () => ({ mode: 'circular', cards: new Map() }),
    writeSessionDocumentFile,
    serializeCanonicalFiles: (_snapshot, _catalog, activeDrawing) => serializeActiveRenderFiles(state, activeDrawing),
    adoptCanonicalRenderArtifacts: () => {},
    validateAnnotationTargets: () => '',
    prepareCandidateCommit: ({ generationResponse, catalogAdmission, featureColorOverrides }) => ({
      results: structuredClone(generationResponse.results),
      featureState: featureStateFromCatalog(catalogAdmission.catalog),
      featureColorOverrides: structuredClone(featureColorOverrides)
    }),
    resetPreviewViewport: () => {}
  }));
  const old = createDrawingState(id, 'circular');
  state.addDrawing(old, 'Old', 2);
  old.inputType.value = 'gb';
  setCircularFile(source('old'), old);
  // The state a Load leaves: the list deferred for the current file.
  Object.assign(old.circularRecordDiscovery, {
    status: 'deferred', error: '', inputType: 'gb', sources: [{ uid: 'record-1', primaryFile: old.circularSources[0].gb, pairedFile: null, error: '' }]
  });
  workerHelperResponses.splice(0, workerHelperResponses.length, reply);
  // The region dialog's read: the list, no Record Order (D-15-205).
  await runner.readDrawingRecordIdentities(old);
  assert.deepEqual(old.adv.multi_record_positions, []);
  state.removeDrawing(id, { takeSlot: () => null });
  const fresh = createDrawingState(id, 'circular');
  state.addDrawing(fresh, 'Fresh', 2);
  try {
    installShownDrawing(id);
    fresh.inputType.value = 'gb';
    fresh.form.prefix = 'circular-reused';
    fresh.form.multi_record_canvas = false;
    fresh.circularConservation.enabled = false;
    setCircularFile(source('fresh'), fresh);
    fresh.circularRecordList.value = [{ uid: 'record-1', k: 1, identity: 'record-1', selector: '#1', record_id: 'reused', record_length: 8, detectedTopology: 'circular' }];
    Object.assign(fresh.circularRecordDiscovery, {
      status: 'ready', error: '', inputType: 'gb', sources: [{ uid: 'record-1', primaryFile: fresh.circularSources[0].gb, pairedFile: null, error: '' }]
    });
    const generated = result('circular-reused.svg', 'circular-reused');
    workerResponses.push(response(generated, validCatalog(generated.name)));
    assert.deepEqual(await runner.runAnalysis(), { status: 'ok' }, JSON.stringify(state.errorLog.value));
    assert.deepEqual(fresh.adv.multi_record_positions, [], 'the old drawing\'s read is not applied to the new one');
  } finally {
    installShownDrawing('circular');
    state.removeDrawing(id, { takeSlot: () => null });
  }
});

// PD-OI-121 scenario revision 3: Generate passes a gzip input as chosen and the Worker's Python expands it, so the
// ring plan's reference is staged as Generate's request names it (its gzip bytes and their size), and the plan reads
// the expanded record (its length is the search space's; the gzip case of circular-ring-evalue-scale.playwright.spec.js).
test('a gzip Circular reference stages for the LOSAT ring plan as Generate names it, its gzip bytes that expand to the record', async () => {
  const { gzipSync, gunzipSync } = await import('node:zlib');
  const { requestRecordSources } = await import('../../gbdraw/web/js/services/session-request.js');
  const { createDiagramResourceTransport } = await import('../../gbdraw/web/js/services/diagram-resource-staging.js');
  const text = 'LOCUS       REF1                    2000 bp    DNA     circular UNK 01-JAN-1980\n' +
    `ORIGIN\n        1 ${'acgt'.repeat(500)}\n//\n`;
  const chosen = gzipSync(Buffer.from(text));
  const file = new File([chosen], 'ref.gb.gz', { lastModified: 1 });
  // The ring planner stages the reference under the ids the Generate's request names it by (planRingScales).
  // Its descriptors are those Generate's request writes (`serializeFile`), as the page passes them.
  const reference = await requestRecordSources([{ gb: file, gff: null, fasta: null }], 'genbank', serializeFile);
  assert.ok(reference);
  assert.equal(reference.resources['record-1-genbank'].size, chosen.length);
  const prepared = await createDiagramResourceTransport().prepare({
    request: reference.sources, resources: reference.resources, retain: true
  });
  assert.deepEqual(prepared.resourceManifest.map(({ resourceId, size }) => [resourceId, size]),
    [['record-1-genbank', chosen.length]]);
  const [staged] = prepared.stagedResources;
  const stagedBytes = Buffer.from(await staged.blob.arrayBuffer());
  assert.deepEqual(stagedBytes, chosen);
  assert.equal(gunzipSync(stagedBytes).toString('utf8'), text);
});
