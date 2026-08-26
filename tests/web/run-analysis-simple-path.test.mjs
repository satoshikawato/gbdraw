import assert from 'node:assert/strict';
import { webcrypto } from 'node:crypto';
import test from 'node:test';
import { installFakeSvgDom } from './fake-svg-dom.mjs';

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
globalThis.document = {
  body: { appendChild: (node) => { node.parentNode = globalThis.document.body; } },
  createElement: () => ({
    addEventListener() {},
    click() { downloadedFilename = this.download; },
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
let inactiveFileReads = 0;

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
const workerMessages = [];

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
      queueMicrotask(() => this.emit('message', {
        type: 'run',
        requestId: message.requestId,
        ok: true,
        results: structuredClone(response)
      }));
      return;
    }
    if (message.type === 'helper') {
      const response = workerHelperResponses.shift();
      if (!response) throw new Error('missing helper worker response');
      Promise.resolve(response).then((payload) => this.emit('message', {
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

const { createRunAnalysis } = await import('../../gbdraw/web/js/app/run-analysis.js');
const {
  resolveLinearComparisonPlan
} = await import('../../gbdraw/web/js/app/linear-comparisons.js');
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
  serializeResults,
  SESSION_VERSION
} = await import('../../gbdraw/web/js/services/config.js');
const { createHistoryFileStore } = await import(
  '../../gbdraw/web/js/services/history-files.js'
);
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
const { normalizeLinearSeqList, state } = await import('../../gbdraw/web/js/state.js');

const artifactSnapshots = createHistorySnapshotService({
  state,
  fileStore: createHistoryFileStore(),
  nextTick: window.Vue.nextTick,
  normalizeLinearSeqList,
  buildUiStateData,
  applyUiStateData,
  buildFeatureStateData,
  applyFeatureStateData,
  buildEditorStateData,
  applyEditorStateData,
  buildOrthogroupStateData,
  applyOrthogroupStateData,
  serializeResults,
  applyResultsData,
  buildRunStateData,
  applyRunStateData
});
const generatedArtifactHandleOptions = {
  captureGeneratedArtifactHandle: artifactSnapshots.captureGeneratedArtifactHandle,
  restoreGeneratedArtifactHandle: artifactSnapshots.restoreGeneratedArtifactHandle,
  setGeneratedArtifactIdentity: artifactSnapshots.setGeneratedArtifactIdentity
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
  schema: 3,
  items: [{
    resultIndex: 0,
    resultName: name,
    recordKeys: ['record-1'],
    features: [{
      svgId: 'rendered-feature-1',
      recordKey: 'record-1',
      biologicalFeatureId: 'biological-feature-1',
      fillColor: '#abcdef'
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
      qualifiers: { product: ['audit protein'] }
    }],
    orthogroups: [],
    annotations: [],
    comparisonMatches: []
  }]
});

const response = (logicalResult, featureCatalog) => ({
  results: [logicalResult],
  metadata: featureCatalog === undefined ? {} : { featureCatalog }
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
      text: decoder.decode(bytes.subarray(dataStart, dataStart + size))
    });
    offset = dataStart + size;
  }
  return entries;
};

const committedFeatureState = () => structuredClone({
  results: state.results.value,
  selectedResultIndex: state.selectedResultIndex.value,
  featureCatalog: state.featureCatalog.value,
  extractedFeatures: state.extractedFeatures.value,
  biologicalFeatures: state.biologicalFeatures.value,
  featureSelectorSafetyScope: state.featureSelectorSafetyScope.value,
  selectedFeatureIds: [...state.selectedFeatureIds.value],
  selectedFeatureAnchorId: state.selectedFeatureAnchorId.value,
  featureSelectionStatus: state.featureSelectionStatus.value,
  featureColorOverrides: state.featureColorOverrides,
  featureStrokeOverrides: state.featureStrokeOverrides,
  legendEntries: state.legendEntries.value,
  featureEditorStatus: state.featureEditorStatus,
  featureExtractionPending: state.featureExtractionPending.value,
  featureExtractionError: state.featureExtractionError.value,
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

test('audit-5 owner: direct simple createRunAnalysis path is worker-only and catalog-transactional', async () => {
  const structuralMetrics = {};
  globalThis.__GBDRAW_TEST_HOOKS__ = {
    onStructuralMetric(metric) {
      structuralMetrics[metric.name] = (
        Number(structuralMetrics[metric.name] || 0) + Number(metric.value || 0)
      );
    }
  };
  const primary = new AuditFile(['LOCUS audit\nORIGIN\n//\n'], 'active.gb', {
    type: 'text/plain',
    lastModified: 7
  });
  const fallbackPrimary = new AuditFile(['record parser fallback'], 'fallback.gb', {
    type: 'application/octet-stream',
    lastModified: 8
  });
  const inactive = new AuditInactiveFile(['unused'], 'inactive.dat');

  state.mode.value = 'circular';
  state.cInputType.value = 'gb';
  state.files.c_gb = primary;
  state.files.c_gff = null;
  state.files.c_fasta = null;
  state.files.c_depth = inactive;
  state.files.c_conservation_blasts = [];
  state.files.c_conservation_fastas = [];
  state.files.c_conservation_sequence_sources = [];
  state.files.d_color = inactive;
  state.files.t_color = inactive;
  state.files.blacklist = inactive;
  state.files.whitelist = inactive;
  state.files.qualifier_priority = inactive;
  state.form.prefix = 'audit-simple';
  state.form.show_depth = false;
  state.form.multi_record_canvas = false;
  state.adv.circular_track_slots_enabled = false;
  state.circularConservation.enabled = false;
  state.circularConservation.source = 'upload';
  state.annotationSets.splice(0);
  state.circularRecordList.value = [{ selector: '#1', record_id: 'audit' }];
  Object.assign(state.circularRecordDiscovery, {
    status: 'ready',
    error: '',
    inputType: 'gb',
    primaryFile: primary,
    pairedFile: null
  });
  state.linearSeqs.splice(0, state.linearSeqs.length, {
    uid: 'inactive-linear-row',
    gb: inactive,
    gff: null,
    fasta: null,
    depth: inactive,
    blast: null,
    losat_gencode: 1,
    losat_filename: '',
    definition: '',
    record_subtitle: '',
    region_record_id: '',
    region_start: null,
    region_end: null,
    region_reverse: false
  });

  let adoptedArtifacts = 0;
  let failArtifactAdoption = false;
  let cancelDuringCandidate = false;
  let runner;
  runner = wireGeneratedArtifactRuntimeOwner(createRunAnalysis({
    ...generatedArtifactHandleOptions,
    state,
    serializeCanonicalFiles: () => serializeActiveRenderFiles(state.mode.value, state),
    canonicalSessionVersion: SESSION_VERSION,
    adoptCanonicalRenderArtifacts: () => {
      if (failArtifactAdoption) {
        throw new Error('forced late canonical artifact adoption failure');
      }
      adoptedArtifacts += 1;
    },
    prepareCandidateCommit: ({
      results,
      catalog,
      featureColorOverrides,
      featureStrokeOverrides
    }) => {
      const candidate = {
        results: structuredClone(results),
        featureState: featureStateFromCatalog(catalog),
        featureColorOverrides: structuredClone(featureColorOverrides),
        featureStrokeOverrides: structuredClone(featureStrokeOverrides)
      };
      if (cancelDuringCandidate) runner.cancelRunAnalysis();
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

  assert.deepEqual(
    await runner.runAnalysis(),
    { status: 'ok' },
    JSON.stringify(state.errorLog.value)
  );
  assert.deepEqual(state.results.value, [committedResult]);
  assert.deepEqual(state.featureCatalog.value, committedCatalog);
  assert.equal(state.extractedFeatures.value.length, 1);
  assert.equal(state.biologicalFeatures.value.length, 1);
  assert.equal(structuralMetrics.canonicalReplayFullSerializationCount || 0, 0);
  runner.downloadCliHelperFiles();
  assert.equal(structuralMetrics.canonicalReplayFullSerializationCount, 1);
  assert.equal(downloadedFilename, 'audit-simple-cli-files.zip');
  const helperEntries = await storedZipEntries(downloadedBlob);
  const replayEntry = helperEntries.find(({ name }) => name.endsWith('.gbdraw-session.json'));
  assert.ok(replayEntry);
  const replay = JSON.parse(replayEntry.text);
  assert.equal(replay.format, 'gbdraw-session');
  assert.equal(replay.version, SESSION_VERSION);
  assert.equal(replay.renderRequest.schema, 6);
  assert.ok(Object.keys(replay.resources).length > 0);

  const firstRunPayload = workerMessages.find(({ type }) => type === 'run').payload;
  assert.equal(Object.hasOwn(firstRunPayload, 'resources'), false);
  assert.ok(firstRunPayload.resourceManifest.length > 0);
  assert.ok(firstRunPayload.stagedResources.every(({ bytes }) => bytes instanceof ArrayBuffer));

  const committedState = committedFeatureState();
  const committedExtractedFeatureIdentity = state.extractedFeatures.value;
  const committedBiologicalFeatureIdentity = state.biologicalFeatures.value;

  workerResponses.push(response(result('missing.svg', 'missing'), undefined));
  assert.deepEqual(await runner.runAnalysis(), { status: 'error' });
  assert.match(
    state.errorLog.value?.summary || '',
    /incompatible feature metadata/
  );
  assert.deepEqual(committedFeatureState(), committedState);

  workerResponses.push(response(
    result('malformed.svg', 'malformed'),
    { schema: 2, items: [] }
  ));
  assert.deepEqual(await runner.runAnalysis(), { status: 'error' });
  assert.match(
    state.errorLog.value?.summary || '',
    /incompatible feature metadata/
  );
  assert.deepEqual(committedFeatureState(), committedState);

  state.zoom.value = 1.7;
  state.canvasPan.x = 31;
  state.canvasPan.y = -12;
  const lateFailureState = committedFeatureState();
  failArtifactAdoption = true;
  const lateFailureResult = result('late-failure.svg', 'late-failure');
  workerResponses.push(response(
    lateFailureResult,
    validCatalog(lateFailureResult.name)
  ));
  assert.deepEqual(await runner.runAnalysis(), { status: 'error' });
  failArtifactAdoption = false;
  assert.match(
    state.errorLog.value?.summary || '',
    /forced late canonical artifact adoption failure/
  );
  assert.deepEqual(committedFeatureState(), lateFailureState);
  assert.equal(state.extractedFeatures.value, committedExtractedFeatureIdentity);
  assert.equal(state.biologicalFeatures.value, committedBiologicalFeatureIdentity);
  assert.equal(state.zoom.value, 1.7);
  assert.deepEqual(state.canvasPan, { x: 31, y: -12 });

  const canceledState = committedFeatureState();
  const canceledResultIdentity = state.results.value;
  const canceledResult = result('canceled.svg', 'canceled');
  workerResponses.push(response(canceledResult, validCatalog(canceledResult.name)));
  cancelDuringCandidate = true;
  assert.deepEqual(await runner.runAnalysis(), { status: 'canceled' });
  cancelDuringCandidate = false;
  assert.deepEqual(committedFeatureState(), canceledState);
  assert.notEqual(state.results.value, canceledResultIdentity);
  assert.deepEqual(state.results.value, canceledResultIdentity);
  assert.equal(state.results.value[0], canceledResultIdentity[0]);
  assert.equal(state.extractedFeatures.value, committedExtractedFeatureIdentity);
  assert.equal(state.biologicalFeatures.value, committedBiologicalFeatureIdentity);
  assert.equal(state.processingStatus.value, 'Canceled.');

  assert.equal(activePrimaryReads, 1);
  assert.equal(inactiveFileReads, 0);
  assert.equal(adoptedArtifacts, 1);
  assert.equal(workerMessages.filter(({ type }) => type === 'run').length, 5);
  assert.equal(workerMessages.filter(({ type }) => type === 'feature-extraction').length, 0);

  state.form.multi_record_canvas = true;
  state.circularRecordList.value = [{ selector: '#1', record_id: 'audit' }];
  Object.assign(state.circularRecordDiscovery, {
    status: 'ready',
    error: '',
    inputType: 'gb',
    primaryFile: primary,
    pairedFile: null
  });
  const readyResult = result('ready.svg', 'ready');
  workerResponses.push(response(readyResult, validCatalog(readyResult.name)));
  assert.deepEqual(await runner.runAnalysis(), { status: 'ok' });
  assert.equal(activePrimaryReads, 2);

  Object.assign(state.circularRecordDiscovery, {
    status: 'idle',
    error: '',
    inputType: '',
    primaryFile: null,
    pairedFile: null
  });
  state.files.c_gb = fallbackPrimary;
  const workerRunCountBeforeDiscoveryFailure = workerMessages
    .filter(({ type }) => type === 'run')
    .length;
  const workerHelperCountBeforeDiscoveryFailure = workerMessages
    .filter(({ type }) => type === 'helper')
    .length;
  workerHelperResponses.push({
    ok: false,
    error: { name: 'Error', message: 'injected record discovery helper failure' }
  });
  assert.deepEqual(await runner.runAnalysis(), { status: 'error' });
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
    /Could not read records from the circular input file/
  );
  assert.doesNotMatch(
    JSON.stringify(state.errorLog.value),
    /Traceback|simple render must not initialize Pyodide/
  );

  state.form.multi_record_canvas = false;
  state.form.show_depth = false;
  state.files.c_depth = [[inactive]];
  state.adv.circular_track_slots_enabled = true;
  state.adv.circular_track_slots = [{
    id: 'custom-depth',
    renderer: 'depth',
    enabled: true,
    side: 'inside',
    params: { track_index: 0 }
  }];
  Object.assign(state.circularRecordDiscovery, {
    status: 'idle',
    error: '',
    inputType: '',
    primaryFile: null,
    pairedFile: null
  });
  workerHelperResponses.push({
    ok: false,
    error: { name: 'Error', message: 'injected depth record discovery failure' }
  });
  assert.deepEqual(await runner.runAnalysis(), { status: 'error' });
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
  Object.assign(state.circularRecordDiscovery, {
    status: 'idle', error: '', inputType: '', primaryFile: null, pairedFile: null
  });
  const firstCoalescedDiscovery = runner.refreshCircularRecordOrder();
  const secondCoalescedDiscovery = runner.refreshCircularRecordOrder();
  assert.equal(firstCoalescedDiscovery, secondCoalescedDiscovery);
  await Promise.resolve();
  releaseCoalescedDiscovery();
  await Promise.all([firstCoalescedDiscovery, secondCoalescedDiscovery]);

  state.form.multi_record_canvas = true;
  state.files.c_depth = null;
  state.adv.circular_track_slots_enabled = false;
  state.adv.circular_track_slots = [];
  state.circularRecordList.value = [];
  Object.assign(state.circularRecordDiscovery, {
    status: 'idle',
    error: '',
    inputType: '',
    primaryFile: null,
    pairedFile: null
  });
  let releaseDiscovery;
  workerHelperResponses.push(new Promise((resolve) => {
    releaseDiscovery = () => resolve({
      ok: true,
      result: {
        records: [{ selector: '#1', record_id: 'audit', record_length: 10 }]
      }
    });
  }));
  const canceledRun = runner.runAnalysis();
  await Promise.resolve();
  await runner.cancelRunAnalysis();
  releaseDiscovery();
  assert.deepEqual(await canceledRun, { status: 'canceled' });
  assert.equal(state.processing.value, false);
  assert.equal(state.generationCancelRequested.value, false);
  assert.equal(
    workerMessages.filter(({ type }) => type === 'run').length,
    workerRunCountBeforeDiscoveryFailure
  );

  state.form.multi_record_canvas = true;
  state.files.c_gb = fallbackPrimary;
  state.adv.multi_record_positions.splice(
    0,
    state.adv.multi_record_positions.length,
    { selector: '#1', row: 2 }
  );
  state.circularRecordList.value = [{
    selector: '#1',
    record_id: 'preserved-record',
    record_length: 123
  }];
  state.semanticFileWatchersSuppressed.value = false;
  state.sessionImportRollbackInProgress.value = false;
  let releaseStaleDiscovery;
  workerHelperResponses.push(new Promise((resolve) => {
    releaseStaleDiscovery = () => resolve({
      ok: false,
      error: { name: 'Error', message: 'injected restored-file reparse failure' }
    });
  }));
  const staleDiscovery = runner.refreshCircularRecordOrder();
  await Promise.resolve();
  const discoveryStateBeforeRollback = {
    records: structuredClone(state.circularRecordList.value),
    positions: structuredClone(state.adv.multi_record_positions),
    discovery: { ...state.circularRecordDiscovery }
  };
  state.semanticFileWatchersSuppressed.value = true;
  await runner.refreshCircularRecordOrder({ suppress: true });
  releaseStaleDiscovery();
  await staleDiscovery;
  assert.deepEqual(state.circularRecordList.value, discoveryStateBeforeRollback.records);
  assert.deepEqual(state.adv.multi_record_positions, discoveryStateBeforeRollback.positions);
  assert.deepEqual(state.circularRecordDiscovery, discoveryStateBeforeRollback.discovery);
  state.semanticFileWatchersSuppressed.value = false;
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
  const comparisonTexts = [
    '>comparison-b\nTTTTAAAA\n',
    '>comparison-c\nAAAATTTT\n',
    '>comparison-d\nCCCCGGGG\n'
  ];
  const resourceTable = adoptCurrentSessionResources({
    'reference-a': encodedResource('genbank', 'reference-a.gb', referenceText),
    'comparison-b': encodedResource(
      'conservation-fasta-file',
      'comparison-b.fasta',
      comparisonTexts[0]
    ),
    'comparison-c': encodedResource(
      'conservation-fasta-file',
      'comparison-c.fasta',
      comparisonTexts[1]
    ),
    'comparison-d': encodedResource(
      'conservation-fasta-file',
      'comparison-d.fasta',
      comparisonTexts[2]
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

  state.mode.value = 'circular';
  state.cInputType.value = 'gb';
  state.form.prefix = 'lazy-circular';
  state.form.show_depth = false;
  state.form.multi_record_canvas = false;
  state.adv.circular_track_slots_enabled = false;
  state.adv.circular_track_slots = [];
  Object.assign(state.adv, {
    min_bitscore: 0,
    evalue: '1e-2',
    identity: 0,
    alignment_length: 0
  });
  Object.assign(state.circularConservation, {
    enabled: true,
    source: 'losat',
    losat_program: 'blastn',
    subject_gencode: 1,
    reference: 'auto',
    labels: '',
    ring_width: 5,
    ring_gap: 2
  });
  state.circularConservation.series.splice(
    0,
    state.circularConservation.series.length,
    ...comparisonViews.map((file, index) => ({
      sourceKey: `${file.name}|${file.size}|0|${index}`,
      fileName: file.name,
      sourceIndex: index,
      label: ['comparison-b', 'comparison-c', 'comparison-d'][index],
      color: ['#4e79a7', '#e15759', '#59a14f'][index],
      losat_gencode: 1
    }))
  );
  state.files.c_gb = referenceView;
  state.files.c_gff = null;
  state.files.c_fasta = null;
  state.files.c_depth = null;
  state.files.c_conservation_blasts = [];
  state.files.c_conservation_fastas = comparisonViews;
  state.files.c_conservation_sequence_sources = [];
  state.annotationSets.splice(0);
  state.circularRecordList.value = [];
  Object.assign(state.circularRecordDiscovery, {
    status: 'idle', error: '', inputType: '', primaryFile: null, pairedFile: null
  });

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
      subjectCanonicalHash
    }]);
  }
  state.losatCache.value = new Map(cacheEntries);
  state.losatDerivedCache.value = new Map();

  const runner = wireGeneratedArtifactRuntimeOwner(createRunAnalysis({
    ...generatedArtifactHandleOptions,
    state,
    serializeCanonicalFiles: () => serializeActiveRenderFiles(state.mode.value, state),
    canonicalSessionVersion: SESSION_VERSION,
    adoptCanonicalRenderArtifacts: () => {},
    prepareLinearRecordCatalog: async () => ({
      catalog: { mode: 'linear', status: 'ready', records: [] },
      error: ''
    }),
    prepareCandidateCommit: ({
      results,
      catalog,
      featureColorOverrides,
      featureStrokeOverrides
    }) => ({
      results: structuredClone(results),
      featureState: featureStateFromCatalog(catalog),
      featureColorOverrides: structuredClone(featureColorOverrides),
      featureStrokeOverrides: structuredClone(featureStrokeOverrides)
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
      text: 'SECOND\tTHIRD\t100\t4\t0\t0\t1\t4\t1\t4\t1e-10\t20\n'
    }));
  };

  try {
    const circularResult = result('lazy-circular.svg', 'lazy-circular');
    workerResponses.push(response(circularResult, validCatalog(circularResult.name)));
    assert.deepEqual(
      await runner.runAnalysis(),
      { status: 'ok' },
      JSON.stringify(state.errorLog.value)
    );
    assert.equal(losatCalls, 0, 'verified circular cache entries must prevent LOSAT execution');
    assert.deepEqual(
      resourceMetrics.filter(({ name }) => name === 'resourceTextReadCount')
        .map(({ resourceId }) => resourceId),
      ['reference-a', 'comparison-b', 'comparison-c', 'comparison-d']
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

    const invalidFile = Object.freeze({
      name: 'missing.gb', type: 'text/plain', size: 0, lastModified: 0
    });
    state.files.c_gb = invalidFile;
    state.circularRecordList.value = [
      { selector: '#1', record_id: 'missing', record_length: 1 }
    ];
    Object.assign(state.circularRecordDiscovery, {
      status: 'ready',
      error: '',
      inputType: 'gb',
      primaryFile: invalidFile,
      pairedFile: null
    });
    assert.deepEqual(await runner.runAnalysis(), { status: 'error' });
    assert.match(
      state.errorLog.value?.summary || '',
      /A File-like object with arrayBuffer\(\) or text\(\) is required/
    );

    resourceMetrics.splice(0);
    const linearTable = adoptCurrentSessionResources({
      multi: encodedResource('genbank', 'multi.gb', [
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
      third: encodedResource('genbank', 'third.gb', [
        'LOCUS       THIRD 8 bp DNA',
        'ORIGIN',
        '        1 ttttaaaa',
        '//',
        ''
      ].join('\n'))
    });
    const multiView = createSessionResourceFileView(linearTable, 'multi');
    const thirdView = createSessionResourceFileView(linearTable, 'third');
    state.mode.value = 'linear';
    state.lInputType.value = 'gb';
    state.form.prefix = 'lazy-linear';
    state.form.show_depth = false;
    state.form.linear_track_layout = 'middle';
    state.adv.linear_track_slots_enabled = false;
    state.adv.linear_track_slots = [];
    state.adv.comparison_height = null;
    state.linearRecordLayoutEnabled.value = false;
    state.linearRecordRows.splice(0);
    state.linearSeqs.splice(0, state.linearSeqs.length, {
      uid: 'multi',
      gb: multiView,
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
      region_reverse: true
    }, {
      uid: 'third',
      gb: thirdView,
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
      region_reverse: false
    });
    Object.assign(state.linearComparisonPlan, {
      mode: 'selected',
      defaultSource: 'losat'
    });
    state.linearComparisonPlan.edges.splice(0, state.linearComparisonPlan.edges.length, {
      id: 'lazy-linear-edge',
      queryUid: 'multi',
      subjectUid: 'third',
      included: true,
      source: 'losat',
      file: null,
      fileActive: false,
      losatFilename: '',
      losatFilenameActive: false
    });
    state.files.linearCanonicalComparisons = [];
    state.losatProgram.value = 'blastn';
    state.losat.blastn.task = 'megablast';
    state.losat.executionMode = 'serial';
    state.losatCache.value = new Map();

    const comparisonPlanSnapshot = resolveLinearComparisonPlan({
      plan: state.linearComparisonPlan,
      sequences: state.linearSeqs,
      layout: [],
      losatProgram: state.losatProgram.value,
      blastpMode: state.losat.blastp.mode
    });
    const linearResult = result('lazy-linear.svg', 'lazy-linear');
    workerHelperResponses.push({
      ok: true,
      result: {
        tsv: 'SECOND\tTHIRD\t100\t4\t0\t0\t1\t4\t1\t4\t1e-10\t20\n'
      }
    });
    workerResponses.push(response(linearResult, validCatalog(linearResult.name)));
    assert.deepEqual(
      await runner.runAnalysis(comparisonPlanSnapshot),
      { status: 'ok' },
      JSON.stringify({ error: state.errorLog.value, lastWorker: workerMessages.at(-1) })
    );
    assert.equal(losatCalls, 1);
    assert.ok(capturedSequences.includes('>SECOND\nACCG\n'));
    assert.ok(capturedSequences.includes('>THIRD\nTTTTAAAA\n'));
    const conversionRequest = workerMessages.findLast((message) => (
      message.type === 'helper' &&
      message.operation === 'convertLosatNucleotideToDisplayTsv'
    ));
    assert.deepEqual(conversionRequest?.payload.queryViewTransform, {
      length: 4,
      reverse: true
    });
    assert.equal(
      resourceMetrics.filter(({ name }) => name === 'resourceByteReadCount').length,
      2
    );
    assert.equal(
      resourceMetrics.filter(({ name }) => name === 'resourceTextReadCount').length,
      0,
      'the explicit staged text fast path must not materialize file text again'
    );

    const warmCacheInfo = [{
      key: 'warm-cache', edgeKey: 'multi->third', queryUid: 'multi', subjectUid: 'third',
      queryIndex: 0, subjectIndex: 1, ordinal: 0, display: true
    }];
    state.losatCacheInfo.value = structuredClone(warmCacheInfo);
    state.losatProgram.value = 'blastp';
    state.losat.blastp.mode = 'pairwise';
    state.files.linearCanonicalComparisons = [{
      kind: 'generatedProteinComparison', mode: 'none', pairs: [], settings: {}
    }];
    const warmComparisonPlan = resolveLinearComparisonPlan({
      plan: state.linearComparisonPlan,
      sequences: state.linearSeqs,
      layout: [],
      losatProgram: state.losatProgram.value,
      blastpMode: state.losat.blastp.mode
    });
    const warmResult = result('lazy-linear-warm.svg', 'lazy-linear-warm');
    workerResponses.push(response(warmResult, validCatalog(warmResult.name)));
    assert.deepEqual(
      await runner.runAnalysis(warmComparisonPlan),
      { status: 'ok' },
      JSON.stringify(state.errorLog.value)
    );
    assert.equal(losatCalls, 1, 'resolved protein artifacts must bypass LOSAT execution');
    assert.deepEqual(state.losatCacheInfo.value, warmCacheInfo);
  } finally {
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

test('Linear mode none ignores dormant comparison state while active depth and annotations render', async () => {
  const primaryReadsBefore = activePrimaryReads;
  const inactiveReadsBefore = inactiveFileReads;
  const workerRunCountBefore = workerMessages.filter(({ type }) => type === 'run').length;
  const first = new AuditFile(['LOCUS first\nORIGIN\n        1 acgtacgt\n//\n'], 'first.gb');
  const second = new AuditFile(['LOCUS second\nORIGIN\n        1 acgtacgt\n//\n'], 'second.gb');
  const depth = new AuditFile(['position\tdepth\n1\t5\n'], 'first-depth.tsv');
  const dormantUpload = new AuditInactiveFile(['unused upload'], 'dormant-upload.tsv');
  const staleCanonical = new AuditInactiveFile(['unused canonical'], 'stale-canonical.tsv');

  state.mode.value = 'linear';
  state.lInputType.value = 'gb';
  state.form.prefix = 'linear-none';
  state.form.show_depth = true;
  state.form.multi_record_canvas = false;
  state.adv.linear_track_slots_enabled = false;
  state.adv.linear_track_slots = [];
  state.adv.depth_tracks = [{
    label: 'Coverage',
    color: '#4A90E2',
    height: null,
    large_tick_interval: null,
    small_tick_interval: null,
    tick_font_size: null
  }];
  state.linearRecordLayoutEnabled.value = false;
  state.linearRecordRows.splice(0);
  state.linearSeqs.splice(0, state.linearSeqs.length,
    {
      uid: 'linear-none-first',
      gb: first,
      gff: null,
      fasta: null,
      depth,
      blast: dormantUpload,
      losat_gencode: 1,
      definition: '',
      record_subtitle: '',
      region_record_id: '',
      region_start: null,
      region_end: null,
      region_reverse: false
    },
    {
      uid: 'linear-none-second',
      gb: second,
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
      region_reverse: false
    }
  );
  state.linearComparisonPlan.mode = 'none';
  state.linearComparisonPlan.defaultSource = 'losat';
  state.linearComparisonPlan.edges.splice(0, state.linearComparisonPlan.edges.length, {
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
  Object.assign(state.adv, {
    comparison_height: 'dormant-invalid-height',
    min_bitscore: 'dormant-invalid-bitscore',
    evalue: 'dormant-invalid-evalue',
    identity: 'dormant-invalid-identity',
    alignment_length: 'dormant-invalid-alignment',
    pairwise_match_style: 'dormant-invalid-style'
  });
  state.losatProgram.value = 'dormant-invalid-program';
  Object.assign(state.losat, {
    totalThreadBudget: 'dormant-invalid-budget',
    threadsPerJob: 'dormant-invalid-threads',
    parallelWorkers: 'dormant-invalid-workers'
  });
  Object.assign(state.losat.blastp, {
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
      comparison_height: state.adv.comparison_height,
      min_bitscore: state.adv.min_bitscore,
      evalue: state.adv.evalue,
      identity: state.adv.identity,
      alignment_length: state.adv.alignment_length,
      pairwise_match_style: state.adv.pairwise_match_style
    },
    losatProgram: state.losatProgram.value,
    losat: state.losat
  });
  const dormantRawCache = new AuditCacheMap([['dormant-raw', { stale: true }]]);
  const dormantDerivedCache = new AuditCacheMap([['dormant-derived', { stale: true }]]);
  state.losatCache.value = dormantRawCache;
  state.losatDerivedCache.value = dormantDerivedCache;
  state.losatCacheInfo.value = [{ key: 'dormant-display-cache' }];
  state.files.linearCanonicalComparisons = [{
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
  state.annotationSets.splice(0, state.annotationSets.length, {
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
    defaultStyle: null,
    legendLabel: null
  });

  const comparisonPlanSnapshot = resolveLinearComparisonPlan({
    plan: state.linearComparisonPlan,
    sequences: state.linearSeqs,
    layout: [],
    losatProgram: state.losatProgram.value,
    blastpMode: state.losat.blastp.mode
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
    serializeCanonicalFiles: (snapshot, recordCatalog) => {
      serializeCalls += 1;
      serializedSnapshot = snapshot;
      serializedRecordCatalog = recordCatalog;
      return serializeActiveRenderFiles(state.mode.value, state, snapshot);
    },
    prepareLinearRecordCatalog: (...args) => prepareLinearRecordCatalogImpl(...args),
    canonicalSessionVersion: SESSION_VERSION,
    adoptCanonicalRenderArtifacts: () => {},
    validateAnnotationTargets: ({ loadComparison }) => {
      annotationValidationCalls += 1;
      assert.equal(loadComparison, false);
      return annotationValidationError;
    },
    prepareCandidateCommit: ({
      results,
      catalog,
      featureColorOverrides,
      featureStrokeOverrides
    }) => ({
      results: structuredClone(results),
      featureState: featureStateFromCatalog(catalog),
      featureColorOverrides: structuredClone(featureColorOverrides),
      featureStrokeOverrides: structuredClone(featureStrokeOverrides)
    }),
    resetPreviewViewport: ({ pan = null } = {}) => {
      state.canvasPan.x = Number(pan?.x) || 0;
      state.canvasPan.y = Number(pan?.y) || 0;
      if (!pan) state.zoom.value = 1;
    }
  }));

  const linearResult = result('linear-none.svg', 'linear-none');
  workerResponses.push(response(linearResult, validCatalog(linearResult.name)));
  try {
    assert.deepEqual(
      await runner.runAnalysis(comparisonPlanSnapshot),
      { status: 'ok' },
      JSON.stringify(state.errorLog.value)
    );
  } finally {
    if (previousLosatExecutor === undefined) {
      delete globalThis.__GBDRAW_LOSAT_EXECUTOR__;
    } else {
      globalThis.__GBDRAW_LOSAT_EXECUTOR__ = previousLosatExecutor;
    }
  }

  const workerRunMessages = workerMessages.filter(({ type }) => type === 'run');
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
      comparison_height: state.adv.comparison_height,
      min_bitscore: state.adv.min_bitscore,
      evalue: state.adv.evalue,
      identity: state.adv.identity,
      alignment_length: state.adv.alignment_length,
      pairwise_match_style: state.adv.pairwise_match_style
    },
    losatProgram: state.losatProgram.value,
    losat: state.losat
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
  assert.deepEqual(await runner.runAnalysis(comparisonPlanSnapshot), { status: 'error' });
  assert.match(state.errorLog.value?.summary || '', /injected annotation target failure/);
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
  const canceledRun = runner.runAnalysis(comparisonPlanSnapshot);
  await Promise.resolve();
  await runner.cancelRunAnalysis();
  releaseRecordCatalog({ catalog: preparedRecordCatalog, error: '' });
  assert.deepEqual(await canceledRun, { status: 'canceled' });
  assert.notEqual(state.results.value, committedResults);
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
  assert.deepEqual(await runner.runAnalysis(comparisonPlanSnapshot), { status: 'error' });
  assert.match(state.errorLog.value?.summary || '', /injected record catalog failure/);
  assert.equal(state.processing.value, false);

  prepareLinearRecordCatalogImpl = async () => ({
    catalog: preparedRecordCatalog,
    error: ''
  });
  const retryResult = result('linear-none-retry.svg', 'linear-none-retry');
  workerResponses.push(response(retryResult, validCatalog(retryResult.name)));
  assert.deepEqual(await runner.runAnalysis(comparisonPlanSnapshot), { status: 'ok' });
});
