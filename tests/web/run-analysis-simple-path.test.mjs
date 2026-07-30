import assert from 'node:assert/strict';
import { webcrypto } from 'node:crypto';
import test from 'node:test';

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
globalThis.document = {};
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

globalThis.File = AuditFile;

const { EXPECTED_WEB_RUNTIME_CAPABILITIES } = await import(
  '../../gbdraw/web/js/services/runtime-capabilities.js'
);

const workerResponses = [];
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
    throw new Error(`unexpected worker message: ${String(message.type)}`);
  }

  terminate() {}
}

globalThis.Worker = AuditSimplePathWorker;

const { createRunAnalysis } = await import('../../gbdraw/web/js/app/run-analysis.js');
const {
  serializeActiveRenderFiles,
  SESSION_VERSION
} = await import('../../gbdraw/web/js/services/config.js');
const {
  featureStateFromCatalog
} = await import('../../gbdraw/web/js/services/feature-catalog.js');
const { state } = await import('../../gbdraw/web/js/state.js');

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

test('audit-5 owner: direct simple createRunAnalysis path is worker-only and catalog-transactional', async () => {
  const primary = new AuditFile(['LOCUS audit\nORIGIN\n//\n'], 'active.gb', {
    type: 'text/plain',
    lastModified: 7
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

  let ensurePyodideCalls = 0;
  let ensurePyodideImpl = async () => {
    throw new Error('simple render must not initialize Pyodide');
  };
  let adoptedArtifacts = 0;
  let failArtifactAdoption = false;
  const runner = createRunAnalysis({
    state,
    getPyodide: () => null,
    ensurePyodide: async () => {
      ensurePyodideCalls += 1;
      return ensurePyodideImpl();
    },
    writeFileToFs: async () => {
      throw new Error('simple render must not write through Pyodide');
    },
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
  });

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

  const committedState = committedFeatureState();

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
  assert.equal(state.zoom.value, 1.7);
  assert.deepEqual(state.canvasPan, { x: 31, y: -12 });

  assert.equal(ensurePyodideCalls, 0);
  assert.equal(activePrimaryReads, 1);
  assert.equal(inactiveFileReads, 0);
  assert.equal(adoptedArtifacts, 1);
  assert.equal(workerMessages.filter(({ type }) => type === 'run').length, 4);
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
  assert.equal(ensurePyodideCalls, 0);
  assert.equal(activePrimaryReads, 1);

  Object.assign(state.circularRecordDiscovery, {
    status: 'idle',
    error: '',
    inputType: '',
    primaryFile: null,
    pairedFile: null
  });
  const workerRunCountBeforeDiscoveryFailure = workerMessages
    .filter(({ type }) => type === 'run')
    .length;
  assert.deepEqual(await runner.runAnalysis(), { status: 'error' });
  assert.equal(ensurePyodideCalls, 1);
  assert.equal(
    workerMessages.filter(({ type }) => type === 'run').length,
    workerRunCountBeforeDiscoveryFailure
  );
  assert.match(
    state.errorLog.value?.summary || '',
    /Could not start the record discovery helper/
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
  assert.deepEqual(await runner.runAnalysis(), { status: 'error' });
  assert.equal(ensurePyodideCalls, 2);
  assert.equal(
    workerMessages.filter(({ type }) => type === 'run').length,
    workerRunCountBeforeDiscoveryFailure
  );

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
  ensurePyodideImpl = () => new Promise((resolve) => {
    releaseDiscovery = resolve;
  });
  const canceledRun = runner.runAnalysis();
  assert.equal(ensurePyodideCalls, 3);
  await runner.cancelRunAnalysis();
  releaseDiscovery();
  assert.deepEqual(await canceledRun, { status: 'canceled' });
  assert.equal(state.processing.value, false);
  assert.equal(state.generationCancelRequested.value, false);
  assert.equal(
    workerMessages.filter(({ type }) => type === 'run').length,
    workerRunCountBeforeDiscoveryFailure
  );
});
