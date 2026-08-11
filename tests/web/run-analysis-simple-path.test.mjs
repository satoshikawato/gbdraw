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
const generatedArtifactSnapshotOptions = {
  buildGeneratedArtifactSnapshot: artifactSnapshots.buildGeneratedArtifactSnapshot,
  applyGeneratedArtifactSnapshot: artifactSnapshots.applyGeneratedArtifactSnapshot
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

  let ensurePyodideCalls = 0;
  let ensurePyodideImpl = async () => {
    throw new Error('simple render must not initialize Pyodide');
  };
  let adoptedArtifacts = 0;
  let failArtifactAdoption = false;
  let cancelDuringCandidate = false;
  let runner;
  runner = createRunAnalysis({
    ...generatedArtifactSnapshotOptions,
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
  assert.equal(state.results.value, canceledResultIdentity);
  assert.equal(state.extractedFeatures.value, committedExtractedFeatureIdentity);
  assert.equal(state.biologicalFeatures.value, committedBiologicalFeatureIdentity);
  assert.equal(state.processingStatus.value, 'Canceled.');

  assert.equal(ensurePyodideCalls, 0);
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
  assert.equal(ensurePyodideCalls, 0);
  assert.equal(activePrimaryReads, 1);

  Object.assign(state.circularRecordDiscovery, {
    status: 'idle',
    error: '',
    inputType: '',
    primaryFile: null,
    pairedFile: null
  });
  const primaryText = primary.text.bind(primary);
  Object.defineProperty(primary, 'text', { value: undefined, configurable: true });
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

  let releaseCoalescedDiscovery;
  ensurePyodideImpl = () => new Promise((resolve) => {
    releaseCoalescedDiscovery = resolve;
  });
  Object.assign(state.circularRecordDiscovery, {
    status: 'idle', error: '', inputType: '', primaryFile: null, pairedFile: null
  });
  const firstCoalescedDiscovery = runner.refreshCircularRecordOrder();
  const secondCoalescedDiscovery = runner.refreshCircularRecordOrder();
  assert.equal(firstCoalescedDiscovery, secondCoalescedDiscovery);
  assert.equal(ensurePyodideCalls, 3);
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
  ensurePyodideImpl = () => new Promise((resolve) => {
    releaseDiscovery = resolve;
  });
  const canceledRun = runner.runAnalysis();
  assert.equal(ensurePyodideCalls, 4);
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
  state.files.c_gb = primary;
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
  let rejectStaleDiscovery;
  ensurePyodideImpl = () => new Promise((_resolve, reject) => {
    rejectStaleDiscovery = reject;
  });
  const staleDiscovery = runner.refreshCircularRecordOrder();
  await Promise.resolve();
  const discoveryStateBeforeRollback = {
    records: structuredClone(state.circularRecordList.value),
    positions: structuredClone(state.adv.multi_record_positions),
    discovery: { ...state.circularRecordDiscovery }
  };
  state.semanticFileWatchersSuppressed.value = true;
  await runner.refreshCircularRecordOrder({ suppress: true });
  rejectStaleDiscovery(new Error('injected restored-file reparse failure'));
  await staleDiscovery;
  assert.deepEqual(state.circularRecordList.value, discoveryStateBeforeRollback.records);
  assert.deepEqual(state.adv.multi_record_positions, discoveryStateBeforeRollback.positions);
  assert.deepEqual(state.circularRecordDiscovery, discoveryStateBeforeRollback.discovery);
  Object.defineProperty(primary, 'text', { value: primaryText, configurable: true });
  state.semanticFileWatchersSuppressed.value = false;
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

  let ensurePyodideCalls = 0;
  let pyodideWrites = 0;
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

  const runner = createRunAnalysis({
    ...generatedArtifactSnapshotOptions,
    state,
    getPyodide: () => null,
    ensurePyodide: async () => {
      ensurePyodideCalls += 1;
      throw new Error('mode none must not initialize Pyodide');
    },
    writeFileToFs: async () => {
      pyodideWrites += 1;
      throw new Error('mode none must not stage LOSAT input through Pyodide');
    },
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
  });

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
  assert.equal(ensurePyodideCalls, 0);
  assert.equal(pyodideWrites, 0);
  assert.equal(losatCalls, 0);
  assert.equal(dormantRawCache.probes, 0);
  assert.equal(dormantDerivedCache.probes, 0);
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
    Object.values(payload.resources).some((resource) => (
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
  assert.equal(state.results.value, committedResults);
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
