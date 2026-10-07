// E1: each diagram mode keeps its own generated artifact. The History snapshot
// service captures and installs one mode's artifact slot by reference, key by
// key (`ARTIFACT_SLOT_KEYS`). Draft and project state stays in place: the
// specific color rules, the file Legend captions, the Linear record
// orientations, the pending palette, the LOSAT caches and manifests, the
// alignment plan and record translations, and every editor intent.
import assert from 'node:assert/strict';
import { test } from 'node:test';

globalThis.window = { Vue: {
  ref: value => ({ value }), reactive: value => value, toRaw: value => value,
  computed: getter => ({ get value() { return getter(); } }), nextTick: async () => {}
} };
const { state } = await import('../../gbdraw/web/js/state.js');
const { adoptCanonicalRenderArtifacts, canonicalRenderArtifactOwner } = await import('../../gbdraw/web/js/services/config.js');
const { createHistorySnapshotService } = await import('../../gbdraw/web/js/services/history-snapshot.js');
const { ARTIFACT_SLOT_KEYS } = await import('../../gbdraw/web/js/services/artifact-slot.js');
const { createHistoryFileStore } = await import('../../gbdraw/web/js/services/history-files.js');
const { readFile } = await import('node:fs/promises');

const snapshots = createHistorySnapshotService({ state, fileStore: createHistoryFileStore() });
snapshots.setGeneratedArtifactRuntimeOwner({
  capture: () => ({ canonical: canonicalRenderArtifactOwner.capture(), cli: 'circular helper files' }),
  restore: (runtime) => { canonicalRenderArtifactOwner.restore(runtime?.canonical ?? null); }
});

const SWAPPED = [
  'results', 'selectedResultIndex', 'featureCatalog', 'extractedFeatures', 'biologicalFeatures',
  'featureRecordIds', 'orthogroups', 'featureOrthogroupIndex', 'collinearGroups', 'trackSlotResolvedGeometry',
  'annotationWarnings', 'featureIdentityNotices', 'comparisonWarnings', 'lastRunInfo', 'pairwiseMatchFactors',
  'editableLabels', 'generatedLegendPosition', 'generatedMode', 'generatedMultiRecordCanvas',
  'generatedCircularPlotTitlePosition', 'appliedPaletteName', 'appliedPaletteColors',
  'similarityAlignmentResetReceipt', 'originalLegendColors', 'originalSvgStroke'
];
// Draft and project state that sits in today's owner set or beside it.
const KEPT = [
  'pendingPaletteName', 'pendingPaletteColors', 'losatCache', 'losatDerivedCache', 'losatCacheInfo',
  'proteinIdentityManifest', 'legacyProteinRawCandidates', 'legacyProteinDerivedEvidence',
  'similarityAlignmentPlan', 'linearRecordTranslations', 'legendEntries', 'deletedLegendEntries',
  // The Legend owner keeps each Result's inventory; the slot carries it as `legendInventory`.
  'originalLegendOrder'
];

const fillCircularArtifact = async () => {
  const canonical = JSON.parse(await readFile('gbdraw/web/gallery/sessions/HmmtDNA_basic_circular.gbdraw-session.json', 'utf8'));
  adoptCanonicalRenderArtifacts(canonical, { adoptOwnedRequest: true });
  const values = {
    results: [{ name: 'a.svg', content: '<svg id="a"/>' }, { name: 'b.svg', content: '<svg id="b"/>' }],
    selectedResultIndex: 1, featureCatalog: { schema: 5, items: [] }, extractedFeatures: [{ id: 'f1' }],
    biologicalFeatures: [{ id: 'b1' }], featureRecordIds: ['r1'], orthogroups: [{ id: 'og1' }],
    featureOrthogroupIndex: new Map([['f1', 'og1']]), collinearGroups: [{ id: 'c1' }],
    trackSlotResolvedGeometry: { slots: [] }, annotationWarnings: [{ code: 'w' }],
    featureIdentityNotices: [{ status: 'dormant' }], comparisonWarnings: [{ code: 'c' }],
    lastRunInfo: { invocation: { mode: 'circular' } }, pairwiseMatchFactors: { p1: 1.5 },
    editableLabels: [{ key: 'l1' }], generatedLegendPosition: 'right', generatedMode: 'circular',
    generatedMultiRecordCanvas: true, generatedCircularPlotTitlePosition: 'top',
    appliedPaletteName: 'tableau', appliedPaletteColors: { CDS: '#111111' },
    similarityAlignmentResetReceipt: null,
    originalLegendColors: { CDS: '#111111' }, originalSvgStroke: { color: '#000000', width: 1 }
  };
  Object.entries(values).forEach(([key, value]) => { state[key].value = value; });
  const kept = {
    pendingPaletteName: 'pending', pendingPaletteColors: { CDS: '#222222' }, losatCache: new Map([['k', {}]]),
    losatDerivedCache: new Map([['d', {}]]), losatCacheInfo: [{ edgeKey: 'e' }],
    proteinIdentityManifest: { schema: 2, proteinSets: {} }, legacyProteinRawCandidates: { schema: 1, entries: [] },
    legacyProteinDerivedEvidence: { schema: 1, entries: [] }, similarityAlignmentPlan: { groupId: 'g' },
    linearRecordTranslations: [{ recordKey: 'r', x: 1, y: 2 }], legendEntries: [{ caption: 'CDS' }],
    deletedLegendEntries: [{ caption: 'gene' }], originalLegendOrder: ['CDS', 'tRNA']
  };
  Object.entries(kept).forEach(([key, value]) => { state[key].value = value; });
  state.manualSpecificRules.splice(0, state.manualSpecificRules.length, { qual: 'product', val: 'x', color: '#333333', cap: 'X' });
  state.fileLegendCaptions.value = new Set(['From file']);
  state.linearSeqs.splice(0, state.linearSeqs.length, { uid: 'seq-1', region_reverse: true });
  state.matchSequenceRegistry.register({ origin: 'linear-record', sourceKey: 'k1', recordId: 'R1', sequence: 'ACGT' });
  return { values, kept };
};

test('the artifact slot names exactly the displayed-Result part of the generated artifact', () => {
  assert.deepEqual([...ARTIFACT_SLOT_KEYS].sort(), [...SWAPPED].sort());
  KEPT.forEach((key) => assert.ok(!ARTIFACT_SLOT_KEYS.includes(key), `${key} is draft or project state`));
});

test('an empty slot clears only the artifact, and the captured slot comes back by reference', async () => {
  const { values, kept } = await fillCircularArtifact();
  const committed = canonicalRenderArtifactOwner.capture();
  const registryOwner = state.matchSequenceRegistry.captureTrustedOwner();
  const rules = [...state.manualSpecificRules];

  const circular = snapshots.captureArtifactSlot();
  assert.ok(Object.isFrozen(circular));
  assert.deepEqual(circular.legendInventory, ['CDS', 'tRNA']);
  snapshots.installArtifactSlot(null, { mode: 'linear' });

  // The arriving empty Linear slot: no Result, its own mode, no committed Session.
  assert.deepEqual(state.results.value, []);
  assert.equal(state.selectedResultIndex.value, 0);
  assert.equal(state.featureCatalog.value, null);
  assert.equal(state.generatedMode.value, 'linear');
  assert.equal(canonicalRenderArtifactOwner.capture().committedCanonicalSession, null);
  assert.equal(state.matchSequenceRegistry.values().length, 0);
  // An empty slot keeps the applied palette: it has no Result to project.
  assert.equal(state.appliedPaletteName.value, 'tableau');
  assert.strictEqual(state.appliedPaletteColors.value, values.appliedPaletteColors);
  // Draft and project state stays.
  Object.entries(kept).forEach(([key, value]) => assert.strictEqual(state[key].value, value, key));
  assert.deepEqual(state.manualSpecificRules, rules);
  assert.deepEqual([...state.fileLegendCaptions.value], ['From file']);
  assert.equal(state.linearSeqs[0].region_reverse, true);

  const linear = snapshots.captureArtifactSlot();
  snapshots.installArtifactSlot(circular, { mode: 'circular' });
  SWAPPED.forEach((key) => assert.strictEqual(state[key].value, values[key], key));
  assert.strictEqual(canonicalRenderArtifactOwner.capture().committedCanonicalSession, committed.committedCanonicalSession);
  assert.strictEqual(canonicalRenderArtifactOwner.capture().activeSessionResourceTable, committed.activeSessionResourceTable);
  assert.strictEqual(state.matchSequenceRegistry.captureTrustedOwner(), registryOwner);
  Object.entries(kept).forEach(([key, value]) => assert.strictEqual(state[key].value, value, key));

  // The Linear slot captured while empty installs as empty again.
  snapshots.installArtifactSlot(linear, { mode: 'linear' });
  assert.deepEqual(state.results.value, []);
  assert.equal(state.generatedMode.value, 'linear');
});

test('installing a slot of another mode is refused', async () => {
  await fillCircularArtifact();
  const circular = snapshots.captureArtifactSlot();
  assert.throws(() => snapshots.installArtifactSlot(circular, { mode: 'linear' }), /circular artifact.*linear/);
});

// OV-104 (PD-OI-052): decoration continuity throws on a mode change
// (`tests/web/composition-layout.test.mjs`). The app reads the displayed mode's
// Results for it (`legend-layout.js`), so after a switch to a mode without a
// Result a Generate finds no other mode's Result to carry moves from.
test('a Generate after a switch reads no decoration of the other mode\'s Result', async () => {
  const { captureDecorationContinuity } = await import('../../gbdraw/web/js/app/legend-layout/decoration-continuity.js');
  const { values } = await fillCircularArtifact();
  const circular = snapshots.captureArtifactSlot();
  snapshots.installArtifactSlot(null, { mode: 'linear' });
  const capture = () => captureDecorationContinuity({
    canonical: { renderRequest: { mode: 'linear', records: [] }, resources: {} },
    results: state.results.value,
    catalog: state.featureCatalog.value,
    mountedSvg: null,
    selectedResultIndex: state.selectedResultIndex.value,
    projectRecordIdentity: () => assert.fail('no Result of the other mode reaches continuity')
  });
  assert.equal(capture(), null);
  snapshots.installArtifactSlot(circular, { mode: 'circular' });
  assert.strictEqual(state.results.value, values.results);
});

// REVIEW-1 P3: a History step that switches the mode reads its alignment Reset
// receipt against the committed Session of the step's mode before anything is
// applied, so a receipt that does not bind it leaves the shown mode and its
// artifact as they were.
test('a mode-switch step whose alignment receipt does not bind fails before the switch', async () => {
  const { values } = await fillCircularArtifact();
  state.mode.value = 'circular';
  /** @type {string[]} */
  const transitions = [];
  snapshots.registerModeTransition((mode) => { transitions.push(mode); });
  snapshots.registerModeCommittedSession?.((mode) => (
    mode === 'circular' ? canonicalRenderArtifactOwner.capture().committedCanonicalSession : null
  ));
  try {
    await assert.rejects(snapshots.applyHistoryIntent(
      { ui: { mode: 'linear' }, alignmentState: { receipt: { binding: 'stale', directions: [], referenceDeltaX: null }, plan: null } },
      { changes: [{ path: ['ui', 'mode'] }, { path: ['alignmentState', 'receipt'] }] }
    ), /ALIGNMENT_RESET_EVIDENCE|receipt/i);
    assert.deepEqual(transitions, []);
    assert.equal(state.mode.value, 'circular');
    assert.strictEqual(state.results.value, values.results);
  } finally {
    snapshots.registerModeTransition(null);
    snapshots.registerModeCommittedSession?.(null);
  }
});
