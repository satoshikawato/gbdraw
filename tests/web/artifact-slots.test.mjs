// E1: each diagram mode keeps its own generated artifact. The History snapshot
// service captures and installs one mode's artifact slot by reference, key by
// key (`ARTIFACT_SLOT_KEYS`). Draft and project state stays in place: the
// specific color rules, the file Legend captions, the Linear record
// orientations, the pending palette, the LOSAT caches and manifests, the
// alignment plan and record translations, and every editor intent.
import assert from 'node:assert/strict';
import { test } from 'node:test';
import { createHash } from 'node:crypto';
import { savedDrawingView } from './helpers/saved-drawing-view.mjs';
import { modeDrawing } from './helpers/drawing-state.mjs';

globalThis.window = { Vue: {
  ref: value => ({ value }), reactive: value => value, toRaw: value => value,
  computed: getter => ({ get value() { return getter(); } }), nextTick: async () => {}
} };
const { state } = await import('../../gbdraw/web/js/state.js');
await (await import('./helpers/drawing-state.mjs')).addStartLinearDrawing();
const { adoptCanonicalRenderArtifacts, canonicalRenderArtifactOwner, installShownDrawing } = await import('../../gbdraw/web/js/services/config.js');
const { createHistorySnapshotService } = await import('../../gbdraw/web/js/services/history-snapshot.js');
const { ARTIFACT_SLOT_KEYS, createArtifactSlotStore } = await import('../../gbdraw/web/js/services/artifact-slot.js');
const { createInputFileStore } = await import('../../gbdraw/web/js/services/history-files.js');
const { readFile } = await import('node:fs/promises');

const snapshots = createHistorySnapshotService({ state, fileStore: createInputFileStore(), currentResults: () => state.results.value });
snapshots.setGeneratedArtifactRuntimeOwner({
  capture: () => ({ canonical: canonicalRenderArtifactOwner.capture(), cli: 'circular helper files' }),
  restore: (runtime) => { canonicalRenderArtifactOwner.restore(runtime?.canonical ?? null); }
});

const SWAPPED = [
  'results', 'selectedResultIndex', 'featureCatalog', 'extractedFeatures', 'featureEditorStatus', 'biologicalFeatures',
  'featureRecordIds', 'orthogroups', 'featureOrthogroupIndex', 'collinearGroups', 'trackSlotResolvedGeometry',
  'annotationWarnings', 'featureIdentityNotices', 'comparisonWarnings', 'labelFitReports', 'inputSummaries', 'lastRunInfo',
  'pairwiseMatchFactors',
  'editableLabels', 'generatedLegendPosition', 'generatedMultiRecordCanvas',
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

// The kept members of the app state and of the shown mode's drawing.
const memberOf = (key) => (Object.hasOwn(state, key) ? state : state.activeDrawing())[key];

const fillCircularArtifact = async () => {
  const saved = JSON.parse(await readFile('gbdraw/web/gallery/sessions/HmmtDNA_basic_circular.gbdraw-session.json', 'utf8'));
  const canonical = savedDrawingView(saved, saved.ui.activeDrawingId);
  adoptCanonicalRenderArtifacts(canonical, { adoptOwnedRequest: true, conservationPairs: new Set(['["s","q"]']) });
  const values = {
    results: [{ name: 'a.svg', content: '<svg id="a"/>' }, { name: 'b.svg', content: '<svg id="b"/>' }],
    selectedResultIndex: 1, featureCatalog: { schema: 5, items: [] }, extractedFeatures: [{ id: 'f1' }],
    featureEditorStatus: { status: 'failed', error: 'recovery failed' },
    biologicalFeatures: [{ id: 'b1' }], featureRecordIds: ['r1'], orthogroups: [{ id: 'og1' }],
    featureOrthogroupIndex: new Map([['f1', 'og1']]), collinearGroups: [{ id: 'c1' }],
    trackSlotResolvedGeometry: { slots: [] }, annotationWarnings: [{ code: 'w' }],
    featureIdentityNotices: [{ status: 'dormant' }], comparisonWarnings: [{ code: 'c' }],
    labelFitReports: [{ resultIndex: 0, resultName: 'a.svg', report: { drawn: true }, rowToken: null }],
    inputSummaries: [{ resultIndex: 0, resultName: 'a.svg', lines: ['Record r1 (a.gb), record 1 of 1, 8 bp: drew 0 features.'], warnings: [] }],
    lastRunInfo: { invocation: { mode: 'circular' } }, pairwiseMatchFactors: { p1: 1.5 },
    editableLabels: [{ key: 'l1' }], generatedLegendPosition: 'right',
    generatedMultiRecordCanvas: true, generatedCircularPlotTitlePosition: 'top',
    appliedPaletteName: 'tableau', appliedPaletteColors: { CDS: '#111111' },
    similarityAlignmentResetReceipt: null,
    originalLegendColors: { CDS: '#111111' }, originalSvgStroke: { color: '#000000', width: 1 }
  };
  Object.entries(values).forEach(([key, value]) => { state[key].value = value; });
  state.generatedDrawingId.value = modeDrawing(state, 'circular').id; // the displayed artifact is Circular's
  const kept = {
    pendingPaletteName: 'pending', pendingPaletteColors: { CDS: '#222222' }, losatCache: new Map([['k', {}]]),
    losatDerivedCache: new Map([['d', {}]]), losatCacheInfo: [{ edgeKey: 'e' }],
    proteinIdentityManifest: { schema: 2, proteinSets: {} }, legacyProteinRawCandidates: { schema: 1, entries: [] },
    legacyProteinDerivedEvidence: { schema: 1, entries: [] }, similarityAlignmentPlan: { groupId: 'g' },
    linearRecordTranslations: [{ recordKey: 'r', x: 1, y: 2 }], legendEntries: [{ caption: 'CDS' }],
    deletedLegendEntries: [{ caption: 'gene' }], originalLegendOrder: ['CDS', 'tRNA']
  };
  Object.entries(kept).forEach(([key, value]) => { memberOf(key).value = value; });
  state.activeDrawing().manualSpecificRules.splice(0, state.activeDrawing().manualSpecificRules.length, { qual: 'product', val: 'x', color: '#333333', cap: 'X' });
  modeDrawing(state, 'linear').linearSeqs.splice(0, modeDrawing(state, 'linear').linearSeqs.length, { uid: 'seq-1', region_reverse: true });
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
  const rules = [...state.activeDrawing().manualSpecificRules];

  const circular = snapshots.captureArtifactSlot();
  assert.ok(Object.isFrozen(circular));
  assert.deepEqual(circular.legendInventory, ['CDS', 'tRNA']);
  snapshots.installArtifactSlot(null, { drawing: modeDrawing(state, 'linear') });

  // The arriving empty Linear slot: no Result, its own mode, no committed Session.
  assert.deepEqual(state.results.value, []);
  assert.equal(state.selectedResultIndex.value, 0);
  assert.equal(state.featureCatalog.value, null);
  assert.equal(state.generatedDrawingId.value, 'linear');
  assert.equal(canonicalRenderArtifactOwner.capture().committedCanonicalSession, null);
  assert.equal(canonicalRenderArtifactOwner.capture().conservationPairs.size, 0);
  assert.equal(state.matchSequenceRegistry.values().length, 0);
  // An empty slot keeps the applied palette: it has no Result to project.
  assert.equal(state.appliedPaletteName.value, 'tableau');
  assert.strictEqual(state.appliedPaletteColors.value, values.appliedPaletteColors);
  // Draft and project state stays.
  Object.entries(kept).forEach(([key, value]) => assert.strictEqual(memberOf(key).value, value, key));
  assert.deepEqual(state.activeDrawing().manualSpecificRules, rules);
  assert.equal(modeDrawing(state, 'linear').linearSeqs[0].region_reverse, true);

  // The displayed artifact is the Linear drawing's now (its mode is linear; the test installs it without a switch).
  const linear = snapshots.captureArtifactSlot();
  snapshots.installArtifactSlot(circular, { drawing: modeDrawing(state, 'circular') });
  SWAPPED.forEach((key) => assert.strictEqual(state[key].value, values[key], key));
  assert.strictEqual(canonicalRenderArtifactOwner.capture().committedCanonicalSession, committed.committedCanonicalSession);
  assert.strictEqual(canonicalRenderArtifactOwner.capture().activeSessionResourceTable, committed.activeSessionResourceTable);
  // D-15-124: the request's conservation hash pairs wait in its slot with it.
  assert.strictEqual(canonicalRenderArtifactOwner.capture().conservationPairs, committed.conservationPairs);
  assert.deepEqual([...committed.conservationPairs], ['["s","q"]']);
  assert.strictEqual(state.matchSequenceRegistry.captureTrustedOwner(), registryOwner);
  Object.entries(kept).forEach(([key, value]) => assert.strictEqual(memberOf(key).value, value, key));

  // The Linear slot captured while empty installs as empty again.
  snapshots.installArtifactSlot(linear, { drawing: modeDrawing(state, 'linear') });
  assert.deepEqual(state.results.value, []);
  assert.equal(state.generatedDrawingId.value, 'linear');
});

test('installing a slot of another drawing is refused', async () => {
  await fillCircularArtifact();
  const circular = snapshots.captureArtifactSlot();
  assert.equal(circular.drawingId, state.activeDrawing().id, 'the slot is of the drawing it was captured for');
  assert.throws(() => snapshots.installArtifactSlot(circular, { drawing: modeDrawing(state, 'linear') }),
    (error) => error.code === 'RESULT_INVALID' && error.context.reason === 'DRAWING_MISMATCH'
      && error.context.artifactDrawingId === 'circular' && error.context.drawingId === 'linear');
});

// OV-104 (PD-OI-052): decoration continuity throws on a mode change
// (`tests/web/composition-layout.test.mjs`). The app reads the displayed mode's
// Results for it (`legend-layout.js`), so after a switch to a mode without a
// Result a Generate finds no other mode's Result to carry moves from.
test('a Generate after a switch reads no decoration of the other mode\'s Result', async () => {
  const { captureDecorationContinuity } = await import('../../gbdraw/web/js/app/legend-layout/decoration-continuity.js');
  const { values } = await fillCircularArtifact();
  const circular = snapshots.captureArtifactSlot();
  snapshots.installArtifactSlot(null, { drawing: modeDrawing(state, 'linear') });
  const capture = () => captureDecorationContinuity({
    canonical: { renderRequest: { mode: 'linear', records: [] }, resources: {} },
    results: state.results.value,
    catalog: state.featureCatalog.value,
    mountedSvg: null,
    selectedResultIndex: state.selectedResultIndex.value,
    projectRecordIdentity: () => assert.fail('no Result of the other mode reaches continuity')
  });
  assert.equal(capture(), null);
  snapshots.installArtifactSlot(circular, { drawing: modeDrawing(state, 'circular') });
  assert.strictEqual(state.results.value, values.results);
});

// REVIEW-1 P3: a History step reads its alignment Reset receipt against the
// committed Session of its drawing's artifact before anything is applied, so a
// receipt that does not bind it leaves the drawing and its artifact as they
// were. A step is always of the shown drawing (PD-OI-090), and a step of a
// drawing that is not shown is refused.
test('a step whose alignment receipt does not bind fails before anything is applied', async () => {
  const { values } = await fillCircularArtifact();
  installShownDrawing('circular');
  const prefix = state.activeDrawing().form.prefix;
  await assert.rejects(snapshots.drawingHistory('circular').applyIntent(
    {
      config: { form: { ...state.activeDrawing().form, prefix: 'not applied' } },
      alignmentState: { receipt: { binding: 'stale', directions: [], referenceDeltaX: null }, plan: null }
    },
    { changes: [{ path: ['config', 'form', 'prefix'] }, { path: ['alignmentState', 'receipt'] }] }
  ), /ALIGNMENT_RESET_EVIDENCE|receipt/i);
  assert.equal(state.activeDrawing().form.prefix, prefix);
  assert.strictEqual(state.results.value, values.results);
});

// Review alignment compares its mode's own Result whichever mode is shown (OV-622, D-45).
test('a value of a drawing\'s artifact is read from the shown one, or from the slot that keeps it', async () => {
  const { createArtifactSlot } = await import('../../gbdraw/web/js/services/artifact-slot.js');
  const linearResult = { name: 'linear.svg' };
  const otherLinearResult = { name: 'other.svg' };
  const circularResult = { name: 'circular.svg' };
  const kept = createArtifactSlot({ drawing: modeDrawing(state, 'linear'), values: { results: [otherLinearResult, linearResult], selectedResultIndex: 1 } });
  assert.equal(kept.drawingId, 'linear');
  assert.equal(Object.hasOwn(kept.values, 'generatedMode'), false, 'the slot names its drawing by ID only');
  const slots = createArtifactSlotStore();
  slots.replaceLoaded([kept]);
  /** @param {'circular' | 'linear'} id @param {'circular' | 'linear'} shownId @param {unknown} shown */
  const selected = (id, shownId, shown) => slots.readOfDrawing(id, shownId, () => shown,
    (slot) => slot.values.results[slot.values.selectedResultIndex]);
  assert.equal(selected('linear', 'linear', linearResult), linearResult);
  assert.equal(selected('linear', 'circular', circularResult), linearResult);
  assert.equal(selected('circular', 'linear', linearResult), undefined);
});

// REVIEW-1 L1: a drawing with no slot captures a null applied palette (the
// empty slot's "keep the displayed palette"). Restored while that drawing is
// shown, the null keeps the displayed palette, as it does while it is hidden;
// it does not fall back to the drawing's draft palette.
test('a null applied palette restored in the shown drawing keeps the displayed palette', async () => {
  const { applyUiStateData } = await import('../../gbdraw/web/js/services/config.js');
  const service = createHistorySnapshotService({ state, fileStore: createInputFileStore(), currentResults: () => state.results.value, applyUiStateData });
  installShownDrawing('circular');
  const drawing = modeDrawing(state, 'circular');
  const draftPalette = drawing.selectedPalette.value;
  drawing.selectedPalette.value = 'draft-palette';
  state.appliedPaletteName.value = 'displayed-palette';
  try {
    await service.drawingHistory('circular').applyIntent(
      { ui: { appliedPaletteName: null, appliedPaletteColors: null } },
      { changes: [{ path: ['ui', 'appliedPaletteName'] }] }
    );
    assert.equal(state.appliedPaletteName.value, 'displayed-palette');
    await service.drawingHistory('circular').applyIntent(
      { ui: { appliedPaletteName: 'kept-palette', appliedPaletteColors: null } },
      { changes: [{ path: ['ui', 'appliedPaletteName'] }] }
    );
    assert.equal(state.appliedPaletteName.value, 'kept-palette', 'a captured palette is restored');
  } finally {
    drawing.selectedPalette.value = draftPalette;
  }
});

// REVIEW-1 M1: one store keeps the hidden drawings' artifacts, and only its
// named operations change it; an entry's key is its slot's drawing.
test('the slot store keeps a hidden drawing\'s artifact through its named operations', async () => {
  const { createArtifactSlot } = await import('../../gbdraw/web/js/services/artifact-slot.js');
  const circular = modeDrawing(state, 'circular');
  const linear = modeDrawing(state, 'linear');
  const store = createArtifactSlotStore();
  const circularSlot = createArtifactSlot({ drawing: circular, values: { results: [{ name: 'c.svg' }] } });
  assert.equal(store.swap(circularSlot, 'linear'), null, 'a drawing with no artifact yet has no slot');
  assert.deepEqual(store.kept(), [circularSlot]);
  const linearSlot = createArtifactSlot({ drawing: linear, values: { results: [{ name: 'l.svg' }] } });
  assert.strictEqual(store.swap(linearSlot, 'circular'), circularSlot, 'the arriving slot, by reference');
  assert.deepEqual(store.kept(), [linearSlot], 'the shown drawing has no slot');
  const before = store.snapshot();
  store.writeValues(linear, { appliedPaletteName: 'Pastel' });
  const written = store.kept()[0];
  assert.equal(written.drawingId, 'linear');
  assert.equal(written.values.appliedPaletteName, 'Pastel');
  assert.strictEqual(written.values.results, linearSlot.values.results, 'the other values stay by reference');
  store.writeValues(circular, { appliedPaletteName: 'Vivid' });
  assert.deepEqual(store.kept().map((slot) => slot.drawingId), ['linear', 'circular'], 'an empty slot takes the values');
  store.replaceLoaded([]);
  assert.deepEqual(store.kept(), []);
  store.restore(before);
  assert.deepEqual(store.kept(), [linearSlot]);
  assert.ok(Object.isFrozen(store.kept()));
});

// OV-111: a label rerender installs through the artifact owner of the drawing
// whose Result it draws again, which may differ from the drawing shown (a Load
// can show Circular with the Linear artifact displayed); the artifact of
// another drawing is refused and keeps its values.
test('a rerender lands in the artifact of its Result\'s drawing', async () => {
  await fillCircularArtifact();
  const before = { results: state.results.value, featureCatalog: state.featureCatalog.value };
  const redrawn = {
    featureCatalog: { schema: 5, items: [{ id: 'redrawn' }] }, extractedFeatures: [{ id: 'f2' }],
    biologicalFeatures: [{ id: 'b2' }], results: [{ name: 'a.svg', content: '<svg id="a2"/>' }],
    labelFitReports: [],
    // The input summary the rerender drew replaces the Generate's (it draws the live Feature visibility).
    inputSummaries: [{ resultIndex: 0, resultName: 'a.svg', lines: ['redrawn'], warnings: [], drawnVisibility: '[[],[]]' }]
  };
  /** @type {unknown[][]} */
  const installed = [];
  const installResults = (/** @type {unknown[]} */ results) => { installed.push(results); state.results.value = results; };
  assert.throws(() => snapshots.installRedrawnArtifact('linear', redrawn, { installResults }),
    (error) => error.code === 'RESULT_INVALID' && error.context.reason === 'DRAWING_MISMATCH'
      && error.context.drawingId === 'linear' && error.context.artifactDrawingId === 'circular');
  assert.deepEqual(installed, []);
  assert.strictEqual(state.results.value, before.results);
  assert.strictEqual(state.featureCatalog.value, before.featureCatalog);

  snapshots.installRedrawnArtifact('circular', redrawn, { installResults });
  assert.deepEqual(installed, [redrawn.results], 'the Results commit through the caller');
  for (const key of ['featureCatalog', 'extractedFeatures', 'biologicalFeatures', 'labelFitReports', 'inputSummaries']) {
    assert.strictEqual(state[key].value, redrawn[key], key);
  }
});

// U7 (R1): each capture of the Results (an artifact slot, a generated
// artifact's owner set, a History checkpoint) reads them through the preview
// owner's `currentResults`, which writes the displayed Result's edits first.
test('the slot, owner-set and checkpoint captures read the Results through currentResults', async () => {
  const written = [{ name: 'shown.svg', content: '<svg id="written"></svg>' }];
  let reads = 0;
  const service = createHistorySnapshotService({
    state, fileStore: createInputFileStore(), currentResults: () => { reads += 1; return written; }
  });
  assert.equal(service.captureArtifactSlot().values.results, written);
  assert.equal(service.captureGeneratedArtifactOwnerSet(state.activeDrawing()).results, written);
  const checkpoint = await service.drawingHistory(state.activeDrawingId.value).buildCheckpoint();
  assert.deepEqual(checkpoint.results, [{ name: 'shown.svg', content: '<svg id="written"></svg>' }]);
  assert.equal(reads, 3);
});

// OV-2210: the project's LOSAT cache-info rows sit in every drawing's owner set,
// and the transport identity names them by reference. A pass over them that
// changes no row (another drawing's Linear file change filters its pairs out
// of an empty list) keeps the array, so the shown drawing's Result keeps its
// fingerprint and an unchanged Generate there records no step.
test('a pass over the LOSAT rows that changes none keeps the shown artifact\'s transport identity', async () => {
  const { writeLosatCacheInfo } = await import('../../gbdraw/web/js/services/config.js');
  await fillCircularArtifact();
  state.losatCacheInfo.value = [{ key: 'k', filename: 'pair.tsv', display: true }];
  snapshots.setGeneratedArtifactIdentity({ schema: 1, algorithm: 'SHA-256', fingerprint: 'b'.repeat(64), retainedBytes: 10 },
    { results: state.results.value });
  const fingerprint = () => snapshots.captureGeneratedArtifactHandle(state.activeDrawing()).identity.fingerprint;
  assert.equal(fingerprint(), 'b'.repeat(64));
  const rows = state.losatCacheInfo.value;
  writeLosatCacheInfo(rows.filter((entry) => !entry?.edgeKey));
  assert.strictEqual(state.losatCacheInfo.value, rows, 'no row changed: the array stays');
  assert.equal(fingerprint(), 'b'.repeat(64), 'the Result keeps its identity');
  writeLosatCacheInfo([]);
  assert.deepEqual(state.losatCacheInfo.value, [], 'a pass that drops a row writes');
  // OV-2212: the rows are project state, not the drawing's Result.
  assert.equal(fingerprint(), 'b'.repeat(64), 'changed project rows leave the Result its identity');
});

// OV-2212: a drawing's transport identity matches only what its Result is: its
// artifact slot values. Each Generate installs its own project LOSAT and protein
// values (new objects), so a Generate in another drawing replaces them; the
// Result's fingerprint covers what of them it drew, and the compact signature
// covers the alignment plan and record translations by value. So an unchanged
// Generate after a Generate in another drawing is the same artifact (no step),
// and a Generate whose Result differs is still a step.
test('an unchanged Generate in a drawing records no step after a Generate in another drawing', async () => {
  const { createHistoryManager } = await import('../../gbdraw/web/js/services/history.js');
  await fillCircularArtifact();
  installShownDrawing('circular');
  const histories = Object.fromEntries(['circular', 'linear'].map((id) => [id, createHistoryManager({
    ...snapshots.drawingHistory(id),
    captureGeneratedArtifactHandle: () => snapshots.captureGeneratedArtifactHandle(modeDrawing(state, id)),
    restoreGeneratedArtifactHandle: snapshots.restoreGeneratedArtifactHandle,
    compareGeneratedArtifactHandles: snapshots.compareGeneratedArtifactHandles
  })]));
  const slots = new Map();
  const baselines = new Set();
  /** @param {'circular' | 'linear'} id */
  const show = async (id) => {
    slots.set(state.generatedDrawingId.value, snapshots.captureArtifactSlot());
    installShownDrawing(id);
    snapshots.installArtifactSlot(slots.get(id) ?? null, { drawing: modeDrawing(state, id) });
    if (!baselines.has(id)) await histories[id].initializeIntentBaseline();
    baselines.add(id);
  };
  /** @param {'circular' | 'linear'} id @param {string} content */
  const generate = (id, content) => histories[id].runUndoableArtifactReplacement('Generate diagram', () => {
    const drawing = modeDrawing(state, id);
    const results = [{ name: `${id}.svg`, content }];
    snapshots.installGeneratedArtifactOwnerSet({
      ...snapshots.captureGeneratedArtifactOwnerSet(drawing),
      results,
      losatCache: new Map(state.losatCache.value),
      losatDerivedCache: new Map(state.losatDerivedCache.value),
      losatCacheInfo: [...state.losatCacheInfo.value],
      proteinIdentityManifest: { ...state.proteinIdentityManifest.value }
    }, { drawing, installResults: (next) => { state.results.value = next; } });
    snapshots.setGeneratedArtifactIdentity({
      schema: 1, algorithm: 'SHA-256', fingerprint: createHash('sha256').update(content).digest('hex'), retainedBytes: 10
    }, { results });
  });

  await show('circular');
  await generate('circular', '<svg id="c1"/>');
  assert.equal(histories.circular.getUndoCount(), 1);
  await show('linear');
  await generate('linear', '<svg id="l1"/>');
  assert.equal(histories.linear.getUndoCount(), 1);
  await show('circular');
  await generate('circular', '<svg id="c1"/>');
  assert.equal(histories.circular.getUndoCount(), 1, 'the unchanged Generate is the same artifact');
  await generate('circular', '<svg id="c2"/>');
  assert.equal(histories.circular.getUndoCount(), 2, 'a changed Result is a step');
  assert.equal(histories.linear.getUndoCount(), 1);
});
