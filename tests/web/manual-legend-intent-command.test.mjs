import assert from 'node:assert/strict';
import test from 'node:test';

import {
  applyManualLegendIntent,
  buildManualLegendIntentCommand,
  normalizeManualLegendEntries,
  orderedLegendCaptions
} from '../../gbdraw/web/js/app/legend/manual-intent-command.js';
import { createHistoryManager } from '../../gbdraw/web/js/services/history.js';
import { styleFingerprint } from '../../gbdraw/web/js/services/style-revision.js';

const ref = (value) => ({ value });
const clone = (value) => structuredClone(value);

const catalog = { schema: 3, items: [{
  resultKey: 'result-a', resultIndex: 0, resultName: 'a.svg', recordKeys: ['a'],
  biologicalFeatures: [], features: []
}, {
  resultKey: 'result-b', resultIndex: 1, resultName: 'b.svg', recordKeys: ['b'],
  biologicalFeatures: [], features: []
}] };

const fingerprint = styleFingerprint({
  rules: [], appliedPaletteName: 'default', appliedPaletteColors: { CDS: '#cccccc' }
});

const createState = () => ({
  manualSpecificRules: [],
  appliedPaletteName: ref('default'),
  appliedPaletteColors: ref({ CDS: '#cccccc' }),
  manualLegendEntries: ref([{ caption: 'Notes', color: '#111111', owner: 'manual', featureIds: [] }]),
  legendOrderIntent: ref([]),
  legendEntries: ref([{ caption: 'Notes', color: '#111111', owner: 'manual', featureIds: [] }]),
  results: ref([{ name: 'a.svg', content: 'a:Notes:#111111' }, {
    name: 'b.svg', content: 'b:Notes:#111111'
  }]),
  validatedStyleFingerprintByResultKey: ref({
    'result-a': fingerprint, 'result-b': fingerprint
  }),
  semanticStyleRevision: ref(4),
  semanticStyleFingerprint: ref(fingerprint),
  documentEpoch: ref(2),
  resultGenerationKey: ref(3),
  selectedResultIndex: ref(0)
});

const intent = (kind, values = {}) => ({
  kind,
  documentEpoch: 2,
  resultGenerationKey: 3,
  semanticStyleRevision: 4,
  styleFingerprint: fingerprint,
  ...values
});

const projectionFactory = ({ fail = false, commitFailure = false } = {}) => async ({
  results, manualLegendEntries, direction
}) => {
  if (fail) throw new Error('inactive Result projection failed');
  const projectionText = manualLegendEntries.map((entry) => `${entry.caption}:${entry.color}`).join('|');
  const nextResults = results.map((result, index) => ({
    ...result,
    content: `${index ? 'b' : 'a'}:${projectionText}`
  }));
  return {
    nextResults,
    affectedResultKeys: ['result-a', 'result-b'],
    selectedLegendEntries: clone(manualLegendEntries),
    counters: { affectedResults: 2, detachedPasses: 2, changedResults: 2 },
    ...(commitFailure ? { commitFailure: direction } : {})
  };
};

const orderProjectionFactory = ({ fail = false } = {}) => async ({
  results, legendOrderIntent, orderByResult
}) => {
  if (fail) throw new Error('inactive Result order projection failed');
  const keys = ['result-a', 'result-b'];
  const nextResults = results.map((result, index) => {
    const current = result.content.split(':').at(-1).split('|').filter(Boolean);
    const requested = orderByResult?.[keys[index]] || legendOrderIntent;
    return { ...result, content: `${index ? 'b' : 'a'}:${orderedLegendCaptions(current, requested).join('|')}` };
  });
  const selectedCaptions = nextResults[0].content.split(':').at(-1).split('|').filter(Boolean);
  return {
    nextResults,
    affectedResultKeys: keys,
    selectedLegendEntries: selectedCaptions.map((caption) => ({ caption, color: '#111111' })),
    sourceOrderByResult: Object.fromEntries(results.map((result, index) => [
      keys[index], result.content.split(':').at(-1).split('|').filter(Boolean)
    ])),
    counters: {
      affectedResults: 2,
      detachedPasses: 2,
      changedResults: nextResults.filter((result, index) => result.content !== results[index].content).length,
      movedEntries: 4
    }
  };
};

test('manual intent add, rename, color, stroke, and remove are immutable and preserve style fields', () => {
  const original = [{ caption: 'Notes', color: '#111111', showStroke: true }];
  const added = applyManualLegendIntent({
    entries: original,
    intent: { kind: 'add', caption: 'Review', color: '#abcdef' }
  });
  assert.deepEqual(original, [{ caption: 'Notes', color: '#111111', showStroke: true }]);
  assert.equal(added[1].owner, 'manual');
  assert.deepEqual(added[1].featureIds, []);

  const renamed = applyManualLegendIntent({
    entries: added,
    intent: { kind: 'rename', caption: 'Review', newCaption: 'Curated' }
  });
  assert.equal(renamed[1].caption, 'Curated');
  const recolored = applyManualLegendIntent({
    entries: renamed,
    intent: { kind: 'color', caption: 'Curated', color: '#fedcba' }
  });
  assert.equal(recolored[1].color, '#fedcba');
  const stroked = applyManualLegendIntent({
    entries: recolored,
    intent: {
      kind: 'stroke',
      caption: 'Curated',
      value: { kind: 'stroke', strokeColor: '#123456', strokeWidth: 1.25 }
    }
  });
  assert.equal(stroked[1].strokeColor, '#123456');
  assert.equal(stroked[1].strokeWidth, 1.25);
  const inheritedStroke = applyManualLegendIntent({
    entries: stroked,
    intent: { kind: 'stroke', caption: 'Curated', value: { kind: 'inherit' } }
  });
  assert.equal(inheritedStroke[1].strokeColor, null);
  assert.equal(inheritedStroke[1].strokeWidth, null);
  const removed = applyManualLegendIntent({
    entries: inheritedStroke,
    intent: { kind: 'remove', caption: 'Curated' }
  });
  assert.deepEqual(removed, normalizeManualLegendEntries(original));
});

test('manual intents reject duplicate captions and Feature-derived targets', () => {
  assert.throws(() => applyManualLegendIntent({
    entries: [{ caption: 'Notes', color: '#111111' }],
    intent: { kind: 'add', caption: 'notes', color: '#222222' }
  }), /already exists/);
  assert.throws(() => applyManualLegendIntent({
    entries: [{ caption: 'Core', color: '#111111' }],
    selectedEntries: [{ caption: 'Core', color: '#111111', owner: 'feature', featureIds: ['f1'] }],
    intent: { kind: 'color', caption: 'Core', color: '#222222' }
  }), /Feature-derived/);
  assert.throws(() => applyManualLegendIntent({
    entries: [],
    featureCaptionKeys: new Set(['core']),
    intent: { kind: 'add', caption: 'Core', color: '#222222' }
  }), /Feature-derived/);
});

test('one command publishes manual intent, selected view, and every Result atomically', async () => {
  const state = createState();
  const originalResults = state.results.value;
  const command = await buildManualLegendIntentCommand({
    intent: intent('color', { caption: 'Notes', color: '#abcdef' }),
    state,
    catalog,
    prepareProjection: projectionFactory()
  });
  assert.equal(state.results.value, originalResults);
  assert.equal(state.manualLegendEntries.value[0].color, '#111111');
  assert.equal(await command.apply(), true);
  assert.equal(state.manualLegendEntries.value[0].color, '#abcdef');
  assert.equal(state.legendEntries.value[0].color, '#abcdef');
  assert.deepEqual(state.results.value.map((result) => result.content), [
    'a:Notes:#abcdef', 'b:Notes:#abcdef'
  ]);
  assert.equal(state.semanticStyleRevision.value, 5);
  assert.equal(command.counters.affectedResults, 2);
  assert.equal(command.counters.resultArraySwaps, 1);
  assert.ok(command.estimateBytes() > 0);
});

test('inactive Result preparation failure and stale tokens mutate nothing', async () => {
  const state = createState();
  const before = clone({
    manual: state.manualLegendEntries.value,
    legend: state.legendEntries.value,
    results: state.results.value,
    revision: state.semanticStyleRevision.value,
    ledger: state.validatedStyleFingerprintByResultKey.value
  });
  await assert.rejects(() => buildManualLegendIntentCommand({
    intent: intent('remove', { caption: 'Notes' }),
    state,
    catalog,
    prepareProjection: projectionFactory({ fail: true })
  }), /inactive Result/);
  assert.deepEqual({
    manual: state.manualLegendEntries.value,
    legend: state.legendEntries.value,
    results: state.results.value,
    revision: state.semanticStyleRevision.value,
    ledger: state.validatedStyleFingerprintByResultKey.value
  }, before);
  await assert.rejects(() => buildManualLegendIntentCommand({
    intent: intent('remove', { caption: 'Notes', documentEpoch: 99 }),
    state,
    catalog,
    prepareProjection: projectionFactory()
  }), /stale/);
});

test('commit failure restores exact refs, manual intent, Results, revision, and runtime', async () => {
  const state = createState();
  const refs = {
    manual: state.manualLegendEntries.value,
    legend: state.legendEntries.value,
    results: state.results.value,
    ledger: state.validatedStyleFingerprintByResultKey.value
  };
  const runtime = { dirty: false, generation: 8 };
  const command = await buildManualLegendIntentCommand({
    intent: intent('rename', { caption: 'Notes', newCaption: 'Curated' }),
    state,
    catalog,
    prepareProjection: projectionFactory(),
    captureRuntimeState: () => clone(runtime),
    restoreRuntimeState: (snapshot) => {
      Object.assign(runtime, snapshot);
      return true;
    },
    reconcile: () => {
      runtime.dirty = true;
      throw new Error('commit seam failed');
    }
  });
  await assert.rejects(() => command.apply(), /commit seam failed/);
  assert.equal(state.manualLegendEntries.value, refs.manual);
  assert.equal(state.legendEntries.value, refs.legend);
  assert.equal(state.results.value, refs.results);
  assert.equal(state.validatedStyleFingerprintByResultKey.value, refs.ledger);
  assert.equal(state.semanticStyleRevision.value, 4);
  assert.deepEqual(runtime, { dirty: false, generation: 8 });
});

test('History undo and redo reproject all Results as one command', async () => {
  const state = createState();
  const history = createHistoryManager({
    buildIntent: () => ({
      manual: clone(state.manualLegendEntries.value),
      results: state.results.value.map((result) => result.content),
      revision: state.semanticStyleRevision.value
    }),
    applyIntent: () => { throw new Error('command History must not use intent restoration'); },
    buildCheckpoint: () => ({}),
    applyCheckpoint: () => {}
  });
  await history.captureBaseline('initial');
  assert.equal(await history.runUndoableCommand('Add manual legend', () => (
    buildManualLegendIntentCommand({
      intent: intent('add', { caption: 'Review', color: '#abcdef' }),
      state,
      catalog,
      prepareProjection: projectionFactory()
    })
  )), true);
  assert.equal(history.getUndoCount(), 1);
  assert.deepEqual(state.results.value.map((result) => result.content), [
    'a:Notes:#111111|Review:#abcdef', 'b:Notes:#111111|Review:#abcdef'
  ]);
  assert.equal(await history.undo(), true);
  assert.deepEqual(state.results.value.map((result) => result.content), [
    'a:Notes:#111111', 'b:Notes:#111111'
  ]);
  assert.equal(await history.redo(), true);
  assert.deepEqual(state.results.value.map((result) => result.content), [
    'a:Notes:#111111|Review:#abcdef', 'b:Notes:#111111|Review:#abcdef'
  ]);
});

test('History post-apply capture failure compensates the command exactly', async () => {
  const state = createState();
  const refs = {
    manual: state.manualLegendEntries.value,
    legend: state.legendEntries.value,
    results: state.results.value,
    ledger: state.validatedStyleFingerprintByResultKey.value
  };
  let captures = 0;
  const history = createHistoryManager({
    buildIntent: () => {
      captures += 1;
      if (captures > 1) throw new Error('post-apply capture failed');
      return { manual: clone(state.manualLegendEntries.value) };
    },
    applyIntent: () => {},
    buildCheckpoint: () => ({}),
    applyCheckpoint: () => {}
  });
  await history.captureBaseline('initial');
  await assert.rejects(() => history.runUndoableCommand('Remove manual legend', () => (
    buildManualLegendIntentCommand({
      intent: intent('remove', { caption: 'Notes' }),
      state,
      catalog,
      prepareProjection: projectionFactory()
    })
  )), /post-apply capture failed/);
  assert.equal(state.manualLegendEntries.value, refs.manual);
  assert.equal(state.legendEntries.value, refs.legend);
  assert.equal(state.results.value, refs.results);
  assert.equal(state.validatedStyleFingerprintByResultKey.value, refs.ledger);
  assert.equal(state.semanticStyleRevision.value, 4);
  assert.equal(history.getUndoCount(), 0);
  assert.equal(history.getRedoCount(), 0);
});

test('document-global order intent stages heterogeneous Results and undo/redo restores exact order', async () => {
  const state = createState();
  state.results.value = [
    { name: 'a.svg', content: 'a:Alpha|Beta|Only A' },
    { name: 'b.svg', content: 'b:Gamma|Beta' }
  ];
  state.legendEntries.value = [
    { caption: 'Alpha', color: '#111111' },
    { caption: 'Beta', color: '#222222' },
    { caption: 'Only A', color: '#333333' }
  ];
  const command = await buildManualLegendIntentCommand({
    intent: intent('order', { order: ['Beta', 'Alpha', 'Only A'] }),
    state,
    catalog,
    prepareOrderProjection: orderProjectionFactory()
  });

  assert.deepEqual(state.legendOrderIntent.value, []);
  assert.equal(await command.apply(), true);
  assert.deepEqual(state.legendOrderIntent.value, ['Beta', 'Alpha', 'Only A']);
  assert.deepEqual(state.results.value.map((result) => result.content), [
    'a:Beta|Alpha|Only A',
    'b:Beta|Gamma'
  ]);
  assert.deepEqual(state.legendEntries.value.map((entry) => entry.caption), [
    'Beta', 'Alpha', 'Only A'
  ]);
  assert.equal(command.counters.affectedResults, 2);
  assert.equal(command.counters.resultArraySwaps, 1);
  assert.equal(command.counters.legendOrderMoves, 4);

  assert.equal(await command.revert(), true);
  assert.deepEqual(state.legendOrderIntent.value, []);
  assert.deepEqual(state.results.value.map((result) => result.content), [
    'a:Alpha|Beta|Only A',
    'b:Gamma|Beta'
  ]);
  assert.equal(await command.apply({ direction: 'redo' }), true);
  assert.deepEqual(state.results.value.map((result) => result.content), [
    'a:Beta|Alpha|Only A',
    'b:Beta|Gamma'
  ]);
});

test('legend order preflight and commit failures leave order, Results, and revision exact', async () => {
  const state = createState();
  state.results.value = [
    { name: 'a.svg', content: 'a:Alpha|Beta' },
    { name: 'b.svg', content: 'b:Gamma|Beta' }
  ];
  state.legendEntries.value = [
    { caption: 'Alpha', color: '#111111' },
    { caption: 'Beta', color: '#222222' }
  ];
  const before = {
    order: state.legendOrderIntent.value,
    legend: state.legendEntries.value,
    results: state.results.value,
    ledger: state.validatedStyleFingerprintByResultKey.value,
    revision: state.semanticStyleRevision.value
  };
  await assert.rejects(() => buildManualLegendIntentCommand({
    intent: intent('order', { order: ['Beta', 'Alpha'] }),
    state,
    catalog,
    prepareOrderProjection: orderProjectionFactory({ fail: true })
  }), /inactive Result order projection failed/);
  assert.equal(state.legendOrderIntent.value, before.order);
  assert.equal(state.legendEntries.value, before.legend);
  assert.equal(state.results.value, before.results);
  assert.equal(state.semanticStyleRevision.value, before.revision);

  const command = await buildManualLegendIntentCommand({
    intent: intent('order', { order: ['Beta', 'Alpha'] }),
    state,
    catalog,
    prepareOrderProjection: orderProjectionFactory(),
    reconcile: () => { throw new Error('order commit failed'); }
  });
  await assert.rejects(() => command.apply(), /order commit failed/);
  assert.equal(state.legendOrderIntent.value, before.order);
  assert.equal(state.legendEntries.value, before.legend);
  assert.equal(state.results.value, before.results);
  assert.equal(state.validatedStyleFingerprintByResultKey.value, before.ledger);
  assert.equal(state.semanticStyleRevision.value, before.revision);
});
