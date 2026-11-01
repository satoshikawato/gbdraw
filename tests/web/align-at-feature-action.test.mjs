// The Align at feature owner (gbdraw/web/js/app/feature-anchor.js) with fake ports:
// the reply check, the preview, Apply over an active Similarity alignment, the
// stale refusal, a Python refusal through the error model, and the notice Undo retires.
import assert from 'node:assert/strict';
import test from 'node:test';

const { createFeatureAnchorAction, validateFeatureAnchorReply } = await import('../../gbdraw/web/js/app/feature-anchor.js');

const KEYS = ['record-1', 'record-2'];
const committedSession = (layout = null) => ({
  renderRequest: { mode: 'linear', layout,
    records: KEYS.map((recordKey) => ({ recordKey, presentation: { reverseComplement: false } })) },
  resources: {}
});
const reply = () => ({
  selector: 'gene=COX1',
  records: [
    { recordKey: 'record-1', recordId: 'NC_012920.1', status: 'reference', beforeX: 0, afterX: 0 },
    { recordKey: 'record-2', recordId: 'NC_002333.2', status: 'aligned', beforeX: 12, afterX: -50.7862 }
  ],
  recordTranslations: [{ recordKey: 'record-1', x: 0, y: 0 }, { recordKey: 'record-2', x: -50.7862, y: 3 }]
});

// A fixture page: one committed Session, the drawn translations, and the calls the owner makes.
const fixture = ({ plan = false, helper = async () => ({ result: reply() }) } = {}) => {
  const page = { committed: committedSession(), plan, drawn: [], runs: [], replaced: [], helperCalls: [] };
  const action = createFeatureAnchorAction({
    materializedSession: () => {
      if (!page.committed) throw new Error('Generate a Linear diagram before aligning.');
      return page.committed;
    },
    hasAlignmentPlan: () => page.plan,
    runHelperOperation: async (operation, payload) => { page.helperCalls.push({ operation, payload }); return helper(); },
    resolveOperation: 'resolveFeatureAnchor',
    runPlacement: async (run) => {
      page.runs.push(run);
      // As the History step does: the candidate is committed, the plan is gone, its translations are drawn.
      page.committed = run.canonical;
      page.plan = false;
      page.drawn = run.canonical.renderRequest.layout.recordTranslations.map((entry) => ({ ...entry }));
      return { status: 'ok' };
    },
    notePlanReplaced: (reason) => page.replaced.push(reason),
    projectPlacement: ({ committed, plan, translations }) => ({ ...committed,
      renderRequest: { ...committed.renderRequest, layout: { similarityAlignment: plan, recordTranslations: translations } } }),
    recordTranslations: () => page.drawn
  });
  return { page, action };
};

const reasonOf = (error) => error?.context?.reason;

test('validateFeatureAnchorReply admits one entry per record in request order and refuses anything else', () => {
  assert.deepEqual(validateFeatureAnchorReply(reply(), KEYS), reply());
  const refusals = [
    { ...reply(), records: reply().records.slice().reverse() },
    { ...reply(), records: reply().records.slice(0, 1) },
    { ...reply(), recordTranslations: [{ recordKey: 'record-1', x: 0, y: 0 }, { recordKey: 'record-2', x: Number.NaN, y: 0 }] },
    { ...reply(), records: [{ ...reply().records[0], status: 'moved' }, reply().records[1]] },
    { ...reply(), selector: 1 },
    null
  ];
  for (const value of refusals) {
    assert.throws(() => validateFeatureAnchorReply(value, KEYS), (error) => error.code === 'RESULT_INVALID'
      && error.operation === 'resolveFeatureAnchor', JSON.stringify(value));
  }
});

test('Apply over an active Similarity alignment writes the translations without the plan, in one step, and announces the clearing', async () => {
  const { page, action } = fixture({ plan: true });
  action.selector.value = 'gene=COX1';
  assert.deepEqual(await action.resolve(), { status: 'reviewing' });
  assert.equal(page.helperCalls[0].operation, 'resolveFeatureAnchor');
  assert.equal(page.helperCalls[0].payload.selector, 'gene=COX1');
  assert.equal(action.preview.value.clearsAlignment, true);
  assert.deepEqual(action.previewLines.value, [
    '#1 NC_012920.1: keeps its place (first record with the feature)',
    '#2 NC_002333.2: moves -62.8 px'
  ]);
  assert.equal(page.runs.length, 0, 'nothing changes before Apply');

  assert.deepEqual(await action.apply(), { status: 'ok' });
  assert.equal(page.runs.length, 1);
  assert.equal(page.runs[0].label, 'Align at feature');
  assert.deepEqual(page.runs[0].canonical.renderRequest.layout, { similarityAlignment: null, recordTranslations: reply().recordTranslations });
  assert.deepEqual(page.replaced, ['records aligned at a feature.']);
  assert.equal(action.notice.value, 'Records aligned at gene=COX1.');
  assert.equal(action.preview.value, null);
  assert.equal(action.status.value, 'idle');

  // Undo puts the earlier translations back: the notice goes with them.
  page.drawn = [{ recordKey: 'record-1', x: 0, y: 0 }, { recordKey: 'record-2', x: 12, y: 3 }];
  assert.equal(action.notice.value, '');
});

test('Apply without an active plan clears nothing', async () => {
  const { page, action } = fixture();
  action.selector.value = 'gene=COX1';
  await action.resolve();
  assert.equal(action.preview.value.clearsAlignment, false);
  assert.deepEqual(await action.apply(), { status: 'ok' });
  assert.deepEqual(page.replaced, []);
});

test('Apply after the Result changed, or is gone, is refused with the reason and returns to idle', async () => {
  for (const change of [
    (page) => { page.committed = committedSession({ recordTranslations: [{ recordKey: 'record-2', x: 9, y: 0 }] }); },
    (page) => { page.committed = null; }
  ]) {
    const { page, action } = fixture();
    action.selector.value = 'gene=COX1';
    await action.resolve();
    change(page);
    assert.deepEqual(await action.apply(), { status: 'refused' });
    assert.equal(page.runs.length, 0);
    assert.equal(action.error.value.code, 'INPUT_INVALID');
    assert.equal(reasonOf(action.error.value), 'ANCHOR_STALE');
    assert.match(action.error.value.summary, /The diagram changed after the preview\. Align at feature again\./);
    assert.equal(action.status.value, 'idle');
    assert.equal(action.preview.value, null);
  }
});

test('a Python refusal reaches the page through the error model, with its field and reason', async () => {
  const refusal = Object.assign(new Error('Worker error'), {
    code: 'INPUT_INVALID', operation: 'resolveFeatureAnchor', stage: 'helper',
    context: { field: 'anchor_feature', reason: 'SELECTOR_NO_MATCH' }
  });
  const { page, action } = fixture({ helper: async () => { throw refusal; } });
  action.selector.value = 'gene=nonexistent';
  assert.deepEqual(await action.resolve(), { status: 'error' });
  assert.equal(action.status.value, 'idle');
  assert.equal(action.preview.value, null);
  assert.equal(action.error.value.operation, 'resolveFeatureAnchor');
  assert.equal(reasonOf(action.error.value), 'SELECTOR_NO_MATCH');
  assert.match(action.error.value.summary, /No record has a feature that this selector names\./);
  assert.equal(page.runs.length, 0);
});
