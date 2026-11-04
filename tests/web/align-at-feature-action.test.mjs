// The Align at feature owner (gbdraw/web/js/app/align-at-feature.js) with fake ports:
// the reply check, the preview, Apply over an active Similarity alignment, the
// stale refusal, a Python refusal through the error model, and the notice Undo retires.
import assert from 'node:assert/strict';
import test from 'node:test';

const { createFeatureAnchorAction, validateFeatureAnchorReply } = await import('../../gbdraw/web/js/app/align-at-feature.js');
const { ref } = await import('../../gbdraw/web/js/utils/reactivity.js');

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
  const page = { committed: committedSession(), plan, drawn: [], runs: [], replaced: [], helperCalls: [], sessions: ref(0), steps: [] };
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
      // As `runCommittedCanonicalCandidate` and its History step do: capture the intent before, commit the candidate
      // (the plan is gone, its translations are drawn), commit the intent, capture it after. Undo puts both back.
      const before = { committed: page.committed, plan: page.plan, drawn: page.drawn, intent: run.intent?.capture() };
      page.committed = run.canonical;
      page.plan = false;
      page.drawn = run.canonical.renderRequest.layout.recordTranslations.map((entry) => ({ ...entry }));
      run.intent?.commit();
      const after = { committed: page.committed, plan: page.plan, drawn: page.drawn, intent: run.intent?.capture() };
      /** @param {typeof before} side */
      const install = (side) => {
        Object.assign(page, { committed: side.committed, plan: side.plan, drawn: side.drawn });
        if (run.intent) run.intent.restore(side.intent);
      };
      page.steps.push({ undo: () => install(before), redo: () => install(after) });
      return { status: 'ok' };
    },
    notePlanReplaced: (reason) => page.replaced.push(reason),
    projectPlacement: ({ committed, plan, translations }) => ({ ...committed,
      renderRequest: { ...committed.renderRequest, layout: { similarityAlignment: plan, recordTranslations: translations } } }),
    recordTranslations: () => page.drawn,
    sessionReplacements: () => page.sessions.value
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
    '#1 NC_012920.1 keeps its place (first record with the feature)',
    '#2 NC_002333.2 moves -62.8 px'
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

test('applied holds while the records keep the positions Apply wrote, which Wrap Rows reads as Align at feature', async () => {
  const { page, action } = fixture();
  assert.equal(action.applied.value, false, 'nothing applied yet');
  action.selector.value = 'gene=COX1';
  await action.resolve();
  assert.equal(action.applied.value, false, 'a preview moves nothing');
  await action.apply();
  assert.equal(action.applied.value, true);
  // Undo puts other positions back: Align at feature no longer holds.
  page.drawn = [{ recordKey: 'record-1', x: 0, y: 0 }, { recordKey: 'record-2', x: 12, y: 3 }];
  assert.equal(action.applied.value, false);
  assert.equal(action.notice.value, '');
});

// Review 3 F4 (D-54): a Session Load replaces the document. A Session that holds the positions Apply wrote (the one
// saved after Apply, loaded again in the same page) is the same Session as in a fresh page, where Align at feature
// is not known: Wrap Rows stays enabled and the notice does not come back.
test('a Session Load retires Align at feature even when the loaded positions are the ones Apply wrote', async () => {
  const { page, action } = fixture();
  action.selector.value = 'gene=COX1';
  await action.resolve();
  await action.apply();
  assert.equal(action.applied.value, true);
  const written = page.drawn.map((entry) => ({ ...entry }));
  // Load another Session, then the saved one: the drawn positions are those Apply wrote.
  page.drawn = [];
  page.sessions.value += 1;
  page.drawn = written;
  page.sessions.value += 1;
  assert.equal(action.applied.value, false);
  assert.equal(action.notice.value, '');
  // A new Apply in the loaded Session holds again.
  page.committed = committedSession();
  await action.resolve();
  await action.apply();
  assert.equal(action.applied.value, true);
});

// Review 2 M1 (D-74, OV-2471): the hold is a value of Apply's History step. Its Undo ends it also when Apply moved
// no record (the positions Undo puts back are the ones Apply wrote), and its Redo brings it back.
test('the Undo of an Apply that moved no record ends Align at feature, and its Redo brings it back', async () => {
  const { page, action } = fixture({ helper: async () => ({ result: { ...reply(),
    records: reply().records.map((record) => ({ ...record, status: 'unmatched', afterX: record.beforeX })),
    recordTranslations: [{ recordKey: 'record-1', x: 0, y: 0 }, { recordKey: 'record-2', x: 0, y: 0 }] } }) });
  page.drawn = [{ recordKey: 'record-1', x: 0, y: 0 }, { recordKey: 'record-2', x: 0, y: 0 }];
  action.selector.value = 'gene=COX1';
  await action.resolve();
  assert.deepEqual(await action.apply(), { status: 'ok' });
  assert.equal(action.applied.value, true);
  assert.equal(page.steps.length, 1);
  page.steps[0].undo();
  assert.deepEqual(page.drawn, [{ recordKey: 'record-1', x: 0, y: 0 }, { recordKey: 'record-2', x: 0, y: 0 }]);
  assert.equal(action.applied.value, false, 'the same positions, but Align at feature no longer holds');
  assert.equal(action.notice.value, '');
  page.steps[0].redo();
  assert.equal(action.applied.value, true);
  assert.equal(action.notice.value, 'Records aligned at gene=COX1.');
});

test('a later move of the records retires the hold, and putting the positions back restores it with its notice', async () => {
  const { page, action } = fixture();
  action.selector.value = 'gene=COX1';
  await action.resolve();
  await action.apply();
  const written = page.drawn;
  // Reset record positions writes x 0 and keeps the hold value; the positions alone retire it.
  page.drawn = written.map((entry) => ({ ...entry, x: 0 }));
  assert.equal(action.applied.value, false);
  // Its Undo puts Apply's positions back: the notice is back too.
  page.drawn = written;
  assert.equal(action.applied.value, true);
  assert.equal(action.notice.value, 'Records aligned at gene=COX1.');
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
