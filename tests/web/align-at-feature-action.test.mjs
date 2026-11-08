// The Align at feature owner (gbdraw/web/js/app/align-at-feature.js) with fake ports:
// the reply check, the preview, Apply over an active Similarity alignment, the
// starts Apply writes for the records Python rotated, the stale refusal, a Python
// refusal through the error model, and the notice Undo retires.
import assert from 'node:assert/strict';
import test from 'node:test';

const { createFeatureAnchorAction, validateFeatureAnchorReply } = await import('../../gbdraw/web/js/app/align-at-feature.js');
const { ref } = await import('../../gbdraw/web/js/utils/reactivity.js');

const KEYS = ['record-1', 'record-2'];
const RECORDS = KEYS.map((recordKey) => ({ recordKey, cardinality: 'exactly_one' }));
const committedSession = (layout = {}, starts = [null, null]) => ({
  renderRequest: { mode: 'linear', layout,
    records: KEYS.map((recordKey, index) => ({ recordKey, cardinality: 'exactly_one',
      display: { isCircular: null, startCoordinate: starts[index] }, presentation: { reverseComplement: false }, windows: [] })) },
  resources: {}
});
// A linear human record moves below the rotated zebrafish record's left end.
const reply = () => ({
  selector: 'gene=COX1',
  records: [
    { recordKey: 'record-1', recordId: 'NC_012920.1', status: 'aligned', startCoordinate: null, beforeX: 0, afterX: -50.7862 },
    { recordKey: 'record-2', recordId: 'NC_002333.2', status: 'rotated', startCoordinate: 6425, beforeX: 12, afterX: 0 }
  ],
  recordTranslations: [{ recordKey: 'record-1', x: -50.7862, y: 0 }, { recordKey: 'record-2', x: 0, y: 3 }]
});
// Both circular records rotated: nothing moves sideways.
const rotatedReply = () => ({
  selector: 'gene=COX1',
  records: [
    { recordKey: 'record-1', recordId: 'NC_012920.1', status: 'rotated', startCoordinate: 5904, beforeX: 0, afterX: 0 },
    { recordKey: 'record-2', recordId: 'NC_002333.2', status: 'rotated', startCoordinate: 6425, beforeX: 0, afterX: 0 }
  ],
  recordTranslations: [{ recordKey: 'record-1', x: 0, y: 0 }, { recordKey: 'record-2', x: 0, y: 0 }]
});

// A fixture page: one committed Session, and the calls the owner makes. The candidate writes the translations and
// the starts it is given, as `projectCommittedSimilarityAlignment` does.
const fixture = ({ plan = false, helper = async () => ({ result: reply() }), committed = committedSession() } = {}) => {
  const page = { committed, plan, runs: [], projections: [], replaced: [], helperCalls: [], sessions: ref(0), revision: ref(0) };
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
      // As the History step does: the candidate is committed and the plan is gone.
      page.before = page.committed;
      page.committed = run.canonical;
      page.revision.value += 1;
      page.plan = false;
      return { status: 'ok' };
    },
    notePlanReplaced: (reason) => page.replaced.push(reason),
    projectPlacement: (input) => {
      page.projections.push(input);
      const { committed: base, plan: alignment, translations, starts } = input;
      const startOf = new Map(starts.map(({ recordKey, startCoordinate }) => [recordKey, startCoordinate]));
      return { ...base, renderRequest: { ...base.renderRequest,
        records: base.renderRequest.records.map((record) => (startOf.has(record.recordKey)
          ? { ...record, display: { ...record.display, startCoordinate: startOf.get(record.recordKey) } } : record)),
        layout: { similarityAlignment: alignment, recordTranslations: translations } } };
    },
    committedRequest: () => { void page.revision.value; return page.committed?.renderRequest ?? null; },
    sessionReplacements: () => page.sessions.value
  });
  // Undo: the committed Session before the last Apply comes back.
  const undo = () => { page.committed = page.before; page.revision.value += 1; };
  return { page, action, undo };
};

const reasonOf = (error) => error?.context?.reason;

test('validateFeatureAnchorReply admits one entry per record in request order and refuses anything else', () => {
  assert.deepEqual(validateFeatureAnchorReply(reply(), RECORDS), reply());
  assert.deepEqual(validateFeatureAnchorReply(rotatedReply(), RECORDS), rotatedReply());
  const withRecord = (index, patch) => ({ ...reply(), records: reply().records.map((record, at) => (at === index ? { ...record, ...patch } : record)) });
  const refusals = [
    { ...reply(), records: reply().records.slice().reverse() },
    { ...reply(), records: reply().records.slice(0, 1) },
    { ...reply(), recordTranslations: [{ recordKey: 'record-1', x: 0, y: 0 }, { recordKey: 'record-2', x: Number.NaN, y: 0 }] },
    withRecord(0, { status: 'moved' }),
    // A rotated record names its new start, a positive whole base; no other start is read as one.
    withRecord(1, { startCoordinate: null }),
    withRecord(1, { startCoordinate: 0 }),
    withRecord(1, { startCoordinate: 12.5 }),
    withRecord(0, { startCoordinate: '7' }),
    // An exactly-one record is not expanded into members.
    { ...reply(), records: [{ ...reply().records[0], recordKey: 'record-1:1' }, reply().records[1]],
      recordTranslations: [{ recordKey: 'record-1:1', x: 0, y: 0 }, reply().recordTranslations[1]] },
    { ...reply(), selector: 1 },
    null
  ];
  for (const value of refusals) {
    assert.throws(() => validateFeatureAnchorReply(value, RECORDS), (error) => error.code === 'RESULT_INVALID'
      && error.operation === 'resolveFeatureAnchor', JSON.stringify(value));
  }
});

// OV-2474: Python draws a File card that draws several records as one record per member (`<card>:<k>`), and replies
// per drawn record: its keys are checked against the card's request record, not copied from it.
test('validateFeatureAnchorReply admits the members of a File card that draws several records, in order', () => {
  const records = [{ recordKey: 'card', cardinality: 'all' }, { recordKey: 'other', cardinality: 'exactly_one' }];
  const entry = (recordKey, status, startCoordinate) => ({ recordKey, recordId: recordKey, status, startCoordinate, beforeX: 0, afterX: 0 });
  const value = (keys) => ({ selector: 'gene=COX1',
    records: keys.map((key) => entry(key, 'rotated', 10)), recordTranslations: keys.map((recordKey) => ({ recordKey, x: 0, y: 0 })) });
  assert.deepEqual(validateFeatureAnchorReply(value(['card:1', 'card:2', 'other']), records), value(['card:1', 'card:2', 'other']));
  // A card of one record keeps its key.
  assert.deepEqual(validateFeatureAnchorReply(value(['card', 'other']), records), value(['card', 'other']));
  for (const keys of [['card:2', 'card:1', 'other'], ['other', 'card:1'], ['card:1', 'card:1', 'other'], ['card:1', 'card:x', 'other'], ['card:1']]) {
    assert.throws(() => validateFeatureAnchorReply(value(keys), records), (error) => error.code === 'RESULT_INVALID', keys.join());
  }
});

// OV-2474: Apply writes each record's start on its File card's committed row; the members of one card that draws
// several records have no row of their own, so Resolve refuses such a reply with the reason.
test('a reply for the members of a File card that draws several records is refused at Resolve with the reason', async () => {
  const committed = { renderRequest: { mode: 'linear', layout: {},
    records: [{ recordKey: 'card', cardinality: 'all', display: { isCircular: null, startCoordinate: null }, presentation: { reverseComplement: false }, windows: [] }] }, resources: {} };
  const members = { selector: 'gene=COX1',
    records: [1, 2].map((k) => ({ recordKey: `card:${k}`, recordId: `R${k}`, status: 'rotated', startCoordinate: 10 * k, beforeX: 0, afterX: 0 })),
    recordTranslations: [1, 2].map((k) => ({ recordKey: `card:${k}`, x: 0, y: 0 })) };
  const { page, action } = fixture({ committed, helper: async () => ({ result: members }) });
  action.selector.value = 'gene=COX1';
  assert.deepEqual(await action.resolve(), { status: 'error' });
  assert.equal(action.error.value.code, 'INPUT_INVALID');
  assert.equal(reasonOf(action.error.value), 'ANCHOR_RECORD_CARDS');
  assert.match(action.error.value.summary, /Click Generate Diagram, which gives each record its own File card, then use Align at feature again\./);
  assert.equal(action.preview.value, null);
  assert.equal(action.status.value, 'idle');
  assert.equal(page.runs.length, 0);
});

// Align at feature reads the drawn request: new Record Layout rows count after a Generate (OV-2501).
test('a refusal for records that share a row names the Generate that makes new rows count', async () => {
  const refusal = Object.assign(new Error('Worker error'), {
    code: 'INPUT_INVALID', operation: 'resolveFeatureAnchor', stage: 'helper',
    context: { field: 'anchor_feature', reason: 'ANCHOR_SHARED_ROW' }
  });
  const { action } = fixture({ helper: async () => { throw refusal; } });
  action.selector.value = 'gene=COX1';
  assert.deepEqual(await action.resolve(), { status: 'error' });
  assert.equal(reasonOf(action.error.value), 'ANCHOR_SHARED_ROW');
  assert.match(action.error.value.summary,
    /Give each record its own row in Record Layout, click Generate Diagram, then use Align at feature again\./);
});

test('Apply over an active Similarity alignment writes the translations and starts without the plan, in one step, and announces the clearing', async () => {
  const { page, action, undo } = fixture({ plan: true });
  action.selector.value = 'gene=COX1';
  assert.deepEqual(await action.resolve(), { status: 'reviewing' });
  assert.equal(page.helperCalls[0].operation, 'resolveFeatureAnchor');
  assert.equal(page.helperCalls[0].payload.selector, 'gene=COX1');
  assert.equal(action.preview.value.clearsAlignment, true);
  assert.deepEqual(action.previewLines.value, [
    '#1 NC_012920.1 moves -50.8 px',
    '#2 NC_002333.2 starts at base 6425 (rotated)'
  ]);
  assert.equal(page.runs.length, 0, 'nothing changes before Apply');

  assert.deepEqual(await action.apply(), { status: 'ok' });
  assert.equal(page.runs.length, 1);
  assert.equal(page.runs[0].label, 'Align at feature');
  // Only the rotated record's start is written: Python's value, not one the page computes.
  assert.deepEqual(page.runs[0].starts, [{ recordKey: 'record-2', startCoordinate: 6425 }]);
  assert.deepEqual(page.projections[0].starts, page.runs[0].starts);
  assert.deepEqual(page.runs[0].canonical.renderRequest.layout, { similarityAlignment: null, recordTranslations: reply().recordTranslations });
  assert.deepEqual(page.runs[0].canonical.renderRequest.records.map(({ display }) => display.startCoordinate), [null, 6425]);
  assert.deepEqual(page.replaced, ['records aligned at a feature.']);
  assert.equal(action.notice.value, 'Records aligned at gene=COX1.');
  assert.equal(action.preview.value, null);
  assert.equal(action.status.value, 'idle');

  // Undo puts the earlier Session back: the notice goes with it.
  undo();
  assert.equal(action.notice.value, '');
});

test('the notice holds while the records keep the positions Apply wrote', async () => {
  const { page, action, undo } = fixture();
  assert.equal(action.notice.value, '', 'nothing applied yet');
  action.selector.value = 'gene=COX1';
  await action.resolve();
  assert.equal(action.notice.value, '', 'a preview moves nothing');
  await action.apply();
  assert.equal(action.notice.value, 'Records aligned at gene=COX1.');
  // A drag or Reset record positions moves a record: Align at feature no longer holds.
  page.committed = { ...page.committed, renderRequest: { ...page.committed.renderRequest,
    layout: { similarityAlignment: null, recordTranslations: [{ recordKey: 'record-1', x: 0, y: 0 }, { recordKey: 'record-2', x: 0, y: 3 }] } } };
  page.revision.value += 1;
  assert.equal(action.notice.value, '');
  undo();
  assert.equal(action.notice.value, '');
});

// An Apply that only rotated writes x 0 for every record, which the records may already have had: the notice holds
// while the starts Apply wrote are the drawn ones too, so the Undo that puts the earlier starts back retires it.
test('the notice goes after the Undo of an Apply that only rotated', async () => {
  const atZero = { similarityAlignment: null, recordTranslations: rotatedReply().recordTranslations };
  const { page, action, undo } = fixture({ helper: async () => ({ result: rotatedReply() }), committed: committedSession(atZero, [1, null]) });
  action.selector.value = 'gene=COX1';
  await action.resolve();
  assert.deepEqual(action.previewLines.value, ['#1 NC_012920.1 starts at base 5904 (rotated)', '#2 NC_002333.2 starts at base 6425 (rotated)']);
  await action.apply();
  assert.deepEqual(page.runs[0].starts, [{ recordKey: 'record-1', startCoordinate: 5904 }, { recordKey: 'record-2', startCoordinate: 6425 }]);
  assert.equal(action.notice.value, 'Records aligned at gene=COX1.');
  undo();
  assert.deepEqual(page.committed.renderRequest.layout, atZero, 'the positions before Apply are the ones it wrote');
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
  assert.equal(action.notice.value, 'Records aligned at gene=COX1.');
  const written = page.committed;
  // Load another Session, then the saved one: the drawn positions are those Apply wrote.
  page.committed = committedSession();
  page.sessions.value += 1;
  page.committed = written;
  page.sessions.value += 1;
  assert.equal(action.notice.value, '');
  // A new Apply in the loaded Session holds again.
  page.committed = committedSession();
  await action.resolve();
  await action.apply();
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
