// Reset record positions (gbdraw/web/js/app/record-position-reset.js) with fake ports: one press asks for one
// placement run (one History step under its own label, the channel Align at feature's Apply uses) that writes x 0
// for every record of the committed request and keeps each y and record direction. The Undo of that step is the
// run's (`runCommittedCanonicalCandidate`), checked in the browser by tests/web/linear-wrap.playwright.spec.js.
import assert from 'node:assert/strict';
import test from 'node:test';

const { createRecordPositionReset, recordTranslationsAtX0 } = await import('../../gbdraw/web/js/app/record-position-reset.js');

const committedSession = (recordTranslations) => ({
  renderRequest: { mode: 'linear', layout: { similarityAlignment: null, recordTranslations },
    records: ['record-1', 'record-2', 'record-3'].map((recordKey, index) => ({ recordKey,
      presentation: { reverseComplement: index === 1 } })) },
  resources: {}
});
const MOVED = [{ recordKey: 'record-1', x: 40, y: 0 }, { recordKey: 'record-2', x: -12.5, y: 7 }];

const fixture = ({ plan = false, committed = committedSession(MOVED) } = {}) => {
  const page = { committed, plan, steps: [], errors: [], projections: [], replaced: [] };
  const action = createRecordPositionReset({
    materializedSession: () => {
      if (!page.committed) throw new Error('Generate a Linear diagram first.');
      return page.committed;
    },
    hasAlignmentPlan: () => page.plan,
    notePlanReplaced: (reason) => page.replaced.push(reason),
    runPlacement: async (run) => {
      // As the candidate run's History step does: the candidate is committed, with no plan; Undo puts both back.
      page.steps.push({ label: run.label, before: { committed: page.committed, plan: page.plan } });
      page.committed = run.canonical;
      page.plan = run.canonical.renderRequest.layout.similarityAlignment !== null;
      return { status: 'ok' };
    },
    projectPlacement: (input) => {
      page.projections.push(input);
      return { ...input.committed, renderRequest: { ...input.committed.renderRequest,
        layout: { similarityAlignment: input.plan, recordTranslations: input.translations } } };
    },
    onError: (error) => page.errors.push(error)
  });
  return { page, action };
};

test('recordTranslationsAtX0 gives every record of the request x 0 and keeps its y, also a record with no entry', () => {
  assert.deepEqual(recordTranslationsAtX0(committedSession(MOVED).renderRequest), [
    { recordKey: 'record-1', x: 0, y: 0 }, { recordKey: 'record-2', x: 0, y: 7 }, { recordKey: 'record-3', x: 0, y: 0 }]);
  assert.deepEqual(recordTranslationsAtX0({ ...committedSession([]).renderRequest, layout: {} }).map(({ x }) => x), [0, 0, 0]);
});

test('one press asks for one placement run under its own label with every record at x 0, y and directions kept', async () => {
  const { page, action } = fixture();
  const before = page.committed;
  assert.deepEqual(await action.reset(), { status: 'ok' });
  assert.equal(page.steps.length, 1);
  assert.equal(page.steps[0].label, 'Reset record positions');
  assert.equal(page.steps[0].before.committed, before, 'the run starts from the committed Session, which its Undo restores');
  assert.deepEqual(page.committed.renderRequest.layout.recordTranslations.map(({ x }) => x), [0, 0, 0]);
  assert.deepEqual(page.committed.renderRequest.layout.recordTranslations.map(({ y }) => y), [0, 7, 0]);
  assert.deepEqual(page.projections[0].orientations, [
    { recordKey: 'record-1', reverseComplement: false }, { recordKey: 'record-2', reverseComplement: true },
    { recordKey: 'record-3', reverseComplement: false }]);
  assert.equal(page.projections[0].plan, null);
  assert.deepEqual(page.replaced, [], 'no alignment to end');
  assert.equal(action.busy.value, false);
});

// Q5 A (2): a Similarity alignment's positions are moved positions; the press ends the alignment in its one step.
test('a press under a Similarity alignment drops the plan in the same step and says so', async () => {
  const { page, action } = fixture({ plan: true });
  assert.deepEqual(await action.reset(), { status: 'ok' });
  assert.equal(page.steps.length, 1);
  assert.equal(page.steps[0].before.plan, true, 'the step that its Undo reverts held the plan');
  assert.equal(page.plan, false);
  assert.deepEqual(page.committed.renderRequest.layout.recordTranslations.map(({ x }) => x), [0, 0, 0]);
  assert.deepEqual(page.replaced, ['record positions reset.']);
});

test('without a committed Linear Result the press runs nothing and reports the reason', async () => {
  const { page, action } = fixture({ committed: null });
  assert.deepEqual(await action.reset(), { status: 'refused' });
  assert.equal(page.steps.length, 0);
  assert.equal(page.errors.length, 1);
});
