// Reset record positions (gbdraw/web/js/app/record-position-reset.js) with fake ports: one press is one placement
// run (one History step under its own label) that writes x 0 for every record of the committed request, keeps each
// y and record direction, and drops a Similarity alignment plan; the Undo of that step puts the committed Session
// back (the fake History below stands in for the checkpoint of `runCommittedCanonicalCandidate`). A press that
// draws retires Align at feature.
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
  const page = { committed, plan, steps: [], replaced: [], errors: [], projections: [], retired: 0 };
  const action = createRecordPositionReset({
    materializedSession: () => {
      if (!page.committed) throw new Error('Generate a Linear diagram first.');
      return page.committed;
    },
    hasAlignmentPlan: () => page.plan,
    runPlacement: async ({ canonical, label }) => {
      page.steps.push({ label, before: { committed: page.committed, plan: page.plan } });
      page.committed = canonical;
      page.plan = false;
      return { status: 'ok' };
    },
    notePlanReplaced: (reason) => page.replaced.push(reason),
    projectPlacement: (input) => {
      page.projections.push(input);
      return { ...input.committed, renderRequest: { ...input.committed.renderRequest,
        layout: { similarityAlignment: input.plan, recordTranslations: input.translations } } };
    },
    retireFeatureAlignment: () => { page.retired += 1; },
    onError: (error) => page.errors.push(error)
  });
  const undo = () => { const step = page.steps.pop(); page.committed = step.before.committed; page.plan = step.before.plan; };
  return { page, action, undo };
};

test('recordTranslationsAtX0 gives every record of the request x 0 and keeps its y, also a record with no entry', () => {
  assert.deepEqual(recordTranslationsAtX0(committedSession(MOVED).renderRequest), [
    { recordKey: 'record-1', x: 0, y: 0 }, { recordKey: 'record-2', x: 0, y: 7 }, { recordKey: 'record-3', x: 0, y: 0 }]);
  assert.deepEqual(recordTranslationsAtX0({ ...committedSession([]).renderRequest, layout: {} }).map(({ x }) => x), [0, 0, 0]);
});

test('one press is one History step that puts every record at x 0, and its Undo restores the positions', async () => {
  const { page, action, undo } = fixture();
  const before = page.committed;
  assert.deepEqual(await action.reset(), { status: 'ok' });
  assert.equal(page.steps.length, 1);
  assert.equal(page.steps[0].label, 'Reset record positions');
  assert.deepEqual(page.committed.renderRequest.layout.recordTranslations.map(({ x }) => x), [0, 0, 0]);
  assert.deepEqual(page.committed.renderRequest.layout.recordTranslations.map(({ y }) => y), [0, 7, 0]);
  assert.deepEqual(page.projections[0].orientations, [
    { recordKey: 'record-1', reverseComplement: false }, { recordKey: 'record-2', reverseComplement: true },
    { recordKey: 'record-3', reverseComplement: false }]);
  assert.equal(page.projections[0].plan, null);
  assert.deepEqual(page.replaced, []);
  assert.equal(page.retired, 1);
  assert.equal(action.busy.value, false);
  undo();
  assert.equal(page.committed, before);
  assert.deepEqual(page.committed.renderRequest.layout.recordTranslations, MOVED);
});

test('a press over a Similarity alignment drops the plan in the same step and says so', async () => {
  const { page, action } = fixture({ plan: true });
  assert.deepEqual(await action.reset(), { status: 'ok' });
  assert.equal(page.steps.length, 1);
  assert.equal(page.plan, false);
  assert.deepEqual(page.replaced, ['record positions reset.']);
});

test('without a committed Linear Result the press runs nothing and reports the reason', async () => {
  const { page, action } = fixture({ committed: null });
  assert.deepEqual(await action.reset(), { status: 'refused' });
  assert.equal(page.steps.length, 0);
  assert.equal(page.errors.length, 1);
  assert.equal(page.retired, 0);
});
