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

const fixture = ({ committed = committedSession(MOVED) } = {}) => {
  const page = { committed, steps: [], errors: [], projections: [] };
  const action = createRecordPositionReset({
    materializedSession: () => {
      if (!page.committed) throw new Error('Generate a Linear diagram first.');
      return page.committed;
    },
    runPlacement: async (run) => {
      page.steps.push({ label: run.label, intent: run.intent, before: page.committed });
      page.committed = run.canonical;
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

test('one press asks for one placement run under its own label with every record at x 0, y and directions kept, and no other value', async () => {
  const { page, action } = fixture();
  const before = page.committed;
  assert.deepEqual(await action.reset(), { status: 'ok' });
  assert.equal(page.steps.length, 1);
  assert.equal(page.steps[0].label, 'Reset record positions');
  assert.equal(page.steps[0].before, before, 'the run starts from the committed Session, which its Undo restores');
  assert.equal(page.steps[0].intent, undefined, 'no value besides the positions (D-74: no plan drop, no retire)');
  assert.deepEqual(page.committed.renderRequest.layout.recordTranslations.map(({ x }) => x), [0, 0, 0]);
  assert.deepEqual(page.committed.renderRequest.layout.recordTranslations.map(({ y }) => y), [0, 7, 0]);
  assert.deepEqual(page.projections[0].orientations, [
    { recordKey: 'record-1', reverseComplement: false }, { recordKey: 'record-2', reverseComplement: true },
    { recordKey: 'record-3', reverseComplement: false }]);
  assert.equal(page.projections[0].plan, null);
  assert.equal(action.busy.value, false);
});

test('without a committed Linear Result the press runs nothing and reports the reason', async () => {
  const { page, action } = fixture({ committed: null });
  assert.deepEqual(await action.reset(), { status: 'refused' });
  assert.equal(page.steps.length, 0);
  assert.equal(page.errors.length, 1);
});
