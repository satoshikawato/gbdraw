// The step arguments the random walk draws (`nextStep` in
// tests/web/helpers/live-random-walk.mjs) and the check of a recorded step
// before a replay applies it (`whyNotApplicable`), over a stub state instead of
// a browser.
import assert from 'node:assert/strict';
import test from 'node:test';
import { nextStep, whyNotApplicable } from './helpers/live-random-walk.mjs';
import { createRandom } from './helpers/live-random-walk-plan.mjs';

const LEGEND_ROW_KINDS = ['stroke', 'legendRename', 'legendColor', 'legendDelete'];
const context = { fixture: { mode: 'linear', results: 'single' }, canUndo: false };
const walkState = ({ legend, deletedLegend }) => ({
  drawing: 'linear',
  features: [{ key: 'record 0, feature f38ba7c3f', label: 'CDS|1|300|FL1|', span: '1|300', type: 'CDS', tag: 'FL1', product: '' }],
  resultCount: 1,
  resultIndex: 0,
  visibilityRules: 0,
  colorRules: 0,
  legend,
  deletedLegend,
  palettes: ['default'],
  depth: false,
  undoDepth: 0,
  canUndo: false
});
// Draws of `kind` that name a Legend row (a stroke draws the popup target half the time).
const legendRowSteps = (kind, state) => Array.from({ length: 40 }, (_, seed) => nextStep(createRandom(seed + 1), state, kind, context))
  .filter((step) => step.kind !== 'stroke' || step.target === 'legend');

// OV-2263: seed 40904959235, walk 4 deleted the last listed row of the linear
// drawing and then drew a Legend stroke step with no row.
test('a Legend-row step on an empty Legend restores the deleted rows and names one of them', () => {
  const state = walkState({ legend: [], deletedLegend: ['repeat_region', 'CDS'] });
  for (const kind of LEGEND_ROW_KINDS) {
    const steps = legendRowSteps(kind, state);
    assert.ok(steps.length > 0, kind);
    for (const step of steps) {
      assert.ok(state.deletedLegend.includes(step.row), `${kind}: row ${step.row}`);
      assert.deepEqual(step.prepare, { kind: 'legendRestore' }, kind);
      assert.equal(whyNotApplicable(step, state, context), '', `${kind}: the replay applies it`);
    }
  }
});

test('a Legend-row step on a listed Legend names a listed row and prepares nothing', () => {
  const state = walkState({ legend: ['CDS', 'repeat_region'], deletedLegend: [] });
  for (const kind of LEGEND_ROW_KINDS) {
    for (const step of legendRowSteps(kind, state)) {
      assert.ok(state.legend.includes(step.row), `${kind}: row ${step.row}`);
      assert.equal(step.prepare, undefined, kind);
    }
  }
});

test('a Legend-row step with no row listed or deleted fails with the reason', () => {
  const state = walkState({ legend: [], deletedLegend: [] });
  for (const kind of [...LEGEND_ROW_KINDS, 'legendRestore']) {
    assert.throws(() => legendRowSteps(kind, state), /lists no row and has none deleted/, kind);
  }
});

test('a replay refuses a restoring step when no row is deleted', () => {
  const step = { kind: 'legendColor', prepare: { kind: 'legendRestore' }, row: 'CDS', color: '#123456' };
  assert.equal(whyNotApplicable(step, walkState({ legend: ['CDS'], deletedLegend: [] }), context),
    'its prepare step: no deleted Legend row');
});
