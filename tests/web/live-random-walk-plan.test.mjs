// The kind schedule of the live-vs-Generate random walk
// (tests/web/helpers/live-random-walk-plan.mjs): the pull-request plan and the
// promotion plans cover every target of their level within the step budgets.
import assert from 'node:assert/strict';
import test from 'node:test';
import {
  BARRIER_KINDS, DEFAULT_WALK_LENGTH, FIXTURES, PR_SEED, REPLACED_MATRIX_CASES, WALK_KINDS,
  allowedIn, coverageOf, coverageTargets, planWalk, stepCount
} from './helpers/live-random-walk-plan.mjs';

// About 1.1 s per step in CI: about 3 minutes for the pull-request walk and
// 20 minutes for the promotion walk.
const PR_STEP_BUDGET = 150;
const PROMOTION_STEP_BUDGET = 1000;
// The seeds of 100 promotion runs: small numbers and run-ID-sized ones (the
// promotion job seeds the walk with the GitHub run ID).
const PROMOTION_SEEDS = [...Array.from({ length: 94 }, (_, offset) => offset + 1),
  40077991168, 39506600666, 18890421337, 4294967295, 4294967296, 20261006];

test('the pull-request plan covers every pull-request target and holds the replaced matrix cases', () => {
  const plan = planWalk(PR_SEED, { level: 'pr', required: REPLACED_MATRIX_CASES });
  assert.deepEqual(coverageOf(plan, 'pr', { required: REPLACED_MATRIX_CASES }).missing, []);
  for (const { fixture, kind, popup, reflow } of REPLACED_MATRIX_CASES) {
    const fixtureIndex = FIXTURES.findIndex(({ name }) => name === fixture);
    assert.ok(plan.walks.some((walk) => walk.fixtureIndex === fixtureIndex
      && walk.steps.some((step) => step.kind === kind && step.popup === popup && step.reflow === reflow)), `${fixture}: ${kind}`);
  }
  assert.ok(stepCount(plan) <= PR_STEP_BUDGET, `${stepCount(plan)} steps`);
  for (const walk of plan.walks) assert.ok(walk.steps.length <= DEFAULT_WALK_LENGTH.pr);
});

test('every promotion plan covers every promotion target within the step budget', () => {
  for (const seed of PROMOTION_SEEDS) {
    const plan = planWalk(seed, { level: 'promotion' });
    assert.deepEqual(coverageOf(plan, 'promotion').missing, [], `seed ${seed}`);
    assert.ok(stepCount(plan) <= PROMOTION_STEP_BUDGET, `seed ${seed}: ${stepCount(plan)} steps`);
    for (const walk of plan.walks) assert.ok(walk.steps.length <= DEFAULT_WALK_LENGTH.promotion, `seed ${seed}`);
  }
});

test('a plan is the same for the same seed and differs between seeds', () => {
  assert.deepEqual(planWalk(PR_SEED, { level: 'pr', required: REPLACED_MATRIX_CASES }),
    planWalk(PR_SEED, { level: 'pr', required: REPLACED_MATRIX_CASES }));
  assert.deepEqual(planWalk(7, { level: 'promotion' }), planWalk(7, { level: 'promotion' }));
  assert.notDeepEqual(planWalk(7, { level: 'promotion' }), planWalk(8, { level: 'promotion' }));
});

test('every kind is drawn, only in a fixture that allows it', () => {
  for (const plan of [planWalk(PR_SEED, { level: 'pr', required: REPLACED_MATRIX_CASES }), planWalk(3, { level: 'promotion' })]) {
    const drawn = new Set(plan.walks.flatMap(({ steps }) => steps.map(({ kind }) => kind)));
    assert.deepEqual(WALK_KINDS.filter((kind) => !drawn.has(kind)), []);
    for (const { fixtureIndex, steps } of plan.walks) {
      for (const { kind } of steps) assert.ok(allowedIn(kind, FIXTURES[fixtureIndex]), `${kind} in ${FIXTURES[fixtureIndex].name}`);
    }
  }
});

test('a walk length from the environment caps every walk and still covers the targets', () => {
  const plan = planWalk(11, { level: 'promotion', walkLength: 20 });
  assert.ok(plan.walks.every(({ steps }) => steps.length <= 20));
  assert.deepEqual(coverageOf(plan, 'promotion').missing, []);
});

test('coverage counts a pair across a barrier only when the barrier is one of the pair', () => {
  const walk = (kinds) => ({ walks: [{ fixtureIndex: 0, steps: kinds.map((kind) => ({ kind, popup: false, reflow: false })) }] });
  const { missing } = coverageOf(walk(['fill', 'reset', 'stroke']), 'pr');
  assert.ok(missing.includes('sequence fill>stroke'));
  assert.ok(!missing.includes('sequence fill>reset'));
  assert.ok(!missing.includes('sequence reset>stroke'));
  assert.ok(missing.includes('sequence fill>fill'), 'a self pair needs the kind twice');
  assert.ok(!coverageOf(walk(['fill', 'stroke', 'fill']), 'pr').missing.includes('sequence fill>fill'));
  // Adjacency counts only neighbours.
  const promotion = coverageOf(walk(['fill', 'stroke', 'labelText']), 'promotion').missing;
  assert.ok(!promotion.includes('adjacent fill>stroke'));
  assert.ok(promotion.includes('adjacent fill>labelText'));
  assert.ok(!promotion.includes('sequence fill>labelText'));
});

test('the targets leave out the pairs no fixture can run', () => {
  const { sequence } = coverageTargets('pr');
  assert.ok(!sequence.has('resultSwitch>reset'));
  assert.ok(!sequence.has('drawingSwitch>resultSwitch'));
  assert.ok(sequence.has('resultSwitch>load'));
  assert.equal(sequence.size, WALK_KINDS.length ** 2 - 4);
  assert.deepEqual([...BARRIER_KINDS].sort(), ['load', 'reset']);
});
