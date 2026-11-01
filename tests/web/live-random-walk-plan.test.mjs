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
  for (const { fixture, kind, popup } of REPLACED_MATRIX_CASES) {
    const fixtureIndex = FIXTURES.findIndex(({ name }) => name === fixture);
    assert.ok(plan.walks.some((walk) => walk.fixtureIndex === fixtureIndex
      && walk.steps.some((step) => step.kind === kind && step.popup === popup)), `${fixture}: ${kind}`);
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
  const walk = (kinds) => ({ walks: [{ fixtureIndex: 0, steps: kinds.map((kind) => ({ kind, popup: false })) }] });
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

// The targets of a level, built here from the kinds, the fixtures and
// `allowedIn` with loops of their own: a test oracle for `coverageTargets`,
// which the planner and `coverageOf` share, so that narrowing a target group
// fails here instead of passing in both.
const POPUP_STATES = [false, true];
const flagKey = (kind, popup) => `${kind}|${popup ? 'popup' : 'none'}`;
const expectedTargets = (level, required = []) => {
  const pairs = new Set();
  const flags = new Set();
  const fixtureFlags = new Set();
  for (const first of WALK_KINDS) {
    for (const second of WALK_KINDS) {
      if (FIXTURES.some((fixture) => allowedIn(first, fixture) && allowedIn(second, fixture))) pairs.add(`${first}>${second}`);
    }
    for (const popup of POPUP_STATES) {
      flags.add(flagKey(first, popup));
      for (const fixture of FIXTURES) if (allowedIn(first, fixture)) fixtureFlags.add(`${fixture.name}|${flagKey(first, popup)}`);
    }
  }
  const targets = {
    fixtures: new Set(FIXTURES.map(({ name }) => name)),
    sequence: pairs,
    flags,
    required: new Set(required.map(({ fixture, kind, popup }) => `${fixture}|${flagKey(kind, popup)}`))
  };
  return level === 'pr' ? targets : { ...targets, adjacent: new Set(pairs), fixtureFlags };
};

test('the targets of each level are every feasible pair, kind state, and fixture state', () => {
  for (const [level, required] of [['pr', REPLACED_MATRIX_CASES], ['promotion', []]]) {
    const actual = coverageTargets(level, required);
    const expected = expectedTargets(level, required);
    assert.deepEqual(Object.keys(actual).sort(), Object.keys(expected).sort(), level);
    for (const group of Object.keys(expected)) assert.deepEqual([...actual[group]].sort(), [...expected[group]].sort(), `${level} ${group}`);
  }
  const promotion = coverageTargets('promotion');
  // 26 kinds: 26 x 26 ordered pairs less the 4 no fixture runs; 26 kinds x 2
  // popup states; (25 + 25 + 24) kinds that the three fixtures run x 2 popup states.
  assert.deepEqual(Object.fromEntries(Object.entries(promotion).map(([group, keys]) => [group, keys.size])),
    { fixtures: 3, sequence: 672, flags: 52, required: 0, adjacent: 672, fixtureFlags: 148 });
  // A pair with a barrier is a target in both pair groups: a barrier may be A or B.
  for (const barrier of BARRIER_KINDS) {
    for (const kind of ['fill', barrier]) {
      for (const pair of [`${barrier}>${kind}`, `${kind}>${barrier}`]) {
        assert.ok(promotion.sequence.has(pair) && promotion.adjacent.has(pair), pair);
      }
    }
  }
});
