// The PR random walk (PD-OI-066, gbdraw/web/CLAUDE.md R3): the seeded walk of
// tests/web/helpers/live-random-walk.mjs at the pull-request level of
// tests/web/helpers/live-random-walk-plan.mjs, with a fixed seed: every
// feasible ordered pair of edit kinds, and every kind in each popup and Auto Reflow
// state. Every step compares the live Result with a fresh Generate from the
// same draft and checks that the drawing it did not edit keeps its slices
// (M6). A failing walk is shrunk to a minimal reproduction within
// GBDRAW_RANDOM_WALK_SHRINK_REPLAYS replays (default 10 here), for the first
// GBDRAW_RANDOM_WALK_SHRINK_WALKS failing walks of the run (default 2). The
// promotion spec live-generate-random-walk.promotion.spec.js runs the
// promotion level.
const { expect, test } = require('@playwright/test');
const { ALLOWLIST_WALK_KINDS, randomWalkTests, recordReplayTest, shrinkLimitsFromEnv } = require('./helpers/live-random-walk.mjs');
const { FIXTURES, PR_SEED, REPLACED_MATRIX_CASES, WALK_KINDS, coverageOf, planWalk } = require('./helpers/live-random-walk-plan.mjs');
const { WORK_ALLOWLIST } = require('./helpers/live-work-allowlist.mjs');

test.describe.configure({ retries: 0 });

const PLAN = { level: 'pr', required: REPLACED_MATRIX_CASES };
// At most 2 x 10 replays of at most 40 steps: about 10 minutes more for this
// spec's shard when every walk fails (its job limit is 45 minutes).
const SHRINK = { replays: 10, walks: 2 };

// Each row of the R3 work allowlist names the walk kind that draws its path
// (or null), and the walk draws every kind, every feasible ordered pair of kinds, and
// every kind in each popup and Auto Reflow state on every run.
test('the walk draws every kind and names a kind for each work allowlist row', () => {
  expect(Object.keys(ALLOWLIST_WALK_KINDS)).toEqual(WORK_ALLOWLIST.map(({ kind }) => kind));
  expect(Object.values(ALLOWLIST_WALK_KINDS).filter((kind) => kind !== null && !WALK_KINDS.includes(kind))).toEqual([]);
  const plan = planWalk(PR_SEED, PLAN);
  const drawn = new Set(plan.walks.flatMap(({ steps }) => steps.map(({ kind }) => kind)));
  expect(WALK_KINDS.filter((kind) => !drawn.has(kind))).toEqual([]);
  expect(coverageOf(plan, 'pr', { required: REPLACED_MATRIX_CASES }).missing).toEqual([]);
});

// The case-matrix cases of live-generate-parity.playwright.spec.js that this
// walk replaced (R3, REPLACED_MATRIX_CASES): each is a step of this plan, in
// the fixture, popup state, and Auto Reflow state the case pinned.
test('the walk draws each matrix case it replaced in that case\'s states', () => {
  const { walks } = planWalk(PR_SEED, PLAN);
  const missing = REPLACED_MATRIX_CASES.filter(({ kind, fixture, popup, reflow }) => !walks.some((walk) => FIXTURES[walk.fixtureIndex].name === fixture
    && walk.steps.some((step) => step.kind === kind && step.popup === popup && step.reflow === reflow))).map(({ edit }) => edit);
  expect(missing).toEqual([]);
});

const { title, body } = recordReplayTest({ seed: PR_SEED });
test(title, body);

for (const walk of randomWalkTests({ seed: PR_SEED, ...PLAN, shrink: shrinkLimitsFromEnv(SHRINK) })) test(walk.title, walk.body);
