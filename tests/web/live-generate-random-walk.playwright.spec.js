// The PR random walk (PD-OI-066, gbdraw/web/CLAUDE.md R3): the seeded walk of
// tests/web/helpers/live-random-walk.mjs with a fixed seed and 10 steps per
// fixture. Every step compares the live Result with a fresh Generate from the
// same draft and checks that the drawing it did not edit keeps its slices
// (M6). The promotion spec live-generate-random-walk.promotion.spec.js runs
// the same walk with a longer budget.
const { expect, test } = require('@playwright/test');
const { ALLOWLIST_WALK_KINDS, randomWalkTests, planWalk, WALK_KINDS } = require('./helpers/live-random-walk.mjs');
const { WORK_ALLOWLIST } = require('./helpers/live-work-allowlist.mjs');

test.describe.configure({ retries: 0 });

const SEED = 20261018;
const STEPS = 10;

// Each row of the R3 work allowlist names the walk kind that draws its path
// (or null), and the walk draws every kind at least once per run.
test('the walk draws every kind and names a kind for each work allowlist row', () => {
  expect(Object.keys(ALLOWLIST_WALK_KINDS)).toEqual(WORK_ALLOWLIST.map(({ kind }) => kind));
  expect(Object.values(ALLOWLIST_WALK_KINDS).filter((kind) => kind !== null && !WALK_KINDS.includes(kind))).toEqual([]);
  const drawn = new Set(planWalk(SEED, STEPS).flat().map(({ kind }) => kind));
  expect(WALK_KINDS.filter((kind) => !drawn.has(kind))).toEqual([]);
});

for (const { title, body } of randomWalkTests({ seed: SEED, steps: STEPS })) test(title, body);
