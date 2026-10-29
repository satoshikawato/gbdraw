// Promotion check, not a PR or dev check (docs/internal/WEB_PERIODIC_AUDIT.md,
// Promotion PR checklist): the seeded random walk of
// tests/web/helpers/live-random-walk.mjs at the promotion level of
// tests/web/helpers/live-random-walk-plan.mjs (every ordered pair of edit
// kinds as adjacent steps and within a walk, every fixture x kind x popup x
// Auto Reflow; about 800 steps in about 17 walks).
//
// Run (one browser, fixed seed):
//   GBDRAW_RANDOM_WALK_SEED=20261006 npx playwright test --config=tests/web/playwright/promotion.config.js --workers=1
// GBDRAW_RANDOM_WALK_SEED (default 20261006) seeds the plan and the steps;
// GBDRAW_RANDOM_WALK_STEPS (default 50) caps the steps of each walk; the same
// seed and cap replay the same walks, and `-g "walk 3 of"` replays one.
// GBDRAW_RANDOM_WALK_SHRINK_REPLAYS (default 40, 0 disables) bounds the
// replays that shrink a failing walk to a minimal reproduction.
// GBDRAW_RANDOM_WALK_REPLAY=<walk-steps.json or walk-shrunk.json> replays that
// recorded run instead of planning. The file suffix is `.promotion.spec.js`
// so that no PR or dev Playwright configuration (testMatch
// `.playwright.spec.js`) collects it.
const { test } = require('@playwright/test');
const { randomWalkTests, replayWalkTest, shrinkReplaysFromEnv } = require('./helpers/live-random-walk.mjs');
const { DEFAULT_WALK_LENGTH } = require('./helpers/live-random-walk-plan.mjs');

test.describe.configure({ retries: 0 });

const SEED = Number.parseInt(process.env.GBDRAW_RANDOM_WALK_SEED || '20261006', 10);
const STEPS = Number.parseInt(process.env.GBDRAW_RANDOM_WALK_STEPS || `${DEFAULT_WALK_LENGTH.promotion}`, 10);
const REPLAY = process.env.GBDRAW_RANDOM_WALK_REPLAY || '';
const SHRINK_REPLAYS = shrinkReplaysFromEnv(40);
if (!Number.isInteger(SEED) || SEED < 0) throw new Error('GBDRAW_RANDOM_WALK_SEED must be a non-negative integer');
if (!Number.isInteger(STEPS) || STEPS < 2) throw new Error('GBDRAW_RANDOM_WALK_STEPS must be an integer of at least 2');

test.beforeAll(() => {
  // eslint-disable-next-line no-console
  console.log(REPLAY ? `live-generate random walk: replay of ${REPLAY}`
    : `live-generate random walk: seed ${SEED}, at most ${STEPS} steps per walk`);
});

const tests = REPLAY ? [replayWalkTest(REPLAY, { shrinkReplays: SHRINK_REPLAYS })]
  : randomWalkTests({ seed: SEED, level: 'promotion', walkLength: STEPS, shrinkReplays: SHRINK_REPLAYS });
for (const { title, body } of tests) test(title, body);
