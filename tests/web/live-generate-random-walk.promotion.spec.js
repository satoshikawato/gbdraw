// Promotion check, not a PR or dev check (docs/internal/WEB_PERIODIC_AUDIT.md,
// Promotion PR checklist): the seeded random walk of
// tests/web/helpers/live-random-walk.mjs with a longer budget than the PR walk
// (tests/web/live-generate-random-walk.playwright.spec.js).
//
// Run (one browser, fixed seed, default budget):
//   GBDRAW_RANDOM_WALK_SEED=20261006 npx playwright test --config=tests/web/playwright/promotion.config.js --workers=1
// GBDRAW_RANDOM_WALK_SEED (default 20261006) seeds the walk;
// GBDRAW_RANDOM_WALK_STEPS (default 20) is the number of steps per fixture.
// The same seed and budget replay the same walk. The file suffix is
// `.promotion.spec.js` so that no PR or dev Playwright configuration
// (testMatch `.playwright.spec.js`) collects it.
const { test } = require('@playwright/test');
const { randomWalkTests } = require('./helpers/live-random-walk.mjs');

test.describe.configure({ retries: 0 });

const SEED = Number.parseInt(process.env.GBDRAW_RANDOM_WALK_SEED || '20261006', 10);
const STEPS = Number.parseInt(process.env.GBDRAW_RANDOM_WALK_STEPS || '20', 10);
if (!Number.isInteger(SEED) || SEED < 0) throw new Error('GBDRAW_RANDOM_WALK_SEED must be a non-negative integer');
if (!Number.isInteger(STEPS) || STEPS < 1) throw new Error('GBDRAW_RANDOM_WALK_STEPS must be a positive integer');

test.beforeAll(() => {
  // eslint-disable-next-line no-console
  console.log(`live-generate random walk: seed ${SEED}, ${STEPS} steps per fixture`);
});

for (const { title, body } of randomWalkTests({ seed: SEED, steps: STEPS })) test(title, body);
