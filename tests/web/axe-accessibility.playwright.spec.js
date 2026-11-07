// axe-core, the devDependency that package.json pins and node_modules holds,
// checks the start page and the page after one Circular and one Linear Generate.
// AXE_VIOLATIONS holds, per state and axe rule, the number of elements that
// violate the rule today. It is the allowlist of the registered rule
// AXE-VIOLATIONS (tools/web-design-rule-guards.json, kind count-map): adding
// or raising an entry is an authority-only change that needs the Owner's
// approval; lowering or removing one is free, and the test names the entries
// that may be lowered.
const { test, expect } = require('@playwright/test');
const { readFileSync } = require('node:fs');
const { openApp } = require('./helpers/app-lifecycle.cjs');
const { seeds, load, generate } = require('./helpers/mode-transition.cjs');

const AXE_SOURCE = readFileSync(require.resolve('axe-core/axe.min.js'), 'utf8');

const AXE_VIOLATIONS = {
  'start color-contrast': 1,
  'start label-content-name-mismatch': 1,
  'start nested-interactive': 13,
  'start region': 2,
  'circular color-contrast': 1,
  'circular label-content-name-mismatch': 4,
  'circular nested-interactive': 13,
  'circular region': 2,
  'linear label-content-name-mismatch': 4,
  'linear nested-interactive': 12,
  'linear region': 2
};

// Runs axe on the page and compares the elements each rule finds with the
// entries of `state` in AXE_VIOLATIONS.
const expectAxeBaseline = async (page, state, testInfo) => {
  // axe reads the page at rest: a transition still running (the Generate bar
  // after a Generate) would show its intermediate colors.
  await page.evaluate(() => document.getAnimations()
    .filter((animation) => Number.isFinite(animation.effect?.getComputedTiming().endTime))
    .forEach((animation) => animation.finish()));
  await page.evaluate(AXE_SOURCE);
  const violations = await page.evaluate(async () => (await window.axe.run(document, { resultTypes: ['violations'] }))
    .violations.map(({ id, help, nodes }) => ({ id, help, count: nodes.length, targets: nodes.slice(0, 3).map(({ target }) => target.join(' ')) })));
  const { compareCountBaseline } = await import('./helpers/count-baseline.mjs');
  const recorded = Object.fromEntries(Object.entries(AXE_VIOLATIONS).filter(([key]) => key.startsWith(`${state} `)));
  const observed = Object.fromEntries(violations.map(({ id, count }) => [`${state} ${id}`, count]));
  const { above, lowerable } = compareCountBaseline(observed, recorded);
  if (lowerable.length) testInfo.annotations.push({ type: 'AXE_VIOLATIONS may be lowered', description: lowerable.join(', ') });
  expect(above.map(({ key, observed: count, recorded: allowed }) => {
    const { help, targets } = violations.find(({ id }) => `${state} ${id}` === key);
    return `${key}: ${count} elements (AXE_VIOLATIONS ${allowed}); ${help}; ${targets.join(' | ')}`;
  }), 'axe violations above AXE_VIOLATIONS').toEqual([]);
};

test('axe finds no violation above AXE_VIOLATIONS on the start page', async ({ page }, testInfo) => {
  await openApp(page);
  await expectAxeBaseline(page, 'start', testInfo);
});

for (const mode of ['circular', 'linear']) {
  test(`axe finds no violation above AXE_VIOLATIONS after a ${mode} Generate`, async ({ browser }, testInfo) => {
    test.setTimeout(240000);
    const page = await load(browser, seeds[mode]);
    try {
      await generate(page);
      await expectAxeBaseline(page, mode, testInfo);
    } finally {
      await page.context().close();
    }
  });
}
