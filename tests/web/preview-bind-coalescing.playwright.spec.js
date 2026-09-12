const { test, expect } = require('@playwright/test');
const { load, seeds, generate, switchMode } = require('./helpers/mode-transition.cjs');
const { capture, assertCoherent } = require('./helpers/visual-state.cjs');

for (const mode of ['circular', 'linear']) {
  test(`${mode} replacement previews bind once and keep stroke editing after mode switches`, async ({ browser }, info) => {
    test.setTimeout(240000);
    const page = await load(browser, seeds[mode]);
    try {
      await page.evaluate(() => {
        window.__BIND_METRICS__ = {};
        window.__GBDRAW_TEST_HOOKS__.onStructuralMetric = ({ name, value }) => {
          window.__BIND_METRICS__[name] = (window.__BIND_METRICS__[name] || 0) + value;
        };
      });
      for (let run = 1; run <= 2; run++) {
        await page.evaluate(() => { window.__BIND_METRICS__ = {}; });
        await generate(page);
        const metrics = await page.evaluate(() => {
          const count = name => window.__BIND_METRICS__[name] || 0;
          return {
            mounts: count('previewMountAdoptionCount'),
            binds: count('previewBinderInvocationCount'),
            duplicates: count('previewDuplicateBindRejectedCount'),
            accepted: count('previewReadyReceiptAcceptedCount')
          };
        });
        expect.soft(metrics, `Generate ${run} must complete one fresh bind without duplicate requests`)
          .toEqual({ mounts: 1, binds: 1, duplicates: 0, accepted: 1 });
        await assertCoherent(await capture(page, info, `generate-${run}`), `Generate ${run}`, true);
      }
      const before = await page.evaluate(() => window.__GBDRAW_APP__.results[0].content);
      await switchMode(page, mode === 'circular' ? 'linear' : 'circular');
      await switchMode(page, mode);
      expect(await page.evaluate(() => window.__GBDRAW_APP__.results[0].content)).toBe(before);
      await page.locator('summary[aria-label="Features"]').click();
      await page.getByLabel('Block Stroke Width', { exact: true }).fill('3');
      await expect.poll(() => page.evaluate(() => {
        const app = window.__GBDRAW_APP__;
        const stored = new DOMParser().parseFromString(app.results[0].content, 'image/svg+xml');
        const mounted = app.svgContainer.querySelector('svg');
        return [stored, mounted].map(svg => {
          const blocks = [...svg.querySelectorAll('[data-gbdraw-feature-part="block"]')];
          return blocks.length > 0 && blocks.every(el => el.getAttribute('stroke-width') === '3');
        });
      })).toEqual([true, true]);
      await assertCoherent(await capture(page, info, 'stroke-after-mode-return'), 'live stroke edit', true);
      expect(page.externalRequests).toEqual([]);
    } finally { await page.context().close(); }
  });
}
