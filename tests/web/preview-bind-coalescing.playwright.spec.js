const { test, expect } = require('@playwright/test');
const { load, seeds, generate, switchMode } = require('./helpers/mode-transition.cjs');
const { capture, assertCoherent } = require('./helpers/visual-state.cjs');

for (const mode of ['circular', 'linear']) {
  test(`${mode} replacement previews bind once and apply stroke edits on Generate after mode switches`, async ({ browser }, info) => {
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
      await page.locator('summary[aria-label="Features"]').press('Enter');
      // PD-OI-060 (D-05): the global block width applies on Generate.
      await page.getByLabel('Block Stroke Width', { exact: true }).fill('3');
      await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.adv.block_stroke_width)).toBe(3);
      expect(await page.evaluate(() => window.__GBDRAW_APP__.results[0].content)).toBe(before);
      await generate(page);
      expect(await page.evaluate(() => {
        const app = window.__GBDRAW_APP__;
        const stored = new DOMParser().parseFromString(app.results[0].content, 'image/svg+xml');
        const mounted = app.svgContainer.querySelector('svg');
        return [stored, mounted].map(svg => {
          const blocks = [...svg.querySelectorAll('[data-gbdraw-feature-part="block"]')];
          return blocks.length > 0 && blocks.every(el => Number(el.getAttribute('stroke-width')) === 3);
        });
      })).toEqual([true, true]);
      await assertCoherent(await capture(page, info, 'stroke-after-mode-return'), 'generated stroke edit', true);
      expect(page.externalRequests).toEqual([]);
    } finally { await page.context().close(); }
  });
}

// OV-1030: once a script has read `labels` of the page's controls (Playwright's
// `getByLabel` does), Chromium updates every such live list for each node the
// HTML parser inserts into the page document: a large Result took seconds to
// mount instead of tens of milliseconds. The preview parses a Result outside
// the page document and inserts the parsed tree once.
test('a Generate parses the preview Result outside the page document (OV-1030)', async ({ browser }) => {
  test.setTimeout(240000);
  const page = await load(browser, seeds.linear);
  try {
    await generate(page);
    await page.evaluate(() => {
      for (const control of document.querySelectorAll('input, select, textarea, button')) void control.labels?.length;
      window.__PAGE_DOCUMENT_SVG_PARSES__ = 0;
      // Count every way to parse SVG markup into the page document.
      const count = () => { window.__PAGE_DOCUMENT_SVG_PARSES__ += 1; };
      const hasSvg = markup => String(markup).includes('<svg');
      const inPage = element => element.ownerDocument === document && !(element instanceof HTMLTemplateElement);
      const wrapSetter = name => {
        const property = /** @type {PropertyDescriptor} */ (Object.getOwnPropertyDescriptor(Element.prototype, name));
        Object.defineProperty(Element.prototype, name, {
          ...property,
          set(markup) {
            if (inPage(this) && hasSvg(markup)) count();
            property.set.call(this, markup);
          }
        });
      };
      wrapSetter('innerHTML');
      wrapSetter('outerHTML');
      const wrapMethod = (owner, name, counts) => {
        const original = owner[name];
        if (typeof original !== 'function') return;
        owner[name] = function (...args) {
          const result = original.apply(this, args);
          if (counts(this, args, result)) count();
          return result;
        };
      };
      wrapMethod(Element.prototype, 'insertAdjacentHTML', (el, [, markup]) => inPage(el) && hasSvg(markup));
      wrapMethod(Element.prototype, 'setHTMLUnsafe', (el, [markup]) => inPage(el) && hasSvg(markup));
      wrapMethod(Range.prototype, 'createContextualFragment', (range, [markup], fragment) =>
        fragment.ownerDocument === document && hasSvg(markup));
      wrapMethod(Document.prototype, 'importNode', (doc, [node], result) =>
        doc === document && node.ownerDocument !== document && result.ownerDocument === document
        && (node.localName === 'svg' || Boolean(node.querySelector?.('svg'))));
    });
    await generate(page);
    const facts = await page.evaluate(() => {
      const app = window.__GBDRAW_APP__;
      const parses = window.__PAGE_DOCUMENT_SVG_PARSES__;
      // What `innerHTML` makes of the Result, parsed in another document.
      const reference = document.implementation.createHTMLDocument('').createElement('div');
      reference.innerHTML = app.results[0].content;
      const mounted = app.svgContainer;
      return {
        parses,
        elements: mounted.getElementsByTagName('*').length === reference.getElementsByTagName('*').length,
        rootInPage: mounted.firstElementChild?.ownerDocument === document && mounted.isConnected
      };
    });
    expect(facts).toEqual({ parses: 0, elements: true, rootInPage: true });
    expect(page.externalRequests).toEqual([]);
  } finally { await page.context().close(); }
});
