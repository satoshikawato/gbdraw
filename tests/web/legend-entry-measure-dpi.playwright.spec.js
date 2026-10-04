const { test, expect } = require('@playwright/test');
const { load, generate } = require('./helpers/mode-transition.cjs');
const { evaluateWithRetainedPromise } = require('./helpers/app-lifecycle.cjs');

// Python lays the legend out at the DPI its config resolves to (96 by default);
// Add entry must measure the new caption at that DPI, not at 72. The expected
// widths are gbdraw.core.text.calculate_bbox_dimensions('Added legend entry',
// 'Arial', 20, dpi): the HmmtDNA circular legend draws 20 px captions.
const CAPTION = 'Added legend entry';
const WIDTH_AT_RENDER_DPI = 179.44010416666669;
const WIDTH_AT_72_DPI = 134.580078125;

test('Add entry measures the caption at the committed Result render DPI', async ({ browser }) => {
  test.setTimeout(600_000);
  const page = await load(browser);
  page.setDefaultTimeout(180_000);
  try {
    await page.evaluate(() => {
      window.__MEASURE_EXCHANGES__ = [];
      const original = Worker.prototype.postMessage;
      Worker.prototype.postMessage = function (message, ...rest) {
        if (message?.type === 'helper' && message.operation === 'measureLegendText') {
          this.addEventListener('message', (event) => {
            if (event.data?.result?.width !== undefined) {
              window.__MEASURE_EXCHANGES__.push({ payload: message.payload, width: event.data.result.width });
            }
          });
        }
        return original.call(this, message, ...rest);
      };
    });
    await generate(page);
    const committedOptions = await page.evaluate(async () => {
      const { getCommittedCanonicalRenderRequest } = await import('./js/services/config.js');
      const options = getCommittedCanonicalRenderRequest().diagramOptions;
      return JSON.parse(JSON.stringify({ config: options.config ?? null, configOverrides: options.configOverrides ?? {} }));
    });
    await evaluateWithRetainedPromise(page, async (caption) => {
      const app = window.__GBDRAW_APP__;
      app.newLegendCaption = caption;
      app.newLegendColor = '#884422';
      await app.addNewLegendEntry();
    }, CAPTION);
    await expect.poll(() => page.evaluate(() => window.__MEASURE_EXCHANGES__.length)).toBe(1);
    const [exchange] = await page.evaluate(() => window.__MEASURE_EXCHANGES__);
    expect(exchange.payload.caption).toBe(CAPTION);
    expect(exchange.payload.fontSize).toBe(20);
    expect({ config: exchange.payload.config, configOverrides: exchange.payload.configOverrides }).toEqual(committedOptions);
    expect(Math.abs(exchange.width - WIDTH_AT_RENDER_DPI)).toBeLessThan(0.5);
    expect(Math.abs(exchange.width - WIDTH_AT_72_DPI)).toBeGreaterThan(1);
  } finally {
    await page.context().close();
  }
});
