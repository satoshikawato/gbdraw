const { test, expect } = require('@playwright/test');
const { readFileSync, writeFileSync } = require('node:fs');
const { join } = require('node:path');
const { openApp } = require('./helpers/app-lifecycle.cjs');

const source = readFileSync(join(__dirname, '../test_inputs/HmmtDNA.gbk'), 'utf8');
const status = page => page.locator('p[role="status"][aria-live="polite"][aria-atomic="true"]');
const generate = page => page.getByRole('button', { name: 'Generate Diagram', exact: true });
const prepare = async (page, baseURL) => {
  const externalRequests = [];
  await page.context().route(/^https?:\/\//, async route => {
    if (new URL(route.request().url()).origin !== new URL(baseURL).origin) {
      externalRequests.push(route.request().url());
      return route.abort('blockedbyclient');
    }
    return route.continue();
  });
  await page.addInitScript(() => {
    window.feedbackEvents = [];
    window.__GBDRAW_TEST_HOOKS__ = {
      onSessionLifecycleEvent(event) {
        window.feedbackEvents.push({ ...event, observedAt: performance.now() });
      }
    };
  });
  await openApp(page);
  await page.evaluate(() => {
    window.feedbackStatuses = [];
    window.Vue.watch(() => window.__GBDRAW_APP__.processingStatus, value => {
      window.feedbackStatuses.push({ value, at: performance.now() });
    }, { flush: 'sync' });
  });
  return externalRequests;
};
const resetObservations = page => page.evaluate(() => {
  window.feedbackStatuses = [];
  window.feedbackEvents = [];
});
const awaitSuccess = async page => {
  await expect.poll(() => page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    return { processing: app.processing, error: app.errorLog?.summary || '', results: app.results.length };
  }), { timeout: 180_000 }).toEqual({ processing: false, error: '', results: 1 });
  await expect(status(page)).toHaveCount(0);
  await expect(generate(page)).toBeEnabled();
  await expect(page.locator('.gbdraw-preview-surface svg').first()).toBeVisible();
};
const record = async (page, testInfo, label) => {
  const evidence = await page.evaluate(() => ({
    statuses: window.feedbackStatuses,
    events: window.feedbackEvents,
    workers: window.__GBDRAW_DIAGRAM_WORKER_ACTIVITY__
  }));
  const path = testInfo.outputPath(`${label}.json`);
  writeFileSync(path, JSON.stringify(evidence, null, 2));
  await testInfo.attach(label, { path, contentType: 'application/json' });
  return evidence;
};
const screenshot = async (page, testInfo, label) => {
  const path = testInfo.outputPath(`${label}.png`);
  await page.screenshot({ path });
  await testInfo.attach(label, { path, contentType: 'image/png' });
};
const expectRenderStages = evidence => {
  const values = evidence.statuses.map(entry => entry.value);
  const expected = [
    'Preparing render inputs and session...',
    'Preparing diagram input resources...',
    'Rendering diagram...',
    'Finalizing diagram results...',
    'Preparing preview...',
    ''
  ];
  let previous = -1;
  for (const value of expected) {
    const index = values.indexOf(value, previous + 1);
    expect(index, JSON.stringify(values)).toBeGreaterThan(previous);
    previous = index;
  }
  const names = evidence.events.map(event => event.name);
  expect(names.indexOf('generate.completed')).toBeGreaterThan(names.indexOf('preview.post-bind-frame-completed'));
  expect(values.every(value => !value.includes('%'))).toBe(true);
};

test('Generate reports real cold/warm stages and preserves the successful Result through cancel, error, and retry', async ({ page, baseURL }, testInfo) => {
  test.setTimeout(180_000);
  const external = await prepare(page, baseURL);
  const setInput = text => page.evaluate(async content => {
    const app = window.__GBDRAW_APP__;
    app.files.c_gb = new File([content], 'HmmtDNA.gbk', { type: 'text/plain' });
    await window.Vue.nextTick();
  }, text);
  await setInput(source);
  await generate(page).click();
  await expect(status(page)).toHaveText('Preparing diagram runtime (first use)...');
  await expect(generate(page)).toBeDisabled();
  await screenshot(page, testInfo, 'cold-runtime-visible');
  await awaitSuccess(page);
  const cold = await record(page, testInfo, 'cold');
  expectRenderStages(cold);
  expect(cold.workers.constructions).toBe(1);

  await resetObservations(page);
  await generate(page).click();
  await awaitSuccess(page);
  const warm = await record(page, testInfo, 'warm');
  expectRenderStages(warm);
  expect(warm.statuses.map(entry => entry.value)).not.toContain('Preparing diagram runtime (first use)...');
  expect(warm.workers.constructions).toBe(1);

  await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    window.successfulResults = app.results;
    window.__GBDRAW_TEST_HOOKS__.beforeDiagramGenerationResponse = () => new Promise(resolve => {
      window.releaseFeedbackResponse = resolve;
    });
  });
  await generate(page).click();
  await page.waitForFunction(() => Boolean(window.releaseFeedbackResponse));
  await expect(status(page)).toHaveText('Finalizing diagram results...');
  await expect(generate(page)).toBeDisabled();
  expect(await page.evaluate(() => window.__GBDRAW_APP__.results === window.successfulResults)).toBe(true);
  await screenshot(page, testInfo, 'finalization-visible');
  await page.getByRole('button', { name: /Cancel$/ }).click();
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.processing)).toBe(false);
  await page.evaluate(() => {
    delete window.__GBDRAW_TEST_HOOKS__.beforeDiagramGenerationResponse;
    window.releaseFeedbackResponse();
  });
  expect(await page.evaluate(() => ({
    status: window.__GBDRAW_APP__.processingStatus,
    sameResult: window.__GBDRAW_APP__.results === window.successfulResults,
    preserved: window.__GBDRAW_APP__.failedGeneratePreservedResult
  }))).toEqual({ status: 'Canceled.', sameResult: true, preserved: true });
  await record(page, testInfo, 'canceled');

  await setInput('invalid GenBank');
  await generate(page).click();
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.processing), { timeout: 180_000 }).toBe(false);
  await expect(page.getByRole('alert', { name: 'Generation Error' })).toBeVisible();
  await expect(generate(page)).toBeEnabled();
  expect(await page.evaluate(() => window.__GBDRAW_APP__.results === window.successfulResults)).toBe(true);
  await record(page, testInfo, 'error');

  await setInput(source);
  await page.setViewportSize({ width: 390, height: 844 });
  await resetObservations(page);
  await generate(page).click();
  await awaitSuccess(page);
  expectRenderStages(await record(page, testInfo, 'retry-mobile'));
  await screenshot(page, testInfo, 'ready-mobile');
  expect(external).toEqual([]);
});

test('Linear comparison preparation, real search counts, cache reuse, and no-comparison generation use the same status', async ({ page, baseURL }, testInfo) => {
  test.setTimeout(180_000);
  const external = await prepare(page, baseURL);
  await page.evaluate(async content => {
    const app = window.__GBDRAW_APP__;
    app.mode = 'linear';
    await window.Vue.nextTick();
    while (app.linearSeqs.length < 2) app.addLinearSeq();
    for (let i = 0; i < 2; i += 1) {
      app.setLinearSeqPrimaryFile(i, 'gb', new File([content], `record-${i}.gbk`, { type: 'text/plain' }));
    }
    app.setLinearComparisonGlobalAction('losat');
    await window.Vue.nextTick();
  }, source);
  await generate(page).click();
  await awaitSuccess(page);
  const search = await record(page, testInfo, 'linear-search');
  expectRenderStages(search);
  expect(search.statuses.map(entry => entry.value)).toEqual(expect.arrayContaining([
    'Preparing comparisons...', 'Preparing LOSAT jobs...',
    'Running LOSAT: 0/1 LOSAT jobs complete', 'Running LOSAT: 1/1 LOSAT jobs complete',
    'Preparing nucleotide comparison results...'
  ]));

  await resetObservations(page);
  await generate(page).click();
  await awaitSuccess(page);
  const cached = await record(page, testInfo, 'linear-cached');
  expectRenderStages(cached);
  expect(cached.statuses.map(entry => entry.value)).toContain('Using cached LOSAT results...');
  expect(cached.statuses.some(entry => entry.value.startsWith('Running LOSAT:'))).toBe(false);

  await page.evaluate(() => window.__GBDRAW_APP__.setLinearComparisonGlobalAction('none'));
  await resetObservations(page);
  await generate(page).click();
  await awaitSuccess(page);
  const none = await record(page, testInfo, 'linear-none');
  expectRenderStages(none);
  expect(none.statuses.some(entry => /comparison|LOSAT/.test(entry.value))).toBe(false);
  expect(external).toEqual([]);
});
