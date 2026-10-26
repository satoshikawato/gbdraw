const { test, expect } = require('@playwright/test');
const { readFileSync } = require('node:fs');
const {
  assertDiagramWorkerIdle,
  assertSessionLoadLeftWorkerIdle,
  assertWorkerReuseAcrossHelperAndRender,
  openApp,
  evaluateWithRetainedPromise
} = require('./helpers/app-lifecycle.cjs');

test('uncached protein LOSAT helpers and render share one lazy Worker runtime', async ({ page }) => {
  test.setTimeout(240000);
  page.on('dialog', (dialog) => dialog.accept());
  await page.addInitScript(() => {
    window.__GBDRAW_LOSAT_EXECUTOR_CALLS__ = 0;
    window.__GBDRAW_LOSAT_EXECUTOR__ = async (jobs, options) => {
      window.__GBDRAW_LOSAT_EXECUTOR_CALLS__ += 1;
      const firstFastaId = (text) => String(text || '').match(/^>([^\s]+)/m)?.[1] || '';
      return jobs.map((job) => ({
        cacheKey: job.cacheKey,
        text: [
          firstFastaId(options.sequences.get(job.querySequenceKey)),
          firstFastaId(options.sequences.get(job.subjectSequenceKey)),
          '100', '30', '0', '0', '1', '30', '1', '30', '1e-20', '100'
        ].join('\t') + '\n'
      }));
    };
  });
  await openApp(page);
  await page.evaluate(() => new Promise((resolve) => requestAnimationFrame(() => {
    requestAnimationFrame(resolve);
  })));
  await assertDiagramWorkerIdle(page);
  expect(await page.evaluate(() => ({
    mainLoaderPresent: typeof window.loadPyodide === 'function',
    mainRuntimeStatePresent: Object.prototype.hasOwnProperty.call(
      window.__GBDRAW_APP__,
      'pyodideReady'
    )
  }))).toEqual({ mainLoaderPresent: false, mainRuntimeStatePresent: false });

  const imported = await evaluateWithRetainedPromise(page, async () => {
    const response = await fetch('/gbdraw/web/gallery/sessions/BGC0000708-BGC0000713.gbdraw-session.json');
    const file = new File(
      [await response.text()],
      'BGC0000708-BGC0000713.gbdraw-session.json',
      { type: 'application/json' }
    );
    return window.__GBDRAW_APP__.importSession({
      target: { files: [file], value: '' }
    });
  });

  expect(imported?.status).toBe('ok');
  await assertSessionLoadLeftWorkerIdle(page);

  const proteinMode = await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    const { state } = await import('/gbdraw/web/js/state.js');
    state.losatCache.value.clear();
    state.losatDerivedCache.value.clear();
    await app.setLinearComparisonGlobalAction('losat');
    app.setLinearComparisonLosatMode('blastp');
    app.setLinearComparisonLosatpMode('pairwise');
    return {
      program: app.losatProgram,
      blastpMode: app.losat.blastp.mode,
      cacheEntries: state.losatCache.value.size,
      activeAlignmentPlan: state.similarityAlignmentPlan.value
    };
  });
  expect(proteinMode).toEqual({
    program: 'blastp',
    blastpMode: 'pairwise',
    cacheEntries: 0,
    activeAlignmentPlan: null
  });

  const generated = await evaluateWithRetainedPromise(page, async () => {
    const app = window.__GBDRAW_APP__;
    const result = await app.runAnalysis();
    return {
      result,
      errorLog: app.errorLog,
      resultCount: app.results.length,
      executorCalls: window.__GBDRAW_LOSAT_EXECUTOR_CALLS__,
      cacheMisses: app.lastRunInfo?.losatTelemetry?.cacheMisses
    };
  });

  expect(generated.result).toEqual({ status: 'ok' });
  expect(generated.errorLog).toBeNull();
  expect(generated.resultCount).toBeGreaterThan(0);
  expect(generated.executorCalls).toBe(1);
  expect(generated.cacheMisses).toBeGreaterThan(0);

  const hydrationDownloadPromise = page.waitForEvent('download', { timeout: 120000 });
  await page.evaluate(() => window.__GBDRAW_APP__.downloadLosatPair(0));
  const hydrationDownload = await hydrationDownloadPromise;
  const hydratedTsv = readFileSync(await hydrationDownload.path(), 'utf8');
  expect(hydratedTsv.trim().split('\t')).toHaveLength(12);
  expect(hydratedTsv).not.toMatch(/\bh_[a-z2-7]{26}\b/);

  const repeated = await evaluateWithRetainedPromise(page, async () => ({
    result: await window.__GBDRAW_APP__.runAnalysis(),
    errorLog: window.__GBDRAW_APP__.errorLog,
    executorCalls: window.__GBDRAW_LOSAT_EXECUTOR_CALLS__
  }));
  expect(repeated).toEqual({
    result: { status: 'ok' },
    errorLog: null,
    executorCalls: 1
  });

  const activity = await assertWorkerReuseAcrossHelperAndRender(page);
  expect(activity.instances[0].initializations).toBe(1);
  expect(activity.instances[0].runs).toHaveLength(2);
  expect(activity.instances[0].helpers.length).toBeGreaterThan(0);
  expect(activity.instances[0].helpers.some(({ operation }) => (
    operation === 'hydrateProteinLosatTsv'
  ))).toBe(true);
  expect(activity.instances[0].helpers.some((helper) => (
    helper.transferCount > 0 && helper.transferredBytes > 0
  ))).toBe(true);
  expect(activity.instances[0].terminated).toBe(false);
  expect(await page.evaluate(() => ({
    mainLoaderPresent: typeof window.loadPyodide === 'function',
    mainRuntimeStatePresent: Object.prototype.hasOwnProperty.call(
      window.__GBDRAW_APP__,
      'pyodideReady'
    )
  }))).toEqual({ mainLoaderPresent: false, mainRuntimeStatePresent: false });
});

test('Gallery session colors, record labels, and feature labels survive regeneration', async ({ page }) => {
  test.setTimeout(240000);
  page.on('dialog', (dialog) => dialog.accept());
  await openApp(page, { waitForPalette: false });

  const imported = await evaluateWithRetainedPromise(page, async () => {
    const inspectSettings = (app) => ({
      cdsColor: app.currentColors.CDS,
      labels: app.form.show_labels_linear,
      // The BGC Session names each record in the record's own Definition.
      recordLabels: app.linearSeqs.map((record) => record.definition),
      rules: app.manualSpecificRules.map(({ feat, qual, val, color, cap }) => ({
        feat, qual, val, color, cap
      }))
    });
    const response = await fetch('/gbdraw/web/gallery/sessions/BGC0000708-BGC0000713.gbdraw-session.json');
    const text = await response.text();
    const file = new File([text], 'BGC0000708-BGC0000713.gbdraw-session.json', {
      type: 'application/json'
    });
    const result = await window.__GBDRAW_APP__.importSession({
      target: { files: [file], value: '' }
    });
    return {
      result: result?.status,
      settings: inspectSettings(window.__GBDRAW_APP__)
    };
  });

  expect(imported.result).toBe('ok');
  expect(await page.evaluate(() => Object.prototype.hasOwnProperty.call(
    window.__GBDRAW_APP__,
    'pyodideReady'
  ))).toBe(false);
  expect(imported.settings.cdsColor).toBe('#808080');
  expect(imported.settings.labels).toBe('first');
  expect(imported.settings.recordLabels[0]).toContain('Streptomyces lividus');
  expect(imported.settings.rules).toHaveLength(4);

  const regenerated = await evaluateWithRetainedPromise(page, async () => {
    const app = window.__GBDRAW_APP__;
    const inspectSettings = () => ({
      cdsColor: app.currentColors.CDS,
      labels: app.form.show_labels_linear,
      // The BGC Session names each record in the record's own Definition.
      recordLabels: app.linearSeqs.map((record) => record.definition),
      rules: app.manualSpecificRules.map(({ feat, qual, val, color, cap }) => ({
        feat, qual, val, color, cap
      }))
    });
    const result = await app.runAnalysis();
    const content = app.currentResults()[0]?.content || '';
    return {
      result,
      errorLog: app.errorLog,
      settings: inspectSettings(),
      hasCoreColor: content.toLowerCase().includes('#810e15'),
      hasCustomRecordLabel: content.includes('Streptomyces lividus')
    };
  });

  expect(regenerated.result).toEqual({ status: 'ok' });
  expect(regenerated.errorLog).toBeNull();
  expect(regenerated.settings).toEqual(imported.settings);
  expect(regenerated.hasCoreColor).toBe(true);
  expect(regenerated.hasCustomRecordLabel).toBe(true);
});

// OV-1049: a fresh protein analysis must not make the next unrelated Generate rebuild the comparison stages.
// The analysis writes the pairwise tables that the next request reads from the committed Session, so each table
// has to keep one resource identity (resource id and Worker cache token) from the analysis to the next request.
test('the Generate after a fresh protein analysis rebuilds only what its own change needs', async ({ page }) => {
  test.setTimeout(480000);
  page.on('dialog', (dialog) => dialog.accept());
  await page.addInitScript(() => {
    window.__OV1049_DIAGNOSTICS__ = [];
    window.__GBDRAW_TEST_HOOKS__ = {
      ...(window.__GBDRAW_TEST_HOOKS__ || {}),
      onSessionLifecycleEvent: (event) => {
        if (event?.name === 'python-diagnostics') window.__OV1049_DIAGNOSTICS__.push(event.metrics || {});
      }
    };
  });
  await openApp(page, { waitForPalette: false });
  const imported = await evaluateWithRetainedPromise(page, async () => {
    const response = await fetch('/gbdraw/web/gallery/sessions/BGC0000708-BGC0000713.gbdraw-session.json');
    const file = new File([await response.text()], 'BGC0000708-BGC0000713.gbdraw-session.json', {
      type: 'application/json'
    });
    return (await window.__GBDRAW_APP__.importSession({ target: { files: [file], value: '' } }))?.status;
  });
  expect(imported).toBe('ok');

  const setInput = async (selector, value) => {
    const input = page.locator(selector).first();
    await input.evaluate((element) => {
      for (let details = element.closest('details'); details; details = details.parentElement?.closest('details')) {
        details.open = true;
      }
    });
    await input.fill(String(value));
    await input.dispatchEvent('change');
  };
  // The render-stage builds of one Generate, by stage name.
  const generate = () => evaluateWithRetainedPromise(page, async () => {
    window.__OV1049_DIAGNOSTICS__.length = 0;
    const app = window.__GBDRAW_APP__;
    const result = await app.runAnalysis();
    const metrics = window.__OV1049_DIAGNOSTICS__.at(-1) || {};
    return {
      result,
      errorLog: app.errorLog,
      builds: Object.fromEntries(Object.entries(metrics)
        .filter(([name]) => name.startsWith('renderStageBuild:'))
        .map(([name, count]) => [name.slice('renderStageBuild:'.length), count]))
    };
  });
  const ok = (run) => {
    expect(run.result).toEqual({ status: 'ok' });
    expect(run.errorLog).toBeNull();
    return run.builds;
  };

  ok(await generate());
  await setInput('input[aria-label="Plot Title"]', 'Title A');
  const withoutAnalysis = ok(await generate());
  await setInput('input[aria-label="Linear comparison minimum identity"]', '50');
  const freshAnalysis = ok(await generate());
  // The identity filter draws other comparison rows, so this Generate rebuilds the comparison stages.
  expect(freshAnalysis.comparison_rows).toBeGreaterThan(0);
  await setInput('input[aria-label="Plot Title"]', 'Title B');
  const afterAnalysis = ok(await generate());
  expect(afterAnalysis).toEqual(withoutAnalysis);
});
