const { test, expect } = require('@playwright/test');
const { readFileSync } = require('node:fs');
const {
  assertDiagramWorkerIdle,
  assertSessionLoadLeftWorkerIdle,
  assertWorkerReuseAcrossHelperAndRender,
  openApp
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

  const imported = await page.evaluate(async () => {
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
    app.setLinearComparisonGlobalAction('losat');
    app.setLinearComparisonLosatMode('blastp');
    app.setLinearComparisonLosatpMode('pairwise');
    return {
      program: app.losatProgram,
      blastpMode: app.losat.blastp.mode,
      cacheEntries: state.losatCache.value.size
    };
  });
  expect(proteinMode).toEqual({
    program: 'blastp',
    blastpMode: 'pairwise',
    cacheEntries: 0
  });

  const generated = await page.evaluate(async () => {
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

  const repeated = await page.evaluate(async () => ({
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

  const imported = await page.evaluate(async () => {
    const inspectSettings = (app) => ({
      cdsColor: app.currentColors.CDS,
      labels: app.form.show_labels_linear,
      recordLabels: app.linearSeqs.map((record) => record.definition),
      rules: app.manualSpecificRules.map(({ feat, qual, val, color, cap }) => ({
        feat, qual, val, color, cap
      }))
    });
    const response = await fetch('/gbdraw/web/gallery/sessions/BGC0000708-BGC0000713.gbdraw-session.json');
    const text = await response.text();
    const legacySession = JSON.parse(text);
    legacySession.version = 39;
    legacySession.renderRequest.schema = 2;
    legacySession.renderRequest.records.forEach((record) => {
      record.presentation = {
        label: null,
        subtitle: null,
        reverseComplement: false,
        gridRow: null,
        gridColumn: null
      };
    });
    const legacyOptions = legacySession.renderRequest.diagramOptions;
    legacyOptions.output.outputPrefix = legacySession.renderRequest.output.prefix;
    legacyOptions.colors = {
      colorTable: null,
      colorTableFile: {
        resourceId: 'colors-color-table-file',
        representation: 'file'
      },
      defaultColors: null,
      defaultColorsPalette: 'default',
      defaultColorsFile: {
        resourceId: 'colors-default-colors-file',
        representation: 'file'
      }
    };
    const file = new File([JSON.stringify(legacySession)], 'BGC0000708-BGC0000713.gbdraw-session.json', {
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
  expect(imported.settings.cdsColor).toBe('#dddddd');
  expect(imported.settings.labels).toBe('first');
  expect(imported.settings.recordLabels[0]).toContain('Streptomyces lividus');
  expect(imported.settings.rules).toHaveLength(4);

  const regenerated = await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    const inspectSettings = () => ({
      cdsColor: app.currentColors.CDS,
      labels: app.form.show_labels_linear,
      recordLabels: app.linearSeqs.map((record) => record.definition),
      rules: app.manualSpecificRules.map(({ feat, qual, val, color, cap }) => ({
        feat, qual, val, color, cap
      }))
    });
    const result = await app.runAnalysis();
    const content = app.results?.[0]?.content || '';
    return {
      result,
      errorLog: app.errorLog,
      settings: inspectSettings(),
      hasCoreColor: content.toLowerCase().includes('#d03535'),
      hasCustomRecordLabel: content.includes('Streptomyces lividus')
    };
  });

  expect(regenerated.result).toEqual({ status: 'ok' });
  expect(regenerated.errorLog).toBeNull();
  expect(regenerated.settings).toEqual(imported.settings);
  expect(regenerated.hasCoreColor).toBe(true);
  expect(regenerated.hasCustomRecordLabel).toBe(true);
});
