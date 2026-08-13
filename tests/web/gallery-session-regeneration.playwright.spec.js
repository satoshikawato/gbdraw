const { test, expect } = require('@playwright/test');

const installDiagramWorkerActivityTracker = async (page) => {
  await page.addInitScript(() => {
    window.__GBDRAW_DIAGRAM_WORKER_ACTIVITY__ = {
      constructions: 0,
      initializations: 0,
      runs: 0
    };
    const NativeWorker = window.Worker;
    window.Worker = new Proxy(NativeWorker, {
      construct(target, args) {
        const worker = Reflect.construct(target, args, target);
        if (!String(args[0] || '').includes('diagram-generation-worker.js')) {
          return worker;
        }
        window.__GBDRAW_DIAGRAM_WORKER_ACTIVITY__.constructions += 1;
        const nativePostMessage = worker.postMessage.bind(worker);
        worker.postMessage = (message, transfer) => {
          if (message?.type === 'init') {
            window.__GBDRAW_DIAGRAM_WORKER_ACTIVITY__.initializations += 1;
          } else if (message?.type === 'run') {
            window.__GBDRAW_DIAGRAM_WORKER_ACTIVITY__.runs += 1;
          }
          if (transfer === undefined) return nativePostMessage(message);
          return nativePostMessage(message, transfer);
        };
        return worker;
      }
    });
  });
};

test('Diagram worker stays idle through session load and starts once on Generate', async ({ page }) => {
  test.setTimeout(240000);
  page.on('dialog', (dialog) => dialog.accept());
  await installDiagramWorkerActivityTracker(page);
  await page.goto('/gbdraw/web/index.html', { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);
  await page.evaluate(() => new Promise((resolve) => requestAnimationFrame(() => {
    requestAnimationFrame(resolve);
  })));
  expect(await page.evaluate(() => window.__GBDRAW_DIAGRAM_WORKER_ACTIVITY__)).toEqual({
    constructions: 0,
    initializations: 0,
    runs: 0
  });

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
  expect(await page.evaluate(() => window.__GBDRAW_DIAGRAM_WORKER_ACTIVITY__)).toEqual({
    constructions: 0,
    initializations: 0,
    runs: 0
  });

  const generated = await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    const result = await app.runAnalysis();
    return {
      result,
      errorLog: app.errorLog,
      resultCount: app.results.length
    };
  });

  expect(generated.result).toEqual({ status: 'ok' });
  expect(generated.errorLog).toBeNull();
  expect(generated.resultCount).toBeGreaterThan(0);
  expect(await page.evaluate(() => window.__GBDRAW_DIAGRAM_WORKER_ACTIVITY__)).toEqual({
    constructions: 1,
    initializations: 1,
    runs: 1
  });
});

test('Gallery session colors, record labels, and feature labels survive regeneration', async ({ page }) => {
  test.setTimeout(240000);
  page.on('dialog', (dialog) => dialog.accept());
  await page.goto('/gbdraw/web/index.html', { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);

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
  expect(await page.evaluate(() => window.__GBDRAW_APP__?.pyodideReady)).toBe(false);
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
