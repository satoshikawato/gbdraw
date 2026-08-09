const { test, expect } = require('@playwright/test');

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
  await page.waitForFunction(
    () => window.__GBDRAW_APP__?.pyodideReady === true,
    null,
    { timeout: 180000 }
  );
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
