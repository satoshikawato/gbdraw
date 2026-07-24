const { test, expect } = require('@playwright/test');
const { createReadStream, existsSync } = require('node:fs');
const { createServer } = require('node:http');
const { extname, resolve, sep } = require('node:path');

const repoRoot = resolve(process.env.GBDRAW_REPO || process.cwd());
const contentTypes = {
  '.html': 'text/html; charset=utf-8',
  '.js': 'text/javascript; charset=utf-8',
  '.mjs': 'text/javascript; charset=utf-8',
  '.json': 'application/json; charset=utf-8',
  '.svg': 'image/svg+xml',
  '.wasm': 'application/wasm',
  '.whl': 'application/octet-stream'
};

let server;
let baseUrl;

test.beforeAll(async () => {
  await new Promise((resolveServer, rejectServer) => {
    server = createServer((request, response) => {
      const url = new URL(request.url || '/', 'http://127.0.0.1');
      const filePath = resolve(repoRoot, decodeURIComponent(url.pathname).replace(/^[/\\]+/, ''));
      if ((!filePath.startsWith(`${repoRoot}${sep}`) && filePath !== repoRoot) || !existsSync(filePath)) {
        response.writeHead(404);
        response.end('Not found');
        return;
      }
      response.writeHead(200, {
        'Content-Type': contentTypes[extname(filePath)] || 'application/octet-stream'
      });
      createReadStream(filePath).pipe(response);
    });
    server.once('error', rejectServer);
    server.listen(0, '127.0.0.1', () => {
      baseUrl = `http://127.0.0.1:${server.address().port}`;
      resolveServer();
    });
  });
});

test.afterAll(async () => {
  await new Promise((resolveClose) => server.close(resolveClose));
});

test('Gallery session colors, record labels, and feature labels survive regeneration', async ({ page }) => {
  test.setTimeout(240000);
  page.on('dialog', (dialog) => dialog.accept());
  await page.goto(`${baseUrl}/gbdraw/web/index.html`, { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);
  await page.waitForFunction(() => window.__GBDRAW_APP__?.pyodideReady === true, null, { timeout: 180000 });

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
    delete legacyOptions.configOverrides;
    delete legacyOptions.featureShapes;
    delete legacyOptions.plotTitle;
    delete legacyOptions.qualifierPriorityFile;
    legacyOptions.colors = {
      colorTable: null,
      colorTableFile: null,
      defaultColors: {
        resourceId: 'colors-default-colors',
        representation: 'canonicalTsv'
      },
      defaultColorsPalette: 'default',
      defaultColorsFile: null
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
