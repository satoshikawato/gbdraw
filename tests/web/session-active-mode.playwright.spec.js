const { test, expect } = require('@playwright/test');
const fs = require('node:fs/promises');
const path = require('node:path');
const { gunzipSync } = require('node:zlib');
const { execFile } = require('node:child_process');
const { createHash } = require('node:crypto');
const { promisify } = require('node:util');
const {
  openApp, assertSessionLoadLeftWorkerIdle, getDiagramWorkerActivity,
  generateAndWaitForResult, evaluateWithRetainedPromise
} = require('./helpers/app-lifecycle.cjs');

const fixture = path.join(process.cwd(), 'gbdraw/web/gallery/sessions/tobacco-chloroplast.gbdraw-session.json');
const divergentFixture = path.join(process.cwd(), 'gbdraw/web/gallery/sessions/HmmtDNA_basic_circular.gbdraw-session.json');
const linearSource = path.join(process.cwd(), 'tests/fixtures/sessions/cli-web-mito.gb');
const sessionInput = 'input[accept^=".json,"]';
const readSession = async file => JSON.parse(gunzipSync(await fs.readFile(file)));

const load = async (page, file) => {
  const dialogs = [];
  const accept = async dialog => { dialogs.push(dialog.message()); await dialog.accept(); };
  page.on('dialog', accept);
  try {
    await page.locator(sessionInput).setInputFiles(file);
    await page.waitForFunction(() => !window.__GBDRAW_APP__.sessionImportPending, null, { timeout: 180_000 });
    const error = await page.evaluate(() => window.__GBDRAW_APP__.errorLog);
    expect(error, `Load ${path.basename(file)} failed`).toBeNull();
    expect(dialogs).toContain('Session loaded successfully!');
  } finally {
    page.off('dialog', accept);
  }
};

const save = async (page, file) => {
  const downloadPromise = page.waitForEvent('download', { timeout: 180_000 });
  await page.getByRole('banner').getByRole('button', { name: 'Save Session', exact: true }).click();
  await (await downloadPromise).saveAs(file);
  return readSession(file);
};

const snapshot = page => page.evaluate(async () => {
  const { state } = await import('./js/state.js');
  const { buildConfigData, getCommittedCanonicalRenderRequest } = await import('./js/services/config.js');
  const { readFileBytes } = await import('./js/services/file-content-cache.js');
  const digest = async bytes => [...new Uint8Array(await crypto.subtle.digest('SHA-256', bytes))]
    .map(byte => byte.toString(16).padStart(2, '0')).join('');
  const fileInfo = async file => file ? {
    name: file.name, type: file.type, lastModified: file.lastModified,
    bytes: await digest(await readFileBytes(file))
  } : null;
  return {
    mode: state.mode.value,
    generatedMode: state.generatedMode.value,
    config: buildConfigData(),
    request: getCommittedCanonicalRenderRequest(),
    results: state.results.value.map(result => ({ name: result.name, content: result.content })),
    selectedResultIndex: state.selectedResultIndex.value,
    featurePanelTab: state.featurePanelTab.value,
    downloadDpi: state.downloadDpi.value,
    circularSource: await fileInfo(state.files.c_gb),
    linearSource: await fileInfo(state.linearSeqs[0]?.gb),
    history: [window.__GBDRAW_HISTORY__.getUndoCount(), window.__GBDRAW_HISTORY__.getRedoCount()]
  };
});

test('current biological Save, fresh Load, and re-save keep a Linear draft beside a Circular Result', async ({ page, browser }, info) => {
  test.setTimeout(360_000);
  await openApp(page);
  await load(page, divergentFixture);
  await page.getByRole('button', { name: 'Linear', exact: true }).click();
  await page.getByTestId('linear-genbank-1').setInputFiles(linearSource);
  await page.evaluate(async () => {
    const { state } = await import('./js/state.js');
    await window.__GBDRAW_HISTORY__.runUndoable('Pending Linear scale size', () => {
      state.adv.scale_font_size = 19;
    });
    state.featurePanelTab.value = 'labels';
  });
  const before = await snapshot(page);
  expect(before.mode).toBe('linear');
  expect(before.generatedMode).toBe('circular');
  expect(before.request.mode).toBe('circular');
  expect(before.config.adv.scale_font_size).toBe(19);
  expect(before.linearSource).toBeTruthy();
  expect(before.history[0]).toBeGreaterThan(0);

  const firstFile = info.outputPath('linear-draft.gbdraw-session.json.gz');
  const first = await save(page, firstFile);
  expect(first.version).toBe(44);
  expect(first.ui.mode).toBe('linear');
  expect(first.config.modeProfiles.activeMode).toBe('linear');
  expect(first.renderRequest.mode).toBe('circular');
  expect(first.config.adv.scale_font_size).toBe(19);

  const context = await browser.newContext({
    baseURL: `http://127.0.0.1:${process.env.GBDRAW_WEB_TEST_PORT || 4173}`,
    acceptDownloads: true
  });
  try {
    const fresh = await context.newPage();
    await openApp(fresh);
    await load(fresh, firstFile);
    const loadWorker = await getDiagramWorkerActivity(fresh);
    expect(loadWorker.runs).toBe(0);
    expect(loadWorker.constructions).toBeLessThanOrEqual(1);
    expect(loadWorker.instances.flatMap(instance => instance.helpers.map(helper => helper.operation)))
      .toEqual(loadWorker.helpers ? ['validateConfigOverrides'] : []);
    const restored = await snapshot(fresh);
    expect(restored.mode).toBe('linear');
    expect(restored.generatedMode).toBe('circular');
    expect(restored.config.modeProfiles).toEqual(before.config.modeProfiles);
    expect(restored.config.adv.scale_font_size).toBe(19);
    expect(restored.circularSource).toEqual(before.circularSource);
    expect(restored.linearSource).toEqual(before.linearSource);
    expect(restored.request).toEqual(before.request);
    expect(restored.results).toEqual(before.results);
    expect(restored.selectedResultIndex).toBe(before.selectedResultIndex);
    expect(restored.featurePanelTab).toBe(before.featurePanelTab);
    expect(restored.downloadDpi).toBe(before.downloadDpi);
    expect(restored.history).toEqual([0, 0]);

    const second = await save(fresh, info.outputPath('linear-draft-resaved.gbdraw-session.json.gz'));
    expect(second.ui.mode).toBe('linear');
    expect(second.config.modeProfiles).toEqual(first.config.modeProfiles);
    expect(second.config.adv.scale_font_size).toBe(19);
    expect(second.renderRequest).toEqual(first.renderRequest);
    expect(second.results).toEqual(first.results);
    expect(second.resources).toEqual(first.resources);
    expect(second.webFiles.bindings.c_gb).toEqual(first.webFiles.bindings.c_gb);
    expect(second.webFiles.bindings.linearSeqs[0].gb).toEqual(first.webFiles.bindings.linearSeqs[0].gb);
    expect(second.webFiles.resourceOriginalNames).toEqual(first.webFiles.resourceOriginalNames);
    expect(second.ui.selectedResultIndex).toBe(first.ui.selectedResultIndex);
    expect(second.ui.featurePanelTab).toBe(first.ui.featurePanelTab);

    await generateAndWaitForResult(fresh);
    const generated = await snapshot(fresh);
    expect(generated.mode).toBe('linear');
    expect(generated.generatedMode).toBe('linear');
    expect(generated.request.mode).toBe('linear');
    expect(generated.request.diagramOptions.configOverrides['objects.scale.font_size.short']).toBe(19);
    expect(generated.request.diagramOptions.configOverrides['objects.scale.font_size.long']).toBe(19);
    expect(generated.results[0].content).not.toBe(restored.results[0].content);
    expect(generated.linearSource).toEqual(before.linearSource);
  } finally {
    await context.close();
  }
});

test('each mode keeps its own title and fonts while missing Linear layout starts fresh', async ({ page, browser }, info) => {
  test.setTimeout(360_000);
  await openApp(page);
  await load(page, fixture);
  const titles = page.locator('summary[aria-label="Titles and Record Labels"]');
  if (!(await titles.evaluate(summary => summary.parentElement.open))) await titles.press('Enter');
  const plotTitle = page.getByRole('textbox', { name: 'Plot Title', exact: true });
  await plotTitle.fill('Circular title');
  await plotTitle.press('Tab');
  const draft = target => target.evaluate(async () => {
    const { state } = await import('./js/state.js');
    const size = value => (value === null || value === undefined || value === '' ? null : Number(value));
    return {
      title: state.form.plot_title, titleFont: size(state.adv.plot_title_font_size),
      definitionFont: size(state.adv.def_font_size),
      rows: state.linearRecordLayoutEnabled.value, replicon: state.adv.linear_show_replicon
    };
  });
  const circular = await draft(page);
  expect(circular).toMatchObject({ title: 'Circular title', titleFont: 32, definitionFont: 28 });

  // The fixture omits config.linearRecordLayout; Linear starts from fresh defaults.
  await page.getByRole('button', { name: 'Linear', exact: true }).click();
  expect(await draft(page)).toEqual({ title: '', titleFont: null, definitionFont: null, rows: true, replicon: false });
  await expect(plotTitle).toHaveValue('');
  await plotTitle.fill('Linear title');
  await plotTitle.press('Tab');
  await page.getByRole('button', { name: 'Circular', exact: true }).click();
  expect(await draft(page)).toMatchObject({ title: 'Circular title', titleFont: 32, definitionFont: 28 });

  const file = info.outputPath('per-mode-title.gbdraw-session.json.gz');
  const saved = await save(page, file);
  expect(saved.config.modeProfiles.profiles.circular.values.plot_title).toBe('Circular title');
  expect(saved.config.modeProfiles.profiles.linear.values.plot_title).toBe('Linear title');
  const context = await browser.newContext({
    baseURL: `http://127.0.0.1:${process.env.GBDRAW_WEB_TEST_PORT || 4173}`
  });
  try {
    const fresh = await context.newPage();
    await openApp(fresh);
    await load(fresh, file);
    expect(await draft(fresh)).toMatchObject({ title: 'Circular title', titleFont: 32, definitionFont: 28 });
    await fresh.getByRole('button', { name: 'Linear', exact: true }).click();
    expect(await draft(fresh)).toEqual({ title: 'Linear title', titleFont: null, definitionFont: null, rows: true, replicon: false });
  } finally {
    await context.close();
  }
});

test('matching current, historical, CLI-origin, and settings-only modes keep their Load behavior', async ({ page, browser }, info) => {
  test.setTimeout(360_000);
  await openApp(page);
  await load(page, fixture);
  let state = await snapshot(page);
  expect([state.mode, state.generatedMode, state.request.mode]).toEqual(['circular', 'circular', 'circular']);
  const previewWorker = await getDiagramWorkerActivity(page);
  expect(previewWorker.runs).toBe(0);
  expect(previewWorker.constructions).toBeLessThanOrEqual(1);

  const current = await save(page, info.outputPath('matching.gbdraw-session.json.gz'));
  expect(current.ui.mode).toBe('circular');
  await load(page, info.outputPath('matching.gbdraw-session.json.gz'));
  state = await snapshot(page);
  expect([state.mode, state.generatedMode, state.request.mode]).toEqual(['circular', 'circular', 'circular']);

  await load(page, path.join(process.cwd(), 'tests/fixtures/sessions/single.v41-bindings1.json'));
  state = await snapshot(page);
  expect([state.mode, state.generatedMode, state.request.mode]).toEqual(['circular', 'circular', 'circular']);

  const cliPrefix = info.outputPath('cli-origin');
  const cliFile = `${cliPrefix}.gbdraw-session.json`;
  await promisify(execFile)('python', ['-m', 'gbdraw.cli', 'circular', '--gbk', linearSource,
    '-o', cliPrefix, '--format', 'svg', '--save_session'], {
    cwd: info.outputDir, env: { ...process.env, PYTHONPATH: process.cwd() }, timeout: 180_000
  });
  const cli = JSON.parse(await fs.readFile(cliFile, 'utf8'));
  expect(cli.config).toBeUndefined();
  await load(page, cliFile);
  state = await snapshot(page);
  expect([state.mode, state.generatedMode, state.request.mode]).toEqual(['circular', 'circular', 'circular']);
  delete cli.ui.mode;
  const missingModeFile = info.outputPath('cli-origin-missing-ui-mode.gbdraw-session.json');
  await fs.writeFile(missingModeFile, JSON.stringify(cli));
  await load(page, missingModeFile);
  state = await snapshot(page);
  expect([state.mode, state.generatedMode, state.request.mode]).toEqual(['circular', 'circular', 'circular']);

  const settingsFile = path.join(process.cwd(), 'tests/fixtures/sessions/settings-only.v42.json.gz');
  const settings = await readSession(settingsFile);
  const context = await browser.newContext({
    baseURL: `http://127.0.0.1:${process.env.GBDRAW_WEB_TEST_PORT || 4173}`
  });
  try {
    const fresh = await context.newPage();
    await openApp(fresh);
    await load(fresh, settingsFile);
    state = await snapshot(fresh);
    expect(state.mode).toBe(settings.ui.mode);
    expect(state.request).toBeNull();
    expect(state.results).toEqual([]);
    expect(state.config.modeProfiles.activeMode).toBe(settings.config.modeProfiles.activeMode);
    for (const mode of ['circular', 'linear']) {
      expect(state.config.modeProfiles.profiles[mode].values.axis_stroke_color)
        .toBe(settings.config.modeProfiles.profiles[mode].values.axis_stroke_color);
    }
    await assertSessionLoadLeftWorkerIdle(fresh);
  } finally {
    await context.close();
  }
});

test('failed import after draft application restores the prior mode, profiles, sources, and Result', async ({ page }, info) => {
  test.setTimeout(240_000);
  await openApp(page);
  await load(page, fixture);
  await page.getByRole('button', { name: 'Linear', exact: true }).click();
  await page.getByTestId('linear-genbank-1').setInputFiles(linearSource);
  const before = await snapshot(page);
  const failure = await evaluateWithRetainedPromise(page, async () => {
    const name = 'tobacco-chloroplast.gbdraw-session.json';
    const file = new File([await (await fetch(`/gbdraw/web/gallery/sessions/${name}`)).arrayBuffer()], name);
    const { importSession } = await import('./js/services/config.js');
    const result = await importSession({ target: { files: [file], value: 'selected' } }, {
      beforePreviewMount: () => { throw new Error('S02 forced rollback'); }
    });
    return { status: result.status, error: result.error };
  });
  expect(failure).toMatchObject({ status: 'error', error: { code: 'UNKNOWN', stage: 'request-validation' } });
  expect(JSON.stringify(failure)).not.toContain('S02 forced rollback');
  const restored = await snapshot(page);
  for (const key of ['mode', 'generatedMode', 'config', 'request', 'selectedResultIndex',
    'featurePanelTab', 'downloadDpi', 'circularSource', 'linearSource', 'history']) {
    expect(restored[key], `rollback ${key}`).toEqual(before[key]);
  }
  const canonicalSvgHashes = async results => page.evaluate(async rows => {
    const { sanitizeSvgContent } = await import('./js/services/svg-sanitization.js');
    const { serializeCleanSvg } = await import('./js/services/svg-serialization.js');
    return rows.map(({ content }) => {
      const svg = new DOMParser().parseFromString(
        sanitizeSvgContent(content, window.DOMPurify), 'image/svg+xml'
      ).documentElement;
      svg.removeAttribute('xmlns');
      svg.removeAttribute('xmlns:xlink');
      return serializeCleanSvg(svg);
    });
  }, results);
  const digest = content => createHash('sha256').update(content).digest('hex');
  expect((await canonicalSvgHashes(restored.results)).map(digest))
    .toEqual((await canonicalSvgHashes(before.results)).map(digest));
  const saved = await save(page, info.outputPath('after-rollback.gbdraw-session.json.gz'));
  expect(saved.ui.mode).toBe('linear');
  expect(saved.renderRequest.mode).toBe('circular');
});
