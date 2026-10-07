const { test, expect } = require('@playwright/test');
const fs = require('node:fs/promises');
const path = require('node:path');
const { gunzipSync } = require('node:zlib');
const { promoteRequest } = require('./helpers/request-schema.cjs');
const { execFile } = require('node:child_process');
const { createHash } = require('node:crypto');
const { promisify } = require('node:util');
const {
  CURRENT_SESSION_VERSION,
  assertSessionLoadLeftWorkerIdle,
  evaluateWithRetainedPromise,
  generateAndWaitForResult,
  getDiagramWorkerActivity,
  openApp
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
    config: buildConfigData(state.activeDrawing()),
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
  const circular = await snapshot(page);
  await page.getByRole('button', { name: 'Linear', exact: true }).click();
  await page.getByTestId('linear-genbank-1').setInputFiles(linearSource);
  await page.evaluate(async () => {
    const { state } = await import('./js/state.js');
    await window.__GBDRAW_HISTORY__.runUndoable('Pending Linear scale size', () => {
      state.activeDrawing().adv.scale_font_size = 19;
    });
    state.featurePanelTab.value = 'labels';
  });
  const before = await snapshot(page);
  // E1: Linear shows its own empty Result slot; the Circular Result waits in its slot.
  expect(before.mode).toBe('linear');
  expect(before.generatedMode).toBe('linear');
  expect(before.request).toBeNull();
  expect(before.results).toEqual([]);
  expect(before.config.adv.scale_font_size).toBe(19);
  expect(before.linearSource).toBeTruthy();
  expect(before.history[0]).toBeGreaterThan(0);

  const firstFile = info.outputPath('linear-draft.gbdraw-session.json.gz');
  const first = await save(page, firstFile);
  expect(first.version).toBe(CURRENT_SESSION_VERSION);
  expect(first.ui.mode).toBe('linear');
  expect(first.renderRequest.mode).toBe('circular');
  // The Linear drawing is the Linear slice (Session 46); the Circular slice keeps its own.
  expect(first.modes.linear.config.adv.scale_font_size).toBe(19);
  expect(first.modes.circular.config.adv.scale_font_size).not.toBe(19);

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
    // E1 (contract): a Session with a Result opens on that Result's mode.
    const restored = await snapshot(fresh);
    expect(restored.mode).toBe('circular');
    expect(restored.generatedMode).toBe('circular');
    expect(restored.circularSource).toEqual(before.circularSource);
    expect(restored.linearSource).toEqual(before.linearSource);
    expect(restored.request).toEqual(await promoteRequest(fresh, circular.request));
    expect(restored.results).toEqual(circular.results);
    expect(restored.selectedResultIndex).toBe(circular.selectedResultIndex);
    expect(restored.featurePanelTab).toBe(before.featurePanelTab);
    expect(restored.downloadDpi).toBe(before.downloadDpi);
    expect(restored.history).toEqual([0, 0]);
    // The Linear draft waits in the Linear drawing.
    await fresh.getByRole('button', { name: 'Linear', exact: true }).click();
    const linearDraft = await snapshot(fresh);
    expect(linearDraft.mode).toBe('linear');
    expect(linearDraft.results).toEqual([]);
    expect(linearDraft.config).toEqual(before.config);
    expect(linearDraft.config.adv.scale_font_size).toBe(19);

    const second = await save(fresh, info.outputPath('linear-draft-resaved.gbdraw-session.json.gz'));
    expect(second.ui.mode).toBe('linear');
    expect(second.modes).toEqual(first.modes);
    expect(second.modes.linear.config.adv.scale_font_size).toBe(19);
    expect(second.renderRequest).toEqual(first.renderRequest);
    expect(second.results).toEqual(first.results);
    expect(second.resources).toEqual(first.resources);
    expect(second.webFiles.bindings.c_gb).toEqual(first.webFiles.bindings.c_gb);
    expect(second.webFiles.bindings.linearSeqs[0].gb).toEqual(first.webFiles.bindings.linearSeqs[0].gb);
    expect(second.webFiles.resourceOriginalNames).toEqual(first.webFiles.resourceOriginalNames);
    expect(second.ui.selectedResultIndex).toBe(first.ui.selectedResultIndex);

    await generateAndWaitForResult(fresh);
    const generated = await snapshot(fresh);
    expect(generated.mode).toBe('linear');
    expect(generated.generatedMode).toBe('linear');
    expect(generated.request.mode).toBe('linear');
    expect(generated.request.diagramOptions.configOverrides['objects.scale.font_size.short']).toBe(19);
    expect(generated.request.diagramOptions.configOverrides['objects.scale.font_size.long']).toBe(19);
    expect(generated.results[0].content).not.toBe(circular.results[0].content);
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
      title: state.activeDrawing().form.plot_title, titleFont: size(state.activeDrawing().adv.plot_title_font_size),
      definitionFont: size(state.activeDrawing().adv.def_font_size),
      rows: state.activeDrawing().linearRecordLayoutEnabled.value, replicon: state.activeDrawing().adv.linear_show_replicon
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
  // Each mode's drawing is its slice (Session 46).
  expect(saved.modes.circular.config.form.plot_title).toBe('Circular title');
  expect(saved.modes.linear.config.form.plot_title).toBe('Linear title');
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
    // The shown mode takes the saved flat value and the other mode its saved
    // profile value, each in its own drawing.
    const axisColors = await fresh.evaluate(async () => {
      const { state: live } = await import('./js/state.js');
      return { circular: live.drawings.circular.adv.axis_stroke_color, linear: live.drawings.linear.adv.axis_stroke_color };
    });
    for (const mode of ['circular', 'linear']) {
      expect(axisColors[mode]).toBe(mode === settings.ui.mode
        ? settings.config.adv.axis_stroke_color
        : settings.config.modeProfiles.profiles[mode].values.axis_stroke_color);
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
