const { test, expect } = require('@playwright/test');
const fs = require('node:fs/promises');
const { gunzipSync } = require('node:zlib');
const { openApp, assertDiagramWorkerIdle, generateAndWaitForResult, getDiagramWorkerActivity } = require('./helpers/app-lifecycle.cjs');
const { capture, assertCoherent, settle } = require('./helpers/visual-state.cjs');
const path = require('node:path');

const snapshot = page => page.evaluate(async () => {
  const { state: s } = await import('./js/state.js');
  const { buildConfigData, getCommittedCanonicalSession } = await import('./js/services/config.js');
  return {
    mode: s.mode.value, config: buildConfigData(),
    circularSources: ['c_gb', 'c_gff', 'c_fasta'].map(key => Boolean(s.files[key])),
    linearSources: s.linearSeqs.map(row => ['gb', 'gff', 'fasta'].map(key => Boolean(row[key]))),
    results: s.results.value.length, catalog: s.featureCatalog.value,
    committed: getCommittedCanonicalSession(), error: s.errorLog.value
  };
});

const assertSourceFree = state => {
  expect(state.circularSources).toEqual([false, false, false]);
  expect(state.linearSources.flat().some(Boolean)).toBe(false);
  expect(state.results).toBe(0);
  expect(state.catalog).toBeNull();
  expect(state.committed).toBeNull();
};

const reveal = async locator => {
  for (const details of await locator.locator('xpath=ancestor::details').all()) {
    if (await details.getAttribute('open') === null) await details.locator(':scope > summary').click();
  }
  return locator;
};

const save = async (page, testInfo, label) => {
  let downloaded;
  const listener = download => { downloaded = download; };
  page.on('download', listener);
  try {
    await page.getByRole('button', { name: 'Save Session', exact: true }).click();
    // Await operation settlement, then assert once. Preserve the product error
    // as the red witness instead of reporting only a download timeout.
    await expect.poll(async () => Boolean(downloaded || (await snapshot(page)).error), { timeout: 180_000 }).toBe(true);
    const boundary = await snapshot(page);
    await fs.writeFile(testInfo.outputPath(`${label}-save-boundary.json`), JSON.stringify(boundary, null, 2));
    expect(boundary.error, 'Ordinary source-free Save must produce a Session download').toBeNull();
    expect(downloaded).toBeTruthy();
    const file = testInfo.outputPath(`${label}.gbdraw-session.json.gz`);
    await downloaded.saveAs(file);
    return { file, document: JSON.parse(gunzipSync(await fs.readFile(file))) };
  } finally {
    page.off('download', listener);
  }
};

const loadFile = async (page, file, message = 'Session loaded successfully!') => {
  const dialog = page.waitForEvent('dialog');
  await page.locator('input[accept^=".json,"]').setInputFiles(file);
  const actual = (await dialog).message();
  await page.waitForFunction(() => !window.__GBDRAW_APP__.sessionImportPending);
  if (typeof message === 'string') expect(actual).toBe(message);
  else expect(actual).toMatch(message);
};

test('settings-only Session preserves non-default Circular and Linear profiles through real Save and fresh Load', async ({ browser, viewport }, testInfo) => {
  test.setTimeout(1_800_000);
  const contexts = [];
  const externalRequests = [];
  const fresh = async () => {
    const context = await browser.newContext({ viewport });
    contexts.push(context);
    await context.route('**/*', route => {
      if (new URL(route.request().url()).hostname === '127.0.0.1') return route.continue();
      externalRequests.push(route.request().url());
      return route.abort();
    });
    const page = await context.newPage();
    page.setDefaultTimeout(20000);
    page.on('dialog', dialog => dialog.accept(dialog.type() === 'prompt' ? 'J08 settings' : undefined));
    await openApp(page);
    return page;
  };
  try {
    const page = await fresh();
    assertSourceFree(await snapshot(page));
    await (await reveal(page.locator('#circular-label-mode'))).selectOption('both');
    await (await reveal(page.locator('#circular-track-preset'))).selectOption('middle');
    await page.getByRole('button', { name: 'Linear', exact: true }).click();
    await (await reveal(page.locator('#linear-show-labels'))).selectOption('all');
    await page.getByRole('button', { name: 'Circular', exact: true }).click();
    const before = await snapshot(page);
    expect(before.config.form.labels_mode).toBe('both');
    expect(before.config.form.track_type).toBe('middle');
    const saved = await save(page, testInfo, 'first');
    assertSourceFree(await snapshot(page));
    await assertDiagramWorkerIdle(page);

    const restored = await fresh();
    const loaded = restored.waitForEvent('dialog');
    await restored.locator('input[accept^=".json,"]').setInputFiles(saved.file);
    expect((await loaded).message()).toBe('Session loaded successfully!');
    await restored.waitForFunction(() => !window.__GBDRAW_APP__.sessionImportPending);
    const after = await snapshot(restored);
    assertSourceFree(after);
    expect(after.mode).toBe('circular');
    expect(after.config.form.labels_mode).toBe('both');
    expect(after.config.form.track_type).toBe('middle');
    expect(after.config).toEqual(before.config);
    await assertDiagramWorkerIdle(restored);
    const again = await save(restored, testInfo, 'second');
    expect(again.document.config).toEqual(saved.document.config);
    expect(again.document.resources).toEqual(saved.document.resources);
    await restored.getByRole('button', { name: 'Linear', exact: true }).click();
    expect((await snapshot(restored)).config.form.show_labels_linear).toBe('all');
    await restored.getByRole('button', { name: 'Circular', exact: true }).click();
    expect((await snapshot(restored)).config.form.labels_mode).toBe('both');
    await assertDiagramWorkerIdle(restored);
    // Original J08's real source/Generate/Undo/Redo suffix, with restored
    // non-default settings checked against the resulting request and output.
    const missingSource = await generateAndWaitForResult(restored, { expectedStatus: 'error' });
    expect(missingSource.errorSummary).toMatch(/input|GenBank|file/i);
    assertSourceFree(await snapshot(restored));
    await restored.getByLabel('GenBank/DDBJ File', { exact: true }).setInputFiles(
      path.join(process.cwd(), 'tests/fixtures/sessions/cli-web-mito.gb'));
    await restored.waitForFunction(async () =>
      (await import('./js/state.js')).state.circularRecordDiscovery.status !== 'loading');
    await generateAndWaitForResult(restored);
    const generated = await snapshot(restored);
    expect(generated.committed.renderRequest.diagramOptions.configOverrides['canvas.circular.track_type']).toBe('middle');
    expect(generated.committed.renderRequest.diagramOptions.configOverrides['labels.circular.scope']).toBe('both');
    const visual = await capture(restored, testInfo, 'generated-settings');
    await assertCoherent(visual, 'settings applied by Generate');
    expect(await restored.locator('.origin-top svg text[data-label-feature-id]').count()).toBeGreaterThan(0);
    await restored.getByRole('button', { name: 'Undo', exact: true }).click();
    await settle(restored);
    const undone = await snapshot(restored);
    expect(undone.results).toBe(0);
    expect(undone.committed).toBeNull();
    expect(undone.config.form.labels_mode).toBe('both');
    await restored.getByRole('button', { name: 'Redo', exact: true }).click();
    await settle(restored);
    expect((await snapshot(restored)).committed.renderRequest).toEqual(generated.committed.renderRequest);
    const redone = await capture(restored, testInfo, 'redone-settings');
    await assertCoherent(redone, 'Redo Generate');
    expect(redone.completed.selected).toEqual(visual.completed.selected);
    await generateAndWaitForResult(restored);
    const repeated = await capture(restored, testInfo, 'repeated-settings');
    await assertCoherent(repeated, 'repeated Generate');
    expect(repeated.completed.selected).toEqual(visual.completed.selected);
    expect(externalRequests).toEqual([]);
  } finally {
    for (const context of contexts) await context.close();
  }
});

test('settings-only Load replaces existing work and rejected candidates preserve the complete prior Session', async ({ page }, testInfo) => {
  test.setTimeout(1_800_000);
  page.setDefaultTimeout(20000);
  page.on('dialog', dialog => dialog.accept(dialog.type() === 'prompt' ? 'Full control' : undefined));
  await openApp(page);
  await loadFile(page, path.join(process.cwd(), 'gbdraw/web/gallery/sessions/HmmtDNA_basic_circular.gbdraw-session.json'));
  const committed = (await snapshot(page)).committed.renderRequest;
  await page.getByRole('button', { name: 'Linear', exact: true }).click();
  await settle(page);
  const inactive = await snapshot(page);
  expect(inactive.circularSources[0]).toBe(true);
  expect(inactive.linearSources.flat().some(Boolean)).toBe(false);
  await assertCoherent(await capture(page, testInfo, 'inactive-source'), 'mode switch preserves the full Session');
  const inactiveSave = await save(page, testInfo, 'inactive-source');
  expect(inactiveSave.document.renderRequest).toEqual(committed);
  await loadFile(page, inactiveSave.file);
  expect((await snapshot(page)).circularSources[0]).toBe(true);
  expect((await snapshot(page)).committed.renderRequest).toEqual(committed);
  const full = await save(page, testInfo, 'full-control');
  const before = await snapshot(page);
  const visual = await capture(page, testInfo, 'before-rejected-loads');
  await assertCoherent(visual, 'prior valid Session');
  const cases = [
    ['request', d => { delete d.renderRequest; }, /canonical render request/],
    ['null-request-with-result', d => { d.renderRequest = null; }, /biological sources|committed render artifacts/],
    ['resource', d => { d.resources = {}; }, /[Mm]issing.*resource|resource.*missing/],
    ['config', d => { d.config.form = []; }, /active form/],
    ['binding', d => { d.webFiles.bindings.c_gb.resourceId = 'absent'; }, /[Mm]issing.*resource/],
    ['version', d => { d.version = 999; }, /newer than/],
    // This valid-shaped profile reaches the existing apply transaction before
    // its managed flag is rejected, exercising rollback after reset has begun.
    ['profile-rollback', d => { d.config.modeProfiles.profiles.circular.managed.axis_stroke_color = 'invalid'; }, /Invalid managed flag/]
  ];
  for (const [name, mutate, error] of cases) {
    const candidate = structuredClone(full.document);
    mutate(candidate);
    const file = testInfo.outputPath(`${name}.json`);
    await fs.writeFile(file, JSON.stringify(candidate));
    await loadFile(page, file, error);
    expect(await snapshot(page), `${name}: complete state at failed Load completion`).toEqual(before);
    const after = await capture(page, testInfo, `rejected-${name}`);
    await assertCoherent(after, `${name}: rejected Load`);
    expect(after.completed.selected).toEqual(visual.completed.selected);
  }
  const priorActivity = await getDiagramWorkerActivity(page);
  await loadFile(page, path.join(process.cwd(), 'tests/fixtures/sessions/settings-only.v42.json.gz'));
  const replaced = await snapshot(page);
  assertSourceFree(replaced);
  expect(replaced.config.form.labels_mode).toBe('both');
  expect(replaced.config.form.track_type).toBe('middle');
  expect(await page.locator('.origin-top svg').count()).toBe(0);
  const replacedActivity = await getDiagramWorkerActivity(page);
  for (const field of ['constructions', 'initializations', 'helpers', 'runs']) {
    expect(replacedActivity[field], `settings-only Load starts no Worker work: ${field}`).toBe(priorActivity[field]);
  }
  const saved = await save(page, testInfo, 'replacement');
  expect(saved.document.renderRequest).toBeNull();
  expect(saved.document.resources).toEqual({});
});

test('both mode profiles and auxiliary file bytes survive settings-only Save Load Save', async ({ page, browser }, testInfo) => {
  test.setTimeout(1_800_000);
  page.setDefaultTimeout(20000);
  page.on('dialog', dialog => dialog.accept(dialog.type() === 'prompt' ? 'Profiles and priority' : undefined));
  await openApp(page);
  // Exercise the existing mode-profile owner with distinct explicit values.
  await page.evaluate(() => { window.__GBDRAW_APP__.adv.axis_stroke_color = '#13579b'; });
  await page.getByRole('button', { name: 'Linear', exact: true }).click();
  await (await reveal(page.locator('#linear-show-labels'))).selectOption('all');
  await page.evaluate(() => { window.__GBDRAW_APP__.adv.axis_stroke_color = '#2468ac'; });
  // Empty input rows can carry auxiliary depth files and retain their identity.
  await page.getByRole('button', { name: 'Add sequence', exact: true }).first().click();
  await page.evaluate(() => {
    window.__GBDRAW_APP__.linearSeqs[1].depth = new File(['chr1\t1\t10\n'], 'depth.tsv', { lastModified: 123 });
  });
  const bytes = Buffer.from('CDS\tgene,product,locus_tag\r\n');
  await page.getByLabel('Priority File (TSV)', { exact: true }).setInputFiles({
    name: 'priority.tsv', mimeType: 'text/tab-separated-values', buffer: bytes
  });
  const original = await save(page, testInfo, 'profiles');
  await assertDiagramWorkerIdle(page);
  const binding = original.document.webFiles.bindings.qualifier_priority;
  expect(Buffer.from(original.document.resources[binding.resourceId].data, 'base64')).toEqual(bytes);
  expect(original.document.config.modeProfiles.profiles.circular.values.axis_stroke_color).toBe('#13579b');
  expect(original.document.config.modeProfiles.profiles.linear.values.axis_stroke_color).toBe('#2468ac');
  const context = await browser.newContext();
  try {
    const restored = await context.newPage();
    restored.on('dialog', dialog => dialog.accept());
    await openApp(restored);
    await loadFile(restored, original.file);
    const state = await snapshot(restored);
    assertSourceFree(state);
    expect(state.mode).toBe('linear');
    expect(state.config.adv.axis_stroke_color).toBe('#2468ac');
    expect(state.config.modeProfiles).toEqual(original.document.config.modeProfiles);
    const again = await save(restored, testInfo, 'profiles-again');
    expect(again.document.config).toEqual(original.document.config);
    expect(again.document.webFiles).toEqual(original.document.webFiles);
    expect(again.document.resources).toEqual(original.document.resources);
    await restored.getByRole('button', { name: 'Circular', exact: true }).click();
    expect((await snapshot(restored)).config.adv.axis_stroke_color).toBe('#13579b');
    await restored.getByRole('button', { name: 'Linear', exact: true }).click();
    expect((await snapshot(restored)).config.adv.axis_stroke_color).toBe('#2468ac');
    await assertDiagramWorkerIdle(restored);
  } finally { await context.close(); }
});
