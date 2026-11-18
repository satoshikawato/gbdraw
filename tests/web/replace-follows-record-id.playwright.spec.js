const { test, expect } = require('@playwright/test');
const { readFileSync } = require('node:fs');
const { join, resolve } = require('node:path');
const { openApp, evaluateWithRetainedPromise, waitForAppShell, expectSessionLoadedWithoutWarnings } = require('./helpers/app-lifecycle.cjs');

// OV-3846 (D-05): a Replace in the same File card continues each record by
// its record ID when every ID of both files is unique, else by place. A file
// with the same two records in the other order keeps each record's Depth file
// and feature edit with that record: through a Session saved and loaded before
// the next Generate (the Depth rows move at the read after the Replace, the
// edits at that Generate, from the committed request's own file), and through a
// Replace and Generate in one sitting. A's Depth is never drawn on B.
const repoRoot = resolve(process.env.GBDRAW_REPO || process.cwd());
const genbankA = readFileSync(join(repoRoot, 'tests/test_inputs/BGC0000708.gbk'), 'utf8');
const genbankB = readFileSync(join(repoRoot, 'tests/test_inputs/BGC0000709.gbk'), 'utf8');
const depthTsv = (recordId, length, depth) => ['reference_name\tposition\tdepth',
  ...Array.from({ length: Math.ceil(length / 1000) }, (_, index) => `${recordId}\t${Math.min(length, index * 1000 + 1)}\t${depth + (index % 3)}`),
  `${recordId}\t${length}\t${depth}`].join('\n');

// What the shown drawing draws and holds: the record each Depth series is
// drawn on in the Result, the Depth matrix by record ID, and the record of
// each feature edit by record ID.
const drawn = (page) => page.evaluate(async () => {
  const app = window.__GBDRAW_APP__;
  const { state } = await import('./js/state.js');
  const drawing = state.activeDrawing();
  const records = app.circularRecordList.map((record) => [record.identity, record.record_id]);
  const idOf = new Map(records);
  // The record each Depth row of the Result is drawn on (`depth_1` is sample-a's series).
  const svg = new DOMParser().parseFromString(app.currentResults()[0]?.content || '', 'image/svg+xml');
  const drawnOn = Object.fromEntries([...svg.querySelectorAll('g[data-gbdraw-slot-renderer="depth"]')].map((group) => [
    group.getAttribute('data-gbdraw-slot-id'),
    (group.closest('[data-gbdraw-record-id]') ?? group.parentElement?.querySelector('[data-gbdraw-record-id]'))?.getAttribute('data-gbdraw-record-id') ?? null
  ]));
  return {
    recordIds: records.map(([, id]) => id),
    depth: Object.fromEntries(drawing.files.c_depth.map((row, index) => [records[index]?.[1], row.map((ref) => ref?.name ?? null)])),
    drawnOn,
    edits: Object.values(drawing.featureOverrides).map((row) => [idOf.get(row.recordKey) ?? row.recordKey, row.labelText])
  };
});

// Replaces the first File card's GenBank file through its file control, as a pick does.
const replaceCardFile = async (page, name, text) => {
  const card = page.locator('[data-circular-source-card]').first();
  const replace = card.getByRole('button', { name: 'Replace GenBank/DDBJ File' });
  if (!(await replace.isVisible())) await card.locator('summary').click({ position: { x: 8, y: 8 } });
  await expect(replace).toBeVisible();
  await card.locator('input[type=file]').setInputFiles({ name, mimeType: 'text/plain', buffer: Buffer.from(text) });
};
const recordOrder = (page, ids) => page.waitForFunction((joined) => window.__GBDRAW_APP__.circularRecordList
  ?.map((record) => record.record_id).join() === joined, ids.join(), { timeout: 180_000 });

test('a reordered Replace keeps each record\'s Depth file and feature edit through a Save/Load before its Generate and in one sitting', async ({ page }, testInfo) => {
  test.setTimeout(900_000);
  await openApp(page, { waitForPalette: false });
  await page.evaluate(async ({ a, b, depthA, depthB }) => {
    const inputs = await import('./js/services/drawing/inputs.js');
    const { state } = await import('./js/state.js');
    const app = window.__GBDRAW_APP__;
    app.selectDrawing('circular');
    app.setCircularInputType('gb');
    await app.addCircularSourceFiles([new File([`${a}\n${b}`], 'ab.gbk', { type: 'text/plain', lastModified: 1 })]);
    app.files.c_depth = inputs.registerInputFileTree(state.inputFiles, [
      [new File([depthA], 'sample-a.tsv', { type: 'text/tab-separated-values', lastModified: 2 }), null],
      [null, new File([depthB], 'sample-b.tsv', { type: 'text/tab-separated-values', lastModified: 3 })]
    ]);
    Object.assign(app.form, { multi_record_canvas: true, show_depth: true, suppress_gc: true, suppress_skew: true, legend: 'none' });
    app.adv.depth_tracks.splice(0, app.adv.depth_tracks.length,
      { label: 'Sample A', color: '#1f77b4', large_tick_interval: null, small_tick_interval: null, tick_font_size: null },
      { label: 'Sample B', color: '#d62728', large_tick_interval: null, small_tick_interval: null, tick_font_size: null });
    app.resetCircularTrackSlotsFromSimpleControls();
    app.setCircularTrackSlotsEnabled(true);
    app.sessionTitle = 'replace-follows-record-id';
  }, { a: genbankA, b: genbankB, depthA: depthTsv('BGC0000708', 40579, 10), depthB: depthTsv('BGC0000709', 50466, 50) });
  expect(await evaluateWithRetainedPromise(page, () => window.__GBDRAW_APP__.runAnalysis())).toEqual({ status: 'ok' });
  // A label edit on a CDS of BGC0000708 (record A).
  await evaluateWithRetainedPromise(page, async () => {
    const app = window.__GBDRAW_APP__;
    const { state } = await import('./js/state.js');
    const { setFeatureOverrideRow } = await import('./js/services/drawing/feature-edits.js');
    const feature = app.extractedFeatures.find((entry) => entry.record_id === 'BGC0000708' && entry.type === 'CDS');
    const row = { recordKey: feature.record_key, biologicalFeatureId: feature.biological_feature_id,
      featureVisibility: null, labelVisibility: 'on', labelText: 'EDIT_ON_A', labelSourceText: null };
    await window.__GBDRAW_HISTORY__.runUndoable('Label edit', () => {
      setFeatureOverrideRow(state.activeDrawing().slices.featureEdits, JSON.stringify([row.recordKey, row.biologicalFeatureId]), row);
    });
  });
  const before = await drawn(page);
  expect(before.recordIds).toEqual(['BGC0000708', 'BGC0000709']);
  expect(before.edits).toEqual([['BGC0000708', 'EDIT_ON_A']]);

  // Replace the card's file with the same records in the other order, then save and load before any Generate: the
  // loaded Session keeps no identities of the committed request.
  await replaceCardFile(page, 'ba.gbk', `${genbankB}\n${genbankA}`);
  await expect(page.locator('[data-choice-prompt]')).toHaveCount(0);
  await recordOrder(page, ['BGC0000709', 'BGC0000708']);
  const download = page.waitForEvent('download', { timeout: 120_000 });
  await evaluateWithRetainedPromise(page, () => window.__GBDRAW_APP__.saveSessionWithTitle());
  const saved = await (await download).path();
  await page.reload({ waitUntil: 'domcontentloaded' });
  await waitForAppShell(page, { waitForPalette: false });
  await page.locator('input[accept*=".gbdraw"]').first().setInputFiles(saved);
  await expectSessionLoadedWithoutWarnings(page);
  expect(await evaluateWithRetainedPromise(page, () => window.__GBDRAW_APP__.runAnalysis())).toEqual({ status: 'ok' });
  await recordOrder(page, ['BGC0000709', 'BGC0000708']);
  const loaded = await drawn(page);
  // The figure after the Replace, for review images (GBDRAW_REPLACE_SVG names the file).
  if (process.env.GBDRAW_REPLACE_SVG) {
    require('node:fs').writeFileSync(process.env.GBDRAW_REPLACE_SVG, await page.evaluate(() => window.__GBDRAW_APP__.currentResults()[0].content));
  }
  expect(loaded).toEqual({
    recordIds: ['BGC0000709', 'BGC0000708'],
    depth: { BGC0000709: [null, 'sample-b.tsv'], BGC0000708: ['sample-a.tsv', null] },
    drawnOn: { depth_1: 'BGC0000708', depth_2: 'BGC0000709' },
    edits: [['BGC0000708', 'EDIT_ON_A']]
  });
  expect(await page.evaluate(() => window.__GBDRAW_APP__.currentResults()[0].content)).toContain('EDIT_ON_A');

  // A Replace back to the first order and a Generate in one sitting move them again, from the identities the
  // Generate kept.
  await replaceCardFile(page, 'ab.gbk', `${genbankA}\n${genbankB}`);
  await expect(page.locator('[data-choice-prompt]')).toHaveCount(0);
  await recordOrder(page, ['BGC0000708', 'BGC0000709']);
  expect(await evaluateWithRetainedPromise(page, () => window.__GBDRAW_APP__.runAnalysis())).toEqual({ status: 'ok' });
  const back = await drawn(page);
  expect(back).toEqual({
    recordIds: ['BGC0000708', 'BGC0000709'],
    depth: { BGC0000708: ['sample-a.tsv', null], BGC0000709: [null, 'sample-b.tsv'] },
    drawnOn: { depth_1: 'BGC0000708', depth_2: 'BGC0000709' },
    edits: [['BGC0000708', 'EDIT_ON_A']]
  });
  expect(await page.evaluate(() => window.__GBDRAW_APP__.currentResults()[0].content)).toContain('EDIT_ON_A');
  await testInfo.attach('drawn.json', { body: JSON.stringify({ before, loaded, back }, null, 2), contentType: 'application/json' });
});
