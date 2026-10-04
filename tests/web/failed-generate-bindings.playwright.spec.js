const { test, expect } = require('@playwright/test');
const fs = require('node:fs/promises');
const path = require('node:path');
const { createHash } = require('node:crypto');
const { gunzipSync } = require('node:zlib');
const { openApp, getDiagramWorkerActivity, evaluateWithRetainedPromise } = require('./helpers/app-lifecycle.cjs');

const fixture = path.join(process.cwd(), 'gbdraw/web/gallery/sessions/tobacco-chloroplast.gbdraw-session.json');
const digest = bytes => createHash('sha256').update(bytes).digest('hex');
const input = 'input[accept^=".json,"]';

const load = async (page, file) => {
  page.on('dialog', dialog => dialog.accept());
  await openApp(page);
  await page.locator(input).setInputFiles(file);
  await page.waitForFunction(() => !window.__GBDRAW_APP__.sessionImportPending
    && window.__GBDRAW_APP__.results.length === 1, null, { timeout: 180_000 });
  expect(await page.evaluate(() => window.__GBDRAW_APP__.errorLog)).toBeNull();
};

const save = async (page, file) => {
  const pending = page.waitForEvent('download', { timeout: 120_000 });
  await page.getByRole('banner').getByRole('button', { name: 'Save Session', exact: true }).click();
  await (await pending).saveAs(file);
  return JSON.parse(gunzipSync(await fs.readFile(file)));
};

const exportSvg = async page => {
  const pending = page.waitForEvent('download', { timeout: 120_000 });
  await page.evaluate(() => window.__GBDRAW_APP__.downloadSVG());
  return digest(await fs.readFile(await (await pending).path()));
};

const snapshot = page => page.evaluate(async () => {
  const { state } = await import('./js/state.js');
  const { getCommittedCanonicalRenderRequest } = await import('./js/services/config.js');
  return {
    width: JSON.parse(JSON.stringify(state.adv.circular_track_slots.find(row => row.id === 'gc_content').width)),
    positions: JSON.parse(JSON.stringify(state.adv.multi_record_positions)),
    annotations: state.annotationSets.flatMap(set => set.annotations.map(annotation => ({
      id: annotation.id,
      record: JSON.parse(JSON.stringify(annotation.target?.record)),
      key: annotation.metadata?._gbdraw_web_target_record_key || null
    }))),
    request: getCommittedCanonicalRenderRequest(),
    results: state.results.value.map(result => ({ name: result.name, content: result.content })),
    history: [window.__GBDRAW_HISTORY__.getUndoCount(), window.__GBDRAW_HISTORY__.getRedoCount()]
  };
});

const resourceHashes = document => Object.fromEntries(Object.entries(document.resources || {})
  .map(([id, resource]) => [id, digest(Buffer.from(resource.data, 'base64'))]));

test('failed Generate retains only valid same-source bindings in Save and fresh Load', async ({ page, browser, baseURL }, info) => {
  test.setTimeout(360_000);
  await page.addInitScript(() => { window.__GBDRAW_TEST_HOOKS__ = {}; });
  await load(page, fixture);
  await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    const slot = app.adv.circular_track_slots.find(row => row.id === 'gc_content');
    await window.__GBDRAW_HISTORY__.runUndoable('Pending width', () =>
      app.updateCircularTrackSlotMeasure(slot, 'width', { value: '0.09', unit: 'factor' }));
  });
  const before = await snapshot(page);
  expect(before.positions).toEqual([]);
  expect(before.annotations).toHaveLength(4);
  expect(before.annotations.every(annotation => annotation.key === null)).toBe(true);
  const beforeFile = info.outputPath('before.gbdraw-session.json.gz');
  const beforeSaved = await save(page, beforeFile);
  const oldExport = await exportSvg(page);

  try {
    await page.evaluate(() => {
      window.__GBDRAW_TEST_HOOKS__.beforeDiagramGenerationResponse = () => {
        throw new Error('S04 post-response fault');
      };
    });
    await page.getByRole('button', { name: 'Generate Diagram', exact: true }).click();
    await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.processing),
      { timeout: 180_000 }).toBe(false);
  } finally {
    await page.evaluate(() => { delete window.__GBDRAW_TEST_HOOKS__.beforeDiagramGenerationResponse; });
  }
  expect(await page.evaluate(() => window.__GBDRAW_APP__.errorLog)).toMatchObject({
    code: 'UNKNOWN', operation: 'generate', stage: 'render'
  });
  const failed = await snapshot(page);
  expect(failed.width).toEqual(before.width);
  expect(failed.request).toEqual(before.request);
  expect(failed.results).toEqual(before.results);
  expect(failed.history).toEqual(before.history);
  expect(failed.positions).toEqual([{ selector: '#1', row: 1 }]);
  expect(failed.annotations.map(annotation => annotation.record)).toEqual(
    before.annotations.map(annotation => annotation.record));
  const keys = failed.annotations.map(annotation => annotation.key);
  expect(keys.every(key => typeof key === 'string' && key.endsWith('::[0,"NC_001879.2",155943]'))).toBe(true);
  expect(new Set(keys).size).toBe(1);

  const workerAfterFailure = await getDiagramWorkerActivity(page);
  const failedFile = info.outputPath('failed.gbdraw-session.json.gz');
  const failedSaved = await save(page, failedFile);
  expect(await exportSvg(page)).toBe(oldExport);
  const workerAfterSaveAndExport = await getDiagramWorkerActivity(page);
  expect(workerAfterSaveAndExport.runs).toBe(workerAfterFailure.runs);
  expect(workerAfterSaveAndExport.settledHelpers).toBe(workerAfterSaveAndExport.helpers);
  expect(workerAfterSaveAndExport.settledRuns).toBe(workerAfterSaveAndExport.runs);
  expect(resourceHashes(failedSaved)).toEqual(resourceHashes(beforeSaved));
  expect(failedSaved.renderRequest).toEqual(beforeSaved.renderRequest);
  expect(failedSaved.results).toEqual(beforeSaved.results);
  expect(failedSaved.config.adv.circular_track_slots.find(row => row.id === 'gc_content').width)
    .toEqual(beforeSaved.config.adv.circular_track_slots.find(row => row.id === 'gc_content').width);
  expect(failedSaved.config.adv.multi_record_positions).toEqual(failed.positions);
  expect(failedSaved.config.annotationSets.flatMap(set => set.annotations.map(annotation =>
    annotation.metadata?._gbdraw_web_target_record_key))).toEqual(keys);

  const freshContext = await browser.newContext({ baseURL, acceptDownloads: true });
  try {
    const fresh = await freshContext.newPage();
    await load(fresh, failedFile);
    const restored = await snapshot(fresh);
    expect(restored.width).toEqual(before.width);
    expect(restored.request).toEqual(before.request);
    expect(restored.results).toEqual(before.results);
    expect(restored.positions).toEqual(failed.positions);
    expect(restored.annotations).toEqual(failed.annotations);
    expect(await exportSvg(fresh)).toBe(oldExport);
    const resaved = await save(fresh, info.outputPath('reloaded.gbdraw-session.json.gz'));
    expect(resourceHashes(resaved)).toEqual(resourceHashes(beforeSaved));
    expect(resaved.renderRequest).toEqual(beforeSaved.renderRequest);
    expect(resaved.results).toEqual(beforeSaved.results);
  } finally {
    await freshContext.close();
  }

  await page.getByRole('button', { name: 'Generate Diagram', exact: true }).click();
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.processing),
    { timeout: 180_000 }).toBe(false);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.errorLog)).toBeNull();
  expect((await snapshot(page)).history[0]).toBe(before.history[0] + 1);
});

test('canceled post-response Generate keeps valid bindings and old artifact through stale completion', async ({ page, browser, baseURL }, info) => {
  test.setTimeout(300_000);
  await page.addInitScript(() => { window.__GBDRAW_TEST_HOOKS__ = {}; });
  await load(page, fixture);
  const before = await snapshot(page);
  const oldExport = await exportSvg(page);
  await page.evaluate(() => {
    window.__GBDRAW_TEST_HOOKS__.beforeDiagramGenerationResponse = () => new Promise(resolve => {
      window.releaseS04Response = resolve;
    });
  });
  try {
    await page.getByRole('button', { name: 'Generate Diagram', exact: true }).click();
    await page.waitForFunction(() => Boolean(window.releaseS04Response), null, { timeout: 180_000 });
    await page.getByRole('button', { name: /Cancel$/ }).click();
    await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.processing)).toBe(false);
  } finally {
    await page.evaluate(() => {
      delete window.__GBDRAW_TEST_HOOKS__.beforeDiagramGenerationResponse;
      window.releaseS04Response?.();
    });
  }
  const canceled = await snapshot(page);
  expect(canceled.request).toEqual(before.request);
  expect(canceled.results).toEqual(before.results);
  expect(canceled.width).toEqual(before.width);
  expect(canceled.history).toEqual(before.history);
  expect(canceled.positions).toEqual([{ selector: '#1', row: 1 }]);
  expect(canceled.annotations.every(annotation => annotation.record?.value === 'NC_001879.2'
    && annotation.key?.endsWith('::[0,"NC_001879.2",155943]'))).toBe(true);
  const canceledFile = info.outputPath('canceled.gbdraw-session.json.gz');
  const saved = await save(page, canceledFile);
  expect(saved.renderRequest).toEqual(before.request);
  expect(saved.config.adv.multi_record_positions).toEqual(canceled.positions);
  expect(await exportSvg(page)).toBe(oldExport);
  const freshContext = await browser.newContext({ baseURL, acceptDownloads: true });
  try {
    const fresh = await freshContext.newPage();
    await load(fresh, canceledFile);
    const restored = await snapshot(fresh);
    expect(restored.request).toEqual(before.request);
    expect(restored.results).toEqual(before.results);
    expect(restored.positions).toEqual(canceled.positions);
    expect(restored.annotations).toEqual(canceled.annotations);
  } finally {
    await freshContext.close();
  }
  await expect.poll(async () => {
    const activity = await getDiagramWorkerActivity(page);
    return [activity.settledHelpers === activity.helpers, activity.settledRuns === activity.runs];
  }).toEqual([true, true]);
});

const genbankRecord = id => {
  const sequence = 'atg'.repeat(100);
  const origin = sequence.match(/.{1,60}/g).map((line, index) =>
    `${String(index * 60 + 1).padStart(9)} ${line.match(/.{1,10}/g).join(' ')}`).join('\n');
  return `LOCUS       ${id.padEnd(24)} 300 bp    DNA     linear   UNA 01-JAN-2000
DEFINITION  S04 public synthetic record.
ACCESSION   ${id}
VERSION     ${id}
KEYWORDS    .
SOURCE      synthetic construct
  ORGANISM  synthetic construct
            .
FEATURES             Location/Qualifiers
     CDS             1..90
                     /product="test protein"
ORIGIN
${origin}
//
`;
};

test('failed multi-record Generate retains a validated binding to the selected second record', async ({ page, browser, baseURL }, info) => {
  test.setTimeout(360_000);
  page.on('dialog', dialog => dialog.accept(dialog.type() === 'prompt' ? 'S04 multi' : undefined));
  await page.addInitScript(() => { window.__GBDRAW_TEST_HOOKS__ = {}; });
  await openApp(page);
  await page.evaluate(content => {
    const app = window.__GBDRAW_APP__;
    app.files.c_gb = new File([content], 'two-records.gb', { type: 'text/plain', lastModified: 1 });
    Object.assign(app.form, { multi_record_canvas: true, suppress_gc: true, suppress_skew: true, labels_mode: 'none' });
    const set = app.addAnnotationSet('s04');
    const annotation = app.addCoordinateAnnotation(set, { start: 1, end: 10 });
    annotation.target.record = { kind: 'recordId', value: 'RecB' };
  }, genbankRecord('RecA') + genbankRecord('RecB'));
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.circularRecordList.length),
    { timeout: 60_000 }).toBe(2);
  expect(await evaluateWithRetainedPromise(page, () => window.__GBDRAW_APP__.runAnalysis())).toEqual({ status: 'ok' });
  await page.evaluate(() => {
    delete window.__GBDRAW_APP__.annotationSets[0].annotations[0].metadata._gbdraw_web_target_record_key;
  });
  const beforeSaved = await save(page, info.outputPath('multi-before.gbdraw-session.json.gz'));
  const before = await snapshot(page);
  expect(before.positions.map(position => position.selector)).toEqual(['#1', '#2']);
  expect(before.annotations[0].record).toEqual({ kind: 'recordId', value: 'RecB' });
  expect(before.annotations[0].key).toBeNull();
  try {
    await page.evaluate(() => {
      window.__GBDRAW_TEST_HOOKS__.beforeDiagramGenerationResponse = () => {
        throw new Error('S04 multi-record post-response fault');
      };
    });
    await page.getByRole('button', { name: 'Generate Diagram', exact: true }).click();
    await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.processing),
      { timeout: 180_000 }).toBe(false);
  } finally {
    await page.evaluate(() => { delete window.__GBDRAW_TEST_HOOKS__.beforeDiagramGenerationResponse; });
  }
  const failed = await snapshot(page);
  expect(failed.request).toEqual(before.request);
  expect(failed.results).toEqual(before.results);
  expect(failed.history).toEqual(before.history);
  expect(failed.positions.map(position => position.selector)).toEqual(['#1', '#2']);
  expect(failed.annotations[0].record).toEqual(before.annotations[0].record);
  expect(failed.annotations[0].key).toMatch(/::\[1,"RecB",300\]$/);
  const failedFile = info.outputPath('multi-failed.gbdraw-session.json.gz');
  const failedSaved = await save(page, failedFile);
  expect(resourceHashes(failedSaved)).toEqual(resourceHashes(beforeSaved));
  expect(failedSaved.renderRequest).toEqual(beforeSaved.renderRequest);
  expect(failedSaved.results).toEqual(beforeSaved.results);
  const freshContext = await browser.newContext({ baseURL, acceptDownloads: true });
  try {
    const fresh = await freshContext.newPage();
    await load(fresh, failedFile);
    const restored = await snapshot(fresh);
    expect(restored.request).toEqual(before.request);
    expect(restored.results).toEqual(before.results);
    expect(restored.positions).toEqual(failed.positions);
    expect(restored.annotations).toEqual(failed.annotations);
  } finally {
    await freshContext.close();
  }
});
