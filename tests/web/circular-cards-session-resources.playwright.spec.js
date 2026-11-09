const { test, expect } = require('@playwright/test');
const fs = require('node:fs/promises');
const { createHash } = require('node:crypto');
const { gunzipSync } = require('node:zlib');
const { execFile } = require('node:child_process');
const { promisify } = require('node:util');
const path = require('node:path');
const { CURRENT_REQUEST_SCHEMA, CURRENT_SESSION_VERSION, answerChoicePrompt, answerLoadRepairs, openApp, reveal } = require('./helpers/app-lifecycle.cjs');
const { shownSavedDrawing } = require('./helpers/saved-session.cjs');

// A CLI Session of three single-record GenBank files (`--gbk a b c`). The
// Session 44 bound them as one composite `c_gb`; Load makes one Circular File
// card of each file, and every later Save binds one resource per card
// (`webFiles.bindings.circularSources`).
const seed = 'tests/fixtures/sessions/composite-circular-three-files.v44-schema8.gbdraw-session.json.gz';
const recordCount = 3;
const loadTimeout = 300_000;
const generateTimeout = 600_000;
const hash = bytes => createHash('sha256').update(bytes).digest('hex');
const resourceIdentity = descriptor => ({
  size: descriptor.size,
  sha256: hash(Buffer.from(descriptor.data, 'base64'))
});
const cardLocators = page => page.locator('[data-circular-source-card]');
// A card's file as the page shows it, without the resource id a Save may rename.
const withoutResourceId = cards => cards.map(({ resourceId, ...card }) => card);
// A saved card binding without the resource id of each of its files.
const savedCardShape = card => Object.fromEntries(Object.entries(card).map(([key, value]) => [key,
  value && typeof value === 'object' ? (({ resourceId, ...metadata }) => metadata)(value) : value]));

const load = async (browser, file, contexts, viewport = { width: 1600, height: 1000 }) => {
  const context = await browser.newContext({ viewport });
  contexts.push(context);
  const external = [];
  await context.route('**/*', route => {
    if (new URL(route.request().url()).hostname === '127.0.0.1') return route.continue();
    external.push(route.request().url());
    return route.abort();
  });
  const page = await context.newPage();
  page.setDefaultTimeout(loadTimeout);
  // Accept both the large-session download confirmation and repeated filename confirmation.
  page.on('dialog', dialog => dialog.accept());
  await page.addInitScript(() => {
    window.__CARD_METRICS__ = [];
    window.__CARD_PREFLIGHT_METRICS__ = null;
    window.__GBDRAW_TEST_HOOKS__ = {
      onStructuralMetric: metric => window.__CARD_METRICS__.push(metric),
      onSessionLifecycleEvent: event => {
        if (event.name === 'current-session-preflight-end') {
          window.__CARD_PREFLIGHT_METRICS__ = [...window.__CARD_METRICS__];
        }
      }
    };
  });
  await openApp(page);
  await page.locator('input[accept^=".json,"]').setInputFiles(file);
  await expect.poll(() => page.evaluate(() => !window.__GBDRAW_APP__.sessionImportPending
    && window.__GBDRAW_APP__.extractedFeatures.length > 0), { timeout: loadTimeout }).toBe(true);
  expect(external).toEqual([]);
  const genomeReads = await page.evaluate(async () => {
    const { state } = await import('./js/state.js');
    const { getSessionResourceSource } = await import('./js/services/file-content-cache.js');
    const ids = state.activeDrawing().circularSources
      .map(card => getSessionResourceSource(state.inputFiles.fileOf(card.gb))?.resourceId).filter(Boolean);
    // The later definition-preview helper legitimately consumes source bytes.
    // Observe metadata projection at its boundary, independent of trace overhead.
    return window.__CARD_PREFLIGHT_METRICS__?.filter(metric => metric.name === 'base64DecodeCount' && ids.includes(metric.resourceId));
  });
  expect(genomeReads).toEqual([]);
  return { page, external };
};

const snapshot = page => page.evaluate(async () => {
  const { state } = await import('./js/state.js');
  const { getSessionResourceSource, readFileBytes } = await import('./js/services/file-content-cache.js');
  const { getCommittedCanonicalSession } = await import('./js/services/config.js');
  const digest = async bytes => [...new Uint8Array(await crypto.subtle.digest('SHA-256', bytes))]
    .map(byte => byte.toString(16).padStart(2, '0')).join('');
  // A committed resource's bytes, read by their owner: a Session's text, or the File Generate described.
  const backing = await import('./js/services/session-resource-backing.js');
  const identity = async descriptor => ({ size: descriptor.size, sha256: await digest(await backing.readSessionResourceBytes(
    backing.createSessionResourceFileView(backing.adoptCurrentSessionResources({ resource: descriptor }), 'resource'))) });
  const drawing = state.activeDrawing();
  const committed = getCommittedCanonicalSession();
  return {
    // One File card per file, each with its own GenBank File and no other file.
    cards: await Promise.all(drawing.circularSources.map(async card => {
      const file = state.inputFiles.fileOf(card.gb);
      return { uid: card.uid, gff: card.gff, fasta: card.fasta, name: file.name, size: file.size, type: file.type,
        lastModified: file.lastModified, sha256: await digest(await readFileBytes(file)),
        resourceId: getSessionResourceSource(file)?.resourceId ?? null };
    })),
    records: drawing.circularRecordList.value.map(record => ({ uid: record.uid, identity: record.identity, recordId: record.record_id })),
    request: committed.renderRequest,
    resources: Object.fromEntries(await Promise.all(Object.entries(committed.resources)
      .map(async ([id, descriptor]) => [id, await identity(descriptor)]))),
    results: window.__GBDRAW_APP__.currentResults().map(result => result.content),
    resultNames: state.results.value.map(result => result.name),
    error: state.errorLog.value
  };
});

const generate = async page => {
  const key = await page.evaluate(async () => (await import('./js/state.js')).state.resultGenerationKey.value);
  await page.getByRole('button', { name: 'Generate Diagram', exact: true }).click();
  await expect.poll(() => page.evaluate(async () => {
    const { state } = await import('./js/state.js');
    return { key: state.resultGenerationKey.value, processing: state.processing.value, error: state.errorLog.value };
  }), { timeout: generateTimeout }).toEqual({ key: key + 1, processing: false, error: null });
};

// `again`: the page already downloaded this Session file name, so the Save asks first.
const save = async (page, testInfo, label, { again = false } = {}) => {
  let download = null;
  const onDownload = value => { download = value; };
  page.on('download', onDownload);
  try {
    await page.evaluate(() => { window.__CARD_METRICS__ = []; });
    const started = performance.now();
    await page.getByRole('button', { name: 'Save Session', exact: true }).click();
    if (again) await answerChoicePrompt(page, 'Download the Session again?', 'download');
    await expect.poll(async () => download ? 'download' : await page.evaluate(() =>
      window.__GBDRAW_APP__.errorLog || 'pending'), { timeout: loadTimeout }).not.toBe('pending');
    expect(await page.evaluate(() => window.__GBDRAW_APP__.errorLog)).toBeNull();
    expect(download).not.toBeNull();
    const saved = testInfo.outputPath(`${label}.gbdraw-session.json.gz`);
    await download.saveAs(saved);
    expect(await download.failure()).toBeNull();
    const bytes = await fs.readFile(saved);
    const session = JSON.parse((bytes[0] === 0x1f ? gunzipSync(bytes) : bytes).toString());
    const metrics = await page.evaluate(() => window.__CARD_METRICS__);
    await fs.writeFile(testInfo.outputPath(`${label}-metrics.json`), JSON.stringify({
      milliseconds: performance.now() - started, gzipBytes: bytes.length,
      jsonBytes: Buffer.byteLength(JSON.stringify(session)),
      resourceBytes: Object.values(session.resources).reduce((sum, r) => sum + r.size, 0), metrics
    }, null, 2));
    return { saved, session, metrics };
  } finally {
    page.off('download', onDownload);
  }
};

const readSeed = async () => {
  const bytes = await fs.readFile(seed);
  expect(hash(bytes)).toBe('87fef00880b1f9f2cd81851f84049e75076a8be1db2e2d0523dadb592e2e7312');
  const session = JSON.parse(gunzipSync(bytes));
  const components = session.webFiles.bindings.c_gb.components;
  expect(components.map(component => component.name))
    .toEqual(['SARS-CoV-1.gbk', 'SARS-CoV-2.gbk', 'MERS-CoV.gbk']);
  return { session, components, componentIds: components.map(component => component.resourceId) };
};

// The saved File cards: one row per card, each GenBank file a leaf whose
// resource holds that card's bytes.
const savedCards = (session, cards) => {
  const rows = shownSavedDrawing(session).webFiles.bindings.circularSources;
  expect(rows.map(row => row.uid)).toEqual(cards.map(card => card.uid));
  expect(rows.map(row => [row.gff, row.fasta])).toEqual(cards.map(() => [null, null]));
  expect(rows.map(row => resourceIdentity(session.resources[row.gb.resourceId])))
    .toEqual(cards.map(({ size, sha256 }) => ({ size, sha256 })));
  expect(rows.map(({ gb: { resourceId, ...metadata } }) => metadata))
    .toEqual(cards.map(({ name, type, lastModified }) => ({ name, type, lastModified })));
  return rows;
};

for (const journey of ['minimal', 'grid-batch-grid']) {
  test(`Circular File cards keep one resource per file through ${journey} Save, fresh Load and Generate`, async ({ browser }, testInfo) => {
    test.setTimeout(6 * generateTimeout + 6 * loadTimeout);
    const contexts = [];
    try {
      const { session: seedSession, components, componentIds: seedIds } = await readSeed();
      const { page, external } = await load(browser, seed, contexts);
      const original = await snapshot(page);
      expect(original.request.records).toHaveLength(recordCount);
      // Load makes one card of each file the Session bound, in its order, with its name and bytes.
      expect(original.cards.map(({ name, type, lastModified, gff, fasta }) => ({ name, type, lastModified, gff, fasta })))
        .toEqual(components.map(({ name, type, lastModified }) => ({ name, type, lastModified, gff: null, fasta: null })));
      expect(original.cards.map(({ size, sha256 }) => ({ size, sha256 })))
        .toEqual(seedIds.map(id => resourceIdentity(seedSession.resources[id])));
      expect(new Set(original.cards.map(card => card.uid)).size).toBe(recordCount);
      await expect(cardLocators(page)).toHaveCount(recordCount);
      if (journey === 'minimal') {
        const control = await save(page, testInfo, 'before-generate');
        expect(shownSavedDrawing(control.session).results).toHaveLength(1);
        savedCards(control.session, original.cards);
      }
      await generate(page);
      if (journey === 'grid-batch-grid') {
        const multiRecordCanvas = await reveal(page.getByRole('checkbox', { name: 'Multi-Record Canvas', exact: true, includeHidden: true }));
        await multiRecordCanvas.uncheck();
        await generate(page);
        expect((await snapshot(page)).results).toHaveLength(recordCount);
        await multiRecordCanvas.check();
        await generate(page);
      }
      const generated = await snapshot(page);
      // One record of each card, keyed by the card (a one-record file's identity is its card's uid).
      expect(generated.records.map(record => record.identity)).toEqual(original.cards.map(card => card.uid));
      expect(generated.request.records.map(record => record.recordKey)).toEqual(original.cards.map(card => card.uid));
      expect(generated.results).toHaveLength(1);
      expect(generated.cards).toEqual(original.cards);
      const { saved, session, metrics } = await save(page, testInfo, 'generated', { again: journey === 'minimal' });
      expect([session.version, shownSavedDrawing(session).webFiles.bindings.schema, shownSavedDrawing(session).renderRequest.schema]).toEqual([CURRENT_SESSION_VERSION, 2, CURRENT_REQUEST_SCHEMA]);
      const bindings = shownSavedDrawing(session).webFiles.bindings;
      expect(['c_gb', 'c_gff', 'c_fasta'].filter(slot => Object.hasOwn(bindings, slot))).toEqual([]);
      const rows = savedCards(session, original.cards);
      // Save decodes nothing, and encodes no File Generate described more than once (OV-2612).
      expect(metrics.filter(metric => metric.name === 'base64DecodeCount')).toEqual([]);
      const encoded = metrics.filter(metric => metric.name === 'base64EncodeCount').map(metric => metric.resourceName);
      expect(new Set(encoded).size).toBe(encoded.length);
      const resources = Object.fromEntries(Object.entries(session.resources)
        .map(([id, descriptor]) => [id, resourceIdentity(descriptor)]));
      for (const [id, identity] of Object.entries(generated.resources)) expect(resources[id]).toEqual(identity);
      expect(shownSavedDrawing(session).renderRequest).toEqual(generated.request);
      await fs.writeFile(testInfo.outputPath('resource-inspection.json'), JSON.stringify({
        original: { ...original, results: original.results.map(hash) },
        generated: { ...generated, results: generated.results.map(hash) },
        resources, webFiles: shownSavedDrawing(session).webFiles
      }, null, 2));
      expect(external).toEqual([]);
      await page.context().close();

      const restored = await load(browser, saved, contexts);
      const loaded = await snapshot(restored.page);
      expect(withoutResourceId(loaded.cards)).toEqual(withoutResourceId(original.cards));
      expect(loaded.request).toEqual(generated.request);
      expect(loaded.resultNames).toEqual(generated.resultNames);
      expect(loaded.results).toEqual(shownSavedDrawing(session).results.map(result => result.content));
      const resaved = await save(restored.page, testInfo, 'restored');
      expect(shownSavedDrawing(resaved.session).webFiles.bindings.circularSources.map(savedCardShape))
        .toEqual(rows.map(savedCardShape));
      savedCards(resaved.session, original.cards);
      await generate(restored.page);
      const regenerated = await snapshot(restored.page);
      expect(regenerated.records).toEqual(generated.records);
      expect(regenerated.results).toEqual(generated.results);
      expect(regenerated.cards).toEqual(loaded.cards);
      await restored.page.screenshot({ path: testInfo.outputPath('regenerated.png') });
      if (journey === 'minimal') {
        // Replace the first card's file: that card keeps its place and uid, the others are unchanged.
        const first = cardLocators(restored.page).first();
        await first.locator('summary').click({ position: { x: 8, y: 8 } });
        await first.locator('input[type=file]').setInputFiles({
          name: 'replacement.gb', mimeType: 'text/plain',
          buffer: Buffer.from(seedSession.resources[seedIds[0]].data, 'base64')
        });
        const replaced = await save(restored.page, testInfo, 'replaced', { again: true });
        const replacedRows = shownSavedDrawing(replaced.session).webFiles.bindings.circularSources;
        expect(replacedRows.map(row => row.uid)).toEqual(original.cards.map(card => card.uid));
        expect(replacedRows[0].gb.name).toBe('replacement.gb');
        expect(resourceIdentity(replaced.session.resources[replacedRows[0].gb.resourceId]))
          .toEqual(resourceIdentity(seedSession.resources[seedIds[0]]));
        expect(replacedRows.slice(1).map(savedCardShape)).toEqual(rows.slice(1).map(savedCardShape));
        const mobile = await load(browser, saved, contexts, { width: 390, height: 844 });
        const mobileState = await snapshot(mobile.page);
        expect(mobileState.cards).toEqual(loaded.cards);
        expect(mobileState.request).toEqual(generated.request);
        expect(mobileState.results).toEqual(loaded.results);
        await mobile.page.screenshot({ path: testInfo.outputPath('mobile-restored.png') });
        expect(mobile.external).toEqual([]);
      }
      expect(restored.external).toEqual([]);
    } finally {
      for (const context of contexts) await context.close();
    }
  });
}

test('Session CLI sidecar keeps the Circular File cards for fresh Web Load and Generate', async ({ browser }, testInfo) => {
  test.setTimeout(4 * generateTimeout);
  const contexts = [];
  try {
    await readSeed();
    const { page } = await load(browser, seed, contexts);
    await generate(page);
    const original = await snapshot(page);
    const { saved, session } = await save(page, testInfo, 'cli-source');
    await page.context().close();
    const prefix = testInfo.outputPath('cli-replayed');
    const sidecar = `${prefix}.json.gz`;
    const { stdout, stderr } = await promisify(execFile)('python', [
      '-m', 'gbdraw.cli', 'circular', '--session', saved,
      '-o', prefix, '--session_output', sidecar
    ], {
      cwd: testInfo.outputDir,
      env: { ...process.env, PYTHONPATH: process.cwd() },
      timeout: generateTimeout, maxBuffer: 1_000_000
    });
    await fs.writeFile(testInfo.outputPath('cli.log'), stdout + stderr);
    const replayed = JSON.parse(gunzipSync(await fs.readFile(sidecar)));
    expect([replayed.version, shownSavedDrawing(replayed).webFiles.bindings.schema, shownSavedDrawing(replayed).renderRequest.schema]).toEqual([CURRENT_SESSION_VERSION, 2, CURRENT_REQUEST_SCHEMA]);
    // The replay keeps each card, its uid and its file's name and bytes.
    const expected = savedCards(session, original.cards);
    const actual = savedCards(replayed, original.cards);
    expect(actual.map(savedCardShape)).toEqual(expected.map(savedCardShape));
    expect(shownSavedDrawing(replayed).results[0].content).toBe(await fs.readFile(`${prefix}.svg`, 'utf8'));
    expect(shownSavedDrawing(replayed).results[0].name).toBe(path.basename(prefix));
    expect(shownSavedDrawing(replayed).editorState.featureCatalog.items[0].resultName).toBe(path.basename(prefix));
    const restored = await load(browser, sidecar, contexts);
    const loaded = await snapshot(restored.page);
    expect(withoutResourceId(loaded.cards)).toEqual(withoutResourceId(original.cards));
    const admittedContents = await restored.page.evaluate(async results => {
      const { sanitizeSvgContent } = await import('./js/services/svg-sanitization.js');
      return results.map(result => sanitizeSvgContent(result.content, window.DOMPurify));
    }, shownSavedDrawing(replayed).results);
    expect(loaded.results.map(hash)).toEqual(admittedContents.map(hash));
    expect(loaded.request).toEqual(shownSavedDrawing(replayed).renderRequest);
    await generate(restored.page);
    expect((await snapshot(restored.page)).records).toEqual(original.records);
    expect(restored.external).toEqual([]);
    await restored.page.screenshot({ path: testInfo.outputPath('cli-regenerated.png') });
  } finally {
    for (const context of contexts) await context.close();
  }
});

// BA-D25 (receipt 〔確認〕): a malformed binding is not a refusal class; Load lists
// what it drops before it changes anything, and Cancel keeps the current work.
test('Session import lists malformed File card bindings and Cancel keeps existing work @pr-smoke', async ({ browser }, testInfo) => {
  test.setTimeout(6 * loadTimeout);
  const contexts = [];
  try {
    const seedPath = 'gbdraw/web/gallery/sessions/HmmtDNA_basic_circular.gbdraw-session.json';
    const { page, external } = await load(browser, seedPath, contexts);
    const before = await snapshot(page);
    const beforeHistory = await page.evaluate(()=>window.__GBDRAW_HISTORY__.getUndoCount());
    const seedSession = JSON.parse(await fs.readFile(seedPath));
    const seedDrawing = shownSavedDrawing(seedSession);
    const resourceId = seedDrawing.renderRequest.records[0].source.resourceId;
    const at = `drawings[${seedSession.drawings.indexOf(seedDrawing)}]`;
    const leaf = { resourceId, name: 'part.gb', type: 'text/plain', lastModified: 0 };
    const card = { uid: seedDrawing.renderRequest.records[0].recordKey, gb: leaf, gff: null, fasta: null };

    const inventoryDropped = [`Session file: the field ${at}.webFiles.bindings could not be read and was left out.`];
    const entryDropped = place => [`Input files: 1 entry of ${at}.webFiles.bindings.circularSources${place} could not be read and was dropped.`];
    const cards = s => shownSavedDrawing(s).webFiles.bindings.circularSources;
    const mutations = {
      schema: s => { shownSavedDrawing(s).webFiles.bindings.schema = 99; },
      missing: s => { delete shownSavedDrawing(s).webFiles.bindings.circularSources; },
      list: s => { shownSavedDrawing(s).webFiles.bindings.circularSources = { 0: structuredClone(card) }; },
      field: s => { cards(s)[0].depth = null; },
      uid: s => { cards(s)[0].uid = ` ${card.uid}`; },
      repeated: s => { cards(s).push(structuredClone(card)); },
      null: s => { cards(s)[0] = null; },
      composite: s => { cards(s)[0].gb = { kind: 'composite', components: [structuredClone(leaf), structuredClone(leaf)],
        name: 'logical.gb', type: 'text/plain', lastModified: 0 }; },
      dangling: s => { cards(s)[0].gb.resourceId = 'missing'; },
      metadata: s => { cards(s)[0].gb.lastModified = -1; },
      // The slot the Circular input had before the File cards is no field of a Session 51 drawing.
      retired: s => { shownSavedDrawing(s).webFiles.bindings.c_gb = structuredClone(leaf); },
      payload: s => { s.resources[resourceId].encoding = 'unsupported'; }
    };
    const rows = {
      schema: inventoryDropped, missing: inventoryDropped, list: inventoryDropped,
      field: entryDropped('[0]'), uid: entryDropped('[0]'), repeated: entryDropped('[1]'), null: entryDropped('[0]'),
      retired: [`Session file: the field ${at}.webFiles.bindings.c_gb could not be read and was left out.`],
      payload: [
        `Saved diagram: the settings of the saved diagram (${at}.renderRequest) could not be read. The file loads as settings only; the saved Results, the GenBank, GFF and FASTA inputs, the feature list and the caches are left out.`,
        `Input files: 1 entry of resources.${resourceId} could not be read and was dropped.`,
        ...entryDropped('[0].gb')
      ]
    };
    const canceled = [];
    for (const [label, mutate] of Object.entries(mutations)) {
      const invalid = structuredClone(seedSession);
      const invalidDrawing = shownSavedDrawing(invalid);
      invalidDrawing.webFiles ||= {};
      invalidDrawing.webFiles.bindings = { schema: 2, circularSources: [structuredClone(card)] };
      mutate(invalid);
      await page.locator('input[accept^=".json,"]').setInputFiles({
        name: `${label}.gbdraw-session.json`, mimeType: 'application/json', buffer: Buffer.from(JSON.stringify(invalid))
      });
      // The rows name places in the file, never the stored file names.
      await answerLoadRepairs(page, rows[label] || entryDropped('[0].gb'), 'cancel', { timeout: loadTimeout });
      await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.sessionImportPending), { timeout: loadTimeout }).toBe(false);
      expect(await snapshot(page)).toEqual(before);
      expect(await page.evaluate(()=>window.__GBDRAW_HISTORY__.getUndoCount())).toBe(beforeHistory);
      canceled.push(label);
    }
    expect(external).toEqual([]);
    await fs.writeFile(testInfo.outputPath('transactional-cancels.json'), JSON.stringify(canceled, null, 2));
  } finally {
    for (const context of contexts) await context.close();
  }
});
