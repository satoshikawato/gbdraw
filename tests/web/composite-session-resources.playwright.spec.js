const { test, expect } = require('@playwright/test');
const fs = require('node:fs/promises');
const { createHash } = require('node:crypto');
const { gunzipSync } = require('node:zlib');
const { execFile } = require('node:child_process');
const { promisify } = require('node:util');
const path = require('node:path');
const { openApp } = require('./helpers/app-lifecycle.cjs');

const seed = 'gbdraw/web/gallery/sessions/Vnig_TUMSAT-TG-2018.gbdraw-session.json.gz';
const loadTimeout = 300_000;
const generateTimeout = 1_800_000;
const hash = bytes => createHash('sha256').update(bytes).digest('hex');
const resourceIdentity = descriptor => ({
  size: descriptor.size,
  sha256: hash(Buffer.from(descriptor.data, 'base64'))
});

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
    window.__COMPOSITE_METRICS__ = [];
    window.__COMPOSITE_PREFLIGHT_METRICS__ = null;
    window.__GBDRAW_TEST_HOOKS__ = {
      onStructuralMetric: metric => window.__COMPOSITE_METRICS__.push(metric),
      onSessionLifecycleEvent: event => {
        if (event.name === 'current-session-preflight-end') {
          window.__COMPOSITE_PREFLIGHT_METRICS__ = [...window.__COMPOSITE_METRICS__];
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
    const source = getSessionResourceSource(state.files.c_gb);
    const ids = (source?.descriptors || (source?.descriptor ? [source] : [])).map(part => part.resourceId);
    // The later definition-preview helper legitimately consumes source bytes.
    // Observe metadata projection at its boundary, independent of trace overhead.
    return window.__COMPOSITE_PREFLIGHT_METRICS__?.filter(metric => metric.name === 'base64DecodeCount' && ids.includes(metric.resourceId));
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
  const identity = async descriptor => ({ size: descriptor.size,
    sha256: await digest(Uint8Array.from(atob(descriptor.data), char => char.charCodeAt(0))) });
  const source = getSessionResourceSource(state.files.c_gb);
  const components = source?.descriptors || (source?.descriptor ? [source] : []);
  const committed = getCommittedCanonicalSession();
  return {
    components: await Promise.all(components.map(async ({ resourceId, descriptor, name, type, lastModified }) => ({
      resourceId, name, type, lastModified, ...await identity(descriptor)
    }))),
    source: { name: state.files.c_gb.name, size: state.files.c_gb.size,
      type: state.files.c_gb.type, lastModified: state.files.c_gb.lastModified,
      isArray: Array.isArray(state.files.c_gb),
      sha256: await digest(await readFileBytes(state.files.c_gb)) },
    records: state.circularRecordList.value.map(record => ({ id: record.id, recordKey: record.recordKey })),
    request: committed.renderRequest,
    resources: Object.fromEntries(await Promise.all(Object.entries(committed.resources)
      .map(async ([id, descriptor]) => [id, await identity(descriptor)]))),
    results: state.results.value.map(result => result.content),
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

const save = async (page, testInfo, label) => {
  let download = null;
  const onDownload = value => { download = value; };
  page.on('download', onDownload);
  try {
    await page.evaluate(() => { window.__COMPOSITE_METRICS__ = []; });
    const started = performance.now();
    await page.getByRole('button', { name: 'Save Session', exact: true }).click();
    await expect.poll(async () => download ? 'download' : await page.evaluate(() =>
      window.__GBDRAW_APP__.errorLog || 'pending'), { timeout: loadTimeout }).not.toBe('pending');
    expect(await page.evaluate(() => window.__GBDRAW_APP__.errorLog)).toBeNull();
    expect(download).not.toBeNull();
    const saved = testInfo.outputPath(`${label}.gbdraw-session.json.gz`);
    await download.saveAs(saved);
    expect(await download.failure()).toBeNull();
    const bytes = await fs.readFile(saved);
    const session = JSON.parse((bytes[0] === 0x1f ? gunzipSync(bytes) : bytes).toString());
    const metrics = await page.evaluate(() => window.__COMPOSITE_METRICS__);
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

const prepareSeed = async (journey, testInfo) => {
  const bytes = await fs.readFile(seed);
  expect(hash(bytes)).toBe('11009dd28a1d98afa79e8b5d22b94217f61d9f53cb6353010095c983c3a9d90e');
  const session = JSON.parse(gunzipSync(bytes));
  if (journey !== 'minimal') return { file: seed, session };

  // Three real plasmids exercise composite ordering and persistence in PR smoke.
  // Keep the six-replicon Gallery fixture for grid/batch and CLI acceptance.
  const records = session.renderRequest.records.slice(2, 5)
    .map((record, index) => ({ ...record, recordKey: `record-${index + 1}` }));
  const resourceIds = new Set(records.map(record => record.source.resourceId));
  const title = 'Three Vibrio plasmids';
  const input = testInfo.outputPath('three-plasmids-input.json');
  const file = testInfo.outputPath('three-plasmids.gbdraw-session.json.gz');
  await fs.writeFile(input, JSON.stringify({
    format: session.format, version: session.version,
    results: [], editorState: { featureCatalog: null },
    ui: session.ui,
    config: { ...session.config,
      form: { ...session.config.form, plot_title: title },
      adv: { ...session.config.adv, multi_record_positions: [] } },
    renderRequest: { ...session.renderRequest, records,
      layout: { ...session.renderRequest.layout, multiRecordPositions: [] },
      diagramOptions: { ...session.renderRequest.diagramOptions, plotTitle: title } },
    resources: Object.fromEntries(Object.entries(session.resources)
      .filter(([id, resource]) => resource.kind !== 'genbank' || resourceIds.has(id)))
  }));
  // Materialize the matching SVG and catalog through the normal CLI writer.
  const { stdout, stderr } = await promisify(execFile)('python', [
    '-m', 'gbdraw.cli', 'circular', '--session', input,
    '-o', testInfo.outputPath('three-plasmids'), '--session_output', file
  ], { cwd: testInfo.outputDir, env: { ...process.env, PYTHONPATH: process.cwd() },
    timeout: generateTimeout, maxBuffer: 1_000_000 });
  await fs.writeFile(testInfo.outputPath('three-plasmids-cli.log'), stdout + stderr);
  return { file, session: JSON.parse(gunzipSync(await fs.readFile(file))) };
};

for (const journey of ['minimal', 'grid-batch-grid']) {
  test(`Session export Vibrio composite resources survive ${journey} Save, fresh Load and Generate`, async ({ browser }, testInfo) => {
    test.setTimeout(6 * generateTimeout + 6 * loadTimeout);
    const contexts = [];
    try {
      const { file, session: seedSession } = await prepareSeed(journey, testInfo);
      const recordCount = journey === 'minimal' ? 3 : 6;
      const { page, external } = await load(browser, file, contexts);
      const original = await snapshot(page);
      expect(original.request.records).toHaveLength(recordCount);
      expect(original.components).toHaveLength(recordCount);
      const seedIds = [...new Set(seedSession.renderRequest.records.map(r => r.source.resourceId))];
      expect(original.components.map(({ size, sha256 }) => ({ size, sha256 })))
        .toEqual(seedIds.map(id => resourceIdentity(seedSession.resources[id])));
      if (journey === 'minimal') {
        const control = await save(page, testInfo, 'before-generate');
        expect(control.session.results).toHaveLength(1);
      }
      await generate(page);
      if (journey === 'grid-batch-grid') {
        await page.getByRole('checkbox', { name: 'Multi-Record Canvas', exact: true }).uncheck();
        await generate(page);
        expect((await snapshot(page)).results).toHaveLength(recordCount);
        await page.getByRole('checkbox', { name: 'Multi-Record Canvas', exact: true }).check();
        await generate(page);
      }
      const generated = await snapshot(page);
      expect(generated.records).toHaveLength(recordCount);
      expect(generated.results).toHaveLength(1);
      expect(generated.components).toEqual(original.components);
      const { saved, session, metrics } = await save(page, testInfo, 'generated');
      expect([session.version, session.webFiles.bindings.schema, session.renderRequest.schema]).toEqual([42, 2, 7]);
      const composite = session.webFiles.bindings.c_gb;
      expect(composite.kind).toBe('composite');
      expect(composite.components).toHaveLength(recordCount);
      expect(composite.components.map(part => resourceIdentity(session.resources[part.resourceId])))
        .toEqual(original.components.map(({ size, sha256 }) => ({ size, sha256 })));
      expect(composite.components.map(({ resourceId, ...metadata }) => metadata))
        .toEqual(original.components.map(({ resourceId, size, sha256, ...metadata }) => metadata));
      expect(metrics.filter(metric => ['base64EncodeCount', 'base64DecodeCount'].includes(metric.name)))
        .toEqual([]);
      const genomicSizes = Object.values(session.resources).filter(r => r.kind === 'genbank' || composite.components.some(c => session.resources[c.resourceId] === r)).map(r => r.size);
      expect(genomicSizes.reduce((sum, n) => sum + n, 0))
        .toBe(2 * original.components.reduce((sum, component) => sum + component.size, 0));
      const resources = Object.fromEntries(Object.entries(session.resources)
        .map(([id, descriptor]) => [id, resourceIdentity(descriptor)]));
      for (const [id, identity] of Object.entries(generated.resources)) expect(resources[id]).toEqual(identity);
      for (const component of original.components) {
        expect(Object.values(resources)).toContainEqual({ size: component.size, sha256: component.sha256 });
      }
      expect(session.renderRequest).toEqual(generated.request);
      await fs.writeFile(testInfo.outputPath('resource-inspection.json'), JSON.stringify({
        original: { ...original, results: original.results.map(hash) },
        generated: { ...generated, results: generated.results.map(hash) },
        resources, webFiles: session.webFiles
      }, null, 2));
      expect(external).toEqual([]);
      await page.context().close();

      const restored = await load(browser, saved, contexts);
      const loaded = await snapshot(restored.page);
      expect(loaded.source).toEqual(original.source);
      expect(loaded.components.map(({ resourceId, ...identity }) => identity))
        .toEqual(original.components.map(({ resourceId, ...identity }) => identity));
      expect(loaded.request).toEqual(generated.request);
      expect(loaded.resultNames).toEqual(generated.resultNames);
      expect(loaded.source.isArray).toBe(false);
      expect(loaded.results).toEqual(session.results.map(result => result.content));
      const resaved = await save(restored.page, testInfo, 'restored');
      expect(resaved.session.webFiles.bindings.c_gb).toEqual(composite);
      await generate(restored.page);
      const regenerated = await snapshot(restored.page);
      expect(regenerated.records).toEqual(generated.records);
      expect(regenerated.results).toEqual(generated.results);
      expect(regenerated.components).toEqual(loaded.components);
      await restored.page.screenshot({ path: testInfo.outputPath('regenerated.png') });
      if (journey === 'minimal') {
        await restored.page.getByLabel('GenBank/DDBJ File', { exact: true }).setInputFiles({
          name: 'replacement.gb', mimeType: 'text/plain',
          buffer: Buffer.from(seedSession.resources[seedIds[0]].data, 'base64')
        });
        const replaced = await save(restored.page, testInfo, 'replaced');
        expect(replaced.session.webFiles.bindings.c_gb.kind).toBeUndefined();
        expect(replaced.session.webFiles.bindings.c_gb.components).toBeUndefined();
        expect(resourceIdentity(replaced.session.resources[replaced.session.webFiles.bindings.c_gb.resourceId]))
          .toEqual(resourceIdentity(seedSession.resources[seedIds[0]]));
        const mobile = await load(browser, saved, contexts, { width: 390, height: 844 });
        const mobileState = await snapshot(mobile.page);
        expect(mobileState.components).toEqual(loaded.components);
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

test('Session CLI sidecar preserves the six-source draft for fresh Web Load and Generate', async ({ browser }, testInfo) => {
  test.setTimeout(4 * generateTimeout);
  const contexts = [];
  try {
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
    expect([replayed.version, replayed.webFiles.bindings.schema, replayed.renderRequest.schema]).toEqual([42, 2, 7]);
    const expected = session.webFiles.bindings.c_gb;
    const actual = replayed.webFiles.bindings.c_gb;
    expect(actual.components).toHaveLength(6);
    expect(actual.components.map(({ resourceId, ...metadata }) => metadata))
      .toEqual(expected.components.map(({ resourceId, ...metadata }) => metadata));
    expect(actual.components.map(part => resourceIdentity(replayed.resources[part.resourceId])))
      .toEqual(expected.components.map(part => resourceIdentity(session.resources[part.resourceId])));
    expect({ ...actual, components: [] }).toEqual({ ...expected, components: [] });
    expect(replayed.results[0].content).toBe(await fs.readFile(`${prefix}.svg`, 'utf8'));
    expect(replayed.results[0].name).toBe(path.basename(prefix));
    expect(replayed.editorState.featureCatalog.items[0].resultName).toBe(path.basename(prefix));
    const restored = await load(browser, sidecar, contexts);
    const loaded = await snapshot(restored.page);
    expect(loaded.source).toEqual(original.source);
    expect(loaded.source.isArray).toBe(false);
    expect(loaded.components.map(({ resourceId, ...identity }) => identity))
      .toEqual(original.components.map(({ resourceId, ...identity }) => identity));
    const admittedContents = await restored.page.evaluate(async results => {
      const { sanitizeSvgContent } = await import('./js/services/svg-sanitization.js');
      return results.map(result => sanitizeSvgContent(result.content, window.DOMPurify));
    }, replayed.results);
    expect(loaded.results.map(hash)).toEqual(admittedContents.map(hash));
    expect(loaded.request).toEqual(replayed.renderRequest);
    await generate(restored.page);
    expect((await snapshot(restored.page)).records).toEqual(original.records);
    expect(restored.external).toEqual([]);
    await restored.page.screenshot({ path: testInfo.outputPath('cli-regenerated.png') });
  } finally {
    for (const context of contexts) await context.close();
  }
});

test('Session import rejects malformed composite bindings without replacing existing work @pr-smoke', async ({ browser }, testInfo) => {
  test.setTimeout(6 * loadTimeout);
  const contexts = [];
  try {
    const seedPath = 'gbdraw/web/gallery/sessions/HmmtDNA_basic_circular.gbdraw-session.json';
    const { page, external } = await load(browser, seedPath, contexts);
    const before = await snapshot(page);
    const seedSession = JSON.parse(await fs.readFile(seedPath));
    const resourceId = seedSession.renderRequest.records[0].source.resourceId;
    const leaf = { resourceId, name: 'part.gb', type: 'text/plain', lastModified: 0 };
    const dialogs = [];
    page.on('dialog', dialog => { dialogs.push(dialog.message()); });
    const mutations = {
      schema: s => { s.webFiles.bindings.schema = 99; },
      kind: s => { s.webFiles.bindings.c_gb.kind = 'unknown'; },
      nested: s => { s.webFiles.bindings.c_gb.components[0] = structuredClone(s.webFiles.bindings.c_gb); },
      empty: s => { s.webFiles.bindings.c_gb.components = []; },
      singleton: s => { s.webFiles.bindings.c_gb.components.length = 1; },
      missing: s => { delete s.webFiles.bindings.c_gb; },
      dangling: s => { s.webFiles.bindings.c_gb.components[0].resourceId = 'missing'; },
      null: s => { s.webFiles.bindings.c_gb.components[0] = null; },
      metadata: s => { s.webFiles.bindings.c_gb.lastModified = -1; },
      mixed: s => { s.webFiles.bindings.c_gb.resourceId = resourceId; },
      slot: s => { s.webFiles.bindings.c_fasta = s.webFiles.bindings.c_gb; },
      payload: s => { s.resources[resourceId].encoding = 'unsupported'; }
    };
    const rejected = [];
    for (const [label, mutate] of Object.entries(mutations)) {
      const invalid = structuredClone(seedSession);
      invalid.webFiles ||= {};
      invalid.webFiles.bindings = { schema: 2, c_gb: { kind: 'composite',
        components: [structuredClone(leaf), structuredClone(leaf)],
        name: 'logical.gb', type: 'text/plain', lastModified: 0 } };
      mutate(invalid);
      const previous = dialogs.length;
      await page.locator('input[accept^=".json,"]').setInputFiles({
        name: `${label}.gbdraw-session.json`, mimeType: 'application/json', buffer: Buffer.from(JSON.stringify(invalid))
      });
      await expect.poll(() => dialogs.length, { timeout: loadTimeout }).toBe(previous + 1);
      await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.sessionImportPending), { timeout: loadTimeout }).toBe(false);
      expect(dialogs.at(-1)).toMatch(/^Failed to load session:/);
      expect(await snapshot(page)).toEqual(before);
      rejected.push({ label, error: dialogs.at(-1) });
    }
    expect(external).toEqual([]);
    await fs.writeFile(testInfo.outputPath('transactional-rejections.json'), JSON.stringify(rejected, null, 2));
  } finally {
    for (const context of contexts) await context.close();
  }
});
