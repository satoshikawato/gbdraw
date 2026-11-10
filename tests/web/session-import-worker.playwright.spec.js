const { test, expect } = require('@playwright/test');
const { readFileSync } = require('node:fs');
const { gunzipSync, gzipSync } = require('node:zlib');
const path = require('node:path');
const { pathToFileURL } = require('node:url');
const { openApp } = require('./helpers/app-lifecycle.cjs');
const { installImportReadGate } = require('./helpers/session-import-gate.cjs');

const inputSelector = 'input[type="file"][accept*=".gbdraw"]';
// A current Session file container, and a JSON Session an older gbdraw saved.
let baseline;
let container;
const legacyGzip = readFileSync('tests/fixtures/sessions/HmmtDNA-wrap-5kb.v50-web0145.gbdraw-session.json.gz');
test.beforeAll(async () => {
  const { gallerySessionFile } = await import(pathToFileURL(path.resolve('tests/web/helpers/written-resources.mjs')).href);
  baseline = Buffer.from(await (await gallerySessionFile('HmmtDNA_basic_circular')).arrayBuffer());
  container = await import(pathToFileURL(path.resolve('gbdraw/web/js/services/session-container.js')).href);
});
const installProbe = page => page.addInitScript(() => {
  const NativeWorker = Worker;
  const probe = window.__importWorkers = { entries: [], events: [], urlsCreated: 0, urlsRevoked: 0 };
  window.__GBDRAW_TEST_HOOKS__ = { onSessionLifecycleEvent: event => probe.events.push(event) };
  window.Worker = class extends NativeWorker {
    constructor(...args) {
      super(...args);
      this.entry = { url: String(args[0]), listeners: 0, terminated: 0 };
      probe.entries.push(this.entry);
    }
    addEventListener(...args) { this.entry.listeners += 1; return super.addEventListener(...args); }
    removeEventListener(...args) { this.entry.listeners -= 1; return super.removeEventListener(...args); }
    terminate() { this.entry.terminated += 1; return super.terminate(); }
  };
  const create = URL.createObjectURL, revoke = URL.revokeObjectURL;
  URL.createObjectURL = (...args) => { probe.urlsCreated += 1; return create(...args); };
  URL.revokeObjectURL = (...args) => { probe.urlsRevoked += 1; return revoke(...args); };
});
// Load opens a Session by its content, not by its name or type.
const load = async (page, buffer, name = 'session.gbdraw') => {
  await page.locator(inputSelector).setInputFiles({ name, mimeType: 'application/octet-stream', buffer });
  await page.waitForFunction(() => !window.__GBDRAW_APP__.sessionImportPending);
};
// The baseline container with its document's version changed.
const containerOfVersion = async (version) => {
  const blob = new Blob([baseline]);
  const read = await container.readSessionContainer(blob);
  const document = { ...read.document, version };
  const written = await container.writeSessionContainer({
    document,
    documentBytes: () => new Blob([JSON.stringify(document)]).stream(),
    readResource: (id) => container.openSessionContainerEntry(blob, read.resources[id])
  });
  return Buffer.from(await written.arrayBuffer());
};
// The files this origin holds in browser storage (the input copies).
const storedFiles = page => page.evaluate(async () => {
  const count = async (folder) => {
    let files = 0;
    for await (const handle of folder.values()) files += handle.kind === 'file' ? 1 : await count(handle);
    return files;
  };
  return count(await navigator.storage.getDirectory());
});
const retained = page => page.evaluate(async () => {
  const { state } = await import('/gbdraw/web/js/state.js');
  const service = await import('/gbdraw/web/js/services/config.js');
  const history = window.__GBDRAW_HISTORY__;
  window.__beforeImport = { results: state.results.value, canonical: service.getCommittedCanonicalSession(),
    files: state.activeDrawing().circularSources.map(card => card.gb), cache: state.losatCache.value, undo: history.getUndoCount(), redo: history.getRedoCount() };
});
const expectRetained = async page => expect(await page.evaluate(async () => {
  const { state } = await import('/gbdraw/web/js/state.js');
  const service = await import('/gbdraw/web/js/services/config.js');
  const history = window.__GBDRAW_HISTORY__, before = window.__beforeImport;
  return before.results === state.results.value && before.canonical === service.getCommittedCanonicalSession()
    && before.files.length === state.activeDrawing().circularSources.length
    && before.files.every((file, index) => file === state.activeDrawing().circularSources[index].gb) && before.cache === state.losatCache.value
    && before.undo === history.getUndoCount() && before.redo === history.getRedoCount();
})).toBe(true);

test('real import Worker container, JSON and gzip loads terminate before preflight with no transport resources retained', async ({ page }) => {
  await installProbe(page);
  page.on('dialog', dialog => dialog.accept());
  await openApp(page);
  const files = [[baseline, 'gbdraw'], [baseline, 'gbdraw'], [legacyGzip, 'json.gz'], [gunzipSync(legacyGzip), 'json'], [baseline, 'gbdraw']];
  for (const [i, [buffer, extension]] of files.entries()) {
    await load(page, buffer, `repeat-${i}.${extension}`);
    expect(await page.evaluate(() => window.__GBDRAW_APP__.errorLog), `repeat-${i}`).toBeNull();
    // A current Session loads with no diagram Worker; an older JSON Session's
    // conversion may start one.
    if (i === 1) {
      expect((await page.evaluate(() => window.__importWorkers)).entries
        .filter(e => e.url.includes('diagram-generation-worker'))).toHaveLength(0);
    }
  }
  const probe = await page.evaluate(() => window.__importWorkers);
  const imports = probe.entries.filter(e => e.url.includes('session-import-worker'));
  expect(imports).toHaveLength(5);
  for (const entry of imports) expect([entry.terminated, entry.listeners]).toEqual([1, 0]);
  expect(probe.urlsCreated).toBe(probe.urlsRevoked);
  const names = probe.events.map(e => e.name);
  let start = 0;
  for (let i = 0; i < 5; i += 1) {
    const end = names.indexOf('session-import-worker-terminated', start);
    const preflight = names.indexOf('current-session-preflight-start', start);
    expect(end).toBeGreaterThanOrEqual(start);
    expect(preflight).toBeGreaterThan(end);
    start = preflight + 1;
  }
});

test('unsafe keys, malformed JSON/gzip, fatal UTF-8, refused containers, storage refusal and Worker crash preserve the old document and History', async ({ page, context }) => {
  await installProbe(page);
  const dialogs = [];
  page.on('dialog', dialog => { dialogs.push(dialog.message()); dialog.accept(); });
  await openApp(page);
  await load(page, baseline);
  await expect(page.locator('[data-app-notice]')).toContainText('Session loaded: session.gbdraw', { timeout: 120_000 });
  await expect(page.locator('[data-app-notice] li')).toHaveCount(0);
  await retained(page);
  for (const [index, buffer] of [Buffer.from('{"nested":{"__proto__":{}}}'), Buffer.from('{private broken'),
    Buffer.from([0xff]), Buffer.from([0x1f, 0x8b, 1, 2, 3]), gzipSync(Buffer.from([0xc3]))].entries()) {
    await load(page, buffer);
    const alert = page.getByRole('alert', { name: 'Session not loaded' });
    await expect(alert).toBeVisible();
    const error = await page.evaluate(() => window.__GBDRAW_APP__.errorLog);
    // X-01 (SE-09): the Worker's final diagnostic reaches the panel unchanged.
    expect(error).toMatchObject(index === 0
      // OV-952: an unsafe key is the Session-schema diagnostic.
      ? { code: 'INPUT_INVALID', operation: 'session-load', stage: 'request-validation', context: { field: 'schema', reason: 'FIELDS' } }
      : index === 1
        ? { code: 'INPUT_INVALID', operation: 'session-load', stage: 'parse', context: { field: 'schema', reason: 'JSON_FORMAT' } }
        // Bytes that are not UTF-8 or not valid gzip fail while the Worker reads them.
        : { code: 'INPUT_UNREADABLE', operation: 'session-load', stage: 'read', context: { field: 'schema' } });
    expect(JSON.stringify(error)).not.toMatch(/private broken|__proto__|Traceback/);
    expect(dialogs).toEqual([]);
    await expect(page.locator('[data-app-notice]', { hasText: 'Session loaded' })).toHaveCount(0);
    await expectRetained(page);
  }
  // A refused container, a Session 51 JSON file, a key written twice (OV-2694)
  // and a container whose inputs do not fit in browser storage each leave the
  // loaded Session, its Result and History as they were (PD-OI-016).
  const read = await container.readSessionContainer(new Blob([baseline]));
  const entry = Object.values(read.resources)[0];
  const damaged = Buffer.from(baseline);
  damaged[entry.end - 1] ^= 0xff;
  const refusals = [
    [damaged, { code: 'INPUT_INVALID', stage: 'read', context: { field: 'schema', reason: 'CONTAINER_DAMAGED' } }],
    [Buffer.from(JSON.stringify(read.document)),
      { code: 'INPUT_INVALID', stage: 'parse', context: { field: 'schema', reason: 'SESSION_VERSION_UNSUPPORTED' } }],
    [Buffer.from('{"version":50,"version":50}'),
      { code: 'INPUT_INVALID', stage: 'parse', context: { field: 'schema', reason: 'JSON_FORMAT' } }],
    [Buffer.from(gunzipSync(legacyGzip).toString('utf8').replace('{', '{"version":50,')),
      { code: 'INPUT_INVALID', stage: 'parse', context: { field: 'schema', reason: 'JSON_FORMAT' } }],
    // A JSON Session of a newer version says a newer gbdraw saved it; a
    // container holds a Session 51 only, and one of a newer gbdraw is refused
    // before any of its inputs is copied.
    [Buffer.from(JSON.stringify({ ...read.document, version: 52 })),
      { code: 'INPUT_INVALID', context: { field: 'schema', reason: 'SESSION_TOO_NEW', sessionVersion: 52 } }],
    [await containerOfVersion(50),
      { code: 'INPUT_INVALID', stage: 'parse', context: { field: 'schema', reason: 'SESSION_VERSION_UNSUPPORTED' } }],
    [await containerOfVersion(52),
      { code: 'INPUT_INVALID', context: { field: 'schema', reason: 'SESSION_TOO_NEW', sessionVersion: 52 } }],
  ];
  for (const [buffer, expected] of refusals) {
    const copiesBeforeRefusal = await storedFiles(page);
    await load(page, buffer);
    await expect(page.getByRole('alert', { name: 'Session not loaded' })).toBeVisible();
    expect(await page.evaluate(() => window.__GBDRAW_APP__.errorLog)).toMatchObject({ operation: 'session-load', ...expected });
    await expect(page.locator('[data-app-notice]', { hasText: 'Session loaded' })).toHaveCount(0);
    if (expected.context.reason.startsWith('SESSION_')) expect(await storedFiles(page)).toBe(copiesBeforeRefusal);
    await expectRetained(page);
  }
  // The storage check runs before any input is written.
  await page.evaluate(() => {
    window.__GBDRAW_ESTIMATES__ = 0;
    navigator.storage.estimate = async () => { window.__GBDRAW_ESTIMATES__ += 1; return { quota: 1_000_000, usage: 999_000 }; };
  });
  const copiesBefore = await storedFiles(page);
  await load(page, baseline);
  await expect(page.getByRole('alert', { name: 'Session not loaded' })).toContainText(
    `The Session was not loaded: its inputs need ${(read.inputBytes / 1048576).toFixed(1)} MB of browser storage and 0.0 MB is free.`);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.errorLog)).toMatchObject({
    code: 'SESSION_SIZE_LIMIT', operation: 'session-load', stage: 'read',
    context: { storageNeedBytes: read.inputBytes, storageFreeBytes: 1000 }
  });
  expect(await storedFiles(page)).toBe(copiesBefore);
  expect(await page.evaluate(() => window.__GBDRAW_ESTIMATES__)).toBe(1);
  await expectRetained(page);
  await context.route('**/workers/session-import-worker.js', route => route.fulfill({
    contentType: 'text/javascript', body: 'throw new Error("controlled Worker crash")'
  }));
  await load(page, baseline);
  await expect(page.getByRole('alert', { name: 'Session not loaded' })).toBeVisible();
  // A crashed import Worker leaves no import path: a transport failure.
  expect(await page.evaluate(() => window.__GBDRAW_APP__.errorLog)).toMatchObject({
    code: 'SESSION_IMPORT_UNAVAILABLE', operation: 'session-load', stage: 'transport'
  });
  expect(dialogs).toEqual([]);
  await expect(page.locator('[data-app-notice]', { hasText: 'Session loaded' })).toHaveCount(0);
  await expectRetained(page);
  // Every import Worker ended; the input copy Worker lives on with the tab.
  const entries = await page.evaluate(() => window.__importWorkers.entries);
  expect(entries.filter(e => !e.url.includes('session-import-worker')).map(e => e.url.replace(/^.*\//, '')))
    .toEqual(['input-copy-worker.js']);
  for (const entry of entries.filter(e => e.url.includes('session-import-worker'))) {
    expect([entry.terminated, entry.listeners]).toEqual([1, 0]);
  }
});

test('teardown cancels the actual pending Worker before adoption and permits a new load', async ({ page }) => {
  await installProbe(page);
  page.on('dialog', dialog => dialog.accept());
  await openApp(page);
  await load(page, baseline);
  await retained(page);
  await installImportReadGate(page, 'pending.gbdraw');
  await page.locator(inputSelector).setInputFiles({ name: 'pending.gbdraw', mimeType: 'application/octet-stream', buffer: baseline });
  await page.waitForFunction(() => window.__GBDRAW_SESSION_IMPORT_GATE__.streamInvocations === 1);
  const outcome = await page.evaluate(async () => {
    const service = await import('/gbdraw/web/js/services/config.js');
    service.disposeSessionOperations();
    return window.__GBDRAW_APP__.sessionImportPending;
  });
  expect(outcome).toBe(false);
  await expectRetained(page);
  await page.evaluate(() => window.__GBDRAW_SESSION_IMPORT_GATE__.release());
  await load(page, baseline, 'after-teardown.gbdraw');
  // The input copy Worker a container Load starts lives on with the tab.
  const entries = (await page.evaluate(() => window.__importWorkers.entries)).filter(e => e.url.includes('session-import-worker'));
  expect(entries).toHaveLength(3);
  for (const entry of entries) expect([entry.terminated, entry.listeners]).toEqual([1, 0]);
});
