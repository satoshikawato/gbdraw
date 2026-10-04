const { test, expect } = require('@playwright/test');
const { readFileSync } = require('node:fs');
const { gzipSync } = require('node:zlib');
const { openApp } = require('./helpers/app-lifecycle.cjs');
const { installImportReadGate } = require('./helpers/session-import-gate.cjs');

const inputSelector = 'input[type="file"][accept*="application/json"][accept*="application/gzip"]';
const baseline = readFileSync('gbdraw/web/gallery/sessions/HmmtDNA_basic_circular.gbdraw-session.json');
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
const load = async (page, buffer, name = 'session.json') => {
  await page.locator(inputSelector).setInputFiles({ name, mimeType: 'application/json', buffer });
  await page.waitForFunction(() => !window.__GBDRAW_APP__.sessionImportPending);
};
const retained = page => page.evaluate(async () => {
  const { state } = await import('/gbdraw/web/js/state.js');
  const service = await import('/gbdraw/web/js/services/config.js');
  const history = window.__GBDRAW_HISTORY__;
  window.__beforeImport = { results: state.results.value, canonical: service.getCommittedCanonicalSession(),
    file: state.files.c_gb, cache: state.losatCache.value, undo: history.getUndoCount(), redo: history.getRedoCount() };
});
const expectRetained = async page => expect(await page.evaluate(async () => {
  const { state } = await import('/gbdraw/web/js/state.js');
  const service = await import('/gbdraw/web/js/services/config.js');
  const history = window.__GBDRAW_HISTORY__, before = window.__beforeImport;
  return before.results === state.results.value && before.canonical === service.getCommittedCanonicalSession()
    && before.file === state.files.c_gb && before.cache === state.losatCache.value
    && before.undo === history.getUndoCount() && before.redo === history.getRedoCount();
})).toBe(true);

test('real import Worker JSON/gzip and repeated loads terminate before preflight with no transport resources retained', async ({ page }) => {
  await installProbe(page);
  page.on('dialog', dialog => dialog.accept());
  await openApp(page);
  for (let i = 0; i < 5; i += 1) await load(page, i % 2 ? gzipSync(baseline) : baseline, `repeat-${i}.json`);
  const probe = await page.evaluate(() => window.__importWorkers);
  expect(probe.entries.filter(e => e.url.includes('diagram-generation-worker'))).toHaveLength(0);
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

test('unsafe keys, malformed JSON/gzip, fatal UTF-8 and Worker crash preserve the old document and History', async ({ page, context }) => {
  await installProbe(page);
  const dialogs = [];
  page.on('dialog', dialog => { dialogs.push(dialog.message()); dialog.accept(); });
  await openApp(page);
  await load(page, baseline);
  await retained(page);
  for (const [index, buffer] of [Buffer.from('{"nested":{"__proto__":{}}}'), Buffer.from('{private broken'),
    Buffer.from([0xff]), Buffer.from([0x1f, 0x8b, 1, 2, 3]), gzipSync(Buffer.from([0xc3]))].entries()) {
    await load(page, buffer);
    const alert = page.getByRole('alert', { name: 'Operation error' });
    await expect(alert).toBeVisible();
    const error = await page.evaluate(() => window.__GBDRAW_APP__.errorLog);
    // X-01 (SE-09): the Worker's final diagnostic reaches the panel unchanged.
    expect(error).toMatchObject(index === 0
      ? { code: 'UNKNOWN', operation: 'unknown', stage: 'request-validation' }
      : index === 1
        ? { code: 'INPUT_INVALID', operation: 'unknown', stage: 'parse', context: { field: 'schema', reason: 'JSON_FORMAT' } }
        // Bytes that are not UTF-8 or not valid gzip fail while the Worker reads them.
        : { code: 'INPUT_UNREADABLE', operation: 'unknown', stage: 'read', context: { field: 'schema' } });
    expect(JSON.stringify(error)).not.toMatch(/private broken|__proto__|Traceback/);
    expect(dialogs).toEqual(['Session loaded successfully!']);
    await expectRetained(page);
  }
  await context.route('**/workers/session-import-worker.js', route => route.fulfill({
    contentType: 'text/javascript', body: 'throw new Error("controlled Worker crash")'
  }));
  await load(page, baseline);
  await expect(page.getByRole('alert', { name: 'Operation error' })).toBeVisible();
  // A crashed import Worker leaves no import path: a transport failure.
  expect(await page.evaluate(() => window.__GBDRAW_APP__.errorLog)).toMatchObject({
    code: 'SESSION_IMPORT_UNAVAILABLE', operation: 'unknown', stage: 'transport'
  });
  expect(dialogs).toEqual(['Session loaded successfully!']);
  await expectRetained(page);
  const entries = await page.evaluate(() => window.__importWorkers.entries);
  for (const entry of entries) expect([entry.terminated, entry.listeners]).toEqual([1, 0]);
});

test('teardown cancels the actual pending Worker before adoption and permits a new load', async ({ page }) => {
  await installProbe(page);
  page.on('dialog', dialog => dialog.accept());
  await openApp(page);
  await load(page, baseline);
  await retained(page);
  await installImportReadGate(page, 'pending.json');
  await page.locator(inputSelector).setInputFiles({ name: 'pending.json', mimeType: 'application/json', buffer: baseline });
  await page.waitForFunction(() => window.__GBDRAW_SESSION_IMPORT_GATE__.streamInvocations === 1);
  const outcome = await page.evaluate(async () => {
    const service = await import('/gbdraw/web/js/services/config.js');
    service.disposeSessionOperations();
    return window.__GBDRAW_APP__.sessionImportPending;
  });
  expect(outcome).toBe(false);
  await expectRetained(page);
  await page.evaluate(() => window.__GBDRAW_SESSION_IMPORT_GATE__.release());
  await load(page, baseline, 'after-teardown.json');
  const entries = await page.evaluate(() => window.__importWorkers.entries);
  expect(entries).toHaveLength(3);
  for (const entry of entries) expect([entry.terminated, entry.listeners]).toEqual([1, 0]);
});
