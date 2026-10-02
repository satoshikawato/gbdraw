// Shared helpers for the manual audit sweeps in tools/audit/.
// Evidence goes to GBDRAW_AUDIT_OUT (default: <os tmpdir>/gbdraw-audit), never into the repository.
const { execFileSync } = require('node:child_process');
const { mkdirSync, readFileSync, readdirSync, writeFileSync } = require('node:fs');
const { tmpdir } = require('node:os');
const { join, resolve } = require('node:path');
const { gunzipSync } = require('node:zlib');

const REPO_ROOT = resolve(__dirname, '..', '..', '..');
const AUDIT_OUT = resolve(process.env.GBDRAW_AUDIT_OUT || join(tmpdir(), 'gbdraw-audit'));

const helpers = () => require(join(REPO_ROOT, 'tests/web/helpers/app-lifecycle.cjs'));

const outDir = (...parts) => {
  const dir = join(AUDIT_OUT, ...parts);
  mkdirSync(dir, { recursive: true });
  return dir;
};

const writeEvidence = (dir, name, data) => {
  writeFileSync(join(dir, name), typeof data === 'string' ? data : JSON.stringify(data, null, 2));
};

// Page errors, console errors, and warnings seen during a sweep.
const collect = (page) => {
  const log = { pageErrors: [], consoleErrors: [], console: [] };
  page.on('pageerror', (e) => log.pageErrors.push(String(e?.stack || e?.message || e)));
  page.on('console', (m) => {
    if (m.type() === 'error') log.consoleErrors.push(m.text());
    if (m.type() === 'warning' || m.type() === 'error') log.console.push(`${m.type()}: ${m.text()}`);
  });
  return log;
};

const setFile = (page, key, path, name) => page.evaluate(async ({ key, content, name }) => {
  const app = window.__GBDRAW_APP__;
  app.files[key] = new File([content], name, { type: 'text/plain' });
  await window.Vue.nextTick();
}, { key, content: readFileSync(path, 'utf8'), name: name || path.split('/').pop() });

// Generate and return the Result plus the Source recipe the GUI offers for it.
const runGenerate = (page) => helpers().evaluateWithRetainedPromise(page, async () => {
  const app = window.__GBDRAW_APP__;
  const started = performance.now();
  let result;
  try {
    result = await app.runAnalysis();
  } catch (error) {
    result = { status: 'threw', message: String(error?.message || error) };
  }
  const recipe = app.lastRunInfo?.sourceRecipe || null;
  return {
    status: result?.status,
    elapsedMs: Math.round(performance.now() - started),
    errorSummary: String(app.errorLog?.summary || ''),
    errorDetails: Array.isArray(app.errorLog?.details)
      ? app.errorLog.details.map((d) => (d && typeof d === 'object') ? `${d.label}: ${d.text}` : String(d))
      : [],
    svg: String(app.results?.[0]?.content || ''),
    command: recipe?.command || app.lastRunInfo?.command || '',
    sourceRecipe: recipe
      ? { level: recipe.level, reason: recipe.reason || recipe.unavailableReason || '' }
      : null,
    exactReplay: app.lastRunInfo?.exactReplay?.command || ''
  };
});

// Download the Run info helper files (tables the Source recipe refers to) into dir.
const captureHelpers = async (page, dir) => {
  const has = await page.evaluate(() => {
    const files = window.__GBDRAW_APP__.lastRunInfo?.helperFiles;
    return Array.isArray(files) && files.length > 0;
  });
  if (!has) return [];
  mkdirSync(dir, { recursive: true });
  const [download] = await Promise.all([
    page.waitForEvent('download', { timeout: 30_000 }),
    page.evaluate(() => window.__GBDRAW_APP__.downloadCliHelperFiles())
  ]);
  const zipPath = join(dir, 'helpers.zip');
  await download.saveAs(zipPath);
  execFileSync('python3', ['-c', 'import zipfile,sys; zipfile.ZipFile(sys.argv[1]).extractall(sys.argv[2])', zipPath, dir]);
  return readdirSync(dir);
};

const readSession = (file) => {
  const buffer = readFileSync(file);
  const gzip = buffer[0] === 0x1f && buffer[1] === 0x8b;
  return JSON.parse((gzip ? gunzipSync(buffer) : buffer).toString('utf8'));
};

const SESSION_INPUT = 'input[accept^=".json,"]';

const loadSession = async (page, file) => {
  const dialogs = [];
  const accept = async (dialog) => { dialogs.push(dialog.message()); await dialog.accept(); };
  page.on('dialog', accept);
  const started = Date.now();
  try {
    await page.locator(SESSION_INPUT).setInputFiles(file);
    await page.waitForFunction(() => !window.__GBDRAW_APP__.sessionImportPending, null, { timeout: 300_000 });
    const error = await page.evaluate(() => {
      const log = window.__GBDRAW_APP__.errorLog;
      return log ? JSON.parse(JSON.stringify(log)) : null;
    });
    return { dialogs, error, ms: Date.now() - started };
  } finally {
    page.off('dialog', accept);
  }
};

const saveSession = async (page, file) => {
  const dialogs = [];
  const accept = async (dialog) => { dialogs.push(dialog.message()); await dialog.accept(); };
  page.on('dialog', accept);
  try {
    const pending = page.waitForEvent('download', { timeout: 300_000 });
    await page.getByRole('banner').getByRole('button', { name: 'Save Session', exact: true }).click();
    const download = await pending;
    await download.saveAs(file);
    await page.waitForFunction(() => !window.__GBDRAW_APP__.sessionSavePending, null, { timeout: 120_000 });
    return { doc: readSession(file), dialogs, suggested: download.suggestedFilename() };
  } finally {
    page.off('dialog', accept);
  }
};

// Elements that extend past the viewport, for the viewport-width sweeps.
const measureOverflow = (page) => page.evaluate(() => {
  const doc = document.documentElement;
  const vw = window.innerWidth;
  const offenders = [];
  for (const el of document.querySelectorAll('body *')) {
    const r = el.getBoundingClientRect();
    if (r.width > 0 && r.right > vw + 1 && getComputedStyle(el).position !== 'fixed') {
      offenders.push({
        tag: el.tagName,
        id: el.id,
        cls: String(el.className).slice(0, 80),
        right: Math.round(r.right),
        text: (el.textContent || '').trim().slice(0, 40)
      });
      if (offenders.length > 15) break;
    }
  }
  return {
    vw,
    scrollWidth: doc.scrollWidth,
    clientWidth: doc.clientWidth,
    overflowX: doc.scrollWidth > doc.clientWidth + 1,
    offenders
  };
});

module.exports = {
  AUDIT_OUT,
  REPO_ROOT,
  captureHelpers,
  collect,
  helpers,
  loadSession,
  measureOverflow,
  outDir,
  readSession,
  runGenerate,
  saveSession,
  setFile,
  writeEvidence
};
