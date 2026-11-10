// The temporary copies of the chosen input files (PD-OI-121 in
// docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md): the page keeps the
// content of a chosen file in a copy in the browser's private site storage, so
// a file changed, touched, moved, or deleted on disk after it was chosen still
// Generates and Saves as chosen. A copy is deleted when no state uses its
// File, and the folders of closed tabs are deleted when the page next opens.
const { test, expect } = require('@playwright/test');
const fs = require('node:fs');
const path = require('node:path');
const { gzipSync } = require('node:zlib');
const { evaluateWithRetainedPromise, routeCrossOriginIsolated } = require('./helpers/app-lifecycle.cjs');
const { loadSessionFile, openFresh } = require('./helpers/audit-browser.cjs');
const { generate } = require('./helpers/mode-transition.cjs');
const { circularFileInput } = require('./helpers/circular-inputs.cjs');
const { showDrawing } = require('./helpers/drawing-bar.cjs');
const { readSessionFileEntries } = require('./helpers/session-file.cjs');

const SOURCE = path.resolve('tests/test_inputs/NC_001416.gb');
const OTHER = path.resolve('tests/test_inputs/AP027078.gb');

// Every copy file of every tab folder: [{ tab, name, size }]. The app writes and removes
// copies while this lists them, so an entry that is gone when read is left out.
const copyFiles = (page) => page.evaluate(async () => {
  const base = await (await navigator.storage.getDirectory()).getDirectoryHandle('gbdraw-input-copies', { create: true });
  const files = [];
  const skipGone = (error) => {
    if (error.name !== 'NotFoundError') throw error;
  };
  for await (const [tab, folder] of base.entries()) {
    try {
      for await (const [name, handle] of folder.entries()) {
        try {
          files.push({ tab, name, size: (await handle.getFile()).size });
        } catch (error) {
          skipGone(error);
        }
      }
    } catch (error) {
      skipGone(error);
    }
  }
  return files;
});
const copySizes = async (page) => (await copyFiles(page)).map(({ size }) => size).sort((a, b) => a - b);

const placeCopy = (info, from, name) => {
  const to = info.outputPath(name);
  fs.mkdirSync(path.dirname(to), { recursive: true });
  fs.copyFileSync(from, to);
  return to;
};
const changes = {
  modified: (file) => fs.appendFileSync(file, '\n'),
  touched: (file) => { const later = new Date(Date.now() + 5000); fs.utimesSync(file, later, later); },
  moved: (file) => fs.renameSync(file, `${file}.moved`),
  deleted: (file) => fs.unlinkSync(file)
};

const chooseGenBank = async (page, file) => {
  await circularFileInput(page, 'gb').setInputFiles(file);
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.circularRecordDiscoveryState.status),
    { timeout: 60_000 }).toBe('ready');
};
// Saves the Session; returns the saved file and each GenBank resource as the file stores it:
// [{ size (its declared size), method (0 stored, 8 deflated), bytes }].
const saveSessionFile = async (page, title) => {
  const pending = page.waitForEvent('download', { timeout: 60_000 });
  await evaluateWithRetainedPromise(page, async (sessionTitle) => {
    const app = window.__GBDRAW_APP__;
    app.sessionTitle = sessionTitle;
    await app.saveSessionWithTitle();
  }, title);
  const file = await (await pending).path();
  const { document, resources } = await readSessionFileEntries(file);
  const genbank = Object.entries(document.resources).filter(([, descriptor]) => descriptor.kind === 'genbank')
    .map(([id, descriptor]) => ({ size: descriptor.size, ...resources[id] }));
  return { file, genbank };
};
const saveSession = async (page, title) => (await saveSessionFile(page, title)).genbank.map(({ bytes }) => bytes);

test('a file changed, touched, moved, or deleted after it was chosen Saves and Generates as chosen', async ({ page }, info) => {
  test.setTimeout(480_000);
  const chosenBytes = fs.readFileSync(SOURCE);
  await openFresh(page);
  for (const [kind, change] of Object.entries(changes)) {
    const chosen = placeCopy(info, SOURCE, `${kind}/chosen.gb`);
    await chooseGenBank(page, chosen);
    await expect.poll(() => copySizes(page), { timeout: 30_000 }).toContain(chosenBytes.byteLength);
    change(chosen);
    const saved = await saveSession(page, `copy-${kind}`);
    expect(saved, kind).toHaveLength(1);
    expect(saved[0].equals(chosenBytes), `${kind}: Save writes the content as chosen`).toBe(true);
    await generate(page);
    expect(await page.evaluate(() => window.__GBDRAW_APP__.results.length), kind).toBe(1);
  }
});

test('a replaced file\'s copy is deleted once History no longer holds the File', async ({ page }, info) => {
  test.setTimeout(240_000);
  await openFresh(page);
  const first = fs.statSync(SOURCE).size;
  const second = fs.statSync(OTHER).size;
  await chooseGenBank(page, placeCopy(info, SOURCE, 'first.gb'));
  await chooseGenBank(page, placeCopy(info, OTHER, 'second.gb'));
  // Undo can bring the first file back, so its copy stays.
  await expect.poll(() => copySizes(page), { timeout: 30_000 }).toEqual([first, second].sort((a, b) => a - b));
  await page.evaluate(() => window.__GBDRAW_HISTORY__.captureBaseline('No more Undo'));
  await expect.poll(() => copySizes(page), { timeout: 30_000 }).toEqual([second]);
});

test('the page deletes the folders of closed tabs when it opens, and keeps an open tab\'s', async ({ page }, info) => {
  test.setTimeout(240_000);
  await openFresh(page);
  await chooseGenBank(page, placeCopy(info, SOURCE, 'live.gb'));
  await expect.poll(() => copyFiles(page), { timeout: 30_000 }).toHaveLength(1);
  const [{ tab: liveTab }] = await copyFiles(page);
  // A tab that ended without deleting its folder (a crash): nobody holds its lock.
  await page.evaluate(async () => {
    const base = await (await navigator.storage.getDirectory()).getDirectoryHandle('gbdraw-input-copies');
    const folder = await base.getDirectoryHandle('tab-closed', { create: true });
    const writable = await (await folder.getFileHandle('copy-1.bin', { create: true })).createWritable();
    await writable.write(new Uint8Array(64));
    await writable.close();
  });
  const later = await page.context().newPage();
  await openFresh(later);
  await expect.poll(async () => (await copyFiles(later)).map(({ tab }) => tab)).toEqual([liveTab]);
  await later.close();
});

test('a file whose copy could not be written is read itself, and refuses once it changed', async ({ page }, info) => {
  test.setTimeout(240_000);
  await page.addInitScript(() => {
    window.__GBDRAW_TEST_HOOKS__ = { ...(window.__GBDRAW_TEST_HOOKS__ || {}), failInputCopyWrites: true };
  });
  await openFresh(page);
  const chosen = placeCopy(info, SOURCE, 'unwritten.gb');
  await chooseGenBank(page, chosen);
  await page.waitForTimeout(1000);
  expect(await copyFiles(page)).toEqual([]);
  changes.modified(chosen);
  await page.getByRole('button', { name: 'Generate Diagram', exact: true }).click();
  await expect.poll(() => page.evaluate(async () => (await import('./js/state.js')).state.errorLog.value?.code ?? null),
    { timeout: 180_000 }).toBe('INPUT_UNREADABLE');
  await expect(page.getByText('An input could not be read. Replace or reselect it.')).toBeVisible();
});

// Receipt session.file-container: a gzip input goes into the Session file as it is.
test('a gzip input\'s copy holds its gzip bytes, Save stores them as they are, and the loaded Session draws', async ({ page }, info) => {
  test.setTimeout(480_000);
  const plain = fs.readFileSync(SOURCE);
  const gzip = gzipSync(plain);
  const gzipped = info.outputPath('NC_001416.gb.gz');
  fs.writeFileSync(gzipped, gzip);
  await openFresh(page);
  await chooseGenBank(page, gzipped);
  await expect.poll(() => copySizes(page), { timeout: 30_000 }).toEqual([gzip.byteLength]);
  const copied = await page.evaluate(async () => {
    const base = await (await navigator.storage.getDirectory()).getDirectoryHandle('gbdraw-input-copies');
    for await (const folder of base.values()) {
      for await (const handle of folder.values()) return Array.from(new Uint8Array(await (await handle.getFile()).arrayBuffer()));
    }
    return null;
  });
  expect(Buffer.from(copied).equals(gzip)).toBe(true);
  await generate(page);
  const saved = await saveSessionFile(page, 'gzip input');
  expect(saved.genbank.map(({ size, method }) => ({ size, method }))).toEqual([{ size: gzip.byteLength, method: 0 }]);
  expect(saved.genbank[0].bytes.equals(gzip), 'the entry is the gzip file as chosen').toBe(true);
  // A fresh page loads it and draws it again from the stored gzip bytes.
  const later = await page.context().newPage();
  await openFresh(later);
  await loadSessionFile(later, saved.file);
  await generate(later);
  expect(await later.evaluate(() => window.__GBDRAW_APP__.results.length)).toBe(1);
  await later.close();
});

// OV-3364: a LOSATP search stages a card's file under the id the Generate's request names. A gzip input staged with
// another size than the request counts fails the Generate (a cache miss in the Worker).
test('a Linear LOSATP Generate draws with a gzip GenBank input', async ({ page }, info) => {
  test.setTimeout(480_000);
  const gzipped = info.outputPath('NC_001416.gb.gz');
  fs.writeFileSync(gzipped, gzipSync(fs.readFileSync(SOURCE)));
  await routeCrossOriginIsolated(page.context());
  await openFresh(page);
  await showDrawing(page, { kind: 'linear' });
  await page.getByTestId('linear-genbank-1').setInputFiles(SOURCE);
  await page.evaluate(() => window.__GBDRAW_APP__.addLinearSeq());
  await page.getByTestId('linear-genbank-2').setInputFiles(gzipped);
  await page.getByRole('button', { name: 'Run LOSAT for all adjacent pairs', exact: true }).click();
  await page.evaluate(() => window.__GBDRAW_APP__.setLinearComparisonLosatMode('blastp'));
  const outcome = await evaluateWithRetainedPromise(page, async () => {
    const result = await window.__GBDRAW_APP__.runAnalysis();
    return { status: result?.status, code: result?.error?.code ?? null };
  });
  expect(outcome).toEqual({ status: 'ok', code: null });
  expect(await page.evaluate(() => window.__GBDRAW_APP__.results.length)).toBe(1);
});
