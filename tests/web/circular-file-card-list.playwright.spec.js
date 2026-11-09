const { test, expect } = require('@playwright/test');
const { execFile } = require('node:child_process');
const { promisify } = require('node:util');
const path = require('node:path');
const { readFile } = require('node:fs/promises');
const { openApp, evaluateWithRetainedPromise } = require('./helpers/app-lifecycle.cjs');
const { circularFileInput } = require('./helpers/circular-inputs.cjs');

// The Circular input is a list of File cards in the shell the Linear File
// cards use (circular-multi-file D-01, D-10): no card and the Add files zone
// while there is no file; a pick or a drop of several files on the zone adds a
// card per file (a GFF3 + FASTA pair by file-name stem) in one History step;
// each card has Replace, its moves and Remove, which asks first. The pure parts
// are in tests/web/circular-file-card-list.test.mjs.
const FIXTURES = path.resolve('tests/fixtures');
const twoRecords = path.join(FIXTURES, 'web_batch_two_records.gb');
const mito = path.join(FIXTURES, 'sessions', 'cli-web-mito.gb');
const small = path.join(FIXTURES, 'feature_location_search.gb');
const regex = path.join(FIXTURES, 'regex_rules.gb');

const cards = (page) => page.locator('[data-circular-source-card]');
const addZone = (page) => page.locator('[data-circular-file-add]');
const cardNames = (page) => cards(page).locator('summary').evaluateAll((summaries) => summaries.map((summary) => summary.getAttribute('aria-label')));

test('Circular File cards: empty state, add, replace, move, and remove with its dialog', async ({ page }) => {
  test.setTimeout(300_000);
  await openApp(page);
  const section = page.locator('[data-setting=inputs]');
  // Empty: no card, the Add files zone, no input Generate can read.
  await expect(cards(page)).toHaveCount(0);
  await expect(section.getByRole('button', { name: 'Choose GenBank/DDBJ files' })).toBeVisible();
  expect(await page.evaluate(() => window.__GBDRAW_APP__.shownInputComplete('circular'))).toBe(false);

  await section.locator('input[type=file]').first().setInputFiles(small);
  await expect(cards(page)).toHaveCount(1);
  for (const file of [twoRecords, mito]) await addZone(page).locator('input[type=file]').setInputFiles(file);
  await expect.poll(() => cardNames(page)).toEqual([
    'File 1: feature_location_search.gb', 'File 2: web_batch_two_records.gb', 'File 3: cli-web-mito.gb'
  ]);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.shownInputComplete('circular'))).toBe(true);
  // A card shows folded with its file chip; once read, the open card says what it holds.
  await expect(page.locator('[data-circular-discovery-status]')).toContainText('4 source records inspected', { timeout: 180_000 });
  const second = cards(page).nth(1);
  await expect(second.locator('.file-card-chip')).toHaveText('web_batch_two_records.gb');
  await second.locator('summary').click({ position: { x: 8, y: 8 } });
  await expect(second.locator('[data-circular-file-records]')).toHaveText(/^2 records · [\d,]+ bp$/);
  // The Topology card names its file and, for a file of several records, its #k.
  await expect(page.locator('[data-record-rotation]').nth(2)).toContainText('web_batch_two_records.gb #2');

  // Replace: the open card's slot takes another file; the card keeps its place.
  await expect(second.getByRole('button', { name: 'Replace GenBank/DDBJ File' })).toBeVisible();
  await second.locator('input[type=file]').setInputFiles(regex);
  await expect.poll(() => cardNames(page)).toEqual([
    'File 1: feature_location_search.gb', 'File 2: regex_rules.gb', 'File 3: cli-web-mito.gb'
  ]);

  // Move File 3 up; File 1 cannot move up.
  await expect(cards(page).nth(0).getByRole('button', { name: 'Move File 1 up' })).toBeDisabled();
  await cards(page).nth(2).getByRole('button', { name: 'Move File 3 up' }).click();
  await expect.poll(() => cardNames(page)).toEqual([
    'File 1: feature_location_search.gb', 'File 2: cli-web-mito.gb', 'File 3: regex_rules.gb'
  ]);

  // Remove asks first; Cancel keeps the card, Remove file removes it.
  await cards(page).nth(1).getByRole('button', { name: 'Remove File 2' }).click();
  const dialog = page.locator('[data-circular-source-removal-dialog] [role=dialog]');
  await expect(dialog).toContainText('Remove cli-web-mito.gb?');
  await dialog.getByRole('button', { name: 'Cancel' }).click();
  await expect(dialog).toHaveCount(0);
  await expect(cards(page)).toHaveCount(3);
  await cards(page).nth(1).getByRole('button', { name: 'Remove File 2' }).click();
  await dialog.getByRole('button', { name: 'Remove file' }).click();
  await expect.poll(() => cardNames(page)).toEqual(['File 1: feature_location_search.gb', 'File 2: regex_rules.gb']);
  // The removal is one History step.
  await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.undo());
  await expect.poll(() => cardNames(page)).toEqual([
    'File 1: feature_location_search.gb', 'File 2: cli-web-mito.gb', 'File 3: regex_rules.gb'
  ]);

  // Choose records… groups the records by file.
  await page.locator('[data-record-list-open=circular]').click();
  await expect(page.locator('[data-record-list-group]')).toHaveText([
    'feature_location_search.gb', 'cli-web-mito.gb', 'regex_rules.gb'
  ]);
});

// D-07: one drop of three GenBank files is three cards in the order dropped and one History step; one pick of
// GFF3 and FASTA files pairs them by file-name stem, an unpaired GFF3 shows its empty FASTA slot; the zone keeps
// the focus (OV-2962); a FASTA among GenBank files fails its own card and the read error names it.
test('Add files: a drop of several files is one card each and one undo; GFF3 + FASTA pairs by stem', async ({ page }) => {
  test.setTimeout(300_000);
  await openApp(page);
  const zone = addZone(page).getByRole('button', { name: 'Choose GenBank/DDBJ files' });
  await expect(zone).toContainText('Add GenBank/DDBJ files: choose several or drop them here');
  const dropped = await page.evaluateHandle((files) => {
    const transfer = new DataTransfer();
    files.forEach(([name, text]) => transfer.items.add(new File([text], name, { type: 'text/plain' })));
    return transfer;
  }, await Promise.all([small, twoRecords, mito].map(async (file) => [path.basename(file), await readFile(file, 'utf8')])));
  const undoBefore = await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount());
  await zone.dispatchEvent('dragover', { dataTransfer: dropped });
  await expect(zone).toHaveClass(/drag-over/);
  await zone.dispatchEvent('drop', { dataTransfer: dropped });
  await expect.poll(() => cardNames(page)).toEqual([
    'File 1: feature_location_search.gb', 'File 2: web_batch_two_records.gb', 'File 3: cli-web-mito.gb'
  ]);
  await expect(zone).not.toHaveClass(/drag-over/);
  await expect(page.locator('[data-circular-discovery-status]')).toContainText('4 source records inspected', { timeout: 180_000 });
  expect(await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount())).toBe(undoBefore + 1);
  await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.undo());
  await expect(cards(page)).toHaveCount(0);
  await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.redo());
  await expect(cards(page)).toHaveCount(3);

  // A FASTA picked as GenBank fails its own card; the other cards still list and the error names the file.
  await addZone(page).locator('input[type=file]').setInputFiles('tests/test_inputs/NC_013668.fasta');
  await expect(page.locator('[data-circular-discovery-status]')).toContainText('NC_013668.fasta: No records were found.', { timeout: 180_000 });
  await expect(cards(page).nth(3).locator('.file-card-chip.failed')).toHaveText('NC_013668.fasta');
  await expect(cards(page).locator('.file-card-chip:not(.failed)')).toHaveCount(3);

  // GFF3 + FASTA: one pick pairs by stem; the unpaired GFF3 shows open with its empty FASTA slot.
  for (let n = 4; n > 0; n -= 1) {
    await cards(page).first().getByRole('button', { name: 'Remove File 1' }).click();
    await page.locator('[data-circular-source-removal-dialog] [role=dialog]').getByRole('button', { name: 'Remove file' }).click();
    await expect(cards(page)).toHaveCount(n - 1);
  }
  await page.getByLabel('GFF3 + FASTA', { exact: true }).check();
  const gffZone = addZone(page).getByRole('button', { name: 'Choose GFF3 and FASTA files' });
  await gffZone.focus();
  const gff = await readFile('tests/test_inputs/NC_013668.gff3');
  const fasta = await readFile('tests/test_inputs/NC_013668.fasta');
  await addZone(page).locator('input[type=file]').setInputFiles([
    { name: 'NC_013668.fasta', mimeType: 'text/plain', buffer: fasta },
    { name: 'NC_013668.gff3', mimeType: 'text/plain', buffer: gff },
    { name: 'NC_013668-reannotated.gff3', mimeType: 'text/plain', buffer: gff }
  ]);
  await expect.poll(() => cardNames(page)).toEqual([
    'File 1: NC_013668.gff3 + NC_013668.fasta', 'File 2: NC_013668-reannotated.gff3 + no FASTA'
  ]);
  await expect(cards(page).nth(1).locator('[data-circular-source-fold]')).toHaveAttribute('open', '');
  await expect(cards(page).nth(1).getByRole('button', { name: 'Choose FASTA File', exact: true })).toBeVisible();
  await expect(gffZone).toBeFocused();
});

// OV-2840: a CLI Session of several --gbk files shows one card per file, each named after its own file.
test('a CLI Session of several --gbk files loads as one File card per file', async ({ page }, testInfo) => {
  test.setTimeout(600_000);
  const output = testInfo.outputPath('cli.gbdraw-session.json.gz');
  await promisify(execFile)('python', ['-m', 'gbdraw.cli', 'circular', '--gbk', twoRecords, mito, small,
    '--multi_record_canvas', '-o', testInfo.outputPath('cli'), '--session_output', output],
  { cwd: testInfo.outputDir, env: { ...process.env, PYTHONPATH: process.cwd() }, timeout: 300_000 });
  page.on('dialog', (dialog) => dialog.accept());
  await openApp(page);
  await page.locator('input[accept^=".json,"]').setInputFiles(output);
  await expect.poll(() => page.evaluate(() => !window.__GBDRAW_APP__.sessionImportPending), { timeout: 180_000 }).toBe(true);
  await expect.poll(() => cardNames(page)).toEqual([
    'File 1: web_batch_two_records.gb', 'File 2: cli-web-mito.gb', 'File 3: feature_location_search.gb'
  ]);
  await expect(cards(page).locator('.file-card-chip')).toHaveText([
    'web_batch_two_records.gb', 'cli-web-mito.gb', 'feature_location_search.gb'
  ]);
});
