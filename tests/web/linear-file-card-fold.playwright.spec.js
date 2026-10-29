// A Linear File card folds to its header line: File N, the file names, the
// problems the card marks, the move buttons, the caret. The header line is the
// summary of a native disclosure that holds the whole card, so folding keeps
// every node. The folded cards are view state keyed by the File's uid: a
// loaded Session shows the cards that hold a file folded, and no Session,
// History step, Undo, Reset Settings, upload, or mode switch holds or changes
// the state.
const { test, expect } = require('@playwright/test');
const { readFileSync } = require('node:fs');
const { join } = require('node:path');
const { gunzipSync } = require('node:zlib');
const { answerSessionTitle, evaluateWithRetainedPromise, openApp } = require('./helpers/app-lifecycle.cjs');

test.use({ actionTimeout: 20_000 });

const FIVE_FILES = 'BGC0000708-BGC0000713.gbdraw-session.json';
const ONE_FILE = 'lambda_basic_linear.gbdraw-session.json';
const TWO_RECORD_FILES = 'vibrio-harveyi-group-collinear.gbdraw-session.json.gz';
const genbank = readFileSync(join(__dirname, '../test_inputs/BGC0000708.gbk'));

const cards = (page) => page.locator('[data-linear-source-card]');
const header = (card) => card.locator(':scope > details > summary');
const chips = (card) => header(card).locator('.linear-file-chip');
const undoCount = (page) => page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount());
const sessionSettled = (page) => page.waitForFunction(() => !window.__GBDRAW_APP__.sessionImportPending, null, { timeout: 180_000 });

const loadSession = async (page, path) => {
  await page.locator('input[accept^=".json,"]').setInputFiles(path);
  await sessionSettled(page);
  await page.waitForFunction(() => window.__GBDRAW_APP__.results.length > 0, null, { timeout: 180_000 });
};
const openWithSession = async (page, name = FIVE_FILES) => {
  await page.setViewportSize({ width: 1600, height: 1000 });
  page.on('dialog', (dialog) => dialog.accept());
  await openApp(page);
  await loadSession(page, join(process.cwd(), 'gbdraw/web/gallery/sessions', name));
};
const openLinear = async (page) => {
  await openApp(page);
  await page.getByRole('button', { name: 'Linear', exact: true }).click();
  await expect(cards(page)).toHaveCount(1);
};
const chooseGenBank = (card, name, content = genbank) => card.locator('input[type="file"]').first()
  .setInputFiles({ name, mimeType: 'text/plain', buffer: Buffer.from(content) });

// Each card in list order: its File uid, whether it is folded, and its height.
// `folded` is the disclosure's state where the owner state agrees with it, and
// a text no expectation matches while the toggle event has not arrived yet.
const foldState = (page) => page.evaluate(async () => {
  const { state } = await import('./js/state.js');
  return Array.from(document.querySelectorAll('[data-linear-source-card]'), (card) => {
    const closed = !card.querySelector(':scope > details').open;
    const owned = state.linearFoldedSources.has(card.dataset.linearSourceCard);
    return {
      uid: card.dataset.linearSourceCard,
      folded: owned === closed ? closed : `disclosure ${closed ? 'closed' : 'open'}, owner state ${owned ? 'folded' : 'open'}`,
      height: Math.round(card.getBoundingClientRect().height)
    };
  });
});
const folded = async (page) => (await foldState(page)).map((card) => card.folded);
// The layout a width change asks for is in place once the pane or the window has that width.
const paneWidth = (page) => page.locator('.settings-pane').evaluate((pane) => Math.round(pane.getBoundingClientRect().width));
// The toggle event that reports a disclosure's state arrives after the click.
const expectFolded = (page, expected) => expect.poll(() => folded(page)).toEqual(expected);
const listHeight = (page) => page.locator('[data-linear-record-list]').evaluate((list) => Math.round(list.getBoundingClientRect().height));
const noSidewaysScroll = (page) => page.evaluate(() => {
  const area = document.querySelector('.settings-scroll');
  const wide = Array.from(document.querySelectorAll('[data-linear-source-card]'))
    .filter((card) => card.scrollWidth > card.clientWidth + 1 || card.getBoundingClientRect().right > area.getBoundingClientRect().right + 1);
  return { area: area.scrollWidth <= area.clientWidth + 1, page: document.documentElement.scrollWidth <= window.innerWidth + 1, wideCards: wide.length };
});

test('a loaded Session shows its File cards folded to one line each, and the header line folds and unfolds a card without rebuilding it', async ({ page }) => {
  test.setTimeout(420_000);
  await openWithSession(page);
  await expect(cards(page)).toHaveCount(5);
  const loaded = await foldState(page);
  expect(loaded.map((card) => card.folded)).toEqual([true, true, true, true, true]);
  for (const card of loaded) expect(card.height, 'a folded card is one line').toBeLessThanOrEqual(40);
  const foldedList = await listHeight(page);
  expect(foldedList).toBeLessThan(260);

  // The folded line: File N, the file name in the uploader's ready look, the move buttons.
  const first = cards(page).first();
  await expect(header(first)).toHaveAttribute('aria-label', 'File 1: BGC0000708.gbk');
  await expect(header(first).locator(':scope > span').first()).toHaveText('File 1');
  await expect(chips(first)).toHaveText(['BGC0000708.gbk']);
  await expect(chips(first)).toHaveAttribute('title', 'BGC0000708.gbk');
  await expect(chips(first)).not.toHaveClass(/failed/);
  await expect(chips(first).locator('i')).toHaveClass(/ph-check-circle/);
  expect(await chips(first).evaluate((chip) => {
    const zone = document.createElement('div');
    zone.className = 'upload-zone ready';
    document.body.append(zone);
    const look = (element) => [getComputedStyle(element).borderTopColor, getComputedStyle(element).backgroundColor];
    const same = JSON.stringify(look(chip)) === JSON.stringify(look(zone));
    zone.remove();
    return same;
  }), 'the chip has the border and background of a ready uploader').toBe(true);
  await expect(first.getByRole('button', { name: 'Move File 1 up', exact: true })).toBeDisabled();
  await expect(first.getByRole('button', { name: 'Move File 1 down', exact: true })).toBeEnabled();
  // The content of a folded card is neither shown, nor a Tab stop, nor in the accessibility tree.
  const snapshot = await first.ariaSnapshot();
  expect(snapshot).not.toMatch(/GenBank \/ DDBJ File|Definition for sequence|Depth tracks|Record options/);
  expect(await first.evaluate((card) => Array.from(card.querySelectorAll('input, select, button, summary, [tabindex]'))
    .filter((element) => element.checkVisibility() && element.tabIndex >= 0 && !element.disabled)
    .map((element) => element.getAttribute('aria-label')))).toEqual(['File 1: BGC0000708.gbk', 'Move File 1 down']);
  await header(first).focus();
  const stops = [];
  for (let index = 0; index < 4; index += 1) {
    await page.keyboard.press('Tab');
    stops.push(await page.evaluate(() => document.activeElement.getAttribute('aria-label')));
  }
  expect(stops).toEqual(['Move File 1 down', 'File 2: BGC0000709.gbk', 'Move File 2 up', 'Move File 2 down']);

  // A press on the header line opens the card and begins no History step.
  const steps = await undoCount(page);
  await header(first).locator(':scope > span').first().click();
  await expectFolded(page, [false, true, true, true, true]);
  await expect(header(first)).toHaveAttribute('aria-label', 'File 1: BGC0000708.gbk');
  await expect(chips(first)).toBeHidden();
  const name = first.getByLabel('Definition for sequence 1', { exact: true });
  await expect(name).toBeVisible();
  expect(await undoCount(page)).toBe(steps);

  // Folding keeps the nodes: a typed value and the inner disclosures survive.
  await name.fill('Typed and kept');
  await name.evaluate((field) => { field.dataset.foldProbe = 'same-node'; });
  const depth = first.locator('[data-linear-source-depth]');
  const options = first.locator('[data-linear-record-options]');
  await depth.locator(':scope > summary').click();
  await options.locator(':scope > summary').click();
  await expect(options).toHaveAttribute('open', '');
  const afterEdit = await undoCount(page);
  await header(first).focus();
  await page.keyboard.press('Enter');
  await expectFolded(page, [true, true, true, true, true]);
  await expect(name).toBeHidden();
  await page.keyboard.press('Space');
  await expectFolded(page, [false, true, true, true, true]);
  await expect(name).toHaveValue('Typed and kept');
  await expect(name).toHaveAttribute('data-fold-probe', 'same-node');
  await expect(depth).toHaveAttribute('open', '');
  await expect(options).toHaveAttribute('open', '');
  expect(await undoCount(page), 'the header line begins no History step').toBe(afterEdit);

  // A move button moves its card and folds or unfolds nothing: the state follows the File.
  await cards(page).nth(1).getByRole('button', { name: 'Move File 2 up', exact: true }).click();
  await expect.poll(async () => (await foldState(page)).map((card) => card.uid)).toEqual(['record-2', 'record-1', 'record-3', 'record-4', 'record-5']);
  await expectFolded(page, [true, false, true, true, true]);
  await cards(page).first().getByRole('button', { name: 'Move File 1 down', exact: true }).focus();
  await page.keyboard.press('Enter');
  await expect.poll(async () => (await foldState(page)).map((card) => card.uid)).toEqual(['record-1', 'record-2', 'record-3', 'record-4', 'record-5']);
  await expectFolded(page, [false, true, true, true, true]);
  // Space on a move button of a folded card moves the card as well, and unfolds nothing.
  for (const [name, order] of [
    ['Move File 2 down', ['record-1', 'record-3', 'record-2', 'record-4', 'record-5']],
    ['Move File 3 up', ['record-1', 'record-2', 'record-3', 'record-4', 'record-5']]
  ]) {
    await page.getByRole('button', { name, exact: true }).focus();
    await page.keyboard.press('Space');
    await expect.poll(async () => (await foldState(page)).map((card) => card.uid)).toEqual(order);
    await expectFolded(page, [false, true, true, true, true]);
  }
  // A press on the disabled move button, or between the two, folds nothing either.
  await first.getByRole('button', { name: 'Move File 1 up', exact: true }).click({ force: true });
  const gap = await first.locator('[data-linear-source-move="up"]').evaluate((button) => {
    const box = button.getBoundingClientRect();
    return { x: box.right + 2, y: box.top + box.height / 2 };
  });
  await page.mouse.click(gap.x, gap.y);
  // A click that folded would have closed the disclosure before it returned.
  await expect(first.locator(':scope > details')).toHaveJSProperty('open', true);
  await expectFolded(page, [false, true, true, true, true]);

  // Every card opened by hand: the list is as tall as it was before cards folded.
  for (const index of [1, 2, 3, 4]) await header(cards(page).nth(index)).press('Enter');
  await expectFolded(page, [false, false, false, false, false]);
  expect(await listHeight(page)).toBeGreaterThan(foldedList * 5);
});

test('Add sequence, an upload, Undo, Redo, a mode switch, and Reset Settings leave the folded cards as they are', async ({ page }) => {
  test.setTimeout(420_000);
  await openWithSession(page);
  await header(cards(page).nth(2)).click();
  await expectFolded(page, [true, true, false, true, true]);

  // A card the user adds is open, and choosing its file does not fold it.
  await page.locator('[data-linear-file-add]').click();
  await expect(cards(page)).toHaveCount(6);
  await expectFolded(page, [true, true, false, true, true, false]);
  const added = cards(page).nth(5);
  await expect(header(added)).toHaveAttribute('aria-label', 'File 6: No file');
  await chooseGenBank(added, 'added-file.gbk');
  await expect(header(added)).toHaveAttribute('aria-label', 'File 6: added-file.gbk');
  await expect(added.locator('.upload-zone.ready')).toBeVisible({ timeout: 120_000 });
  await expectFolded(page, [true, true, false, true, true, false]);

  // Undo and Redo of the upload and the added card change no card that stays.
  await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.undo());
  await expect(header(cards(page).nth(5))).toHaveAttribute('aria-label', 'File 6: No file');
  await expectFolded(page, [true, true, false, true, true, false]);
  await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.undo());
  await expect(cards(page)).toHaveCount(5);
  await expectFolded(page, [true, true, false, true, true]);
  await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.redo());
  await expect(cards(page)).toHaveCount(6);
  await expectFolded(page, [true, true, false, true, true, false]);
  // A card folded by hand comes back folded through Undo and Redo of the step that added it.
  await header(cards(page).nth(5)).click();
  await expectFolded(page, [true, true, false, true, true, true]);
  await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.undo());
  await expect(cards(page)).toHaveCount(5);
  await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.redo());
  await expect(cards(page)).toHaveCount(6);
  await expectFolded(page, [true, true, false, true, true, true]);
  await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.undo());
  await expect(cards(page)).toHaveCount(5);

  // A mode round trip and Reset Settings.
  await page.getByRole('button', { name: 'Circular', exact: true }).click();
  await expect(cards(page)).toHaveCount(0);
  await page.getByRole('button', { name: 'Linear', exact: true }).click();
  await expect(cards(page)).toHaveCount(5);
  await expectFolded(page, [true, true, false, true, true]);
  await page.getByRole('button', { name: 'Reset Settings', exact: true }).click();
  await page.locator('[data-choice-prompt] [data-choice-prompt-action="reset"]').click();
  await expect.poll(() => page.evaluate(() => !window.__GBDRAW_HISTORY__.mutationPending())).toBe(true);
  await expectFolded(page, [true, true, false, true, true]);
});

test('Save Session holds no fold state, Load Session folds the loaded cards, and a failed load folds nothing', async ({ page }, testInfo) => {
  test.setTimeout(420_000);
  await openWithSession(page);
  for (const index of [0, 1, 2, 3, 4]) await header(cards(page).nth(index)).press('Enter');
  await expectFolded(page, [false, false, false, false, false]);

  const pending = page.waitForEvent('download');
  await page.getByRole('banner').getByRole('button', { name: 'Save Session', exact: true }).click();
  const download = await pending;
  const saved = testInfo.outputPath(download.suggestedFilename());
  await download.saveAs(saved);
  const bytes = readFileSync(saved);
  const text = (bytes[0] === 0x1f && bytes[1] === 0x8b ? gunzipSync(bytes) : bytes).toString('utf8');
  const keys = new Set();
  JSON.parse(text, (key, value) => { keys.add(key); return value; });
  expect([...keys].filter((key) => /fold/i.test(key)), 'the Session holds no fold state').toEqual([]);
  await expectFolded(page, [false, false, false, false, false]);

  // A load that fails leaves every card as the user had it.
  await header(cards(page).nth(1)).press('Enter');
  await expectFolded(page, [false, true, false, false, false]);
  const broken = testInfo.outputPath('broken.gbdraw-session.json');
  require('node:fs').writeFileSync(broken, text.slice(0, Math.floor(text.length / 2)));
  await page.locator('input[accept^=".json,"]').setInputFiles(broken);
  await sessionSettled(page);
  await expect(page.locator('[data-session-import-status]')).toHaveCount(0);
  expect(await page.evaluate(() => Boolean(window.__GBDRAW_APP__.errorLog)), 'the broken Session is refused').toBe(true);
  await expect(cards(page)).toHaveCount(5);
  await expectFolded(page, [false, true, false, false, false]);
  // A load that fails after it replaced the cards is rolled back to the same state. The
  // one File of this Session has the uid of the first card, which the user had left open.
  await page.evaluate(() => { window.__GBDRAW_APP__.errorLog = null; });
  const rolledBack = await evaluateWithRetainedPromise(page, async (name) => {
    const file = new File([await (await fetch(`/gbdraw/web/gallery/sessions/${name}`)).arrayBuffer()], name);
    const { importSession } = (await import('/tests/web/helpers/session-caller.mjs')).callerSession(await import('./js/services/config.js'), (await import('./js/state.js')).state);
    const { state } = await import('./js/state.js');
    let foldedInside = null;
    const result = await importSession({ target: { files: [file], value: 'selected' } }, {
      beforePreviewMount: () => { foldedInside = [...state.linearFoldedSources]; throw new Error('forced rollback'); }
    });
    return { status: result.status, foldedInside };
  }, ONE_FILE);
  expect(rolledBack).toEqual({ status: 'error', foldedInside: ['record-1'] });
  await expect(cards(page)).toHaveCount(5);
  await expectFolded(page, [false, true, false, false, false]);

  // The Session saved with every card open loads with its cards folded.
  await page.evaluate(() => { window.__GBDRAW_APP__.errorLog = null; });
  await loadSession(page, saved);
  await expect(cards(page)).toHaveCount(5);
  await expectFolded(page, [true, true, true, true, true]);
});

test('Load example folds the five files of BGC0000708 and the one file of lambda_basic_linear', async ({ page }) => {
  test.setTimeout(420_000);
  await page.setViewportSize({ width: 1600, height: 1000 });
  await openApp(page);
  const chooser = page.getByRole('dialog', { name: 'Choose an example', exact: true });
  await page.getByRole('button', { name: 'Load an example', exact: true }).click();
  await chooser.locator('[data-gallery-example="BGC0000708-BGC0000713"]').click();
  await page.waitForFunction(() => window.__GBDRAW_APP__.results.length > 0 && !window.__GBDRAW_APP__.sessionImportPending, null, { timeout: 180_000 });
  await expect(cards(page)).toHaveCount(5);
  await expectFolded(page, [true, true, true, true, true]);

  await page.locator('header.app-header').getByRole('button', { name: 'Load example', exact: true }).click();
  await chooser.locator('[data-gallery-example="lambda_basic_linear"]').click();
  await expect(cards(page)).toHaveCount(1, { timeout: 180_000 });
  await sessionSettled(page);
  await expectFolded(page, [true]);
  // A single card has no move buttons: File N, the name, the caret.
  const only = cards(page).first();
  await expect(only.locator('[data-linear-source-move]')).toHaveCount(0);
  await expect(chips(only)).toHaveCount(1);
  expect((await foldState(page))[0].height).toBeLessThanOrEqual(40);
});

test('Show in settings of a record opens its folded card and focuses the name field, in the stacked display and in the rail', async ({ page }) => {
  test.setTimeout(420_000);
  await openWithSession(page);
  await expectFolded(page, [true, true, true, true, true]);
  const showInSettings = async () => {
    const definition = page.locator('.gbdraw-preview-surface [data-gbdraw-role^="record-definition"] text').first();
    const box = await definition.boundingBox();
    await page.mouse.click(box.x + box.width / 2, box.y + box.height / 2);
    await page.getByRole('button', { name: 'Show in settings', exact: true }).click();
  };
  const field = page.getByLabel('Definition for sequence 1', { exact: true });
  await showInSettings();
  await expect(field).toBeFocused();
  await expectFolded(page, [false, true, true, true, true]);

  await header(cards(page).first()).click();
  await expectFolded(page, [true, true, true, true, true]);
  await page.getByRole('button', { name: 'Collapse settings to icons', exact: true }).click();
  await page.getByRole('button', { name: 'Close section', exact: true }).click();
  await showInSettings();
  await expect(field).toBeFocused();
  await expect(page.locator('.settings-pane')).toHaveAttribute('data-settings-docked', 'input-genomes');
  await expectFolded(page, [false, true, true, true, true]);
});

test('a File with several records keeps its record list as it was through a fold, and Show in settings opens into it', async ({ page }) => {
  test.setTimeout(420_000);
  await openWithSession(page, TWO_RECORD_FILES);
  await expectFolded(page, [true, true]);
  const card = cards(page).first();
  const records = card.locator('[data-linear-source-records]');
  await header(card).click();
  await expectFolded(page, [false, true]);
  await expect(records.locator(':scope > summary')).toContainText('Number of records: 2');
  await expect(records).toHaveJSProperty('open', false);
  // Opened by hand, the list is the same node and still open after the card was folded and opened.
  await records.locator(':scope > summary').click();
  await expect(records).toHaveJSProperty('open', true);
  await records.evaluate((node) => { node.dataset.foldProbe = 'kept'; });
  await header(card).click();
  await expectFolded(page, [true, true]);
  await expect(records).toBeHidden();
  await header(card).click();
  await expectFolded(page, [false, true]);
  await expect(records).toHaveAttribute('data-fold-probe', 'kept');
  await expect(records).toHaveJSProperty('open', true);
  // Closed by hand, it stays closed through a fold, and Show in settings opens the card and the list.
  await records.locator(':scope > summary').click();
  await expect(records).toHaveJSProperty('open', false);
  await header(card).click();
  await expectFolded(page, [true, true]);
  const definition = page.locator('.gbdraw-preview-surface [data-gbdraw-role^="record-definition"] text').first();
  const box = await definition.boundingBox();
  await page.mouse.click(box.x + box.width / 2, box.y + box.height / 2);
  await page.getByRole('button', { name: 'Show in settings', exact: true }).click();
  await expectFolded(page, [false, true]);
  await expect(records).toHaveJSProperty('open', true);
  await expect(records.locator('[data-linear-record-names] input:focus')).toHaveCount(1);
});

// Save Session refuses a GFF3 without its FASTA, so no saved Session holds such a card; the
// load rule for one is pinned in tests/web/linear-sources.test.mjs.
test('a loaded GFF3 + FASTA Session folds a card that holds both of its files', async ({ page }, testInfo) => {
  test.setTimeout(420_000);
  await page.setViewportSize({ width: 1600, height: 1000 });
  page.on('dialog', (dialog) => dialog.accept());
  await openLinear(page);
  await page.getByRole('radio', { name: 'GFF3 + FASTA' }).check();
  const card = cards(page).first();
  await card.locator('input[type="file"]').nth(0).setInputFiles(join(__dirname, '../test_inputs/NC_013668.gff3'));
  await expect(header(card)).toHaveAttribute('aria-label', 'File 1: NC_013668.gff3 + no FASTA');
  await card.locator('input[type="file"]').nth(1).setInputFiles(join(__dirname, '../test_inputs/NC_013668.fasta'));
  await expect(header(card)).toHaveAttribute('aria-label', 'File 1: NC_013668.gff3 + NC_013668.fasta');
  await expectFolded(page, [false]);

  const pending = page.waitForEvent('download');
  await page.getByRole('banner').getByRole('button', { name: 'Save Session', exact: true }).click();
  await answerSessionTitle(page, '');
  const download = await pending;
  const saved = testInfo.outputPath(download.suggestedFilename());
  await download.saveAs(saved);
  // A fresh page is on GenBank: the load rule reads the input type of the Session it loads.
  await openLinear(page);
  await expect(page.getByRole('radio', { name: 'GenBank', exact: true }).first()).toBeChecked();
  await page.locator('input[accept^=".json,"]').setInputFiles(saved);
  await sessionSettled(page);
  await expect(cards(page)).toHaveCount(1);
  await expectFolded(page, [true]);
  await expect(chips(card)).toHaveText(['NC_013668.gff3 + NC_013668.fasta']);
});

test('removing a File moves the focus to the header line of a folded neighbour and leaves it folded', async ({ page }) => {
  test.setTimeout(420_000);
  await openWithSession(page);
  const dialog = page.locator('[data-linear-source-removal-dialog]');
  await header(cards(page).nth(1)).click();
  await expectFolded(page, [true, false, true, true, true]);

  // Cancel returns to the Remove button of the open card.
  const remove = cards(page).nth(1).getByRole('button', { name: 'Remove GenBank / DDBJ File', exact: true });
  await remove.click();
  await dialog.getByRole('button', { name: 'Cancel', exact: true }).click();
  await expect(remove).toBeFocused();

  // Remove file: the next card takes the removed card's place, folded, with the focus on its header line.
  await remove.click();
  await dialog.getByRole('button', { name: 'Remove file', exact: true }).click();
  await expect(cards(page)).toHaveCount(4);
  await expect(header(cards(page).nth(1))).toBeFocused();
  await expect(header(cards(page).nth(1))).toHaveAttribute('aria-label', 'File 2: BGC0000711.gbk');
  await expectFolded(page, [true, true, true, true]);

  // The Remove button under the list removes the last File: the focus goes to the new last card.
  await page.getByRole('button', { name: 'Remove last sequence', exact: true }).click();
  await dialog.getByRole('button', { name: 'Remove file', exact: true }).click();
  await expect(cards(page)).toHaveCount(3);
  await expect(header(cards(page).nth(2))).toBeFocused();
  await expectFolded(page, [true, true, true]);
  // Undo brings the removed File back as it was: folded.
  await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.undo());
  await expect(cards(page)).toHaveCount(4);
  await expect(header(cards(page).nth(3))).toHaveAttribute('aria-label', 'File 4: BGC0000713.gbk');
  await expectFolded(page, [true, true, true, true]);
  await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.redo());
  await expect(cards(page)).toHaveCount(3);
  await expectFolded(page, [true, true, true]);

  // Empty this slot: the card's place holds a new, open card, and its uploader takes the focus.
  await header(cards(page).first()).click();
  await cards(page).first().getByRole('button', { name: 'Remove GenBank / DDBJ File', exact: true }).click();
  await dialog.getByRole('button', { name: 'Empty this slot', exact: true }).click();
  await expect(cards(page)).toHaveCount(3);
  await expectFolded(page, [false, true, true]);
  await expect(cards(page).first().getByRole('button', { name: 'Choose GenBank / DDBJ File', exact: true })).toBeFocused();
});

test('the folded line names both files of a GFF3 + FASTA card, cuts a long name, and says when a card has no file', async ({ page }) => {
  test.setTimeout(300_000);
  await page.setViewportSize({ width: 1600, height: 1000 });
  await openLinear(page);
  const card = cards(page).first();
  // The first card of a fresh page is open; folded by hand it says that it holds no file.
  await expectFolded(page, [false]);
  await header(card).click();
  await expectFolded(page, [true]);
  await expect(header(card)).toHaveAttribute('aria-label', 'File 1: No file');
  await expect(header(card).locator('.linear-file-fold-names')).toHaveText('No file');
  await expect(chips(card)).toHaveCount(0);
  await header(card).click();

  // A long name is cut with an ellipsis, stays in the title, and never widens the pane.
  const longName = `${'a-very-long-genbank-file-name-'.repeat(6)}.gbk`;
  await chooseGenBank(card, longName);
  await expect(card.locator('.upload-zone.ready')).toBeVisible({ timeout: 120_000 });
  await expectFolded(page, [false]);
  await header(card).click();
  await expectFolded(page, [true]);
  await expect(chips(card)).toHaveAttribute('title', longName);
  await expect(header(card)).toHaveAttribute('aria-label', `File 1: ${longName}`);
  expect(await chips(card).locator('span').evaluate((text) => text.scrollWidth > text.clientWidth), 'the name is cut').toBe(true);
  await page.locator('[data-linear-file-add]').click();
  await header(cards(page).nth(1)).click();
  await expectFolded(page, [true, true]);
  for (const width of [240, 320, 500]) {
    await page.evaluate((value) => { window.__GBDRAW_APP__.sidebarWidth = value; }, width);
    await expect.poll(() => paneWidth(page)).toBe(width);
    expect(await noSidewaysScroll(page), `pane width ${width}`).toEqual({ area: true, page: true, wideCards: 0 });
    expect((await foldState(page)).map((entry) => entry.height <= 40), `one line at ${width}`).toEqual([true, true]);
  }
  await page.setViewportSize({ width: 390, height: 844 });
  await expect.poll(() => paneWidth(page)).toBeLessThanOrEqual(390);
  // The header takes its two-row layout in the page's resize event, which can run after setViewportSize returns.
  await expect.poll(() => noSidewaysScroll(page), { message: 'phone width' }).toEqual({ area: true, page: true, wideCards: 0 });
  await page.setViewportSize({ width: 1600, height: 1000 });

  // GFF3 + FASTA. A card that holds only a file of the other input type reads as a card without a file.
  await page.evaluate(() => { window.__GBDRAW_APP__.sidebarWidth = 400; });
  await page.getByRole('radio', { name: 'GFF3 + FASTA' }).check();
  const pair = cards(page).first();
  await expect(header(pair)).toHaveAttribute('aria-label', 'File 1: No file');
  await expect(chips(pair)).toHaveCount(0);
  // One file of the pair: the line names the missing one and does not look ready.
  await header(pair).click();
  await pair.locator('input[type="file"]').nth(0).setInputFiles({ name: 'genome.gff3', mimeType: 'text/plain', buffer: Buffer.from('not a gff\n') });
  await expect(pair.locator('.upload-zone.ready')).toHaveCount(1);
  await header(pair).click();
  await expectFolded(page, [true, true]);
  await expect(header(pair)).toHaveAttribute('aria-label', 'File 1: genome.gff3 + no FASTA');
  await expect(header(pair).locator('.linear-file-fold-names')).toHaveText('genome.gff3 + no FASTA');
  await expect(chips(pair)).toHaveCount(0);
  await expect(header(pair).locator('.ph-check-circle')).toHaveCount(0);
  // Both files, and the failed look of a file whose inspection failed.
  await header(pair).click();
  await pair.locator('input[type="file"]').nth(1).setInputFiles({ name: 'genome.fasta', mimeType: 'text/plain', buffer: Buffer.from('not a fasta\n') });
  await expect(pair.locator('.upload-zone.failed')).toHaveCount(2, { timeout: 120_000 });
  await header(pair).click();
  await expect(chips(pair)).toHaveText(['genome.gff3 + genome.fasta']);
  // The failed state is also a description of the header line, not only a colour and an icon.
  expect(await header(pair).evaluate((summary) => summary.getAttribute('aria-describedby').split(' ')
    .map((id) => document.getElementById(id)?.textContent))).toEqual(['This file could not be read.']);
  await expect(header(cards(page).nth(1))).not.toHaveAttribute('aria-describedby', /.+/);
  await expect(chips(pair)).toHaveAttribute('title', 'genome.gff3 + genome.fasta');
  await expect(header(pair)).toHaveAttribute('aria-label', 'File 1: genome.gff3 + genome.fasta');
  for (const chip of await chips(pair).all()) {
    await expect(chip).toHaveClass(/failed/);
    await expect(chip.locator('i')).toHaveClass(/ph-warning-circle/);
  }
  expect(await chips(pair).first().evaluate((chip) => {
    const zone = document.createElement('div');
    zone.className = 'upload-zone failed';
    document.body.append(zone);
    const same = getComputedStyle(chip).borderTopColor === getComputedStyle(zone).borderTopColor
      && getComputedStyle(chip).backgroundColor === getComputedStyle(zone).backgroundColor;
    zone.remove();
    return same;
  }), 'the chip has the look of a failed uploader').toBe(true);
});

test('a folded card that TLOSATX cannot search says Gencode required, as its open card does', async ({ page }) => {
  test.setTimeout(420_000);
  await openWithSession(page);
  const card = cards(page).first();
  await expect(header(card)).not.toContainText('Gencode required');
  await expect(header(card)).not.toHaveAttribute('aria-describedby', /.+/);
  await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    app.losatProgram = 'tblastx';
    app.linearSeqs[0].losat_gencode = '';
  });
  const mark = header(card).getByText('Gencode required', { exact: true });
  await expect(mark).toBeVisible();
  expect(await header(card).evaluate((summary) => document.getElementById(summary.getAttribute('aria-describedby'))?.textContent))
    .toBe('Gencode required');
  await expect(header(cards(page).nth(1))).not.toContainText('Gencode required');
  // The mark is whole: beside the name where the line has the room, else on a second line.
  expect(await mark.evaluate((element) => element.scrollWidth <= element.clientWidth + 1)).toBe(true);
  expect(await chips(card).evaluate((chip) => chip.getBoundingClientRect().width), 'the name keeps its room').toBeGreaterThan(90);
  expect((await foldState(page))[0].height).toBeLessThanOrEqual(56);
  // On the narrowest pane the mark wraps inside the line and stays clear of the move buttons (OV-602).
  const paneWidth = await page.evaluate(() => window.__GBDRAW_APP__.sidebarWidth);
  await page.evaluate(() => { window.__GBDRAW_APP__.sidebarWidth = 240; });
  await expect.poll(() => header(card).evaluate((summary) => summary.querySelector('[id^="linear-file-gencode-"]').getBoundingClientRect().right
    <= summary.querySelector('[data-linear-source-move]').getBoundingClientRect().left)).toBe(true);
  await page.evaluate((value) => { window.__GBDRAW_APP__.sidebarWidth = value; }, paneWidth);
  // The open card says it on its record options, from the same predicate, and the header line no longer does.
  await header(card).click();
  await expect(mark).toBeHidden();
  await expect(card.locator('[data-linear-record-options] > summary')).toContainText('Gencode required');
  await page.evaluate(() => { window.__GBDRAW_APP__.linearSeqs[0].losat_gencode = 11; });
  await expect(card.locator('[data-linear-record-options] > summary')).not.toContainText('Gencode required');
  await header(card).click();
  await expect(header(card)).not.toContainText('Gencode required');
});
