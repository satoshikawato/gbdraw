// The Circular Ring colors setting (PD-OI-110 scenario revision 3): a CLI
// Session opens drawing the colours the command line drew (OV-3054), a ring's
// own colour wins on every scheme, and re-assigning changes scheme colours only,
// in the drawn order.
const { test, expect } = require('@playwright/test');
const { execFileSync } = require('node:child_process');
const { readFileSync, writeFileSync } = require('node:fs');
const path = require('node:path');
const { generateAndWaitForResult, reveal } = require('./helpers/app-lifecycle.cjs');
const { openFresh, loadSessionFile, settle } = require('./helpers/audit-browser.cjs');
const { download } = require('./helpers/mode-transition.cjs');
const { shownSavedDrawing } = require('./helpers/saved-session.cjs');
const { readSessionFile } = require('./helpers/session-file.cjs');

const root = path.resolve(process.env.GBDRAW_REPO || process.cwd());
const RING_FIXTURE = 'tests/fixtures/forced_label_underlay.gb';

// Three rings on the single-record fixture (subject FORCEDLBL), each with three hits.
const ringTable = (ring) => [40 + 10 * ring, 60 + 10 * ring, 95].map((identity, index) => [
  `RING${ring}`, 'FORCEDLBL', identity.toFixed(1), '300', '0', '0', String(1 + 300 * index), String(300 * (index + 1)),
  String(1 + 300 * index), String(300 * (index + 1)), '1e-50', '200'
].join('\t')).join('\n') + '\n';

// The Session and the ring colours (by source index) the command line writes with `--conservation_color_scheme`.
const cliRingSession = (testInfo, scheme) => {
  const tables = [1, 2, 3].map((ring) => {
    const file = testInfo.outputPath(`ring${ring}.blast`);
    writeFileSync(file, ringTable(ring));
    return file;
  });
  const session = testInfo.outputPath(`rings-${scheme}.gbdraw`);
  execFileSync(process.env.PYTHON || 'python', [
    '-m', 'gbdraw.cli', 'circular', '--gbk', path.resolve(root, RING_FIXTURE), '--conservation_blast', ...tables,
    '--conservation_color_scheme', scheme, '-o', testInfo.outputPath(`rings-${scheme}`), '-f', 'svg', '--session_output', session
  ], { cwd: root, stdio: 'pipe', env: { ...process.env, PYTHONPATH: root } });
  return { session, svg: readFileSync(testInfo.outputPath(`rings-${scheme}.svg`), 'utf8') };
};

// Each drawn ring, outermost first: its colour (`data-track-color`, absent for a ring without a colour of its own)
// and the fills of its hits; of the SVG text given, else of the shown Result.
const ringView = (text) => {
  const svg = text === null
    ? window.__GBDRAW_APP__.svgContainer.querySelector('svg')
    : new DOMParser().parseFromString(text, 'image/svg+xml').documentElement;
  return [...svg.querySelectorAll('[data-gbdraw-slot-renderer="sequence_conservation"]')]
    .sort((a, b) => Number(a.dataset.sourceIndex) - Number(b.dataset.sourceIndex))
    .map((slot) => ({
      color: slot.getAttribute('data-track-color'),
      fills: [...new Set([...slot.querySelectorAll('[fill]')].map((node) => node.getAttribute('fill')))].sort()
    }));
};
const ringViewOf = (page, svgText) => page.evaluate(ringView, svgText);
const generatedRingView = (page) => page.evaluate(ringView, null);
const generatedRingColors = async (page) => (await generatedRingView(page)).map((ring) => ring.color);

const schemes = (page) => page.evaluate(async () => {
  const module = await import('./js/utils/conservation-colors.generated.js');
  return { cycles: module.CONSERVATION_COLOR_SCHEMES, same: module.CONSERVATION_SAME_COLOR };
});
const swatches = (page) => page.locator('input[type=color][aria-label^="Comparison ring color for"]')
  .evaluateAll((inputs) => inputs.map((input) => input.value));
const ringFiles = (page) => page.evaluate(() => window.__GBDRAW_APP__.circularConservationSeriesRows.map((row) => row.filename));
const ownColors = (page) => page.evaluate(() => window.__GBDRAW_APP__.circularConservation.series.map((entry) => entry.color));
const ringControls = (page) => reveal(page.getByTestId('circular-ring-colors'));
const pickScheme = async (page, scheme) => {
  await ringControls(page);
  await page.locator(`input[name="circular-ring-color-scheme"][value="${scheme}"]`).check();
  await settle(page);
};
const historyStep = async (page, step) => {
  await page.evaluate((name) => window.__GBDRAW_HISTORY__[name](), step);
  await settle(page);
};
const savedRings = async (page, testInfo, name) => {
  const file = testInfo.outputPath(`${name}.gbdraw`);
  // A Session opened from a file has its title, so Save asks for none.
  const titled = await page.evaluate(async () => Boolean((await import('./js/state.js')).state.sessionTitle?.value));
  await download(page, 'Save Session', file, titled ? {} : { title: name });
  const saved = await readSessionFile(file);
  return shownSavedDrawing(saved).config.circularConservation;
};

for (const scheme of ['ten', 'same']) {
  test(`ring colors: a CLI Session saved with ${scheme} opens and generates the colours the CLI drew, and saves none of its own (OV-3054)`, async ({ page }, testInfo) => {
    test.setTimeout(300_000);
    const cli = cliRingSession(testInfo, scheme);
    await openFresh(page);
    const drawn = await ringViewOf(page, cli.svg);
    const { same } = await schemes(page);
    await loadSessionFile(page, cli.session);
    await ringControls(page);
    await expect(page.locator(`input[name="circular-ring-color-scheme"][value="${scheme}"]`)).toBeChecked();
    expect(await swatches(page), 'the panel shows the drawn colours')
      .toEqual(drawn.map((ring) => ring.color ?? same));
    await expect(page.getByTestId('circular-ring-use-scheme-color')).toHaveCount(0);
    await generateAndWaitForResult(page);
    expect(await generatedRingView(page), 'Generate draws what the CLI drew').toEqual(drawn);
    const saved = await savedRings(page, testInfo, `saved-${scheme}`);
    expect(saved.colorScheme).toBe(scheme);
    expect(saved.series.map((entry) => entry.color)).toEqual([null, null, null]);
  });
}

test('ring colors: a hand-set colour wins on every scheme; re-assigning follows the drawn order, by the button or on a move with the checkbox on', async ({ page }, testInfo) => {
  test.setTimeout(420_000);
  const cli = cliRingSession(testInfo, 'ten');
  await openFresh(page);
  await loadSessionFile(page, cli.session);
  const { cycles, same } = await schemes(page);
  const files = await ringFiles(page);
  const ten = (...indexes) => indexes.map((index) => cycles.ten[index]);
  // The scheme index of each ring, by file; and the colour each ring draws after a Generate, by file.
  const schemeIndexes = () => page.evaluate(() => Object.fromEntries(window.__GBDRAW_APP__.circularConservation.series
    .map((entry) => [entry.fileName, entry.schemeIndex])));
  const drawnByFile = async () => {
    await generateAndWaitForResult(page);
    const view = await generatedRingView(page);
    const rows = await page.evaluate(() => window.__GBDRAW_APP__.circularConservation.series.map((entry) => entry.fileName));
    return Object.fromEntries(rows.map((file, index) => [file, view[index]?.color]));
  };

  // A colour set by hand on the first ring stays through every scheme; the others take the scheme's.
  const color = page.getByLabel(`Comparison ring color for ${files[0]}`, { exact: true });
  await reveal(color);
  await color.fill('#123456');
  await settle(page);
  expect(await swatches(page)).toEqual(['#123456', cycles.ten[1], cycles.ten[2]]);
  await expect(page.getByTestId('circular-ring-use-scheme-color')).toHaveCount(1);
  await pickScheme(page, 'same');
  expect(await swatches(page)).toEqual(['#123456', same, same]);
  await expect(page.getByTestId('circular-ring-colors-reassign')).toBeDisabled();
  await expect(page.getByTestId('circular-ring-colors-reassign-on-reorder')).toBeDisabled();
  await pickScheme(page, 'per_five');
  expect(await swatches(page)).toEqual(['#123456', cycles.per_five[1], cycles.per_five[2]]);
  await historyStep(page, 'undo');
  await expect(page.locator('input[name="circular-ring-color-scheme"][value="same"]'), 'a scheme change is one History step').toBeChecked();
  expect(await swatches(page)).toEqual(['#123456', same, same]);
  await pickScheme(page, 'ten');

  // With the checkbox off, a move keeps each ring's colour (and its index); the button re-assigns in the new order.
  const outward = page.getByRole('button', { name: 'Move series outward', exact: true });
  await outward.nth(2).click();
  await settle(page);
  expect(await ringFiles(page)).toEqual([files[0], files[2], files[1]]);
  expect(await schemeIndexes()).toEqual({ [files[0]]: 0, [files[1]]: 1, [files[2]]: 2 });
  expect(await swatches(page)).toEqual(['#123456', ...ten(2, 1)]);
  expect(await drawnByFile(), 'Generate draws the colour that stayed with each ring')
    .toEqual({ [files[0]]: '#123456', [files[1]]: cycles.ten[1], [files[2]]: cycles.ten[2] });
  await page.getByTestId('circular-ring-colors-reassign').click();
  await settle(page);
  expect(await schemeIndexes()).toEqual({ [files[0]]: 0, [files[1]]: 2, [files[2]]: 1 });
  expect(await swatches(page)).toEqual(['#123456', ...ten(1, 2)]);
  // Undo before the Generate: a Generate is a History step of its own.
  await historyStep(page, 'undo');
  expect(await schemeIndexes(), 'the button is one History step').toEqual({ [files[0]]: 0, [files[1]]: 1, [files[2]]: 2 });
  await historyStep(page, 'redo');
  expect(await schemeIndexes()).toEqual({ [files[0]]: 0, [files[1]]: 2, [files[2]]: 1 });
  expect(await drawnByFile()).toEqual({ [files[0]]: '#123456', [files[1]]: cycles.ten[2], [files[2]]: cycles.ten[1] });

  // With the checkbox on, a move re-assigns in the move's step; the hand-set colour stays.
  await page.getByTestId('circular-ring-colors-reassign-on-reorder').check();
  await settle(page);
  await historyStep(page, 'undo');
  await expect(page.getByTestId('circular-ring-colors-reassign-on-reorder'), 'the checkbox is one History step').not.toBeChecked();
  await historyStep(page, 'redo');
  await expect(page.getByTestId('circular-ring-colors-reassign-on-reorder')).toBeChecked();
  await outward.nth(2).click();
  await settle(page);
  expect(await ringFiles(page)).toEqual([files[0], files[1], files[2]]);
  expect(await schemeIndexes()).toEqual({ [files[0]]: 0, [files[1]]: 1, [files[2]]: 2 });
  expect(await swatches(page)).toEqual(['#123456', ...ten(1, 2)]);
  await historyStep(page, 'undo');
  expect(await ringFiles(page), 'the move and its re-assignment are one History step').toEqual([files[0], files[2], files[1]]);
  expect(await schemeIndexes()).toEqual({ [files[0]]: 0, [files[1]]: 2, [files[2]]: 1 });
  await historyStep(page, 'redo');
  expect(await ringFiles(page)).toEqual([files[0], files[1], files[2]]);
  expect(await drawnByFile()).toEqual({ [files[0]]: '#123456', [files[1]]: cycles.ten[1], [files[2]]: cycles.ten[2] });

  // Under Same color the checkbox has no effect (D-21): a move there leaves the indexes, so 10 colors comes back
  // with those of before the move.
  await pickScheme(page, 'same');
  await outward.nth(2).click();
  await settle(page);
  expect(await ringFiles(page)).toEqual([files[0], files[2], files[1]]);
  expect(await schemeIndexes()).toEqual({ [files[0]]: 0, [files[1]]: 1, [files[2]]: 2 });
  await pickScheme(page, 'ten');
  expect(await swatches(page)).toEqual(['#123456', ...ten(2, 1)]);

  // Use the scheme color drops the ring's own colour, as one History step.
  await page.getByTestId('circular-ring-use-scheme-color').click();
  await settle(page);
  expect(await ownColors(page)).toEqual([null, null, null]);
  expect(await swatches(page)).toEqual(ten(0, 2, 1));
  await expect(page.getByTestId('circular-ring-use-scheme-color')).toHaveCount(0);
  await historyStep(page, 'undo');
  expect(await ownColors(page)).toEqual(['#123456', null, null]);

  // Generate draws the panel's colours (one owner for both).
  await generateAndWaitForResult(page);
  expect(await generatedRingColors(page)).toEqual(await swatches(page));
  const saved = await savedRings(page, testInfo, 'saved-hand-set');
  expect(saved).toMatchObject({ colorScheme: 'ten', reassignOnReorder: true });
  expect(saved.series.map((entry) => [entry.fileName, entry.color, entry.schemeIndex]))
    .toEqual([[files[0], '#123456', 0], [files[2], null, 2], [files[1], null, 1]]);
});
