// The export and the stored Result of a Result whose labels the editor bound
// and edited. OV-3848: an export never carries the label binder's marks
// (data-label-editable, data-label-key, data-label-source-text), and a
// Result's text (Save, History) keeps only the source text of a renamed label,
// so neither depends on whether a label was clicked, and a renamed label's old
// text does not leave the page. Reset shows the source text after a remount
// and a Load. OV-2438: a label text edit and a label hide show their result on
// the mounted Result (the projection that reads each label's parts from one
// index of the SVG).
const { expect } = require('@playwright/test');
const fs = require('node:fs/promises');
const { evaluateWithRetainedPromise, warmAppTest: test } = require('./helpers/app-lifecycle.cjs');
const { generate, legendRowColor, open, popupEdit } = require('./helpers/live-generate-parity-steps.cjs');
const { settleLive, showResult } = require('./helpers/live-generate-parity.cjs');
const { download, load } = require('./helpers/mode-transition.cjs');
const { readSessionFile, readSessionFileJson, sessionFileBytesOf } = require('./helpers/session-file.cjs');

test.describe.configure({ retries: 0 });

const BINDER_MARK = /data-label-(?:editable|key|source-text)=/;

const exportSvg = async (page, testInfo, name) => {
  const pending = page.waitForEvent('download');
  await page.evaluate(async () => { await window.__GBDRAW_APP__.downloadSVG(); });
  const file = testInfo.outputPath(`${name}.svg`);
  await (await pending).saveAs(file);
  return fs.readFile(file, 'utf8');
};

// The mounted Result's labels of a feature (by locus_tag): each part's text and display.
const labelParts = (page, locusTag) => page.evaluate(async (tag) => {
  const { state } = await import('./js/state.js');
  const feature = window.__GBDRAW_APP__.filteredFeatures.find((item) => item.locus_tag === tag);
  const svg = state.svgContainer.value.querySelector('svg');
  return [...svg.querySelectorAll(`[data-label-feature-id="${CSS.escape(feature.svg_id)}"]`)]
    .map((part) => ({ tag: part.localName, text: part.localName === 'text' ? part.textContent : null, display: part.getAttribute('display') }));
}, locusTag);

const boundLabelCount = (page) => page.evaluate(async () => {
  const { state } = await import('./js/state.js');
  return state.svgContainer.value.querySelectorAll('svg text[data-label-editable="true"]').length;
});

test('an SVG export is the same before and after a popup binds the labels', async ({ page }, testInfo) => {
  test.setTimeout(240_000);
  await open(page, { mode: 'linear', results: 'single', reflow: 'off' });
  expect(await boundLabelCount(page)).toBe(0);
  const unbound = await exportSvg(page, testInfo, 'unbound');
  expect(unbound).not.toMatch(BINDER_MARK);

  // A popup opened and closed binds the labels of the mounted Result.
  await popupEdit(page, 'FL2');
  expect(await boundLabelCount(page)).toBeGreaterThan(0);
  expect(await exportSvg(page, testInfo, 'bound')).toBe(unbound);
});

test('a label text edit and a label hide show on the Result, and the export holds neither the old text nor the marks', async ({ page }, testInfo) => {
  test.setTimeout(240_000);
  await open(page, { mode: 'linear', results: 'single', reflow: 'off' });
  const before = await labelParts(page, 'FL1');
  expect(before.find((part) => part.tag === 'text')?.text).toBe('alpha protein');

  await popupEdit(page, 'FL1', { labelText: 'renamed label' });
  await popupEdit(page, 'FL2', { labelVisibility: 'off' });

  const renamed = await labelParts(page, 'FL1');
  expect(renamed.filter((part) => part.tag === 'text').map((part) => part.text)).toEqual(['renamed label']);
  const hidden = await labelParts(page, 'FL2');
  expect(hidden.length).toBeGreaterThan(0);
  expect(hidden.map((part) => part.display)).toEqual(hidden.map(() => 'none'));
  await fs.writeFile(testInfo.outputPath('mounted-after-edits.svg'), await page.evaluate(async () => {
    const { state } = await import('./js/state.js');
    return state.svgContainer.value.querySelector('svg').outerHTML;
  }));

  const exported = await exportSvg(page, testInfo, 'edited');
  expect(exported).not.toMatch(BINDER_MARK);
  expect(exported).toContain('renamed label');
  expect(exported).not.toContain('alpha protein');
});

// The text of each feature's label (by locus_tag) on the displayed Result.
const labelTexts = (page, tags) => page.evaluate(async (wanted) => {
  const { state } = await import('./js/state.js');
  const app = window.__GBDRAW_APP__;
  const svg = state.svgContainer.value.querySelector('svg');
  return Object.fromEntries(wanted.map((tag) => {
    const feature = app.filteredFeatures.find((item) => item.locus_tag === tag);
    const text = feature && svg.querySelector(`text[data-label-feature-id="${CSS.escape(feature.svg_id)}"]`);
    return [tag, text ? text.textContent : null];
  }));
}, tags);

test('renamed labels show their original text after Reset label edits, also after a Result switch and a Save and Load', async ({ page, browser }, testInfo) => {
  test.setTimeout(240_000);
  await open(page, { mode: 'circular', results: 'batch', reflow: 'off' });
  const tags = ['TESTA_0001', 'TESTA_0002', 'TESTA_0003'];
  const original = await labelTexts(page, tags);
  expect(original).toEqual({ TESTA_0001: 'duplicate protein', TESTA_0002: 'duplicate protein', TESTA_0003: 'origin spanning protein' });

  // One label renamed in its popup, and every label of one text renamed (the label list's "all labels" choice).
  await popupEdit(page, 'TESTA_0003', { labelText: 'renamed single' });
  await evaluateWithRetainedPromise(page, async () => {
    const app = window.__GBDRAW_APP__;
    const feature = app.filteredFeatures.find((item) => item.locus_tag === 'TESTA_0001');
    const entry = app.editableLabels.find((item) => item.featureId === feature.svg_id);
    // The app's actions wait for the edit queue; each is awaited.
    await app.requestLabelTextChangeByKey(entry.key, 'renamed bulk');
    await app.handleLabelTextScopeChoice('all');
  });
  await settleLive(page);
  const renamed = { TESTA_0001: 'renamed bulk', TESTA_0002: 'renamed bulk', TESTA_0003: 'renamed single' };
  expect(await labelTexts(page, tags)).toEqual(renamed);

  // The Result remounts from its stored text, which holds no binder mark but
  // the source text of each renamed label (the record of what Python drew).
  await showResult(page, 1);
  await showResult(page, 0);
  expect(await labelTexts(page, tags)).toEqual(renamed);
  const stored = await page.evaluate(() => window.__GBDRAW_APP__.currentResults().map((result) => result.content).join(''));
  expect(stored).not.toMatch(/data-label-(editable|key)=/);
  expect([...stored.matchAll(/data-label-source-text="([^"]*)"/g)].map((match) => match[1]).sort())
    .toEqual(['duplicate protein', 'duplicate protein', 'origin spanning protein']);

  const saved = testInfo.outputPath('renamed-labels.gbdraw');
  await download(page, 'Save Session', saved, { title: '' });
  const loaded = await load(browser, saved);
  try {
    for (const shown of [page, loaded]) {
      expect(await labelTexts(shown, tags)).toEqual(renamed);
      await evaluateWithRetainedPromise(shown, async () => { await window.__GBDRAW_APP__.resetAllLabelTextOverrides(); });
      await settleLive(shown);
      expect(await labelTexts(shown, tags)).toEqual(original);
    }
  } finally {
    await loaded.context().close();
  }
});

// The Results of a saved Session: the SVG text of each drawing's Results.
const savedResultSvgs = async (file) => (await readSessionFile(file)).drawings
  .flatMap((drawing) => (drawing.results || []).map((result) => String(result.content)));
const sourceTextRecords = (svgs) => svgs.flatMap((svg) => [...svg.matchAll(/data-label-source-text="([^"]*)"/g)].map((match) => match[1]));

// The color of the first Legend row, set to `color` (the Legend editor's control; no popup, so no label binding).
const colorCdsRow = async (page, color) => {
  const row = await page.evaluate(() => {
    const [entry] = window.__GBDRAW_APP__.legendEntries;
    return entry ? { caption: entry.caption, color: window.__GBDRAW_APP__.legendEntryColor(entry) } : null;
  });
  expect(row?.color, 'the first Legend row').toBeTruthy();
  await legendRowColor(page, row.caption, color);
  return row.color;
};

// OV-3848: the record of a renamed label's source text is written when the
// Result enters the page (Generate, Load), from the feature row, so a saved
// Session does not depend on whether a label was clicked or the displayed
// Result was written again.
test('a Session saved after a Generate of a label edit is the same with and without a label click', async ({ page, browser }, testInfo) => {
  test.setTimeout(240_000);
  await open(page, { mode: 'linear', results: 'single', reflow: 'off' });
  await popupEdit(page, 'FL1', { labelText: 'renamed label' });
  await generate(page);

  // Saved at once, and after the displayed Result was written again (a color changed and back).
  const atGenerate = testInfo.outputPath('at-generate.gbdraw');
  await download(page, 'Save Session', atGenerate, { title: '' });
  const original = await colorCdsRow(page, '#7b2cbf');
  await colorCdsRow(page, original);
  const rewritten = testInfo.outputPath('rewritten.gbdraw');
  await download(page, 'Save Session', rewritten, { title: '' });
  expect(sourceTextRecords(await savedResultSvgs(atGenerate))).toEqual(['alpha protein']);
  expect(sourceTextRecords(await savedResultSvgs(rewritten))).toEqual(['alpha protein']);

  // Loaded twice: one page clicks a label (a popup opened and closed binds the
  // labels), the other does not; both change one color and save.
  const saves = [];
  for (const click of [true, false]) {
    const loaded = await load(browser, atGenerate);
    try {
      if (click) await popupEdit(loaded, 'FL2');
      await colorCdsRow(loaded, '#2a9d8f');
      const file = testInfo.outputPath(`loaded-${click ? 'clicked' : 'not-clicked'}.gbdraw`);
      await download(loaded, 'Save Session', file, { title: '' });
      saves.push(await savedResultSvgs(file));
    } finally {
      await loaded.context().close();
    }
  }
  expect(sourceTextRecords(saves[0])).toEqual(['alpha protein']);
  expect(saves[0]).toEqual(saves[1]);
});

// A Result saved before its labels named their features (Sessions 40 to 43) has no
// `data-label-feature-id`, so Load cannot tie its labels to their rows; the label
// projection ties them when the Result is shown and writes the record then. It
// writes it into the Result, so a Save right after Load holds the record, as a
// Save after the shown Result was written again does.
test('a loaded Result whose labels name no feature records the renamed label source text without a rewrite', async ({ page, browser }, testInfo) => {
  test.setTimeout(240_000);
  await open(page, { mode: 'linear', results: 'single', reflow: 'off' });
  await popupEdit(page, 'FL1', { labelText: 'renamed label' });
  await generate(page);
  const atGenerate = testInfo.outputPath('at-generate.gbdraw');
  await download(page, 'Save Session', atGenerate, { title: '' });

  // The same Session with Results as a version before label feature IDs saved them: no label names its feature, no record.
  const document = await readSessionFileJson(atGenerate);
  let stripped = 0;
  document.drawings.forEach((drawing) => (drawing.results || []).forEach((result) => {
    const content = String(result.content);
    result.content = content.replace(/ data-label-(?:feature-id|source-text)="[^"]*"/g, '');
    stripped += content.length - result.content.length;
  }));
  expect(stripped, 'label feature IDs removed').toBeGreaterThan(0);
  const unnamed = testInfo.outputPath('unnamed.gbdraw');
  await fs.writeFile(unnamed, await sessionFileBytesOf(document));

  const saves = [];
  for (const rewrite of [false, true]) {
    const loaded = await load(browser, unnamed);
    try {
      if (rewrite) {
        await colorCdsRow(loaded, await colorCdsRow(loaded, '#2a9d8f'));
      } else {
        // Writing the record into the Result is no edit: no History step, no unsaved work.
        expect(await loaded.evaluate(async () => ({
          canUndo: window.__GBDRAW_APP__.canUndoHistory,
          unsaved: await window.__GBDRAW_HISTORY__.hasChangesSinceSavePoint()
        }))).toEqual({ canUndo: false, unsaved: false });
      }
      const file = testInfo.outputPath(`loaded-${rewrite ? 'rewritten' : 'at-once'}.gbdraw`);
      await download(loaded, 'Save Session', file, { title: '' });
      saves.push(await savedResultSvgs(file));
    } finally {
      await loaded.context().close();
    }
  }
  expect(sourceTextRecords(saves[0])).toEqual(['alpha protein']);
  // A color changed and back keeps the record of Python's swatch fill (`data-gbdraw-base-fill`, OV-3951), which
  // this case does not judge.
  const withoutPaintBases = (svgs) => svgs.map((svg) => svg.replace(/ data-gbdraw-base-[a-z-]+="[^"]*"/g, ''));
  expect(withoutPaintBases(saves[0])).toEqual(withoutPaintBases(saves[1]));
});
