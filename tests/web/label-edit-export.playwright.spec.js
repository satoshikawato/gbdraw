// The SVG export of a Result whose labels the editor bound and edited.
// OV-3848: the label binder's marks (data-label-editable, data-label-key,
// data-label-source-text) stay on the page; an export never carries them, so
// its bytes do not depend on whether a label was clicked, and a renamed
// label's old text does not leave the page. OV-2438: a label text edit and a
// label hide show their result on the mounted Result (the projection that
// reads each label's parts from one index of the SVG).
const { expect } = require('@playwright/test');
const fs = require('node:fs/promises');
const { warmAppTest: test } = require('./helpers/app-lifecycle.cjs');
const { open, popupEdit } = require('./helpers/live-generate-parity-steps.cjs');

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
