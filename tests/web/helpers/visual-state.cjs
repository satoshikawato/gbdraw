const { expect } = require('@playwright/test');
const fs = require('node:fs/promises');
const { load: loadSeed, popup } = require('./mode-transition.cjs');

const semantics = (page, content) => page.evaluate(async content => {
  const { svgVisualSemantics } = await import('/tests/web/helpers/svg-visual-semantics.mjs');
  return svgVisualSemantics(content);
}, content);

const capture = async (page, testInfo, name) => {
  // No Save, Generate, flush or editor close. Capture action completion before
  // invoking the public export handler, then capture its adjacent state.
  const pending = page.waitForEvent('download');
  const state = await page.evaluate(async () => {
    const { state: s } = await import('./js/state.js');
    const cap = () => ({ selected: s.results.value[s.selectedResultIndex.value]?.content || '',
      mounted: s.svgContainer.value?.querySelector('svg')?.outerHTML || '',
      catalogIds: s.featureCatalog.value?.items?.find(item => item.resultIndex === s.selectedResultIndex.value)
        ?.features?.map(feature => feature.svgId) || [] });
    const completed = cap();
    await window.__GBDRAW_APP__.downloadSVG();
    return { completed, boundary: cap() };
  });
  const file = testInfo.outputPath(`${name}.svg`);
  await (await pending).saveAs(file);
  const exported = await semantics(page, await fs.readFile(file, 'utf8'));
  for (const phase of Object.values(state)) {
    phase.selected = await semantics(page, phase.selected);
    phase.mounted = await semantics(page, phase.mounted);
  }
  const observation = { ...state, exported };
  await fs.writeFile(testInfo.outputPath(`${name}.json`), JSON.stringify(observation));
  return observation;
};
const assertCoherent = async (observation, message, soft = false) => {
  const { compareVisualSemantics } = await import('./svg-visual-semantics.mjs');
  const check = soft ? expect.soft : expect;
  for (const [phase, state] of Object.entries(observation).filter(([key]) => key !== 'exported')) {
    const options = { catalogIds: state.catalogIds };
    const differences = compareVisualSemantics(state.selected, state.mounted, options);
    check(differences.slice(0, 5), `${message}: ${phase} selected/mounted`).toEqual([]);
    if (phase === 'boundary') {
      check(compareVisualSemantics(state.selected, observation.exported, options).slice(0, 5),
        `${message}: download boundary selected/export`).toEqual([]);
      check(compareVisualSemantics(state.mounted, observation.exported, options).slice(0, 5),
        `${message}: download boundary mounted/export`).toEqual([]);
    }
  }
};
const assertNonTargetsPreserved = async (before, after, targetIds, message, soft = false) => {
  const { nonTargetFeatures, compareVisualSemantics } = await import('./svg-visual-semantics.mjs');
  (soft ? expect.soft : expect)(compareVisualSemantics(nonTargetFeatures(before, targetIds),
    nonTargetFeatures(after, targetIds)).slice(0, 5), message).toEqual([]);
};
const load = async (...args) => { const page = await loadSeed(...args); page.setDefaultTimeout(20000); return page; };

const settle = page => page.waitForFunction(() => !window.__GBDRAW_HISTORY__.restoring.value
  && !window.__GBDRAW_HISTORY__.capturing.value && !window.__GBDRAW_APP__.processing && !window.__GBDRAW_APP__.featureStyleScopeDialog.show);
const agree = async (page, info, name) => {
  await settle(page);
  const observation = await capture(page, info, name);
  await assertCoherent(observation, name, true);
  return observation;
};
const reveal = async locator => {
  for (const details of await locator.locator('xpath=ancestor::details').all()) {
    if (await details.getAttribute('open') === null) await details.locator(':scope > summary').click();
  }
  return locator;
};
const check = async (page, name, enabled) => {
  await (await reveal(page.getByLabel(name, { exact: true }))).setChecked(enabled);
  await settle(page);
};
const color = async page => {
  const target = await popup(page);
  await page.getByLabel('Feature fill color', { exact: true }).first().evaluate(element => {
    element.value = '#c83366'; element.dispatchEvent(new Event('change', { bubbles: true }));
  });
  await page.getByText('This feature only', { exact: true }).click();
  await settle(page);
  return target;
};
const label = async (page, text) => {
  const input = page.locator('.feature-popup input[placeholder="Edit label text"]');
  await input.fill(text);
  await page.getByRole('button', { name: 'Apply Label', exact: true }).click();
  if (await page.getByRole('heading', { name: 'Enable Labels', exact: true }).isVisible()) {
    await page.getByRole('button', { name: /Show all labels/ }).click();
  }
  await settle(page);
  expect(Object.values(await page.evaluate(async () => (await import('./js/state.js')).state.labelTextFeatureOverrides))).toContain(text);
};

module.exports = { semantics, capture, assertCoherent, assertNonTargetsPreserved, load, settle, agree, reveal, check, color, label };
