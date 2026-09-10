const { test, expect } = require('@playwright/test');
const { seeds, load, generate, switchMode, popup, closeEditor, download, snapshot } = require('./helpers/mode-transition.cjs');

const rename = async (page, text) => {
  const target = await popup(page);
  await page.locator('.feature-popup input[placeholder="Edit label text"]').fill(text);
  await page.getByRole('button', { name: 'Apply Label', exact: true }).click();
  await closeEditor(page);
  return target;
};

const labelIntent = state => ({ labels: state.labels, bulk: state.bulkLabels,
  sources: state.labelSources, visibility: state.visibility });

for (const mode of ['circular', 'linear']) {
  test(`@pr-smoke ${mode} label intent survives mode Undo/Redo and Save/Load/Generate`, async ({ browser }, testInfo) => {
    test.setTimeout(360000);
    const page = await load(browser, seeds[mode]);
    let fresh;
    const text = 'MODE_RETAINED_LABEL';
    try {
      await generate(page);
      await rename(page, text);
      const edited = await snapshot(page);
      expect(Object.values(edited.labels)).toContain(text);
      expect(edited.mounted).toContain(text);
      const temporary = mode === 'circular' ? 'linear' : 'circular';
      await switchMode(page, temporary);
      expect.soft(labelIntent(await snapshot(page))).toEqual(labelIntent(edited));
      expect((await snapshot(page)).featureCount).toBe(0);
      await switchMode(page, mode);
      const check = async name => {
        const returned = await snapshot(page);
        await testInfo.attach(name, { body: JSON.stringify(returned), contentType: 'application/json' });
        expect.soft(labelIntent(returned)).toEqual(labelIntent(edited));
        expect.soft(returned.mounted).toContain(text);
        expect.soft(returned.result).toContain(text);
      };
      await check('returned');
      for (const action of ['Undo', 'Redo']) {
        for (const expectedMode of [temporary, mode]) {
          await page.getByRole('button', { name: action, exact: true }).click();
          await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.mode)).toBe(expectedMode);
          expect.soft(labelIntent(await snapshot(page))).toEqual(labelIntent(edited));
        }
        await check(action);
      }
      const file = testInfo.outputPath('roundtrip.gbdraw-session.json.gz');
      await download(page, 'Save Session', file);
      fresh = await load(browser, file);
      expect.soft(labelIntent(await snapshot(fresh))).toEqual(labelIntent(edited));
      await generate(fresh);
      expect.soft((await snapshot(fresh)).mounted).toContain(text);
      await generate(page);
      expect.soft((await snapshot(page)).mounted).toContain(text);
      expect(page.externalRequests).toEqual([]);
      expect(fresh.externalRequests).toEqual([]);
    } finally {
      await page.context().close();
      if (fresh) await fresh.context().close();
    }
  });
}

test('label visibility intent survives mode inactivity and regeneration', async ({ browser }, testInfo) => {
  test.setTimeout(240000);
  const page = await load(browser);
  try {
    await generate(page);
    const target = await popup(page);
    await expect(page.locator(`.origin-top text[data-label-feature-id="${target.featureId}"]`)).toBeVisible();
    await page.locator('.feature-popup select').filter({ has: page.locator('option[value="on"]')
      .filter({ hasText: 'On (force show, bypass filters)' }) }).selectOption('off');
    await page.getByRole('button', { name: 'Apply Label', exact: true }).click();
    await closeEditor(page);
    const edited = await snapshot(page);
    expect(Object.values(edited.visibility)).toEqual(['off']);
    const visible = () => page.locator(`.origin-top text[data-label-feature-id="${target.featureId}"]`).isVisible();
    expect(await visible()).toBe(false);
    await switchMode(page, 'linear');
    await switchMode(page, 'circular');
    expect.soft((await snapshot(page)).visibility).toEqual(edited.visibility);
    expect.soft(await visible()).toBe(false);
    await generate(page);
    expect.soft((await snapshot(page)).visibility).toEqual(edited.visibility);
    expect.soft(await visible()).toBe(false);
    await testInfo.attach('visibility-state', { body: JSON.stringify(await snapshot(page)), contentType: 'application/json' });
  } finally { await page.context().close(); }
});

test('an applicable bulk label override remains dormant across mode navigation', async ({ browser }, testInfo) => {
  test.setTimeout(240000);
  const page = await load(browser);
  try {
    await generate(page);
    await popup(page);
    await closeEditor(page);
    // Seed the existing bulk representation through the Session feature-state owner.
    // Normal feature scope UI resolves matching labels into feature-specific keys.
    await page.evaluate(async () => {
      const { buildFeatureStateData, applyFeatureStateData } = await import('./js/services/config.js');
      applyFeatureStateData({ ...buildFeatureStateData(),
        labelTextBulkOverrides: { 'tRNA-Phe': 'BULK_RETAINED_LABEL' } });
      window.__GBDRAW_APP__.syncLabelEditor();
    });
    const before = await snapshot(page);
    expect(before.bulkLabels).toEqual({ 'tRNA-Phe': 'BULK_RETAINED_LABEL' });
    expect(before.mounted).toContain('BULK_RETAINED_LABEL');
    await switchMode(page, 'linear');
    expect.soft((await snapshot(page)).bulkLabels).toEqual(before.bulkLabels);
    await switchMode(page, 'circular');
    expect.soft((await snapshot(page)).bulkLabels).toEqual(before.bulkLabels);
    expect.soft((await snapshot(page)).mounted).toContain('BULK_RETAINED_LABEL');
    await generate(page);
    expect.soft((await snapshot(page)).mounted).toContain('BULK_RETAINED_LABEL');
    await testInfo.attach('bulk-state', { body: JSON.stringify(await snapshot(page)), contentType: 'application/json' });
  } finally { await page.context().close(); }
});

test('incompatible source replacement reconciles old label intent after a mode round trip', async ({ browser }) => {
  test.setTimeout(240000);
  const page = await load(browser);
  try {
    await generate(page);
    await rename(page, 'SOURCE_A_ONLY_LABEL');
    const edited = await snapshot(page);
    await switchMode(page, 'linear');
    await switchMode(page, 'circular');
    expect.soft((await snapshot(page)).labels).toEqual(edited.labels);
    await page.getByLabel('GenBank/DDBJ File', { exact: true }).setInputFiles('examples/MellatMJNV.gb');
    await expect.poll(() => page.evaluate(async () =>
      (await import('./js/state.js')).state.circularRecordDiscovery.status)).toBe('ready');
    await generate(page);
    const replaced = await snapshot(page);
    expect(replaced.labels).toEqual({});
    expect(replaced.mounted).not.toContain('SOURCE_A_ONLY_LABEL');
    expect(replaced.result).not.toContain('SOURCE_A_ONLY_LABEL');
  } finally { await page.context().close(); }
});

test('feature visibility overrides and qualifier rules survive mode navigation', async ({ browser }) => {
  test.setTimeout(240000);
  const page = await load(browser);
  try {
    await generate(page);
    await popup(page);
    await page.getByLabel('Feature visibility', { exact: true }).selectOption('off');
    if (await page.getByRole('heading', { name: 'Feature Visibility Scope', exact: true }).isVisible()) {
      await page.getByText('This feature', { exact: true }).click();
    }
    await closeEditor(page);
    await popup(page, 1);
    await page.getByLabel('Feature visibility', { exact: true }).selectOption('off');
    await page.getByText(/^Exact product:/).click();
    await closeEditor(page);
    const before = await snapshot(page);
    expect(Object.values(before.featureVisibility)).toEqual(['off']);
    expect(before.visibilityRules).toHaveLength(1);
    await switchMode(page, 'linear');
    await switchMode(page, 'circular');
    const after = await snapshot(page);
    expect.soft(after.featureVisibility).toEqual(before.featureVisibility);
    expect.soft(after.visibilityRules).toEqual(before.visibilityRules);
    await generate(page);
    const generated = await snapshot(page);
    expect.soft(generated.featureVisibility).toEqual(before.featureVisibility);
    expect.soft(generated.visibilityRules).toEqual(before.visibilityRules);
  } finally { await page.context().close(); }
});
