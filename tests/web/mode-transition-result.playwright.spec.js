const { test, expect } = require('@playwright/test');
const { seeds, load, generate, switchMode, popup, closeEditor, download, snapshot, semantics } = require('./helpers/mode-transition.cjs');

const color = async (page, value, index = 0) => {
  const target = await popup(page, index);
  await page.getByLabel('Feature fill color', { exact: true }).first().evaluate((element, value) => {
    element.value = value;
    element.dispatchEvent(new Event('change', { bubbles: true }));
  }, value);
  await page.getByText('This feature only', { exact: true }).click();
  await closeEditor(page);
  return target;
};

const agree = async (page, expected, testInfo, name) => {
  const state = await snapshot(page);
  const exported = (await download(page, 'SVG', testInfo.outputPath(`${name}.svg`))).toString();
  await testInfo.attach(`${name}-state`, { body: JSON.stringify(state), contentType: 'application/json' });
  expect.soft(state.colors).toEqual(expected.colors);
  expect.soft(state.rules).toEqual(expected.rules);
  expect.soft(await semantics(page, state.result)).toEqual(await semantics(page, expected.result));
  expect.soft(await semantics(page, state.mounted)).toEqual(await semantics(page, expected.result));
  expect.soft(await semantics(page, exported)).toEqual(await semantics(page, expected.result));
  expect.soft(state.markedMounted).toBe(true);
  return state;
};

for (const mode of ['circular', 'linear']) {
  test(`@pr-smoke ${mode} color Result agrees with a fresh mode root and export before regeneration`, async ({ browser }, testInfo) => {
    test.setTimeout(240000);
    const page = await load(browser, seeds[mode]);
    try {
      await generate(page);
      const original = await snapshot(page);
      await color(page, '#c83366');
      const edited = await snapshot(page);
      expect(edited.result).toContain('#c83366');
      expect(edited.mounted).toContain('#c83366');
      expect(edited.resultIdentity).toBe(original.resultIdentity);
      await page.evaluate(() => { window.__MODE_EDITED_ROOT__ = document.querySelector('.origin-top svg'); });
      await switchMode(page, mode === 'circular' ? 'linear' : 'circular');
      if (mode === 'circular') {
        const scale = page.getByLabel('Show Coordinate Scale (Linear)', { exact: true });
        for (const details of await scale.locator('xpath=ancestor::details').all()) {
          if (await details.getAttribute('open') === null) await details.locator(':scope > summary').click();
        }
        await scale.uncheck();
        expect((await snapshot(page)).request).toEqual(edited.request);
      }
      await agree(page, edited, testInfo, 'inactive');
      await switchMode(page, mode);
      const returned = await agree(page, edited, testInfo, 'returned');
      expect(returned.sameRoot).toBe(false);
      expect(returned.resultIdentity).toBe(edited.resultIdentity);
      expect(returned.generation).toBe(edited.generation);
      expect.soft(returned.mountEvents.at(-1).rootGeneration).toBeGreaterThan(edited.mountEvents.at(-1).rootGeneration);
      expect.soft(returned.mountEvents.at(-1).bindSequence).toBeGreaterThan(edited.mountEvents.at(-1).bindSequence);
      expect(page.externalRequests).toEqual([]);
    } finally { await page.context().close(); }
  });
}

test('incremental color edits survive mode history, dirty draft, and fresh Save/Load/export', async ({ browser }, testInfo) => {
  test.setTimeout(360000);
  const page = await load(browser);
  let fresh;
  try {
    await generate(page);
    await color(page, '#c83366');
    await color(page, '#267b91', 1);
    const edited = await snapshot(page);
    expect(edited.result).toContain('#c83366');
    expect(edited.result).toContain('#267b91');
    const savedBefore = testInfo.outputPath('before.gbdraw-session.json.gz');
    await download(page, 'Save Session', savedBefore);
    await page.getByRole('button', { name: 'Undo', exact: true }).click();
    const undone = await snapshot(page);
    expect(undone.colors).not.toEqual(edited.colors);
    await switchMode(page, 'linear');
    await switchMode(page, 'circular');
    await agree(page, undone, testInfo, 'after-edit-undo');
    for (const mode of ['linear', 'circular']) {
      await page.getByRole('button', { name: 'Undo', exact: true }).click();
      await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.mode)).toBe(mode);
    }
    await agree(page, undone, testInfo, 'after-mode-undo');
    // Reload the two-edit checkpoint; Generate must not be needed to restore it.
    fresh = await load(browser, savedBefore);
    await agree(fresh, edited, testInfo, 'saved-before-loaded');
    await switchMode(fresh, 'linear');
    await switchMode(fresh, 'circular');
    await agree(fresh, edited, testInfo, 'two-edits-returned');
    const savedAfter = testInfo.outputPath('after.gbdraw-session.json.gz');
    await download(fresh, 'Save Session', savedAfter);
    await fresh.context().close();
    fresh = await load(browser, savedAfter);
    await agree(fresh, edited, testInfo, 'saved-after-loaded');
    await generate(fresh);
    await agree(fresh, edited, testInfo, 'regenerated');
    expect(fresh.externalRequests).toEqual([]);
  } finally {
    await page.context().close();
    if (fresh) await fresh.context().close();
  }
});

test('unedited Result and same-mode navigation keep feature presentation unchanged', async ({ browser }, testInfo) => {
  test.setTimeout(240000);
  for (const mode of ['circular', 'linear']) {
    const page = await load(browser, seeds[mode]);
    try {
      await generate(page);
      const original = await snapshot(page);
      await page.evaluate(() => { window.__MODE_EDITED_ROOT__ = document.querySelector('.origin-top svg'); });
      await switchMode(page, mode);
      expect((await snapshot(page)).sameRoot).toBe(true);
      await switchMode(page, mode === 'circular' ? 'linear' : 'circular');
      await switchMode(page, mode);
      const returned = await agree(page, original, testInfo, `unedited-${mode}`);
      expect(returned.sameRoot).toBe(false);
      expect(returned.resultIdentity).toBe(original.resultIdentity);
      expect(returned.generation).toBe(original.generation);
    } finally { await page.context().close(); }
  }
});
