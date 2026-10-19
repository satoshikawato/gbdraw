// A click (or Enter) on a Legend row or a record definition in the preview
// opens its popup while Layout edit is off. The popup edits through the Legend
// panel's operations and the sidebar's record settings, one History step each;
// Escape returns focus to the diagram item. With Layout edit on, no popup opens.
const { test, expect } = require('@playwright/test');
const { HMMT_SESSION, loadSessionFile, openFresh, settle } = require('./helpers/audit-browser.cjs');
const { evaluateWithRetainedPromise } = require('./helpers/app-lifecycle.cjs');

test.describe.configure({ retries: 0 });

const LINEAR_SESSION = 'gbdraw/web/gallery/sessions/lambda_basic_linear.gbdraw-session.json';
const undoCount = (page) => page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount());
const popup = (page) => page.locator('.preview-item-popup[role="dialog"]');
const legendRow = (page, key) => page.locator(`.gbdraw-preview-surface g[data-legend-key="${key}"]`).first();

// Clicks the centre of `locator` and returns the clicked point.
const clickCentre = async (locator) => {
  await locator.scrollIntoViewIfNeeded();
  const box = await locator.boundingBox();
  await locator.click({ position: { x: box.width / 2, y: box.height / 2 } });
  return { x: box.x + box.width / 2, y: box.y + box.height / 2 };
};

// The popups open as the feature popup does (OV-371, OV-384): inside the
// window, below the app header, and beside the clicked point.
const expectPlacedBeside = async (page, point) => {
  await expect.poll(async () => {
    const box = await popup(page).boundingBox();
    const limits = await page.evaluate(() => ({
      top: document.querySelector('header.app-header').getBoundingClientRect().bottom,
      bottom: Math.min(window.innerHeight, document.querySelector('[data-app-footer]')?.getBoundingClientRect().top ?? Infinity),
      right: window.innerWidth
    }));
    const covers = point.x >= box.x && point.x <= box.x + box.width && point.y >= box.y && point.y <= box.y + box.height;
    return {
      belowHeader: box.y >= limits.top - 0.5,
      aboveBottom: box.y + box.height <= limits.bottom + 0.5,
      insideWidth: box.x >= -0.5 && box.x + box.width <= limits.right + 0.5,
      coversPoint: covers
    };
  }).toEqual({ belowHeader: true, aboveBottom: true, insideWidth: true, coversPoint: false });
};

test('a Legend row opens its popup; rename, remove, keyboard, and Layout edit', async ({ page }) => {
  test.setTimeout(300_000);
  await openFresh(page);
  await loadSessionFile(page, HMMT_SESSION);

  const legendPoint = await clickCentre(legendRow(page, 'tRNA').locator('path').first());
  const dialog = page.getByRole('dialog', { name: 'Legend item: tRNA', exact: true });
  await expect(dialog).toBeVisible();
  await expectPlacedBeside(page, legendPoint);
  await expect(dialog.locator('[data-preview-legend-color-source]')).toHaveText(/^Sets the tRNA default color/);

  // Rename: one History step, shown in the Legend panel's rows.
  const before = await undoCount(page);
  await dialog.getByLabel('Legend name', { exact: true }).fill('Transfer RNA');
  await dialog.getByLabel('Legend name', { exact: true }).press('Enter');
  await settle(page);
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.legendEntries.map((entry) => entry.caption)))
    .toContain('Transfer RNA');
  expect(await undoCount(page)).toBe(before + 1);
  await expect(page.getByRole('dialog', { name: 'Legend item: Transfer RNA', exact: true })).toBeVisible();

  // Remove from legend: the row moves to the deleted list and the popup closes.
  await popup(page).getByRole('button', { name: 'Remove from legend', exact: true }).click();
  await settle(page);
  await expect(popup(page)).toHaveCount(0);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.deletedLegendEntries.map((entry) => entry.caption)))
    .toContain('Transfer RNA');

  // Keyboard: Enter opens the popup with focus inside; Escape returns focus.
  const key = await page.evaluate(() => window.__GBDRAW_APP__.legendEntries[0].caption);
  const row = legendRow(page, key);
  await row.focus();
  await page.keyboard.press('Enter');
  await expect(popup(page)).toBeVisible();
  expect(await popup(page).evaluate((element) => element.contains(document.activeElement))).toBe(true);
  await page.keyboard.press('Escape');
  await expect(popup(page)).toHaveCount(0);
  expect(await row.evaluate((element) => element === document.activeElement)).toBe(true);

  // Layout edit on: a press drags the legend; no popup opens.
  await page.evaluate(() => { window.__GBDRAW_APP__.layoutRepositionMode = true; });
  await legendRow(page, key).locator('path').first().click();
  await expect(popup(page)).toHaveCount(0);

  // The exported SVG carries none of the preview's popup affordances.
  const exported = await page.evaluate(async () => {
    const { serializeCleanSvg } = await import('./js/services/svg-serialization.js');
    return serializeCleanSvg(window.__GBDRAW_APP__.svgContainer.querySelector('svg'));
  });
  expect(exported).not.toContain('data-gbdraw-preview-item');
});

// OV-1272: a GC content or GC skew row names its plot, not a feature count,
// and says that its Color changes the swatch only.
test('a GC content or GC skew Legend row popup names its plot, not features', async ({ page }) => {
  test.setTimeout(300_000);
  await openFresh(page);
  await loadSessionFile(page, HMMT_SESSION);

  for (const [key, subtitle] of [['GC content', 'Content plot'], ['GC skew (+)', 'Skew plot']]) {
    await clickCentre(legendRow(page, key).locator('path').first());
    const dialog = page.getByRole('dialog', { name: `Legend item: ${key}`, exact: true });
    await expect(dialog).toBeVisible();
    await expect(dialog.locator('[data-preview-item-subtitle]')).toHaveText(subtitle);
    await expect(dialog.locator('[data-preview-legend-color-source]'))
      .toHaveText('Changes this legend swatch only. The plot keeps the color set in Settings.');
    await page.keyboard.press('Escape');
    await expect(popup(page)).toHaveCount(0);
  }
});

// OV-1281 (Owner, Q-09 a): on a row that a feature type's default color draws
// (no Specific color rule, no color of its own on a feature), the popup's Color
// sets that type's Default color as the scope dialog's Legend item choice does
// (OIC-029): the features and the legend change together, as one History step.
test('a default-colored Legend row popup Color sets the type\'s default color', async ({ page }) => {
  test.setTimeout(300_000);
  await openFresh(page);
  await loadSessionFile(page, HMMT_SESSION);
  const state = () => page.evaluate(async () => {
    const { state: s } = await import('./js/state.js');
    const app = window.__GBDRAW_APP__;
    const svg = app.svgContainer.querySelector('svg');
    const fills = app.extractedFeatures.filter((feature) => feature.type === 'tRNA').map((feature) => {
      const element = svg.querySelector(`[id="${CSS.escape(feature.svg_id)}"]`);
      const shape = element?.matches('path, rect, polygon') ? element : element?.querySelector('path, rect, polygon');
      return shape?.getAttribute('fill') || null;
    }).filter(Boolean);
    const swatch = svg.querySelector('#legend g[data-legend-key="tRNA"] path[fill]')?.getAttribute('fill') || null;
    return { fills: [...new Set(fills)], swatch, defaultColor: app.currentColors.tRNA, overrides: { ...s.activeDrawing().legendColorOverrides } };
  });
  const before = await state();
  expect(before.fills.length).toBe(1);

  await clickCentre(legendRow(page, 'tRNA').locator('path').first());
  const dialog = page.getByRole('dialog', { name: 'Legend item: tRNA', exact: true });
  await expect(dialog).toBeVisible();
  await expect(dialog.locator('[data-preview-legend-color-source]')).toHaveText(
    'Sets the tRNA default color (tRNA features without their own color or rule, also hidden ones).');
  const steps = await undoCount(page);
  await dialog.getByLabel('Color of tRNA', { exact: true }).fill('#c026d3');
  await settle(page);
  await expect.poll(state).toEqual({ fills: ['#c026d3'], swatch: '#c026d3', defaultColor: '#c026d3', overrides: before.overrides });
  expect(await undoCount(page), 'one History step').toBe(steps + 1);

  await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.undo());
  await settle(page);
  await expect.poll(state).toEqual(before);
  await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.redo());
  await settle(page);
  await expect.poll(state).toEqual({ fills: ['#c026d3'], swatch: '#c026d3', defaultColor: '#c026d3', overrides: before.overrides });
  await evaluateWithRetainedPromise(page, () => window.__GBDRAW_APP__.runAnalysis());
  await settle(page);
  await expect.poll(state, 'Generate draws the same').toEqual({ fills: ['#c026d3'], swatch: '#c026d3', defaultColor: '#c026d3', overrides: before.overrides });
});

test('a record definition opens its popup; its edits write the record settings', async ({ page }) => {
  test.setTimeout(300_000);
  await openFresh(page);
  await loadSessionFile(page, LINEAR_SESSION);

  const recordPoint = await clickCentre(page.locator('.gbdraw-preview-surface [data-gbdraw-role^="record-definition"] text').first());
  const dialog = popup(page);
  await expect(dialog).toBeVisible();
  await expectPlacedBeside(page, recordPoint);
  await expect(dialog).toHaveAttribute('aria-label', /^Record /);
  await expect(dialog.locator('dl')).toContainText('Accession');
  await expect(dialog.locator('dl')).toContainText('Whole record');

  const before = await undoCount(page);
  await dialog.getByRole('checkbox', { name: 'Reverse complement', exact: true }).check();
  await settle(page);
  expect(await undoCount(page)).toBe(before + 1);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.linearSeqs[0].region_reverse)).toBe(true);
  // Enter commits a text field and keeps focus in it; leaving the field afterwards commits nothing more.
  const name = dialog.getByLabel('Organism / strain', { exact: true });
  await name.fill('Popup organism');
  await name.press('Enter');
  await settle(page);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.linearSeqs[0].definition)).toBe('Popup organism');
  await expect(name).toBeFocused();
  expect(await undoCount(page)).toBe(before + 2);
  await dialog.getByLabel('Subtitle', { exact: true }).fill('Popup subtitle');
  await dialog.getByLabel('Subtitle', { exact: true }).press('Tab');
  await settle(page);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.linearSeqs[0].record_subtitle)).toBe('Popup subtitle');
  expect(await undoCount(page)).toBe(before + 3);
  // Show in settings focuses the record's own name field. This File has one
  // record, so the field is on its File card, not in Record options (OV-508).
  await dialog.getByRole('button', { name: 'Show in settings', exact: true }).click();
  await expect(page.getByLabel('Definition for sequence 1', { exact: true })).toBeFocused();
  await expect(dialog).toBeHidden();
  await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.undo());
  await settle(page);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.linearSeqs[0].record_subtitle || '')).not.toBe('Popup subtitle');
});

// A Circular record selected on its own draws its Record label, so the popup
// edits that field, and the popup leaves the whole center
// definition uncovered, not only the clicked point.
test('a Circular center definition opens beside the definition and edits the Record label', async ({ page }) => {
  test.setTimeout(300_000);
  await openFresh(page);
  await loadSessionFile(page, HMMT_SESSION);
  // A loaded Session's records are not inspected; the Record label fields open once they are.
  await page.locator('[data-circular-inspect]').click();
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.circularSingleRecordPresentationEnabled)).toBe(true);

  const definition = page.locator('.gbdraw-preview-surface [data-gbdraw-role^="record-definition"]').first();
  const recordPoint = await clickCentre(definition.locator('text').first());
  const dialog = popup(page);
  await expect(dialog).toBeVisible();
  await expectPlacedBeside(page, recordPoint);
  await expect.poll(async () => {
    const [a, b] = [await dialog.boundingBox(), await definition.boundingBox()];
    return a.x >= b.x + b.width || a.x + a.width <= b.x || a.y >= b.y + b.height || a.y + a.height <= b.y;
  }).toBe(true);

  const before = await undoCount(page);
  await dialog.getByLabel('Record label', { exact: true }).fill('Popup label');
  await dialog.getByLabel('Record label', { exact: true }).press('Tab');
  await settle(page);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.form.circular_record_label)).toBe('Popup label');
  expect(await undoCount(page)).toBe(before + 1);
  const subtitle = dialog.getByLabel('Subtitle', { exact: true });
  await subtitle.fill('Popup subtitle');
  await subtitle.press('Enter');
  await settle(page);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.form.circular_record_subtitle)).toBe('Popup subtitle');
  await expect(subtitle).toBeFocused();
});

// A rename onto another row's name asks first; the popup keeps its row until
// the choice, stays on it after Cancel, and follows it after Keep both.
// OV-431: the Enter that commits the name does not also press the dialog's focused choice.
test('a Legend popup rename onto another name follows the row through the choice', async ({ page }) => {
  test.setTimeout(300_000);
  await openFresh(page);
  await loadSessionFile(page, HMMT_SESSION);
  const target = await page.evaluate(() => window.__GBDRAW_APP__.legendEntries
    .map((entry) => entry.caption).find((caption) => caption !== 'tRNA'));
  await clickCentre(legendRow(page, 'tRNA').locator('path').first());
  const tRNA = page.getByRole('dialog', { name: 'Legend item: tRNA', exact: true });
  const choice = page.getByRole('dialog').filter({ has: page.locator('#legend-rename-title') });

  // The pointer still rests on the row: the "Turn on Layout edit" hover hint
  // does not show over the popup the click opened.
  await expect(tRNA).toBeVisible();
  await expect(page.locator('.preview-layout-hint')).toBeHidden();

  // The "Session loaded" notice stays readable above the popup, and its
  // Dismiss closes the notice only.
  const notice = page.locator('[data-app-notice]');
  await expect(notice).toContainText('Session loaded');
  expect((await tRNA.boundingBox()).y).toBeGreaterThanOrEqual((await notice.boundingBox()).y + (await notice.boundingBox()).height - 0.5);
  await notice.getByRole('button', { name: 'Dismiss notice', exact: true }).click();
  await expect(notice).toHaveCount(0);
  await expect(tRNA).toBeVisible();
  expect(await page.evaluate(() => document.activeElement?.classList.contains('result-pane'))).toBe(true);

  const rename = async () => {
    await tRNA.getByLabel('Legend name', { exact: true }).fill(target);
    await tRNA.getByLabel('Legend name', { exact: true }).press('Enter');
    await expect(choice).toBeVisible();
    await expect(tRNA).toBeVisible();
  };

  await rename();
  await choice.getByRole('button', { name: 'Cancel', exact: true }).click();
  await expect(choice).toHaveCount(0);
  await expect(tRNA).toBeVisible();
  // Focus returns to the field the name was typed in, which shows the row's name again.
  await expect(tRNA.getByLabel('Legend name', { exact: true })).toBeFocused();
  await expect(tRNA.getByLabel('Legend name', { exact: true })).toHaveValue('tRNA');

  await rename();
  const keepBoth = choice.getByRole('button', { name: /^Keep both: name it / });
  const suffix = (await keepBoth.innerText()).match(/"(.+)"/)[1];
  await keepBoth.click();
  await settle(page);
  await expect(page.getByRole('dialog', { name: `Legend item: ${suffix}`, exact: true })).toBeVisible();
});

// PD-OI-042 scenario revision 2 (D-01): a caption written in two colors has two
// Legend rows, whose keys carry a color ID. No popup, list, or dialog shows a
// key, and no edit writes one as a name (OV-640).
test('two color rules of one caption: the popups and the deleted list show names, and a color edit writes no row key', async ({ page }) => {
  test.setTimeout(300_000);
  await openFresh(page);
  await loadSessionFile(page, HMMT_SESSION);
  await evaluateWithRetainedPromise(page, async () => {
    const app = window.__GBDRAW_APP__;
    for (const [val, color] of [['NADH', '#112233'], ['(?i)cytochrome', '#445566']]) {
      Object.assign(app.newSpecRule, { feat: 'CDS', qual: 'product', val, color, cap: 'Shared' });
      await app.addSpecificRule();
    }
  });
  await settle(page);
  const keys = ['Shared [#112233]', 'Shared [#445566]'];
  const drawnNames = () => page.evaluate(() => Object.fromEntries(
    [...document.querySelectorAll('.gbdraw-preview-surface g[data-legend-key]')]
      .map((row) => [row.getAttribute('data-legend-key'), row.querySelector('text')?.textContent])
  ));
  const names = await drawnNames();
  expect(keys.map((key) => names[key]).sort()).toEqual(['Shared (1)', 'Shared (2)']);

  // A feature's hover summary names its row as the Legend draws it, never by the key.
  const hoverFeature = (type) => page.evaluate((eventType) => {
    const id = window.__GBDRAW_APP__.extractedFeatures.find((feature) => /NADH/.test(feature.product)).svg_id;
    const element = [...document.querySelectorAll('.gbdraw-preview-surface :is(path, polygon, rect)')]
      .find((candidate) => candidate.getAttribute('data-gbdraw-feature-id') === id || candidate.id === id);
    const box = element.getBoundingClientRect();
    element.dispatchEvent(new MouseEvent(eventType, { bubbles: true, clientX: box.x + box.width / 2, clientY: box.y + box.height / 2 }));
  }, type);
  await hoverFeature('mouseover');
  const hoverSummary = page.locator('.feature-hover-summary');
  await expect(hoverSummary).toBeVisible();
  await expect.soft(hoverSummary.locator('.feature-hover-summary-row').filter({ hasText: /^\s*Legend/ })
    .locator('.feature-hover-summary-value')).toHaveText(names[keys[0]]);
  await hoverFeature('mouseout');

  // The Legend item popup is titled and filled with the drawn name.
  await clickCentre(legendRow(page, keys[0]).locator('path').first());
  const legendDialog = page.getByRole('dialog', { name: `Legend item: ${names[keys[0]]}`, exact: true });
  await expect(legendDialog).toBeVisible();
  await expect(legendDialog.getByLabel('Legend name', { exact: true })).toHaveValue(names[keys[0]]);
  await page.keyboard.press('Escape');
  await expect(popup(page)).toHaveCount(0);

  // The feature popup shows the drawn name of the feature's row.
  await evaluateWithRetainedPromise(page, async () => {
    const app = window.__GBDRAW_APP__;
    await app.openFeatureEditorFromList(app.extractedFeatures.find((feature) => /NADH/.test(feature.product)), null);
    await window.Vue.nextTick();
  });
  const nameField = page.getByLabel('Feature legend name', { exact: true });
  await expect(nameField).toHaveValue(names[keys[0]]);
  // D-11 (OV-587): the field suggests no group; the list under it is the only place that offers one.
  expect.soft(await nameField.getAttribute('list')).toBeNull();
  expect.soft(await page.locator('#specific-rule-legend-name-options').count()).toBe(0);
  // The list shows each group as the Legend shows its row: by the drawn name, never by the key.
  const offered = () => page.getByLabel('Legend names of color rules', { exact: true }).locator('option').allTextContents()
    .then((texts) => texts.map((text) => text.trim()).slice(1).filter((text) => text.startsWith('Shared')).sort());
  expect.soft(await offered()).toEqual(['Shared (1)', 'Shared (2)']);

  // OV-640: a color for this feature only names its row by the feature, never by the key of the row it leaves.
  await evaluateWithRetainedPromise(page, async () => {
    const app = window.__GBDRAW_APP__;
    await app.updateClickedFeatureColor('#c83366');
    if (!app.featureStyleScopeDialog.show) throw new Error('no Color Change Scope dialog');
  });
  await expect.soft(page.getByRole('button', { name: `Legend item "${names[keys[0]]}"`, exact: false })).toBeVisible();
  await evaluateWithRetainedPromise(page, () => window.__GBDRAW_APP__.handleFeatureStyleScopeChoice('single'));
  await settle(page);
  const written = await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    return {
      captions: app.manualSpecificRules.map((rule) => rule.cap),
      edits: Object.values(app.featureColorOverrides).map((edit) => edit?.caption).filter(Boolean)
    };
  });
  expect.soft(written.captions.filter((caption) => caption.includes(' [#'))).toEqual([]);
  expect.soft(written.edits.filter((caption) => caption.includes(' [#'))).toEqual([]);
  expect.soft(written.captions.filter((caption) => caption === 'Shared')).toHaveLength(2);
  expect.soft(written.captions).toHaveLength(3);
  expect.soft(Object.values(await drawnNames()).filter((name) => name.includes(' [#'))).toEqual([]);
  await page.evaluate(() => { window.__GBDRAW_APP__.clickedFeature = null; });

  // A deleted row the Result no longer holds is listed by its name (D-15), never by its key.
  await clickCentre(legendRow(page, keys[1]).locator('path').first());
  await popup(page).getByRole('button', { name: 'Remove from legend', exact: true }).click();
  await settle(page);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.deletedLegendEntries.map((entry) => entry.caption))).toContain(keys[1]);
  if (!await page.evaluate(() => window.__GBDRAW_APP__.showRightDrawer)) await page.locator('.drawer-toggle').click();
  await page.locator('.right-drawer').getByRole('button', { name: 'Legend' }).click();
  await expect.soft(page.locator('.right-drawer').getByRole('button', { name: 'Restore Shared', exact: true })).toBeVisible();
  await expect.soft(page.locator('.right-drawer').getByText(' [#')).toHaveCount(0);

  // The Legend numbers only the rows it draws: one row is left, drawn without a number, and the group
  // of the deleted row draws none, so the list shows that group by its name as written. The two
  // groups read the same here; how to tell them apart is an open Owner question (P-4).
  expect.soft(await page.evaluate(() => [...document.querySelectorAll('.gbdraw-preview-surface g[data-legend-key]')]
    .filter((row) => row.getAttribute('display') !== 'none').map((row) => row.querySelector('text')?.textContent)
    .filter((name) => name.startsWith('Shared')))).toEqual(['Shared']);
  await evaluateWithRetainedPromise(page, async () => {
    const app = window.__GBDRAW_APP__;
    await app.openFeatureEditorFromList(app.extractedFeatures.find((feature) => /cytochrome/i.test(feature.product)), null);
    await window.Vue.nextTick();
  });
  expect.soft(await offered()).toEqual(['Shared', 'Shared']);
});
