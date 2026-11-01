// The two displays of the settings pane. One button turns the stacked
// settings into an icon rail (a 60 px column of icon and short name) with at
// most one section docked beside it (320 px), and the first button of the
// rail turns it back. Every section and every control exists once and stays
// mounted in both displays. The browser remembers the display; the Session,
// History, the settings, and the Result do not know it.
const { test, expect } = require('@playwright/test');
const { readFileSync } = require('node:fs');
const { join } = require('node:path');
const { gunzipSync } = require('node:zlib');
const { openApp } = require('./helpers/app-lifecycle.cjs');
const { showDrawing } = require('./helpers/drawing-bar.cjs');

// A click on a button that is not there fails at once instead of at the test's own limit.
test.use({ actionTimeout: 20_000 });

const STORAGE_KEY = 'gbdraw.settings-pane-display';
const LINEAR_SESSION = 'BGC0000708-BGC0000713.gbdraw-session.json';

// Short name and section title of every rail entry, per mode, in rail order.
const RAIL = {
  circular: [
    ['Input', 'Input Genomes'], ['Regions', 'Region Annotations'], ['Compare', 'Pairwise Comparisons'],
    ['Layout', 'Layout'], ['Font', 'Font'], ['Titles', 'Titles & Record Labels'], ['Legend', 'Legend'], ['Colors', 'Colors'],
    ['Features', 'Features'], ['Labels', 'Labels'], ['Axis', 'Axis & Scale'],
    ['GC / skew', 'Dinucleotide content/skew'], ['About', 'About & Citation']
  ],
  linear: [
    ['Input', 'Input Genomes'], ['Compare', 'Comparison'], ['Advanced', 'Advanced comparison and layout'],
    ['Regions', 'Region Annotations'], ['Layout', 'Layout'], ['Font', 'Font'], ['Titles', 'Titles & Record Labels'],
    ['Legend', 'Legend'], ['Colors', 'Colors'], ['Features', 'Features'], ['Labels', 'Labels'],
    ['Axis', 'Axis & Scale'], ['GC / skew', 'Dinucleotide content/skew'], ['About', 'About & Citation']
  ]
};

const index = (page) => page.getByRole('navigation', { name: 'Settings sections' });
const entry = (page, title) => index(page).getByRole('button', { name: `Go to ${title}`, exact: true });
const toRail = (page) => page.getByRole('button', { name: 'Collapse settings to icons', exact: true });
const toStacked = (page) => page.getByRole('button', { name: 'Expand settings', exact: true });
const closeSection = (page) => page.getByRole('button', { name: 'Close section', exact: true });
const generate = (page) => page.getByRole('button', { name: 'Generate Diagram', exact: true });
const setMode = (page, mode) => showDrawing(page, { kind: mode });
const undoCount = (page) => page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount());

const loadSession = async (page, name) => {
  page.on('dialog', (dialog) => dialog.accept());
  await page.locator('input[accept^=".json,"]').setInputFiles(join(process.cwd(), 'gbdraw/web/gallery/sessions', name));
  await page.waitForFunction(() => !window.__GBDRAW_APP__.sessionImportPending
    && window.__GBDRAW_APP__.results.length > 0, null, { timeout: 180000 });
};

// What the pane shows: its display, the docked entry, the widths, the
// sections on screen, the marked entries, and the Tab stops in order.
const paneState = (page) => page.evaluate((storageKey) => {
  const pane = document.querySelector('.settings-pane');
  const area = pane.querySelector('.settings-scroll');
  const visible = (element) => element.checkVisibility({ visibilityProperty: true, opacityProperty: false });
  const box = (element) => element.getBoundingClientRect();
  const stops = Array.from(pane.querySelectorAll('a[href], button, input, select, textarea, summary, [tabindex]'))
    .filter((element) => visible(element) && !element.disabled && element.tabIndex >= 0);
  const place = (element) => {
    if (element.closest('.generate-bar')) return 'generate';
    if (element.closest('nav[aria-label="Settings sections"]')) return 'entry';
    if (element.closest('.settings-index')) return 'switch';
    if (element.closest('.settings-scroll')) return 'section';
    return element.getAttribute('aria-label') || element.outerHTML.slice(0, 60);
  };
  const bar = pane.querySelector('.generate-bar');
  let stored = null;
  try { stored = localStorage.getItem(storageKey); } catch { stored = 'unreadable'; }
  return {
    display: pane.dataset.settingsDisplay,
    docked: pane.dataset.settingsDocked ?? null,
    paneWidth: Math.round(box(pane).width),
    resultWidth: Math.round(box(document.querySelector('.result-pane')).width),
    areaShown: visible(area),
    areaWidth: visible(area) ? Math.round(box(area).width) : 0,
    shown: Array.from(area.querySelectorAll('[data-settings-section]')).filter(visible)
      .map((section) => section.dataset.settingsSection),
    marked: Array.from(pane.querySelectorAll('nav[aria-label="Settings sections"] [aria-current="true"]'),
      (button) => button.getAttribute('aria-label')),
    stops: stops.map(place).filter((name, order, names) => name !== names[order - 1]),
    lastStop: stops.at(-1)?.getAttribute('aria-label') || '',
    hiddenStops: stops.filter((element) => box(element).width === 0 || box(element).height === 0).length,
    generate: { position: getComputedStyle(bar).position, width: Math.round(box(bar).width),
      height: Math.round(box(bar).height), left: Math.round(box(bar).left), bottomGap: Math.round(innerHeight - box(bar).bottom) },
    generateButton: (({ width, height, bottom }) => ({ size: `${Math.round(width)}x${Math.round(height)}`,
      bottomGap: Math.round(innerHeight - bottom) }))(box(bar.querySelector(':scope > button'))),
    focus: document.activeElement?.getAttribute('aria-label') || document.activeElement?.id || document.activeElement?.tagName,
    resizeHandle: visible(document.querySelector('.resize-handle')),
    stored
  };
}, STORAGE_KEY);

// Where a section is inside the docked panel or the stacked settings.
const sectionState = (page, key) => page.evaluate((key) => {
  const area = document.querySelector('.settings-scroll');
  const section = area.querySelector(`[data-settings-section="${key}"]`);
  const details = section.querySelector(':scope > details:first-child');
  const heading = details ? details.querySelector(':scope > summary') : section.querySelector(':scope > [tabindex="-1"]');
  return {
    open: details ? details.open : null,
    offset: Math.round(section.getBoundingClientRect().top - area.getBoundingClientRect().top),
    atEnd: area.scrollTop + area.clientHeight >= area.scrollHeight - 1,
    headingFocused: document.activeElement === heading
  };
}, key);

// Whether an element lies inside what the settings scroll area shows.
const inSettingsView = (locator) => locator.evaluate((element) => {
  const area = document.querySelector('.settings-scroll').getBoundingClientRect();
  const box = element.getBoundingClientRect();
  return box.height > 0 && box.top >= area.top - 1 && box.bottom <= area.bottom + 1;
});

// The entries of the rail in the column that scrolls them.
const railFit = (page) => page.evaluate(() => {
  const nav = document.querySelector('nav[aria-label="Settings sections"]');
  const view = nav.getBoundingClientRect();
  const bar = document.querySelector('.generate-bar').getBoundingClientRect();
  const entries = Array.from(nav.querySelectorAll('button'));
  const boxes = entries.map((button) => button.getBoundingClientRect());
  const names = entries.map((button) => button.querySelector('span'));
  const outside = (inner, outer) => inner.left < outer.left - 0.5 || inner.right > outer.right + 0.5
    || inner.top < outer.top - 0.5 || inner.bottom > outer.bottom + 0.5;
  return {
    scrolls: nav.scrollHeight > nav.clientHeight,
    scrollbar: nav.offsetWidth - nav.clientWidth,
    // Entries the column does not show in full, by short name.
    cut: entries.filter((button, order) => outside(boxes[order], view)).map((button) => button.textContent.trim()),
    underBar: boxes.filter((box) => box.bottom > bar.top).length,
    overlap: boxes.filter((box, order) => order > 0 && box.top < boxes[order - 1].bottom - 0.5).length,
    smallestName: Math.min(...names.map((name) => parseFloat(getComputedStyle(name).fontSize))),
    cutNames: names.filter((name, order) => outside(name.getBoundingClientRect(), boxes[order])).length,
    minHeight: Math.round(Math.min(...boxes.map((box) => box.height)))
  };
});

const near = (actual, expected, label) => expect(Math.abs(actual - expected), `${label}: ${actual} px, expected about ${expected}`)
  .toBeLessThanOrEqual(3);

for (const mode of ['circular', 'linear']) {
  test(`${mode}: the switch turns the pane into the rail and back, and controls keep their node and draft`, async ({ page }) => {
    test.setTimeout(240_000);
    await page.setViewportSize({ width: 1600, height: 1000 });
    await openApp(page);
    await setMode(page, mode);
    const rail = RAIL[mode];
    const rowTitles = rail.filter(([short]) => short !== 'About').map(([, title]) => `Go to ${title}`);

    // The stacked display is the default: the row, the switch at its right end, and the resize handle.
    const stacked = await paneState(page);
    expect(stacked).toMatchObject({ display: 'stacked', docked: null, paneWidth: 320, areaShown: true, resizeHandle: true,
      lastStop: 'Generate Diagram', stored: null });
    near(stacked.resultWidth, 1228, 'figure area, stacked');
    expect(stacked.stops.slice(0, 3)).toEqual(['entry', 'switch', 'section']);
    await expect(toStacked(page)).toHaveCount(0);
    await expect(closeSection(page)).toHaveCount(0);
    expect(await index(page).getByRole('button').evaluateAll((buttons) => buttons.map((button) => button.getAttribute('aria-label'))))
      .toEqual(rowTitles);
    const row = await page.evaluate(() => {
      const wrap = document.querySelector('.settings-index').getBoundingClientRect();
      const toggle = document.querySelector('.settings-index [aria-label="Collapse settings to icons"]').getBoundingClientRect();
      const first = document.querySelector('nav[aria-label="Settings sections"] button').getBoundingClientRect();
      return { right: Math.round(wrap.right - toggle.right), sameRow: Math.abs(toggle.top - first.top) <= 1,
        size: `${Math.round(toggle.width)}x${Math.round(toggle.height)}` };
    });
    expect(row.sameRow, 'the switch is on the first row').toBe(true);
    expect(row.size).toBe('26x26');
    expect(row.right, 'the switch ends the row').toBeLessThanOrEqual(8);

    // A typed draft and a marked node, in a section that the rail will hide.
    const prefix = page.locator('#output-prefix');
    await prefix.fill('rail-draft');
    await prefix.evaluate((element) => { element.dataset.probe = 'kept'; });
    // The typed text is one History step once the focus leaves the box.
    await prefix.blur();
    await page.waitForTimeout(150);
    const history = await undoCount(page);
    await entry(page, 'Colors').click();
    expect(await sectionState(page, 'colors')).toMatchObject({ open: true, headingFocused: true });

    // To the rail: the marked entry is docked, the focus is on the switch back.
    await toRail(page).click();
    await expect(page.locator('.settings-pane')).toHaveAttribute('data-settings-display', 'rail');
    await expect(toStacked(page)).toBeFocused();
    await expect(toRail(page)).toHaveCount(0);
    const docked = await paneState(page);
    expect(docked).toMatchObject({ display: 'rail', docked: 'colors', paneWidth: 388, areaShown: true, areaWidth: 320,
      shown: ['colors'], marked: ['Go to Colors'], resizeHandle: false, lastStop: 'Generate Diagram', hiddenStops: 0,
      stored: 'rail' });
    near(docked.resultWidth, 1166, 'figure area, rail with a section');
    // Tab order follows the visual order: the switch, the entries, the close button, the section, Generate.
    expect(docked.stops).toEqual(['switch', 'entry', 'Close section', 'section', 'generate']);
    expect(docked.generate).toMatchObject({ position: 'fixed', width: 388, left: 16, bottomGap: 32 });
    expect(await sectionState(page, 'colors')).toMatchObject({ open: true, offset: 0 });
    expect(await index(page).getByRole('button').evaluateAll((buttons) => buttons.map((button) => [
      button.textContent.trim(), button.getAttribute('aria-label'), button.getAttribute('title'),
      button.querySelector('i')?.getAttribute('aria-hidden')
    ]))).toEqual(rail.map(([short, title]) => [short, `Go to ${title}`, null, 'true']));
    const column = await page.evaluate(() => {
      const wrap = document.querySelector('.settings-index').getBoundingClientRect();
      const buttons = Array.from(document.querySelectorAll('.settings-index button'), (button) => button.getBoundingClientRect());
      const about = buttons.at(-1);
      return { width: Math.round(wrap.width), first: document.querySelector('.settings-index button').getAttribute('aria-label'),
        inside: buttons.every((box) => box.left >= wrap.left && box.right <= wrap.right && box.top >= wrap.top && box.bottom <= wrap.bottom),
        aboutAtBottom: Math.round(wrap.bottom - about.bottom) <= 8, entryHeight: buttons[1].height };
    });
    // Below the drawing bar a docked Linear rail is about 9 px short of 48 px entries at 1600x1000 (OV-1785 A:
    // only windows at most 800 px high are compacted), so its entries shrink a little; Circular's keep 48 px.
    expect(column).toMatchObject({ width: 60, first: 'Expand settings', inside: true, aboutAtBottom: true });
    if (mode === 'circular') expect(Math.round(column.entryHeight), 'entry height').toBe(48);
    else expect(column.entryHeight, 'entry height').toBeGreaterThanOrEqual(47);
    // The hidden control is the same node with its draft.
    expect(await prefix.evaluate((element) => [element.dataset.probe, element.value, element.checkVisibility()]))
      .toEqual(['kept', 'rail-draft', false]);

    // An entry docks its section: open, at the top of the panel, focus on its heading. Pressing it again keeps it.
    for (const press of [1, 2]) {
      await entry(page, 'Labels').click();
      await expect(page.locator('.settings-pane')).toHaveAttribute('data-settings-docked', 'labels');
      expect(await sectionState(page, 'labels'), `press ${press}`).toMatchObject({ open: true, offset: 0, headingFocused: true });
      expect(await paneState(page)).toMatchObject({ docked: 'labels', shown: ['labels'], marked: ['Go to Labels'] });
    }
    // The section's own disclosure keeps working inside the panel.
    await page.locator('[data-settings-section="labels"] > details > summary').press('Enter');
    expect((await sectionState(page, 'labels')).open).toBe(false);
    await entry(page, 'Labels').press('Enter');
    expect(await sectionState(page, 'labels')).toMatchObject({ open: true, headingFocused: true });

    // The hover text of a rail entry is the full title, to the right of the entry.
    await entry(page, 'Titles & Record Labels').hover();
    const tip = page.locator('.settings-index-tip');
    await expect(tip).toHaveText('Titles & Record Labels');
    const [target, shown] = [await entry(page, 'Titles & Record Labels').boundingBox(), await tip.boundingBox()];
    expect(Math.abs(shown.x - (target.x + target.width) - 6)).toBeLessThanOrEqual(1);
    expect(Math.abs((shown.y + shown.height / 2) - (target.y + target.height / 2))).toBeLessThanOrEqual(1);
    await page.mouse.move(900, 500);
    await expect(tip).toBeHidden();

    // Close: the rail alone, the focus back on the entry that was docked, a compact Generate.
    await closeSection(page).click();
    await expect(entry(page, 'Labels')).toBeFocused();
    const alone = await paneState(page);
    expect(alone).toMatchObject({ display: 'rail', docked: null, paneWidth: 60, areaShown: false, shown: [], marked: [],
      stops: ['switch', 'entry', 'generate'], lastStop: 'Generate Diagram', hiddenStops: 0, resizeHandle: false, stored: 'rail' });
    near(alone.resultWidth, 1494, 'figure area, rail alone');
    expect(alone.generate).toMatchObject({ position: 'fixed', width: 60, left: 16, bottomGap: 32 });
    // The compact button is 60 px square and ends where the stacked button ends.
    expect(alone.generateButton).toEqual({ size: '60x60', bottomGap: stacked.generateButton.bottomGap });
    await expect(closeSection(page)).toHaveCount(0);
    await expect(generate(page)).toBeVisible();
    await expect(generate(page)).toHaveText('Generate', { useInnerText: true });
    // Nothing of the rail is under the bar.
    expect(await page.evaluate(() => Math.round(document.querySelector('.generate-bar').getBoundingClientRect().top
      - document.querySelector('.settings-index').getBoundingClientRect().bottom))).toBeGreaterThanOrEqual(4);

    // From the rail alone the stacked settings come back where they were (Colors at the top).
    await toStacked(page).click();
    await expect(toRail(page)).toBeFocused();
    const back = await paneState(page);
    expect(back).toMatchObject({ display: 'stacked', docked: null, paneWidth: 320, areaShown: true, resizeHandle: true,
      marked: ['Go to Colors'], stored: 'stacked' });
    near(back.resultWidth, 1228, 'figure area, stacked again');
    const colors = await sectionState(page, 'colors');
    expect(Math.abs(colors.offset) <= 2 || colors.atEnd, `Colors is where it was (offset ${colors.offset})`).toBe(true);
    expect(await prefix.evaluate((element) => [element.dataset.probe, element.value])).toEqual(['kept', 'rail-draft']);

    // With a section docked, the stacked display opens it at the top as its row entry would.
    await toRail(page).click();
    await entry(page, 'Features').click();
    await page.locator('[data-settings-section="features"] > details > summary').press('Enter');
    await toStacked(page).click();
    await expect(toRail(page)).toBeFocused();
    const features = await sectionState(page, 'features');
    expect(features.open).toBe(true);
    expect(Math.abs(features.offset) <= 2 || features.atEnd, `Features at the top (offset ${features.offset})`).toBe(true);
    expect(await paneState(page)).toMatchObject({ display: 'stacked', marked: ['Go to Features'], stored: 'stacked' });
    // The dragged width returns with the stacked display.
    await page.evaluate(() => { window.__GBDRAW_APP__.sidebarWidth = 400; });
    await toRail(page).click();
    expect((await paneState(page)).paneWidth).toBe(388);
    await toStacked(page).click();
    expect((await paneState(page)).paneWidth).toBe(400);

    // The switch, an entry, and the close button begin no History step (a click commits on a timer).
    await page.waitForTimeout(150);
    expect(await undoCount(page)).toBe(history);
  });
}

test('Linear: Input docks without Comparison, Compare docks alone, and a mode change keeps a docked entry valid', async ({ page }) => {
  test.setTimeout(240_000);
  await page.setViewportSize({ width: 1600, height: 1000 });
  await openApp(page);
  await setMode(page, 'linear');
  await toRail(page).click();
  // The top of the stacked settings is Input: its card, Basic, and no Comparison.
  expect(await paneState(page)).toMatchObject({ docked: 'input-genomes', shown: ['input-genomes', 'basic'], marked: ['Go to Input Genomes'] });
  await expect(page.locator('#output-prefix')).toBeVisible();

  // Compare: the Comparison section alone, without the card around it.
  await entry(page, 'Comparison').click();
  const compare = await page.evaluate(() => {
    const section = document.querySelector('[data-settings-section="comparison"]');
    const card = document.querySelector('[data-settings-section="input-genomes"]');
    const style = getComputedStyle(card);
    const visible = (element) => element.checkVisibility({ visibilityProperty: true, opacityProperty: false });
    const area = document.querySelector('.settings-scroll').getBoundingClientRect();
    return {
      section: visible(section), headingFocused: document.activeElement === document.getElementById('linear-comparison-heading'),
      cardHeader: visible(card.querySelector(':scope > .card-header')),
      frame: [style.borderTopWidth, style.paddingTop, style.paddingLeft, style.boxShadow, style.backgroundColor],
      otherContent: Array.from(card.querySelectorAll('input, select, textarea, button, a[href]'))
        .filter((control) => visible(control) && !section.contains(control)).length,
      offset: Math.round(section.getBoundingClientRect().top - area.top),
      width: Math.round(section.getBoundingClientRect().width),
      // Nothing beside the section takes room: the card is as tall as the section.
      roomBeside: Math.round(card.getBoundingClientRect().height - section.getBoundingClientRect().height)
    };
  });
  expect(compare).toEqual({ section: true, headingFocused: true, cardHeader: false,
    frame: ['0px', '0px', '0px', 'none', 'rgba(0, 0, 0, 0)'], otherContent: 0, offset: 0, width: 312, roomBeside: 0 });
  expect(await paneState(page)).toMatchObject({ docked: 'comparison', shown: ['input-genomes', 'comparison'], marked: ['Go to Comparison'] });
  await expect(page.locator('#output-prefix')).toBeHidden();

  // Circular has no Comparison entry: Input is docked in its place, with content on screen.
  await setMode(page, 'circular');
  await expect(page.locator('.settings-pane')).toHaveAttribute('data-settings-docked', 'input-genomes');
  expect(await paneState(page)).toMatchObject({ display: 'rail', shown: ['input-genomes', 'basic'], marked: ['Go to Input Genomes'] });
  await expect(page.locator('#output-prefix')).toBeVisible();
  // The same through a mode written by Undo or a loaded Session.
  await setMode(page, 'linear');
  await entry(page, 'Advanced comparison and layout').click();
  await expect(page.locator('.settings-pane')).toHaveAttribute('data-settings-docked', 'advanced-comparison');
  await page.evaluate(() => { window.__GBDRAW_APP__.selectDrawing('circular'); });
  await expect(page.locator('.settings-pane')).toHaveAttribute('data-settings-docked', 'input-genomes');
  // An entry both modes have stays docked; the rail alone stays alone.
  await entry(page, 'Labels').click();
  await setMode(page, 'linear');
  expect(await paneState(page)).toMatchObject({ docked: 'labels', shown: ['labels'], marked: ['Go to Labels'] });
  await closeSection(page).click();
  await setMode(page, 'circular');
  expect(await paneState(page)).toMatchObject({ display: 'rail', docked: null, areaShown: false, marked: [] });

  // About & Citation is the last entry and docks its card.
  await entry(page, 'About & Citation').click();
  expect(await paneState(page)).toMatchObject({ docked: 'about', shown: ['about'], marked: ['Go to About & Citation'] });
  expect(await sectionState(page, 'about')).toMatchObject({ open: true, headingFocused: true, offset: 0 });

  // Axis & Scale is one card per mode under one title: the docked entry follows the mode, in both directions.
  await entry(page, 'Axis & Scale').click();
  expect(await paneState(page)).toMatchObject({ docked: 'circular-axis', shown: ['circular-axis'], marked: ['Go to Axis & Scale'] });
  await setMode(page, 'linear');
  await expect(page.locator('.settings-pane')).toHaveAttribute('data-settings-docked', 'linear-axis');
  expect(await paneState(page)).toMatchObject({ shown: ['linear-axis'], marked: ['Go to Axis & Scale'] });
  await page.evaluate(() => { window.__GBDRAW_APP__.selectDrawing('circular'); });
  await expect(page.locator('.settings-pane')).toHaveAttribute('data-settings-docked', 'circular-axis');
  expect(await paneState(page)).toMatchObject({ shown: ['circular-axis'], marked: ['Go to Axis & Scale'] });

  // The row has no entry for About & Citation: after Expand settings from it, the row marks the section at the
  // top of the settings at once, by the rule it follows while the settings scroll.
  await entry(page, 'About & Citation').click();
  await toStacked(page).click();
  const row = await page.evaluate(() => {
    const area = document.querySelector('.settings-scroll');
    const entries = Array.from(document.querySelectorAll('nav[aria-label="Settings sections"] button'));
    const reached = entries.filter((button) => area.querySelector(`[data-settings-section="${button.dataset.settingsEntry}"]`)
      .getBoundingClientRect().top <= area.getBoundingClientRect().top + 1);
    return { marked: entries.filter((button) => button.getAttribute('aria-current') === 'true').map((button) => button.getAttribute('aria-label')),
      atTop: reached.at(-1).getAttribute('aria-label'), aboutShown: area.querySelector('[data-settings-section="about"] > details').open };
  });
  expect(row.marked).toEqual([row.atTop]);
  expect(row.aboutShown).toBe(true);
});

test('the browser remembers the display and nothing else; storage that fails gives the stacked display', async ({ page, browser, baseURL }) => {
  test.setTimeout(300_000);
  await page.setViewportSize({ width: 1280, height: 800 });
  const errors = [];
  page.on('pageerror', (error) => errors.push(String(error)));
  page.on('console', (message) => { if (message.type() === 'error') errors.push(message.text()); });
  await openApp(page);
  expect(await page.evaluate(() => Object.keys(localStorage))).toEqual([]);
  await entry(page, 'Labels').click();
  await toRail(page).click();
  await expect(page.locator('.settings-pane')).toHaveAttribute('data-settings-docked', 'labels');
  await page.evaluate(() => { window.__GBDRAW_APP__.sidebarWidth = 420; });
  expect(await page.evaluate(() => Object.entries(localStorage))).toEqual([[STORAGE_KEY, 'rail']]);

  // A reload opens in the rail with Input docked: the docked entry and the width are not stored.
  // The first frame the app paints is already the rail.
  await page.addInitScript(() => {
    window.__firstDisplays = [];
    const read = () => {
      const pane = document.querySelector('.settings-pane');
      if (pane?.dataset.settingsDisplay) window.__firstDisplays.push(pane.dataset.settingsDisplay);
      if (window.__firstDisplays.length < 3) requestAnimationFrame(read);
    };
    requestAnimationFrame(read);
  });
  await openApp(page);
  expect(await paneState(page)).toMatchObject({ display: 'rail', docked: 'input-genomes', paneWidth: 388, marked: ['Go to Input Genomes'], stored: 'rail' });
  expect(await page.evaluate(() => window.__firstDisplays)).toEqual(['rail', 'rail', 'rail']);
  await toStacked(page).click();
  expect(await paneState(page)).toMatchObject({ display: 'stacked', paneWidth: 320, stored: 'stacked' });
  expect(await page.evaluate(() => Object.keys(localStorage))).toEqual([STORAGE_KEY]);

  // A stored value the app does not know gives the default and is left alone until the next switch.
  await page.evaluate((key) => localStorage.setItem(key, 'icons'), STORAGE_KEY);
  await openApp(page);
  expect(await paneState(page)).toMatchObject({ display: 'stacked', docked: null, stored: 'icons' });
  expect(errors).toEqual([]);

  // Storage that throws on every access: the default, no error, and the switch still works for this page.
  const context = await browser.newContext({ baseURL, viewport: { width: 1280, height: 800 } });
  try {
    const blocked = await context.newPage();
    const blockedErrors = [];
    blocked.on('pageerror', (error) => blockedErrors.push(String(error)));
    blocked.on('console', (message) => { if (message.type() === 'error') blockedErrors.push(message.text()); });
    await blocked.addInitScript(() => {
      Object.defineProperty(window, 'localStorage', {
        configurable: true,
        get() { throw new DOMException('The operation is insecure.', 'SecurityError'); }
      });
    });
    await openApp(blocked);
    expect(await paneState(blocked)).toMatchObject({ display: 'stacked', docked: null, stored: 'unreadable' });
    await toRail(blocked).click();
    expect(await paneState(blocked)).toMatchObject({ display: 'rail', docked: 'input-genomes', paneWidth: 388 });
    await toStacked(blocked).click();
    expect(await paneState(blocked)).toMatchObject({ display: 'stacked', paneWidth: 320 });
    await toRail(blocked).click();
    await openApp(blocked);
    expect(await paneState(blocked)).toMatchObject({ display: 'stacked', docked: null });
    expect(blockedErrors).toEqual([]);

    // A storage that reads but refuses a write (quota, private mode): the display still switches.
    const full = await context.newPage();
    await full.addInitScript(() => {
      Storage.prototype.setItem = () => { throw new DOMException('The quota has been exceeded.', 'QuotaExceededError'); };
    });
    full.on('pageerror', (error) => blockedErrors.push(String(error)));
    full.on('console', (message) => { if (message.type() === 'error') blockedErrors.push(message.text()); });
    await openApp(full);
    await toRail(full).click();
    expect(await paneState(full)).toMatchObject({ display: 'rail', docked: 'input-genomes', stored: null });
    expect(blockedErrors).toEqual([]);
  } finally {
    await context.close();
  }
});

test('a path that leads to a setting docks its section when the rail is alone', async ({ page }) => {
  test.setTimeout(420_000);
  await page.setViewportSize({ width: 1600, height: 1000 });
  await openApp(page);
  await loadSession(page, LINEAR_SESSION);
  await toRail(page).click();
  const alone = async () => {
    if (await closeSection(page).count()) await closeSection(page).click();
    expect(await paneState(page)).toMatchObject({ display: 'rail', docked: null, areaShown: false });
  };
  const inView = (locator) => inSettingsView(locator);
  const comparisonAtTop = () => page.evaluate(() => {
    const card = document.querySelector('[data-linear-comparison-card]');
    // Rounded to whole pixels; `|| 0` reads a sub-pixel -0.4 as 0, not -0.
    return { offset: Math.round(card.getBoundingClientRect().top - document.querySelector('.settings-scroll').getBoundingClientRect().top) || 0,
      focusInside: card.contains(document.activeElement) && document.activeElement.matches('button') };
  });

  // The record popup's Show in settings: the record's field, in Input.
  await alone();
  const definition = page.locator('.gbdraw-preview-surface [data-gbdraw-role^="record-definition"] text').first();
  const box = await definition.boundingBox();
  await page.mouse.click(box.x + box.width / 2, box.y + box.height / 2);
  await page.getByRole('button', { name: 'Show in settings', exact: true }).click();
  const field = page.getByLabel('Definition for sequence 1', { exact: true });
  await expect(field).toBeFocused();
  expect(await paneState(page)).toMatchObject({ docked: 'input-genomes', marked: ['Go to Input Genomes'] });
  expect(await inView(field), 'the field is inside the panel').toBe(true);

  // The cause of a comparison error: the Execution setting, in Compare. Its disclosure is closed first.
  await alone();
  await page.evaluate(() => {
    document.querySelector('[data-linear-comparison-disclosure="settings"]').open = false;
    window.__GBDRAW_APP__.errorLog = { code: 'LOSAT_RUNTIME', stage: 'comparison' };
  });
  await page.getByRole('button', { name: 'Show in Comparison', exact: true }).click();
  const execution = page.getByRole('combobox', { name: 'LOSAT execution' });
  await expect(execution).toBeFocused();
  expect(await paneState(page)).toMatchObject({ docked: 'comparison', marked: ['Go to Comparison'] });
  expect(await inView(execution), 'the Execution setting is inside the panel').toBe(true);
  // A comparison error whose cause has no control: the Comparison section at the top, its first button focused.
  await alone();
  await page.evaluate(() => { window.__GBDRAW_APP__.errorLog = { code: 'COMPARISON_IDENTITY', stage: 'comparison' }; });
  await page.getByRole('button', { name: 'Show in Comparison', exact: true }).click();
  expect(await paneState(page)).toMatchObject({ docked: 'comparison', marked: ['Go to Comparison'] });
  expect(await comparisonAtTop()).toEqual({ offset: 0, focusInside: true });
  await page.evaluate(() => { window.__GBDRAW_APP__.errorLog = null; });

  // Record Labels visibility, reached from the Layout note: the field's setting, in Titles.
  await alone();
  const visibility = page.locator('select[id^="linear-label-visibility-"]').first();
  const key = (await visibility.getAttribute('id')).replace('linear-label-visibility-', '');
  await page.evaluate((key) => window.__GBDRAW_APP__.focusLinearLabelVisibility(key), key);
  await expect(visibility).toBeFocused();
  expect(await paneState(page)).toMatchObject({ docked: 'titles', marked: ['Go to Titles & Record Labels'] });
  expect(await inView(visibility), 'the visibility setting is inside the panel').toBe(true);

  // The Circular Multi-Record Canvas setting, reached from Input: in Layout.
  await setMode(page, 'circular');
  await alone();
  await page.evaluate(() => { window.__GBDRAW_APP__.form.multi_record_canvas = true; });
  expect(await page.evaluate(async () => (await window.__GBDRAW_APP__.showCircularCanvasSetting()).status)).toBe('ok');
  const canvas = page.locator('[data-circular-canvas-setting]');
  await expect(canvas).toBeFocused();
  const state = await paneState(page);
  expect(state.docked).toBe(await canvas.evaluate((element) => element.closest('[data-settings-section]').dataset.settingsSection));
  expect(await inView(canvas), 'the canvas setting is inside the panel').toBe(true);

  // The stacked display: the same paths open the section and bring the setting into view.
  await toStacked(page).click();
  await page.evaluate(() => {
    document.querySelectorAll('.settings-scroll > [data-settings-section] > details').forEach((details) => { details.open = false; });
    document.querySelector('.settings-scroll').scrollTop = 0;
  });
  expect(await page.evaluate(async () => (await window.__GBDRAW_APP__.showCircularCanvasSetting()).status)).toBe('ok');
  await expect(canvas).toBeFocused();
  expect(await inView(canvas), 'stacked: the canvas setting is in view').toBe(true);
  expect(await page.evaluate(() => [document.querySelector('.app-main').scrollTop, document.querySelector('.settings-pane').scrollTop,
    document.scrollingElement.scrollTop]), 'only the settings scroll area moved').toEqual([0, 0, 0]);
  // Stacked, the comparison error without a control: the Comparison section goes to the top.
  await setMode(page, 'linear');
  await page.evaluate(() => {
    document.querySelector('.settings-scroll').scrollTop = 0;
    window.__GBDRAW_APP__.errorLog = { code: 'COMPARISON_IDENTITY', stage: 'comparison' };
  });
  await page.getByRole('button', { name: 'Show in Comparison', exact: true }).click();
  expect(await comparisonAtTop()).toEqual({ offset: 0, focusInside: true });
});

test('phone width shows the stacked settings without the switch, whatever is remembered', async ({ page }) => {
  test.setTimeout(240_000);
  await page.setViewportSize({ width: 390, height: 844 });
  await page.addInitScript((key) => { localStorage.setItem(key, 'rail'); }, STORAGE_KEY);
  await openApp(page);
  await setMode(page, 'linear');
  const rowCount = RAIL.linear.length - 1;
  const phone = await paneState(page);
  expect(phone).toMatchObject({ display: 'stacked', docked: null, areaShown: true, stored: 'rail' });
  expect(phone.paneWidth).toBeGreaterThan(350);
  expect(phone.shown).toEqual(expect.arrayContaining(['input-genomes', 'comparison', 'basic', 'labels', 'about']));
  await expect(index(page).getByRole('button')).toHaveCount(rowCount);
  await expect(toRail(page)).toHaveCount(0);
  await expect(toStacked(page)).toHaveCount(0);
  await expect(closeSection(page)).toHaveCount(0);
  expect(await page.evaluate(() => document.documentElement.scrollWidth - innerWidth)).toBeLessThanOrEqual(0);
  // An entry opens its section below the sticky header, as the stacked row does.
  await entry(page, 'Labels').click();
  const labels = await page.evaluate(() => {
    const section = document.querySelector('[data-settings-section="labels"]');
    return { open: section.querySelector(':scope > details').open,
      belowHeader: Math.round(section.getBoundingClientRect().top - document.querySelector('.app-header').getBoundingClientRect().bottom) };
  });
  expect(labels.open).toBe(true);
  expect(Math.abs(labels.belowHeader)).toBeLessThanOrEqual(2);

  // A wide window shows the remembered rail with Input docked; narrow again, every section is back.
  await page.setViewportSize({ width: 1280, height: 800 });
  await expect(page.locator('.settings-pane')).toHaveAttribute('data-settings-display', 'rail');
  expect(await paneState(page)).toMatchObject({ docked: 'input-genomes', paneWidth: 388, shown: ['input-genomes', 'basic'], stored: 'rail' });
  await expect(index(page).getByRole('button')).toHaveCount(RAIL.linear.length);
  await closeSection(page).click();
  await page.setViewportSize({ width: 390, height: 844 });
  await expect(page.locator('.settings-pane')).toHaveAttribute('data-settings-display', 'stacked');
  const narrow = await paneState(page);
  expect(narrow).toMatchObject({ display: 'stacked', docked: null, areaShown: true, stored: 'rail' });
  expect(narrow.shown).toEqual(phone.shown);
  await expect(index(page).getByRole('button')).toHaveCount(rowCount);
  await expect(toStacked(page)).toHaveCount(0);
  // The Generate bar spans the phone width.
  expect(narrow.generate.width).toBe(390);
});

test('Save Session, Load Session, and Reset Settings leave the display alone', async ({ page }, testInfo) => {
  test.setTimeout(420_000);
  await page.setViewportSize({ width: 1600, height: 1000 });
  await openApp(page);
  await loadSession(page, LINEAR_SESSION);
  await toRail(page).click();
  await entry(page, 'Labels').click();
  const before = await paneState(page);
  expect(before).toMatchObject({ display: 'rail', docked: 'labels', stored: 'rail' });

  const pending = page.waitForEvent('download');
  await page.getByRole('banner').getByRole('button', { name: 'Save Session', exact: true }).click();
  const download = await pending;
  const saved = testInfo.outputPath(download.suggestedFilename());
  await download.saveAs(saved);
  const bytes = readFileSync(saved);
  const text = (bytes[0] === 0x1f && bytes[1] === 0x8b ? gunzipSync(bytes) : bytes).toString('utf8');
  const keys = new Set();
  JSON.parse(text, (key, value) => { keys.add(key); return value; });
  expect([...keys].filter((key) => /settings.?pane|sidebar|docked/i.test(key)), 'the Session holds no pane display').toEqual([]);
  expect(await paneState(page)).toMatchObject({ display: 'rail', docked: 'labels', stored: 'rail' });

  // Loading the Session saved in the rail into the stacked display leaves it stacked.
  await toStacked(page).click();
  await page.locator('input[accept^=".json,"]').setInputFiles(saved);
  await page.waitForFunction(() => !window.__GBDRAW_APP__.sessionImportPending, null, { timeout: 180000 });
  await page.waitForTimeout(300);
  expect(await paneState(page)).toMatchObject({ display: 'stacked', docked: null, stored: 'stacked' });
  // And loading it in the rail keeps the rail and the docked entry.
  await toRail(page).click();
  await entry(page, 'Features').click();
  await page.locator('input[accept^=".json,"]').setInputFiles(saved);
  await page.waitForFunction(() => !window.__GBDRAW_APP__.sessionImportPending, null, { timeout: 180000 });
  await page.waitForTimeout(300);
  expect(await paneState(page)).toMatchObject({ display: 'rail', docked: 'features', shown: ['features'], stored: 'rail' });

  await page.getByRole('button', { name: 'Reset Settings', exact: true }).click();
  await page.locator('[data-choice-prompt] [data-choice-prompt-action="reset"]').click();
  await expect.poll(() => page.evaluate(() => !window.__GBDRAW_HISTORY__.mutationPending())).toBe(true);
  expect(await paneState(page)).toMatchObject({ display: 'rail', docked: 'features', shown: ['features'], stored: 'rail' });

  // While a Session operation is pending, the rail's own buttons stay usable and everything else is disabled.
  await page.evaluate(() => { window.__GBDRAW_APP__.sessionSavePending = true; });
  const enabled = await page.evaluate(() => Array.from(document.querySelectorAll('.settings-pane :is(input, select, textarea, button)'))
    .filter((control) => control.checkVisibility() && !control.closest('.help-tip') && !control.disabled
      && control.getAttribute('aria-disabled') !== 'true' && !control.matches('button[aria-expanded], button[aria-controls], button:has(.ph-download-simple)'))
    .map((control) => control.getAttribute('aria-label') || control.textContent.trim()));
  await page.evaluate(() => { window.__GBDRAW_APP__.sessionSavePending = false; });
  expect(enabled).toEqual([]);
});

// The Review alignment panel is docked in the Result: a right column 24rem wide
// when the Result's workspace is at least 56rem wide, below the figure when it
// is narrower. The pane's display changes that width, so the panel takes the
// side its rule gives in each display and the review goes on.
test('an open Review alignment stays docked in the Result through every display of the pane', async ({ page }) => {
  test.setTimeout(420_000);
  await page.setViewportSize({ width: 1600, height: 1000 });
  await openApp(page);
  await loadSession(page, LINEAR_SESSION);
  await page.locator('[data-gbdraw-feature-id="f94b7f30b_record_2"]').click({ force: true });
  await page.locator('.feature-popup [data-similarity-alignment-popup-review]').click();
  const review = page.getByRole('dialog', { name: 'Review alignment' });
  await expect(review.locator('[data-alignment-record-key]')).toHaveCount(4, { timeout: 180000 });
  const chosen = () => review.locator('input:checked').evaluateAll((inputs) => inputs.map((input) => input.value));
  const choices = await chosen();
  expect(choices.length, 'the review has choices made').toBeGreaterThan(0);

  const layout = () => page.evaluate(() => {
    const box = (selector) => document.querySelector(selector).getBoundingClientRect();
    const [panel, figure, workspace, result, pane, bar] = ['[data-similarity-alignment-dialog]', '.preview-viewport',
      '.preview-workspace', '.result-pane', '.settings-pane', '.generate-bar'].map(box);
    const rem = parseFloat(getComputedStyle(document.documentElement).fontSize);
    const beside = Math.abs(panel.left - figure.right) <= 1 && Math.abs(panel.top - figure.top) <= 1 && Math.round(panel.width) === 24 * rem;
    const below = panel.top >= figure.bottom - 1 && Math.abs(panel.left - figure.left) <= 1 && Math.abs(panel.width - figure.width) <= 2;
    return {
      side: beside ? 'beside' : below ? 'below' : 'neither',
      wide: result.width >= 56.125 * rem,
      resultWidth: Math.round(result.width),
      figureRoom: figure.width >= 300 && figure.height >= 200,
      inWorkspace: panel.left >= workspace.left - 1 && panel.right <= workspace.right + 1
        && panel.bottom <= document.querySelector('[data-app-footer]').getBoundingClientRect().top,
      clearOfSettings: Math.max(pane.right, bar.right) <= result.left
    };
  });
  // The side follows the width the display leaves for the Result.
  const expectDocked = async (state, side, resultWidth) => {
    await expect.poll(async () => (await layout()).side, { message: state }).toBe(side);
    const found = await layout();
    expect(found, state).toMatchObject({ wide: side === 'beside', figureRoom: true, inWorkspace: true, clearOfSettings: true });
    near(found.resultWidth, resultWidth, `Result width, ${state}`);
    await expect(review.locator('[data-alignment-record-key]'), state).toHaveCount(4);
    expect(await chosen(), `${state}: the choices of the review`).toEqual(choices);
  };

  // 1600 px: the Result stays wide in every display, so the panel stays beside the figure.
  await expectDocked('1600 stacked', 'beside', 1228);
  await toRail(page).click();
  await expectDocked('1600 rail with a section', 'beside', 1166);
  await closeSection(page).click();
  await expectDocked('1600 rail alone', 'beside', 1494);
  await entry(page, 'Labels').click();
  await expectDocked('1600 rail with Labels', 'beside', 1166);
  await toStacked(page).click();
  await expectDocked('1600 stacked again', 'beside', 1228);

  // 1280 px: the docked section takes 62 px more than the stacked settings, which leaves the
  // Result narrower than 56rem; the panel is below the figure there, and beside it again without the section.
  await page.setViewportSize({ width: 1280, height: 800 });
  await expectDocked('1280 stacked', 'beside', 908);
  await toRail(page).click();
  await expectDocked('1280 rail with a section', 'below', 846);
  await closeSection(page).click();
  await expectDocked('1280 rail alone', 'beside', 1174);
  await toStacked(page).click();
  await expectDocked('1280 stacked again', 'beside', 908);
});

test('at 1280x720 the rail shows every entry without scrolling, with a section docked and alone', async ({ page }) => {
  test.setTimeout(180_000);
  await page.setViewportSize({ width: 1280, height: 720 });
  await page.addInitScript((key) => { localStorage.setItem(key, 'rail'); }, STORAGE_KEY);
  await openApp(page);
  for (const mode of ['circular', 'linear']) {
    await setMode(page, mode);
    for (const closed of [false, true]) {
      if (closed) await closeSection(page).click();
      const fit = await railFit(page);
      expect(fit, `${mode}, ${closed ? 'rail alone' : 'docked'}`).toMatchObject({ scrolls: false, cut: [], underBar: 0, overlap: 0, cutNames: 0 });
      expect(fit.smallestName, `${mode}: name size`).toBeGreaterThanOrEqual(9);
      expect(fit.minHeight, `${mode}: entry height`).toBeGreaterThanOrEqual(30);
      if (closed) await entry(page, 'Input Genomes').click();
    }
  }
});

// A lower window: the column scrolls. Headless Chromium hides scrollbars, so
// this case runs in its own browser with scrollbars shown.
test('a rail too low for its entries scrolls with a scrollbar and brings the focused and the docked entry into view (1366x650)', async ({ playwright, baseURL }) => {
  test.setTimeout(180_000);
  const browser = await playwright.chromium.launch({ ignoreDefaultArgs: ['--hide-scrollbars'] });
  try {
    const page = await browser.newPage({ baseURL, viewport: { width: 1366, height: 650 } });
    await page.addInitScript((key) => { localStorage.setItem(key, 'rail'); }, STORAGE_KEY);
    await openApp(page);
    await setMode(page, 'linear');
    const column = index(page);
    const toTop = () => column.evaluate((nav) => { nav.scrollTop = 0; });
    const tip = page.locator('.settings-index-tip');

    // Linear, Input docked: the last entries are below the column's edge, and a 4 px scrollbar shows it.
    const low = await railFit(page);
    expect(low).toMatchObject({ scrolls: true, scrollbar: 4, overlap: 0, cutNames: 0 });
    expect(low.cut).toContain('About');
    expect(low.smallestName).toBeGreaterThanOrEqual(9);
    expect(low.minHeight).toBeGreaterThanOrEqual(30);

    // The text of a hovered entry goes when the column scrolls under the pointer.
    await entry(page, 'Input Genomes').hover();
    await expect(tip).toHaveText('Input Genomes');
    await page.mouse.wheel(0, 60);
    await expect.poll(() => column.evaluate((nav) => nav.scrollTop)).toBeGreaterThan(0);
    await expect(tip).toHaveCount(0);
    await page.mouse.move(700, 300);

    // The keyboard reaches every entry: the focused one comes into view, and its text stays beside it.
    await toTop();
    await toStacked(page).focus();
    for (let stop = 0; stop < RAIL.linear.length; stop += 1) await page.keyboard.press('Tab');
    await expect(entry(page, 'About & Citation')).toBeFocused();
    expect((await railFit(page)).cut).not.toContain('About');
    await expect(tip).toHaveText('About & Citation');
    await expect.poll(async () => {
      const [target, shown] = [await entry(page, 'About & Citation').boundingBox(), await tip.boundingBox()];
      return Math.abs((shown.y + shown.height / 2) - (target.y + target.height / 2));
    }).toBeLessThanOrEqual(2);

    // Close section returns the focus to the docked entry, in view.
    await entry(page, 'About & Citation').click();
    await toTop();
    await closeSection(page).click();
    await expect(entry(page, 'About & Citation')).toBeFocused();
    expect((await railFit(page)).cut).not.toContain('About');

    // A path that docks an entry below the column's edge brings that entry into view.
    await toTop();
    await page.evaluate(() => window.__GBDRAW_APP__.revealSetting(document.querySelector('[data-settings-section="about"] summary')));
    await expect(page.locator('.settings-pane')).toHaveAttribute('data-settings-docked', 'about');
    await expect.poll(async () => (await railFit(page)).cut).not.toContain('About');
  } finally {
    await browser.close();
  }
});

// A color rule with a pattern Python rejects fails with an error in the Result
// that offers Edit pattern: the pattern field of the Colors section.
test('Edit pattern of a failed color rule leads to the pattern field in every display', async ({ page }) => {
  test.setTimeout(420_000);
  await page.setViewportSize({ width: 1600, height: 1000 });
  await openApp(page);
  await loadSession(page, 'HmmtDNA_basic_circular.gbdraw-session.json');
  await entry(page, 'Colors').click();
  await page.getByLabel('Qualifier name', { exact: true }).fill('product');
  const pattern = page.getByLabel('New color rule pattern', { exact: true });
  await pattern.fill('(?<enzyme>NADH)');
  await page.getByRole('button', { name: 'Add color rule', exact: true }).click();
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.errorLog?.code ?? null), { timeout: 180_000 }).toBe('REGEX_SYNTAX');
  const edit = page.getByRole('button', { name: 'Edit pattern', exact: true });
  const area = page.locator('.settings-scroll');

  // Stacked: the field out of view, then inside a closed Colors section.
  await area.evaluate((element) => { element.scrollTop = 0; });
  expect(await inSettingsView(pattern)).toBe(false);
  await edit.click();
  await expect(pattern).toBeFocused();
  expect(await inSettingsView(pattern), 'stacked: the field is in view').toBe(true);
  await area.evaluate((element) => {
    element.querySelector('[data-settings-section="colors"] > details').open = false;
    element.scrollTop = 0;
  });
  await edit.click();
  await expect(pattern).toBeFocused();
  expect(await inSettingsView(pattern), 'stacked, Colors closed: the field is in view').toBe(true);

  // The rail with another section docked, then the rail alone: Colors is docked and the field focused.
  await toRail(page).click();
  for (const display of ['another section docked', 'rail alone']) {
    if (display === 'rail alone') await closeSection(page).click();
    else await entry(page, 'Layout').click();
    expect((await paneState(page)).docked, display).toBe(display === 'rail alone' ? null : 'layout');
    await edit.click();
    await expect(pattern, display).toBeFocused();
    expect(await paneState(page), display).toMatchObject({ docked: 'colors', marked: ['Go to Colors'] });
    expect(await inSettingsView(pattern), `${display}: the field is inside the panel`).toBe(true);
  }
});

// Beside the rail alone Generate is a 60 px button that ends where the stacked
// button ends, and it is the last Tab stop of the pane in every display.
test('beside the rail alone Generate is compact, also while it is disabled or running', async ({ page }) => {
  test.setTimeout(180_000);
  await page.setViewportSize({ width: 1024, height: 768 });
  await openApp(page);
  const lastStop = () => page.evaluate(() => {
    const visible = (element) => element.checkVisibility({ visibilityProperty: true, opacityProperty: false });
    return Array.from(document.querySelector('.settings-pane').querySelectorAll('a[href], button, input, select, textarea, summary, [tabindex]'))
      .filter((element) => visible(element) && !element.disabled && element.tabIndex >= 0).at(-1)?.getAttribute('aria-label');
  });
  const stackedGap = (await paneState(page)).generateButton.bottomGap;
  expect(await lastStop(), 'stacked').toBe('Generate Diagram');
  await toRail(page).click();
  expect(await lastStop(), 'rail with a section').toBe('Generate Diagram');
  await closeSection(page).click();
  expect(await lastStop(), 'rail alone').toBe('Generate Diagram');
  const compact = await paneState(page);
  expect(compact.generateButton).toEqual({ size: '60x60', bottomGap: stackedGap });
  expect(compact.generate.width).toBe(60);

  // Disabled by a Session operation, and while Generate runs: the same 60 px button, the loader alone.
  await page.evaluate(() => { window.__GBDRAW_APP__.sessionSavePending = true; });
  await expect(generate(page)).toBeDisabled();
  expect((await paneState(page)).generateButton).toEqual({ size: '60x60', bottomGap: stackedGap });
  await page.evaluate(() => { Object.assign(window.__GBDRAW_APP__, { sessionSavePending: false, processing: true }); });
  await expect(generate(page)).toBeDisabled();
  await expect(page.locator('.generate-bar > button .gbdraw-loader')).toBeVisible();
  await expect(page.locator('.generate-bar > button .generate-bar-long')).toBeHidden();
  expect((await paneState(page)).generateButton).toEqual({ size: '60x60', bottomGap: stackedGap });
  await page.evaluate(() => { window.__GBDRAW_APP__.processing = false; });
  await expect(generate(page)).toBeEnabled();
});

// In the phone layout the page scrolls the settings and the Generate bar is
// fixed over its bottom edge: a setting under the bar is not in view.
test('phone width: Record Labels brings its setting into view above the fixed Generate bar', async ({ page }) => {
  test.setTimeout(180_000);
  await page.setViewportSize({ width: 390, height: 844 });
  await openApp(page);
  await setMode(page, 'linear');
  const visibility = page.locator('select[id^="linear-label-visibility-"]').first();
  const key = (await visibility.getAttribute('id')).replace('linear-label-visibility-', '');
  const place = () => visibility.evaluate((select) => {
    const box = select.getBoundingClientRect();
    return { aboveBar: Math.round(document.querySelector('.generate-bar').getBoundingClientRect().top - box.bottom),
      belowHeader: Math.round(box.top - document.querySelector('.app-header').getBoundingClientRect().bottom) };
  });
  // The setting lies in the strip the bar covers.
  await visibility.evaluate((select) => {
    for (let details = select.closest('details'); details; details = details.parentElement?.closest('details')) details.open = true;
    window.scrollBy(0, select.getBoundingClientRect().top - document.querySelector('.generate-bar').getBoundingClientRect().top - 4);
  });
  expect((await place()).aboveBar).toBeLessThan(0);
  await page.evaluate((key) => window.__GBDRAW_APP__.focusLinearLabelVisibility(key), key);
  await expect(visibility).toBeFocused();
  const shown = await place();
  expect(shown.aboveBar, 'the setting ends above the Generate bar').toBeGreaterThanOrEqual(0);
  expect(shown.belowHeader, 'the setting starts below the header').toBeGreaterThanOrEqual(0);
});

// The close button sits over the top right corner of the docked panel and
// stays there while the section scrolls. With the section at its top, where
// it cannot be scrolled further down, no text and no control lies under it.
test('nothing of a docked section lies under the close button while the section is at its top', async ({ page }) => {
  test.setTimeout(240_000);
  await page.setViewportSize({ width: 1600, height: 1000 });
  await openApp(page);
  await loadSession(page, LINEAR_SESSION);
  await toRail(page).click();
  for (const mode of ['linear', 'circular']) {
    await setMode(page, mode);
    for (const [short, title] of RAIL[mode]) {
      await entry(page, title).click();
      const under = await page.evaluate(() => {
        const area = document.querySelector('.settings-scroll');
        const close = document.querySelector('.settings-dock-close').getBoundingClientRect();
        const visible = (element) => element.checkVisibility({ visibilityProperty: true, opacityProperty: false });
        const covered = (box) => box.width > 0 && box.height > 0 && box.right > close.left && box.left < close.right
          && box.bottom > close.top && box.top < close.bottom;
        const found = [];
        const walker = document.createTreeWalker(area, NodeFilter.SHOW_TEXT);
        for (let node = walker.nextNode(); node; node = walker.nextNode()) {
          if (!node.textContent.trim() || !visible(node.parentElement)) continue;
          const range = document.createRange();
          range.selectNodeContents(node);
          if (Array.from(range.getClientRects()).some(covered)) found.push(node.textContent.trim().slice(0, 40));
        }
        for (const control of area.querySelectorAll('button, input, select, textarea, a[href], i')) {
          if (visible(control) && covered(control.getBoundingClientRect())) found.push(control.getAttribute('aria-label') || control.className);
        }
        return { scrollTop: area.scrollTop, found };
      });
      expect(under, `${mode}, ${short}`).toEqual({ scrollTop: 0, found: [] });
    }
  }
});
