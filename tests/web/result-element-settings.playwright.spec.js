// A click on a drawn element of the Result opens the settings section that
// holds its settings, as the section index does for the same section: the
// scale ticks and the Linear scale open Axis & Scale, a dinucleotide track
// opens Dinucleotide content/skew, and the plot title opens Titles & Record
// Labels. It works while Layout edit is off and no Review alignment is open,
// writes no setting, History step, or Result, and leaves no trace in an export.
const { test, expect } = require('@playwright/test');
const { readFileSync } = require('node:fs');
const { join } = require('node:path');
const { openApp } = require('./helpers/app-lifecycle.cjs');

test.use({ actionTimeout: 20_000 });

const CUSTOM_SLOTS_SESSION = 'HmmtDNA_ATskew.gbdraw-session.json';
const DEFAULT_CIRCULAR_SESSION = 'HmmtDNA_basic_circular.gbdraw-session.json';
const LINEAR_SESSION = 'BGC0000708-BGC0000713.gbdraw-session.json';
const SETTINGS_HINT = 'Click to open its settings';
const MOVE_HINT = 'Turn on Layout edit to move this item';
// One line each: the hint box stays as wide as the longer sentence.
const BOTH_HINTS = `${SETTINGS_HINT} ${MOVE_HINT}`;
const hintLines = (page) => hint(page).evaluate((element) => new Set(Array.from(element.querySelectorAll('span'))
  .filter((span) => span.checkVisibility()).map((span) => Math.round(span.getBoundingClientRect().top))).size);
const MARKER = '[data-gbdraw-preview-item="settings"]';

const index = (page) => page.getByRole('navigation', { name: 'Settings sections' });
const entry = (page, title) => index(page).getByRole('button', { name: `Go to ${title}`, exact: true });
const toRail = (page) => page.getByRole('button', { name: 'Collapse settings to icons', exact: true });
const closeSection = (page) => page.getByRole('button', { name: 'Close section', exact: true });
const layoutEdit = (page) => page.getByRole('button', { name: 'Toggle layout edit mode' });
const hint = (page) => page.locator('.preview-layout-hint');
const undoCount = (page) => page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount());
const frames = (page) => page.evaluate(() => new Promise((resolve) => requestAnimationFrame(() => requestAnimationFrame(resolve))));

const sessionPath = (name) => join(process.cwd(), 'gbdraw/web/gallery/sessions', name);
// `file` is a Gallery Session by name, or a Session file as `setInputFiles` takes it.
const loadSession = async (page, file) => {
  page.on('dialog', (dialog) => dialog.accept());
  await page.locator('input[accept^=".json,"]').setInputFiles(typeof file === 'string' ? sessionPath(file) : file);
  await page.waitForFunction(() => !window.__GBDRAW_APP__.sessionImportPending
    && window.__GBDRAW_APP__.results.length > 0, null, { timeout: 180000 });
  await frames(page);
  // The whole figure in the preview, so that every drawn element can be pointed at.
  await page.getByRole('button', { name: 'Fit to preview', exact: true }).click();
};

// The figure has settled once its box stays where it is for two looks (a fit
// moves it, and the pane's display changes the room it has).
const figureSettled = async (page) => {
  await page.evaluate(() => { window.__figureBox = ''; });
  await page.waitForFunction(() => {
    const svg = document.querySelector('.gbdraw-preview-surface svg');
    const box = svg ? JSON.stringify(svg.getBoundingClientRect()) : '';
    const same = Boolean(box) && window.__figureBox === box;
    window.__figureBox = box;
    return same;
  }, null, { polling: 250 });
};

// A point of the window on the ink of the first element that `selector` finds
// in the Result.
const inkPoint = async (page, selector) => {
  await figureSettled(page);
  const point = await page.evaluate((selector) => {
    const surface = document.querySelector('.gbdraw-preview-surface');
    for (const element of surface.querySelectorAll(selector)) {
      const box = element.getBoundingClientRect();
      for (let row = 1; row < 32; row += 1) {
        for (let column = 1; column < 32; column += 1) {
          const x = box.left + (box.width * column) / 32;
          const y = box.top + (box.height * row) / 32;
          const hit = x >= 0 && y >= 0 && x < innerWidth && y < innerHeight ? document.elementFromPoint(x, y) : null;
          if (hit && element.contains(hit)) return { x, y };
        }
      }
    }
    return null;
  }, selector);
  expect(point, `a point on ${selector}`).not.toBeNull();
  return point;
};
const clickInk = async (page, selector, modifier) => {
  const point = await inkPoint(page, selector);
  if (modifier) await page.keyboard.down(modifier);
  await page.mouse.click(point.x, point.y);
  if (modifier) await page.keyboard.up(modifier);
  return point;
};
const hover = async (page, selector) => {
  const point = await inkPoint(page, selector);
  await page.mouse.move(point.x, point.y);
  return page.evaluate(({ x, y }) => {
    const hit = document.elementFromPoint(x, y);
    const target = hit.closest('[data-gbdraw-preview-item], .gbdraw-preview-layout-target') || hit;
    return { cursor: getComputedStyle(hit).cursor, outline: getComputedStyle(target).outlineStyle };
  }, point);
};

// Where a section is: its disclosure, its distance below the top edge of what
// the settings show, the focus, the marked entry, and the pane's display.
const sectionState = (page, key) => page.evaluate((key) => {
  const pane = document.querySelector('.settings-pane');
  const area = pane.querySelector('.settings-scroll');
  const section = area.querySelector(`[data-settings-section="${key}"]`);
  const details = section.querySelector(':scope > details:first-child');
  const scrolls = getComputedStyle(area).overflowY !== 'visible';
  const top = scrolls ? area.getBoundingClientRect().top : document.querySelector('.app-header').getBoundingClientRect().bottom;
  const owner = scrolls ? area : document.scrollingElement;
  return {
    open: details.open,
    shown: section.checkVisibility(),
    offset: Math.round(section.getBoundingClientRect().top - top),
    atEnd: owner.scrollTop + owner.clientHeight >= owner.scrollHeight - 1,
    headingFocused: document.activeElement === details.querySelector(':scope > summary'),
    marked: Array.from(pane.querySelectorAll('nav[aria-label="Settings sections"] [aria-current="true"]'),
      (button) => button.getAttribute('aria-label')),
    display: pane.dataset.settingsDisplay,
    docked: pane.dataset.settingsDocked ?? null
  };
}, key);

// Every section closed, the settings at their top, the focus nowhere.
const closeSections = (page) => page.evaluate(() => {
  document.querySelectorAll('.settings-scroll [data-settings-section] > details').forEach((details) => { details.open = false; });
  document.querySelector('.settings-scroll').scrollTop = 0;
  document.scrollingElement.scrollTop = 0;
  document.activeElement?.blur();
});

const expectOpened = async (page, key, title, label = title) => {
  await expect.poll(async () => (await sectionState(page, key)).headingFocused, { message: `${label}: focus on the heading` }).toBe(true);
  await frames(page);
  const state = await sectionState(page, key);
  expect(state, label).toMatchObject({ open: true, shown: true, marked: [`Go to ${title}`] });
  expect(Math.abs(state.offset) <= 2 || state.atEnd, `${label}: at the top edge or at the end of the scroll (offset ${state.offset})`).toBe(true);
  return state;
};
const expectClosed = async (page, key, label) => {
  await frames(page);
  await page.waitForTimeout(150);
  expect(await sectionState(page, key), label).toMatchObject({ open: false, headingFocused: false });
};

const appState = (page) => page.evaluate(() => {
  const app = window.__GBDRAW_APP__;
  return JSON.stringify([window.__GBDRAW_HISTORY__.getUndoCount(), app.form, app.adv, app.results.map((result) => result.content)]);
});
const exported = (page) => page.evaluate(async () => {
  const { serializeCleanSvg } = await import('./js/services/svg-serialization.js');
  return serializeCleanSvg(window.__GBDRAW_APP__.svgContainer.querySelector('svg'));
});

const TICKS = 'g[data-gbdraw-slot-renderer="ticks"] text';
const GC_CONTENT = 'g[data-gbdraw-slot-id="gc_content"] path';
const GC_SKEW = 'g[data-gbdraw-slot-id="gc_skew"] path';

test('Circular with custom track slots: a tick label opens Axis & Scale and a dinucleotide track its section', async ({ page }) => {
  test.setTimeout(300_000);
  await page.setViewportSize({ width: 1600, height: 1000 });
  await openApp(page);
  await loadSession(page, CUSTOM_SLOTS_SESSION);
  const before = await appState(page);

  // The targets are marked for the pointer and are no buttons: the section index is the keyboard way.
  const marked = await page.locator(`.gbdraw-preview-surface ${MARKER}`).evaluateAll((elements) => elements.map((element) => ({
    slot: element.getAttribute('data-gbdraw-slot-id'), role: element.getAttribute('role'), tabindex: element.getAttribute('tabindex')
  })));
  expect(marked.map((target) => target.slot).sort()).toEqual(['a_skew_2', 'gc_content', 'gc_skew', 'ticks']);
  expect(marked.every((target) => target.role === null && target.tabindex === null)).toBe(true);

  // Hover: the pointer cursor and the hint, without a frame around the ring.
  expect(await hover(page, TICKS)).toEqual({ cursor: 'pointer', outline: 'none' });
  await expect(hint(page)).toBeVisible();
  await expect(hint(page)).toHaveText(SETTINGS_HINT, { useInnerText: true });
  expect(await hover(page, GC_SKEW)).toEqual({ cursor: 'pointer', outline: 'none' });
  await expect(hint(page)).toHaveText(SETTINGS_HINT, { useInnerText: true });
  // The Legend keeps its own hint.
  await hover(page, '#legend text');
  await expect(hint(page)).toHaveText(MOVE_HINT, { useInnerText: true });

  await clickInk(page, TICKS);
  await expectOpened(page, 'circular-axis', 'Axis & Scale');
  for (const [selector, label] of [[GC_SKEW, 'GC skew'], [GC_CONTENT, 'GC content'], ['g[data-gbdraw-slot-id="a_skew_2"] path', 'AT skew']]) {
    await closeSections(page);
    await clickInk(page, selector);
    await expectOpened(page, 'dinucleotide', 'Dinucleotide content/skew', label);
  }
  // Only the settings scroll area moved.
  expect(await page.evaluate(() => [document.querySelector('.app-main').scrollTop, document.querySelector('.settings-pane').scrollTop,
    document.scrollingElement.scrollTop])).toEqual([0, 0, 0]);

  // A press that pans the figure ends on the track it began on, and opens nothing.
  await closeSections(page);
  const pan = await page.evaluate(() => ({ ...window.__GBDRAW_APP__.canvasPan }));
  const grab = await inkPoint(page, GC_CONTENT);
  await page.mouse.move(grab.x, grab.y);
  await page.mouse.down();
  await page.mouse.move(grab.x + 48, grab.y + 30, { steps: 6 });
  // While it is panned, the track under the pointer shows the canvas's cursor and no hint.
  await frames(page);
  expect(await page.evaluate(({ x, y }) => {
    const hit = document.elementFromPoint(x, y);
    return [Boolean(hit.closest('[data-gbdraw-preview-item="settings"]')), getComputedStyle(hit).cursor];
  }, { x: grab.x + 48, y: grab.y + 30 })).toEqual([true, 'grabbing']);
  await expect(hint(page)).toBeHidden();
  await page.mouse.up();
  await expectClosed(page, 'dinucleotide', 'after a pan');
  expect(await page.evaluate(() => window.__GBDRAW_APP__.canvasPan.x)).toBeCloseTo(pan.x + 48, 0);
  // A modifier click keeps its meaning (feature selection).
  await clickInk(page, GC_CONTENT, 'Control');
  await clickInk(page, GC_CONTENT, 'Shift');
  await expectClosed(page, 'dinucleotide', 'after a modifier click');

  // A feature keeps its precedence, also when it is drawn inside a track's group.
  await page.evaluate(() => {
    const surface = document.querySelector('.gbdraw-preview-surface');
    const feature = surface.querySelector('path[data-gbdraw-feature-id]');
    feature.setAttribute('data-test-moved-feature', '');
    surface.querySelector('g[data-gbdraw-slot-id="gc_content"]').appendChild(feature);
  });
  await clickInk(page, '[data-test-moved-feature]');
  await expect(page.locator('.feature-popup').first()).toBeVisible();
  await expectClosed(page, 'dinucleotide', 'a feature inside a track group');
  // A click on a drawn element that closes an open popup does only that, as a click beside
  // the feature does: the settings stay where they are. The next click opens the section.
  const settled = await sectionState(page, 'dinucleotide');
  await clickInk(page, GC_SKEW);
  await expect(page.locator('.feature-popup')).toHaveCount(0);
  await expectClosed(page, 'dinucleotide', 'the click that closed the feature popup');
  expect(await sectionState(page, 'dinucleotide')).toEqual(settled);
  await clickInk(page, GC_SKEW);
  await expectOpened(page, 'dinucleotide', 'Dinucleotide content/skew', 'the click after the popup closed');
  // The same for a feature selection: the first click clears it.
  await closeSections(page);
  await page.keyboard.down('Control');
  await page.locator('.gbdraw-preview-surface path[data-gbdraw-feature-id]').nth(3).click({ force: true });
  await page.keyboard.up('Control');
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.selectedFeatureCount)).toBe(1);
  await clickInk(page, TICKS);
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.selectedFeatureCount)).toBe(0);
  await expectClosed(page, 'circular-axis', 'the click that cleared the selection');
  await clickInk(page, TICKS);
  await expectOpened(page, 'circular-axis', 'Axis & Scale', 'the click after the selection cleared');
  // And for the popup of a Legend row.
  await closeSections(page);
  const rowPopup = page.locator('.preview-item-popup[role="dialog"]');
  await clickInk(page, 'g[data-legend-key] text');
  await expect(rowPopup).toBeVisible();
  await clickInk(page, GC_SKEW);
  await expect(rowPopup).toHaveCount(0);
  await expectClosed(page, 'dinucleotide', 'the click that closed the popup of a Legend row');
  await clickInk(page, GC_SKEW);
  await expectOpened(page, 'dinucleotide', 'Dinucleotide content/skew', 'the click after the popup of a Legend row closed');

  // No setting, History step, or Result changed, and an export carries no marker.
  await page.waitForTimeout(150);
  expect(await appState(page)).toBe(before);
  expect(await exported(page)).not.toContain('data-gbdraw-preview-item');

  // Layout edit on: a click on a tick label or a track opens nothing and adds no History step, and the hint is gone.
  const steps = await undoCount(page);
  await layoutEdit(page).click();
  await expect(layoutEdit(page)).toHaveAttribute('aria-pressed', 'true');
  await closeSections(page);
  await clickInk(page, TICKS);
  await expectClosed(page, 'circular-axis', 'Layout edit on, a tick label');
  await clickInk(page, GC_SKEW);
  await expectClosed(page, 'dinucleotide', 'Layout edit on, a track');
  await expect(hint(page)).toHaveCount(0);
  expect(await undoCount(page)).toBe(steps);
});

test('Circular in the icon rail: the clicked section replaces the docked one, and docks beside the rail alone', async ({ page }) => {
  test.setTimeout(300_000);
  await page.setViewportSize({ width: 1600, height: 1000 });
  await openApp(page);
  await loadSession(page, CUSTOM_SLOTS_SESSION);
  await toRail(page).click();
  await entry(page, 'Layout').click();
  expect(await sectionState(page, 'layout')).toMatchObject({ display: 'rail', docked: 'layout' });

  await clickInk(page, TICKS);
  expect(await expectOpened(page, 'circular-axis', 'Axis & Scale', 'with Layout docked')).toMatchObject({ display: 'rail', docked: 'circular-axis', offset: 0 });
  expect((await sectionState(page, 'layout')).shown, 'Layout left the dock').toBe(false);
  await clickInk(page, GC_SKEW);
  expect(await expectOpened(page, 'dinucleotide', 'Dinucleotide content/skew', 'with Axis & Scale docked')).toMatchObject({ docked: 'dinucleotide', offset: 0 });

  // The rail alone: the click docks the section.
  await closeSection(page).click();
  expect(await sectionState(page, 'dinucleotide')).toMatchObject({ display: 'rail', docked: null, shown: false });
  await clickInk(page, GC_CONTENT);
  expect(await expectOpened(page, 'dinucleotide', 'Dinucleotide content/skew', 'the rail alone')).toMatchObject({ display: 'rail', docked: 'dinucleotide', offset: 0 });
  await closeSection(page).click();
  await clickInk(page, TICKS);
  expect(await expectOpened(page, 'circular-axis', 'Axis & Scale', 'the rail alone')).toMatchObject({ docked: 'circular-axis', offset: 0 });
});

test('Circular, the default layout with a plot title: the title opens Titles & Record Labels and keeps its Layout edit hint', async ({ page }) => {
  test.setTimeout(420_000);
  await page.setViewportSize({ width: 1600, height: 1000 });
  await openApp(page);
  await loadSession(page, DEFAULT_CIRCULAR_SESSION);
  const generated = await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    app.form.plot_title = 'Human mitochondrial genome';
    app.adv.plot_title_position = 'top';
    return (await app.runAnalysis()).status;
  });
  expect(generated).toBe('ok');
  await frames(page);
  const before = await appState(page);

  // The title is a Layout edit target too: its help cursor and frame stay, and the hint reads both.
  expect(await hover(page, '#plot_title text')).toEqual({ cursor: 'help', outline: 'dashed' });
  await expect(hint(page)).toHaveText(BOTH_HINTS, { useInnerText: true });
  await clickInk(page, '#plot_title text');
  await expectOpened(page, 'titles', 'Titles & Record Labels');

  // The groups of the default layout, by their own ids.
  for (const [selector, key, title] of [['#tick text', 'circular-axis', 'Axis & Scale'],
    ['#gc_content path', 'dinucleotide', 'Dinucleotide content/skew'], ['#skew path', 'dinucleotide', 'Dinucleotide content/skew']]) {
    await closeSections(page);
    await clickInk(page, selector);
    await expectOpened(page, key, title, selector);
  }
  await page.waitForTimeout(150);
  expect(await appState(page)).toBe(before);
  const svg = await exported(page);
  expect(svg).toContain('id="plot_title"');
  expect(svg).not.toContain('data-gbdraw-preview-item');
});

test('Linear: the scale opens Axis & Scale and the plot title Titles; a record definition and Layout edit open nothing', async ({ page }) => {
  test.setTimeout(420_000);
  await page.setViewportSize({ width: 1600, height: 1000 });
  await openApp(page);
  // The stored Result, with the definition of its first record named as that of a record called
  // `tick_chr1` is (`<record id>_definition…`): the id reads like a tick group's.
  const stored = readFileSync(sessionPath(LINEAR_SESSION), 'utf8');
  expect(stored.split('BGC0000708_definition_record_1')).toHaveLength(2);
  await loadSession(page, { name: LINEAR_SESSION, mimeType: 'application/json',
    buffer: Buffer.from(stored.replace('BGC0000708_definition_record_1', 'tick_chr1_definition_record_1')) });
  const before = await appState(page);

  // It stays the button of its popup, and only the scale and the title are marked for their settings.
  await expect(page.locator('.gbdraw-preview-surface [id="tick_chr1_definition_record_1"]')).toHaveAttribute('data-gbdraw-preview-item', '');
  expect((await page.locator(`.gbdraw-preview-surface ${MARKER}`).evaluateAll((elements) => elements.map((element) => element.id))).sort())
    .toEqual(['length_bar', 'plot_title']);

  for (const [selector, key, title] of [['#length_bar', 'linear-axis', 'Axis & Scale'], ['#plot_title', 'titles', 'Titles & Record Labels']]) {
    await closeSections(page);
    expect(await hover(page, `${selector} text`), selector).toEqual({ cursor: 'help', outline: 'dashed' });
    await expect(hint(page)).toHaveText(BOTH_HINTS, { useInnerText: true });
    expect(await hintLines(page), 'two lines').toBe(2);
    await clickInk(page, `${selector} text`);
    await expectOpened(page, key, title);
    // A held press pans the figure: the Layout edit sentence stays, alone, and says why the item did not move.
    const held = await inkPoint(page, `${selector} text`);
    await page.mouse.move(held.x, held.y);
    await page.mouse.down();
    await page.mouse.move(held.x + 12, held.y + 8, { steps: 4 });
    await expect(hint(page), `${selector}, held`).toHaveText(MOVE_HINT, { useInnerText: true });
    await page.mouse.up();
  }
  await hover(page, '#legend text');
  await expect(hint(page)).toHaveText(MOVE_HINT, { useInnerText: true });
  expect(await appState(page)).toBe(before);

  // The rail with another section docked.
  await toRail(page).click();
  await entry(page, 'Layout').click();
  await clickInk(page, '#length_bar text');
  expect(await expectOpened(page, 'linear-axis', 'Axis & Scale', 'the rail')).toMatchObject({ display: 'rail', docked: 'linear-axis', offset: 0 });
  await page.getByRole('button', { name: 'Expand settings', exact: true }).click();

  // Layout edit on: a click opens nothing, and a press still drags the title.
  await layoutEdit(page).click();
  await expect(layoutEdit(page)).toHaveAttribute('aria-pressed', 'true');
  await closeSections(page);
  await clickInk(page, '#plot_title text');
  await expectClosed(page, 'titles', 'Layout edit on, the title');
  await clickInk(page, '#length_bar text');
  await expectClosed(page, 'linear-axis', 'Layout edit on, the scale');
  await expect(hint(page)).toHaveCount(0);
  const steps = await undoCount(page);
  const titleAt = () => page.locator('.gbdraw-preview-surface #plot_title').getAttribute('transform');
  const transform = await titleAt();
  const grip = await inkPoint(page, '#plot_title text');
  await page.mouse.move(grip.x, grip.y);
  await page.mouse.down();
  await page.mouse.move(grip.x + 30, grip.y + 14, { steps: 6 });
  await page.mouse.up();
  await expect.poll(() => undoCount(page)).toBe(steps + 1);
  expect(await titleAt()).not.toBe(transform);
  await expectClosed(page, 'titles', 'after a drag of the title');
});

test('Linear with several records: a GC content or GC skew track opens Dinucleotide content/skew; Review alignment offers and opens nothing', async ({ page }) => {
  test.setTimeout(420_000);
  await page.setViewportSize({ width: 1600, height: 1000 });
  await openApp(page);
  await loadSession(page, LINEAR_SESSION);
  // No stored Linear Session draws these tracks: turn them on and generate.
  expect(await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    app.form.show_gc = true;
    app.form.show_skew = true;
    return (await app.runAnalysis()).status;
  })).toBe('ok');
  await frames(page);
  await page.getByRole('button', { name: 'Fit to preview', exact: true }).click();
  const tracks = await page.locator(`.gbdraw-preview-surface ${MARKER}`).evaluateAll((elements) => elements.map((element) => element.id));
  expect(tracks.filter((id) => /^gc_(content|skew)_record_\d+$/.test(id))).toHaveLength(10);
  for (const id of ['gc_content_record_2', 'gc_skew_record_4']) {
    await closeSections(page);
    await clickInk(page, `[id="${id}"] *`);
    await expectOpened(page, 'dinucleotide', 'Dinucleotide content/skew', id);
  }

  // Review alignment open, from the drawer and with no popup held: a track shows the canvas's
  // cursor and no hint, and a click on a drawn element opens nothing.
  await page.locator('[data-gbdraw-feature-id="f94b7f30b_record_2"]').click({ force: true });
  await page.locator('.feature-popup').getByRole('button', { name: 'Show in editor' }).click();
  await page.locator('[data-similarity-alignment-drawer-review]').click();
  const review = page.getByRole('dialog', { name: 'Review alignment' });
  await expect(review.locator('[data-alignment-record-key]')).toHaveCount(4, { timeout: 180000 });
  await expect(page.locator('.feature-popup')).toHaveCount(0);
  await closeSections(page);
  const SKEW_TRACK = '[id="gc_skew_record_4"] *';
  expect(await hover(page, SKEW_TRACK)).toEqual({ cursor: 'grab', outline: 'none' });
  await expect(hint(page)).toBeHidden();
  // What the review leaves as it is: the Layout edit sentence of the title, the pointer of a Legend row.
  expect(await hover(page, '#plot_title text')).toEqual({ cursor: 'help', outline: 'dashed' });
  await expect(hint(page)).toHaveText(MOVE_HINT, { useInnerText: true });
  expect((await hover(page, 'g[data-legend-key] text')).cursor).toBe('pointer');
  for (const [selector, key] of [[SKEW_TRACK, 'dinucleotide'], ['#length_bar text', 'linear-axis'], ['#plot_title text', 'titles']]) {
    await clickInk(page, selector);
    await expectClosed(page, key, `Review alignment open, ${selector}`);
  }
  await expect(review).toBeVisible();
});

test('phone width: a tick label brings Axis & Scale below the header', async ({ page }) => {
  test.setTimeout(300_000);
  await page.setViewportSize({ width: 390, height: 844 });
  await openApp(page);
  await loadSession(page, CUSTOM_SLOTS_SESSION);
  // The Result is below the settings there: the page is scrolled to it.
  await page.locator('.gbdraw-preview-surface').scrollIntoViewIfNeeded();
  await clickInk(page, TICKS);
  await expectOpened(page, 'circular-axis', 'Axis & Scale');
  expect(await page.evaluate(() => document.documentElement.scrollWidth - innerWidth)).toBeLessThanOrEqual(0);
});
