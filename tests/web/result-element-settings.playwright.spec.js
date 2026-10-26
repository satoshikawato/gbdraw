// A click on a drawn element of the Result shows the settings section that
// holds its settings in a popup beside the clicked point: the scale ticks and
// the Linear scale show Axis & Scale, a dinucleotide track shows Dinucleotide
// content/skew, and the plot title shows Titles & Record Labels. The popup
// holds the section itself, moved out of the pane, so its controls keep their
// ids and their History steps. It opens while Layout edit is off and no Review
// alignment is open; opening it writes no setting, History step, or Result,
// and leaves no trace in an export.
const { test, expect } = require('@playwright/test');
const { readFileSync } = require('node:fs');
const { join } = require('node:path');
const { evaluateWithRetainedPromise, generateAndWaitForResult, openApp } = require('./helpers/app-lifecycle.cjs');
const { showDrawing, showPhonePane } = require('./helpers/drawing-bar.cjs');

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
const AXIS = 'Axis & Scale';
const DINUCLEOTIDE = 'Dinucleotide content/skew';
const TITLES = 'Titles & Record Labels';

const index = (page) => page.getByRole('navigation', { name: 'Settings sections' });
const entry = (page, title) => index(page).getByRole('button', { name: `Go to ${title}`, exact: true });
const toRail = (page) => page.getByRole('button', { name: 'Collapse settings to icons', exact: true });
const popup = (page) => page.locator('#settings-popup');
const closePopup = (page) => popup(page).getByRole('button', { name: 'Close section', exact: true });
const closeDock = (page) => page.locator('.settings-pane').getByRole('button', { name: 'Close section', exact: true });
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
  await showPhonePane(page, 'preview');
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
// The box of the drawn element at a point of the window: what the popup of a click there keeps clear of.
const inkBox = (page, point) => page.evaluate(({ x, y }) => {
  const { left, top, right, bottom } = document.elementFromPoint(x, y).getBoundingClientRect();
  return { left, top, right, bottom };
}, point);
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
const appState = (page) => page.evaluate(() => {
  const app = window.__GBDRAW_APP__;
  return JSON.stringify([window.__GBDRAW_HISTORY__.getUndoCount(), app.form, app.adv, app.currentResults().map((result) => result.content)]);
});
const exported = (page) => page.evaluate(async () => {
  const { serializeCleanSvg } = await import('./js/services/svg-serialization.js');
  return serializeCleanSvg(window.__GBDRAW_APP__.svgContainer.querySelector('svg'));
});


// The popup and the section of `key`: what the popup holds, how many elements
// carry the key and how many of them the pane holds, whether the popup lies
// inside the Result's canvas and the window, and the ids that would name two
// elements or none.
const popupState = (page, key) => page.evaluate((key) => {
  const host = document.getElementById('settings-popup');
  const sections = Array.from(document.querySelectorAll(`[data-settings-section="${key}"]`));
  const details = sections[0]?.querySelector(':scope > details:first-child');
  const box = host.getBoundingClientRect();
  const canvas = document.querySelector('.preview-canvas').getBoundingClientRect();
  const ids = Array.from(document.querySelectorAll('[id]'), (element) => element.id);
  return {
    shown: host.checkVisibility(),
    inApp: Boolean(host.closest('#app')),
    held: Array.from(host.querySelectorAll('[data-settings-section]'), (element) => element.dataset.settingsSection),
    copies: sections.length,
    inPane: document.querySelectorAll(`.settings-scroll [data-settings-section="${key}"]`).length,
    open: details?.open ?? null,
    headingFocused: Boolean(details) && document.activeElement === details.querySelector(':scope > summary'),
    inside: box.left >= Math.max(0, canvas.left) && box.right <= Math.min(innerWidth, canvas.right)
      && box.top >= Math.max(0, canvas.top) && box.bottom <= Math.min(innerHeight, canvas.bottom),
    box: [box.left, box.top, box.width, box.height].map(Math.round),
    duplicateIds: ids.filter((id, order) => ids.indexOf(id) !== order),
    dangling: Array.from(document.querySelectorAll('[aria-describedby]'))
      .flatMap((element) => element.getAttribute('aria-describedby').split(/\s+/)).filter((id) => id && !document.getElementById(id)),
    settingsScroll: document.querySelector('.settings-scroll').scrollTop
  };
}, key);

// The popup shows the one section of `key`, with the focus on its heading, inside the Result.
const expectPopup = async (page, key, title, label = title) => {
  await expect.poll(async () => (await popupState(page, key)).headingFocused, { message: `${label}: focus on the popup's heading` }).toBe(true);
  await frames(page);
  const state = await popupState(page, key);
  expect(state, label).toMatchObject({ shown: true, inApp: true, held: [key], copies: 1, inPane: 0, open: true, inside: true });
  await expect(popup(page).locator('summary').first(), label).toContainText(title);
  return state;
};
// No popup: every section is in the pane.
const expectNoPopup = async (page, label) => {
  await frames(page);
  await page.waitForTimeout(150);
  await expect(popup(page), label).toBeHidden();
  expect(await popup(page).locator('[data-settings-section]').count(), label).toBe(0);
};
// A point of the Result's canvas that no drawn element and no popup covers.
const backgroundPoint = async (page) => {
  await figureSettled(page);
  const point = await page.evaluate(() => {
    const canvas = document.querySelector('.preview-canvas');
    const surface = canvas.querySelector('.gbdraw-preview-surface');
    const box = canvas.getBoundingClientRect();
    for (let row = 1; row < 40; row += 1) {
      for (let column = 1; column < 40; column += 1) {
        const x = box.left + (box.width * column) / 40;
        const y = box.top + (box.height * row) / 40;
        const hit = document.elementFromPoint(x, y);
        if (hit && canvas.contains(hit) && (hit.contains(surface) || hit === surface.querySelector('svg'))) return { x, y };
      }
    }
    return null;
  });
  expect(point, 'a point of the background').not.toBeNull();
  return point;
};

const TICKS = 'g[data-gbdraw-slot-renderer="ticks"] text';
const GC_CONTENT = 'g[data-gbdraw-slot-id="gc_content"] path';
const GC_SKEW = 'g[data-gbdraw-slot-id="gc_skew"] path';
const LEGEND_ROW = 'g[data-legend-key] text';
const featurePopup = (page) => page.locator('.feature-popup:not(.pairwise-match-popup)');
const rowPopup = (page) => page.locator('.preview-item-popup[role="dialog"]');

test('Circular with custom track slots: a tick label or a dinucleotide track shows its section in a popup, one popup of the Result at a time', async ({ page }) => {
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
  await hover(page, LEGEND_ROW);
  await expect(hint(page)).toHaveText(MOVE_HINT, { useInnerText: true });

  // In the pane, Dinucleotide content/skew is open and Axis & Scale closed: each returns as it was.
  await closeSections(page);
  await page.locator('.settings-scroll [data-settings-section="dinucleotide"] > details > summary').click();
  await page.evaluate(() => document.activeElement.blur());
  const rest = await popupState(page, 'circular-axis');
  expect(rest).toMatchObject({ shown: false, held: [], copies: 1, inPane: 1, open: false });

  const click = await clickInk(page, TICKS);
  const axis = await expectPopup(page, 'circular-axis', AXIS);
  // Beside the clicked point, as wide as the docked panel, and no id names two elements or none.
  const [left, top, width, height] = axis.box;
  expect(width).toBe(320);
  expect(click.x < left || click.x > left + width || click.y < top || click.y > top + height, 'the clicked point stays uncovered').toBe(true);
  expect([axis.duplicateIds, axis.dangling]).toEqual([rest.duplicateIds, rest.dangling]);
  // The pane neither scrolled nor took the focus.
  expect(axis.settingsScroll).toBe(rest.settingsScroll);
  await expect(closePopup(page)).toBeVisible();
  await expect(page.getByRole('dialog', { name: 'Settings of the clicked element' })).toBeVisible();

  // Another element swaps the section in the one popup; the first is back in the pane, closed as it was.
  for (const [selector, label] of [[GC_SKEW, 'GC skew'], [GC_CONTENT, 'GC content'], ['g[data-gbdraw-slot-id="a_skew_2"] path', 'AT skew']]) {
    await clickInk(page, selector);
    await expectPopup(page, 'dinucleotide', DINUCLEOTIDE, label);
  }
  expect(await popupState(page, 'circular-axis')).toMatchObject({ copies: 1, inPane: 1, open: false });

  // Tab walks the section's controls inside the popup; its heading drags it and does not fold it.
  await page.keyboard.press('Tab');
  expect(await page.evaluate(() => {
    const host = document.getElementById('settings-popup');
    return host.contains(document.activeElement) && !document.activeElement.matches('summary');
  })).toBe(true);
  const heading = popup(page).locator('[data-settings-section] > details > summary');
  await heading.click({ position: { x: 150, y: 14 } });
  await heading.press('Enter');
  const placed = await popupState(page, 'dinucleotide');
  expect(placed, 'the heading folds nothing').toMatchObject({ shown: true, held: ['dinucleotide'], open: true, inside: true });
  const grip = await heading.boundingBox();
  await page.mouse.move(grip.x + 150, grip.y + 14);
  await page.mouse.down();
  await page.mouse.move(grip.x + 110, grip.y + 44, { steps: 5 });
  await page.mouse.up();
  const dragged = await popupState(page, 'dinucleotide');
  expect([dragged.box[0] - placed.box[0], dragged.box[1] - placed.box[1]]).toEqual([-40, 30]);

  // A pan of the figure and a zoom leave it open where it is; a pan that begins on a track opens nothing else.
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
  expect(await page.evaluate(() => window.__GBDRAW_APP__.canvasPan.x)).toBeCloseTo(pan.x + 48, 0);
  const empty = await backgroundPoint(page);
  await page.mouse.move(empty.x, empty.y);
  await page.mouse.down();
  await page.mouse.move(empty.x + 20, empty.y + 12, { steps: 4 });
  await page.mouse.up();
  await page.mouse.wheel(0, -120);
  await frames(page);
  expect(await popupState(page, 'dinucleotide')).toMatchObject({ shown: true, held: ['dinucleotide'], box: dragged.box });
  await page.getByRole('button', { name: 'Fit to preview', exact: true }).click();

  // The close button, Escape, and a click on the background close it; the section is back, open as it was.
  await closePopup(page).click();
  await expectNoPopup(page, 'the close button');
  expect(await popupState(page, 'dinucleotide')).toMatchObject({ copies: 1, inPane: 1, open: true });
  await clickInk(page, TICKS);
  await expectPopup(page, 'circular-axis', AXIS);
  await page.keyboard.press('Escape');
  await expectNoPopup(page, 'Escape');
  expect(await popupState(page, 'circular-axis')).toMatchObject({ copies: 1, inPane: 1, open: false });
  await clickInk(page, TICKS);
  await expectPopup(page, 'circular-axis', AXIS);
  const background = await backgroundPoint(page);
  await page.mouse.click(background.x, background.y);
  await expectNoPopup(page, 'a click on the background');

  // A press that pans the figure ends on the track it began on, and a modifier click selects: neither opens it.
  const held = await inkPoint(page, GC_CONTENT);
  await page.mouse.move(held.x, held.y);
  await page.mouse.down();
  await page.mouse.move(held.x + 30, held.y + 20, { steps: 5 });
  await page.mouse.up();
  await expectNoPopup(page, 'after a pan');
  await clickInk(page, GC_CONTENT, 'Control');
  await clickInk(page, GC_CONTENT, 'Shift');
  await expectNoPopup(page, 'after a modifier click');

  // A feature keeps its precedence, also when it is drawn inside a track's group.
  await page.evaluate(() => {
    const surface = document.querySelector('.gbdraw-preview-surface');
    const feature = surface.querySelector('path[data-gbdraw-feature-id]');
    feature.setAttribute('data-test-moved-feature', '');
    surface.querySelector('g[data-gbdraw-slot-id="gc_content"]').appendChild(feature);
  });
  await clickInk(page, '[data-test-moved-feature]');
  await expect(featurePopup(page)).toBeVisible();
  await expectNoPopup(page, 'a feature inside a track group');
  // The rich feature popup is about half as wide as the Result: drag it by its heading into the
  // lower right corner, so that it leaves the Legend rows uncovered.
  const opened = await featurePopup(page).boundingBox();
  await page.mouse.move(opened.x + opened.width / 2, opened.y + 12);
  await page.mouse.down();
  await page.mouse.move(opened.x + opened.width / 2 + 2000, opened.y + 2000, { steps: 4 });
  await page.mouse.up();
  // The drag moved it and left it open.
  await expect(featurePopup(page)).toBeVisible();
  expect((await featurePopup(page).boundingBox()).x).toBeGreaterThan(opened.x);

  // One popup of the Result at a time, as a click on a Legend row does it: the clicked item's
  // popup takes the place of the open one in the same click.
  await clickInk(page, LEGEND_ROW);
  await expect(rowPopup(page)).toBeVisible();
  await expect(featurePopup(page)).toHaveCount(0);
  await clickInk(page, GC_SKEW);
  await expectPopup(page, 'dinucleotide', DINUCLEOTIDE, 'in place of the popup of a Legend row');
  await expect(rowPopup(page)).toHaveCount(0);
  await clickInk(page, '[data-test-moved-feature]');
  await expect(featurePopup(page)).toBeVisible();
  await expectNoPopup(page, 'a feature popup in its place');
  await clickInk(page, TICKS);
  await expectPopup(page, 'circular-axis', AXIS, 'in place of the feature popup');
  await expect(featurePopup(page)).toHaveCount(0);
  await clickInk(page, LEGEND_ROW);
  await expect(rowPopup(page)).toBeVisible();
  await expectNoPopup(page, 'the popup of a Legend row in its place');
  await page.keyboard.press('Escape');
  // A feature selection stays, as it does under the popup of a Legend row.
  await page.keyboard.down('Control');
  await page.locator('.gbdraw-preview-surface path[data-gbdraw-feature-id]').nth(3).click({ force: true });
  await page.keyboard.up('Control');
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.selectedFeatureCount)).toBe(1);
  await clickInk(page, TICKS);
  await expectPopup(page, 'circular-axis', AXIS, 'with a feature selected');
  expect(await page.evaluate(() => window.__GBDRAW_APP__.selectedFeatureCount)).toBe(1);

  // The index entry of the floating section takes it back into the pane and opens it there.
  await entry(page, AXIS).click();
  await expectNoPopup(page, 'the index entry of the floating section');
  expect(await expectOpened(page, 'circular-axis', AXIS)).toMatchObject({ display: 'stacked' });

  // No setting, History step, or Result changed, and an export carries no marker.
  await page.waitForTimeout(150);
  expect(await appState(page)).toBe(before);
  expect(await exported(page)).not.toContain('data-gbdraw-preview-item');

  // Layout edit: turning it on closes the popup, and a click on a tick label or a track then opens nothing.
  await clickInk(page, GC_SKEW);
  await expectPopup(page, 'dinucleotide', DINUCLEOTIDE);
  const steps = await undoCount(page);
  await layoutEdit(page).click();
  await expect(layoutEdit(page)).toHaveAttribute('aria-pressed', 'true');
  await expectNoPopup(page, 'Layout edit turned on');
  await clickInk(page, TICKS);
  await clickInk(page, GC_SKEW);
  await expectNoPopup(page, 'Layout edit on');
  await expect(hint(page)).toHaveCount(0);
  expect(await undoCount(page)).toBe(steps);
});

test('Circular in the icon rail: the docked section takes the focus where it is, another section floats beside the rail', async ({ page }) => {
  test.setTimeout(300_000);
  await page.setViewportSize({ width: 1600, height: 1000 });
  await openApp(page);
  await loadSession(page, CUSTOM_SLOTS_SESSION);
  await toRail(page).click();
  await entry(page, 'Layout').click();
  expect(await sectionState(page, 'layout')).toMatchObject({ display: 'rail', docked: 'layout' });

  // Another section is docked: the clicked one floats, and the dock keeps its section.
  await clickInk(page, TICKS);
  await expectPopup(page, 'circular-axis', AXIS, 'with Layout docked');
  expect(await sectionState(page, 'layout')).toMatchObject({ display: 'rail', docked: 'layout', shown: true });
  // Its rail entry docks it: the popup closes.
  await entry(page, AXIS).click();
  await expectNoPopup(page, 'the rail entry of the floating section');
  expect(await expectOpened(page, 'circular-axis', AXIS, 'docked')).toMatchObject({ display: 'rail', docked: 'circular-axis', offset: 0 });

  // The docked section's own element: no popup, the focus on its heading.
  await page.evaluate(() => document.activeElement.blur());
  await clickInk(page, TICKS);
  await expectOpened(page, 'circular-axis', AXIS, 'a click on the docked section\'s element');
  await expectNoPopup(page, 'the docked section');
  await clickInk(page, GC_SKEW);
  await expectPopup(page, 'dinucleotide', DINUCLEOTIDE, 'with Axis & Scale docked');
  expect(await sectionState(page, 'circular-axis')).toMatchObject({ docked: 'circular-axis', shown: true });
  // A click on the docked section's element closes the popup of the other.
  await clickInk(page, TICKS);
  await expectNoPopup(page, 'back to the docked section');
  await expectOpened(page, 'circular-axis', AXIS, 'the docked section again');

  // The rail alone: the section floats, and the rail stays alone.
  await closeDock(page).click();
  await clickInk(page, GC_CONTENT);
  await expectPopup(page, 'dinucleotide', DINUCLEOTIDE, 'the rail alone');
  expect(await page.locator('.settings-pane').evaluate((pane) => [pane.dataset.settingsDisplay, pane.dataset.settingsDocked ?? null])).toEqual(['rail', null]);
});

test('Circular, the default layout with a plot title: an edit in the popup is one History step, and Generate keeps the popup open', async ({ page }) => {
  test.setTimeout(480_000);
  await page.setViewportSize({ width: 1600, height: 1000 });
  await openApp(page);
  await loadSession(page, DEFAULT_CIRCULAR_SESSION);
  const generated = await evaluateWithRetainedPromise(page, async () => {
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
  // The groups of the default layout, by their own ids.
  for (const [selector, key, title] of [['#tick text', 'circular-axis', AXIS], ['#gc_content path', 'dinucleotide', DINUCLEOTIDE],
    ['#skew path', 'dinucleotide', DINUCLEOTIDE], ['#plot_title text', 'titles', TITLES]]) {
    await clickInk(page, selector);
    await expectPopup(page, key, title, selector);
  }
  await page.waitForTimeout(150);
  expect(await appState(page)).toBe(before);
  const svg = await exported(page);
  expect(svg).toContain('id="plot_title"');
  expect(svg).not.toContain('data-gbdraw-preview-item');

  // An edit in the popup writes the draft as one History step; Undo and Redo restore it and leave the popup open.
  const steps = await undoCount(page);
  const title = popup(page).getByRole('textbox', { name: 'Plot Title', exact: true });
  await title.fill('Edited in the popup');
  await title.press('Tab');
  const draft = () => page.evaluate(() => window.__GBDRAW_APP__.form.plot_title);
  expect(await draft()).toBe('Edited in the popup');
  await expect.poll(() => undoCount(page)).toBe(steps + 1);
  await page.getByRole('button', { name: 'Undo', exact: true }).click();
  await expect.poll(draft).toBe('Human mitochondrial genome');
  expect(await undoCount(page)).toBe(steps);
  expect(await popupState(page, 'titles')).toMatchObject({ shown: true, held: ['titles'], copies: 1, inPane: 0 });
  await page.getByRole('button', { name: 'Redo', exact: true }).click();
  await expect.poll(draft).toBe('Edited in the popup');

  // Generate keeps the popup open, with its section, where it is; the figure then shows the edit.
  const open = await popupState(page, 'titles');
  await page.getByRole('button', { name: 'Generate Diagram', exact: true }).click();
  await expect(page.locator('.gbdraw-preview-surface #plot_title')).toContainText('Edited in the popup', { timeout: 180_000 });
  await page.waitForFunction(() => !window.__GBDRAW_APP__.processing);
  expect(await popupState(page, 'titles')).toMatchObject({ shown: true, held: ['titles'], copies: 1, inPane: 0, open: true, box: open.box });
  await expect(title).toHaveValue('Edited in the popup');

  // Reset Settings closes it.
  await page.getByRole('button', { name: 'Reset Settings', exact: true }).click();
  await page.getByRole('button', { name: 'Reset all settings', exact: true }).click();
  await expectNoPopup(page, 'Reset Settings');
  expect(await popupState(page, 'titles')).toMatchObject({ copies: 1, inPane: 1 });
});

test('Linear at 1280 x 720: the scale shows Axis & Scale and the plot title Titles; a record definition, a mode switch, a Session, and Layout edit', async ({ page }) => {
  test.setTimeout(420_000);
  await page.setViewportSize({ width: 1280, height: 720 });
  await openApp(page);
  // The stored Result, with the definition of its first record named as that of a record called
  // `tick_chr1` is (`<record id>_definition…`): the id reads like a tick group's.
  const stored = readFileSync(sessionPath(LINEAR_SESSION), 'utf8');
  expect(stored.split('BGC0000708_definition_record_1')).toHaveLength(2);
  const session = { name: LINEAR_SESSION, mimeType: 'application/json',
    buffer: Buffer.from(stored.replace('BGC0000708_definition_record_1', 'tick_chr1_definition_record_1')) };
  await loadSession(page, session);
  const before = await appState(page);

  // It stays the button of its popup, and only the scale and the title are marked for their settings.
  const definition = page.locator('.gbdraw-preview-surface [id="tick_chr1_definition_record_1"]');
  await expect(definition).toHaveAttribute('data-gbdraw-preview-item', '');
  expect((await page.locator(`.gbdraw-preview-surface ${MARKER}`).evaluateAll((elements) => elements.map((element) => element.id))).sort())
    .toEqual(['length_bar', 'plot_title']);

  for (const [selector, key, title] of [['#length_bar', 'linear-axis', AXIS], ['#plot_title', 'titles', TITLES]]) {
    expect(await hover(page, `${selector} text`), selector).toEqual({ cursor: 'help', outline: 'dashed' });
    await expect(hint(page)).toHaveText(BOTH_HINTS, { useInnerText: true });
    expect(await hintLines(page), 'two lines').toBe(2);
    await clickInk(page, `${selector} text`);
    await expectPopup(page, key, title);
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

  // The record definition opens its own popup, in place of the settings popup.
  await clickInk(page, '[id="tick_chr1_definition_record_1"] text');
  await expect(rowPopup(page)).toBeVisible();
  await expectNoPopup(page, 'a record definition');
  await page.keyboard.press('Escape');

  // A mode switch closes the popup: Circular has its own Axis & Scale section.
  await clickInk(page, '#length_bar text');
  await expectPopup(page, 'linear-axis', AXIS);
  await showDrawing(page, { kind: 'circular' });
  await expectNoPopup(page, 'a mode switch');
  await showDrawing(page, { kind: 'linear' });
  await page.getByRole('button', { name: 'Fit to preview', exact: true }).click();
  expect(await popupState(page, 'linear-axis')).toMatchObject({ copies: 1, inPane: 1 });
  // A loaded Session closes it.
  await clickInk(page, '#plot_title text');
  await expectPopup(page, 'titles', TITLES);
  await loadSession(page, session);
  await expectNoPopup(page, 'a loaded Session');

  // Layout edit on: a click opens nothing, and a press still drags the title.
  await layoutEdit(page).click();
  await expect(layoutEdit(page)).toHaveAttribute('aria-pressed', 'true');
  await clickInk(page, '#plot_title text');
  await clickInk(page, '#length_bar text');
  await expectNoPopup(page, 'Layout edit on');
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
  await expectNoPopup(page, 'after a drag of the title');
});

test('Linear with several records: a GC content or GC skew track shows Dinucleotide content/skew; Review alignment closes the popup and offers nothing', async ({ page }) => {
  test.setTimeout(420_000);
  await page.setViewportSize({ width: 1600, height: 1000 });
  await openApp(page);
  await loadSession(page, LINEAR_SESSION);
  // No stored Linear Session draws these tracks: turn them on and generate.
  expect(await evaluateWithRetainedPromise(page, async () => {
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
    const point = await inkPoint(page, `[id="${id}"] *`);
    const track = await inkBox(page, point);
    await page.mouse.click(point.x, point.y);
    const { box: [left, top, width, height] } = await expectPopup(page, 'dinucleotide', DINUCLEOTIDE, id);
    // A track is as wide as the Result and flat: no room beside it, and room below or above it.
    expect(left >= track.right || left + width <= track.left || top >= track.bottom || top + height <= track.top,
      `${id}: the popup stands clear of the track`).toBe(true);
  }

  // Review alignment, opened from the drawer while the popup shows, closes it. The popup of the
  // last track covers the feature that opens the drawer: close it first.
  await closePopup(page).click();
  await expectNoPopup(page, 'Close section');
  await clickInk(page, '[data-gbdraw-feature-id="f94b7f30b_record_2"]');
  await page.locator('.feature-popup').getByRole('button', { name: 'Show in editor' }).click();
  await clickInk(page, '[id="gc_content_record_2"] *');
  await expectPopup(page, 'dinucleotide', DINUCLEOTIDE, 'beside the drawer');
  // Escape closes the top layer alone: the popup goes, and the drawer stays open.
  await page.keyboard.press('Escape');
  await expectNoPopup(page, 'Escape beside the drawer');
  await expect(page.locator('.right-drawer')).toBeVisible();
  await clickInk(page, '[id="gc_content_record_2"] *');
  await expectPopup(page, 'dinucleotide', DINUCLEOTIDE, 'beside the drawer again');
  await page.locator('[data-similarity-alignment-drawer-review]').click();
  const review = page.getByRole('dialog', { name: 'Review alignment' });
  await expect(review.locator('[data-alignment-record-key]')).toHaveCount(4, { timeout: 180000 });
  await expectNoPopup(page, 'Review alignment opened');
  await expect(page.locator('.feature-popup')).toHaveCount(0);
  // While it is open, a track shows the canvas's cursor and no hint, and a click on a drawn element opens nothing.
  const SKEW_TRACK = '[id="gc_skew_record_4"] *';
  expect(await hover(page, SKEW_TRACK)).toEqual({ cursor: 'grab', outline: 'none' });
  await expect(hint(page)).toBeHidden();
  // What the review leaves as it is: the Layout edit sentence of the title, the pointer of a Legend row.
  expect(await hover(page, '#plot_title text')).toEqual({ cursor: 'help', outline: 'dashed' });
  await expect(hint(page)).toHaveText(MOVE_HINT, { useInnerText: true });
  expect((await hover(page, LEGEND_ROW)).cursor).toBe('pointer');
  for (const selector of [SKEW_TRACK, '#length_bar text', '#plot_title text']) await clickInk(page, selector);
  await expectNoPopup(page, 'Review alignment open');
  await expect(review).toBeVisible();
});

test('Circular, a first Result: the popup gives the focus back, follows a window that changes size, and closes when the Result is not shown', async ({ page }) => {
  test.setTimeout(420_000);
  await page.setViewportSize({ width: 1600, height: 1000 });
  await openApp(page);
  await page.getByLabel('GenBank/DDBJ File', { exact: true }).setInputFiles('tests/test_inputs/HmmtDNA.gbk');
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.circularRecordList.length), { timeout: 60_000 }).toBeGreaterThan(0);
  await generateAndWaitForResult(page);
  await frames(page);
  await page.getByRole('button', { name: 'Fit to preview', exact: true }).click();
  const open = async (label) => {
    await clickInk(page, TICKS);
    return expectPopup(page, 'circular-axis', AXIS, label);
  };
  const inPane = () => page.locator('.settings-scroll [data-settings-section="circular-axis"]').count();

  // The focus the popup took goes back to what had it before, here a control of the pane.
  const generate = page.getByRole('button', { name: 'Generate Diagram', exact: true });
  await generate.focus();
  await open('from the Generate button');
  await closePopup(page).click();
  await expectNoPopup(page, 'the close button');
  await expect(generate).toBeFocused();
  // By the keyboard: the close button precedes the heading, so Shift+Tab reaches it from there.
  await open('again');
  await page.keyboard.press('Shift+Tab');
  await expect(closePopup(page)).toBeFocused();
  await page.keyboard.press('Enter');
  await expectNoPopup(page, 'Enter on the close button');
  await expect(generate).toBeFocused();

  // A window that gets smaller moves the popup the least distance back inside what a drag allows; it stays open.
  await open('before the window changes');
  const heading = popup(page).locator('[data-settings-section] > details > summary');
  const grip = await heading.boundingBox();
  await page.mouse.move(grip.x + 150, grip.y + 14);
  await page.mouse.down();
  await page.mouse.move(1590, 990, { steps: 6 });
  await page.mouse.up();
  const dragged = await popupState(page, 'circular-axis');
  expect(dragged.box[0] + dragged.box[2]).toBe(1600 - 12);
  await page.setViewportSize({ width: 1000, height: 640 });
  await frames(page);
  const bounds = await page.evaluate(async () => (await import('./js/utils/popup-bounds.js')).popupViewportBounds());
  const moved = await popupState(page, 'circular-axis');
  expect(moved).toMatchObject({ shown: true, held: ['circular-axis'], copies: 1, inPane: 0 });
  const [left, top, width, height] = moved.box;
  expect([left + width, top + height], 'at the right and bottom edges a drag allows').toEqual([1000 - 12, Math.round(bounds.bottom) - 12]);
  expect(top).toBeGreaterThanOrEqual(Math.round(bounds.top) + 12);
  await page.setViewportSize({ width: 1600, height: 1000 });
  await frames(page);
  expect((await popupState(page, 'circular-axis')).box.slice(0, 2), 'a larger window moves nothing back').toEqual([left, top]);

  // The popup belongs to the shown Result: another tab of the Result panel closes it, and its section is in the pane again.
  await page.getByRole('button', { name: 'Run info', exact: true }).click();
  await expectNoPopup(page, 'the Run info tab');
  expect(await inPane()).toBe(1);
  await page.getByRole('button', { name: 'Preview', exact: true }).click();
  await page.getByRole('button', { name: 'Fit to preview', exact: true }).click();
  // So does an Undo back to before the first Generate, which leaves no Result.
  await open('before the Undo');
  await page.getByRole('button', { name: 'Undo', exact: true }).click();
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.results.length)).toBe(0);
  await expect(popup(page), 'no Result').toBeHidden();
  expect(await popup(page).locator('[data-settings-section]').count()).toBe(0);
  expect(await inPane()).toBe(1);
});

test('phone width: a tick label shows Axis & Scale in a popup inside the Result', async ({ page }) => {
  test.setTimeout(300_000);
  await page.setViewportSize({ width: 390, height: 844 });
  await openApp(page);
  await loadSession(page, CUSTOM_SLOTS_SESSION);
  // The Result is below the settings there: the page is scrolled to it.
  await page.locator('.gbdraw-preview-surface').scrollIntoViewIfNeeded();
  await clickInk(page, TICKS);
  const state = await expectPopup(page, 'circular-axis', AXIS);
  expect(state.box[2]).toBeLessThanOrEqual(390 - 32);
  expect(await page.evaluate(() => document.documentElement.scrollWidth - innerWidth)).toBeLessThanOrEqual(0);
  await closePopup(page).click();
  await expectNoPopup(page, 'the close button');
});
