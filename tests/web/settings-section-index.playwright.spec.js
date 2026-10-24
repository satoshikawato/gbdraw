// The section index of the settings pane: one pinned row of icon buttons above
// the settings, one per section of the shown mode, in card order, with the
// switch to the icon rail at its right end. An entry
// opens its section, brings it to the top of the settings, and moves the
// keyboard focus to its heading; the entry of the section at the top is
// marked; the section title shows on hover and on keyboard focus. The index
// changes only what the pane shows: no setting, History step, or Result.
const { test, expect } = require('@playwright/test');
const { join } = require('node:path');
const { openApp, openLinearFileCards } = require('./helpers/app-lifecycle.cjs');
const { showDrawing } = require('./helpers/drawing-bar.cjs');

// Section key and title of every entry, per mode, in card order. Basic and
// About & Citation have no entry.
const SECTIONS = {
  circular: [
    ['input-genomes', 'Input Genomes'],
    ['region-annotations', 'Region Annotations'],
    ['pairwise-comparisons', 'Pairwise Comparisons'],
    ['layout', 'Layout'],
    ['font', 'Font'],
    ['titles', 'Titles & Record Labels'],
    ['legend', 'Legend'],
    ['colors', 'Colors'],
    ['features', 'Features'],
    ['labels', 'Labels'],
    ['circular-axis', 'Axis & Scale'],
    ['dinucleotide', 'Dinucleotide content/skew']
  ],
  linear: [
    ['input-genomes', 'Input Genomes'],
    ['comparison', 'Comparison'],
    ['advanced-comparison', 'Advanced comparison and layout'],
    ['region-annotations', 'Region Annotations'],
    ['layout', 'Layout'],
    ['font', 'Font'],
    ['titles', 'Titles & Record Labels'],
    ['legend', 'Legend'],
    ['colors', 'Colors'],
    ['features', 'Features'],
    ['labels', 'Labels'],
    ['linear-axis', 'Axis & Scale'],
    ['dinucleotide', 'Dinucleotide content/skew']
  ]
};

// An entry is named for what it does; the section's own disclosure keeps the bare title.
const entryName = (title) => `Go to ${title}`;
const index = (page) => page.getByRole('navigation', { name: 'Settings sections' });
const entry = (page, title) => index(page).getByRole('button', { name: entryName(title), exact: true });
const tip = (page) => page.locator('.settings-index-tip');
const setMode = (page, mode) => showDrawing(page, { kind: mode });

const loadSession = async (page, name) => {
  page.on('dialog', (dialog) => dialog.accept());
  await page.locator('input[accept^=".json,"]').setInputFiles(
    join(process.cwd(), 'gbdraw/web/gallery/sessions', name));
  await page.waitForFunction(() => !window.__GBDRAW_APP__.sessionImportPending
    && window.__GBDRAW_APP__.results.length > 0, null, { timeout: 180000 });
};

// Where one section is: its distance below the top edge of the visible
// settings, its disclosure, the focus, the marked entries, and every scroller
// outside the settings scroll area.
const sectionState = (page, key) => page.evaluate((key) => {
  const area = document.querySelector('.settings-scroll');
  const section = area.querySelector(`[data-settings-section="${key}"]`);
  // A card that is one disclosure is headed by its summary; Input Genomes and
  // Comparison (which holds disclosures of its own) have a heading that takes the focus.
  const details = section.querySelector(':scope > details:first-child');
  const heading = details ? details.querySelector(':scope > summary') : section.querySelector(':scope > [tabindex="-1"]');
  return {
    open: details ? details.open : null,
    offset: Math.round(section.getBoundingClientRect().top - area.getBoundingClientRect().top),
    atEnd: area.scrollTop + area.clientHeight >= area.scrollHeight - 1,
    headingFocused: document.activeElement === heading,
    focusInside: section.contains(document.activeElement) && document.activeElement !== heading
      && Boolean(heading.compareDocumentPosition(document.activeElement) & Node.DOCUMENT_POSITION_FOLLOWING),
    marked: Array.from(document.querySelectorAll('nav[aria-label="Settings sections"] [aria-current]'),
      (button) => button.getAttribute('aria-label')),
    outside: [document.querySelector('.app-main').scrollTop, document.querySelector('.settings-pane').scrollTop,
      document.scrollingElement.scrollTop, window.scrollY]
  };
}, key);

const marked = (page) => index(page).locator('[aria-current="true"]');
const expectMarked = (page, title) => expect(marked(page)).toHaveAttribute('aria-label', entryName(title));

test('the index lists the sections of each mode in card order, before the settings in Tab order', async ({ page }) => {
  test.setTimeout(180_000);
  await page.setViewportSize({ width: 1280, height: 800 });
  await openApp(page);
  for (const mode of ['circular', 'linear']) {
    await setMode(page, mode);
    const found = await page.evaluate(() => {
      const pane = document.querySelector('.settings-pane');
      const nav = pane.querySelector('nav[aria-label="Settings sections"]');
      const area = pane.querySelector('.settings-scroll');
      const iconOf = (element) => Array.from(element?.classList || []).filter((name) => name.startsWith('ph-')).join(' ');
      const visible = (element) => element.checkVisibility({ visibilityProperty: true, opacityProperty: false });
      const stops = Array.from(pane.querySelectorAll('a[href], button, input, select, textarea, summary, [tabindex]'))
        .filter((element) => visible(element) && !element.disabled && element.tabIndex >= 0);
      const buttons = Array.from(nav?.querySelectorAll('button') || []);
      return {
        // The row holds the list of entries and, after it, the switch to the rail.
        first: pane.firstElementChild === nav.parentElement && nav.parentElement.nextElementSibling === area,
        historyIgnore: Boolean(nav.closest('[data-history-ignore]')),
        controls: buttons.every((button) => button.getAttribute('aria-controls') === area.id && area.id),
        firstStops: stops.slice(0, buttons.length).every((stop, order) => stop === buttons[order]),
        entries: buttons.map((button) => ({
          title: button.getAttribute('aria-label'),
          type: button.type,
          nativeTitle: button.getAttribute('title'),
          text: button.textContent.trim(),
          icon: iconOf(button.querySelector('i')),
          iconHidden: button.querySelector('i')?.getAttribute('aria-hidden')
        })),
        sections: Array.from(area.querySelectorAll('[data-settings-section]')).filter(visible).map((section) => {
          const heading = section.querySelector(':scope > details:first-child > summary, :scope > .card-header, :scope > [id$="-heading"]');
          return {
            key: section.dataset.settingsSection,
            heading: heading.textContent.replace(/\s+/g, ' ').trim(),
            icon: iconOf(heading.querySelector(':scope > i')),
            focusable: heading.matches('summary') || heading.tabIndex === -1
          };
        }),
        marked: buttons.filter((button) => button.hasAttribute('aria-current'))
          .map((button) => `${button.getAttribute('aria-label')}=${button.getAttribute('aria-current')}`)
      };
    });
    const expected = SECTIONS[mode];
    expect(found.first, `${mode}: the index is the first child of the pane, before the scroll area`).toBe(true);
    expect(found.historyIgnore, `${mode}: a click on an entry begins no History step`).toBe(true);
    expect(found.controls, `${mode}: each entry names the scroll area it moves`).toBeTruthy();
    expect(found.firstStops, `${mode}: the entries are the first Tab stops of the pane`).toBe(true);
    expect(found.entries.map((item) => item.title), mode).toEqual(expected.map(([, title]) => entryName(title)));
    // Two cards have no row entry: Basic and About & Citation (the rail reaches them).
    const sections = found.sections.filter((section) => !['basic', 'about'].includes(section.key));
    expect(found.sections.length - sections.length, `${mode}: the cards without a row entry`).toBe(2);
    expect(sections.map((section) => section.key), `${mode}: one section per entry, same order`)
      .toEqual(expected.map(([key]) => key));
    for (const [order, item] of found.entries.entries()) {
      const section = sections[order];
      // The name is the aria-label alone: no title attribute, no text, a decorative icon.
      expect(item, `${mode} ${item.title}`).toMatchObject({ type: 'button', nativeTitle: null, text: '', iconHidden: 'true' });
      expect(item.icon, `${mode} ${item.title}: an icon`).toMatch(/^ph-[a-z-]+$/);
      // The card heading shows the same title and icon as the entry.
      expect(section.icon, `${mode} ${item.title}: heading icon`).toBe(item.icon);
      const title = expected[order][1];
      expect(section.heading.startsWith(title), `${mode}: "${section.heading}" starts with "${title}"`).toBe(true);
      expect(section.focusable, `${mode} ${item.title}: heading can take focus`).toBe(true);
    }
    expect(found.marked, `${mode}: Input Genomes is marked before any scroll`).toEqual([`${entryName('Input Genomes')}=true`]);
  }
});

test('an entry opens its section, brings it to the top, focuses its heading, and is marked (Linear, 5 records)', async ({ page }) => {
  test.setTimeout(300_000);
  await page.setViewportSize({ width: 1600, height: 1000 });
  await openApp(page);
  await loadSession(page, 'BGC0000708-BGC0000713.gbdraw-session.json');
  // The test measures scroll distances over Input Genomes with its cards open.
  await openLinearFileCards(page);
  await expect(index(page)).toBeVisible();
  const snapshot = () => page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    return JSON.stringify([app.canUndoHistory, app.form, app.adv, app.results.length]);
  });
  const before = await snapshot();
  const pinnedAt = await index(page).evaluate((nav) => Math.round(nav.getBoundingClientRect().top));

  // Pointer: Labels starts closed and far below. Opened, it is shorter than the
  // settings view at this height, so it comes up as far as the settings scroll.
  const start = await sectionState(page, 'labels');
  expect(start.open).toBe(false);
  expect(start.offset).toBeGreaterThan(1500);
  await entry(page, 'Labels').click();
  const labels = await sectionState(page, 'labels');
  expect(labels, 'at once, without smooth scrolling')
    .toMatchObject({ open: true, headingFocused: true, marked: [entryName('Labels')], outside: [0, 0, 0, 0] });
  expect(Math.abs(labels.offset) <= 2 || labels.atEnd, `at the top edge or at the end of the scroll (offset ${labels.offset})`).toBe(true);
  await expect(marked(page)).toHaveCount(1);
  // The row stays where it was: only the settings scroll area moved.
  expect(await index(page).evaluate((nav) => Math.round(nav.getBoundingClientRect().top))).toBe(pinnedAt);

  // An entry never closes a section. In a lower window the same section reaches the top edge.
  await page.setViewportSize({ width: 1600, height: 640 });
  await entry(page, 'Labels').click();
  expect(await sectionState(page, 'labels')).toMatchObject({ open: true, headingFocused: true, marked: [entryName('Labels')], outside: [0, 0, 0, 0] });
  expect(Math.abs((await sectionState(page, 'labels')).offset)).toBeLessThanOrEqual(2);
  await page.setViewportSize({ width: 1600, height: 1000 });

  // Keyboard: Enter on Comparison, a section without a summary; the next Tab goes into it.
  await entry(page, 'Comparison').focus();
  await page.keyboard.press('Enter');
  const comparison = await sectionState(page, 'comparison');
  expect(comparison).toMatchObject({ open: null, headingFocused: true, marked: [entryName('Comparison')], outside: [0, 0, 0, 0] });
  expect(Math.abs(comparison.offset)).toBeLessThanOrEqual(2);
  await page.keyboard.press('Tab');
  expect((await sectionState(page, 'comparison')).focusInside, 'Tab continues inside Comparison').toBe(true);

  // Space on a summary section; Tab continues inside it too.
  await entry(page, 'Features').focus();
  await page.keyboard.press('Space');
  const features = await sectionState(page, 'features');
  expect(features).toMatchObject({ open: true, headingFocused: true, marked: [entryName('Features')], outside: [0, 0, 0, 0] });
  expect(Math.abs(features.offset)).toBeLessThanOrEqual(2);
  await page.keyboard.press('Tab');
  expect((await sectionState(page, 'features')).focusInside, 'Tab continues inside Features').toBe(true);

  // Enter on a closed summary section leaves it open (the key does not toggle the focused summary).
  await entry(page, 'Colors').focus();
  await page.keyboard.press('Enter');
  expect(await sectionState(page, 'colors')).toMatchObject({ open: true, headingFocused: true, marked: [entryName('Colors')] });

  // Scrolling marks the section at the top edge; the top of the settings is Input Genomes.
  const scrollTo = (key) => page.evaluate((key) => {
    const area = document.querySelector('.settings-scroll');
    const section = key ? area.querySelector(`[data-settings-section="${key}"]`) : null;
    area.scrollTop = section ? area.scrollTop + section.getBoundingClientRect().top - area.getBoundingClientRect().top + 30 : 0;
  }, key);
  await scrollTo('colors');
  await expectMarked(page, 'Colors');
  await scrollTo('layout');
  await expectMarked(page, 'Layout');
  await scrollTo(null);
  await expectMarked(page, 'Input Genomes');
  await expect(marked(page)).toHaveCount(1);

  // The last section is too short to reach the top: its entry stays marked
  // until the scroll position changes again.
  await page.evaluate(() => document.querySelectorAll('.settings-scroll > .card > details')
    .forEach((details) => { details.open = false; }));
  await entry(page, 'Dinucleotide content/skew').click();
  const last = await sectionState(page, 'dinucleotide');
  expect(last).toMatchObject({ open: true, headingFocused: true, outside: [0, 0, 0, 0] });
  expect(last.offset, 'the last section cannot reach the top').toBeGreaterThan(20);
  // Two frames: the scroll event of the activation has been handled.
  await page.evaluate(() => new Promise((resolve) => requestAnimationFrame(() => requestAnimationFrame(resolve))));
  await expectMarked(page, 'Dinucleotide content/skew');
  await page.evaluate(() => { document.querySelector('.settings-scroll').scrollTop -= 10; });
  await expect(marked(page)).not.toHaveAttribute('aria-label', entryName('Dinucleotide content/skew'));
  await expect(marked(page)).toHaveCount(1);

  // No setting, History step, or Result changed (a button click commits on a timer).
  await page.waitForTimeout(100);
  expect(await snapshot()).toBe(before);

  // A mode without the marked section marks the section at its own top edge:
  // the one the geometry gives after the switch, which is not the first here.
  await entry(page, 'Features').click();
  await entry(page, 'Comparison').click();
  await expectMarked(page, 'Comparison');
  await setMode(page, 'circular');
  const atTopEdge = () => page.evaluate(() => {
    const area = document.querySelector('.settings-scroll');
    const edge = area.getBoundingClientRect().top;
    const sections = Array.from(area.querySelectorAll('[data-settings-section]'))
      .filter((section) => section.checkVisibility());
    const reached = sections.filter((section) => section.getBoundingClientRect().top <= edge + 1);
    return (reached.at(-1) ?? sections[0]).dataset.settingsSection;
  });
  const titleOf = Object.fromEntries(SECTIONS.circular);
  expect(await atTopEdge(), 'the settings stay scrolled past Input Genomes').not.toBe('input-genomes');
  await expect.poll(async () => await marked(page).getAttribute('aria-label') === entryName(titleOf[await atTopEdge()])).toBe(true);
  await expect(index(page).locator('[aria-current]')).toHaveCount(1);
});

test('the section title shows below an entry on hover and on keyboard focus, inside the viewport', async ({ page }) => {
  test.setTimeout(180_000);
  await page.setViewportSize({ width: 1280, height: 800 });
  await openApp(page);
  await setMode(page, 'linear');
  await expect(index(page)).toBeVisible();
  await expect(tip(page)).toBeHidden();

  // Hover: below the button by 6 px and centred on it.
  await entry(page, 'Labels').hover();
  await expect(tip(page)).toBeVisible();
  await expect(tip(page)).toHaveText('Labels');
  await expect(tip(page)).toHaveAttribute('aria-hidden', 'true');
  const [button, shown] = [await entry(page, 'Labels').boundingBox(), await tip(page).boundingBox()];
  expect(Math.abs(shown.y - (button.y + button.height) - 6)).toBeLessThanOrEqual(1);
  expect(Math.abs((shown.x + shown.width / 2) - (button.x + button.width / 2))).toBeLessThanOrEqual(1);
  await page.mouse.move(800, 500);
  await expect(tip(page)).toBeHidden();

  // An activation hides it while the pointer is still over the entry.
  await entry(page, 'Labels').hover();
  await expect(tip(page)).toBeVisible();
  await entry(page, 'Labels').click();
  await expect(tip(page)).toBeHidden();
  await page.mouse.move(800, 500);

  // Keyboard focus shows it, Escape hides it and keeps the focus, the next entry shows its own.
  await entry(page, 'Comparison').focus();
  await page.keyboard.press('Shift+Tab');
  await expect(entry(page, 'Input Genomes')).toBeFocused();
  await expect(tip(page)).toHaveText('Input Genomes');
  // Escape hides a shown text and goes no further; with no text shown it reaches the page as usual.
  await page.evaluate(() => {
    window.__escapes = 0;
    document.addEventListener('keydown', (event) => { if (event.key === 'Escape') window.__escapes += 1; });
  });
  await page.keyboard.press('Escape');
  await expect(tip(page)).toBeHidden();
  await expect(entry(page, 'Input Genomes')).toBeFocused();
  expect(await page.evaluate(() => window.__escapes)).toBe(0);
  await page.keyboard.press('Escape');
  expect(await page.evaluate(() => window.__escapes)).toBe(1);
  await expect(entry(page, 'Input Genomes')).toBeFocused();
  await page.keyboard.press('Tab');
  await expect(tip(page)).toHaveText('Comparison');
  await page.keyboard.press('Enter');
  await expect(tip(page)).toBeHidden();
  // Focus that leaves the index hides it.
  await entry(page, 'Layout').focus();
  await page.keyboard.press('Shift+Tab');
  await expect(tip(page)).toBeVisible();
  await page.getByRole('tab', { name: 'Linear', exact: true }).focus();
  await expect(tip(page)).toBeHidden();
  // A mode change (Undo, a loaded Session) that removes the focused entry removes its text.
  await entry(page, 'Input Genomes').focus();
  await page.keyboard.press('Tab');
  await expect(tip(page)).toHaveText('Comparison');
  await page.evaluate(() => { window.__GBDRAW_APP__.selectDrawing('circular'); });
  await expect(index(page).getByRole('button')).toHaveCount(SECTIONS.circular.length);
  await expect(tip(page)).toBeHidden();
  await page.evaluate(() => { window.__GBDRAW_APP__.selectDrawing('linear'); });
  await expect(index(page).getByRole('button')).toHaveCount(SECTIONS.linear.length);

  // At every pane width the text of the first and of the last entry stays 8 px inside the viewport,
  // and the row neither overflows nor scrolls sideways. The switch takes the right end of the first row, so
  // the twelve Linear entries share what is left: 22 px each when there is room, down to 18 px before the
  // list wraps onto a second row.
  for (const width of [240, 320, 500]) {
    await page.evaluate((width) => { window.__GBDRAW_APP__.sidebarWidth = width; }, width);
    const layout = await page.evaluate(() => {
      const pane = document.querySelector('.settings-pane');
      const nav = pane.querySelector('nav[aria-label="Settings sections"]');
      const box = nav.parentElement.getBoundingClientRect();
      const buttons = Array.from(nav.querySelectorAll('button'), (button) => button.getBoundingClientRect());
      const card = getComputedStyle(pane.querySelector('.settings-scroll > .card'));
      const row = getComputedStyle(nav.parentElement);
      const firstRow = buttons.filter((rect) => Math.round(rect.top) === Math.round(buttons[0].top));
      return {
        surface: ['backgroundColor', 'borderTopColor', 'borderTopWidth', 'borderTopLeftRadius', 'boxShadow']
          .filter((property) => row[property] !== card[property]),
        paneWidth: Math.round(pane.getBoundingClientRect().width),
        overflow: Math.max(nav.parentElement.scrollWidth - nav.parentElement.clientWidth, pane.scrollWidth - pane.clientWidth),
        outside: buttons.filter((rect) => rect.left < box.left || rect.right > box.right || rect.top < box.top || rect.bottom > box.bottom).length,
        rows: new Set(buttons.map((rect) => Math.round(rect.top))).size,
        widths: [Math.min(...firstRow.map((rect) => rect.width)), Math.max(...firstRow.map((rect) => rect.width))],
        heights: [...new Set(buttons.map((rect) => Math.round(rect.height)))],
        height: Math.round(box.height),
        gapBelow: Math.round(pane.querySelector('.settings-scroll').getBoundingClientRect().top - box.bottom),
        rightEdge: Math.round(box.right - pane.querySelector('.settings-scroll > .card').getBoundingClientRect().right)
      };
    });
    expect(layout, `pane ${width} px`).toMatchObject({
      surface: [], paneWidth: width, overflow: 0, outside: 0, rows: width === 240 ? 2 : 1, heights: [26], gapBelow: 8, rightEdge: 0
    });
    // The entries of a row are equally wide.
    expect(layout.widths[1] - layout.widths[0], `entry widths at ${width} px`).toBeLessThanOrEqual(0.5);
    expect(layout.widths[0], `entry width at ${width} px`).toBeGreaterThanOrEqual(18);
    expect(layout.widths[1], `entry width at ${width} px`).toBeLessThanOrEqual(22);
    if (width === 500) expect(layout.widths[0]).toBe(22);
    if (width !== 240) expect(layout.height, `row height at ${width} px`).toBe(36);
    for (const title of ['Input Genomes', 'Dinucleotide content/skew']) {
      await entry(page, title).hover();
      await expect(tip(page)).toHaveText(title);
      const [target, box] = [await entry(page, title).boundingBox(), await tip(page).boundingBox()];
      expect(box.x, `${title} at ${width} px: left edge`).toBeGreaterThanOrEqual(8);
      expect(box.x + box.width, `${title} at ${width} px: right edge`).toBeLessThanOrEqual(1280 - 8);
      expect(Math.abs(box.y - (target.y + target.height) - 6), `${title} at ${width} px: below the entry`).toBeLessThanOrEqual(1);
    }
    await page.mouse.move(800, 500);
  }
});

test.describe('390 px touch', () => {
  test.use({ viewport: { width: 390, height: 844 }, hasTouch: true, isMobile: true });

  test('the row is the first row of the settings, a tap opens a section without hover text', async ({ page }) => {
    test.setTimeout(180_000);
    await openApp(page);
    await setMode(page, 'linear');
    await expect(index(page).getByRole('button')).toHaveCount(SECTIONS.linear.length);
    const top = () => index(page).evaluate((nav) => nav.getBoundingClientRect().top);
    // Not pinned: the row scrolls away with the page.
    const before = await top();
    await page.evaluate(() => window.scrollBy(0, 120));
    expect(Math.round(before - await top())).toBe(120);
    await page.evaluate(() => window.scrollTo(0, 0));

    await entry(page, 'Labels').tap();
    const state = await page.evaluate(() => {
      const section = document.querySelector('[data-settings-section="labels"]');
      const summary = section.querySelector(':scope > details > summary');
      return {
        open: section.querySelector(':scope > details').open,
        focused: document.activeElement === summary,
        belowHeader: Math.round(section.getBoundingClientRect().top - document.querySelector('.app-header').getBoundingClientRect().bottom),
        pageOverflow: document.documentElement.scrollWidth - innerWidth
      };
    });
    expect(state).toMatchObject({ open: true, focused: true });
    expect(Math.abs(state.belowHeader), 'the section top is just below the sticky header').toBeLessThanOrEqual(2);
    expect(state.pageOverflow, 'no horizontal page scroll').toBeLessThanOrEqual(0);
    await expect(tip(page)).toBeHidden();
    await expectMarked(page, 'Labels');
    // Back at the top of the page, the first section is the marked one again.
    await page.evaluate(() => window.scrollTo(0, 0));
    await expectMarked(page, 'Input Genomes');
  });
});

// Chromium hides scrollbars in headless runs; a desktop browser shows the
// settings scrollbar beside the cards, and none when the settings fit. This
// case runs in its own browser with scrollbars shown.
test('with classic scrollbars the row ends on the right edge of the cards', async ({ playwright, baseURL }) => {
  test.setTimeout(180_000);
  const browser = await playwright.chromium.launch({ ignoreDefaultArgs: ['--hide-scrollbars'] });
  try {
    const page = await browser.newPage({ baseURL, viewport: { width: 1280, height: 800 } });
    await openApp(page);
    const edges = () => page.evaluate(() => {
      const area = document.querySelector('.settings-scroll');
      return {
        scrollbar: area.offsetWidth - area.clientWidth,
        rightEdge: Math.round(document.querySelector('.settings-index').getBoundingClientRect().right
          - area.querySelector(':scope > .card').getBoundingClientRect().right)
      };
    });
    await expect.poll(edges).toEqual({ scrollbar: 8, rightEdge: 0 });
    await page.setViewportSize({ width: 1280, height: 2400 });
    await expect.poll(edges).toEqual({ scrollbar: 0, rightEdge: 0 });
  } finally {
    await browser.close();
  }
});
