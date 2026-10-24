const { test, expect } = require('@playwright/test');
const { load, popup: openFromFeatureList, seeds } = require('./helpers/mode-transition.cjs');
const { showPhonePane } = require('./helpers/drawing-bar.cjs');

// Small UI fixes from the 2026-10-05 re-audit: panel wording that states what a
// preset or the rule order does (FL-11, FL-12), and two overlays that must not
// cover the header controls (Reset alignment) or the footer (match popup), whose
// text also reaches 4.5:1 (GX-15).
// The Session load is the only setup; each overlay is opened through its state.

const open = browser => load(browser, seeds.circular, { width: 1600, height: 1000 });

test('FL-11 the preset panel says what Bakta changes and that the palette is unchanged', async ({ browser }) => {
  test.setTimeout(300000);
  const page = await open(browser);
  try {
    const effects = page.locator('[data-preset-scheme-effects]');
    await expect(effects).toHaveCount(1);
    await effects.evaluate(element => { for (let node = element.closest('details'); node; node = node.parentElement?.closest('details')) node.open = true; });
    await expect(effects).toContainText('the palette is unchanged');
    // The Bakta side effects show when Bakta is the chosen preset.
    await expect(effects).not.toContainText('Bakta');
    await page.getByLabel('Preset scheme', { exact: true }).selectOption('bakta');
    await expect(effects).toContainText('Bakta also sets the CDS default color to #cccccc');
    await expect(effects).toContainText('Legend box size and font size to 12');
    await expect(effects).not.toContainText('update the palette');
  } finally {
    await page.context().close();
  }
});

test('FL-12 two or more specific rules show one line on what Move rule up/down orders', async ({ browser }) => {
  test.setTimeout(300000);
  const page = await open(browser);
  try {
    const hint = page.locator('[data-specific-rule-order-hint]');
    await expect(hint).toHaveCount(0);
    await page.evaluate(async () => {
      const { state } = await import('./js/state.js');
      for (const [i, color] of ['#112233', '#445566'].entries()) {
        state.activeDrawing().manualSpecificRules.push({ feat: 'CDS', qual: 'product', val: `rule${i}`, color, cap: '' });
      }
    });
    await expect(hint).toHaveCount(1);
    await hint.evaluate(element => { for (let node = element.closest('details'); node; node = node.parentElement?.closest('details')) node.open = true; });
    await expect(hint).toBeVisible();
    await expect(hint).toContainText('Rules on the same qualifier: the first match wins.');
  } finally {
    await page.context().close();
  }
});

test('Reset alignment opens as a centred modal that covers no header control', async ({ browser }, info) => {
  test.setTimeout(300000);
  const page = await open(browser);
  try {
    await page.evaluate(() => { window.__GBDRAW_APP__.similarityAlignmentResetDialogOpen = true; });
    const dialog = page.getByRole('dialog', { name: 'Reset alignment', exact: true });
    await expect(dialog).toBeVisible();
    await expect(dialog).toHaveAttribute('aria-modal', 'true');
    const header = await page.locator('header.app-header').boundingBox();
    const box = await dialog.boundingBox();
    expect(box.y).toBeGreaterThanOrEqual(header.y + header.height);
    expect(Math.abs((box.x + box.width / 2) - 800)).toBeLessThan(2);
    await info.attach('reset-alignment', { body: await page.screenshot(), contentType: 'image/png' });
    await page.evaluate(() => { window.__GBDRAW_APP__.similarityAlignmentResetDialogOpen = false; });
    await expect(dialog).toBeHidden();
  } finally {
    await page.context().close();
  }
});

// GX-15: text in the match popup under 4.5:1 on its own ground (a hovered feature row
// has the row's hover ground).
const faintText = popup => popup.evaluate(root => {
  const rgb = value => (value.match(/[\d.]+/g) || []).map(Number);
  const luminance = color => color.slice(0, 3).map(v => v / 255)
    .map(v => (v <= 0.04045 ? v / 12.92 : ((v + 0.055) / 1.055) ** 2.4))
    .reduce((sum, v, i) => sum + v * [0.2126, 0.7152, 0.0722][i], 0);
  const ground = element => {
    for (let node = element; node; node = node.parentElement) {
      const color = rgb(getComputedStyle(node).backgroundColor);
      if (color.length === 3 || color[3] > 0) return color;
    }
    return [255, 255, 255];
  };
  return Array.from(root.querySelectorAll('*'))
    .filter(element => Array.from(element.childNodes).some(node => node.nodeType === 3 && node.textContent.trim()))
    .map(element => {
      const [a, b] = [luminance(rgb(getComputedStyle(element).color)), luminance(ground(element))];
      return { text: element.textContent.trim().slice(0, 40), ratio: Math.round(((Math.max(a, b) + 0.05) / (Math.min(a, b) + 0.05)) * 100) / 100 };
    })
    .filter(({ ratio }) => ratio < 4.5);
});

test('the match popup stays above the footer and its text reaches 4.5:1', async ({ browser }, info) => {
  test.setTimeout(300000);
  const page = await open(browser);
  try {
    await page.evaluate(async () => {
      const { state } = await import('./js/state.js');
      state.clickedPairwiseMatch.value = {
        title: 'Pairwise match',
        subtitle: 'comparison1_match1',
        fill: '#f87171',
        sections: [{
          title: 'Summary',
          rows: Array.from({ length: 24 }, (_, i) => ({ label: `Row ${i + 1}`, value: `value ${i + 1}` }))
        }, {
          title: 'Features',
          featureRows: [{
            key: 'f1', label: 'nad1', subLabel: 'CDS 3307..4262', record: 'NC_012920.1', location: '3307..4262',
            product: 'NADH dehydrogenase subunit 1', canOpen: true, copyText: 'nad1'
          }]
        }]
      };
      state.clickedPairwiseMatchPos.x = 884;
      state.clickedPairwiseMatchPos.y = 552;
    });
    const popup = page.getByRole('dialog', { name: 'Pairwise match details', exact: true });
    await expect(popup).toBeVisible();
    const footer = await page.locator('[data-app-footer]').boundingBox();
    const popupBox = await popup.boundingBox();
    expect(popupBox.y + popupBox.height).toBeLessThanOrEqual(footer.y + 0.5);
    expect(await faintText(popup), 'popup text under 4.5:1').toEqual([]);
    await popup.locator('.pairwise-match-feature-row').first().hover();
    await expect(popup.locator('.pairwise-match-feature-row').first()).toHaveCSS('background-color', 'rgb(239, 246, 255)');
    expect(await faintText(popup), 'popup text under 4.5:1 with a feature row hovered').toEqual([]);
    await info.attach('match-popup', { body: await page.screenshot(), contentType: 'image/png' });
  } finally {
    await page.context().close();
  }
});

// #941 residual: the feature popup keeps the same footer clearance as the match
// popup when it opens near the bottom, when it is dragged down, and when the
// match popup is dragged down. A drag released past the clamp (here on the
// footer) is not a click outside, so the popup stays open.
test('feature and match popups stay above the footer when opened low or dragged down', async ({ browser }) => {
  test.setTimeout(300000);
  const page = await open(browser);
  try {
    const footerTop = async () => (await page.locator('[data-app-footer]').boundingBox()).y;
    const bottomOf = async (locator) => {
      await expect(locator).toBeVisible();
      const box = await locator.boundingBox();
      return box.y + box.height;
    };
    const dragDown = async (locator) => {
      const box = await locator.boundingBox();
      await page.mouse.move(box.x + box.width / 2, box.y + 8);
      await page.mouse.down();
      await page.mouse.move(box.x + box.width / 2, 990, { steps: 4 });
      await page.mouse.up();
    };
    await page.evaluate(() => {
      const app = window.__GBDRAW_APP__;
      app.openFeatureEditorFromList(app.extractedFeatures.find((feature) => feature.type === 'CDS'), { clientX: 700, clientY: 990 });
    });
    const featurePopup = page.locator('.feature-popup').first();
    await expect(featurePopup).toBeVisible();
    expect(await bottomOf(featurePopup), 'opened low').toBeLessThanOrEqual(await footerTop() + 0.5);
    await dragDown(featurePopup);
    expect(await bottomOf(featurePopup), 'dragged down').toBeLessThanOrEqual(await footerTop() + 0.5);
    await page.keyboard.press('Escape');
    await expect(featurePopup).toHaveCount(0);

    await page.evaluate(async () => {
      const { state } = await import('./js/state.js');
      state.clickedPairwiseMatch.value = {
        title: 'Pairwise match', subtitle: 'comparison1_match1', fill: '#f87171',
        sections: [{ title: 'Summary', rows: [{ label: 'Identity', value: '99%' }] }]
      };
      state.clickedPairwiseMatchPos.x = 884;
      state.clickedPairwiseMatchPos.y = 300;
    });
    const matchPopup = page.getByRole('dialog', { name: 'Pairwise match details', exact: true });
    await expect(matchPopup).toBeVisible();
    await dragDown(matchPopup);
    expect(await bottomOf(matchPopup), 'match popup dragged down').toBeLessThanOrEqual(await footerTop() + 0.5);

    // OV-626: a window that gets narrower moves an open popup the least distance back inside it.
    // The footer may wrap to a second line in the narrower window (OV-1427); a popup that was
    // dragged down to it then rises just enough to keep its margin above it (#941). The margin
    // is placePopupWithin's default in gbdraw/web/js/utils/popup-bounds.js.
    const popupMargin = 12;
    const insideNarrowerWindow = async (locator, name) => {
      const before = await locator.boundingBox();
      expect(before.x + before.width, `${name} starts right of the narrower window`).toBeGreaterThan(1100);
      await page.setViewportSize({ width: 1100, height: 1000 });
      await expect.poll(async () => { const box = await locator.boundingBox(); return box && box.x + box.width; }, { message: `${name} inside the window` })
        .toBeLessThanOrEqual(1100);
      await expect(locator, `${name} stays open`).toBeVisible();
      const after = await locator.boundingBox();
      expect(after.x, `${name} inside the window`).toBeGreaterThanOrEqual(0);
      const footerY = await footerTop();
      expect(after.y + after.height, `${name} stays above the footer`).toBeLessThanOrEqual(footerY - popupMargin + 0.5);
      expect(after.y, `${name} keeps its height on the page unless the footer pushes it up`)
        .toBeCloseTo(Math.min(before.y, footerY - popupMargin - after.height), 0);
      await page.setViewportSize({ width: 1600, height: 1000 });
    };
    await insideNarrowerWindow(matchPopup, 'match popup');
    await page.keyboard.press('Escape');
    await expect(matchPopup).toHaveCount(0);
    await page.evaluate(() => {
      const app = window.__GBDRAW_APP__;
      app.openFeatureEditorFromList(app.extractedFeatures.find((feature) => feature.type === 'CDS'), { clientX: 700, clientY: 400 });
    });
    await expect(featurePopup).toBeVisible();
    await insideNarrowerWindow(featurePopup, 'feature popup');
  } finally {
    await page.context().close();
  }
});

// OV-371: a feature or ribbon clicked at 1600x1000 opened its popup at the
// lower right, 520 px tall with its content scrolled, over the clicked item.
// OV-384: the popup (z-50) could lie over the header (z-40) and cover its
// buttons, such as Undo, when opened near the top-left or dragged up.
const BGC = 'gbdraw/web/gallery/sessions/BGC0000708-BGC0000713.gbdraw-session.json';
const LIVE = 'fed46a3a6_record_1';

// The first visible point of an SVG element, where a pointer would reach it.
const visiblePoint = (page, selector) => page.locator(selector).first().evaluate((element) => {
  const box = element.getBoundingClientRect();
  const steps = [0.5, 0.4, 0.6, 0.3, 0.7, 0.2, 0.8, 0.1, 0.9];
  for (const fy of steps) {
    for (const fx of steps) {
      const [x, y] = [box.left + box.width * fx, box.top + box.height * fy];
      const hit = document.elementFromPoint(x, y);
      if (hit === element || element.contains(hit)) return { x, y };
    }
  }
  throw new Error('no visible point');
});

const clickShape = async (page, selector) => {
  const point = await visiblePoint(page, selector);
  await page.locator(selector).first().dispatchEvent('click', { clientX: point.x, clientY: point.y });
  return point;
};

const popupBounds = (page, popup, bodySelector) => popup.evaluate((element, selector) => {
  const box = element.getBoundingClientRect();
  const body = element.querySelector(selector);
  return {
    top: box.top, bottom: box.bottom, left: box.left, right: box.right,
    header: document.querySelector('header.app-header').getBoundingClientRect().bottom,
    footer: document.querySelector('[data-app-footer]').getBoundingClientRect().top,
    width: innerWidth,
    scrolls: body.scrollHeight > body.clientHeight + 1
  };
}, bodySelector);

const expectPlaced = (placed, point, name) => {
  expect(placed.top, `${name} below the header`).toBeGreaterThanOrEqual(placed.header);
  expect(placed.bottom, `${name} above the footer`).toBeLessThanOrEqual(placed.footer + 0.5);
  expect(placed.left, `${name} inside the window`).toBeGreaterThanOrEqual(0);
  expect(placed.right, `${name} inside the window`).toBeLessThanOrEqual(placed.width);
  const covers = point.x >= placed.left && point.x <= placed.right && point.y >= placed.top && point.y <= placed.bottom;
  expect(covers, `${name} leaves the clicked point visible`).toBe(false);
  // The popup opens as tall as its content needs, up to the space it has.
  if (placed.scrolls) expect(placed.bottom - placed.top, `${name} uses the height it has`).toBeGreaterThan(placed.footer - placed.header - 48);
};

test('OV-371 feature and ribbon popups open inside the window beside the clicked item', async ({ browser }) => {
  test.setTimeout(300000);
  const page = await load(browser, BGC, { width: 1600, height: 1000 });
  try {
    await page.getByRole('button', { name: 'Fit to preview', exact: true }).click();
    const feature = await clickShape(page, `[data-gbdraw-feature-id="${LIVE}"]`);
    const featurePopup = page.locator('.feature-popup[role="dialog"]').first();
    // The first popup after a Session load waits for the Python runtime.
    await expect(featurePopup).toBeVisible({ timeout: 180000 });
    expectPlaced(await popupBounds(page, featurePopup, '.feature-popup-body'), feature, 'feature popup');
    await page.keyboard.press('Escape');
    await expect(featurePopup).toHaveCount(0);

    const ribbon = await clickShape(page, 'path[data-match-kind]');
    const matchPopup = page.getByRole('dialog', { name: 'Pairwise match details', exact: true });
    await expect(matchPopup).toBeVisible();
    expectPlaced(await popupBounds(page, matchPopup, '.pairwise-match-popup-body'), ribbon, 'match popup');
  } finally {
    await page.context().close();
  }
});

test('OV-384 the feature popup never covers the header buttons', async ({ browser }) => {
  test.setTimeout(300000);
  const page = await load(browser, BGC, { width: 1600, height: 1000 });
  try {
    const popup = page.locator('.feature-popup[role="dialog"]').first();
    const undoReachable = () => page.getByRole('button', { name: 'Undo', exact: true }).evaluate((button) => {
      const box = button.getBoundingClientRect();
      return button.contains(document.elementFromPoint(box.left + box.width / 2, box.top + box.height / 2));
    });
    const headerBottom = async () => (await page.locator('header.app-header').boundingBox()).height;
    // Opened at a point above the window (the Gallery test of #986 did this).
    await page.locator(`[data-gbdraw-feature-id="${LIVE}"]`).first().dispatchEvent('click', { clientX: -100, clientY: -100 });
    await expect(popup).toBeVisible({ timeout: 180000 });
    expect((await popup.boundingBox()).y, 'opened near the top-left').toBeGreaterThanOrEqual(await headerBottom());
    expect(await undoReachable(), 'Undo reachable with the popup opened near the top-left').toBe(true);
    // Dragged up past the header.
    const box = await popup.boundingBox();
    await page.mouse.move(box.x + box.width / 2, box.y + 6);
    await page.mouse.down();
    await page.mouse.move(box.x + box.width / 2 + 300, 0, { steps: 6 });
    await page.mouse.up();
    expect((await popup.boundingBox()).y, 'dragged up').toBeGreaterThanOrEqual(await headerBottom());
    expect(await undoReachable(), 'Undo reachable with the popup dragged up').toBe(true);
  } finally {
    await page.context().close();
  }
});

const expectInsideWindow = async (page, popup, name) => {
  const box = await popup.boundingBox();
  const header = await page.locator('header.app-header').boundingBox();
  const { width, height } = page.viewportSize();
  expect(box.y, `${name}: below the header`).toBeGreaterThanOrEqual(header.y + header.height);
  expect(box.y + box.height, `${name}: above the window bottom`).toBeLessThanOrEqual(height);
  expect(box.x, `${name}: inside the window`).toBeGreaterThanOrEqual(0);
  expect(box.x + box.width, `${name}: inside the window`).toBeLessThanOrEqual(width);
};

// OV-496: on a one-column page the notice banner lies below the sidebar, far
// under the window while the page is scrolled above the Result. It was still
// read as the popups' top limit, so a popup dragged or opened there left the
// window.
test('OV-496 at 390x844 a popup stays inside the window, also with the Result scrolled away', async ({ browser }) => {
  test.setTimeout(300000);
  const page = await load(browser, BGC, { width: 390, height: 844 });
  try {
    await showPhonePane(page, 'preview');
    const featureSelector = `[data-gbdraw-feature-id="${LIVE}"]`;
    const expectInside = (popup, name) => expectInsideWindow(page, popup, name);
    const scrollAboveResult = async () => {
      await page.mouse.move(195, 40);
      await page.mouse.wheel(0, -100000);
      await expect.poll(() => page.evaluate(() => window.scrollY)).toBe(0);
      expect((await page.locator('[data-app-notice]').boundingBox()).y, 'the notice is below the window').toBeGreaterThan(844);
    };

    // A tap on a feature that is in view.
    await page.locator(featureSelector).first().evaluate(element => element.scrollIntoView({ block: 'center', inline: 'center' }));
    const tapped = await visiblePoint(page, featureSelector);
    await page.mouse.click(tapped.x, tapped.y);
    const featurePopup = page.locator('.feature-popup[role="dialog"]').first();
    await expect(featurePopup).toBeVisible({ timeout: 180000 });
    await expectInside(featurePopup, 'feature popup opened by a tap');

    // The page scrolled back to the sidebar, then the popup dragged by its top edge.
    await scrollAboveResult();
    const box = await featurePopup.boundingBox();
    await page.mouse.move(box.x + box.width / 2, box.y + 6);
    await page.mouse.down();
    await page.mouse.move(box.x + box.width / 2, box.y + 26, { steps: 4 });
    await page.mouse.up();
    await expectInside(featurePopup, 'feature popup dragged with the Result scrolled away');
    await featurePopup.getByRole('button', { name: 'Close feature popup', exact: true }).click();
    await expect(featurePopup).toHaveCount(0);

    // Enter on a focused Legend row after the page was scrolled away from it.
    await page.locator('.gbdraw-preview-surface g[data-legend-key="Core biosynthetic genes"]').first().focus();
    await scrollAboveResult();
    await page.keyboard.press('Enter');
    const legendPopup = page.locator('.pairwise-match-popup[role="dialog"]');
    await expect(legendPopup).toBeVisible();
    await expectInside(legendPopup, 'Legend row popup opened by Enter');
    await legendPopup.getByRole('button', { name: 'Close popup', exact: true }).click();
    await expect(legendPopup).toHaveCount(0);
  } finally {
    await page.context().close();
  }
});

// OV-503: Edit in the Editor drawer's Features list opened the feature popup
// beside the clicked button, over the drawer's title, tabs, list search and
// upper rows, so another tab or the next row was out of reach. With room
// beside the drawer the popup opens there. A docked Editor (one column) has no
// such room: the popup then opens as it does at any clicked point.
for (const viewport of [{ width: 1600, height: 1000 }, { width: 1280, height: 800 }]) {
  test(`OV-503 at ${viewport.width}x${viewport.height} a feature popup opened from the Features list leaves the Editor drawer uncovered`, async ({ browser }) => {
    test.setTimeout(300000);
    const page = await load(browser, seeds.circular, viewport);
    try {
      const drawer = page.locator('.right-drawer');
      const featurePopup = page.locator('.feature-popup[role="dialog"]').first();
      const expectBesideDrawer = async (name) => {
        const [popupBox, drawerBox] = [await featurePopup.boundingBox(), await drawer.boundingBox()];
        expect(popupBox.x + popupBox.width, `${name}: the popup ends left of the drawer`).toBeLessThanOrEqual(drawerBox.x);
        await expectInsideWindow(page, featurePopup, name);
        return popupBox;
      };
      const first = await openFromFeatureList(page);
      await expectBesideDrawer('first row');
      expect(await drawer.locator('button[aria-controls^="right-drawer-panel-"]').evaluateAll(tabs => tabs.map((tab) => {
        const box = tab.getBoundingClientRect();
        return tab.contains(document.elementFromPoint(box.left + box.width / 2, box.top + box.height / 2));
      })), 'the three drawer tabs are reachable').toEqual([true, true, true]);
      // OV-505: so is the drawer's toggle, which sticks out of the drawer's left edge.
      expect(await page.locator('.drawer-toggle').evaluate((toggle) => {
        const box = toggle.getBoundingClientRect();
        return toggle.contains(document.elementFromPoint(box.left + 4, box.top + box.height / 2));
      }), 'the drawer toggle is reachable').toBe(true);

      // The next row's Edit is reachable: a click opens that feature, beside the drawer again.
      await drawer.getByRole('button', { name: 'Edit', exact: true }).nth(1).click();
      await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.clickedFeature?.id)).not.toBe(first.id);
      const opened = await expectBesideDrawer('next row');

      // A popup the user drags over the drawer stays where it was put.
      await page.mouse.move(opened.x + opened.width / 2, opened.y + 6);
      await page.mouse.down();
      await page.mouse.move(opened.x + opened.width / 2 + 120, opened.y + 6, { steps: 4 });
      await page.mouse.up();
      await featurePopup.getByRole('button', { name: 'Details', exact: true }).click();
      expect((await featurePopup.boundingBox()).x, 'dragged popup').toBeCloseTo(opened.x + 120, 0);
    } finally {
      await page.context().close();
    }
  });
}

test('OV-503 at 390x844 a feature popup opened from the docked Features list stays inside the window', async ({ browser }) => {
  test.setTimeout(300000);
  const page = await load(browser, seeds.circular, { width: 390, height: 844 });
  try {
    await openFromFeatureList(page);
    await expectInsideWindow(page, page.locator('.feature-popup[role="dialog"]').first(), 'feature popup opened from the docked list');
  } finally {
    await page.context().close();
  }
});

// OV-430: the closed Editor drawer lies right of the Result card. The card is
// not a scroll container, so no scroll (a click's scroll-into-view, a focus)
// moves the Result sideways out of reach.
test('OV-430 the Result card cannot be scrolled sideways by the closed Editor drawer', async ({ browser }) => {
  test.setTimeout(300000);
  const page = await open(browser);
  try {
    const card = page.getByRole('region', { name: 'Result Preview', exact: true });
    await expect(page.locator('.right-drawer[aria-hidden="true"]')).toHaveCount(1);
    expect(await card.evaluate((element) => {
      element.scrollLeft = 400;
      return element.scrollLeft;
    })).toBe(0);
  } finally {
    await page.context().close();
  }
});
