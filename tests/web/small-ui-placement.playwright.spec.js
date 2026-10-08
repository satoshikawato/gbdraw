const { test, expect } = require('@playwright/test');
const { load, seeds } = require('./helpers/mode-transition.cjs');

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
        state.manualSpecificRules.push({ feat: 'CDS', qual: 'product', val: `rule${i}`, color, cap: '' });
      }
    });
    await expect(hint).toHaveCount(1);
    await hint.evaluate(element => { for (let node = element.closest('details'); node; node = node.parentElement?.closest('details')) node.open = true; });
    await expect(hint).toBeVisible();
    await expect(hint).toContainText('same qualifier key');
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
