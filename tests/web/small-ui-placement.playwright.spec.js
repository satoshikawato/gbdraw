const { test, expect } = require('@playwright/test');
const { load, seeds } = require('./helpers/mode-transition.cjs');

// Small UI fixes from the 2026-10-05 re-audit: panel wording that states what a
// preset or the rule order does (FL-11, FL-12), and two overlays that must not
// cover the header controls (Reset alignment) or the footer (match popup).
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

test('the match popup stays above the footer when it opens low in the window', async ({ browser }, info) => {
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
    await info.attach('match-popup', { body: await page.screenshot(), contentType: 'image/png' });
  } finally {
    await page.context().close();
  }
});
