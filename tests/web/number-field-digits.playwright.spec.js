// @ts-check
const { test, expect } = require('@playwright/test');
const { openApp } = require('./helpers/app-lifecycle.cjs');

// A number field reads Persian, Arabic-Indic, Bengali and Devanagari digits and the Arabic
// decimal separator (U+066B), typed on a keyboard layout that writes them or pasted, as 0-9
// and '.' (OV-2780). Chromium drops these characters from a number field without help.

/**
 * Presses one key per character, as a keyboard layout that writes `text` does: each key
 * event carries its character, so the browser fires the same events as for a typed digit.
 * @param {import('@playwright/test').Page} page
 * @param {string} text
 */
const pressCharacters = async (page, text) => {
  const cdp = await page.context().newCDPSession(page);
  for (const character of text) {
    const key = { key: character, code: 'Digit1', windowsVirtualKeyCode: 49 };
    await cdp.send('Input.dispatchKeyEvent', { type: 'keyDown', text: character, unmodifiedText: character, ...key });
    await cdp.send('Input.dispatchKeyEvent', { type: 'keyUp', ...key });
  }
  await cdp.detach();
};

/** @param {import('@playwright/test').Page} page */
const firstNumberField = (page) => page.locator('input[type="number"]:visible:enabled').first();

test('a number field reads digits of other scripts and the Arabic decimal separator typed into it', async ({ page }) => {
  await openApp(page);
  const field = firstNumberField(page);
  for (const [typed, read] of [['۱۲', '12'], ['١٢', '12'], ['১২', '12'], ['१२', '12'], ['۲٫۵', '2.5']]) {
    await field.fill('');
    await field.focus();
    await pressCharacters(page, typed);
    await expect(field, typed).toHaveValue(read);
  }
  await field.fill('1');
  await field.press('Home');
  await pressCharacters(page, '۲');
  await expect(field, 'inserted at the caret').toHaveValue('21');
  await field.press('ControlOrMeta+a');
  await pressCharacters(page, '٣');
  await expect(field, 'replaces the selection').toHaveValue('3');
  await field.press('Tab');
  await expect(field, 'the setting keeps the value').toHaveValue('3');
});

test('a number field reads digits of other scripts in pasted text', async ({ page }) => {
  await openApp(page);
  await page.context().grantPermissions(['clipboard-read', 'clipboard-write']);
  const field = firstNumberField(page);
  await field.fill('');
  await field.focus();
  await page.evaluate(() => navigator.clipboard.writeText('۳٫۷۵'));
  await field.press('ControlOrMeta+v');
  await expect(field).toHaveValue('3.75');
});

test('a numeric text field reads other digits; any other text field keeps them', async ({ page }) => {
  await openApp(page);
  await page.evaluate(() => {
    for (const [id, mode] of [['digits-numeric', 'numeric'], ['digits-name', '']]) {
      const input = document.createElement('input');
      input.type = 'text';
      input.id = id;
      if (mode) input.inputMode = mode;
      document.body.prepend(input);
    }
  });
  await page.locator('#digits-numeric').focus();
  await pressCharacters(page, '۱۲۳');
  await expect(page.locator('#digits-numeric')).toHaveValue('123');
  await page.locator('#digits-name').focus();
  await pressCharacters(page, '۱۲۳');
  await expect(page.locator('#digits-name')).toHaveValue('۱۲۳');
});
