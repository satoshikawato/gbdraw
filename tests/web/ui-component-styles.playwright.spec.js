// UI-01: the component rules of index.html (.card, .form-input, .btn, ...) take
// effect. Tailwind Play compiles @apply only inside <style type="text/tailwindcss">;
// in a plain <style> every rule was dropped and inputs had no frame. The rules that
// refine .form-input must follow it in the same layer, or the compact fields clip
// their text and the "(auto)" values drawn behind transparent inputs disappear.
// Checks of the Web GUI audit 2026-10-05 (PLAN §2.3) at a desktop and a phone width,
// with HmmtDNA loaded, a custom track stack, Label Mode Out, and every disclosure open.
const { test, expect } = require('@playwright/test');
const { readFileSync } = require('node:fs');
const { join } = require('node:path');
const { openApp, reveal } = require('./helpers/app-lifecycle.cjs');

const genbank = readFileSync(join(__dirname, '../test_inputs/HmmtDNA.gbk'));

const openAllDetails = (page) => page.evaluate(() => {
  document.querySelectorAll('details').forEach((details) => { details.open = true; });
});

const prepare = async (page) => {
  await openApp(page);
  await page.getByLabel('GenBank/DDBJ File', { exact: true })
    .setInputFiles({ name: 'HmmtDNA.gbk', mimeType: 'text/plain', buffer: genbank });
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.circularRecordList.length),
    { timeout: 120_000 }).toBeGreaterThan(0);
  const panelButton = page.locator('button[aria-controls="circular-custom-track-slots-panel"]');
  await reveal(panelButton);
  if (await panelButton.getAttribute('aria-expanded') !== 'true') await panelButton.click();
  const useStack = page.locator('#circular-custom-track-slots-panel').locator('..')
    .getByLabel('Use custom stack', { exact: true });
  if (!await useStack.isChecked()) await useStack.check();
  await openAllDetails(page);
  await page.locator('#circular-label-mode').selectOption('out');
  await openAllDetails(page);
  await expect(page.locator('[data-capture^="circular-track-slot-"]').first()).toBeVisible();
};

// Computed styles of the visible controls and text levels.
const measure = (page) => page.evaluate(() => {
  const NON_TEXT = new Set(['checkbox', 'radio', 'range', 'color', 'file', 'button', 'submit', 'reset', 'image']);
  const visible = (element) => element.checkVisibility({ visibilityProperty: true })
    && element.getClientRects().length > 0;
  const px = (value) => Number.parseFloat(value);
  const name = (element) => element.getAttribute('aria-label') || element.id || element.className;
  const all = (selector) => Array.from(document.querySelectorAll(selector)).filter(visible);
  const level = (selector) => [...new Set(all(selector).map((element) => {
    const style = getComputedStyle(element);
    return `${style.fontSize} ${style.fontWeight} ${style.textTransform}`;
  }))].sort();
  const textControls = all('input, select').filter((element) => !NON_TEXT.has(element.type));
  const overlays = all('.auto-value-input, .track-slot-geometry-input');
  return {
    lang: document.documentElement.lang,
    formInputBorder: getComputedStyle(document.querySelector('.form-input')).borderTopWidth,
    cardRadius: px(getComputedStyle(document.querySelector('.card')).borderTopLeftRadius),
    textControls: textControls.length,
    clipped: textControls.filter((element) => {
      const style = getComputedStyle(element);
      return element.clientHeight - px(style.paddingTop) - px(style.paddingBottom) < px(style.fontSize) * 1.05;
    }).map(name),
    overlays: overlays.length,
    opaqueOverlays: overlays
      .filter((element) => getComputedStyle(element).backgroundColor !== 'rgba(0, 0, 0, 0)').map(name),
    cardTitles: level('.settings-scroll .card > .card-header, .settings-scroll .card > details > summary'),
    groupHeadings: level('.settings-scroll h4'),
    fieldLabels: level('.settings-scroll .input-label'),
    hints: all('.settings-scroll .ui-hint').map((element) => {
      const style = getComputedStyle(element);
      return {
        size: px(style.fontSize),
        lineHeight: px(style.lineHeight) / px(style.fontSize),
        color: style.color,
        transform: style.textTransform
      };
    })
  };
});

for (const viewport of [{ width: 1280, height: 800 }, { width: 390, height: 844 }]) {
  test.describe(() => {
    test.use({ viewport });

    test(`at ${viewport.width} px the component rules apply, no control clips its text, and the type levels hold (UI-01)`, async ({ page }) => {
      test.setTimeout(240_000);
      await prepare(page);
      const result = await measure(page);

      expect.soft(result.lang, 'the page language (UI-05)').toBe('en');
      // 1. The component layer is compiled.
      expect.soft(result.formInputBorder, '.form-input has its 1 px frame').toBe('1px');
      expect.soft(result.cardRadius, '.card has rounded corners').toBeGreaterThan(0);
      // 2. Every visible text input and select shows a full line of text.
      expect(result.textControls).toBeGreaterThan(50);
      expect.soft(result.clipped, 'controls whose content box is shorter than their text').toEqual([]);
      // 3. The "(auto)" values behind the inputs stay visible.
      expect(result.overlays).toBeGreaterThan(0);
      expect.soft(result.opaqueOverlays, 'auto-value inputs with an opaque background').toEqual([]);
      // 4. Type levels: card title 14 px; group heading 11 px uppercase; field label
      // 12 px, or 11 px where a compact label keeps text-[10px] under the 11 px floor;
      // hint 11 px or more in the muted color.
      expect.soft(result.cardTitles, 'card titles').toEqual(['14px 600 none']);
      expect.soft(result.groupHeadings, 'group headings').toEqual(['11px 600 uppercase']);
      expect.soft(result.fieldLabels, 'field labels').toContain('12px 500 none');
      expect.soft(result.fieldLabels.filter((value) => !['11px 500 none', '12px 500 none'].includes(value)),
        'field labels other than 11-12 px, medium, sentence case').toEqual([]);
      expect(result.hints.length).toBeGreaterThan(0);
      expect.soft(result.hints.filter((hint) => hint.size < 11 || hint.transform !== 'none'
        || hint.color !== 'rgb(100, 116, 139)' || Math.abs(hint.lineHeight - 1.45) > 0.01), 'hints').toEqual([]);
    });
  });
}
