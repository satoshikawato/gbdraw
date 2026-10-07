// UI-01: the component rules of index.html (.card, .form-input, .btn, ...) take
// effect. Tailwind Play compiles @apply only inside <style type="text/tailwindcss">;
// in a plain <style> every rule was dropped and inputs had no frame. The rules that
// refine .form-input must follow it in the same layer, or the compact fields clip
// their text and the "(auto)" values drawn behind transparent inputs disappear.
// Checks of the Web GUI audit 2026-10-05 (PLAN §2.3) at a desktop and a phone width,
// with HmmtDNA and a Depth TSV loaded, a custom track stack, Label Mode Out, and every
// disclosure open.
// G02 (UI-06): hints use .ui-hint; every visible text in the settings pane is 11 px or
// more and reaches 4.5:1 on its own background (icons of icon-only buttons 3:1).
const { test, expect } = require('@playwright/test');
const { readFileSync } = require('node:fs');
const { join } = require('node:path');
const { openApp, reveal } = require('./helpers/app-lifecycle.cjs');

const genbank = readFileSync(join(__dirname, '../test_inputs/HmmtDNA.gbk'));
const depthTsv = ['reference_name\tposition\tdepth', 'NC_012920.1\t1\t10', 'NC_012920.1\t8000\t20', ''].join('\n');

// Text that may stay under the 11 px floor, with the reason. Inactive text (a disabled
// control, its label, or a dimmed group whose controls are all disabled) is exempt
// from the contrast floor, as WCAG 1.4.3 exempts inactive components.
const SMALL_TEXT_EXCEPTIONS = [
  // The "(auto)" value drawn behind a compact input stands for the input's own value,
  // so it keeps the compact input text size (10 px).
  '.auto-value-placeholder'
];

const openAllDetails = (page) => page.evaluate(() => {
  document.querySelectorAll('details').forEach((details) => { details.open = true; });
});

const prepare = async (page) => {
  await openApp(page);
  await page.evaluate((text) => {
    window.__GBDRAW_APP__.setCircularDepthFile(0, new File([text], 'sample.depth.tsv', { type: 'text/plain' }));
  }, depthTsv);
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

// A custom stack with Depth, skew, and annotation rows (Linear: Features, GC, skew, Depth).
const stackWithEveryRowKind = async (page, mode) => {
  await openApp(page);
  await page.evaluate(({ mode, genbank, depthTsv }) => {
    const app = window.__GBDRAW_APP__;
    const file = (text, name) => new File([text], name, { type: 'text/plain' });
    if (mode === 'linear') {
      app.mode = 'linear';
      app.lInputType = 'gb';
      app.setLinearSeqPrimaryFile(0, 'gb', file(genbank, 'HmmtDNA.gbk'));
      app.setLinearDepthFile(app.linearSeqs[0], 0, file(depthTsv, 'sample.depth.tsv'));
    } else {
      app.setCircularDepthFile(0, file(depthTsv, 'sample.depth.tsv'));
    }
  }, { mode, genbank: genbank.toString('utf8'), depthTsv });
  if (mode === 'circular') {
    await page.getByLabel('GenBank/DDBJ File', { exact: true })
      .setInputFiles({ name: 'HmmtDNA.gbk', mimeType: 'text/plain', buffer: genbank });
    await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.circularRecordList.length),
      { timeout: 120_000 }).toBeGreaterThan(0);
  }
  await openAllDetails(page);
  const annotations = page.locator('details').filter({ has: page.locator('summary[aria-label="Region Annotations"]') });
  await annotations.getByRole('button', { name: /Add set/ }).click();
  const panelButton = page.locator(`button[aria-controls="${mode}-custom-track-slots-panel"]`);
  await reveal(panelButton);
  if (await panelButton.getAttribute('aria-expanded') !== 'true') await panelButton.click();
  const useStack = page.locator(`#${mode}-custom-track-slots-panel`).locator('..')
    .getByLabel('Use custom stack', { exact: true });
  if (!await useStack.isChecked()) await useStack.check();
  if (mode === 'linear') {
    // The Linear stack starts from Features only; Reset copies GC, skew, and Depth into it.
    await page.evaluate(() => Object.assign(window.__GBDRAW_APP__.form, { show_gc: true, show_skew: true, show_depth: true }));
    page.once('dialog', (dialog) => dialog.accept());
    await useStack.locator('xpath=ancestor::div[1]').locator('button', { hasText: 'Reset' }).click();
  }
  const panel = page.locator(`#${mode}-custom-track-slots-panel`);
  await panel.getByLabel(`New ${mode} track renderer`, { exact: true }).selectOption('annotations');
  await panel.getByRole('button', { name: 'Add track' }).click();
  await expect(page.locator(`[data-capture^="${mode}-track-slot-"]`)).not.toHaveCount(0);
};

// Computed styles of the visible controls and text levels, after finite CSS
// transitions end (.btn fades its opacity when a row button becomes enabled).
const measure = async (page) => {
  await page.waitForFunction(() => document.getAnimations()
    .every((animation) => animation.playState !== 'running' || animation.effect?.getTiming().iterations === Infinity));
  return page.evaluate((smallTextExceptions) => {
    const auditSettingsPane = () => {
      const pane = document.querySelector('.settings-pane');
      const shown = (element) => element.checkVisibility({ visibilityProperty: true, opacityProperty: true })
        && element.getClientRects().length > 0;
      const rgba = (value) => {
        const parts = (value.match(/rgba?\(([^)]+)\)/)?.[1] || '0 0 0 0').split(/[\s,/]+/).filter(Boolean).map(Number);
        return { r: parts[0], g: parts[1], b: parts[2], a: parts.length > 3 ? parts[3] : 1 };
      };
      const over = (top, bottom) => ({
        r: top.r * top.a + bottom.r * (1 - top.a),
        g: top.g * top.a + bottom.g * (1 - top.a),
        b: top.b * top.a + bottom.b * (1 - top.a),
        a: 1
      });
      const luminance = ({ r, g, b }) => [r, g, b].map((v) => v / 255)
        .map((v) => (v <= 0.03928 ? v / 12.92 : ((v + 0.055) / 1.055) ** 2.4))
        .reduce((sum, v, i) => sum + v * [0.2126, 0.7152, 0.0722][i], 0);
      const ratio = (a, b) => {
        const [hi, lo] = [luminance(a), luminance(b)].sort((x, y) => y - x);
        return (hi + 0.05) / (lo + 0.05);
      };
      // The color painted behind an element: its own and its ancestors' backgrounds over white.
      const backdrop = (element) => {
        const layers = [];
        for (let node = element; node; node = node.parentElement) layers.push(rgba(getComputedStyle(node).backgroundColor));
        return layers.reverse().reduce((under, layer) => (layer.a > 0 ? over(layer, under) : under), { r: 255, g: 255, b: 255, a: 1 });
      };
      const opacity = (element) => {
        let value = 1;
        for (let node = element; node; node = node.parentElement) value *= Number(getComputedStyle(node).opacity);
        return value;
      };
      const contrast = (element) => {
        const style = getComputedStyle(element);
        const color = rgba(style.color);
        const behind = backdrop(element);
        return ratio(over({ ...color, a: color.a * opacity(element) }, behind), behind);
      };
      const CONTROLS = 'input:not([type="hidden"]), select, textarea, button';
      const inactive = (element) => {
        if (element.closest(':disabled, [aria-disabled="true"]')) return true;
        const label = element.closest('label');
        const labelled = label && (label.control || label.querySelector(CONTROLS));
        if (labelled?.disabled) return true;
        // A dimmed group whose controls are all disabled (its help tips still open; the
        // Circular single-record group while Multi-Record Canvas is on, for example).
        for (let node = element; node && node !== pane; node = node.parentElement) {
          if (Number(getComputedStyle(node).opacity) < 1) {
            const controls = Array.from(node.querySelectorAll(CONTROLS)).filter((control) => shown(control) && !control.closest('.help-tip'));
            return controls.length > 0 && controls.every((control) => control.disabled);
          }
        }
        return false;
      };
      const describe = (element, text) => `${text.slice(0, 48)} <${element.tagName.toLowerCase()} class="${(element.getAttribute('class') || '').slice(0, 70)}">`;
      const small = [];
      const faint = [];
      const walker = document.createTreeWalker(pane, NodeFilter.SHOW_TEXT);
      for (let node = walker.nextNode(); node; node = walker.nextNode()) {
        const text = node.textContent.replace(/\s+/g, ' ').trim();
        const element = node.parentElement;
        if (!text || !element || element.closest('select, textarea, .sr-only') || !shown(element)) continue;
        const range = document.createRange();
        range.selectNodeContents(node);
        if (!Array.from(range.getClientRects()).some((rect) => rect.width > 1 && rect.height > 1)) continue;
        const size = Number.parseFloat(getComputedStyle(element).fontSize);
        if (size < 11 && !smallTextExceptions.some((selector) => element.closest(selector))) small.push(`${size}px ${describe(element, text)}`);
        const value = contrast(element);
        if (value < 4.5 && !inactive(element)) faint.push(`${value.toFixed(2)} ${describe(element, text)}`);
      }
      // Icon-only buttons: the glyph is non-text content (WCAG 1.4.11, 3:1).
      const faintIcons = Array.from(pane.querySelectorAll('button')).filter(shown)
        .filter((button) => !inactive(button) && !button.textContent.trim() && button.querySelector('i.ph'))
        .map((button) => ({ button, value: contrast(button.querySelector('i.ph')) }))
        .filter(({ value }) => value < 3)
        .map(({ button, value }) => `${value.toFixed(2)} ${button.getAttribute('aria-label') || button.className}`);
      return { smallText: small, faintText: faint, faintIcons };
    };
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
      ...auditSettingsPane(),
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
  }, SMALL_TEXT_EXCEPTIONS);
};

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
      // 5. Text size and contrast in the settings pane (UI-06).
      expect.soft(result.smallText, 'text under 11 px').toEqual([]);
      expect.soft(result.faintText, 'text under 4.5:1').toEqual([]);
      expect.soft(result.faintIcons, 'icon-only buttons under 3:1').toEqual([]);
    });
  });
}

test('linear: every visible settings text is 11 px or more and 4.5:1 (UI-06)', async ({ page }) => {
  test.setTimeout(240_000);
  await stackWithEveryRowKind(page, 'linear');
  await openAllDetails(page);
  const result = await measure(page);
  expect(result.textControls).toBeGreaterThan(30);
  expect.soft(result.smallText, 'text under 11 px').toEqual([]);
  expect.soft(result.faintText, 'text under 4.5:1').toEqual([]);
  expect.soft(result.faintIcons, 'icon-only buttons under 3:1').toEqual([]);
});
