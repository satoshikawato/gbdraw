// TR-10(a): every visible form control has an author-provided accessible name
// (placeholder text and state-dependent titles are not names).
// TR-10(b), PD-OI-057: every help tip is a keyboard- and tap-openable
// disclosure button outside <label>, described by its own text and referenced
// by the control it explains.
// TR-11: custom-stack rows keep the slot id and renderer readable.
const { test, expect } = require('@playwright/test');
const { readFileSync } = require('node:fs');
const { join } = require('node:path');
const { openApp, reveal } = require('./helpers/app-lifecycle.cjs');

const genbank = readFileSync(join(__dirname, '../test_inputs/HmmtDNA.gbk'), 'utf8');
const depthTsv = ['reference_name\tposition\tdepth', 'NC_012920.1\t1\t10', 'NC_012920.1\t8000\t20', ''].join('\n');

const loadMode = async (page, mode) => {
  await openApp(page);
  await page.evaluate(async ({ mode, genbank, depthTsv }) => {
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
  }, { mode, genbank, depthTsv });
  if (mode === 'circular') {
    await page.getByLabel('GenBank/DDBJ File', { exact: true })
      .setInputFiles({ name: 'HmmtDNA.gbk', mimeType: 'text/plain', buffer: Buffer.from(genbank) });
    await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.circularRecordList.length),
      { timeout: 120_000 }).toBeGreaterThan(0);
  }
  const panelButton = page.locator(`button[aria-controls="${mode}-custom-track-slots-panel"]`);
  await reveal(panelButton);
  if (await panelButton.getAttribute('aria-expanded') !== 'true') await panelButton.click();
  const useStack = page.locator(`#${mode}-custom-track-slots-panel`).locator('..')
    .getByLabel('Use custom stack', { exact: true });
  if (!await useStack.isChecked()) await useStack.check();
  if (mode === 'linear') {
    // The Linear stack starts from Features only; Reset copies the simple
    // GC, skew, and Depth tracks into it so every row type is present.
    await page.evaluate(() => Object.assign(window.__GBDRAW_APP__.form, { show_gc: true, show_skew: true, show_depth: true }));
    page.once('dialog', (dialog) => dialog.accept());
    await useStack.locator('xpath=ancestor::div[1]').locator('button', { hasText: 'Reset' }).click();
  }
  await expect(page.locator(`[data-capture^="${mode}-track-slot-"]`).first()).toBeVisible();
  await page.evaluate(() => document.querySelectorAll('details').forEach((details) => { details.open = true; }));
  const annotations = page.locator('details').filter({ has: page.locator('summary[aria-label="Region Annotations"]') });
  await annotations.getByRole('button', { name: /Add set/ }).click();
  await annotations.getByRole('button', { name: /Coordinates/ }).click();
};

// Returns the problems found among visible controls and help tips.
const auditPage = (page) => page.evaluate(() => {
  const CONTROLS = 'input:not([type="hidden"]):not([type="file"]), select, textarea';
  const visible = (element) => element.checkVisibility({ visibilityProperty: true, opacityProperty: false });
  const text = (value) => String(value || '').replace(/\s+/g, ' ').trim();
  const describe = (element) => text(element.outerHTML).slice(0, 180);
  const authorName = (element) => text((element.getAttribute('aria-labelledby') || '').split(/\s+/)
    .filter(Boolean).map((id) => document.getElementById(id)?.textContent || '').join(' '))
    || text(element.getAttribute('aria-label'))
    || text(Array.from(element.labels || []).map((label) => label.textContent).join(' '));
  const controls = Array.from(document.querySelectorAll(CONTROLS)).filter(visible);
  const unnamed = controls.filter((element) => !authorName(element)).map(describe);
  const visibleTips = Array.from(document.querySelectorAll('.help-tip')).filter(visible);
  const tips = [];
  for (const tip of visibleTips) {
    const button = tip.querySelector(':scope > button');
    const descriptionId = button?.getAttribute('aria-describedby') || '';
    const description = descriptionId ? document.getElementById(descriptionId) : null;
    const problems = [];
    if (!button || button.disabled || button.tabIndex < 0) problems.push('not focusable');
    if (!text(description?.textContent)) problems.push('no description');
    // The description is referenced, not read in place: the text around the
    // tip stays its visible label.
    if (description && !description.hidden) problems.push('description is in the reading order');
    if (description && tip.contains(description)) problems.push('description inside the tip wrapper');
    if (tip.closest('label')) problems.push('inside label');
    if (tip.parentElement.closest('button, [role="button"]:not(summary)')) problems.push('inside a button');
    // A heading or summary computes its name from content, so one that holds a
    // tip keeps its visible name only through aria-label.
    const named = tip.parentElement.closest('summary, h1, h2, h3, h4, h5, h6');
    if (named && !text(named.getAttribute('aria-label'))) problems.push('changes a heading or summary name');
    // The nearest ancestor holding another interactive element is the tip's
    // field; when that field holds exactly one form control, the control is
    // the one the tip explains.
    const others = (element) => Array.from(element.querySelectorAll(`${CONTROLS}, button`))
      .filter((candidate) => candidate !== button && visible(candidate));
    let field = tip.parentElement;
    while (field && others(field).length === 0) field = field.parentElement;
    const interactive = field ? others(field) : [];
    const fieldControls = interactive.filter((element) => element.matches(CONTROLS));
    if (descriptionId && fieldControls.length === 1 && interactive.length === 1
      && !(fieldControls[0].getAttribute('aria-describedby') || '').split(/\s+/).includes(descriptionId)) {
      problems.push(`control does not reference the tip: ${describe(fieldControls[0])}`);
    }
    if (problems.length) tips.push({ tip: text(description?.textContent || tip.textContent).slice(0, 80), problems });
  }
  const ids = Array.from(document.querySelectorAll('[id]'), (element) => element.id);
  const duplicateIds = [...new Set(ids.filter((id, index) => ids.indexOf(id) !== index))];
  const danglingDescriptions = Array.from(document.querySelectorAll('[aria-describedby]'))
    .filter(visible)
    .flatMap((element) => element.getAttribute('aria-describedby').split(/\s+/).filter(Boolean)
      .filter((id) => !document.getElementById(id)).map((id) => `${id} <- ${describe(element)}`));
  const counts = { controls: controls.length, tips: visibleTips.length };
  return { counts, unnamed, tips, duplicateIds, danglingDescriptions };
});

for (const mode of ['circular', 'linear']) {
  test(`${mode} controls have names and every help tip is a reachable disclosure`, async ({ page }) => {
    test.setTimeout(240_000);
    await page.setViewportSize({ width: 1280, height: 900 });
    await loadMode(page, mode);
    const audit = await auditPage(page);
    expect(audit.counts.controls).toBeGreaterThan(80);
    expect(audit.counts.tips).toBeGreaterThan(70);
    expect.soft(audit.unnamed, 'controls without an author-provided name').toEqual([]);
    expect.soft(audit.tips, 'help tips that are not reachable disclosures').toEqual([]);
    expect.soft(audit.duplicateIds, 'duplicate ids').toEqual([]);
    expect(audit.danglingDescriptions, 'aria-describedby without a target').toEqual([]);

    const enable = page.locator(`[data-capture^="${mode}-track-slot-"] input[type="checkbox"]`).first();
    const name = await enable.evaluate((element) => element.getAttribute('aria-label'));
    await enable.click();
    await expect(enable).toHaveAccessibleName(name);
    await enable.click();
    await expect(enable).toHaveAccessibleName(name);
  });
}

const windowTip = async (page) => {
  const section = page.locator('details').filter({ has: page.locator('summary[aria-label="Dinucleotide content/skew"]') });
  if (await section.getAttribute('open') === null) await section.locator(':scope > summary').click();
  const input = section.getByRole('spinbutton', { name: 'Window', exact: true });
  const descriptionId = await input.getAttribute('aria-describedby', { timeout: 30_000 });
  expect(descriptionId).toBeTruthy();
  return { section, input, button: page.locator(`.help-tip > button[aria-describedby="${descriptionId}"]`) };
};

test('help tips open from the keyboard and on hover, and Escape closes them', async ({ page }) => {
  test.setTimeout(180_000);
  await openApp(page);
  const { section, input, button } = await windowTip(page);
  const tooltip = page.locator('[role="tooltip"]');
  // The visible label text stays exactly the label; the tip text is only a
  // referenced description.
  await expect(section.getByText('Window', { exact: true })).toHaveCount(1);
  await expect(input).toHaveAccessibleDescription(/Window size for GC content\/skew/);
  await expect(button).toHaveAccessibleName('Help');
  await expect(button).toHaveAccessibleDescription(/Window size for GC content\/skew/);
  await input.focus();
  await page.keyboard.press('Shift+Tab');
  await expect(button).toBeFocused();
  await expect(button).toHaveAttribute('aria-expanded', 'true');
  await expect(tooltip).toContainText('Window size for GC content/skew');
  await page.keyboard.press('Escape');
  await expect(tooltip).toHaveCount(0);
  await expect(button).toHaveAttribute('aria-expanded', 'false');
  await page.keyboard.press('Enter');
  await expect(tooltip).toContainText('Window size for GC content/skew');
  await page.keyboard.press('Enter');
  await expect(tooltip).toHaveCount(0);
  await page.keyboard.press('Tab');
  await expect(input).toBeFocused();
  await button.hover();
  await expect(tooltip).toBeVisible();
  await page.mouse.move(0, 0);
  await expect(tooltip).toHaveCount(0);
});

test.describe('390 px touch', () => {
  test.use({ viewport: { width: 390, height: 844 }, hasTouch: true, isMobile: true });

  test('a tap opens and closes a help tip without horizontal page scroll', async ({ page }) => {
    test.setTimeout(180_000);
    await openApp(page);
    const { button } = await windowTip(page);
    const tooltip = page.locator('[role="tooltip"]');
    await button.tap();
    await expect(tooltip).toBeVisible();
    await expect(tooltip).toContainText('Window size for GC content/skew');
    const box = await tooltip.boundingBox();
    expect(box.x).toBeGreaterThanOrEqual(0);
    expect(box.x + box.width).toBeLessThanOrEqual(390);
    await button.tap();
    await expect(tooltip).toHaveCount(0);
    await button.tap();
    await expect(tooltip).toBeVisible();
    await page.touchscreen.tap(5, 5);
    await expect(tooltip).toHaveCount(0);
    expect(await page.evaluate(() => document.documentElement.scrollWidth)).toBeLessThanOrEqual(390);
  });
});

for (const width of [390, 1280, 1920]) {
  for (const mode of ['circular', 'linear']) {
    test(`${mode} stack rows keep slot id and renderer readable at ${width} px`, async ({ page }) => {
      test.setTimeout(240_000);
      await page.setViewportSize({ width, height: 900 });
      await loadMode(page, mode);
      const rows = await page.locator(`[data-capture^="${mode}-track-slot-"]`).evaluateAll((elements) => (
        elements.map((row) => {
          const rect = (element) => element.getBoundingClientRect();
          const [id, renderer] = [row.querySelector('input:not([type="checkbox"])'), row.querySelector('select')];
          const rowBox = rect(row);
          return {
            slot: row.dataset.capture,
            id: Math.round(rect(id).width),
            renderer: Math.round(rect(renderer).width),
            buttons: Array.from(row.querySelector('.track-slot-row-actions')?.querySelectorAll('button') || [])
              .filter((button) => rect(button).left < rowBox.left || rect(button).right > rowBox.right + 0.5).length
          };
        })
      ));
      expect(rows.length).toBeGreaterThan(2);
      for (const row of rows) {
        expect(row.id, `${row.slot} slot id width`).toBeGreaterThanOrEqual(96);
        expect(row.renderer, `${row.slot} renderer width`).toBeGreaterThanOrEqual(96);
        expect(row.buttons, `${row.slot} clipped action buttons`).toBe(0);
      }
      await expect(page.locator(`[data-capture^="${mode}-track-slot-"] .track-slot-row-actions`))
        .toHaveCount(rows.length);
      expect(await page.evaluate(() => document.documentElement.scrollWidth)).toBeLessThanOrEqual(width);
    });
  }
}
