// TR-10(a): every visible form control has an author-provided accessible name
// (placeholder text and state-dependent titles are not names).
// TR-10(b), PD-OI-057: every help tip is a keyboard- and tap-openable
// disclosure button outside <label>, described by its own text and referenced
// by the control it explains.
// TR-11: custom-stack rows keep the slot id and renderer readable.
// OV-30 (R12): every visible button has a stable author-provided name (an
// aria-label, aria-labelledby, or visible text), never a Phosphor glyph.
const { test, expect } = require('@playwright/test');
const { readFileSync } = require('node:fs');
const { join } = require('node:path');
const { openApp, reveal } = require('./helpers/app-lifecycle.cjs');
const { seeds, load, popup } = require('./helpers/mode-transition.cjs');

const genbank = readFileSync(join(__dirname, '../test_inputs/HmmtDNA.gbk'), 'utf8');
const depthTsv = ['reference_name\tposition\tdepth', 'NC_012920.1\t1\t10', 'NC_012920.1\t8000\t20', ''].join('\n');

const loadMode = async (page, mode) => {
  await openApp(page);
  await page.evaluate(async ({ mode, genbank, depthTsv }) => {
    const app = window.__GBDRAW_APP__;
    const file = (text, name) => new File([text], name, { type: 'text/plain' });
    if (mode === 'linear') {
      app.setDiagramMode('linear');
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

const PRIVATE_USE = /[\uE000-\uF8FF]/;

// Visible buttons whose name is not author-provided or holds an icon glyph.
// The page's own computation skips aria-hidden subtrees; the browser's name,
// read from the aria snapshot, is checked for glyphs and empty buttons too.
const auditButtons = async (page) => {
  const inPage = await page.evaluate(() => {
    const visible = (element) => element.checkVisibility({ visibilityProperty: true, opacityProperty: false });
    const text = (value) => String(value || '').replace(/\s+/g, ' ').trim();
    const content = (node) => {
      if (node.nodeType === Node.TEXT_NODE) return node.textContent;
      if (node.nodeType !== Node.ELEMENT_NODE || node.getAttribute('aria-hidden') === 'true') return '';
      return Array.from(node.childNodes, content).join(' ');
    };
    const authorName = (element) => text((element.getAttribute('aria-labelledby') || '').split(/\s+/)
      .filter(Boolean).map((id) => content(document.getElementById(id) || document.createTextNode(''))).join(' '))
      || text(element.getAttribute('aria-label'))
      || text(content(element));
    const buttons = Array.from(document.querySelectorAll('button, [role="button"]')).filter(visible);
    const names = new Map(buttons.map((button) => [button, authorName(button)]));
    const problems = [];
    for (const [button, name] of names) {
      const html = text(button.outerHTML).slice(0, 160);
      if (!name) problems.push(`no author-provided name (title only or empty): ${html}`);
      else if (/[\uE000-\uF8FF]/.test(name)) problems.push(`private-use glyph in name "${name}": ${html}`);
    }
    return { count: buttons.length, problems };
  });
  const snapshot = await page.locator('body').ariaSnapshot();
  const browser = snapshot.split('\n').flatMap((line) => {
    const match = line.match(/^\s*- button(?: "((?:[^"\\]|\\.)*)")?/);
    if (!match) return [];
    if (!match[1]) return [`browser reports an unnamed button: ${line.trim()}`];
    return PRIVATE_USE.test(match[1]) ? [`browser name holds a private-use glyph: ${line.trim()}`] : [];
  });
  return { count: inPage.count, problems: [...inPage.problems, ...browser] };
};

// A name must not depend on state: toggle every aria-expanded/aria-pressed
// button once and compare its name.
const auditStateNames = async (page) => {
  const stamped = await page.evaluate(() => {
    const visible = (element) => element.checkVisibility({ visibilityProperty: true, opacityProperty: false });
    const buttons = Array.from(document.querySelectorAll('button[aria-expanded], button[aria-pressed]'))
      .filter((button) => visible(button) && !button.disabled && !button.closest('.help-tip, .app-mode-button')
        && !button.classList.contains('app-mode-button'));
    buttons.forEach((button, index) => button.setAttribute('data-name-probe', String(index)));
    return buttons.length;
  });
  const changed = [];
  for (let index = 0; index < stamped; index += 1) {
    const button = page.locator(`[data-name-probe="${index}"]`);
    if (await button.count() !== 1 || !await button.isVisible()) continue;
    const read = () => button.evaluate((element) => element.getAttribute('aria-label')
      || element.textContent.replace(/\s+/g, ' ').trim());
    const before = await read();
    await button.click({ timeout: 5_000 }).catch(() => {});
    if (await button.count() === 1) {
      const after = await read();
      if (after !== before) changed.push(`${before} -> ${after}`);
      await button.click({ timeout: 5_000 }).catch(() => {});
    }
  }
  return changed;
};

for (const mode of ['circular', 'linear']) {
  test(`${mode} custom stack and Advanced buttons have stable names without icon glyphs`, async ({ page }) => {
    test.setTimeout(240_000);
    await page.setViewportSize({ width: 1280, height: 900 });
    await loadMode(page, mode);
    const audit = await auditButtons(page);
    expect(audit.count).toBeGreaterThan(40);
    expect(audit.problems, 'buttons without a stable author-provided name').toEqual([]);
    const move = mode === 'circular' ? 'Move outside Axis' : 'Move above Axis';
    await expect(page.getByRole('button', { name: new RegExp(`^${mode} track slot \\S+ ${move}$`, 'i') }).first()).toBeVisible();
    expect(await auditStateNames(page), 'button names that change with state').toEqual([]);
  });

  test(`${mode} Result, drawer tabs, and feature popup buttons have stable names`, async ({ browser }) => {
    test.setTimeout(300_000);
    const page = await load(browser, mode === 'circular' ? seeds.circular : seeds.linear);
    try {
      await popup(page);
      const problems = [];
      const record = async (label) => {
        const audit = await auditButtons(page);
        expect(audit.count, `${label} buttons`).toBeGreaterThan(5);
        problems.push(...audit.problems.map((problem) => `${label}: ${problem}`));
      };
      await record('page with drawer and popup');
      for (const tab of ['Details', 'Qualifiers', 'Sequence', 'Edit']) {
        const button = page.locator('.feature-popup').getByRole('button', { name: tab, exact: true });
        if (await button.count()) { await button.click(); await record(`popup ${tab}`); }
      }
      for (const tab of ['Legend', 'Similarity groups', 'Features']) {
        const button = page.locator('.right-drawer').getByRole('button', { name: tab, exact: true });
        if (await button.isEnabled()) { await button.click(); await record(`drawer ${tab}`); }
      }
      expect(problems, 'buttons without a stable author-provided name').toEqual([]);
      expect(await auditStateNames(page), 'button names that change with state').toEqual([]);
    } finally { await page.context().close(); }
  });
}

// Each Legend editor row names its caption field and its icon buttons by its
// caption, so two rows never share a name (R12, gui-fix UI-04).
test('Legend editor rows name their controls by the row caption', async ({ browser }) => {
  test.setTimeout(300_000);
  const page = await load(browser, seeds.circular);
  try {
    await page.locator('.drawer-toggle').click();
    await page.evaluate(() => window.__GBDRAW_APP__.openRightDrawerTab('legend'));
    const drawer = page.locator('.right-drawer');
    const captions = (await page.evaluate(() => window.__GBDRAW_APP__.legendEntries.map((entry) => entry.caption))).slice(0, 2);
    expect(captions).toHaveLength(2);
    const controls = [
      ['textbox', (caption) => `Legend entry ${caption} name`],
      ['button', (caption) => `Move ${caption} up`],
      ['button', (caption) => `Move ${caption} down`],
      ['button', (caption) => `Stroke options for ${caption}`],
      ['button', (caption) => `Remove ${caption}`]
    ];
    for (const [role, nameOf] of controls) {
      const names = captions.map(nameOf);
      expect(new Set(names).size, names.join(' / ')).toBe(2);
      for (const name of names) await expect(drawer.getByRole(role, { name, exact: true })).toHaveCount(1);
    }
    // The controls inside Stroke options carry the caption too.
    for (const caption of captions) {
      await drawer.getByRole('button', { name: `Stroke options for ${caption}`, exact: true }).click();
    }
    const strokeControls = [
      ['combobox', (caption) => `Legend stroke color for ${caption} mode`],
      ['spinbutton', (caption) => `Legend stroke width for ${caption}`],
      ['button', (caption) => `Reset stroke of ${caption} to default`]
    ];
    for (const [role, nameOf] of strokeControls) {
      const names = captions.map(nameOf);
      expect(new Set(names).size, names.join(' / ')).toBe(2);
      for (const name of names) await expect(drawer.getByRole(role, { name, exact: true })).toHaveCount(1);
    }
    for (const caption of captions) {
      await expect(drawer.getByLabel(`Legend stroke color for ${caption}`, { exact: true })).toHaveCount(1);
    }
  } finally { await page.context().close(); }
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
  if (await section.getAttribute('open') === null) await section.locator(':scope > summary').press('Enter');
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
  await expect(button).toHaveAccessibleName('Help: Window');
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

// The author-provided accessible name, as auditButtons and auditPage read it.
const NAME_OF = (element) => {
  const text = (value) => String(value || '').replace(/\s+/g, ' ').trim();
  const content = (node) => {
    if (node.nodeType === Node.TEXT_NODE) return node.textContent;
    if (node.nodeType !== Node.ELEMENT_NODE || node.getAttribute('aria-hidden') === 'true') return '';
    return Array.from(node.childNodes, content).join(' ');
  };
  return text((element.getAttribute('aria-labelledby') || '').split(/\s+/).filter(Boolean)
    .map((id) => content(document.getElementById(id) || document.createTextNode(''))).join(' '))
    || text(element.getAttribute('aria-label'))
    || text(Array.from(element.labels || []).map((label) => content(label)).join(' '))
    || (element.matches('button, [role="button"]') ? text(content(element)) : '');
};

// UI-04, TK-16 (R12): a stack row's controls and help tips carry the slot id,
// so no two controls in one Custom Track Slots panel share a name, and every
// help tip on the page is named after the label it explains.
for (const mode of ['circular', 'linear']) {
  test(`${mode} Custom Track Slots names are unique and help tips name their labels`, async ({ page }) => {
    test.setTimeout(240_000);
    await page.setViewportSize({ width: 1280, height: 900 });
    await loadMode(page, mode);
    const names = await page.evaluate(({ mode, nameOf }) => {
      const name = new Function(`return (${nameOf})`)();
      const visible = (element) => element.checkVisibility({ visibilityProperty: true, opacityProperty: false });
      const panel = document.getElementById(`${mode}-custom-track-slots-panel`);
      const controls = Array.from(panel?.querySelectorAll('input:not([type="file"]), select, textarea, button') || [])
        .filter(visible).map(name);
      const tips = Array.from(document.querySelectorAll('.help-tip > button')).filter(visible).map(name);
      return { controls, tips };
    }, { mode, nameOf: NAME_OF.toString() });
    expect(names.controls.length).toBeGreaterThan(30);
    const repeated = names.controls.filter((name, index) => names.controls.indexOf(name) !== index);
    expect([...new Set(repeated)], 'names repeated inside one Custom Track Slots panel').toEqual([]);
    expect(names.tips.length).toBeGreaterThan(60);
    expect(names.tips.filter((name) => !/^Help: \S/.test(name)), 'help tips not named after a label').toEqual([]);
    await expect(page.getByRole('button', { name: 'Help: Window', exact: true })).toHaveCount(1);
  });
}

// UI-12 (Owner decision 2026-10-07): Generate Diagram is the last stop of the
// settings panel in Tab order, and the bar stays fixed at the bottom.
test('the Generate bar is the last settings Tab stop in both modes and stays fixed', async ({ page }) => {
  test.setTimeout(180_000);
  await page.setViewportSize({ width: 1280, height: 800 });
  await openApp(page);
  for (const mode of ['circular', 'linear']) {
    await page.getByRole('button', { name: mode === 'circular' ? 'Circular' : 'Linear', exact: true }).click();
    const order = await page.evaluate(() => {
      const visible = (element) => element.checkVisibility({ visibilityProperty: true, opacityProperty: false });
      const stops = Array.from(document.querySelector('.settings-pane').querySelectorAll(
        'a[href], button, input, select, textarea, summary, [tabindex]'
      )).filter((element) => visible(element) && !element.disabled && element.tabIndex >= 0);
      const bar = document.querySelector('.generate-bar');
      const box = bar.getBoundingClientRect();
      return {
        last: stops.at(-1)?.getAttribute('aria-label') || '',
        position: getComputedStyle(bar).position,
        bottomGap: Math.round(innerHeight - box.bottom)
      };
    });
    expect(order, mode).toEqual({ last: 'Generate Diagram', position: 'fixed', bottomGap: 32 });
  }
});

// UI-09: a control disabled for a reason the section does not show names it
// in one visible line that the control references.
test('disabled Single-record and LOSAT cache controls reference a visible reason', async ({ page }) => {
  test.setTimeout(180_000);
  await openApp(page);
  await page.getByLabel('GenBank/DDBJ File', { exact: true })
    .setInputFiles({ name: 'HmmtDNA.gbk', mimeType: 'text/plain', buffer: Buffer.from(genbank) });
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.circularRecordList.length),
    { timeout: 120_000 }).toBeGreaterThan(0);
  await page.evaluate(() => document.querySelectorAll('details').forEach((details) => { details.open = true; }));
  for (const name of ['Circular record label', 'Circular record subtitle', 'Circular region start', 'Circular region end',
    'Circular reverse complement', 'Save Raw LOSAT TSV', 'Clear Cache']) {
    const control = page.getByLabel(name, { exact: true }).or(page.getByRole('button', { name, exact: true })).first();
    await expect(control, name).toBeDisabled();
    const reason = await control.evaluate((element) => (element.getAttribute('aria-describedby') || '').split(/\s+/)
      .map((id) => document.getElementById(id)).filter((target) => target?.checkVisibility()
        && /^Available /.test(target.textContent.trim())).map((target) => target.textContent.trim())[0] || '');
    expect(reason, `${name} names why it is disabled`).toMatch(/^Available /);
  }
});

// The visible choice dialogs and the element that has focus.
const focusState = (page) => page.evaluate(() => {
  const active = document.activeElement;
  return {
    inDialog: Boolean(active?.closest('[role="dialog"][aria-modal="true"]')),
    inPopup: Boolean(active?.closest('.feature-popup')),
    name: active?.getAttribute('aria-label') || active?.textContent?.replace(/\s+/g, ' ').trim().slice(0, 40) || active?.tagName
  };
});

const CHOICE_DIALOGS = [
  ['featureVisibilityScopeDialog', { scopes: [{ id: 'product', label: 'All with this product', description: 'Product' }] }, 'Feature Visibility Scope'],
  ['labelTextScopeDialog', { sourceText: 'ND1', matchingCount: 2, featureId: 'f1' }, 'Label Text Scope'],
  ['hiddenLabelTextDialog', { featureId: '', reason: '' }, 'Label Not Shown'],
  ['labelOnDialog', { reason: 'hidden', featureType: 'CDS' }, 'Feature Is Hidden'],
  ['featureStyleScopeDialog', { kind: 'fill' }, 'Color Change Scope'],
  ['legendRenameDialog', { mode: 'scope', oldCaption: 'CDS', newCaption: 'Coding', siblingCount: 1 }, 'Legend Name Scope'],
  ['resetColorDialog', { siblingCount: 1, caption: 'CDS' }, 'Reset Fill Color']
];

// UI-02, UI-03, UI-04: each choice dialog is a named modal that takes focus,
// keeps Tab inside, cancels on Escape without closing the feature popup, and
// returns focus to its opener; the popup takes focus from the keyboard opener
// and returns it on Escape and Close; the drawer tabs announce the shown tab.
test('choice dialogs and the feature popup move focus in and back, and drawer tabs announce the shown tab', async ({ browser }) => {
  test.setTimeout(300_000);
  const page = await load(browser, seeds.circular);
  try {
    await popup(page);
    const featurePopup = page.locator('.feature-popup');
    const opener = featurePopup.getByRole('button', { name: 'Close feature popup', exact: true });
    for (const [dialogState, fields, heading] of CHOICE_DIALOGS) {
      await opener.focus();
      await page.evaluate(({ dialogState, fields }) => {
        Object.assign(window.__GBDRAW_APP__[dialogState], fields, { show: true });
      }, { dialogState, fields });
      const dialog = page.getByRole('dialog', { name: heading });
      await expect(dialog, heading).toHaveAttribute('aria-modal', 'true');
      await expect.poll(() => focusState(page), heading).toMatchObject({ inDialog: true });
      const stops = await dialog.locator('button:enabled').count();
      for (const key of [...Array(stops + 1).fill('Tab'), ...Array(stops + 1).fill('Shift+Tab')]) {
        await page.keyboard.press(key);
        expect((await focusState(page)).inDialog, `${heading}: ${key} stays inside`).toBe(true);
      }
      await page.keyboard.press('Escape');
      // A Cancel handler commits through History, which may first capture the
      // intent; the dialog closes when that returns.
      await expect(dialog, `${heading} closes on Escape`).toHaveCount(0, { timeout: 60_000 });
      await expect(featurePopup, `${heading}: Escape keeps the popup`).toBeVisible();
      await expect(opener, `${heading}: focus returns`).toBeFocused();
    }

    await opener.click();
    await expect(featurePopup).toHaveCount(0);
    const search = page.getByRole('searchbox', { name: 'Search features' });
    await search.fill('ND1');
    await search.press('Enter');
    const open = page.getByRole('button', { name: 'Open active feature', exact: true });
    for (const close of ['Escape', 'Close feature popup']) {
      await open.focus();
      await page.keyboard.press('Enter');
      await expect(featurePopup).toBeVisible();
      await expect.poll(() => focusState(page), 'focus moves into the popup').toMatchObject({ inPopup: true });
      if (close === 'Escape') await page.keyboard.press('Escape');
      else { await opener.focus(); await page.keyboard.press('Enter'); }
      await expect(featurePopup).toHaveCount(0);
      await expect(open, `${close} returns focus to the opener`).toBeFocused();
    }

    if (!await page.locator('.right-drawer').isVisible()) await page.locator('.drawer-toggle').click();
    const drawer = page.locator('.right-drawer');
    for (const tab of ['Legend', 'Features']) {
      await drawer.getByRole('button', { name: tab, exact: true }).click();
      const tabs = await drawer.locator('button[aria-controls^="right-drawer-panel-"]').evaluateAll((buttons) => buttons
        .map((button) => ({ name: button.textContent.trim(), pressed: button.getAttribute('aria-pressed'),
          panel: Boolean(document.getElementById(button.getAttribute('aria-controls'))?.checkVisibility()) })));
      expect(tabs, tab).toEqual(['Legend', 'Features', 'Similarity groups'].map((name) => ({
        name, pressed: String(name === tab), panel: name === tab })));
    }
    for (const placeholder of ['Search similarity groups...']) {
      await expect(drawer.locator(`input[placeholder="${placeholder}"]`)).toHaveAttribute('aria-label', 'Search similarity groups');
    }
    await expect(drawer.locator('select:has(option[value="member_count"])')).toHaveAttribute('aria-label', 'Sort similarity groups');
  } finally { await page.context().close(); }
});

// #941 residual: Escape closes only the top layer. A feature popup opened from
// the Editor drawer closes on Escape and Close with the drawer still open, and
// focus returns to the drawer's Edit button; a second Escape closes the drawer.
test('a feature popup opened from the Editor drawer returns focus to its Edit button', async ({ browser }) => {
  test.setTimeout(300_000);
  const page = await load(browser, seeds.circular);
  try {
    if (!await page.locator('.right-drawer').isVisible()) await page.locator('.drawer-toggle').click();
    const drawer = page.locator('.right-drawer');
    const featurePopup = page.locator('.feature-popup');
    const edit = drawer.getByRole('button', { name: 'Edit', exact: true }).first();
    for (const close of ['Escape', 'Close feature popup']) {
      await edit.focus();
      await page.keyboard.press('Enter');
      await expect(featurePopup).toBeVisible();
      await expect.poll(() => focusState(page), 'focus moves into the popup').toMatchObject({ inPopup: true });
      if (close === 'Escape') await page.keyboard.press('Escape');
      else await featurePopup.getByRole('button', { name: close, exact: true }).press('Enter');
      await expect(featurePopup).toHaveCount(0);
      await expect(drawer, `${close} keeps the drawer open`).toBeVisible();
      await expect(edit, `${close} returns focus to the drawer's Edit button`).toBeFocused();
    }
    await page.keyboard.press('Escape');
    await expect(page.locator('.drawer-toggle')).toHaveAttribute('aria-expanded', 'false');
  } finally { await page.context().close(); }
});

// UI-02 (Owner Q1, Q2 pattern): Label Not Shown asks before anything is
// written, so Cancel applies nothing and records no History step; the popup
// stays open and focus returns to the label text.
test('Label Not Shown Cancel applies nothing and returns focus to the label text', async ({ browser }) => {
  test.setTimeout(300_000);
  const { evaluateWithRetainedPromise } = require('./helpers/app-lifecycle.cjs');
  const { generate } = require('./helpers/mode-transition.cjs');
  const page = await load(browser, seeds.circular);
  try {
    // Show Labels None draws no label, so a text edit leaves its feature unlabeled.
    await page.evaluate(() => { window.__GBDRAW_APP__.form.labels_mode = 'none'; });
    await generate(page);
    const before = await evaluateWithRetainedPromise(page, async () => {
      const app = window.__GBDRAW_APP__;
      const feature = app.extractedFeatures.find((item) => item.type === 'CDS' && !app.getEditableLabelByFeatureId(item.svg_id));
      if (!feature) return null;
      await app.openFeatureEditorFromList(feature, null);
      return { overrides: JSON.stringify(app.featureOverrides), undo: window.__GBDRAW_HISTORY__.getUndoCount() };
    });
    expect(before, 'a feature without a drawn label').not.toBeNull();
    const text = page.locator('.feature-popup input[placeholder="Edit label text"]');
    await text.fill('UNSHOWN_LABEL');
    await text.press('Enter');
    const dialog = page.getByRole('dialog', { name: 'Label Not Shown' });
    await expect(dialog.getByRole('button')).toHaveText([/Show this label/, /Keep hidden/, /Cancel/]);
    await expect.poll(() => focusState(page)).toMatchObject({ inDialog: true });
    await dialog.getByRole('button', { name: 'Cancel', exact: true }).click();
    await expect(dialog).toHaveCount(0);
    await expect(page.locator('.feature-popup')).toBeVisible();
    await expect(text).toBeFocused();
    await expect.poll(() => page.evaluate(() => ({ overrides: JSON.stringify(window.__GBDRAW_APP__.featureOverrides),
      undo: window.__GBDRAW_HISTORY__.getUndoCount() }))).toEqual(before);
  } finally { await page.context().close(); }
});
