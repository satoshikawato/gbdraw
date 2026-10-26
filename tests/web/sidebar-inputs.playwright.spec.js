// Sidebar input cards on loaded Gallery Sessions: record readouts before and
// after the source records are inspected.
const { test, expect } = require('@playwright/test');
const { join } = require('node:path');
const { evaluateWithRetainedPromise, generateAndWaitForResult, openApp, openLinearFileCards, waitForAppShell } = require('./helpers/app-lifecycle.cjs');
const { showDrawing, showPhonePane } = require('./helpers/drawing-bar.cjs');

const loadSession = async (page, name) => {
  page.on('dialog', (dialog) => dialog.accept());
  await page.locator('input[accept^=".json,"]').setInputFiles(
    join(process.cwd(), 'gbdraw/web/gallery/sessions', name));
  await page.waitForFunction(() => !window.__GBDRAW_APP__.sessionImportPending
    && window.__GBDRAW_APP__.results.length > 0, null, { timeout: 180000 });
};

// OV-374: after Inspect, the Display start field shows the start the record
// is drawn from, the same value the readout used to show beside an empty field.
test('OV-374: Display start shows the current start after Inspect', async ({ page }) => {
  test.setTimeout(240000);
  await openApp(page);
  await loadSession(page, 'HmmtDNA_ATskew.gbdraw-session.json');
  await page.locator('[data-circular-inspect]').click();
  const start = page.getByLabel(/^Display start NC_012920\.1/);
  await expect(start).toBeVisible({ timeout: 180000 });
  await expect(start).toHaveValue('1');

  // Reset start leaves with its change; the keyboard focus moves to the field.
  await start.fill('100');
  await start.press('Tab');
  const reset = page.getByRole('button', { name: 'Reset start', exact: true });
  await reset.focus();
  await page.keyboard.press('Enter');
  await expect(start).toBeFocused();
  await expect(start).toHaveValue('1');

  // Clearing the field requests no start, also when none was requested, so
  // it shows the current start again.
  await start.fill('');
  await start.press('Tab');
  await expect(start).toHaveValue('1');

  // With Circular record off the start field is disabled; Reset start then
  // moves the keyboard focus to that checkbox instead of dropping it.
  await start.fill('100');
  await start.press('Tab');
  const topology = page.getByLabel(/^Circular record NC_012920\.1/);
  await topology.uncheck();
  await expect(start).toBeDisabled();
  await reset.focus();
  await page.keyboard.press('Enter');
  await expect(reset).toHaveCount(0);
  await expect(topology).toBeFocused();
});

// OV-373: a single-record Circular Session shows its annotation records as a
// plain record name before Inspect, never "target record unavailable".
test('OV-373: annotation rows name their record before Inspect', async ({ page }) => {
  test.setTimeout(240000);
  await openApp(page);
  await loadSession(page, 'tobacco-chloroplast.gbdraw-session.json');
  const panel = page.locator('details').filter({ has: page.locator('summary[aria-label="Region Annotations"]') });
  await panel.locator('summary').first().click();
  await expect(panel.locator('[data-annotation-record-caption]').first()).toHaveText('Record: NC_001879.2');
  await expect(panel).not.toContainText('unavailable');
});

// OV-513: the closed Action select of a Feature Visibility rule shows the
// whole text of every action, at the default sidebar width and at phone width.
// The width a select needs is the width of a copy that holds only that option.
for (const viewport of [{ width: 1600, height: 1000 }, { width: 390, height: 844 }]) {
  test(`OV-513 at ${viewport.width}x${viewport.height} a Feature Visibility rule shows its whole action`, async ({ page }) => {
    test.setTimeout(240000);
    await page.setViewportSize(viewport);
    await openApp(page);
    const card = page.locator('summary[aria-label="Features"]');
    await card.scrollIntoViewIfNeeded();
    await card.press('Enter');
    await page.getByRole('button', { name: 'Add feature visibility rule', exact: true }).click();
    const action = page.getByLabel('Visibility rule 1 action', { exact: true });
    await action.scrollIntoViewIfNeeded();
    const missingWidth = () => action.evaluate((select) => {
      const copy = /** @type {HTMLSelectElement} */ (select.cloneNode(false));
      copy.removeAttribute('aria-label');
      copy.append(select.selectedOptions[0].cloneNode(true));
      copy.style.cssText = 'position: absolute; visibility: hidden; width: auto; min-width: 0; max-width: none;';
      select.after(copy);
      const needed = copy.getBoundingClientRect().width;
      copy.remove();
      return Math.max(0, Math.ceil(needed - select.getBoundingClientRect().width));
    });
    for (const [value, text] of [['show', 'Show'], ['off', 'Off'], ['exclude_matching', 'Leave out of protein search']]) {
      await action.selectOption(value);
      await expect(action.locator('option:checked')).toHaveText(text);
      expect(await missingWidth(), `px missing to show "${text}"`).toBe(0);
    }
    // Record ID, Feature Type and Qualifier share the first line in thirds.
    const shares = await page.evaluate(() => ['record ID', 'feature type', 'qualifier'].map((name) => {
      const field = /** @type {HTMLElement} */ (document.querySelector(`[aria-label="Visibility rule 1 ${name}"]`));
      return field.getBoundingClientRect().width / /** @type {HTMLElement} */ (field.parentElement).getBoundingClientRect().width;
    }));
    for (const share of shares) expect(share, 'a third of the line').toBeGreaterThan(0.3);
  });
}

// Scrolls the selected-feature bar with the mouse wheel, as a user does, until
// the control is in the bar's visible part, then checks that the control is in
// the window and that nothing covers it.
const expectReachableInBar = async (page, control, name) => {
  const delta = await control.evaluate((element) => {
    const bar = /** @type {HTMLElement} */ (element.closest('.feature-selection-toolbar'));
    const [outer, inner] = [bar.getBoundingClientRect(), element.getBoundingClientRect()];
    return inner.top >= outer.top && inner.bottom <= outer.bottom ? 0
      : inner.top - outer.top - (bar.clientHeight - inner.height) / 2;
  });
  if (delta) {
    const box = await page.locator('.feature-selection-toolbar').boundingBox();
    await page.mouse.move(box.x + 6, box.y + 6);
    await page.mouse.wheel(0, delta);
  }
  await expect(control, name).toBeInViewport({ ratio: 1 });
  await expect.poll(() => control.evaluate((element) => {
    const rect = element.getBoundingClientRect();
    return [rect.left + 3, rect.left + rect.width / 2, rect.right - 3]
      .every((x) => element.contains(document.elementFromPoint(x, rect.top + rect.height / 2)));
  }), `${name} is not covered`).toBe(true);
};

// The selected-feature bar names the action as a Feature Visibility rule does.
// At every width the list shows that option in full beside its Apply button,
// the legend name stays readable, no two controls overlap and every Apply
// button stays in the window, with the default and the widest sidebar.
for (const viewport of [{ width: 1024, height: 768 }, { width: 1024, height: 768, sidebar: 500 },
  { width: 1440, height: 900 }, { width: 1440, height: 900, sidebar: 500 }, { width: 390, height: 844 }]) {
  const sidebar = viewport.sidebar ? ` with a ${viewport.sidebar} px sidebar` : '';
  test(`at ${viewport.width}x${viewport.height}${sidebar} the selected-feature visibility list keeps its Apply button in view`, async ({ page }, testInfo) => {
    test.setTimeout(240000);
    await page.setViewportSize({ width: viewport.width, height: viewport.height });
    await openApp(page);
    await loadSession(page, 'HmmtDNA_ATskew.gbdraw-session.json');
    await page.evaluate((width) => {
      const app = window.__GBDRAW_APP__;
      if (width) app.sidebarWidth = width;
      const id = app.extractedFeatures.map((feature) => String(feature?.svg_id || '')).find(Boolean);
      app.selectedFeatureIds = new Set([id]);
      app.selectedFeatureAnchorId = id;
    }, viewport.sidebar || 0);
    // A phone user taps the Preview pane tab (a Load keeps Settings shown) and
    // scrolls the page to the Preview; hovering its controls does that without
    // scrolling the Preview itself.
    await showPhonePane(page, 'preview');
    await page.locator('button[aria-label="Fit to preview"]').hover();
    // Set without a Playwright action, which would scroll the preview into
    // place where a user cannot scroll it.
    const list = page.locator('select[aria-label="Selected feature visibility"]');
    await list.evaluate((select) => {
      select.value = 'exclude_matching';
      select.dispatchEvent(new Event('change'));
    });
    await expect(list.locator('option:checked')).toHaveText('Leave out of protein search');
    await page.screenshot({ path: testInfo.outputPath('selected-feature-bar.png') });
    // The list shows its longest option in full: the text, its padding and the
    // arrow fit the list's own width.
    const fit = await list.evaluate((select) => {
      const style = getComputedStyle(select);
      const probe = document.createElement('span');
      probe.style.font = style.font;
      probe.style.whiteSpace = 'pre';
      probe.style.position = 'absolute';
      probe.textContent = [...select.options].map((option) => option.text).sort((a, b) => b.length - a.length)[0];
      document.body.append(probe);
      const text = probe.getBoundingClientRect().width;
      probe.remove();
      return { need: text + parseFloat(style.paddingLeft) + parseFloat(style.paddingRight) + 16, width: select.clientWidth,
        scrollWidth: select.scrollWidth };
    });
    expect(fit.scrollWidth, 'no hidden overflow').toBeLessThanOrEqual(fit.width);
    expect(fit.need, 'the longest option and the arrow fit').toBeLessThanOrEqual(fit.width);
    const name = page.locator('input[aria-label="Selected feature legend name"]');
    await expect(name).toHaveValue('Selected features');
    expect(await name.evaluate((input) => input.scrollWidth <= input.clientWidth), 'the legend name fits').toBe(true);
    const overlaps = await list.evaluate((select) => {
      const bar = /** @type {HTMLElement} */ (select.closest('.feature-selection-toolbar'));
      const boxes = [...bar.querySelectorAll('input, select, button')].map((control) => ({
        name: control.getAttribute('aria-label') || control.textContent?.trim() || control.tagName,
        rect: control.getBoundingClientRect() })).filter(({ rect }) => rect.width > 0 && rect.height > 0);
      return boxes.flatMap((a, i) => boxes.slice(i + 1)
        .filter((b) => Math.min(a.rect.right, b.rect.right) - Math.max(a.rect.left, b.rect.left) > 1
          && Math.min(a.rect.bottom, b.rect.bottom) - Math.max(a.rect.top, b.rect.top) > 1)
        .map((b) => `${a.name} / ${b.name}`));
    });
    expect(overlaps, 'no two controls overlap').toEqual([]);
    for (const label of ['Apply color and legend name to selected features', 'Apply visibility to selected features',
      'Apply stroke to selected features']) {
      // Nothing above the bar, such as the Editor drawer toggle, covers it.
      await expectReachableInBar(page, page.locator(`button[aria-label="${label}"]`), label);
    }
    await expectReachableInBar(page, list, 'the visibility list');
    if (viewport.width >= 1440 && !viewport.sidebar) {
      // With room for them, the three edits share one row.
      const rows = await page.locator('.feature-selection-toolbar button.btn-sm').evaluateAll((buttons) => buttons
        .filter((button) => /^Apply /.test(button.getAttribute('aria-label') || ''))
        .map((button) => Math.round(button.getBoundingClientRect().top)));
      expect(Math.max(...rows) - Math.min(...rows), `one row of Apply buttons (${rows})`).toBeLessThan(4);
    }
    const apply = page.locator('button[aria-label="Apply visibility to selected features"]');
    const [listBox, applyBox] = [await list.boundingBox(), await apply.boundingBox()];
    expect(Math.abs((listBox.y + listBox.height / 2) - (applyBox.y + applyBox.height / 2)), 'one row').toBeLessThan(2);
  });
}

// OV-1098: on a phone, a fitted diagram leaves a short preview, and the bar
// keeps every control reachable inside it instead of running under the
// preview controls.
test('OV-1098 at 390x844 the selected-feature bar stays reachable in a short preview', async ({ page }, testInfo) => {
  test.setTimeout(240000);
  await page.setViewportSize({ width: 390, height: 844 });
  await openApp(page);
  await loadSession(page, 'HmmtDNA_basic_circular.gbdraw-session.json');
  await showPhonePane(page, 'preview');
  await page.locator('button[aria-label="Fit to preview"]').click();
  await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    const id = app.extractedFeatures.map((feature) => String(feature?.svg_id || '')).find(Boolean);
    app.selectedFeatureIds = new Set([id]);
    app.selectedFeatureAnchorId = id;
  });
  for (const selector of ['select[aria-label="Selected feature visibility"]',
    'button[aria-label="Apply stroke to selected features"]']) {
    await expectReachableInBar(page, page.locator(selector), selector);
  }
  await page.screenshot({ path: testInfo.outputPath('selected-feature-bar-short.png') });
});

// A default color edited under "More colors" stays in that grid and keeps
// focus while it is edited; the grouping does not follow each edit.
test('a color edited under More colors stays there and keeps focus', async ({ page }) => {
  test.setTimeout(240000);
  await openApp(page);
  await loadSession(page, 'HmmtDNA_ATskew.gbdraw-session.json');
  await page.locator('summary[aria-label="Colors"]').click();
  const more = page.locator('[data-default-colors-more]');
  await more.locator('summary').click();
  const tmrna = more.getByLabel('tmRNA feature color mode');
  await expect(tmrna).toBeVisible();
  await expect(tmrna).toHaveValue('color');
  await tmrna.focus();
  await page.keyboard.press('ArrowUp');
  await expect(more.getByLabel('tmRNA feature color mode')).toHaveValue('none');
  await expect(more.getByLabel('tmRNA feature color mode')).toBeFocused();
});

// A Default colors file replaces the color table, so a type it colors joins
// the user's colors in the first list, as a color added by hand does.
test('a color from a Default colors file leaves More colors', async ({ page }) => {
  test.setTimeout(240000);
  await openApp(page);
  await loadSession(page, 'HmmtDNA_ATskew.gbdraw-session.json');
  await page.locator('summary[aria-label="Colors"]').click();
  const more = page.locator('[data-default-colors-more]');
  await more.locator('summary').click();
  await expect(more.getByLabel('tmRNA feature color mode')).toBeVisible();
  await page.getByLabel('Default colors file', { exact: true }).setInputFiles({
    name: 'colors.tsv', mimeType: 'text/plain', buffer: Buffer.from('tmRNA\t#ff0000\n')
  });
  await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    await app.waitForAuxiliaryFileImport(app.cardFiles.d_color);
  });
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.currentColors.tmRNA)).toBe('#ff0000');
  await expect(more.getByLabel('tmRNA feature color mode')).toHaveCount(0);
  await expect(page.getByLabel('tmRNA feature color mode')).toBeVisible();
});

// OV-379: a loaded Linear Session has not read its records; the record select
// says so instead of a load that never runs.
test('OV-379: Linear record selects say Not inspected after Session load', async ({ page }) => {
  test.setTimeout(240000);
  await openApp(page);
  await showDrawing(page, { kind: 'linear' });
  await loadSession(page, 'BGC0000708-BGC0000713.gbdraw-session.json');
  const labels = await page.locator('select[data-record-selector-uid]').evaluateAll((selects) => (
    selects.map((select) => select.options[select.selectedIndex]?.text || '')));
  expect(labels.length).toBeGreaterThan(0);
  expect(new Set(labels)).toEqual(new Set(['Not inspected']));
});

// Linear File cards: which name fields a card shows follows the File's rows
// and stored File values (web.linear.file-default-definition).
const namedGenbank = (id, organism) => {
  const sequence = 'atg'.repeat(100);
  const origin = sequence.match(/.{1,60}/g).map((chunk, index) => (
    `${String(index * 60 + 1).padStart(9)} ${chunk.match(/.{1,10}/g).join(' ')}`)).join('\n');
  return `LOCUS       ${id.padEnd(24)} ${sequence.length} bp    DNA     linear   UNA 01-JAN-2000
DEFINITION  ${organism} test record.
ACCESSION   ${id}
VERSION     ${id}
KEYWORDS    .
SOURCE      ${organism}
  ORGANISM  ${organism}
            .
FEATURES             Location/Qualifiers
     source          1..${sequence.length}
                     /organism="${organism}"
                     /strain="${id}"
     CDS             1..90
                     /product="test protein"
ORIGIN
${origin}
//
`;
};

const openLinearWithFile = async (page, name, content, rows) => {
  await openApp(page);
  await showDrawing(page, { kind: 'linear' });
  await page.locator('[data-linear-source-card] input[type="file"]').first().setInputFiles({
    name, mimeType: 'text/plain', buffer: Buffer.from(content)
  });
  await expect.poll(() => page.evaluate(() => {
    const seqs = window.__GBDRAW_APP__.linearSeqs;
    return seqs.every((seq) => seq.inferred_definition) ? seqs.length : 0;
  }), { timeout: 180000 }).toBe(rows);
  await page.evaluate(() => Object.assign(window.__GBDRAW_APP__.form, {
    legend: 'none', show_gc: false, show_skew: false, show_labels_linear: 'none'
  }));
};

const saveAndReloadSession = async (page) => {
  const download = page.waitForEvent('download');
  await evaluateWithRetainedPromise(page, async () => {
    window.__GBDRAW_APP__.sessionTitle = 'file-card-names';
    await window.__GBDRAW_APP__.saveSessionWithTitle();
  });
  const saved = await (await download).path();
  await page.reload();
  await waitForAppShell(page);
  await page.locator('input[accept^=".json,"]').first().setInputFiles(saved);
  await page.waitForFunction(() => !window.__GBDRAW_APP__.sessionImportPending
    && window.__GBDRAW_APP__.results.length > 0, null, { timeout: 180000 });
  await openLinearFileCards(page);
};

const drawnDefinitionLines = (page) => page.evaluate(() => {
  const doc = new DOMParser().parseFromString(window.__GBDRAW_APP__.currentResults()[0]?.content || '', 'image/svg+xml');
  return [...doc.querySelectorAll('text[data-definition-line-kind]')].map((text) => text.textContent);
});

test('a one-record File shows its record names on the File card', async ({ page }) => {
  test.setTimeout(240000);
  await openApp(page);
  await showDrawing(page, { kind: 'linear' });
  const card = page.locator('[data-linear-source-card]').first();
  const fileNames = card.locator('[data-linear-source-defaults]');
  const definition = page.getByLabel('Definition for sequence 1', { exact: true });
  const subtitle = page.getByLabel('Subtitle / title for sequence 1', { exact: true });
  const recordOptions = card.locator('[data-linear-record-options]');

  // An empty File card has one row: its record names are on the card.
  await expect(definition).toBeVisible();
  await expect(fileNames).toBeHidden();

  await openLinearWithFile(page, 'one.gb', namedGenbank('ONE_REC', 'Aeromonas hydrophila'), 1);
  const inferred = await page.evaluate(() => window.__GBDRAW_APP__.linearSeqs[0].inferred_definition);
  await expect(recordOptions).not.toHaveAttribute('open', '');
  await expect(card.locator(':scope > details > [data-linear-record-names]')).toHaveCount(1);
  await expect(recordOptions.locator('[data-linear-record-names]')).toHaveCount(0);
  await expect(definition).toBeVisible();
  await expect(subtitle).toBeVisible();
  await expect(definition).toHaveAttribute('placeholder', inferred);
  await expect(fileNames).toBeHidden();
  await expect(card.getByText('Names for all records in this file')).toBeHidden();

  // The card fields are the record's own: they write the record's values,
  // and the next Generate draws them.
  await generateAndWaitForResult(page);
  const drawnBefore = await drawnDefinitionLines(page);
  await definition.fill('Card organism');
  await subtitle.fill('Card subtitle');
  await subtitle.press('Tab');
  expect(await page.evaluate(() => {
    const [seq] = window.__GBDRAW_APP__.linearSeqs;
    return [seq.definition, seq.record_subtitle, seq.file_definition, seq.file_subtitle];
  })).toEqual(['Card organism', 'Card subtitle', '', '']);
  expect(await drawnDefinitionLines(page)).toEqual(drawnBefore);
  await generateAndWaitForResult(page);
  expect(await drawnDefinitionLines(page)).toEqual(expect.arrayContaining(['Card organism', 'Card subtitle']));

  // A cleared record name shows the inferred one as its placeholder again.
  await definition.fill('');
  await expect(definition).toHaveAttribute('placeholder', inferred);
  await definition.fill('Card organism');
  await definition.press('Tab');

  // Session round trip: the record keeps its values, and the card, whose
  // records are not inspected after the load, shows the same fields.
  await saveAndReloadSession(page);
  await expect(definition).toBeVisible();
  await expect(definition).toHaveValue('Card organism');
  await expect(subtitle).toHaveValue('Card subtitle');
  await expect(fileNames).toBeHidden();
  expect(await page.evaluate(() => window.__GBDRAW_APP__.linearRecordOptions(window.__GBDRAW_APP__.linearSeqs[0])[0].label))
    .toBe('Not inspected');

  // Inspect leaves the name fields where they are in the settings. The click
  // on Inspect records may scroll the settings, so the place is read from the
  // top of the scrolled content.
  const placeInSettings = async () => {
    const box = await definition.boundingBox();
    const scrolled = await page.locator('#settings-sections').evaluate((area) => area.scrollTop);
    return { ...box, y: box.y + scrolled };
  };
  const before = await placeInSettings();
  await card.getByRole('button', { name: 'Record options for sequence 1', exact: true }).click();
  await card.getByRole('button', { name: 'Inspect records', exact: true }).click({ timeout: 15000 });
  await expect(card.locator('[data-record-rotation]')).toHaveCount(1, { timeout: 180000 });
  await expect(definition).toBeVisible();
  expect(await placeInSettings()).toEqual(before);
  await expect(fileNames).toBeHidden();
});

test('a one-record File that holds a File names value keeps the File names block alone', async ({ page }) => {
  test.setTimeout(240000);
  await openLinearWithFile(page, 'one.gb', namedGenbank('ONE_REC', 'Aeromonas hydrophila'), 1);
  // A one-record card has no field for a File value; a Session can hold one.
  await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    app.setLinearSourceDefaultDefinition(app.linearSourceGroups[0], 'File organism');
  });
  await generateAndWaitForResult(page);
  expect(await drawnDefinitionLines(page)).toContain('File organism');
  await saveAndReloadSession(page);

  const card = page.locator('[data-linear-source-card]').first();
  const fileNames = card.locator('[data-linear-source-defaults]');
  const fileDefinition = page.getByLabel('Default definition for file 1', { exact: true });
  const cardNames = card.locator(':scope > details > [data-linear-record-names]');
  const names = card.locator('[data-linear-record-options] [data-linear-record-names]');
  const definition = page.getByLabel('Definition for sequence 1', { exact: true });
  await expect(fileNames).toBeVisible();
  await expect(fileNames).toContainText('Names for all records in this file');
  await expect(fileDefinition).toHaveValue('File organism');
  // One names block: the record's own fields stay in Record options.
  await expect(cardNames).toHaveCount(0);
  await expect(definition).toBeHidden();
  await card.getByRole('button', { name: 'Record options for sequence 1', exact: true }).click();
  await expect(definition).toBeVisible();
  await expect(definition).toHaveValue('');
  await expect(definition).toHaveAttribute('placeholder', 'File organism');
  await expect(names.getByText('Using file default')).toHaveCount(1);

  // Record, then File, then inferred; Reset returns to the File value.
  await definition.fill('Own organism');
  await definition.press('Tab');
  await names.getByRole('button', { name: 'Reset to default', exact: true }).click();
  await expect(definition).toHaveValue('');
  expect(await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    return app.resolveLinearRecordEffectiveDefinition(app.linearSeqs[0], app.linearSourceGroups[0]);
  })).toBe('File organism');

  // Clearing the last File value keeps the block, alone, while the focus is in it.
  await fileDefinition.fill('');
  await expect(fileDefinition).toBeFocused();
  await expect(fileNames).toBeVisible();
  await expect(definition).toBeHidden();
  // When the focus leaves, the record's own fields take its place on the card.
  await fileDefinition.blur();
  await expect(fileNames).toBeHidden();
  await expect(cardNames).toHaveCount(1);
  await expect(names).toHaveCount(0);
  await expect(definition).toBeVisible();
  const inferred = await page.evaluate(() => window.__GBDRAW_APP__.linearSeqs[0].inferred_definition);
  await expect(definition).toHaveAttribute('placeholder', inferred || 'e.g., Escherichia coli O157:H7 Sakai');
  // Undo brings the value back, and the File block with it.
  await page.getByRole('button', { name: 'Undo', exact: true }).click();
  await expect(fileNames).toBeVisible();
  await expect(fileDefinition).toHaveValue('File organism');
  await expect(cardNames).toHaveCount(0);
});

// OV-486: an empty card shows the record's own name fields. A name typed there
// stays the record's when the file has one record, and becomes the File's name
// when the upload expands the card to several records.
const gffText = (ids) => `##gff-version 3\n${ids.map((id) => [
  `${id}\ttest\tregion\t1\t300\t.\t+\t.\tID=${id}`,
  `${id}\ttest\tgene\t1\t90\t.\t+\t.\tID=gene_${id}`,
  `${id}\ttest\tCDS\t1\t90\t.\t+\t0\tID=cds_${id};Parent=gene_${id};product=test protein`
].join('\n')).join('\n')}\n`;
const fastaText = (ids) => ids.map((id) => `>${id}\n${'atg'.repeat(100)}\n`).join('');
const chooseFiles = async (card, files) => {
  for (const [index, [name, content]] of files.entries()) {
    await card.locator('input[type="file"]').nth(index).setInputFiles({ name, mimeType: 'text/plain', buffer: Buffer.from(content) });
  }
};
const typedFirstUploads = [
  { kind: 'a one-record GenBank file', files: [['one.gb', namedGenbank('ONE_REC', 'Aeromonas hydrophila')]],
    names: [['Typed first', 'Typed subtitle', '', '']] },
  { kind: 'a two-record GenBank file',
    files: [['two.gb', namedGenbank('TWO_A', 'Aeromonas hydrophila') + namedGenbank('TWO_B', 'Aeromonas caviae')]],
    names: [['', '', 'Typed first', 'Typed subtitle'], ['', '', 'Typed first', 'Typed subtitle']] },
  { kind: 'two-record GFF3 + FASTA files', gff: true,
    files: [['two.gff3', gffText(['recA', 'recB'])], ['two.fasta', fastaText(['recA', 'recB'])]],
    names: [['', '', 'Typed first', 'Typed subtitle'], ['', '', 'Typed first', 'Typed subtitle']] }
];
for (const { kind, files, gff = false, names } of typedFirstUploads) {
  test(`a name typed on an empty File card is kept after ${kind} is chosen`, async ({ page }) => {
    test.setTimeout(240000);
    await openApp(page);
    await showDrawing(page, { kind: 'linear' });
    if (gff) await page.locator('input[type="radio"][value="gff"]:visible').check();
    const card = page.locator('[data-linear-source-card]').first();
    await page.getByLabel('Definition for sequence 1', { exact: true }).fill('Typed first');
    await page.getByLabel('Subtitle / title for sequence 1', { exact: true }).fill('Typed subtitle');
    await page.getByLabel('Subtitle / title for sequence 1', { exact: true }).press('Tab');
    await chooseFiles(card, files);
    // A GenBank upload is read when every record has its inferred name; GFF3 + FASTA infer none.
    await expect.poll(() => page.evaluate((isGff) => {
      const seqs = window.__GBDRAW_APP__.linearSeqs;
      return (isGff ? seqs.length > 1 : seqs.every((seq) => seq.inferred_definition)) ? seqs.length : 0;
    }, gff), { timeout: 180000 }).toBe(names.length);
    expect(await page.evaluate(() => window.__GBDRAW_APP__.linearSeqs.map((seq) => (
      [seq.definition, seq.record_subtitle, seq.file_definition, seq.file_subtitle])))).toEqual(names);
    // The card shows the name where it now lives, and Generate draws it.
    await expect(names.length > 1 ? page.getByLabel('Default definition for file 1', { exact: true })
      : page.getByLabel('Definition for sequence 1', { exact: true })).toHaveValue('Typed first');
    await page.evaluate(() => Object.assign(window.__GBDRAW_APP__.form, {
      legend: 'none', show_gc: false, show_skew: false, show_labels_linear: 'none'
    }));
    await generateAndWaitForResult(page);
    const drawn = await drawnDefinitionLines(page);
    expect(drawn).toEqual(expect.arrayContaining(['Typed first', 'Typed subtitle']));
    expect(drawn.filter((text) => text.includes('Aeromonas'))).toEqual([]);
  });
}

// A one-record file that replaces a one-record file keeps the record's own
// name: the card still shows it, and Generate draws it on the new record.
const oneRecordReplacements = [
  { kind: 'GenBank file', drawnId: 'OTHER_REC',
    first: [['one.gb', namedGenbank('ONE_REC', 'Aeromonas hydrophila')]],
    second: [['other.gb', namedGenbank('OTHER_REC', 'Aeromonas caviae')]] },
  { kind: 'GFF3 + FASTA pair', gff: true, drawnId: 'recZ',
    first: [['one.gff3', gffText(['recA'])], ['one.fasta', fastaText(['recA'])]],
    second: [['other.gff3', gffText(['recZ'])], ['other.fasta', fastaText(['recZ'])]] }
];
for (const { kind, first, second, drawnId, gff = false } of oneRecordReplacements) {
  test(`a name typed on a one-record File card stays the record's name when another one-record ${kind} replaces it`, async ({ page }) => {
    test.setTimeout(240000);
    await openApp(page);
    await showDrawing(page, { kind: 'linear' });
    if (gff) await page.locator('input[type="radio"][value="gff"]:visible').check();
    const card = page.locator('[data-linear-source-card]').first();
    const definition = page.getByLabel('Definition for sequence 1', { exact: true });
    await chooseFiles(card, first);
    await definition.fill('Typed first');
    await definition.press('Tab');
    await chooseFiles(card, second);
    await page.evaluate(() => Object.assign(window.__GBDRAW_APP__.form, {
      legend: 'none', show_gc: false, show_skew: false, show_labels_linear: 'none'
    }));
    await generateAndWaitForResult(page);
    expect(await page.evaluate(() => window.__GBDRAW_APP__.linearSeqs.map((seq) => [seq.definition, seq.file_definition])))
      .toEqual([['Typed first', '']]);
    await expect(definition).toHaveValue('Typed first');
    expect(await drawnDefinitionLines(page)).toEqual(expect.arrayContaining(['Typed first', drawnId]));
  });
}

// Unchanged behaviour: this case passes before and after the one-record change.
test('a File with two records keeps its File names block and per-record names', async ({ page }) => {
  test.setTimeout(240000);
  await openLinearWithFile(page, 'two.gb',
    namedGenbank('TWO_A', 'Aeromonas hydrophila') + namedGenbank('TWO_B', 'Aeromonas hydrophila'), 2);
  const card = page.locator('[data-linear-source-card]').first();
  const fileNames = card.locator('[data-linear-source-defaults]');
  const fileDefinition = page.getByLabel('Default definition for file 1', { exact: true });
  const expectTwoRecordCard = async () => {
    await expect(fileNames).toBeVisible();
    await expect(card.locator(':scope > details > [data-linear-record-names]')).toHaveCount(0);
    await expect(card.locator('[data-linear-record-options]').getByLabel(/^Definition for sequence [12]$/)).toHaveCount(2);
    await expect(page.getByLabel('Definition for sequence 1', { exact: true })).toBeHidden();
  };
  await expectTwoRecordCard();
  await expect(fileDefinition).toHaveValue('');

  await fileDefinition.fill('Shared organism');
  await fileDefinition.press('Tab');
  await generateAndWaitForResult(page);
  // Both records share a row and now one name, drawn once; no inferred name is left.
  const drawn = await drawnDefinitionLines(page);
  expect(drawn).toContain('Shared organism');
  expect(drawn.filter((text) => text.includes('Aeromonas'))).toEqual([]);
  await saveAndReloadSession(page);
  await expectTwoRecordCard();
  await expect(fileDefinition).toHaveValue('Shared organism');
  expect(await page.evaluate(() => window.__GBDRAW_APP__.linearSeqs.map((seq) => [seq.definition, seq.file_definition])))
    .toEqual([['', 'Shared organism'], ['', 'Shared organism']]);
});
