// The drawing bar (PD-OI-090 〔drawing bar と切り替え〕): one tab per drawing,
// named by the drawing, under the header and above both panes. It works at
// 390 px with the compact Preview and Editor conditions of PD-OI-038 (OIC-025),
// and its tab menu renames a drawing as one History step.
const { test, expect } = require('@playwright/test');
const path = require('node:path');
const { readFileSync } = require('node:fs');
const { gunzipSync } = require('node:zlib');
const { openApp, evaluateWithRetainedPromise } = require('./helpers/app-lifecycle.cjs');
const { drawingIds, drawingTab, offerTab, shownDrawing, showDrawing, takeLinearOffer } = require('./helpers/drawing-bar.cjs');

const loadGallerySession = (page, filename) => evaluateWithRetainedPromise(page, async (path) => {
  const name = path.split('/').pop();
  const response = await fetch(path);
  if (!response.ok) throw new Error(`Could not load ${name}: ${response.status}`);
  const file = new File([await response.arrayBuffer()], name, { type: 'application/json' });
  return window.__GBDRAW_APP__.importSession({ target: { files: [file], value: 'selected' } });
}, `/gbdraw/web/gallery/sessions/${filename}`);

const tabName = (page, drawingId) => page.evaluate(async (id) => (
  (await import('./js/state.js')).state.drawingTabs.value.find((tab) => tab.id === id).name), drawingId);

test.describe('390 px', () => {
  test.use({ viewport: { width: 390, height: 844 }, hasTouch: true });

  test('the phone layout shows one picker row and one pane at a time, switched by reveal, Generate and the Editor', async ({ page }) => {
    test.setTimeout(240_000);
    await openApp(page);
    expect((await loadGallerySession(page, 'HmmtDNA_ATskew.gbdraw-session.json')).status).toBe('ok');
    // The Session holds one Circular drawing; the list below needs a second one, made through the offer.
    await takeLinearOffer(page);
    await showDrawing(page, { kind: 'circular' });
    const ids = await drawingIds(page);
    const shown = (await shownDrawing(page)).id;
    const other = ids.find((id) => id !== shown);
    const picker = page.locator('.drawing-bar .drawing-picker');
    const paneTab = (pane) => page.locator(`.phone-pane-tabs [role="tab"][data-phone-pane="${pane}"]`);
    const settings = page.locator('#settings-sections');
    const result = page.locator('#result-pane');
    await expect(page.locator('.drawing-bar-tabs')).toBeHidden();
    await expect(picker).toBeVisible();
    expect(await page.evaluate(() => document.documentElement.scrollWidth)).toBeLessThanOrEqual(390);

    // The picker names the shown drawing and lists every drawing in tab order; a tap shows another.
    await picker.locator('summary').tap();
    const items = picker.locator('.action-menu-item');
    await expect(items).toHaveCount(ids.length);
    await items.nth(ids.indexOf(other)).tap();
    await expect.poll(async () => (await shownDrawing(page)).id).toBe(other);
    // The keyboard reaches it too: Arrow Down opens the list on its first item.
    await picker.locator('summary').focus();
    await page.keyboard.press('ArrowDown');
    for (let index = 0; index < ids.indexOf(shown); index += 1) await page.keyboard.press('ArrowDown');
    await page.keyboard.press('Enter');
    await expect.poll(async () => (await shownDrawing(page)).id).toBe(shown);
    const menuBox = await page.locator('.drawing-bar-menu').evaluate(async (menu) => {
      menu.open = true;
      await new Promise(requestAnimationFrame);
      const box = menu.querySelector('.action-menu-items').getBoundingClientRect();
      menu.open = false;
      return { left: box.left, right: box.right };
    });
    expect(menuBox.left, 'the Drawing menu stays on screen').toBeGreaterThanOrEqual(0);
    expect(menuBox.right).toBeLessThanOrEqual(390);

    // Settings | Preview: one pane at a time, the ARIA tab pattern with manual activation. A Load with a
    // Result shows the Preview, as Generate does, and a drawing switch keeps the pane.
    await expect(paneTab('preview')).toHaveAttribute('aria-selected', 'true');
    await expect(settings).toBeHidden();
    await paneTab('settings').tap();
    await expect(result).toBeHidden();
    await paneTab('preview').tap();
    await expect(settings).toBeHidden();
    await expect(result).toBeVisible();
    await paneTab('preview').focus();
    await page.keyboard.press('ArrowLeft');
    await expect(paneTab('settings')).toBeFocused();
    await page.keyboard.press('Enter');
    await expect(paneTab('settings')).toHaveAttribute('aria-selected', 'true');
    await expect(settings).toBeVisible();

    // A path that reveals a setting shows Settings; Generate and opening the Editor show the Preview.
    await paneTab('preview').tap();
    await page.evaluate(() => window.__GBDRAW_APP__.revealSettingsSection('input-genomes'));
    await expect(paneTab('settings')).toHaveAttribute('aria-selected', 'true');
    await expect(settings).toBeVisible();
    await page.locator('.generate-bar button').click();
    await expect(paneTab('preview')).toHaveAttribute('aria-selected', 'true');
    await page.waitForFunction(() => !window.__GBDRAW_APP__.processing, null, { timeout: 180_000 });
    await expect(page.getByRole('region', { name: 'Result Preview' })).toBeVisible();
    // The approved images (Owner R17-10 A): the Result comes in with its workspace at the pane's top, and the
    // canvas with the controls below it fills the pane; the Result header and the search lie above it.
    const workspaceFit = () => page.evaluate(() => {
      const box = (selector) => document.querySelector(selector).getBoundingClientRect();
      return { top: Math.round(box('.preview-workspace').top - box('#result-pane').top),
        end: Math.round(box('#result-pane').bottom - box('.preview-controls').bottom) };
    });
    await expect.poll(workspaceFit).toEqual({ top: 0, end: 0 });
    await paneTab('settings').tap();
    await page.evaluate(() => window.__GBDRAW_APP__.openRightDrawerTab('features'));
    await expect(paneTab('preview')).toHaveAttribute('aria-selected', 'true');
    expect(await page.evaluate(() => document.documentElement.scrollWidth)).toBeLessThanOrEqual(390);

    // PD-OI-038 at 390 x 844: with the Editor open, the canvas spans the Result card and shows at least
    // 200 px between the header and the Generate bar, the Editor docks under it, and its toggle takes a tap.
    await expect(page.locator('.right-drawer')).toBeVisible();
    const geometry = await page.getByRole('region', { name: 'Result Preview' }).evaluate(async (preview) => {
      await Promise.all(document.getAnimations().map((animation) => animation.finished.catch(() => {})));
      const canvas = preview.querySelector('.preview-canvas');
      const box = (element) => element.getBoundingClientRect();
      const style = getComputedStyle(preview);
      const generate = document.querySelector('.generate-bar');
      const bottom = getComputedStyle(generate).position === 'fixed' ? box(generate).top : innerHeight;
      const toggle = box(document.querySelector('.drawer-toggle'));
      return {
        canvasLeft: box(canvas).left - box(preview).left - parseFloat(style.borderLeftWidth),
        canvasRight: box(preview).right - parseFloat(style.borderRightWidth) - box(canvas).right,
        visibleHeight: Math.min(box(canvas).bottom, bottom) - Math.max(box(canvas).top, box(document.querySelector('.phone-pane-tabs')).bottom),
        drawerBelow: box(document.querySelector('.right-drawer')).top >= box(canvas).bottom - 0.5,
        toggleHit: document.querySelector('.drawer-toggle').contains(document.elementFromPoint(
          toggle.left + toggle.width / 2, toggle.top + toggle.height / 2)),
        pageWidth: document.documentElement.scrollWidth,
        // Opening the Editor brings the workspace to the pane's top: the canvas keeps 13rem, the whole dock
        // ends above the Generate bar, and its feature rows start in view.
        workspaceTop: Math.round(box(preview.querySelector('.preview-workspace')).top - box(document.querySelector('#result-pane')).top),
        canvasHeight: Math.round(box(canvas).height),
        dockAboveGenerate: box(document.querySelector('.right-drawer')).bottom <= bottom,
        rowsInView: [...document.querySelectorAll('.right-drawer [data-feature-list-row]')]
          .filter((row) => box(row).top >= box(canvas).bottom && box(row).bottom <= bottom).length
      };
    });
    expect(geometry).toMatchObject({ canvasLeft: 0, canvasRight: 0, drawerBelow: true, toggleHit: true,
      workspaceTop: 0, canvasHeight: 13 * 16, dockAboveGenerate: true });
    expect(geometry.rowsInView).toBeGreaterThanOrEqual(2);
    expect(geometry.visibleHeight).toBeGreaterThanOrEqual(200);
    expect(geometry.pageWidth).toBeLessThanOrEqual(390);
  });

  // A Load whose shown drawing has a Result shows the Preview pane, as Generate does;
  // a Load of a Session without a Result keeps the pane it started from.
  test('a Session Load shows the Preview pane when the shown drawing has a Result and keeps Settings when it has none', async ({ page }) => {
    test.setTimeout(240_000);
    await openApp(page);
    const tab = (pane) => page.locator(`.phone-pane-tabs [role="tab"][data-phone-pane="${pane}"]`);
    await expect(tab('settings')).toHaveAttribute('aria-selected', 'true');
    await page.locator('input[accept^=".json,"]').setInputFiles(path.join(process.cwd(), 'tests/fixtures/sessions/settings-only.v42.json.gz'));
    await expect(page.locator('[data-app-notice]')).toContainText('Session loaded');
    await page.waitForFunction(() => !window.__GBDRAW_APP__.sessionImportPending);
    expect(await page.locator('.origin-top svg').count()).toBe(0);
    await expect(tab('settings')).toHaveAttribute('aria-selected', 'true');
    // The pane after a Load depends only on the loaded Session: with the Preview pane shown, a Load without a
    // Result switches to Settings.
    await tab('preview').click();
    await expect(tab('preview')).toHaveAttribute('aria-selected', 'true');
    await page.locator('input[accept^=".json,"]').setInputFiles(path.join(process.cwd(), 'tests/fixtures/sessions/settings-only.v42.json.gz'));
    await expect(tab('settings')).toHaveAttribute('aria-selected', 'true');
    expect((await loadGallerySession(page, 'HmmtDNA_ATskew.gbdraw-session.json')).status).toBe('ok');
    await expect(tab('preview')).toHaveAttribute('aria-selected', 'true');
    await expect(page.locator('.origin-top svg')).toBeVisible();
  });
});

// 195 x 422 is a 390 x 844 phone at 200 % zoom (Owner R17-12, OV-1802 and OV-1772): the page scrolls as a whole,
// so the header, the drawing bar and the tabs scroll away and the figure of a Result can be brought into view
// above the fixed Generate bar, whose short label keeps one line.
test.describe('390 px offer', () => {
  test.use({ viewport: { width: 390, height: 844 }, hasTouch: true });

  test('a new Session has one Circular drawing; the Linear offer in the picker makes the Linear drawing as one step', async ({ page }) => {
    test.setTimeout(180_000);
    await openApp(page);
    expect(await drawingIds(page), 'a new Session starts with one Circular drawing').toEqual(['circular']);
    const picker = page.locator('.drawing-bar .drawing-picker');
    await picker.locator('summary').click();
    const items = picker.locator('.action-menu-item');
    await expect(items).toHaveText(['Circular', 'Linear']);
    await expect(items.nth(1).locator('.ph-plus'), 'the offer item carries the plus').toHaveCount(1);
    await items.nth(1).click();
    await expect.poll(() => drawingIds(page)).toEqual(['circular', 'linear']);
    expect((await shownDrawing(page)).id, 'the new Linear drawing is shown').toBe('linear');
    expect(await page.evaluate(() => window.__GBDRAW_HISTORY__.undoLabel())).toBe('New drawing');
    await picker.locator('summary').click();
    await expect(items, 'no offer once a Linear drawing exists').toHaveText(['Circular', 'Linear']);
    await expect(items.locator('.ph-plus')).toHaveCount(0);
    await page.keyboard.press('Escape');

    await page.getByRole('button', { name: 'Undo', exact: true }).first().click();
    await expect.poll(() => drawingIds(page)).toEqual(['circular']);
    expect((await shownDrawing(page)).id, 'Undo shows the drawing shown before').toBe('circular');
  });
});

test.describe('195 px', () => {
  test.use({ viewport: { width: 195, height: 422 }, hasTouch: true });

  test('the Settings pane without a Result does not run past the window width (OV-1810)', async ({ page }) => {
    test.setTimeout(120_000);
    await openApp(page);
    const width = await page.evaluate(() => ({ scroll: document.documentElement.scrollWidth, window: document.documentElement.clientWidth }));
    expect(width.scroll).toBeLessThanOrEqual(width.window);
  });

  test('the page scrolls as a whole, the Preview pane is at least 70 % of the window, and Generate keeps one line', async ({ page }) => {
    test.setTimeout(240_000);
    await openApp(page);
    const label = () => page.locator('.generate-bar button').evaluate((button) => {
      const shown = [...button.querySelectorAll('span > span')].filter((span) => span.checkVisibility());
      return { text: shown.map((span) => span.textContent.trim()), lines: shown.map((span) => span.getClientRects().length) };
    });
    expect(await label()).toEqual({ text: ['Generate'], lines: [1] });
    await expect(page.locator('.generate-bar button')).toHaveAccessibleName('Generate Diagram');
    expect((await loadGallerySession(page, 'HmmtDNA_ATskew.gbdraw-session.json')).status).toBe('ok');
    const previewTab = page.locator('.phone-pane-tabs [role="tab"][data-phone-pane="preview"]');
    if (await previewTab.getAttribute('aria-selected') !== 'true') await previewTab.tap();
    await expect(page.getByRole('region', { name: 'Result Preview' })).toBeVisible();
    expect(await label()).toEqual({ text: ['Generate'], lines: [1] });
    await page.locator('.preview-canvas').evaluate((canvas) => canvas.scrollIntoView({ block: 'start', behavior: 'instant' }));
    const view = await page.evaluate(() => {
      const box = (selector) => document.querySelector(selector).getBoundingClientRect();
      const canvas = box('.preview-canvas');
      const generateTop = box('.generate-bar').top;
      const visible = Math.min(canvas.bottom, generateTop) - Math.max(canvas.top, 0);
      const probe = document.elementFromPoint(canvas.left + canvas.width / 2, Math.max(canvas.top, 0) + Math.min(visible, 200) / 2);
      return { paneShare: box('#result-pane').height / innerHeight, headerTop: box('.app-header').bottom <= 0,
        visible: Math.round(visible), figureHit: Boolean(probe?.closest('.preview-canvas')) };
    });
    expect(view.paneShare, 'the Preview pane is at least 70 % of the window').toBeGreaterThanOrEqual(0.7);
    expect(view.headerTop, 'the header scrolled away with the page').toBe(true);
    expect(view.visible, 'the canvas in view above the Generate bar').toBeGreaterThanOrEqual(200);
    expect(view.figureHit, 'the figure takes the point in view').toBe(true);
  });
});

test('Rename in the tab menu is one History step that Undo restores; a repeated name is refused with its reason', async ({ page }) => {
  test.setTimeout(180_000);
  await page.setViewportSize({ width: 1280, height: 900 });
  await openApp(page);
  const [first] = await drawingIds(page);
  const second = await takeLinearOffer(page);
  const name = await tabName(page, second);
  const firstName = await tabName(page, first);
  const menu = page.locator('.drawing-bar-menu');
  const dialog = page.getByRole('dialog', { name: 'Rename drawing' });
  const rename = async (value) => {
    await menu.locator('summary').click();
    await menu.getByRole('button', { name: 'Rename…', exact: true }).click();
    await expect(dialog.getByRole('textbox', { name: 'Drawing name' })).toHaveValue(name);
    await dialog.getByRole('textbox', { name: 'Drawing name' }).fill(value);
    await dialog.getByRole('button', { name: 'Rename', exact: true }).click();
  };

  await rename(firstName.toUpperCase());
  await expect(dialog).toContainText('Another drawing has this name.');
  await expect(dialog.getByRole('textbox', { name: 'Drawing name' })).toHaveValue(firstName.toUpperCase());
  await dialog.getByRole('button', { name: 'Cancel', exact: true }).click();
  await expect(dialog).toHaveCount(0);
  expect(await tabName(page, second), 'a refused rename changes nothing').toBe(name);

  await rename('Phage comparison');
  await expect(dialog).toHaveCount(0);
  await expect(drawingTab(page, second)).toHaveAccessibleName('Phage comparison');
  await expect(page.getByRole('tab', { name: 'Phage comparison', exact: true })).toHaveAttribute('aria-selected', 'true');
  await expect(menu.locator('summary')).toHaveAttribute('aria-label', 'Drawing menu for Phage comparison');
  expect(await page.evaluate(() => window.__GBDRAW_HISTORY__.undoLabel())).toBe('Rename drawing');

  await page.getByRole('button', { name: 'Undo', exact: true }).click();
  await expect(drawingTab(page, second)).toHaveAccessibleName(name);
  expect(await tabName(page, second)).toBe(name);
  expect((await shownDrawing(page)).id, 'Undo of a rename shows the same drawing').toBe(second);
});

test('Delete in the tab menu removes the shown drawing as one step that Undo restores; the last drawing cannot be deleted', async ({ page }) => {
  test.setTimeout(180_000);
  await page.setViewportSize({ width: 1440, height: 900 });
  await openApp(page);
  const [first] = await drawingIds(page);
  const second = await takeLinearOffer(page);
  const menu = page.locator('.drawing-bar-menu');
  const deleteItem = menu.getByRole('button', { name: 'Delete drawing', exact: true });

  await menu.locator('summary').click();
  await deleteItem.click();
  await expect(page.getByRole('dialog')).toHaveCount(0);
  expect(await drawingIds(page), 'a drawing without Result or edits is deleted without asking').toEqual([first]);
  expect((await shownDrawing(page)).id).toBe(first);
  expect(await page.evaluate(() => window.__GBDRAW_HISTORY__.undoLabel())).toBe('Delete drawing');

  await menu.locator('summary').click();
  await expect(deleteItem, 'the last drawing has no Delete').toBeDisabled();
  await page.keyboard.press('Escape');

  await page.getByRole('button', { name: 'Undo', exact: true }).click();
  expect(await drawingIds(page)).toEqual([first, second]);
  expect((await shownDrawing(page)).id, 'Undo shows the restored drawing').toBe(second);
});

// PD-OI-090 〔作成〕, R16-2: Duplicate is one project step; its Undo shows the source.
test('Duplicate in the tab menu shows a named copy as one step that Undo takes back to the source', async ({ page }) => {
  test.setTimeout(180_000);
  await page.setViewportSize({ width: 1440, height: 900 });
  await openApp(page);
  const [first] = await drawingIds(page);
  const name = await tabName(page, first);
  const menu = page.locator('.drawing-bar-menu');

  await menu.locator('summary').click();
  await menu.getByRole('button', { name: 'Duplicate drawing', exact: true }).click();
  const ids = await drawingIds(page);
  expect(ids).toHaveLength(2);
  const copy = ids[1];
  expect(ids[0]).toBe(first);
  await expect(drawingTab(page, copy)).toHaveAccessibleName(`${name} copy`);
  expect((await shownDrawing(page)).id, 'the copy is shown').toBe(copy);
  expect(await page.evaluate(() => window.__GBDRAW_HISTORY__.undoLabel())).toBe('Duplicate drawing');

  await page.getByRole('button', { name: 'Undo', exact: true }).click();
  expect(await drawingIds(page)).toEqual([first]);
  expect((await shownDrawing(page)).id, 'Undo shows the source').toBe(first);
});

// PD-OI-090 〔入力ファイル〕, T-B4-6: another drawing chooses a file this Session
// already holds without an upload. The Linear drawing's GenBank uploader offers
// the Circular drawing's file; choosing it binds the same file ID as one step.
test('an uploader offers the files in this Session and binds the chosen one without an upload', async ({ page }) => {
  test.setTimeout(180_000);
  await page.setViewportSize({ width: 1440, height: 900 });
  await openApp(page);
  const [circular] = await drawingIds(page);
  await page.getByLabel('GenBank/DDBJ File', { exact: true }).setInputFiles('tests/test_inputs/NC_001416.gb');
  const boundId = (id, read) => page.evaluate(async ([drawingId, field]) => {
    const drawing = (await import('./js/state.js')).state.drawingOf(drawingId);
    const ref = field === 'c_gb' ? drawing.files.c_gb : drawing.linearSeqs[0].gb;
    return ref ? ref.fileId : null;
  }, [id, read]);
  await expect.poll(() => boundId(circular, 'c_gb')).not.toBeNull();
  const fileId = await boundId(circular, 'c_gb');

  const linear = await takeLinearOffer(page);
  const picker = page.getByRole('combobox', { name: 'Choose GenBank / DDBJ File from this Session', exact: true });
  await expect(picker.locator('option')).toHaveText(['Files in this Session', 'NC_001416.gb']);
  // The store names its files file-1, file-2, ...: the highest stored ID.
  const lastStoredId = () => page.evaluate(async () => {
    const { inputFiles } = (await import('./js/state.js')).state;
    let last = 0;
    for (let id = 1; id <= 500; id += 1) if (inputFiles.has(`file-${id}`)) last = id;
    return last;
  });
  const before = await lastStoredId();
  await picker.selectOption({ label: 'NC_001416.gb' });
  await expect.poll(() => boundId(linear, 'gb'), 'the Linear drawing binds the same file').toBe(fileId);
  expect(await lastStoredId(), 'choosing stored no new file').toBe(before);
  expect(await page.evaluate(() => window.__GBDRAW_HISTORY__.undoLabel())).toBe('Change uploaded file');
  await expect(picker, 'the uploader offers no file it already shows').toHaveCount(0);
});

// D-15-120: a Load of a one-mode Session makes no drawing of the other mode. The bar lists the
// Session's drawings; the Linear offer tab, or the menu's New Circular drawing, makes the other
// mode's drawing at its defaults.
test('a one-mode Session Load lists only its own drawings; the offer tab or New makes the other mode at its defaults', async ({ page }) => {
  test.setTimeout(240_000);
  await page.setViewportSize({ width: 1440, height: 900 });
  await openApp(page);
  const atDefaults = (id) => page.evaluate(async (drawingId) => {
    const { state, createDrawingState } = await import('./js/state.js');
    const drawing = state.drawingOf(drawingId);
    const fresh = createDrawingState(drawing.id, drawing.mode);
    return JSON.stringify([drawing.form, drawing.adv]) === JSON.stringify([fresh.form, fresh.adv]);
  }, id);

  const linearOnly = await loadGallerySession(page, 'lambda_basic_linear.gbdraw-session.json');
  expect(linearOnly.status, JSON.stringify(linearOnly.error)).toBe('ok');
  expect(await drawingIds(page)).toEqual(['linear']);
  await expect(page.locator('.drawing-bar [role="tab"]'), 'the bar lists the Linear drawing only').toHaveCount(1);
  await expect(offerTab(page)).toHaveCount(0);
  const menu = page.locator('.drawing-bar-menu');
  await menu.locator('summary').click();
  await menu.getByRole('button', { name: 'New Circular drawing', exact: true }).click();
  await expect.poll(() => drawingIds(page)).toEqual(['linear', 'circular']);
  expect(await shownDrawing(page)).toEqual({ id: 'circular', kind: 'circular' });
  expect(await atDefaults('circular'), 'the new Circular drawing is at its defaults').toBe(true);

  const circularOnly = await loadGallerySession(page, 'HmmtDNA_basic_circular.gbdraw-session.json');
  expect(circularOnly.status, JSON.stringify(circularOnly.error)).toBe('ok');
  expect(await drawingIds(page)).toEqual(['circular']);
  await expect(offerTab(page)).toHaveAccessibleName('Linear');
  await takeLinearOffer(page);
  expect(await drawingIds(page)).toEqual(['circular', 'linear']);
  expect(await atDefaults('linear'), 'the offer tab makes the Linear drawing at its defaults').toBe(true);
});

test('Save then a fresh Load keeps the drawings, their names, and their tab order', async ({ page }) => {
  test.setTimeout(240_000);
  await page.setViewportSize({ width: 1440, height: 900 });
  await openApp(page);
  await expect(offerTab(page)).toHaveAccessibleName('Linear');
  await takeLinearOffer(page);
  await expect(offerTab(page)).toHaveCount(0);
  const menu = page.locator('.drawing-bar-menu');
  await menu.locator('summary').click();
  await menu.getByRole('button', { name: 'New Circular drawing', exact: true }).click();
  await expect.poll(() => drawingIds(page)).toEqual(['circular', 'linear', 'circular-2']);
  expect((await shownDrawing(page)).id).toBe('circular-2');
  // Save keeps a drawing with a changed setting (PD-OI-090): each drawing gets a title.
  await page.evaluate(async () => {
    const { state } = await import('./js/state.js');
    for (const drawing of state.drawingList()) drawing.form.plot_title = `Title of ${drawing.id}`;
  });
  const dialog = page.getByRole('dialog', { name: 'Rename drawing' });
  await menu.locator('summary').click();
  await menu.getByRole('button', { name: 'Rename…', exact: true }).click();
  await dialog.getByRole('textbox', { name: 'Drawing name' }).fill('Plasmid map');
  await dialog.getByRole('button', { name: 'Rename', exact: true }).click();
  await menu.locator('summary').click();
  await menu.getByRole('button', { name: 'Move left', exact: true }).click();
  const tabsBefore = await page.evaluate(async () => (await import('./js/state.js')).state.drawingTabs.value
    .map(({ id, name }) => ({ id, name })));
  expect(tabsBefore).toEqual([
    { id: 'circular', name: 'Circular' }, { id: 'circular-2', name: 'Plasmid map' }, { id: 'linear', name: 'Linear' }]);

  // A titled Session saves without the title question.
  await page.evaluate(() => { window.__GBDRAW_APP__.sessionTitle = 'drawings'; });
  const downloadPromise = page.waitForEvent('download', { timeout: 120_000 });
  expect(await evaluateWithRetainedPromise(page, () => window.__GBDRAW_APP__.saveSessionWithTitle()))
    .toMatchObject({ status: 'saved' });
  const saved = gunzipSync(readFileSync(await (await downloadPromise).path())).toString('utf8');
  const fresh = await page.context().newPage();
  await openApp(fresh);
  expect(await drawingIds(fresh)).toEqual(['circular']);
  const outcome = await evaluateWithRetainedPromise(fresh, async (text) => {
    const file = new File([text], 'drawings.gbdraw-session.json', { type: 'application/json' });
    return window.__GBDRAW_APP__.importSession({ target: { files: [file], value: 'selected' } });
  }, saved);
  expect(outcome.status).toBe('ok');
  expect(await fresh.evaluate(async () => (await import('./js/state.js')).state.drawingTabs.value
    .map(({ id, name }) => ({ id, name })))).toEqual(tabsBefore);
  expect((await shownDrawing(fresh)).id, 'the drawing shown at Save is shown').toBe('circular-2');
  await expect(drawingTab(fresh, 'circular-2')).toHaveAccessibleName('Plasmid map');
});
