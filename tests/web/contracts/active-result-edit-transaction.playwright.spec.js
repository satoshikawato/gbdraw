const { test, expect } = require('@playwright/test');
const { readFileSync } = require('node:fs');
const { join, resolve } = require('node:path');
const { gunzipSync } = require('node:zlib');
const { openApp } = require('../helpers/app-lifecycle.cjs');

const repoRoot = resolve(process.env.GBDRAW_REPO || process.cwd());
const sourceSessionPath = join(
  repoRoot,
  'gbdraw',
  'web',
  'gallery',
  'sessions',
  'HmmtDNA_basic_circular.gbdraw-session.json'
);
const BEFORE_COLOR = '#e8b441';
const AFTER_COLOR = '#2f7fc1';
const CONTINUATION_COLOR = '#6a3d9a';
const TARGET_CAPTION = 'tRNA';

test.describe('active Result Feature fill transaction', () => {
  test.describe.configure({ retries: 0 });

  const installLifecycleEvidence = async (page) => {
    await page.addInitScript(() => {
      const events = [];
      window.__GBDRAW_C1_LIFECYCLE__ = events;
      window.__GBDRAW_TEST_HOOKS__ = {
        onSessionLifecycleEvent(event) {
          events.push({ ...event });
        }
      };
    });
  };

  const observeErrors = (page) => {
    const diagnostics = { pageErrors: [], consoleErrors: [], dialogs: [] };
    page.on('pageerror', (error) => {
      diagnostics.pageErrors.push(String(error?.message || error));
    });
    page.on('console', (message) => {
      if (message.type() === 'error') diagnostics.consoleErrors.push(message.text());
    });
    page.on('dialog', async (dialog) => {
      diagnostics.dialogs.push({ type: dialog.type(), message: dialog.message() });
      await dialog.accept();
    });
    return diagnostics;
  };

  const openObservedApp = async (page) => {
    const diagnostics = observeErrors(page);
    await installLifecycleEvidence(page);
    await openApp(page);
    return diagnostics;
  };

  const loadSessionThroughUi = async (page, filePath) => {
    await page.evaluate(() => {
      window.__GBDRAW_C1_LIFECYCLE__.length = 0;
    });
    const chooserPromise = page.waitForEvent('filechooser');
    await page.getByRole('button', { name: 'Load Session', exact: true }).click();
    const chooser = await chooserPromise;
    await chooser.setFiles(filePath);
    await page.waitForFunction(() => window.__GBDRAW_C1_LIFECYCLE__.some(
      (event) => event.name === 'interactiveReady' && event.status === 'success'
    ), null, { timeout: 180_000 });
    await expect(page.locator('.origin-top svg')).toBeVisible();
  };

  const collectTargetInventory = (page) => page.evaluate(async ({ caption }) => {
    const app = window.__GBDRAW_APP__;
    const svg = document.querySelector('.origin-top svg');
    const { getFeatureFillElements } = await import(
      '/gbdraw/web/js/app/feature-editor/svg-actions.js'
    );
    const { featureOverrideKey } = await import(
      '/gbdraw/web/js/services/feature-override-identity.js'
    );
    const rendered = (feature) => {
      const id = String(feature?.svg_id || '').trim();
      return id && getFeatureFillElements(svg, id).length > 0;
    };
    const targets = app.extractedFeatures.filter(
      (feature) => feature.type === caption && rendered(feature)
    );
    const nonTarget = app.extractedFeatures.find(
      (feature) => feature.type !== caption && rendered(feature)
    );
    if (targets.length < 2 || !nonTarget) {
      throw new Error('The fixture did not expose the required tRNA scope and non-target.');
    }
    const inventoryEntry = (feature) => ({
      id: String(feature.svg_id),
      key: featureOverrideKey(feature),
      type: String(feature.type),
      start: Number(feature.start),
      end: Number(feature.end)
    });
    return {
      targets: targets.map(inventoryEntry),
      nonTarget: inventoryEntry(nonTarget)
    };
  }, { caption: TARGET_CAPTION });

  const inspectState = (page, inventory) => page.evaluate(async ({ targetInventory, caption }) => {
    const app = window.__GBDRAW_APP__;
    const { state } = await import('/gbdraw/web/js/state.js');
    const mountedSvg = document.querySelector('.origin-top svg');
    const selectedResult = app.results[app.selectedResultIndex] || null;
    const selectedSvg = new DOMParser().parseFromString(
      String(selectedResult?.content || ''),
      'image/svg+xml'
    ).documentElement;
    const { getFeatureFillElements } = await import(
      '/gbdraw/web/js/app/feature-editor/svg-actions.js'
    );
    const normalizedFills = (root, id) => getFeatureFillElements(root, id)
      .map((element) => String(element.getAttribute('fill') || '').toLowerCase());
    const legendFill = (root) => {
      const group = root.querySelector(`g[data-legend-key="${CSS.escape(caption)}"]`);
      return String(Array.from(group?.querySelectorAll('path') || []).find((path) => {
        const fill = path.getAttribute('fill');
        return fill && fill !== 'none' && !fill.startsWith('url(');
      })?.getAttribute('fill') || '').toLowerCase();
    };
    const inspectSvg = (root) => ({
      targetFills: Object.fromEntries(targetInventory.targets.map(({ id }) => [
        id,
        normalizedFills(root, id)
      ])),
      nonTargetFills: normalizedFills(root, targetInventory.nonTarget.id),
      legendFill: legendFill(root)
    });
    const plain = (value) => JSON.parse(JSON.stringify(value ?? null));
    const history = window.__GBDRAW_HISTORY__;
    return {
      mounted: inspectSvg(mountedSvg),
      selectedResult: inspectSvg(selectedSvg),
      canonical: {
        targetOverrides: Object.fromEntries(targetInventory.targets.map(({ key }) => [
          key,
          plain(app.featureColorOverrides[key])
        ])),
        targetRules: plain(app.manualSpecificRules.filter((rule) => (
          rule.feat === caption && String(rule.cap || '') === caption
        ))),
        legendOverride: String(state.legendColorOverrides[caption] || '').toLowerCase(),
        legendEntryColor: String(
          app.legendEntries.find((entry) => entry.caption === caption)?.color || ''
        ).toLowerCase()
      },
      history: {
        undoCount: history.getUndoCount(),
        redoCount: history.getRedoCount(),
        undoLabel: String(history.undoLabel?.() || ''),
        redoLabel: String(history.redoLabel?.() || '')
      },
      feedback: app.errorLog ? plain(app.errorLog) : null
    };
  }, { targetInventory: inventory, caption: TARGET_CAPTION });

  const expectFillState = (state, inventory, color, nonTargetBaseline) => {
    for (const { id } of inventory.targets) {
      expect(state.mounted.targetFills[id], `mounted ${id}`).not.toHaveLength(0);
      expect(state.mounted.targetFills[id].every((fill) => fill === color)).toBe(true);
      expect(state.selectedResult.targetFills[id], `selected Result ${id}`).toEqual(
        state.mounted.targetFills[id]
      );
    }
    expect(state.mounted.nonTargetFills).toEqual(nonTargetBaseline);
    expect(state.selectedResult.nonTargetFills).toEqual(nonTargetBaseline);
    expect(state.mounted.legendFill).toBe(color);
    expect(state.selectedResult.legendFill).toBe(color);
    expect(state.canonical.legendEntryColor).toBe(color);
    expect(state.feedback).toBeNull();
  };

  const expectAcceptedCanonicalState = (state, inventory) => {
    for (const { key } of inventory.targets) {
      expect(state.canonical.targetOverrides[key]).toEqual({
        color: AFTER_COLOR,
        caption: TARGET_CAPTION
      });
    }
    expect(state.canonical.targetRules).toHaveLength(inventory.targets.length);
    expect(state.canonical.targetRules.every((rule) => (
      String(rule.qual || '').toLowerCase() === 'hash'
      && String(rule.color || '').toLowerCase() === AFTER_COLOR
      && rule.cap === TARGET_CAPTION
    ))).toBe(true);
    expect(state.canonical.legendOverride).toBe(AFTER_COLOR);
  };

  const inspectSvgText = (page, svgText, inventory) => page.evaluate(
    async ({ text, targetInventory, caption }) => {
      const root = new DOMParser().parseFromString(text, 'image/svg+xml').documentElement;
      const { getFeatureFillElements } = await import(
        '/gbdraw/web/js/app/feature-editor/svg-actions.js'
      );
      const fills = (id) => getFeatureFillElements(root, id)
        .map((element) => String(element.getAttribute('fill') || '').toLowerCase());
      const legendGroup = root.querySelector(
        `g[data-legend-key="${CSS.escape(caption)}"]`
      );
      const legendFill = String(Array.from(legendGroup?.querySelectorAll('path') || [])
        .find((path) => {
          const fill = path.getAttribute('fill');
          return fill && fill !== 'none' && !fill.startsWith('url(');
        })?.getAttribute('fill') || '').toLowerCase();
      return {
        targetFills: Object.fromEntries(targetInventory.targets.map(({ id }) => [id, fills(id)])),
        nonTargetFills: fills(targetInventory.nonTarget.id),
        legendFill
      };
    },
    { text: svgText, targetInventory: inventory, caption: TARGET_CAPTION }
  );

  const readSavedSession = (filePath) => {
    const bytes = readFileSync(filePath);
    const text = bytes[0] === 0x1f && bytes[1] === 0x8b
      ? gunzipSync(bytes).toString('utf8')
      : bytes.toString('utf8');
    return JSON.parse(text);
  };

  const assertNoUnexpectedErrors = async (page, diagnostics) => {
    expect(diagnostics.pageErrors).toEqual([]);
    expect(diagnostics.consoleErrors).toEqual([]);
    expect(diagnostics.dialogs.every(({ message }) => (
      message === 'Session loaded successfully!'
    ))).toBe(true);
    expect(await page.evaluate(() => window.__GBDRAW_APP__.errorLog || null)).toBeNull();
  };

  test('visible tRNA fill scope remains coherent through History, Session, Generate, and SVG export', async ({
    page,
    browser
  }, testInfo) => {
    test.setTimeout(600_000);
    expect(testInfo.retry).toBe(0);

    const initialDiagnostics = await openObservedApp(page);
    await loadSessionThroughUi(page, sourceSessionPath);
    const inventory = await collectTargetInventory(page);
    const before = await inspectState(page, inventory);
    const nonTargetBaseline = before.mounted.nonTargetFills;
    expectFillState(before, inventory, BEFORE_COLOR, nonTargetBaseline);
    expect(before.history).toMatchObject({ undoCount: 0, redoCount: 0 });
    expect(before.canonical.targetRules).toEqual([]);
    expect(Object.values(before.canonical.targetOverrides).every((value) => value === null))
      .toBe(true);

    await page.locator('.drawer-toggle').click();
    const drawer = page.locator('.right-drawer');
    await expect(drawer.getByText(/^Features \(\d+\)$/)).toBeVisible();
    await drawer.getByPlaceholder('Search by feature or annotation...').fill(TARGET_CAPTION);
    await expect(drawer.getByText(`Features (${inventory.targets.length})`, { exact: true }))
      .toBeVisible();
    const selectedTarget = inventory.targets[0];
    const targetRow = drawer
      .locator(`span[title="${selectedTarget.start}..${selectedTarget.end}"]`)
      .locator('..');
    await targetRow.getByRole('button', { name: 'Edit', exact: true }).click();

    const popup = page.getByRole('dialog', { name: /Feature details:/ });
    await expect(popup).toBeVisible();
    const fillPicker = popup.getByLabel('Feature fill color', { exact: true }).first();
    await expect(fillPicker).toHaveValue(BEFORE_COLOR);
    await fillPicker.fill(AFTER_COLOR);
    const scopeDialog = page.getByRole('heading', { name: 'Color Change Scope' }).locator('..');
    await expect(scopeDialog).toBeVisible();
    const applyScopeButton = scopeDialog.getByRole('button').filter({
      hasText: `Apply to all "${TARGET_CAPTION}" (${inventory.targets.length})`
    }).last();
    await expect(applyScopeButton).toBeVisible();
    expect(await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount())).toBe(0);
    await applyScopeButton.click();
    await expect.poll(
      () => page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount())
    ).toBe(1);

    const applied = await inspectState(page, inventory);
    expectFillState(applied, inventory, AFTER_COLOR, nonTargetBaseline);
    expectAcceptedCanonicalState(applied, inventory);
    expect(applied.history).toMatchObject({
      undoCount: 1,
      redoCount: 0,
      undoLabel: 'Change feature color'
    });

    await page.getByRole('button', { name: 'Undo', exact: true }).click();
    await expect.poll(
      () => page.evaluate(() => window.__GBDRAW_HISTORY__.getRedoCount())
    ).toBe(1);
    const undone = await inspectState(page, inventory);
    expectFillState(undone, inventory, BEFORE_COLOR, nonTargetBaseline);
    expect(undone.history).toMatchObject({
      undoCount: 0,
      redoCount: 1,
      redoLabel: 'Change feature color'
    });
    expect(undone.canonical.targetRules).toEqual([]);
    expect(Object.values(undone.canonical.targetOverrides).every((value) => value === null))
      .toBe(true);
    expect(undone.canonical.legendOverride).toBe('');

    await page.getByRole('button', { name: 'Redo', exact: true }).click();
    await expect.poll(
      () => page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount())
    ).toBe(1);
    const redone = await inspectState(page, inventory);
    expectFillState(redone, inventory, AFTER_COLOR, nonTargetBaseline);
    expectAcceptedCanonicalState(redone, inventory);
    expect(redone.history).toMatchObject({ undoCount: 1, redoCount: 0 });

    const sessionDownloadPromise = page.waitForEvent('download', { timeout: 180_000 });
    await page.getByRole('button', { name: 'Save Session', exact: true }).click();
    const sessionDownload = await sessionDownloadPromise;
    const savedSessionPath = await sessionDownload.path();
    expect(savedSessionPath).toBeTruthy();
    const savedSession = readSavedSession(savedSessionPath);
    expect(savedSession).toMatchObject({
      format: 'gbdraw-session',
      version: 40,
      renderRequest: { schema: 6 },
      editorState: {
        legend: { colorOverrides: { [TARGET_CAPTION]: AFTER_COLOR } }
      }
    });
    for (const { key } of inventory.targets) {
      expect(savedSession.features.featureColorOverrides[key]).toEqual({
        color: AFTER_COLOR,
        caption: TARGET_CAPTION
      });
    }
    const savedRules = savedSession.config.rules.filter((rule) => (
      rule.feat === TARGET_CAPTION && rule.cap === TARGET_CAPTION
    ));
    expect(savedRules).toHaveLength(inventory.targets.length);
    expect(savedSession.results[savedSession.ui.selectedResultIndex].content)
      .toContain(AFTER_COLOR);
    await assertNoUnexpectedErrors(page, initialDiagnostics);

    const freshContext = await browser.newContext({ acceptDownloads: true });
    const freshPage = await freshContext.newPage();
    try {
      const freshDiagnostics = await openObservedApp(freshPage);
      await loadSessionThroughUi(freshPage, savedSessionPath);
      const freshlyLoaded = await inspectState(freshPage, inventory);
      expectFillState(freshlyLoaded, inventory, AFTER_COLOR, nonTargetBaseline);
      expectAcceptedCanonicalState(freshlyLoaded, inventory);
      expect(freshlyLoaded.history).toMatchObject({ undoCount: 0, redoCount: 0 });

      await freshPage.evaluate(() => {
        window.__GBDRAW_C1_LIFECYCLE__.length = 0;
      });
      const generateButton = freshPage.getByRole('button', {
        name: 'Generate Diagram',
        exact: true
      });
      await expect(generateButton).toBeEnabled();
      await generateButton.click();
      await freshPage.waitForFunction(() => (
        window.__GBDRAW_APP__?.processing === false
        && window.__GBDRAW_C1_LIFECYCLE__.some(
          (event) => event.name === 'preview-result-commit-end'
        )
      ), null, { timeout: 300_000 });
      const generated = await inspectState(freshPage, inventory);
      expectFillState(generated, inventory, AFTER_COLOR, nonTargetBaseline);
      expectAcceptedCanonicalState(generated, inventory);

      const svgDownloadPromise = freshPage.waitForEvent('download', { timeout: 180_000 });
      await freshPage.getByRole('button', { name: 'SVG', exact: true }).click();
      const svgDownload = await svgDownloadPromise;
      const exportedSvg = readFileSync(await svgDownload.path(), 'utf8');
      const exported = await inspectSvgText(freshPage, exportedSvg, inventory);
      for (const { id } of inventory.targets) {
        expect(exported.targetFills[id]).not.toHaveLength(0);
        expect(exported.targetFills[id].every((fill) => fill === AFTER_COLOR)).toBe(true);
      }
      expect(exported.nonTargetFills).toEqual(nonTargetBaseline);
      expect(exported.legendFill).toBe(AFTER_COLOR);

      const historyBeforeContinuation = generated.history;
      await freshPage.locator('.drawer-toggle').click();
      const freshDrawer = freshPage.locator('.right-drawer');
      await freshDrawer.getByPlaceholder('Search by feature or annotation...').fill(TARGET_CAPTION);
      const freshTargetRow = freshDrawer
        .locator(`span[title="${selectedTarget.start}..${selectedTarget.end}"]`)
        .locator('..');
      await freshTargetRow.getByRole('button', { name: 'Edit', exact: true }).click();
      const continuationPicker = freshPage
        .getByRole('dialog', { name: /Feature details:/ })
        .getByLabel('Feature fill color', { exact: true })
        .first();
      await expect(continuationPicker).toHaveValue(AFTER_COLOR);
      await continuationPicker.fill(CONTINUATION_COLOR);
      const continuationScope = freshPage
        .getByRole('heading', { name: 'Color Change Scope' })
        .locator('..');
      await expect(continuationScope).toBeVisible();
      await continuationScope.getByRole('button', { name: 'Cancel', exact: true }).click();
      const continued = await inspectState(freshPage, inventory);
      expectFillState(continued, inventory, AFTER_COLOR, nonTargetBaseline);
      expectAcceptedCanonicalState(continued, inventory);
      expect(continued.history).toEqual(historyBeforeContinuation);
      await assertNoUnexpectedErrors(freshPage, freshDiagnostics);

      await testInfo.attach('active-result-edit-transaction.json', {
        body: Buffer.from(JSON.stringify({
          targetCaption: TARGET_CAPTION,
          targetCount: inventory.targets.length,
          nonTarget: inventory.nonTarget,
          before: { history: before.history, canonical: before.canonical },
          applied: { history: applied.history, canonical: applied.canonical },
          undone: { history: undone.history, canonical: undone.canonical },
          redone: { history: redone.history, canonical: redone.canonical },
          freshlyLoaded: {
            history: freshlyLoaded.history,
            canonical: freshlyLoaded.canonical
          },
          generated: { history: generated.history, canonical: generated.canonical },
          exported
        }, null, 2)),
        contentType: 'application/json'
      });
    } finally {
      await freshContext.close();
    }
  });
});
