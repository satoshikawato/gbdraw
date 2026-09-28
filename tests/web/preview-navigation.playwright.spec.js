const { test, expect } = require('@playwright/test');
const { join } = require('node:path');
const { writeFile } = require('node:fs/promises');
const { openApp, getDiagramWorkerActivity } = require('./helpers/app-lifecycle.cjs');

const geometry = (page) => page.evaluate(() => {
  const app = window.__GBDRAW_APP__;
  const wrapper = app.svgContainer;
  const svg = wrapper.querySelector('svg');
  const container = app.canvasContainerRef;
  const matrix = svg.getScreenCTM();
  const point = (x, y) => {
    const p = new DOMPoint(x, y).matrixTransform(matrix);
    return { x: p.x, y: p.y };
  };
  const bounds = wrapper.getBoundingClientRect();
  return {
    zoom: app.zoom, pan: { ...app.canvasPan },
    scroll: { x: container.scrollLeft, y: container.scrollTop },
    points: [point(0, 0), point(100, 100)],
    anchor: { x: bounds.x + bounds.width / 2, y: bounds.y },
    panning: app.isPanning
  };
});

const settle = async (page) => {
  await expect.poll(() => page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    const transform = new DOMMatrix(getComputedStyle(app.svgContainer).transform);
    return Math.max(Math.abs(transform.a - app.zoom),
      Math.abs(transform.e - app.canvasPan.x), Math.abs(transform.f - app.canvasPan.y));
  })).toBeLessThan(0.001);
};

const backgroundPoint = (page) => page.evaluate(() => {
  const app = window.__GBDRAW_APP__;
  const container = app.canvasContainerRef;
  const rect = container.getBoundingClientRect();
  const svg = app.svgContainer.querySelector('svg');
  for (let y = rect.bottom - 160; y > rect.y + 30; y -= 30) {
    for (let x = rect.x + 50; x < rect.right - 40; x += 60) {
      const target = document.elementFromPoint(x, y);
      if (target === container || target === app.svgContainer || target === svg) return { x, y };
    }
  }
  throw new Error('No visible non-editing background point.');
});

const move = async (page, point, dx, dy) => {
  await page.mouse.move(point.x, point.y);
  await page.mouse.down();
  await page.mouse.move(point.x + dx, point.y + dy, { steps: 10 });
  await page.mouse.up();
  await settle(page);
};

const expectTranslation = (before, after, dx, dy) => {
  for (let i = 0; i < before.points.length; i += 1) {
    expect(after.points[i].x - before.points[i].x).toBeCloseTo(dx, 1);
    expect(after.points[i].y - before.points[i].y).toBeCloseTo(dy, 1);
  }
  expect(after.panning).toBe(false);
};

test('background pan preserves screen displacement after text selection and zoom', async ({ page }, testInfo) => {
  test.setTimeout(180000);
  await openApp(page);
  page.on('dialog', (dialog) => dialog.dismiss());
  await page.locator('input[type=file][accept*="application/json"][accept*="application/gzip"]')
    .setInputFiles(join(process.cwd(), 'gbdraw/web/gallery/sessions/BGC0000708-BGC0000713.gbdraw-session.json'));
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.results.length)).toBe(1);
  await settle(page);
  const committed = await page.evaluate(async () => {
    const { getCommittedCanonicalSession } = await import('./js/services/config.js');
    return JSON.stringify(getCommittedCanonicalSession());
  });
  const evidence = [];
  for (const zoom of [0.3, 1, 2.8]) {
    await page.getByRole('button', { name: 'Reset layout', exact: true }).click();
    const direction = zoom < 1 ? 'Zoom out' : 'Zoom in';
    for (let i = 0; i < Math.round(Math.abs(zoom - 1) * 10); i += 1) {
      await page.getByRole('button', { name: direction, exact: true }).click();
    }
    await settle(page);
    const before = await geometry(page);
    const point = await backgroundPoint(page);
    await move(page, point, 60, 80);
    const after = await geometry(page);
    evidence.push({ zoom, point, before, after });
    await writeFile(testInfo.outputPath('navigation.json'), JSON.stringify(evidence, null, 2));
    expectTranslation(before, after, 60, 80);
    // Preserve the existing top-center zoom anchor, including after panning.
    await page.mouse.move(point.x, point.y);
    await page.mouse.wheel(0, -100);
    await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.zoom)).toBeCloseTo(zoom + 0.1, 5);
    await settle(page);
    const zoomed = await geometry(page);
    expect(zoomed.anchor.x).toBeCloseTo(after.anchor.x, 1);
    expect(zoomed.anchor.y).toBeCloseTo(after.anchor.y, 1);
    await page.mouse.wheel(0, 100);
    await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.zoom)).toBeCloseTo(zoom, 5);
    await settle(page);
    expectTranslation(after, await geometry(page), 0, 0);
  }
  await page.getByRole('button', { name: 'Reset layout', exact: true }).click();
  await settle(page);
  expect(await geometry(page)).toMatchObject({ zoom: 1, pan: { x: 0, y: 0 }, panning: false });
  expect(await page.evaluate(async () => {
    const { getCommittedCanonicalSession } = await import('./js/services/config.js');
    return JSON.stringify(getCommittedCanonicalSession());
  })).toBe(committed);
  expect(await getDiagramWorkerActivity(page)).toMatchObject({ constructions: 0 });
});

test('preview pan leaves feature and match gestures available', async ({ page }) => {
  test.setTimeout(180000);
  // Keep the wide comparison's editing targets visible at its existing zoom anchor.
  await page.setViewportSize({ width: 2520, height: 1327 });
  await openApp(page);
  page.on('dialog', (dialog) => dialog.dismiss());
  await page.locator('input[type=file][accept*="application/json"][accept*="application/gzip"]')
    .setInputFiles(join(process.cwd(), 'gbdraw/web/gallery/sessions/BGC0000708-BGC0000713.gbdraw-session.json'));
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.results.length)).toBe(1);
  for (let i = 0; i < 7; i += 1) await page.getByRole('button', { name: 'Zoom out', exact: true }).click();
  await settle(page);
  expect(await getDiagramWorkerActivity(page)).toMatchObject({ constructions: 0 });

  for (const enabled of [false, true]) {
    await page.evaluate((enabled) => { window.__GBDRAW_APP__.layoutRepositionMode = enabled; }, enabled);
    for (const selector of ['[data-gbdraw-feature-id]', '[data-gbdraw-pairwise-match-id]']) {
      const point = await page.evaluate((selector) => {
        for (const element of window.__GBDRAW_APP__.svgContainer.querySelectorAll(selector)) {
          const bounds = element.getBoundingClientRect();
          const x = bounds.x + bounds.width / 2;
          const y = bounds.y + bounds.height / 2;
          if (document.elementFromPoint(x, y)?.closest(selector)) return { x, y };
        }
        throw new Error(`No visible gesture target: ${selector}`);
      }, selector);
      const before = await geometry(page);
      await move(page, point, 10, 10);
      expectTranslation(before, await geometry(page), 0, 0);
      await page.mouse.click(point.x, point.y);
      await page.waitForFunction(() => !window.__GBDRAW_APP__.ruleMatchingPending);
      const selectionKey = selector.includes('pairwise') ? 'clickedPairwiseMatch' : 'clickedFeature';
      await expect.poll(() => page.evaluate((key) => Boolean(window.__GBDRAW_APP__[key]), selectionKey)).toBe(true);
      await page.keyboard.press('Escape');
      if (selector === '[data-gbdraw-feature-id]') {
        const id = await page.evaluate(({ x, y }) => document.elementFromPoint(x, y)
          .closest('[data-gbdraw-feature-id]').getAttribute('data-gbdraw-rendered-feature-id')
          || document.elementFromPoint(x, y).closest('[data-gbdraw-feature-id]').id, point);
        for (const modifier of ['Control', 'Shift']) {
          await page.keyboard.down(modifier);
          await page.mouse.click(point.x, point.y);
          await page.keyboard.up(modifier);
          expect(await page.evaluate((id) => window.__GBDRAW_APP__.selectedFeatureIds.has(id), id)).toBe(true);
          expectTranslation(before, await geometry(page), 0, 0);
        }
        await page.evaluate(() => window.__GBDRAW_APP__.clearFeatureSelection());
      }
    }
  }
  await page.evaluate(() => { window.__GBDRAW_APP__.layoutRepositionMode = false; });

  await page.locator('.drawer-toggle').click();
  await expect(page.locator('.right-drawer')).toBeVisible();
  const before = await geometry(page);
  await move(page, await backgroundPoint(page), 30, 20);
  expectTranslation(before, await geometry(page), 30, 20);
  await page.getByRole('button', { name: 'Reset zoom', exact: true }).focus();
  await page.keyboard.press('Enter');
  await settle(page);
  expect(await geometry(page)).toMatchObject({ zoom: 1, pan: { x: 30, y: 20 } });
  await page.getByRole('button', { name: 'Reset layout', exact: true }).click();
  await settle(page);
  expect(await geometry(page)).toMatchObject({ zoom: 1, pan: { x: 0, y: 0 }, panning: false });
  // Feature editing prepares Python rule matches; later navigation reuses that worker.
  expect(await getDiagramWorkerActivity(page)).toMatchObject({ constructions: 1 });
});

test('docked search and controls remain reachable across Preview sizes and Editor transitions', async ({ page }, testInfo) => {
  test.setTimeout(120000);
  await page.setViewportSize({ width: 1440, height: 844 });
  await openApp(page);
  await page.locator('input[accept^=".json,"]').setInputFiles(
    join(process.cwd(), 'gbdraw/web/gallery/sessions/HmmtDNA_basic_circular.gbdraw-session.json')
  );
  await page.waitForFunction(() => window.__GBDRAW_APP__.results.length && window.__GBDRAW_APP__.extractedFeatures.length);
  const search = page.getByRole('searchbox', { name: 'Search features', exact: true });
  await search.fill('tRNA');
  await search.press('Enter');
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.previewFeatureSearchMatches.length)).toBeGreaterThan(1);
  await page.getByRole('button', { name: 'Next match' }).click();
  const active = await page.evaluate(() => window.__GBDRAW_APP__.previewFeatureSearchActiveIndex);
  await search.focus();
  await page.evaluate(() => {
    const preview = document.querySelector('[aria-label="Result Preview"]');
    window.__s03Nodes = {
      search: preview.querySelector('.preview-feature-search'),
      canvas: preview.querySelector('.preview-canvas'),
      svg: window.__GBDRAW_APP__.svgContainer,
      editor: preview.querySelector('.right-drawer'),
      controls: preview.querySelector('.preview-controls')
    };
  });
  const observations = [];
  const matrix = [
    [1440, 844], [1024, 740], [900, 844], [768, 740], [390, 844],
    [390, 740], [390, 480], [1024, 480], [900, 740], [768, 480]
  ];
  for (const [width, height] of matrix) {
    await page.setViewportSize({ width, height });
    const operations = await page.locator(
      '.preview-feature-search input:not(:disabled), .preview-feature-search select:not(:disabled), '
      + '.preview-feature-search button:not(:disabled), .preview-controls button:not(:disabled)'
    ).all();
    for (const operation of operations) {
      await operation.evaluate((element) => element.scrollIntoView({ block: 'center' }));
      const result = await operation.evaluate((element) => {
        const box = element.getBoundingClientRect();
        return element.contains(document.elementFromPoint(box.left + box.width / 2, box.top + box.height / 2));
      });
      expect(result, `${width}x${height}: ${await operation.getAttribute('aria-label') || await operation.textContent()} unreachable`).toBe(true);
    }
    const record = await page.evaluate(() => {
      const preview = document.querySelector('[aria-label="Result Preview"]');
      const { search, canvas, svg, editor, controls } = window.__s03Nodes;
      const workspace = preview.querySelector('.preview-workspace');
      const box = (element) => {
        const { left, top, right, bottom, height } = element.getBoundingClientRect();
        return { left, top, right, bottom, height };
      };
      return {
        size: [innerWidth, innerHeight], containerWidth: preview.clientWidth,
        search: box(search), workspace: box(workspace), controls: box(controls),
        scroll: { page: scrollY, pane: document.querySelector('.result-pane').scrollTop },
        sameNodes: search === preview.querySelector('.preview-feature-search')
          && canvas === preview.querySelector('.preview-canvas')
          && svg === window.__GBDRAW_APP__.svgContainer
          && editor === preview.querySelector('.right-drawer')
          && controls === preview.querySelector('.preview-controls'),
        rowParents: search.parentElement === workspace.parentElement
          && controls.parentElement === workspace.parentElement
          && editor.parentElement === workspace,
        query: window.__GBDRAW_APP__.previewFeatureSearchQuery,
        active: window.__GBDRAW_APP__.previewFeatureSearchActiveIndex
      };
    });
    record.reachableOperations = operations.length;
    observations.push(record);
    expect(record.sameNodes).toBe(true);
    expect(record.rowParents).toBe(true);
    expect(record.search.bottom).toBeLessThanOrEqual(record.workspace.top + 1);
    expect(record.workspace.bottom).toBeLessThanOrEqual(record.controls.top + 1);
    expect(record.workspace.height).toBeGreaterThanOrEqual(200);
    expect(record.query).toBe('tRNA');
    expect(record.active).toBe(active);
  }
  await page.setViewportSize({ width: 1024, height: 740 });
  const beforeResize = await page.locator('[aria-label="Result Preview"]').evaluate((element) => element.clientWidth);
  const handle = await page.locator('.resize-handle').boundingBox();
  await page.mouse.move(handle.x + handle.width / 2, handle.y + handle.height / 2);
  await page.mouse.down();
  await page.mouse.move(handle.x + handle.width / 2 + 45, handle.y + handle.height / 2, { steps: 5 });
  await page.mouse.up();
  const afterResize = await page.locator('[aria-label="Result Preview"]').evaluate((element) => element.clientWidth);
  expect(beforeResize).toBeGreaterThan(640);
  expect(afterResize).toBeLessThanOrEqual(640);
  observations.push({ settingsResize: [beforeResize, afterResize] });
  await search.focus();
  await page.evaluate(() => { window.__GBDRAW_APP__.sidebarWidth -= 45; });
  expect(await search.evaluate((element) => document.activeElement === element)).toBe(true);
  await page.evaluate(() => { window.__GBDRAW_APP__.sidebarWidth += 45; });
  expect(await search.evaluate((element) => document.activeElement === element)).toBe(true);
  await page.evaluate(() => window.__GBDRAW_APP__.openRightDrawerTab('legend'));
  await expect(page.locator('.right-drawer')).toHaveAttribute('aria-hidden', 'false');
  expect(await search.evaluate((element) => document.activeElement === element)).toBe(true);
  await page.locator('.right-drawer').getByRole('button', { name: /Features/ }).click();
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.rightDrawerTab)).toBe('features');
  await page.locator('.right-drawer').getByRole('button', { name: 'Close editor', exact: true }).click();
  await expect(page.locator('.right-drawer')).toHaveAttribute('aria-hidden', 'true');
  await search.focus();
  await page.evaluate(() => window.__GBDRAW_APP__.openRightDrawerTab('legend'));
  await page.keyboard.press('Escape');
  await expect(page.locator('.right-drawer')).toHaveAttribute('aria-hidden', 'true');
  expect(await search.evaluate((element) => document.activeElement === element)).toBe(true);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.previewFeatureSearchActiveIndex)).toBe(active);
  await page.getByRole('button', { name: 'Zoom in' }).click();
  await expect(page.getByRole('button', { name: 'Reset zoom' })).toContainText('110%');
  await page.getByRole('button', { name: 'Reset zoom' }).focus();
  await page.keyboard.press('Enter');
  await expect(page.getByRole('button', { name: 'Reset zoom' })).toContainText('100%');
  await page.getByRole('button', { name: 'Open active feature' }).click();
  await expect(page.locator('.feature-popup')).toBeVisible();
  await page.getByRole('button', { name: 'Close feature popup', exact: true }).click();
  await page.evaluate(() => { window.__GBDRAW_APP__.resultPanelTab = 'run-info'; });
  await expect(search).toBeHidden();
  await page.locator('.result-tabs button').first().click();
  await expect(search).toBeVisible();
  expect(await page.evaluate(() => window.__s03Nodes.search === document.querySelector('.preview-feature-search'))).toBe(true);
  await page.getByRole('combobox', { name: 'Search field' }).selectOption('label');
  await search.fill('tRNA');
  await search.press('Enter');
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.previewFeatureSearchMatches.length)).toBeGreaterThan(0);
  await page.getByRole('checkbox', { name: /Regex/ }).check();
  await search.fill('tRNA|rRNA');
  await search.press('Enter');
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.previewFeatureSearchMatches.length)).toBeGreaterThan(1);
  await page.getByRole('button', { name: 'Previous match' }).click();
  await page.getByRole('button', { name: 'Next match' }).click();
  await page.setViewportSize({ width: 390, height: 740 });
  const cdp = await page.context().newCDPSession(page);
  await cdp.send('Emulation.setPageScaleFactor', { pageScaleFactor: 2 });
  const scaled = await page.evaluate(() => ({ scale: visualViewport.scale, width: visualViewport.width,
    height: visualViewport.height }));
  expect(scaled.scale).toBeGreaterThanOrEqual(1.9);
  for (const target of [search, page.getByRole('button', { name: 'Zoom in' })]) {
    await target.evaluate((element) => element.scrollIntoView({ block: 'center' }));
    const point = await target.evaluate((element) => {
      const box = element.getBoundingClientRect();
      return { x: box.left + box.width / 2 - visualViewport.offsetLeft,
        y: box.top + box.height / 2 - visualViewport.offsetTop };
    });
    await cdp.send('Input.dispatchMouseEvent', { type: 'mousePressed', ...point, button: 'left', clickCount: 1 });
    await cdp.send('Input.dispatchMouseEvent', { type: 'mouseReleased', ...point, button: 'left', clickCount: 1 });
    expect(await target.evaluate((element) => document.activeElement === element)).toBe(true);
  }
  let reachableAt200Percent = 0;
  for (const operation of await page.locator(
    '.preview-feature-search input:not(:disabled), .preview-feature-search select:not(:disabled), '
    + '.preview-feature-search button:not(:disabled), .preview-controls button:not(:disabled)'
  ).all()) {
    const access = await operation.evaluate((element) => {
      element.focus();
      element.scrollIntoView({ block: 'center', inline: 'center' });
      const box = element.getBoundingClientRect();
      const view = visualViewport;
      const x = box.left + box.width / 2;
      const y = box.top + box.height / 2;
      return { focused: document.activeElement === element,
        visible: x >= view.offsetLeft && x <= view.offsetLeft + view.width
          && y >= view.offsetTop && y <= view.offsetTop + view.height,
        hit: element.contains(document.elementFromPoint(x, y)) };
    });
    expect(access, `${await operation.getAttribute('aria-label') || await operation.textContent()} at 200%`).toEqual({
      focused: true, visible: true, hit: true
    });
    reachableAt200Percent += 1;
  }
  observations.push({ visualViewport200Percent: scaled, reachableAt200Percent,
    zoomAfterVisualPointer: await page.evaluate(() => window.__GBDRAW_APP__.zoom) });
  await cdp.send('Emulation.setPageScaleFactor', { pageScaleFactor: 1 });
  await search.focus();
  await cdp.send('Emulation.setDeviceMetricsOverride', { width: 390, height: 480,
    screenWidth: 390, screenHeight: 740, deviceScaleFactor: 1, mobile: true });
  const keyboardViewport = await page.evaluate(() => ({ width: visualViewport.width,
    height: visualViewport.height, focused: document.activeElement?.getAttribute('aria-label') }));
  expect(keyboardViewport.height).toBeLessThanOrEqual(480);
  expect(keyboardViewport.focused).toBe('Search features');
  await search.press('Enter');
  const zoomDuringKeyboard = page.getByRole('button', { name: 'Zoom in' });
  await zoomDuringKeyboard.evaluate((element) => element.scrollIntoView({ block: 'center' }));
  const keyboardPoint = await zoomDuringKeyboard.evaluate((element) => {
    const box = element.getBoundingClientRect();
    return { x: box.left + box.width / 2 - visualViewport.offsetLeft,
      y: box.top + box.height / 2 - visualViewport.offsetTop };
  });
  await cdp.send('Input.dispatchMouseEvent', { type: 'mousePressed', ...keyboardPoint, button: 'left', clickCount: 1 });
  await cdp.send('Input.dispatchMouseEvent', { type: 'mouseReleased', ...keyboardPoint, button: 'left', clickCount: 1 });
  expect(await zoomDuringKeyboard.evaluate((element) => document.activeElement === element)).toBe(true);
  observations.push({ keyboardEquivalentViewport: keyboardViewport,
    zoomAfterKeyboardViewportPointer: await page.evaluate(() => window.__GBDRAW_APP__.zoom) });
  await cdp.send('Emulation.clearDeviceMetricsOverride');
  await writeFile(testInfo.outputPath('docked-preview-matrix.json'), JSON.stringify(observations, null, 2));
});
