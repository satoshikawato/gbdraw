const { test, expect } = require('@playwright/test');
const { writeFileSync } = require('node:fs');
const { semantics, switchMode, snapshot: readSnapshot } = require('./helpers/mode-transition.cjs');
const { semantics: svgSemantics } = require('./helpers/visual-state.cjs');

const snapshot = async (page) => {
  // Opening an imported Editor hydrates transient DOM label bindings and its
  // derived context key. Compare every drawing primitive/text with the existing
  // oracle, while keeping exact Result, payload, request and override checks.
  const { context, mounted, ...artifact } = await readSnapshot(page);
  return { ...artifact, mounted: await svgSemantics(page, mounted) };
};

const loadGallerySession = async (page, filename) => page.evaluate(async (name) => {
  const response = await fetch(`/gbdraw/web/gallery/sessions/${name}`);
  if (!response.ok) throw new Error(`Could not load ${name}: ${response.status}`);
  const bytes = await response.arrayBuffer();
  const file = new File([bytes], name, {
    type: name.endsWith('.gz') ? 'application/gzip' : 'application/json'
  });
  return window.__GBDRAW_APP__.importSession({
    target: { files: [file], value: 'selected' }
  });
}, filename);

const ISSUE_461_VIEWPORTS = [
  { id: 'issue-desktop', width: 2541, height: 1409 },
  { id: 'desktop-common', width: 1920, height: 1080 },
  { id: 'laptop', width: 1366, height: 768 },
  { id: 'narrow-supported', width: 1024, height: 768 }
];

const waitForSettledDrawer = async (page) => {
  await page.locator('.right-drawer').evaluate(async (drawer) => {
    await Promise.all(
      drawer.getAnimations().map((animation) => animation.finished.catch(() => {}))
    );
  });
  await page.locator('.drawer-toggle').evaluate(async (toggle) => {
    await Promise.all(
      toggle.getAnimations().map((animation) => animation.finished.catch(() => {}))
    );
  });
};

const readSettledPreviewSurface = async (page) => page.locator(
  '.gbdraw-preview-surface'
).evaluate(async (surface) => {
  await Promise.all(
    surface.getAnimations().map((animation) => animation.finished.catch(() => {}))
  );
  const box = surface.getBoundingClientRect();
  return {
    left: box.left,
    top: box.top,
    width: box.width,
    height: box.height,
    transform: surface.style.transform
  };
});

const readOverlayGeometry = (page) => page.evaluate(() => {
  const preview = document.querySelector('[role="region"][aria-label="Result Preview"]');
  const toggle = document.querySelector('.drawer-toggle');
  const controls = document.querySelector('.preview-controls');
  const zoomReset = document.querySelector('[aria-label="Reset zoom"]');
  const zoomControls = zoomReset?.parentElement;
  const drawer = document.querySelector('.right-drawer');
  const rect = (element) => {
    const box = element.getBoundingClientRect();
    return {
      left: box.left,
      top: box.top,
      right: box.right,
      bottom: box.bottom,
      width: box.width,
      height: box.height
    };
  };
  const receivesPointerAtCenter = (element) => {
    const box = element.getBoundingClientRect();
    const hit = document.elementFromPoint(
      box.left + box.width / 2,
      box.top + box.height / 2
    );
    return Boolean(hit && element.contains(hit));
  };
  if (!preview || !toggle || !controls || !zoomControls || !drawer) {
    throw new Error('Issue #461 overlay controls are not mounted.');
  }
  const overlays = [toggle, drawer, controls, zoomControls, ...controls.querySelectorAll('button')];
  return {
    viewport: { width: innerWidth, height: innerHeight, pageWidth: document.documentElement.scrollWidth },
    media: { max760: matchMedia('(max-width: 760px)').matches, mobile: matchMedia('(max-width: 768px)').matches,
      compactPreview: preview.parentElement.clientWidth <= 640 },
    hitTargets: overlays.map((element) => {
      const box = element.getBoundingClientRect();
      const hit = document.elementFromPoint(box.left + box.width / 2, box.top + box.height / 2);
      return { name: element.getAttribute('aria-label') || element.title || element.className,
        zIndex: getComputedStyle(element).zIndex,
        hit: hit ? { tag: hit.tagName, className: hit.getAttribute('class'), text: hit.textContent.trim().slice(0, 80) } : null };
    }),
    preview: rect(preview),
    toggle: rect(toggle),
    controls: rect(controls),
    zoomControls: rect(zoomControls),
    drawer: rect(drawer),
    canvas: rect(document.querySelector('.preview-canvas')),
    header: rect(document.querySelector('.app-header')),
    generate: getComputedStyle(document.querySelector('.generate-bar')).position === 'fixed'
      ? rect(document.querySelector('.generate-bar')) : null,
    controlButtons: Array.from(controls.querySelectorAll('button')).map((button) => ({
      name: button.getAttribute('aria-label') || button.getAttribute('title') || '',
      box: rect(button),
      receivesPointerAtCenter: receivesPointerAtCenter(button)
    })),
    toggleReceivesPointerAtCenter: receivesPointerAtCenter(toggle)
  };
});

const boxOverlap = (first, second) => {
  const x = Math.min(first.right, second.right) - Math.max(first.left, second.left);
  const y = Math.min(first.bottom, second.bottom) - Math.max(first.top, second.top);
  return { x, y, intersects: x > 0 && y > 0 };
};

const logOverlayGeometry = (viewport, drawerState, geometry) => console.log(
  `[issue-461] ${JSON.stringify({
    viewport,
    drawerState,
    editorToggle: geometry.toggle,
    zoomControls: geometry.zoomControls,
    drawer: geometry.drawer,
    toggleZoomOverlap: boxOverlap(geometry.toggle, geometry.zoomControls),
    drawerZoomOverlap: boxOverlap(geometry.drawer, geometry.zoomControls)
  })}`
);

const assertOverlayGeometry = (geometry, label, drawerOpen) => {
  const toggleZoomOverlap = boxOverlap(geometry.toggle, geometry.zoomControls);
  expect(
    toggleZoomOverlap.intersects,
    `${label}: editor toggle intersects zoom controls `
      + `(overlapX=${toggleZoomOverlap.x}, overlapY=${toggleZoomOverlap.y})`
  ).toBe(false);

  if (drawerOpen) {
    const drawerZoomOverlap = boxOverlap(geometry.drawer, geometry.zoomControls);
    expect(
      drawerZoomOverlap.intersects,
      `${label}: open Editor drawer intersects zoom controls `
        + `(overlapX=${drawerZoomOverlap.x}, overlapY=${drawerZoomOverlap.y})`
    ).toBe(false);
  }

  const insidePreview = (box) => (
    box.left >= geometry.preview.left
    && box.top >= geometry.preview.top
    && box.right <= geometry.preview.right
    && box.bottom <= geometry.preview.bottom
  );
  for (const [name, box] of [
    ['editor toggle', geometry.toggle],
    ['preview controls', geometry.controls],
    ['zoom controls', geometry.zoomControls]
  ]) {
    expect(insidePreview(box), `${label}: ${name} is outside Preview bounds`).toBe(true);
  }
  expect(geometry.toggleReceivesPointerAtCenter, `${label}: editor toggle is blocked`).toBe(true);
  for (const button of geometry.controlButtons) {
    expect(
      insidePreview(button.box),
      `${label}: ${button.name || 'unnamed control'} is outside Preview bounds`
    ).toBe(true);
    expect(
      button.receivesPointerAtCenter,
      `${label}: ${button.name || 'unnamed control'} is blocked`
    ).toBe(true);
  }
};

const exercisePreviewControls = async (page) => {
  const zoomIn = page.getByRole('button', { name: 'Zoom in', exact: true });
  const resetZoom = page.getByRole('button', { name: 'Reset zoom', exact: true });
  const zoomOut = page.getByRole('button', { name: 'Zoom out', exact: true });

  await resetZoom.click();
  await expect(resetZoom).toContainText('100%');
  await zoomIn.click();
  await expect(resetZoom).toContainText('110%');
  await resetZoom.click();
  await expect(resetZoom).toContainText('100%');
  await zoomOut.click();
  await expect(resetZoom).toContainText('90%');
  await resetZoom.click();
  await expect(resetZoom).toContainText('100%');

  const canvasPadding = page.getByRole('button', {
    name: 'Toggle canvas padding controls',
    exact: true
  });
  await canvasPadding.click();
  await expect(page.getByText('Canvas Padding', { exact: true })).toBeVisible();
  await canvasPadding.click();
  await expect(page.getByText('Canvas Padding', { exact: true })).toBeHidden();
};

test.beforeEach(async ({ page }) => {
  page.overlayDiagnostics = { pageErrors: [], consoleErrors: [], externalRequests: [] };
  page.on('pageerror', (error) => page.overlayDiagnostics.pageErrors.push(error.message));
  page.on('console', (message) => {
    if (message.type() === 'error') page.overlayDiagnostics.consoleErrors.push(message.text());
  });
  await page.route('**/*', (route) => {
    if (new URL(route.request().url()).hostname === '127.0.0.1') return route.continue();
    page.overlayDiagnostics.externalRequests.push(route.request().url());
    return route.abort();
  });
  page.on('dialog', (dialog) => dialog.accept());
  await page.addInitScript(() => {
    window.__MODE_EVENTS__ = [];
    window.__GBDRAW_TEST_HOOKS__ = {
      onSessionLifecycleEvent: (event) => window.__MODE_EVENTS__.push(event)
    };
  });
  await page.goto('/gbdraw/web/index.html', { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);
});

for (const viewport of ISSUE_461_VIEWPORTS) {
  test(`zoom controls do not overlap the Editor toggle or drawer at ${viewport.id}`, async ({
    page
  }) => {
    await page.setViewportSize({ width: viewport.width, height: viewport.height });
    const imported = await loadGallerySession(
      page,
      'majanivirus_orthogroup.gbdraw-session.json.gz'
    );
    expect(imported.status).toBe('ok');

    const toggle = page.locator('.drawer-toggle');
    const drawer = page.locator('.right-drawer');
    await expect(toggle).toHaveAttribute('aria-expanded', 'false');
    await waitForSettledDrawer(page);
    const initialSurface = await readSettledPreviewSurface(page);

    const closedGeometry = await readOverlayGeometry(page);
    logOverlayGeometry(viewport, 'closed', closedGeometry);
    assertOverlayGeometry(closedGeometry, `${viewport.id} drawer closed`, false);
    await exercisePreviewControls(page);

    await toggle.click();
    await expect(toggle).toHaveAttribute('aria-expanded', 'true');
    await expect(drawer).toHaveAttribute('aria-hidden', 'false');
    await waitForSettledDrawer(page);

    const openGeometry = await readOverlayGeometry(page);
    logOverlayGeometry(viewport, 'open', openGeometry);
    assertOverlayGeometry(openGeometry, `${viewport.id} drawer open`, true);
    await exercisePreviewControls(page);
    expect(await readSettledPreviewSurface(page)).toEqual(initialSurface);

    await toggle.click();
    await expect(toggle).toHaveAttribute('aria-expanded', 'false');
    await expect(drawer).toHaveAttribute('aria-hidden', 'true');
    await waitForSettledDrawer(page);
    expect(await readSettledPreviewSurface(page)).toEqual(initialSurface);
  });
}

const MOBILE_VIEWPORTS = [
  { id: 'MOB1', width: 390, height: 844 },
  { id: 'MOB2', width: 375, height: 667 },
  { id: 'M01-740', width: 390, height: 740 },
  { id: 'MOB3', width: 430, height: 932 },
  { id: 'MOB4', width: 844, height: 390 }
];

const centerPreview = async (page) => {
  await waitForSettledDrawer(page);
  await page.getByRole('region', { name: 'Result Preview' }).evaluate(async (preview) => {
    await document.fonts.ready;
    const canvas = preview.querySelector('.preview-canvas');
    const compact = getComputedStyle(preview.querySelector('.preview-editor-layout')).display === 'grid';
    const target = compact ? canvas : preview;
    target.scrollIntoView({ block: 'center' });
    const header = document.querySelector('.app-header').getBoundingClientRect();
    const generate = document.querySelector('.generate-bar');
    const bottom = getComputedStyle(generate).position === 'fixed'
      ? generate.getBoundingClientRect().top : innerHeight;
    const box = target.getBoundingClientRect();
    // Scroll the actual canvas/toolbar/dock into the available screen. The
    // search and export controls precede it and stay reachable by page scroll.
    window.scrollBy(0, compact ? box.top - header.bottom
      : (box.top + box.bottom - header.bottom - bottom) / 2);
  });
};

const assertMobileGeometry = (geometry, drawerOpen) => {
  const { preview, canvas, toggle, drawer, viewport, header, generate } = geometry;
  expect(viewport.pageWidth).toBeLessThanOrEqual(viewport.width);
  expect(toggle.left).toBeGreaterThanOrEqual(Math.max(0, preview.left));
  expect(toggle.right).toBeLessThanOrEqual(Math.min(viewport.width, preview.right));
  // Search and its toggle precede the canvas and may be above the screen after
  // scrolling to the dock; their hit tests run when scrolled into view.
  for (const button of geometry.controlButtons) {
    expect(boxOverlap(toggle, button.box).intersects, button.name).toBe(false);
    expect(boxOverlap(canvas, button.box).intersects, button.name).toBe(false);
    expect(boxOverlap(button.box, drawer).intersects && drawerOpen, button.name).toBe(false);
    expect(button.box.left).toBeGreaterThanOrEqual(canvas.left);
    expect(button.box.right).toBeLessThanOrEqual(canvas.right);
    expect(button.receivesPointerAtCenter, `${button.name} blocked`).toBe(true);
  }
  expect(canvas.left).toBe(preview.left + 2);
  expect(canvas.right).toBe(preview.right - 2);
  if (drawerOpen) {
    expect(drawer.left).toBe(canvas.left);
    expect(drawer.right).toBe(canvas.right);
    expect(drawer.top).toBeGreaterThanOrEqual(canvas.bottom);
    expect(boxOverlap(toggle, drawer).intersects).toBe(false);
    if (viewport.width === 390 && [844, 740].includes(viewport.height)) {
      const visibleHeight = Math.min(canvas.bottom, generate?.top ?? viewport.height)
        - Math.max(canvas.top, header.bottom);
      expect(visibleHeight, 'M01: actual visible canvas below header and above Generate').toBeGreaterThanOrEqual(200);
      expect(drawer.bottom, 'M01: complete dock above Generate').toBeLessThanOrEqual(generate.top);
    }
  }
};

const assertNoOverlayErrors = (page) => expect(page.overlayDiagnostics).toEqual({
  pageErrors: [], consoleErrors: [], externalRequests: []
});

for (const viewport of MOBILE_VIEWPORTS) {
  for (const [journey, filename] of [
    ['J43', 'HmmtDNA_ATskew.gbdraw-session.json'],
    ['J44', 'lambda_basic_linear.gbdraw-session.json']
  ]) {
    test(`mobile Editor overlays ${viewport.id} ${journey}`, async ({ page }, testInfo) => {
      test.setTimeout(180000);
      await page.setViewportSize(viewport);
      await expect(page.locator('.drawer-toggle')).toHaveCount(0);
      expect((await loadGallerySession(page, filename)).status).toBe('ok');
      const toggle = page.locator('.drawer-toggle');
      const drawer = page.locator('.right-drawer');
      const initialAccessibility = await toggle.ariaSnapshot();
      const identity = () => page.evaluate(() => ({
        index: window.__GBDRAW_APP__.selectedResultIndex,
        names: window.__GBDRAW_APP__.results.map((result) => result.name)
      }));
      const content = () => page.locator('.gbdraw-preview-surface svg').evaluate((svg) => svg.outerHTML);
      const beforeIdentity = await identity();
      const historyCounts = () => page.evaluate(() => [
        window.__GBDRAW_HISTORY__.getUndoCount(), window.__GBDRAW_HISTORY__.getRedoCount()
      ]);
      const beforeSemantics = await semantics(page, await content());
      const text = () => page.locator('.gbdraw-preview-surface svg text').allTextContents();
      const beforeText = await text();
      const capture = async (state) => {
        await centerPreview(page);
        const geometry = await readOverlayGeometry(page);
        const path = testInfo.outputPath(`${state}-geometry.json`);
        writeFileSync(path, JSON.stringify(geometry, null, 2));
        await testInfo.attach(`${state}-geometry`, { path, contentType: 'application/json' });
        await page.screenshot({ path: testInfo.outputPath(`${viewport.id}-${journey}-${state}.png`) });
        assertMobileGeometry(geometry, state === 'open');
      };
      await expect(toggle).toHaveAttribute('aria-expanded', 'false');
      await expect(toggle).toHaveAttribute('title', 'Open editor');
      await capture('closed');
      await exercisePreviewControls(page);
      await page.getByRole('button', { name: 'Toggle layout edit mode', exact: true }).click();
      await expect(page.getByRole('button', { name: 'Toggle layout edit mode', exact: true })).toHaveAttribute('aria-pressed', 'true');
      await page.getByRole('button', { name: 'Toggle layout edit mode', exact: true }).click();
      await page.getByRole('button', { name: 'Zoom in', exact: true }).click();
      await page.getByRole('button', { name: 'Reset layout', exact: true }).click();
      await expect(page.getByRole('button', { name: 'Reset zoom', exact: true })).toHaveText('100%');
      expect(await identity()).toEqual(beforeIdentity);
      expect(await semantics(page, await content())).toEqual(beforeSemantics);
      // Reset Layout above is an explicit artifact action. Start the pure
      // open/close baseline after it, so its serialization is not attributed
      // to the drawer's presentation transitions.
      const beforeArtifact = await snapshot(page);
      const beforeHistory = await historyCounts();
      await centerPreview(page);
      await toggle.click();
      await expect(toggle).toHaveAttribute('aria-expanded', 'true');
      await expect(drawer).toHaveAttribute('aria-hidden', 'false');
      await expect(toggle).toHaveAttribute('title', 'Close editor');
      await capture('open');
      await exercisePreviewControls(page);
      await centerPreview(page);
      // Exercise the real scrollable Features list, not just its overflow style.
      const scroll = await drawer.locator('.overflow-y-auto').evaluateAll((elements) => {
        const list = elements.find((element) => element.scrollHeight > element.clientHeight);
        if (!list) throw new Error('Expected a scrollable Editor feature list');
        list.scrollTop = list.scrollHeight;
        return list.scrollTop;
      });
      expect(scroll).toBeGreaterThan(0);
      await drawer.getByRole('button', { name: 'Close editor', exact: true }).click();
      await expect(drawer).toHaveAttribute('aria-hidden', 'true');
      await toggle.focus();
      await expect(toggle).toBeFocused();
      await page.keyboard.press('Enter');
      await expect(toggle).toHaveAttribute('aria-expanded', 'true');
      await page.keyboard.press('Escape');
      await expect(toggle).toHaveAttribute('aria-expanded', 'false');
      await toggle.focus();
      await page.keyboard.press('Space');
      await expect(toggle).toHaveAttribute('aria-expanded', 'true');
      await waitForSettledDrawer(page);
      await toggle.click();
      await expect(drawer).toHaveAttribute('aria-hidden', 'true');
      await centerPreview(page);
      assertMobileGeometry(await readOverlayGeometry(page), false);
      expect(await identity()).toEqual(beforeIdentity);
      expect(await semantics(page, await content())).toEqual(beforeSemantics);
      expect(await text()).toEqual(beforeText);
      expect(await snapshot(page)).toEqual(beforeArtifact);
      expect(await historyCounts()).toEqual(beforeHistory);
      expect(await toggle.ariaSnapshot()).toEqual(initialAccessibility);
      writeFileSync(testInfo.outputPath('diagnostics.json'), JSON.stringify(page.overlayDiagnostics, null, 2));
      assertNoOverlayErrors(page);
    });
  }
}

test('mobile overlays preserve mode, source replacement, and resize behavior', async ({ page }) => {
  test.setTimeout(180000);
  await page.setViewportSize({ width: 1366, height: 768 });
  await expect(page.locator('.drawer-toggle')).toHaveCount(0);
  expect((await loadGallerySession(page, 'HmmtDNA_ATskew.gbdraw-session.json')).status).toBe('ok');
  const toggle = page.locator('.drawer-toggle');
  for (const viewport of [{ width: 1366, height: 768 }, MOBILE_VIEWPORTS[0], { width: 1366, height: 768 }]) {
    await page.setViewportSize(viewport);
    await expect(toggle).toHaveCount(1);
    await expect(page.locator('.right-drawer')).toHaveCount(1);
    for (const open of [false, true]) {
      if (open) await toggle.click();
      await centerPreview(page);
      const geometry = await readOverlayGeometry(page);
      if (viewport.width === 390) assertMobileGeometry(geometry, open);
      else assertOverlayGeometry(geometry, 'resized desktop', open);
    }
    await toggle.click();
    await expect(toggle).toHaveAttribute('aria-expanded', 'false');
    expect(await toggle.getAttribute('style')).toBeNull();
  }
  await page.setViewportSize(MOBILE_VIEWPORTS[0]);
  for (const mode of ['linear', 'circular']) {
    await switchMode(page, mode);
    await expect(toggle).toHaveCount(1);
    await centerPreview(page);
    assertMobileGeometry(await readOverlayGeometry(page), false);
    await toggle.click();
    await expect(toggle).toHaveAttribute('aria-expanded', 'true');
    await waitForSettledDrawer(page);
    await toggle.click();
  }
  // Replace the actual input bytes through the upload control while retaining
  // the established last-Result availability semantics.
  const replacement = await page.evaluate(async () => {
    const session = await (await fetch('/gbdraw/web/gallery/sessions/lambda_basic_linear.gbdraw-session.json')).json();
    const resource = session.resources[session.renderRequest.records[0].source.resourceId];
    return atob(resource.data);
  });
  await page.getByLabel('GenBank/DDBJ File', { exact: true }).setInputFiles({
    name: 'lambda-replacement.gb', mimeType: 'text/plain', buffer: Buffer.from(replacement)
  });
  await expect.poll(() => page.evaluate(async () => (await import('./js/state.js')).state.circularRecordDiscovery.status)).not.toBe('loading');
  await expect(toggle).toHaveCount(1);
  await expect(page.locator('.right-drawer')).toHaveCount(1);
  await centerPreview(page);
  assertMobileGeometry(await readOverlayGeometry(page), false);
  await toggle.click();
  await expect(toggle).toHaveAttribute('aria-expanded', 'true');
  await page.locator('.right-drawer').getByRole('button', { name: 'Close editor', exact: true }).click();
  await expect(toggle).toHaveAttribute('aria-expanded', 'false');
  assertNoOverlayErrors(page);
});

test('a remembered unavailable Similarity groups tab reopens as Features', async ({
  page
}) => {
  const orthogroupImport = await loadGallerySession(
    page,
    'majanivirus_orthogroup.gbdraw-session.json.gz'
  );
  expect(orthogroupImport.status).toBe('ok');
  await expect.poll(
    () => page.evaluate(() => window.__GBDRAW_APP__.orthogroupCount)
  ).toBeGreaterThan(0);

  const toggle = page.locator('.drawer-toggle');
  const drawer = page.locator('.right-drawer');
  await toggle.click();
  await drawer.getByRole('button', { name: 'Similarity groups' }).click();
  await expect.poll(
    () => page.evaluate(() => window.__GBDRAW_APP__.rightDrawerTab)
  ).toBe('orthogroups');
  await toggle.click();
  await expect(toggle).toHaveAttribute('aria-expanded', 'false');

  const groupFreeImport = await loadGallerySession(
    page,
    'HmmtDNA_basic_circular.gbdraw-session.json'
  );
  expect(groupFreeImport.status).toBe('ok');
  await expect.poll(() => page.evaluate(() => ({
    open: window.__GBDRAW_APP__.showRightDrawer,
    tab: window.__GBDRAW_APP__.rightDrawerTab,
    groups: window.__GBDRAW_APP__.orthogroupCount
  }))).toEqual({ open: false, tab: 'features', groups: 0 });

  await toggle.click();
  await expect(toggle).toHaveAttribute('aria-expanded', 'true');
  await expect(drawer).toHaveAttribute('aria-hidden', 'false');
  await expect(drawer.getByText(/^Features \(\d+\)$/)).toBeVisible();
  await expect.poll(
    () => page.evaluate(() => window.__GBDRAW_APP__.rightDrawerTab)
  ).toBe('features');
});

test('preview similarity-group copy actions report isolated accessible outcomes', async ({
  page
}) => {
  const imported = await loadGallerySession(
    page,
    'majanivirus_orthogroup.gbdraw-session.json.gz'
  );
  expect(imported.status).toBe('ok');
  await expect.poll(
    () => page.evaluate(() => window.__GBDRAW_APP__.orthogroupCount)
  ).toBeGreaterThan(0);
  await page.clock.install();
  await page.evaluate(() => {
    window.__copiedText = '';
    Object.defineProperty(navigator, 'clipboard', {
      configurable: true,
      value: {
        writeText: async (value) => { window.__copiedText = String(value); }
      }
    });
  });

  const toggle = page.locator('.drawer-toggle');
  const drawer = page.locator('.right-drawer');
  await toggle.click();
  await drawer.getByRole('button', { name: 'Similarity groups' }).click();
  const groupNtCopy = drawer.locator(
    'button[title="Copy all member nucleotide FASTA"]'
  );
  const groupAaCopy = drawer.locator(
    'button[title="Copy all member amino-acid FASTA"]'
  );
  await expect(groupNtCopy).toHaveText(/Copy nt/);
  await expect(groupAaCopy).toHaveText(/Copy aa/);
  await groupAaCopy.click();
  await expect(groupAaCopy).toHaveText(/Copied!/);
  await expect(groupAaCopy).toHaveAccessibleName('Copied!');
  await expect(groupAaCopy).toHaveAttribute('aria-live', 'polite');
  await expect(groupAaCopy).toHaveAttribute('aria-atomic', 'true');
  await expect(groupNtCopy).toHaveText(/Copy nt/);
  const groupAminoAcidPayload = await page.evaluate(() => window.__copiedText);
  expect(groupAminoAcidPayload).toMatch(/^>/);
  expect(groupAminoAcidPayload).not.toMatch(/h_[a-z2-7]{26}/i);
  await page.clock.fastForward(1500);
  await expect(groupAaCopy).toHaveText(/Copy aa/);

  const firstMemberAaCopy = drawer.locator(
    'button[title="Copy amino-acid FASTA"]'
  ).first();
  const secondMemberAaCopy = drawer.locator(
    'button[title="Copy amino-acid FASTA"]'
  ).nth(1);
  await firstMemberAaCopy.click();
  await expect(firstMemberAaCopy).toHaveText(/Copied!/);
  await expect(secondMemberAaCopy).toHaveText(/^\s*aa\s*$/);
  const memberAminoAcidPayload = await page.evaluate(() => window.__copiedText);
  expect(memberAminoAcidPayload).toMatch(/^>/);
  expect(groupAminoAcidPayload).toContain(memberAminoAcidPayload.trim());
  await page.clock.fastForward(1500);
  await expect(firstMemberAaCopy).toHaveText(/^\s*aa\s*$/);

  await page.evaluate(() => {
    Object.defineProperty(navigator, 'clipboard', {
      configurable: true,
      value: { writeText: async () => { throw new Error('denied'); } }
    });
    document.execCommand = () => false;
  });
  await groupAaCopy.click();
  await expect(groupAaCopy).toHaveText(/Copy failed/);
  await expect(groupAaCopy).not.toHaveText(/Copied!/);
  await expect(groupNtCopy).toHaveText(/Copy nt/);
  await page.clock.fastForward(1500);
  await expect(groupAaCopy).toHaveText(/Copy aa/);

  await page.evaluate(() => {
    Object.defineProperty(navigator, 'clipboard', {
      configurable: true,
      value: {}
    });
    document.execCommand = () => true;
  });
  await groupAaCopy.click();
  await expect(groupAaCopy).toHaveText(/Copied!/);
  await page.clock.fastForward(1500);

  await page.evaluate(() => {
    Object.defineProperty(navigator, 'clipboard', {
      configurable: true,
      value: {
        writeText: async (value) => { window.__copiedText = String(value); }
      }
    });
  });
  await groupAaCopy.click();
  await expect(groupAaCopy).toHaveText(/Copied!/);
  await toggle.click();
  await toggle.click();
  await expect(groupAaCopy).toHaveText(/Copy aa/);

  const openedFeaturePopup = await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    const feature = app.extractedFeatures.find((candidate) => (
      candidate?.orthogroupId && candidate?.type === 'CDS'
    ));
    if (!feature) return false;
    await app.openFeatureEditorFromList(feature, null);
    return Boolean(app.clickedFeature);
  });
  expect(openedFeaturePopup).toBe(true);
  const featurePopup = page.locator('.feature-popup');
  const popupGroupAaCopy = featurePopup.locator(
    'button[title="Copy all member amino-acid FASTA"]'
  );
  const popupMemberAaCopy = featurePopup.locator(
    'button[title="Copy amino-acid FASTA"]'
  ).first();
  await popupGroupAaCopy.click();
  await expect(popupGroupAaCopy).toHaveText(/Copied!/);
  await expect(popupMemberAaCopy).toHaveText(/^\s*aa\s*$/);
  await featurePopup.getByRole('button', { name: 'Close feature popup' }).click();

  expect(await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    const feature = app.extractedFeatures.find((candidate) => (
      candidate?.orthogroupId && candidate?.type === 'CDS'
    ));
    if (!feature) return false;
    await app.openFeatureEditorFromList(feature, null);
    return Boolean(app.clickedFeature);
  })).toBe(true);
  await expect(featurePopup.locator(
    'button[title="Copy all member amino-acid FASTA"]'
  )).toHaveText(/Copy aa/);
});

for (const liveCase of [
  { id: 'desktop', session: 'HmmtDNA_basic_circular.gbdraw-session.json' },
  { id: 'compact Circular', session: 'HmmtDNA_basic_circular.gbdraw-session.json', viewport: { width: 390, height: 740 } },
  { id: 'compact Linear', session: 'lambda_basic_linear.gbdraw-session.json', viewport: { width: 390, height: 844 } }
]) {
test(`individual Feature, Label, and Legend edits update the mounted SVG: ${liveCase.id}`, { tag: liveCase.viewport ? [] : '@pr-smoke' }, async ({
  page
}) => {
  test.setTimeout(180000);
  if (liveCase.viewport) await page.setViewportSize(liveCase.viewport);
  const imported = await loadGallerySession(page, liveCase.session);
  expect(imported.status).toBe('ok');
  await page.locator('.drawer-toggle').click();
  await expect(page.locator('.right-drawer').getByText(/^Features \(\d+\)$/)).toBeVisible();
  expect(await page.evaluate(() => window.__GBDRAW_APP__.runAnalysis())).toEqual({ status: 'ok' });

  await page.evaluate(() => window.__GBDRAW_APP__.openRightDrawerTab('features'));
  if (liveCase.viewport) {
    await centerPreview(page);
    assertMobileGeometry(await readOverlayGeometry(page), true);
  }
  const result = await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    const svg = document.querySelector('.origin-top svg');
    const { getFeatureFillElements } = await import(
      '/gbdraw/web/js/app/feature-editor/svg-actions.js'
    );
    const feature = app.extractedFeatures.find((candidate) => {
      const featureId = String(candidate?.svg_id || '').trim();
      return featureId && getFeatureFillElements(svg, featureId).length > 0;
    });
    if (!feature) throw new Error('No editable rendered feature was found.');

    const featureColor = '#1234ab';
    await app.setFeatureColorValue(feature, featureColor);
    const featureResultContent = app.results[app.selectedResultIndex]?.content || '';

    const labeledFeature = app.extractedFeatures.find((candidate) => (
      app.getEditableLabelByFeatureId(candidate?.svg_id)
    ));
    if (!labeledFeature) throw new Error('No editable feature label was found.');
    const labelEntry = app.getEditableLabelByFeatureId(labeledFeature.svg_id);
    const labelText = 'Drawer live label';
    await app.requestLabelTextChangeByFeatureId(labeledFeature.svg_id, labelText);
    await app.handleLabelTextScopeChoice('single');
    let currentSvg = document.querySelector('.origin-top svg');
    let labelElement = currentSvg.querySelector(
      `text[data-label-key="${CSS.escape(labelEntry.key)}"]`
    );
    const textResultContent = app.results[app.selectedResultIndex]?.content || '';
    app.openFeatureEditorFromList(labeledFeature, null);
    app.clickedFeature.labelVisibility = 'off';
    await app.updateClickedFeatureLabelText();
    for (let attempt = 0; attempt < 1200; attempt += 1) {
      currentSvg = document.querySelector('.origin-top svg');
      const currentText = Array.from(currentSvg.querySelectorAll('text[data-label-feature-id]'))
        .find((element) => element.getAttribute('data-label-feature-id') === labeledFeature.svg_id);
      const currentParts = Array.from(currentSvg.querySelectorAll('[data-label-feature-id]'))
        .filter((element) => element.getAttribute('data-label-feature-id') === labeledFeature.svg_id);
      if (
        currentText?.getAttribute('data-gbdraw-label-binding-schema') === '1'
        && currentParts.length > 0
        && currentParts.every((element) => element.getAttribute('display') === 'none')
      ) break;
      await new Promise((resolve) => window.setTimeout(resolve, 50));
    }
    currentSvg = document.querySelector('.origin-top svg');
    labelElement = Array.from(currentSvg.querySelectorAll('text[data-label-feature-id]'))
      .find((element) => element.getAttribute('data-label-feature-id') === labeledFeature.svg_id);
    const labelParts = Array.from(currentSvg.querySelectorAll('[data-label-feature-id]'))
      .filter((element) => element.getAttribute('data-label-feature-id') === labeledFeature.svg_id);
    const visibilityResultContent = app.results[app.selectedResultIndex]?.content || '';
    const resultSvg = new DOMParser().parseFromString(
      visibilityResultContent,
      'image/svg+xml'
    );
    const resultLabelParts = Array.from(resultSvg.querySelectorAll('[data-label-feature-id]'))
      .filter((element) => element.getAttribute('data-label-feature-id') === labeledFeature.svg_id);

    const legendIndex = app.legendEntries.findIndex((entry) => {
      if (!entry?.caption) return false;
      const group = currentSvg.querySelector(
        `g[data-legend-key="${CSS.escape(entry.caption)}"]`
      );
      return Boolean(group?.querySelector('path[fill]:not([fill="none"])'));
    });
    if (legendIndex < 0) throw new Error('No editable legend entry was found.');
    const legendCaption = app.legendEntries[legendIndex].caption;
    const legendColor = '#ab3412';
    app.updateLegendEntryColor(legendIndex, legendColor);
    currentSvg = document.querySelector('.origin-top svg');
    const legendFill = currentSvg.querySelector(
      `g[data-legend-key="${CSS.escape(legendCaption)}"] path[fill]:not([fill="none"])`
    )?.getAttribute('fill');

    return {
      featureFills: getFeatureFillElements(currentSvg, feature.svg_id)
        .map((element) => element.getAttribute('fill')),
      featureOverrideColors: Object.values(app.featureColorOverrides)
        .map((override) => override?.color),
      labelText: labelElement?.textContent || '',
      labelDisplay: labelElement?.getAttribute('display') || '',
      labelPartCount: labelParts.length,
      labelLeaderCount: labelParts.filter((element) => element.localName === 'line').length,
      labelPartsHidden: labelParts.every((element) => element.getAttribute('display') === 'none'),
      resultLabelPartsHidden: resultLabelParts.length === labelParts.length
        && resultLabelParts.every((element) => element.getAttribute('display') === 'none'),
      labelVisibilityOverride: app.labelVisibilityOverrides[labeledFeature.svg_id],
      labelOverride: app.labelTextFeatureOverrides[labeledFeature.svg_id],
      legendFill,
      legendEntryColor: app.legendEntries[legendIndex]?.color,
      featureResultContent,
      textResultContent,
      visibilityResultContent,
      legendResultContent: app.results[app.selectedResultIndex]?.content || ''
    };
  });

  expect(result.featureFills.length).toBeGreaterThan(0);
  expect(result.featureFills.every((fill) => fill === '#1234ab')).toBe(true);
  expect(result.featureOverrideColors).toContain('#1234ab');
  expect(result.featureResultContent).toContain('#1234ab');
  expect(result.labelText).toBe('Drawer live label');
  expect(result.labelDisplay).toBe('none');
  expect(result.labelPartCount).toBeGreaterThan(1);
  expect(result.labelLeaderCount).toBeGreaterThan(0);
  expect(result.labelPartsHidden).toBe(true);
  expect(result.resultLabelPartsHidden).toBe(true);
  expect(result.labelVisibilityOverride).toBe('off');
  expect(result.labelOverride).toBe('Drawer live label');
  expect(result.legendFill).toBe('#ab3412');
  expect(result.legendEntryColor).toBe('#ab3412');
  expect(result.textResultContent).toContain('Drawer live label');
  expect(result.visibilityResultContent).toContain(
    'data-gbdraw-label-visibility-preview="off"'
  );
  expect(result.visibilityResultContent).toContain('display="none"');
  expect(result.legendResultContent).toContain('#ab3412');
  if (liveCase.viewport) {
    // Isolate retry recovery from the preceding Legend override scenario.
    // Reuse a fresh Gallery request as in the existing S04 failure contract.
    expect((await loadGallerySession(page, liveCase.session)).status).toBe('ok');
    expect(await page.evaluate(() => window.__GBDRAW_APP__.runAnalysis())).toEqual({ status: 'ok' });
    await page.evaluate(() => window.__GBDRAW_APP__.openRightDrawerTab('features'));
    await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.labelReflowProcessing)).toBe(false);
    const target = await page.evaluate(async () => {
      const app = window.__GBDRAW_APP__;
      const { state } = await import('./js/state.js');
      const feature = app.extractedFeatures.find((candidate) => app.getEditableLabelByFeatureId(candidate.svg_id));
      state.autoLabelReflowEnabled.value = true;
      window.__GBDRAW_TEST_HOOKS__.beforeDiagramGenerationResponse = () => {
        throw new Error('S05 forced live rerender failure');
      };
      app.openFeatureEditorFromList(feature, null);
      return feature.svg_id;
    });
    const input = page.locator('.feature-popup input[placeholder="Edit label text"]');
    await input.fill('S05 direct label retained');
    await page.getByRole('button', { name: 'Apply Label', exact: true }).click();
    await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.labelReflowLastError), { timeout: 180000 })
      .toContain('S05 forced live rerender failure');
    await expect(page.locator('[data-live-application-feedback]')).toContainText('Live edit failed');
    expect(await page.evaluate((id) => window.__GBDRAW_APP__.labelTextFeatureOverrides[id], target)).toBe('S05 direct label retained');
    const retained = await snapshot(page);
    expect(retained.result).toContain('S05 direct label retained');
    await page.keyboard.press('Escape');
    await expect(page.locator('.right-drawer')).toHaveAttribute('aria-hidden', 'true');
    expect(await snapshot(page)).toEqual(retained);
    await page.locator('.drawer-toggle').click();
    await page.evaluate((id) => {
      delete window.__GBDRAW_TEST_HOOKS__.beforeDiagramGenerationResponse;
      const app = window.__GBDRAW_APP__;
      app.openFeatureEditorFromList(app.extractedFeatures.find((feature) => feature.svg_id === id), null);
    }, target);
    await input.fill('S05 retry label succeeds');
    await page.getByRole('button', { name: 'Apply Label', exact: true }).click();
    await expect.poll(() => page.evaluate(() => ({
      busy: window.__GBDRAW_APP__.labelReflowProcessing,
      error: window.__GBDRAW_APP__.labelReflowLastError
    })), { timeout: 180000 }).toEqual({ busy: false, error: null });
    await expect(page.locator('[data-live-application-feedback]')).toHaveCount(0);
    expect(await page.evaluate((id) => window.__GBDRAW_APP__.labelTextFeatureOverrides[id], target)).toBe('S05 retry label succeeds');
    expect((await snapshot(page)).result).toContain('S05 retry label succeeds');
    await expect(page.locator('.right-drawer')).toHaveAttribute('aria-hidden', 'false');
  }
});

}

test('adjacent Collinear mixed groups remain selectable after current-session save/load', async ({ page }, testInfo) => {
  test.setTimeout(180000);
  const { execFileSync } = require('node:child_process');
  const { mkdirSync, readFileSync } = require('node:fs');
  const { join } = require('node:path');
  const { gunzipSync } = require('node:zlib');
  const directory = testInfo.outputPath('mixed');
  mkdirSync(directory, { recursive: true });
  execFileSync('python', ['-c',
    'import sys; from pathlib import Path; from tests.test_api_session import _record_local_collinear_session; _record_local_collinear_session(Path(sys.argv[1]))',
    directory
  ], { cwd: process.cwd(), stdio: 'pipe' });
  const typedSession = readFileSync(join(directory, 'mixed.gbdraw-session.json'), 'utf8');
  await page.setViewportSize({ width: 1600, height: 1000 });
  await page.goto('/gbdraw/web/index.html', { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);
  // The Python writer has no active Web draft; derive it from this request.
  const session = Buffer.from(await page.evaluate(async (source) => {
    const { projectCanonicalSessionRequest } = await import('./js/services/session-request.js');
    const document = JSON.parse(source);
    document.config = projectCanonicalSessionRequest({
      renderRequest: document.renderRequest, resources: document.resources,
      webFiles: document.webFiles
    }).config;
    document.config.linearComparisonPlan = { mode: 'none', defaultSource: 'losat', edges: [] };
    return JSON.stringify(document);
  }, typedSession));
  const importSession = async (bytes, name) => page.evaluate(async ({ bytes, name }) => {
    const file = new File([new Uint8Array(bytes)], name);
    const result = await window.__GBDRAW_APP__.importSession({ target: { files: [file], value: 'selected' } });
    if (result.status !== 'ok') throw new Error(result.error?.stack || result.status);
    return { status: result.status };
  }, { bytes: [...bytes], name });
  expect(await importSession(session, 'mixed.gbdraw-session.json')).toMatchObject({ status: 'ok' });

  const verifyGroups = async () => {
    const state = await page.evaluate(async () => {
      const app = window.__GBDRAW_APP__;
      const { state } = await import('./js/state.js');
      return {
        semantic: app.orthogroups.map((group) => group.id),
        presentation: state.collinearGroups.value.map((group) => group.id),
        features: app.extractedFeatures.map((feature) => ({
          protein: feature.protein_id, group: feature.orthogroupId || null,
          role: feature.orthogroupMember?.role || null
        }))
      };
    });
    expect(state.semantic).toEqual(['og_1', 'og_2']);
    expect(state.presentation).toEqual(['og_1']);
    expect(Object.fromEntries(state.features.map((feature) => [feature.protein, feature.group]))).toEqual({
      a0: 'og_1', a1: 'og_1', a3: 'og_2', a4: 'og_2', a2: null, b0: 'og_1', b1: null
    });
    expect(state.features.find((feature) => feature.protein === 'a1').role).toBe('inparalog');
    for (const id of ['a3', 'a4']) {
      expect(state.features.find((feature) => feature.protein === id).role).toBe('local_paralog');
    }
    if (!await page.evaluate(() => window.__GBDRAW_APP__.showRightDrawer)) {
      await page.locator('.drawer-toggle').click();
    }
    const drawer = page.locator('.right-drawer');
    const groupsTab = drawer.getByRole('button', { name: 'Similarity groups' });
    await expect(groupsTab).toBeVisible();
    await groupsTab.click();
    const paths = page.locator('.gbdraw-preview-surface path[data-match-kind="collinear"]');
    await expect(paths).toHaveCount(1);
    const geometry = await paths.evaluateAll((elements) => elements.map((element) => element.getAttribute('d')));
    for (const [id, expectedMembers] of [['og_1', ['a0', 'a1', 'b0']], ['og_2', ['a3', 'a4']]]) {
      const entry = drawer.locator('button').filter({ has: page.locator('.font-mono', { hasText: new RegExp(`^${id}$`) }) });
      await expect(entry).toHaveCount(1);
      await entry.click();
      expect(await page.evaluate(() => window.__GBDRAW_APP__.selectedOrthogroupId)).toBe(id);
      expect(await page.evaluate(() => window.__GBDRAW_APP__.selectedOrthogroup.members.map((member) => member.proteinId).sort())).toEqual(expectedMembers);
      const highlight = drawer.getByRole('button', { name: /Highlight$/ });
      await expect(highlight).toBeVisible();
      await highlight.click();
      const highlighted = await page.evaluate(() => {
        const app = window.__GBDRAW_APP__;
        return app.extractedFeatures.filter((feature) => [...document.querySelectorAll('[data-gbdraw-feature-id]')].some((element) =>
          element.getAttribute('data-gbdraw-feature-id') === feature.svg_id && element.getAttribute('stroke') === '#2563eb'
        )).map((feature) => feature.protein_id).sort();
      });
      expect(highlighted).toEqual(expectedMembers);
      await expect(paths).toHaveCount(1);
      await expect(paths).toHaveAttribute('data-orthogroup-id', 'og_1');
      expect(await paths.evaluateAll((elements) => elements.map((element) => element.getAttribute('d')))).toEqual(geometry);
      await expect(page.locator('.gbdraw-preview-surface path[data-orthogroup-id="og_2"]')).toHaveCount(0);
    }
    await page.screenshot({ path: testInfo.outputPath('mixed-groups.png'), fullPage: true });
    await paths.dispatchEvent('click', { clientX: 500, clientY: 350 });
    expect(await page.evaluate(() => {
      const match = window.__GBDRAW_APP__.clickedPairwiseMatch;
      return match?.blockOrthogroups.map((group) => group.memberRows.length);
    })).toEqual([3]);
    await page.keyboard.press('Escape');
    return geometry;
  };
  const before = await verifyGroups();
  const downloadPromise = page.waitForEvent('download');
  await page.evaluate(async () => {
    window.__GBDRAW_APP__.sessionTitle = 'issue460-round-trip';
    await window.__GBDRAW_APP__.saveSessionWithTitle();
  });
  const download = await downloadPromise;
  const saved = readFileSync(await download.path());
  const document = JSON.parse(gunzipSync(saved).toString('utf8'));
  expect(document.version).toBe(44);
  expect(document.editorState.featureCatalog.items[0].orthogroups.map((group) => group.id)).toEqual(['og_1', 'og_2']);
  expect(await importSession(saved, download.suggestedFilename())).toMatchObject({ status: 'ok' });
  expect(await verifyGroups()).toEqual(before);
  const rollback = await page.evaluate(async (bytes) => {
    const { state } = await import('./js/state.js');
    const { importSession } = await import('./js/services/config.js');
    const previous = state.collinearGroups.value;
    const result = await importSession({ target: {
      files: [new File([new Uint8Array(bytes)], 'failed.gbdraw-session.json.gz')], value: 'selected'
    } }, { beforePreviewMount: () => { throw new Error('Issue 460 rollback probe'); } });
    return {
      status: result.status, message: result.error?.message,
      samePresentation: state.collinearGroups.value === previous
    };
  }, [...saved]);
  expect(rollback).toEqual({ status: 'error', message: 'Issue 460 rollback probe', samePresentation: true });
  expect(await verifyGroups()).toEqual(before);
});


// Scroll recovery checks keep the same center hit oracle as M01. A short
// screen need not show every operation at once, but must expose each in turn.
const exposeControl = async (page, control) => {
  await control.evaluate((element) => {
    element.scrollIntoView({ block: 'center', inline: 'nearest' });
    const header = document.querySelector('.app-header').getBoundingClientRect();
    const generate = document.querySelector('.generate-bar');
    const bottom = getComputedStyle(generate).position === 'fixed'
      ? generate.getBoundingClientRect().top : innerHeight;
    const box = element.getBoundingClientRect();
    window.scrollBy(0, (box.top + box.bottom - header.bottom - bottom) / 2);
  });
  expect(await control.evaluate((element) => {
    const box = element.getBoundingClientRect();
    return element.contains(document.elementFromPoint(
      box.left + box.width / 2, box.top + box.height / 2
    ));
  }), 'reachable control must receive its actual center pointer hit').toBe(true);
};

for (const viewport of [
  { width: 390, height: 500 },
  { width: 320, height: 740 },
  { width: 844, height: 390 },
  // 200% browser zoom has half the CSS width/height of a 390 x 844 screen.
  // This tests that layout size; it does not operate a physical browser's UI.
  { width: 195, height: 422 }
]) {
  test(`compact Editor scroll recovery at ${viewport.width}x${viewport.height}`, async ({ page }, info) => {
    await page.setViewportSize(viewport);
    expect((await loadGallerySession(page, 'majanivirus_orthogroup.gbdraw-session.json.gz')).status).toBe('ok');
    const initialPageWidth = await page.evaluate(() => document.documentElement.scrollWidth);
    const toggle = page.locator('.drawer-toggle');
    await exposeControl(page, toggle);
    await toggle.focus();
    await page.keyboard.press('Enter');
    const drawer = page.locator('.right-drawer');
    await expect(drawer).toHaveAttribute('aria-hidden', 'false');
    for (const name of ['Legend', 'Features', 'Similarity groups']) {
      const tab = drawer.getByRole('button', { name: new RegExp(`${name}$`) });
      await exposeControl(page, tab);
      await tab.focus();
      await page.keyboard.press('Space');
      await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.rightDrawerTab))
        .toBe({ Legend: 'legend', Features: 'features', 'Similarity groups': 'orthogroups' }[name]);
    }
    const close = drawer.getByRole('button', { name: 'Close editor', exact: true });
    await exposeControl(page, close);
    await close.focus();
    await page.keyboard.press('Tab');
    await expect(drawer.getByRole('button', { name: /Legend$/ })).toBeFocused();
    await page.keyboard.press('Tab');
    await expect(drawer.getByRole('button', { name: /Features$/ })).toBeFocused();
    await page.keyboard.press('Tab');
    await expect(drawer.getByRole('button', { name: /Similarity groups$/ })).toBeFocused();
    for (const name of ['Reset layout', 'Toggle layout edit mode', 'Toggle canvas padding controls', 'Zoom in', 'Reset zoom', 'Zoom out']) {
      const control = page.getByRole('button', { name, exact: true });
      await exposeControl(page, control);
      await control.focus();
      await expect(control).toBeFocused();
    }
    await exposeControl(page, page.getByRole('button', { name: 'Zoom in', exact: true }));
    await page.getByRole('button', { name: 'Zoom in', exact: true }).click();
    await expect(page.getByRole('button', { name: 'Reset zoom', exact: true })).toHaveText('110%');
    await exposeControl(page, drawer.getByRole('button', { name: /Features$/ }));
    await drawer.getByRole('button', { name: /Features$/ }).click();
    const input = drawer.getByPlaceholder('Search by feature or annotation...');
    await exposeControl(page, input);
    await input.focus();
    // Keep native input focus while reducing the viewport, as with a software
    // keyboard. The physical OS keyboard is outside this automated evidence.
    await page.setViewportSize({ width: viewport.width, height: Math.min(500, viewport.height) });
    await expect(input).toBeFocused();
    await exposeControl(page, input);
    await input.fill('protein');
    await exposeControl(page, close);
    await page.keyboard.press('Escape');
    await expect(drawer).toHaveAttribute('aria-hidden', 'true');
    expect(await page.evaluate(() => window.__GBDRAW_APP__.rightDrawerTab)).toBe('features');
    await exposeControl(page, toggle);
    await toggle.click();
    await exposeControl(page, close);
    await close.click();
    await expect(drawer).toHaveAttribute('aria-hidden', 'true');
    await exposeControl(page, page.getByRole('button', { name: 'Generate Diagram', exact: true }));
    const geometry = await readOverlayGeometry(page);
    // Below 320 CSS px, pre-existing settings fields can overflow the page.
    // The Editor must add none; all its controls pass exact center hit checks.
    expect(geometry.viewport.pageWidth).toBeLessThanOrEqual(Math.max(viewport.width, initialPageWidth));
    writeFileSync(info.outputPath('recovery-geometry.json'), JSON.stringify(geometry, null, 2));
    assertNoOverlayErrors(page);
  });
}

for (const [mode, session] of [
  ['circular', 'HmmtDNA_basic_circular.gbdraw-session.json'],
  ['linear', 'lambda_basic_linear.gbdraw-session.json']
]) {
  test(`compact ${mode} pan zoom and responsive transitions preserve artifact and History`, async ({ page }, info) => {
    await page.setViewportSize({ width: 390, height: 740 });
    expect((await loadGallerySession(page, session)).status).toBe('ok');
    const before = await snapshot(page);
    const beforeHistory = await page.evaluate(() => [window.__GBDRAW_HISTORY__.getUndoCount(), window.__GBDRAW_HISTORY__.getRedoCount()]);
    await page.evaluate(() => {
      window.s05Root = document.querySelector('.gbdraw-preview-surface svg');
      window.s05WorkerCalls = 0;
      const post = Worker.prototype.postMessage;
      Worker.prototype.postMessage = function (...args) {
        window.s05WorkerCalls += 1;
        return Reflect.apply(post, this, args);
      };
    });
    await page.locator('.drawer-toggle').click();
    await centerPreview(page);
    const canvas = page.locator('.preview-canvas');
    const box = await canvas.boundingBox();
    const pan = await page.evaluate(() => ({ ...window.__GBDRAW_APP__.canvasPan }));
    await page.mouse.move(box.x + 4, box.y + 4);
    await page.mouse.down();
    await page.mouse.move(box.x + 34, box.y + 16, { steps: 5 });
    await page.mouse.up();
    await expect.poll(() => page.evaluate(() => ({ ...window.__GBDRAW_APP__.canvasPan }))).toEqual({ x: pan.x + 30, y: pan.y + 12 });
    await page.mouse.move(box.x + 4, box.y + 4);
    await page.mouse.wheel(0, -100);
    await expect(page.getByRole('button', { name: 'Reset zoom', exact: true })).toHaveText('110%');
    for (const viewport of [{ width: 1600, height: 1000 }, { width: 390, height: 844 }, { width: 390, height: 740 }]) {
      await page.setViewportSize(viewport);
      await waitForSettledDrawer(page);
      expect(await page.evaluate(() => window.__GBDRAW_APP__.showRightDrawer)).toBe(true);
      expect(await snapshot(page)).toEqual(before);
    }
    await page.evaluate(() => window.__GBDRAW_APP__.closeRightDrawer());
    await page.evaluate(() => window.__GBDRAW_APP__.openRightDrawerTab('legend'));
    await page.keyboard.press('Escape');
    expect(await page.evaluate(() => window.__GBDRAW_APP__.rightDrawerTab)).toBe('legend');
    expect(await snapshot(page)).toEqual(before);
    expect(await page.evaluate(() => [window.__GBDRAW_HISTORY__.getUndoCount(), window.__GBDRAW_HISTORY__.getRedoCount()])).toEqual(beforeHistory);
    expect(await page.evaluate(() => ({ same: document.querySelector('.gbdraw-preview-surface svg') === window.s05Root, workers: window.s05WorkerCalls })))
      .toEqual({ same: true, workers: 0 });
    await page.evaluate(() => window.__GBDRAW_APP__.openRightDrawerTab('features'));
    await centerPreview(page);
    assertMobileGeometry(await readOverlayGeometry(page), true);
    writeFileSync(info.outputPath('final-geometry.json'), JSON.stringify(await readOverlayGeometry(page), null, 2));
    // Reproducible evidence camera only: the real pointer pan/wheel above proves
    // operation. Frame the same Gallery SVG for readable final UI inspection.
    for (const viewport of [{ width: 1600, height: 1000 }, { width: 390, height: 844 }, { width: 390, height: 740 }]) {
      await page.setViewportSize(viewport);
      await centerPreview(page);
      await page.evaluate(async () => {
        const app = window.__GBDRAW_APP__;
        const canvas = document.querySelector('.preview-canvas').getBoundingClientRect();
        const svg = document.querySelector('.gbdraw-preview-surface svg');
        const box = svg.getBoundingClientRect();
        const drawer = document.querySelector('.right-drawer');
        const right = getComputedStyle(drawer).position === 'absolute'
          ? drawer.getBoundingClientRect().left : canvas.right;
        app.zoom = Math.min((right - canvas.left - 16) / (box.width / app.zoom),
          (canvas.height - 16) / (box.height / app.zoom));
        await window.Vue.nextTick();
      });
      await readSettledPreviewSurface(page);
      await page.evaluate(async () => {
        const app = window.__GBDRAW_APP__;
        const canvas = document.querySelector('.preview-canvas').getBoundingClientRect();
        const svg = document.querySelector('.gbdraw-preview-surface svg').getBoundingClientRect();
        const drawer = document.querySelector('.right-drawer');
        const right = getComputedStyle(drawer).position === 'absolute'
          ? drawer.getBoundingClientRect().left : canvas.right;
        app.canvasPan.x += (canvas.left + right - svg.left - svg.right) / 2;
        app.canvasPan.y += (canvas.top + canvas.bottom - svg.top - svg.bottom) / 2;
        const pane = drawer.querySelector(':scope > .flex-1:not([style*="display: none"])');
        const list = pane.querySelector('.overflow-y-auto');
        pane.scrollTop = list.offsetTop - pane.offsetTop;
        await window.Vue.nextTick();
      });
      await readSettledPreviewSurface(page);
      const geometry = await readOverlayGeometry(page);
      writeFileSync(info.outputPath(`${viewport.width}x${viewport.height}-geometry.json`), JSON.stringify(geometry, null, 2));
      await page.screenshot({ path: info.outputPath(`${mode}-${viewport.width}x${viewport.height}-editor.png`) });
      expect(await snapshot(page)).toEqual(before);
      expect(await page.evaluate(() => window.s05WorkerCalls)).toBe(0);
    }
    assertNoOverlayErrors(page);
  });
}


test('compact Editor reconciles availability and restores its tab after failed Session replacement', async ({ page }) => {
  await page.setViewportSize({ width: 390, height: 740 });
  expect((await loadGallerySession(page, 'majanivirus_orthogroup.gbdraw-session.json.gz')).status).toBe('ok');
  await page.evaluate(() => window.__GBDRAW_APP__.openRightDrawerTab('orthogroups'));
  const before = await snapshot(page);
  const failure = await page.evaluate(async () => {
    const name = 'HmmtDNA_basic_circular.gbdraw-session.json';
    const file = new File([await (await fetch(`/gbdraw/web/gallery/sessions/${name}`)).arrayBuffer()], name);
    const { importSession } = await import('./js/services/config.js');
    const result = await importSession({ target: { files: [file], value: 'selected' } }, {
      beforePreviewMount: () => { throw new Error('S05 Session rollback probe'); }
    });
    return { status: result.status, message: result.error?.message };
  });
  expect(failure).toEqual({ status: 'error', message: 'S05 Session rollback probe' });
  expect(await page.evaluate(() => ({
    open: window.__GBDRAW_APP__.showRightDrawer, tab: window.__GBDRAW_APP__.rightDrawerTab
  }))).toEqual({ open: true, tab: 'orthogroups' });
  const restored = await snapshot(page);
  // Session rollback serializes the restored SVG and normalizes xmlns/attribute
  // order; verify complete drawing semantics plus exact canonical artifacts.
  expect(await svgSemantics(page, restored.result)).toEqual(await svgSemantics(page, before.result));
  for (const key of ['result', 'payload']) {
    const normalized = await page.evaluate(async ({ actual, expected }) => {
      const { sanitizeSvgContent } = await import('./js/services/svg-sanitization.js');
      const { serializeCleanSvg } = await import('./js/services/svg-serialization.js');
      return [actual, expected].map(content => {
        const svg = new DOMParser().parseFromString(sanitizeSvgContent(content, window.DOMPurify), 'image/svg+xml').documentElement;
        // XMLSerializer may inject these namespaces at the beginning on restore.
        // Reinsert only the namespace declarations in one deterministic order.
        svg.removeAttribute('xmlns'); svg.removeAttribute('xmlns:xlink');
        return serializeCleanSvg(svg);
      });
    }, { actual: restored[key], expected: before[key] });
    expect(normalized[0], `sanitized rollback ${key}`).toBe(normalized[1]);
  }
  for (const key of ['request', 'labels', 'visibility', 'colors']) {
    expect(restored[key], `rollback ${key}`).toEqual(before[key]);
  }
  const reconcile = await page.evaluate(async () => {
    const { state } = await import('./js/state.js');
    const { captureRightDrawerState, restoreRightDrawerState, resetRightDrawerState } = await import('./js/app/right-drawer.js');
    const saved = captureRightDrawerState(state);
    state.orthogroups.value = [];
    state.similarityAlignmentPlan.value = null;
    // No nextTick: the availability watcher and restore use the same sync owner.
    const lost = captureRightDrawerState(state);
    restoreRightDrawerState(state, saved);
    const fallback = captureRightDrawerState(state);
    resetRightDrawerState(state);
    return { lost, fallback, reset: captureRightDrawerState(state) };
  });
  expect(reconcile).toEqual({
    lost: { showRightDrawer: true, rightDrawerTab: 'features' },
    fallback: { showRightDrawer: true, rightDrawerTab: 'features' },
    reset: { showRightDrawer: false, rightDrawerTab: 'features' }
  });
  expect((await loadGallerySession(page, 'lambda_basic_linear.gbdraw-session.json')).status).toBe('ok');
  expect(await page.evaluate(() => ({
    open: window.__GBDRAW_APP__.showRightDrawer, tab: window.__GBDRAW_APP__.rightDrawerTab
  }))).toEqual({ open: false, tab: 'features' });
  await page.locator('.drawer-toggle').click();
  await centerPreview(page);
  assertMobileGeometry(await readOverlayGeometry(page), true);
  expect(page.overlayDiagnostics).toEqual({
    pageErrors: [], externalRequests: [],
    consoleErrors: [expect.stringMatching(/^Error: S05 Session rollback probe\n/)]
  });
});
