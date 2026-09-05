const { test, expect } = require('@playwright/test');

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
  return {
    preview: rect(preview),
    toggle: rect(toggle),
    controls: rect(controls),
    zoomControls: rect(zoomControls),
    drawer: rect(drawer),
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
  page.on('dialog', (dialog) => dialog.accept());
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

  const openedFeaturePopup = await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    const feature = app.extractedFeatures.find((candidate) => (
      candidate?.orthogroupId && candidate?.type === 'CDS'
    ));
    if (!feature) return false;
    app.openFeatureEditorFromList(feature, null);
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

  expect(await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    const feature = app.extractedFeatures.find((candidate) => (
      candidate?.orthogroupId && candidate?.type === 'CDS'
    ));
    if (!feature) return false;
    app.openFeatureEditorFromList(feature, null);
    return Boolean(app.clickedFeature);
  })).toBe(true);
  await expect(featurePopup.locator(
    'button[title="Copy all member amino-acid FASTA"]'
  )).toHaveText(/Copy aa/);
});

test('individual Feature, Label, and Legend edits update the mounted SVG', { tag: '@pr-smoke' }, async ({
  page
}) => {
  const imported = await loadGallerySession(
    page,
    'HmmtDNA_basic_circular.gbdraw-session.json'
  );
  expect(imported.status).toBe('ok');
  await page.locator('.drawer-toggle').click();
  await expect(page.locator('.right-drawer').getByText(/^Features \(\d+\)$/)).toBeVisible();

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
    currentSvg = document.querySelector('.origin-top svg');
    labelElement = currentSvg.querySelector(
      `text[data-label-key="${CSS.escape(labelEntry.key)}"]`
    );
    const visibilityResultContent = app.results[app.selectedResultIndex]?.content || '';

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
});

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
  expect(document.version).toBe(40);
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
