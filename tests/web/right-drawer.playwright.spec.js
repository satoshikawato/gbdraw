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

test.beforeEach(async ({ page }) => {
  page.on('dialog', (dialog) => dialog.accept());
  await page.goto('/gbdraw/web/index.html', { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);
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

test('individual Feature, Label, and Legend edits update the mounted SVG', async ({
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
