const { test, expect } = require('@playwright/test');
const { readFileSync } = require('node:fs');
const { resolve } = require('node:path');
const { gunzipSync } = require('node:zlib');

const SESSION = 'HmmtDNA_basic_circular.gbdraw-session.json';

const loadSession = async (page) => {
  page.on('dialog', (dialog) => dialog.accept());
  await page.goto('/gbdraw/web/index.html', { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);
  const status = await page.evaluate(async (name) => {
    const response = await fetch(`/gbdraw/web/gallery/sessions/${name}`);
    const file = new File([await response.arrayBuffer()], name, {
      type: 'application/json'
    });
    const imported = await window.__GBDRAW_APP__.importSession({
      target: { files: [file], value: 'selected' }
    });
    return imported?.status;
  }, SESSION);
  expect(status).toBe('ok');
  await page.locator('.drawer-toggle').click();
};

const documentSnapshot = async (page) => page.evaluate(() => {
  const app = window.__GBDRAW_APP__;
  const hash = (value) => {
    let result = 2166136261;
    for (let index = 0; index < value.length; index += 1) {
      result ^= value.charCodeAt(index);
      result = Math.imul(result, 16777619);
    }
    return result >>> 0;
  };
  return {
    rules: JSON.stringify(app.manualSpecificRules),
    overrides: JSON.stringify(app.featureColorOverrides),
    strokeOverrides: JSON.stringify(app.featureStrokeOverrides),
    strokeModel: JSON.stringify(
      app.clickedFeature ? app.getFeatureStrokeViewModel(app.clickedFeature.feat) : null
    ),
    labelOverrides: JSON.stringify({
      feature: app.labelTextFeatureOverrides,
      bulk: app.labelTextBulkOverrides
    }),
    legend: JSON.stringify(app.legendEntries),
    results: app.results.map((entry) => [entry.name, entry.content.length, hash(entry.content)]),
    history: window.__GBDRAW_HISTORY__?.getUndoCount?.() ?? null
  };
});

const firstFeatureFillPicker = (page) => (
  page.locator('.right-drawer input[type="color"][aria-label^="Fill color for "]').first()
);

const changePicker = async (picker, color) => {
  await picker.evaluate((element, value) => {
    element.value = value;
    element.dispatchEvent(new Event('change', { bubbles: true }));
  }, color);
};

const scopeDialog = (page) => (
  page.getByRole('heading', { name: 'Feature Fill Scope' }).locator('xpath=..')
);

const strokeScopeDialog = (page) => (
  page.getByRole('heading', { name: 'Feature Stroke Scope' }).locator('xpath=..')
);

const labelScopeDialog = (page) => (
  page.getByRole('heading', { name: 'Feature Label Scope' }).locator('xpath=..')
);

test('inherited Feature fill stays editable and Cancel is mutation-free', async ({ page }) => {
  const errors = [];
  page.on('pageerror', (error) => errors.push(error.message));
  page.on('console', (message) => {
    if (message.type() === 'error') errors.push(message.text());
  });
  await loadSession(page);
  const picker = firstFeatureFillPicker(page);
  await expect(picker).toBeVisible();
  await expect(picker).toBeEnabled();
  await expect(picker.locator('xpath=../..').locator('[data-feature-fill-origin]'))
    .toContainText('Inherited from palette');

  const before = await documentSnapshot(page);
  await changePicker(picker, '#ff00ff');
  const dialog = scopeDialog(page);
  await expect(dialog).toBeVisible();
  await expect(dialog.getByRole('button', { name: /Update all .* features/ })).toBeVisible();
  await expect(dialog.getByRole('button', { name: /Update legend caption/ })).toBeVisible();
  await expect(dialog.getByRole('button', { name: /Update only this feature/ })).toBeVisible();
  expect(await documentSnapshot(page)).toEqual(before);

  await dialog.getByRole('button', { name: 'Cancel' }).click();
  await expect(dialog).toBeHidden();
  expect(await documentSnapshot(page)).toEqual(before);
  expect(errors).toEqual([]);
});

test('real majanivirus Feature click exposes applicable type, source-label, Similarity, and single choices', async ({ page }) => {
  test.setTimeout(180000);
  page.on('dialog', (dialog) => dialog.accept());
  await page.goto('/gbdraw/web/index.html', { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);
  const status = await page.evaluate(async () => {
    const name = 'majanivirus_orthogroup.gbdraw-session.json.gz';
    const response = await fetch(`/gbdraw/web/gallery/sessions/${name}`);
    const file = new File([await response.arrayBuffer()], name, { type: 'application/gzip' });
    return (await window.__GBDRAW_APP__.importSession({
      target: { files: [file], value: 'selected' }
    }))?.status;
  });
  expect(status).toBe('ok');
  const target = await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    const counts = (values) => values.reduce((result, value) => {
      const key = String(value || '').trim();
      if (key) result.set(key, (result.get(key) || 0) + 1);
      return result;
    }, new Map());
    const sourceLabel = (feature) => (
      feature.product || feature.gene || feature.locus_tag || feature.note || ''
    );
    const typeCounts = counts(app.extractedFeatures.map((feature) => feature.type));
    const sourceLabelCounts = counts(app.extractedFeatures.map(sourceLabel));
    const groupedSvgIds = new Set(app.orthogroups.flatMap((group) => (
      (group.members || []).map((member) => String(member.renderedFeatureSvgId || ''))
    )).filter(Boolean));
    return app.extractedFeatures.map((feature) => {
      const svgId = String(feature.svg_id || feature.renderedFeatureSvgId || '');
      const annotationLabel = sourceLabel(feature);
      const element = svgId
        ? document.querySelector(`[data-gbdraw-feature-id="${CSS.escape(svgId)}"]`)
        : null;
      return { svgId, annotationLabel, type: feature.type, element: Boolean(element) };
    }).find((feature) => (
      feature.element
      && groupedSvgIds.has(feature.svgId)
      && feature.annotationLabel
      && (typeCounts.get(String(feature.type || '').trim()) || 0) > 1
      && (sourceLabelCounts.get(feature.annotationLabel) || 0) > 1
    )) || null;
  });
  expect(target).not.toBeNull();

  const targetElement = page.locator(
    `[data-gbdraw-feature-id="${target.svgId}"]`
  ).first();
  await expect(targetElement).toBeVisible();
  await page.waitForFunction((svgId) => (
    document.querySelector(`[data-gbdraw-feature-id="${CSS.escape(svgId)}"]`)?.style?.cursor
      === 'pointer'
  ), target.svgId);
  await targetElement.dispatchEvent('click', { clientX: 500, clientY: 300 });
  const popup = page.getByRole('dialog', { name: /Feature details:/ });
  await expect(popup).toBeVisible();
  const picker = popup.getByLabel('Feature fill color', { exact: true }).first();
  await expect(picker).toBeEnabled();
  await changePicker(picker, '#ff00ff');

  const scoped = await page.evaluate(() => (
    window.__GBDRAW_APP__.pendingFeatureFillPlan?.candidates?.map(
      ({ semanticScope, targetCount }) => [semanticScope, targetCount]
    ) || []
  ));
  const scopes = scoped.map(([scope]) => scope);
  expect(scopes).toContain('feature-type');
  expect(scopes).toContain('source-annotation-label');
  expect(scopes).toContain('similarity-group');
  expect(scopes).toContain('single');
  expect(scoped.filter(([, targetCount]) => targetCount > 1).length).toBeGreaterThanOrEqual(3);
  const dialog = scopeDialog(page);
  await expect(dialog).toBeVisible();
  for (const label of [
    /Update all CDS features/,
    /Update legend caption/,
    /Update source label/,
    /Update similarity group/,
    /Update only this feature/
  ]) {
    await expect(dialog.getByRole('button', { name: label })).toBeVisible();
  }
  await expect(dialog.getByRole('button', { name: /Update matching rule/ })).toHaveCount(0);

  const similarityPlan = await page.evaluate(() => {
    const plan = window.__GBDRAW_APP__.pendingFeatureFillPlan;
    const candidate = plan?.candidates?.find(
      (entry) => entry.semanticScope === 'similarity-group'
    );
    return candidate ? {
      groupId: candidate.durableRuleIntent?.groupId,
      targetCount: candidate.targetCount
    } : null;
  });
  expect(similarityPlan).toEqual(expect.objectContaining({
    groupId: expect.any(String),
    targetCount: expect.any(Number)
  }));
  expect(similarityPlan.targetCount).toBeGreaterThan(1);

  await dialog.getByRole('button', { name: /Update similarity group/ }).first().click();
  if (await page.evaluate(() => window.__GBDRAW_APP__.featureFillPlanStatus === 'conflict')) {
    await dialog.getByRole('button', { name: /^Create “/ }).click();
  }
  await page.waitForFunction(() => (
    ['committed', 'failed', 'stale'].includes(window.__GBDRAW_APP__.featureFillPlanStatus)
  ));

  const literal = `fs1:similarity-group:${encodeURIComponent(similarityPlan.groupId)}`;
  const inspectSimilarityResult = () => page.evaluate(({ selector, targetCount }) => {
    const app = window.__GBDRAW_APP__;
    const rules = app.manualSpecificRules.filter((rule) => (
      rule.qual === '__gbdraw_semantic_scope__' && rule.val === selector
    ));
    const magentaFeatures = new Set();
    let magentaLegendEntries = 0;
    let magentaElementCount = 0;
    const fillSamples = new Set();
    (app.results || []).forEach((result, resultIndex) => {
      const svg = new DOMParser().parseFromString(result.content || '', 'image/svg+xml');
      svg.querySelectorAll('[fill]').forEach((element) => {
        const fill = String(element.getAttribute('fill') || '').toLowerCase();
        if (fill) fillSamples.add(fill);
        if (fill === '#ff00ff') magentaElementCount += 1;
      });
      svg.querySelectorAll('[data-gbdraw-feature-id]').forEach((root) => {
        const colored = [root, ...root.querySelectorAll('[fill]')].some(
          (element) => String(element.getAttribute('fill') || '').toLowerCase() === '#ff00ff'
        );
        if (colored) {
          magentaFeatures.add(`${resultIndex}:${root.getAttribute('data-gbdraw-feature-id')}`);
        }
      });
      magentaLegendEntries += [...svg.querySelectorAll('[data-legend-key]')].filter((entry) => (
        [entry, ...entry.querySelectorAll('[fill]')].some(
          (element) => String(element.getAttribute('fill') || '').toLowerCase() === '#ff00ff'
        )
      )).length;
    });
    return {
      status: app.featureFillPlanStatus,
      rules: rules.map(({ feat, qual, val, color, cap }) => ({ feat, qual, val, color, cap })),
      exactRuleCount: app.manualSpecificRules.filter(
        (rule) => rule.qual === '__gbdraw_instance_hash__'
      ).length,
      magentaFeatureCount: magentaFeatures.size,
      magentaLegendEntries,
      magentaElementCount,
      fillSamples: [...fillSamples].slice(0, 30),
      expectedTargetCount: targetCount
    };
  }, { selector: literal, targetCount: similarityPlan.targetCount });

  const applied = await inspectSimilarityResult();
  expect(applied.status).toBe('committed');
  expect(applied.rules).toHaveLength(1);
  expect(applied.rules[0]).toEqual(expect.objectContaining({
    feat: '*',
    qual: '__gbdraw_semantic_scope__',
    val: literal,
    color: '#ff00ff'
  }));
  expect(applied.exactRuleCount).toBe(0);
  expect(applied.magentaFeatureCount, JSON.stringify(applied, null, 2))
    .toBeGreaterThanOrEqual(similarityPlan.targetCount);
  expect(applied.magentaLegendEntries).toBeGreaterThan(0);

  const downloadPromise = page.waitForEvent('download', { timeout: 120000 });
  expect((await page.evaluate(() => window.__GBDRAW_APP__.saveSessionWithTitle())).status)
    .toBe('saved');
  const download = await downloadPromise;
  const savedPath = await download.path();
  expect(savedPath).toBeTruthy();
  const savedBuffer = readFileSync(savedPath);
  const savedSession = JSON.parse(gunzipSync(savedBuffer).toString('utf8'));
  expect(JSON.stringify(savedSession)).toContain(literal);

  await page.reload({ waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);
  const sessionInput = page.locator(
    'input[type="file"][accept*="application/json"][accept*="application/gzip"]'
  );
  await sessionInput.setInputFiles({
    name: download.suggestedFilename(),
    mimeType: 'application/gzip',
    buffer: savedBuffer
  });
  await page.waitForFunction((selector) => (
    window.__GBDRAW_APP__.manualSpecificRules.some((rule) => (
      rule.qual === '__gbdraw_semantic_scope__' && rule.val === selector
    ))
  ), literal, { timeout: 120000 });

  const regeneratedStatus = await page.evaluate(() => window.__GBDRAW_APP__.runAnalysis());
  expect(regeneratedStatus).toEqual({ status: 'ok' });
  const regenerated = await inspectSimilarityResult();
  expect(regenerated.rules).toHaveLength(1);
  expect(regenerated.exactRuleCount).toBe(0);
  expect(regenerated.magentaFeatureCount).toBeGreaterThanOrEqual(similarityPlan.targetCount);
  expect(regenerated.magentaLegendEntries).toBeGreaterThan(0);
});

test('Reset Settings invalidates an open real Feature fill plan', async ({ page }) => {
  await loadSession(page);
  const target = page.locator('[aria-label="Result Preview"] [data-gbdraw-feature-id]').first();
  await expect(target).toBeVisible();
  await expect.poll(() => target.evaluate((element) => element.style.cursor)).toBe('pointer');
  await target.dispatchEvent('click', { clientX: 500, clientY: 300 });
  const popup = page.getByRole('dialog', { name: /Feature details:/ });
  await expect(popup).toBeVisible();
  await changePicker(popup.getByLabel('Feature fill color', { exact: true }).first(), '#ff00ff');
  const dialog = scopeDialog(page);
  await expect(dialog).toBeVisible();
  const before = await page.evaluate(() => ({
    token: window.__GBDRAW_APP__.pendingFeatureFillPlan?.token,
    candidateId: window.__GBDRAW_APP__.pendingFeatureFillPlan?.candidates?.[0]?.id,
    results: JSON.stringify(window.__GBDRAW_APP__.results),
    rules: JSON.stringify(window.__GBDRAW_APP__.manualSpecificRules)
  }));
  expect(before.token).toBeTruthy();
  expect(before.candidateId).toBeTruthy();

  const resetButton = page.getByRole('button', { name: 'Reset Settings' });
  await expect(resetButton).toBeVisible();
  await resetButton.dispatchEvent('click');
  await expect(dialog).toBeHidden();
  const reset = await page.evaluate(async ({ token, candidateId }) => {
    const app = window.__GBDRAW_APP__;
    const beforeStaleResolve = {
      results: JSON.stringify(app.results),
      rules: JSON.stringify(app.manualSpecificRules)
    };
    const staleAccepted = await app.resolveFeatureFillScope(token, candidateId);
    return {
      staleAccepted,
      pending: app.pendingFeatureFillPlan,
      status: app.featureFillPlanStatus,
      results: app.results.length,
      catalogCleared: app.featureCatalog == null && app.extractedFeatures.length === 0,
      undoCount: window.__GBDRAW_HISTORY__?.getUndoCount?.() ?? null,
      redoCount: window.__GBDRAW_HISTORY__?.getRedoCount?.() ?? null,
      beforeStaleResolve,
      afterStaleResolve: {
        results: JSON.stringify(app.results),
        rules: JSON.stringify(app.manualSpecificRules)
      }
    };
  }, before);
  expect(reset.pending).toBeNull();
  expect(reset.status).toBe('idle');
  expect(reset.results).toBe(0);
  expect(reset.catalogCleared).toBe(true);
  expect(reset.undoCount).toBe(0);
  expect(reset.redoCount).toBe(0);
  expect(reset.staleAccepted).toBe(false);
  expect(reset.afterStaleResolve).toEqual(reset.beforeStaleResolve);
});

test('confirmed legend-caption fill commits one converged History command', async ({ page }) => {
  const errors = [];
  page.on('pageerror', (error) => errors.push(error.message));
  await loadSession(page);
  const picker = firstFeatureFillPicker(page);
  const before = await documentSnapshot(page);

  await changePicker(picker, '#ff00ff');
  const dialog = scopeDialog(page);
  await expect(dialog).toBeVisible();
  await dialog.getByRole('button', { name: /Update legend caption/ }).click();
  await page.waitForFunction(() => (
    ['committed', 'failed', 'stale'].includes(window.__GBDRAW_APP__.featureFillPlanStatus)
  ));

  const after = await documentSnapshot(page);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.featureFillPlanStatus)).toBe('committed');
  expect(after.rules).not.toBe(before.rules);
  expect(after.overrides).not.toBe(before.overrides);
  expect(after.legend).not.toBe(before.legend);
  expect(after.results).not.toEqual(before.results);
  expect(after.history).toBe(1);
  expect(await page.evaluate(() => (
    window.__GBDRAW_APP__.legendEntries.find((entry) => entry.caption === 'tRNA')?.color
  ))).toBe('#ff00ff');
  expect(errors).toEqual([]);
});

test('inherited Feature stroke is editable and the explicit group command is atomic', async ({ page }) => {
  const errors = [];
  page.on('pageerror', (error) => errors.push(error.message));
  await loadSession(page);
  await page.getByRole('button', { name: 'Edit', exact: true }).first().click();
  const picker = page.getByLabel('Feature stroke color');
  await expect(picker).toBeVisible();
  await expect(picker).toBeEnabled();
  await expect(page.locator('[data-feature-stroke-origin]')).toContainText('Inherited');
  const before = await documentSnapshot(page);

  await changePicker(picker, '#123456');
  const dialog = strokeScopeDialog(page);
  await expect(dialog).toBeVisible();
  await dialog.getByRole('button', { name: /Update legend caption/ }).click();
  await page.waitForFunction(() => (
    ['committed', 'failed', 'stale'].includes(window.__GBDRAW_APP__.featureStrokePlanStatus)
  ));

  const after = await documentSnapshot(page);
  expect(
    await page.evaluate(() => ({
      status: window.__GBDRAW_APP__.featureStrokePlanStatus,
      progress: window.__GBDRAW_APP__.featureStrokePlanProgress
    })),
    JSON.stringify(errors)
  ).toEqual(expect.objectContaining({ status: 'committed' }));
  expect(after.strokeModel).not.toBe(before.strokeModel);
  expect(after.results).not.toEqual(before.results);
  expect(after.history).toBe(before.history + 1);
  expect(errors).toEqual([]);
});

test('popup label text uses a planned all-Result command', async ({ page }) => {
  const errors = [];
  page.on('pageerror', (error) => errors.push(error.message));
  await loadSession(page);
  const opened = await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    const feature = app.extractedFeatures.find((entry) => (
      String(entry.product || entry.gene || entry.locus_tag || '') === 'tRNA-Leu'
    ));
    if (!feature) return false;
    app.openFeatureEditorFromList(feature, null);
    return true;
  });
  expect(opened).toBe(true);
  const input = page.locator('input[placeholder="Edit label text"]');
  await expect(input).toBeVisible();
  const before = await documentSnapshot(page);
  await input.fill('Leucine transfer RNA');
  await page.getByRole('button', { name: 'Apply Label' }).click();
  const dialog = labelScopeDialog(page);
  await expect(dialog).toBeVisible();
  await dialog.getByRole('button', { name: /Update rendered label/ }).click();
  await page.waitForFunction(() => (
    ['committed', 'failed', 'stale'].includes(window.__GBDRAW_APP__.featureLabelPlanStatus)
  ));
  const after = await documentSnapshot(page);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.featureLabelPlanStatus)).toBe('committed');
  expect(after.labelOverrides).not.toBe(before.labelOverrides);
  expect(after.results).not.toEqual(before.results);
  expect(after.history).toBe(before.history + 1);
  expect(errors).toEqual([]);
});

test('Circular batch group fill stages every Result before publishing', async ({ page }) => {
  test.setTimeout(120000);
  const errors = [];
  page.on('pageerror', (error) => errors.push(error.message));
  page.on('dialog', (dialog) => dialog.accept());
  await page.goto('/gbdraw/web/index.html', { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);

  const source = readFileSync(resolve('tests/test_inputs/NC_001416.gb'));
  const secondSource = readFileSync(resolve('tests/test_inputs/NC_001454.1.gbk'));
  await page.getByLabel('GenBank/DDBJ File', { exact: true }).setInputFiles({
    name: 'two-distinct-phages.gb',
    mimeType: 'text/plain',
    buffer: Buffer.concat([source, Buffer.from('\n'), secondSource])
  });
  await page.getByRole('button', { name: 'Generate Diagram' }).click();
  await page.waitForFunction(() => window.__GBDRAW_APP__.processing === true);
  await page.waitForFunction(() => window.__GBDRAW_APP__.processing === false, null, {
    timeout: 90000
  });
  const generation = await page.evaluate(() => ({
    results: window.__GBDRAW_APP__.results.length,
    error: window.__GBDRAW_APP__.errorLog
  }));
  expect(generation.results, JSON.stringify(generation.error)).toBe(2);

  await page.locator('.drawer-toggle').click();
  const picker = firstFeatureFillPicker(page);
  await expect(picker).toBeEnabled();
  const before = await documentSnapshot(page);
  await changePicker(picker, '#ff00ff');
  const dialog = scopeDialog(page);
  await expect(dialog).toBeVisible();
  await expect(dialog).toContainText('2 Result(s)');
  await dialog.getByRole('button', { name: /Update legend caption/ }).click();
  await page.waitForFunction(() => (
    ['committed', 'failed', 'stale'].includes(window.__GBDRAW_APP__.featureFillPlanStatus)
  ));

  const after = await documentSnapshot(page);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.featureFillPlanStatus)).toBe('committed');
  expect(after.results).not.toEqual(before.results);
  expect(after.results.every((entry, index) => entry[2] !== before.results[index][2])).toBe(true);
  const serializedFills = await page.evaluate(() => window.__GBDRAW_APP__.results.map((result) => ({
    name: result.name,
    magentaCount: (String(result.content || '').match(/#ff00ff/g) || []).length
  })));
  expect(serializedFills.every((entry) => entry.magentaCount > 0), JSON.stringify(serializedFills)).toBe(true);
  expect(after.history).toBe(before.history + 1);
  await page.getByLabel('Preview Result').selectOption('1');
  await page.waitForFunction(() => Number(window.__GBDRAW_APP__.selectedResultIndex) === 1);
  await page.waitForFunction(() => {
    const app = window.__GBDRAW_APP__;
    return app.legendEntries.some((entry) => (
      String(entry.owner || '').startsWith('feature')
      && Array.isArray(entry.featureIds)
      && entry.featureIds.length > 0
    ));
  });
  const switched = await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    const selectedResultIndex = Number(app.selectedResultIndex || 0);
    const selectedFeatureIds = new Set(app.extractedFeatures.filter((feature) => (
      Number(feature.resultIndex ?? feature.result_index) === selectedResultIndex
    )).map((feature) => String(
      feature.svgId
      ?? feature.svg_id
      ?? feature.renderedFeatureSvgId
      ?? feature.rendered_feature_svg_id
      ?? ''
    )).filter(Boolean));
    const featureLegendEntries = app.legendEntries.filter((entry) => (
      String(entry.owner || '').startsWith('feature')
    ));
    return {
      selectedResultIndex,
      svgContentMagentaCount: (String(app.svgContent || '').match(/#ff00ff/g) || []).length,
      mountedMagentaCount: (String(document.querySelector('[aria-label="Result Preview"] svg')?.outerHTML || '')
        .match(/#ff00ff/g) || []).length,
      featureLegendIdCount: featureLegendEntries.reduce(
        (total, entry) => total + (entry.featureIds || []).length,
        0
      ),
      invalidFeatureLegendIds: featureLegendEntries.flatMap((entry) => (
        (entry.featureIds || []).filter((featureId) => !selectedFeatureIds.has(String(featureId)))
      ))
    };
  });
  expect(switched.selectedResultIndex, JSON.stringify(switched)).toBe(1);
  expect(switched.svgContentMagentaCount, JSON.stringify(switched)).toBeGreaterThan(0);
  expect(switched.featureLegendIdCount, JSON.stringify(switched)).toBeGreaterThan(0);
  expect(switched.invalidFeatureLegendIds, JSON.stringify(switched)).toEqual([]);
  await expect.poll(() => page.locator('[aria-label="Result Preview"] svg [fill="#ff00ff"]').count())
    .toBeGreaterThan(0);

  const historyBeforeLegend = await page.evaluate(() => (
    window.__GBDRAW_HISTORY__?.getUndoCount?.() ?? null
  ));
  const legendDispatch = await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    const index = app.legendEntries.findIndex((entry) => entry.caption === 'CDS');
    const accepted = index >= 0
      ? await app.updateLegendEntryColor(index, '#00ffff')
      : false;
    return {
      accepted,
      status: app.featureFillPlanStatus,
      progress: app.featureFillPlanProgress,
      entry: index >= 0 ? app.legendEntries[index] : null
    };
  });
  expect(legendDispatch.accepted, JSON.stringify(legendDispatch)).toBe(true);
  await page.waitForFunction(() => (
    ['committed', 'failed', 'stale'].includes(window.__GBDRAW_APP__.featureFillPlanStatus)
  ));
  const legendEdit = await page.evaluate(() => ({
    status: window.__GBDRAW_APP__.featureFillPlanStatus,
    cyanCounts: window.__GBDRAW_APP__.results.map((result) => (
      (String(result.content || '').match(/#00ffff/g) || []).length
    )),
    history: window.__GBDRAW_HISTORY__?.getUndoCount?.() ?? null
  }));
  expect(legendEdit.status, JSON.stringify(legendEdit)).toBe('committed');
  expect(legendEdit.cyanCounts.every((count) => count > 0), JSON.stringify(legendEdit)).toBe(true);
  expect(legendEdit.history).toBe(historyBeforeLegend + 1);
  expect(errors).toEqual([]);
});
