const { test, expect } = require('@playwright/test');
const { readFileSync } = require('node:fs');
const { join, resolve } = require('node:path');
const { gunzipSync } = require('node:zlib');

const repoRoot = resolve(process.env.GBDRAW_REPO || process.cwd());
const bgcSessionPath = join(repoRoot, 'tests/test_inputs/BGC0000708-BGC0000713.gbdraw-session.json');
const hmmtDnaPath = join(repoRoot, 'tests/test_inputs/HmmtDNA.gbk');
const repeatRegionGenbankPath = join(repoRoot, 'tests/test_inputs/NC_001454.1.gbk');
const sparseGenbankAPath = join(repoRoot, 'tests/test_inputs/BGC0000708.gbk');
const sparseGenbankBPath = join(repoRoot, 'tests/test_inputs/BGC0000709.gbk');

test('diagram worker startup leaves the main helper runtime lazy', async ({ page }) => {
  test.setTimeout(180000);
  const genbank = readFileSync(hmmtDnaPath, 'utf8');
  await page.goto('/gbdraw/web/index.html', { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);
  await page.waitForFunction(
    () => (
      window.__GBDRAW_APP__?.diagramGenerationWorkerReady === true &&
      Object.keys(window.__GBDRAW_APP__?.paletteDefinitions || {}).length > 0
    ),
    null,
    { timeout: 180000 }
  );
  expect(await page.evaluate(() => ({
    pyodideReady: window.__GBDRAW_APP__.pyodideReady,
    repeatRegion: window.__GBDRAW_APP__.currentColors.repeat_region
  }))).toEqual({
    pyodideReady: false,
    repeatRegion: '#d3d3d3'
  });

  await page.evaluate((genbankText) => {
    const originalLoadPyodide = window.loadPyodide;
    if (typeof originalLoadPyodide !== 'function') {
      throw new Error('Main-thread Pyodide loader is unavailable.');
    }
    window.__GBDRAW_MAIN_PYODIDE_LOADS__ = 0;
    window.loadPyodide = (...args) => {
      window.__GBDRAW_MAIN_PYODIDE_LOADS__ += 1;
      return originalLoadPyodide(...args);
    };
    const app = window.__GBDRAW_APP__;
    app.mode = 'circular';
    app.cInputType = 'gb';
    app.files.c_gb = new File([genbankText], 'simple-worker-only.gbk', {
      type: 'text/plain',
      lastModified: 17
    });
  }, genbank);
  await page.waitForTimeout(250);
  expect(await page.evaluate(() => ({
    loadCalls: window.__GBDRAW_MAIN_PYODIDE_LOADS__,
    pyodideReady: window.__GBDRAW_APP__.pyodideReady
  }))).toEqual({
    loadCalls: 0,
    pyodideReady: false
  });

  expect(await runDiagramWithDiagnostics(page)).toEqual({
    result: { status: 'ok' },
    errorSummary: '',
    errorDetails: []
  });
  expect(await page.evaluate(() => ({
    loadCalls: window.__GBDRAW_MAIN_PYODIDE_LOADS__,
    pyodideReady: window.__GBDRAW_APP__.pyodideReady
  }))).toEqual({
    loadCalls: 0,
    pyodideReady: false
  });
});

const inspectSparseDepthResult = async (page) => page.evaluate(() => {
  const app = window.__GBDRAW_APP__;
  const content = app.results?.[0]?.content || '';
  const svg = new DOMParser().parseFromString(content, 'image/svg+xml').documentElement;
  const depthGroups = (slotId) => [
    ...svg.querySelectorAll(
      `g[data-gbdraw-slot-id="${slotId}"][data-gbdraw-slot-renderer="depth"]`
    )
  ];
  const recordIndexForGroup = (group) => {
    for (
      let sibling = group?.nextElementSibling;
      sibling;
      sibling = sibling.nextElementSibling
    ) {
      if (sibling.matches(
        'g[data-gbdraw-role="record-definition"][data-gbdraw-record-index]'
      )) {
        return Number.parseInt(sibling.getAttribute('data-gbdraw-record-index'), 10);
      }
    }
    return null;
  };
  const depthGroup = (slotId, recordIndex) => depthGroups(slotId)
    .find((group) => recordIndexForGroup(group) === recordIndex);
  const hasAxis = (group) => Boolean(group?.querySelector(':scope > g'));
  const groupFills = (group) => {
    if (!group) return [];
    return [group, ...group.querySelectorAll('[fill]')]
      .map((element) => String(element.getAttribute('fill') || '').toLowerCase())
      .filter(Boolean);
  };
  const depthARecord1 = depthGroup('depth_a', 0);
  const depthARecord2 = depthGroup('depth_a', 1);
  const depthBRecord1 = depthGroup('depth_b', 0);
  const depthBRecord2 = depthGroup('depth_b', 1);
  const args = Array.isArray(app.lastRunInfo?.invocation?.args)
    ? app.lastRunInfo.invocation.args
    : [];
  const depthArgs = [];
  args.forEach((arg, index) => {
    if (arg === '--depth_track') depthArgs.push(args.slice(index + 1, index + 3));
  });
  return {
    resultCount: app.results?.length || 0,
    groups: {
      depthARecord1: Boolean(depthARecord1),
      depthARecord1Axis: hasAxis(depthARecord1),
      depthARecord2: Boolean(depthARecord2),
      depthBRecord1: Boolean(depthBRecord1),
      depthBRecord2: Boolean(depthBRecord2),
      depthBRecord2Axis: hasAxis(depthBRecord2)
    },
    depthAFills: groupFills(depthARecord1),
    depthBFills: groupFills(depthBRecord2),
    depthArgs
  };
});

const inspectCircularSparseDepthResult = async (page) => page.evaluate(() => {
  const app = window.__GBDRAW_APP__;
  const content = app.results?.[0]?.content || '';
  const svg = new DOMParser().parseFromString(content, 'image/svg+xml').documentElement;
  const depthGroups = (slotId) => [
    ...svg.querySelectorAll(
      `g[data-gbdraw-slot-id="${slotId}"][data-gbdraw-slot-renderer="depth"]`
    )
  ];
  const recordIndexForGroup = (group) => {
    const featureGroup = group?.parentElement?.querySelector(
      'g[data-gbdraw-slot-renderer="features"][data-gbdraw-record-index]'
    );
    return featureGroup
      ? Number.parseInt(featureGroup.getAttribute('data-gbdraw-record-index'), 10)
      : null;
  };
  const depthGroup = (slotId, recordIndex) => depthGroups(slotId)
    .find((group) => recordIndexForGroup(group) === recordIndex);
  const hasAxis = (group) => Boolean(group?.querySelector(':scope > g'));
  const groupFills = (group) => {
    if (!group) return [];
    return [group, ...group.querySelectorAll('[fill]')]
      .map((element) => String(element.getAttribute('fill') || '').toLowerCase())
      .filter(Boolean);
  };
  const depthARecord1 = depthGroup('depth_1', 0);
  const depthARecord2 = depthGroup('depth_1', 1);
  const depthBRecord1 = depthGroup('depth_2', 0);
  const depthBRecord2 = depthGroup('depth_2', 1);
  const args = Array.isArray(app.lastRunInfo?.invocation?.args)
    ? app.lastRunInfo.invocation.args
    : [];
  const depthArgs = [];
  args.forEach((arg, index) => {
    if (arg === '--depth_track') depthArgs.push(args.slice(index + 1, index + 3));
  });
  return {
    groups: {
      depthARecord1: Boolean(depthARecord1),
      depthARecord1Axis: hasAxis(depthARecord1),
      depthARecord2: Boolean(depthARecord2),
      depthBRecord1: Boolean(depthBRecord1),
      depthBRecord2: Boolean(depthBRecord2),
      depthBRecord2Axis: hasAxis(depthBRecord2)
    },
    depthAFills: groupFills(depthARecord1),
    depthBFills: groupFills(depthBRecord2),
    depthArgs
  };
});

const runDiagramWithDiagnostics = async (page) => page.evaluate(async () => {
  const app = window.__GBDRAW_APP__;
  const result = await app.runAnalysis();
  return {
    result,
    errorSummary: String(app.errorLog?.summary || ''),
    errorDetails: Array.isArray(app.errorLog?.details)
      ? app.errorLog.details.map((detail) => String(detail))
      : []
  };
});

test('Circular GFF3 mode exposes one annotation and one FASTA uploader', async ({ page }) => {
  await page.goto('/gbdraw/web/index.html', { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);

  await expect(page.getByLabel('GenBank/DDBJ File', { exact: true })).toHaveCount(1);
  await page.getByLabel('GFF3 + FASTA', { exact: true }).check();

  await expect(page.getByLabel('GenBank/DDBJ File', { exact: true })).toHaveCount(0);
  await expect(page.getByLabel('GFF3 File', { exact: true })).toHaveCount(1);
  await expect(page.getByLabel('FASTA File', { exact: true })).toHaveCount(1);
});

test('Show Depth stays disabled until a depth TSV is uploaded', async ({ page }) => {
  await page.goto('/gbdraw/web/index.html', { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);

  const showDepthCheckbox = page.locator('label:has-text("Show Depth") input[type="checkbox"]').first();
  await expect(showDepthCheckbox).toBeDisabled();
  await expect(showDepthCheckbox).not.toBeChecked();

  await page.evaluate(() => {
    window.__GBDRAW_APP__.addCircularDepthTrack();
  });
  await expect(showDepthCheckbox).toBeDisabled();
  await expect(showDepthCheckbox).not.toBeChecked();

  await page.evaluate(() => {
    const file = new File(['position\tdepth\n1\t12\n'], 'depth.tsv', {
      type: 'text/tab-separated-values'
    });
    window.__GBDRAW_APP__.setCircularDepthFile(0, file);
  });
  await expect(showDepthCheckbox).toBeEnabled();
  await expect(showDepthCheckbox).toBeChecked();

  await page.evaluate(() => {
    window.__GBDRAW_APP__.setCircularDepthFile(0, null);
  });
  await expect(showDepthCheckbox).toBeDisabled();
  await expect(showDepthCheckbox).not.toBeChecked();
});

test('Linear Depth above Features generates with repeat_region underlays', async ({ page }) => {
  test.setTimeout(300000);
  const genbank = readFileSync(repeatRegionGenbankPath, 'utf8');
  await page.goto('/gbdraw/web/index.html', { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);
  await page.waitForFunction(() => (
    window.__GBDRAW_APP__?.diagramGenerationWorkerReady === true
  ), null, { timeout: 180000 });

  const outcome = await page.evaluate(async (genbankText) => {
    const app = window.__GBDRAW_APP__;
    const { state } = await import('./js/state.js');
    app.mode = 'linear';
    app.lInputType = 'gb';
    app.setLinearSeqPrimaryFile(0, 'gb', new File([genbankText], 'repeat-region.gbk', {
      type: 'text/plain',
      lastModified: 19
    }));
    app.setLinearDepthFile(app.linearSeqs[0], 0, new File([
      'reference_name\tposition\tdepth\n' +
      'NC_001454.1\t1\t12\n' +
      'NC_001454.1\t80\t18\n' +
      'NC_001454.1\t160\t9\n'
    ], 'repeat-region-depth.tsv', {
      type: 'text/tab-separated-values',
      lastModified: 20
    }));
    Object.assign(app.form, {
      show_gc: false,
      show_skew: false,
      show_depth: true,
      show_labels_linear: 'none',
      legend: 'none'
    });
    if (!app.adv.features.includes('repeat_region')) {
      app.adv.features.push('repeat_region');
    }
    app.adv.feature_shapes = {
      ...(app.adv.feature_shapes || {}),
      repeat_region: 'underlay'
    };
    app.adv.linear_track_slots.splice(
      0,
      app.adv.linear_track_slots.length,
      {
        id: 'depth_3',
        renderer: 'depth',
        enabled: true,
        side: 'above',
        z: 0,
        params: { track_index: 0 }
      },
      {
        id: 'features',
        renderer: 'features',
        enabled: true,
        side: 'overlay',
        z: 0,
        params: {}
      }
    );
    app.adv.linear_track_slots_axis_index = 1;
    app.adv.linear_track_slots_enabled = true;
    app.results.splice(0, app.results.length);
    state.trackSlotResolvedGeometry.value = null;

    const result = await app.runAnalysis();
    const geometry = state.trackSlotResolvedGeometry.value;
    const firstRecord = geometry?.records?.find((record) => (
      Number(record?.resultIndex ?? 0) === 0 &&
      Number(record?.recordIndex ?? 0) === 0
    ));
    const slots = firstRecord?.slots || [];
    const content = app.results?.[0]?.content || '';
    const svg = new DOMParser().parseFromString(content, 'image/svg+xml');
    const underlayPaintIndex = content.indexOf(
      'data-gbdraw-slot-id="__gbdraw_auto_feature_underlay_slot__"'
    );
    const foregroundPaintIndex = content.indexOf(
      'data-gbdraw-slot-id="features"'
    );
    return {
      result,
      errorSummary: String(app.errorLog?.summary || ''),
      errorDetails: Array.isArray(app.errorLog?.details)
        ? app.errorLog.details.map((detail) => String(detail))
        : [],
      resultCount: app.results?.length || 0,
      mode: geometry?.mode || null,
      slotIds: slots.map((slot) => slot.slotId),
      depthSide: slots.find((slot) => slot.slotId === 'depth_3')?.side || null,
      featureSide: slots.find((slot) => slot.slotId === 'features')?.side || null,
      underlayCount: svg.querySelectorAll(
        '[data-gbdraw-auto-feature-underlay="true"]'
      ).length,
      underlayPaintsFirst: (
        underlayPaintIndex >= 0 &&
        foregroundPaintIndex >= 0 &&
        underlayPaintIndex < foregroundPaintIndex
      )
    };
  }, genbank);

  expect(outcome).toMatchObject({
    result: { status: 'ok' },
    errorSummary: '',
    errorDetails: [],
    resultCount: 1,
    mode: 'linear',
    depthSide: 'above',
    featureSide: 'overlay',
    underlayPaintsFirst: true
  });
  expect(outcome.slotIds).toEqual([
    '__gbdraw_auto_feature_underlay_slot__',
    'depth_3',
    'features'
  ]);
  expect(outcome.underlayCount).toBeGreaterThan(0);
});

test('Linear depth add, clear, and remove keep global sparse columns aligned', async ({ page }) => {
  await page.goto('/gbdraw/web/index.html', { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);

  const result = await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    app.mode = 'linear';
    app.addLinearSeq();
    const first = new File(['position\tdepth\n1\t10\n'], 'sample-a.tsv', {
      type: 'text/tab-separated-values'
    });
    const second = new File(['position\tdepth\n1\t20\n'], 'sample-b.tsv', {
      type: 'text/tab-separated-values'
    });
    app.setLinearDepthFile(app.linearSeqs[0], 0, first);
    app.resetLinearTrackSlotsFromSimpleControls();
    app.setLinearTrackSlotsEnabled(true);
    const clonePlain = (value) => JSON.parse(JSON.stringify(value));
    const originalSlots = clonePlain(app.adv.linear_track_slots);
    const originalAxisIndex = app.adv.linear_track_slots_axis_index;
    const originalFeatureSlot = clonePlain(
      app.adv.linear_track_slots.find((slot) => slot.renderer === 'features')
    );
    const duplicateManagedSlot = clonePlain(
      app.adv.linear_track_slots.find((slot) => (
        slot.renderer === 'depth' && slot.params.track_index === 0
      ))
    );
    app.adv.linear_track_slots.splice(
      0,
      app.adv.linear_track_slots.length,
      {
        id: 'manual_depth', renderer: 'depth', enabled: true, side: 'above',
        params: { track_index: 0, custom: 'manual' }
      },
      duplicateManagedSlot,
      originalFeatureSlot
    );
    app.adv.linear_track_slots_axis_index = 2;
    app.ensureLinearTrackDepthSlots();
    const deduplicatedSlots = {
      ids: app.adv.linear_track_slots.map((slot) => slot.id),
      axisIndex: app.adv.linear_track_slots_axis_index
    };
    app.adv.linear_track_slots.splice(0, app.adv.linear_track_slots.length, ...originalSlots);
    app.adv.linear_track_slots_axis_index = originalAxisIndex;
    app.addLinearDepthTrack();
    app.setLinearDepthFile(app.linearSeqs[1], 1, second);
    await new Promise((resolve) => requestAnimationFrame(() => requestAnimationFrame(resolve)));

    const beforeClear = {
      rows: app.linearSeqs.map((seq) => seq.depth.map((file) => file?.name || null)),
      labels: app.adv.depth_tracks.map((track) => track.label),
      slotIndexes: app.adv.linear_track_slots
        .filter((slot) => slot.renderer === 'depth' && slot.enabled !== false)
        .map((slot) => slot.params.track_index)
    };
    app.setLinearDepthFile(app.linearSeqs[0], 0, null);
    const afterClear = {
      rows: app.linearSeqs.map((seq) => seq.depth.map((file) => file?.name || null)),
      labels: app.adv.depth_tracks.map((track) => track.label),
      slotIndexes: app.adv.linear_track_slots
        .filter((slot) => slot.renderer === 'depth' && slot.enabled !== false)
        .map((slot) => slot.params.track_index)
    };
    const featureSlot = app.adv.linear_track_slots.find((slot) => slot.renderer === 'features');
    const trackZeroSlot = app.adv.linear_track_slots.find((slot) => (
      slot.renderer === 'depth' && slot.params.track_index === 0
    ));
    const trackOneSlot = app.adv.linear_track_slots.find((slot) => (
      slot.renderer === 'depth' && slot.params.track_index === 1
    ));
    const manualSlotSource = {
      id: 'custom_depth',
      renderer: 'depth',
      enabled: true,
      side: 'above',
      params: { track_index: 0, custom: 'keep-manual' }
    };
    app.adv.linear_track_slots.splice(
      0,
      app.adv.linear_track_slots.length,
      manualSlotSource,
      trackZeroSlot,
      featureSlot,
      trackOneSlot
    );
    app.adv.linear_track_slots_axis_index = 2;
    app.removeLinearDepthTrack(app.linearSeqs[0], 0);
    const manualSlot = app.adv.linear_track_slots.find((slot) => slot.id === 'custom_depth');
    const { linearTrackAxisIndexForEnabledSlots } = await import(
      new URL('./js/app/linear-track-slots.js', window.location.href).href
    );
    const afterRemove = {
      rows: app.linearSeqs.map((seq) => seq.depth.map((file) => file?.name || null)),
      labels: app.adv.depth_tracks.map((track) => track.label),
      slotIndexes: app.adv.linear_track_slots
        .filter((slot) => slot.renderer === 'depth' && slot.enabled !== false)
        .map((slot) => slot.params.track_index),
      manualSlot: {
        enabled: manualSlot.enabled,
        trackIndex: manualSlot.params.track_index ?? null,
        error: manualSlot.depth_binding_error
      },
      fullAxisIndex: app.adv.linear_track_slots_axis_index,
      emittedAxisIndex: linearTrackAxisIndexForEnabledSlots(
        app.adv.linear_track_slots,
        app.adv.linear_track_slots_axis_index
      ),
      emittedSlotIds: app.adv.linear_track_slots
        .filter((slot) => slot.enabled !== false)
        .map((slot) => slot.id)
    };
    manualSlot.params.track_index = 0;
    app.syncDepthTrackSlotLabel(manualSlot);
    const afterRepair = {
      enabled: manualSlot.enabled,
      trackIndex: manualSlot.params.track_index,
      error: manualSlot.depth_binding_error ?? null
    };
    return { deduplicatedSlots, beforeClear, afterClear, afterRemove, afterRepair };
  });

  expect(result.deduplicatedSlots).toEqual({ ids: ['manual_depth', 'features'], axisIndex: 1 });
  expect(result.beforeClear.rows).toEqual([
    ['sample-a.tsv', null],
    [null, 'sample-b.tsv']
  ]);
  expect(result.beforeClear.slotIndexes).toEqual([0, 1]);
  expect(result.afterClear.rows).toEqual([
    [null, null],
    [null, 'sample-b.tsv']
  ]);
  expect(result.afterClear.labels).toEqual(result.beforeClear.labels);
  expect(result.afterClear.slotIndexes).toEqual(result.beforeClear.slotIndexes);
  expect(result.afterRemove.rows).toEqual([
    [null],
    ['sample-b.tsv']
  ]);
  expect(result.afterRemove.labels).toHaveLength(1);
  expect(result.afterRemove.slotIndexes).toEqual([0]);
  expect(result.afterRemove.manualSlot.enabled).toBe(false);
  expect(result.afterRemove.manualSlot.trackIndex).toBeNull();
  expect(result.afterRemove.manualSlot.error).toContain('logical track index 0');
  expect(result.afterRemove.fullAxisIndex).toBe(1);
  expect(result.afterRemove.emittedAxisIndex).toBe(0);
  expect(result.afterRemove.emittedSlotIds).toEqual(['features', 'depth_2']);
  expect(result.afterRepair).toEqual({ enabled: false, trackIndex: 0, error: null });
});

test('Linear custom slot panel and enable state preserve the explicit stack', async ({ page }) => {
  await page.goto('/gbdraw/web/index.html', { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);

  const result = await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    app.mode = 'linear';
    app.adv.linear_track_slots.splice(
      0,
      app.adv.linear_track_slots.length,
      {
        id: 'custom_gc', renderer: 'dinucleotide_content', enabled: true,
        side: 'above', height: '23px', spacing: '4px', params: { nt: 'AT' }
      },
      {
        id: 'features', renderer: 'features', enabled: true,
        side: 'overlay', height: '41px', spacing: '6px', params: {}
      }
    );
    app.adv.linear_track_slots_axis_index = 1;
    app.setLinearTrackSlotsEnabled(true);
    const snapshot = () => JSON.stringify({
      slots: app.adv.linear_track_slots,
      axis: app.adv.linear_track_slots_axis_index,
      specs: app.adv.linear_track_slots.map((slot) => app.linearTrackSlotCliSpec(slot))
    });
    const initial = snapshot();

    app.toggleLinearTrackSlotsPanel();
    const panelOpened = app.linearTrackSlotsPanelOpen;
    const afterOpen = snapshot();
    app.toggleLinearTrackSlotsPanel();
    const panelClosed = app.linearTrackSlotsPanelOpen;
    const afterClose = snapshot();

    app.setLinearTrackSlotsEnabled(false);
    const afterDisable = snapshot();
    app.form.show_gc = false;
    app.form.show_skew = true;
    app.form.show_depth = false;
    app.form.linear_track_layout = 'below';
    const afterSimpleChanges = snapshot();
    app.setLinearTrackSlotsEnabled(true);
    const afterReenable = snapshot();
    app.resetLinearTrackSlotsFromSimpleControls();
    const afterReset = snapshot();

    return {
      initial,
      panelOpened,
      panelClosed,
      afterOpen,
      afterClose,
      afterDisable,
      afterSimpleChanges,
      afterReenable,
      afterReset
    };
  });

  expect(result.panelOpened).toBe(true);
  expect(result.panelClosed).toBe(false);
  expect(result.afterOpen).toBe(result.initial);
  expect(result.afterClose).toBe(result.initial);
  expect(result.afterDisable).toBe(result.initial);
  expect(result.afterSimpleChanges).toBe(result.initial);
  expect(result.afterReenable).toBe(result.initial);
  expect(result.afterReset).not.toBe(result.initial);
});

test('Custom Track disclosure and editable IDs preserve transient row identity in both modes', async ({ page }) => {
  await page.goto('/gbdraw/web/index.html', { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);

  const result = await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    app.mode = 'linear';
    app.adv.linear_track_slots.splice(
      0,
      app.adv.linear_track_slots.length,
      {
        id: 'features',
        renderer: 'features',
        enabled: true,
        side: 'overlay',
        height: '',
        spacing: '',
        z: 0,
        params: {}
      },
      {
        id: 'linear_space_a',
        renderer: 'spacer',
        enabled: true,
        side: 'below',
        height: '8px',
        spacing: '',
        z: 0,
        params: {}
      },
      {
        id: 'linear_space_b',
        renderer: 'spacer',
        enabled: true,
        side: 'below',
        height: '9px',
        spacing: '',
        z: 0,
        params: {}
      }
    );
    app.adv.linear_track_slots_axis_index = 0;
    app.setLinearTrackSlotsEnabled(true);
    if (!app.linearTrackSlotsPanelOpen) app.toggleLinearTrackSlotsPanel();
    await window.Vue.nextTick();
    const linearKeyBefore = app.linearTrackSlotEditorKey(app.adv.linear_track_slots[0]);
    const linearInput = document.querySelector(
      '#linear-custom-track-slots-panel input[placeholder="slot_id"]'
    );
    linearInput.focus();
    linearInput.value = 'custom_features';
    linearInput.dispatchEvent(new Event('input', { bubbles: true }));
    await window.Vue.nextTick();
    const linearKeyAfter = app.linearTrackSlotEditorKey(app.adv.linear_track_slots[0]);
    const linearFeatures = app.adv.linear_track_slots[0];
    const linearSpaceA = app.adv.linear_track_slots[1];
    const linearSpaceB = app.adv.linear_track_slots[2];
    const linearOriginalKeys = new Map(
      [linearFeatures, linearSpaceA, linearSpaceB].map((slot) => [
        slot,
        app.linearTrackSlotEditorKey(slot)
      ])
    );
    app.moveLinearTrackSlot(2, 1);
    linearSpaceA.renderer = 'dinucleotide_content';
    app.updateLinearTrackSlotRenderer(linearSpaceA);
    app.duplicateLinearTrackSlot(1);
    const linearDuplicate = app.adv.linear_track_slots[2];
    const linearDuplicateKey = app.linearTrackSlotEditorKey(linearDuplicate);
    app.removeLinearTrackSlot(2);
    const linearLifecycleKeysPreserved = (
      app.adv.linear_track_slots.includes(linearFeatures) &&
      app.adv.linear_track_slots.includes(linearSpaceA) &&
      app.adv.linear_track_slots.includes(linearSpaceB) &&
      [linearFeatures, linearSpaceA, linearSpaceB].every(
        (slot) => app.linearTrackSlotEditorKey(slot) === linearOriginalKeys.get(slot)
      ) &&
      ![...linearOriginalKeys.values()].includes(linearDuplicateKey)
    );

    app.mode = 'circular';
    app.adv.circular_track_slots.splice(
      0,
      app.adv.circular_track_slots.length,
      {
        id: 'features',
        renderer: 'features',
        enabled: true,
        side: 'inside',
        width: null,
        radius: null,
        inner_gap_px: null,
        outer_gap_px: null,
        z: 0,
        params: { lane_direction: 'inside' }
      },
      {
        id: 'circular_space_a',
        renderer: 'spacer',
        enabled: true,
        side: 'inside',
        width: '8px',
        radius: null,
        inner_gap_px: null,
        outer_gap_px: null,
        z: 0,
        params: {}
      },
      {
        id: 'circular_space_b',
        renderer: 'spacer',
        enabled: true,
        side: 'inside',
        width: '9px',
        radius: null,
        inner_gap_px: null,
        outer_gap_px: null,
        z: 0,
        params: {}
      }
    );
    app.adv.circular_track_slots_axis_index = 0;
    app.setCircularTrackSlotsEnabled(true);
    if (!app.circularTrackSlotsPanelOpen) app.toggleCircularTrackSlotsPanel();
    await window.Vue.nextTick();
    const circularKeyBefore = app.circularTrackSlotEditorKey(
      app.adv.circular_track_slots[0]
    );
    const beforeDisclosure = JSON.stringify(app.adv.circular_track_slots);
    for (let index = 0; index < 10; index += 1) {
      app.toggleCircularTrackSlotsPanel();
      app.toggleCircularTrackSlotsPanel();
    }
    app.setCircularTrackSlotsEnabled(false);
    app.setCircularTrackSlotsEnabled(true);
    await window.Vue.nextTick();
    const disclosurePreserved =
      beforeDisclosure === JSON.stringify(app.adv.circular_track_slots);
    const circularInput = document.querySelector(
      '#circular-custom-track-slots-panel input[placeholder="slot_id"]'
    );
    circularInput.focus();
    circularInput.value = 'custom_circular_features';
    circularInput.dispatchEvent(new Event('input', { bubbles: true }));
    await window.Vue.nextTick();
    const circularKeyAfter = app.circularTrackSlotEditorKey(
      app.adv.circular_track_slots[0]
    );
    const circularFeatures = app.adv.circular_track_slots[0];
    const circularSpaceA = app.adv.circular_track_slots[1];
    const circularSpaceB = app.adv.circular_track_slots[2];
    const circularOriginalKeys = new Map(
      [circularFeatures, circularSpaceA, circularSpaceB].map((slot) => [
        slot,
        app.circularTrackSlotEditorKey(slot)
      ])
    );
    app.moveCircularTrackSlot(2, 1);
    app.updateCircularTrackSlotRenderer(circularSpaceA, 'dinucleotide_content');
    app.duplicateCircularTrackSlot(1);
    const circularDuplicate = app.adv.circular_track_slots[2];
    const circularDuplicateKey = app.circularTrackSlotEditorKey(circularDuplicate);
    app.removeCircularTrackSlot(2);
    const circularLifecycleKeysPreserved = (
      app.adv.circular_track_slots.includes(circularFeatures) &&
      app.adv.circular_track_slots.includes(circularSpaceA) &&
      app.adv.circular_track_slots.includes(circularSpaceB) &&
      [circularFeatures, circularSpaceA, circularSpaceB].every(
        (slot) => (
          app.circularTrackSlotEditorKey(slot) === circularOriginalKeys.get(slot)
        )
      ) &&
      ![...circularOriginalKeys.values()].includes(circularDuplicateKey)
    );

    return {
      linearId: app.adv.linear_track_slots[0].id,
      linearKeyBefore,
      linearKeyAfter,
      linearLifecycleKeysPreserved,
      circularId: app.adv.circular_track_slots[0].id,
      circularKeyBefore,
      circularKeyAfter,
      circularLifecycleKeysPreserved,
      disclosurePreserved
    };
  });

  expect(result.linearId).toBe('custom_features');
  expect(result.linearKeyAfter).toBe(result.linearKeyBefore);
  expect(result.linearLifecycleKeysPreserved).toBe(true);
  expect(result.circularId).toBe('custom_circular_features');
  expect(result.circularKeyAfter).toBe(result.circularKeyBefore);
  expect(result.circularLifecycleKeysPreserved).toBe(true);
  expect(result.disclosurePreserved).toBe(true);
});

test('Color value controls preserve Auto, None, and named colors', async ({ page }) => {
  await page.goto('/gbdraw/web/index.html', { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);

  const result = await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    app.currentColors.repeat_region = 'red';
    app.adv.block_stroke_color = 'none';
    app.adv.line_stroke_color = null;
    await window.Vue.nextTick();

    const select = (label) => document.querySelector(`select[aria-label="${label} mode"]`);
    const picker = (label) => document.querySelector(`input[aria-label="${label}"]`);
    const changeMode = async (label, value) => {
      const control = select(label);
      control.value = value;
      control.dispatchEvent(new Event('change', { bubbles: true }));
      await window.Vue.nextTick();
    };
    const snapshot = (label) => ({
      mode: select(label)?.value || '',
      picker: picker(label)?.value || '',
      disabled: Boolean(picker(label)?.disabled)
    });

    const initial = {
      repeat: snapshot('repeat_region feature color'),
      block: snapshot('Block stroke color'),
      line: snapshot('Line stroke color')
    };

    await changeMode('repeat_region feature color', 'none');
    const repeatNone = app.currentColors.repeat_region;
    await changeMode('repeat_region feature color', 'auto');
    const repeatAuto = app.currentColors.repeat_region;

    app.adv.block_stroke_color = 'navy';
    await window.Vue.nextTick();
    const blockNamed = snapshot('Block stroke color');
    await changeMode('Block stroke color', 'none');
    const blockNone = app.adv.block_stroke_color;

    await changeMode('Line stroke color', 'color');
    const lineColor = app.adv.line_stroke_color;
    await changeMode('Line stroke color', 'auto');
    const lineAuto = app.adv.line_stroke_color;

    return {
      initial,
      repeatNone,
      repeatAuto,
      blockNamed,
      blockNone,
      lineColor,
      lineAuto
    };
  });

  expect(result.initial.repeat).toEqual({
    mode: 'color',
    picker: '#ff0000',
    disabled: false
  });
  expect(result.initial.block).toEqual({
    mode: 'none',
    picker: '#000000',
    disabled: true
  });
  expect(result.initial.line).toEqual({
    mode: 'auto',
    picker: '#000000',
    disabled: true
  });
  expect(result.repeatNone).toBe('none');
  expect(result.repeatAuto).toBeNull();
  expect(result.blockNamed).toEqual({
    mode: 'color',
    picker: '#000080',
    disabled: false
  });
  expect(result.blockNone).toBe('none');
  expect(result.lineColor).toBe('#000000');
  expect(result.lineAuto).toBeNull();
});

test('Definition line colors preserve raw values and present named, Auto, and None states', async ({ page }) => {
  await page.goto('/gbdraw/web/index.html', { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);

  const result = await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    app.mode = 'linear';
    app.adv.linear_definition_line_styles.name.fill = 'red';
    app.adv.linear_definition_line_styles.subtitle.fill = null;
    app.adv.linear_definition_line_styles.replicon.fill = 'none';
    await window.Vue.nextTick();

    const picker = (label) => document.querySelector(
      `input[aria-label="${label} definition line color"]`
    );
    const state = (label) => document.querySelector(
      `[aria-label="${label} definition line color state"]`
    );
    const initial = {
      namePicker: picker('Name / Species')?.value || '',
      nameState: state('Name / Species')?.textContent?.trim() || '',
      nameRaw: app.adv.linear_definition_line_styles.name.fill,
      subtitlePicker: Boolean(picker('Subtitle')),
      subtitleState: state('Subtitle')?.textContent?.trim() || '',
      subtitleRaw: app.adv.linear_definition_line_styles.subtitle.fill,
      repliconPicker: Boolean(picker('Replicon')),
      repliconState: state('Replicon')?.textContent?.trim() || '',
      repliconRaw: app.adv.linear_definition_line_styles.replicon.fill
    };

    const namePicker = picker('Name / Species');
    namePicker.value = '#00ff00';
    namePicker.dispatchEvent(new Event('input', { bubbles: true }));
    await window.Vue.nextTick();

    return {
      initial,
      editedNameRaw: app.adv.linear_definition_line_styles.name.fill
    };
  });

  expect(result.initial).toEqual({
    namePicker: '#ff0000',
    nameState: '',
    nameRaw: 'red',
    subtitlePicker: false,
    subtitleState: 'Auto',
    subtitleRaw: null,
    repliconPicker: false,
    repliconState: 'None',
    repliconRaw: 'none'
  });
  expect(result.editedNameRaw).toBe('#00ff00');
});

test('Invalid Annotation slot is rejected before worker startup and preserves committed state', async ({ page }) => {
  await page.addInitScript(() => {
    const nativePostMessage = Worker.prototype.postMessage;
    window.__GBDRAW_DIAGRAM_RUN_MESSAGES__ = 0;
    Worker.prototype.postMessage = function instrumentDiagramWorker(message, ...args) {
      if (
        message?.type === 'run' &&
        String(this.__gbdrawWorkerUrl || '').includes('diagram-generation-worker')
      ) {
        window.__GBDRAW_DIAGRAM_RUN_MESSAGES__ += 1;
      }
      return nativePostMessage.call(this, message, ...args);
    };
    const NativeWorker = Worker;
    window.Worker = class InstrumentedWorker extends NativeWorker {
      constructor(url, options) {
        super(url, options);
        this.__gbdrawWorkerUrl = String(url || '');
      }
    };
  });
  page.on('dialog', (dialog) => dialog.accept());
  await page.goto('/gbdraw/web/index.html', { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);
  await page.evaluate(() => {
    window.__GBDRAW_APP__.pyodideReady = true;
  });

  await page.locator('input[accept^=".json,"]').first().setInputFiles(bgcSessionPath);
  await page.waitForFunction(() => window.__GBDRAW_APP__?.results?.length > 0);

  const outcome = await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    const featureId = String(app.extractedFeatures?.[0]?.svg_id || 'transaction-feature');
    app.selectedFeatureIds = new Set([featureId]);
    app.selectedFeatureAnchorId = featureId;
    app.featureColorOverrides.__transaction_probe = {
      color: '#123456',
      caption: 'Transaction probe'
    };
    app.featureStrokeOverrides.__transaction_probe = {
      color: '#654321',
      width: 2
    };
    app.legendEntries.push({
      caption: 'Transaction probe',
      color: '#abcdef',
      showStroke: true,
      strokeColor: '#654321',
      strokeWidth: 2
    });
    app.annotationSets.splice(0, app.annotationSets.length);
    app.linearComparisonPlan.mode = 'none';
    app.mode = 'linear';
    app.adv.linear_track_slots.splice(
      0,
      app.adv.linear_track_slots.length,
      {
        id: 'features',
        renderer: 'features',
        enabled: true,
        side: 'overlay',
        height: '',
        spacing: '',
        z: 0,
        params: {}
      },
      {
        id: 'invalid_annotation',
        renderer: 'annotations',
        enabled: true,
        side: 'overlay',
        height: '',
        spacing: '',
        z: 1,
        params: {
          set_id: 'missing',
          anchor_slot: 'features',
          layer: 'foreground'
        }
      }
    );
    app.adv.linear_track_slots_axis_index = 0;
    app.setLinearTrackSlotsEnabled(true);
    await window.Vue.nextTick();
    app.editableLabels = [{
      key: 'transaction-label',
      idx: 1,
      text: 'Keep this label',
      sourceText: 'Original label',
      featureId,
      kind: 'regular',
      draftText: 'Keep this label'
    }];

    const snapshot = () => ({
      results: app.results.map((entry) => ({
        name: entry.name,
        type: entry.type,
        content: entry.content
      })),
      selectedResultIndex: app.selectedResultIndex,
      selectedFeatureIds: [...app.selectedFeatureIds].sort(),
      selectedFeatureAnchorId: app.selectedFeatureAnchorId,
      featureColorOverrides: JSON.parse(JSON.stringify(app.featureColorOverrides)),
      featureStrokeOverrides: JSON.parse(JSON.stringify(app.featureStrokeOverrides)),
      legendEntries: JSON.parse(JSON.stringify(app.legendEntries)),
      editableLabels: JSON.parse(JSON.stringify(app.editableLabels)),
      liveLegend: (
        document.querySelector('svg #legend, svg #feature_legend, svg [id*="legend" i]')
          ?.outerHTML || null
      )
    });
    const before = snapshot();
    const workerMessagesBefore = window.__GBDRAW_DIAGRAM_RUN_MESSAGES__;
    const result = await app.runAnalysis();
    const after = snapshot();
    const same = (field) => JSON.stringify(after[field]) === JSON.stringify(before[field]);

    return {
      result,
      errorSummary: String(app.errorLog?.summary || ''),
      diagramWorkerMessages:
        window.__GBDRAW_DIAGRAM_RUN_MESSAGES__ - workerMessagesBefore,
      beforeResultCount: before.results.length,
      resultsPreserved: same('results'),
      resultSelectionPreserved: same('selectedResultIndex'),
      featureSelectionPreserved:
        same('selectedFeatureIds') && same('selectedFeatureAnchorId'),
      fillOverridesPreserved: same('featureColorOverrides'),
      strokeOverridesPreserved: same('featureStrokeOverrides'),
      legendPreserved: same('legendEntries') && same('liveLegend'),
      editableLabelsPreserved: same('editableLabels')
    };
  });

  expect(outcome.result?.status).toBe('error');
  expect(outcome.errorSummary).toContain("references unknown set 'missing'");
  expect(outcome.diagramWorkerMessages).toBe(0);
  expect(outcome.beforeResultCount).toBeGreaterThan(0);
  expect(outcome.resultsPreserved).toBe(true);
  expect(outcome.resultSelectionPreserved).toBe(true);
  expect(outcome.featureSelectionPreserved).toBe(true);
  expect(outcome.fillOverridesPreserved).toBe(true);
  expect(outcome.strokeOverridesPreserved).toBe(true);
  expect(outcome.legendPreserved).toBe(true);
  expect(outcome.editableLabelsPreserved).toBe(true);
});

test('Session preflight rejects invalid canonical data without resetting live state', async ({ page }) => {
  await page.goto('/gbdraw/web/index.html', { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);

  await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    const { state } = await import('./js/state.js');
    app.mode = 'linear';
    await window.Vue.nextTick();
    app.sessionTitle = 'keep-live-state';
    app.adv.linear_track_slots.splice(
      0,
      app.adv.linear_track_slots.length,
      {
        id: 'keep_gc', renderer: 'dinucleotide_content', enabled: true,
        side: 'above', height: '23px', spacing: '4px', params: { nt: 'AT' }
      },
      {
        id: 'keep_features', renderer: 'features', enabled: true,
        side: 'overlay', height: '41px', spacing: '6px', params: {}
      }
    );
    app.adv.linear_track_slots_axis_index = 1;
    app.setLinearTrackSlotsEnabled(true);
    state.suppressCircularMultiRecordDefaults.value = true;
  });

  const snapshot = () => page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    const { state } = await import('./js/state.js');
    return {
      mode: app.mode,
      title: app.sessionTitle,
      enabled: app.adv.linear_track_slots_enabled,
      axisIndex: app.adv.linear_track_slots_axis_index,
      slots: JSON.parse(JSON.stringify(app.adv.linear_track_slots)),
      suppressDefaults: state.suppressCircularMultiRecordDefaults.value
    };
  });
  const before = await snapshot();
  const dialogs = [];
  page.on('dialog', async (dialog) => {
    dialogs.push(dialog.message());
    await dialog.accept();
  });
  const input = page.locator('input[accept^=".json,"]').first();

  const invalidCanonicalSession = {
    format: 'gbdraw-session',
    version: 32,
    renderRequest: {
      schema: 2,
      mode: 'linear',
      records: [],
      diagramOptions: {},
      output: { prefix: 'invalid' }
    },
    resources: {}
  };
  await input.setInputFiles({
    name: 'invalid-canonical.gbdraw-session.json',
    mimeType: 'application/json',
    buffer: Buffer.from(JSON.stringify(invalidCanonicalSession))
  });
  await expect.poll(() => dialogs.length).toBe(1);
  expect(dialogs[0]).toContain('Canonical renderRequest records are required.');
  expect(await snapshot()).toEqual(before);

  const invalidSlotSchemaSession = {
    ...invalidCanonicalSession,
    version: 33,
    config: {
      adv: {
        linear_track_slots_schema_version: '1',
        linear_track_slots: [{ id: 'features', renderer: 'features', params: {} }]
      }
    }
  };
  await input.setInputFiles({
    name: 'invalid-slot-schema.gbdraw-session.json',
    mimeType: 'application/json',
    buffer: Buffer.from(JSON.stringify(invalidSlotSchemaSession))
  });
  await expect.poll(() => dialogs.length).toBe(2);
  expect(dialogs[1]).toContain('Custom Track Slots use an obsolete schema.');
  expect(await snapshot()).toEqual(before);

  const recordText = `LOCUS       PREFLIGHT                   4 bp    DNA     linear   UNK 01-JAN-1980
DEFINITION  Session preflight fixture.
ACCESSION   PREFLIGHT
FEATURES             Location/Qualifiers
ORIGIN
        1 atgc
//
`;
  const unknownAuthoritySession = {
    format: 'gbdraw-session',
    version: 33,
    renderRequest: {
      schema: 2,
      mode: 'linear',
      records: [{
        recordKey: 'record-1',
        source: { kind: 'genbank', resourceId: 'record-1-genbank' },
        selector: null,
        region: null,
        presentation: {
          label: null, subtitle: null, reverseComplement: false,
          gridRow: null, gridColumn: null
        }
      }],
      diagramOptions: {
        configOverrides: {},
        tracks: {
          linearTrackSlots: ['features:features@side=overlay'],
          linearTrackAxisIndex: 0
        }
      },
      layout: {},
      comparisons: [],
      output: {
        prefix: 'preflight', formats: ['interactive_svg'], overwrite: true,
        interactiveMetadataPolicy: 'auto'
      }
    },
    resources: {
      'record-1-genbank': {
        kind: 'genbank', name: 'record.gb', type: 'text/plain',
        size: Buffer.byteLength(recordText), lastModified: 0,
        encoding: 'base64', data: Buffer.from(recordText).toString('base64')
      }
    },
    unknownSemanticState: { shouldNotBeAccepted: true },
    config: {
      adv: {
        linear_track_slots_enabled: true,
        linear_track_slots_schema_version: 2,
        linear_track_slots_axis_index: 1,
        linear_track_slots: [
          {
            id: 'invalid_depth', renderer: 'depth', enabled: true,
            side: 'above', z: 0, params: { track_index: 0.5 }
          },
          {
            id: 'features', renderer: 'features', enabled: true,
            side: 'overlay', z: 0, params: {}
          }
        ]
      }
    }
  };
  await input.setInputFiles({
    name: 'unknown-authority.gbdraw-session.json',
    mimeType: 'application/json',
    buffer: Buffer.from(JSON.stringify(unknownAuthoritySession))
  });
  await expect.poll(() => dialogs.length).toBe(3);
  expect(dialogs[2]).toContain('unclassified top-level field');
  expect(await snapshot()).toEqual(before);

  const invalidLegacyConfig = {
    form: { multi_record_canvas: true },
    adv: {
      circular_track_slots: [{ id: 'features', renderer: 'features', params: {} }]
    }
  };
  await input.setInputFiles({
    name: 'invalid-legacy-config.json',
    mimeType: 'application/json',
    buffer: Buffer.from(JSON.stringify(invalidLegacyConfig))
  });
  await expect.poll(() => dialogs.length).toBe(4);
  expect(dialogs[3]).toContain('Custom Track Slots use an obsolete schema.');
  expect(await snapshot()).toEqual(before);

  const repairedCanonicalSession = {
    ...unknownAuthoritySession,
    unknownSemanticState: undefined,
    renderRequest: {
      ...unknownAuthoritySession.renderRequest,
      diagramOptions: {
        configOverrides: { comparison_height: -2 },
        output: {
          outputPrefix: 'preflight',
          legend: 'bottom',
          plotTitlePosition: 'top'
        },
        tracks: {}
      }
    },
    config: {
      form: { legend: 'right' },
      adv: { comparison_height: 99, plot_title_position: 'bottom' }
    },
    ui: { mode: 'circular', legend: 'left', linearLegendPosition: 'right' }
  };
  expect(repairedCanonicalSession.renderRequest.diagramOptions.output.legend).toBe('bottom');
  await input.setInputFiles({
    name: 'repaired-canonical.gbdraw-session.json',
    mimeType: 'application/json',
    buffer: Buffer.from(JSON.stringify(repairedCanonicalSession))
  });
  await expect.poll(() => dialogs.length).toBe(5);
  expect(dialogs[4]).toBe('Session loaded successfully!');
  expect(await page.evaluate(() => ({
    mode: window.__GBDRAW_APP__.mode,
    comparisonHeight: window.__GBDRAW_APP__.adv.comparison_height,
    legend: window.__GBDRAW_APP__.form.legend,
    plotTitlePosition: window.__GBDRAW_APP__.adv.plot_title_position
  }))).toEqual({
    mode: 'linear',
    comparisonHeight: null,
    legend: 'right',
    plotTitlePosition: 'top'
  });
});

test('v40 layout preferences and mode profiles survive fresh-page export and import', async ({ page }) => {
  test.setTimeout(120000);
  const expectedPreferences = {
    circular: {
      single: { legend: 'right', plotTitlePosition: 'top' },
      multi: { legend: 'upper_left', plotTitlePosition: 'bottom' }
    },
    linear: { legend: 'top', plotTitlePosition: 'center' }
  };
  const recordText = `LOCUS       LAYOUT_PREFS                4 bp    DNA     linear   UNK 01-JAN-1980
DEFINITION  Layout preference session fixture.
ACCESSION   LAYOUT_PREFS
FEATURES             Location/Qualifiers
ORIGIN
        1 atgc
//
`;
  const activeProjection = async (mode, multiRecord) => page.evaluate(
    async ({ nextMode, nextMultiRecord }) => {
      const app = window.__GBDRAW_APP__;
      const { state } = await import('./js/state.js');
      app.mode = nextMode;
      if (nextMode === 'circular') app.form.multi_record_canvas = nextMultiRecord;
      await window.Vue.nextTick();
      await window.Vue.nextTick();
      return {
        active: JSON.parse(JSON.stringify(state.activeLayoutPreferences.value)),
        formLegend: app.form.legend,
        plotTitlePosition: app.adv.plot_title_position
      };
    },
    { nextMode: mode, nextMultiRecord: multiRecord }
  );
  const expectActiveProjection = async (mode, multiRecord, expected) => {
    expect(await activeProjection(mode, multiRecord)).toEqual({
      active: expected,
      formLegend: expected.legend,
      plotTitlePosition: expected.plotTitlePosition
    });
  };
  const activeModeProfile = async (mode) => page.evaluate(async (nextMode) => {
    const app = window.__GBDRAW_APP__;
    app.mode = nextMode;
    await window.Vue.nextTick();
    await window.Vue.nextTick();
    return {
      identity: app.adv.identity,
      axisStrokeColor: app.adv.axis_stroke_color
    };
  }, mode);
  const layoutPreferenceTree = () => page.evaluate(async () => {
    const { state } = await import('./js/state.js');
    return JSON.parse(JSON.stringify(state.layoutPreferences));
  });

  await page.goto('/gbdraw/web/index.html', { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);
  await page.evaluate(async ({ genbankText }) => {
    const app = window.__GBDRAW_APP__;
    app.mode = 'circular';
    app.form.multi_record_canvas = false;
    await window.Vue.nextTick();
    app.form.legend = 'right';
    app.adv.plot_title_position = 'top';

    app.form.multi_record_canvas = true;
    await window.Vue.nextTick();
    app.form.legend = 'upper_left';
    app.adv.plot_title_position = 'bottom';
    app.adv.identity = 88;
    app.adv.axis_stroke_color = '#123456';

    app.mode = 'linear';
    await window.Vue.nextTick();
    app.form.legend = 'top';
    app.adv.plot_title_position = 'center';
    app.adv.identity = 77;
    app.adv.axis_stroke_color = '#654321';
    await window.Vue.nextTick();

    app.cInputType = 'gb';
    app.files.c_gb = new File([genbankText], 'layout-preferences.gbk', {
      type: 'text/plain',
      lastModified: 1
    });
    app.sessionTitle = 'layout-preferences-v40';
  }, { genbankText: recordText });

  expect(await layoutPreferenceTree()).toEqual(expectedPreferences);
  await expectActiveProjection('circular', false, expectedPreferences.circular.single);
  await expectActiveProjection('circular', true, expectedPreferences.circular.multi);
  await expectActiveProjection('linear', false, expectedPreferences.linear);

  await activeProjection('circular', false);
  await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    const { state } = await import('./js/state.js');
    state.biologicalFeatures.value = [{ nucleotide_sequence: 'ATGC' }];
    state.linearSeqs[0].gb = new File(
      [await app.files.c_gb.arrayBuffer()],
      'inactive-linear-layout.gbk',
      { type: 'application/genbank', lastModified: 22 }
    );
    state.linearSeqs[0].definition = 'Inactive Linear draft';
    app.adv.linear_track_slots.splice(
      0,
      app.adv.linear_track_slots.length,
      {
        id: 'inactive_spacer',
        renderer: 'spacer',
        enabled: false,
        side: 'below',
        height: '17px',
        spacing: '3px',
        z: 2,
        params: {}
      }
    );
    app.adv.linear_track_slots_enabled = false;
    app.adv.linear_track_slots_axis_index = 1;
  });
  const downloadPromise = page.waitForEvent('download', { timeout: 60000 });
  await page.evaluate(async () => window.__GBDRAW_APP__.saveSessionWithTitle());
  const download = await downloadPromise;
  const savedSessionPath = await download.path();
  expect(savedSessionPath).toBeTruthy();

  const exportedSession = JSON.parse(
    gunzipSync(readFileSync(savedSessionPath)).toString('utf8')
  );
  const legacyLayoutFields = [
    'legend',
    'circularLegendPosition',
    'linearLegendPosition',
    'circularPlotTitlePosition',
    'linearPlotTitlePosition',
    'circularSingleRecordLegendPosition',
    'circularSingleRecordPlotTitlePosition',
    'circularMultiRecordLegendPosition',
    'circularMultiRecordPlotTitlePosition'
  ];
  expect(exportedSession.version).toBe(41);
  expect(exportedSession).not.toHaveProperty('files');
  expect(exportedSession.webFiles).toEqual(expect.any(Object));
  expect(exportedSession.webFiles.bindings.schema).toBe(1);
  expect(exportedSession.webFiles.bindings.c_gb.name).toBe('layout-preferences.gbk');
  expect(exportedSession.webFiles.bindings.linearSeqs[0].gb.name).toBe(
    'inactive-linear-layout.gbk'
  );
  expect(exportedSession.webFiles.bindings.c_gb.resourceId).toBe(
    exportedSession.webFiles.bindings.linearSeqs[0].gb.resourceId
  );
  expect(exportedSession.editorState).toEqual(expect.any(Object));
  expect(exportedSession.editorState.featureCatalog).toBeNull();
  expect(exportedSession.features).not.toHaveProperty('extractedFeatures');
  expect(exportedSession.features).not.toHaveProperty('biologicalFeatures');
  expect(exportedSession.orthogroupState).not.toHaveProperty('groups');
  expect(exportedSession.config.adv.circular_track_slots).toEqual(expect.any(Array));
  expect(exportedSession.config.adv.linear_track_slots).toEqual(expect.any(Array));
  expect(exportedSession.config.adv.linear_track_slots).toEqual([
    expect.objectContaining({
      id: 'inactive_spacer',
      renderer: 'spacer',
      enabled: false,
      height: '17px',
      spacing: '3px'
    })
  ]);
  expect(exportedSession.config.adv).toHaveProperty('circular_track_slots_enabled');
  expect(exportedSession.config.adv).toHaveProperty('linear_track_slots_enabled');
  expect(exportedSession.config.modeProfiles.schema).toBe(1);
  expect(Object.keys(exportedSession.config.modeProfiles.profiles).sort()).toEqual([
    'circular',
    'linear'
  ]);
  expect(exportedSession.config.modeProfiles.profiles.circular).toEqual(
    expect.objectContaining({
      values: expect.objectContaining({
        identity: 88,
        axis_stroke_color: '#123456'
      }),
      managed: expect.objectContaining({
        identity: false,
        axis_stroke_color: false
      })
    })
  );
  expect(exportedSession.config.modeProfiles.profiles.linear).toEqual(
    expect.objectContaining({
      values: expect.objectContaining({
        identity: 77,
        axis_stroke_color: '#654321'
      }),
      managed: expect.objectContaining({
        identity: false,
        axis_stroke_color: false
      })
    })
  );
  expect(exportedSession.ui.layoutPreferences).toEqual(expectedPreferences);
  expect(legacyLayoutFields.filter((field) => (
    Object.prototype.hasOwnProperty.call(exportedSession.ui, field)
  ))).toEqual([]);
  expect(exportedSession.config.form).not.toHaveProperty('legend');
  expect(exportedSession.config.adv).not.toHaveProperty('plot_title_position');

  await page.reload({ waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);
  const dialogPromise = page.waitForEvent('dialog', { timeout: 120000 });
  await page.locator('input[accept^=".json,"]').first().setInputFiles(savedSessionPath);
  const dialog = await dialogPromise;
  expect(dialog.message()).toBe('Session loaded successfully!');
  await dialog.accept();

  expect(await layoutPreferenceTree()).toEqual(expectedPreferences);
  await expectActiveProjection('circular', false, expectedPreferences.circular.single);
  await expectActiveProjection('circular', true, expectedPreferences.circular.multi);
  await expectActiveProjection('linear', false, expectedPreferences.linear);
  expect(await activeModeProfile('circular')).toEqual({
    identity: 88,
    axisStrokeColor: '#123456'
  });
  expect(await activeModeProfile('linear')).toEqual({
    identity: 77,
    axisStrokeColor: '#654321'
  });
  expect(await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    const { state } = await import('./js/state.js');
    return {
      linearFileName: state.linearSeqs[0].gb?.name,
      inactiveSlot: JSON.parse(JSON.stringify(
        app.adv.linear_track_slots.find((slot) => slot.id === 'inactive_spacer')
      )),
      axisIndex: app.adv.linear_track_slots_axis_index,
      enabled: app.adv.linear_track_slots_enabled
    };
  })).toEqual({
    linearFileName: 'inactive-linear-layout.gbk',
    inactiveSlot: expect.objectContaining({
      id: 'inactive_spacer',
      renderer: 'spacer',
      enabled: false,
      height: '17px',
      spacing: '3px'
    }),
    axisIndex: 1,
    enabled: false
  });
});

test('P3 Custom Track drafts survive fresh-page session re-save and Reset history', async ({ page }) => {
  test.setTimeout(180000);
  const recordText = readFileSync(hmmtDnaPath, 'utf8');
  const styleOverride = {
    stroke: '#112233',
    strokeWidth: 2,
    hatch: { pattern: 'diagonal', spacing: 4 },
    label: { color: '#445566', fontSize: 9 }
  };
  const p3Draft = (session) => ({
    circularEnabled: session.config.adv.circular_track_slots_enabled,
    circularAxis: session.config.adv.circular_track_slots_axis_index,
    circularSlots: session.config.adv.circular_track_slots,
    linearEnabled: session.config.adv.linear_track_slots_enabled,
    linearAxis: session.config.adv.linear_track_slots_axis_index,
    linearSlots: session.config.adv.linear_track_slots
  });
  const readSessionDownload = async (download) => {
    const filePath = await download.path();
    expect(filePath).toBeTruthy();
    return JSON.parse(gunzipSync(readFileSync(filePath)).toString('utf8'));
  };
  const loadSession = async (filePath) => {
    const dialogPromise = page.waitForEvent('dialog', { timeout: 120000 });
    await page.locator('input[accept^=".json,"]').first().setInputFiles(filePath);
    const dialog = await dialogPromise;
    expect(dialog.message()).toBe('Session loaded successfully!');
    await dialog.accept();
    await page.waitForFunction(() => {
      const app = window.__GBDRAW_APP__;
      return (
        app.mode === 'circular' &&
        app.adv?.circular_track_slots?.some((slot) => slot.id === 'review_overlay') &&
        app.adv?.linear_track_slots?.some((slot) => slot.id === 'inactive_overlay')
      );
    }, null, { timeout: 120000 });
  };
  const browserDraft = () => page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    return JSON.parse(JSON.stringify({
      circularEnabled: app.adv.circular_track_slots_enabled,
      circularAxis: app.adv.circular_track_slots_axis_index,
      circularSlots: app.adv.circular_track_slots,
      linearEnabled: app.adv.linear_track_slots_enabled,
      linearAxis: app.adv.linear_track_slots_axis_index,
      linearSlots: app.adv.linear_track_slots
    }));
  });

  await page.goto('/gbdraw/web/index.html', { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);
  await page.evaluate(({ genbankText, nestedStyle }) => {
    const app = window.__GBDRAW_APP__;
    app.mode = 'circular';
    app.cInputType = 'gb';
    app.files.c_gb = new File([genbankText], 'p3-session.gbk', {
      type: 'application/genbank',
      lastModified: 31
    });
    app.annotationSets.splice(0, app.annotationSets.length);
    const annotationSet = app.addAnnotationSet('review');
    const annotation = app.addCoordinateAnnotation(annotationSet, { start: 1, end: 20 });
    annotation.label = 'Review';
    annotation.mark = 'line';

    app.adv.circular_track_slots_enabled = true;
    app.adv.circular_track_slots_axis_index = 1;
    app.adv.circular_track_slots.splice(
      0,
      app.adv.circular_track_slots.length,
      {
        id: 'disabled_outer_space',
        renderer: 'spacer',
        enabled: false,
        side: 'outside',
        width: '13px',
        radius: null,
        inner_gap_px: '1px',
        outer_gap_px: '2px',
        z: 0,
        params: {}
      },
      {
        id: 'features',
        renderer: 'features',
        enabled: true,
        side: 'inside',
        width: null,
        radius: null,
        inner_gap_px: null,
        outer_gap_px: null,
        z: 0,
        params: { lane_direction: 'inside' }
      },
      {
        id: 'review_overlay',
        renderer: 'annotations',
        enabled: true,
        side: 'overlay',
        width: null,
        radius: null,
        inner_gap_px: null,
        outer_gap_px: null,
        z: 2,
        params: {
          set_id: 'review',
          marks: ['line', 'band'],
          lane_gap_px: 5,
          padding_px: 0,
          overflow: 'error',
          show_labels: true,
          anchor_slot: 'features',
          layer: 'foreground',
          cover_anchor: true,
          style_override: structuredClone(nestedStyle)
        }
      },
      {
        id: 'inner_space',
        renderer: 'spacer',
        enabled: true,
        side: 'inside',
        width: '9px',
        radius: null,
        inner_gap_px: '2px',
        outer_gap_px: '3px',
        z: 0,
        params: {}
      },
      {
        id: 'disabled_annotation',
        renderer: 'annotations',
        enabled: false,
        side: 'inside',
        width: '18px',
        radius: null,
        inner_gap_px: null,
        outer_gap_px: null,
        z: 1,
        params: {
          set_id: 'review',
          marks: ['bracket', 'highlight'],
          lane_gap_px: 7,
          padding_px: 1,
          overflow: 'clip',
          show_labels: false,
          layer: 'underlay',
          style_override: structuredClone(nestedStyle)
        }
      }
    );

    app.adv.linear_track_slots_enabled = false;
    app.adv.linear_track_slots_axis_index = 2;
    app.adv.linear_track_slots.splice(
      0,
      app.adv.linear_track_slots.length,
      {
        id: 'inactive_above_space',
        renderer: 'spacer',
        enabled: true,
        side: 'above',
        height: '14px',
        spacing: '3px',
        z: 0,
        params: {}
      },
      {
        id: 'inactive_disabled_annotation',
        renderer: 'annotations',
        enabled: false,
        side: 'above',
        height: '17px',
        spacing: '4px',
        z: 1,
        params: {
          set_id: 'review',
          marks: ['highlight'],
          lane_gap_px: 8,
          padding_px: 1,
          overflow: 'compress',
          show_labels: false,
          layer: 'underlay',
          style_override: structuredClone(nestedStyle)
        }
      },
      {
        id: 'linear_features',
        renderer: 'features',
        enabled: true,
        side: 'below',
        height: '',
        spacing: '',
        z: 0,
        params: {}
      },
      {
        id: 'inactive_overlay',
        renderer: 'annotations',
        enabled: true,
        side: 'overlay',
        height: '',
        spacing: '',
        z: 2,
        params: {
          set_id: 'review',
          marks: ['line', 'bracket', 'band', 'highlight'],
          lane_gap_px: 6,
          padding_px: 0,
          overflow: 'error',
          show_labels: true,
          anchor_slot: 'linear_features',
          layer: 'foreground',
          cover_anchor: true,
          style_override: structuredClone(nestedStyle)
        }
      },
      {
        id: 'inactive_disabled_space',
        renderer: 'spacer',
        enabled: false,
        side: 'below',
        height: '11px',
        spacing: '2px',
        z: 0,
        params: {}
      }
    );
    app.sessionTitle = 'p3-drafts-v40';
  }, { genbankText: recordText, nestedStyle: styleOverride });

  await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    if (!app.circularTrackSlotsPanelOpen) app.toggleCircularTrackSlotsPanel();
    await window.Vue.nextTick();
  });
  const circularAnnotationRow = page.locator(
    '[data-capture="circular-track-slot-review_overlay"]'
  );
  await page.locator(
    '[data-capture="circular-track-slot-disabled_outer_space"] input[title="Width"]'
  ).fill('15px');
  await circularAnnotationRow
    .locator('.track-slot-field')
    .filter({ hasText: 'padding' })
    .locator('input[type="number"]')
    .fill('4');
  await circularAnnotationRow
    .locator('label')
    .filter({ hasText: 'Cover anchor' })
    .locator('input[type="checkbox"]')
    .click();
  await circularAnnotationRow
    .getByRole('checkbox', { name: 'band', exact: true })
    .click();

  await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    app.mode = 'linear';
    if (!app.linearTrackSlotsPanelOpen) app.toggleLinearTrackSlotsPanel();
    await window.Vue.nextTick();
    await window.Vue.nextTick();
  });
  const linearAnnotationRow = page.locator(
    '[data-capture="linear-track-slot-inactive_overlay"]'
  );
  await page.locator(
    '[data-capture="linear-track-slot-inactive_above_space"] input[title="Height"]'
  ).fill('16px');
  await linearAnnotationRow
    .locator('.track-slot-field')
    .filter({ hasText: 'lane gap' })
    .locator('input[type="number"]')
    .fill('7');
  await linearAnnotationRow
    .locator('.track-slot-field')
    .filter({ hasText: 'padding' })
    .locator('input[type="number"]')
    .fill('3');
  await linearAnnotationRow
    .locator('label')
    .filter({ hasText: 'Cover anchor' })
    .locator('input[type="checkbox"]')
    .click();
  await linearAnnotationRow
    .getByRole('checkbox', { name: 'highlight', exact: true })
    .click();
  await page.evaluate(async () => {
    window.__GBDRAW_APP__.mode = 'circular';
    await window.Vue.nextTick();
    await window.Vue.nextTick();
  });
  expect(await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    const circularSpacer = app.adv.circular_track_slots.find(
      (slot) => slot.id === 'disabled_outer_space'
    );
    const circularAnnotation = app.adv.circular_track_slots.find(
      (slot) => slot.id === 'review_overlay'
    );
    const linearSpacer = app.adv.linear_track_slots.find(
      (slot) => slot.id === 'inactive_above_space'
    );
    const linearAnnotation = app.adv.linear_track_slots.find(
      (slot) => slot.id === 'inactive_overlay'
    );
    return {
      circularSpacerWidth: circularSpacer?.width,
      circularAnnotation: {
        marks: circularAnnotation?.params?.marks,
        padding: circularAnnotation?.params?.padding_px,
        coverAnchor: circularAnnotation?.params?.cover_anchor ?? false
      },
      linearSpacerHeight: linearSpacer?.height,
      linearAnnotation: {
        marks: linearAnnotation?.params?.marks,
        laneGap: linearAnnotation?.params?.lane_gap_px,
        padding: linearAnnotation?.params?.padding_px,
        coverAnchor: linearAnnotation?.params?.cover_anchor ?? false
      }
    };
  })).toEqual({
    circularSpacerWidth: '15px',
    circularAnnotation: {
      marks: ['line'],
      padding: 4,
      coverAnchor: false
    },
    linearSpacerHeight: '16px',
    linearAnnotation: {
      marks: ['line', 'bracket', 'band'],
      laneGap: 7,
      padding: 3,
      coverAnchor: false
    }
  });

  const initialDownloadPromise = page.waitForEvent('download', { timeout: 60000 });
  await page.evaluate(async () => window.__GBDRAW_APP__.saveSessionWithTitle());
  const initialDownload = await initialDownloadPromise;
  const initialPath = await initialDownload.path();
  expect(initialPath).toBeTruthy();
  const initialSession = JSON.parse(
    gunzipSync(readFileSync(initialPath)).toString('utf8')
  );
  expect(initialSession.version).toBe(41);
  const expectedDraft = p3Draft(initialSession);
  expect(expectedDraft.circularEnabled).toBe(true);
  expect(expectedDraft.linearEnabled).toBe(false);
  expect(expectedDraft.circularSlots.map((slot) => [slot.id, slot.enabled])).toEqual([
    ['disabled_outer_space', false],
    ['features', true],
    ['review_overlay', true],
    ['inner_space', true],
    ['disabled_annotation', false]
  ]);
  expect(expectedDraft.linearSlots.map((slot) => [slot.id, slot.enabled])).toEqual([
    ['inactive_above_space', true],
    ['inactive_disabled_annotation', false],
    ['linear_features', true],
    ['inactive_overlay', true],
    ['inactive_disabled_space', false]
  ]);
  expect(
    expectedDraft.circularSlots.find((slot) => slot.id === 'review_overlay')
      .params.style_override
  ).toEqual(styleOverride);
  expect(
    expectedDraft.linearSlots.find((slot) => slot.id === 'inactive_overlay')
      .params.style_override
  ).toEqual(styleOverride);

  await page.reload({ waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);
  await loadSession(initialPath);
  expect(await browserDraft()).toEqual(expectedDraft);

  await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    if (!app.circularTrackSlotsPanelOpen) app.toggleCircularTrackSlotsPanel();
    await window.Vue.nextTick();
  });
  const annotationRow = page.locator(
    '[data-capture="circular-track-slot-review_overlay"]'
  );
  const laneGapInput = annotationRow
    .locator('.track-slot-field')
    .filter({ hasText: 'lane gap' })
    .locator('input[type="number"]');
  await laneGapInput.fill('9');
  await laneGapInput.evaluate((input) => input.blur());
  await page.waitForFunction(() => (
    window.__GBDRAW_HISTORY__?.canUndo() &&
    window.__GBDRAW_APP__.adv.circular_track_slots
      .find((slot) => slot.id === 'review_overlay')?.params.lane_gap_px === 9
  ));
  expect(
    (await browserDraft()).circularSlots
      .find((slot) => slot.id === 'review_overlay').params.style_override
  ).toEqual(styleOverride);

  await page.evaluate(async () => window.__GBDRAW_APP__.undoHistory());
  await page.waitForFunction(() => (
    window.__GBDRAW_HISTORY__?.canRedo() &&
    window.__GBDRAW_APP__.adv.circular_track_slots
      .find((slot) => slot.id === 'review_overlay')?.params.lane_gap_px === 5
  ));
  expect(await browserDraft()).toEqual(expectedDraft);
  await page.evaluate(async () => window.__GBDRAW_APP__.redoHistory());
  await page.waitForFunction(() => (
    window.__GBDRAW_APP__.adv.circular_track_slots
      .find((slot) => slot.id === 'review_overlay')
      ?.params.lane_gap_px === 9
  ));
  expect(
    (await browserDraft()).circularSlots
      .find((slot) => slot.id === 'review_overlay').params.style_override
  ).toEqual(styleOverride);
  await page.evaluate(async () => window.__GBDRAW_APP__.undoHistory());
  await page.waitForFunction(() => (
    window.__GBDRAW_APP__.adv.circular_track_slots
      .find((slot) => slot.id === 'review_overlay')
      ?.params.lane_gap_px === 5
  ));

  await page.evaluate(() => {
    window.__GBDRAW_APP__.sessionTitle = 'p3-drafts-resaved';
  });
  const secondDownloadPromise = page.waitForEvent('download', { timeout: 60000 });
  await page.evaluate(async () => window.__GBDRAW_APP__.saveSessionWithTitle());
  const secondDownload = await secondDownloadPromise;
  const secondPath = await secondDownload.path();
  expect(secondPath).toBeTruthy();
  const secondSession = JSON.parse(
    gunzipSync(readFileSync(secondPath)).toString('utf8')
  );
  expect(p3Draft(secondSession)).toEqual(expectedDraft);

  await page.reload({ waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);
  await loadSession(secondPath);
  expect(await browserDraft()).toEqual(expectedDraft);
  await page.evaluate(() => {
    window.__GBDRAW_APP__.sessionTitle = 'p3-drafts-fresh-resave';
  });
  const thirdDownloadPromise = page.waitForEvent('download', { timeout: 60000 });
  await page.evaluate(async () => window.__GBDRAW_APP__.saveSessionWithTitle());
  const thirdSession = await readSessionDownload(await thirdDownloadPromise);
  expect(p3Draft(thirdSession)).toEqual(expectedDraft);

  await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    if (!app.circularTrackSlotsPanelOpen) app.toggleCircularTrackSlotsPanel();
    await window.Vue.nextTick();
  });
  await page.locator(
    'button[title="Replace the current custom stack from the simple controls; this is the only action that regenerates it"]:visible'
  ).click();
  await page.waitForFunction(() => (
    !window.__GBDRAW_APP__.adv.circular_track_slots
      .some((slot) => slot.id === 'review_overlay')
  ));
  const resetDraft = await browserDraft();
  expect(resetDraft.linearSlots).toEqual(expectedDraft.linearSlots);
  expect(resetDraft.linearAxis).toBe(expectedDraft.linearAxis);
  expect(resetDraft.linearEnabled).toBe(expectedDraft.linearEnabled);

  await page.evaluate(async () => window.__GBDRAW_APP__.undoHistory());
  await page.waitForFunction(() => (
    window.__GBDRAW_APP__.adv.circular_track_slots
      .some((slot) => slot.id === 'review_overlay')
  ));
  expect(await browserDraft()).toEqual(expectedDraft);
  await page.evaluate(async () => window.__GBDRAW_APP__.redoHistory());
  await page.waitForFunction(() => (
    !window.__GBDRAW_APP__.adv.circular_track_slots
      .some((slot) => slot.id === 'review_overlay')
  ));
  expect(await browserDraft()).toEqual(resetDraft);
});

test('Session commit failure restores the pre-import state', async ({ page }) => {
  await page.goto('/gbdraw/web/index.html', { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);

  const recordText = `LOCUS       ROLLBACK                    4 bp    DNA     linear   UNK 01-JAN-1980
DEFINITION  Session rollback fixture.
ACCESSION   ROLLBACK
FEATURES             Location/Qualifiers
ORIGIN
        1 atgc
//
`;
  const session = {
    format: 'gbdraw-session',
    version: 33,
    title: 'replacement',
    renderRequest: {
      schema: 2,
      mode: 'linear',
      records: [{
        recordKey: 'record-1',
        source: { kind: 'genbank', resourceId: 'record-1-genbank' },
        selector: null,
        region: null,
        presentation: {
          label: null, subtitle: null, reverseComplement: false,
          gridRow: null, gridColumn: null
        }
      }],
      diagramOptions: { configOverrides: {}, tracks: {} },
      layout: {},
      comparisons: [],
      output: {
        prefix: 'rollback', formats: ['interactive_svg'], overwrite: true,
        interactiveMetadataPolicy: 'auto'
      }
    },
    resources: {
      'record-1-genbank': {
        kind: 'genbank', name: 'rollback.gb', type: 'text/plain',
        size: Buffer.byteLength(recordText), lastModified: 0,
        encoding: 'base64', data: Buffer.from(recordText).toString('base64')
      }
    }
  };

  await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    const { state } = await import('./js/state.js');
    app.mode = 'circular';
    app.sessionTitle = 'keep-after-rollback';
    app.form.legend = 'left';
    app.adv.comparison_height = 37;
    const original = state.normalizePaletteColors;
    let injected = false;
    state.normalizePaletteColors = (...args) => {
      if (!injected) {
        injected = true;
        throw new Error('Injected session commit failure');
      }
      return original(...args);
    };
  });
  const dialogs = [];
  page.on('dialog', async (dialog) => {
    dialogs.push(dialog.message());
    await dialog.accept();
  });
  await page.locator('input[accept^=".json,"]').first().setInputFiles({
    name: 'rollback.gbdraw-session.json',
    mimeType: 'application/json',
    buffer: Buffer.from(JSON.stringify(session))
  });
  await expect.poll(() => dialogs.length).toBe(1);
  expect(dialogs[0]).toContain('Injected session commit failure');
  expect(await page.evaluate(() => ({
    mode: window.__GBDRAW_APP__.mode,
    title: window.__GBDRAW_APP__.sessionTitle,
    legend: window.__GBDRAW_APP__.form.legend,
    comparisonHeight: window.__GBDRAW_APP__.adv.comparison_height
  }))).toEqual({
    mode: 'circular',
    title: 'keep-after-rollback',
    legend: 'left',
    comparisonHeight: 37
  });
});

test('HmmtDNA middle overlap layout keeps feature, GC, and skew bands disjoint', async ({ page }) => {
  test.setTimeout(240000);
  const genbank = readFileSync(hmmtDnaPath, 'utf8');

  await page.goto('/gbdraw/web/index.html', { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);
  await page.evaluate((genbankText) => {
    const app = window.__GBDRAW_APP__;
    app.mode = 'linear';
    app.lInputType = 'gb';
    app.setLinearSeqPrimaryFile(0, 'gb', new File([genbankText], 'HmmtDNA.gbk', {
      type: 'text/plain', lastModified: 1
    }));
    Object.assign(app.form, {
      linear_track_layout: 'middle',
      separate_strands: true,
      show_gc: true,
      show_skew: true,
      show_depth: false,
      show_labels_linear: 'none',
      legend: 'none'
    });
    app.adv.resolve_overlaps = true;
    app.setLinearTrackSlotsEnabled(false);
  }, genbank);

  const generated = await runDiagramWithDiagnostics(page);
  expect(generated).toEqual({
    result: { status: 'ok' },
    errorSummary: '',
    errorDetails: []
  });
  await page.waitForFunction(() => {
    const feature = document.querySelector(
      '[data-gbdraw-feature-id][data-gbdraw-record-id="NC_012920.1"]'
    );
    const svg = feature?.ownerSVGElement;
    return Boolean(svg?.getElementById('gc_content') && svg?.getElementById('gc_skew'));
  });

  const geometry = await page.evaluate(() => {
    const featureSelector = '[data-gbdraw-feature-id][data-gbdraw-record-id="NC_012920.1"]';
    const feature = document.querySelector(featureSelector);
    const svg = feature?.ownerSVGElement;
    if (!svg) throw new Error('Live HmmtDNA SVG was not found.');
    const features = Array.from(svg.querySelectorAll(featureSelector));
    const gc = svg.getElementById('gc_content');
    const skew = svg.getElementById('gc_skew');
    if (features.length === 0 || !gc || !skew) {
      throw new Error('Expected live feature, GC content, and GC skew geometry.');
    }

    const clientBand = (elements) => {
      const rects = elements.map((element) => element.getBoundingClientRect());
      return {
        top: Math.min(...rects.map((rect) => rect.top)),
        bottom: Math.max(...rects.map((rect) => rect.bottom))
      };
    };
    const bboxBand = (elements) => {
      let top = Number.POSITIVE_INFINITY;
      let bottom = Number.NEGATIVE_INFINITY;
      elements.forEach((element) => {
        const box = element.getBBox();
        const matrix = element.getCTM();
        if (!matrix) throw new Error(`Missing SVG transform for ${element.id || element.tagName}.`);
        [
          [box.x, box.y],
          [box.x + box.width, box.y],
          [box.x, box.y + box.height],
          [box.x + box.width, box.y + box.height]
        ].forEach(([x, y]) => {
          const point = new DOMPoint(x, y).matrixTransform(matrix);
          top = Math.min(top, point.y);
          bottom = Math.max(bottom, point.y);
        });
      });
      return { top, bottom };
    };
    const args = Array.isArray(window.__GBDRAW_APP__.lastRunInfo?.invocation?.args)
      ? window.__GBDRAW_APP__.lastRunInfo.invocation.args
      : [];
    return {
      client: {
        features: clientBand(features),
        gc: clientBand([gc]),
        skew: clientBand([skew])
      },
      bbox: {
        features: bboxBand(features),
        gc: bboxBand([gc]),
        skew: bboxBand([skew])
      },
      args,
      settings: {
        separateStrands: window.__GBDRAW_APP__.form.separate_strands,
        resolveOverlaps: window.__GBDRAW_APP__.adv.resolve_overlaps,
        showGc: window.__GBDRAW_APP__.form.show_gc,
        showSkew: window.__GBDRAW_APP__.form.show_skew
      },
      reproducibilityLevel: window.__GBDRAW_APP__.lastRunInfo?.reproducibility?.level
    };
  });

  expect(geometry.settings).toEqual({
    separateStrands: true,
    resolveOverlaps: true,
    showGc: true,
    showSkew: true
  });
  expect(geometry.args).toEqual(['--session', 'out.gbdraw-session.json']);
  expect(geometry.reproducibilityLevel).toBe('canonical-request');
  for (const bands of [geometry.client, geometry.bbox]) {
    for (const band of Object.values(bands)) {
      expect(Number.isFinite(band.top)).toBe(true);
      expect(Number.isFinite(band.bottom)).toBe(true);
      expect(band.bottom).toBeGreaterThan(band.top);
    }
    expect(bands.features.bottom).toBeLessThanOrEqual(bands.gc.top + 0.5);
    expect(bands.gc.bottom).toBeLessThanOrEqual(bands.skew.top + 0.5);
  }
});

test('Linear sparse diagonal depth generates and survives a session round trip', async ({ page }) => {
  test.setTimeout(240000);
  const genbankA = readFileSync(sparseGenbankAPath, 'utf8');
  const genbankB = readFileSync(sparseGenbankBPath, 'utf8');
  const makeDepthTsv = (recordId, length, depth) => [
    'reference_name\tposition\tdepth',
    ...Array.from(
      { length: Math.ceil(length / 1000) },
      (_, index) => `${recordId}\t${Math.min(length, index * 1000 + 1)}\t${depth + (index % 3)}`
    ),
    `${recordId}\t${length}\t${depth}`
  ].join('\n');
  const depthA = makeDepthTsv('BGC0000711', 30837, 10);
  const depthB = makeDepthTsv('BGC0000713', 31892, 50);

  await page.goto('/gbdraw/web/index.html', { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);

  await page.evaluate(({ genbankAText, genbankBText, depthAText, depthBText }) => {
    const app = window.__GBDRAW_APP__;
    app.mode = 'linear';
    app.lInputType = 'gb';
    app.addLinearSeq();
    app.setLinearSeqPrimaryFile(0, 'gb', new File([genbankAText], 'BGC0000711.gbk', {
      type: 'text/plain', lastModified: 1
    }));
    app.setLinearSeqPrimaryFile(1, 'gb', new File([genbankBText], 'BGC0000713.gbk', {
      type: 'text/plain', lastModified: 2
    }));
    app.setLinearDepthFile(app.linearSeqs[0], 0, new File([depthAText], 'sample-a.tsv', {
      type: 'text/tab-separated-values', lastModified: 3
    }));
    app.addLinearDepthTrack();
    app.setLinearDepthFile(app.linearSeqs[1], 1, new File([depthBText], 'sample-b.tsv', {
      type: 'text/tab-separated-values', lastModified: 4
    }));
    Object.assign(app.adv.depth_tracks[0], {
      label: 'Sample A', color: '#112233', height: 12
    });
    Object.assign(app.adv.depth_tracks[1], {
      label: 'Sample B', color: '#445566', height: 18
    });
    app.form.show_gc = false;
    app.form.show_skew = false;
    app.form.legend = 'none';
    app.setLinearTrackSlotsEnabled(true);
    app.adv.linear_track_slots.splice(
      0,
      app.adv.linear_track_slots.length,
      {
        id: 'depth_a', renderer: 'depth', enabled: true, side: 'above', height: '12px',
        params: { track_index: 0, legend_label: 'Sample A' }
      },
      {
        id: 'depth_b', renderer: 'depth', enabled: true, side: 'above', height: '18px',
        params: { track_index: 1, legend_label: 'Sample B' }
      },
      { id: 'features', renderer: 'features', enabled: true, side: 'below', params: {} }
    );
    app.adv.linear_track_slots_axis_index = 2;
    app.sessionTitle = 'sparse-depth-e2e';
  }, {
    genbankAText: genbankA,
    genbankBText: genbankB,
    depthAText: depthA,
    depthBText: depthB
  });

  const firstRun = await runDiagramWithDiagnostics(page);
  expect(firstRun).toEqual({
    result: { status: 'ok' },
    errorSummary: '',
    errorDetails: []
  });
  const firstResult = await inspectSparseDepthResult(page);
  expect(firstResult.resultCount).toBe(1);
  expect(firstResult.groups).toEqual({
    depthARecord1: true,
    depthARecord1Axis: true,
    depthARecord2: false,
    depthBRecord1: false,
    depthBRecord2: true,
    depthBRecord2Axis: true
  });
  expect(firstResult.depthAFills).toContain('#112233');
  expect(firstResult.depthBFills).toContain('#445566');
  expect(firstResult.depthArgs).toEqual([]);

  const downloadPromise = page.waitForEvent('download', { timeout: 60000 });
  await page.evaluate(async () => window.__GBDRAW_APP__.saveSessionWithTitle());
  const download = await downloadPromise;
  const savedSessionPath = await download.path();
  expect(savedSessionPath).toBeTruthy();

  await page.reload({ waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);
  const dialogPromise = page.waitForEvent('dialog', { timeout: 120000 });
  await page.locator('input[accept^=".json,"]').first().setInputFiles(savedSessionPath);
  const dialog = await dialogPromise;
  expect(dialog.message()).toBe('Session loaded successfully!');
  await dialog.accept();
  await page.waitForFunction(() => {
    const app = window.__GBDRAW_APP__;
    return app.mode === 'linear' &&
      app.linearSeqs?.length === 2 &&
      app.adv?.linear_track_slots_axis_index === 2;
  }, null, { timeout: 120000 });

  const restoredState = await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    return {
      sparseRows: app.linearSeqs.map((seq) => seq.depth.map(Boolean)),
      labels: app.adv.depth_tracks.map((track) => track.label),
      colors: app.adv.depth_tracks.map((track) => track.color.toLowerCase()),
      slotIds: app.adv.linear_track_slots.map((slot) => slot.id),
      slotIndexes: app.adv.linear_track_slots
        .filter((slot) => slot.renderer === 'depth')
        .map((slot) => slot.params.track_index),
      axisIndex: app.adv.linear_track_slots_axis_index
    };
  });
  expect(restoredState).toEqual({
    sparseRows: [[true, false], [false, true]],
    labels: ['Sample A', 'Sample B'],
    colors: ['#112233', '#445566'],
    slotIds: ['depth_a', 'depth_b', 'features'],
    slotIndexes: [0, 1],
    axisIndex: 2
  });

  const secondRun = await runDiagramWithDiagnostics(page);
  expect(secondRun).toEqual({
    result: { status: 'ok' },
    errorSummary: '',
    errorDetails: []
  });
  const secondResult = await inspectSparseDepthResult(page);
  expect(secondResult.groups).toEqual(firstResult.groups);
  expect(secondResult.depthAFills).toContain('#112233');
  expect(secondResult.depthBFills).toContain('#445566');
  expect(secondResult.depthArgs).toEqual([]);
});

test('Circular sparse diagonal depth survives a session round trip and track removal', async ({ page }) => {
  test.setTimeout(240000);
  const genbankA = readFileSync(sparseGenbankAPath, 'utf8');
  const genbankB = readFileSync(sparseGenbankBPath, 'utf8');
  const makeDepthTsv = (recordId, length, depth) => [
    'reference_name\tposition\tdepth',
    ...Array.from(
      { length: Math.ceil(length / 1000) },
      (_, index) => `${recordId}\t${Math.min(length, index * 1000 + 1)}\t${depth + (index % 3)}`
    ),
    `${recordId}\t${length}\t${depth}`
  ].join('\n');
  const depthA = makeDepthTsv('BGC0000708', 40579, 10);
  const depthB = makeDepthTsv('BGC0000709', 50466, 50);

  await page.goto('/gbdraw/web/index.html', { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);
  await page.evaluate(({ genbankAText, genbankBText, depthAText, depthBText }) => {
    const app = window.__GBDRAW_APP__;
    app.mode = 'circular';
    app.cInputType = 'gb';
    const combinedFile = new File([`${genbankAText}\n${genbankBText}`], 'combined.gbk', {
      type: 'text/plain', lastModified: 1
    });
    const originalArrayBuffer = combinedFile.arrayBuffer.bind(combinedFile);
    let primaryReadCount = 0;
    Object.defineProperty(combinedFile, 'arrayBuffer', {
      value: () => {
        primaryReadCount += 1;
        return originalArrayBuffer();
      }
    });
    window.__GBDRAW_PRIMARY_READ_COUNT__ = () => primaryReadCount;
    app.files.c_gb = combinedFile;
    app.files.c_depth = [
      [new File([depthAText], 'sample-a.tsv', {
        type: 'text/tab-separated-values', lastModified: 2
      }), null],
      [null, new File([depthBText], 'sample-b.tsv', {
        type: 'text/tab-separated-values', lastModified: 3
      })]
    ];
    Object.assign(app.form, {
      multi_record_canvas: true,
      show_depth: true,
      suppress_gc: true,
      suppress_skew: true,
      labels_mode: 'none',
      legend: 'none'
    });
    app.adv.depth_tracks.splice(
      0,
      app.adv.depth_tracks.length,
      { label: 'Sample A', color: '#112233', large_tick_interval: null, small_tick_interval: null, tick_font_size: null },
      { label: 'Sample B', color: '#445566', large_tick_interval: null, small_tick_interval: null, tick_font_size: null }
    );
    app.resetCircularTrackSlotsFromSimpleControls();
    app.setCircularTrackSlotsEnabled(true);
    app.sessionTitle = 'circular-sparse-depth-e2e';
  }, {
    genbankAText: genbankA,
    genbankBText: genbankB,
    depthAText: depthA,
    depthBText: depthB
  });

  const firstRun = await runDiagramWithDiagnostics(page);
  expect(firstRun).toEqual({
    result: { status: 'ok' },
    errorSummary: '',
    errorDetails: []
  });
  await page.waitForFunction(() => window.__GBDRAW_APP__?.circularRecordList?.length === 2, null, {
    timeout: 120000
  });
  expect(await page.evaluate(() => window.__GBDRAW_PRIMARY_READ_COUNT__())).toBe(1);
  const firstResult = await inspectCircularSparseDepthResult(page);
  expect(firstResult.groups).toEqual({
    depthARecord1: true,
    depthARecord1Axis: true,
    depthARecord2: false,
    depthBRecord1: false,
    depthBRecord2: true,
    depthBRecord2Axis: true
  });
  expect(firstResult.depthAFills).toContain('#112233');
  expect(firstResult.depthBFills).toContain('#445566');
  expect(firstResult.depthArgs).toEqual([]);

  const downloadPromise = page.waitForEvent('download', { timeout: 60000 });
  await page.evaluate(async () => window.__GBDRAW_APP__.saveSessionWithTitle());
  const download = await downloadPromise;
  const savedSessionPath = await download.path();
  expect(savedSessionPath).toBeTruthy();

  await page.reload({ waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);
  const dialogPromise = page.waitForEvent('dialog', { timeout: 120000 });
  await page.locator('input[accept^=".json,"]').first().setInputFiles(savedSessionPath);
  const dialog = await dialogPromise;
  expect(dialog.message()).toBe('Session loaded successfully!');
  await dialog.accept();
  await page.waitForFunction(() => {
    const app = window.__GBDRAW_APP__;
    return app.mode === 'circular' &&
      app.circularRecordList?.length === 2 &&
      Array.isArray(app.files?.c_depth) &&
      app.files.c_depth.length === 2;
  }, null, { timeout: 120000 });

  const restoredState = await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    return {
      sparseRows: app.files.c_depth.map((row) => row.map(Boolean)),
      labels: app.adv.depth_tracks.map((track) => track.label),
      colors: app.adv.depth_tracks.map((track) => track.color.toLowerCase())
    };
  });
  expect(restoredState).toEqual({
    sparseRows: [[true, false], [false, true]],
    labels: ['Sample A', 'Sample B'],
    colors: ['#112233', '#445566']
  });

  const secondRun = await runDiagramWithDiagnostics(page);
  expect(secondRun).toEqual({
    result: { status: 'ok' },
    errorSummary: '',
    errorDetails: []
  });
  const secondResult = await inspectCircularSparseDepthResult(page);
  expect(secondResult.groups).toEqual(firstResult.groups);
  expect(secondResult.depthAFills).toContain('#112233');
  expect(secondResult.depthBFills).toContain('#445566');
  expect(secondResult.depthArgs).toEqual([]);

  await page.evaluate(() => window.__GBDRAW_APP__.removeCircularDepthTrack(1));
  await page.waitForFunction(() => {
    const app = window.__GBDRAW_APP__;
    const remainingFiles = app.files.c_depth.flat().filter(Boolean);
    return remainingFiles.length === 1 &&
      app.adv.circular_track_slots
        .filter((slot) => slot?.renderer === 'depth' && slot.enabled !== false)
        .every((slot) => Number(slot.params?.track_index) < 1);
  });

  const afterRemoval = await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    return {
      names: app.files.c_depth.flat().filter(Boolean).map((file) => file.name),
      labels: app.adv.depth_tracks.map((track) => track?.label || ''),
      slots: app.adv.circular_track_slots
        .filter((slot) => slot?.renderer === 'depth')
        .map((slot) => ({
          id: slot.id,
          enabled: slot.enabled !== false,
          trackIndex: slot.params?.track_index,
          legendLabel: slot.params?.legend_label
        }))
    };
  });

  expect(afterRemoval.names).toEqual(['sample-a.tsv']);
  expect(afterRemoval.labels).toEqual(['Sample A']);
  expect(afterRemoval.slots.filter((slot) => slot.enabled).map((slot) => slot.trackIndex)).toEqual([0]);
  expect(afterRemoval.slots.some((slot) => Number(slot.trackIndex) >= 1)).toBe(false);
  expect(afterRemoval.slots.filter((slot) => slot.enabled).map((slot) => slot.legendLabel)).toEqual([
    'Sample A'
  ]);
});

test('BGC session keeps restored feature metadata selectable in the preview', async ({ page }) => {
  page.on('dialog', (dialog) => dialog.accept());
  await page.goto('/gbdraw/web/index.html', { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);
  await page.evaluate(() => {
    window.__GBDRAW_APP__.pyodideReady = true;
  });

  await page.locator('input[accept^=".json,"]').first().setInputFiles(bgcSessionPath);
  await page.waitForFunction(() => window.__GBDRAW_APP__?.results?.length > 0);
  await page.waitForTimeout(250);

  const summary = await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    const svg =
      document.querySelector('[data-gbdraw-feature-id]')?.ownerSVGElement ||
      document.querySelector('svg');
    const renderedIds = Array.from(
      svg?.querySelectorAll('[data-gbdraw-feature-id], path[id^="f"], polygon[id^="f"], rect[id^="f"]') || []
    )
      .map((el) => el.getAttribute('data-gbdraw-feature-id') || el.id || '')
      .filter(Boolean);
    const uniqueRenderedIds = Array.from(new Set(renderedIds));
    const featureIds = new Set(
      (Array.isArray(app.extractedFeatures) ? app.extractedFeatures : [])
        .map((feature) => String(feature?.svg_id || '').trim())
        .filter(Boolean)
    );
    return {
      extractedCount: Array.isArray(app.extractedFeatures) ? app.extractedFeatures.length : 0,
      uniqueRenderedCount: uniqueRenderedIds.length,
      matchedCount: uniqueRenderedIds.filter((id) => featureIds.has(id)).length
    };
  });

  expect(summary.extractedCount).toBeGreaterThan(0);
  expect(summary.uniqueRenderedCount).toBeGreaterThan(0);
  expect(summary.matchedCount).toBe(summary.uniqueRenderedCount);

  const target = await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    const svg =
      document.querySelector('[data-gbdraw-feature-id]')?.ownerSVGElement ||
      document.querySelector('svg');
    const features = Array.isArray(app.extractedFeatures) ? app.extractedFeatures : [];
    const targetFeature = features.find((feature) => String(feature?.orthogroupId || '').trim() === 'og_1');
    const targetFeatureId = String(targetFeature?.svg_id || '').trim();
    const element = Array.from(
      svg.querySelectorAll('[data-gbdraw-feature-id], path[id^="f"], polygon[id^="f"], rect[id^="f"]')
    ).find((candidate) => {
      const candidateId = candidate.getAttribute('data-gbdraw-feature-id') || candidate.id || '';
      return targetFeatureId && candidateId === targetFeatureId;
    });
    if (!element) throw new Error('Rendered og_1 feature was not found');
    element.scrollIntoView({ block: 'center', inline: 'center' });
    await new Promise((resolve) => requestAnimationFrame(() => requestAnimationFrame(resolve)));
    const rect = element.getBoundingClientRect();
    return {
      id: element.getAttribute('data-gbdraw-feature-id') || element.id || '',
      x: rect.left + rect.width / 2,
      y: rect.top + rect.height / 2
    };
  });

  await page.mouse.click(target.x, target.y);
  await page.waitForFunction(
    (featureId) => {
      const feature = window.__GBDRAW_APP__?.clickedFeature;
      return feature?.svg_id === featureId && feature?.orthogroupId === 'og_1';
    },
    target.id
  );
  await expect(page.locator('.feature-popup')).toBeVisible();
  await expect(page.locator('.feature-popup').getByRole('button', { name: /Align/ })).toBeVisible();
});

test('BGC session selected feature Hide undo redo keeps visibility and legend stable', async ({ page }) => {
  page.on('dialog', (dialog) => dialog.accept());
  await page.goto('/gbdraw/web/index.html', { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);
  await page.evaluate(() => {
    window.__GBDRAW_APP__.pyodideReady = true;
  });

  await page.locator('input[accept^=".json,"]').first().setInputFiles(bgcSessionPath);
  await page.waitForFunction(() => window.__GBDRAW_APP__?.results?.length > 0);
  await page.waitForTimeout(250);

  const featureIds = await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    const svg = document.querySelector('[data-gbdraw-feature-id]')?.ownerSVGElement || document.querySelector('svg');
    const renderedIds = new Set(
      Array.from(svg?.querySelectorAll('[data-gbdraw-feature-id], path[id^="f"], polygon[id^="f"], rect[id^="f"]') || [])
        .map((el) => String(el.getAttribute('data-gbdraw-feature-id') || el.id || '').replace(/__part\d+$/, ''))
        .filter(Boolean)
    );
    return (Array.isArray(app.extractedFeatures) ? app.extractedFeatures : [])
      .map((feature) => String(feature?.svg_id || '').trim())
      .filter((id) => renderedIds.has(id))
      .slice(0, 2);
  });
  expect(featureIds.length).toBeGreaterThanOrEqual(2);

  const legendTransformBefore = await page.evaluate(() => {
    const svg = document.querySelector('svg');
    const legend = svg?.querySelector('#legend, #feature_legend, [id*="legend" i]');
    return legend ? legend.getAttribute('transform') || '' : null;
  });
  expect(legendTransformBefore).not.toBeNull();

  const getStates = async () => page.evaluate((ids) => {
    const collect = (root) => ids.map((id) => {
      const elements = Array.from(
        root?.querySelectorAll?.('[data-gbdraw-feature-id], path[id^="f"], polygon[id^="f"], rect[id^="f"]') || []
      ).filter((el) => String(el.getAttribute('data-gbdraw-feature-id') || el.id || '').replace(/__part\d+$/, '') === id);
      return {
        id,
        count: elements.length,
        hidden: elements.length > 0 && elements.every((el) => el.getAttribute('display') === 'none')
      };
    });

    const app = window.__GBDRAW_APP__;
    const liveSvg = document.querySelector('svg');
    const content = app.results?.[app.selectedResultIndex]?.content || '';
    const parsedSvg = new DOMParser().parseFromString(content, 'image/svg+xml').querySelector('svg');
    return {
      live: collect(liveSvg),
      serialized: collect(parsedSvg),
      serializedContent: content
    };
  }, featureIds);

  await page.evaluate(async (ids) => {
    const app = window.__GBDRAW_APP__;
    app.selectedFeatureBulkVisibility = 'off';
    app.selectedFeatureIds = new Set(ids);
    app.selectedFeatureAnchorId = ids[0];
    await app.applySelectedFeatureVisibility();
    await new Promise((resolve) => requestAnimationFrame(() => requestAnimationFrame(resolve)));
  }, featureIds);
  let states = await getStates();
  expect(states.live.every((state) => state.count > 0 && state.hidden)).toBe(true);
  expect(states.serialized.every((state) => state.count > 0 && state.hidden)).toBe(true);
  expect(states.serializedContent).not.toContain('gbdraw-feature-selected');

  await page.evaluate(async () => {
    await window.__GBDRAW_APP__.undoHistory();
    await new Promise((resolve) => requestAnimationFrame(() => requestAnimationFrame(resolve)));
  });
  states = await getStates();
  expect(states.live.every((state) => state.count > 0 && !state.hidden)).toBe(true);
  expect(states.serialized.every((state) => state.count > 0 && !state.hidden)).toBe(true);

  const legendTransformAfterUndo = await page.evaluate(() => {
    const svg = document.querySelector('svg');
    const legend = svg?.querySelector('#legend, #feature_legend, [id*="legend" i]');
    return legend ? legend.getAttribute('transform') || '' : null;
  });
  expect(legendTransformAfterUndo).toBe(legendTransformBefore);

  await page.evaluate(async () => {
    await window.__GBDRAW_APP__.redoHistory();
    await new Promise((resolve) => requestAnimationFrame(() => requestAnimationFrame(resolve)));
  });
  states = await getStates();
  expect(states.live.every((state) => state.count > 0 && state.hidden)).toBe(true);
  expect(states.serialized.every((state) => state.count > 0 && state.hidden)).toBe(true);
  expect(states.serializedContent).not.toContain('gbdraw-feature-selected');
});
