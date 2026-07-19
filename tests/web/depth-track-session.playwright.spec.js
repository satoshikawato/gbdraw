const { test, expect } = require('@playwright/test');
const { createReadStream, existsSync, readFileSync } = require('node:fs');
const { createServer } = require('node:http');
const { extname, join, normalize, resolve, sep } = require('node:path');

const repoRoot = resolve(process.env.GBDRAW_REPO || process.cwd());
const sessionPath = join(repoRoot, 'tests/test_inputs/2026-06-16_wssv.gbdraw-session.json');
const bgcSessionPath = join(repoRoot, 'tests/test_inputs/BGC0000708-BGC0000713.gbdraw-session.json');
const sparseGenbankAPath = join(repoRoot, 'tests/test_inputs/BGC0000711.gbk');
const sparseGenbankBPath = join(repoRoot, 'tests/test_inputs/BGC0000713.gbk');

const contentTypes = {
  '.html': 'text/html; charset=utf-8',
  '.js': 'text/javascript; charset=utf-8',
  '.mjs': 'text/javascript; charset=utf-8',
  '.css': 'text/css; charset=utf-8',
  '.json': 'application/json; charset=utf-8',
  '.svg': 'image/svg+xml',
  '.wasm': 'application/wasm',
  '.whl': 'application/octet-stream',
  '.data': 'application/octet-stream'
};

let server;
let baseUrl;

test.beforeAll(async () => {
  await new Promise((resolveServer, rejectServer) => {
    server = createServer((request, response) => {
      const url = new URL(request.url || '/', 'http://127.0.0.1');
      const requestedPath = normalize(decodeURIComponent(url.pathname)).replace(/^(\.\.(?:\/|\\|$))+/, '');
      const filePath = resolve(repoRoot, requestedPath.replace(/^[/\\]+/, ''));
      if (!filePath.startsWith(`${repoRoot}${sep}`) && filePath !== repoRoot) {
        response.writeHead(403);
        response.end('Forbidden');
        return;
      }
      const finalPath = filePath === repoRoot ? join(repoRoot, 'gbdraw/web/index.html') : filePath;
      if (!existsSync(finalPath)) {
        response.writeHead(404);
        response.end('Not found');
        return;
      }
      response.writeHead(200, {
        'Content-Type': contentTypes[extname(finalPath)] || 'application/octet-stream'
      });
      createReadStream(finalPath).pipe(response);
    });
    server.once('error', rejectServer);
    server.listen(0, '127.0.0.1', () => {
      const address = server.address();
      baseUrl = `http://127.0.0.1:${address.port}`;
      resolveServer();
    });
  });
});

test.afterAll(async () => {
  await new Promise((resolveClose) => server.close(resolveClose));
});

const inspectSparseDepthResult = async (page) => page.evaluate(() => {
  const app = window.__GBDRAW_APP__;
  const content = app.results?.[0]?.content || '';
  const svg = new DOMParser().parseFromString(content, 'image/svg+xml').documentElement;
  const hasGroup = (id) => Boolean(svg.querySelector(`[id="${id}"]`));
  const groupFills = (id) => {
    const group = svg.querySelector(`[id="${id}"]`);
    if (!group) return [];
    return [group, ...group.querySelectorAll('[fill]')]
      .map((element) => String(element.getAttribute('fill') || '').toLowerCase())
      .filter(Boolean);
  };
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
      depthARecord1: hasGroup('depth_a_record_1'),
      depthARecord1Axis: hasGroup('depth_a_record_1_axis'),
      depthARecord2: hasGroup('depth_a_record_2'),
      depthBRecord1: hasGroup('depth_b_record_1'),
      depthBRecord2: hasGroup('depth_b_record_2'),
      depthBRecord2Axis: hasGroup('depth_b_record_2_axis')
    },
    depthAFills: groupFills('depth_a_record_1'),
    depthBFills: groupFills('depth_b_record_2'),
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

test('Show Depth stays disabled until a depth TSV is uploaded', async ({ page }) => {
  await page.goto(`${baseUrl}/gbdraw/web/index.html`, { waitUntil: 'domcontentloaded' });
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

test('Linear depth add, clear, and remove keep global sparse columns aligned', async ({ page }) => {
  await page.goto(`${baseUrl}/gbdraw/web/index.html`, { waitUntil: 'domcontentloaded' });
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

  await page.goto(`${baseUrl}/gbdraw/web/index.html`, { waitUntil: 'domcontentloaded' });
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
  expect(firstResult.depthArgs).toEqual([
    ['sample-a.tsv', ''],
    ['', 'sample-b.tsv']
  ]);

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
  expect(secondResult.depthArgs).toEqual([
    ['depth-track-files-1-1-sample-a.tsv', ''],
    ['', 'depth-track-files-2-2-sample-b.tsv']
  ]);
});

test('WSSV depth session removes stale circular depth metadata and slots', async ({ page }) => {
  page.on('dialog', (dialog) => dialog.accept());
  await page.goto(`${baseUrl}/gbdraw/web/index.html`, { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);

  await page.locator('input[accept^=".json,"]').first().setInputFiles(sessionPath);
  await page.waitForFunction(() => {
    const app = window.__GBDRAW_APP__;
    return Array.isArray(app?.files?.c_depth) && app.files.c_depth.length === 2;
  });

  const loaded = await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    return {
      names: app.files.c_depth.map((file) => file?.name || ''),
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

  expect(loaded.names).toHaveLength(2);
  expect(loaded.labels).toHaveLength(2);
  expect(loaded.labels[1]).toContain('SRR14027721');
  expect(loaded.slots.filter((slot) => slot.enabled).map((slot) => slot.trackIndex)).toEqual([0, 1]);
  expect(loaded.slots.filter((slot) => slot.enabled).map((slot) => slot.legendLabel)).toEqual(loaded.labels);

  await page.evaluate(() => window.__GBDRAW_APP__.removeCircularDepthTrack(1));
  await page.waitForFunction(() => {
    const app = window.__GBDRAW_APP__;
    return Array.isArray(app.files.c_depth) &&
      app.files.c_depth.length === 1 &&
      app.adv.circular_track_slots
        .filter((slot) => slot?.renderer === 'depth' && slot.enabled !== false)
        .every((slot) => Number(slot.params?.track_index) < 1);
  });

  const afterRemoval = await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    return {
      names: app.files.c_depth.map((file) => file?.name || ''),
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

  expect(afterRemoval.names).toHaveLength(1);
  expect(afterRemoval.labels[0]).toContain('SRR14027705');
  expect(afterRemoval.slots.filter((slot) => slot.enabled).map((slot) => slot.trackIndex)).toEqual([0]);
  expect(afterRemoval.slots.some((slot) => Number(slot.trackIndex) >= 1)).toBe(false);
  expect(afterRemoval.slots.some((slot) => String(slot.legendLabel || '').includes('SRR14027712'))).toBe(false);
});

test('BGC session keeps restored feature metadata selectable in the preview', async ({ page }) => {
  page.on('dialog', (dialog) => dialog.accept());
  await page.goto(`${baseUrl}/gbdraw/web/index.html`, { waitUntil: 'domcontentloaded' });
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
  await page.goto(`${baseUrl}/gbdraw/web/index.html`, { waitUntil: 'domcontentloaded' });
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
