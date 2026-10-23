const { test, expect } = require('@playwright/test');
const { readFileSync } = require('node:fs');
const { join, resolve } = require('node:path');
const { gunzipSync } = require('node:zlib');
const { evaluateWithRetainedPromise, openApp, expectSessionLoadedWithoutWarnings } = require('./helpers/app-lifecycle.cjs');
const { loadEditorLegendRows } = require('./helpers/mode-transition.cjs');

const repoRoot = resolve(process.env.GBDRAW_REPO || process.cwd());
const genbankPath = join(repoRoot, 'tests', 'test_inputs', 'HmmtDNA.gbk');
const renamedLegendCaption = [
  'WP5 renamed browser composition legend entry with deliberately wide text',
  'that forces the live legend to be measured and reflowed'
].join(' ');
// The key of the editor row the test adds; renamed, the row keeps its key
// (D-25) and draws `renamedLegendCaption`.
const addedLegendKey = 'WP5 legend entry';
// A generated row deleted after the drags: another measured reflow that keeps
// the manual composition deltas.
const postDragDeletedCaption = 'rRNA';
const legendColor = '#336699';

const blockExternalHttpRequests = async (page, baseURL) => {
  const appOrigin = new URL(baseURL).origin;
  const externalRequests = [];
  await page.context().route(/^https?:\/\//, async (route) => {
    const url = route.request().url();
    if (new URL(url).origin !== appOrigin) {
      externalRequests.push(url);
      await route.abort('blockedbyclient');
      return;
    }
    await route.continue();
  });
  return externalRequests;
};

const renderRealDiagram = async (
  page,
  mode,
  genbankText,
  { legendSide = 'right' } = {}
) => {
  const outcome = await evaluateWithRetainedPromise(page, async ({ diagramMode, source, requestedLegendSide }) => {
    const bound = (await import('./js/services/drawing/inputs.js')).registerInputFileTree.bind(null, (await import('./js/state.js')).state.inputFiles);
    const app = window.__GBDRAW_APP__;
    app.selectDrawing(diagramMode);
    await window.Vue.nextTick();
    await window.Vue.nextTick();
    app.form.prefix = `wp5-real-${diagramMode}`;
    app.form.plot_title = `WP5 real ${diagramMode} composition`;
    app.form.multi_record_canvas = false;
    app.form.legend = requestedLegendSide;
    app.form.labels_mode = 'none';
    app.form.show_labels_linear = 'none';
    app.form.suppress_gc = true;
    app.form.suppress_skew = true;
    app.form.show_gc = false;
    app.form.show_skew = false;
    app.form.show_depth = false;
    app.adv.plot_title_position = 'top';
    app.sessionTitle = `WP5 real ${diagramMode} composition`;

    const input = new File([source], `wp5-${diagramMode}-HmmtDNA.gbk`, {
      type: 'text/plain',
      lastModified: diagramMode === 'circular' ? 51 : 52
    });
    if (diagramMode === 'circular') {
      app.cInputType = 'gb';
      app.files.c_gb = bound(input);
      app.setCircularTrackSlotsEnabled(false);
    } else {
      app.lInputType = 'gb';
      app.setLinearSeqPrimaryFile(0, 'gb', input);
      await app.setLinearComparisonGlobalAction('none');
      app.setLinearTrackSlotsEnabled(false);
    }
    await window.Vue.nextTick();
    const result = await app.runAnalysis();
    return {
      result,
      errorSummary: String(app.errorLog?.summary || ''),
      errorDetails: Array.isArray(app.errorLog?.details)
        ? app.errorLog.details.map((detail) => String(detail))
        : [],
      resultCount: app.results?.length || 0
    };
  }, { diagramMode: mode, source: genbankText, requestedLegendSide: legendSide });

  expect(outcome, JSON.stringify(outcome, null, 2)).toEqual({
    result: { status: 'ok' },
    errorSummary: '',
    errorDetails: [],
    resultCount: 1
  });
  await page.waitForFunction((expectedLegendSide) => {
    const svg = document.querySelector(
      'svg[data-gbdraw-composition-schema="1"][data-gbdraw-composition]'
    );
    if (!svg) return false;
    const metadata = JSON.parse(svg.getAttribute('data-gbdraw-composition'));
    const hasExpectedLegend = expectedLegendSide === 'none'
      ? metadata?.legend === null &&
        metadata?.legendReflow === null &&
        !svg.querySelector('[data-gbdraw-composition-role="legend"]')
      : metadata?.legend &&
        metadata?.legendReflow &&
        svg.querySelector('[data-gbdraw-composition-role="legend"]');
    return Boolean(
      metadata?.legendSide === expectedLegendSide &&
      metadata?.primary &&
      metadata?.title &&
      hasExpectedLegend &&
      svg.querySelector('[data-gbdraw-composition-role="primary"]') &&
      svg.querySelector('[data-gbdraw-composition-role="title"]')
    );
  }, legendSide, { timeout: 120000 });
};

const inspectLiveComposition = async (page, key = addedLegendKey) => (
  page.evaluate(async (expectedKey) => {
    const app = window.__GBDRAW_APP__;
    const composition = await import('./js/app/legend-layout/composition-actions.js');
    const svg = app.svgContainer?.querySelector?.('svg') || document.querySelector(
      'svg[data-gbdraw-composition-schema="1"]'
    );
    if (!svg) return { error: 'The live composition SVG is unavailable.' };

    const binding = composition.bindCompositionMetadata(svg);
    const deltas = composition.compositionUserDeltas(svg);
    // Zero shift: the Legend's recorded local bounds are the bounds the Legend
    // layout owner gives its rows (Python's layout), not its painted box.
    const { createLegendLayoutActions } = await import('./js/app/legend/layout-actions.js');
    const legendLayout = createLegendLayoutActions({ readRegisteredFontFaces: () => [] });
    await legendLayout.prepareLegendLayout();
    const laidOut = svg.getElementById('legend')
      ? legendLayout.layOutLegend(new DOMParser().parseFromString(
        new XMLSerializer().serializeToString(svg), 'image/svg+xml'
      ).documentElement)
      : null;
    const layoutLegendBounds = laidOut && {
      x: laidOut.minX, y: laidOut.minY, width: laidOut.maxX - laidOut.minX, height: laidOut.maxY - laidOut.minY
    };
    const transformPoint = (matrix, x, y) => ({
      x: matrix.a * x + matrix.c * y + matrix.e,
      y: matrix.b * x + matrix.d * y + matrix.f
    });
    const boundsInRootSvg = (target) => {
      const bounds = target.getBBox();
      const matrix = target.getCTM();
      if (!matrix || !['a', 'b', 'c', 'd', 'e', 'f'].every(
        (key) => Number.isFinite(matrix[key])
      )) {
        throw new Error('A composition target has no finite root-SVG transform matrix.');
      }
      const corners = [
        transformPoint(matrix, bounds.x, bounds.y),
        transformPoint(matrix, bounds.x + bounds.width, bounds.y),
        transformPoint(matrix, bounds.x, bounds.y + bounds.height),
        transformPoint(matrix, bounds.x + bounds.width, bounds.y + bounds.height)
      ];
      const x = Math.min(...corners.map((point) => point.x));
      const y = Math.min(...corners.map((point) => point.y));
      const right = Math.max(...corners.map((point) => point.x));
      const bottom = Math.max(...corners.map((point) => point.y));
      return { x, y, width: right - x, height: bottom - y };
    };
    const translated = (bounds, dx, dy) => bounds && ({
      x: bounds.x + dx,
      y: bounds.y + dy,
      width: bounds.width,
      height: bounds.height
    });
    const union = (boundsList) => {
      const present = boundsList.filter(Boolean);
      if (!present.length) return null;
      const x = Math.min(...present.map((bounds) => bounds.x));
      const y = Math.min(...present.map((bounds) => bounds.y));
      const right = Math.max(...present.map((bounds) => bounds.x + bounds.width));
      const bottom = Math.max(...present.map((bounds) => bounds.y + bounds.height));
      return { x, y, width: right - x, height: bottom - y };
    };
    const primaryTargetBounds = binding.primary.targets.map(boundsInRootSvg);
    const actualBounds = {
      primary: union(primaryTargetBounds),
      legend: binding.legend.targets[0]
        ? boundsInRootSvg(binding.legend.targets[0])
        : null,
      title: binding.title.targets[0]
        ? boundsInRootSvg(binding.title.targets[0])
        : null
    };
    const automaticBounds = {
      primary: union(primaryTargetBounds.map((bounds, index) => translated(
        bounds,
        -(deltas.primary[index]?.[0] || 0),
        -(deltas.primary[index]?.[1] || 0)
      ))),
      legend: actualBounds.legend && translated(
        actualBounds.legend,
        -(deltas.legend?.[0] || 0),
        -(deltas.legend?.[1] || 0)
      ),
      title: actualBounds.title && translated(
        actualBounds.title,
        -(deltas.title?.[0] || 0),
        -(deltas.title?.[1] || 0)
      )
    };
    const localBounds = {
      legend: actualBounds.legend && translated(
        actualBounds.legend,
        -(binding.metadata.legend.automaticTranslation[0] + (deltas.legend?.[0] || 0)),
        -(binding.metadata.legend.automaticTranslation[1] + (deltas.legend?.[1] || 0))
      ),
      title: actualBounds.title && translated(
        actualBounds.title,
        -(binding.metadata.title.automaticTranslation[0] + (deltas.title?.[0] || 0)),
        -(binding.metadata.title.automaticTranslation[1] + (deltas.title?.[1] || 0))
      )
    };
    const boundsError = (actual, expected) => {
      if (!actual || !expected) return actual === expected ? 0 : Number.POSITIVE_INFINITY;
      return Math.max(...['x', 'y', 'width', 'height'].map(
        (key) => Math.abs(actual[key] - expected[key])
      ));
    };
    const containmentError = (inner, outer) => {
      if (!inner || !outer) return inner === outer ? 0 : Number.POSITIVE_INFINITY;
      return Math.max(
        0,
        outer.x - inner.x,
        outer.y - inner.y,
        inner.x + inner.width - (outer.x + outer.width),
        inner.y + inner.height - (outer.y + outer.height)
      );
    };
    const dockGap = (() => {
      const primary = automaticBounds.primary;
      const legend = automaticBounds.legend;
      if (!primary || !legend) return null;
      if (binding.metadata.legendSide === 'left') {
        return primary.x - (legend.x + legend.width);
      }
      if (binding.metadata.legendSide === 'right') {
        return legend.x - (primary.x + primary.width);
      }
      if (binding.metadata.legendSide === 'top') {
        return primary.y - (legend.y + legend.height);
      }
      if (binding.metadata.legendSide === 'bottom') {
        return legend.y - (primary.y + primary.height);
      }
      return null;
    })();
    const resultText = String(app.results?.[app.selectedResultIndex]?.content || '');
    const resultDocument = new DOMParser().parseFromString(resultText, 'image/svg+xml');
    const resultSvg = resultDocument.documentElement;
    const resultDeltas = resultDocument.querySelector('parsererror')
      ? null
      : composition.compositionUserDeltas(resultSvg);
    const legendTarget = binding.legend.targets[0] || null;
    const entry = Array.from(
      legendTarget?.querySelectorAll('g[data-legend-key]') || []
    ).find((group) => group.getAttribute('data-legend-key') === expectedKey) || null;
    const entryText = entry?.querySelector('text') || null;
    const entryFill = Array.from(entry?.querySelectorAll('path[fill]') || [])
      .map((path) => String(path.getAttribute('fill') || '').toLowerCase())
      .find((fill) => fill && fill !== 'none' && !fill.startsWith('url(')) || '';
    const viewBox = String(svg.getAttribute('viewBox') || '').trim().split(/\s+/).map(Number);
    const paintedBounds = union(Object.values(actualBounds));
    const viewBoxContainmentError = paintedBounds && viewBox.length === 4
      ? Math.max(
          0,
          viewBox[0] - paintedBounds.x,
          viewBox[1] - paintedBounds.y,
          paintedBounds.x + paintedBounds.width - (viewBox[0] + viewBox[2]),
          paintedBounds.y + paintedBounds.height - (viewBox[1] + viewBox[3])
        )
      : Number.POSITIVE_INFINITY;
    const eventAttributes = [svg, ...svg.querySelectorAll('*')].flatMap((element) => (
      element.getAttributeNames().filter((name) => /^on/i.test(name))
    ));
    return {
      error: '',
      schema: svg.getAttribute(composition.COMPOSITION_SCHEMA_ATTRIBUTE),
      metadata: binding.metadata,
      deltas,
      resultDeltas,
      resultSchema: resultSvg.getAttribute(composition.COMPOSITION_SCHEMA_ATTRIBUTE),
      actualBounds,
      automaticBounds,
      dockGap,
      legendBox: actualBounds.legend,
      metadataBoundsErrors: {
        primary: boundsError(automaticBounds.primary, binding.metadata.primary.finalBounds),
        legend: boundsError(localBounds.legend, binding.metadata.legend?.localBounds || null),
        legendLayout: boundsError(layoutLegendBounds, binding.metadata.legend?.localBounds || null),
        title: boundsError(localBounds.title, binding.metadata.title?.localBounds || null)
      },
      primaryMetadataContainmentError: containmentError(
        automaticBounds.primary,
        binding.metadata.primary.finalBounds
      ),
      entry: entry
        ? {
            caption: entry.getAttribute('data-legend-key'),
            fill: entryFill,
            text: entryText?.textContent || '',
            textWidth: entryText?.getBBox?.().width || 0
          }
        : null,
      roleCounts: {
        primary: binding.primary.targets.length,
        legend: binding.legend.targets.length,
        title: binding.title.targets.length
      },
      viewBox,
      viewBoxContainmentError,
      safePreview: !svg.querySelector('script') && eventAttributes.length === 0
    };
  }, key)
);

const deletePostDragLegendEntry = async (page) => {
  await evaluateWithRetainedPromise(page, async (caption) => {
    const app = window.__GBDRAW_APP__;
    await app.deleteLegendEntry(app.legendEntries.findIndex((entry) => entry.caption === caption));
  }, postDragDeletedCaption);
  await page.waitForFunction((caption) => {
    const app = window.__GBDRAW_APP__;
    const svg = document.querySelector('svg[data-gbdraw-composition-schema="1"]');
    return !app.legendEntries.some((entry) => entry.caption === caption) &&
      Array.from(svg?.querySelectorAll('g[data-legend-key]') || [])
        .filter((entry) => entry.getAttribute('data-legend-key') === caption)
        .every((entry) => entry.getAttribute('display') === 'none');
  }, postDragDeletedCaption, { timeout: 180000 });
};

const loadAndRenameLegendEntry = async (page) => {
  const addedCaption = addedLegendKey;
  // An editor row from a Session (R15-2 retired Add legend item).
  await loadEditorLegendRows(page, [[addedCaption, legendColor]]);
  await page.waitForFunction((caption) => (
    window.__GBDRAW_APP__.legendEntries.some((entry) => entry.caption === caption) &&
    Array.from(document.querySelectorAll('g[data-legend-key]'))
      .some((entry) => entry.getAttribute('data-legend-key') === caption)
  ), addedCaption, { timeout: 180000 });
  const added = await inspectLiveComposition(page, addedCaption);

  const addedIndex = await page.evaluate((caption) => (
    window.__GBDRAW_APP__.legendEntries.findIndex((entry) => entry.caption === caption)
  ), addedCaption);
  expect(addedIndex).toBeGreaterThanOrEqual(0);
  await page.evaluate(async ({ index, caption }) => {
    await window.__GBDRAW_APP__.renameLegendEntry(index, caption);
  }, { index: addedIndex, caption: renamedLegendCaption });
  await page.waitForFunction(({ key, caption }) => (
    window.__GBDRAW_APP__.legendEntries.some((entry) => entry.caption === caption) &&
    Array.from(document.querySelectorAll('g[data-legend-key]'))
      .some((entry) => entry.getAttribute('data-legend-key') === key && entry.querySelector('text')?.textContent === caption)
  ), { key: addedCaption, caption: renamedLegendCaption }, { timeout: 180000 });
  const renamed = await inspectLiveComposition(page);

  expect(added.entry).toMatchObject({
    caption: addedCaption,
    fill: legendColor,
    text: addedCaption
  });
  expect(renamed.entry).toMatchObject({
    caption: addedLegendKey,
    fill: legendColor,
    text: renamedLegendCaption
  });
  expect(renamed.entry.textWidth).toBeGreaterThan(added.entry.textWidth);
  expect(
    renamed.legendBox.width !== added.legendBox.width ||
    renamed.legendBox.height !== added.legendBox.height
  ).toBe(true);
  expect(renamed.metadataBoundsErrors.legendLayout).toBeLessThan(1e-9);
  expect(renamed.metadata.legend.localBounds).not.toEqual(added.metadata.legend.localBounds);
  return { added, renamed };
};

const dragCompositionRole = async (page, role, dx, dy) => {
  await page.evaluate(async ({ targetRole, deltaX, deltaY }) => {
    const app = window.__GBDRAW_APP__;
    app.layoutRepositionMode = true;
    await window.Vue.nextTick();
    const svg = app.svgContainer.querySelector('svg');
    const candidates = Array.from(
      svg.querySelectorAll(`[data-gbdraw-composition-role="${targetRole}"]`)
    );
    const target = targetRole === 'primary'
      ? candidates.find((element) => element.id !== 'length_bar') || candidates[0]
      : candidates[0];
    if (!target) throw new Error(`No ${targetRole} composition drag target was found.`);

    const bounds = target.getBoundingClientRect();
    const startX = Number.isFinite(bounds.left) ? bounds.left + Math.min(bounds.width / 2, 8) : 100;
    const startY = Number.isFinite(bounds.top) ? bounds.top + Math.min(bounds.height / 2, 8) : 100;
    const event = (type, x, y, buttons) => new MouseEvent(type, {
      bubbles: true,
      cancelable: true,
      clientX: x,
      clientY: y,
      buttons,
      view: window
    });

    target.dispatchEvent(event('mousedown', startX, startY, 1));
    await new Promise((resolveFrame) => requestAnimationFrame(resolveFrame));
    const moveTarget = targetRole === 'legend' ? svg : document;
    moveTarget.dispatchEvent(event('mousemove', startX + deltaX, startY + deltaY, 1));
    await new Promise((resolveFrame) => requestAnimationFrame(resolveFrame));
    moveTarget.dispatchEvent(event('mouseup', startX + deltaX, startY + deltaY, 0));
    await new Promise((resolveWait) => setTimeout(resolveWait, 100));
  }, { targetRole: role, deltaX: dx, deltaY: dy });

  const expectedLabel = role === 'primary' ? 'Move diagram' : `Move ${role === 'title' ? 'plot title' : role}`;
  await expect.poll(
    () => page.evaluate(() => window.__GBDRAW_HISTORY__?.undoLabel?.() || ''),
    { timeout: 30000 }
  ).toBe(expectedLabel);
};

const saveSession = async (page) => {
  const downloadPromise = page.waitForEvent('download', { timeout: 120000 });
  const outcome = await evaluateWithRetainedPromise(page, async () => ({
    result: await window.__GBDRAW_APP__.saveSessionWithTitle(),
    errorLog: window.__GBDRAW_APP__.errorLog
  }));
  expect(outcome, JSON.stringify(outcome, null, 2)).toMatchObject({
    result: { status: 'saved' },
    errorLog: null
  });
  const download = await downloadPromise;
  const path = await download.path();
  expect(path).toBeTruthy();
  const buffer = readFileSync(path);
  const session = JSON.parse(gunzipSync(buffer).toString('utf8'));
  return {
    file: {
      name: download.suggestedFilename(),
      mimeType: 'application/gzip',
      buffer
    },
    session
  };
};

const reloadSession = async (page, file, mode, legendSide) => {
  await page.reload({ waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);
  await page.evaluate(() => {
    const history = window.__GBDRAW_HISTORY__;
    if (!history?.initializeIntentBaseline) {
      throw new Error('The session-import history boundary is unavailable.');
    }
    const originalInitializeIntentBaseline = history.initializeIntentBaseline;
    window.__GBDRAW_SESSION_IMPORT_COMPLETE__ = null;
    history.initializeIntentBaseline = async (label, ...args) => {
      try {
        const result = await originalInitializeIntentBaseline(label, ...args);
        if (label === 'Loaded session') {
          history.initializeIntentBaseline = originalInitializeIntentBaseline;
          window.__GBDRAW_SESSION_IMPORT_COMPLETE__ = { status: 'ok', label };
        }
        return result;
      } catch (error) {
        history.initializeIntentBaseline = originalInitializeIntentBaseline;
        window.__GBDRAW_SESSION_IMPORT_COMPLETE__ = {
          status: 'error',
          message: String(error?.message || error)
        };
        throw error;
      }
    };
  });
  const sessionInput = page.locator(
    'input[type="file"][accept*="application/json"][accept*="application/gzip"]'
  );
  await expect(sessionInput).toHaveCount(1);
  await sessionInput.setInputFiles(file);
  await expectSessionLoadedWithoutWarnings(page, { timeout: 180000 });
  await page.waitForFunction(
    () => window.__GBDRAW_SESSION_IMPORT_COMPLETE__ !== null,
    null,
    { timeout: 180000 }
  );
  const importCompletion = await page.evaluate(() => window.__GBDRAW_SESSION_IMPORT_COMPLETE__);
  expect(importCompletion).toEqual({ status: 'ok', label: 'Loaded session' });
  await page.waitForFunction(({ expectedMode, expectedLegendSide, key, caption }) => {
    const app = window.__GBDRAW_APP__;
    const svg = document.querySelector('svg[data-gbdraw-composition-schema="1"]');
    if (!svg) return false;
    const liveMetadata = JSON.parse(svg.getAttribute('data-gbdraw-composition'));
    const resultText = String(app?.results?.[app.selectedResultIndex]?.content || '');
    const resultDocument = new DOMParser().parseFromString(resultText, 'image/svg+xml');
    const resultSvg = resultDocument.documentElement;
    const resultMetadata = resultDocument.querySelector('parsererror')
      ? null
      : JSON.parse(resultSvg.getAttribute('data-gbdraw-composition'));
    return Boolean(
      app?.mode === expectedMode &&
      app?.results?.length === 1 &&
      app.form.legend === expectedLegendSide &&
      liveMetadata.legendSide === expectedLegendSide &&
      resultMetadata?.legendSide === expectedLegendSide &&
      Array.from(document.querySelectorAll('g[data-legend-key]'))
        .some((entry) => entry.getAttribute('data-legend-key') === key && entry.querySelector('text')?.textContent === caption)
    );
  }, {
    expectedMode: mode,
    expectedLegendSide: legendSide,
    key: addedLegendKey,
    caption: renamedLegendCaption
  }, { timeout: 180000 });
};

const downloadPlainSvg = async (page) => {
  const downloadPromise = page.waitForEvent('download', { timeout: 120000 });
  await page.evaluate(() => window.__GBDRAW_APP__.downloadSVG());
  const download = await downloadPromise;
  expect(download.suggestedFilename()).toMatch(/\.svg$/i);
  const path = await download.path();
  expect(path).toBeTruthy();
  return readFileSync(path, 'utf8');
};

const inspectExportedComposition = async (page, svgText) => (
  page.evaluate(async ({ source, key }) => {
    const composition = await import('./js/app/legend-layout/composition-actions.js');
    const documentNode = new DOMParser().parseFromString(source, 'image/svg+xml');
    const svg = documentNode.documentElement;
    if (documentNode.querySelector('parsererror')) return { error: 'Exported SVG parse error.' };
    const binding = composition.bindCompositionMetadata(svg);
    const entry = Array.from(svg.querySelectorAll('g[data-legend-key]'))
      .find((group) => group.getAttribute('data-legend-key') === key) || null;
    const fill = Array.from(entry?.querySelectorAll('path[fill]') || [])
      .map((path) => String(path.getAttribute('fill') || '').toLowerCase())
      .find((value) => value && value !== 'none' && !value.startsWith('url(')) || '';
    const eventAttributeCount = [svg, ...svg.querySelectorAll('*')]
      .flatMap((element) => element.getAttributeNames())
      .filter((name) => /^on/i.test(name)).length;
    return {
      error: '',
      schema: svg.getAttribute(composition.COMPOSITION_SCHEMA_ATTRIBUTE),
      metadata: binding.metadata,
      deltas: composition.compositionUserDeltas(svg),
      entry: entry ? { caption: entry.getAttribute('data-legend-key'), text: entry.querySelector('text')?.textContent || '', fill } : null,
      hasScript: Boolean(svg.querySelector('script')),
      eventAttributeCount,
      viewBox: String(svg.getAttribute('viewBox') || '').trim().split(/\s+/).map(Number)
    };
  }, { source: svgText, key: addedLegendKey })
);

const expectValidComposition = (
  snapshot,
  { legendSide = 'bottom', hasLegend = legendSide !== 'none' } = {}
) => {
  expect(snapshot.error || '').toBe('');
  expect(snapshot.schema).toBe('1');
  expect(snapshot.resultSchema).toBe('1');
  expect(snapshot.metadata).toMatchObject({
    legendSide,
    titleSide: 'top',
    spacing: {
      dockGapPx: expect.any(Number),
      edgePaddingPx: expect.any(Number),
      overlayClearancePx: expect.any(Number),
      stackGapPx: expect.any(Number),
      titleGapPx: expect.any(Number)
    }
  });
  if (hasLegend) {
    expect(snapshot.metadata.legend).toEqual(expect.any(Object));
    // Python's Legend layout inputs survive every edit, History step, Session
    // round trip, and export below.
    expect(snapshot.metadata.legendReflow).toMatchObject({
      colorRectSize: expect.any(Number),
      dpi: 96,
      fontFile: 'LiberationSans-Regular',
      fontSize: expect.any(Number),
      lineHeight: expect.any(Number),
      primaryLocalBounds: {
        minX: expect.any(Number),
        minY: expect.any(Number),
        maxX: expect.any(Number),
        maxY: expect.any(Number)
      },
      textXOffset: expect.any(Number),
      wrapWidth: expect.any(Number)
    });
    expect(snapshot.roleCounts.legend).toBe(1);
    if (['left', 'right', 'top', 'bottom'].includes(legendSide)) {
      expect(snapshot.dockGap).toBeGreaterThanOrEqual(snapshot.metadata.spacing.dockGapPx);
    }
  } else {
    expect(snapshot.metadata.legend).toBeNull();
    expect(snapshot.metadata.legendReflow).toBeNull();
    expect(snapshot.roleCounts.legend).toBe(0);
    expect(snapshot.dockGap).toBeNull();
    expect(snapshot.metadataBoundsErrors.legend).toBe(0);
  }
  expect(snapshot.roleCounts.primary).toBeGreaterThan(0);
  expect(snapshot.roleCounts.title).toBe(1);
  expect(snapshot.actualBounds.primary.width).toBeGreaterThan(0);
  expect(snapshot.actualBounds.primary.height).toBeGreaterThan(0);
  expect(snapshot.primaryMetadataContainmentError).toBeLessThan(1);
  expect(snapshot.viewBox).toHaveLength(4);
  expect(snapshot.viewBox.every(Number.isFinite)).toBe(true);
  expect(snapshot.viewBox[2]).toBeGreaterThan(0);
  expect(snapshot.viewBox[3]).toBeGreaterThan(0);
  expect(snapshot.viewBoxContainmentError).toBeLessThan(1);
  expect(snapshot.safePreview).toBe(true);
};

const expectDeltasClose = (actual, expected) => {
  expect(actual.primary).toHaveLength(expected.primary.length);
  actual.primary.forEach((vector, index) => {
    expect(vector[0]).toBeCloseTo(expected.primary[index][0], 5);
    expect(vector[1]).toBeCloseTo(expected.primary[index][1], 5);
  });
  for (const role of ['legend', 'title']) {
    if (expected[role] === null) {
      expect(actual[role]).toBeNull();
      continue;
    }
    expect(actual[role]).toHaveLength(2);
    expect(actual[role][0]).toBeCloseTo(expected[role][0], 5);
    expect(actual[role][1]).toBeCloseTo(expected[role][1], 5);
  }
};

const expectZeroDeltas = (deltas) => {
  expect(deltas.primary.length).toBeGreaterThan(0);
  deltas.primary.forEach((vector) => {
    expect(vector[0]).toBeCloseTo(0, 8);
    expect(vector[1]).toBeCloseTo(0, 8);
  });
  for (const role of ['legend', 'title']) {
    if (deltas[role] === null) continue;
    expect(deltas[role][0]).toBeCloseTo(0, 8);
    expect(deltas[role][1]).toBeCloseTo(0, 8);
  }
};

const expectRenamedEntry = (snapshot) => {
  expect(snapshot.entry).toMatchObject({
    caption: addedLegendKey,
    fill: legendColor,
    text: renamedLegendCaption
  });
};

for (const mode of ['circular', 'linear']) {
  test(`${mode} real render preserves composition edits through history, session, reset, and plain SVG export`, async ({ page, baseURL }) => {
    test.setTimeout(600000);
    const externalRequests = await blockExternalHttpRequests(page, baseURL);
    const genbankText = readFileSync(genbankPath, 'utf8');

    await openApp(page);
    // A legend side applies on Generate in both modes (D-30, PD-OI-084).
    const initialLegendSide = 'bottom';
    await renderRealDiagram(page, mode, genbankText, { legendSide: initialLegendSide });
    const fresh = await inspectLiveComposition(page, '');
    expectValidComposition(fresh, { legendSide: initialLegendSide });
    expectZeroDeltas(fresh.deltas);
    const baselineDeltas = fresh.deltas;

    const { renamed } = await loadAndRenameLegendEntry(page);
    expectValidComposition(renamed);

    await dragCompositionRole(page, 'legend', 8, -7);
    await dragCompositionRole(page, 'primary', 6, 5);
    const beforeTitleDrag = await inspectLiveComposition(page);
    await dragCompositionRole(page, 'title', 5, 6);
    const moved = await inspectLiveComposition(page);
    expectValidComposition(moved);
    expectRenamedEntry(moved);
    expect(moved.deltas.legend.some((value) => Math.abs(value) > 0.1)).toBe(true);
    expect(moved.deltas.primary.every(
      (vector) => vector.some((value) => Math.abs(value) > 0.1)
    )).toBe(true);
    expect(moved.deltas.title.some((value) => Math.abs(value) > 0.1)).toBe(true);

    await page.evaluate(async () => window.__GBDRAW_APP__.undoHistory());
    await expect.poll(
      async () => {
        const title = (await inspectLiveComposition(page)).deltas.title;
        return Math.max(...title.map((value, index) => (
          Math.abs(value - beforeTitleDrag.deltas.title[index])
        )));
      },
      { timeout: 30000 }
    ).toBeLessThan(1e-6);
    const undone = await inspectLiveComposition(page);
    expectValidComposition(undone);
    expectDeltasClose(undone.deltas, beforeTitleDrag.deltas);
    expectRenamedEntry(undone);

    await page.evaluate(async () => window.__GBDRAW_APP__.redoHistory());
    await expect.poll(
      async () => {
        const title = (await inspectLiveComposition(page)).deltas.title;
        return Math.max(...title.map((value, index) => (
          Math.abs(value - moved.deltas.title[index])
        )));
      },
      { timeout: 30000 }
    ).toBeLessThan(1e-6);
    const redone = await inspectLiveComposition(page);
    expectValidComposition(redone);
    expectDeltasClose(redone.deltas, moved.deltas);
    expectRenamedEntry(redone);

    await deletePostDragLegendEntry(page);
    const reflowedAfterDrag = await inspectLiveComposition(page);
    expectValidComposition(reflowedAfterDrag);
    expectDeltasClose(reflowedAfterDrag.deltas, moved.deltas);
    expectRenamedEntry(reflowedAfterDrag);
    expect(reflowedAfterDrag.metadataBoundsErrors.legendLayout).toBeLessThan(1e-9);

    const saved = await saveSession(page);
    expect(saved.session.format).toBe('gbdraw-session');
    expect(saved.session.renderRequest.diagramOptions.output.legend).toBe(initialLegendSide);
    expect(saved.session.ui.generatedLegendPosition).toBe('bottom');
    expect(saved.session.results).toHaveLength(1);
    expect(saved.session.results[0].content).toContain('data-gbdraw-composition-schema="1"');
    expect(saved.session.results[0].content).toContain(renamedLegendCaption);
    const savedComposition = await inspectExportedComposition(
      page,
      saved.session.results[0].content
    );
    expect(savedComposition.error || '').toBe('');
    expect(savedComposition.metadata.legendSide).toBe('bottom');
    expectDeltasClose(savedComposition.deltas, moved.deltas);

    await reloadSession(page, saved.file, mode, 'bottom');
    const loaded = await inspectLiveComposition(page);
    expectValidComposition(loaded);
    expectDeltasClose(loaded.deltas, moved.deltas);
    expectRenamedEntry(loaded);

    const exportedText = await downloadPlainSvg(page);
    const exported = await inspectExportedComposition(page, exportedText);
    expect(exported.error || '').toBe('');
    expect(exported).toMatchObject({
      schema: '1',
      entry: { caption: addedLegendKey, text: renamedLegendCaption, fill: legendColor },
      hasScript: false,
      eventAttributeCount: 0
    });
    expect(exported.viewBox).toEqual(loaded.viewBox);
    expect(exported.metadata).toEqual(loaded.metadata);
    expectDeltasClose(exported.deltas, moved.deltas);

    await page.evaluate(async () => {
      window.__GBDRAW_APP__.resetLayout();
      await window.Vue.nextTick();
    });
    const reset = await inspectLiveComposition(page);
    expectValidComposition(reset);
    expectDeltasClose(reset.deltas, baselineDeltas);
    expectRenamedEntry(reset);

    expect(externalRequests).toEqual([]);
  });
}

test('real Circular render without a legend serializes diagram/title drags and Reset', async ({ page, baseURL }) => {
  test.setTimeout(600000);
  const externalRequests = await blockExternalHttpRequests(page, baseURL);
  const genbankText = readFileSync(genbankPath, 'utf8');

  await openApp(page);
  await renderRealDiagram(page, 'circular', genbankText, { legendSide: 'none' });
  const fresh = await inspectLiveComposition(page, '');
  expectValidComposition(fresh, { legendSide: 'none', hasLegend: false });
  expectZeroDeltas(fresh.deltas);

  await dragCompositionRole(page, 'primary', 6, 5);
  await dragCompositionRole(page, 'title', 5, 6);
  const moved = await inspectLiveComposition(page, '');
  expectValidComposition(moved, { legendSide: 'none', hasLegend: false });
  expect(moved.deltas.primary.every(
    (vector) => vector.some((value) => Math.abs(value) > 0.1)
  )).toBe(true);
  expect(moved.deltas.legend).toBeNull();
  expect(moved.deltas.title.some((value) => Math.abs(value) > 0.1)).toBe(true);
  expectDeltasClose(moved.resultDeltas, moved.deltas);

  await page.evaluate(async () => {
    window.__GBDRAW_APP__.resetLayout();
    await window.Vue.nextTick();
  });
  const reset = await inspectLiveComposition(page, '');
  expectValidComposition(reset, { legendSide: 'none', hasLegend: false });
  expectZeroDeltas(reset.deltas);
  expectZeroDeltas(reset.resultDeltas);
  expect(reset.metadata.legend).toBeNull();
  expect(reset.metadata.legendReflow).toBeNull();
  expect(externalRequests).toEqual([]);
});

// OV-1605 (R3): a row rename that widens an Upper Left Legend makes the
// planner grow the canvas above the record, onto the top title's side. The
// live replan keeps the title outside the Legend, and Generate draws the same
// composition.
test('real Circular render keeps a top title outside an Upper Left Legend a renamed row moves above the record', async ({ page, baseURL }) => {
  test.setTimeout(600000);
  const externalRequests = await blockExternalHttpRequests(page, baseURL);
  await openApp(page);
  await renderRealDiagram(page, 'circular', readFileSync(genbankPath, 'utf8'), { legendSide: 'upper_left' });
  const { renamed: live } = await loadAndRenameLegendEntry(page);
  expectValidComposition(live, { legendSide: 'upper_left' });
  const overlaps = (first, second) => (
    first.x < second.x + second.width && second.x < first.x + first.width &&
    first.y < second.y + second.height && second.y < first.y + first.height
  );
  const { primary, legend, title } = live.automaticBounds;
  expect(legend.y + legend.height, 'the Legend grew the canvas above the record').toBeLessThanOrEqual(primary.y);
  expect(overlaps(title, legend), 'the live title clears the Legend').toBe(false);
  expect(title.y + title.height).toBeLessThanOrEqual(legend.y);

  const outcome = await evaluateWithRetainedPromise(page, async () => window.__GBDRAW_APP__.runAnalysis());
  expect(outcome).toEqual({ status: 'ok' });
  const generated = await inspectLiveComposition(page);
  expectValidComposition(generated, { legendSide: 'upper_left' });
  expectRenamedEntry(generated);
  for (const role of ['legend', 'title']) {
    generated.metadata[role].automaticTranslation.forEach((value, index) => {
      expect(value, `${role} translation ${index}: live equals Generate`)
        .toBeCloseTo(live.metadata[role].automaticTranslation[index], 6);
    });
  }
  expect(overlaps(generated.automaticBounds.title, generated.automaticBounds.legend)).toBe(false);
  expect(externalRequests).toEqual([]);
});
