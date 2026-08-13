const { test, expect } = require('@playwright/test');
const { join, resolve } = require('node:path');

const repoRoot = resolve(process.env.GBDRAW_REPO || process.cwd());
const wssvSessionPath = join(
  repoRoot,
  'gbdraw',
  'web',
  'gallery',
  'sessions',
  'WSSV_genome_comparison.gbdraw-session.json'
);
const FEATURE_SELECTOR = [
  'path[data-gbdraw-feature-id]',
  'polygon[data-gbdraw-feature-id]',
  'rect[data-gbdraw-feature-id]',
  'path[id^="f"]',
  'polygon[id^="f"]',
  'rect[id^="f"]'
].join(', ');

const percentile = (values, fraction) => {
  const ordered = [...values].sort((a, b) => a - b);
  if (ordered.length === 0) return 0;
  return ordered[Math.min(ordered.length - 1, Math.ceil(ordered.length * fraction) - 1)];
};

const median = (values) => percentile(values, 0.5);

const openApp = async (page) => {
  await page.addInitScript(() => {
    window.__GBDRAW_TEST_HOOKS__ = {
      historyDiagnostics: [],
      onHistoryDiagnostic(detail) {
        window.__GBDRAW_TEST_HOOKS__.historyDiagnostics.push({ ...(detail || {}) });
      }
    };
  });
  await page.goto('/gbdraw/web/index.html', { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);
  await page.waitForFunction(
    () => Object.keys(window.__GBDRAW_APP__?.paletteDefinitions || {}).length > 0,
    null,
    { timeout: 180_000 }
  );
};

const installBrowserProbe = async (page) => page.evaluate((featureSelector) => {
  const largeSvgThreshold = 1_000_000;
  const rootIds = new WeakMap();
  let nextRootId = 1;
  const probe = {
    sanitizeCalls: [],
    parseCalls: [],
    serializeCalls: [],
    resultContentUpdates: [],
    previewMounts: [],
    featureIndexBuilds: [],
    featureHandlerSetups: [],
    mainThreadPyodideLoaderPresent: typeof window.loadPyodide === 'function',
    longTasks: [],
    reset() {
      this.sanitizeCalls.length = 0;
      this.parseCalls.length = 0;
      this.serializeCalls.length = 0;
      this.resultContentUpdates.length = 0;
      this.previewMounts.length = 0;
      this.featureIndexBuilds.length = 0;
      this.featureHandlerSetups.length = 0;
      this.longTasks.length = 0;
      if (window.__GBDRAW_TEST_HOOKS__) {
        window.__GBDRAW_TEST_HOOKS__.historyDiagnostics.length = 0;
      }
    }
  };

  const rootId = (root) => {
    if (!rootIds.has(root)) {
      rootIds.set(root, nextRootId);
      nextRootId += 1;
    }
    return rootIds.get(root);
  };
  const isPreviewRoot = (node) => Boolean(
    node instanceof SVGSVGElement
    && node.parentElement?.classList?.contains('origin-top')
    && node.parentElement?.classList?.contains('shadow-xl')
  );
  const svgFingerprint = (value) => {
    const text = String(value || '');
    return `${text.length}:${text.slice(0, 48)}:${text.slice(-48)}`;
  };

  const originalSanitize = window.DOMPurify.sanitize;
  window.DOMPurify.sanitize = function instrumentedSanitize(value, ...args) {
    const source = String(value || '');
    if (source.length >= largeSvgThreshold && /<svg[\s>]/i.test(source)) {
      probe.sanitizeCalls.push({ identity: svgFingerprint(source), length: source.length });
    }
    return originalSanitize.call(this, value, ...args);
  };

  const originalParseFromString = DOMParser.prototype.parseFromString;
  DOMParser.prototype.parseFromString = function instrumentedParse(value, type) {
    const source = String(value || '');
    if (
      type === 'image/svg+xml'
      && source.length >= largeSvgThreshold
      && /<svg[\s>]/i.test(source)
    ) {
      probe.parseCalls.push({ identity: svgFingerprint(source), length: source.length });
    }
    return originalParseFromString.call(this, value, type);
  };

  const originalSerializeToString = XMLSerializer.prototype.serializeToString;
  XMLSerializer.prototype.serializeToString = function instrumentedSerialize(node) {
    if (node instanceof SVGSVGElement || node?.documentElement instanceof SVGSVGElement) {
      const root = node instanceof SVGSVGElement ? node : node.documentElement;
      probe.serializeCalls.push({ rootId: rootId(root) });
    }
    return originalSerializeToString.call(this, node);
  };

  const originalQuerySelectorAll = Element.prototype.querySelectorAll;
  Element.prototype.querySelectorAll = function instrumentedQuerySelectorAll(selector) {
    if (isPreviewRoot(this) && String(selector) === featureSelector) {
      probe.featureIndexBuilds.push({ rootId: rootId(this) });
    }
    return originalQuerySelectorAll.call(this, selector);
  };

  const originalAddEventListener = SVGSVGElement.prototype.addEventListener;
  SVGSVGElement.prototype.addEventListener = function instrumentedAddEventListener(type, listener, options) {
    if (isPreviewRoot(this) && type === 'mouseover') {
      probe.featureHandlerSetups.push({ rootId: rootId(this) });
    }
    return originalAddEventListener.call(this, type, listener, options);
  };

  const recordPreviewRoot = (root) => {
    if (!isPreviewRoot(root)) return;
    probe.previewMounts.push({ rootId: rootId(root) });
  };
  const observer = new MutationObserver((records) => {
    records.forEach((record) => {
      record.addedNodes.forEach((node) => {
        if (!(node instanceof Element)) return;
        recordPreviewRoot(node);
        node.querySelectorAll?.('svg').forEach(recordPreviewRoot);
      });
    });
  });
  observer.observe(document.body, { childList: true, subtree: true });

  window.Vue.watch(
    () => window.__GBDRAW_APP__.results.map((result) => String(result?.content || '')),
    (next, previous = []) => {
      const count = Math.max(next.length, previous.length);
      for (let index = 0; index < count; index += 1) {
        if (next[index] === previous[index]) continue;
        probe.resultContentUpdates.push({
          index,
          previousIdentity: svgFingerprint(previous[index]),
          nextIdentity: svgFingerprint(next[index])
        });
      }
    },
    { flush: 'sync' }
  );

  if (PerformanceObserver.supportedEntryTypes?.includes('longtask')) {
    const longTaskObserver = new PerformanceObserver((list) => {
      list.getEntries().forEach((entry) => probe.longTasks.push(entry.duration));
    });
    longTaskObserver.observe({ type: 'longtask', buffered: true });
  }

  window.__GBDRAW_WEBAPP_PERF_PROBE__ = probe;
}, FEATURE_SELECTOR);

const loadWssvSession = async (page) => {
  await page.evaluate(() => {
    const history = window.__GBDRAW_HISTORY__;
    if (!history?.captureBaseline) {
      throw new Error('The session-import History boundary is unavailable.');
    }
    const originalCaptureBaseline = history.captureBaseline;
    window.__GBDRAW_WSSV_IMPORT_COMPLETE__ = null;
    history.captureBaseline = async (label, ...args) => {
      try {
        const result = await originalCaptureBaseline(label, ...args);
        if (label === 'Loaded session') {
          history.captureBaseline = originalCaptureBaseline;
          window.__GBDRAW_WSSV_IMPORT_COMPLETE__ = { status: 'ok', label };
        }
        return result;
      } catch (error) {
        history.captureBaseline = originalCaptureBaseline;
        window.__GBDRAW_WSSV_IMPORT_COMPLETE__ = {
          status: 'error',
          message: String(error?.message || error)
        };
        throw error;
      }
    };
  });

  const input = page.locator(
    'input[type="file"][accept*="application/json"][accept*="application/gzip"]'
  );
  await expect(input).toHaveCount(1);
  const startedAt = Date.now();
  const dialogPromise = page.waitForEvent('dialog', { timeout: 180_000 });
  await input.setInputFiles(wssvSessionPath);
  const dialog = await dialogPromise;
  expect(dialog.message()).toBe('Session loaded successfully!');
  await dialog.accept();
  await page.waitForFunction(
    () => window.__GBDRAW_WSSV_IMPORT_COMPLETE__ !== null,
    null,
    { timeout: 180_000 }
  );
  const completion = await page.evaluate(() => window.__GBDRAW_WSSV_IMPORT_COMPLETE__);
  expect(completion).toEqual({ status: 'ok', label: 'Loaded session' });
  await page.waitForFunction(() => {
    const app = window.__GBDRAW_APP__;
    return Boolean(
      app?.results?.length === 1
      && document.querySelector('.shadow-xl.origin-top > svg')
    );
  }, null, { timeout: 180_000 });
  await page.evaluate(() => new Promise((resolve) => (
    requestAnimationFrame(() => requestAnimationFrame(resolve))
  )));
  return Date.now() - startedAt;
};

const getProbeSnapshot = (page) => page.evaluate(() => {
  const probe = window.__GBDRAW_WEBAPP_PERF_PROBE__;
  return {
    sanitizeCalls: [...probe.sanitizeCalls],
    parseCalls: [...probe.parseCalls],
    serializeCalls: [...probe.serializeCalls],
    resultContentUpdates: [...probe.resultContentUpdates],
    previewMounts: [...probe.previewMounts],
    featureIndexBuilds: [...probe.featureIndexBuilds],
    featureHandlerSetups: [...probe.featureHandlerSetups],
    mainThreadPyodideLoaderPresent: probe.mainThreadPyodideLoaderPresent,
    longTasks: [...probe.longTasks],
    historyDiagnostics: [
      ...(window.__GBDRAW_TEST_HOOKS__?.historyDiagnostics || [])
    ]
  };
});

const resetProbe = (page) => page.evaluate(() => window.__GBDRAW_WEBAPP_PERF_PROBE__.reset());

test.describe.configure({ mode: 'serial' });

test('WSSV restore and ordinary edits keep History and SVG work bounded', async ({ page }, testInfo) => {
  test.setTimeout(300_000);
  await openApp(page);
  await installBrowserProbe(page);
  await resetProbe(page);

  const restoreMs = await loadWssvSession(page);
  const restoreProbe = await getProbeSnapshot(page);
  expect(restoreProbe.sanitizeCalls, 'WSSV must cross one sanitizer boundary').toHaveLength(1);
  expect(restoreProbe.parseCalls, 'WSSV must cross one full SVG parse boundary').toHaveLength(1);
  expect(restoreProbe.serializeCalls.length, 'WSSV commit may serialize at most once').toBeLessThanOrEqual(1);
  expect(restoreProbe.resultContentUpdates, 'WSSV must have one final Result commit').toHaveLength(1);
  expect(restoreProbe.previewMounts, 'the committed WSSV Result must mount once').toHaveLength(1);
  expect(restoreProbe.featureIndexBuilds, 'the feature DOM index must be built once per root').toHaveLength(1);
  expect(restoreProbe.featureHandlerSetups, 'delegated handlers must be installed once per root').toHaveLength(1);
  expect(restoreProbe.mainThreadPyodideLoaderPresent).toBe(false);
  expect(await page.evaluate(() => Object.prototype.hasOwnProperty.call(
    window.__GBDRAW_APP__,
    'pyodideReady'
  ))).toBe(false);

  await resetProbe(page);
  const timing = await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    const history = window.__GBDRAW_HISTORY__;
    const beforeDiagnostics = history.getDiagnostics();
    const beginMs = [];
    const commitMs = [];
    const originalSpecies = String(app.form.species || '');
    for (let index = 0; index < 10; index += 1) {
      let startedAt = performance.now();
      const transaction = await history.begin(`WSSV setting ${index + 1}`);
      beginMs.push(performance.now() - startedAt);
      app.form.species = index % 2 === 0 ? `${originalSpecies} ${index}` : originalSpecies;
      startedAt = performance.now();
      await history.commit(transaction);
      commitMs.push(performance.now() - startedAt);
    }
    const afterDiagnostics = history.getDiagnostics();
    return {
      beginMs,
      commitMs,
      diagnosticsDelta: Object.fromEntries(
        Object.keys(afterDiagnostics).map((key) => [
          key,
          Number(afterDiagnostics[key] || 0) - Number(beforeDiagnostics[key] || 0)
        ])
      )
    };
  });
  const editProbe = await getProbeSnapshot(page);
  const historyEventTypes = editProbe.historyDiagnostics.map((event) => (
    String(event.type || event.kind || event.name || '').toLowerCase()
  ));
  expect(
    historyEventTypes.some((type) => type.includes('begin')),
    'History must emit test-only structural diagnostics for the measured edits'
  ).toBe(true);
  expect(historyEventTypes.filter((type) => /artifact.*snapshot.*build/.test(type))).toHaveLength(0);
  expect(historyEventTypes.filter((type) => /history.*svg.*(clone|serial)/.test(type))).toHaveLength(0);
  expect(historyEventTypes.filter((type) => /history.*result.*serial/.test(type))).toHaveLength(0);
  expect(historyEventTypes.filter((type) => /previous.*(re-?sign|sign)/.test(type))).toHaveLength(0);
  expect(historyEventTypes.filter((type) => /previous.*(re-?siz|siz)/.test(type))).toHaveLength(0);
  expect(timing.diagnosticsDelta.artifactCheckpointBuilds).toBe(0);
  expect(timing.diagnosticsDelta.historySvgBytes).toBe(0);
  expect(timing.diagnosticsDelta.signatureComputations).toBe(20);
  expect(timing.diagnosticsDelta.byteEstimateComputations).toBe(10);
  expect(editProbe.sanitizeCalls).toHaveLength(0);
  expect(editProbe.parseCalls).toHaveLength(0);
  expect(editProbe.serializeCalls).toHaveLength(0);
  expect(editProbe.resultContentUpdates).toHaveLength(0);
  expect(editProbe.previewMounts).toHaveLength(0);
  expect(percentile(timing.beginMs, 0.95)).toBeLessThan(50);
  expect(percentile(timing.commitMs, 0.95)).toBeLessThan(50);
  const combinedEditTimes = timing.beginMs.map((value, index) => value + timing.commitMs[index]);
  expect(median(combinedEditTimes.slice(5))).toBeLessThanOrEqual(
    median(combinedEditTimes.slice(0, 5)) * 1.2
  );

  await resetProbe(page);
  const noOpCounts = await page.evaluate(async () => {
    const history = window.__GBDRAW_HISTORY__;
    const before = history.getUndoCount();
    const beforeDiagnostics = history.getDiagnostics();
    await history.runUndoable('WSSV no-op', () => {});
    const afterDiagnostics = history.getDiagnostics();
    return {
      before,
      after: history.getUndoCount(),
      diagnosticsDelta: Object.fromEntries(
        Object.keys(afterDiagnostics).map((key) => [
          key,
          Number(afterDiagnostics[key] || 0) - Number(beforeDiagnostics[key] || 0)
        ])
      )
    };
  });
  expect(noOpCounts.after).toBe(noOpCounts.before);
  expect(noOpCounts.diagnosticsDelta.artifactCheckpointBuilds).toBe(0);
  expect(noOpCounts.diagnosticsDelta.historySvgBytes).toBe(0);
  expect(noOpCounts.diagnosticsDelta.byteEstimateComputations).toBe(0);
  const noOpProbe = await getProbeSnapshot(page);
  expect(noOpProbe.sanitizeCalls).toHaveLength(0);
  expect(noOpProbe.parseCalls).toHaveLength(0);
  expect(noOpProbe.serializeCalls).toHaveLength(0);
  expect(noOpProbe.resultContentUpdates).toHaveLength(0);
  expect(noOpProbe.previewMounts).toHaveLength(0);

  const featureId = await page.evaluate(() => (
    document.querySelector('[data-gbdraw-feature-id]')?.getAttribute('data-gbdraw-feature-id') || ''
  ));
  expect(featureId).not.toBe('');
  await resetProbe(page);
  const visibilityChanged = await page.evaluate(async (targetId) => {
    const app = window.__GBDRAW_APP__;
    app.selectedFeatureBulkVisibility = 'off';
    app.selectedFeatureIds = new Set([targetId]);
    app.selectedFeatureAnchorId = targetId;
    const changed = await app.applySelectedFeatureVisibility();
    await new Promise((resolve) => requestAnimationFrame(() => requestAnimationFrame(resolve)));
    const elements = Array.from(document.querySelectorAll('[data-gbdraw-feature-id]'))
      .filter((element) => element.getAttribute('data-gbdraw-feature-id') === targetId);
    return Boolean(changed && elements.length > 0 && elements.every(
      (element) => element.getAttribute('display') === 'none'
    ));
  }, featureId);
  expect(visibilityChanged).toBe(true);
  const featureEditProbe = await getProbeSnapshot(page);
  expect(featureEditProbe.sanitizeCalls).toHaveLength(0);
  expect(featureEditProbe.parseCalls).toHaveLength(0);
  expect(featureEditProbe.serializeCalls.length).toBeLessThanOrEqual(1);
  expect(featureEditProbe.resultContentUpdates.length).toBeLessThanOrEqual(1);
  expect(featureEditProbe.previewMounts, 'a live feature edit must not remount the SVG').toHaveLength(0);
  expect(featureEditProbe.featureIndexBuilds).toHaveLength(0);
  expect(featureEditProbe.featureHandlerSetups).toHaveLength(0);

  const evidence = {
    restoreMs,
    maxLongTaskMs: Math.max(0, ...restoreProbe.longTasks),
    restore: {
      sanitizes: restoreProbe.sanitizeCalls.length,
      parses: restoreProbe.parseCalls.length,
      serializes: restoreProbe.serializeCalls.length,
      resultContentUpdates: restoreProbe.resultContentUpdates.length,
      mounts: restoreProbe.previewMounts.length,
      featureIndexBuilds: restoreProbe.featureIndexBuilds.length,
      featureHandlerSetups: restoreProbe.featureHandlerSetups.length,
      mainThreadPyodideLoaderPresent: restoreProbe.mainThreadPyodideLoaderPresent
    },
    edits: {
      beginMs: timing.beginMs,
      commitMs: timing.commitMs,
      beginP95Ms: percentile(timing.beginMs, 0.95),
      commitP95Ms: percentile(timing.commitMs, 0.95),
      earlyMedianMs: median(combinedEditTimes.slice(0, 5)),
      laterMedianMs: median(combinedEditTimes.slice(5))
    }
  };
  console.info(`WSSV performance evidence: ${JSON.stringify(evidence)}`);
  await testInfo.attach('wssv-performance-evidence.json', {
    body: Buffer.from(JSON.stringify(evidence, null, 2)),
    contentType: 'application/json'
  });
});

test('the browser SVG sanitizer retains its malicious-input security profile', async ({ page }) => {
  await page.goto('/gbdraw/web/index.html', { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.DOMPurify?.sanitize);
  const sanitized = await page.evaluate(async () => {
    const { sanitizeSvgContent } = await import('/gbdraw/web/js/services/svg-sanitization.js');
    return sanitizeSvgContent([
      '<svg xmlns="http://www.w3.org/2000/svg" onload="window.__svg_xss = 1">',
      '<script>window.__svg_xss = 2</script>',
      '<foreignObject><div xmlns="http://www.w3.org/1999/xhtml">unsafe</div></foreignObject>',
      '<a href="javascript:window.__svg_xss = 3"><rect width="10" height="10" /></a>',
      '<rect id="safe-feature" data-gbdraw-feature-id="safe-feature" ',
      'width="10" height="10" fill="#54bcf8" onclick="window.__svg_xss = 4" />',
      '</svg>'
    ].join(''));
  });
  expect(sanitized).not.toMatch(/<script/i);
  expect(sanitized).not.toMatch(/foreignObject/i);
  expect(sanitized).not.toMatch(/\sonload=/i);
  expect(sanitized).not.toMatch(/\sonclick=/i);
  expect(sanitized).not.toMatch(/javascript:/i);
  expect(sanitized).toContain('data-gbdraw-feature-id="safe-feature"');
  expect(await page.evaluate(() => window.__svg_xss)).toBeUndefined();
});
