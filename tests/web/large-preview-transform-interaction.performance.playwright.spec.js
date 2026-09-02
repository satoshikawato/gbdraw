const { test, expect } = require('@playwright/test');
const { createHash } = require('node:crypto');
const { readFileSync } = require('node:fs');
const os = require('node:os');
const { join, resolve } = require('node:path');
const {
  getDiagramWorkerActivity,
  openApp
} = require('./helpers/app-lifecycle.cjs');

const repoRoot = resolve(process.env.GBDRAW_REPO || process.cwd());
const fixturePath = join(
  repoRoot,
  'gbdraw',
  'web',
  'gallery',
  'sessions',
  'hepatoplasmataceae_collinear.gbdraw-session.json.gz'
);
const EXPECTED_FIXTURE_SHA = '2afba7111520b9dd7b00dffd351a1f4a21d4e0d9bfbaebc25ba7b5a937288577';
const fixtureSha = createHash('sha256').update(readFileSync(fixturePath)).digest('hex');
const OPERATION_TIMEOUT_MS = 900_000;
const INTERACTION_TIMEOUT_MS = 60_000;
const PAN_STEPS = 150;
const PAN_INTERVAL_MS = 16;
const WHEEL_STEPS = 30;
const WHEEL_INTERVAL_MS = 80;
const POST_INTERACTION_SETTLE_MS = 300;
const installInteractionProbe = (page) => page.addInitScript(() => {
  const FEATURE_SELECTOR = [
    'path[data-gbdraw-feature-id]',
    'polygon[data-gbdraw-feature-id]',
    'rect[data-gbdraw-feature-id]',
    'path[id^="f"]',
    'polygon[id^="f"]',
    'rect[id^="f"]'
  ].join(', ');
  const MATCH_SELECTOR = [
    'path[data-gbdraw-match-id]',
    'path[data-gbdraw-pairwise-match-id]',
    'path[data-match-kind]',
    'path[data-pairwise-match-style]'
  ].join(', ');
  const metrics = Object.create(null);
  const details = [];
  const lifecycle = [];
  const longTasks = [];
  let measurement = null;
  let transformInteractionActive = false;
  let heartbeatLast = performance.now();
  const nativeRequestAnimationFrame = window.requestAnimationFrame.bind(window);
  const nativeCancelAnimationFrame = window.cancelAnimationFrame.bind(window);
  const animationFrameOwners = new Map();

  const currentWrapper = () => document.querySelector('.gbdraw-preview-surface');
  const currentSvg = () => currentWrapper()?.querySelector?.('svg') || null;
  const copyMetrics = () => Object.fromEntries(
    Object.entries(metrics).map(([name, value]) => [name, Number(value || 0)])
  );
  const deltaMetrics = (before, after = metrics) => {
    const names = new Set([...Object.keys(before || {}), ...Object.keys(after || {})]);
    return Object.fromEntries([...names].sort().map((name) => [
      name,
      Number(after?.[name] || 0) - Number(before?.[name] || 0)
    ]).filter(([, value]) => value !== 0));
  };
  const quantile = (values, probability) => {
    if (values.length === 0) return 0;
    const sorted = [...values].sort((left, right) => left - right);
    const index = Math.min(
      sorted.length - 1,
      Math.max(0, Math.ceil(probability * sorted.length) - 1)
    );
    return sorted[index];
  };
  const summarizeIntervals = (values) => ({
    count: values.length,
    p50Ms: quantile(values, 0.5),
    p90Ms: quantile(values, 0.9),
    p95Ms: quantile(values, 0.95),
    p99Ms: quantile(values, 0.99),
    maximumGapMs: values.length ? Math.max(...values) : 0,
    over33ms: values.filter((value) => value > 33).length,
    over50ms: values.filter((value) => value > 50).length,
    over100ms: values.filter((value) => value > 100).length
  });
  const emptyMutationCounts = () => ({
    featureMatchStyle: 0,
    featureMatchOpacityFilter: 0,
    featureMatchHoverState: 0,
    featureMatchCursorStyle: 0,
    wrapperStyle: 0
  });
  const emptyEventCounts = () => ({
    mouseover: 0,
    mouseout: 0,
    mousemove: 0,
    pointermove: 0,
    wheel: 0
  });
  const isFeatureOrMatch = (element) => Boolean(
    element?.matches?.(FEATURE_SELECTOR) || element?.matches?.(MATCH_SELECTOR)
  );
  const styleTouches = (record, names) => {
    const pattern = new RegExp(`(?:^|;)\\s*(?:${names.join('|')})\\s*:`, 'i');
    return pattern.test(String(record.oldValue || ''))
      || pattern.test(String(record.target?.getAttribute?.('style') || ''));
  };
  const shouldCountActive = () => Boolean(
    measurement?.sustainedStarted && transformInteractionActive
  );

  const processSvgMutations = (records) => {
    if (!measurement) return;
    for (const record of records) {
      if (!isFeatureOrMatch(record.target)) continue;
      const name = String(record.attributeName || '').toLowerCase();
      const active = shouldCountActive();
      if (name === 'style') {
        measurement.totalMutations.featureMatchStyle += 1;
        if (active) measurement.activeMutations.featureMatchStyle += 1;
        if (styleTouches(record, ['opacity', 'filter'])) {
          measurement.totalMutations.featureMatchOpacityFilter += 1;
          if (active) measurement.activeMutations.featureMatchOpacityFilter += 1;
        }
        if (styleTouches(record, ['cursor'])) {
          measurement.totalMutations.featureMatchCursorStyle += 1;
          if (active) measurement.activeMutations.featureMatchCursorStyle += 1;
        }
      }
      if (name === 'opacity' || name === 'filter') {
        measurement.totalMutations.featureMatchOpacityFilter += 1;
        if (active) measurement.activeMutations.featureMatchOpacityFilter += 1;
      }
      if (name === 'data-gbdraw-hover-opacity' || name === 'data-gbdraw-hover-filter') {
        measurement.totalMutations.featureMatchHoverState += 1;
        if (active) measurement.activeMutations.featureMatchHoverState += 1;
      }
    }
  };
  const processWrapperMutations = (records) => {
    if (!measurement) return;
    for (const record of records) {
      if (record.attributeName !== 'style') continue;
      measurement.totalMutations.wrapperStyle += 1;
      if (shouldCountActive()) measurement.activeMutations.wrapperStyle += 1;
    }
  };

  let longTaskSupported = false;
  try {
    const observer = new PerformanceObserver((list) => {
      for (const entry of list.getEntries()) {
        longTasks.push({
          startTime: Number(entry.startTime || 0),
          duration: Number(entry.duration || 0)
        });
      }
    });
    observer.observe({ type: 'longtask', buffered: true });
    longTaskSupported = true;
  } catch (_error) {
    longTaskSupported = false;
  }

  setInterval(() => {
    const now = performance.now();
    const gap = now - heartbeatLast;
    heartbeatLast = now;
    if (!measurement) return;
    measurement.heartbeatCount += 1;
    measurement.maximumHeartbeatGapMs = Math.max(measurement.maximumHeartbeatGapMs, gap);
    if (performance.memory) {
      measurement.memoryHighWaterBytes = Math.max(
        measurement.memoryHighWaterBytes,
        Number(performance.memory.usedJSHeapSize || 0)
      );
    }
  }, 100);

  const nativeQuerySelectorAll = Element.prototype.querySelectorAll;
  Element.prototype.querySelectorAll = function measuredQuerySelectorAll(selector) {
    const active = shouldCountActive() && this === currentSvg();
    const result = nativeQuerySelectorAll.call(this, selector);
    if (active) {
      const key = String(selector);
      measurement.activeSvgQueries[key] = Number(measurement.activeSvgQueries[key] || 0) + 1;
    }
    return result;
  };

  ['mouseover', 'mouseout', 'mousemove', 'pointermove', 'wheel'].forEach((type) => {
    window.addEventListener(type, (event) => {
      if (!measurement) return;
      const svg = currentSvg();
      if (!svg || (!svg.contains(event.target) && !svg.contains(event.relatedTarget))) return;
      measurement.totalEvents[type] += 1;
      if (shouldCountActive()) measurement.activeEvents[type] += 1;
    }, true);
  });

  ['mousedown', 'mouseup', 'mouseleave', 'pointercancel', 'lostpointercapture'].forEach((type) => {
    window.addEventListener(type, (event) => {
      if (!measurement || event.target !== window.__GBDRAW_APP__?.canvasContainerRef) return;
      measurement.panBoundaryEvents.push({
        type,
        timestamp: performance.now(),
        clientX: Number(event.clientX || 0),
        clientY: Number(event.clientY || 0),
        buttons: Number(event.buttons || 0),
        pointerId: Number(event.pointerId || 0),
        relatedTarget: event.relatedTarget?.tagName || null,
        relatedTargetClass: event.relatedTarget?.getAttribute?.('class') || '',
        interactionActive: transformInteractionActive
      });
    }, true);
  });

  const nativeAddEventListener = EventTarget.prototype.addEventListener;
  const nativeRemoveEventListener = EventTarget.prototype.removeEventListener;
  EventTarget.prototype.addEventListener = function measuredAddEventListener(type, listener, options) {
    const result = nativeAddEventListener.call(this, type, listener, options);
    if (measurement && this === currentWrapper() && type === 'transitionend') {
      measurement.listeners.transitionAdds += 1;
      measurement.listeners.activeTransitionListeners.add(listener);
      measurement.listeners.maximumTransitionListeners = Math.max(
        measurement.listeners.maximumTransitionListeners,
        measurement.listeners.activeTransitionListeners.size
      );
    }
    if (
      measurement
      && this === currentSvg()
      && ['mouseover', 'mouseout', 'mousemove', 'click', 'keydown'].includes(type)
    ) {
      measurement.listeners.delegatedAdds += 1;
    }
    return result;
  };
  EventTarget.prototype.removeEventListener = function measuredRemoveEventListener(type, listener, options) {
    const result = nativeRemoveEventListener.call(this, type, listener, options);
    if (measurement && this === currentWrapper() && type === 'transitionend') {
      measurement.listeners.transitionRemoves += 1;
      measurement.listeners.activeTransitionListeners.delete(listener);
    }
    if (
      measurement
      && this === currentSvg()
      && ['mouseover', 'mouseout', 'mousemove', 'click', 'keydown'].includes(type)
    ) {
      measurement.listeners.delegatedRemoves += 1;
    }
    return result;
  };

  window.requestAnimationFrame = (callback) => {
    let frameId = null;
    const wrapped = (timestamp) => {
      const owner = animationFrameOwners.get(frameId);
      owner?.frames.live.delete(frameId);
      animationFrameOwners.delete(frameId);
      callback(timestamp);
    };
    frameId = nativeRequestAnimationFrame(wrapped);
    if (measurement) {
      animationFrameOwners.set(frameId, measurement);
      measurement.frames.requests += 1;
      measurement.frames.live.add(frameId);
      measurement.frames.maximumLive = Math.max(
        measurement.frames.maximumLive,
        measurement.frames.live.size
      );
    }
    return frameId;
  };
  window.cancelAnimationFrame = (frameId) => {
    const owner = animationFrameOwners.get(frameId);
    if (owner) {
      owner.frames.cancellations += 1;
      owner.frames.live.delete(frameId);
      animationFrameOwners.delete(frameId);
    }
    return nativeCancelAnimationFrame(frameId);
  };

  const nativeSetTimeout = window.setTimeout.bind(window);
  const nativeClearTimeout = window.clearTimeout.bind(window);
  window.setTimeout = (callback, delay = 0, ...args) => {
    const numericDelay = Number(delay || 0);
    const tracked = Boolean(measurement && (numericDelay === 220 || numericDelay === 260));
    let timerId = null;
    const wrapped = (...callbackArgs) => {
      if (tracked && measurement) {
        measurement.timers.liveByDelay[numericDelay].delete(timerId);
      }
      if (typeof callback === 'function') return callback(...callbackArgs);
      return undefined;
    };
    timerId = nativeSetTimeout(wrapped, delay, ...args);
    if (tracked) {
      measurement.timers.setsByDelay[numericDelay] += 1;
      measurement.timers.liveByDelay[numericDelay].add(timerId);
      measurement.timers.maximumLiveByDelay[numericDelay] = Math.max(
        measurement.timers.maximumLiveByDelay[numericDelay],
        measurement.timers.liveByDelay[numericDelay].size
      );
    }
    return timerId;
  };
  window.clearTimeout = (timerId) => {
    if (measurement) {
      for (const delay of [220, 260]) {
        if (measurement.timers.liveByDelay[delay].delete(timerId)) {
          measurement.timers.clearsByDelay[delay] += 1;
        }
      }
    }
    return nativeClearTimeout(timerId);
  };

  const beginMeasurement = (label) => {
    if (measurement) throw new Error('An interaction measurement is already active.');
    const wrapper = currentWrapper();
    const svg = currentSvg();
    if (!wrapper || !svg) throw new Error('The exact-fixture preview is not mounted.');
    const result = window.__GBDRAW_APP__?.results?.[
      Number(window.__GBDRAW_APP__?.selectedResultIndex || 0)
    ];
    const startedAt = performance.now();
    measurement = {
      label,
      startedAt,
      metricBefore: copyMetrics(),
      activeMetricBefore: null,
      activeMetricAfter: null,
      detailStart: details.length,
      lifecycleStart: lifecycle.length,
      resultReferenceBefore: result,
      resultBefore: String(result?.content || ''),
      resultCountBefore: Number(window.__GBDRAW_APP__?.results?.length || 0),
      selectedResultIndexBefore: Number(window.__GBDRAW_APP__?.selectedResultIndex || 0),
      sustainedStarted: false,
      totalEvents: emptyEventCounts(),
      activeEvents: emptyEventCounts(),
      totalMutations: emptyMutationCounts(),
      activeMutations: emptyMutationCounts(),
      activeSvgQueries: Object.create(null),
      panBoundaryEvents: [],
      rafIntervals: [],
      rafLast: null,
      heartbeatCount: 0,
      maximumHeartbeatGapMs: 0,
      memoryStartBytes: Number(performance.memory?.usedJSHeapSize || 0),
      memoryHighWaterBytes: Number(performance.memory?.usedJSHeapSize || 0),
      listeners: {
        transitionAdds: 0,
        transitionRemoves: 0,
        maximumTransitionListeners: 0,
        activeTransitionListeners: new Set(),
        delegatedAdds: 0,
        delegatedRemoves: 0
      },
      frames: {
        requests: 0,
        cancellations: 0,
        maximumLive: 0,
        live: new Set()
      },
      timers: {
        setsByDelay: { 220: 0, 260: 0 },
        clearsByDelay: { 220: 0, 260: 0 },
        maximumLiveByDelay: { 220: 0, 260: 0 },
        liveByDelay: { 220: new Set(), 260: new Set() }
      },
      svgObserver: new MutationObserver(processSvgMutations),
      wrapperObserver: new MutationObserver(processWrapperMutations)
    };
    measurement.svgObserver.observe(svg, {
      subtree: true,
      attributes: true,
      attributeOldValue: true
    });
    measurement.wrapperObserver.observe(wrapper, {
      attributes: true,
      attributeFilter: ['style'],
      attributeOldValue: true
    });
    heartbeatLast = startedAt;
    const activeMeasurement = measurement;
    const rafTick = (timestamp) => {
      if (measurement !== activeMeasurement) return;
      if (activeMeasurement.rafLast !== null) {
        activeMeasurement.rafIntervals.push(timestamp - activeMeasurement.rafLast);
      }
      activeMeasurement.rafLast = timestamp;
      if (performance.memory) {
        activeMeasurement.memoryHighWaterBytes = Math.max(
          activeMeasurement.memoryHighWaterBytes,
          Number(performance.memory.usedJSHeapSize || 0)
        );
      }
      nativeRequestAnimationFrame(rafTick);
    };
    nativeRequestAnimationFrame(rafTick);
    const style = getComputedStyle(wrapper);
    return {
      metricBefore: { ...measurement.metricBefore },
      transition: {
        property: style.transitionProperty,
        duration: style.transitionDuration,
        timingFunction: style.transitionTimingFunction,
        delay: style.transitionDelay
      }
    };
  };

  const markSustainedWindow = () => {
    if (!measurement || !transformInteractionActive) {
      throw new Error('Cannot mark a sustained window before transform interaction activation.');
    }
    processSvgMutations(measurement.svgObserver.takeRecords());
    processWrapperMutations(measurement.wrapperObserver.takeRecords());
    measurement.activeEvents = emptyEventCounts();
    measurement.activeMutations = emptyMutationCounts();
    measurement.activeSvgQueries = Object.create(null);
    measurement.activeMetricBefore = copyMetrics();
    measurement.sustainedStarted = true;
  };

  const endMeasurement = async () => {
    if (!measurement) throw new Error('No interaction measurement is active.');
    const activeMeasurement = measurement;
    await new Promise((resolve) => nativeRequestAnimationFrame(() => (
      nativeRequestAnimationFrame(resolve)
    )));
    processSvgMutations(activeMeasurement.svgObserver.takeRecords());
    processWrapperMutations(activeMeasurement.wrapperObserver.takeRecords());
    const endedAt = performance.now();
    activeMeasurement.svgObserver.disconnect();
    activeMeasurement.wrapperObserver.disconnect();
    const result = window.__GBDRAW_APP__?.results?.[
      Number(window.__GBDRAW_APP__?.selectedResultIndex || 0)
    ];
    const phaseLongTasks = longTasks.filter((entry) => (
      entry.startTime >= activeMeasurement.startedAt && entry.startTime <= endedAt
    ));
    const wrapper = currentWrapper();
    const style = wrapper ? getComputedStyle(wrapper) : null;
    const summary = {
      label: activeMeasurement.label,
      wallDurationMs: endedAt - activeMeasurement.startedAt,
      interactionActiveAtEnd: transformInteractionActive,
      raf: summarizeIntervals(activeMeasurement.rafIntervals),
      heartbeat: {
        count: activeMeasurement.heartbeatCount,
        maximumGapMs: activeMeasurement.maximumHeartbeatGapMs
      },
      longTasks: {
        supported: longTaskSupported,
        count: phaseLongTasks.length,
        totalDurationMs: phaseLongTasks.reduce((total, entry) => total + entry.duration, 0),
        longestDurationMs: phaseLongTasks.length
          ? Math.max(...phaseLongTasks.map((entry) => entry.duration))
          : 0
      },
      memory: {
        supported: Boolean(performance.memory),
        startBytes: activeMeasurement.memoryStartBytes,
        highWaterBytes: activeMeasurement.memoryHighWaterBytes,
        highWaterDeltaBytes:
          activeMeasurement.memoryHighWaterBytes - activeMeasurement.memoryStartBytes
      },
      totalEvents: { ...activeMeasurement.totalEvents },
      activeEvents: { ...activeMeasurement.activeEvents },
      totalMutations: { ...activeMeasurement.totalMutations },
      activeMutations: { ...activeMeasurement.activeMutations },
      activeSvgQueries: { ...activeMeasurement.activeSvgQueries },
      panBoundaryEvents: activeMeasurement.panBoundaryEvents.map((entry) => ({ ...entry })),
      structuralMetricDelta: deltaMetrics(activeMeasurement.metricBefore),
      activeStructuralMetricDelta: deltaMetrics(
        activeMeasurement.activeMetricBefore || {},
        activeMeasurement.activeMetricAfter || copyMetrics()
      ),
      structuralDetails: details.slice(activeMeasurement.detailStart).map((entry) => ({ ...entry })),
      lifecycleEvents: lifecycle.slice(activeMeasurement.lifecycleStart).map((entry) => ({ ...entry })),
      listeners: {
        transitionAdds: activeMeasurement.listeners.transitionAdds,
        transitionRemoves: activeMeasurement.listeners.transitionRemoves,
        maximumTransitionListeners: activeMeasurement.listeners.maximumTransitionListeners,
        activeTransitionListeners: activeMeasurement.listeners.activeTransitionListeners.size,
        delegatedAdds: activeMeasurement.listeners.delegatedAdds,
        delegatedRemoves: activeMeasurement.listeners.delegatedRemoves
      },
      frames: {
        requests: activeMeasurement.frames.requests,
        cancellations: activeMeasurement.frames.cancellations,
        maximumLive: activeMeasurement.frames.maximumLive,
        live: activeMeasurement.frames.live.size
      },
      timers: {
        setsByDelay: { ...activeMeasurement.timers.setsByDelay },
        clearsByDelay: { ...activeMeasurement.timers.clearsByDelay },
        maximumLiveByDelay: { ...activeMeasurement.timers.maximumLiveByDelay },
        liveByDelay: {
          220: activeMeasurement.timers.liveByDelay[220].size,
          260: activeMeasurement.timers.liveByDelay[260].size
        }
      },
      result: {
        sameObject: activeMeasurement.resultReferenceBefore === result,
        unchanged: activeMeasurement.resultBefore === String(result?.content || ''),
        sameCount:
          activeMeasurement.resultCountBefore
            === Number(window.__GBDRAW_APP__?.results?.length || 0),
        sameSelectedIndex:
          activeMeasurement.selectedResultIndexBefore
            === Number(window.__GBDRAW_APP__?.selectedResultIndex || 0),
        beforeLength: activeMeasurement.resultBefore.length,
        afterLength: String(result?.content || '').length
      },
      transition: style ? {
        property: style.transitionProperty,
        duration: style.transitionDuration,
        timingFunction: style.transitionTimingFunction,
        delay: style.transitionDelay
      } : null
    };
    measurement = null;
    return summary;
  };

  window.__GBDRAW_TEST_HOOKS__ = {
    onStructuralMetric(metric) {
      const name = String(metric?.name || '');
      metrics[name] = Number(metrics[name] || 0) + Number(metric?.value || 0);
      details.push({ ...metric, timestamp: performance.now() });
      if (name === 'previewTransformInteractionStartCount') {
        transformInteractionActive = true;
      } else if (name === 'previewTransformHoverCleanupCount') {
        if (measurement && transformInteractionActive && !measurement.sustainedStarted) {
          markSustainedWindow();
        }
      } else if (name === 'previewTransformInteractionEndCount') {
        if (measurement?.sustainedStarted) measurement.activeMetricAfter = copyMetrics();
        transformInteractionActive = false;
      }
    },
    onSessionLifecycleEvent(event) {
      lifecycle.push({ ...event, timestamp: performance.now() });
    }
  };

  window.__GBDRAW_LARGE_PREVIEW_INTERACTION_PROBE__ = {
    beginMeasurement,
    endMeasurement,
    snapshot() {
      return {
        metrics: copyMetrics(),
        lifecycle: lifecycle.map((entry) => ({ ...entry })),
        transformInteractionActive
      };
    }
  };
});

const loadExactSession = async (page, dialogs) => {
  const chooserPromise = page.waitForEvent('filechooser', { timeout: OPERATION_TIMEOUT_MS });
  await page.getByRole('button', { name: 'Load Session', exact: true }).click();
  const chooser = await chooserPromise;
  await chooser.setFiles(fixturePath);
  await page.waitForFunction(() => (
    window.__GBDRAW_LARGE_PREVIEW_INTERACTION_PROBE__?.snapshot?.().lifecycle.some(
      (event) => event.name === 'interactiveReady'
    )
  ), null, { timeout: OPERATION_TIMEOUT_MS });
  await page.evaluate(() => new Promise((resolveFrame) => (
    requestAnimationFrame(() => requestAnimationFrame(resolveFrame))
  )));
  expect(dialogs).toContain('Session loaded successfully!');
};

const ensureDenseRegionVisible = (page) => page.evaluate(() => {
  const wrapper = document.querySelector('.gbdraw-preview-surface');
  const feature = wrapper?.querySelector?.(
    '[data-gbdraw-feature-id]:not([data-gbdraw-auto-feature-underlay="true"])'
  );
  feature?.scrollIntoView?.({ block: 'center', inline: 'center' });
  return new Promise((resolveFrame) => requestAnimationFrame(() => requestAnimationFrame(resolveFrame)));
});

const readGeometry = (page) => page.evaluate(() => {
  const wrapper = document.querySelector('.gbdraw-preview-surface');
  const svg = wrapper?.querySelector?.('svg');
  const container = window.__GBDRAW_APP__?.canvasContainerRef;
  if (!wrapper || !svg || !container) throw new Error('Preview geometry is unavailable.');
  const containerRect = container.getBoundingClientRect();
  const wrapperRect = wrapper.getBoundingClientRect();
  const viewport = {
    left: Math.max(0, containerRect.left + 8),
    right: Math.min(window.innerWidth, containerRect.right - 8),
    top: Math.max(0, containerRect.top + 8),
    bottom: Math.min(window.innerHeight, containerRect.bottom - 8)
  };
  const featureSelector =
    '[data-gbdraw-feature-id]:not([data-gbdraw-auto-feature-underlay="true"])';
  const feature = svg.querySelector(featureSelector);
  const featureRect = feature?.getBoundingClientRect?.();
  if (!featureRect || featureRect.width <= 0 || featureRect.height <= 0) {
    throw new Error('The selected rendered Feature has no visible geometry.');
  }
  const zoomPoint = {
    x: Math.max(
      viewport.left + 4,
      Math.min(viewport.right - 4, (featureRect.left + featureRect.right) / 2)
    ),
    y: Math.max(
      viewport.top + 4,
      Math.min(viewport.bottom - 4, (featureRect.top + featureRect.bottom) / 2)
    )
  };
  const editingSelector = [
    '[data-gbdraw-feature-id]',
    'path[id^="f"]',
    'polygon[id^="f"]',
    'rect[id^="f"]',
    '[data-gbdraw-pairwise-match-id]',
    '[data-match-kind]',
    '[data-pairwise-match-style]',
    '[data-collinearity-block-id]',
    '[data-collinear-group-scope]',
    'text[data-label-editable="true"]'
  ].join(',');
  const offsets = [0, -5, 5, -10, 10, -16, 16, -24, 24, -32, 32, -45, 45, -60, 60];
  const visibleWrapperCenterY = (
    Math.max(viewport.top, wrapperRect.top) + Math.min(viewport.bottom, wrapperRect.bottom)
  ) / 2;
  const yAnchors = [visibleWrapperCenterY, zoomPoint.y];
  const centerX = (viewport.left + viewport.right) / 2;
  const xCandidates = [0, -40, 40, -80, 80, -120, 120, -170, 170, -220, 220]
    .map((delta) => centerX + delta);
  let panPoint = null;
  for (const yAnchor of yAnchors) {
    for (const yOffset of offsets) {
      const y = Math.max(viewport.top + 20, Math.min(viewport.bottom - 20, yAnchor + yOffset));
      for (const rawX of xCandidates) {
        const x = Math.max(viewport.left + 30, Math.min(viewport.right - 30, rawX));
        const target = document.elementFromPoint(x, y);
        if (!target || !wrapper.contains(target)) continue;
        if (target.closest?.(editingSelector)) continue;
        if (target.closest?.('button,input,textarea,select,a,[role="button"]')) continue;
        panPoint = { x, y, tag: target.tagName, id: target.id || '' };
        break;
      }
      if (panPoint) break;
    }
    if (panPoint) break;
  }
  if (!panPoint) throw new Error('Could not find a non-editing SVG background point.');
  const panAmplitude = Math.max(40, Math.min(
    180,
    panPoint.x - viewport.left - 20,
    viewport.right - panPoint.x - 20
  ));
  return {
    viewport,
    panPoint,
    panAmplitude,
    zoomPoint,
    renderedFeatureElementCount: svg.querySelectorAll(featureSelector).length,
    wrapperRect: wrapperRect.toJSON(),
    svgRect: svg.getBoundingClientRect().toJSON()
  };
});

const metricDelta = (before, after) => {
  const names = new Set([...Object.keys(before || {}), ...Object.keys(after || {})]);
  return Object.fromEntries([...names].sort().map((name) => [
    name,
    Number(after?.[name] || 0) - Number(before?.[name] || 0)
  ]).filter(([, value]) => value !== 0));
};

const prewarmHoverIndexes = async (page, geometry) => {
  const before = await page.evaluate(() => (
    window.__GBDRAW_LARGE_PREVIEW_INTERACTION_PROBE__.snapshot().metrics
  ));
  const feature = page.locator(
    '.gbdraw-preview-surface > svg '
      + '[data-gbdraw-feature-id]:not([data-gbdraw-auto-feature-underlay="true"])'
  ).first();
  await expect(feature).toBeVisible({ timeout: INTERACTION_TIMEOUT_MS });
  await feature.hover();
  await expect(feature).toHaveAttribute('data-gbdraw-hover-opacity', /.*/);
  await page.mouse.move(geometry.panPoint.x, geometry.panPoint.y);
  await page.waitForFunction(() => (
    document.querySelector('.gbdraw-preview-surface > svg')
      ?.querySelectorAll('[data-gbdraw-hover-opacity]').length === 0
  ));
  await page.waitForTimeout(100);
  const after = await page.evaluate(() => (
    window.__GBDRAW_LARGE_PREVIEW_INTERACTION_PROBE__.snapshot().metrics
  ));
  return metricDelta(before, after);
};

const startMeasurement = (page, label) => page.evaluate((phase) => (
  window.__GBDRAW_LARGE_PREVIEW_INTERACTION_PROBE__.beginMeasurement(phase)
), label);

const waitForInteractionStart = async (page, beforeMetrics) => {
  await page.waitForFunction((before) => {
    const metrics = window.__GBDRAW_LARGE_PREVIEW_INTERACTION_PROBE__.snapshot().metrics;
    return (
      Number(metrics.previewTransformInteractionStartCount || 0)
        > Number(before.previewTransformInteractionStartCount || 0)
      && Number(metrics.previewTransformHoverCleanupCount || 0)
        > Number(before.previewTransformHoverCleanupCount || 0)
    );
  }, beforeMetrics, { timeout: INTERACTION_TIMEOUT_MS });
};

const waitForInteractionEnd = (page, beforeMetrics) => page.waitForFunction((before) => {
  const snapshot = window.__GBDRAW_LARGE_PREVIEW_INTERACTION_PROBE__.snapshot();
  return (
    !snapshot.transformInteractionActive
    && Number(snapshot.metrics.previewTransformInteractionEndCount || 0)
      > Number(before.previewTransformInteractionEndCount || 0)
  );
}, beforeMetrics, { timeout: INTERACTION_TIMEOUT_MS });

const runPan = async (page, geometry, beforeMetrics) => {
  const { x, y } = geometry.panPoint;
  await page.mouse.move(x, y);
  await page.mouse.down({ button: 'left' });
  await waitForInteractionStart(page, beforeMetrics);
  for (let index = 1; index <= PAN_STEPS; index += 1) {
    const progress = index / PAN_STEPS;
    await page.mouse.move(
      x + geometry.panAmplitude * Math.sin(progress * Math.PI * 8),
      y + 18 * Math.sin(progress * Math.PI * 16)
    );
    await page.waitForTimeout(PAN_INTERVAL_MS);
  }
  await page.mouse.move(x, y);
  await page.mouse.up({ button: 'left' });
};

const runWheel = async (page, geometry, beforeMetrics) => {
  await page.mouse.move(geometry.zoomPoint.x, geometry.zoomPoint.y);
  await page.waitForTimeout(40);
  await page.evaluate(({ x, y, steps, intervalMs }) => new Promise((resolve, reject) => {
    let index = 0;
    const dispatchWheel = () => {
      const target = document.elementFromPoint(x, y);
      if (!target) {
        reject(new Error('The wheel target left the viewport.'));
        return;
      }
      target.dispatchEvent(new WheelEvent('wheel', {
        bubbles: true,
        cancelable: true,
        clientX: x,
        clientY: y,
        deltaY: index < 15 ? -100 : 100,
        view: window
      }));
      index += 1;
      if (index < steps) {
        window.setTimeout(dispatchWheel, intervalMs);
      } else {
        resolve();
      }
    };
    dispatchWheel();
  }), {
    x: geometry.zoomPoint.x,
    y: geometry.zoomPoint.y,
    steps: WHEEL_STEPS,
    intervalMs: WHEEL_INTERVAL_MS
  });
  await waitForInteractionStart(page, beforeMetrics);
};

const workerTotals = (activity) => ({
  constructions: Number(activity.constructions || 0),
  initializations: Number(activity.initializations || 0),
  helpers: Number(activity.helpers || 0),
  runs: Number(activity.runs || 0)
});

const subtractWorkerTotals = (before, after) => Object.fromEntries(
  Object.keys(before).map((name) => [name, Number(after[name] || 0) - Number(before[name] || 0)])
);

const assertPostSettleInteraction = async (page) => {
  const feature = page.locator(
    '.gbdraw-preview-surface > svg '
      + '[data-gbdraw-feature-id]:not([data-gbdraw-auto-feature-underlay="true"])'
  ).first();
  await page.mouse.move(5, 5);
  await feature.hover();
  await expect(feature).toHaveCSS('cursor', 'pointer');
  expect(await feature.evaluate((element) => element.style.cursor)).toBe('');
  await expect(feature).toHaveAttribute('data-gbdraw-hover-opacity', /.*/);
  await feature.click();
  const dialog = page.getByRole('dialog', { name: /Feature details:/ });
  await expect(dialog).toBeVisible({ timeout: INTERACTION_TIMEOUT_MS });
  await dialog.getByRole('button', { name: 'Close feature popup' }).click();
  await expect(dialog).toHaveCount(0);
};

const assertStructuralContract = ({ interaction, summary, workerDelta }) => {
  const allowedStructuralMetrics = new Set([
    'featureDomFullScanCount',
    'previewTransformHoverCleanupCount',
    'previewTransformHoverReconcileCount',
    'previewTransformInteractionEndCount',
    'previewTransformInteractionStartCount'
  ]);
  expect(summary.interactionActiveAtEnd).toBe(false);
  expect(summary.result.sameObject).toBe(true);
  expect(summary.result.unchanged).toBe(true);
  expect(summary.result.sameCount).toBe(true);
  expect(summary.result.sameSelectedIndex).toBe(true);
  expect(workerDelta).toEqual({ constructions: 0, initializations: 0, helpers: 0, runs: 0 });
  expect(summary.lifecycleEvents).toEqual([]);
  expect(
    Object.keys(summary.structuralMetricDelta)
      .filter((name) => !allowedStructuralMetrics.has(name))
  ).toEqual([]);
  expect(Number(summary.structuralMetricDelta.previewTransformInteractionStartCount || 0)).toBe(1);
  expect(Number(summary.structuralMetricDelta.previewTransformInteractionEndCount || 0)).toBe(1);
  expect(Number(summary.structuralMetricDelta.previewTransformHoverCleanupCount || 0)).toBe(1);
  expect(Number(summary.activeStructuralMetricDelta.featureDomFullScanCount || 0)).toBe(0);
  expect(Number(summary.activeStructuralMetricDelta.comparisonDomFullScanCount || 0)).toBe(0);
  expect(Number(summary.activeStructuralMetricDelta.featureSearchIndexBuildCount || 0)).toBe(0);
  expect(summary.activeMutations.featureMatchStyle).toBe(0);
  expect(summary.activeMutations.featureMatchOpacityFilter).toBe(0);
  expect(summary.activeMutations.featureMatchHoverState).toBe(0);
  expect(summary.activeMutations.featureMatchCursorStyle).toBe(0);
  expect(Object.values(summary.activeSvgQueries).reduce((total, count) => total + count, 0)).toBe(0);
  expect(summary.activeMutations.wrapperStyle).toBeGreaterThan(0);
  expect(summary.listeners.delegatedAdds).toBe(0);
  expect(summary.listeners.delegatedRemoves).toBe(0);
  expect(summary.listeners.activeTransitionListeners).toBe(0);
  expect(summary.frames.live).toBe(0);
  expect(summary.timers.liveByDelay).toEqual({ 220: 0, 260: 0 });
  expect(summary.raf.count).toBeGreaterThan(0);
  expect(summary.transition).toEqual({
    property: 'transform',
    duration: '0.2s',
    timingFunction: 'ease',
    delay: '0s'
  });
  if (interaction === 'wheel') {
    expect(summary.activeMutations.wrapperStyle).toBeGreaterThanOrEqual(WHEEL_STEPS - 1);
    expect(summary.listeners.transitionAdds).toBe(1);
    expect(summary.listeners.transitionRemoves).toBe(1);
    expect(summary.listeners.maximumTransitionListeners).toBe(1);
    expect(summary.timers.setsByDelay).toEqual({ 220: WHEEL_STEPS, 260: WHEEL_STEPS });
    expect(summary.timers.maximumLiveByDelay).toEqual({ 220: 1, 260: 1 });
  } else {
    expect(summary.activeMutations.wrapperStyle).toBeGreaterThanOrEqual(PAN_STEPS);
    expect(summary.listeners.transitionAdds).toBe(0);
    expect(summary.listeners.transitionRemoves).toBe(0);
    expect(summary.timers.setsByDelay).toEqual({ 220: 0, 260: 0 });
  }
};

test.describe.configure({ mode: 'serial' });

for (const interaction of ['pan', 'wheel']) {
  for (const thermalState of ['cold', 'warm']) {
    test(`${thermalState} sustained ${interaction} preserves the exact-fixture interaction contract`, async ({
      page,
      browser
    }, testInfo) => {
      test.setTimeout(1_800_000);
      page.setDefaultTimeout(INTERACTION_TIMEOUT_MS);
      await page.setViewportSize({ width: 1440, height: 1000 });
      expect(fixtureSha).toBe(EXPECTED_FIXTURE_SHA);

      const dialogs = [];
      const pageErrors = [];
      const consoleErrors = [];
      page.on('dialog', async (dialog) => {
        dialogs.push(dialog.message());
        await dialog.accept();
      });
      page.on('pageerror', (error) => pageErrors.push(String(error?.message || error)));
      page.on('console', (message) => {
        if (message.type() === 'error') consoleErrors.push(message.text());
      });

      await installInteractionProbe(page);
      await openApp(page);
      await loadExactSession(page, dialogs);
      await ensureDenseRegionVisible(page);
      const geometry = await readGeometry(page);
      const loadWorker = workerTotals(await getDiagramWorkerActivity(page));
      expect(loadWorker).toEqual({ constructions: 0, initializations: 0, helpers: 0, runs: 0 });

      const prewarmMetricDelta = thermalState === 'warm'
        ? await prewarmHoverIndexes(page, geometry)
        : {};
      const begin = await startMeasurement(page, `${thermalState}-${interaction}`);
      expect(begin.transition).toEqual({
        property: 'transform',
        duration: '0.2s',
        timingFunction: 'ease',
        delay: '0s'
      });

      const gestureStartedAt = Date.now();
      if (interaction === 'pan') {
        await runPan(page, geometry, begin.metricBefore);
      } else {
        await runWheel(page, geometry, begin.metricBefore);
      }
      await waitForInteractionEnd(page, begin.metricBefore);
      const gestureWallDurationMs = Date.now() - gestureStartedAt;
      await page.waitForTimeout(POST_INTERACTION_SETTLE_MS);
      const summary = await page.evaluate(() => (
        window.__GBDRAW_LARGE_PREVIEW_INTERACTION_PROBE__.endMeasurement()
      ));
      const afterWorker = workerTotals(await getDiagramWorkerActivity(page));
      const workerDelta = subtractWorkerTotals(loadWorker, afterWorker);
      console.log(`GBDRAW_LARGE_PREVIEW_INTERACTION_PREASSERT ${JSON.stringify({
        profile: { thermalState, interaction },
        geometry,
        summary,
        workerDelta
      })}`);
      assertStructuralContract({ interaction, summary, workerDelta });
      await assertPostSettleInteraction(page);

      expect(pageErrors).toEqual([]);
      expect(consoleErrors).toEqual([]);
      const evidence = {
        fixtureSha,
        browser: {
          project: testInfo.project.name,
          version: browser.version(),
          platform: process.platform,
          node: process.version,
          osRelease: os.release()
        },
        profile: { thermalState, interaction },
        schedule: {
          viewport: { width: 1440, height: 1000 },
          panSteps: PAN_STEPS,
          panIntervalMs: PAN_INTERVAL_MS,
          wheelSteps: WHEEL_STEPS,
          wheelIntervalMs: WHEEL_INTERVAL_MS,
          wheelDelivery: 'browser-timer-dispatched WheelEvent',
          postInteractionSettleMs: POST_INTERACTION_SETTLE_MS
        },
        geometry,
        prewarmMetricDelta,
        gestureWallDurationMs,
        workerBefore: loadWorker,
        workerAfter: afterWorker,
        workerDelta,
        summary,
        pageErrors,
        consoleErrors
      };
      console.log(`GBDRAW_LARGE_PREVIEW_INTERACTION_EVIDENCE ${JSON.stringify(evidence)}`);
      await testInfo.attach(`${thermalState}-${interaction}-evidence.json`, {
        body: Buffer.from(JSON.stringify(evidence, null, 2)),
        contentType: 'application/json'
      });
    });
  }
}
