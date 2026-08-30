const { test, expect } = require('@playwright/test');
const { spawnSync } = require('node:child_process');
const { createHash } = require('node:crypto');
const { mkdtempSync, readFileSync, writeFileSync } = require('node:fs');
const { tmpdir } = require('node:os');
const { join, resolve } = require('node:path');
const { gunzipSync } = require('node:zlib');
const {
  getDiagramWorkerActivity,
  openApp
} = require('../helpers/app-lifecycle.cjs');

const repoRoot = resolve(process.env.GBDRAW_REPO || process.cwd());
const svgComparisonRoot = mkdtempSync(join(tmpdir(), 'gbdraw-svg-compare-'));
let svgComparisonIndex = 0;
const sourceSessionPath = join(
  repoRoot,
  'gbdraw',
  'web',
  'gallery',
  'sessions',
  'HmmtDNA_basic_circular.gbdraw-session.json'
);
const sourceSession = JSON.parse(readFileSync(sourceSessionPath, 'utf8'));
const sessionInputSelector =
  'input[type="file"][accept*="application/json"][accept*="application/gzip"]';
const DIRECT_FILL = '#c026d3';
const DIRECT_STROKE = '#059669';
const DIRECT_STROKE_WIDTH = 5;
const DIRECT_LABEL_TEXT = 'Loaded preview direct label';
const DIRECT_LEGEND_COLOR = '#f43f5e';
const DIRECT_LEGEND_STROKE = '#1d4ed8';
const DIRECT_LEGEND_STROKE_WIDTH = 3;
const DIRECT_REGENERATED_TITLE = 'Loaded preview direct edit regeneration';

const ACTIVE_CONFIG_DOMAINS = Object.freeze([
  'form',
  'adv',
  'losat',
  'cliOptions',
  'colors',
  'palette',
  'paletteInstantPreviewEnabled',
  'rules',
  'qualifierPriorityRules',
  'filterMode',
  'whitelist',
  'blacklistText',
  'losatProgram',
  'circularConservation',
  'annotationSets',
  'modeProfiles',
  'linearRecordLayout',
  'linearComparisonPlan',
  'importedComparisonResolution',
  'webEdits'
]);

const stableValue = (value) => {
  if (Array.isArray(value)) return value.map((item) => stableValue(item));
  if (!value || typeof value !== 'object') return value;
  return Object.fromEntries(
    Object.keys(value).sort().map((key) => [key, stableValue(value[key])])
  );
};

const jsonDigest = (value) => createHash('sha256')
  .update(JSON.stringify(stableValue(value)))
  .digest('hex');

const sessionCanonicalEvidence = (session) => {
  const active = {
    palette: session.config?.palette,
    colors: session.config?.colors,
    rules: (session.config?.rules || []).map((rule) => ({
      ...rule,
      fromFile: Boolean(rule.fromFile)
    })),
    qualifierPriorityRules: session.config?.qualifierPriorityRules,
    filterMode: session.config?.filterMode,
    whitelist: session.config?.whitelist,
    blacklistText: session.config?.blacklistText,
    form: {
      plot_title: session.config?.form?.plot_title,
      labels_mode: session.config?.form?.labels_mode,
      show_scale: session.config?.form?.show_scale
    },
    adv: {
      axis_stroke_width: session.config?.adv?.axis_stroke_width,
      label_font_size: session.config?.adv?.label_font_size,
      feature_width_circular: session.config?.adv?.feature_width_circular
    },
    annotationSetIds: (session.config?.annotationSets || []).map(({ id }) => id)
  };
  return {
    renderRequestSha256: jsonDigest(session.renderRequest),
    activeConfigSha256: jsonDigest(active),
    resourcePayloadsSha256: jsonDigest(Object.fromEntries(
      Object.entries(session.resources || {}).map(([resourceId, resource]) => [
        resourceId,
        {
          kind: resource.kind,
          encoding: resource.encoding,
          payloadSha256: createHash('sha256')
            .update(Buffer.from(String(resource.data || ''), 'base64'))
            .digest('hex')
        }
      ])
    )),
    featureCatalogSha256: jsonDigest(session.editorState?.featureCatalog || null),
    mode: session.renderRequest?.mode,
    requestSchema: session.renderRequest?.schema,
    recordCount: session.renderRequest?.records?.length || 0,
    comparisonCount: session.renderRequest?.comparisons?.length || 0,
    active
  };
};

const installIntentProbe = (page) => page.addInitScript(() => {
  const metrics = {};
  const details = [];
  const lifecycle = [];
  window.__GBDRAW_SESSION_INTENT_PROBE__ = {
    metrics,
    details,
    lifecycle,
    reset() {
      Object.keys(metrics).forEach((name) => delete metrics[name]);
      details.length = 0;
      lifecycle.length = 0;
    },
    snapshot() {
      return {
        metrics: { ...metrics },
        details: details.map((entry) => ({ ...entry })),
        lifecycle: lifecycle.map((entry) => ({ ...entry }))
      };
    }
  };
  window.__GBDRAW_TEST_HOOKS__ = {
    onStructuralMetric(metric) {
      const name = String(metric?.name || '');
      metrics[name] = Number(metrics[name] || 0) + Number(metric?.value || 0);
      details.push({ ...metric });
    },
    onSessionLifecycleEvent(event) {
      lifecycle.push({ ...event });
    }
  };
});

const openIntentApp = async (page) => {
  page.on('dialog', (dialog) => dialog.accept());
  await installIntentProbe(page);
  await openApp(page);
};

const loadCurrentSession = async (page, filePath, session) => {
  await page.evaluate(() => window.__GBDRAW_SESSION_INTENT_PROBE__.reset());
  await page.locator(sessionInputSelector).setInputFiles(filePath);
  await page.waitForFunction(() => (
    window.__GBDRAW_SESSION_INTENT_PROBE__?.lifecycle?.some(
      (event) => event.name === 'interactiveReady'
    )
  ), null, { timeout: 180_000 });
  await page.evaluate(() => new Promise((resolveFrame) => (
    requestAnimationFrame(() => requestAnimationFrame(resolveFrame))
  )));

  const probe = await page.evaluate(() => (
    window.__GBDRAW_SESSION_INTENT_PROBE__.snapshot()
  ));
  const ready = probe.lifecycle.find((event) => event.name === 'interactiveReady');
  expect(ready).toMatchObject({ status: 'success', degradedRecovery: false });
  expect(probe.metrics.currentWriterActiveConfigRestoreCount).toBe(1);
  expect(probe.metrics.activeConfigCanonicalOverwriteCount || 0).toBe(0);

  const suppliedDomains = ACTIVE_CONFIG_DOMAINS.filter((domain) => (
    Object.prototype.hasOwnProperty.call(session.config, domain)
  ));
  const restored = probe.details.find(
    (entry) => entry.name === 'currentWriterActiveConfigRestoreCount'
  );
  const noCanonicalOverwrite = probe.details.find(
    (entry) => entry.name === 'activeConfigCanonicalOverwriteCount'
  );
  expect(restored?.domains).toEqual(suppliedDomains);
  expect(noCanonicalOverwrite).toMatchObject({ value: 0, domains: suppliedDomains });
  return probe;
};

const saveCurrentSession = async (page, title) => {
  await page.evaluate((nextTitle) => {
    window.__GBDRAW_APP__.sessionTitle = nextTitle;
  }, title);
  const downloadPromise = page.waitForEvent('download', { timeout: 180_000 });
  const result = await page.evaluate(() => (
    window.__GBDRAW_APP__.saveSessionWithTitle()
  ));
  expect(result).toMatchObject({ status: 'saved' });
  const download = await downloadPromise;
  const path = await download.path();
  expect(path).toBeTruthy();
  const session = JSON.parse(gunzipSync(readFileSync(path)).toString('utf8'));
  expect(session).toMatchObject({
    format: 'gbdraw-session',
    version: 40,
    renderRequest: { schema: 6 }
  });
  return { path, session };
};

const settleMountedDom = (page) => page.evaluate(() => new Promise((resolveFrame) => (
  requestAnimationFrame(() => requestAnimationFrame(resolveFrame))
)));

const diagramWorkerResourceBytes = (activity) => (
  (activity.instances || []).reduce((total, instance) => (
    total
    + (instance.helpers || []).reduce(
      (helperTotal, helper) => helperTotal + Number(helper.transferredBytes || 0),
      0
    )
    + (instance.runs || []).reduce(
      (runTotal, run) => runTotal + Number(run.transferredBytes || 0),
      0
    )
  ), 0)
);

const captureIdleWorkerStage = async (page, label, evidence) => {
  const activity = await getDiagramWorkerActivity(page);
  const stage = {
    label,
    constructions: activity.constructions,
    initializations: activity.initializations,
    helpers: activity.helpers,
    runs: activity.runs,
    transferredResourceBytes: diagramWorkerResourceBytes(activity)
  };
  expect(stage, `${label}: diagram Worker must remain idle`).toEqual({
    label,
    constructions: 0,
    initializations: 0,
    helpers: 0,
    runs: 0,
    transferredResourceBytes: 0
  });
  evidence.push(stage);
  return activity;
};

const prepareLoadedPreviewDirectEditTarget = (page, expected = null) => page.evaluate(
  async (expectedTarget) => {
    const { state } = await import('/gbdraw/web/js/state.js');
    const { featureOverrideKey } = await import(
      '/gbdraw/web/js/services/feature-override-identity.js'
    );
    const { serializeCleanSvg } = await import(
      '/gbdraw/web/js/services/svg-serialization.js'
    );
    const app = window.__GBDRAW_APP__;
    const svg = state.svgContainer.value?.querySelector?.('svg');
    if (!svg) throw new Error('The loaded session did not mount its committed SVG preview.');

    const featureId = (feature) => String(
      feature?.rendered_svg_id
      || feature?.renderedSvgId
      || feature?.rendered_feature_svg_id
      || feature?.renderedFeatureSvgId
      || feature?.svg_id
      || feature?.id
      || ''
    ).trim();
    const hasFeatureElement = (id) => Boolean(
      svg.querySelector(`[data-gbdraw-feature-id="${CSS.escape(id)}"]`)
      || svg.querySelector(`[data-gbdraw-rendered-feature-id="${CSS.escape(id)}"]`)
      || svg.querySelector(`#${CSS.escape(id)}`)
      || svg.querySelector(`[id^="${CSS.escape(id)}__"]`)
    );
    const hasLegendEntry = (caption) => Boolean(
      svg.querySelector(`g[data-legend-key="${CSS.escape(caption)}"]`)
    );

    let feature = null;
    let fillLegendEntry = null;
    if (expectedTarget) {
      feature = app.extractedFeatures.find(
        (candidate) => featureId(candidate) === expectedTarget.featureId
      );
      fillLegendEntry = app.legendEntries.find(
        (entry) => entry.caption === expectedTarget.fillLegendCaption
      );
    } else {
      const groupedEntries = app.legendEntries.filter((entry) => (
        Array.isArray(entry?.featureIds)
        && entry.featureIds.length > 1
        && hasLegendEntry(String(entry.caption || ''))
      ));
      for (const entry of groupedEntries) {
        feature = app.extractedFeatures.find((candidate) => (
          entry.featureIds.includes(featureId(candidate))
          && app.canEditFeatureColor(candidate)
          && hasFeatureElement(featureId(candidate))
        ));
        if (feature) {
          fillLegendEntry = entry;
          break;
        }
      }
      if (!feature) {
        fillLegendEntry = app.legendEntries.find(
          (entry) => hasLegendEntry(String(entry?.caption || ''))
        ) || null;
        feature = app.extractedFeatures.find((candidate) => (
          app.canEditFeatureColor(candidate) && hasFeatureElement(featureId(candidate))
        ));
      }
    }
    if (!feature || !fillLegendEntry) {
      throw new Error('No loaded Feature with a supported legend entry is editable.');
    }

    const labelFeature = expectedTarget
      ? app.extractedFeatures.find(
          (candidate) => featureId(candidate) === expectedTarget.labelFeatureId
        ) || null
      : [...app.extractedFeatures].reverse().find((candidate) => (
          hasFeatureElement(featureId(candidate))
          && app.getEditableLabelByFeatureId(featureId(candidate))
        )) || null;
    const labelEntry = labelFeature
      ? app.getEditableLabelByFeatureId(featureId(labelFeature))
      : null;
    const labelVisibilityFeature = expectedTarget
      ? app.extractedFeatures.find(
          (candidate) => featureId(candidate) === expectedTarget.labelVisibilityFeatureId
        ) || null
      : null;
    const labelVisibilityEntry = labelVisibilityFeature
      ? app.getEditableLabelByFeatureId(featureId(labelVisibilityFeature))
      : null;

    const visibilityFeature = expectedTarget
      ? app.extractedFeatures.find(
          (candidate) => featureId(candidate) === expectedTarget.visibilityFeatureId
        ) || null
      : app.extractedFeatures.find((candidate) => (
          featureId(candidate) !== featureId(feature)
          && hasFeatureElement(featureId(candidate))
        )) || null;
    if (!visibilityFeature) {
      throw new Error('No second loaded Feature is available for visibility editing.');
    }

    const legendEntry = expectedTarget
      ? app.legendEntries.find((entry) => entry.caption === expectedTarget.legendCaption)
      : app.legendEntries.find((entry) => (
          entry.caption !== fillLegendEntry.caption
          && hasLegendEntry(String(entry?.caption || ''))
        )) || fillLegendEntry;
    if (!legendEntry) throw new Error('No editable loaded legend entry is available.');

    const target = {
      feature,
      featureId: featureId(feature),
      featureOverrideKey: featureOverrideKey(feature),
      featureRecordKey: String(feature.record_key || feature.recordKey || ''),
      biologicalFeatureId: String(
        feature.biological_feature_id || feature.biologicalFeatureId || ''
      ),
      visibilityFeature,
      visibilityFeatureId: featureId(visibilityFeature),
      visibilityFeatureOverrideKey: featureOverrideKey(visibilityFeature),
      fillLegendCaption: String(fillLegendEntry.caption || ''),
      labelFeature,
      labelFeatureId: labelFeature ? featureId(labelFeature) : String(
        expectedTarget?.labelFeatureId || ''
      ),
      labelKey: String(labelEntry?.key || expectedTarget?.labelKey || ''),
      labelVisibilityFeature,
      labelVisibilityFeatureId: labelVisibilityFeature
        ? featureId(labelVisibilityFeature)
        : String(expectedTarget?.labelVisibilityFeatureId || ''),
      labelVisibilityKey: String(
        labelVisibilityEntry?.key || expectedTarget?.labelVisibilityKey || ''
      ),
      legendCaption: String(legendEntry.caption || ''),
      baselineMountedContent: serializeCleanSvg(svg),
      loadedResultContent: String(
        state.results.value[state.selectedResultIndex.value]?.content || ''
      )
    };
    if (!target.featureOverrideKey) {
      throw new Error('The loaded Feature has no stable compound override identity.');
    }
    window.__GBDRAW_LOADED_PREVIEW_DIRECT_EDIT_TARGET__ = target;
    return {
      featureId: target.featureId,
      featureOverrideKey: target.featureOverrideKey,
      featureRecordKey: target.featureRecordKey,
      biologicalFeatureId: target.biologicalFeatureId,
      visibilityFeatureId: target.visibilityFeatureId,
      visibilityFeatureOverrideKey: target.visibilityFeatureOverrideKey,
      fillLegendCaption: target.fillLegendCaption,
      labelFeatureId: target.labelFeatureId,
      labelKey: target.labelKey,
      labelVisibilityFeatureId: target.labelVisibilityFeatureId,
      labelVisibilityKey: target.labelVisibilityKey,
      legendCaption: target.legendCaption
    };
  },
  expected
);

const bindLoadedPreviewDirectEditLabel = (page, expected = null) => page.evaluate(
  async (expectedTarget) => {
    const { state } = await import('/gbdraw/web/js/state.js');
    const { serializeCleanSvg } = await import(
      '/gbdraw/web/js/services/svg-serialization.js'
    );
    const app = window.__GBDRAW_APP__;
    const target = window.__GBDRAW_LOADED_PREVIEW_DIRECT_EDIT_TARGET__;
    const svg = state.svgContainer.value?.querySelector?.('svg');
    if (!target || !svg) throw new Error('The loaded direct-edit preview is unavailable.');
    const featureId = (feature) => String(feature?.svg_id || feature?.id || '').trim();
    const resolveTargets = () => {
      const candidates = app.extractedFeatures.filter((candidate) => (
        featureId(candidate) !== target.visibilityFeatureId
        && app.getEditableLabelByFeatureId(featureId(candidate))
      ));
      const labelFeature = expectedTarget?.labelFeatureId
        ? candidates.find(
            (candidate) => featureId(candidate) === expectedTarget.labelFeatureId
          )
        : candidates[0];
      const labelVisibilityFeature = expectedTarget?.labelVisibilityFeatureId
        ? candidates.find(
            (candidate) => featureId(candidate) === expectedTarget.labelVisibilityFeatureId
          )
        : [...candidates].reverse().find(
            (candidate) => featureId(candidate) !== featureId(labelFeature)
          );
      return {
        labelFeature,
        labelEntry: labelFeature
          ? app.getEditableLabelByFeatureId(featureId(labelFeature))
          : null,
        labelVisibilityFeature,
        labelVisibilityEntry: labelVisibilityFeature
          ? app.getEditableLabelByFeatureId(featureId(labelVisibilityFeature))
          : null
      };
    };
    let resolved = resolveTargets();
    if (!resolved.labelEntry || !resolved.labelVisibilityEntry) {
      app.openFeatureEditorFromList(target.feature, { clientX: 200, clientY: 200 });
      resolved = resolveTargets();
    }
    if (
      !resolved.labelFeature || !resolved.labelEntry?.key ||
      !resolved.labelVisibilityFeature || !resolved.labelVisibilityEntry?.key
    ) {
      throw new Error('Opening the loaded Feature did not initialize two editable label bindings.');
    }
    target.labelFeature = resolved.labelFeature;
    target.labelFeatureId = featureId(resolved.labelFeature);
    target.labelKey = String(resolved.labelEntry.key || '');
    target.labelVisibilityFeature = resolved.labelVisibilityFeature;
    target.labelVisibilityFeatureId = featureId(resolved.labelVisibilityFeature);
    target.labelVisibilityKey = String(resolved.labelVisibilityEntry.key || '');
    target.baselineMountedContent = serializeCleanSvg(svg);
    target.loadedResultContent = String(
      state.results.value[state.selectedResultIndex.value]?.content || ''
    );
    return {
      labelFeatureId: target.labelFeatureId,
      labelKey: target.labelKey,
      labelVisibilityFeatureId: target.labelVisibilityFeatureId,
      labelVisibilityKey: target.labelVisibilityKey
    };
  },
  expected
);

const startLoadedPreviewDirectEditProbe = (page) => page.evaluate(async () => {
  const { state } = await import('/gbdraw/web/js/state.js');
  const owned = (value) => (
    typeof window.Vue.toRaw === 'function' ? window.Vue.toRaw(value) : value
  );
  const counters = {
    directEditResultFlushCount: 0,
    directEditSvgSerializationCount: 0,
    directEditFeatureCatalogReplacementCount: 0,
    directEditExtractedFeatureReplacementCount: 0,
    directEditBiologicalFeatureReplacementCount: 0,
    directEditOrthogroupReplacementCount: 0
  };
  const ownerCounters = {
    featureCatalog: 'directEditFeatureCatalogReplacementCount',
    extractedFeatures: 'directEditExtractedFeatureReplacementCount',
    biologicalFeatures: 'directEditBiologicalFeatureReplacementCount',
    orthogroups: 'directEditOrthogroupReplacementCount'
  };
  const ownerRefs = {
    featureCatalog: state.featureCatalog,
    extractedFeatures: state.extractedFeatures,
    biologicalFeatures: state.biologicalFeatures,
    orthogroups: state.orthogroups
  };
  let active = true;
  const stopOwnerWatches = Object.entries(ownerCounters).map(([owner, counter]) => (
    window.Vue.watch(
      () => ownerRefs[owner].value,
      (value, previousValue) => {
        if (active && owned(value) !== owned(previousValue)) counters[counter] += 1;
      },
      { flush: 'sync' }
    )
  ));

  const nativeSerializeToString = XMLSerializer.prototype.serializeToString;
  const trackedSerializeToString = function trackedSerializeToString(...args) {
    if (active && String(args[0]?.localName || '').toLowerCase() === 'svg') {
      counters.directEditSvgSerializationCount += 1;
    }
    return nativeSerializeToString.apply(this, args);
  };
  XMLSerializer.prototype.serializeToString = trackedSerializeToString;

  const stopResultWatch = window.Vue.watch(
    () => state.results.value[state.selectedResultIndex.value]?.content,
    (content, previousContent) => {
      if (active && content !== previousContent) counters.directEditResultFlushCount += 1;
    },
    { flush: 'sync' }
  );
  window.__GBDRAW_LOADED_PREVIEW_DIRECT_EDIT_PROBE__ = {
    snapshot() {
      return { ...counters };
    },
    stop() {
      active = false;
      stopResultWatch();
      stopOwnerWatches.forEach((stopWatch) => stopWatch());
      if (XMLSerializer.prototype.serializeToString === trackedSerializeToString) {
        XMLSerializer.prototype.serializeToString = nativeSerializeToString;
      }
      return { ...counters };
    }
  };
  return window.__GBDRAW_LOADED_PREVIEW_DIRECT_EDIT_PROBE__.snapshot();
});

const captureLoadedPreviewDirectEditState = (page) => page.evaluate(async () => {
  const { state } = await import('/gbdraw/web/js/state.js');
  const { featureOverrideKey } = await import(
    '/gbdraw/web/js/services/feature-override-identity.js'
  );
  const target = window.__GBDRAW_LOADED_PREVIEW_DIRECT_EDIT_TARGET__;
  if (!target) throw new Error('The loaded-preview direct-edit target is not configured.');
  const svg = state.svgContainer.value?.querySelector?.('svg');
  if (!svg) throw new Error('The loaded SVG preview is no longer mounted.');
  const resultContent = String(
    state.results.value[state.selectedResultIndex.value]?.content || ''
  );
  const resultDocument = new DOMParser().parseFromString(resultContent, 'image/svg+xml');
  const resultSvg = resultDocument.documentElement;
  if (!resultSvg || resultSvg.localName === 'parsererror') {
    throw new Error('The selected Result no longer contains valid SVG.');
  }
  const digest = async (value) => {
    const bytes = new TextEncoder().encode(String(value || ''));
    const hash = new Uint8Array(await crypto.subtle.digest('SHA-256', bytes));
    return [...hash].map((entry) => entry.toString(16).padStart(2, '0')).join('');
  };
  const renderedFeatureId = (feature) => String(
    feature?.rendered_svg_id
    || feature?.renderedSvgId
    || feature?.rendered_feature_svg_id
    || feature?.renderedFeatureSvgId
    || feature?.svg_id
    || feature?.id
    || ''
  ).trim();
  const currentFeature = state.extractedFeatures.value.find(
    (feature) => featureOverrideKey(feature) === target.featureOverrideKey
  );
  const currentVisibilityFeature = state.extractedFeatures.value.find(
    (feature) => featureOverrideKey(feature) === target.visibilityFeatureOverrideKey
  );
  const featureIds = Array.from(new Set([
    target.featureId,
    renderedFeatureId(currentFeature)
  ].filter(Boolean)));
  const visibilityFeatureIds = Array.from(new Set([
    target.visibilityFeatureId,
    renderedFeatureId(currentVisibilityFeature)
  ].filter(Boolean)));
  const collectFeatureElements = (root, ids) => {
    const elements = ids.flatMap((featureId) => {
      const escaped = CSS.escape(featureId);
      return [
        ...root.querySelectorAll(`[data-gbdraw-feature-id="${escaped}"]`),
        ...root.querySelectorAll(`[data-gbdraw-rendered-feature-id="${escaped}"]`),
        ...root.querySelectorAll(`#${escaped}`),
        ...root.querySelectorAll(`[id^="${escaped}__"]`)
      ];
    });
    return Array.from(new Set(elements));
  };
  const featureAttributes = (root, ids) => collectFeatureElements(root, ids).map(
    (element) => ({
      id: String(element.id || ''),
      part: String(element.getAttribute('data-gbdraw-feature-part') || ''),
      fill: element.getAttribute('fill'),
      stroke: element.getAttribute('stroke'),
      strokeWidth: element.getAttribute('stroke-width'),
      display: element.getAttribute('display')
    })
  );
  const labelAttributes = (root, featureId, expectedText = '') => {
    const element = root.querySelector(
      `text[data-label-feature-id="${CSS.escape(featureId)}"]`
    ) || (expectedText
      ? Array.from(root.querySelectorAll('text')).find(
          (candidate) => String(candidate.textContent || '') === expectedText
        )
      : null);
    return element
      ? { text: String(element.textContent || ''), display: element.getAttribute('display') }
      : null;
  };
  const legendAttributes = (root) => {
    const group = root.querySelector(
      `g[data-legend-key="${CSS.escape(target.legendCaption)}"]`
    );
    const swatch = Array.from(group?.querySelectorAll?.('path') || []).find((path) => {
      const fill = path.getAttribute('fill');
      return fill && fill !== 'none' && !fill.startsWith('url(');
    });
    return swatch
      ? {
          fill: swatch.getAttribute('fill'),
          stroke: swatch.getAttribute('stroke'),
          strokeWidth: swatch.getAttribute('stroke-width')
        }
      : null;
  };
  const plain = (value) => JSON.parse(JSON.stringify(value ?? null));
  const historyDiagnostics = window.__GBDRAW_HISTORY__.getDiagnostics();
  const baselineDifferenceIndex = (() => {
    const baseline = String(target.baselineMountedContent || '');
    const length = Math.min(baseline.length, resultContent.length);
    for (let index = 0; index < length; index += 1) {
      if (baseline[index] !== resultContent[index]) return index;
    }
    return baseline.length === resultContent.length ? -1 : length;
  })();
  return {
    identity: {
      featureOverrideKey: target.featureOverrideKey,
      originalFeatureId: target.featureId,
      currentFeatureId: renderedFeatureId(currentFeature),
      featureIds,
      visibilityFeatureOverrideKey: target.visibilityFeatureOverrideKey,
      originalVisibilityFeatureId: target.visibilityFeatureId,
      currentVisibilityFeatureId: renderedFeatureId(currentVisibilityFeature),
      visibilityFeatureIds
    },
    result: {
      characters: resultContent.length,
      sha256: await digest(resultContent),
      baselineCharacters: String(target.baselineMountedContent || '').length,
      baselineSha256: await digest(target.baselineMountedContent),
      baselineDifferenceIndex,
      baselineDifferenceContext: baselineDifferenceIndex < 0
        ? ''
        : {
            baseline: String(target.baselineMountedContent || '').slice(
              Math.max(0, baselineDifferenceIndex - 80),
              baselineDifferenceIndex + 160
            ),
            result: resultContent.slice(
              Math.max(0, baselineDifferenceIndex - 80),
              baselineDifferenceIndex + 160
            )
          },
      equalsBaselineMounted: resultContent === target.baselineMountedContent,
      equalsLoadedResult: resultContent === target.loadedResultContent
    },
    dom: {
      feature: featureAttributes(svg, featureIds),
      visibilityFeature: featureAttributes(svg, visibilityFeatureIds),
      label: labelAttributes(
        svg,
        target.labelFeatureId,
        String(state.labelTextFeatureOverrides[target.labelFeatureId] || '')
      ),
      labelVisibility: labelAttributes(svg, target.labelVisibilityFeatureId),
      legend: legendAttributes(svg)
    },
    resultDom: {
      feature: featureAttributes(resultSvg, featureIds),
      visibilityFeature: featureAttributes(resultSvg, visibilityFeatureIds),
      label: labelAttributes(
        resultSvg,
        target.labelFeatureId,
        String(state.labelTextFeatureOverrides[target.labelFeatureId] || '')
      ),
      labelVisibility: labelAttributes(resultSvg, target.labelVisibilityFeatureId),
      legend: legendAttributes(resultSvg)
    },
    overrides: {
      fill: plain(state.featureColorOverrides[target.featureOverrideKey]),
      stroke: plain(state.featureStrokeOverrides[target.featureOverrideKey]),
      visibility: state.featureVisibilityOverrides[target.visibilityFeatureId] ?? null,
      labelText: state.labelTextFeatureOverrides[target.labelFeatureId] ?? null,
      labelVisibility:
        state.labelVisibilityOverrides[target.labelVisibilityFeatureId] ?? null,
      legendColor: state.legendColorOverrides[target.legendCaption] ?? null,
      legendStroke: plain(state.legendStrokeOverrides[target.legendCaption])
    },
    history: {
      undoCount: window.__GBDRAW_HISTORY__.getUndoCount(),
      redoCount: window.__GBDRAW_HISTORY__.getRedoCount(),
      currentCheckpointAbsent: window.__GBDRAW_HISTORY__.getCurrentCheckpoint() === null,
      artifactReplacementHistoryEntryCount: Number(
        historyDiagnostics.artifactReplacementHistoryEntryCount || 0
      ),
      artifactCheckpointBuilds: Number(historyDiagnostics.artifactCheckpointBuilds || 0),
      generatedArtifactFullCloneCount: Number(
        historyDiagnostics.generatedArtifactFullCloneCount || 0
      ),
      generatedArtifactFullSerializationCount: Number(
        historyDiagnostics.generatedArtifactFullSerializationCount || 0
      )
    },
    autoLabelReflowEnabled: state.autoLabelReflowEnabled.value,
    structural: window.__GBDRAW_LOADED_PREVIEW_DIRECT_EDIT_PROBE__?.snapshot?.() || null
  };
});

const expectBiologicalOwnersUnchanged = (snapshot) => {
  expect(snapshot.structural).toMatchObject({
    directEditFeatureCatalogReplacementCount: 0,
    directEditExtractedFeatureReplacementCount: 0,
    directEditBiologicalFeatureReplacementCount: 0,
    directEditOrthogroupReplacementCount: 0
  });
};

const expectDirectEditFlushed = (before, after) => {
  expect(after.result.sha256).not.toBe(before.result.sha256);
  expect(after.structural.directEditResultFlushCount).toBeGreaterThan(
    before.structural.directEditResultFlushCount
  );
  expect(after.structural.directEditSvgSerializationCount).toBeGreaterThan(
    before.structural.directEditSvgSerializationCount
  );
  expect(after.history.artifactReplacementHistoryEntryCount).toBe(
    before.history.artifactReplacementHistoryEntryCount
  );
  expect(after.history.artifactCheckpointBuilds).toBe(before.history.artifactCheckpointBuilds);
  expect(after.history.generatedArtifactFullCloneCount).toBe(
    before.history.generatedArtifactFullCloneCount
  );
  expect(after.history.generatedArtifactFullSerializationCount).toBe(
    before.history.generatedArtifactFullSerializationCount
  );
  expectBiologicalOwnersUnchanged(after);
};

const capturePageEvidence = (page) => page.evaluate(async () => {
  const { state } = await import('/gbdraw/web/js/state.js');
  const { resolveColorToHex } = await import('/gbdraw/web/js/app/color-utils.js');
  const { serializeCleanSvg } = await import('/gbdraw/web/js/services/svg-serialization.js');
  const app = window.__GBDRAW_APP__;
  const history = window.__GBDRAW_HISTORY__;
  const plain = (value) => JSON.parse(JSON.stringify(value ?? null));
  const sortedObject = (value) => Object.fromEntries(
    Object.entries(plain(value) || {}).sort(([left], [right]) => left.localeCompare(right))
  );
  const hashText = (text) => {
    let hash = 2166136261;
    for (let index = 0; index < text.length; index += 1) {
      hash ^= text.charCodeAt(index);
      hash = Math.imul(hash, 16777619);
    }
    return { length: text.length, hash: hash >>> 0 };
  };
  const mountedSvg = state.svgContainer.value?.querySelector?.('svg');
  const svgContent = mountedSvg
    ? serializeCleanSvg(mountedSvg)
    : String(app.results?.[app.selectedResultIndex]?.content || '');
  const svgDocument = new DOMParser().parseFromString(svgContent, 'image/svg+xml');
  const root = svgDocument.documentElement;
  if (!root || root.localName === 'parsererror') {
    throw new Error('The committed preview is not valid SVG.');
  }
  const canonicalSvg = new XMLSerializer().serializeToString(root);
  const semanticAttributes = [
    'id',
    'viewBox',
    'd',
    'points',
    'x',
    'y',
    'x1',
    'y1',
    'x2',
    'y2',
    'cx',
    'cy',
    'r',
    'rx',
    'ry',
    'transform',
    'fill',
    'fill-opacity',
    'stroke',
    'stroke-width',
    'stroke-dasharray',
    'display',
    'visibility',
    'font-size'
  ];
  const semanticRows = [...root.querySelectorAll('*')].map((element) => ({
    tag: element.localName,
    attributes: Object.fromEntries(semanticAttributes
      .filter((name) => element.hasAttribute(name))
      .map((name) => {
        const raw = String(element.getAttribute(name) || '');
        if (['fill', 'stroke'].includes(name) && !/^(none|url\(|currentcolor)/i.test(raw)) {
          return [name, String(resolveColorToHex(raw) || raw).toLowerCase()];
        }
        if (['stroke-width', 'fill-opacity', 'font-size'].includes(name)) {
          const numeric = Number(raw);
          if (Number.isFinite(numeric)) return [name, String(numeric)];
        }
        return [name, raw];
      })),
    text: element.children.length === 0 ? String(element.textContent || '') : ''
  }));
  const svgBlobUrl = URL.createObjectURL(new Blob([svgContent], { type: 'image/svg+xml' }));
  const svgImage = new Image();
  svgImage.src = svgBlobUrl;
  await new Promise((resolveImage, rejectImage) => {
    svgImage.onload = resolveImage;
    svgImage.onerror = () => rejectImage(new Error('Could not rasterize the committed SVG.'));
  });
  const rasterWidth = 600;
  const rasterHeight = Math.max(
    1,
    Math.round(rasterWidth * Number(svgImage.naturalHeight || 1) / Number(svgImage.naturalWidth || 1))
  );
  const rasterCanvas = document.createElement('canvas');
  rasterCanvas.width = rasterWidth;
  rasterCanvas.height = rasterHeight;
  const rasterContext = rasterCanvas.getContext('2d', { willReadFrequently: true });
  rasterContext.drawImage(svgImage, 0, 0, rasterWidth, rasterHeight);
  const rasterPixels = rasterContext.getImageData(0, 0, rasterWidth, rasterHeight).data;
  const rasterDigest = new Uint8Array(await crypto.subtle.digest('SHA-256', rasterPixels));
  URL.revokeObjectURL(svgBlobUrl);
  const catalogJson = JSON.stringify(state.featureCatalog.value || null);
  const catalogBytes = new TextEncoder().encode(catalogJson);
  const catalogDigest = new Uint8Array(await crypto.subtle.digest('SHA-256', catalogBytes));
  const historyDiagnostics = history.getDiagnostics();
  const normalizeStrokeOverrides = (value) => Object.fromEntries(
    Object.entries(plain(value) || {})
      .sort(([left], [right]) => left.localeCompare(right))
      .map(([key, override]) => [key, {
        ...override,
        ...(override?.strokeColor == null
          ? {}
          : { strokeColor: resolveColorToHex(String(override.strokeColor)) }),
        ...(override?.originalStrokeColor == null
          ? {}
          : { originalStrokeColor: resolveColorToHex(String(override.originalStrokeColor)) })
      }])
  );

  return {
    active: {
      modeAndInput: {
        mode: state.mode.value,
        inputType: state.mode.value === 'circular'
          ? state.cInputType.value
          : state.lInputType.value
      },
      palette: {
        selected: state.selectedPalette.value,
        currentColors: sortedObject(state.currentColors.value),
        selectedDefaults: sortedObject(
          state.paletteDefinitions.value?.[state.selectedPalette.value] || {}
        ),
        instantPreview: state.paletteInstantPreviewEnabled.value,
        appliedName: state.appliedPaletteName.value,
        appliedColors: sortedObject(state.appliedPaletteColors.value),
        pendingName: state.pendingPaletteName.value,
        pendingColors: sortedObject(state.pendingPaletteColors.value)
      },
      specificRules: state.manualSpecificRules.map((rule) => ({
        ...plain(rule),
        fromFile: Boolean(rule.fromFile)
      })),
      qualifierPriorityRules: plain(state.manualPriorityRules),
      filters: {
        mode: state.filterMode.value,
        whitelist: plain(state.manualWhitelist),
        blacklistText: state.manualBlacklist.value
      },
      form: {
        plot_title: state.form.plot_title,
        labels_mode: state.form.labels_mode,
        show_scale: state.form.show_scale,
        legend: state.form.legend,
        track_type: state.form.track_type
      },
      adv: {
        axis_stroke_width: state.adv.axis_stroke_width,
        label_font_size: state.adv.label_font_size,
        feature_width_circular: state.adv.feature_width_circular,
        plot_title_position: state.adv.plot_title_position
      },
      annotations: state.annotationSets.map((set) => ({
        id: set.id,
        annotations: set.annotations.map((annotation) => ({
          id: annotation.id,
          mark: annotation.mark,
          target: {
            kind: annotation.target?.kind,
            record: annotation.target?.record ?? null,
            start: annotation.target?.start,
            end: annotation.target?.end
          }
        }))
      })),
      tracks: {
        enabled: state.adv.circular_track_slots_enabled,
        axisIndex: state.adv.circular_track_slots_axis_index,
        slots: state.adv.circular_track_slots.map((slot) => ({
          id: slot.id,
          renderer: slot.renderer,
          enabled: slot.enabled,
          side: slot.side,
          width: slot.width,
          radius: slot.radius
        }))
      },
      layout: {
        preferences: plain(state.layoutPreferences),
        linearRecordLayout: {
          enabled: state.linearRecordLayoutEnabled.value,
          recordGap: state.linearRecordGap.value,
          rows: plain(state.linearRecordRows)
        }
      },
      editorOverrides: {
        fills: sortedObject(state.featureColorOverrides),
        strokes: normalizeStrokeOverrides(state.featureStrokeOverrides),
        visibility: sortedObject(state.featureVisibilityOverrides),
        labelText: sortedObject(state.labelTextFeatureOverrides),
        labelVisibility: sortedObject(state.labelVisibilityOverrides),
        legendColors: sortedObject(state.legendColorOverrides),
        legendStrokes: normalizeStrokeOverrides(state.legendStrokeOverrides),
        orthogroupNames: sortedObject(state.orthogroupNameOverrides),
        orthogroupDescriptions: sortedObject(state.orthogroupDescriptionOverrides)
      }
    },
    svg: {
      raw: svgContent,
      source: hashText(svgContent),
      canonical: hashText(canonicalSvg),
      semantic: hashText(JSON.stringify(semanticRows)),
      visual: {
        width: rasterWidth,
        height: rasterHeight,
        sha256: [...rasterDigest]
          .map((value) => value.toString(16).padStart(2, '0'))
          .join('')
      },
      name: String(app.results?.[app.selectedResultIndex]?.name || '')
    },
    featureCatalog: {
      sha256: [...catalogDigest]
        .map((value) => value.toString(16).padStart(2, '0'))
        .join(''),
      jsonBytes: catalogBytes.byteLength,
      schema: Number(state.featureCatalog.value?.schema || 0),
      itemCount: state.featureCatalog.value?.items?.length || 0
    },
    history: {
      undoCount: history.getUndoCount(),
      redoCount: history.getRedoCount(),
      currentCheckpointAbsent: history.getCurrentCheckpoint() === null,
      artifactCheckpointBuilds: Number(historyDiagnostics.artifactCheckpointBuilds || 0),
      artifactCheckpointSignatureComputations: Number(
        historyDiagnostics.artifactCheckpointSignatureComputations || 0
      ),
      generatedArtifactFullCloneCount: Number(
        historyDiagnostics.generatedArtifactFullCloneCount || 0
      ),
      generatedArtifactFullSerializationCount: Number(
        historyDiagnostics.generatedArtifactFullSerializationCount || 0
      )
    }
  };
});

const runGenerate = async (page) => {
  const result = await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    const outcome = await app.runAnalysis();
    return {
      outcome,
      errorSummary: String(app.errorLog?.summary || ''),
      resultCount: app.results?.length || 0
    };
  });
  expect(result).toMatchObject({
    outcome: { status: 'ok' },
    errorSummary: ''
  });
  expect(result.resultCount).toBeGreaterThan(0);
  return capturePageEvidence(page);
};

const svgEquivalence = (left, right) => ({
  exact: left.source.length === right.source.length && left.source.hash === right.source.hash,
  canonical:
    left.canonical.length === right.canonical.length &&
    left.canonical.hash === right.canonical.hash,
  semantic:
    left.semantic.length === right.semantic.length &&
    left.semantic.hash === right.semantic.hash,
  visual:
    left.visual.width === right.visual.width &&
    left.visual.height === right.visual.height &&
    left.visual.sha256 === right.visual.sha256
});

const expectSvgEquivalent = (left, right, label) => {
  const comparison = svgEquivalence(left, right);
  svgComparisonIndex += 1;
  const leftPath = join(svgComparisonRoot, `${svgComparisonIndex}-left.svg`);
  const rightPath = join(svgComparisonRoot, `${svgComparisonIndex}-right.svg`);
  writeFileSync(leftPath, left.raw, 'utf8');
  writeFileSync(rightPath, right.raw, 'utf8');
  const command = [
    'import sys',
    'from tests.utils.svg_compare import compare_svgs',
    'result = compare_svgs(sys.argv[1], sys.argv[2])',
    'print(result.message)',
    'print("\\n".join(result.differences))',
    'raise SystemExit(0 if result.equal else 1)'
  ].join(';');
  const semanticComparison = spawnSync(
    process.env.GBDRAW_PYTHON || 'python',
    ['-c', command, leftPath, rightPath],
    { cwd: repoRoot, encoding: 'utf8' }
  );
  expect(
    semanticComparison.status,
    `${label}: ${semanticComparison.stdout}\n${semanticComparison.stderr}`
  ).toBe(0);
  return comparison;
};

const applyDivergentDraft = async (page) => page.evaluate(async () => {
  const { state } = await import('/gbdraw/web/js/state.js');
  const { createAnnotationSet } = await import('/gbdraw/web/js/app/annotations/state.js');
  const app = window.__GBDRAW_APP__;
  const history = window.__GBDRAW_HISTORY__;
  state.autoLabelReflowEnabled.value = false;

  const orange = state.paletteDefinitions.value.orange || state.currentColors.value;
  const draftColors = state.normalizePaletteColors({
    ...orange,
    CDS: '#0b4f6c',
    tRNA: '#f59e0b'
  });
  state.selectedPalette.value = 'orange';
  state.currentColors.value = draftColors;
  state.paletteInstantPreviewEnabled.value = false;
  state.pendingPaletteName.value = 'orange';
  state.pendingPaletteColors.value = { ...draftColors };

  state.manualSpecificRules.splice(0, state.manualSpecificRules.length, {
    feat: 'CDS',
    qual: 'gene',
    val: '^ND2$',
    color: '#dc2626',
    cap: 'Saved ND2 rule',
    fromFile: false
  });
  state.manualPriorityRules.splice(0, state.manualPriorityRules.length, {
    feat: 'CDS',
    order: 'product,gene'
  });
  state.filterMode.value = 'Whitelist';
  state.manualWhitelist.splice(0, state.manualWhitelist.length, {
    feat: 'CDS',
    qual: 'gene',
    key: 'ND1'
  });
  state.manualBlacklist.value = 'hypothetical, draft-only';
  Object.assign(state.form, {
    plot_title: 'Saved active draft',
    labels_mode: 'both',
    show_scale: false
  });
  Object.assign(state.adv, {
    axis_stroke_width: 6,
    label_font_size: 27,
    feature_width_circular: 18
  });
  state.annotationSets.splice(0, state.annotationSets.length, createAnnotationSet({
    id: 'saved-regenerate-annotations',
    legendLabel: 'Saved annotation',
    annotations: [{
      id: 'saved-window',
      target: {
        kind: 'coordinateSpan',
        record: null,
        start: 100,
        end: 300,
        coordinateSpace: 'source'
      },
      label: 'Saved window',
      mark: 'band'
    }]
  }));

  const feature = app.extractedFeatures.find((candidate) => (
    candidate.type === 'CDS' && app.canEditFeatureColor(candidate)
  )) || app.extractedFeatures.find((candidate) => app.canEditFeatureColor(candidate));
  if (!feature) throw new Error('No editable Feature is available for the round-trip contract.');
  const featureId = String(feature.svg_id || feature.id || '');
  const fillApplied = await app.setFeatureColorValue(
    feature,
    '#7c3aed',
    'Saved direct fill'
  );
  app.openFeatureEditorFromList(feature, { clientX: 200, clientY: 200 });
  const strokeApplied = await app.updateClickedFeatureStroke('#111827', 4);

  const labelFeature = app.extractedFeatures.find((candidate) => (
    app.getEditableLabelByFeatureId(String(candidate.svg_id || candidate.id || ''))
  ));
  if (!labelFeature) throw new Error('No editable Feature label is available for the contract.');
  const labelFeatureId = String(labelFeature.svg_id || labelFeature.id || '');
  const labelRequested = await app.requestLabelTextChangeByFeatureId(
    labelFeatureId,
    'Saved direct label'
  );
  await app.handleLabelTextScopeChoice('single');
  app.openFeatureEditorFromList(labelFeature, { clientX: 220, clientY: 220 });
  if (!app.clickedFeature?.hasEditableLabel) {
    throw new Error('The selected Feature label could not be opened for direct editing.');
  }
  app.clickedFeature.labelText = 'Saved direct label';
  app.clickedFeature.labelVisibility = 'off';
  await app.updateClickedFeatureLabelText();
  const labelVisibilityApplied = state.labelVisibilityOverrides[labelFeatureId] === 'off';
  const visibilityApplied = await app.setFeatureVisibility(
    feature,
    'off',
    { triggerReflow: false, scope: { id: 'feature' } }
  );
  state.legendColorOverrides['Saved direct fill'] = '#7c3aed';
  state.legendStrokeOverrides['Saved direct fill'] = {
    strokeColor: '#111827',
    strokeWidth: 4
  };
  state.clickedFeature.value = null;

  await new Promise((resolveFrame) => (
    requestAnimationFrame(() => requestAnimationFrame(resolveFrame))
  ));
  const fillOverrideKey = Object.keys(state.featureColorOverrides).find((key) => (
    state.featureColorOverrides[key]?.color === '#7c3aed'
  ));
  const strokeOverrideKey = Object.keys(state.featureStrokeOverrides).find((key) => (
    state.featureStrokeOverrides[key]?.strokeColor === '#111827'
  ));
  return {
    featureId,
    labelFeatureId,
    fillOverrideKey,
    strokeOverrideKey,
    fillApplied,
    strokeApplied,
    labelRequested,
    labelVisibilityApplied,
    visibilityApplied,
    undoCount: history.getUndoCount()
  };
});

test.describe.configure({ mode: 'serial' });

test('no-draft session preserves preview, active config, canonical request, and catalog', async ({
  page,
  context
}, testInfo) => {
  test.setTimeout(420_000);
  await openIntentApp(page);
  await loadCurrentSession(page, sourceSessionPath, sourceSession);
  const initiallyLoaded = await capturePageEvidence(page);

  const generatedA = await runGenerate(page);
  const sourceFixtureToSeed = svgEquivalence(initiallyLoaded.svg, generatedA.svg);
  const firstSave = await saveCurrentSession(page, 'intent-no-draft-first');
  const firstCanonical = sessionCanonicalEvidence(firstSave.session);

  const freshPage = await context.newPage();
  await openIntentApp(freshPage);
  await loadCurrentSession(freshPage, firstSave.path, firstSave.session);
  const reloadedA = await capturePageEvidence(freshPage);
  expect(reloadedA.active).toEqual(generatedA.active);
  expect(reloadedA.featureCatalog).toEqual(generatedA.featureCatalog);
  const savedToLoaded = expectSvgEquivalent(
    generatedA.svg,
    reloadedA.svg,
    'saved generated preview versus fresh load'
  );

  const generatedAPrime = await runGenerate(freshPage);
  const firstToRegenerated = expectSvgEquivalent(
    generatedA.svg,
    generatedAPrime.svg,
    'no-draft Generate A versus A prime'
  );
  const secondSave = await saveCurrentSession(freshPage, 'intent-no-draft-second');
  const secondCanonical = sessionCanonicalEvidence(secondSave.session);
  expect(secondCanonical).toEqual(firstCanonical);
  expect(generatedAPrime.featureCatalog).toEqual(generatedA.featureCatalog);

  await testInfo.attach('session-regenerate-no-draft.json', {
    body: Buffer.from(JSON.stringify({
      sourceFixtureToSeed,
      savedToLoaded,
      firstToRegenerated,
      firstCanonical,
      secondCanonical,
      activeIntentSha256: jsonDigest(generatedA.active),
      featureCatalog: generatedA.featureCatalog
    }, null, 2)),
    contentType: 'application/json'
  });
  await freshPage.close();
});

test('loaded current preview supports direct edits before the first Generate', async ({
  page,
  context
}, testInfo) => {
  test.setTimeout(600_000);
  const idleWorkerStages = [];
  const featureBlocks = (snapshot, owner = 'dom') => snapshot[owner].feature.filter(
    (element) => element.part !== 'connector' && !element.id.includes('__line')
  );
  const expectFeatureFill = (snapshot, color) => {
    const domBlocks = featureBlocks(snapshot);
    const resultBlocks = featureBlocks(snapshot, 'resultDom');
    expect(domBlocks.length, JSON.stringify(snapshot.identity)).toBeGreaterThan(0);
    expect(resultBlocks.length).toBeGreaterThan(0);
    expect(domBlocks.every((element) => element.fill === color)).toBe(true);
    expect(resultBlocks.every((element) => element.fill === color)).toBe(true);
  };
  const expectFeatureStroke = (snapshot) => {
    expect(snapshot.dom.feature.length).toBeGreaterThan(0);
    expect(snapshot.resultDom.feature.length).toBeGreaterThan(0);
    expect(snapshot.dom.feature.every((element) => (
      element.stroke === DIRECT_STROKE
      && Number(element.strokeWidth) === DIRECT_STROKE_WIDTH
    ))).toBe(true);
    expect(snapshot.resultDom.feature.every((element) => (
      element.stroke === DIRECT_STROKE
      && Number(element.strokeWidth) === DIRECT_STROKE_WIDTH
    ))).toBe(true);
  };
  const expectFeatureHidden = (snapshot, { generated = false } = {}) => {
    const domFeatures = snapshot.dom.visibilityFeature;
    const resultFeatures = snapshot.resultDom.visibilityFeature;
    if (!generated) {
      expect(domFeatures.length).toBeGreaterThan(0);
      expect(resultFeatures.length).toBeGreaterThan(0);
    }
    expect(domFeatures.every((element) => element.display === 'none')).toBe(true);
    expect(resultFeatures.every((element) => element.display === 'none')).toBe(true);
  };
  const expectLabelEdited = (snapshot, { generated = false } = {}) => {
    expect(
      snapshot.dom.label,
      JSON.stringify({ identity: snapshot.identity, overrides: snapshot.overrides })
    ).toEqual({ text: DIRECT_LABEL_TEXT, display: null });
    expect(snapshot.resultDom.label).toEqual({ text: DIRECT_LABEL_TEXT, display: null });
    if (!generated) {
      expect(snapshot.dom.labelVisibility?.display).toBe('none');
      expect(snapshot.resultDom.labelVisibility?.display).toBe('none');
    } else {
      expect(
        snapshot.dom.labelVisibility === null ||
          snapshot.dom.labelVisibility.display === 'none'
      ).toBe(true);
      expect(
        snapshot.resultDom.labelVisibility === null ||
          snapshot.resultDom.labelVisibility.display === 'none'
      ).toBe(true);
    }
  };
  const expectLegendEdited = (snapshot) => {
    expect(snapshot.dom.legend?.fill).toBe(DIRECT_LEGEND_COLOR);
    expect(snapshot.resultDom.legend?.fill).toBe(DIRECT_LEGEND_COLOR);
    expect(snapshot.dom.legend).toMatchObject({
      stroke: DIRECT_LEGEND_STROKE,
      strokeWidth: String(DIRECT_LEGEND_STROKE_WIDTH)
    });
    expect(snapshot.resultDom.legend).toMatchObject({
      stroke: DIRECT_LEGEND_STROKE,
      strokeWidth: String(DIRECT_LEGEND_STROKE_WIDTH)
    });
  };
  const expectAllDirectEdits = (snapshot, options = {}) => {
    expectFeatureFill(snapshot, DIRECT_FILL);
    expectFeatureStroke(snapshot);
    expectFeatureHidden(snapshot, options);
    expectLabelEdited(snapshot, options);
    expectLegendEdited(snapshot);
    expect(snapshot.overrides).toMatchObject({
      fill: { color: DIRECT_FILL },
      stroke: {
        strokeColor: DIRECT_STROKE,
        strokeWidth: DIRECT_STROKE_WIDTH
      },
      visibility: 'off',
      labelText: DIRECT_LABEL_TEXT,
      labelVisibility: 'off',
      legendColor: DIRECT_LEGEND_COLOR,
      legendStroke: {
        strokeColor: DIRECT_LEGEND_STROKE,
        strokeWidth: DIRECT_LEGEND_STROKE_WIDTH
      }
    });
  };

  await openIntentApp(page);
  await loadCurrentSession(page, sourceSessionPath, sourceSession);
  await captureIdleWorkerStage(page, 'after Load', idleWorkerStages);
  const target = await prepareLoadedPreviewDirectEditTarget(page);
  expect(target.featureRecordKey).toBeTruthy();
  expect(target.biologicalFeatureId).toBeTruthy();
  expect(target.featureOverrideKey).toBe(
    `${target.featureRecordKey}\u0000${target.biologicalFeatureId}`
  );
  expect(target.visibilityFeatureOverrideKey).toBeTruthy();
  expect(target.visibilityFeatureId).not.toBe(target.featureId);
  await captureIdleWorkerStage(page, 'before direct edits', idleWorkerStages);

  const featureSelector = [
    `[data-gbdraw-feature-id="${target.featureId}"]`,
    `[data-gbdraw-rendered-feature-id="${target.featureId}"]`,
    `#${target.featureId}`,
    `[id^="${target.featureId}__"]`
  ].join(', ');
  const featureLocator = page.locator(featureSelector).first();
  await featureLocator.scrollIntoViewIfNeeded();
  const featureClickPoint = await featureLocator.evaluate((element, featureId) => {
    const svg = element.closest('svg');
    const featureSelector = [
      'path[data-gbdraw-feature-id]',
      'polygon[data-gbdraw-feature-id]',
      'rect[data-gbdraw-feature-id]',
      'path[id^="f"]',
      'polygon[id^="f"]',
      'rect[id^="f"]'
    ].join(', ');
    const identity = (candidate) => String(
      candidate?.getAttribute?.('data-gbdraw-rendered-feature-id')
      || candidate?.getAttribute?.('data-gbdraw-feature-id')
      || candidate?.id
      || ''
    ).replace(/__(?:part|line)\d+$/, '');
    const pointIsTarget = (x, y) => {
      if (x < 0 || y < 0 || x > innerWidth || y > innerHeight) return false;
      const topFeature = document.elementsFromPoint(x, y)
        .map((candidate) => (
          candidate.matches?.(featureSelector)
            ? candidate
            : candidate.closest?.(featureSelector)
        ))
        .find((candidate) => candidate && svg.contains(candidate));
      return identity(topFeature) === featureId;
    };
    const screenPoint = (point) => {
      const matrix = element.getScreenCTM();
      if (!matrix) return null;
      const transformed = new DOMPoint(point.x, point.y).matrixTransform(matrix);
      return { x: transformed.x, y: transformed.y };
    };
    const candidates = [];
    if (typeof element.getTotalLength === 'function') {
      const length = element.getTotalLength();
      for (let index = 1; index < 40; index += 1) {
        candidates.push(screenPoint(element.getPointAtLength(length * index / 40)));
      }
    }
    const bounds = element.getBBox?.();
    if (bounds) {
      for (const xFactor of [0.25, 0.5, 0.75]) {
        for (const yFactor of [0.25, 0.5, 0.75]) {
          candidates.push(screenPoint({
            x: bounds.x + bounds.width * xFactor,
            y: bounds.y + bounds.height * yFactor
          }));
        }
      }
    }
    for (const candidate of candidates.filter(Boolean)) {
      for (const dx of [0, -1, 1]) {
        for (const dy of [0, -1, 1]) {
          const x = candidate.x + dx;
          const y = candidate.y + dy;
          if (pointIsTarget(x, y)) return { x, y };
        }
      }
    }
    throw new Error(`No painted point is available for loaded Feature ${featureId}.`);
  }, target.featureId);
  await page.mouse.click(featureClickPoint.x, featureClickPoint.y);
  const featureClickDiagnostics = await page.evaluate((featureId) => {
    const svg = window.__GBDRAW_APP__.svgContainer?.querySelector?.('svg');
    const escaped = CSS.escape(featureId);
    const elements = svg
      ? Array.from(svg.querySelectorAll([
          `[data-gbdraw-feature-id="${escaped}"]`,
          `[data-gbdraw-rendered-feature-id="${escaped}"]`,
          `#${escaped}`,
          `[id^="${escaped}__"]`
        ].join(', ')))
      : [];
    return {
      clickedFeatureId: String(window.__GBDRAW_APP__.clickedFeature?.svg_id || ''),
      elements: elements.map((element) => ({
        id: String(element.id || ''),
        cursor: String(element.style?.cursor || ''),
        display: element.getAttribute('display'),
        pointerEvents: getComputedStyle(element).pointerEvents
      }))
    };
  }, target.featureId);
  const featureDialog = page.getByRole('dialog', { name: /Feature details:/ });
  await expect(featureDialog, JSON.stringify(featureClickDiagnostics, null, 2)).toBeVisible();
  Object.assign(target, await bindLoadedPreviewDirectEditLabel(page));
  await page.evaluate(() => window.__GBDRAW_SESSION_INTENT_PROBE__.reset());
  await startLoadedPreviewDirectEditProbe(page);
  const initiallyLoaded = await captureLoadedPreviewDirectEditState(page);
  expect(initiallyLoaded.autoLabelReflowEnabled).toBe(false);
  expect(initiallyLoaded.structural).toMatchObject({
    directEditResultFlushCount: 0,
    directEditSvgSerializationCount: 0
  });
  expect(initiallyLoaded.history).toMatchObject({
    artifactReplacementHistoryEntryCount: 0,
    generatedArtifactFullCloneCount: 0,
    generatedArtifactFullSerializationCount: 0
  });
  expectBiologicalOwnersUnchanged(initiallyLoaded);
  await featureDialog.getByLabel('Feature fill color').first().fill(DIRECT_FILL);
  const colorScopeHeading = page.getByRole('heading', { name: 'Color Change Scope' });
  await expect(colorScopeHeading).toBeVisible();
  const colorScopeDialog = colorScopeHeading.locator('..');
  const captionScopeButton = colorScopeDialog.getByRole('button').filter({
    hasText: `Apply to all "${target.fillLegendCaption}"`
  }).last();
  await expect(captionScopeButton).toBeVisible();
  await captionScopeButton.click();
  await page.waitForFunction(({ overrideKey, color }) => {
    const app = window.__GBDRAW_APP__;
    return app.featureColorOverrides?.[overrideKey]?.color === color
      && String(app.results?.[app.selectedResultIndex]?.content || '').includes(color);
  }, { overrideKey: target.featureOverrideKey, color: DIRECT_FILL });
  await settleMountedDom(page);
  const afterFill = await captureLoadedPreviewDirectEditState(page);
  expectDirectEditFlushed(initiallyLoaded, afterFill);
  expectFeatureFill(afterFill, DIRECT_FILL);
  expect(afterFill.overrides.fill).toMatchObject({ color: DIRECT_FILL });
  expect(afterFill.history.undoCount).toBeGreaterThan(initiallyLoaded.history.undoCount);
  await captureIdleWorkerStage(page, 'after feature fill', idleWorkerStages);

  expect(await page.evaluate(() => window.__GBDRAW_HISTORY__.undo())).toBe(true);
  await settleMountedDom(page);
  const afterFillUndo = await captureLoadedPreviewDirectEditState(page);
  expectDirectEditFlushed(afterFill, afterFillUndo);
  expectFeatureFill(afterFillUndo, featureBlocks(initiallyLoaded)[0].fill);
  expect(afterFillUndo.overrides.fill).toBeNull();
  await captureIdleWorkerStage(page, 'after direct-edit Undo', idleWorkerStages);

  expect(await page.evaluate(() => window.__GBDRAW_HISTORY__.redo())).toBe(true);
  await settleMountedDom(page);
  const afterFillRedo = await captureLoadedPreviewDirectEditState(page);
  expectDirectEditFlushed(afterFillUndo, afterFillRedo);
  expectFeatureFill(afterFillRedo, DIRECT_FILL);
  expect(afterFillRedo.overrides.fill).toMatchObject({ color: DIRECT_FILL });
  await captureIdleWorkerStage(page, 'after direct-edit Redo', idleWorkerStages);

  const strokeApplied = await page.evaluate(async ({ color, width }) => {
    const app = window.__GBDRAW_APP__;
    const targetState = window.__GBDRAW_LOADED_PREVIEW_DIRECT_EDIT_TARGET__;
    app.openFeatureEditorFromList(targetState.feature, { clientX: 220, clientY: 220 });
    return app.updateClickedFeatureStroke(color, width);
  }, { color: DIRECT_STROKE, width: DIRECT_STROKE_WIDTH });
  expect(strokeApplied).toBe(true);
  await settleMountedDom(page);
  const afterStroke = await captureLoadedPreviewDirectEditState(page);
  expectDirectEditFlushed(afterFillRedo, afterStroke);
  expectFeatureStroke(afterStroke);
  expect(afterStroke.overrides.stroke).toMatchObject({
    strokeColor: DIRECT_STROKE,
    strokeWidth: DIRECT_STROKE_WIDTH
  });
  await captureIdleWorkerStage(page, 'after feature stroke', idleWorkerStages);

  const visibilityApplied = await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    const targetState = window.__GBDRAW_LOADED_PREVIEW_DIRECT_EDIT_TARGET__;
    return app.setFeatureVisibility(targetState.visibilityFeature, 'off', {
      triggerReflow: false,
      scope: { id: 'feature' }
    });
  });
  expect(visibilityApplied).toBe(true);
  await settleMountedDom(page);
  const afterVisibility = await captureLoadedPreviewDirectEditState(page);
  expectDirectEditFlushed(afterStroke, afterVisibility);
  expectFeatureHidden(afterVisibility);
  expect(afterVisibility.overrides.visibility).toBe('off');
  await captureIdleWorkerStage(page, 'after feature visibility', idleWorkerStages);

  const labelApplied = await page.evaluate(async ({ text }) => {
    const app = window.__GBDRAW_APP__;
    const targetState = window.__GBDRAW_LOADED_PREVIEW_DIRECT_EDIT_TARGET__;
    const requested = await app.requestLabelTextChangeByFeatureId(
      targetState.labelFeatureId,
      text
    );
    await app.handleLabelTextScopeChoice('single');
    app.openFeatureEditorFromList(targetState.labelFeature, { clientX: 240, clientY: 240 });
    if (!app.clickedFeature?.hasEditableLabel) {
      throw new Error('The loaded Feature label could not be reopened for direct editing.');
    }
    app.clickedFeature.labelText = text;
    await app.updateClickedFeatureLabelText();
    app.openFeatureEditorFromList(
      targetState.labelVisibilityFeature,
      { clientX: 260, clientY: 260 }
    );
    if (!app.clickedFeature?.hasEditableLabel) {
      throw new Error('The second loaded Feature label could not be opened.');
    }
    app.clickedFeature.labelVisibility = 'off';
    await app.updateClickedFeatureLabelText();
    return {
      requested,
      textOverride: app.labelTextFeatureOverrides[targetState.labelFeatureId],
      visibilityOverride:
        app.labelVisibilityOverrides[targetState.labelVisibilityFeatureId]
    };
  }, { text: DIRECT_LABEL_TEXT });
  expect(labelApplied).toEqual({
    requested: true,
    textOverride: DIRECT_LABEL_TEXT,
    visibilityOverride: 'off'
  });
  await settleMountedDom(page);
  const afterLabel = await captureLoadedPreviewDirectEditState(page);
  expectDirectEditFlushed(afterVisibility, afterLabel);
  expectLabelEdited(afterLabel);
  expect(afterLabel.overrides).toMatchObject({
    labelText: DIRECT_LABEL_TEXT,
    labelVisibility: 'off'
  });
  await captureIdleWorkerStage(page, 'after label text and visibility', idleWorkerStages);

  const legendApplied = await page.evaluate(async ({ color, stroke, strokeWidth }) => {
    const app = window.__GBDRAW_APP__;
    const targetState = window.__GBDRAW_LOADED_PREVIEW_DIRECT_EDIT_TARGET__;
    const index = app.legendEntries.findIndex(
      (entry) => entry.caption === targetState.legendCaption
    );
    if (index < 0) throw new Error('The loaded legend entry is no longer available.');
    return {
      color: app.updateLegendEntryColor(index, color),
      stroke: await app.setLegendEntryStrokeColorValue(index, stroke),
      strokeWidth: app.updateLegendEntryStrokeWidth(index, strokeWidth)
    };
  }, {
    color: DIRECT_LEGEND_COLOR,
    stroke: DIRECT_LEGEND_STROKE,
    strokeWidth: DIRECT_LEGEND_STROKE_WIDTH
  });
  expect(legendApplied).toEqual({ color: true, stroke: true, strokeWidth: true });
  await settleMountedDom(page);
  const afterLegend = await captureLoadedPreviewDirectEditState(page);
  expectDirectEditFlushed(afterLabel, afterLegend);
  expectLegendEdited(afterLegend);
  expect(afterLegend.overrides.legendColor).toBe(DIRECT_LEGEND_COLOR);
  expect(afterLegend.overrides.legendStroke).toMatchObject({
    strokeColor: DIRECT_LEGEND_STROKE,
    strokeWidth: DIRECT_LEGEND_STROKE_WIDTH
  });
  await captureIdleWorkerStage(page, 'after legend color', idleWorkerStages);

  const titleSummary = page.locator('summary[aria-label="Title & Legend"]');
  const titleDetails = titleSummary.locator('..');
  if ((await titleDetails.getAttribute('open')) === null) await titleSummary.click();
  await page.getByLabel('Plot Title', { exact: true }).fill(DIRECT_REGENERATED_TITLE);
  await page.getByLabel('Plot Title Position').selectOption('top');
  await settleMountedDom(page);
  const afterDraftConfig = await captureLoadedPreviewDirectEditState(page);
  expect(afterDraftConfig.result.sha256).toBe(afterLegend.result.sha256);
  expectBiologicalOwnersUnchanged(afterDraftConfig);
  await captureIdleWorkerStage(page, 'after divergent active config', idleWorkerStages);

  const saved = await saveCurrentSession(page, 'loaded-preview-direct-edit');
  const afterSave = await captureLoadedPreviewDirectEditState(page);
  expect(afterSave.result.sha256).toBe(afterLegend.result.sha256);
  expectAllDirectEdits(afterSave);
  expectBiologicalOwnersUnchanged(afterSave);
  await captureIdleWorkerStage(page, 'after Save without Generate', idleWorkerStages);
  expect(
    saved.session.results[saved.session.ui.selectedResultIndex].content
  ).toBe(await page.evaluate(() => (
    window.__GBDRAW_APP__.results[window.__GBDRAW_APP__.selectedResultIndex].content
  )));
  expect(saved.session.config).toMatchObject({
    form: { plot_title: DIRECT_REGENERATED_TITLE }
  });
  expect(saved.session.ui).toMatchObject({
    layoutPreferences: {
      circular: { single: { plotTitlePosition: 'top' } }
    }
  });
  expect(saved.session.features).toMatchObject({
    featureColorOverrides: {
      [target.featureOverrideKey]: { color: DIRECT_FILL }
    },
    featureVisibilityOverrides: { [target.visibilityFeatureId]: 'off' },
    labelTextFeatureOverrides: { [target.labelFeatureId]: DIRECT_LABEL_TEXT },
    labelVisibilityOverrides: { [target.labelVisibilityFeatureId]: 'off' }
  });
  expect(saved.session.editorState).toMatchObject({
    legend: {
      colorOverrides: { [target.legendCaption]: DIRECT_LEGEND_COLOR },
      strokeOverrides: {
        [target.legendCaption]: {
          strokeColor: DIRECT_LEGEND_STROKE,
          strokeWidth: DIRECT_LEGEND_STROKE_WIDTH
        }
      }
    },
    featureStrokes: {
      overrides: {
        [target.featureOverrideKey]: {
          strokeColor: DIRECT_STROKE,
          strokeWidth: DIRECT_STROKE_WIDTH
        }
      }
    }
  });
  const directMetrics = await page.evaluate(() => (
    window.__GBDRAW_LOADED_PREVIEW_DIRECT_EDIT_PROBE__.stop()
  ));
  expect(directMetrics.directEditResultFlushCount).toBeGreaterThanOrEqual(7);
  expect(directMetrics.directEditSvgSerializationCount).toBeGreaterThanOrEqual(
    directMetrics.directEditResultFlushCount
  );
  expect(directMetrics).toMatchObject({
    directEditFeatureCatalogReplacementCount: 0,
    directEditExtractedFeatureReplacementCount: 0,
    directEditBiologicalFeatureReplacementCount: 0,
    directEditOrthogroupReplacementCount: 0
  });

  const freshPage = await context.newPage();
  await openIntentApp(freshPage);
  await loadCurrentSession(freshPage, saved.path, saved.session);
  await captureIdleWorkerStage(freshPage, 'after fresh Load', idleWorkerStages);
  const freshTarget = await prepareLoadedPreviewDirectEditTarget(freshPage, target);
  Object.assign(
    freshTarget,
    await bindLoadedPreviewDirectEditLabel(freshPage, target)
  );
  expect(freshTarget).toEqual(target);
  await startLoadedPreviewDirectEditProbe(freshPage);
  const freshlyLoaded = await captureLoadedPreviewDirectEditState(freshPage);
  expectAllDirectEdits(freshlyLoaded);
  expect(freshlyLoaded.result.equalsLoadedResult).toBe(true);
  expect(freshlyLoaded.history).toMatchObject({
    undoCount: 0,
    redoCount: 0,
    currentCheckpointAbsent: true,
    artifactReplacementHistoryEntryCount: 0,
    generatedArtifactFullCloneCount: 0,
    generatedArtifactFullSerializationCount: 0
  });
  expect(freshlyLoaded.autoLabelReflowEnabled).toBe(false);
  expectBiologicalOwnersUnchanged(freshlyLoaded);
  await captureIdleWorkerStage(
    freshPage,
    'before first Generate on fresh Load',
    idleWorkerStages
  );
  const freshLoadMetrics = await freshPage.evaluate(() => (
    window.__GBDRAW_LOADED_PREVIEW_DIRECT_EDIT_PROBE__.stop()
  ));
  expect(freshLoadMetrics).toMatchObject({
    directEditResultFlushCount: 0,
    directEditSvgSerializationCount: 0,
    directEditFeatureCatalogReplacementCount: 0,
    directEditExtractedFeatureReplacementCount: 0,
    directEditBiologicalFeatureReplacementCount: 0,
    directEditOrthogroupReplacementCount: 0
  });

  await runGenerate(freshPage);
  await settleMountedDom(freshPage);
  const regenerated = await captureLoadedPreviewDirectEditState(freshPage);
  expect(regenerated.result.sha256).not.toBe(freshlyLoaded.result.sha256);
  expectAllDirectEdits(regenerated, { generated: true });
  expect(regenerated.history).toMatchObject({
    undoCount: 1,
    redoCount: 0,
    currentCheckpointAbsent: true,
    artifactReplacementHistoryEntryCount: 1,
    generatedArtifactFullCloneCount: 0,
    generatedArtifactFullSerializationCount: 0
  });
  expect(regenerated.history.artifactCheckpointBuilds).toBe(
    freshlyLoaded.history.artifactCheckpointBuilds
  );
  const generatedWorker = await getDiagramWorkerActivity(freshPage);
  expect({
    constructions: generatedWorker.constructions,
    initializations: generatedWorker.initializations,
    helpers: generatedWorker.helpers,
    runs: generatedWorker.runs
  }).toEqual({ constructions: 1, initializations: 1, helpers: 0, runs: 1 });

  expect(await freshPage.evaluate(() => window.__GBDRAW_HISTORY__.undo())).toBe(true);
  await settleMountedDom(freshPage);
  const generatedUndo = await captureLoadedPreviewDirectEditState(freshPage);
  expect(generatedUndo.result.equalsLoadedResult).toBe(true);
  expectAllDirectEdits(generatedUndo);
  expect(generatedUndo.history).toMatchObject({
    undoCount: 0,
    redoCount: 1,
    artifactReplacementHistoryEntryCount: 1
  });
  const workerAfterGenerateUndo = await getDiagramWorkerActivity(freshPage);
  expect({
    constructions: workerAfterGenerateUndo.constructions,
    initializations: workerAfterGenerateUndo.initializations,
    helpers: workerAfterGenerateUndo.helpers,
    runs: workerAfterGenerateUndo.runs
  }).toEqual({ constructions: 1, initializations: 1, helpers: 0, runs: 1 });

  expect(await freshPage.evaluate(() => window.__GBDRAW_HISTORY__.redo())).toBe(true);
  await settleMountedDom(freshPage);
  const generatedRedo = await captureLoadedPreviewDirectEditState(freshPage);
  expect(generatedRedo.result.sha256).toBe(regenerated.result.sha256);
  expectAllDirectEdits(generatedRedo, { generated: true });
  expect(generatedRedo.history).toMatchObject({
    undoCount: 1,
    redoCount: 0,
    artifactReplacementHistoryEntryCount: 1
  });
  const workerAfterGenerateRedo = await getDiagramWorkerActivity(freshPage);
  expect({
    constructions: workerAfterGenerateRedo.constructions,
    initializations: workerAfterGenerateRedo.initializations,
    helpers: workerAfterGenerateRedo.helpers,
    runs: workerAfterGenerateRedo.runs
  }).toEqual({ constructions: 1, initializations: 1, helpers: 0, runs: 1 });

  await testInfo.attach('loaded-preview-direct-edit.json', {
    body: Buffer.from(JSON.stringify({
      target,
      actualUiOperation: 'SVG Feature click -> Feature details -> Feature fill color',
      idleWorkerStages,
      directMetrics: {
        ...directMetrics,
        directEditWorkerConstructionDelta: 0,
        directEditWorkerInitializationDelta: 0,
        directEditWorkerRunDelta: 0,
        directEditResourceTransferBytes: 0
      },
      domEvidence: {
        before: initiallyLoaded.dom,
        afterFill: afterFill.dom,
        afterStroke: afterStroke.dom,
        afterVisibility: afterVisibility.dom,
        afterLabel: afterLabel.dom,
        afterLegend: afterLegend.dom,
        freshLoad: freshlyLoaded.dom,
        regenerated: regenerated.dom
      },
      resultEvidence: {
        before: initiallyLoaded.result,
        afterFill: afterFill.result,
        afterStroke: afterStroke.result,
        afterVisibility: afterVisibility.result,
        afterLabel: afterLabel.result,
        afterLegend: afterLegend.result,
        freshLoad: freshlyLoaded.result,
        regenerated: regenerated.result,
        generatedUndo: generatedUndo.result,
        generatedRedo: generatedRedo.result
      },
      workerAfterGenerate: {
        constructions: generatedWorker.constructions,
        initializations: generatedWorker.initializations,
        helpers: generatedWorker.helpers,
        runs: generatedWorker.runs,
        transferredResourceBytes: diagramWorkerResourceBytes(generatedWorker)
      }
    }, null, 2)),
    contentType: 'application/json'
  });
  await freshPage.close();
});

test('divergent draft and direct editor overrides survive repeated Save, Load, and Generate', async ({
  page,
  context
}, testInfo) => {
  test.setTimeout(600_000);
  await openIntentApp(page);
  await loadCurrentSession(page, sourceSessionPath, sourceSession);
  const generatedA = await runGenerate(page);
  const actions = await applyDivergentDraft(page);
  expect(actions).toMatchObject({
    fillApplied: true,
    strokeApplied: true,
    labelRequested: true,
    labelVisibilityApplied: true,
    visibilityApplied: true
  });
  expect(actions.undoCount).toBeGreaterThan(0);

  const savedDraftIntent = await capturePageEvidence(page);
  expect(savedDraftIntent.active.palette).toMatchObject({
    selected: 'orange',
    instantPreview: false,
    pendingName: 'orange'
  });
  expect(savedDraftIntent.active.palette.currentColors).toMatchObject({
    CDS: '#0b4f6c',
    tRNA: '#f59e0b'
  });
  expect(savedDraftIntent.active.specificRules).toEqual(expect.arrayContaining([
    expect.objectContaining({ val: '^ND2$', color: '#dc2626' })
  ]));
  expect(savedDraftIntent.active.annotations.map(({ id }) => id)).toContain(
    'saved-regenerate-annotations'
  );
  expect(savedDraftIntent.active.editorOverrides).toMatchObject({
    fills: { [actions.fillOverrideKey]: expect.objectContaining({ color: '#7c3aed' }) },
    strokes: {
      [actions.strokeOverrideKey]: expect.objectContaining({ strokeColor: '#111827' })
    },
    visibility: { [actions.featureId]: 'off' },
    labelText: { [actions.labelFeatureId]: 'Saved direct label' },
    labelVisibility: { [actions.labelFeatureId]: 'off' }
  });

  const divergentSave = await saveCurrentSession(page, 'intent-divergent-draft');
  expect(divergentSave.session.config.palette).toBe('orange');
  expect(divergentSave.session.config.colors).toMatchObject({
    CDS: '#0b4f6c',
    tRNA: '#f59e0b'
  });
  expect(divergentSave.session.renderRequest.diagramOptions.colors.defaultColorsPalette).toBe(
    sourceSession.renderRequest.diagramOptions.colors.defaultColorsPalette
  );

  const expectedB = await runGenerate(page);
  expect(svgEquivalence(expectedB.svg, generatedA.svg).semantic).toBe(false);
  const expectedBSave = await saveCurrentSession(page, 'intent-expected-b');
  const expectedBCanonical = sessionCanonicalEvidence(expectedBSave.session);

  const freshPage = await context.newPage();
  await openIntentApp(freshPage);
  await loadCurrentSession(freshPage, divergentSave.path, divergentSave.session);
  const loadedDraft = await capturePageEvidence(freshPage);
  expect(loadedDraft.active).toEqual(savedDraftIntent.active);
  expectSvgEquivalent(
    loadedDraft.svg,
    savedDraftIntent.svg,
    'divergent session must retain the saved committed preview'
  );
  expect(svgEquivalence(loadedDraft.svg, expectedB.svg).visual).toBe(false);
  expect(loadedDraft.history).toMatchObject({
    undoCount: 0,
    redoCount: 0,
    currentCheckpointAbsent: true,
    generatedArtifactFullCloneCount: 0,
    generatedArtifactFullSerializationCount: 0
  });
  await freshPage.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    window.__GBDRAW_LOADED_DRAFT_CONTENT__ = String(
      app.results?.[app.selectedResultIndex]?.content || ''
    );
  });

  const beforeGenerateHistory = await freshPage.evaluate(() => (
    window.__GBDRAW_HISTORY__.getDiagnostics()
  ));
  const actualBPrime = await runGenerate(freshPage);
  const expectedToActual = expectSvgEquivalent(
    expectedB.svg,
    actualBPrime.svg,
    'expected B versus actual B prime after loading the divergent draft'
  );
  expect(actualBPrime.active.editorOverrides).toEqual(
    expectedB.active.editorOverrides
  );
  const actualBSave = await saveCurrentSession(freshPage, 'intent-actual-b');
  const actualBCanonical = sessionCanonicalEvidence(actualBSave.session);
  expect(actualBCanonical).toEqual(expectedBCanonical);
  const afterGenerateHistory = await freshPage.evaluate(() => ({
    diagnostics: window.__GBDRAW_HISTORY__.getDiagnostics(),
    undoCount: window.__GBDRAW_HISTORY__.getUndoCount(),
    redoCount: window.__GBDRAW_HISTORY__.getRedoCount()
  }));
  expect(afterGenerateHistory.undoCount).toBe(1);
  expect(afterGenerateHistory.redoCount).toBe(0);
  expect(
    Number(afterGenerateHistory.diagnostics.artifactCheckpointBuilds || 0) -
      Number(beforeGenerateHistory.artifactCheckpointBuilds || 0)
  ).toBe(0);
  expect(
    Number(afterGenerateHistory.diagnostics.generatedArtifactFullCloneCount || 0) -
      Number(beforeGenerateHistory.generatedArtifactFullCloneCount || 0)
  ).toBe(0);

  const undoRedo = await freshPage.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    const history = window.__GBDRAW_HISTORY__;
    const content = () => String(app.results?.[app.selectedResultIndex]?.content || '');
    const generated = content();
    const undo = await history.undo();
    const loaded = content();
    const redo = await history.redo();
    return {
      undo,
      redo,
      loadedRestored: loaded === window.__GBDRAW_LOADED_DRAFT_CONTENT__,
      generatedRestored: content() === generated,
      undoCount: history.getUndoCount(),
      redoCount: history.getRedoCount()
    };
  });
  expect(undoRedo).toMatchObject({
    undo: true,
    redo: true,
    loadedRestored: true,
    generatedRestored: true,
    undoCount: 1,
    redoCount: 0
  });

  const undoCountBeforeWarmGenerate = await freshPage.evaluate(() => (
    window.__GBDRAW_HISTORY__.getUndoCount()
  ));
  const warmBPrime = await runGenerate(freshPage);
  expectSvgEquivalent(actualBPrime.svg, warmBPrime.svg, 'unchanged warm Generate');
  expect(await freshPage.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount())).toBe(
    undoCountBeforeWarmGenerate
  );
  const worker = await getDiagramWorkerActivity(freshPage);
  expect(worker.constructions).toBe(1);
  expect(worker.initializations).toBe(1);
  expect(worker.runs).toBe(2);
  expect(worker.instances[0].runs).toHaveLength(2);
  expect(worker.instances[0].runs[1].stagedResourceBytes).toBeLessThan(1_024);
  expect(worker.instances[0].runs[1].stagedResourceBytes).toBeLessThan(
    worker.instances[0].runs[1].referencedDeclaredBytes
  );
  expect(worker.instances[0].runs[1].transferredBytes).toBe(
    worker.instances[0].runs[1].stagedResourceBytes
  );

  const repeatedSave = await saveCurrentSession(freshPage, 'intent-repeated-b');
  const activeBeforeRepeatedLoad = await capturePageEvidence(freshPage);

  const repeatedPage = await context.newPage();
  await openIntentApp(repeatedPage);
  await loadCurrentSession(repeatedPage, repeatedSave.path, repeatedSave.session);
  const repeatedLoad = await capturePageEvidence(repeatedPage);
  expect(repeatedLoad.active).toEqual(activeBeforeRepeatedLoad.active);
  expect(repeatedLoad.featureCatalog).toEqual(activeBeforeRepeatedLoad.featureCatalog);
  expectSvgEquivalent(
    repeatedLoad.svg,
    activeBeforeRepeatedLoad.svg,
    'second saved preview versus second load'
  );
  const repeatedGenerate = await runGenerate(repeatedPage);
  const repeatedComparison = expectSvgEquivalent(
    actualBPrime.svg,
    repeatedGenerate.svg,
    'repeated Save, Load, and Generate'
  );
  expect(repeatedGenerate.active.editorOverrides).toEqual(
    activeBeforeRepeatedLoad.active.editorOverrides
  );

  await testInfo.attach('session-regenerate-divergent.json', {
    body: Buffer.from(JSON.stringify({
      actions,
      expectedToActual,
      repeatedComparison,
      savedDraftIntentSha256: jsonDigest(savedDraftIntent.active),
      loadedDraftIntentSha256: jsonDigest(loadedDraft.active),
      expectedBCanonical,
      actualBCanonical,
      worker: {
        constructions: worker.constructions,
        initializations: worker.initializations,
        runs: worker.instances[0].runs
      }
    }, null, 2)),
    contentType: 'application/json'
  });
  await freshPage.close();
  await repeatedPage.close();
});

test('bare legacy configuration drives the next canonical request and SVG', async ({ page }) => {
  test.setTimeout(300_000);
  await openIntentApp(page);
  await loadCurrentSession(page, sourceSessionPath, sourceSession);
  const baseline = await capturePageEvidence(page);
  const legacyImport = await page.evaluate(async () => {
    const legacyConfig = {
      form: {
        plot_title: 'Legacy JSON Generate',
        labels_mode: 'both',
        show_scale: false
      },
      adv: {
        axis_stroke_width: 8,
        label_font_size: 25,
        feature_width_circular: 22
      },
      colors: { CDS: '#0b4f6c', tRNA: '#f59e0b' },
      palette: 'orange',
      paletteInstantPreviewEnabled: true,
      rules: [{
        feat: 'CDS',
        qual: 'gene',
        val: '^ND3$',
        color: '#be123c',
        cap: 'Legacy ND3',
        fromFile: false
      }],
      qualifierPriorityRules: [{ feat: 'CDS', order: 'product,gene' }],
      filterMode: 'Blacklist',
      whitelist: [{ feat: 'CDS', qual: 'gene', key: 'ND3' }],
      blacklistText: 'legacy-hidden'
    };
    const result = await window.__GBDRAW_APP__.importSession({
      target: {
        files: [new File(
          [JSON.stringify(legacyConfig)],
          'legacy-regenerate-config.json',
          { type: 'application/json' }
        )],
        value: 'selected'
      }
    });
    return { status: result?.status };
  });
  expect(legacyImport).toEqual({ status: 'legacy' });
  const loadedLegacyIntent = await capturePageEvidence(page);
  expect(loadedLegacyIntent.active).toMatchObject({
    palette: {
      selected: 'orange',
      currentColors: { CDS: '#0b4f6c', tRNA: '#f59e0b' },
      appliedName: 'orange',
      pendingName: ''
    },
    filters: {
      mode: 'Blacklist',
      blacklistText: 'legacy-hidden'
    },
    form: {
      plot_title: 'Legacy JSON Generate',
      labels_mode: 'both',
      show_scale: false
    },
    adv: {
      axis_stroke_width: 8,
      label_font_size: 25,
      feature_width_circular: 22
    }
  });

  const generated = await runGenerate(page);
  expect(svgEquivalence(generated.svg, baseline.svg).visual).toBe(false);
  const saved = await saveCurrentSession(page, 'legacy-config-generated');
  expect(saved.session.config).toMatchObject({
    palette: 'orange',
    colors: { CDS: '#0b4f6c', tRNA: '#f59e0b' },
    filterMode: 'Blacklist',
    blacklistText: 'legacy-hidden',
    form: {
      plot_title: 'Legacy JSON Generate',
      labels_mode: 'both',
      show_scale: false
    },
    adv: {
      axis_stroke_width: 8,
      label_font_size: 25,
      feature_width_circular: 22
    }
  });
  expect(saved.session.config.rules).toEqual([
    expect.objectContaining({ val: '^ND3$', color: '#be123c' })
  ]);
  expect(saved.session.config.qualifierPriorityRules).toEqual([
    { feat: 'CDS', order: 'product,gene' }
  ]);
  const colors = saved.session.renderRequest.diagramOptions.colors;
  const defaultColorRef = colors.defaultColorsFile || colors.defaultColors;
  expect(defaultColorRef?.resourceId).toBeTruthy();
  const defaultColorResource = saved.session.resources[defaultColorRef.resourceId];
  const defaultColorText = Buffer.from(defaultColorResource.data, 'base64').toString('utf8');
  expect(defaultColorText).toContain('CDS\t#0b4f6c');
  expect(defaultColorText).toContain('tRNA\t#f59e0b');
});
