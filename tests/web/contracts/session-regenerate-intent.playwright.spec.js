const { test, expect } = require('@playwright/test');
const { createHash } = require('node:crypto');
const { readFileSync } = require('node:fs');
const { join, resolve } = require('node:path');
const { gunzipSync } = require('node:zlib');
const {
  getDiagramWorkerActivity,
  openApp
} = require('../helpers/app-lifecycle.cjs');

const repoRoot = resolve(process.env.GBDRAW_REPO || process.cwd());
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
    renderRequest: { schema: 5 }
  });
  return { path, session };
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
  expect(
    comparison.exact || comparison.canonical || comparison.semantic || comparison.visual,
    `${label}: ${JSON.stringify({ left, right, comparison }, null, 2)}`
  ).toBe(true);
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
