import assert from 'node:assert/strict';
import test from 'node:test';

import {
  applyLabelVisibility,
  applyLegendColorOverridesToSvg,
  applyStrokeOverridesToSvg,
  createEditorSvgProjection,
  setLabelText
} from '../../gbdraw/web/js/app/editor-svg-projection.js';
import { createHistoryManager } from '../../gbdraw/web/js/services/history.js';
import {
  FakeSvgElement,
  appendEditableLabel,
  appendFeature,
  appendFeatureLegend
} from './helpers/editor-svg-fixture.mjs';

globalThis.CSS = { escape: (value) => String(value) };

const stableKey = 'record-a\u0000biological-a';
const feature = {
  id: 'feature-a',
  recordKey: 'record-a',
  biologicalFeatureId: 'biological-a',
  svg_id: 'rendered-a',
  type: 'CDS',
  product: 'Unrelated feature caption',
  qualifiers: { gene: ['projection-target'] }
};

const buildSvg = ({ labels = 1, includeFeature = true, includeSuppressedLegend = true } = {}) => {
  const svg = new FakeSvgElement('svg');
  const featureElement = includeFeature ? appendFeature(svg, 'rendered-a') : null;
  const labelElements = [];
  for (let index = 0; index < labels; index += 1) {
    labelElements.push(appendEditableLabel(svg, 'rendered-a', 'source label'));
  }
  const legend = appendFeatureLegend(svg, 'Legend A');
  const comparisonLegend = includeSuppressedLegend
    ? svg.appendChild(new FakeSvgElement('g', {
        'data-gbdraw-role': 'comparison-legend',
        'data-gbdraw-orientation': 'h'
      }))
    : null;
  return { comparisonLegend, featureElement, labelElements, legend, svg };
};

const completeProjection = () => createEditorSvgProjection({
  features: [feature],
  featureColorOverrides: {
    [stableKey]: { color: '#abcdef', caption: 'Projected feature' }
  },
  featureStrokeOverrides: {
    [stableKey]: { strokeColor: '#fedcba', strokeWidth: 2 }
  },
  featureVisibilityOverrides: { 'rendered-a': 'off' },
  labelTextFeatureOverrides: { 'rendered-a': 'edited label' },
  labelTextFeatureOverrideSources: { 'rendered-a': 'source label' },
  labelVisibilityOverrides: { 'rendered-a': 'off' },
  legendEntries: [{
    originalCaption: 'Legend A',
    caption: 'Renamed legend',
    color: '#123456'
  }],
  legendColorOverrides: { 'Renamed legend': '#234567' },
  legendStrokeOverrides: {
    'Renamed legend': { strokeColor: '#345678', strokeWidth: 4 }
  },
  suppressPairwiseIdentityLegend: true
});

test('base label visibility provenance transitions remain stable', () => {
  // Independent transition table observed at PR base 3a098d8f in
  // createLabelActions.applyLabelVisibilityPreview(). The expected state is
  // specified here rather than derived from the projection implementation.
  const cases = [
    {
      name: 'renderer-hidden default',
      display: 'none',
      editorMarker: null,
      mode: 'default',
      expectedDisplay: 'none',
      expectedMarker: null,
      expectedChanged: false
    },
    {
      name: 'editor-hidden default',
      display: 'none',
      editorMarker: 'off',
      mode: 'default',
      expectedDisplay: null,
      expectedMarker: null,
      expectedChanged: true
    },
    {
      name: 'explicit on overrides renderer hiding',
      display: 'none',
      editorMarker: null,
      mode: 'on',
      expectedDisplay: null,
      expectedMarker: null,
      expectedChanged: true
    },
    {
      name: 'explicit off owns hiding',
      display: null,
      editorMarker: null,
      mode: 'off',
      markPreview: true,
      expectedDisplay: 'none',
      expectedMarker: 'off',
      expectedChanged: true
    }
  ];

  cases.forEach((entry) => {
    const label = new FakeSvgElement('text');
    label.appendChild(new FakeSvgElement('textPath', {}, 'path label'));
    if (entry.display !== null) label.setAttribute('display', entry.display);
    if (entry.editorMarker !== null) {
      label.setAttribute('data-gbdraw-label-visibility-preview', entry.editorMarker);
    }
    assert.equal(
      applyLabelVisibility(label, entry.mode, { markPreview: Boolean(entry.markPreview) }),
      entry.expectedChanged,
      entry.name
    );
    assert.equal(label.getAttribute('display'), entry.expectedDisplay, entry.name);
    assert.equal(
      label.getAttribute('data-gbdraw-label-visibility-preview'),
      entry.expectedMarker,
      entry.name
    );
  });
});

test('the projection owner applies every persisted editor domain to a fresh SVG', () => {
  const fixture = buildSvg();
  const projection = completeProjection();
  const result = projection.project(fixture.svg);

  assert.equal(result.changed, true);
  assert.equal(fixture.featureElement.getAttribute('fill'), '#abcdef');
  assert.equal(fixture.featureElement.getAttribute('stroke'), '#fedcba');
  assert.equal(fixture.featureElement.getAttribute('stroke-width'), '2');
  assert.equal(fixture.featureElement.getAttribute('display'), 'none');
  assert.equal(fixture.labelElements[0].textContent, 'edited label');
  assert.equal(fixture.labelElements[0].getAttribute('display'), 'none');
  assert.equal(fixture.legend.entry.getAttribute('data-legend-key'), 'Renamed legend');
  assert.equal(fixture.legend.entry.querySelector('text').textContent, 'Renamed legend');
  assert.equal(fixture.legend.swatch.getAttribute('fill'), '#234567');
  assert.equal(fixture.legend.swatch.getAttribute('stroke'), '#345678');
  assert.equal(fixture.legend.swatch.getAttribute('stroke-width'), '4');
  assert.equal(fixture.comparisonLegend.parentElement, null);
  assert.deepEqual(result.unresolvedLabelFeatureIds, []);

  assert.deepEqual(projection.project(fixture.svg), {
    changed: false,
    featureFillCount: 0,
    featureStrokeCount: 0,
    featureVisibilityCount: 0,
    labelTextCount: 0,
    labelVisibilityCount: 0,
    legendCount: 0,
    suppressionCount: 0,
    unresolvedLabelFeatureIds: []
  });
});

const equivalenceCases = [
  {
    name: 'Feature fill',
    input: { features: [feature], featureColorOverrides: { [stableKey]: { color: '#abcdef' } } },
    direct: (fixture) => fixture.featureElement.setAttribute('fill', '#abcdef'),
    read: (fixture) => fixture.featureElement.getAttribute('fill')
  },
  {
    name: 'Feature stroke',
    input: {
      features: [feature],
      featureStrokeOverrides: { [stableKey]: { strokeColor: '#abcdef', strokeWidth: 3 } }
    },
    direct: (fixture) => {
      fixture.featureElement.setAttribute('stroke', '#abcdef');
      fixture.featureElement.setAttribute('stroke-width', '3');
    },
    read: (fixture) => [
      fixture.featureElement.getAttribute('stroke'),
      fixture.featureElement.getAttribute('stroke-width')
    ]
  },
  {
    name: 'Feature visibility',
    input: { features: [feature], featureVisibilityOverrides: { 'rendered-a': 'off' } },
    direct: (fixture) => fixture.featureElement.setAttribute('display', 'none'),
    read: (fixture) => fixture.featureElement.getAttribute('display')
  },
  {
    name: 'Label text',
    input: {
      labelTextFeatureOverrides: { 'rendered-a': 'edited label' },
      labelTextFeatureOverrideSources: { 'rendered-a': 'source label' }
    },
    direct: (fixture) => setLabelText(fixture.labelElements[0], 'edited label'),
    read: (fixture) => fixture.labelElements[0].textContent
  },
  {
    name: 'Label visibility',
    input: { labelVisibilityOverrides: { 'rendered-a': 'off' } },
    direct: (fixture) => applyLabelVisibility(fixture.labelElements[0], 'off'),
    read: (fixture) => fixture.labelElements[0].getAttribute('display')
  },
  {
    name: 'Legend fill',
    input: { legendColorOverrides: { 'Legend A': '#abcdef' } },
    direct: (fixture) => applyLegendColorOverridesToSvg({
      svg: fixture.svg,
      legendColorOverrides: { 'Legend A': '#abcdef' }
    }),
    read: (fixture) => fixture.legend.swatch.getAttribute('fill')
  },
  {
    name: 'Legend stroke',
    input: {
      legendStrokeOverrides: { 'Legend A': { strokeColor: '#abcdef', strokeWidth: 3 } }
    },
    direct: (fixture) => applyStrokeOverridesToSvg({
      svg: fixture.svg,
      legendStrokeOverrides: { 'Legend A': { strokeColor: '#abcdef', strokeWidth: 3 } }
    }),
    read: (fixture) => [
      fixture.legend.swatch.getAttribute('stroke'),
      fixture.legend.swatch.getAttribute('stroke-width')
    ]
  }
];

for (const equivalenceCase of equivalenceCases) {
  test(`${equivalenceCase.name} direct and fresh-SVG projection paths are equivalent`, () => {
    const mounted = buildSvg({ includeSuppressedLegend: false });
    const fresh = buildSvg({ includeSuppressedLegend: false });
    const canonicalBefore = JSON.stringify(equivalenceCase.input);
    equivalenceCase.direct(mounted);
    createEditorSvgProjection(equivalenceCase.input).project(fresh.svg);
    assert.deepEqual(equivalenceCase.read(fresh), equivalenceCase.read(mounted));
    assert.equal(JSON.stringify(equivalenceCase.input), canonicalBefore);
  });
}

test('specific-rule-derived fills use the same Feature identity and fill targets', () => {
  const fixture = buildSvg({ includeSuppressedLegend: false });
  const projection = createEditorSvgProjection({
    features: [feature],
    manualSpecificRules: [{
      feat: 'CDS',
      qual: 'gene',
      val: '^projection-target$',
      color: '#778899',
      cap: 'Rule-derived'
    }]
  });
  assert.equal(projection.project(fixture.svg).changed, true);
  assert.equal(fixture.featureElement.getAttribute('fill'), '#778899');
  assert.deepEqual(projection.featureColorOverrides[stableKey], {
    color: '#778899',
    caption: 'Rule-derived'
  });
});

test('legend rename replay preserves a fresh renderer palette without an explicit color override', () => {
  const fixture = buildSvg({ includeSuppressedLegend: false });
  fixture.legend.swatch.setAttribute('fill', '#f59e0b');
  const projection = createEditorSvgProjection({
    legendEntries: [{
      originalCaption: 'Legend A',
      caption: 'Renamed legend',
      color: '#e8b441'
    }]
  });
  assert.equal(projection.project(fixture.svg).changed, true);
  assert.equal(fixture.legend.entry.getAttribute('data-legend-key'), 'Renamed legend');
  assert.equal(fixture.legend.swatch.getAttribute('fill'), '#f59e0b');
});

test('fresh-SVG projection replays Legend add, delete, rename, duplicate, color, and order state', () => {
  const fixture = buildSvg({ includeSuppressedLegend: false });
  const appendEntry = (caption, color, y) => {
    const entry = fixture.legend.featureLegend.appendChild(
      new FakeSvgElement('g', { 'data-legend-key': caption })
    );
    entry.appendChild(new FakeSvgElement('path', {
      fill: color,
      transform: `translate(0, ${y})`
    }));
    entry.appendChild(new FakeSvgElement('text', {
      transform: `translate(20, ${y})`
    }, caption));
    return entry;
  };
  fixture.legend.swatch.setAttribute('transform', 'translate(0, 0)');
  fixture.legend.entry.querySelector('text').setAttribute('transform', 'translate(20, 0)');
  appendEntry('Legend B', '#222222', 20);
  appendEntry('Legend C', '#333333', 40);

  const projection = createEditorSvgProjection({
    legendEntries: [
      { caption: 'Legend B', originalCaption: 'Legend B', color: '#222222' },
      { caption: 'Renamed A', originalCaption: 'Legend A', color: '#abcdef' },
      { caption: 'Renamed A (1)', originalCaption: 'Renamed A (1)', color: '#fedcba' }
    ],
    deletedLegendEntries: [
      { caption: 'Legend C', originalCaption: 'Legend C', color: '#333333' }
    ],
    originalLegendOrder: ['Renamed A', 'Legend B', 'Legend C'],
    addedLegendCaptions: ['Renamed A (1)'],
    legendColorOverrides: {
      'Renamed A': '#abcdef',
      'Renamed A (1)': '#fedcba'
    }
  });
  assert.equal(projection.project(fixture.svg).changed, true);

  const projected = fixture.legend.featureLegend.children;
  assert.deepEqual(
    projected.map((entry) => entry.getAttribute('data-legend-key')),
    ['Legend B', 'Renamed A', 'Renamed A (1)']
  );
  assert.deepEqual(
    projected.map((entry) => entry.querySelector('text').textContent),
    ['Legend B', 'Renamed A', 'Renamed A (1)']
  );
  assert.deepEqual(
    projected.map((entry) => entry.querySelector('path').getAttribute('fill')),
    ['#222222', '#abcdef', '#fedcba']
  );
  assert.deepEqual(
    projected.map((entry) => entry.querySelector('text').getAttribute('transform')),
    ['translate(20, 0)', 'translate(20, 20)', 'translate(20, 40)']
  );
});

test('a late Legend projection failure rolls back earlier Legend mutations', () => {
  const fixture = buildSvg({ includeSuppressedLegend: false });
  const secondEntry = fixture.legend.featureLegend.appendChild(
    new FakeSvgElement('g', { 'data-legend-key': 'Legend B' })
  );
  const secondSwatch = secondEntry.appendChild(new FakeSvgElement('path', { fill: '#222222' }));
  secondEntry.appendChild(new FakeSvgElement('text', {}, 'Legend B'));
  const originalSetter = secondSwatch.setAttribute.bind(secondSwatch);
  secondSwatch.setAttribute = (name, value) => {
    if (name === 'fill' && value === '#bbbbbb') throw new Error('late Legend mutation failed');
    originalSetter(name, value);
  };
  const projection = createEditorSvgProjection({
    legendColorOverrides: {
      'Legend A': '#aaaaaa',
      'Legend B': '#bbbbbb'
    }
  });

  assert.throws(() => projection.project(fixture.svg), /late Legend mutation failed/);
  assert.equal(fixture.legend.swatch.getAttribute('fill'), '#111111');
  assert.equal(secondSwatch.getAttribute('fill'), '#222222');
});

test('projection borrows the validated Feature owner without cloning or serializing it', () => {
  const metrics = [];
  const previousHooks = globalThis.__GBDRAW_TEST_HOOKS__;
  const previousStructuredClone = globalThis.structuredClone;
  const previousStringify = JSON.stringify;
  let featureCloneCalls = 0;
  let featureSerializationCalls = 0;
  globalThis.__GBDRAW_TEST_HOOKS__ = {
    onStructuralMetric: (metric) => metrics.push(metric)
  };
  try {
    const features = [feature];
    globalThis.structuredClone = (value, ...args) => {
      if (value === features) featureCloneCalls += 1;
      return previousStructuredClone(value, ...args);
    };
    JSON.stringify = (value, ...args) => {
      if (value === features) featureSerializationCalls += 1;
      return previousStringify(value, ...args);
    };
    const projection = createEditorSvgProjection({ features });
    projection.project(buildSvg({ includeSuppressedLegend: false }).svg);
    const borrowed = metrics.find((metric) => metric.name === 'editorProjectionBorrowedFeatureOwnerCount');
    const accessed = metrics.find((metric) => metric.name === 'editorProjectionFeatureBindingAccessCount');
    assert.equal(borrowed.featureOwner, features);
    assert.equal(accessed.featureOwner, features);
    assert.equal(accessed.value, 1);
    assert.equal(featureCloneCalls, 0);
    assert.equal(featureSerializationCalls, 0);
  } finally {
    globalThis.structuredClone = previousStructuredClone;
    JSON.stringify = previousStringify;
    if (previousHooks === undefined) delete globalThis.__GBDRAW_TEST_HOOKS__;
    else globalThis.__GBDRAW_TEST_HOOKS__ = previousHooks;
  }
});

test('fresh projection applies ordered product rules while direct Feature overrides win', () => {
  const fixture = buildSvg({ includeSuppressedLegend: false });
  const secondFeatureElement = appendFeature(fixture.svg, 'rendered-b');
  const secondFeature = {
    recordKey: 'record-b',
    biologicalFeatureId: 'biological-b',
    svg_id: 'rendered-b',
    type: 'CDS',
    product: 'shared product',
    qualifiers: { protein_id: ['WP_B.1'] }
  };
  const firstFeature = {
    ...feature,
    product: 'shared product',
    qualifiers: { protein_id: ['WP_A.1'] }
  };
  const projection = createEditorSvgProjection({
    features: [firstFeature, secondFeature],
    featureVisibilityOverrides: { 'rendered-b': 'on' },
    featureVisibilityManualRules: [{
      recordId: '*',
      featureType: 'CDS',
      qualifier: 'product',
      value: '^shared product$',
      action: 'off'
    }]
  });
  projection.project(fixture.svg);
  assert.equal(fixture.featureElement.getAttribute('display'), 'none');
  assert.equal(secondFeatureElement.getAttribute('display'), null);
});

test('catalog-backed admission fails when a required rendered Feature binding is missing', () => {
  const fixture = buildSvg({ includeFeature: false, labels: 0, includeSuppressedLegend: false });
  const projection = createEditorSvgProjection({ features: [feature] });
  assert.throws(() => projection.project(fixture.svg, {
    item: {
      features: [{
        recordKey: 'record-a',
        biologicalFeatureId: 'biological-a',
        svgId: 'rendered-a'
      }]
    },
    requireFeatureBindings: true
  }), /missing a rendered feature binding/);
});

test('label replay retains stable targets across geometry changes and fails closed otherwise', () => {
  const projection = createEditorSvgProjection({
    labelTextFeatureOverrides: { 'rendered-a': 'retained edit' },
    labelTextFeatureOverrideSources: { 'rendered-a': 'source label' },
    labelVisibilityOverrides: { 'rendered-a': 'off' }
  });

  const changedGeometry = buildSvg({ includeSuppressedLegend: false });
  changedGeometry.labelElements[0].setAttribute('transform', 'translate(999,777)');
  const retained = projection.project(changedGeometry.svg);
  assert.equal(changedGeometry.labelElements[0].textContent, 'retained edit');
  assert.equal(changedGeometry.labelElements[0].getAttribute('display'), 'none');
  assert.deepEqual(retained.unresolvedLabelFeatureIds, []);

  const removedTarget = buildSvg({ labels: 0, includeSuppressedLegend: false });
  const removed = projection.project(removedTarget.svg);
  assert.deepEqual(removed.unresolvedLabelFeatureIds, ['rendered-a']);

  const ambiguousTarget = buildSvg({ labels: 2, includeSuppressedLegend: false });
  const ambiguous = projection.project(ambiguousTarget.svg);
  assert.deepEqual(ambiguous.unresolvedLabelFeatureIds, ['rendered-a']);
  assert.deepEqual(
    ambiguousTarget.labelElements.map((label) => label.textContent),
    ['source label', 'source label']
  );
});

test('History replay can reset removed label overrides before applying current state', () => {
  const fixture = buildSvg({ includeSuppressedLegend: false });
  fixture.labelElements[0].textContent = 'stale edit';
  fixture.labelElements[0].setAttribute('display', 'none');
  fixture.labelElements[0].setAttribute('data-gbdraw-label-visibility-preview', 'off');
  fixture.labelElements[0].setAttribute('data-gbdraw-label-renderer-display', '');
  const projection = createEditorSvgProjection({});
  const result = projection.project(fixture.svg, { resetLabelState: true });
  assert.equal(result.changed, true);
  assert.equal(fixture.labelElements[0].textContent, 'source label');
  assert.equal(fixture.labelElements[0].getAttribute('display'), null);
});

test('History replay resets removed Feature visibility without changing fresh renderer visibility', () => {
  const fresh = buildSvg({ includeSuppressedLegend: false });
  fresh.featureElement.setAttribute('display', 'none');
  const projection = createEditorSvgProjection({ features: [feature] });
  assert.equal(projection.project(fresh.svg).changed, false);
  assert.equal(fresh.featureElement.getAttribute('display'), 'none');

  const history = buildSvg({ includeSuppressedLegend: false });
  history.featureElement.setAttribute('display', 'none');
  const result = projection.project(history.svg, { resetFeatureVisibility: true });
  assert.equal(result.changed, false);
  assert.equal(history.featureElement.getAttribute('display'), 'none');
});

test('renderer-hidden direct on and default round-trip while an absent Feature is a no-op', () => {
  const fixture = buildSvg({ includeSuppressedLegend: false });
  fixture.featureElement.setAttribute('display', 'none');
  const directOn = createEditorSvgProjection({
    features: [feature],
    featureVisibilityOverrides: { 'rendered-a': 'on' }
  });
  assert.equal(directOn.project(fixture.svg).changed, true);
  assert.equal(fixture.featureElement.getAttribute('display'), null);

  const directDefault = createEditorSvgProjection({
    features: [feature],
    featureVisibilityOverrides: { 'rendered-a': 'default' }
  });
  assert.equal(directDefault.project(fixture.svg, { resetFeatureVisibility: true }).changed, true);
  assert.equal(fixture.featureElement.getAttribute('display'), 'none');

  const absent = buildSvg({ includeFeature: false, labels: 0, includeSuppressedLegend: false });
  assert.equal(directDefault.project(absent.svg).changed, false);
});

test('History manager replays renderer visibility provenance through undo and redo', async () => {
  const fixture = buildSvg({ includeSuppressedLegend: false });
  fixture.featureElement.setAttribute('display', 'none');
  let mode = 'default';
  const apply = (nextMode) => createEditorSvgProjection({
    features: [feature],
    featureVisibilityOverrides: nextMode === 'default' ? {} : { 'rendered-a': nextMode }
  }).project(fixture.svg, { resetFeatureVisibility: true });
  const history = createHistoryManager({
    buildIntent: () => ({ mode }),
    applyIntent: (intent) => {
      mode = intent.mode;
      apply(mode);
    },
    buildCheckpoint: () => ({ mode }),
    applyCheckpoint: () => {}
  });
  await history.captureBaseline('renderer visibility baseline');
  await history.runUndoable('Force Feature visible', () => {
    mode = 'on';
    apply(mode);
  });
  assert.equal(fixture.featureElement.getAttribute('display'), null);
  await history.undo();
  assert.equal(fixture.featureElement.getAttribute('display'), 'none');
  await history.redo();
  assert.equal(fixture.featureElement.getAttribute('display'), null);
});

test('Feature visibility provenance survives History, Save/Load, fresh Generate, and reflow projection', () => {
  const metrics = [];
  const previousHooks = globalThis.__GBDRAW_TEST_HOOKS__;
  globalThis.__GBDRAW_TEST_HOOKS__ = {
    onStructuralMetric: (metric) => metrics.push(metric)
  };
  try {
    const manualOn = createEditorSvgProjection({
      features: [feature],
      featureVisibilityManualRules: [{
        recordId: '*',
        featureType: 'CDS',
        qualifier: 'product',
        value: '^Unrelated feature caption$',
        action: 'show'
      }]
    });
    const noEditorDecision = createEditorSvgProjection({ features: [feature] });

    const historyFixture = buildSvg({ includeSuppressedLegend: false });
    historyFixture.featureElement.setAttribute('display', 'none');
    assert.equal(manualOn.project(historyFixture.svg).changed, true);
    assert.equal(historyFixture.featureElement.getAttribute('display'), null);
    assert.equal(
      historyFixture.featureElement.getAttribute('data-gbdraw-feature-renderer-display'),
      'none'
    );
    assert.equal(
      historyFixture.featureElement.getAttribute('data-gbdraw-feature-visibility-preview'),
      'on'
    );

    const undo = noEditorDecision.project(historyFixture.svg, { resetFeatureVisibility: true });
    assert.equal(undo.changed, true);
    assert.equal(historyFixture.featureElement.getAttribute('display'), 'none');
    assert.equal(
      historyFixture.featureElement.getAttribute('data-gbdraw-feature-visibility-preview'),
      null
    );
    assert.equal(manualOn.project(historyFixture.svg).changed, true);
    assert.equal(historyFixture.featureElement.getAttribute('display'), null);

    for (const boundary of ['Save/Load', 'fresh Generate', 'reflow']) {
      const fresh = buildSvg({ includeSuppressedLegend: false });
      fresh.featureElement.setAttribute('display', 'none');
      assert.equal(manualOn.project(fresh.svg).changed, true, boundary);
      assert.equal(fresh.featureElement.getAttribute('display'), null, boundary);
      assert.equal(
        fresh.featureElement.getAttribute('data-gbdraw-feature-renderer-display'),
        'none',
        boundary
      );
    }

    const visibleWithManualOff = buildSvg({ includeSuppressedLegend: false });
    const manualOff = createEditorSvgProjection({
      features: [feature],
      featureVisibilityOverrides: { 'rendered-a': 'default' },
      featureVisibilityManualRules: [{
        recordId: '*',
        featureType: 'CDS',
        qualifier: 'product',
        value: '^Unrelated feature caption$',
        action: 'off'
      }]
    });
    assert.equal(manualOff.project(visibleWithManualOff.svg).changed, true);
    assert.equal(visibleWithManualOff.featureElement.getAttribute('display'), 'none');
    assert.equal(
      noEditorDecision.project(visibleWithManualOff.svg, { resetFeatureVisibility: true }).changed,
      true
    );
    assert.equal(visibleWithManualOff.featureElement.getAttribute('display'), null);

    assert.ok(
      metrics.filter(({ name }) => name === 'featureRendererBaselineCaptureCount').length >= 6
    );
  } finally {
    if (previousHooks === undefined) delete globalThis.__GBDRAW_TEST_HOOKS__;
    else globalThis.__GBDRAW_TEST_HOOKS__ = previousHooks;
  }
});

const appendPositionedLegendEntry = (parent, caption, x, y) => {
  const entry = parent.appendChild(new FakeSvgElement('g', { 'data-legend-key': caption }));
  entry.appendChild(new FakeSvgElement('path', {
    fill: '#111111',
    transform: `translate(${x},${y})`
  }));
  entry.appendChild(new FakeSvgElement('text', {
    transform: `translate(${x},${y})`
  }, caption));
  return entry;
};

const buildPositionedLegendSvg = ({ orientation, anchors }) => {
  const svg = new FakeSvgElement('svg');
  const legend = svg.appendChild(new FakeSvgElement('g', { id: 'legend' }));
  const featureLegend = legend.appendChild(new FakeSvgElement('g', {
    id: orientation === 'horizontal' ? 'feature_legend_h' : 'feature_legend_v'
  }));
  anchors.forEach(({ caption, x, y }) => appendPositionedLegendEntry(featureLegend, caption, x, y));
  return { featureLegend, svg };
};

const projectAddedLegendEntries = ({ orientation, anchors, captions }) => {
  const fixture = buildPositionedLegendSvg({ orientation, anchors });
  const originalCaptions = anchors.map((entry) => entry.caption);
  const projection = createEditorSvgProjection({
    legendEntries: captions.map((caption) => ({ caption, color: '#111111' })),
    originalLegendOrder: originalCaptions,
    addedLegendCaptions: captions.filter((caption) => !originalCaptions.includes(caption))
  });
  projection.project(fixture.svg);
  return Object.fromEntries(
    fixture.featureLegend.querySelectorAll('g[data-legend-key]').map((entry) => [
      entry.getAttribute('data-legend-key'),
      entry.querySelector('text').getAttribute('transform')
    ])
  );
};

test('Legend replay extrapolates renderer-owned horizontal and vertical slot topology', () => {
  assert.deepEqual(projectAddedLegendEntries({
    orientation: 'horizontal',
    anchors: [
      { caption: 'A', x: 20, y: 0 },
      { caption: 'B', x: 120, y: 0 }
    ],
    captions: ['A', 'B', 'C', 'D']
  }), {
    A: 'translate(20,0)',
    B: 'translate(120,0)',
    C: 'translate(220, 0)',
    D: 'translate(320, 0)'
  });
  assert.deepEqual(projectAddedLegendEntries({
    orientation: 'vertical',
    anchors: [
      { caption: 'A', x: 0, y: 10 },
      { caption: 'B', x: 0, y: 30 }
    ],
    captions: ['A', 'B', 'C']
  }), {
    A: 'translate(0,10)',
    B: 'translate(0,30)',
    C: 'translate(0, 50)'
  });
});

test('Legend replay continues wrapped renderer rows before starting the next row', () => {
  assert.deepEqual(projectAddedLegendEntries({
    orientation: 'horizontal',
    anchors: [
      { caption: 'A', x: 0, y: 0 },
      { caption: 'B', x: 100, y: 0 },
      { caption: 'C', x: 0, y: 20 }
    ],
    captions: ['A', 'B', 'C', 'D', 'E']
  }), {
    A: 'translate(0,0)',
    B: 'translate(100,0)',
    C: 'translate(0,20)',
    D: 'translate(100, 20)',
    E: 'translate(0, 40)'
  });
});

test('Legend replay reuses a removed renderer slot before extrapolating another slot', () => {
  const fixture = buildPositionedLegendSvg({
    orientation: 'horizontal',
    anchors: [
      { caption: 'A', x: 20, y: 0 },
      { caption: 'B', x: 120, y: 0 }
    ]
  });
  const projection = createEditorSvgProjection({
    legendEntries: [
      { caption: 'A', color: '#111111' },
      { caption: 'C', color: '#111111' }
    ],
    deletedLegendEntries: [{ caption: 'B', originalCaption: 'B' }],
    originalLegendOrder: ['A', 'B'],
    addedLegendCaptions: ['C']
  });

  projection.project(fixture.svg);
  assert.deepEqual(
    fixture.featureLegend.querySelectorAll('g[data-legend-key]').map((entry) => [
      entry.getAttribute('data-legend-key'),
      entry.querySelector('text').getAttribute('transform')
    ]),
    [
      ['A', 'translate(20,0)'],
      ['C', 'translate(120, 0)']
    ]
  );
});

test('horizontal Legend topology is stable across undo/redo and fresh admission boundaries', () => {
  const metrics = [];
  const previousHooks = globalThis.__GBDRAW_TEST_HOOKS__;
  globalThis.__GBDRAW_TEST_HOOKS__ = {
    onStructuralMetric: (metric) => metrics.push(metric)
  };
  const fixture = buildPositionedLegendSvg({
    orientation: 'horizontal',
    anchors: [
      { caption: 'A', x: 20, y: 0 },
      { caption: 'B', x: 120, y: 0 }
    ]
  });
  const withAddedEntry = createEditorSvgProjection({
    legendEntries: [
      { caption: 'A', color: '#111111' },
      { caption: 'B', color: '#111111' },
      { caption: 'C', color: '#111111' }
    ],
    originalLegendOrder: ['A', 'B'],
    addedLegendCaptions: ['C']
  });
  const withoutAddedEntry = createEditorSvgProjection({
    legendEntries: [
      { caption: 'A', color: '#111111' },
      { caption: 'B', color: '#111111' }
    ],
    originalLegendOrder: ['A', 'B']
  });

  withAddedEntry.project(fixture.svg);
  assert.equal(
    fixture.featureLegend.querySelector('g[data-legend-key="C"]')
      .querySelector('text').getAttribute('transform'),
    'translate(220, 0)'
  );
  withoutAddedEntry.project(fixture.svg);
  assert.deepEqual(
    fixture.featureLegend.querySelectorAll('g[data-legend-key]')
      .map((entry) => entry.getAttribute('data-legend-key')),
    ['A', 'B']
  );
  withAddedEntry.project(fixture.svg);
  assert.equal(
    fixture.featureLegend.querySelector('g[data-legend-key="C"]')
      .querySelector('text').getAttribute('transform'),
    'translate(220, 0)'
  );

  for (const boundary of ['Save/Load', 'fresh Generate', 'reflow']) {
    assert.equal(projectAddedLegendEntries({
      orientation: 'horizontal',
      anchors: [
        { caption: 'A', x: 20, y: 0 },
        { caption: 'B', x: 120, y: 0 }
      ],
      captions: ['A', 'B', 'C']
    }).C, 'translate(220, 0)', boundary);
  }
  assert.ok(metrics.some(({ name }) => name === 'legendTopologyResolutionCount'));
  if (previousHooks === undefined) delete globalThis.__GBDRAW_TEST_HOOKS__;
  else globalThis.__GBDRAW_TEST_HOOKS__ = previousHooks;
});

test('History manager restores and reapplies a horizontal renderer Legend slot', async () => {
  const fixture = buildPositionedLegendSvg({
    orientation: 'horizontal',
    anchors: [
      { caption: 'A', x: 20, y: 0 },
      { caption: 'B', x: 120, y: 0 }
    ]
  });
  let captions = ['A', 'B'];
  const apply = () => createEditorSvgProjection({
    legendEntries: captions.map((caption) => ({ caption, color: '#111111' })),
    originalLegendOrder: ['A', 'B'],
    addedLegendCaptions: captions.filter((caption) => !['A', 'B'].includes(caption))
  }).project(fixture.svg);
  const history = createHistoryManager({
    buildIntent: () => ({ captions: [...captions] }),
    applyIntent: (intent) => {
      captions = [...intent.captions];
      apply();
    },
    buildCheckpoint: () => ({ captions: [...captions] }),
    applyCheckpoint: () => {}
  });
  await history.captureBaseline('horizontal Legend baseline');
  await history.runUndoable('Add horizontal Legend entry', () => {
    captions.push('C');
    apply();
  });
  assert.equal(
    fixture.featureLegend.querySelector('g[data-legend-key="C"]')
      .querySelector('text').getAttribute('transform'),
    'translate(220, 0)'
  );
  await history.undo();
  assert.deepEqual(
    fixture.featureLegend.querySelectorAll('g[data-legend-key]')
      .map((entry) => entry.getAttribute('data-legend-key')),
    ['A', 'B']
  );
  await history.redo();
  assert.equal(
    fixture.featureLegend.querySelector('g[data-legend-key="C"]')
      .querySelector('text').getAttribute('transform'),
    'translate(220, 0)'
  );
});

test('invalid persisted paint fails before an invalid value is projected', () => {
  const fixture = buildSvg({ includeSuppressedLegend: false });
  const projection = createEditorSvgProjection({
    features: [feature],
    featureColorOverrides: { [stableKey]: { color: 'url(javascript:bad)' } }
  });
  assert.throws(() => projection.project(fixture.svg), /Invalid feature fill override/);
  assert.equal(fixture.featureElement.getAttribute('fill'), '#111111');
});

test('projection rolls back an early DOM mutation when a later target fails', () => {
  const fixture = buildSvg({ includeSuppressedLegend: false });
  const secondElement = appendFeature(fixture.svg, 'rendered-b');
  const secondFeature = {
    id: 'feature-b',
    recordKey: 'record-b',
    biologicalFeatureId: 'biological-b',
    svg_id: 'rendered-b',
    type: 'CDS',
    qualifiers: {}
  };
  const originalSetAttribute = secondElement.setAttribute.bind(secondElement);
  secondElement.setAttribute = (name, value) => {
    if (name === 'fill' && value === '#fedcba') throw new Error('late projection target failed');
    originalSetAttribute(name, value);
  };
  const projection = createEditorSvgProjection({
    features: [feature, secondFeature],
    featureColorOverrides: {
      [stableKey]: { color: '#abcdef' },
      ['record-b\u0000biological-b']: { color: '#fedcba' }
    }
  });

  assert.throws(() => projection.project(fixture.svg), /late projection target failed/);
  assert.equal(fixture.featureElement.getAttribute('fill'), '#111111');
  assert.equal(secondElement.getAttribute('fill'), '#111111');
});
