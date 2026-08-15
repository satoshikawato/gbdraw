import assert from 'node:assert/strict';
import test from 'node:test';

import {
  applyLabelVisibility,
  applyLegendColorOverridesToSvg,
  applyStrokeOverridesToSvg,
  createEditorSvgProjection,
  setLabelText
} from '../../gbdraw/web/js/app/editor-svg-projection.js';
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
  assert.equal(result.changed, true);
  assert.equal(history.featureElement.getAttribute('display'), null);
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
