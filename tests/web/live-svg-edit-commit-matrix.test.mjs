import assert from 'node:assert/strict';
import test from 'node:test';

import { createPreviewRuntime } from '../../gbdraw/web/js/app/preview-runtime.js';
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
  appendFeatureLegend,
  serializeFakeSvg
} from './helpers/editor-svg-fixture.mjs';

globalThis.CSS = { escape: (value) => String(value) };

const ref = (value) => ({ value });

const createFixture = ({ includeSecondFeature = false } = {}) => {
  const svg = new FakeSvgElement('svg');
  const feature = appendFeature(svg, 'feature-a');
  const secondFeature = includeSecondFeature ? appendFeature(svg, 'feature-b') : null;
  const label = appendEditableLabel(svg, 'feature-a', 'source label');
  const { swatch } = appendFeatureLegend(svg, 'Legend A');
  let replacementCount = 0;
  const initialContent = serializeFakeSvg(svg);
  const initialResult = { name: 'diagram.svg', content: initialContent };
  const results = new Proxy([initialResult], {
    set(target, property, value) {
      if (/^\d+$/.test(String(property))) replacementCount += 1;
      target[property] = value;
      return true;
    }
  });
  const state = {
    results: ref(results),
    selectedResultIndex: ref(0),
    skipCaptureBaseConfig: ref(false),
    svgContainer: ref({ querySelector: (selector) => selector === 'svg' ? svg : null })
  };
  let serializationCount = 0;
  const structuralMetrics = [];
  globalThis.__GBDRAW_TEST_HOOKS__ = {
    onStructuralMetric(metric) {
      structuralMetrics.push(metric);
    }
  };
  const runtime = createPreviewRuntime({
    state,
    serializeSvg: (targetSvg) => {
      serializationCount += 1;
      return serializeFakeSvg(targetSvg);
    }
  });
  runtime.mountResultSvg(0, svg);
  return {
    feature,
    initialContent,
    initialResult,
    label,
    runtime,
    secondFeature,
    state,
    swatch,
    structuralMetrics,
    counts: {
      replacements: () => replacementCount,
      serializations: () => serializationCount
    }
  };
};

const directEditCases = [
  {
    name: 'Feature fill',
    apply: (fixture, canonical) => {
      canonical.featureFills['feature-a'] = '#abcdef';
      return fixture.runtime.applyFeatureFillChanges([
        { featureId: 'feature-a', color: '#abcdef' }
      ]);
    },
    assertDom: ({ feature }) => assert.equal(feature.getAttribute('fill'), '#abcdef')
  },
  {
    name: 'Feature stroke',
    apply: (fixture, canonical) => {
      canonical.featureStrokes['feature-a'] = { strokeColor: '#abcdef', strokeWidth: 3 };
      return fixture.runtime.applyFeatureStrokeChanges([
        { featureId: 'feature-a', strokeColor: '#abcdef', strokeWidth: 3 }
      ]);
    },
    assertDom: ({ feature }) => {
      assert.equal(feature.getAttribute('stroke'), '#abcdef');
      assert.equal(feature.getAttribute('stroke-width'), '3');
    }
  },
  {
    name: 'Feature visibility',
    apply: (fixture, canonical) => {
      canonical.featureVisibility['feature-a'] = 'off';
      return fixture.runtime.applyFeatureVisibilityChanges([
        { featureId: 'feature-a', mode: 'off' }
      ]);
    },
    assertDom: ({ feature }) => assert.equal(feature.getAttribute('display'), 'none')
  },
  {
    name: 'Label text',
    apply: (fixture, canonical) => {
      canonical.labelText['feature-a'] = 'edited label';
      return fixture.runtime.commitDomEdit({
        reason: 'label-text',
        mutate: ({ mutation }) => setLabelText(
          fixture.label,
          'edited label',
          { mutation }
        )
      }).changed;
    },
    assertDom: ({ label }) => assert.equal(label.textContent, 'edited label')
  },
  {
    name: 'Label visibility',
    apply: (fixture, canonical) => {
      canonical.labelVisibility['feature-a'] = 'off';
      return fixture.runtime.commitDomEdit({
        reason: 'label-visibility',
        mutate: ({ mutation }) => applyLabelVisibility(
          fixture.label,
          'off',
          { markPreview: true, mutation }
        )
      }).changed;
    },
    assertDom: ({ label }) => assert.equal(label.getAttribute('display'), 'none')
  },
  {
    name: 'Legend fill',
    apply: (fixture, canonical) => {
      canonical.legendFills['Legend A'] = '#abcdef';
      return fixture.runtime.commitDomEdit({
        reason: 'legend-fill',
        invalidateIndexes: ['legend'],
        mutate: ({ svg, mutation }) => applyLegendColorOverridesToSvg({
          svg,
          legendColorOverrides: canonical.legendFills,
          mutation
        })
      }).changed;
    },
    assertDom: ({ swatch }) => assert.equal(swatch.getAttribute('fill'), '#abcdef')
  },
  {
    name: 'Legend stroke',
    apply: (fixture, canonical) => {
      canonical.legendStrokes['Legend A'] = { strokeColor: '#abcdef', strokeWidth: 3 };
      return fixture.runtime.commitDomEdit({
        reason: 'legend-stroke',
        invalidateIndexes: ['legend'],
        mutate: ({ svg, mutation }) => applyStrokeOverridesToSvg({
          svg,
          legendStrokeOverrides: canonical.legendStrokes,
          mutation
        })
      }).changed;
    },
    assertDom: ({ swatch }) => {
      assert.equal(swatch.getAttribute('stroke'), '#abcdef');
      assert.equal(swatch.getAttribute('stroke-width'), '3');
    }
  }
];

for (const directEditCase of directEditCases) {
  test(`${directEditCase.name} commits once and repeats as a no-op`, () => {
    const fixture = createFixture();
    const canonical = {
      featureFills: {},
      featureStrokes: {},
      featureVisibility: {},
      labelText: {},
      labelVisibility: {},
      legendFills: {},
      legendStrokes: {}
    };
    const canonicalBefore = JSON.stringify(canonical);

    assert.equal(directEditCase.apply(fixture, canonical), true);
    assert.notEqual(JSON.stringify(canonical), canonicalBefore);
    directEditCase.assertDom(fixture);
    assert.notEqual(fixture.state.results.value[0].content, fixture.initialContent);
    assert.equal(fixture.counts.serializations(), 1);
    assert.equal(fixture.counts.replacements(), 1);
    const replacementMetric = fixture.structuralMetrics.find(
      ({ name }) => name === 'domEditResultReplacementCount'
    );
    assert.ok(replacementMetric, 'the assertion observes the production commit boundary');
    assert.equal(replacementMetric.resultOwnerBefore, fixture.initialResult);
    assert.equal(replacementMetric.resultOwnerAfter, fixture.state.results.value[0]);
    assert.equal(replacementMetric.resultsOwner, fixture.state.results.value);
    assert.equal(replacementMetric.svgOwner, fixture.state.svgContainer.value.querySelector('svg'));
    assert.deepEqual(
      fixture.structuralMetrics.filter(({ name }) => [
        'workerConstructionCount',
        'workerInitializationCount',
        'pythonHelperRequestCount'
      ].includes(name)),
      [],
      'the shared production hook observed no Worker or Python helper lifecycle event'
    );

    const resultAfterFirstEdit = fixture.state.results.value[0];
    assert.equal(directEditCase.apply(fixture, canonical), false);
    assert.equal(fixture.state.results.value[0], resultAfterFirstEdit);
    assert.equal(fixture.counts.serializations(), 1);
    assert.equal(fixture.counts.replacements(), 1);
  });
}

test('bulk Feature fill changes multiple targets with one commit', () => {
  const fixture = createFixture({ includeSecondFeature: true });
  const canonical = {
    'feature-a': '#abcdef',
    'feature-b': '#fedcba'
  };
  assert.equal(fixture.runtime.applyFeatureFillChanges(
    Object.entries(canonical).map(([featureId, color]) => ({ featureId, color })),
    { reason: 'bulk-feature-fill' }
  ), true);
  assert.equal(fixture.feature.getAttribute('fill'), '#abcdef');
  assert.equal(fixture.secondFeature.getAttribute('fill'), '#fedcba');
  assert.equal(fixture.counts.serializations(), 1);
  assert.equal(fixture.counts.replacements(), 1);
});

test('a failed optional reflow does not roll back an already committed direct edit', async () => {
  const fixture = createFixture();
  const canonical = { labelText: { 'feature-a': 'edited before reflow' } };
  const committed = fixture.runtime.commitDomEdit({
    reason: 'label-text-before-reflow',
    mutate: ({ mutation }) => setLabelText(
      fixture.label,
      canonical.labelText['feature-a'],
      { mutation }
    )
  });
  const committedResult = fixture.state.results.value[0];
  assert.deepEqual(committed, {
    changed: true,
    flushed: true,
    resultIndex: 0,
    reason: 'label-text-before-reflow'
  });

  await assert.rejects(
    Promise.reject(new Error('optional reflow failed')),
    /optional reflow failed/
  );
  assert.equal(fixture.label.textContent, 'edited before reflow');
  assert.equal(canonical.labelText['feature-a'], 'edited before reflow');
  assert.equal(fixture.state.results.value[0], committedResult);
  assert.equal(fixture.counts.serializations(), 1);
  assert.equal(fixture.counts.replacements(), 1);
});

test('Linear direct edit and fresh-render replay converge through the shared owner', () => {
  const mounted = createFixture();
  mounted.state.svgContainer.value.querySelector('svg').setAttribute('data-diagram-mode', 'linear');
  assert.equal(mounted.runtime.applyFeatureFillChanges([
    { featureId: 'feature-a', color: '#0f766e' }
  ]), true);

  const freshSvg = new FakeSvgElement('svg', { 'data-diagram-mode': 'linear' });
  const freshFeature = appendFeature(freshSvg, 'feature-a');
  const features = [{
    svg_id: 'feature-a',
    recordKey: 'linear-record',
    biologicalFeatureId: 'linear-feature',
    type: 'CDS'
  }];
  createEditorSvgProjection({
    features,
    featureColorOverrides: {
      ['linear-record\u0000linear-feature']: { color: '#0f766e' }
    }
  }).project(freshSvg);

  assert.equal(mounted.feature.getAttribute('fill'), '#0f766e');
  assert.equal(freshFeature.getAttribute('fill'), '#0f766e');
});

test('multi-Result edit isolates Result 2 and its saved selection round-trips all artifacts', () => {
  const firstSvg = new FakeSvgElement('svg');
  const secondSvg = new FakeSvgElement('svg');
  appendFeature(firstSvg, 'feature-a');
  const secondFeature = appendFeature(secondSvg, 'feature-a');
  const firstResult = { name: 'record-1.svg', content: serializeFakeSvg(firstSvg) };
  const secondResult = { name: 'record-2.svg', content: serializeFakeSvg(secondSvg) };
  const results = [firstResult, secondResult];
  let mountedSvg = firstSvg;
  const state = {
    results: ref(results),
    selectedResultIndex: ref(0),
    skipCaptureBaseConfig: ref(false),
    svgContainer: ref({ querySelector: (selector) => selector === 'svg' ? mountedSvg : null })
  };
  const runtime = createPreviewRuntime({ state, serializeSvg: serializeFakeSvg });
  runtime.mountResultSvg(0, firstSvg);
  runtime.selectResult(1);
  mountedSvg = secondSvg;
  runtime.mountResultSvg(1, secondSvg);

  assert.equal(runtime.applyFeatureFillChanges([
    { featureId: 'feature-a', color: '#7c3aed' }
  ]), true);
  assert.equal(state.results.value[0], firstResult);
  assert.equal(state.results.value[0].content, firstResult.content);
  assert.notEqual(state.results.value[1], secondResult);
  assert.equal(secondFeature.getAttribute('fill'), '#7c3aed');

  const saved = JSON.stringify({
    results: state.results.value,
    ui: { selectedResultIndex: state.selectedResultIndex.value }
  });
  const loaded = JSON.parse(saved);
  assert.deepEqual(loaded.results.map(({ name }) => name), ['record-1.svg', 'record-2.svg']);
  assert.equal(loaded.results[0].content, firstResult.content);
  assert.equal(loaded.results[1].content, state.results.value[1].content);
  assert.equal(loaded.ui.selectedResultIndex, 1);
});
