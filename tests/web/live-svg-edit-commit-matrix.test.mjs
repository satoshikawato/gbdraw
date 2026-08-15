import assert from 'node:assert/strict';
import test from 'node:test';

import { createPreviewRuntime } from '../../gbdraw/web/js/app/preview-runtime.js';
import {
  applyLabelVisibility,
  applyLegendColorOverridesToSvg,
  applyStrokeOverridesToSvg,
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
  const results = new Proxy([{ name: 'diagram.svg', content: initialContent }], {
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
    label,
    runtime,
    secondFeature,
    state,
    swatch,
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
        mutate: () => setLabelText(fixture.label, 'edited label')
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
        mutate: () => applyLabelVisibility(fixture.label, 'off', { markPreview: true })
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
        mutate: () => applyLegendColorOverridesToSvg({
          svg: fixture.state.svgContainer.value.querySelector('svg'),
          legendColorOverrides: canonical.legendFills
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
        mutate: () => applyStrokeOverridesToSvg({
          svg: fixture.state.svgContainer.value.querySelector('svg'),
          legendStrokeOverrides: canonical.legendStrokes
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
    const owners = {
      featureCatalog: {},
      extractedFeatures: [],
      biologicalFeatures: [],
      orthogroups: []
    };
    const ownerReferences = { ...owners };
    const execution = {
      workerConstruction: 0,
      workerInitialization: 0,
      workerRuns: 0,
      pythonCalls: 0,
      historyEntries: 0
    };
    const canonicalBefore = JSON.stringify(canonical);

    assert.equal(directEditCase.apply(fixture, canonical), true);
    assert.notEqual(JSON.stringify(canonical), canonicalBefore);
    directEditCase.assertDom(fixture);
    assert.notEqual(fixture.state.results.value[0].content, fixture.initialContent);
    assert.equal(fixture.counts.serializations(), 1);
    assert.equal(fixture.counts.replacements(), 1);
    assert.deepEqual(execution, {
      workerConstruction: 0,
      workerInitialization: 0,
      workerRuns: 0,
      pythonCalls: 0,
      historyEntries: 0
    });
    Object.entries(ownerReferences).forEach(([key, owner]) => assert.equal(owners[key], owner));

    const resultAfterFirstEdit = fixture.state.results.value[0];
    assert.equal(directEditCase.apply(fixture, canonical), false);
    assert.equal(fixture.state.results.value[0], resultAfterFirstEdit);
    assert.equal(fixture.counts.serializations(), 1);
    assert.equal(fixture.counts.replacements(), 1);
    assert.equal(execution.historyEntries, 0);
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
    mutate: () => setLabelText(fixture.label, canonical.labelText['feature-a'])
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
