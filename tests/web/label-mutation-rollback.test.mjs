import assert from 'node:assert/strict';
import test from 'node:test';

import { createFeatureLabelActions } from '../../gbdraw/web/js/app/feature-editor/label-actions.js';
import { createPreviewRuntime } from '../../gbdraw/web/js/app/preview-runtime.js';
import {
  FakeSvgElement,
  appendEditableLabel,
  serializeFakeSvg
} from './helpers/editor-svg-fixture.mjs';

const ref = (value) => ({ value });

test('a failed label Result commit restores DOM, canonical overrides, dirty state, and indexes', () => {
  globalThis.CSS ??= { escape: (value) => String(value) };
  const svg = new FakeSvgElement('svg');
  const label = appendEditableLabel(svg, 'feature-a', 'renderer label');
  label.setAttribute('data-label-key', 'label-1');
  const querySelector = svg.querySelector.bind(svg);
  svg.querySelector = (selector) => (
    selector === 'text[data-label-key="label-1"]' ? label : querySelector(selector)
  );
  const initialResult = { name: 'label.svg', content: serializeFakeSvg(svg) };
  const labelTextFeatureOverrides = {};
  const labelTextFeatureOverrideSources = {};
  const state = {
    autoLabelReflowEnabled: ref(false),
    clickedFeature: ref(null),
    editableLabels: ref([]),
    extractedFeatures: ref([]),
    filterMode: ref('None'),
    form: { labels_mode: 'out', show_labels_linear: 'all' },
    globalLabelModeDialog: {},
    labelOverrideBuildWarning: ref(''),
    labelOverrideContextKey: ref(''),
    labelReflowForceRequestReason: ref(''),
    labelReflowForceRequestSeq: ref(0),
    labelReflowLastError: ref(null),
    labelReflowRequestReason: ref(''),
    labelReflowRequestSeq: ref(0),
    labelTextBulkOverrides: {},
    labelTextFeatureOverrideSources,
    labelTextFeatureOverrides,
    labelTextScopeDialog: {
      featureId: 'feature-a',
      labelKey: 'label-1',
      matchingCount: 1,
      newText: 'edited label',
      show: true,
      sourceText: 'renderer label'
    },
    labelVisibilityOverrides: {},
    manualWhitelist: [],
    mode: ref('circular'),
    results: ref([initialResult]),
    selectedResultIndex: ref(0),
    skipCaptureBaseConfig: ref(false),
    svgContainer: ref({ querySelector: (selector) => selector === 'svg' ? svg : null })
  };
  const runtime = createPreviewRuntime({
    state,
    serializeSvg: () => {
      throw new Error('label serialization failed');
    }
  });
  runtime.mountResultSvg(0, svg);
  const retainedFeatureIndex = { retained: true };
  runtime.getActiveRuntime().indexes.features = retainedFeatureIndex;
  const actions = createFeatureLabelActions({ state, previewRuntime: runtime });

  assert.throws(
    () => actions.handleLabelTextScopeChoice('single'),
    /label serialization failed/
  );

  assert.equal(label.textContent, 'renderer label');
  assert.deepEqual(labelTextFeatureOverrides, {});
  assert.deepEqual(labelTextFeatureOverrideSources, {});
  assert.equal(state.results.value[0], initialResult);
  assert.equal(state.skipCaptureBaseConfig.value, false);
  assert.equal(runtime.getActiveRuntime().dirty, false);
  assert.deepEqual([...runtime.getActiveRuntime().dirtyReasons], []);
  assert.equal(runtime.getActiveRuntime().indexes.features, retainedFeatureIndex);
});
