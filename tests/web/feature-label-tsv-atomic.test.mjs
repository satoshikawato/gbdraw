import assert from 'node:assert/strict';
import test from 'node:test';

import { createFeatureLabelActions } from '../../gbdraw/web/js/app/feature-editor/label-actions.js';

const ref = (value) => ({ value });

class FakeElement {
  constructor(attributes = {}, text = '') {
    this.attributes = new Map(Object.entries(attributes));
    this.textContent = text;
    this.style = {};
    Object.defineProperty(this.style, 'cursor', {
      set: (value) => this.setAttribute('style', `cursor: ${value}`)
    });
  }

  get id() {
    return this.getAttribute('id') || '';
  }

  getAttribute(name) {
    return this.attributes.has(name) ? this.attributes.get(name) : null;
  }

  hasAttribute(name) {
    return this.attributes.has(name);
  }

  setAttribute(name, value) {
    this.attributes.set(name, String(value));
  }

  removeAttribute(name) {
    this.attributes.delete(name);
  }

  querySelector() {
    return null;
  }

  closest(selector) {
    return selector === 'g[id]' ? { id: 'features_1' } : null;
  }

  getBBox() {
    return { x: 10, y: 20, width: 10, height: 10 };
  }

  getCTM() {
    return null;
  }
}

class FakeSvg {
  constructor() {
    this.feature = new FakeElement({ id: 'f1', 'data-gbdraw-feature-id': 'f1' });
    this.label = new FakeElement({
      'dominant-baseline': 'central',
      'data-label-editable': 'true',
      'data-label-key': 'label-1',
      'data-label-source-text': 'Source label',
      'data-label-feature-id': 'f1',
      'data-gbdraw-label-visibility-preview': 'off',
      display: 'none',
      style: 'font-size: 10px'
    }, 'Before label');
  }

  querySelectorAll(selector) {
    if (selector === 'text[dominant-baseline="central"]') return [this.label];
    if (selector === 'text[data-label-editable="true"]') {
      return this.label.getAttribute('data-label-editable') === 'true' ? [this.label] : [];
    }
    if (selector.includes('path[data-gbdraw-feature-id]')) return [this.feature];
    return [];
  }

  createSVGPoint() {
    return { x: 0, y: 0, matrixTransform: () => ({ x: 0, y: 0 }) };
  }
}

const makeState = (svg = new FakeSvg()) => {
  const initialEditableLabels = [{
    key: 'label-1',
    text: 'Before label',
    sourceText: 'Source label',
    featureId: 'f1',
    kind: 'regular',
    draftText: 'Before label'
  }];
  const dialogResolver = () => {};
  return {
    state: {
      mode: ref('linear'),
      form: {},
      filterMode: ref('None'),
      manualWhitelist: [],
      results: ref([{ name: 'diagram', content: 'before-svg' }]),
      selectedResultIndex: ref(0),
      svgContainer: ref(svg ? { querySelector: (selector) => selector === 'svg' ? svg : null } : null),
      skipCaptureBaseConfig: ref(false),
      editableLabels: ref(initialEditableLabels),
      extractedFeatures: ref([{
        svg_id: 'f1',
        record_id: 'Rec1',
        type: 'CDS',
        qualifiers: {}
      }]),
      clickedFeature: ref(null),
      labelTextFeatureOverrides: { f1: 'Before label' },
      labelTextBulkOverrides: { 'Source label': 'Before all' },
      labelTextFeatureOverrideSources: { f1: 'Source label' },
      labelVisibilityOverrides: { f1: 'off' },
      labelOverrideContextKey: ref('linear:f1'),
      labelOverrideBuildWarning: ref('existing warning'),
      globalLabelModeDialog: {
        show: true,
        featureId: 'f1',
        featureType: 'CDS',
        resolve: dialogResolver
      },
      autoLabelReflowEnabled: ref(true),
      labelReflowRequestSeq: ref(0),
      labelReflowRequestReason: ref(''),
      labelReflowForceRequestSeq: ref(0),
      labelReflowForceRequestReason: ref(''),
      labelReflowLastError: ref(null)
    },
    svg,
    initialEditableLabels,
    dialogResolver
  };
};

const inputEvent = (text) => ({
  target: {
    files: [{ text: async () => text }],
    value: 'labels.tsv'
  }
});

const validTsv = 'Rec1\tCDS\thash\t^f1$\tImported label\n';

test('Label TSV preflight failure does not mutate the current label state', async () => {
  const { state, svg, initialEditableLabels, dialogResolver } = makeState();
  const originalResults = state.results.value;
  const originalResult = originalResults[0];
  const originalAttributes = Object.fromEntries(svg.label.attributes);
  const alerts = [];
  globalThis.window = { alert: (message) => alerts.push(message) };
  const originalConsoleError = console.error;
  console.error = () => {};
  try {
    const actions = createFeatureLabelActions({
      state,
      serializeSvg: () => { throw new Error('serialization failed'); }
    });
    const event = inputEvent(validTsv);
    await actions.loadLabelOverrideTable(event);

    assert.equal(event.target.value, '');
    assert.equal(svg.label.textContent, 'Before label');
    assert.deepEqual(Object.fromEntries(svg.label.attributes), originalAttributes);
    assert.deepEqual(state.labelTextFeatureOverrides, { f1: 'Before label' });
    assert.deepEqual(state.labelTextBulkOverrides, { 'Source label': 'Before all' });
    assert.deepEqual(state.labelTextFeatureOverrideSources, { f1: 'Source label' });
    assert.deepEqual(state.labelVisibilityOverrides, { f1: 'off' });
    assert.equal(state.labelOverrideContextKey.value, 'linear:f1');
    assert.equal(state.labelOverrideBuildWarning.value, 'existing warning');
    assert.equal(state.editableLabels.value, initialEditableLabels);
    assert.equal(state.results.value, originalResults);
    assert.equal(state.results.value[0], originalResult);
    assert.equal(state.skipCaptureBaseConfig.value, false);
    assert.equal(state.globalLabelModeDialog.show, true);
    assert.equal(state.globalLabelModeDialog.resolve, dialogResolver);
    assert.equal(state.labelReflowRequestSeq.value, 0);
    assert.match(alerts.at(-1), /^Failed to load label TSV\. serialization failed$/);
  } finally {
    console.error = originalConsoleError;
  }
});

test('Label TSV commit publishes the prepared label and override state together', async () => {
  const { state, svg } = makeState();
  const alerts = [];
  globalThis.window = { alert: (message) => alerts.push(message) };
  const actions = createFeatureLabelActions({
    state,
    serializeSvg: () => '<svg>imported</svg>'
  });
  const event = inputEvent(validTsv);
  await actions.loadLabelOverrideTable(event);

  assert.equal(svg.label.textContent, 'Imported label');
  assert.equal(svg.label.hasAttribute('data-gbdraw-label-visibility-preview'), false);
  assert.equal(svg.label.hasAttribute('display'), false);
  assert.deepEqual(state.labelTextFeatureOverrides, { f1: 'Imported label' });
  assert.deepEqual(state.labelTextBulkOverrides, {});
  assert.deepEqual(state.labelTextFeatureOverrideSources, { f1: 'Source label' });
  assert.deepEqual(state.labelVisibilityOverrides, {});
  assert.equal(state.results.value[0].content, '<svg>imported</svg>');
  assert.equal(state.skipCaptureBaseConfig.value, true);
  assert.equal(state.globalLabelModeDialog.show, false);
  assert.equal(state.labelReflowRequestSeq.value, 1);
  assert.match(alerts.at(-1), /^Loaded 1 row\(s\)\. Applied to 1 label\(s\)\.$/);
});

test('Label TSV parsing completes before any state is cleared', async () => {
  const { state, svg } = makeState();
  const originalAttributes = Object.fromEntries(svg.label.attributes);
  const alerts = [];
  let serializeCount = 0;
  globalThis.window = { alert: (message) => alerts.push(message) };
  const originalConsoleError = console.error;
  console.error = () => {};
  try {
    const actions = createFeatureLabelActions({
      state,
      serializeSvg: () => { serializeCount += 1; return 'unexpected'; }
    });
    await actions.loadLabelOverrideTable(inputEvent('Rec1\tCDS\thash\t^f1$\n'));

    assert.equal(serializeCount, 0);
    assert.equal(svg.label.textContent, 'Before label');
    assert.deepEqual(Object.fromEntries(svg.label.attributes), originalAttributes);
    assert.deepEqual(state.labelTextFeatureOverrides, { f1: 'Before label' });
    assert.deepEqual(state.labelVisibilityOverrides, { f1: 'off' });
    assert.match(alerts.at(-1), /^Failed to load label TSV\./);
  } finally {
    console.error = originalConsoleError;
  }
});

test('Label TSV import without a diagram no longer calls the retired scope dialog', async () => {
  const { state } = makeState(null);
  const alerts = [];
  globalThis.window = { alert: (message) => alerts.push(message) };
  const actions = createFeatureLabelActions({ state });
  await actions.loadLabelOverrideTable(inputEvent(validTsv));

  assert.deepEqual(state.labelTextFeatureOverrides, {});
  assert.deepEqual(state.labelTextBulkOverrides, {});
  assert.deepEqual(state.labelTextFeatureOverrideSources, {});
  assert.deepEqual(state.labelVisibilityOverrides, {});
  assert.deepEqual(state.editableLabels.value, []);
  assert.equal(state.labelOverrideContextKey.value, '');
  assert.match(alerts.at(-1), /No diagram is currently displayed\.$/);
});
