import assert from 'node:assert/strict';
import test from 'node:test';

import { createFeatureLabelActions } from '../../gbdraw/web/js/app/feature-editor/label-actions.js';
import { installFakeSvgDom } from './fake-svg-dom.mjs';

installFakeSvgDom();
globalThis.CSS ||= { escape: (value) => String(value) };

const ref = (value) => ({ value });

const labelSvg = ({
  featureId = 'feature:one/[a]',
  marker = '1',
  leaderCount = 1,
  targetKey = 'label-1'
} = {}) => {
  const markerAttribute = marker === null
    ? ''
    : ` data-gbdraw-label-binding-schema="${marker}"`;
  const leaders = Array.from({ length: leaderCount }, (_, index) => (
    `<line data-label-feature-id="${featureId}" data-part="leader-${index}" />`
  )).join('');
  return new DOMParser().parseFromString([
    '<svg>',
    `<text data-label-editable="true" data-label-key="${targetKey}" `,
    `data-label-feature-id="${featureId}"${markerAttribute}></text>`,
    leaders,
    `<line data-label-feature-id="${featureId}-similar" data-part="similar" />`,
    '<line data-label-feature-id="unrelated" data-part="unrelated" />',
    '</svg>'
  ].join(''), 'image/svg+xml').documentElement;
};

const buildHarness = ({
  svg = labelSvg(),
  featureId = 'feature:one/[a]',
  labelKey = 'label-1',
  visibilityOverrides = {}
} = {}) => {
  const mutations = { dirty: 0, flush: 0 };
  const state = {
    mode: ref('linear'),
    generatedMode: ref('linear'),
    form: { show_labels_linear: 'all' },
    filterMode: ref('None'),
    manualWhitelist: [],
    results: ref([{ content: '<svg></svg>' }]),
    selectedResultIndex: ref(0),
    svgContainer: ref({ querySelector: (selector) => (selector === 'svg' ? svg : null) }),
    skipCaptureBaseConfig: ref(false),
    editableLabels: ref([{
      key: labelKey,
      featureId,
      text: '',
      sourceText: ''
    }]),
    extractedFeatures: ref([]),
    clickedFeature: ref({
      svg_id: featureId,
      feat: { type: 'CDS' },
      labelText: '',
      labelSourceText: '',
      labelVisibility: 'default',
      hasEditableLabel: true
    }),
    labelTextScopeDialog: { show: false },
    labelTextFeatureOverrides: {},
    labelTextBulkOverrides: {},
    labelTextFeatureOverrideSources: {},
    labelVisibilityOverrides: { ...visibilityOverrides },
    labelOverrideContextKey: ref(''),
    labelOverrideBuildWarning: ref(''),
    globalLabelModeDialog: { show: false },
    autoLabelReflowEnabled: ref(false),
    labelReflowRequestSeq: ref(0),
    labelReflowRequestReason: ref(''),
    labelReflowForceRequestSeq: ref(0),
    labelReflowForceRequestReason: ref(''),
    labelReflowLastError: ref(null)
  };
  const actions = createFeatureLabelActions({
    state,
    previewRuntime: {
      markActiveResultDirty() {
        mutations.dirty += 1;
        return true;
      },
      flushActiveResult() {
        mutations.flush += 1;
      }
    },
    rulePreparation: {}
  });
  return { actions, mutations, state, svg };
};

const exactParts = (svg, featureId) => svg
  .querySelectorAll('[data-label-feature-id]')
  .filter((element) => element.getAttribute('data-label-feature-id') === featureId);

for (const leaderCount of [0, 1, 2]) {
  test(`visibility Off and On mutate one complete ${leaderCount}-leader visual unit`, async () => {
    const harness = buildHarness({ svg: labelSvg({ leaderCount }) });
    harness.state.clickedFeature.value.labelVisibility = 'off';
    await harness.actions.updateClickedFeatureLabelText();

    const parts = exactParts(harness.svg, harness.state.clickedFeature.value.svg_id);
    assert.equal(parts.length, leaderCount + 1);
    parts.forEach((part) => {
      assert.equal(part.getAttribute('display'), 'none');
      assert.equal(part.getAttribute('data-gbdraw-label-visibility-preview'), 'off');
    });
    assert.equal(harness.svg.querySelector('[data-part="similar"]').getAttribute('display'), null);
    assert.equal(harness.svg.querySelector('[data-part="unrelated"]').getAttribute('display'), null);
    assert.deepEqual(harness.mutations, { dirty: 1, flush: 1 });
    assert.equal(harness.state.labelReflowForceRequestSeq.value, 0);

    harness.state.clickedFeature.value.labelVisibility = 'on';
    await harness.actions.updateClickedFeatureLabelText();
    parts.forEach((part) => {
      assert.equal(part.getAttribute('display'), null);
      assert.equal(part.getAttribute('data-gbdraw-label-visibility-preview'), null);
    });
    assert.deepEqual(harness.mutations, { dirty: 2, flush: 2 });
    assert.equal(harness.state.labelReflowForceRequestSeq.value, 0);
  });
}

test('Default restores only display values owned by the visibility preview', async () => {
  const harness = buildHarness();
  harness.state.clickedFeature.value.labelVisibility = 'off';
  await harness.actions.updateClickedFeatureLabelText();
  harness.state.clickedFeature.value.labelVisibility = 'default';
  await harness.actions.updateClickedFeatureLabelText();
  exactParts(harness.svg, harness.state.clickedFeature.value.svg_id).forEach((part) => {
    assert.equal(part.getAttribute('display'), null);
    assert.equal(part.getAttribute('data-gbdraw-label-visibility-preview'), null);
  });

  const authoredHidden = buildHarness();
  const leader = authoredHidden.svg.querySelector('[data-part="leader-0"]');
  leader.setAttribute('display', 'none');
  authoredHidden.state.clickedFeature.value.labelVisibility = 'on';
  await authoredHidden.actions.updateClickedFeatureLabelText();
  assert.equal(leader.getAttribute('display'), 'none');
  assert.deepEqual(authoredHidden.mutations, { dirty: 0, flush: 0 });
});

test('stored visibility projection uses the same complete visual-unit mutation', () => {
  const featureId = 'feature:one/[a]';
  const harness = buildHarness({ visibilityOverrides: { [featureId]: 'off' } });
  harness.actions.reconcileLabelOverrides();
  exactParts(harness.svg, featureId).forEach((part) => {
    assert.equal(part.getAttribute('display'), 'none');
    assert.equal(part.getAttribute('data-gbdraw-label-visibility-preview'), 'off');
  });
  assert.deepEqual(harness.mutations, { dirty: 1, flush: 1 });
  assert.equal(harness.state.labelReflowForceRequestSeq.value, 0);
});

test('one action serializes combined text and visual-unit visibility changes once', async () => {
  const harness = buildHarness();
  const textElement = harness.svg.querySelector('text[data-label-key="label-1"]');
  textElement.textContent = 'Original label';
  harness.state.editableLabels.value[0].text = 'Original label';
  harness.state.editableLabels.value[0].sourceText = 'Original label';
  Object.assign(harness.state.clickedFeature.value, {
    labelText: 'Renamed label',
    labelSourceText: 'Original label',
    labelVisibility: 'off'
  });
  await harness.actions.updateClickedFeatureLabelText();
  assert.equal(textElement.textContent, 'Renamed label');
  exactParts(harness.svg, harness.state.clickedFeature.value.svg_id).forEach((part) => {
    assert.equal(part.getAttribute('display'), 'none');
  });
  assert.deepEqual(harness.mutations, { dirty: 1, flush: 1 });
});

for (const [name, options, labelKey] of [
  ['missing schema marker', { marker: null }, 'label-1'],
  ['unsupported schema marker', { marker: '2' }, 'label-1'],
  ['empty exact identity', { featureId: '' }, 'label-1'],
  ['missing target text', {}, 'label-missing']
]) {
  test(`${name} fails closed and queues exactly one forced regeneration`, async () => {
    const svg = labelSvg(options);
    const featureId = options.featureId === '' ? 'feature:one/[a]' : 'feature:one/[a]';
    const harness = buildHarness({ svg, featureId, labelKey });
    harness.state.clickedFeature.value.labelVisibility = 'off';
    await harness.actions.updateClickedFeatureLabelText();
    svg.querySelectorAll('[data-label-feature-id]').forEach((part) => {
      assert.equal(part.getAttribute('display'), null);
      assert.equal(part.getAttribute('data-gbdraw-label-visibility-preview'), null);
    });
    assert.deepEqual(harness.mutations, { dirty: 0, flush: 0 });
    assert.equal(harness.state.labelReflowForceRequestSeq.value, 1);
    assert.equal(harness.state.labelReflowForceRequestReason.value, 'label-visibility-apply');
  });
}

test('stored override on a metadata-free Result fails closed without partial mutation', () => {
  const featureId = 'feature:one/[a]';
  const harness = buildHarness({
    svg: labelSvg({ marker: null }),
    visibilityOverrides: { [featureId]: 'off' }
  });
  harness.actions.reconcileLabelOverrides();
  exactParts(harness.svg, featureId).forEach((part) => {
    assert.equal(part.getAttribute('display'), null);
    assert.equal(part.getAttribute('data-gbdraw-label-visibility-preview'), null);
  });
  assert.deepEqual(harness.mutations, { dirty: 0, flush: 0 });
  assert.equal(harness.state.labelReflowForceRequestSeq.value, 1);
  assert.equal(
    harness.state.labelReflowForceRequestReason.value,
    'label-visibility-binding-refresh'
  );
});
