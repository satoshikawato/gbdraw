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
  visibilityOverrides = {},
  rulePreparation = {}
} = {}) => {
  const mutations = { commit: 0 };
  const state = {
    mode: ref('linear'),
    generatedMode: ref('linear'),
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
    hiddenLabelTextDialog: { show: false, featureId: '' },
    labelTextFeatureOverrides: {},
    labelTextBulkOverrides: {},
    labelTextFeatureOverrideSources: {},
    labelVisibilityOverrides: { ...visibilityOverrides },
    labelOverrideBuildWarning: ref(''),
    autoLabelReflowEnabled: ref(false),
    labelReflowRequestSeq: ref(0),
    labelReflowRequestReason: ref(''),
    labelReflowForceRequestSeq: ref(0),
    labelReflowForceRequestReason: ref(''),
    labelReflowLastError: ref(null)
  };
  const actions = createFeatureLabelActions({
    ref, computed: get => ({ get value() { return get(); } }),
    state,
    previewRuntime: {
      commitActiveResultEdit(reason) {
        assert.equal(reason, 'feature-label');
        mutations.commit += 1;
        return true;
      }
    },
    rulePreparation
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
    assert.deepEqual(harness.mutations, { commit: 1 });
    assert.equal(harness.state.labelReflowForceRequestSeq.value, 0);

    harness.state.clickedFeature.value.labelVisibility = 'on';
    await harness.actions.updateClickedFeatureLabelText();
    parts.forEach((part) => {
      assert.equal(part.getAttribute('display'), null);
      assert.equal(part.getAttribute('data-gbdraw-label-visibility-preview'), null);
    });
    assert.deepEqual(harness.mutations, { commit: 2 });
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
  assert.deepEqual(authoredHidden.mutations, { commit: 0 });
});

test('stored visibility projection uses the same complete visual-unit mutation', () => {
  const featureId = 'feature:one/[a]';
  const harness = buildHarness({ visibilityOverrides: { [featureId]: 'off' } });
  harness.actions.reconcileLabelOverrides();
  exactParts(harness.svg, featureId).forEach((part) => {
    assert.equal(part.getAttribute('display'), 'none');
    assert.equal(part.getAttribute('data-gbdraw-label-visibility-preview'), 'off');
  });
  assert.deepEqual(harness.mutations, { commit: 1 });
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
  assert.deepEqual(harness.mutations, { commit: 1 });
});

test('a text edit asks whether to show a label only when the feature has none', async () => {
  const featureId = 'feature:one/[a]';
  const labeled = buildHarness();
  Object.assign(labeled.state.clickedFeature.value, { labelText: 'Renamed', labelSourceText: 'Original' });
  await labeled.actions.updateClickedFeatureLabelText();
  assert.equal(labeled.state.hiddenLabelTextDialog.show, false);

  const harness = buildHarness();
  harness.state.editableLabels.value = [];
  Object.assign(harness.state.clickedFeature.value, { labelText: 'Renamed', labelSourceText: 'Original' });
  await harness.actions.updateClickedFeatureLabelText();
  assert.deepEqual({ ...harness.state.hiddenLabelTextDialog }, { show: true, featureId });
  assert.equal(harness.state.labelTextFeatureOverrides[featureId], 'Renamed');
  harness.actions.handleHiddenLabelTextChoice('text_only');
  assert.deepEqual({ ...harness.state.hiddenLabelTextDialog }, { show: false, featureId: '' });
  assert.deepEqual(harness.state.labelVisibilityOverrides, {});
  assert.equal(harness.state.labelReflowForceRequestSeq.value, 0);

  harness.state.clickedFeature.value.labelText = 'Renamed again';
  await harness.actions.updateClickedFeatureLabelText();
  assert.equal(harness.state.hiddenLabelTextDialog.show, true);
  harness.actions.handleHiddenLabelTextChoice('show');
  assert.equal(harness.state.hiddenLabelTextDialog.show, false);
  assert.deepEqual(harness.state.labelVisibilityOverrides, { [featureId]: 'on' });
  assert.equal(harness.state.clickedFeature.value.labelVisibility, 'on');
  assert.equal(harness.state.labelTextFeatureOverrides[featureId], 'Renamed again');
  assert.equal(harness.state.labelReflowForceRequestSeq.value, 1);
  assert.equal(harness.state.labelReflowForceRequestReason.value, 'label-visibility-apply');
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
    assert.deepEqual(harness.mutations, { commit: 0 });
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
  assert.deepEqual(harness.mutations, { commit: 0 });
  assert.equal(harness.state.labelReflowForceRequestSeq.value, 1);
  assert.equal(
    harness.state.labelReflowForceRequestReason.value,
    'label-visibility-binding-refresh'
  );
});

test('mounting Results with disjoint features keeps every label override (FE-01)', () => {
  const resultSvg = (featureId) => new DOMParser().parseFromString([
    '<svg>',
    `<path data-gbdraw-feature-id="${featureId}" d="M0 0" />`,
    `<text data-label-editable="true" data-label-key="label-1" data-label-feature-id="${featureId}" `,
    'data-gbdraw-label-binding-schema="1">source</text>',
    '</svg>'
  ].join(''), 'image/svg+xml').documentElement;
  const resultA = resultSvg('fa');
  const resultB = resultSvg('fb');
  const harness = buildHarness({ svg: resultA, featureId: 'fa' });
  const { state, actions } = harness;
  let mounted = resultA;
  state.svgContainer.value = { querySelector: (selector) => (selector === 'svg' ? mounted : null) };
  Object.assign(state.labelTextFeatureOverrides, { fb: 'EDITED_B' });
  Object.assign(state.labelTextFeatureOverrideSources, { fb: 'source' });
  Object.assign(state.labelTextBulkOverrides, { source: 'BULK' });
  Object.assign(state.labelVisibilityOverrides, { fb: 'off' });
  const before = JSON.stringify([
    state.labelTextFeatureOverrides, state.labelTextFeatureOverrideSources,
    state.labelTextBulkOverrides, state.labelVisibilityOverrides
  ]);
  for (const next of [resultB, resultA, resultB]) {
    mounted = next;
    actions.syncLabelEditor({ queueIncompleteVisibility: false });
    assert.equal(JSON.stringify([
      state.labelTextFeatureOverrides, state.labelTextFeatureOverrideSources,
      state.labelTextBulkOverrides, state.labelVisibilityOverrides
    ]), before);
  }
});

// F-3: Generate draws no label for a hidden feature. The label projection that
// live edits, History, and Result display share hides it with the feature.
test('a hidden feature hides its label through the stored visibility projection', () => {
  const featureId = 'feature:one/[a]';
  const harness = buildHarness();
  harness.state.featureVisibilityOverrides = { [featureId]: 'off' };
  harness.state.featureVisibilityManualRules = [];
  harness.actions.reconcileLabelOverrides();
  exactParts(harness.svg, featureId).forEach((part) => {
    assert.equal(part.getAttribute('display'), 'none');
  });
  assert.equal(harness.svg.querySelector('[data-part="unrelated"]').getAttribute('display'), null);

  delete harness.state.featureVisibilityOverrides[featureId];
  harness.actions.reconcileLabelOverrides();
  exactParts(harness.svg, featureId).forEach((part) => {
    assert.equal(part.getAttribute('display'), null);
  });
});

test('a feature visibility edit commits its label once and queues the label reflow', () => {
  const featureId = 'feature:one/[a]';
  const harness = buildHarness();
  harness.state.autoLabelReflowEnabled.value = true;
  harness.state.featureVisibilityOverrides = { [featureId]: 'off' };
  assert.equal(harness.actions.applyFeatureVisibilityToLabels('feature-visibility'), true);
  exactParts(harness.svg, featureId).forEach((part) => {
    assert.equal(part.getAttribute('display'), 'none');
  });
  assert.deepEqual(harness.mutations, { commit: 1 });
  assert.equal(harness.state.labelReflowRequestSeq.value, 1);
  assert.equal(harness.state.labelReflowRequestReason.value, 'feature-visibility');
  assert.equal(harness.state.labelReflowForceRequestSeq.value, 0);

  harness.state.featureVisibilityOverrides = {};
  assert.equal(harness.actions.applyFeatureVisibilityToLabels('feature-visibility', { reflow: false }), true);
  exactParts(harness.svg, featureId).forEach((part) => {
    assert.equal(part.getAttribute('display'), null);
  });
  assert.equal(harness.state.labelReflowRequestSeq.value, 1, 'a declined reflow is not queued');
});

// B6 (R2, R3): a Label TSV import writes the label intent of every batch Result
// once. The displayed Result binds its labels as the editor does; another
// Result's committed content carries the renderer's label bindings.
const batchLabel = (featureId, source, extra = '') => [
  `<path data-gbdraw-feature-id="${featureId}" d="M0 0" />`,
  `<text ${extra}data-label-feature-id="${featureId}" data-gbdraw-label-binding-schema="1" `,
  `data-label-source-text="${source}"></text>`
].join('');

test('a Label TSV import writes the label intent of every batch Result once (B6)', async () => {
  const displayed = new DOMParser().parseFromString(
    `<svg>${batchLabel('fa', 'alpha', 'dominant-baseline="central" ')}</svg>`, 'image/svg+xml'
  ).documentElement;
  displayed.querySelector('text').textContent = 'alpha';
  const evaluated = [];
  const harness = buildHarness({
    svg: displayed,
    featureId: 'fa',
    rulePreparation: {
      snapshot: () => ({}),
      isCurrent: () => true,
      async evaluate({ features, rules }) {
        evaluated.push(features.map((feature) => feature.label));
        return { winners: features.map((feature) => rules.findIndex((rule) => (
          new RegExp(rule.valueRegex).test(feature.label)
        ))) };
      }
    }
  });
  const { state, actions } = harness;
  state.errorLog = ref(null);
  state.results.value = [{ content: '<svg></svg>' }, { content: `<svg>${batchLabel('fb', 'beta')}</svg>` }];
  const messages = [];
  globalThis.window = { alert: (message) => messages.push(message) };
  try {
    const tsv = '*\tCDS\tlabel\t^beta$\tIMPORTED_B\n*\t*\tlabel\t^alpha$\tBULK_A\n';
    await actions.loadLabelOverrideTable({ target: { files: [{ text: async () => tsv }], value: 'labels.tsv' } });
  } finally {
    delete globalThis.window;
  }
  assert.deepEqual(evaluated, [['alpha', 'beta']], 'one evaluation covers the labels of both Results');
  assert.deepEqual(state.labelTextFeatureOverrides, { fb: 'IMPORTED_B' });
  assert.deepEqual(state.labelTextBulkOverrides, { alpha: 'BULK_A' });
  assert.deepEqual(state.labelTextFeatureOverrideSources, { fa: 'alpha', fb: 'beta' });
  assert.equal(displayed.querySelector('text').textContent, 'BULK_A');
  assert.deepEqual(harness.mutations, { commit: 1 });
  assert.deepEqual(messages, ['Loaded 2 row(s). Applied to 2 label(s).']);
});

// The display of a Result projects the current label intent, also when an
// Undo or another import removed the intent it was last shown with.
test('displaying a Result shows the current label intent when none remains (B6, R3)', () => {
  const svg = new DOMParser().parseFromString(
    `<svg>${batchLabel('fb', 'beta', 'dominant-baseline="central" ')}</svg>`, 'image/svg+xml'
  ).documentElement;
  svg.querySelector('text').textContent = 'IMPORTED_B';
  const harness = buildHarness({ svg, featureId: 'fb' });
  harness.actions.syncLabelEditor();
  assert.equal(svg.querySelector('text').textContent, 'beta');
  assert.deepEqual(harness.mutations, { commit: 1 });
});
