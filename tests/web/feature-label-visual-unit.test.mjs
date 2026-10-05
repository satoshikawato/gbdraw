import assert from 'node:assert/strict';
import test from 'node:test';

import { createFeatureLabelActions } from '../../gbdraw/web/js/app/feature-editor/label-actions.js';
import { installFakeSvgDom } from './fake-svg-dom.mjs';

installFakeSvgDom();
globalThis.CSS ||= { escape: (value) => String(value) };

const ref = (value) => ({ value });

// Per-feature label and visibility edits are identity rows of the Result's
// mode (design Q4, R2).
const featureFor = (featureId, type = 'CDS') => ({
  svg_id: featureId, type, scope: 'linear', record_key: 'record-1', biological_feature_id: `bio-${featureId}`
});
const keyFor = (featureId) => JSON.stringify(['linear', 'record-1', `bio-${featureId}`]);
const rowFor = (featureId, fields = {}) => ({
  scope: 'linear',
  recordKey: 'record-1',
  biologicalFeatureId: `bio-${featureId}`,
  featureVisibility: null,
  labelVisibility: null,
  labelText: null,
  labelSourceText: null,
  ...fields
});
const setRow = (state, featureId, fields) => {
  state.featureOverrides[keyFor(featureId)] = rowFor(featureId, {
    ...(state.featureOverrides[keyFor(featureId)] || {}), ...fields
  });
};
const rowOf = (state, featureId) => state.featureOverrides[keyFor(featureId)] || null;

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
  rulePreparation = {},
  diagramOptions = null
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
    extractedFeatures: ref([featureFor(featureId)]),
    clickedFeature: ref({
      svg_id: featureId,
      feat: featureFor(featureId),
      labelText: '',
      labelSourceText: '',
      labelVisibility: 'default',
      hasEditableLabel: true
    }),
    labelTextScopeDialog: { show: false },
    hiddenLabelTextDialog: { show: false, featureId: '' },
    labelOnDialog: { show: false, reason: '', featureType: '' },
    featureOverrides: Object.fromEntries(Object.entries(visibilityOverrides).map(([id, mode]) => [
      keyFor(id), rowFor(id, { labelVisibility: mode })
    ])),
    labelTextBulkOverrides: {},
    labelOverrideBuildWarning: ref(''),
    autoLabelReflowEnabled: ref(false),
    labelReflowRequestSeq: ref(0),
    labelReflowForceRequestSeq: ref(0),
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
    rulePreparation,
    getCommittedRequest: () => (diagramOptions ? { diagramOptions } : null)
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
  assert.equal(rowOf(harness.state, featureId).labelText, 'Renamed');
  harness.actions.handleHiddenLabelTextChoice('text_only');
  assert.deepEqual({ ...harness.state.hiddenLabelTextDialog }, { show: false, featureId: '' });
  // "Keep hidden (apply text only)" keeps Label visibility unset (design Q4 6.2).
  assert.equal(rowOf(harness.state, featureId).labelVisibility, null);
  assert.equal(harness.state.labelReflowForceRequestSeq.value, 0);

  harness.state.clickedFeature.value.labelText = 'Renamed again';
  await harness.actions.updateClickedFeatureLabelText();
  assert.equal(harness.state.hiddenLabelTextDialog.show, true);
  harness.actions.handleHiddenLabelTextChoice('show');
  assert.equal(harness.state.hiddenLabelTextDialog.show, false);
  assert.equal(rowOf(harness.state, featureId).labelVisibility, 'on');
  assert.equal(harness.state.clickedFeature.value.labelVisibility, 'on');
  assert.equal(rowOf(harness.state, featureId).labelText, 'Renamed again');
  assert.equal(harness.state.labelReflowForceRequestSeq.value, 1);
  assert.equal(harness.state.labelReflowRequestSeq.value, 0);
});

// Owner decisions Q1 and Q2 (2026-10-04): the popup note names why the
// committed request draws no label, and On only where it can be drawn.
test('the popup note names why the diagram draws no label for the feature', () => {
  const featureId = 'feature:one/[a]';
  const diagramOptions = { featureShapes: { repeat_region: 'underlay' }, configOverrides: { 'labels.rendering': 'auto' } };
  const harness = buildHarness({ diagramOptions });
  const { state, actions } = harness;
  const clicked = state.clickedFeature.value;
  const hint = () => actions.clickedFeatureLabelHint.value;
  const none = 'This feature has no label in the current Result.';
  const onDrawable = ' Its Label visibility "On" applies when the label can be drawn.';

  clicked.labelKey = 'label-1';
  assert.equal(hint(), '');
  clicked.labelKey = '';
  assert.equal(hint(), `${none} Choose On to show it.`);
  clicked.labelVisibility = 'off';
  assert.equal(hint(), '');
  clicked.labelVisibility = 'on';
  setRow(state, featureId, { labelVisibility: 'on' });
  assert.equal(hint(), `${none}${onDrawable}`);

  clicked.feat = featureFor(featureId, 'repeat_region');
  assert.equal(hint(), `${none} Labels are not drawn for features drawn as "Underlay".${onDrawable}`);
  delete state.featureOverrides[keyFor(featureId)];
  clicked.labelVisibility = 'default';
  assert.equal(hint(), `${none} Labels are not drawn for features drawn as "Underlay".`);

  clicked.feat = featureFor(featureId);
  diagramOptions.configOverrides['labels.rendering'] = 'embedded_only';
  assert.equal(hint(), `${none} With "Label Rendering" = "Embedded Only", a label is drawn only when it fits inside its feature.`);

  // A hidden feature hides its label too, so the note shows with a label key.
  setRow(state, featureId, { featureVisibility: 'off', labelVisibility: 'on' });
  clicked.labelKey = 'label-1';
  assert.equal(hint(), `${none} The feature is hidden. Its Label visibility "On" applies when the feature is shown.`);
});

test('Label visibility On for an underlay feature waits for Keep without label or Cancel', async () => {
  const featureId = 'feature:one/[a]';
  const harness = buildHarness({ diagramOptions: { featureShapes: {}, configOverrides: {} } });
  const { state, actions } = harness;
  harness.state.editableLabels.value = [];
  Object.assign(state.clickedFeature.value, {
    feat: featureFor(featureId, 'repeat_region'), labelVisibility: 'on', labelText: 'RPT'
  });
  const canceled = actions.updateClickedFeatureLabelText();
  assert.deepEqual({ ...state.labelOnDialog }, { show: true, reason: 'underlay', featureType: 'repeat_region' });
  actions.handleLabelOnChoice('cancel');
  assert.equal(await canceled, false);
  assert.equal(state.labelOnDialog.show, false);
  assert.deepEqual(state.featureOverrides, {});
  assert.equal(state.clickedFeature.value.labelVisibility, 'default');

  Object.assign(state.clickedFeature.value, { labelVisibility: 'on', labelText: 'RPT' });
  const kept = actions.updateClickedFeatureLabelText();
  actions.handleLabelOnChoice('keep');
  await kept;
  assert.deepEqual(state.featureOverrides, {
    [keyFor(featureId)]: rowFor(featureId, { labelVisibility: 'on', labelText: 'RPT' })
  });
  assert.equal(state.labelReflowForceRequestSeq.value, 1);
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
    assert.equal(harness.state.labelReflowRequestSeq.value, 0);
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
  assert.equal(harness.state.labelReflowRequestSeq.value, 0);
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
  setRow(state, 'fb', { labelText: 'EDITED_B', labelSourceText: 'source', labelVisibility: 'off' });
  Object.assign(state.labelTextBulkOverrides, { source: 'BULK' });
  const before = JSON.stringify([state.featureOverrides, state.labelTextBulkOverrides]);
  for (const next of [resultB, resultA, resultB]) {
    mounted = next;
    actions.syncLabelEditor({ queueIncompleteVisibility: false });
    assert.equal(JSON.stringify([state.featureOverrides, state.labelTextBulkOverrides]), before);
  }
});

// F-3: Generate draws no label for a hidden feature. The label projection that
// live edits, History, and Result display share hides it with the feature.
test('a hidden feature hides its label through the stored visibility projection', () => {
  const featureId = 'feature:one/[a]';
  const harness = buildHarness();
  setRow(harness.state, featureId, { featureVisibility: 'off' });
  harness.state.featureVisibilityManualRules = [];
  harness.actions.reconcileLabelOverrides();
  exactParts(harness.svg, featureId).forEach((part) => {
    assert.equal(part.getAttribute('display'), 'none');
  });
  assert.equal(harness.svg.querySelector('[data-part="unrelated"]').getAttribute('display'), null);

  delete harness.state.featureOverrides[keyFor(featureId)];
  harness.actions.reconcileLabelOverrides();
  exactParts(harness.svg, featureId).forEach((part) => {
    assert.equal(part.getAttribute('display'), null);
  });
});

test('a feature visibility edit commits its label once and queues the label reflow', () => {
  const featureId = 'feature:one/[a]';
  const harness = buildHarness();
  harness.state.autoLabelReflowEnabled.value = true;
  setRow(harness.state, featureId, { featureVisibility: 'off' });
  assert.equal(harness.actions.applyFeatureVisibilityToLabels(), true);
  exactParts(harness.svg, featureId).forEach((part) => {
    assert.equal(part.getAttribute('display'), 'none');
  });
  assert.deepEqual(harness.mutations, { commit: 1 });
  assert.equal(harness.state.labelReflowRequestSeq.value, 1);
  assert.equal(harness.state.labelReflowForceRequestSeq.value, 0);

  delete harness.state.featureOverrides[keyFor(featureId)];
  assert.equal(harness.actions.applyFeatureVisibilityToLabels({ reflow: false }), true);
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
  state.results.value = [
    { name: 'a.svg', content: '<svg></svg>' },
    { name: 'b.svg', content: `<svg>${batchLabel('fb', 'beta')}</svg>` }
  ];
  // The other Result's labels reach their features through its catalog item.
  const item = (resultIndex, featureId) => ({
    resultIndex,
    resultName: `${'ab'[resultIndex]}.svg`,
    recordKeys: ['record-1'],
    features: [{ svgId: featureId, recordKey: 'record-1', biologicalFeatureId: `bio-${featureId}`, fillColor: '', drawnSelector: null }],
    biologicalFeatures: [{
      recordKey: 'record-1', biologicalFeatureId: `bio-${featureId}`, type: 'CDS',
      anchorProfile: { precision: 'exact', operator: 'single', partOrder: 'biological', strand: '+' }
    }],
    orthogroups: [],
    annotations: [],
    comparisonMatches: []
  });
  state.featureCatalog = ref({ schema: 5, items: [item(0, 'fa'), item(1, 'fb')] });
  const messages = [];
  globalThis.window = { alert: (message) => messages.push(message) };
  try {
    const tsv = '*\tCDS\tlabel\t^beta$\tIMPORTED_B\n*\t*\tlabel\t^alpha$\tBULK_A\n';
    await actions.loadLabelOverrideTable({ target: { files: [{ text: async () => tsv }], value: 'labels.tsv' } });
  } finally {
    delete globalThis.window;
  }
  assert.deepEqual(evaluated, [['alpha', 'beta']], 'one evaluation covers the labels of both Results');
  assert.deepEqual(state.featureOverrides, {
    [keyFor('fa')]: rowFor('fa', { labelSourceText: 'alpha' }),
    [keyFor('fb')]: rowFor('fb', { labelText: 'IMPORTED_B', labelSourceText: 'beta' })
  });
  assert.deepEqual(state.labelTextBulkOverrides, { alpha: 'BULK_A' });
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
