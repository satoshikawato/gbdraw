import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import test from 'node:test';

import { LABEL_ABSENCE_REASONS, createFeatureLabelActions } from '../../gbdraw/web/js/app/feature-editor/label-actions.js';
import { createHistoryManager } from '../../gbdraw/web/js/services/history.js';
import { installFakeSvgDom } from './fake-svg-dom.mjs';
import { withDrawings, withModeDrawings } from './helpers/drawing-state.mjs';

installFakeSvgDom();
globalThis.CSS ||= { escape: (value) => String(value) };

const ref = (value) => ({ value });

// Per-feature label and visibility edits are identity rows of the drawing of
// the Result's mode (design Q4, R2; PR-1: the drawing is the mode, so a row
// has no scope).
const featureFor = (featureId, type = 'CDS') => ({
  svg_id: featureId, type, scope: 'linear', record_key: 'record-1', biological_feature_id: `bio-${featureId}`
});
const keyFor = (featureId) => JSON.stringify(['record-1', `bio-${featureId}`]);
const rowFor = (featureId, fields = {}) => ({
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
  evaluateLabelRules = async () => ({ winners: [] }),
  diagramOptions = null,
  showNotice = () => {}
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
    hiddenLabelTextDialog: { show: false, featureId: '', reason: '' },
    labelOnDialog: { show: false, reason: '', featureType: '' },
    featureOverrides: Object.fromEntries(Object.entries(visibilityOverrides).map(([id, mode]) => [
      keyFor(id), rowFor(id, { labelVisibility: mode })
    ])),
    labelTextBulkOverrides: {},
    canonicalLabelOverrideRows: ref([]),
    autoLabelReflowEnabled: ref(false),
    labelReflowRequestSeq: ref(0),
    labelReflowForceRequestSeq: ref(0),
    labelReflowLastError: ref(null)
  };
  const actions = createFeatureLabelActions({
    ref, computed: get => ({ get value() { return get(); } }),
    state: withDrawings(state),
    commitActiveResultEdit(reason) {
      assert.equal(reason, 'feature-label');
      mutations.commit += 1;
      return true;
    },
    rulePreparation,
    evaluateLabelRules,
    showNotice,
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

// Label Not Shown asks before anything is written (Owner Q1, Q2, as On asks):
// the Apply stays open until the choice, and Cancel writes nothing.
test('a text edit asks whether to show a label only when the feature has none', async () => {
  const featureId = 'feature:one/[a]';
  const labeled = buildHarness();
  Object.assign(labeled.state.clickedFeature.value, { labelText: 'Renamed', labelSourceText: 'Original' });
  await labeled.actions.updateClickedFeatureLabelText();
  assert.equal(labeled.state.hiddenLabelTextDialog.show, false);

  const harness = buildHarness();
  harness.state.editableLabels.value = [];
  Object.assign(harness.state.clickedFeature.value, { labelText: 'Renamed', labelSourceText: 'Original' });
  const canceled = harness.actions.updateClickedFeatureLabelText();
  assert.deepEqual({ ...harness.state.hiddenLabelTextDialog }, { show: true, featureId, reason: '' });
  assert.equal(rowOf(harness.state, featureId)?.labelText ?? null, null);
  harness.actions.handleHiddenLabelTextChoice('cancel');
  await canceled;
  assert.equal(harness.state.hiddenLabelTextDialog.show, false);
  assert.equal(rowOf(harness.state, featureId)?.labelText ?? null, null);
  assert.equal(harness.state.labelReflowForceRequestSeq.value + harness.state.labelReflowRequestSeq.value, 0);

  const kept = harness.actions.updateClickedFeatureLabelText();
  assert.equal(harness.state.hiddenLabelTextDialog.show, true);
  harness.actions.handleHiddenLabelTextChoice('text_only');
  await kept;
  assert.deepEqual({ ...harness.state.hiddenLabelTextDialog }, { show: false, featureId: '', reason: '' });
  assert.equal(rowOf(harness.state, featureId).labelText, 'Renamed');
  // "Keep hidden (apply text only)" keeps Label visibility unset (design Q4 6.2).
  assert.equal(rowOf(harness.state, featureId).labelVisibility, null);
  assert.equal(harness.state.labelReflowForceRequestSeq.value, 0);

  harness.state.clickedFeature.value.labelText = 'Renamed again';
  const shown = harness.actions.updateClickedFeatureLabelText();
  assert.equal(harness.state.hiddenLabelTextDialog.show, true);
  harness.actions.handleHiddenLabelTextChoice('show');
  await shown;
  assert.equal(harness.state.hiddenLabelTextDialog.show, false);
  assert.equal(rowOf(harness.state, featureId).labelVisibility, 'on');
  assert.equal(harness.state.clickedFeature.value.labelVisibility, 'on');
  assert.equal(rowOf(harness.state, featureId).labelText, 'Renamed again');
  assert.equal(harness.state.labelReflowForceRequestSeq.value, 1);
  assert.equal(harness.state.labelReflowRequestSeq.value, 0);
});

// R-4 (FINDINGS, residual of #764): the Label Not Shown dialog and the popup
// note state the one reason that decides why the diagram draws no label.
const absenceCases = [
  ['hidden', { rows: { featureVisibility: 'off' } }, ' The feature is hidden.'],
  ['underlay', { type: 'repeat_region' }, ' Labels are not drawn for features drawn as "Underlay".'],
  ['scope_none', { scope: 'none' },
    ' Labels are set to "None" ("Show Labels" or "Label Mode"), so only a feature with Label visibility "On" has one.'],
  ['scope_first', { scope: 'first', recordIdx: 1 },
    ' "Show Labels" is "First Record Only", so labels are drawn only in the first record unless a feature has Label visibility "On".'],
  ['scope_orthogroup_top', { scope: 'orthogroup_top' },
    ' "Show Labels" is "Top Similarity Group Record", so labels are drawn only for the features that setting selects, unless a feature has Label visibility "On".'],
  ['whitelist', { whitelist: 'labelWhitelistFile' }, ' A label whitelist is set, and labels are drawn only for the features it lists.'],
  // The CLI and Python write the whitelist as `labelWhitelistTable`.
  ['whitelist', { whitelist: 'labelWhitelistTable' }, ' A label whitelist is set, and labels are drawn only for the features it lists.'],
  ['blacklist', { blacklist: ['putative'] },
    ' A label blacklist is set, and a label whose text contains one of its keywords is not drawn.'],
  ['embedded_only', { rendering: 'embedded_only' },
    ' With "Label Rendering" = "Embedded Only", a label is drawn only when it fits inside its feature.'],
  ['', {}, '']
];
const absenceHarness = ({ type = 'CDS', rows, scope = 'all', recordIdx = 0, whitelist, blacklist, rendering = 'auto' }) => {
  const featureId = 'feature:one/[a]';
  const feature = { ...featureFor(featureId, type), record_idx: recordIdx };
  const configOverrides = {
    'labels.rendering': rendering,
    'labels.circular.scope': scope,
    ...(blacklist ? { 'labels.filtering.blacklist_keywords': blacklist } : {})
  };
  const harness = buildHarness({
    diagramOptions: { featureShapes: {}, configOverrides, ...(whitelist ? { [whitelist]: { resourceId: 'w' } } : {}) }
  });
  harness.state.featureVisibilityManualRules = [];
  harness.state.editableLabels.value = [];
  if (rows) setRow(harness.state, featureId, rows);
  Object.assign(harness.state.clickedFeature.value, {
    feat: feature, labelText: 'Renamed', labelSourceText: 'Original'
  });
  harness.state.extractedFeatures.value = [feature];
  return { ...harness, feature, featureId };
};

for (const [reason, options, sentence] of absenceCases) {
  test(`the Label Not Shown dialog names the reason "${reason || 'unknown'}" as the popup note does`, async () => {
    const { actions, state } = absenceHarness(options);
    const none = 'This feature has no label in the current Result.';
    // The popup note, before the edit is applied, names the same reason in a few words.
    state.clickedFeature.value.labelText = '';
    const cause = LABEL_ABSENCE_REASONS[reason]?.cause;
    assert.equal(actions.clickedFeatureLabelHint.value.startsWith(`No label shown${cause ? ` (${cause})` : ''}.`), true,
      actions.clickedFeatureLabelHint.value);
    state.clickedFeature.value.labelText = 'Renamed';
    const applied = actions.updateClickedFeatureLabelText();
    assert.equal(state.hiddenLabelTextDialog.show, true);
    assert.equal(state.hiddenLabelTextDialog.reason, reason);
    const message = actions.hiddenLabelTextMessage.value;
    assert.equal(message.startsWith(`${none}${sentence || ' '}`), true, message);
    assert.match(message, /The edited text will not appear/);
    actions.handleHiddenLabelTextChoice('cancel');
    await applied;
  });
}

test('every reason the resolver can return has a sentence, and every sentence a reason', async () => {
  const returned = new Set();
  for (const [, options] of absenceCases) {
    const { actions, state } = absenceHarness(options);
    const applied = actions.updateClickedFeatureLabelText();
    returned.add(state.hiddenLabelTextDialog.reason);
    actions.handleHiddenLabelTextChoice('cancel');
    await applied;
  }
  returned.delete('');
  assert.deepEqual([...returned].sort(), Object.keys(LABEL_ABSENCE_REASONS).sort());
  Object.values(LABEL_ABSENCE_REASONS).forEach(({ sentence, cause }) => {
    assert.equal(typeof sentence, 'string');
    assert.equal(typeof cause, 'string');
  });
});

// The popup's Apply as app-setup.js runs it, in a step of the real History
// manager, up to the open Label Not Shown dialog.
const applyInHistoryStep = async (options = {}) => {
  const harness = buildHarness(options);
  harness.state.editableLabels.value = [];
  Object.assign(harness.state.clickedFeature.value, { labelText: 'Renamed', labelSourceText: 'Original' });
  const history = createHistoryManager({
    buildIntent: () => structuredClone(harness.state.featureOverrides),
    applyIntent: () => {}, buildCheckpoint: () => ({}), applyCheckpoint: () => {}
  });
  await history.initializeIntentBaseline();
  const step = { ended: false };
  // app-setup.js: `undoableAction('Change label text', updateClickedFeatureLabelText)`.
  step.done = history.runUndoable('Change label text', () => harness.actions.updateClickedFeatureLabelText())
    .then(() => { step.ended = true; });
  await new Promise((resolve) => setImmediate(resolve));
  assert.equal(harness.state.hiddenLabelTextDialog.show, true);
  assert.equal(history.mutationPending(), true, 'the step of the Apply is open while the dialog asks');
  return { ...harness, history, step };
};

// OV-462: the dialog can be hidden without a choice. Hiding it answers it as
// Cancel, so the History step of the Apply that opened it ends; before, the
// step stayed open (no Undo, Redo, Save or Load) until the next undoable action.
const loadedLabelSvg = () => new DOMParser().parseFromString([
  '<svg><path data-gbdraw-feature-id="fa" d="M0 0" />',
  '<text dominant-baseline="central" data-label-feature-id="fa" data-gbdraw-label-binding-schema="1" ',
  'data-label-source-text="alpha">alpha</text></svg>'
].join(''), 'image/svg+xml').documentElement;
for (const [hidden, hide, options] of [
  ['Reset label edits', (actions) => actions.resetAllLabelTextOverrides()],
  ['a Load Label TSV that applies', async (actions, state) => {
    await actions.loadLabelOverrideTable({ target: { files: [{ text: async () => '*\tCDS\tlabel\t^alpha$\tLOADED\n' }], value: 'labels.tsv' } });
    assert.equal(rowOf(state, 'fa').labelText, 'LOADED', 'the table was applied');
  }, {
    svg: loadedLabelSvg(), featureId: 'fa', rulePreparation: { snapshot: () => ({}), isCurrent: () => true },
    evaluateLabelRules: async ({ features }) => ({ winners: features.map(() => 0) })
  }],
  // The port a mode switch calls (watchers.js `resetModeTransientUi`).
  ['closing it from another owner', (actions) => actions.closeHiddenLabelTextDialog()]
]) {
  test(`${hidden} answers an open Label Not Shown dialog as Cancel, and its History step ends`, async () => {
    const { actions, state, history, step } = await applyInHistoryStep(options);
    await hide(actions, state);
    await new Promise((resolve) => setImmediate(resolve));
    assert.equal(step.ended, true, 'the Apply that opened the dialog ends');
    assert.deepEqual({ ...state.hiddenLabelTextDialog }, { show: false, featureId: '', reason: '' });
    assert.equal(history.mutationPending(), false);
    assert.notEqual(rowOf(state, state.clickedFeature.value.svg_id)?.labelText, 'Renamed', 'the text is not written');
  });
}

// The choice is a second `runUndoable` call on the step of the Apply. It joins
// that step (`joinsOpenStep` in app-setup.js): one step, under the Apply's name.
test('the Label Not Shown choice and the Apply that opened the dialog are one History step under the Apply\'s name', async () => {
  const wiring = readFileSync(new URL('../../gbdraw/web/js/app/app-setup.js', import.meta.url), 'utf8')
    .match(/const handleHiddenLabelTextChoiceWithHistory = undoableAction\('([^']+)', handleHiddenLabelTextChoice(, \{ joinsOpenStep: true \})?\);/);
  assert.ok(wiring, 'app-setup.js wires the choice through undoableAction');
  const [, label, joins] = wiring;
  const { actions, state, history, step } = await applyInHistoryStep();
  await history.runUndoable(label, () => actions.handleHiddenLabelTextChoice('text_only'), joins ? { joinsOpenStep: true } : {});
  await step.done;
  assert.equal(rowOf(state, 'feature:one/[a]').labelText, 'Renamed');
  assert.equal(history.getUndoCount(), 1);
  assert.equal(history.undoLabel(), 'Change label text');
});

// "Show this label" applies Label visibility On through the Q1 and Q2 dialogs.
for (const [reason, options] of [['underlay', { type: 'repeat_region' }], ['hidden', { rows: { featureVisibility: 'off' } }]]) {
  test(`"Show this label" for the reason "${reason}" asks the Label On question`, async () => {
    const { actions, state, featureId } = absenceHarness(options);
    const applied = actions.updateClickedFeatureLabelText();
    assert.equal(state.hiddenLabelTextDialog.reason, reason);
    const shown = actions.handleHiddenLabelTextChoice('show');
    await Promise.resolve();
    await Promise.resolve();
    assert.equal(state.labelOnDialog.show, true);
    assert.equal(state.labelOnDialog.reason, reason);
    actions.handleLabelOnChoice('cancel');
    assert.equal(await shown, false);
    await applied;
    assert.equal(rowOf(state, featureId)?.labelVisibility ?? null, null);
    assert.equal(rowOf(state, featureId)?.labelText ?? null, null);
  });
}

// Owner decisions Q1 and Q2 (2026-10-04): the popup note names why the
// committed request draws no label, and On only where it can be drawn, in one
// line.
test('the popup note names why the diagram draws no label for the feature', () => {
  const featureId = 'feature:one/[a]';
  const diagramOptions = { featureShapes: { repeat_region: 'underlay' }, configOverrides: { 'labels.rendering': 'auto' } };
  const harness = buildHarness({ diagramOptions });
  const { state, actions } = harness;
  const clicked = state.clickedFeature.value;
  const hint = () => actions.clickedFeatureLabelHint.value;
  const onDrawable = ' On applies when the label can be drawn.';

  clicked.labelKey = 'label-1';
  assert.equal(hint(), '');
  clicked.labelKey = '';
  assert.equal(hint(), 'No label shown. Choose On to show it.');
  clicked.labelVisibility = 'off';
  assert.equal(hint(), '');
  clicked.labelVisibility = 'on';
  setRow(state, featureId, { labelVisibility: 'on' });
  assert.equal(hint(), `No label shown.${onDrawable}`);

  clicked.feat = featureFor(featureId, 'repeat_region');
  assert.equal(hint(), `No label shown (drawn as Underlay).${onDrawable}`);
  delete state.featureOverrides[keyFor(featureId)];
  clicked.labelVisibility = 'default';
  assert.equal(hint(), 'No label shown (drawn as Underlay).');

  clicked.feat = featureFor(featureId);
  diagramOptions.configOverrides['labels.rendering'] = 'embedded_only';
  assert.equal(hint(), 'No label shown (does not fit; Label Rendering = Embedded Only).');

  // A hidden feature hides its label too, so the note shows with a label key.
  setRow(state, featureId, { featureVisibility: 'off', labelVisibility: 'on' });
  clicked.labelKey = 'label-1';
  assert.equal(hint(), 'No label shown (the feature is hidden). On applies when the feature is shown.');
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

// U2BFIX2 review M1 (OV-200): with Embedded Only only Python decides whether a
// label text fits, so a projected text asks for the rerender; a label the
// projection hides needs none. A loaded table projects the feature visibility
// and the labels in one follow, which queues one request.
const EMBEDDED_ONLY = { configOverrides: { 'labels.rendering': 'embedded_only' } };
test('with Embedded Only, a projected label text asks for one rerender and a hidden label for none', () => {
  const featureId = 'feature:one/[a]';
  const hidden = buildHarness({ diagramOptions: EMBEDDED_ONLY });
  setRow(hidden.state, featureId, { labelVisibility: 'off' });
  assert.equal(hidden.actions.reconcileLabelOverrides(), true);
  assert.equal(hidden.state.labelReflowForceRequestSeq.value, 0, 'a hidden label');

  const edited = buildHarness({ diagramOptions: EMBEDDED_ONLY });
  const { state, actions } = edited;
  setRow(state, featureId, { labelText: 'edited' });
  assert.equal(actions.reconcileLabelOverrides(), true);
  assert.equal(state.labelReflowForceRequestSeq.value, 1, 'a changed text');
  state.autoLabelReflowEnabled.value = true;
  setRow(state, featureId, { labelText: 'edited again' });
  assert.equal(actions.applyFeatureVisibilityToLabels({ reflow: true, labels: true }), true);
  assert.deepEqual([state.labelReflowForceRequestSeq.value, state.labelReflowRequestSeq.value], [2, 0],
    'one request for the loaded table');
});

// U2BFIX2 review L3 (OV-237): the labels of a Result no popup bound are bound
// once per mounted root; a later projection (a History step) does not bind
// them again.
test('the label projection binds the labels of a mounted root once', () => {
  const svg = new DOMParser().parseFromString(
    '<svg><path data-gbdraw-feature-id="fa" d="M0 0" /></svg>', 'image/svg+xml'
  ).documentElement;
  const harness = buildHarness({ svg, featureId: 'fa' });
  /** @type {string[]} */
  const scans = [];
  const query = svg.querySelectorAll.bind(svg);
  svg.querySelectorAll = (selector) => { scans.push(selector); return query(selector); };
  harness.actions.reconcileLabelOverrides();
  harness.actions.reconcileLabelOverrides();
  assert.equal(scans.filter((selector) => selector.includes('dominant-baseline')).length, 1);
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
  // The outcome shows in the app's notice banner, not a native alert.
  const messages = [];
  const harness = buildHarness({
    svg: displayed,
    featureId: 'fa',
    showNotice: (notice) => messages.push(notice),
    rulePreparation: { snapshot: () => ({}), isCurrent: () => true },
    async evaluateLabelRules({ features, rules }) {
      evaluated.push(features.map((feature) => feature.label));
      return { winners: features.map((feature) => rules.findIndex((rule) => (
        new RegExp(rule.valueRegex).test(feature.label)
      ))) };
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
  const tsv = '*\tCDS\tlabel\t^beta$\tIMPORTED_B\n*\t*\tlabel\t^alpha$\tBULK_A\n';
  await actions.loadLabelOverrideTable({ target: { files: [{ text: async () => tsv }], value: 'labels.tsv' } });
  assert.deepEqual(evaluated, [['alpha', 'beta']], 'one evaluation covers the labels of both Results');
  assert.deepEqual(state.featureOverrides, {
    [keyFor('fa')]: rowFor('fa', { labelSourceText: 'alpha' }),
    [keyFor('fb')]: rowFor('fb', { labelText: 'IMPORTED_B', labelSourceText: 'beta' })
  });
  assert.deepEqual(state.labelTextBulkOverrides, { alpha: 'BULK_A' });
  assert.equal(displayed.querySelector('text').textContent, 'BULK_A');
  assert.deepEqual(harness.mutations, { commit: 1 });
  assert.deepEqual(messages, [{ text: 'Loaded 2 row(s). Applied to 2 label(s).', items: [] }]);
});

// OV-131: a Label TSV that applies to no label changes nothing. The import
// replaces the label edits only when it applies, so the bulk and per-feature
// label edits stay, no state change reaches History, and the message says why.
const KEPT = 'Loaded 1 row(s). Not applied; the existing label edits were kept.';
for (const [name, body, tsv, message] of [
  ['no row matches a label', batchLabel('fa', 'alpha', 'dominant-baseline="central" '),
    '*\t*\tproduct\t^no such product$\tZZ\n',
    { text: KEPT, items: ['No row matched a label of the diagram.'] }],
  ['the displayed Result has no label', '', '*\t*\tproduct\t.*\tZZ\n',
    { text: KEPT, items: ['No row matched a label of the diagram.'] }],
  ['every match lacks a feature key',
    '<text dominant-baseline="central" data-label-source-text="alpha">alpha</text>',
    '*\tCDS\tlabel\t^alpha$\tZZ\n',
    { text: KEPT, items: ['The 1 matched label(s) lacked a feature key.'] }]
]) {
  test(`a Label TSV import keeps the label edits when ${name} (OV-131)`, async () => {
    const svg = new DOMParser().parseFromString(`<svg>${body}</svg>`, 'image/svg+xml').documentElement;
    const messages = [];
    const harness = buildHarness({
      svg,
      featureId: 'fa',
      showNotice: (notice) => messages.push(notice),
      rulePreparation: { snapshot: () => ({}), isCurrent: () => true },
      async evaluateLabelRules({ features, rules }) {
        return { winners: features.map((feature) => rules.findIndex((rule) => (
          rule.qualifier === 'label' && new RegExp(rule.valueRegex).test(feature.label)
        ))) };
      }
    });
    const { state, actions } = harness;
    state.errorLog = ref(null);
    state.labelTextBulkOverrides.ND1 = 'MyLabel';
    setRow(state, 'fa', { labelText: 'EDITED', labelSourceText: 'alpha' });
    const before = structuredClone({ rows: state.featureOverrides, bulk: state.labelTextBulkOverrides });
    await actions.loadLabelOverrideTable({ target: { files: [{ text: async () => tsv }], value: 'labels.tsv' } });
    assert.deepEqual({ rows: state.featureOverrides, bulk: state.labelTextBulkOverrides }, before);
    assert.deepEqual(harness.mutations, { commit: 0 });
    assert.deepEqual(messages, [message]);
  });
}

// OV-1102: an action that changes the bulk label edits replaces the saved label
// table (the rules a loaded Session sent) itself, in its own History step, so
// Undo returns both. An action that leaves the bulk label edits as they are
// keeps the table: a per-feature label edit applies before it (design Q4).
const SAVED_LABEL_ROW = Object.freeze({
  recordId: '*', featureType: '*', qualifier: 'product', valueRegex: '^alpha$', labelText: 'SAVED'
});
for (const [name, tsv, bulk, rows] of [
  ['a global label row replaces it', '*\t*\tlabel\t^alpha$\tBULK_A\n', { alpha: 'BULK_A' }, []],
  ['a per-feature row keeps it', '*\tCDS\tlabel\t^alpha$\tLOADED\n', {}, [SAVED_LABEL_ROW]]
]) {
  test(`a Label TSV load: ${name} (OV-1102)`, async () => {
    const harness = buildHarness({
      svg: loadedLabelSvg(), featureId: 'fa', rulePreparation: { snapshot: () => ({}), isCurrent: () => true },
      evaluateLabelRules: async ({ features }) => ({ winners: features.map(() => 0) })
    });
    const { state, actions } = harness;
    state.errorLog = ref(null);
    state.canonicalLabelOverrideRows.value = [SAVED_LABEL_ROW];
    await actions.loadLabelOverrideTable({ target: { files: [{ text: async () => tsv }], value: 'labels.tsv' } });
    assert.deepEqual([state.labelTextBulkOverrides, state.canonicalLabelOverrideRows.value], [bulk, rows]);
  });
}
test('Reset label edits and an "all labels" text edit replace the saved label table when they change the bulk label edits (OV-1102)', () => {
  const svg = new DOMParser().parseFromString(
    '<svg><text data-label-editable="true" data-label-source-text="alpha">alpha</text></svg>', 'image/svg+xml'
  ).documentElement;
  const { state, actions } = buildHarness({ svg });
  state.canonicalLabelOverrideRows.value = [SAVED_LABEL_ROW];
  actions.resetAllLabelTextOverrides();
  assert.deepEqual(state.canonicalLabelOverrideRows.value, [SAVED_LABEL_ROW], 'no bulk label edit changed');
  Object.assign(state.labelTextScopeDialog, { show: true, sourceText: 'alpha', newText: 'ALL', labelKey: '', featureId: '' });
  actions.handleLabelTextScopeChoice('all');
  assert.deepEqual([state.labelTextBulkOverrides, state.canonicalLabelOverrideRows.value], [{ alpha: 'ALL' }, []]);
  state.canonicalLabelOverrideRows.value = [SAVED_LABEL_ROW];
  actions.resetAllLabelTextOverrides();
  assert.deepEqual([state.labelTextBulkOverrides, state.canonicalLabelOverrideRows.value], [{}, []]);
});

// OV-447, OV-449 (Owner decisions): a Label TSV load is declined with a
// notice that names the cause when, while the table is being evaluated, a
// label is edited (the edit stays), the displayed diagram is replaced, or
// another input of the rule preparation changes. A label edit in the app also
// makes the rule preparation stale (a new feature row); its own notice is
// shown. A load a newer load replaced ends without a notice.
const heldLabelImport = ({ preparationFollowsRows = false } = {}) => {
  const svg = new DOMParser().parseFromString(
    `<svg>${batchLabel('fa', 'alpha', 'dominant-baseline="central" ')}</svg>`, 'image/svg+xml'
  ).documentElement;
  svg.querySelector('text').textContent = 'alpha';
  const evaluations = [];
  const messages = [];
  const preparation = { inputs: 0 };
  const snapshot = () => `${preparation.inputs}:${preparationFollowsRows ? Object.keys(harness.state.featureOverrides).length : 0}`;
  const harness = buildHarness({
    svg,
    featureId: 'fa',
    showNotice: (notice) => messages.push(notice),
    rulePreparation: { snapshot, isCurrent: (before) => before === snapshot() },
    evaluateLabelRules: ({ features }) => new Promise((resolve) => {
      evaluations.push(() => resolve({ winners: features.map(() => 0) }));
    })
  });
  harness.state.errorLog = ref(null);
  const target = () => ({ files: [{ text: async () => '*\t*\tlabel\t^alpha$\tFROM_TABLE\n' }], value: 'labels.tsv' });
  const load = async (input, during) => {
    const pending = harness.actions.loadLabelOverrideTable({ target: input });
    while (!evaluations.length) await new Promise((resolve) => setTimeout(resolve, 0));
    await during?.();
    evaluations.splice(0).forEach((release) => release());
    await pending;
  };
  return { harness, messages, target, load, preparation };
};
for (const preparationFollowsRows of [false, true]) {
  test(`a Label TSV load is declined with a notice when a label is edited while it is evaluated (OV-447, preparation ${preparationFollowsRows ? 'stale' : 'current'})`, async () => {
    const { harness, messages, target, load } = heldLabelImport({ preparationFollowsRows });
    const { state } = harness;
    const input = target();
    await load(input, () => { setRow(state, 'fa', { labelText: 'EDITED', labelSourceText: 'alpha' }); });
    assert.deepEqual(messages, [{ text: 'Loaded 1 row(s). Not applied: label text, label visibility, or feature visibility changed while '
      + 'the table loaded. The changes were kept. Load the table again.', tone: 'warning' }]);
    assert.deepEqual(state.featureOverrides, { [keyFor('fa')]: rowFor('fa', { labelText: 'EDITED', labelSourceText: 'alpha' }) });
    assert.deepEqual(state.labelTextBulkOverrides, {});
    assert.deepEqual(harness.mutations, { commit: 0 }, 'the declined load reaches no History step');
    assert.equal(input.value, '', 'the same file can be picked again');
  });
}
for (const [cause, change, notice] of [
  ['the displayed diagram is replaced', ({ harness }) => {
    const next = new DOMParser().parseFromString(`<svg>${batchLabel('fa', 'alpha')}</svg>`, 'image/svg+xml').documentElement;
    harness.state.svgContainer.value = { querySelector: (selector) => (selector === 'svg' ? next : null) };
  }, 'the displayed diagram changed while the table loaded'],
  ['a color rule or another rule input changes', ({ preparation }) => { preparation.inputs += 1; },
    "colors, color rules, the Legend, or the diagram's inputs changed while the table loaded"]
]) {
  test(`a Label TSV load is declined with a notice when ${cause} while it is evaluated (OV-449)`, async () => {
    const held = heldLabelImport();
    const { harness, messages, target, load } = held;
    const input = target();
    await load(input, () => change(held));
    assert.deepEqual(messages, [{ text: `Loaded 1 row(s). Not applied: ${notice}. Load the table again.`, tone: 'warning' }]);
    assert.deepEqual([harness.state.featureOverrides, harness.state.labelTextBulkOverrides], [{}, {}]);
    assert.equal(harness.svg.querySelector('text').textContent, 'alpha');
    assert.deepEqual(harness.mutations, { commit: 0 }, 'the declined load reaches no History step');
    assert.equal(input.value, '', 'the same file can be picked again');
  });
}
test('a Label TSV load a newer load replaced ends without a notice, also after a label edit (OV-447)', async () => {
  const { harness, messages, target, load } = heldLabelImport();
  const { state, actions } = harness;
  await load(target(), () => {
    setRow(state, 'fa', { labelText: 'EDITED', labelSourceText: 'alpha' });
    actions.loadLabelOverrideTable({ target: { files: [{ text: () => new Promise(() => {}) }], value: 'newer.tsv' } });
  });
  assert.deepEqual(messages, []);
  assert.deepEqual(harness.mutations, { commit: 0 });
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

// A label's source text is free text, so the bulk text table is read by own
// key: a label that reads `constructor` is not replaced by Object.prototype.
test('a label whose text names an Object.prototype member keeps its text', () => {
  const svg = labelSvg();
  const text = svg.querySelector('text');
  text.textContent = 'constructor';
  const harness = buildHarness({ svg });
  harness.actions.reconcileLabelOverrides();
  assert.equal(text.textContent, 'constructor');
  assert.deepEqual(harness.mutations, { commit: 0 });
});

// D-36, review C S3: a requested redraw is owed by the drawing it was asked
// for (each mode has its own, PD-OI-086), so a rerender shown in one mode does
// not pay the other mode's.
test('a requested redraw stays owed until a rerender of its own drawing is shown', () => {
  const harness = buildHarness();
  const circular = {};
  const linear = {};
  withModeDrawings(harness.state, { circular, linear });
  assert.equal(harness.actions.requestAutomaticRerender(), true);
  assert.equal(harness.actions.isRedrawOwed(linear), true, 'asked in the shown mode');
  assert.equal(harness.actions.isRedrawOwed(circular), false, 'the other mode owes nothing');
  harness.actions.redrawShown(circular);
  assert.equal(harness.actions.isRedrawOwed(linear), true, 'the other drawing\'s rerender');
  harness.actions.redrawShown(linear);
  assert.equal(harness.actions.isRedrawOwed(linear), false);
});

// OV-923: opening the editor or a feature popup binds the labels of the
// mounted Result, measuring every feature. While nothing the binding reads
// has changed on that root since the last binding, the binding is reused; a
// changed text, child, or geometry attribute binds again.
test('opening the editor reuses a current label binding until the root changes (OV-923)', () => {
  const observers = [];
  const previous = globalThis.MutationObserver;
  globalThis.MutationObserver = class {
    constructor(callback) { this.callback = callback; this.pending = []; this.connected = false; observers.push(this); }
    observe(target, options) { this.target = target; this.options = options; this.connected = true; }
    disconnect() { this.connected = false; }
    takeRecords() { const records = this.pending; this.pending = []; return records; }
    mutate(record) { if (this.connected) this.pending.push(record); }
    deliver() { const records = this.takeRecords(); if (records.length) this.callback(records, this); }
  };
  try {
    const svg = new DOMParser().parseFromString([
      '<svg>',
      '<g id="record_1">',
      '<path data-gbdraw-feature-id="feature-one" data-fake-bbox="0 0 10 10" />',
      '<path data-gbdraw-feature-id="feature-two" data-fake-bbox="50 0 10 10" />',
      '<text dominant-baseline="central" data-fake-bbox="0 12 10 4"></text>',
      '</g>',
      '</svg>'
    ].join(''), 'image/svg+xml').documentElement;
    const harness = buildHarness({ svg });
    const proto = Object.getPrototypeOf(svg.querySelector('path'));
    const getBBox = proto.getBBox;
    let measured = 0;
    proto.getBBox = function countedGetBBox() { measured += 1; return getBBox.call(this); };
    try {
      harness.actions.syncLabelEditor({ reuseCurrentBinding: true });
      const firstBinding = measured;
      assert.ok(firstBinding > 0, 'the first opening binds the labels');
      const label = svg.querySelector('text');
      const bound = [label.getAttribute('data-label-key'), label.getAttribute('data-label-feature-id')];
      const observer = observers.at(-1);
      assert.equal(observer.options.subtree, true);
      assert.equal(observer.options.characterData, true);
      assert.equal(observer.options.childList, true);
      assert.ok(observer.options.attributeFilter.includes('transform'));
      assert.equal(observer.options.attributeFilter.includes('style'), false, 'a hover highlight keeps the binding');

      harness.actions.syncLabelEditor({ reuseCurrentBinding: true });
      harness.actions.syncLabelEditor({ reuseCurrentBinding: true });
      assert.equal(measured, firstBinding, 'reopening measures nothing');
      assert.deepEqual([label.getAttribute('data-label-key'), label.getAttribute('data-label-feature-id')], bound);

      // A sync without the reuse (an edit, a Result display) always binds.
      harness.actions.syncLabelEditor();
      assert.equal(measured, firstBinding * 2);

      // A mutation the binding reads, delivered or still queued, binds again.
      observers.at(-1).mutate({ type: 'characterData' });
      harness.actions.syncLabelEditor({ reuseCurrentBinding: true });
      assert.equal(measured, firstBinding * 3);
      observers.at(-1).mutate({ type: 'attributes', attributeName: 'transform' });
      observers.at(-1).deliver();
      harness.actions.syncLabelEditor({ reuseCurrentBinding: true });
      assert.equal(measured, firstBinding * 4);
      harness.actions.syncLabelEditor({ reuseCurrentBinding: true });
      assert.equal(measured, firstBinding * 4);
    } finally {
      proto.getBBox = getBBox;
    }
  } finally {
    globalThis.MutationObserver = previous;
  }
});

// OV-1291: a Linear Result draws other text with the central baseline of a
// feature label: a region label in its annotation's group, the copy drawn in
// front of the leader lines, and the depth and GC content axis ticks. None of
// them is a feature label, so the editor neither lists it nor binds it to the
// unlabeled feature next to it.
test('a Linear region label or axis tick is never an editable feature label (OV-1291)', () => {
  const central = (id, bbox) => `<text id="${id}" dominant-baseline="central" data-fake-bbox="${bbox}"></text>`;
  const svg = new DOMParser().parseFromString([
    '<svg>',
    '<g id="record_group_1">',
    '<path data-gbdraw-feature-id="labeled" data-fake-bbox="0 0 10 10" />',
    '<path data-gbdraw-feature-id="unlabeled" data-fake-bbox="50 0 10 10" />',
    central('feature-label', '0 12 10 4'),
    '</g>',
    '<g id="gbdraw-annotation-track-regions_1-1">',
    '<g id="gbdraw-annotation-regions_1-rrna" data-gbdraw-annotation-id="rrna">',
    '<rect data-fake-bbox="50 20 10 4" />', central('region-label', '50 12 10 4'),
    '</g>',
    '</g>',
    '<g id="gbdraw-annotation-track-regions_1-1-front">',
    '<g data-gbdraw-region-label-for="gbdraw-annotation-regions_1-control">', central('front-region-label', '52 12 6 4'), '</g>',
    '</g>',
    '<g id="depth"><g id="depth_axis">', central('depth-tick', '48 14 2 2'), '</g></g>',
    '<g id="gc_content" data-gbdraw-slot-renderer="dinucleotide_content"><g id="gc_content_axis">',
    central('gc-tick', '48 16 2 2'), '</g></g>',
    '</svg>'
  ].join(''), 'image/svg+xml').documentElement;
  const harness = buildHarness({ svg, featureId: 'labeled' });
  harness.state.extractedFeatures.value = [featureFor('labeled'), featureFor('unlabeled')];
  harness.actions.syncLabelEditor();
  const byId = (id) => svg.querySelector(`#${id}`);
  assert.equal(byId('feature-label').getAttribute('data-label-editable'), 'true');
  assert.equal(byId('feature-label').getAttribute('data-label-feature-id'), 'labeled');
  for (const id of ['region-label', 'front-region-label', 'depth-tick', 'gc-tick']) {
    assert.equal(byId(id).getAttribute('data-label-editable'), null, `${id} is not an editable label`);
    assert.equal(byId(id).getAttribute('data-label-feature-id'), null, `${id} is bound to no feature`);
  }
  assert.deepEqual(harness.state.editableLabels.value.map((entry) => entry.featureId), ['labeled']);
});

test('a curved Circular region label is never an editable feature label (OV-1210, OV-1291)', () => {
  const curved = (id, bbox) => `<text id="${id}" data-fake-bbox="${bbox}"><textPath href="#${id}-path">label</textPath></text>`;
  const svg = new DOMParser().parseFromString([
    '<svg>',
    '<g id="record_group_1">',
    '<path data-gbdraw-feature-id="labeled" data-fake-bbox="0 0 10 10" />',
    '<path data-gbdraw-feature-id="unlabeled" data-fake-bbox="50 0 10 10" />',
    curved('feature-label', '2 3 6 4'),
    '</g>',
    '<g id="gbdraw-annotation-track-regions_1-1">',
    '<g id="gbdraw-annotation-regions_1-rrna" data-gbdraw-annotation-id="rrna">',
    '<path data-fake-bbox="50 20 10 4" />', curved('region-label', '50 12 10 4'),
    '</g>',
    '</g>',
    '</svg>'
  ].join(''), 'image/svg+xml').documentElement;
  const harness = buildHarness({ svg, featureId: 'labeled' });
  harness.state.mode.value = 'circular';
  harness.state.generatedMode.value = 'circular';
  harness.state.extractedFeatures.value = [featureFor('labeled'), featureFor('unlabeled')];
  harness.actions.syncLabelEditor();
  const byId = (id) => svg.querySelector(`#${id}`);
  assert.equal(byId('feature-label').getAttribute('data-label-editable'), 'true');
  assert.equal(byId('feature-label').getAttribute('data-label-feature-id'), 'labeled');
  assert.equal(byId('region-label').getAttribute('data-label-editable'), null);
  assert.equal(byId('region-label').getAttribute('data-label-feature-id'), null);
  assert.deepEqual(harness.state.editableLabels.value.map((entry) => entry.featureId), ['labeled']);
});
