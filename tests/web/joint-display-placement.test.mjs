import assert from 'node:assert/strict';
import test from 'node:test';
import { readFile } from 'node:fs/promises';
import { gunzipSync } from 'node:zlib';
import {
  CANONICAL_REQUEST_SCHEMA, buildCanonicalRequestState, buildCanonicalRenderRequest, canonicalFeaturePlacements,
  projectCanonicalSessionRequest, promoteCanonicalRenderRequestToCurrent
} from '../../gbdraw/web/js/services/session-request.js';
import { createDefaultForm, createDefaultAdv } from '../../gbdraw/web/js/services/session-active-config-contract.js';
import {
  buildRecordDisplayRows, requestedRecordDisplay, validateRecordDisplayDrafts
} from '../../gbdraw/web/js/app/record-display-options.js';
import { createFeaturePlacementActions } from '../../gbdraw/web/js/app/feature-editor/placement-actions.js';
import { nameFeaturePlacementFailure } from '../../gbdraw/web/js/services/feature-placement.js';

import { resolveLinearComparisonPlan } from '../../gbdraw/web/js/app/linear-comparisons.js';

const target = (recordKey = 'card', biologicalFeatureId = 'feature') => ({
  recordKey, biologicalFeatureId, placement: { kind: 'main' }
});
const draft = (row, startCoordinate) => ({ scope: row.scope, sourceUid: row.sourceUid,
  selector: row.selector, recordId: row.recordId, topologyOverride: null, startCoordinate,
  reverseComplementOverride: null, anchorIntent: null });
const stateFor = (mode) => buildCanonicalRequestState({
  session: { renderRequest: { records: [] } },
  projection: { mode, inputType: 'gb', files: {}, config: {} },
  config: { form: createDefaultForm(), adv: createDefaultAdv(mode), colors: {} }
});
const text = 'LOCUS       same 100 bp DNA circular\n//\nLOCUS       same 100 bp DNA circular\n//\n';
const file = { name: 'same.gb', size: text.length, type: 'text/plain', lastModified: 0, data: btoa(text) };

test('schema 7 canonical placement uses code-point ordering and rejects transient or invalid wire intent', () => {
  const rows = [target('z'), target('A'), target('a'), target('\u{10000}'), target('\ue000')];
  assert.deepEqual(canonicalFeaturePlacements(rows, 'linear').map((row) => row.recordKey), ['A', 'a', 'z', '\ue000', '\u{10000}']);
  for (const placement of [{ kind: 'auto' }, { kind: 'lane', side: 'above', level: 2 },
    { kind: 'lane', side: 'outward', level: 1 }, { kind: 'main', feature_track_id: 0 }]) {
    assert.throws(() => canonicalFeaturePlacements([{ ...target(), placement }], 'linear'));
  }
  // A malformed row comes only from a Session file: a typed Session-format diagnostic (R6).
  const sessionFormat = (error) => error.code === 'INPUT_INVALID' && error.context.field === 'schema';
  assert.throws(() => canonicalFeaturePlacements([target(), target()], 'linear'), sessionFormat);
  assert.throws(() => canonicalFeaturePlacements({ 'card\0feature': target() }, 'linear'), sessionFormat);
});

test('inactive rotation retains raw intent but cannot emit an invalid effective display', () => {
  const row = buildRecordDisplayRows({ scope: 'linear', sourceUid: 'card', source: {},
    records: [{ recordId: 'same', recordLength: 100, detectedTopology: 'circular' }] })[0];
  const saved = { ...draft(row, 71), topologyOverride: false };
  validateRecordDisplayDrafts([saved]);
  assert.deepEqual(requestedRecordDisplay(row, saved), { isCircular: false, startCoordinate: null });
  assert.equal(requestedRecordDisplay(row, { ...saved, topologyOverride: null }).startCoordinate, 71);
  assert.equal(requestedRecordDisplay(row, draft(row, 1)).startCoordinate, 1);
  assert.equal(requestedRecordDisplay(row, draft(row, null)).startCoordinate, null);
  assert.throws(() => requestedRecordDisplay(row, draft(row, 101)), /between 1 and 100/);
  assert.throws(() => validateRecordDisplayDrafts([saved, saved]), /Duplicate/);
});

for (const mode of ['circular', 'linear']) {
  test(`${mode} shared writer preserves rotation placement and nondefault tolerance`, () => {
    const state = stateFor(mode);
    const sourceUid = mode === 'circular' ? 'circular' : 'card';
    const rows = buildRecordDisplayRows({ scope: mode, sourceUid, source: {},
      records: [1, 2].map(() => ({ recordId: 'same', recordLength: 100, detectedTopology: 'circular' })) });
    state.circularRecordList.value = rows.map((row) => ({ ...row, record_id: row.recordId }));
    state.recordDisplayDrafts = [draft(rows[0], 1), draft(rows[1], 71)];
    state.adv.feature_overlap_tolerance_bp = 2;
    const recordKey = mode === 'circular' ? 'record-2' : 'card:2';
    const override = target(recordKey);
    const draftRow = { scope: mode, ...override };
    state.featurePlacementOverrides = { [JSON.stringify([mode, recordKey, 'feature'])]: draftRow };
    const filesData = { c_gb: file, linearSeqs: [{ uid: 'card', gb: file }] };
    const result = buildCanonicalRenderRequest({ state, filesData, recordDisplayRows: rows, comparisonPlanSnapshot: mode === 'linear'
      ? resolveLinearComparisonPlan({ plan: state.linearComparisonPlan, sequences: filesData.linearSeqs, layout: [], losatProgram: 'blastn', blastpMode: 'orthogroup' }) : null });
    assert.equal(result.renderRequest.schema, CANONICAL_REQUEST_SCHEMA);
    assert.deepEqual(result.renderRequest.records.map((record) => record.display.startCoordinate), [1, 71]);
    assert.deepEqual(result.renderRequest.diagramOptions.featurePlacements, [override]);
    assert.equal(new Set(result.renderRequest.records.map((record) => record.source.resourceId)).size, 1);
    if (mode === 'linear') {
      assert.deepEqual(result.renderRequest.records.map((record) => record.presentation.gridRow), [1, 1]);
      assert.deepEqual(result.renderRequest.records.map((record) => record.cardinality), ['exactly_one', 'exactly_one']);
    }
    const projected = projectCanonicalSessionRequest(result);
    assert.equal(projected.config.adv.feature_overlap_tolerance_bp, 2);
    assert.deepEqual(Object.values(projected.config.featurePlacementOverrides), [draftRow]);
    const old = structuredClone(result.renderRequest);
    old.schema = 6;
    assert.throws(() => projectCanonicalSessionRequest({ ...result, renderRequest: old }), /schema 7/);
  });
}

// OV-08, R-1 (R2): each draft row names its mode, so the draft keeps both
// modes' rows and a request carries only the rows of its own mode and records,
// also when the other mode uses the same record key.
test('a request carries only its records\' placement rows and the draft keeps the other mode', () => {
  const row = (scope, recordKey, biologicalFeatureId, side) => ({ scope, recordKey, biologicalFeatureId,
    placement: side ? { kind: 'lane', side, level: 1 } : { kind: 'main' } });
  const rows = [row('circular', 'record-1', 'circular', 'outward'), row('linear', 'card', 'above', 'above'),
    row('linear', 'card', 'main'), row('circular', 'card', 'collision', 'inward'), row('circular', 'card', 'main'),
    row('linear', 'removed', 'gone', 'below')];
  const overrides = Object.fromEntries(rows.map((entry) => [
    JSON.stringify([entry.scope, entry.recordKey, entry.biologicalFeatureId]), entry]));
  assert.equal(canonicalFeaturePlacements(overrides).length, rows.length);
  // A draft lane is a side of its own mode.
  assert.throws(() => canonicalFeaturePlacements({ ...overrides,
    [JSON.stringify(['linear', 'card', 'wrong'])]: row('linear', 'card', 'wrong', 'inward') }),
  (error) => error.code === 'INPUT_INVALID' && error.context.field === 'schema');
  const state = stateFor('linear');
  state.featurePlacementOverrides = overrides;
  // Per-feature edits follow the same rule.
  const edit = (scope, featureVisibility) => ({ scope, recordKey: 'card', biologicalFeatureId: 'main',
    featureVisibility, labelVisibility: null, labelText: null, labelSourceText: null });
  state.featureOverrides = { [JSON.stringify(['circular', 'card', 'main'])]: edit('circular', 'off'),
    [JSON.stringify(['linear', 'card', 'main'])]: edit('linear', 'on') };
  const filesData = { linearSeqs: [{ uid: 'card', gb: file }] };
  const { renderRequest } = buildCanonicalRenderRequest({ state, filesData, comparisonPlanSnapshot:
    resolveLinearComparisonPlan({ plan: state.linearComparisonPlan, sequences: filesData.linearSeqs, layout: [], losatProgram: 'blastn', blastpMode: 'orthogroup' }) });
  const requestRow = ({ scope: _scope, ...entry }) => entry;
  assert.deepEqual(renderRequest.diagramOptions.featurePlacements, [requestRow(rows[1]), requestRow(rows[2])]);
  assert.deepEqual(renderRequest.diagramOptions.featureOverrides, [
    { recordKey: 'card', biologicalFeatureId: 'main', featureVisibility: 'on', labelVisibility: null, labelText: null }
  ]);
  assert.equal(Object.keys(state.featurePlacementOverrides).length, rows.length);
  assert.equal(Object.keys(state.featureOverrides).length, 2);
});

test('the Web names the feature of a placement row that Python reports (OV-09)', () => {
  const request = { mode: 'circular', diagramOptions: { featurePlacements: [
    { recordKey: 'record-1', biologicalFeatureId: 'b', placement: { kind: 'lane', side: 'outward', level: 1 } },
    { recordKey: 'record-1', biologicalFeatureId: 'a', placement: { kind: 'main' } }] } };
  const features = [{ record_key: 'record-1', biological_feature_id: 'b', type: 'CDS', product: 'NADH dehydrogenase subunit 1' }];
  const error = { code: 'FEATURE_PLACEMENT', context: { reason: 'SPLIT_LANES', placementIndex: 1 } };
  assert.deepEqual(nameFeaturePlacementFailure(error, request, features).context,
    { reason: 'SPLIT_LANES', placementIndex: 1, featureCaption: 'NADH dehydrogenase subunit 1' });
  for (const unnamed of [{ ...error, context: { placementIndex: 2 } }, { ...error, code: 'RENDER_FAILED' }]) {
    assert.equal(nameFeaturePlacementFailure(unnamed, request, features), unnamed);
  }
  assert.equal(nameFeaturePlacementFailure(error, request, []), error);
});

test('Main, resolved side, bulk Auto and history share one draft owner', async () => {
  const features = ['one', 'two'].map((id) => ({ scope: 'linear', record_key: 'card', biological_feature_id: id }));
  const overrides = {};
  const transactions = [];
  const state = { mode: { value: 'linear' }, selectedResultIndex: { value: 0 }, featurePlacementOverrides: overrides,
    form: { linear_track_layout: 'middle', separate_strands: false }, adv: {},
    featureCatalog: { value: { items: [{ recordKeys: ['card'] }] } },
    trackSlotResolvedGeometry: { value: { mode: 'linear', records: [{ recordIndex: 0, resultIndex: 0,
      featurePlacementTargets: [{ kind: 'main' }, { kind: 'lane', side: 'below', level: 1 }] }] } } };
  const actions = createFeaturePlacementActions({ state, isCurrentFeature: () => true,
    getCommittedRequest: () => ({ mode: 'linear', records: [{ recordKey: 'card' }], grouping: 'single' }),
    history: { runUndoable: async (label, fn) => { transactions.push({ label, before: structuredClone(overrides) }); fn(); } } });
  assert.equal(actions.choices(features).find((choice) => choice.value === 'above').enabled, true);
  await actions.setPlacement(features, 'below');
  assert.equal(Object.keys(overrides).length, 2);
  assert.equal(transactions.length, 1);
  await actions.setPlacement(features, 'auto');
  assert.deepEqual(overrides, {});
  assert.equal(transactions.length, 2);
  assert.equal(Object.keys(transactions[1].before).length, 2);
  state.form.separate_strands = true;
  assert.throws(() => actions.setPlacement(features, 'above'), /Unavailable/);
});

for (const mode of ['circular', 'linear']) {
  test(`${mode} placement admission follows draft slots and preserves artifact geometry`, () => {
    const feature = { scope: mode, record_key: 'record-1', biological_feature_id: 'feature' };
    const state = { mode: { value: mode }, form: createDefaultForm(), adv: createDefaultAdv(mode),
      featurePlacementOverrides: {}, featureCatalog: { value: { items: [{ recordKeys: ['record-1'] }] } },
      trackSlotResolvedGeometry: { value: { mode, records: [] } } };
    const geometry = structuredClone(state.trackSlotResolvedGeometry.value);
    const actions = createFeaturePlacementActions({ state, getCommittedRequest: () => ({ mode }),
      isCurrentFeature: (entry) => entry === feature, history: { runUndoable: (_label, fn) => fn() } });
    const sides = mode === 'circular' ? ['outward', 'inward'] : ['above', 'below'];
    const check = (enabled) => {
      for (const side of sides) {
        assert.equal(actions.choices([feature]).find((choice) => choice.value === side).enabled, enabled);
        if (enabled) {
          actions.setPlacement([feature], side);
          assert.equal(actions.valueFor(feature), side);
        } else {
          const before = structuredClone(state.featurePlacementOverrides);
          assert.throws(() => actions.setPlacement([feature], side), /current draft feature slot/);
          assert.deepEqual(state.featurePlacementOverrides, before);
        }
      }
      assert.deepEqual(state.trackSlotResolvedGeometry.value, geometry);
    };
    const layoutField = mode === 'circular' ? 'track_type' : 'linear_track_layout';
    for (const separate of [true, false]) {
      state.form.separate_strands = separate;
      for (const layout of mode === 'circular' ? ['middle', 'tuckin', 'spreadout'] : ['middle', 'above', 'below']) {
        state.form[layoutField] = layout;
        check(layout === 'middle' && (mode === 'circular' || !separate));
      }
    }
    // A custom slot's resolved direction takes precedence over the preset name.
    state.adv[`${mode}_track_slots_enabled`] = true;
    state.form[layoutField] = mode === 'circular' ? 'tuckin' : 'above';
    state.form.separate_strands = false;
    const slot = { id: 'features', renderer: 'features', enabled: true, side: 'overlay',
      params: mode === 'circular' ? { lane_direction: 'split' } : {} };
    state.adv[`${mode}_track_slots`] = [slot];
    state.adv[`${mode}_track_slots_axis_index`] = 0;
    check(true);
    state.form[layoutField] = 'middle';
    slot.side = mode === 'circular' ? 'inside' : 'below';
    slot.params = {};
    check(false);
    // Omitted side uses the same Axis resolution as request projection.
    delete slot.side;
    check(false);
    slot.enabled = false;
    assert.deepEqual(actions.choices([feature]).filter((entry) => entry.enabled).map((entry) => entry.value), ['auto']);
    assert.throws(() => actions.setPlacement([feature], 'main'), /Unavailable/);
    assert.ok(actions.choices([{ ...feature, record_key: 'unknown' }]).every((entry) => !entry.enabled));
  });
}

// Q3 (Owner, 2026-10-04): a layout edit that would leave lane placements
// undrawable asks first. Reset is one step that applies the edit and removes
// exactly those rows; cancel records nothing (OV-09, R10, R11). The stack
// editors' edits use the same transition (tests/web/track-layout-transition.test.mjs).
test('a layout edit that drops lanes asks before it resets those placements', async () => {
  const lane = (scope, recordKey, side) => ({ scope, recordKey, biologicalFeatureId: 'f',
    placement: { kind: 'lane', side, level: 1 } });
  const overrides = { outward: lane('circular', 'c', 'outward'), above: lane('linear', 'l', 'above'),
    main: { scope: 'circular', recordKey: 'c', biologicalFeatureId: 'g', placement: { kind: 'main' } } };
  const state = { mode: { value: 'circular' }, form: { ...createDefaultForm(), track_type: 'middle' },
    adv: createDefaultAdv('circular'), featurePlacementOverrides: overrides };
  const steps = [];
  const actions = createFeaturePlacementActions({ state, getCommittedRequest: () => null, isCurrentFeature: () => true,
    history: { runUndoable: async (label, fn) => { steps.push(label); fn(); } } });
  const select = (value, label) => ({ target: { type: 'select-one', value, labels: [{ textContent: ` ${label} ` }],
    selectedOptions: [{ text: `${value[0].toUpperCase()}${value.slice(1)}` }] } });
  const ask = (event, field) => {
    assert.equal(actions.changeLayoutSetting(event, field), false);
    return { ...actions.layoutChange };
  };

  assert.deepEqual(ask(select('spreadout', 'Track Preset'), 'track_type'),
    { open: true, count: 1, setting: 'Track Preset', value: 'Spreadout', scope: '' });
  assert.equal(await actions.resolveLayoutChange('cancel'), false);
  assert.equal(actions.layoutChange.open, false);
  assert.equal(state.form.track_type, 'middle');
  assert.deepEqual(Object.keys(overrides), ['outward', 'above', 'main']);
  assert.deepEqual(steps, []);

  ask(select('tuckin', 'Track Preset'), 'track_type');
  await actions.resolveLayoutChange('reset');
  assert.equal(state.form.track_type, 'tuckin');
  // The other mode's lane and the Main row stay (R2).
  assert.deepEqual(Object.keys(overrides), ['above', 'main']);
  assert.deepEqual(steps, ['Change setting and reset Feature placements']);
  // With no drawable lane to lose, the edit applies now and the control's
  // History adapter records it (R11).
  actions.changeLayoutSetting(select('middle', 'Track Preset'), 'track_type');
  assert.equal(actions.layoutChange.open, false);
  assert.equal(state.form.track_type, 'middle');
  assert.equal(steps.length, 1);

  state.mode.value = 'linear';
  state.adv = createDefaultAdv('linear');
  Object.assign(state.form, { linear_track_layout: 'middle', separate_strands: false });
  const strands = { target: { type: 'checkbox', checked: true, getAttribute: () => 'Separate Strands' } };
  assert.deepEqual(ask(strands, 'separate_strands'),
    { open: true, count: 1, setting: 'Separate Strands', value: 'On', scope: '' });
  // The checkbox shows the kept value while the dialog asks.
  assert.equal(strands.target.checked, false);
  await actions.resolveLayoutChange('reset');
  assert.equal(state.form.separate_strands, true);
  assert.deepEqual(Object.keys(overrides), ['main']);
});

test('historical session 40 schema 6 promotes without Generate and preserves cardinality', async () => {
  const bytes = await readFile('tests/fixtures/sessions/test_linear_cli_sidecar_reuses0.v40-schema6.json.gz');
  const session = JSON.parse(gunzipSync(bytes));
  assert.equal(session.version, 40);
  assert.equal(session.renderRequest.schema, 6);
  const before = structuredClone(session.renderRequest);
  const current = promoteCanonicalRenderRequestToCurrent(session.renderRequest);
  assert.equal(current.schema, CANONICAL_REQUEST_SCHEMA);
  assert.deepEqual(current.records.map((record) => record.cardinality), before.records.map((record) => record.cardinality));
  assert.ok(current.records.every((record) => record.display.isCircular === null && record.display.startCoordinate === null));
  assert.deepEqual(current.diagramOptions.featurePlacements, []);
  assert.deepEqual(session.renderRequest, before);
});
