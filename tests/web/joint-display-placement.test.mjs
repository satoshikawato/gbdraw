import assert from 'node:assert/strict';
import test from 'node:test';
import { readFile } from 'node:fs/promises';
import { gunzipSync } from 'node:zlib';
import {
  buildCanonicalRequestState, buildCanonicalRenderRequest, canonicalFeaturePlacements,
  projectCanonicalSessionRequest, promoteCanonicalRenderRequestToCurrent
} from '../../gbdraw/web/js/services/session-request.js';
import { createDefaultForm, createDefaultAdv } from '../../gbdraw/web/js/services/session-active-config-contract.js';
import {
  buildRecordDisplayRows, requestedRecordDisplay, validateRecordDisplayDrafts
} from '../../gbdraw/web/js/app/record-display-options.js';
import { createFeaturePlacementActions } from '../../gbdraw/web/js/app/feature-editor/placement-actions.js';

import { resolveLinearComparisonPlan } from '../../gbdraw/web/js/app/linear-comparisons.js';

const target = (recordKey = 'card', biologicalFeatureId = 'feature') => ({
  recordKey, biologicalFeatureId, placement: { kind: 'main' }
});
const draft = (row, startCoordinate) => ({ scope: row.scope, sourceUid: row.sourceUid,
  selector: row.selector, recordId: row.recordId, topologyOverride: null, startCoordinate });
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
  assert.throws(() => canonicalFeaturePlacements([target(), target()], 'linear'), /Duplicate/);
  assert.throws(() => canonicalFeaturePlacements({ 'card\0feature': target() }, 'linear'), /JSON pair/);
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
    state.recordDisplayRows = { value: rows };
    state.recordDisplayDrafts = [draft(rows[0], 1), draft(rows[1], 71)];
    state.adv.feature_overlap_tolerance_bp = 2;
    const recordKey = mode === 'circular' ? 'record-2' : 'card:2';
    const override = target(recordKey);
    state.featurePlacementOverrides = { [JSON.stringify([recordKey, 'feature'])]: override };
    const filesData = { c_gb: file, linearSeqs: [{ uid: 'card', gb: file }] };
    const result = buildCanonicalRenderRequest({ state, filesData, comparisonPlanSnapshot: mode === 'linear'
      ? resolveLinearComparisonPlan({ plan: state.linearComparisonPlan, sequences: filesData.linearSeqs, layout: [], losatProgram: 'blastn', blastpMode: 'orthogroup' }) : null });
    assert.equal(result.renderRequest.schema, 7);
    assert.deepEqual(result.renderRequest.records.map((record) => record.display.startCoordinate), [1, 71]);
    assert.deepEqual(result.renderRequest.diagramOptions.featurePlacements, [override]);
    assert.equal(new Set(result.renderRequest.records.map((record) => record.source.resourceId)).size, 1);
    if (mode === 'linear') {
      assert.deepEqual(result.renderRequest.records.map((record) => record.presentation.gridRow), [1, 1]);
      assert.deepEqual(result.renderRequest.records.map((record) => record.cardinality), ['exactly_one', 'exactly_one']);
    }
    const projected = projectCanonicalSessionRequest(result);
    assert.equal(projected.config.adv.feature_overlap_tolerance_bp, 2);
    assert.deepEqual(Object.values(projected.config.featurePlacementOverrides), [override]);
    const old = structuredClone(result.renderRequest);
    old.schema = 6;
    assert.throws(() => projectCanonicalSessionRequest({ ...result, renderRequest: old }), /schema 7/);
  });
}

test('Main, resolved side, bulk Auto and history share one draft owner', async () => {
  const features = ['one', 'two'].map((id) => ({ record_key: 'card', biological_feature_id: id }));
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
    const feature = { record_key: 'record-1', biological_feature_id: 'feature' };
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

test('historical session 40 schema 6 promotes without Generate and preserves cardinality', async () => {
  const bytes = await readFile('tests/fixtures/sessions/test_linear_cli_sidecar_reuses0.v40-schema6.json.gz');
  const session = JSON.parse(gunzipSync(bytes));
  assert.equal(session.version, 40);
  assert.equal(session.renderRequest.schema, 6);
  const before = structuredClone(session.renderRequest);
  const current = promoteCanonicalRenderRequestToCurrent(session.renderRequest);
  assert.equal(current.schema, 7);
  assert.deepEqual(current.records.map((record) => record.cardinality), before.records.map((record) => record.cardinality));
  assert.ok(current.records.every((record) => record.display.isCircular === null && record.display.startCoordinate === null));
  assert.deepEqual(current.diagramOptions.featurePlacements, []);
  assert.deepEqual(session.renderRequest, before);
});
