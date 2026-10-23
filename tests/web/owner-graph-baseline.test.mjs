import assert from 'node:assert/strict';
import test from 'node:test';

import { evaluateWebOwnerGraphAt } from '../../tools/report-web-owner-graph.mjs';
import { compareCountBaseline, reportLowerable } from './helpers/count-baseline.mjs';

// R13 baseline (gbdraw/web/CLAUDE.md "Owner layering and ports"). The
// detectors in tools/web-owner-graph-detectors.mjs observe the coupling
// channels the import graph does not: owner objects injected at composition,
// closures bound to owners created later, functions assigned into state,
// factories that take whole owner objects, heavy derived trigger sites,
// identities read from display values, and the call shapes of each projection
// domain outside its owner.
//
// The two literals are registered design-rule allowlists
// (tools/web-design-rule-guards.json, R13, kind `count-map`): an observed
// count may not exceed its entry, and adding or raising an entry is an
// authority-only change. Each entry is `detector|subject` with its count (1
// for a subject that is present; the number of trigger sites for a
// trigger-site module). Lowering or removing an entry whose subject is gone is
// optional; the test names those entries.
//
// Print the current subjects with `node tools/report-web-owner-graph.mjs --at worktree`.
const OWNER_GRAPH_BASELINE = {
  // owner-graph.injection-edge.v1
  'owner-graph.injection-edge.v1|app/app-setup.js|featureActions<-featureSelection': 1,
  'owner-graph.injection-edge.v1|app/app-setup.js|featureActions<-previewTransformInteraction': 1,
  'owner-graph.injection-edge.v1|app/app-setup.js|featureActions<-rulePreparation': 1,
  'owner-graph.injection-edge.v1|app/app-setup.js|featureRecordRotation<-featureRecordRotationAction': 1,
  'owner-graph.injection-edge.v1|app/app-setup.js|historySnapshots<-historyFileStore': 1,
  'owner-graph.injection-edge.v1|app/app-setup.js|legendLayout<-previewRuntime': 1,
  'owner-graph.injection-edge.v1|app/feature-editor.js|colorActions<-ruleActions': 1,
  'owner-graph.injection-edge.v1|app/feature-editor.js|featureSvgActions<-featureSelection': 1,
  'owner-graph.injection-edge.v1|app/feature-editor.js|featureSvgActions<-previewTransformInteraction': 1,
  'owner-graph.injection-edge.v1|app/feature-editor.js|labelActions<-rulePreparation': 1,
  'owner-graph.injection-edge.v1|app/feature-editor.js|ruleActions<-rulePreparation': 1,
  'owner-graph.injection-edge.v1|app/feature-editor.js|visibilityActions<-rulePreparation': 1,
  // owner-graph.whole-object-port.v1
  'owner-graph.whole-object-port.v1|app/feature-editor/color-actions.js|createFeatureColorActions|ruleActions': 1,
  'owner-graph.whole-object-port.v1|app/feature-editor/label-actions.js|createFeatureLabelActions|rulePreparation': 1,
  'owner-graph.whole-object-port.v1|app/feature-editor/rule-actions.js|createFeatureRuleActions|rulePreparation': 1,
  'owner-graph.whole-object-port.v1|app/feature-editor/svg-actions.js|createFeatureSvgActions|featureSelection': 1,
  'owner-graph.whole-object-port.v1|app/feature-editor/svg-actions.js|createFeatureSvgActions|previewTransformInteraction': 1,
  'owner-graph.whole-object-port.v1|app/feature-editor/visibility-actions.js|createFeatureVisibilityActions|rulePreparation': 1,
  'owner-graph.whole-object-port.v1|app/run-analysis.js|createRunAnalysis|previewRuntime': 1,
  'owner-graph.whole-object-port.v1|app/run-analysis.js|createRunAnalysis|rulePreparation': 1,
  'owner-graph.whole-object-port.v1|services/history-snapshot.js|createHistorySnapshotService|fileStore': 1,
  'owner-graph.whole-object-port.v1|services/history.js|createHistoryManager|fileStore': 1,
  // heavy-derived.trigger-site.v2 (sites per module; includes calls of a
  // producer port such as rulePreparation.run handed to an owner)
  'heavy-derived.trigger-site.v2|rulePreparation|app/app-setup.js': 2,
  'heavy-derived.trigger-site.v2|rulePreparation|app/feature-editor/rule-actions.js': 1,
  'heavy-derived.trigger-site.v2|rulePreparation|app/feature-editor/visibility-actions.js': 1,
  'heavy-derived.trigger-site.v2|rulePreparation|app/run-analysis.js': 1,
  // owner-graph.identity-from-display.v1 (a Legend row or a feature joined on
  // a fill or stroke, `paint-join`, or on a value a Result shows, `shown-join`)
  'owner-graph.identity-from-display.v1|app/candidate-render.js|matchingRuleDerivedFill|paint-join': 1,
  'owner-graph.identity-from-display.v1|app/feature-editor/color-actions.js|findLegendRowNamed|paint-join': 1,
  'owner-graph.identity-from-display.v1|app/feature-editor/rule-actions.js|commitOnce|paint-join': 1,
  'owner-graph.identity-from-display.v1|app/legend/entry-actions.js|findLegendRowNamed|paint-join': 1,
  'owner-graph.identity-from-display.v1|app/legend/entry-actions.js|listLegendRows|shown-join': 1,
  'owner-graph.identity-from-display.v1|app/legend/entry-actions.js|namedBy|shown-join': 1,
  'owner-graph.identity-from-display.v1|app/legend/entry-actions.js|prepareFileLegendEntries|paint-join': 1,
  'owner-graph.identity-from-display.v1|services/legend-svg.js|legendRowFeatureIds|paint-join': 1,
  'owner-graph.identity-from-display.v1|services/specific-color-rules.js|captionFills|paint-join': 1,
  'owner-graph.identity-from-display.v1|services/specific-color-rules.js|createRuleLegendCaptions|paint-join': 1,
  'owner-graph.identity-from-display.v1|services/svg-result-ingestion.js|applyLegendOperations|shown-join': 1,
  'owner-graph.identity-from-display.v1|services/svg-result-ingestion.js|legends|shown-join': 1,
  'owner-graph.identity-from-display.v1|services/svg-result-ingestion.js|nameSavedRuleRows|shown-join': 1,
  'owner-graph.identity-from-display.v1|services/svg-result-ingestion.js|recordSavedLegendKeys|shown-join': 1,
  'owner-graph.identity-from-display.v1|services/svg-result-ingestion.js|restoreLegendStructure|shown-join': 1
};

// Distinct call shapes per projection domain outside its owner (R3: one
// projection is one call shape). May not grow.
const PROJECTION_SHAPE_BASELINE = {
  'feature-visibility': 1,
  'feature-visibility-labels': 1,
  'label-intent': 1,
  'legend-order': 1,
  'palette-rules': 1,
  'strokes': 0,
  'composition': 1
};

const SUBJECT_DETECTORS = [
  'owner-graph.injection-edge.v1',
  'owner-graph.forward-closure.v2',
  'owner-graph.state-backdoor.v1',
  'owner-graph.whole-object-port.v1',
  'owner-graph.identity-from-display.v1'
];

const { results } = evaluateWebOwnerGraphAt('worktree');

test('no owner-graph subject in the working tree exceeds its baseline entry (R13)', (t) => {
  const observed = {};
  SUBJECT_DETECTORS.forEach((id) => {
    results[id].subjects.forEach((subject) => { observed[`${id}|${subject}`] = 1; });
  });
  Object.entries(results['heavy-derived.trigger-site.v2'].countsBySubject).forEach(([subject, count]) => {
    observed[`heavy-derived.trigger-site.v2|${subject}`] = count;
  });
  const { above, lowerable } = compareCountBaseline(observed, OWNER_GRAPH_BASELINE);
  reportLowerable(t, 'OWNER_GRAPH_BASELINE', lowerable);
  assert.deepEqual(above.map(({ key, observed: count }) => `${key}: ${count}`), [],
    'new owner-graph subjects (R13): remove the coupling, or register the entry in an authority-only pull request');
});

test('the distinct projection call shapes per domain do not grow (R3, R13)', (t) => {
  const { above, lowerable } = compareCountBaseline(results['projection.call-shape.v1'].shapesByDomain, PROJECTION_SHAPE_BASELINE);
  reportLowerable(t, 'PROJECTION_SHAPE_BASELINE', lowerable);
  assert.deepEqual(above.map(({ key, observed: count }) => `${key}: ${count}`), [],
    `new projection call shapes (R3): ${JSON.stringify(results['projection.call-shape.v1'].observedShapes.map(({ domain, path, line, shape }) => `${domain} ${path}:${line} ${shape}`), null, 1)}`);
});

// R13 layers (`webLayerOf` in tools/web-owner-graph-detectors.mjs): a module
// imports its own layer and lower ones. The layering slices removed every
// upward import, so none is allowed.
test('no module imports from a higher R13 layer', () => {
  const upward = Object.entries(results['layer.import-direction.v1'].countsBySubject)
    .map(([subject, count]) => `${subject}: ${count}`);
  assert.deepEqual(upward, [],
    'upward imports (R13): import from the same layer or a lower one (move the code down, or take a port)');
});
