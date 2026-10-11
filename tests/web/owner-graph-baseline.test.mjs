import assert from 'node:assert/strict';
import test from 'node:test';

import { evaluateWebOwnerGraphAt } from '../../tools/report-web-owner-graph.mjs';

// R13 baseline (gbdraw/web/CLAUDE.md "Owner layering and ports"). The
// detectors in tools/web-owner-graph-detectors.mjs observe the coupling
// channels the import graph does not: owner objects injected at composition,
// closures bound to owners created later, functions assigned into state,
// factories that take whole owner objects, heavy derived trigger sites, and
// the call shapes of each projection domain outside its owner.
//
// The three literals are registered design-rule allowlists
// (tools/web-design-rule-guards.json, R13): a runtime change removes the
// subjects it fixes from the baseline in the same pull request, and an
// addition is an authority-only change. Each entry is `detector|subject`
// with its count (1 for a subject that is present; the number of trigger
// sites for a trigger-site module). Lower an entry or remove it when the
// subject is gone; the assertion names what to change.
//
// Print the current subjects with `node tools/report-web-owner-graph.mjs --at worktree`.
const OWNER_GRAPH_BASELINE = {
  // owner-graph.injection-edge.v1
  'owner-graph.injection-edge.v1|app/app-setup.js|featureActions<-featureSelection': 1,
  'owner-graph.injection-edge.v1|app/app-setup.js|featureActions<-previewTransformInteraction': 1,
  'owner-graph.injection-edge.v1|app/app-setup.js|featureActions<-rulePreparation': 1,
  'owner-graph.injection-edge.v1|app/app-setup.js|featureRecordRotation<-featureRecordRotationAction': 1,
  'owner-graph.injection-edge.v1|app/app-setup.js|history<-historyFileStore': 1,
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
  'heavy-derived.trigger-site.v2|rulePreparation|app/run-analysis.js': 1
};

// Distinct call shapes per projection domain outside its owner (R3: one
// projection is one call shape). May only decrease.
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
  'owner-graph.whole-object-port.v1'
];

const { results } = evaluateWebOwnerGraphAt('worktree');

test('every owner-graph subject in the working tree is in the baseline, and every baseline entry is still observed (R13)', () => {
  const observed = new Map();
  SUBJECT_DETECTORS.forEach((id) => {
    results[id].subjects.forEach((subject) => observed.set(`${id}|${subject}`, 1));
  });
  Object.entries(results['heavy-derived.trigger-site.v2'].countsBySubject).forEach(([subject, count]) => {
    observed.set(`heavy-derived.trigger-site.v2|${subject}`, count);
  });
  const added = [...observed].filter(([key, count]) => !(key in OWNER_GRAPH_BASELINE) || count > OWNER_GRAPH_BASELINE[key]);
  const fixed = Object.entries(OWNER_GRAPH_BASELINE).filter(([key, count]) => !observed.has(key) || observed.get(key) < count);
  assert.deepEqual(added.map(([key, count]) => `${key}: ${count}`), [],
    'new owner-graph subjects (R13): remove the coupling, or register the entry in an authority-only pull request');
  assert.deepEqual(fixed.map(([key, count]) => `${key}: ${count} -> ${observed.get(key) ?? 0}`), [],
    'fixed owner-graph subjects: lower or remove these OWNER_GRAPH_BASELINE entries in this pull request');
});

test('the distinct projection call shapes per domain do not grow (R3, R13)', () => {
  const shapes = results['projection.call-shape.v1'].shapesByDomain;
  const grown = Object.entries(shapes).filter(([domain, count]) => !(domain in PROJECTION_SHAPE_BASELINE) || count > PROJECTION_SHAPE_BASELINE[domain]);
  const shrunk = Object.entries(PROJECTION_SHAPE_BASELINE).filter(([domain, count]) => (shapes[domain] ?? 0) < count);
  assert.deepEqual(grown.map(([domain, count]) => `${domain}: ${count}`), [],
    `new projection call shapes (R3): ${JSON.stringify(results['projection.call-shape.v1'].observedShapes.map(({ domain, path, line, shape }) => `${domain} ${path}:${line} ${shape}`), null, 1)}`);
  assert.deepEqual(shrunk.map(([domain, count]) => `${domain}: ${count} -> ${shapes[domain] ?? 0}`), [],
    'fewer projection call shapes: lower these PROJECTION_SHAPE_BASELINE entries in this pull request');
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
