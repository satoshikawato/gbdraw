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
// Both literals are registered design-rule allowlists
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
  'owner-graph.injection-edge.v1|app/app-setup.js|featureActions<-history': 1,
  'owner-graph.injection-edge.v1|app/app-setup.js|featureActions<-legendActions': 1,
  'owner-graph.injection-edge.v1|app/app-setup.js|featureActions<-previewRuntime': 1,
  'owner-graph.injection-edge.v1|app/app-setup.js|featureActions<-previewTransformInteraction': 1,
  'owner-graph.injection-edge.v1|app/app-setup.js|featureActions<-rulePreparation': 1,
  'owner-graph.injection-edge.v1|app/app-setup.js|featureActions<-svgActions': 1,
  'owner-graph.injection-edge.v1|app/app-setup.js|featureRecordRotation<-featureRecordRotationAction': 1,
  'owner-graph.injection-edge.v1|app/app-setup.js|history<-historyFileStore': 1,
  'owner-graph.injection-edge.v1|app/app-setup.js|historySnapshots<-historyFileStore': 1,
  'owner-graph.injection-edge.v1|app/app-setup.js|legendActions<-previewRuntime': 1,
  'owner-graph.injection-edge.v1|app/app-setup.js|legendActions<-rulePreparation': 1,
  'owner-graph.injection-edge.v1|app/app-setup.js|legendLayout<-previewRuntime': 1,
  'owner-graph.injection-edge.v1|app/app-setup.js|previewFeatureSearch<-previewRuntime': 1,
  'owner-graph.injection-edge.v1|app/app-setup.js|recordDisplayControls<-history': 1,
  'owner-graph.injection-edge.v1|app/app-setup.js|recordDisplayControls<-linearRecordSelector': 1,
  'owner-graph.injection-edge.v1|app/app-setup.js|svgActions<-legendActions': 1,
  'owner-graph.injection-edge.v1|app/app-setup.js|svgActions<-previewRuntime': 1,
  'owner-graph.injection-edge.v1|app/app-setup.js|svgActions<-rulePreparation': 1,
  'owner-graph.injection-edge.v1|app/feature-editor.js|colorActions<-featureSvgActions': 1,
  'owner-graph.injection-edge.v1|app/feature-editor.js|colorActions<-legendActions': 1,
  'owner-graph.injection-edge.v1|app/feature-editor.js|colorActions<-previewRuntime': 1,
  'owner-graph.injection-edge.v1|app/feature-editor.js|colorActions<-ruleActions': 1,
  'owner-graph.injection-edge.v1|app/feature-editor.js|colorActions<-rulePreparation': 1,
  'owner-graph.injection-edge.v1|app/feature-editor.js|colorActions<-svgActions': 1,
  'owner-graph.injection-edge.v1|app/feature-editor.js|featureSvgActions<-featureSelection': 1,
  'owner-graph.injection-edge.v1|app/feature-editor.js|featureSvgActions<-previewRuntime': 1,
  'owner-graph.injection-edge.v1|app/feature-editor.js|featureSvgActions<-previewTransformInteraction': 1,
  'owner-graph.injection-edge.v1|app/feature-editor.js|featureSvgActions<-rulePreparation': 1,
  'owner-graph.injection-edge.v1|app/feature-editor.js|labelActions<-previewRuntime': 1,
  'owner-graph.injection-edge.v1|app/feature-editor.js|labelActions<-rulePreparation': 1,
  'owner-graph.injection-edge.v1|app/feature-editor.js|ruleActions<-history': 1,
  'owner-graph.injection-edge.v1|app/feature-editor.js|ruleActions<-legendActions': 1,
  'owner-graph.injection-edge.v1|app/feature-editor.js|ruleActions<-rulePreparation': 1,
  'owner-graph.injection-edge.v1|app/feature-editor.js|ruleActions<-svgActions': 1,
  'owner-graph.injection-edge.v1|app/feature-editor.js|visibilityActions<-featureSvgActions': 1,
  'owner-graph.injection-edge.v1|app/feature-editor.js|visibilityActions<-previewRuntime': 1,
  'owner-graph.injection-edge.v1|app/feature-editor.js|visibilityActions<-rulePreparation': 1,
  'owner-graph.injection-edge.v1|app/legend-layout.js|canvasActions<-previewRuntime': 1,
  'owner-graph.injection-edge.v1|app/legend-layout.js|diagramActions<-previewRuntime': 1,
  'owner-graph.injection-edge.v1|app/legend-layout.js|repositionActions<-previewRuntime': 1,
  'owner-graph.injection-edge.v1|app/legend.js|dragActions<-previewRuntime': 1,
  'owner-graph.injection-edge.v1|app/legend.js|entryActions<-layoutActions': 1,
  'owner-graph.injection-edge.v1|app/legend.js|entryActions<-previewRuntime': 1,
  'owner-graph.injection-edge.v1|app/legend.js|sortActions<-previewRuntime': 1,
  'owner-graph.injection-edge.v1|app/legend.js|strokeActions<-previewRuntime': 1,
  // owner-graph.forward-closure.v1
  'owner-graph.forward-closure.v1|app/app-setup.js|legendActions.commitSpecificRules->featureActions.commitSpecificRules': 1,
  'owner-graph.forward-closure.v1|app/app-setup.js|previewFeatureSearch.resolveOrthogroups->orthogroupActions.resolveOrthogroupDescription': 1,
  'owner-graph.forward-closure.v1|app/app-setup.js|previewFeatureSearch.resolveOrthogroups->orthogroupActions.resolveOrthogroupName': 1,
  'owner-graph.forward-closure.v1|app/app-setup.js|rightDrawerActions.onClose->featureActions.suspendSpecificRulePatternDrafts': 1,
  // owner-graph.whole-object-port.v1
  'owner-graph.whole-object-port.v1|app/feature-editor/color-actions.js|createFeatureColorActions|featureSvgActions': 1,
  'owner-graph.whole-object-port.v1|app/feature-editor/color-actions.js|createFeatureColorActions|legendActions': 1,
  'owner-graph.whole-object-port.v1|app/feature-editor/color-actions.js|createFeatureColorActions|previewRuntime': 1,
  'owner-graph.whole-object-port.v1|app/feature-editor/color-actions.js|createFeatureColorActions|ruleActions': 1,
  'owner-graph.whole-object-port.v1|app/feature-editor/color-actions.js|createFeatureColorActions|rulePreparation': 1,
  'owner-graph.whole-object-port.v1|app/feature-editor/color-actions.js|createFeatureColorActions|svgActions': 1,
  'owner-graph.whole-object-port.v1|app/feature-editor/label-actions.js|createFeatureLabelActions|previewRuntime': 1,
  'owner-graph.whole-object-port.v1|app/feature-editor/label-actions.js|createFeatureLabelActions|rulePreparation': 1,
  'owner-graph.whole-object-port.v1|app/feature-editor/placement-actions.js|createFeaturePlacementActions|history': 1,
  'owner-graph.whole-object-port.v1|app/feature-editor/rule-actions.js|createFeatureRuleActions|history': 1,
  'owner-graph.whole-object-port.v1|app/feature-editor/rule-actions.js|createFeatureRuleActions|legendActions': 1,
  'owner-graph.whole-object-port.v1|app/feature-editor/rule-actions.js|createFeatureRuleActions|rulePreparation': 1,
  'owner-graph.whole-object-port.v1|app/feature-editor/rule-actions.js|createFeatureRuleActions|svgActions': 1,
  'owner-graph.whole-object-port.v1|app/feature-editor/svg-actions.js|createFeatureSvgActions|featureSelection': 1,
  'owner-graph.whole-object-port.v1|app/feature-editor/svg-actions.js|createFeatureSvgActions|previewRuntime': 1,
  'owner-graph.whole-object-port.v1|app/feature-editor/svg-actions.js|createFeatureSvgActions|previewTransformInteraction': 1,
  'owner-graph.whole-object-port.v1|app/feature-editor/svg-actions.js|createFeatureSvgActions|rulePreparation': 1,
  'owner-graph.whole-object-port.v1|app/feature-editor/visibility-actions.js|createFeatureVisibilityActions|featureSvgActions': 1,
  'owner-graph.whole-object-port.v1|app/feature-editor/visibility-actions.js|createFeatureVisibilityActions|previewRuntime': 1,
  'owner-graph.whole-object-port.v1|app/feature-editor/visibility-actions.js|createFeatureVisibilityActions|rulePreparation': 1,
  'owner-graph.whole-object-port.v1|app/feature-search/preview-actions.js|createPreviewFeatureSearch|previewRuntime': 1,
  'owner-graph.whole-object-port.v1|app/legend-layout/canvas-actions.js|createLegendCanvasActions|previewRuntime': 1,
  'owner-graph.whole-object-port.v1|app/legend-layout/diagram-drag.js|createDiagramDragActions|previewRuntime': 1,
  'owner-graph.whole-object-port.v1|app/legend-layout/reposition-actions.js|createLegendRepositionActions|previewRuntime': 1,
  'owner-graph.whole-object-port.v1|app/legend/drag-actions.js|createLegendDragActions|previewRuntime': 1,
  'owner-graph.whole-object-port.v1|app/legend/entry-actions.js|createLegendEntryActions|layoutActions': 1,
  'owner-graph.whole-object-port.v1|app/legend/entry-actions.js|createLegendEntryActions|previewRuntime': 1,
  'owner-graph.whole-object-port.v1|app/legend/sort-actions.js|createLegendSortActions|previewRuntime': 1,
  'owner-graph.whole-object-port.v1|app/legend/stroke-actions.js|createLegendStrokeActions|previewRuntime': 1,
  'owner-graph.whole-object-port.v1|app/record-display-options.js|createRecordDisplayControls|history': 1,
  'owner-graph.whole-object-port.v1|app/record-display-options.js|createRecordDisplayControls|linearRecordSelector': 1,
  'owner-graph.whole-object-port.v1|app/run-analysis.js|createRunAnalysis|previewRuntime': 1,
  'owner-graph.whole-object-port.v1|app/run-analysis.js|createRunAnalysis|rulePreparation': 1,
  'owner-graph.whole-object-port.v1|app/svg-styles.js|createSvgStyles|legendActions': 1,
  'owner-graph.whole-object-port.v1|app/svg-styles.js|createSvgStyles|previewRuntime': 1,
  'owner-graph.whole-object-port.v1|app/svg-styles.js|createSvgStyles|rulePreparation': 1,
  'owner-graph.whole-object-port.v1|services/history-snapshot.js|createHistorySnapshotService|fileStore': 1,
  'owner-graph.whole-object-port.v1|services/history.js|createHistoryManager|fileStore': 1,
  // heavy-derived.trigger-site.v1 (sites per module)
  'heavy-derived.trigger-site.v1|rulePreparation|app/app-setup.js': 2,
  'heavy-derived.trigger-site.v1|rulePreparation|app/feature-editor/color-actions.js': 2,
  'heavy-derived.trigger-site.v1|rulePreparation|app/feature-editor/label-actions.js': 2,
  'heavy-derived.trigger-site.v1|rulePreparation|app/feature-editor/rule-actions.js': 1,
  'heavy-derived.trigger-site.v1|rulePreparation|app/feature-editor/svg-actions.js': 1,
  'heavy-derived.trigger-site.v1|rulePreparation|app/feature-editor/visibility-actions.js': 1,
  'heavy-derived.trigger-site.v1|rulePreparation|app/legend.js': 1,
  'heavy-derived.trigger-site.v1|rulePreparation|app/run-analysis.js': 2,
  'heavy-derived.trigger-site.v1|rulePreparation|app/svg-styles.js': 1
};

// Distinct call shapes per projection domain outside its owner (R3: one
// projection is one call shape). May only decrease.
const PROJECTION_SHAPE_BASELINE = {
  'feature-visibility': 1,
  'feature-visibility-labels': 1,
  'label-intent': 1,
  'legend-order': 2,
  'palette-rules': 4,
  'strokes': 2,
  'composition': 1
};

const SUBJECT_DETECTORS = [
  'owner-graph.injection-edge.v1',
  'owner-graph.forward-closure.v1',
  'owner-graph.state-backdoor.v1',
  'owner-graph.whole-object-port.v1'
];

const { results } = evaluateWebOwnerGraphAt('worktree');

test('every owner-graph subject in the working tree is in the baseline, and every baseline entry is still observed (R13)', () => {
  const observed = new Map();
  SUBJECT_DETECTORS.forEach((id) => {
    results[id].subjects.forEach((subject) => observed.set(`${id}|${subject}`, 1));
  });
  Object.entries(results['heavy-derived.trigger-site.v1'].countsBySubject).forEach(([subject, count]) => {
    observed.set(`heavy-derived.trigger-site.v1|${subject}`, count);
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
