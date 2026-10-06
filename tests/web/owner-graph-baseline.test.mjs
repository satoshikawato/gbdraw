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
  'heavy-derived.trigger-site.v2|rulePreparation|app/feature-editor/color-actions.js': 2,
  'heavy-derived.trigger-site.v2|rulePreparation|app/feature-editor/label-actions.js': 2,
  'heavy-derived.trigger-site.v2|rulePreparation|app/feature-editor/rule-actions.js': 1,
  'heavy-derived.trigger-site.v2|rulePreparation|app/feature-editor/svg-actions.js': 1,
  'heavy-derived.trigger-site.v2|rulePreparation|app/feature-editor/visibility-actions.js': 1,
  'heavy-derived.trigger-site.v2|rulePreparation|app/run-analysis.js': 1
};

// Distinct call shapes per projection domain outside its owner (R3: one
// projection is one call shape). May only decrease.
const PROJECTION_SHAPE_BASELINE = {
  'feature-visibility': 1,
  'feature-visibility-labels': 1,
  'label-intent': 1,
  'legend-order': 2,
  'palette-rules': 2,
  'strokes': 2,
  'composition': 1
};

// Upward imports under the R13 layers (`webLayerOf` in
// tools/web-owner-graph-detectors.mjs): a module imports its own layer and
// lower ones. Each entry is `layer.import-direction.v1|importer->target` with
// the number of distinct names the importer takes from the target. A slice
// that moves the imported code down removes its entries in the same pull
// request; the list ends empty. May only shrink.
const LAYER_IMPORT_BASELINE = {
  'layer.import-direction.v1|mode-profiles.js->services/error-normalization.js': 1,
  'layer.import-direction.v1|services/config.js->app/annotations/state.js': 1,
  'layer.import-direction.v1|services/config.js->app/circular-track-slots.js': 7,
  'layer.import-direction.v1|services/config.js->app/color-utils.js': 1,
  'layer.import-direction.v1|services/config.js->app/current-option-values.js': 9,
  'layer.import-direction.v1|services/config.js->app/definition-line-style-state.js': 1,
  'layer.import-direction.v1|services/config.js->app/depth-track-state.js': 8,
  'layer.import-direction.v1|services/config.js->app/feature-visibility.js': 3,
  'layer.import-direction.v1|services/config.js->app/layout-preferences.js': 4,
  'layer.import-direction.v1|services/config.js->app/legend-layout/composition-actions.js': 3,
  'layer.import-direction.v1|services/config.js->app/legend/stroke-actions.js': 1,
  'layer.import-direction.v1|services/config.js->app/linear-comparisons.js': 6,
  'layer.import-direction.v1|services/config.js->app/linear-label-visibility.js': 2,
  'layer.import-direction.v1|services/config.js->app/linear-track-slots.js': 7,
  'layer.import-direction.v1|services/config.js->app/linear-typography.js': 1,
  'layer.import-direction.v1|services/config.js->app/losat-cache.js': 13,
  'layer.import-direction.v1|services/config.js->app/losat-normalization.js': 3,
  'layer.import-direction.v1|services/config.js->app/match-sequences.js': 3,
  'layer.import-direction.v1|services/config.js->app/plot-title-position.js': 1,
  'layer.import-direction.v1|services/config.js->app/record-display-options.js': 1,
  'layer.import-direction.v1|services/config.js->app/right-drawer.js': 3,
  'layer.import-direction.v1|services/config.js->app/run-info.js': 1,
  'layer.import-direction.v1|services/config.js->app/session-feature-metadata.js': 4,
  'layer.import-direction.v1|services/config.js->app/specific-color-rules.js': 1,
  'layer.import-direction.v1|services/export.js->app/feature-search/preview-svg.js': 1,
  'layer.import-direction.v1|services/feature-placement.js->app/feature-utils.js': 1,
  'layer.import-direction.v1|services/gallery-session-migration.js->app/circular-track-slots.js': 3,
  'layer.import-direction.v1|services/gallery-session-migration.js->app/current-option-values.js': 4,
  'layer.import-direction.v1|services/gallery-session-migration.js->app/layout-preferences.js': 1,
  'layer.import-direction.v1|services/gallery-session-migration.js->app/linear-comparisons.js': 3,
  'layer.import-direction.v1|services/gallery-session-migration.js->app/linear-label-visibility.js': 1,
  'layer.import-direction.v1|services/gallery-session-publication.js->app/layout-preferences.js': 1,
  'layer.import-direction.v1|services/gallery-session-publication.js->app/linear-label-visibility.js': 1,
  'layer.import-direction.v1|services/gallery-session-publication.js->app/record-display-options.js': 1,
  'layer.import-direction.v1|services/history-snapshot.js->app/feature-visibility.js': 1,
  'layer.import-direction.v1|services/history-snapshot.js->app/layout-preferences.js': 1,
  'layer.import-direction.v1|services/main-session-comparison-frame.js->app/record-discovery.js': 2,
  'layer.import-direction.v1|services/orthogroup-feature-metadata.js->app/losat-normalization.js': 1,
  'layer.import-direction.v1|services/reset.js->app/layout-preferences.js': 1,
  'layer.import-direction.v1|services/session-active-config-contract.js->app/current-option-values.js': 20,
  'layer.import-direction.v1|services/session-active-config-contract.js->app/definition-line-style-state.js': 1,
  'layer.import-direction.v1|services/session-active-config-contract.js->app/linear-label-visibility.js': 1,
  'layer.import-direction.v1|services/session-active-config-contract.js->app/linear-track-slots.js': 4,
  'layer.import-direction.v1|services/session-active-config-contract.js->app/record-display-options.js': 2,
  'layer.import-direction.v1|services/session-authority.js->app/linear-label-visibility.js': 1,
  'layer.import-direction.v1|services/session-authority.js->app/record-display-options.js': 1,
  'layer.import-direction.v1|services/session-request.js->app/annotations/state.js': 2,
  'layer.import-direction.v1|services/session-request.js->app/circular-track-slots.js': 10,
  'layer.import-direction.v1|services/session-request.js->app/circular-track-slots/measure-editor.js': 1,
  'layer.import-direction.v1|services/session-request.js->app/color-utils.js': 2,
  'layer.import-direction.v1|services/session-request.js->app/conservation-series.js': 2,
  'layer.import-direction.v1|services/session-request.js->app/current-option-values.js': 23,
  'layer.import-direction.v1|services/session-request.js->app/depth-track-state.js': 3,
  'layer.import-direction.v1|services/session-request.js->app/feature-editor/label-override-table.js': 3,
  'layer.import-direction.v1|services/session-request.js->app/feature-visibility.js': 2,
  'layer.import-direction.v1|services/session-request.js->app/file-imports.js': 7,
  'layer.import-direction.v1|services/session-request.js->app/genbank-header.js': 1,
  'layer.import-direction.v1|services/session-request.js->app/layout-preferences.js': 2,
  'layer.import-direction.v1|services/session-request.js->app/linear-comparisons.js': 1,
  'layer.import-direction.v1|services/session-request.js->app/linear-label-visibility.js': 1,
  'layer.import-direction.v1|services/session-request.js->app/linear-record-layout.js': 2,
  'layer.import-direction.v1|services/session-request.js->app/linear-sources.js': 2,
  'layer.import-direction.v1|services/session-request.js->app/linear-track-slots.js': 5,
  'layer.import-direction.v1|services/session-request.js->app/losat-normalization.js': 3,
  'layer.import-direction.v1|services/session-request.js->app/record-display-options.js': 3,
  'layer.import-direction.v1|services/session-request.js->app/record-options.js': 3,
  'layer.import-direction.v1|services/session-request.js->app/track-slot-validation.js': 4,
  'layer.import-direction.v1|services/standalone-interactivity.js->app/feature-utils.js': 1,
  'layer.import-direction.v1|services/standalone-interactivity.js->app/record-source-coordinates.js': 2,
  'layer.import-direction.v1|services/svg-result-ingestion.js->app/feature-dom.js': 3,
  'layer.import-direction.v1|services/svg-result-ingestion.js->app/legend/utils.js': 4,
  'layer.import-direction.v1|state.js->app/color-utils.js': 3,
  'layer.import-direction.v1|state.js->app/feature-selector.js': 1,
  'layer.import-direction.v1|state.js->app/feature-visibility.js': 3,
  'layer.import-direction.v1|state.js->app/layout-preferences.js': 3,
  'layer.import-direction.v1|state.js->app/linear-comparisons.js': 2,
  'layer.import-direction.v1|state.js->app/match-sequences.js': 1,
  'layer.import-direction.v1|state.js->app/plot-title-position.js': 1,
  'layer.import-direction.v1|utils/optional-positive-number.js->services/error-normalization.js': 1
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

test('every upward import under the R13 layers is in the baseline, and every baseline entry is still observed', () => {
  const observed = new Map(Object.entries(results['layer.import-direction.v1'].countsBySubject)
    .map(([subject, count]) => [`layer.import-direction.v1|${subject}`, count]));
  const added = [...observed].filter(([key, count]) => !(key in LAYER_IMPORT_BASELINE) || count > LAYER_IMPORT_BASELINE[key]);
  const fixed = Object.entries(LAYER_IMPORT_BASELINE).filter(([key, count]) => !observed.has(key) || observed.get(key) < count);
  assert.deepEqual(added.map(([key, count]) => `${key}: ${count}`), [],
    'new upward imports (R13): import from the same layer or a lower one (move the code down, or take a port), or register the entry in an authority-only pull request');
  assert.deepEqual(fixed.map(([key, count]) => `${key}: ${count} -> ${observed.get(key) ?? 0}`), [],
    'fixed upward imports: lower or remove these LAYER_IMPORT_BASELINE entries in this pull request');
});
