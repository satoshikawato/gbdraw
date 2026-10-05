import assert from 'node:assert/strict';
import test from 'node:test';

import { evaluateWebOwnerGraphAt, revisionExists } from '../../tools/report-web-owner-graph.mjs';

// Untouched-base characterization and historical replay of the owner-graph
// detectors (implementation plan Phase B1). Both read revisions from Git, so
// they run where the history is present (a full checkout or the trusted-base
// policy job) and skip in a shallow checkout.
const BASE = '8f450194';

// First-parent merges that introduced a coupling the audit attributes to the
// named detector; each must be observable as a NEW subject against its first
// parent.
const REPLAY = [
  { pr: '#617', merge: 'f5f86634', detector: 'owner-graph.forward-closure.v1', subject: 'app/app-setup.js|legendActions.commitSpecificRules->featureActions.commitSpecificRules' },
  { pr: '#622', merge: '27939fae', detector: 'owner-graph.forward-closure.v1', subject: 'app/app-setup.js|historySnapshots.buildLegendEntryOwners->legendActions.captureLegendEntryOwners' },
  { pr: '#641', merge: 'b1744142', detector: 'owner-graph.state-backdoor.v1', subject: 'app/app-setup.js|sessionPreparationBusyReason' },
  { pr: '#692', merge: '942fa954', detector: 'owner-graph.injection-edge.v1', subject: 'app/feature-editor.js|visibilityActions<-labelActions' },
  { pr: '#737', merge: 'cafed980', detector: 'owner-graph.forward-closure.v1', subject: 'app/app-setup.js|historySnapshots.buildCompositionIntent->legendLayout.captureCompositionIntent' },
  { pr: '#764', merge: '0bf9034a', detector: 'owner-graph.forward-closure.v1', subject: 'app/feature-editor.js|labelActions.setFeatureVisibility->visibilityActions.setFeatureVisibility' },
  { pr: '#805', merge: '1808d20d', detector: 'owner-graph.whole-object-port.v1', subject: 'app/circular-track-slots.js|createCircularTrackSlotEditor|trackLayoutActions' },
  { pr: '#807', merge: '4cf78f5d', detector: 'owner-graph.state-backdoor.v1', subject: 'app/app-setup.js|committedDiagramOptions' }
];

const historyAvailable = revisionExists(BASE) && REPLAY.every(({ merge }) => revisionExists(`${merge}^1`));

test(`untouched-base characterization at ${BASE}`, { skip: historyAvailable ? false : 'revision history not available in this checkout' }, () => {
  const { results, summary } = evaluateWebOwnerGraphAt(BASE);
  assert.deepEqual(summary, {
    injectionEdges: 57,
    forwardClosures: 7,
    stateBackdoors: 3,
    wholeObjectPorts: 46,
    projectionShapes: {
      'feature-visibility': 2,
      'feature-visibility-labels': 5,
      'label-intent': 2,
      'legend-order': 3,
      'palette-rules': 4,
      strokes: 2,
      composition: 1
    },
    triggerSites: 15,
    triggerModules: 9
  });
  assert.deepEqual(results['owner-graph.forward-closure.v1'].subjects, [
    'app/app-setup.js|historySnapshots.buildCompositionIntent->legendLayout.captureCompositionIntent',
    'app/app-setup.js|historySnapshots.buildLegendEntryOwners->legendActions.captureLegendEntryOwners',
    'app/app-setup.js|legendActions.commitSpecificRules->featureActions.commitSpecificRules',
    'app/app-setup.js|previewFeatureSearch.resolveOrthogroups->orthogroupActions.resolveOrthogroupDescription',
    'app/app-setup.js|previewFeatureSearch.resolveOrthogroups->orthogroupActions.resolveOrthogroupName',
    'app/app-setup.js|rightDrawerActions.onClose->featureActions.suspendSpecificRulePatternDrafts',
    'app/feature-editor.js|labelActions.setFeatureVisibility->visibilityActions.setFeatureVisibility'
  ]);
  assert.deepEqual(results['owner-graph.state-backdoor.v1'].subjects, [
    'app/app-setup.js|committedDiagramOptions',
    'app/app-setup.js|recordDisplayRows',
    'app/app-setup.js|sessionPreparationBusyReason'
  ]);
  assert.deepEqual(results['heavy-derived.trigger-site.v1'].countsBySubject, {
    'rulePreparation|app/app-setup.js': 3,
    'rulePreparation|app/feature-editor/color-actions.js': 2,
    'rulePreparation|app/feature-editor/label-actions.js': 2,
    'rulePreparation|app/feature-editor/rule-actions.js': 1,
    'rulePreparation|app/feature-editor/svg-actions.js': 1,
    'rulePreparation|app/feature-editor/visibility-actions.js': 2,
    'rulePreparation|app/legend.js': 1,
    'rulePreparation|app/run-analysis.js': 2,
    'rulePreparation|app/svg-styles.js': 1
  });
  assert.deepEqual(
    results['owner-graph.injection-edge.v1'].subjects.filter((subject) => subject.startsWith('app/feature-editor.js|visibilityActions')),
    [
      'app/feature-editor.js|visibilityActions<-featureSvgActions',
      'app/feature-editor.js|visibilityActions<-labelActions',
      'app/feature-editor.js|visibilityActions<-previewRuntime',
      'app/feature-editor.js|visibilityActions<-rulePreparation'
    ]
  );
});

test('the 2026-09-01 baseline had no forward closures, state backdoors, or trigger sites', { skip: revisionExists('e11e1d02') ? false : 'revision history not available in this checkout' }, () => {
  const { summary } = evaluateWebOwnerGraphAt('e11e1d02');
  assert.equal(summary.forwardClosures, 0);
  assert.equal(summary.stateBackdoors, 0);
  assert.equal(summary.triggerSites, 0);
  assert.equal(summary.injectionEdges, 32);
});

REPLAY.forEach(({ pr, merge, detector, subject }) => {
  test(`${pr} (${merge}) is observable as a NEW ${detector} subject against its first parent`, {
    skip: historyAvailable ? false : 'revision history not available in this checkout'
  }, () => {
    const before = new Set(evaluateWebOwnerGraphAt(`${merge}^1`).results[detector].subjects);
    const after = new Set(evaluateWebOwnerGraphAt(merge).results[detector].subjects);
    assert.ok(!before.has(subject), `${pr}: ${subject} already present before the merge`);
    assert.ok(after.has(subject), `${pr}: ${subject} not observed at the merge`);
  });
});
