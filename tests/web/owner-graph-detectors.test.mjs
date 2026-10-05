import assert from 'node:assert/strict';
import test from 'node:test';

import {
  detectWebOwnerGraph,
  summarizeWebOwnerGraph,
  WEB_OWNER_GRAPH_DEFAULTS,
  WEB_OWNER_GRAPH_DETECTOR_IDS,
  WEB_OWNER_GRAPH_DETECTORS
} from '../../tools/web-owner-graph-detectors.mjs';

// Synthetic fixtures: a composition root (app/app-setup.js) that creates
// owners and connects them, two owner modules, and the state module.
const ROOT = [
  "import { createHistoryManager } from '../services/history.js';",
  "import { createLegendManager } from './legend.js';",
  "import { createFeatureEditor } from './feature-editor.js';",
  '',
  'export const createAppSetup = () => {',
  '  let featureActions = null;',
  '  const historySnapshots = createHistorySnapshotService({',
  '    state,',
  '    // reads the legend owner that is created below',
  '    buildLegendEntryOwners: () => legendActions.captureLegendEntryOwners(),',
  '    buildFiles: (ids) => fileStore.collect(ids)',
  '  });',
  '  const history = createHistoryManager({ state, snapshots: historySnapshots, maxActions: 50 });',
  '  const legendActions = createLegendManager({',
  '    state,',
  '    history,',
  '    commitSpecificRules: (...args) => featureActions.commitSpecificRules(...args),',
  "    describe: (entry) => `${entry.caption} (${entry.color})`",
  '  });',
  '  featureActions = createFeatureEditor({ state, history, legendActions, nextTick });',
  '  state.recordDisplayRows = recordDisplayControls.allRows;',
  '  state.committedDiagramOptions = () => getCommittedCanonicalRenderRequest()?.diagramOptions || null;',
  '  state.sessionPreparationBusyReason = null;',
  "  state.errorLog = 'reset';",
  '  if (visibility) await projectFeatureVisibility();',
  '  if (visibility) await projectFeatureVisibility();',
  '  return { history, legendActions, featureActions };',
  '};',
  ''
].join('\n');

const LEGEND_OWNER = [
  'export const createLegendManager = ({',
  '  state,',
  '  history,',
  '  commitSpecificRules = null,',
  '  describe = () => "",',
  '  maxActions = 10,',
  '  trackLayout = "middle"',
  '}) => {',
  '  const orderLegendEntries = (group) => group;',
  '  const sync = () => rulePreparation.prepare();',
  '  const again = () => rulePreparation?.prepareDrawn?.();',
  '  return { orderLegendEntries, sync, again };',
  '};',
  ''
].join('\n');

const FEATURE_OWNER = [
  "import { orderLegendEntries } from './legend/utils.js';",
  'export const createFeatureEditor = ({ state, history, legendActions, nextTick = () => {} }) => {',
  '  const apply = () => {',
  '    orderLegendEntries(group, captions);',
  '    orderLegendEntries(group, captions, { keepFollowed: true });',
  '    orderLegendEntries(group,   captions);',
  '  };',
  '  return { apply, orderLegendEntries: legendActions.orderLegendEntries };',
  '};',
  ''
].join('\n');

const STATE = [
  'export const state = { errorLog: null };',
  'state.featureList = () => [];',
  ''
].join('\n');

const REGISTRY = {
  compositionRoots: ['app/app-setup.js'],
  ownerObjectNames: ['history', 'historySnapshots', 'legendActions', 'featureActions', 'rulePreparation', 'fileStore'],
  ownerObjectNamePattern: WEB_OWNER_GRAPH_DEFAULTS.ownerObjectNamePattern,
  stateModule: 'state.js',
  projectionDomains: [
    { name: 'legend-order', owners: ['app/legend.js'], functions: ['orderLegendEntries'] },
    { name: 'feature-visibility', owners: ['app/feature-editor/visibility-actions.js'], functions: ['projectFeatureVisibility'] }
  ],
  heavyProducers: [{ name: 'rulePreparation', owner: 'app/rule-matching.js', methods: ['prepare', 'prepareDrawn'] }]
};

const SOURCES = new Map([
  ['gbdraw/web/js/app/app-setup.js', ROOT],
  ['gbdraw/web/js/app/legend.js', LEGEND_OWNER],
  ['gbdraw/web/js/app/feature-editor.js', FEATURE_OWNER],
  ['gbdraw/web/js/state.js', STATE]
]);

test('every detector exposes the evidence contract and a stable id list', () => {
  assert.deepEqual(WEB_OWNER_GRAPH_DETECTOR_IDS, [
    'owner-graph.injection-edge.v1',
    'owner-graph.forward-closure.v1',
    'owner-graph.state-backdoor.v1',
    'owner-graph.whole-object-port.v1',
    'projection.call-shape.v1',
    'heavy-derived.trigger-site.v1'
  ]);
  WEB_OWNER_GRAPH_DETECTOR_IDS.forEach((id) => {
    const detector = WEB_OWNER_GRAPH_DETECTORS[id];
    assert.equal(typeof detector.detect, 'function', id);
    assert.equal(typeof detector.encodeSubject, 'function', id);
    assert.equal(typeof detector.subjectCategory, 'string', id);
    assert.ok(Object.isFrozen(detector), id);
  });
});

test('injection edges name every owner object passed to another factory in a root', () => {
  const result = WEB_OWNER_GRAPH_DETECTORS['owner-graph.injection-edge.v1'].detect(SOURCES, REGISTRY);
  assert.deepEqual(result.subjects, [
    'app/app-setup.js|featureActions<-history',
    'app/app-setup.js|featureActions<-legendActions',
    'app/app-setup.js|history<-historySnapshots',
    'app/app-setup.js|legendActions<-history'
  ]);
  // A deferred `let featureActions = null; featureActions = create…` is a
  // factory result; `state`, `nextTick`, and a number are not owner objects.
  assert.deepEqual(result.observedEdges.map(({ consumer, provider, line }) => `${consumer}<-${provider}@${line}`), [
    'featureActions<-history@20', 'featureActions<-legendActions@20', 'history<-historySnapshots@13', 'legendActions<-history@14'
  ]);
  assert.equal(WEB_OWNER_GRAPH_DETECTORS['owner-graph.injection-edge.v1'].encodeSubject({
    path: 'gbdraw/web/js/app/app-setup.js', consumer: 'legendActions', provider: 'history'
  }), 'app/app-setup.js|legendActions<-history');
  // The same sources in another order give the same subjects.
  const reversed = new Map([...SOURCES].reverse());
  assert.deepEqual(WEB_OWNER_GRAPH_DETECTORS['owner-graph.injection-edge.v1'].detect(reversed, REGISTRY).subjects, result.subjects);
});

test('forward closures are closures that reach an owner created later in the root', () => {
  const result = WEB_OWNER_GRAPH_DETECTORS['owner-graph.forward-closure.v1'].detect(SOURCES, REGISTRY);
  assert.deepEqual(result.subjects, [
    'app/app-setup.js|historySnapshots.buildLegendEntryOwners->legendActions.captureLegendEntryOwners',
    'app/app-setup.js|legendActions.commitSpecificRules->featureActions.commitSpecificRules'
  ]);
  // `buildFiles` reaches `fileStore`, which is not created in this root, and
  // `describe` reaches no owner: neither is a forward closure.
  assert.deepEqual(result.observedClosures.map(({ consumer, port, provider, method }) => `${consumer}.${port}->${provider}.${method}`), [
    'historySnapshots.buildLegendEntryOwners->legendActions.captureLegendEntryOwners',
    'legendActions.commitSpecificRules->featureActions.commitSpecificRules'
  ]);
  // Reordering creation removes the closure.
  const reordered = new Map(SOURCES);
  reordered.set('gbdraw/web/js/app/app-setup.js', ROOT.replace(
    '  const historySnapshots = createHistorySnapshotService({',
    "  const legendActionsEarly = createLegendManager({ state });\n  const historySnapshots = createHistorySnapshotService({"
  ).replace('legendActions.captureLegendEntryOwners()', 'legendActionsEarly.captureLegendEntryOwners()'));
  assert.deepEqual(WEB_OWNER_GRAPH_DETECTORS['owner-graph.forward-closure.v1'].detect(reordered, REGISTRY).subjects, [
    'app/app-setup.js|legendActions.commitSpecificRules->featureActions.commitSpecificRules'
  ]);
});

test('state backdoors are functions or owner members assigned into state outside the state module', () => {
  const result = WEB_OWNER_GRAPH_DETECTORS['owner-graph.state-backdoor.v1'].detect(SOURCES, REGISTRY);
  assert.deepEqual(result.subjects, [
    'app/app-setup.js|committedDiagramOptions',
    'app/app-setup.js|recordDisplayRows'
  ]);
  assert.deepEqual(result.observedAssignments, [
    { path: 'app/app-setup.js', name: 'committedDiagramOptions', line: 22, kind: 'function' },
    { path: 'app/app-setup.js', name: 'recordDisplayRows', line: 21, kind: 'owner-member' }
  ]);
  // `= null`, a string literal, and state.js's own assignments are not backdoors.
});

test('whole-object ports are owner objects taken by factories outside the roots', () => {
  const result = WEB_OWNER_GRAPH_DETECTORS['owner-graph.whole-object-port.v1'].detect(SOURCES, REGISTRY);
  assert.deepEqual(result.subjects, [
    'app/feature-editor.js|createFeatureEditor|history',
    'app/feature-editor.js|createFeatureEditor|legendActions',
    'app/legend.js|createLegendManager|history'
  ]);
  // `commitSpecificRules` and `describe` are ports (functions); `maxActions`
  // and `trackLayout` match the name pattern but are never factory results.
  assert.deepEqual(result.observedParameters.map(({ factory, parameter }) => `${factory}:${parameter}`), [
    'createFeatureEditor:history', 'createFeatureEditor:legendActions', 'createLegendManager:history'
  ]);
});

test('projection call shapes count distinct call lines outside the domain owners', () => {
  const result = WEB_OWNER_GRAPH_DETECTORS['projection.call-shape.v1'].detect(SOURCES, REGISTRY);
  assert.deepEqual(result.shapesByDomain, { 'legend-order': 2, 'feature-visibility': 1 });
  assert.deepEqual(result.subjects, [
    'feature-visibility|app/app-setup.js|if (visibility) await projectFeatureVisibility();',
    'legend-order|app/feature-editor.js|orderLegendEntries(group, captions);',
    'legend-order|app/feature-editor.js|orderLegendEntries(group, captions, { keepFollowed: true });'
  ]);
  // Whitespace differences, repeated lines, the import, the export mapping
  // (`orderLegendEntries: legendActions.orderLegendEntries`), and the owner's
  // own definition do not add shapes.
  assert.equal(result.observedShapes.filter(({ domain }) => domain === 'legend-order').length, 3);
});

test('heavy-derived trigger sites count producer calls per module outside the owner', () => {
  const result = WEB_OWNER_GRAPH_DETECTORS['heavy-derived.trigger-site.v1'].detect(SOURCES, REGISTRY);
  assert.deepEqual(result.subjects, ['rulePreparation|app/legend.js']);
  assert.deepEqual(result.countsBySubject, { 'rulePreparation|app/legend.js': 2 });
  assert.equal(result.siteCount, 2);
  assert.deepEqual(result.observedSites.map(({ method, line }) => `${method}@${line}`), ['prepare@10', 'prepareDrawn@11']);
  const owner = new Map(SOURCES);
  owner.set('gbdraw/web/js/app/rule-matching.js', 'export const createRulePreparation = () => ({ prepare: () => rulePreparation.prepare() });\n');
  assert.equal(WEB_OWNER_GRAPH_DETECTORS['heavy-derived.trigger-site.v1'].detect(owner, REGISTRY).siteCount, 2);
});

test('detectWebOwnerGraph runs every detector and summarizeWebOwnerGraph counts subjects', () => {
  const results = detectWebOwnerGraph(SOURCES, REGISTRY);
  assert.deepEqual(Object.keys(results), WEB_OWNER_GRAPH_DETECTOR_IDS);
  assert.deepEqual(summarizeWebOwnerGraph(results), {
    injectionEdges: 4,
    forwardClosures: 2,
    stateBackdoors: 2,
    wholeObjectPorts: 3,
    projectionShapes: { 'legend-order': 2, 'feature-visibility': 1 },
    triggerSites: 2,
    triggerModules: 1
  });
  // The default registry applies when none is given.
  const withDefaults = detectWebOwnerGraph(SOURCES);
  assert.ok(Array.isArray(withDefaults['owner-graph.injection-edge.v1'].subjects));
  assert.ok(Object.isFrozen(WEB_OWNER_GRAPH_DEFAULTS));
  assert.deepEqual(WEB_OWNER_GRAPH_DEFAULTS.compositionRoots, [
    'app/app-setup.js', 'app/feature-editor.js', 'app/legend.js', 'app/legend-layout.js'
  ]);
});

test('the detectors read text only: a module that would throw when imported is still observed', () => {
  const sources = new Map(SOURCES);
  sources.set('gbdraw/web/js/app/app-setup.js', `throw new Error('never run');\n${ROOT}`);
  const results = detectWebOwnerGraph(sources, REGISTRY);
  assert.equal(results['owner-graph.forward-closure.v1'].subjects.length, 2);
});
