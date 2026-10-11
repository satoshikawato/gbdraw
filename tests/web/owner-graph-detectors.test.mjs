import assert from 'node:assert/strict';
import test from 'node:test';

import {
  detectWebOwnerGraph,
  summarizeWebOwnerGraph,
  WEB_OWNER_GRAPH_DEFAULTS,
  WEB_OWNER_GRAPH_DETECTOR_IDS,
  WEB_OWNER_GRAPH_DETECTORS,
  webLayerOf
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
    'owner-graph.forward-closure.v2',
    'owner-graph.state-backdoor.v1',
    'owner-graph.whole-object-port.v1',
    'projection.call-shape.v1',
    'heavy-derived.trigger-site.v1',
    'heavy-derived.trigger-site.v2',
    'layer.import-direction.v1',
    'owner-graph.identity-from-display.v1'
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

// A composition root whose closures sit below the top level: a named helper
// that an owner's port reaches, an object literal inside an argument, a
// method, a watcher, and function variables bound late.
const NESTED_ROOT = [
  'export const createAppSetup = () => {',
  '  let blocksEditor = () => false;',
  '  let refreshCanvas = () => {};',
  '  let neverRead = () => {};',
  '  let counter = 0;',
  '  let lazyTimer = null;',
  '  const restoreWithDrafts = async (restore, ...args) => {',
  '    const drafts = featureActions.captureDrafts();',
  '    try { return await restore(...args); } finally { featureActions?.restoreDrafts(drafts); }',
  '  };',
  '  const history = createHistoryManager({',
  '    applyIntent: (...args) => restoreWithDrafts(applyHistoryIntent, ...args),',
  '    describe: () => history.label',
  '  });',
  '  const legendLayout = createLegendLayout({',
  '    lifecycle: {',
  '      beforeDrag: () => alignmentActions.beforeDrag(),',
  '      afterDrag: (options) => alignmentActions.afterDrag(options)',
  '    }',
  '  });',
  '  const drawer = createDrawer({',
  "    getReason: () => blocksEditor() ? 'blocked' : ''",
  '  });',
  '  previewRuntime.configureBinder({',
  '    install(context) {',
  '      refreshCanvas();',
  '    }',
  '  });',
  '  watch(selected, () => refreshCanvas());',
  '  watch(other, () => refreshCanvas());',
  '  const afterOnly = () => featureActions.afterCreation();',
  '  const chained = () => afterOnly();',
  '  const unusedHelper = () => featureActions.neverPassed();',
  '  const sidebar = createSidebar({ unusedHelper: null });',
  '  const shadowed = (featureActions) => featureActions.local();',
  '  const wrapped = createWrapper(shadowed);',
  '  const laterTimer = () => { lazyTimer = () => {}; };',
  '  const readsLater = () => featureActions.tdzCase();',
  '  const bridge = createBridge({ callback: () => viaLater() });',
  '  const featureActions = createFeatureEditor({ history });',
  '  const alignmentActions = createAlignment({ legendLayout });',
  '  const viaLater = () => readsLater();',
  '  const earlier = () => featureActions.readsAnOwnerCreatedEarlier();',
  '  const reuse = createReuse({ earlier });',
  '  blocksEditor = () => Boolean(alignmentActions.dialogOpen);',
  '  refreshCanvas = () => { alignmentActions.refresh(); };',
  '  neverRead = () => {};',
  '  counter = 5;',
  '  return { history, chained, wrapped, sidebar, drawer, laterTimer, reuse };',
  '};',
  ''
].join('\n');
const NESTED_SOURCES = new Map([['gbdraw/web/js/app/app-setup.js', NESTED_ROOT]]);

test('v2 forward closures find closures at any depth that an owner can reach before the owner they read exists', () => {
  const detector = WEB_OWNER_GRAPH_DETECTORS['owner-graph.forward-closure.v2'];
  const result = detector.detect(NESTED_SOURCES, REGISTRY);
  assert.deepEqual(result.subjects, [
    // A named helper that createHistoryManager's port names (reached through the port).
    'app/app-setup.js|restoreWithDrafts->featureActions.captureDrafts',
    'app/app-setup.js|restoreWithDrafts->featureActions.restoreDrafts',
    // An object literal inside an argument, which v1's entry check does not open.
    'app/app-setup.js|beforeDrag->alignmentActions.beforeDrag',
    'app/app-setup.js|afterDrag->alignmentActions.afterDrag',
    // A late-bound function variable read by a port, a method, and watchers.
    'app/app-setup.js|getReason->blocksEditor',
    'app/app-setup.js|install->refreshCanvas',
    'app/app-setup.js|watch(callback)->refreshCanvas'
  ].sort());
  assert.deepEqual(detector.encodeSubject({
    path: 'gbdraw/web/js/app/app-setup.js', consumer: 'getReason', provider: 'blocksEditor', method: ''
  }), 'app/app-setup.js|getReason->blocksEditor');
  assert.deepEqual(detector.encodeSubject({
    path: 'gbdraw/web/js/app/app-setup.js', consumer: 'beforeDrag', provider: 'alignmentActions', method: 'beforeDrag'
  }), 'app/app-setup.js|beforeDrag->alignmentActions.beforeDrag');
  // Both watchers share one subject; the references keep their lines.
  assert.deepEqual(result.observedReferences.filter(({ consumer }) => consumer === 'watch(callback)').map(({ line }) => line), [29, 30]);
  // Not subjects: `afterOnly` and `chained` are reachable only from the returned
  // bindings after every owner exists; `unusedHelper` is a property key, not a
  // use; `shadowed` declares its own `featureActions`; `earlier` reads an
  // owner created before it; `readsLater` is named only by `viaLater`, a `const`
  // function defined after `featureActions` exists, which cannot run before
  // that; `neverRead`, `counter` (not a function), and `lazyTimer` (assigned
  // inside a function) are not late-bound function variables.
  ['afterOnly', 'chained', 'unusedHelper', 'shadowed', 'earlier', 'readsLater', 'tdzCase', 'neverRead', 'counter', 'lazyTimer'].forEach((name) => {
    assert.ok(!result.subjects.some((subject) => subject.includes(name)), name);
  });
});

test('v2 forward closures keep every v1 subject and add nothing when the owner is created first', () => {
  const v1 = WEB_OWNER_GRAPH_DETECTORS['owner-graph.forward-closure.v1'].detect(SOURCES, REGISTRY).subjects;
  const v2 = WEB_OWNER_GRAPH_DETECTORS['owner-graph.forward-closure.v2'].detect(SOURCES, REGISTRY);
  assert.deepEqual(v2.subjects, v1);
  assert.deepEqual(v2.observedReferences, []);
  // The same closures with the late owner created first are not forward.
  const early = NESTED_ROOT
    .replace("  const featureActions = createFeatureEditor({ history });\n", '')
    .replace("  const alignmentActions = createAlignment({ legendLayout });\n", '')
    .replace("  let blocksEditor = () => false;\n", "  const featureActions = createFeatureEditor({ history });\n  const alignmentActions = createAlignment({ legendLayout });\n  let blocksEditor = () => false;\n")
    .replace("  blocksEditor = () => Boolean(alignmentActions.dialogOpen);\n  refreshCanvas = () => { alignmentActions.refresh(); };\n", '');
  const result = WEB_OWNER_GRAPH_DETECTORS['owner-graph.forward-closure.v2'].detect(new Map([['gbdraw/web/js/app/app-setup.js', early]]), REGISTRY);
  assert.deepEqual(result.subjects, []);
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

// A root that hands producer functions to a factory instead of the producer.
const PORT_ROOT = [
  'export const createFeatureEditor = ({ rulePreparation }) => {',
  '  const svg = createSvgOwner({',
  '    state,',
  '    runWithDrawn: rulePreparation.runDrawn,',
  '    runWithRules: rulePreparation.run,',
  '    onChanges: rulePreparation.notifyChanges,',
  '    isReady: rulePreparation.isPrepared',
  '  });',
  '  const other = createOtherOwner({ unrelated: someOtherObject.run, runWithDrawn: rulePreparation.runDrawn });',
  '  const missing = createMissingOwner({ runWithDrawn: rulePreparation.runDrawn });',
  '  return { svg, other, missing };',
  '};',
  ''
].join('\n');
const PORT_OWNER = [
  'export const createSvgOwner = ({',
  '  state,',
  '  runWithDrawn = (commit) => commit(),',
  '  runWithRules: runRules,',
  '  onChanges = () => {},',
  '  isReady = () => true',
  '}) => {',
  '  const open = () => runWithDrawn(() => draw());',
  '  const refresh = () => runRules(rules, () => draw());',
  '  const guarded = () => runWithDrawn?.(() => draw());',
  '  const notify = () => onChanges(candidate);',
  '  return { open, refresh, guarded, notify, ready: isReady() };',
  '};',
  ''
].join('\n');
const OTHER_OWNER = [
  'export const createOtherOwner = ({ unrelated, runWithDrawn }) => {',
  '  const run = () => unrelated();',
  '  return { run, expose: runWithDrawn };',
  '};',
  ''
].join('\n');
const PORT_REGISTRY = {
  ...REGISTRY,
  compositionRoots: ['app/feature-editor.js'],
  heavyProducers: [{
    name: 'rulePreparation', owner: 'app/rule-matching.js', methods: ['prepare', 'prepareDrawn', 'run'], v2Methods: ['runDrawn']
  }]
};
const PORT_SOURCES = new Map([
  ['gbdraw/web/js/app/feature-editor.js', PORT_ROOT],
  ['gbdraw/web/js/app/svg-owner.js', PORT_OWNER],
  ['gbdraw/web/js/app/other-owner.js', OTHER_OWNER],
  ['gbdraw/web/js/app/rule-matching.js', 'export const createRulePreparation = () => ({ run: () => rulePreparation.run() });\n']
]);

test('v2 trigger sites count calls of a producer port in the module that receives it', () => {
  const v1 = WEB_OWNER_GRAPH_DETECTORS['heavy-derived.trigger-site.v1'].detect(PORT_SOURCES, PORT_REGISTRY);
  // v1 sees no call of `rulePreparation.` in any module: the port hides them.
  assert.deepEqual(v1.countsBySubject, {});
  const detector = WEB_OWNER_GRAPH_DETECTORS['heavy-derived.trigger-site.v2'];
  const result = detector.detect(PORT_SOURCES, PORT_REGISTRY);
  assert.deepEqual(result.countsBySubject, { 'rulePreparation|app/svg-owner.js': 3 });
  assert.deepEqual(result.subjects, ['rulePreparation|app/svg-owner.js']);
  assert.equal(result.siteCount, 3);
  assert.deepEqual(result.observedSites.map(({ method, via, line }) => `${method}@${line} via ${via}`), [
    'runDrawn@8 via app/feature-editor.js:createSvgOwner.runWithDrawn',
    'run@9 via app/feature-editor.js:createSvgOwner.runWithRules',
    'runDrawn@10 via app/feature-editor.js:createSvgOwner.runWithDrawn'
  ]);
  assert.equal(detector.encodeSubject({ producer: 'rulePreparation', path: 'app/svg-owner.js' }), 'rulePreparation|app/svg-owner.js');
  // Not sites: the default parameter (`runWithDrawn = (commit) => commit()`),
  // the binding itself, `notifyChanges` and `isPrepared` (not run methods), a
  // port the module only hands on (`other-owner.js`), a factory that is not in
  // the sources, and another object's `run`.
  assert.ok(!result.subjects.some((subject) => /other-owner|rule-matching|feature-editor/.test(subject)));
  // A direct call counts as in v1, and `runDrawn` is a producer method in v2 only.
  const direct = new Map(PORT_SOURCES);
  direct.set('gbdraw/web/js/app/direct.js', 'export const x = () => rulePreparation.runDrawn(() => 1);\nexport const y = () => rulePreparation?.prepareDrawn?.();\n');
  assert.deepEqual(detector.detect(direct, PORT_REGISTRY).countsBySubject, {
    'rulePreparation|app/direct.js': 2, 'rulePreparation|app/svg-owner.js': 3
  });
  assert.deepEqual(WEB_OWNER_GRAPH_DETECTORS['heavy-derived.trigger-site.v1'].detect(direct, PORT_REGISTRY).countsBySubject, {
    'rulePreparation|app/direct.js': 1
  });
  // The default registry names `runDrawn` as a v2 producer method.
  assert.ok(WEB_OWNER_GRAPH_DEFAULTS.heavyProducers.find(({ name }) => name === 'rulePreparation').v2Methods.includes('runDrawn'));
});

const LAYER_SOURCES = new Map([
  ['gbdraw/web/js/utils/leaf.js', [
    "import { helper } from './other-leaf.js';",
    "import { VERSION } from '../config.js';",
    "import { normalize } from '../services/error-normalization.js';",
    "import { service } from '../services/svc.js';",
    ''
  ].join('\n')],
  ['gbdraw/web/js/utils/other-leaf.js', 'export const helper = 1;\n'],
  ['gbdraw/web/js/config.js', "import { helper } from './utils/other-leaf.js';\nexport const VERSION = 1;\n"],
  ['gbdraw/web/js/services/error-normalization.js', 'export const normalize = (value) => value;\n'],
  ['gbdraw/web/js/services/svc.js', [
    "import { helper } from '../utils/other-leaf.js';",
    "import { sibling } from './sibling.js';",
    "import { state } from '../state.js';",
    "import { load } from './config.js';",
    "import { a, b as renamed, default as fallback } from '../app/owner.js';",
    "export { a } from '../app/owner.js';",
    "import * as everything from '../app/ns.js';",
    "import '../app/side.js';",
    "import defaultOwner from '../app/default-owner.js';",
    "import Vue from 'vue';",
    "/** @import { OwnerType } from '../app/typed-only.js' */",
    "// import { skipped } from '../app/line-comment.js';",
    "/* import { skipped } from '../app/block-comment.js'; */",
    "const text = \"import { skipped } from '../app/in-string.js'\";",
    "export const lazy = () => import('../app/dynamic.js');",
    ''
  ].join('\n')],
  ['gbdraw/web/js/services/sibling.js', 'export const sibling = 1;\n'],
  ['gbdraw/web/js/services/config.js', [
    "import { state } from '../state.js';",
    "import { service } from './svc.js';",
    "import { a } from '../app/owner.js';",
    ''
  ].join('\n')],
  ['gbdraw/web/js/state.js', [
    "import { service } from './services/svc.js';",
    "import { helper } from './utils/other-leaf.js';",
    "import { a } from './app/owner.js';",
    ''
  ].join('\n')],
  ['gbdraw/web/js/app/owner.js', [
    "import { state } from '../state.js';",
    "import { load } from '../services/config.js';",
    "import { createAppSetup } from './app-setup.js';",
    "import { other } from './other-owner.js';",
    ''
  ].join('\n')],
  ['gbdraw/web/js/app/other-owner.js', 'export const other = 1;\n'],
  ['gbdraw/web/js/app/app-setup.js', "import { a } from './owner.js';\nimport { createLegendManager } from './legend.js';\nexport const createAppSetup = () => a;\n"],
  ['gbdraw/web/js/app/legend.js', "import { state } from '../state.js';\nexport const createLegendManager = () => state;\n"],
  ['gbdraw/web/js/app.js', "import { createAppSetup } from './app/app-setup.js';\nimport { state } from './state.js';\n"],
  ['gbdraw/web/js/workers/worker.js', "import { a } from '../app/owner.js';\n"],
  // Targets that exist, so a missed mask would show as a subject.
  ...['default-owner', 'dynamic', 'ns', 'side', 'typed-only', 'line-comment', 'block-comment', 'in-string']
    .map((name) => [`gbdraw/web/js/app/${name}.js`, 'export const x = 1;\n'])
]);

test('webLayerOf ranks modules by the R13 layers and leaves workers unranked', () => {
  const ranks = (paths) => paths.map((path) => webLayerOf(path));
  assert.deepEqual(ranks(['utils/zip.js', 'config.js', 'web-ux-profile.js', 'mode-profiles.generated.js', 'mode-profiles.js', 'mode-scoped-settings.generated.js']), [0, 0, 0, 0, 0, 0]);
  assert.deepEqual(ranks(['services/svg-serialization.js', 'services/losat.js', 'services/error-normalization.js']), [1, 1, 1]);
  assert.equal(webLayerOf('state.js'), 2);
  assert.deepEqual(ranks(['services/config.js', 'services/reset.js']), [3, 3]);
  assert.deepEqual(ranks(['app/run-analysis.js', 'app/legend/utils.js', 'app/feature-editor/label-actions.js']), [4, 4, 4]);
  assert.deepEqual(ranks(['app/app-setup.js', 'app/feature-editor.js', 'app/legend.js', 'app/legend-layout.js']), [5, 5, 5, 5]);
  assert.deepEqual(ranks(['app.js', 'components.js']), [6, 6]);
  assert.deepEqual(ranks(['workers/losat-worker.js', 'package.json']), [null, null]);
  // The source prefix is accepted, and every module under the registered roots is ranked.
  assert.equal(webLayerOf('gbdraw/web/js/state.js'), 2);
  assert.equal(webLayerOf('gbdraw/web/js/services/config.js'), 3);
});

test('layer import direction reports a module that imports a higher layer, and only that', () => {
  const detector = WEB_OWNER_GRAPH_DETECTORS['layer.import-direction.v1'];
  assert.equal(detector.subjectCategory, 'layer-import');
  const result = detector.detect(LAYER_SOURCES);
  // Counts are the distinct imported names: `a`, `b`, `default` (a repeated name counts once);
  // a namespace, a side-effect, and a dynamic import count one each.
  assert.deepEqual(result.countsBySubject, {
    'app/owner.js->app/app-setup.js': 1,
    'services/config.js->app/owner.js': 1,
    'services/svc.js->app/default-owner.js': 1,
    'services/svc.js->app/dynamic.js': 1,
    'services/svc.js->app/ns.js': 1,
    'services/svc.js->app/owner.js': 3,
    'services/svc.js->app/side.js': 1,
    'services/svc.js->services/config.js': 1,
    'services/svc.js->state.js': 1,
    'state.js->app/owner.js': 1,
    'utils/leaf.js->services/error-normalization.js': 1,
    'utils/leaf.js->services/svc.js': 1
  });
  assert.deepEqual(result.subjects, Object.keys(result.countsBySubject));
  assert.equal(result.nameCount, 14);
  // A state-free service never imports state.js, an owner never imports a composition root.
  assert.deepEqual(result.observedImports.filter(({ target }) => target === 'state.js' || target === 'app/app-setup.js')
    .map(({ path, target, layers }) => `${path}->${target} (${layers})`), [
    'app/owner.js->app/app-setup.js (owner module -> composition root)',
    'services/svc.js->state.js (state-free service -> state)'
  ]);
  // Not violations: same layer, a lower layer, a leaf importing a leaf, an entry module importing
  // anything, a composition root importing an owner, a bare specifier, a JSDoc `@import`,
  // comments, a string, and a module under workers/.
  const flagged = new Set(result.observedImports.map(({ path }) => path));
  ['config.js', 'utils/other-leaf.js', 'services/sibling.js', 'app/other-owner.js', 'app/app-setup.js', 'app/legend.js', 'app.js', 'workers/worker.js']
    .forEach((path) => assert.ok(!flagged.has(path), path));
  const targets = result.observedImports.map(({ target }) => target);
  ['app/typed-only.js', 'app/line-comment.js', 'app/block-comment.js', 'app/in-string.js', 'config.js']
    .forEach((target) => assert.ok(!targets.includes(target), target));
  assert.equal(detector.encodeSubject({ path: 'gbdraw/web/js/services/svc.js', target: 'app/owner.js' }), 'services/svc.js->app/owner.js');
  // The same sources in another order give the same subjects.
  assert.deepEqual(detector.detect(new Map([...LAYER_SOURCES].reverse())).subjects, result.subjects);
  // A graph that only goes down reports nothing.
  assert.deepEqual(detector.detect(new Map([
    ['app/owner.js', "import { x } from '../services/svc.js';\nimport { y } from '../utils/leaf.js';\n"],
    ['services/svc.js', "import { y } from '../utils/leaf.js';\n"],
    ['utils/leaf.js', 'export const y = 1;\n']
  ])).subjects, []);
  // The summary counts importer->target pairs.
  assert.equal(summarizeWebOwnerGraph(detectWebOwnerGraph(LAYER_SOURCES)).layerImports, 12);
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
    triggerModules: 1,
    layerImports: 0
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

// Identity read from a display value: one module per join form, the three
// known instances (OV-243, OV-288, review M1), and the comparisons that are
// not joins.
const IDENTITY_FROM_DISPLAY = 'owner-graph.identity-from-display.v1';
const IDENTITY_SOURCES = new Map([
  // OV-243: a mounted Legend row matched to an intent entry by the caption
  // and color the Result shows.
  ['gbdraw/web/js/app/legend/entry-actions.js', [
    'export const readMountedLegend = (svg, previousEntries) => {',
    '  const rows = [];',
    "  svg.querySelectorAll('g[data-legend-key]').forEach((entryGroup) => {",
    "    const caption = entryGroup.getAttribute('data-legend-key');",
    "    const color = entryGroup.querySelector('path').getAttribute('fill');",
    '    const existingEntry = previousEntries.find((entry) => (',
    '      entry.caption === caption',
    '      && normalizedColor(entry.color) === normalizedColor(color)',
    '    ));',
    '    rows.push(existingEntry);',
    '  });',
    '  return rows;',
    '};',
    // Review M1: a caption fallback that names a row by its shown key.
    'export const createEntryActions = ({ entries, dormant }) => {',
    '  const intentEntry = (test) => entries.find(test) || dormant.find(test);',
    '  const namedBy = (row) => intentEntry((each) => legendCaption(each) === row.shownKey);',
    '  const untypedNamedBy = (row) => intentEntry((each) => legendCaption(each) === row.key);',
    '  return { namedBy, untypedNamedBy };',
    '};',
    ''
  ].join('\n')],
  // OV-288: stroke membership read from the Legend row color.
  ['gbdraw/web/js/services/legend-svg.js', [
    'export const legendRowFeatureIds = (entry, drawnFills) => {',
    '  const color = paintKey(entry?.color);',
    '  const reached = new Set();',
    '  for (const [id, fill] of drawnFills) {',
    "    if (color && color !== 'none' && paintKey(fill) === color) reached.add(id);",
    '  }',
    '  return [...reached];',
    '};',
    ''
  ].join('\n')],
  // Lookup keys.
  ['gbdraw/web/js/app/legend/stroke-actions.js', [
    'export const strokeOf = (entry, strokesByColor) => strokesByColor.get(entry.color);',
    'export const rowsByColor = (rows) => new Map(rows.map((row) => [row.color, row]));',
    'export const overrideOfShownRow = (overrides, row) => overrides[row.textContent];',
    "export const rowOf = (svg, caption) => svg.querySelector(`[data-legend-key=\"${caption}\"]`);",
    'export const groupCaptions = (entries) => groupBy(entries, legendCaption);',
    ''
  ].join('\n')],
  // Not joins: literals, typeof, an identity in the same conjunction,
  // change detection of one field in a loop, caption dedupe, and a
  // comparison in straight-line code.
  ['gbdraw/web/js/app/legend/layout-actions.js', [
    "export const unpainted = (entries) => entries.filter((entry) => entry.color === 'none');",
    'export const typed = (entries, kind) => entries.filter((entry) => typeof entry.color === kind);',
    'export const sameRow = (entries, key, color) => entries.find((entry) => entry.originalCaption === key && entry.color === color);',
    'export const changed = (entries, previous) => {',
    '  let count = 0;',
    '  entries.forEach((entry, index) => { if (entry.color !== previous[index].color) count += 1; });',
    '  return count;',
    '};',
    'export const uniqueCaptions = (entries) => {',
    '  const seen = new Set();',
    '  return entries.filter((entry) => !seen.has(entry.caption) && seen.add(entry.caption));',
    '};',
    'export const recolor = (entry, nextColor) => {',
    '  if (entry.color !== nextColor) entry.color = nextColor;',
    '  return entry;',
    '};',
    ''
  ].join('\n')],
  // Outside the registered modules.
  ['gbdraw/web/js/app/track-slots.js', 'export const slotOf = (slots, color) => slots.find((slot) => slot.color === color);\n']
]);

test('identity from display reports paint and shown joins per function and inventories caption keys', () => {
  const result = WEB_OWNER_GRAPH_DETECTORS[IDENTITY_FROM_DISPLAY].detect(IDENTITY_SOURCES);
  assert.deepEqual(result.countsBySubject, {
    'app/legend/entry-actions.js|namedBy|shown-join': 1,
    'app/legend/entry-actions.js|readMountedLegend|shown-join': 2,
    'app/legend/stroke-actions.js|overrideOfShownRow|shown-join': 1,
    'app/legend/stroke-actions.js|rowsByColor|paint-join': 1,
    'app/legend/stroke-actions.js|strokeOf|paint-join': 1,
    'services/legend-svg.js|legendRowFeatureIds|paint-join': 1
  });
  assert.deepEqual(result.subjects, Object.keys(result.countsBySubject));
  assert.equal(result.siteCount, 7);
  const site = (path, fn) => result.observedSites
    .filter((record) => record.path === path && record.function === fn)
    .map((record) => `${record.class} ${record.form}:${record.field}`);
  // Search callbacks; the mounted row's caption and color are locals bound
  // from what the Result shows.
  assert.deepEqual(site('app/legend/entry-actions.js', 'readMountedLegend'), [
    'shown-join search:caption',
    'shown-join search:color'
  ]);
  // A predicate passed to another function; `row.key` is an untyped name,
  // so only the shown field makes the M1 form a subject.
  assert.deepEqual(site('app/legend/entry-actions.js', 'namedBy'), ['shown-join predicate:legendCaption()']);
  assert.deepEqual(site('app/legend/entry-actions.js', 'untypedNamedBy'), ['caption-key predicate:legendCaption()']);
  // A loop that joins two different paint fields.
  assert.deepEqual(site('services/legend-svg.js', 'legendRowFeatureIds'), ['paint-join loop-join:fill']);
  assert.deepEqual(result.observedSites.filter((record) => record.path === 'app/legend/stroke-actions.js')
    .map((record) => `${record.function} ${record.class} ${record.form}`), [
    'strokeOf paint-join key-get',
    'rowsByColor paint-join key-entries',
    'overrideOfShownRow shown-join key-index',
    'rowOf caption-key key-selector',
    'groupCaptions caption-key key-function'
  ]);
  // The non-joins: literals, typeof, and straight-line code are no site; the
  // identity-qualified comparison, same-field change detection, and caption
  // dedupe are inventory only.
  assert.deepEqual(result.observedSites.filter((record) => record.path === 'app/legend/layout-actions.js')
    .map((record) => `${record.function} ${record.class} ${record.form}`), [
    'sameRow id-qualified search',
    'changed paint-compare loop-join',
    'uniqueCaptions caption-key key-has',
    'uniqueCaptions caption-key key-add'
  ]);
  assert.deepEqual(result.inventory, {
    'caption-key|app/legend/entry-actions.js': 1,
    'caption-key|app/legend/stroke-actions.js': 2,
    'caption-key|app/legend/layout-actions.js': 2,
    'id-qualified|app/legend/layout-actions.js': 1,
    'paint-compare|app/legend/layout-actions.js': 1
  });
  assert.equal(result.observedSites.some((record) => record.path === 'app/track-slots.js'), false);
  // Subjects encode as path|function|class.
  const record = result.observedSites.find((each) => each.function === 'legendRowFeatureIds');
  assert.equal(WEB_OWNER_GRAPH_DETECTORS[IDENTITY_FROM_DISPLAY].encodeSubject(record), 'services/legend-svg.js|legendRowFeatureIds|paint-join');
});

test('identity from display scopes a shown local to the function that binds it and reads its modules from the registry', () => {
  const sources = new Map([
    ['gbdraw/web/js/app/legend/entry-actions.js', [
      "export const shownCaption = (group) => { const caption = group.getAttribute('data-legend-key'); return caption; };",
      'export const intentOf = (entries, caption) => entries.find((entry) => entry.caption === caption);',
      ''
    ].join('\n')]
  ]);
  const detect = (registry) => WEB_OWNER_GRAPH_DETECTORS[IDENTITY_FROM_DISPLAY].detect(sources, registry);
  // `caption` in intentOf is its own parameter, not the shown local.
  assert.deepEqual(detect().observedSites.map((record) => `${record.function} ${record.class}`), ['intentOf caption-key']);
  assert.deepEqual(detect().subjects, []);
  const outside = { identityFromDisplay: { ...WEB_OWNER_GRAPH_DEFAULTS.identityFromDisplay, modules: ['services/'] } };
  assert.deepEqual(detect(outside).observedSites, []);
  assert.ok(Object.isFrozen(WEB_OWNER_GRAPH_DEFAULTS.identityFromDisplay));
});

// Review fixes of v1 before a baseline references it: optional chaining,
// the reach of a shown local, a comparison that qualified itself, shown
// locals with any name, write-if-different, ternaries, and one-line bodies.
const IDENTITY_REVIEW_SOURCE = [
  // H1: `?.` stays inside the operand on both sides of a comparison.
  'export const hasShownCaption = (texts, finalCaption) => texts.find((text) => text.textContent?.trim() === finalCaption);',
  'export const optionalIdentity = (entries, c, k) => entries.find((e) => e.color === c && e?.originalCaption === k);',
  // M1: a shown local reaches the end of the function that binds it, not a
  // sibling callback, and a nested parameter of the same name shadows it.
  'export const siblings = (groups, captions, intents) => {',
  '  let count = 0;',
  '  groups.forEach((group) => {',
  "    const caption = group.getAttribute('data-legend-key');",
  '    if (intents.some((intent) => intent.caption === caption)) count += 1;',
  '  });',
  '  captions.forEach((caption) => {',
  '    if (intents.some((intent) => intent.caption === caption)) count += 1;',
  '  });',
  '  return count;',
  '};',
  'export const outerName = (groups, caption, intents) => {',
  "  groups.forEach((group) => { const caption = group.getAttribute('data-legend-key'); show(caption); });",
  '  return intents.find((intent) => intent.caption === caption);',
  '};',
  'export const shadowed = (group, intents, captions, seen) => {',
  "  const caption = group.getAttribute('data-legend-key');",
  '  const hit = intents.find((intent) => intent.caption === caption);',
  '  return [hit, captions.filter((caption) => seen.has(caption))];',
  '};',
  // M2: an identity compared with a shown value is the join itself.
  "export const ownIdentity = (entries, group) => entries.find((entry) => entry && entry.originalCaption === group.getAttribute('data-legend-key'));",
  // M3: a local bound from a shown read, whatever its name, as the operand
  // value; not as an argument of another call.
  'export const overrideOfKey = (overrides, group) => {',
  "  const legendKey = group.getAttribute('data-legend-key');",
  '  return overrides[legendKey];',
  '};',
  "const keyOf = (group) => String(group.getAttribute('data-legend-key') || '').trim();",
  'export const rankOf = (rank, group) => rank.get(keyOf(group));',
  'export const sameFeatures = (rows, group, featureRow) => {',
  "  const key = group.getAttribute('data-legend-key');",
  '  return rows.some((row) => drawsFeatures(row, key) === featureRow);',
  '};',
  // L2: write-if-different is change detection.
  'export const syncStroke = (paths, color) => {',
  '  for (const path of paths) {',
  "    if (path.getAttribute('stroke') !== String(color)) path.setAttribute('stroke', String(color));",
  '  }',
  '};',
  'export const syncText = (groups, entries) => groups.forEach((group, index) => {',
  "  const text = group.querySelector('text');",
  "  if (text && String(text.textContent || '') !== entries[index].caption) { text.textContent = entries[index].caption; }",
  '});',
  // L3: an identity across `?`/`:` does not qualify the other branch.
  'export const ternary = (entries, k, c) => entries.find((e) => (k && e.id === k ? true : e.color === c));',
  // L4: a shown local in a one-line body.
  "export const oneLine = (es, g) => { const caption = g.getAttribute('data-legend-key'); return es.find((e) => e.caption === caption); };",
  ''
].join('\n');

test('identity from display keeps optional chains, scopes shown locals, and exempts write-if-different', () => {
  const result = WEB_OWNER_GRAPH_DETECTORS[IDENTITY_FROM_DISPLAY].detect(new Map([
    ['gbdraw/web/js/app/legend/sort-actions.js', IDENTITY_REVIEW_SOURCE]
  ]));
  assert.deepEqual(result.observedSites.map((record) => `${record.function} ${record.class} ${record.form}:${record.field}`), [
    'hasShownCaption shown-join search:textContent',
    'optionalIdentity id-qualified search:color',
    'siblings shown-join search:caption',
    'siblings caption-key search:caption',
    'outerName caption-key search:caption',
    'shadowed shown-join search:caption',
    'shadowed caption-key key-has:caption',
    'ownIdentity shown-join search:@data-legend-key',
    'overrideOfKey shown-join key-index:legendKey',
    'rankOf shown-join key-get:keyOf',
    'syncStroke write-if-different loop-join:@stroke',
    'syncText write-if-different loop-join:textContent',
    'ternary paint-join search:color',
    'oneLine shown-join search:caption'
  ]);
});
