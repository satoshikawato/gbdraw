// Web function size (docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md;
// WEB-FUNCTION-SIZE in tools/web-design-rule-guards.json). Every function of
// the modules under gbdraw/web/js/, except the generated `*.generated.js`
// ones, has at most MAX_FUNCTION_LINES lines. FUNCTION_LINES lists the longer
// functions the rule started with, at their length: a listed function may
// shrink but not grow, and any other function, a new one included, stays
// within the limit. Lowering or removing an entry is optional; the test names
// the entries that may be lowered. MAX_FUNCTION_LINES and FUNCTION_LINES are
// registered design-rule allowlists: raising the limit, and adding or raising
// an entry, are authority-only changes that need the Owner's approval.
//
// What counts. A function is each function declaration, function expression,
// arrow function, method, accessor and constructor, nested ones included. Its
// lines run from its first token (`export`, `async`, `function`, its name or
// its parameters; a JSDoc block above it does not count) to its last, so a
// nested function counts in the functions around it as well.
//
// Names. A function is `<path under gbdraw/web/js/>::<qualified name>`, the
// qualified name being the names of the functions, classes and named object
// literals around it and its own, joined by `.`. A function expression or an
// object literal takes the name it is bound to (a variable, a property, a
// class field, an assignment target, a parameter or destructuring default),
// through the calls and expressions around it: in `const x = computed(() =>
// ...)` the arrow is `x`. A callback bound to no name is named after the call
// it is passed to (`forEach()`, `new Promise()`), and any other function is
// `(anonymous)`. Functions with one name in one module count as the longest
// of them.
//
// Logic in composition roots. A composition root of the Web owner graph
// (WEB_OWNER_GRAPH_DEFAULTS.compositionRoots in
// tools/web-owner-graph-detectors.mjs) should hold wiring: it creates the
// owners, passes them to each other as ports and returns what the app uses.
// Its root function is the function at the top level of the module that calls
// a `create*` factory imported by name from another module (an owner
// creation); each composition root has exactly one. The logic written in the
// root is every function whose nearest enclosing function is the root
// function: a `const f = () => {...}` statement, a function passed inline to a
// factory as a port, a method of an object literal built in the root, a
// computed getter, a watch callback. Each counts with its whole span, as
// above, so the functions nested in it count once, in it.
// COMPOSITION_ROOT_LOGIC_LINES records the total lines of that logic per
// `<path>::<root function>`. A total may shrink but not grow, and every
// composition root has an entry. When it fails, move the new logic into a
// module factory and pass it to the root as a port, or shrink other logic in
// the same root. Adding or raising an entry (a new composition root, a renamed
// root function, more logic) needs the Owner's approval; an entry for a module
// that is no longer a composition root is reported as removable.
import assert from 'node:assert/strict';
import { mkdtempSync, readdirSync, rmSync, writeFileSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { dirname, join, sep } from 'node:path';
import test from 'node:test';
import { fileURLToPath } from 'node:url';

import { WEB_OWNER_GRAPH_DEFAULTS } from '../../tools/web-owner-graph-detectors.mjs';
import { compareCountBaseline, reportLowerable, typescriptApi } from './helpers/count-baseline.mjs';

const MAX_FUNCTION_LINES = 150;

const FUNCTION_LINES = {
  'app/annotations.js::createAnnotationEditor': 332,
  'app/app-setup.js::createAppSetup': 3853,
  'app/candidate-render.js::compilePlanBundle': 409,
  'app/circular-track-slots.js::createCircularTrackSlotEditor': 1456,
  'app/drawing-list.js::createDrawingList': 220,
  'app/feature-editor.js::createFeatureEditor': 252,
  'app/feature-editor/color-actions.js::createFeatureColorActions': 1622,
  'app/feature-editor/color-actions.js::createFeatureColorActions.continueLegendRenameRequest': 167,
  'app/feature-editor/label-actions.js::createFeatureLabelActions': 1243,
  'app/feature-editor/placement-actions.js::createFeaturePlacementActions': 157,
  'app/feature-editor/rule-actions.js::createFeatureRuleActions': 1150,
  'app/feature-editor/svg-actions.js::createFeatureSvgActions': 1921,
  'app/feature-editor/svg-actions.js::createFeatureSvgActions.attachSvgFeatureHandlers': 1100,
  'app/feature-editor/svg-actions.js::createFeatureSvgActions.attachSvgFeatureHandlers.showAlignmentOverlay': 206,
  'app/feature-editor/visibility-actions.js::createFeatureVisibilityActions': 715,
  'app/feature-search/preview-actions.js::createPreviewFeatureSearch': 476,
  'app/feature-selection.js::createFeatureSelection': 513,
  'app/font-settings.js::createFontSettings': 324,
  'app/history-inputs.js::setupHistoryInputs': 181,
  'app/label-fit-dialog.js::createLabelFitDialog': 164,
  'app/legend-layout/diagram-drag.js::createDiagramDragActions': 546,
  'app/legend/drag-actions.js::createLegendDragActions': 223,
  'app/legend/entry-actions.js::createLegendEntryActions': 889,
  'app/linear-record-selector.js::createLinearRecordSelector': 320,
  'app/linear-track-slots.js::createLinearTrackSlotEditor': 1037,
  'app/orthogroups.js::createOrthogroupEditor': 671,
  'app/pairwise-match-popup.js::buildMatchPopupPayload': 217,
  'app/pairwise-match-popup.js::createPairwisePayloadContext': 216,
  'app/preview-item-popup.js::createPreviewItemPopup': 214,
  'app/preview-runtime.js::createPreviewRuntime': 947,
  'app/preview-runtime.js::createPreviewRuntime.bindMountedResult': 159,
  'app/record-display-options.js::createRecordDisplayControls': 514,
  'app/record-display/feature-record-rotation.js::createFeatureRecordRotationWorkflow': 279,
  'app/record-selection.js::createRecordSelection': 183,
  'app/region-drawing.js::createRegionDrawing': 160,
  'app/region-select.js::createRegionSelect': 396,
  'app/results.js::createResultsManager': 285,
  'app/rule-matching.js::createRulePreparation': 224,
  'app/run-analysis.js::createRunAnalysis': 4867,
  'app/run-analysis.js::createRunAnalysis.runAnalysis': 179,
  'app/run-analysis.js::createRunAnalysis.runAnalysisInternal': 2810,
  'app/run-analysis.js::createRunAnalysis.runAnalysisInternal.getSeqEntry': 180,
  'app/run-analysis.js::createRunAnalysis.runAnalysisInternal.runCircularLosatConservation': 157,
  'app/run-analysis.js::createRunAnalysis.runCommittedCanonicalCandidateInternal': 233,
  'app/run-info.js::appendConfigOverrides': 160,
  'app/run-info.js::buildRunInfo': 201,
  'app/run-info.js::validateCurrentSemanticCoverage': 165,
  'app/settings-pane.js::createSettingsPaneController': 178,
  'app/similarity-alignment.js::createSimilarityAlignmentActions': 1112,
  'app/svg-styles.js::createSvgStyles': 186,
  'app/ui.js::createPanZoom': 414,
  'app/watchers.js::setupWatchers': 449,
  'components.js::SettingsSectionIndex.setup': 163,
  'services/config.js::exportSessionDocument': 193,
  'services/config.js::importSessionDocument': 1418,
  'services/diagram-generation.js::runDiagramGeneration': 169,
  'services/diagram-resource-staging.js::createDiagramResourceTransport': 155,
  'services/feature-catalog.js::validateAndProjectCatalogItem': 372,
  'services/feature-override-identity.js::migrateLegacyFeatureOverrides': 152,
  'services/history-snapshot.js::createHistorySnapshotService': 1116,
  'services/history.js::createHistoryManager': 1623,
  'services/history.js::createHistoryRouter': 653,
  'services/legacy-similarity-alignment.js::materializeLegacySimilarityAlignment': 183,
  'services/linear-comparisons.js::buildLinearComparisonTimeline': 183,
  'services/losat.js::runLosatPairsThreaded': 220,
  'services/losat.js::runLosatPairsThreaded.new Promise()': 179,
  'services/losat.js::runLosatPairsWithWorkers': 174,
  'services/match-sequences.js::analyzeCatalogSequenceSourceCoverage': 319,
  'services/match-sequences.js::createSequenceSourceRegistry': 167,
  'services/orthogroup-feature-metadata.js::createOrthogroupFeatureProjection': 214,
  'services/session-feature-recovery.js::alignRecoveredFeatureIdsToRenderedSvg': 181,
  'services/session-feature-recovery.js::buildSessionFeatureRecoveryPlan': 188,
  'services/session-request.js::buildComparisons': 273,
  'services/session-request.js::buildConfigOverrides': 187,
  'services/session-request.js::buildRecords': 151,
  'services/session-request.js::projectCanonicalRenderInput': 399,
  'services/session-request.js::projectCanonicalSessionRequest': 1106,
  'services/session-request.js::trackPlanProjection': 201,
  'services/session-resources.js::buildSessionDrawingResources': 376,
  'services/track-slot-validation.js::collectCustomTrackIssues': 589,
  'services/track-slot-validation.js::collectCustomTrackIssues.forEach()': 351,
  'services/track-slot-validation.js::validateTrackSlotBindingInvariants': 175,
  'workers/diagram-generation-worker.js::runGeneration': 231,
  'workers/losat-threaded-worker.js::runThreadedLosat': 191
};

const COMPOSITION_ROOT_LOGIC_LINES = {
  'app/app-setup.js::createAppSetup': 1256,
  'app/feature-editor.js::createFeatureEditor': 9,
  'app/legend-layout.js::createLegendLayout': 75,
  'app/legend.js::createLegendManager': 21
};

const REPOSITORY_ROOT = join(dirname(fileURLToPath(import.meta.url)), '..', '..');
const WEB_JS = join(REPOSITORY_ROOT, 'gbdraw', 'web', 'js');
const TSCONFIG = join(REPOSITORY_ROOT, 'tests', 'web', 'types', 'tsconfig.json');

// The functions of one parsed source file: [{ name, line, lines, parent }],
// `line` being the first line of the function and `parent` the nearest
// function around it (null at the top level of the module).
const functionsOf = (source, SyntaxKind) => {
  const FUNCTIONS = new Set([
    SyntaxKind.FunctionDeclaration, SyntaxKind.FunctionExpression, SyntaxKind.ArrowFunction,
    SyntaxKind.MethodDeclaration, SyntaxKind.GetAccessor, SyntaxKind.SetAccessor, SyntaxKind.Constructor
  ]);
  const NAMED_SCOPES = new Set([
    ...FUNCTIONS, SyntaxKind.ClassDeclaration, SyntaxKind.ClassExpression, SyntaxKind.ObjectLiteralExpression
  ]);
  // Expressions a function or object literal passes through to the name it is
  // bound to; a call or `new` passes it on only as an argument.
  const PASS_THROUGH = new Set([
    SyntaxKind.ParenthesizedExpression, SyntaxKind.AsExpression, SyntaxKind.SatisfiesExpression,
    SyntaxKind.NonNullExpression, SyntaxKind.AwaitExpression, SyntaxKind.TypeAssertionExpression,
    SyntaxKind.ConditionalExpression, SyntaxKind.BinaryExpression
  ]);
  const CALLS = new Set([SyntaxKind.CallExpression, SyntaxKind.NewExpression]);
  const BINDINGS = new Set([
    SyntaxKind.VariableDeclaration, SyntaxKind.PropertyAssignment, SyntaxKind.PropertyDeclaration,
    SyntaxKind.Parameter, SyntaxKind.BindingElement
  ]);
  const lineOf = (position) => source.getLineAndCharacterOfPosition(position).line + 1;
  const textOf = (node) => source.text.slice(node.getStart(source), node.end).replace(/\s+/g, ' ');
  const nameText = (name) => ([SyntaxKind.Identifier, SyntaxKind.PrivateIdentifier, SyntaxKind.StringLiteral,
    SyntaxKind.NumericLiteral].includes(name.kind) ? name.text : textOf(name));
  const callName = (call) => {
    const callee = call.expression;
    const name = callee.kind === SyntaxKind.Identifier ? callee.text
      : callee.kind === SyntaxKind.PropertyAccessExpression ? callee.name.text : null;
    return name === null ? null : `${call.kind === SyntaxKind.NewExpression ? 'new ' : ''}${name}()`;
  };
  // null: an object literal bound to no name, which adds nothing to the names
  // inside it.
  const nameOf = (node) => {
    if (node.kind === SyntaxKind.Constructor) return 'constructor';
    if (node.kind === SyntaxKind.FunctionDeclaration || node.kind === SyntaxKind.ClassDeclaration) {
      return node.name ? node.name.text : 'default';
    }
    if ([SyntaxKind.MethodDeclaration, SyntaxKind.GetAccessor, SyntaxKind.SetAccessor].includes(node.kind)) {
      return nameText(node.name);
    }
    let call = null;
    for (let child = node, up = node.parent; ; child = up, up = up.parent) {
      if (up.kind === SyntaxKind.BinaryExpression && up.operatorToken.kind === SyntaxKind.EqualsToken) {
        if (up.right === child) return textOf(up.left);
        break;
      }
      if (PASS_THROUGH.has(up.kind)) continue;
      if (CALLS.has(up.kind) && up.expression !== child) {
        call ??= callName(up);
        continue;
      }
      if (BINDINGS.has(up.kind) && up.initializer === child && up.name.kind !== SyntaxKind.ObjectBindingPattern
        && up.name.kind !== SyntaxKind.ArrayBindingPattern) {
        return nameText(up.name);
      }
      if (up.kind === SyntaxKind.ExportAssignment) return 'default';
      break;
    }
    if (node.name) return node.name.text;
    return call ?? (node.kind === SyntaxKind.ObjectLiteralExpression ? null : '(anonymous)');
  };
  const found = [];
  const visit = (node, scope, enclosing) => {
    let inner = scope;
    let innerEnclosing = enclosing;
    if (NAMED_SCOPES.has(node.kind)) {
      const name = nameOf(node);
      if (name !== null) inner = [...scope, name];
      if (FUNCTIONS.has(node.kind)) {
        const line = lineOf(node.getStart(source));
        innerEnclosing = { name: inner.join('.'), line, lines: lineOf(node.end) - line + 1, parent: enclosing };
        found.push(innerEnclosing);
      }
    }
    node.forEachChild((child) => visit(child, inner, innerEnclosing));
  };
  visit(source, [], null);
  return found;
};

// Parses `files` (absolute paths) with the pinned TypeScript; the modules of
// the R14 config in its project, any other file in an inferred one.
// Map(file -> read(source, SyntaxKind)), `read` being functionsOf by default.
const measure = async (files, read = functionsOf) => {
  const { API, ast: { SyntaxKind } } = await typescriptApi();
  const api = new API({ cwd: REPOSITORY_ROOT });
  try {
    const snapshot = api.updateSnapshot({ openProjects: [TSCONFIG], openFiles: files });
    return new Map(files.map((file) => {
      const source = snapshot.getDefaultProjectForFile(file)?.program.getSourceFile(file);
      assert.ok(source, `the TypeScript API returned no source file for ${file}`);
      return [file, read(source, SyntaxKind)];
    }));
  } finally {
    api.close();
  }
};

// The root functions of one parsed composition root module and the logic
// written in it: { roots, logic }, both as functionsOf entries. A root function
// is a top-level function that calls a `create*` factory the module imports by
// name; `logic` lists the functions whose parent is the only root function.
const compositionRootOf = (source, SyntaxKind) => {
  const functions = functionsOf(source, SyntaxKind);
  const factories = new Set(source.statements
    .filter((statement) => statement.kind === SyntaxKind.ImportDeclaration && statement.importClause
      && statement.importClause.namedBindings && statement.importClause.namedBindings.kind === SyntaxKind.NamedImports)
    .flatMap((statement) => statement.importClause.namedBindings.elements.map((element) => element.name.text))
    .filter((name) => /^create[A-Z]/.test(name)));
  const creationLines = [];
  const visit = (node) => {
    if (node.kind === SyntaxKind.CallExpression && node.expression.kind === SyntaxKind.Identifier
      && factories.has(node.expression.text)) {
      creationLines.push(source.getLineAndCharacterOfPosition(node.getStart(source)).line + 1);
    }
    node.forEachChild(visit);
  };
  visit(source);
  const roots = functions.filter(({ parent, line, lines }) => parent === null
    && creationLines.some((creation) => creation >= line && creation < line + lines));
  return { roots, logic: roots.length === 1 ? functions.filter(({ parent }) => parent === roots[0]) : [] };
};

test('a function is measured from its first token to its last and named by its binding', async () => {
  const directory = mkdtempSync(join(tmpdir(), 'web-function-size-'));
  const file = join(directory, 'fixture.js');
  writeFileSync(file, [
    '/** A JSDoc block does not count. */',
    'export async function outer(items) {',
    '  const inner = (value) => {',
    '    return value;',
    '  };',
    '  items.forEach((item) => {',
    '    inner(item);',
    '  });',
    '  const ready = await new Promise((resolve) => {',
    '    resolve(true);',
    '  });',
    '  const doubled = items.map((item) => item * 2);',
    '  const actions = { run() { return 1; }, stop: () => 2 };',
    '  return ready && doubled && actions;',
    '}',
    'class Store {',
    '  constructor() {}',
    '  get size() { return 0; }',
    '  handler = function named() {};',
    '}',
    'window.onload = () => {};',
    'const pick = (choose = () => 1) => choose;',
    '[1].map(function () {});',
    '(() => {})();',
    'export default () => {};',
    ''
  ].join('\n'));
  try {
    const functions = (await measure([file])).get(file);
    assert.deepEqual(functions.map(({ name, line, lines }) => `${line}+${lines} ${name}`), [
      '2+14 outer',
      '3+3 outer.inner',
      '6+3 outer.forEach()',
      '9+3 outer.ready',
      '12+1 outer.doubled',
      '13+1 outer.actions.run',
      '13+1 outer.actions.stop',
      '17+1 Store.constructor',
      '18+1 Store.size',
      '19+1 Store.handler',
      '21+1 window.onload',
      '22+1 pick',
      '22+1 pick.choose',
      '23+1 map()',
      '24+1 (anonymous)',
      '25+1 default'
    ]);
  } finally {
    rmSync(directory, { recursive: true, force: true });
  }
});

test('no Web function is longer than MAX_FUNCTION_LINES or its FUNCTION_LINES entry', async (t) => {
  assert.deepEqual(Object.keys(FUNCTION_LINES), Object.keys(FUNCTION_LINES).sort(), 'keep FUNCTION_LINES sorted by key');
  Object.entries(FUNCTION_LINES).forEach(([key, lines]) => {
    assert.ok(Number.isInteger(lines) && lines > MAX_FUNCTION_LINES,
      `FUNCTION_LINES['${key}'] must be an integer above MAX_FUNCTION_LINES (${MAX_FUNCTION_LINES})`);
  });
  const paths = readdirSync(WEB_JS, { recursive: true })
    .map((path) => path.split(sep).join('/'))
    .filter((path) => path.endsWith('.js') && !path.endsWith('.generated.js'))
    .sort();
  const measured = await measure(paths.map((path) => join(WEB_JS, path)));
  // The longest function of each key.
  const longest = new Map();
  paths.forEach((path) => measured.get(join(WEB_JS, path)).forEach((found) => {
    const key = `${path}::${found.name}`;
    if (!longest.has(key) || found.lines > longest.get(key).lines) longest.set(key, { path, ...found });
  }));
  // A key without an entry counts as 0, so every function above the limit is
  // compared, and every listed one, to name the entries that may be lowered.
  const observed = Object.fromEntries([...longest]
    .filter(([key, { lines }]) => lines > MAX_FUNCTION_LINES || Object.hasOwn(FUNCTION_LINES, key))
    .map(([key, { lines }]) => [key, lines]));
  const { above, lowerable } = compareCountBaseline(observed, FUNCTION_LINES);
  const failures = above.map(({ key, observed: lines, recorded }) => {
    const { path, line, name } = longest.get(key);
    return `gbdraw/web/js/${path}:${line} ${name}: ${lines} lines (limit ${MAX_FUNCTION_LINES}, `
      + `FUNCTION_LINES ${recorded || 'none'})`;
  });
  t.diagnostic(`functions: ${[...measured.values()].reduce((sum, found) => sum + found.length, 0)} in ${paths.length} modules, `
    + `${Object.keys(observed).length} above ${MAX_FUNCTION_LINES} lines`);
  reportLowerable(t, 'FUNCTION_LINES', lowerable.map((entry) => entry.replace(/ -> 0$/, ' -> removed')));
  assert.deepEqual(failures, [], `functions above ${MAX_FUNCTION_LINES} lines, or above their FUNCTION_LINES entry. `
    + 'Split them at an owner boundary (gbdraw/web/CLAUDE.md "Module ownership"); adding or raising an entry '
    + `in FUNCTION_LINES needs the Owner's approval.\n  ${failures.join('\n  ')}`);
});

test('the logic of a composition root is every function written directly in its root function', async () => {
  const directory = mkdtempSync(join(tmpdir(), 'web-root-logic-'));
  const file = join(directory, 'fixture.js');
  writeFileSync(file, [
    "import { createOwner } from './owner.js';",
    "import * as helpers from './helpers.js';",
    'const createLocal = () => ({});',
    'const local = () => createLocal() && helpers.createOther();',
    'export const createState = () => ({ value: 0 });',
    'export const createRoot = ({ watch }) => {',
    '  const owner = createOwner({',
    '    onPick: (item) => {',
    '      return item;',
    '    }',
    '  });',
    '  const run = () => {',
    '    const nested = () => 1;',
    '    return nested();',
    '  };',
    '  watch(() => owner.value, () => {',
    '    run();',
    '  });',
    '  const steps = { first() { return 1; } };',
    '  return { owner, run, steps, get size() { return 0; } };',
    '};',
    ''
  ].join('\n'));
  try {
    const { roots, logic } = (await measure([file], compositionRootOf)).get(file);
    assert.deepEqual(roots.map(({ name, line, lines }) => `${line}+${lines} ${name}`), ['6+16 createRoot']);
    assert.deepEqual(logic.map(({ name, line, lines }) => `${line}+${lines} ${name}`), [
      '8+3 createRoot.owner.onPick',
      '12+4 createRoot.run',
      '16+1 createRoot.watch()',
      '16+3 createRoot.watch()',
      '19+1 createRoot.steps.first',
      '20+1 createRoot.size'
    ]);
  } finally {
    rmSync(directory, { recursive: true, force: true });
  }
});

test('the logic written in each composition root does not grow beyond COMPOSITION_ROOT_LOGIC_LINES', async (t) => {
  const recorded = COMPOSITION_ROOT_LOGIC_LINES;
  assert.deepEqual(Object.keys(recorded), Object.keys(recorded).sort(), 'keep COMPOSITION_ROOT_LOGIC_LINES sorted by key');
  Object.entries(recorded).forEach(([key, lines]) => {
    assert.ok(Number.isInteger(lines) && lines >= 0, `COMPOSITION_ROOT_LOGIC_LINES['${key}'] must be an integer of 0 or more`);
  });
  const paths = [...WEB_OWNER_GRAPH_DEFAULTS.compositionRoots].sort();
  const measured = await measure(paths.map((path) => join(WEB_JS, path)), compositionRootOf);
  const observed = {};
  const rows = new Map();
  const rootProblems = paths.flatMap((path) => {
    const { roots, logic } = measured.get(join(WEB_JS, path));
    if (roots.length !== 1) {
      return [`gbdraw/web/js/${path}: ${roots.length} top-level functions create owners`
        + `${roots.length ? ` (${roots.map(({ name, line }) => `${name} at line ${line}`).join(', ')})` : ''}; `
        + 'a composition root creates its owners in one root function'];
    }
    const key = `${path}::${roots[0].name}`;
    observed[key] = logic.reduce((sum, { lines }) => sum + lines, 0);
    rows.set(key, { path, root: roots[0], functions: logic.length });
    return [];
  });
  assert.deepEqual(rootProblems, [], rootProblems.join('\n'));
  rows.forEach(({ root, functions }, key) => t.diagnostic(
    `${key}: ${observed[key]} lines in ${functions} functions (the root function has ${root.lines} lines)`
  ));
  const stale = Object.keys(recorded).filter((key) => !Object.hasOwn(observed, key));
  const { above, lowerable } = compareCountBaseline(observed, recorded);
  reportLowerable(t, 'COMPOSITION_ROOT_LOGIC_LINES', lowerable.map((entry) => (
    stale.some((key) => entry.startsWith(`${key} `)) ? entry.replace(/ -> 0$/, ' -> removed (no longer a composition root)') : entry
  )));
  const failures = [
    ...Object.keys(observed).filter((key) => !Object.hasOwn(recorded, key)).map((key) => {
      const { path, root } = rows.get(key);
      return `gbdraw/web/js/${path}:${root.line} ${root.name}: ${observed[key]} lines of logic in a composition root `
        + 'without a COMPOSITION_ROOT_LOGIC_LINES entry (a new composition root or a renamed root function)';
    }),
    ...above.filter(({ key }) => Object.hasOwn(recorded, key)).map(({ key, observed: lines, recorded: entry }) => {
      const { path, root, functions } = rows.get(key);
      return `gbdraw/web/js/${path}:${root.line} ${root.name}: ${lines} lines of logic in ${functions} functions `
        + `(COMPOSITION_ROOT_LOGIC_LINES ${entry})`;
    })
  ];
  assert.deepEqual(failures, [], 'the logic written in a composition root grew. Move the new logic into a module '
    + 'factory and pass it to the root as a port, or shrink other logic in the same root '
    + '(docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md "Web function size"); adding or raising an entry in '
    + `COMPOSITION_ROOT_LOGIC_LINES needs the Owner's approval.\n  ${failures.join('\n  ')}`);
});
