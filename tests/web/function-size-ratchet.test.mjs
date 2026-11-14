// Web function size (docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md;
// WEB-FUNCTION-SIZE in tools/web-design-rule-guards.json). Every function of
// the modules under gbdraw/web/js/, except the generated `*.generated.js`
// ones, has at most MAX_FUNCTION_STATEMENTS statements. FUNCTION_STATEMENTS
// lists the larger functions at their size: a listed function may shrink but
// not grow, and any other function, a new one included, stays within the
// limit. Lowering or removing an entry is optional; the test names the entries
// that may be lowered. MAX_FUNCTION_STATEMENTS and FUNCTION_STATEMENTS are
// registered design-rule allowlists: raising the limit, and adding or raising
// an entry, are authority-only changes that need the Owner's approval.
//
// What counts. A function is each function declaration, function expression,
// arrow function, method, accessor and constructor, nested ones included. Its
// statements are the statements in its body, those of the functions and
// classes nested in it included, so a nested function counts in the functions
// around it as well. A block (`{ ... }`) and an empty statement are not
// statements; an arrow function's expression body (`=> value`) counts as one,
// its implicit `return`. Comments, JSDoc, blank lines and line breaks change
// no count.
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
// root is every function of two or more statements whose nearest enclosing
// function is the root function: a `const f = () => {...}` statement, a
// function passed inline to a factory as a port, a method of an object literal
// built in the root, a computed getter, a watch callback. Each counts with all
// its statements, as above, so the functions nested in it count once, in it. A
// function of one statement written directly in the root function, such as
// the port `drawingTabs: () => state.drawingTabs.value`, is wiring and does
// not count; a one-statement function inside a counted one counts with it.
// Statements written directly in the root body and the other top-level
// functions of the module are not counted here; the root's whole size stays
// bounded by MAX_FUNCTION_STATEMENTS and FUNCTION_STATEMENTS.
// COMPOSITION_ROOT_LOGIC_STATEMENTS records the total statements of that logic
// per `<path>::<root function>`. A total may shrink but not grow, and every
// composition root has an entry. When it fails, move the new logic into a
// module factory and pass it to the root as a port, or shrink other logic in
// the same root. Adding or raising an entry (a new composition root, a renamed
// root function, more logic) needs the Owner's approval. An entry whose key
// names no root function now (the module is no longer a composition root, or
// its root function took another name) is reported as removable, whatever its
// value.
import assert from 'node:assert/strict';
import { mkdtempSync, readdirSync, rmSync, writeFileSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { dirname, join, sep } from 'node:path';
import test from 'node:test';
import { fileURLToPath } from 'node:url';

import { WEB_OWNER_GRAPH_DEFAULTS } from '../../tools/web-owner-graph-detectors.mjs';
import { compareCountBaseline, reportLowerable, typescriptApi } from './helpers/count-baseline.mjs';

// 150 lines at the median density of the functions of 100 to 200 lines when
// the limit moved from lines to statements (1.76 lines a statement, 2026-11-10).
const MAX_FUNCTION_STATEMENTS = 85;

const FUNCTION_STATEMENTS = {
  'app/align-at-feature.js::createFeatureAnchorAction': 94,
  'app/annotations.js::createAnnotationEditor': 218,
  'app/app-setup.js::createAppSetup': 1098,
  'app/candidate-render.js::compilePlanBundle': 227,
  'app/circular-record-controls.js::createCircularSourceEdits': 86,
  'app/circular-track-slots.js::createCircularTrackSlotEditor': 857,
  'app/drawing-list.js::createDrawingList': 107,
  'app/feature-editor/color-actions.js::createFeatureColorActions': 924,
  'app/feature-editor/label-actions.js::createFeatureLabelActions': 634,
  'app/feature-editor/placement-actions.js::createFeaturePlacementActions': 100,
  'app/feature-editor/rule-actions.js::createFeatureRuleActions': 646,
  'app/feature-editor/svg-actions.js::createFeatureSvgActions': 1176,
  'app/feature-editor/svg-actions.js::createFeatureSvgActions.attachSvgFeatureHandlers': 773,
  'app/feature-editor/svg-actions.js::createFeatureSvgActions.attachSvgFeatureHandlers.showAlignmentOverlay': 147,
  'app/feature-editor/visibility-actions.js::createFeatureVisibilityActions': 359,
  'app/feature-search/preview-actions.js::createPreviewFeatureSearch': 237,
  'app/feature-selection.js::createFeatureSelection': 321,
  'app/font-settings.js::createFontSettings': 185,
  'app/history-inputs.js::setupHistoryInputs': 127,
  'app/history-list.js::createHistoryControls': 93,
  'app/label-fit-dialog.js::createLabelFitDialog': 99,
  'app/legend-layout/diagram-drag.js::createDiagramDragActions': 338,
  'app/legend/drag-actions.js::createLegendDragActions': 154,
  'app/legend/entry-actions.js::createLegendEntryActions': 485,
  'app/legend/entry-actions.js::createLegendEntryActions.prepareFileLegendEntries': 98,
  'app/linear-comparison-edits.js::createLinearComparisonGapEdits': 98,
  'app/linear-comparison-edits.js::createLinearComparisonPairEdits': 102,
  'app/linear-record-selector.js::createLinearRecordSelector': 154,
  'app/linear-track-slots.js::createLinearTrackSlotEditor': 649,
  'app/orthogroups.js::createOrthogroupEditor': 275,
  'app/preview-runtime.js::createPreviewRuntime': 381,
  'app/record-display-options.js::createRecordDisplayControls': 307,
  'app/record-display/feature-record-rotation.js::createFeatureRecordRotationWorkflow': 165,
  'app/record-selection.js::createRecordSelection': 115,
  'app/region-drawing.js::createRegionDrawing': 94,
  'app/region-select.js::createRegionSelect': 240,
  'app/results.js::createResultsManager': 142,
  'app/rule-matching.js::createRulePreparation': 112,
  'app/run-analysis.js::createRunAnalysis': 2138,
  'app/run-analysis.js::createRunAnalysis.runAnalysis': 111,
  'app/run-analysis.js::createRunAnalysis.runAnalysisInternal': 1175,
  'app/run-analysis.js::createRunAnalysis.runCommittedCanonicalCandidateInternal': 110,
  'app/run-info.js::appendConfigOverrides': 100,
  'app/run-info.js::appendInputArgs': 92,
  'app/sequence-search.js::createSequenceSearch': 89,
  'app/settings-pane.js::createSettingsPaneController': 123,
  'app/similarity-alignment.js::createSimilarityAlignmentActions': 587,
  'app/svg-styles.js::createSvgStyles': 97,
  'app/svg-styles.js::createSvgStyles.applyPaletteToSvg': 87,
  'app/ui.js::createPreviewTransformInteraction': 88,
  'app/watchers.js::setupWatchers': 224,
  'components.js::FileUploader.setup': 120,
  'components.js::SettingsSectionIndex.setup': 115,
  'services/config.js::importSessionDocument': 566,
  'services/diagram-generation.js::runDiagramGeneration': 96,
  'services/diagram-resource-staging.js::createDiagramResourceTransport': 87,
  'services/feature-catalog.js::validateAndProjectCatalogItem': 191,
  'services/feature-override-identity.js::migrateLegacyFeatureOverrides': 94,
  'services/history-snapshot.js::createHistorySnapshotService': 390,
  'services/history.js::createHistoryManager': 897,
  'services/history.js::createHistoryRouter': 449,
  'services/legend-layout.js::buildCircularLegendLayout': 97,
  'services/losat.js::runLosatPairsThreaded': 136,
  'services/losat.js::runLosatPairsThreaded.new Promise()': 116,
  'services/losat.js::runLosatPairsWithWorkers': 115,
  'services/orthogroup-feature-metadata.js::createOrthogroupFeatureProjection': 113,
  'services/session-feature-recovery.js::alignRecoveredFeatureIdsToRenderedSvg': 98,
  'services/session-request.js::buildComparisons': 113,
  'services/session-request.js::projectCanonicalRenderInput': 119,
  'services/session-request.js::projectCanonicalSessionRequest': 332,
  'services/session-request.js::trackPlanProjection': 145,
  'services/session-resources.js::buildSessionDrawingResources': 205,
  'services/svg-result-ingestion.js::applyLegendOperations': 96,
  'services/svg-result-ingestion.js::recordSavedEditBases': 94,
  'services/track-slot-validation.js::collectCustomTrackIssues': 192,
  'services/track-slot-validation.js::collectCustomTrackIssues.forEach()': 104,
  'services/track-slot-validation.js::validateTrackSlotBindingInvariants': 98,
  'utils/contract-repair.js::repairPlace': 121,
  'workers/diagram-generation-worker.js::runGeneration': 96,
  'workers/losat-threaded-worker.js::runThreadedLosat': 104
};

const COMPOSITION_ROOT_LOGIC_STATEMENTS = {
  'app/app-setup.js::createAppSetup': 625,
  'app/feature-editor.js::createFeatureEditor': 0,
  'app/legend-layout.js::createLegendLayout': 56,
  'app/legend.js::createLegendManager': 21
};

const REPOSITORY_ROOT = join(dirname(fileURLToPath(import.meta.url)), '..', '..');
const WEB_JS = join(REPOSITORY_ROOT, 'gbdraw', 'web', 'js');
const TSCONFIG = join(REPOSITORY_ROOT, 'tests', 'web', 'types', 'tsconfig.json');

// The functions of one parsed source file: [{ name, line, start, end,
// statements, parent }], `line` being the first line of the function, `start`
// and `end` its first and last positions in the text, and `parent` the nearest
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
  const STATEMENTS = new Set([
    SyntaxKind.VariableStatement, SyntaxKind.ExpressionStatement, SyntaxKind.IfStatement, SyntaxKind.DoStatement,
    SyntaxKind.WhileStatement, SyntaxKind.ForStatement, SyntaxKind.ForInStatement, SyntaxKind.ForOfStatement,
    SyntaxKind.ContinueStatement, SyntaxKind.BreakStatement, SyntaxKind.ReturnStatement, SyntaxKind.WithStatement,
    SyntaxKind.SwitchStatement, SyntaxKind.LabeledStatement, SyntaxKind.ThrowStatement, SyntaxKind.TryStatement,
    SyntaxKind.DebuggerStatement, SyntaxKind.FunctionDeclaration, SyntaxKind.ClassDeclaration
  ]);
  const expressionBody = (node) => node.kind === SyntaxKind.ArrowFunction && node.body.kind !== SyntaxKind.Block;
  // The statements of a function body, nested ones included; an expression
  // body, the function's own or a nested arrow function's, is one.
  const statementsOf = (body) => {
    let count = body.kind === SyntaxKind.Block ? 0 : 1;
    const walk = (node) => {
      if (STATEMENTS.has(node.kind) || expressionBody(node)) count += 1;
      node.forEachChild(walk);
    };
    walk(body);
    return count;
  };
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
        const start = node.getStart(source);
        innerEnclosing = {
          name: inner.join('.'), line: lineOf(start), start, end: node.end,
          statements: statementsOf(node.body), parent: enclosing
        };
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
// written in it: { roots, logic, wiring }, `roots` and `logic` as functionsOf
// entries. A root function is a top-level function that calls a `create*`
// factory the module imports by name. `logic` lists the functions of two or
// more statements whose parent is the only root function, and `wiring` counts
// its functions of one statement or none.
const compositionRootOf = (source, SyntaxKind) => {
  const functions = functionsOf(source, SyntaxKind);
  const factories = new Set(source.statements
    .filter((statement) => statement.kind === SyntaxKind.ImportDeclaration && statement.importClause
      && statement.importClause.namedBindings && statement.importClause.namedBindings.kind === SyntaxKind.NamedImports)
    .flatMap((statement) => statement.importClause.namedBindings.elements.map((element) => element.name.text))
    .filter((name) => /^create[A-Z]/.test(name)));
  const creations = [];
  const visit = (node) => {
    if (node.kind === SyntaxKind.CallExpression && node.expression.kind === SyntaxKind.Identifier
      && factories.has(node.expression.text)) {
      creations.push(node.getStart(source));
    }
    node.forEachChild(visit);
  };
  visit(source);
  const roots = functions.filter(({ parent, start, end }) => parent === null
    && creations.some((creation) => creation >= start && creation < end));
  const written = roots.length === 1 ? functions.filter(({ parent }) => parent === roots[0]) : [];
  return {
    roots,
    logic: written.filter(({ statements }) => statements > 1),
    wiring: written.filter(({ statements }) => statements <= 1).length
  };
};

// Why a composition root module at `path` (under gbdraw/web/js/) has no single
// root function, or null when `roots` (compositionRootOf) holds exactly one.
const rootFunctionProblem = (path, roots) => (roots.length === 1 ? null
  : `gbdraw/web/js/${path}: ${roots.length} top-level functions create owners`
    + `${roots.length ? ` (${roots.map(({ name, line }) => `${name} at line ${line}`).join(', ')})` : ''}; `
    + 'a composition root creates its owners in one root function');

test('a function counts the statements in its body and is named by its binding', async () => {
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
    assert.deepEqual(functions.map(({ name, line, statements }) => `${line} ${name}: ${statements}`), [
      '2 outer: 12',
      '3 outer.inner: 1',
      '6 outer.forEach(): 1',
      '9 outer.ready: 1',
      '12 outer.doubled: 1',
      '13 outer.actions.run: 1',
      '13 outer.actions.stop: 1',
      '17 Store.constructor: 0',
      '18 Store.size: 1',
      '19 Store.handler: 0',
      '21 window.onload: 0',
      '22 pick: 1',
      '22 pick.choose: 1',
      '23 map(): 0',
      '24 (anonymous): 0',
      '25 default: 0'
    ]);
  } finally {
    rmSync(directory, { recursive: true, force: true });
  }
});

test('comments, JSDoc, blank lines and line breaks change no statement count', async () => {
  const directory = mkdtempSync(join(tmpdir(), 'web-function-statements-'));
  const compact = join(directory, 'compact.js');
  const spread = join(directory, 'spread.js');
  writeFileSync(compact, [
    'export function sum(items) { const total = items.reduce((all, item) => all + item, 0); '
      + 'if (total > 0) { return total; } return items.map((item) => { return item; }); }',
    ''
  ].join('\n'));
  writeFileSync(spread, [
    '/**',
    ' * @param {number[]} items',
    ' */',
    'export function sum(items) {',
    '  // A comment.',
    '  const total = items.reduce(',
    '    (all, item) => all + item,',
    '    0',
    '  );',
    '',
    '  if (total > 0) {',
    '    return total;',
    '  }',
    '  return items.map((item) => {',
    '    return item;',
    '  });',
    '}',
    ''
  ].join('\n'));
  try {
    const measured = await measure([compact, spread]);
    const counts = (file) => measured.get(file).map(({ name, statements }) => `${name}: ${statements}`);
    // total, the expression body of the reduce() callback, if, its return,
    // return, and the return of the map() callback.
    assert.deepEqual(counts(compact), ['sum: 6', 'sum.total: 1', 'sum.map(): 1']);
    assert.deepEqual(counts(spread), counts(compact));
  } finally {
    rmSync(directory, { recursive: true, force: true });
  }
});

test('no Web function has more statements than MAX_FUNCTION_STATEMENTS or its FUNCTION_STATEMENTS entry', async (t) => {
  assert.deepEqual(Object.keys(FUNCTION_STATEMENTS), Object.keys(FUNCTION_STATEMENTS).sort(),
    'keep FUNCTION_STATEMENTS sorted by key');
  Object.entries(FUNCTION_STATEMENTS).forEach(([key, statements]) => {
    assert.ok(Number.isInteger(statements) && statements > MAX_FUNCTION_STATEMENTS,
      `FUNCTION_STATEMENTS['${key}'] must be an integer above MAX_FUNCTION_STATEMENTS (${MAX_FUNCTION_STATEMENTS})`);
  });
  const paths = readdirSync(WEB_JS, { recursive: true })
    .map((path) => path.split(sep).join('/'))
    .filter((path) => path.endsWith('.js') && !path.endsWith('.generated.js'))
    .sort();
  const measured = await measure(paths.map((path) => join(WEB_JS, path)));
  // The largest function of each key.
  const largest = new Map();
  paths.forEach((path) => measured.get(join(WEB_JS, path)).forEach((found) => {
    const key = `${path}::${found.name}`;
    if (!largest.has(key) || found.statements > largest.get(key).statements) largest.set(key, { path, ...found });
  }));
  // A key without an entry counts as 0, so every function above the limit is
  // compared, and every listed one, to name the entries that may be lowered.
  const observed = Object.fromEntries([...largest]
    .filter(([key, { statements }]) => statements > MAX_FUNCTION_STATEMENTS || Object.hasOwn(FUNCTION_STATEMENTS, key))
    .map(([key, { statements }]) => [key, statements]));
  const { above, lowerable } = compareCountBaseline(observed, FUNCTION_STATEMENTS);
  const failures = above.map(({ key, observed: statements, recorded }) => {
    const { path, line, name } = largest.get(key);
    return `gbdraw/web/js/${path}:${line} ${name}: ${statements} statements (limit ${MAX_FUNCTION_STATEMENTS}, `
      + `FUNCTION_STATEMENTS ${recorded || 'none'})`;
  });
  t.diagnostic(`functions: ${[...measured.values()].reduce((sum, found) => sum + found.length, 0)} in ${paths.length} modules, `
    + `${Object.keys(observed).length} above ${MAX_FUNCTION_STATEMENTS} statements`);
  reportLowerable(t, 'FUNCTION_STATEMENTS', lowerable.map((entry) => entry.replace(/ -> 0$/, ' -> removed')));
  assert.deepEqual(failures, [], `functions above ${MAX_FUNCTION_STATEMENTS} statements, or above their FUNCTION_STATEMENTS `
    + 'entry. Split them at an owner boundary (gbdraw/web/CLAUDE.md "Module ownership"); adding or raising an entry '
    + `in FUNCTION_STATEMENTS needs the Owner's approval.\n  ${failures.join('\n  ')}`);
});

test('the logic of a composition root is every function of two or more statements written directly in its root function', async () => {
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
    '      const picked = item;',
    '      return picked;',
    '    },',
    '    getValue: () => 0,',
    '    getNext: () =>',
    '      1',
    '  });',
    '  const run = () => {',
    '    const nested = () => 1;',
    '    return nested();',
    '  };',
    '  const reset = () => { owner.value = 0; run(); };',
    '  watch(() => owner.value, () => {',
    '    run();',
    '  });',
    '  const steps = { first() { return 1; } };',
    '  return { owner, run, reset, steps, get size() { return 0; } };',
    '};',
    ''
  ].join('\n'));
  try {
    const { roots, logic, wiring } = (await measure([file], compositionRootOf)).get(file);
    assert.deepEqual(roots.map(({ name, line }) => `${line} ${name}`), ['6 createRoot']);
    assert.deepEqual(logic.map(({ name, line, statements }) => `${line} ${name}: ${statements}`), [
      '8 createRoot.owner.onPick: 2',
      '16 createRoot.run: 3',
      '20 createRoot.reset: 2'
    ]);
    // getValue, getNext, the getter and the callback of watch(), steps.first
    // and size: one statement each.
    assert.equal(wiring, 6);
  } finally {
    rmSync(directory, { recursive: true, force: true });
  }
});

test('a composition root module without exactly one root function fails by name', async () => {
  const directory = mkdtempSync(join(tmpdir(), 'web-root-function-'));
  const none = join(directory, 'none.js');
  const two = join(directory, 'two.js');
  writeFileSync(none, [
    "import { createOwner } from './owner.js';",
    'const createLocal = () => ({});',
    'export const createRoot = () => createLocal();',
    'export const owner = createOwner;',
    ''
  ].join('\n'));
  writeFileSync(two, [
    "import { createOwner, createOther } from './owner.js';",
    'export const createFirst = () => createOwner();',
    'export const createSecond = () => {',
    '  return createOther();',
    '};',
    ''
  ].join('\n'));
  try {
    const measured = await measure([none, two], compositionRootOf);
    assert.equal(rootFunctionProblem('app/none.js', measured.get(none).roots),
      'gbdraw/web/js/app/none.js: 0 top-level functions create owners; '
      + 'a composition root creates its owners in one root function');
    assert.equal(rootFunctionProblem('app/two.js', measured.get(two).roots),
      'gbdraw/web/js/app/two.js: 2 top-level functions create owners (createFirst at line 2, createSecond at line 3); '
      + 'a composition root creates its owners in one root function');
    assert.deepEqual(measured.get(two).logic, []);
  } finally {
    rmSync(directory, { recursive: true, force: true });
  }
});

test('the logic written in each composition root does not grow beyond COMPOSITION_ROOT_LOGIC_STATEMENTS', async (t) => {
  const recorded = COMPOSITION_ROOT_LOGIC_STATEMENTS;
  assert.deepEqual(Object.keys(recorded), Object.keys(recorded).sort(), 'keep COMPOSITION_ROOT_LOGIC_STATEMENTS sorted by key');
  Object.entries(recorded).forEach(([key, statements]) => {
    assert.ok(Number.isInteger(statements) && statements >= 0,
      `COMPOSITION_ROOT_LOGIC_STATEMENTS['${key}'] must be an integer of 0 or more`);
  });
  const paths = [...WEB_OWNER_GRAPH_DEFAULTS.compositionRoots].sort();
  const measured = await measure(paths.map((path) => join(WEB_JS, path)), compositionRootOf);
  const observed = {};
  const rows = new Map();
  const rootProblems = paths.flatMap((path) => {
    const { roots, logic, wiring } = measured.get(join(WEB_JS, path));
    const problem = rootFunctionProblem(path, roots);
    if (problem !== null) return [problem];
    const key = `${path}::${roots[0].name}`;
    observed[key] = logic.reduce((sum, { statements }) => sum + statements, 0);
    rows.set(key, { path, root: roots[0], functions: logic.length, wiring });
    return [];
  });
  assert.deepEqual(rootProblems, [], rootProblems.join('\n'));
  rows.forEach(({ root, functions, wiring }, key) => t.diagnostic(`${key}: ${observed[key]} statements in ${functions} `
    + `functions of two or more statements; ${wiring} functions of one statement or none not counted `
    + `(the root function has ${root.statements} statements)`));
  // An entry whose key is no root function now, whatever its value (0 too).
  const stale = Object.keys(recorded).filter((key) => !Object.hasOwn(observed, key));
  const rootFunctionOf = new Map([...rows.values()].map(({ path, root }) => [path, root.name]));
  const { above, lowerable } = compareCountBaseline(observed, recorded);
  reportLowerable(t, 'COMPOSITION_ROOT_LOGIC_STATEMENTS', [
    ...lowerable.filter((entry) => !stale.some((key) => entry.startsWith(`${key} `))),
    ...stale.map((key) => {
      const path = key.slice(0, key.indexOf('::'));
      return `${key} ${recorded[key]} -> removed (${rootFunctionOf.has(path)
        ? `the root function of ${path} is now ${rootFunctionOf.get(path)}` : `${path} is no longer a composition root`})`;
    })
  ]);
  const failures = [
    ...Object.keys(observed).filter((key) => !Object.hasOwn(recorded, key)).map((key) => {
      const { path, root } = rows.get(key);
      return `gbdraw/web/js/${path}:${root.line} ${root.name}: ${observed[key]} statements of logic in a composition `
        + 'root without a COMPOSITION_ROOT_LOGIC_STATEMENTS entry (a new composition root or a renamed root function)';
    }),
    ...above.filter(({ key }) => Object.hasOwn(recorded, key)).map(({ key, observed: statements, recorded: entry }) => {
      const { path, root, functions } = rows.get(key);
      return `gbdraw/web/js/${path}:${root.line} ${root.name}: ${statements} statements of logic in ${functions} `
        + `functions of two or more statements (COMPOSITION_ROOT_LOGIC_STATEMENTS ${entry})`;
    })
  ];
  assert.deepEqual(failures, [], 'the logic written in a composition root grew. Move the new logic into a module '
    + 'factory and pass it to the root as a port, or shrink other logic in the same root '
    + '(docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md "Web function size"); adding or raising an entry in '
    + `COMPOSITION_ROOT_LOGIC_STATEMENTS needs the Owner's approval.\n  ${failures.join('\n  ')}`);
});
