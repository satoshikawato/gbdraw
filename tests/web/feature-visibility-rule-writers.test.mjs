// OV-19 (PD-OI-066; gbdraw/web/CLAUDE.md R1, R10, R11): one transition,
// `editFeatureVisibilityRules` in app/feature-editor/visibility-actions.js,
// writes the Feature Visibility rules for an edit and projects them, so a rule
// edit is live as Generate draws it. This guard lists every writer of
// `featureVisibilityManualRules` in the Web modules and the template and fails
// on a writer outside the owners below.
//
// R2: the paint intent of a drawing (its Legend rows and their deleted,
// dormant, and added lists, and the Legend and feature color and stroke
// edits) has the writers `PAINT_INTENT_WRITERS` lists, each by the route that
// reaches it; a live edit goes through `editEditorIntent` in app/app-setup.js.
// The same scan pins both: a write through a name (a mutating method, an
// assignment to it, an element, or a member, a `delete`), or a call that
// receives it, in its enclosing named function. Blind spot (P10): a write
// through an alias is not seen, such as `entry.color = x` on a row read from
// the list, or `delete overrides[key]` on `const overrides = drawing.<store>`;
// the function that made the alias is listed only when it writes by name too.
// Nor is a store handed on inside an object literal under another name, or a
// ref replaced without `.value` (`drawing.legendEntries = ref([])`). A write
// inside a shorthand method is reported under the function name ''.
// A drawing slice (js/services/drawing/) is no owner: it writes its members
// only for its caller, so a call of a slice function that receives the handle
// of a slice holding one of the names (`drawing.slices.<slice>`, or a name
// equal to the handle's, such as a helper's `paint` parameter) is a write of
// the function that calls it, except the slice readers that `SLICE_READERS`
// lists by name (tests/web/helpers/drawing-state.mjs; the M6 table proves
// each changes neither drawing) and the calls `READERS` lists. A handle held
// under another name is not seen (the slices' own guard, tests/web/drawing-slices.test.mjs, covers
// their members).
import assert from 'node:assert/strict';
import { readdirSync, readFileSync } from 'node:fs';
import test from 'node:test';
import { DRAWING_SLICES, SLICE_READERS } from './helpers/drawing-state.mjs';
import { NOT_CALLS as KEYWORDS, enclosingCallee, sliceHandleCalls } from './helpers/drawing-writes.mjs';

const WEB_ROOT = new URL('../../gbdraw/web/', import.meta.url);
const NAME = 'featureVisibilityManualRules';

// The functions that may write the rules, by module.
export const RULE_WRITERS = Object.freeze({
  // The edit transition: the Features panel and the popup's rule scopes.
  'js/app/feature-editor/visibility-actions.js': ['editFeatureVisibilityRules'],
  // Undo and Redo install a captured copy (R11); the History apply projects it.
  'js/services/history-snapshot.js': ['replaceFeatureEditState'],
  // A Session replacement installs the saved rules with the saved Result.
  'js/services/config.js': ['resetSessionBaseline', 'applyDrawingFeatureData'],
  // Reset Settings, which keeps the current Results, and the Session baseline.
  'js/services/reset.js': ['resetEditorDraftState']
});

// The paint intent: the Vue refs (written through `.value`) and the stores.
const PAINT_INTENT_REFS = ['legendEntries', 'deletedLegendEntries', 'dormantLegendEntries', 'addedLegendCaptions'];
const PAINT_INTENT_NAMES = [
  ...PAINT_INTENT_REFS, 'legendColorOverrides', 'legendStrokeOverrides', 'featureColorOverrides', 'featureStrokeOverrides'
];
// The routes that write the paint intent; each holds one History step and
// shows the edit through `showEditorIntent` unless it says otherwise.
const PAINT_ROUTES = new Set([
  // An action `editEditorIntent` wraps (`history.runUndoable`, `showingEditorIntent`).
  'edit',
  // A Change Scope dialog choice (`dialogChoice.withHistory`): the stroke choice
  // shows through `showingEditorIntent`, the color choice through the rule commit,
  // or through the palette when it sets the type's default color.
  'scope',
  // The rule commit (`commitOnce`, R13): the rules and the Legend rows, captions,
  // and colors they draw, shown by `projectPaletteAndRules` with the rule rows.
  'rule',
  // Restore of deleted Legend rows (`restoreLegendItems`, `restoreLegendRows`).
  'restore',
  // Sort and Move of Legend rows (`showLegendStructure`); Sort's step is the
  // History input adapter's button step, Move's `history.runUndoable`.
  'order',
  // Another owner's edit retires the paint edits of what it removed, in that
  // edit's step: track data that no longer names a caption (OV-87), and Delete
  // settings of OFF records, which no Result draws.
  'retire',
  // The Legend name input shows its row's name again after a refused or
  // cancelled rename: a copy of the list with the same rows, outside a step.
  'redisplay',
  // Undo and Redo install a captured step (R11).
  'history',
  // A Session load, replacement, or baseline; a loaded Session's payload.
  'session',
  // Reset Settings.
  'reset',
  // Generate and the Result binding: Python's Legend rows and feature colors,
  // and the edits a source-replacing Generate drops.
  'generate'
]);

// The functions that may write the paint intent, by module, with their route.
export const PAINT_INTENT_WRITERS = {
  'js/app/feature-editor/color-actions.js': {
    // The commit callback of `applyColorToFeatureGroup`.
    afterCommit: 'rule',
    applyStrokeToLegendEntry: 'scope',
    clearFeatureStrokeOverride: 'edit',
    clearLegendRenameDialog: 'redisplay',
    handleColorScopeChoice: 'scope',
    // Also the stroke scope choice (`applyStrokeToSelectedFeatures`).
    recordFeatureStrokeOverride: 'edit',
    renameLegendEntry: 'redisplay',
    setClickedFeatureStrokeColorValue: 'edit'
  },
  'js/app/feature-editor/rule-actions.js': { commitOnce: 'rule', refreshFeatureOverrides: 'rule' },
  'js/app/legend/entry-actions.js': {
    adoptRestoredLegend: 'history',
    // The rule commit's Legend rows (`prepareFileLegendEntries`).
    apply: 'rule',
    deleteLegendEntry: 'edit',
    // Also Reset Settings, after its projection.
    extractLegendEntries: 'generate',
    moveLegendRecords: 'rule',
    restoreDeletedLegendEntries: 'restore',
    // Also the rule commit's callback and the color scope choice.
    writeLegendRowColor: 'edit'
  },
  'js/app/legend/sort-actions.js': { applyLegendEntryOrder: 'order' },
  'js/app/legend/stroke-actions.js': {
    resetAllStrokes: 'edit', resetLegendEntryStroke: 'edit', setLegendEntryStrokeColorValue: 'edit', writeStrokeField: 'edit'
  },
  'js/app/legend/track-data-styles.js': { buildLegendStyleRetirement: 'retire' },
  'js/app/run-analysis.js': { runAnalysisInternal: 'generate' },
  // Also History's checkpoint restore (`applyEditorStateData`, `applyFeatureStateData`).
  'js/services/config.js': { applyDrawingEditorData: 'session', applyDrawingFeatureData: 'session', resetSessionBaseline: 'session' },
  'js/services/feature-visibility.js': { pruneUnmatchedFeatureOverrides: 'generate' },
  'js/services/history-snapshot.js': {
    applyEditorIntentData: 'history', applyFallbackFeatureStateData: 'history', applyFeatureIntentData: 'history'
  },
  'js/services/record-draw-selection.js': { forEachRecordEdit: 'retire' },
  // Also the Session baseline (`resetSessionBaseline`).
  'js/services/reset.js': { resetEditorDraftState: 'reset' },
  'js/services/session-feature-recovery.js': { migrateFeatureOverrideState: 'session' }
};

// The template binds the rule rows to these transition actions only.
const ROW_ACTIONS = new Set([
  'setFeatureVisibilityRuleField', 'moveFeatureVisibilityRuleUp', 'moveFeatureVisibilityRuleDown',
  'removeFeatureVisibilityRule'
]);
const MUTATING_METHODS = new Set([
  'push', 'pop', 'shift', 'unshift', 'splice', 'sort', 'reverse', 'fill', 'copyWithin', 'add', 'delete', 'clear', 'set'
]);
// Calls that read the values they are given; any other call that receives one is a writer.
const READERS = new Set([
  'Array.from', 'Array.isArray', 'Object.entries', 'Object.hasOwn', 'Object.keys', 'clonePlainObject', 'cloneJsonArray',
  'cloneJsonData', 'cloneJsonObject', 'cloneJsonValue', 'cloneFeatureVisibilityRules', 'getFeatureOverride',
  'getGeneratedArtifactRef', 'getRef', 'hasOwn', 'hasOwnEntry', 'isLegendOrderEdited', 'normalizeFeatureVisibilityRulesForSession',
  'normalizedLegendEntries', 'publicationClone', 'readFeatureColorEdits', 'requestFeatureVisibilityRules', 'rewriteOverrideMap',
  'rowOf', 'savedFeatureColorEdits', 'serializeFeatureVisibilityRules'
]);
const ASSIGNMENT = /^\s*(?:[-+*/%&|^]|\*\*|<<|>>>?|&&|\|\||\?\?)?=(?![=>])/;

// Comments and quoted strings blanked, so offsets stay and prose is no writer.
const codeOnly = (source) => source.replace(
  /\/\*[\s\S]*?\*\/|\/\/[^\n]*|'(?:\\.|[^'\\\n])*'|"(?:\\.|[^"\\\n])*"/g,
  (text) => text.replace(/[^\n]/g, ' ')
);

// Whether the code after a name writes through it: a mutating method, or an
// assignment to it, an element, or a member.
const writesThrough = (after) => {
  let rest = after;
  let last = '';
  for (let step; (step = rest.match(/^\s*(?:\??\.\s*([\w$]+)|\[[^\]]*\])/)); rest = rest.slice(step[0].length)) {
    last = step[1] || '';
  }
  return (MUTATING_METHODS.has(last) && /^\s*\(/.test(rest)) || ASSIGNMENT.test(rest);
};

// The index after the bracket group that opens at `start`.
const groupEnd = (code, start) => {
  for (let at = start, depth = 0; at < code.length; at += 1) {
    if ('([{'.includes(code[at])) depth += 1;
    else if (')]}'.includes(code[at]) && --depth === 0) return at + 1;
  }
  return code.length;
};

// Each named function's body as [name, start, end]: a declared function, or
// a function assigned to a name or a property. A callback has no name.
const namedFunctions = (code) => [...code.matchAll(
  /\b(?:(?:const|let|var)\s+([\w$]+)\s*=|([\w$]+)\s*:)\s*(?:async\b\s*)?(?:function\b\s*\*?\s*[\w$]*\s*)?|\bfunction\s*\*?\s*([\w$]+)\s*/g
)].flatMap((match) => {
  let at = match.index + match[0].length;
  if (code[at] === '(') at = groupEnd(code, at);
  else if (/^[\w$]+\s*=>/.test(code.slice(at, at + 64))) at += code.slice(at).match(/^[\w$]+/)[0].length;
  else return [];
  at += code.slice(at).match(/^\s*/)[0].length;
  const arrow = code.startsWith('=>', at);
  if (arrow) at += 2 + code.slice(at + 2).match(/^\s*/)[0].length;
  else if (!/function\b/.test(match[0]) || code[at] !== '{') return [];
  let end = at;
  if (code[at] === '{') end = groupEnd(code, at);
  else while (end < code.length && !',;)]}'.includes(code[end])) end = '([{'.includes(code[end]) ? groupEnd(code, end) : end + 1;
  return [[match.slice(1).find(Boolean), at, end]];
});

// The innermost named function whose body holds `index`.
const enclosingFunction = (functions, index) => functions
  .filter(([, start, end]) => start <= index && index < end).at(-1)?.[0] || '';

// The writer sites of `names` in one module as `function:line`: a write
// through a name (a mutating method, an assignment to it, an element, or a
// member, a `delete`), or a call that receives it. A name assigned without a
// member access is a local binding. A name in `refs` is a Vue ref: its value
// is written through `.value`, and the ref itself only by a call.
export const writerSites = (source, names, refs = []) => {
  const code = codeOnly(source);
  let functions = null;
  const at = (start) => {
    functions ||= namedFunctions(code);
    return `${enclosingFunction(functions, start)}:${code.slice(0, start).split('\n').length}`;
  };
  const handles = Object.entries(DRAWING_SLICES)
    .filter(([, slice]) => slice.members.some((member) => names.includes(member))).map(([handle]) => handle);
  const sliceCalls = sliceHandleCalls(code, handles)
    .filter(({ callee }) => !READERS.has(callee) && !SLICE_READERS.includes(callee)).map(({ index }) => at(index));
  return [...sliceCalls, ...[...code.matchAll(new RegExp(`(?<![\\w$])(?:${names.join('|')})\\b`, 'g'))].flatMap((match) => {
    const start = match.index;
    const before = code.slice(0, start);
    const ref = refs.includes(match[0]);
    const value = ref ? code.slice(start + match[0].length).match(/^\s*\??\.\s*value\b/)?.[0] || '' : '';
    const after = code.slice(start + match[0].length + value.length);
    const direct = ASSIGNMENT.test(after);
    const stepped = /^\s*(?:\??\.|\[)/.test(after);
    const callee = stepped || direct ? '' : enclosingCallee(code, start);
    const handed = Boolean(callee) && !KEYWORDS.has(callee) && !READERS.has(callee);
    const through = direct ? Boolean(value) || /\??\.\s*$/.test(before) : stepped && writesThrough(after);
    const writes = !/\b(?:const|let|var)\s+$/.test(before) && (handed || ((!ref || value)
      && (through || /\bdelete\s+\(?\s*(?:[\w$]+\s*\??\.\s*)*$/.test(before))));
    if (!writes) return [];
    return [at(start)];
  })];
};

export const ruleWriters = (source) => writerSites(source, [NAME]);
export const paintWriters = (source) => writerSites(source, PAINT_INTENT_NAMES, PAINT_INTENT_REFS);

const BINDINGS = /\s(v-model(?:\.[\w.]+)?|@[\w.:-]+|v-on:[\w.:-]+)="([^"]*)"/g;
const lineIn = (html, index) => html.slice(0, index).split('\n').length;

// The template writers of `names`: a `v-model` on one, and a handler that
// writes through one.
const templateWriters = (html, names) => [...html.matchAll(BINDINGS)].flatMap((match) => {
  const [, attribute, expression] = match;
  const code = codeOnly(expression);
  const uses = [...code.matchAll(new RegExp(`(?<![\\w$])(?:${names.join('|')})\\b`, 'g'))];
  return (attribute.startsWith('v-model') ? uses.length > 0
    : uses.some((use) => writesThrough(code.slice(use.index + use[0].length))))
    ? [`${attribute}="${expression}":${lineIn(html, match.index)}`] : [];
});

// The markup of the element whose start tag holds `index`.
const elementAt = (html, index) => {
  const start = html.lastIndexOf('<', index);
  const tag = new RegExp(`<(/?)${html.slice(start + 1).match(/^[\w-]+/)[0]}\\b`, 'g');
  tag.lastIndex = start;
  for (let depth = 0, match; (match = tag.exec(html));) {
    if ((depth += match[1] ? -1 : 1) === 0) return html.slice(start, html.indexOf('>', match.index) + 1);
  }
  return html.slice(start);
};

// The paint writers of the template: a binding above, and in a row of a
// `v-for` over the paint intent, a `v-model` on the row or a handler that
// writes through it.
export const templatePaintWriters = (html) => [
  ...templateWriters(html, PAINT_INTENT_NAMES),
  ...[...html.matchAll(new RegExp(`v-for="\\(?([\\w$]+)[^"]* in (?:${PAINT_INTENT_NAMES.join('|')})"`, 'g'))]
    .flatMap(({ 1: row, index }) => {
      const start = html.lastIndexOf('<', index);
      const rows = elementAt(html, index);
      return [...rows.matchAll(BINDINGS)].filter(([, attribute, expression]) => {
        const code = codeOnly(expression);
        const uses = [...code.matchAll(new RegExp(`(?<![\\w$.])${row}\\b`, 'g'))];
        return attribute.startsWith('v-model') ? uses.length > 0
          : uses.some((use) => writesThrough(code.slice(use.index + row.length)));
      }).map((match) => `${row} row ${match[0].trim()}:${lineIn(html, start + match.index)}`);
    })
];

// The template writers: a `v-model` on the rules or in a rule row, a handler
// that writes through the rules, and a rule-row handler that assigns or calls
// an action other than the transition's.
export const templateRuleWriters = (html) => {
  const writers = templateWriters(html, [NAME]);
  const rowsStart = html.search(new RegExp(`v-for="\\(rule, i\\) in ${NAME}"`));
  const rowsEnd = html.indexOf(`v-if="${NAME}.length === 0"`, rowsStart);
  if (rowsStart < 0 || rowsEnd < 0) return [...writers, 'rule rows: markup not found; update this guard'];
  const rows = html.slice(rowsStart, rowsEnd);
  for (const match of rows.matchAll(BINDINGS)) {
    const [, attribute, expression] = match;
    const action = expression.match(/^\s*([\w$]+)\s*\(/)?.[1];
    const assigns = /(?<![=!<>])=(?![=>])/.test(codeOnly(expression));
    if (attribute.startsWith('v-model') || assigns || !ROW_ACTIONS.has(action)) {
      writers.push(`rule row ${attribute}="${expression}":${lineIn(html, rowsStart + match.index)}`);
    }
  }
  return writers;
};

const listModules = (url, prefix = '') => readdirSync(url, { withFileTypes: true }).flatMap((entry) => {
  const path = `${prefix}${entry.name}`;
  if (entry.isDirectory()) return listModules(new URL(`${entry.name}/`, url), `${path}/`);
  return entry.name.endsWith('.js') ? [path] : [];
});

// The writer sites of the Web modules that have one, by module; the slices
// write for their callers.
const scanModules = (writers) => Object.fromEntries(listModules(new URL('js/', WEB_ROOT), 'js/')
  .filter((file) => !file.startsWith('js/services/drawing/'))
  .map((file) => [file, writers(readFileSync(new URL(file, WEB_ROOT), 'utf8'))])
  .filter(([, sites]) => sites.length));

// The writers outside `owners` (a module's functions, listed or mapped to
// their route), and the listed functions that no longer write.
const writerProblems = (writersByFile, owners, { writes, route, list }) => {
  const listed = (file) => (Array.isArray(owners[file]) ? owners[file] : Object.keys(owners[file] || {}));
  const problems = [];
  for (const [file, writers] of Object.entries(writersByFile)) {
    writers.filter((writer) => !listed(file).includes(writer.split(':')[0]))
      .forEach((writer) => problems.push(`${file}: ${writer} writes ${writes}; route the edit through ${route}`));
  }
  for (const file of Object.keys(owners)) {
    listed(file).filter((name) => !(writersByFile[file] || []).some((writer) => writer.startsWith(`${name}:`)))
      .forEach((name) => problems.push(`${file}: ${name} no longer writes ${writes}; remove it from ${list}`));
  }
  return problems;
};

export const ruleWriterProblems = (writersByFile, owners = RULE_WRITERS) => writerProblems(writersByFile, owners,
  { writes: NAME, route: 'editFeatureVisibilityRules', list: 'RULE_WRITERS' });
export const paintWriterProblems = (writersByFile, owners = PAINT_INTENT_WRITERS) => writerProblems(writersByFile, owners,
  { writes: 'the paint intent', route: 'editEditorIntent (app/app-setup.js)', list: 'PAINT_INTENT_WRITERS' });

test('only the rule edit transition, History restore, Session load, and Reset Settings write the visibility rules (OV-19)', () => {
  const writersByFile = scanModules(ruleWriters);
  assert.deepEqual(writersByFile['js/app/feature-editor/visibility-actions.js']?.length, 1,
    'the transition writes the rules once');
  assert.deepEqual(ruleWriterProblems(writersByFile), []);
  assert.deepEqual(templateRuleWriters(readFileSync(new URL('index.html', WEB_ROOT), 'utf8')), []);
});

test('only editEditorIntent and the routes PAINT_INTENT_WRITERS names write the paint intent of a drawing (R2)', () => {
  assert.deepEqual(paintWriterProblems(scanModules(paintWriters)), []);
  assert.deepEqual(Object.values(PAINT_INTENT_WRITERS).flatMap(Object.values).filter((route) => !PAINT_ROUTES.has(route)), []);
  assert.deepEqual(templatePaintWriters(readFileSync(new URL('index.html', WEB_ROOT), 'utf8')), []);
});

test('the writer guard flags a direct edit and leaves reads alone', () => {
  const writers = ruleWriters([
    'const featureVisibilityManualRules = reactive([]);',
    'const addRule = () => {',
    '  featureVisibilityManualRules.push(rule);',
    '  state.featureVisibilityManualRules[0].value = "x";',
    '  upsertEditorQualifierFeatureVisibilityRule(state.featureVisibilityManualRules, input, mode);',
    '};',
    'const readRules = () => {',
    '  // featureVisibilityManualRules.splice(0) in prose is no writer',
    '  const rows = [...requestFeatureVisibilityRules(state?.featureVisibilityManualRules)];',
    '  if (featureVisibilityManualRules.length === 0) return [...featureVisibilityManualRules];',
    '  return featureVisibilityManualRules.map((rule) => rule.value);',
    '};'
  ].join('\n'));
  assert.deepEqual(writers, ['addRule:3', 'addRule:4', 'addRule:5']);
  // A slice writer that receives the rules' slice handle writes; its capture reads.
  const [handle] = Object.entries(DRAWING_SLICES).find(([, slice]) => slice.members.includes(NAME));
  assert.deepEqual(ruleWriters([
    'const restore = (drawing, rules) => {',
    `  SLICE.install(drawing.slices.${handle}, rules);`,
    `  return VISIBILITY_RULES_SLICE.capture(drawing?.slices.${handle});`,
    '};',
    // A helper's parameter named as the handle is the handle.
    `const prune = (${handle}, key) => {`,
    `  if (Object.keys(SLICE.rulesOf(${handle})).length) removeRule(${handle}, key);`,
    `  return { ${handle}: ${handle}.rules, count: [${handle}].length };`,
    '};'
  ].join('\n')), ['restore:2', 'prune:6', 'prune:6']);
  // A slice reader is one by the shared list, not by its name.
  assert.deepEqual(paintWriters([
    'const clear = (drawing, paint) => {',
    '  clearStrokesOf(drawing.slices.paint);',
    '  captureAndClearStrokes(paint);',
    '  return [paintOf(paint), PAINT_SLICE.capture(paint), captureLegendPaint(drawing.slices.paint)];',
    '};'
  ].join('\n')), ['clear:2', 'clear:3']);
  assert.match(ruleWriterProblems({ 'js/app/new-import.js': ['importRules:9'] }).join('\n'),
    /js\/app\/new-import\.js: importRules:9 writes featureVisibilityManualRules/);
  assert.match(ruleWriterProblems({}).join('\n'), /editFeatureVisibilityRules no longer writes/);

  const rows = (inner) => `<div v-for="(rule, i) in featureVisibilityManualRules">${inner}</div>`
    + '<div v-if="featureVisibilityManualRules.length === 0"></div>';
  assert.deepEqual(templateRuleWriters(rows(
    '<input :value="rule.value" @change="setFeatureVisibilityRuleField(i, \'value\', $event.target.value)">'
  )), []);
  assert.equal(templateRuleWriters(rows('<input v-model="rule.value">')).length, 1);
  assert.equal(templateRuleWriters(rows('<button @click="rule.action = \'off\'"></button>')).length, 1);
  assert.equal(templateRuleWriters(rows('') + '<button @click="featureVisibilityManualRules.splice(0)"></button>').length, 1);

  const paint = paintWriters([
    'const editRows = (drawing) => {',
    '  drawing.legendEntries.value = [];',
    '  delete drawing.legendColorOverrides[key];',
    '  replacePlainObject(drawing.featureStrokeOverrides, next);',
    '  drawing.addedLegendCaptions.value.add(caption);',
    '  drawing?.featureColorOverrides?.clear();',
    '  state.legendStrokeOverrides.set(key, stroke);',
    '  installRows(drawing.dormantLegendEntries);',
    '};',
    'const closeRename = (drawing) => closeAfterChoice(() => {',
    '  drawing.deletedLegendEntries.value.push(entry);',
    '});',
    'const readRows = ({ legendEntries = [], featureColorOverrides = {} } = {}) => {',
    '  // drawing.legendEntries.value = [] in prose is no writer',
    '  const keys = Object.keys(drawing.legendStrokeOverrides);',
    '  built.legendEntries = new Map();',
    '  return drawing.legendEntries.value.filter(Boolean).length + legendEntries.length;',
    '};'
  ].join('\n'));
  assert.deepEqual(paint, [
    'editRows:2', 'editRows:3', 'editRows:4', 'editRows:5', 'editRows:6', 'editRows:7', 'editRows:8', 'closeRename:11'
  ]);
  assert.match(paintWriterProblems({ 'js/app/legend/new-import.js': ['importLegend:4'] }).join('\n'),
    /js\/app\/legend\/new-import\.js: importLegend:4 writes the paint intent; route the edit through editEditorIntent/);
  assert.match(paintWriterProblems({}, { 'js/app/legend/sort-actions.js': { applyLegendEntryOrder: 'order' } }).join('\n'),
    /applyLegendEntryOrder no longer writes the paint intent; remove it from PAINT_INTENT_WRITERS/);

  const legendRows = (inner) => `<div v-for="(entry, idx) in legendEntries" :key="idx"><div>${inner}</div></div>`
    + '<div><input v-model="entry.caption"></div>';
  assert.deepEqual(templatePaintWriters(legendRows(
    '<input :value="entry.caption" @change="renameLegendEntry(idx, $event.target.value)">'
  )), []);
  assert.equal(templatePaintWriters(legendRows('<input v-model="entry.caption">')).length, 1);
  assert.equal(templatePaintWriters(legendRows('<button @click="entry.color = \'#000\'"></button>')).length, 1);
  assert.equal(templatePaintWriters('<button @click="legendColorOverrides[key] = color"></button>').length, 1);
});
