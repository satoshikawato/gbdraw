// OV-19 (PD-OI-066; gbdraw/web/CLAUDE.md R1, R10, R11): one transition,
// `editFeatureVisibilityRules` in app/feature-editor/visibility-actions.js,
// writes the Feature Visibility rules for an edit and projects them, so a rule
// edit is live as Generate draws it. This guard lists every writer of
// `featureVisibilityManualRules` in the Web modules and the template and fails
// on a writer outside the owners below.
import assert from 'node:assert/strict';
import { readdirSync, readFileSync } from 'node:fs';
import test from 'node:test';

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

// The template binds the rule rows to these transition actions only.
const ROW_ACTIONS = new Set([
  'setFeatureVisibilityRuleField', 'moveFeatureVisibilityRuleUp', 'moveFeatureVisibilityRuleDown',
  'removeFeatureVisibilityRule'
]);
const MUTATING_METHODS = new Set(['push', 'pop', 'shift', 'unshift', 'splice', 'sort', 'reverse', 'fill', 'copyWithin']);
// Calls that read the rules they are given; any other call that receives them is a writer.
const READERS = new Set([
  'Array.isArray', 'cloneFeatureVisibilityRules', 'normalizeFeatureVisibilityRulesForSession', 'publicationClone',
  'requestFeatureVisibilityRules', 'serializeFeatureVisibilityRules'
]);
const KEYWORDS = new Set(['if', 'while', 'for', 'switch', 'catch', 'return', 'typeof', 'await', 'void']);
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
  for (let step; (step = rest.match(/^\s*(?:\.\s*([\w$]+)|\[[^\]]*\])/)); rest = rest.slice(step[0].length)) {
    last = step[1] || '';
  }
  return (MUTATING_METHODS.has(last) && /^\s*\(/.test(rest)) || ASSIGNMENT.test(rest);
};

// The call whose arguments hold `index`, or '' outside a call.
const enclosingCallee = (code, index) => {
  const closers = [];
  for (let at = index - 1; at >= 0; at -= 1) {
    const char = code[at];
    if (')]}'.includes(char)) closers.push(char);
    else if ('([{'.includes(char)) {
      if (closers.length) closers.pop();
      else return char === '(' ? (code.slice(0, at).match(/([\w$]+(?:\.[\w$]+)*)\s*$/)?.[1] || '') : '';
    }
  }
  return '';
};

const enclosingFunction = (code, index) => [...code.slice(0, index).matchAll(
  /(?:\b(?:const|let|var)\s+([\w$]+)\s*=\s*(?:async\s*)?(?:\([^)]*\)|[\w$]+)\s*=>|\bfunction\s+([\w$]+)\s*\()/g
)].at(-1)?.slice(1).find(Boolean) || '';

// The writer sites of one module as `function:line`.
export const ruleWriters = (source) => {
  const code = codeOnly(source);
  return [...code.matchAll(new RegExp(`(?<![\\w$])${NAME}\\b`, 'g'))].flatMap((match) => {
    const start = match.index;
    const end = start + NAME.length;
    const after = code.slice(end);
    const declaration = /\b(?:const|let|var)\s+$/.test(code.slice(0, start));
    const bare = !/^\s*(?:\.|\[)/.test(after) && !ASSIGNMENT.test(after);
    const callee = bare ? enclosingCallee(code, start) : '';
    const writes = !declaration && (writesThrough(after)
      || (bare && callee && !KEYWORDS.has(callee) && !READERS.has(callee)));
    return writes ? [`${enclosingFunction(code, start)}:${code.slice(0, start).split('\n').length}`] : [];
  });
};

// The template writers: a `v-model` on the rules or in a rule row, a handler
// that writes through the rules, and a rule-row handler that assigns or calls
// an action other than the transition's.
export const templateRuleWriters = (html) => {
  const writers = [];
  const lineOf = (index) => html.slice(0, index).split('\n').length;
  for (const match of html.matchAll(/\s(v-model(?:\.[\w.]+)?|@[\w.:-]+|v-on:[\w.:-]+)="([^"]*)"/g)) {
    const [, attribute, expression] = match;
    const code = codeOnly(expression);
    const uses = [...code.matchAll(new RegExp(`(?<![\\w$])${NAME}\\b`, 'g'))];
    if (attribute.startsWith('v-model') ? uses.length > 0
      : uses.some((use) => writesThrough(code.slice(use.index + NAME.length)))) {
      writers.push(`${attribute}="${expression}":${lineOf(match.index)}`);
    }
  }
  const rowsStart = html.search(new RegExp(`v-for="\\(rule, i\\) in ${NAME}"`));
  const rowsEnd = html.indexOf(`v-if="${NAME}.length === 0"`, rowsStart);
  if (rowsStart < 0 || rowsEnd < 0) return [...writers, 'rule rows: markup not found; update this guard'];
  const rows = html.slice(rowsStart, rowsEnd);
  for (const match of rows.matchAll(/\s(v-model(?:\.[\w.]+)?|@[\w.:-]+|v-on:[\w.:-]+)="([^"]*)"/g)) {
    const [, attribute, expression] = match;
    const action = expression.match(/^\s*([\w$]+)\s*\(/)?.[1];
    const assigns = /(?<![=!<>])=(?![=>])/.test(codeOnly(expression));
    if (attribute.startsWith('v-model') || assigns || !ROW_ACTIONS.has(action)) {
      writers.push(`rule row ${attribute}="${expression}":${lineOf(rowsStart + match.index)}`);
    }
  }
  return writers;
};

const listModules = (url, prefix = '') => readdirSync(url, { withFileTypes: true }).flatMap((entry) => {
  const path = `${prefix}${entry.name}`;
  if (entry.isDirectory()) return listModules(new URL(`${entry.name}/`, url), `${path}/`);
  return entry.name.endsWith('.js') && !path.endsWith('.generated.js') ? [path] : [];
});

export const ruleWriterProblems = (writersByFile, owners = RULE_WRITERS) => {
  const problems = [];
  for (const [file, writers] of Object.entries(writersByFile)) {
    writers.filter((writer) => !(owners[file] || []).includes(writer.split(':')[0]))
      .forEach((writer) => problems.push(`${file}: ${writer} writes ${NAME}; route the edit through editFeatureVisibilityRules`));
  }
  for (const [file, functions] of Object.entries(owners)) {
    functions.filter((name) => !(writersByFile[file] || []).some((writer) => writer.startsWith(`${name}:`)))
      .forEach((name) => problems.push(`${file}: ${name} no longer writes ${NAME}; remove it from RULE_WRITERS`));
  }
  return problems;
};

test('only the rule edit transition, History restore, Session load, and Reset Settings write the visibility rules (OV-19)', () => {
  const writersByFile = Object.fromEntries(listModules(new URL('js/', WEB_ROOT), 'js/')
    .map((file) => [file, ruleWriters(readFileSync(new URL(file, WEB_ROOT), 'utf8'))])
    .filter(([, writers]) => writers.length));
  assert.deepEqual(writersByFile['js/app/feature-editor/visibility-actions.js']?.length, 1,
    'the transition writes the rules once');
  assert.deepEqual(ruleWriterProblems(writersByFile), []);
  assert.deepEqual(templateRuleWriters(readFileSync(new URL('index.html', WEB_ROOT), 'utf8')), []);
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
});
