// @ts-check
// The Undo and Redo titles translate a step's label where it is shown
// (`undoHistoryTitle` in app/app-setup.js), by the label as the key: a label the
// catalogs do not hold shows English in every language (OV-2415).
import assert from 'node:assert/strict';
import test from 'node:test';

import { skipExpression, sourceCalls, sourceFiles, sourceText, webPath } from '../../tools/i18n/source-keys.mjs';

// The History calls whose first argument names a step.
const LABEL_TAKERS = /(?<![\w$])(?:runUndoable|runUndoableCheckpoint|runUndoableArtifactReplacement|runUndoableCommand|runProjectStep|begin|beginCheckpoint|beginArtifactReplacement|beginHistoryTransaction|undoableAction|editEditorIntent|withHistory|withDialogChoice|scopeChoiceWithHistory|restoreLegendItems)\(\s*/g;
// Two wrappers take the label later: as the second argument, or as `label:` in it.
const SECOND_ARGUMENT_TAKER = /(?<![\w$])commitSpecificRules\(\s*/g;
const LABEL_OPTION_TAKER = /(?<![\w$])askTrackHide\(\s*/g;
// A literal the argument evaluates to: on its own, or the result of an arrow
// function, a ternary, or a fallback. A literal inside `msg()` or a comparison is not.
const RESULT_LITERAL = /(?:^|=>|\?\??|:|\|\||\breturn\s)\s*(['"])((?:\\.|(?!\1)[^\\])*)\1/g;
// A label is stored and translated where it is shown, so it must not be translated here.
const TRANSLATED_NOW = /(?<![\w$.])(?:t|tn|tx|tm)\(/;

/** @param {string} source @param {number} start the first character of an argument */
const argumentAt = (source, start) => source.slice(start, skipExpression(source, start)).trim();

/** @param {string} source @param {RegExpExecArray | RegExpMatchArray} match a taker's opening */
const secondArgumentStart = (source, match) => skipExpression(source, Number(match.index) + match[0].length) + 1;

/**
 * Every argument of `source` that names a History step.
 * @param {string} source
 * @returns {string[]}
 */
const labelArguments = (source) => [
  ...[...source.matchAll(LABEL_TAKERS)].map((match) => argumentAt(source, match.index + match[0].length)),
  ...[...source.matchAll(SECOND_ARGUMENT_TAKER)].map((match) => argumentAt(source, secondArgumentStart(source, match))),
  ...[...source.matchAll(LABEL_OPTION_TAKER)].flatMap((match) => {
    const start = secondArgumentStart(source, match);
    const option = /\blabel:\s*/.exec(source.slice(start, skipExpression(source, start)));
    return option ? [argumentAt(source, start + option.index + option[0].length)] : [];
  })
];

/**
 * The label arguments of `source` that no catalog holds: any template literal, a
 * text joined at run time, a text translated too early, or a literal result that
 * is no key. A `msg()` literal, and a name that forwards a label from another
 * call, pass.
 * @param {string} source @param {ReadonlySet<string>} keys
 * @returns {string[]}
 */
const unkeyedLabels = (source, keys) => labelArguments(source).filter((argument) => (
  argument.includes('`') || /['"]\s*\+|\+\s*['"]/.test(argument) || TRANSLATED_NOW.test(argument)
  || [...argument.matchAll(RESULT_LITERAL)].some((literal) => !keys.has(literal[2]))
));

const catalogKeys = () => new Set(sourceCalls().flatMap((call) => ('key' in call ? [call.key] : [])));

test('a History step label is always a catalog key', () => {
  const keys = catalogKeys();
  const found = sourceFiles()
    .filter((path) => path.endsWith('.js'))
    .flatMap((path) => unkeyedLabels(sourceText(path), keys).map((label) => `${webPath(path)}: ${label}`));
  assert.deepEqual(found, []);
});

test('the check refuses a label built at run time and an unmarked text, and sees the marked ones', () => {
  const keys = new Set(['Change setting']);
  const shapes = [
    'runUndoable(`Reset ${path}`, fn)',
    'history.runUndoable(`${action} comparison`, fn)',
    'withHistory(() => `Fit labels: ${title}`, fn, close)',
    "undoableAction('Not a key', fn, { joinsOpenStep: true })",
    "runUndoable('Reset ' + path, fn)",
    'runUndoable(`Not a key`, fn)',
    "runUndoable(t('Change setting'), fn)",
    "withHistory(() => 'Not a key', fn, close)",
    "runUndoable(flag ? 'Not a key' : msg('Any text'), fn)",
    "runUndoableCommand('Not a key', fn)",
    "withDialogChoice(() => 'Not a key', fn)",
    "scopeChoiceWithHistory(() => 'Not a key', fn, cancel)",
    "restoreLegendItems('Not a key', fn)",
    "commitSpecificRules((current) => current, 'Not a key', { drawing })",
    "askTrackHide(ports, { rows, label: 'Not a key', mode })"
  ];
  for (const source of shapes) assert.equal(unkeyedLabels(source, keys).length, 1, source);
  for (const source of [
    "undoableAction(msg('Any text'), fn)", "runProjectStep('Change setting', fn)", 'runUndoable(label, fn)',
    "withHistory(() => msg('Fit labels'), fn, close)", 'runUndoable(labels[action], fn)',
    "runUndoable(kind === 'stroke' ? msg('Any text') : msg('Other text'), fn)",
    "commitSpecificRules((current) => current, msg('Any text'), { drawing })",
    'commitSpecificRules(build, label, options)',
    "askTrackHide(ports, { rows, label: own[renderer] || msg('Any text') })",
    'askTrackHide(ports, options)'
  ]) assert.deepEqual(unkeyedLabels(source, keys), [], source);
  assert.ok(sourceCalls().filter((call) => 'key' in call && call.key === 'Change setting').length > 0);
});
