// @ts-check
// The Undo and Redo titles translate a step's label where it is shown
// (`undoHistoryTitle` in app/app-setup.js), by the label as the key: a label the
// catalogs do not hold shows English in every language (OV-2415).
import assert from 'node:assert/strict';
import test from 'node:test';

import { skipExpression, sourceCalls, sourceFiles, sourceText, webPath } from '../../tools/i18n/source-keys.mjs';

// The History calls whose first argument names a step.
const LABEL_TAKERS = /(?<![\w$])(?:runUndoable|runUndoableCheckpoint|runUndoableArtifactReplacement|runProjectStep|begin|beginCheckpoint|beginArtifactReplacement|beginHistoryTransaction|undoableAction|editEditorIntent|withHistory)\(\s*/g;
const PLAIN_LITERAL = /^(['"])(?:\\.|(?!\1)[^\\])*\1$/;

/**
 * The label arguments of `source` that no catalog holds: a template literal with
 * a substitution, a text joined at run time, or a plain literal that is no key.
 * A `msg()` literal, and a name that forwards a label from another call, pass.
 * @param {string} source @param {ReadonlySet<string>} keys
 * @returns {string[]}
 */
const unkeyedLabels = (source, keys) => [...source.matchAll(LABEL_TAKERS)].flatMap((match) => {
  const start = match.index + match[0].length;
  const argument = source.slice(start, skipExpression(source, start)).trim();
  if (argument.includes('${') || /['"`]\s*\+|\+\s*['"`]/.test(argument)) return [argument];
  return PLAIN_LITERAL.test(argument) && !keys.has(argument.slice(1, -1)) ? [argument] : [];
});

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
    "runUndoable('Reset ' + path, fn)"
  ];
  for (const source of shapes) assert.equal(unkeyedLabels(source, keys).length, 1, source);
  for (const source of [
    "undoableAction(msg('Any text'), fn)", "runProjectStep('Change setting', fn)", 'runUndoable(label, fn)',
    "withHistory(() => msg('Fit labels'), fn, close)", 'runUndoable(labels[action], fn)'
  ]) assert.deepEqual(unkeyedLabels(source, keys), [], source);
  assert.ok(sourceCalls().filter((call) => 'key' in call && call.key === 'Change setting').length > 0);
});
