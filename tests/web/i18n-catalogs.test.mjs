// The Web UI translation runtime (gbdraw/web/js/utils/i18n.js) and its
// catalogs (gbdraw/web/js/utils/i18n-catalogs/<locale>.js). English is the
// key: the literal first argument of `t()`, the `other` argument of `tn()`,
// and `context + '\u0004' + text` of `tx()`. This test lists those literals in
// index.html and the modules, and fails on a call whose key is not a literal,
// on a catalog key that no call produces, on a translation whose placeholders
// differ from the English, and on a value that is not a string. It prints
// the coverage of each locale without failing.
import assert from 'node:assert/strict';
import { existsSync, readFileSync, readdirSync } from 'node:fs';
import { join, relative, sep } from 'node:path';
import test from 'node:test';
import { pathToFileURL } from 'node:url';

import { CONTEXT_SEPARATOR, SOURCE_LOCALE, setUiLocale, t, tn, tx, uiLocale } from '../../gbdraw/web/js/utils/i18n.js';

const webRoot = join(process.cwd(), 'gbdraw', 'web');
const jsRoot = join(webRoot, 'js');
const catalogRoot = join(jsRoot, 'utils', 'i18n-catalogs');
// The runtime's own source and the catalogs hold no UI text of their own.
const RUNTIME_FILES = new Set(['utils/i18n.js']);

// For each function, the argument positions that are literal keys, and how
// they make the catalog key.
const CALLS = {
  t: { literals: [0], key: ([text]) => text },
  tn: { literals: [1, 2], key: ([, , other]) => other },
  tx: { literals: [0, 1], key: ([context, text]) => context + CONTEXT_SEPARATOR + text }
};
const CALL_START = /(?<![\w$.])(tn|tx|t)\(/g;
const ESCAPES = { n: '\n', t: '\t', r: '\r', b: '\b', f: '\f', v: '\v', 0: '\0' };
const PLACEHOLDER = /\{([A-Za-z_][A-Za-z0-9_]*)\}/g;

/**
 * Reads the string literal that starts at `start` ('…', "…", or `…` without
 * `${`). Returns its value and the index after it, or null for a template
 * literal with a substitution.
 */
const readLiteral = (source, start) => {
  const quote = source[start];
  let value = '';
  for (let at = start + 1; at < source.length; at += 1) {
    const char = source[at];
    if (char === quote) return { value, end: at + 1 };
    if (quote === '`' && char === '$' && source[at + 1] === '{') return null;
    if (char !== '\\') {
      value += char;
      continue;
    }
    const next = source[at + 1];
    if (next === 'u' && source[at + 2] === '{') {
      const close = source.indexOf('}', at);
      value += String.fromCodePoint(parseInt(source.slice(at + 3, close), 16));
      at = close;
    } else if (next === 'u') {
      value += String.fromCharCode(parseInt(source.slice(at + 2, at + 6), 16));
      at += 5;
    } else if (next === 'x') {
      value += String.fromCharCode(parseInt(source.slice(at + 2, at + 4), 16));
      at += 3;
    } else {
      value += ESCAPES[next] ?? next;
      at += 1;
    }
  }
  throw new Error(`unterminated literal at ${start}`);
};

/** Skips one argument expression; returns the index of the `,` or `)` that ends it. */
const skipExpression = (source, start) => {
  let depth = 0;
  for (let at = start; at < source.length; at += 1) {
    const char = source[at];
    if (char === '\'' || char === '"' || char === '`') {
      const literal = readLiteral(source, at);
      if (!literal) {
        // A template literal with substitutions: skip to its closing backtick.
        at = source.indexOf('`', source.indexOf('}', at));
        continue;
      }
      at = literal.end - 1;
    } else if ('([{'.includes(char)) depth += 1;
    else if (')]}'.includes(char)) {
      if (depth === 0) return at;
      depth -= 1;
    } else if (char === ',' && depth === 0) return at;
  }
  throw new Error(`unterminated call at ${start}`);
};

/**
 * The calls of `t`, `tn`, `tx` in `source`: each with its key, or with the
 * problem that keeps the key from being read.
 * @param {string} source
 * @param {string} file the name that problems report
 */
const extractCalls = (source, file) => {
  const calls = [];
  for (const match of source.matchAll(CALL_START)) {
    const name = match[1];
    const line = source.slice(0, match.index).split('\n').length;
    const where = `${file}:${line}`;
    const args = [];
    let at = match.index + match[0].length;
    let problem = '';
    while (args.length <= Math.max(...CALLS[name].literals)) {
      while (/\s/.test(source[at])) at += 1;
      const index = args.length;
      const literal = '\'"`'.includes(source[at]) ? readLiteral(source, at) : null;
      let end = literal ? literal.end : at;
      while (/\s/.test(source[end])) end += 1;
      // A literal that is only part of the argument (`'Go to ' + name`) is no key.
      const whole = literal && (source[end] === ',' || source[end] === ')');
      args.push(whole ? literal.value : null);
      if (!whole) end = skipExpression(source, at);
      if (CALLS[name].literals.includes(index) && args[index] === null) {
        problem = `${where}: argument ${index + 1} of ${name}() is not a string literal`;
        break;
      }
      if (source[end] === ')' && args.length <= Math.max(...CALLS[name].literals)) {
        problem = `${where}: ${name}() has too few arguments`;
        break;
      }
      at = end + 1;
    }
    if (problem) {
      calls.push({ where, problem });
      continue;
    }
    const key = CALLS[name].key(args);
    // An in-DOM template sees an entity decoded, an x-template does not: write the character.
    if (/&(?:[A-Za-z]+|#\d+|#x[0-9A-Fa-f]+);/.test(args.join(''))) {
      calls.push({ where, problem: `${where}: ${name}() key holds an HTML entity; write the character` });
      continue;
    }
    calls.push({ where, key, english: name === 'tx' ? args[1] : key });
  }
  return calls;
};

const sourceFiles = () => [
  join(webRoot, 'index.html'),
  ...readdirSync(jsRoot, { recursive: true })
    .map((path) => String(path))
    .filter((path) => path.endsWith('.js'))
    .filter((path) => !RUNTIME_FILES.has(path.split(sep).join('/')))
    .filter((path) => !path.split(sep).join('/').startsWith('utils/i18n-catalogs/'))
    .map((path) => join(jsRoot, path))
];

const sourceCalls = () => sourceFiles().flatMap((path) => (
  extractCalls(readFileSync(path, 'utf8'), relative(webRoot, path).split(sep).join('/'))
));

const placeholdersOf = (text) => [...new Set([...text.matchAll(PLACEHOLDER)].map((match) => match[1]))].sort();

/**
 * The problems of one catalog against the English keys of the sources.
 * @param {Map<string, string>} english the key of each call and its English text
 * @param {Record<string, unknown>} catalog
 */
const catalogProblems = (english, catalog) => Object.entries(catalog).flatMap(([key, value]) => {
  if (typeof value !== 'string') return [`${JSON.stringify(key)}: the value is not a string`];
  if (!english.has(key)) return [`${JSON.stringify(key)}: no t(), tn() or tx() call has this key`];
  const expected = placeholdersOf(english.get(key));
  const actual = placeholdersOf(value);
  return expected.join() === actual.join() ? []
    : [`${JSON.stringify(key)}: placeholders {${actual.join('}, {')}} differ from the English {${expected.join('}, {')}}`];
});

const catalogFiles = () => (existsSync(catalogRoot)
  ? readdirSync(catalogRoot).filter((name) => name.endsWith('.js')).sort()
  : []);

test('extraction reads the quoting forms of templates and modules', () => {
  const source = [
    '<button :title="t(\'Undo\')" aria-hidden="true">{{ t(\'Don\\\'t save\') }}</button>',
    '{{ t("Say \\"hi\\"") }} {{ t(`Plain`) }} {{ t(\'Rows {from}–{to}\', { from: a, to: b }) }}',
    'const label = tn(rows.filter((row) => row.on).length, \'{count} row\', \'{count} rows\');',
    'tx(\'Linear (size scaling)\', \'Linear\'); foo.t(notAKey); split(\'x\');'
  ].join('\n');
  assert.deepEqual(extractCalls(source, 'sample').map((call) => call.key ?? call.problem), [
    'Undo', 'Don\'t save', 'Say "hi"', 'Plain', 'Rows {from}–{to}', '{count} rows',
    `Linear (size scaling)${CONTEXT_SEPARATOR}Linear`
  ]);
});

test('extraction refuses a key that is not a literal', () => {
  const source = 't(entry.title); t(`Go to ${name}`); t(\'Help: \' + name); tn(count, one, \'{count} rows\'); t(\'A &amp; B\')';
  assert.deepEqual(extractCalls(source, 'sample').map((call) => call.problem), [
    'sample:1: argument 1 of t() is not a string literal',
    'sample:1: argument 1 of t() is not a string literal',
    'sample:1: argument 1 of t() is not a string literal',
    'sample:1: argument 2 of tn() is not a string literal',
    'sample:1: t() key holds an HTML entity; write the character'
  ]);
});

test('every t(), tn() and tx() call in the Web UI has a literal key', () => {
  const problems = sourceCalls().filter((call) => call.problem).map((call) => call.problem);
  assert.deepEqual(problems, []);
});

test('catalog check: stale keys, placeholder sets and non-string values', () => {
  const english = new Map([['Save Session', 'Save Session'], ['{count} rows', '{count} rows']]);
  assert.deepEqual(catalogProblems(english, {
    'Save Session': 'セッションを保存',
    '{count} rows': '{n} 行',
    'Load Session': 'セッションを読み込む',
    'Undo': 1
  }), [
    '"{count} rows": placeholders {n} differ from the English {count}',
    '"Load Session": no t(), tn() or tx() call has this key',
    '"Undo": the value is not a string'
  ]);
});

test('every catalog matches the keys and placeholders of the sources', async () => {
  const english = new Map(sourceCalls().filter((call) => call.key).map((call) => [call.key, call.english]));
  for (const file of catalogFiles()) {
    const { default: catalog } = await import(pathToFileURL(join(catalogRoot, file)).href);
    assert.ok(Object.isFrozen(catalog), `${file} exports a frozen object`);
    assert.deepEqual(catalogProblems(english, catalog), [], file);
    const missing = [...english.keys()].filter((key) => !Object.hasOwn(catalog, key));
    console.log(`${file}: ${english.size - missing.length} of ${english.size} keys translated`);
  }
  console.log(`${english.size} keys in the sources; catalogs: ${catalogFiles().join(', ') || 'none'}`);
});

test('English passes through with its placeholders filled', () => {
  assert.equal(uiLocale.value, SOURCE_LOCALE);
  assert.equal(t('Save Session'), 'Save Session');
  assert.equal(t('Rows {from}–{to}', { from: 3, to: 4 }), 'Rows 3–4');
  assert.equal(t('Help: {name}'), 'Help: {name}', 'a placeholder without a value stays');
  assert.equal(t('{a} and {b}', { a: '{b}', b: 'x' }), '{b} and x', 'a value is not substituted again');
  assert.equal(tn(1, '{count} row', '{count} rows'), '1 row');
  assert.equal(tn(2, '{count} row', '{count} rows'), '2 rows');
  assert.equal(tn(0, 'one record', '{count} records of {file}', { file: 'a.gb' }), '0 records of a.gb');
  assert.equal(tx('Linear (size scaling)', 'Linear'), 'Linear');
});

test('a locale shows its catalog and falls back to English', async () => {
  const catalog = Object.freeze({
    'Save Session': 'セッションを保存',
    'Rows {from}–{to}': '{from}–{to} 行',
    '{count} rows': '{count} 行',
    [`Linear (size scaling)${CONTEXT_SEPARATOR}Linear`]: '線形（サイズ比例）'
  });
  const loaded = [];
  try {
    await setUiLocale('ja', async (locale) => {
      loaded.push(locale);
      return { default: catalog };
    });
    assert.deepEqual(loaded, ['ja']);
    assert.equal(uiLocale.value, 'ja');
    assert.equal(t('Save Session'), 'セッションを保存');
    assert.equal(t('Rows {from}–{to}', { from: 3, to: 4 }), '3–4 行');
    assert.equal(tn(1, '{count} row', '{count} rows'), '1 行');
    assert.equal(tx('Linear (size scaling)', 'Linear'), '線形（サイズ比例）');
    assert.equal(t('Linear'), 'Linear', 'tx context keys are not t keys');
    assert.equal(t('Load Session'), 'Load Session', 'a missing translation shows the English');
    assert.equal(tn(1, 'one record', '{count} records'), 'one record');
    await assert.rejects(setUiLocale('../x'), /Unknown UI locale/);
    assert.equal(uiLocale.value, 'ja', 'a refused locale changes nothing');
  } finally {
    await setUiLocale(SOURCE_LOCALE);
  }
  assert.equal(t('Save Session'), 'Save Session');
});
