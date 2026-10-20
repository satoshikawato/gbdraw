// The Web UI translation runtime (gbdraw/web/js/utils/i18n.js) and its
// catalogs (gbdraw/web/js/utils/i18n-catalogs/<locale>.js). English is the
// key: the literal first argument of `t()` and of `msg()` (a text kept until
// it is shown, which `tm()` translates there), the `other` argument of `tn()`,
// `context + '\u0004' + text` of `tx()`, and the table values of the error
// wording owner (tests/web/helpers/error-text-keys.mjs). This test lists those
// literals in index.html and the modules, and fails on a call whose key is not
// a literal,
// on a `tm()` given a literal (that is a `t()` call),
// on a catalog key that no call produces, on a translation whose placeholders
// differ from the English, and on a value that is not a string. It prints
// the coverage of each locale without failing.
import assert from 'node:assert/strict';
import { existsSync, readFileSync, readdirSync } from 'node:fs';
import { join, relative, sep } from 'node:path';
import test from 'node:test';
import { pathToFileURL } from 'node:url';

import {
  CONTEXT_SEPARATOR, SOURCE_LOCALE, UI_LOCALES, joinSentences, listText, msg, sentenceSeparator, setUiLocale, t, textParts, tm, tn, tx, uiLocale, uiLocaleOf
} from '../../gbdraw/web/js/utils/i18n.js';
import { errorTextKeys } from './helpers/error-text-keys.mjs';

const webRoot = join(process.cwd(), 'gbdraw', 'web');
const jsRoot = join(webRoot, 'js');
const catalogRoot = join(jsRoot, 'utils', 'i18n-catalogs');
// The runtime's own source and the catalogs hold no UI text of their own.
const RUNTIME_FILES = new Set(['utils/i18n.js']);

// For each function, the argument positions that are literal keys, and how
// they make the catalog key.
const CALLS = {
  t: { literals: [0], key: ([text]) => text },
  msg: { literals: [0], key: ([text]) => text },
  tn: { literals: [1, 2], key: ([, , other]) => other },
  tx: { literals: [0, 1], key: ([context, text]) => context + CONTEXT_SEPARATOR + text }
};
const CALL_START = /(?<![\w$.])(msg|tn|tx|t)\(/g;
// `tm()` translates a value that `msg()` marked; a literal there is a `t()` call.
const TM_LITERAL = /(?<![\w$.])tm\(\s*['"`]/g;
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
  for (const match of source.matchAll(TM_LITERAL)) {
    const where = `${file}:${source.slice(0, match.index).split('\n').length}`;
    calls.push({ where, problem: `${where}: tm() takes a value that msg() marked; write t() for a literal` });
  }
  return calls;
};

const NAMED_ENTITIES = { quot: '"', amp: '&', lt: '<', gt: '>', apos: '\'', '#39': '\'', nbsp: ' ' };
const ENTITY = /&(#x[0-9A-Fa-f]+|#\d+|[A-Za-z][A-Za-z0-9]*);/g;
const X_TEMPLATE = /<script\b[^>]*type="text\/x-template"[^>]*>[\s\S]*?<\/script>/g;

/**
 * Decodes the entities of the list above and numeric references the way the
 * browser does; any other named entity stays as written, so the check below
 * can refuse it.
 */
const decodeEntities = (text) => text.replace(ENTITY, (whole, name) => {
  if (name[0] === '#' && name[1] !== 'x' && name[1] !== 'X') return String.fromCodePoint(Number(name.slice(1)));
  if (name[0] === '#' && /^#[xX]/.test(name)) return String.fromCodePoint(parseInt(name.slice(2), 16));
  return NAMED_ENTITIES[name] ?? whole;
});

/**
 * index.html as Vue compiles it: Vue decodes entities itself, in attributes
 * and in {{ }}, in the in-DOM root and in an x-template alike.
 */
const compiledHtml = (html) => decodeEntities(html);

const VOID_TAGS = new Set(['area', 'base', 'br', 'col', 'embed', 'hr', 'img', 'input', 'link', 'meta', 'source', 'track', 'wbr']);
const TRANSLATED = /(?<![\w$.])(?:tn|tx|t)\(|<i18n-text\b/;

/**
 * The elements with a `v-memo` whose subtree shows UI text (a t(), tn() or
 * tx() call, or an <i18n-text>) and whose dependency list lacks
 * `uiLocale.value`: a locale switch would not update the memoized subtree.
 * @returns {string[]} one `file:line: ...` problem per element
 */
const vMemoProblems = (source, file) => {
  const problems = [];
  for (const start of source.matchAll(/<([A-Za-z][\w-]*)\b(?:"[^"]*"|'[^']*'|[^'">])*?\sv-memo=(?:"([^"]*)"|'([^']*)')(?:"[^"]*"|'[^']*'|[^'">])*>/g)) {
    const tag = start[1].toLowerCase();
    const deps = start[2] ?? start[3];
    // Find the matching end tag by counting the same tag name.
    const tags = new RegExp(`<(/?)${tag}\\b(?:"[^"]*"|'[^']*'|[^'">])*?(/?)>`, 'gi');
    tags.lastIndex = start.index + start[0].length;
    let depth = 1;
    let end = source.length;
    for (let next = tags.exec(source); next; next = tags.exec(source)) {
      if (next[1]) depth -= 1;
      else if (!next[2] && !VOID_TAGS.has(tag)) depth += 1;
      if (depth === 0) { end = next.index; break; }
    }
    const body = source.slice(start.index, end);
    if (TRANSLATED.test(body) && !deps.includes('uiLocale.value')) {
      problems.push(`${file}:${source.slice(0, start.index).split('\n').length}: v-memo="${deps}" shows UI text but does not list uiLocale.value`);
    }
  }
  return problems;
};

/**
 * The start tags of `html` that break on a raw `"` inside a double-quoted
 * attribute value and that hold a t(), tn() or tx() call. Such a quote ends
 * the value early (write `&quot;` in index.html outside x-templates).
 * @returns {string[]} one `file:line: ...` problem per tag
 */
const rawQuoteProblems = (html, file) => {
  const problems = [];
  const tags = /<[A-Za-z][^\s/>]*(?:"[^"]*"|'[^']*'|[^'">])*>?/g;
  for (const tag of html.matchAll(tags)) {
    const text = tag[0];
    if (!/(?<![\w$.])(?:tn|tx|t)\(/.test(text)) continue;
    // Walk the attributes the way an HTML parser does.
    let at = text.search(/\s/);
    let broken = false;
    while (at >= 0 && at < text.length && !broken) {
      while (/\s/.test(text[at] ?? '')) at += 1;
      if (at >= text.length || text[at] === '>' || text[at] === '/') break;
      const name = /^[^\s=>/]+/.exec(text.slice(at))?.[0] ?? '';
      if (!name || /["'()<]/.test(name)) { broken = true; break; }
      at += name.length;
      while (/\s/.test(text[at] ?? '')) at += 1;
      if (text[at] !== '=') continue;
      at += 1;
      while (/\s/.test(text[at] ?? '')) at += 1;
      const quote = text[at];
      if (quote === '"' || quote === '\'') {
        const close = text.indexOf(quote, at + 1);
        if (close < 0) { broken = true; break; }
        at = close + 1;
        if (at < text.length && !/[\s>/]/.test(text[at])) broken = true;
      } else at += (/^[^\s>]*/.exec(text.slice(at))?.[0] ?? '').length;
    }
    if (broken) {
      const line = html.slice(0, tag.index).split('\n').length;
      problems.push(`${file}:${line}: a double-quoted attribute holds a raw "; write &quot; (${text.slice(0, 60).replace(/\s+/g, ' ')}...)`);
    }
  }
  return problems;
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

const sourceText = (path) => {
  const source = readFileSync(path, 'utf8');
  return path.endsWith('.html') ? compiledHtml(source) : source;
};

const sourceCalls = () => sourceFiles().flatMap((path) => (
  extractCalls(sourceText(path), relative(webRoot, path).split(sep).join('/'))
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
  const source = 't(entry.title); t(`Go to ${name}`); t(\'Help: \' + name); tn(count, one, \'{count} rows\'); t(\'A &amp; B\'); msg(title)';
  assert.deepEqual(extractCalls(source, 'sample').map((call) => call.problem), [
    'sample:1: argument 1 of t() is not a string literal',
    'sample:1: argument 1 of t() is not a string literal',
    'sample:1: argument 1 of t() is not a string literal',
    'sample:1: argument 2 of tn() is not a string literal',
    'sample:1: t() key holds an HTML entity; write the character',
    'sample:1: argument 1 of msg() is not a string literal'
  ]);
});

test('msg() literals are keys; tm() reads a marked value and refuses a literal', () => {
  const source = [
    'const TITLES = { canvas: msg(\'Canvas\'), font: msg("Font") };',
    '{{ tm(section.title) }} {{ tm(note.text, { n: 2 }) }} item.msg(\'not a key\'); notmsg(\'x\');',
    'const shown = tm(\'Canvas\') + tm(`Font`);'
  ].join('\n');
  assert.deepEqual(extractCalls(source, 'sample').map((call) => call.key ?? call.problem), [
    'Canvas', 'Font',
    'sample:3: tm() takes a value that msg() marked; write t() for a literal',
    'sample:3: tm() takes a value that msg() marked; write t() for a literal'
  ]);
});

test('index.html keys are read after Vue decodes entities, in the root and in x-templates', () => {
  const html = [
    '<div :title="t(\'Change &quot;{from}&quot; to {to} &amp; &#x41; &#66;\')"></div>',
    '<script type="text/x-template" id="a">{{ t(\'Keep &quot;{name}&quot; &amp; "{other}"\') }}</script>',
    '<p>{{ t(\'Fish &amp; chips\') }} {{ t(\'Wait&nbsp;now &copy;\') }}</p>'
  ].join('\n');
  assert.deepEqual(extractCalls(compiledHtml(html), 'sample').map((call) => call.key ?? call.problem), [
    'Change "{from}" to {to} & A B',
    'Keep "{name}" & "{other}"',
    'Fish & chips',
    'sample:3: t() key holds an HTML entity; write the character'
  ]);
  // A module's key is read as written: an entity there is refused.
  assert.deepEqual(extractCalls('t(\'A &amp; B\')', 'module.js').map((call) => call.problem), [
    'module.js:1: t() key holds an HTML entity; write the character'
  ]);
});

test('a v-memo that shows UI text lists uiLocale.value in its dependencies', () => {
  const html = [
    '<li v-for="r in rows" v-memo="[r.id]"><b>{{ t(\'Edited\') }}</b><i>x</i></li>',
    '<li v-for="r in rows" v-memo="[r.id, uiLocale.value]"><b>{{ t(\'Edited\') }}</b></li>',
    '<li v-for="r in rows" v-memo="[r.id]"><b>{{ r.name }}</b></li>',
    '<ul v-memo="[x]"><li><i18n-text :text="t(\'a {b}\')"></i18n-text></li></ul><p>{{ t(\'Outside\') }}</p>',
    "<script>const template = `<div v-memo='[a]'><span>${'{{ t(\\'In JS\\') }}'}</span></div>`;</script>"
  ].join('\n');
  assert.deepEqual(vMemoProblems(html, 'sample').map((problem) => problem.split(' shows')[0]), [
    'sample:1: v-memo="[r.id]"',
    'sample:4: v-memo="[x]"',
    'sample:5: v-memo="[a]"'
  ]);
});

/**
 * Vue drops a whitespace-only text node that is the first or last child of an
 * element, so a space written before a leading `{{` or after a trailing `}}`
 * would vanish from the page (the English it replaced showed it). Write
 * `{{ ' ' }}` there.
 * @returns {string[]} one `file:line: ...` problem per place
 */
const edgeSpaceProblems = (html, file) => {
  const problems = [];
  const edges = /<([A-Za-z][\w-]*)(?:"[^"]*"|'[^']*'|[^'">])*(?<!\/)>[ \t]+\{\{|\}\}[ \t]+<\//g;
  for (const match of html.matchAll(edges)) {
    if (match[1] && VOID_TAGS.has(match[1].toLowerCase())) continue;
    problems.push(`${file}:${html.slice(0, match.index).split('\n').length}: a space at the start or end of an element beside {{ }} is dropped; write {{ ' ' }}`);
  }
  return problems;
};

test('a space at the edge of an element beside a mustache is written as {{ \' \' }}', () => {
  const html = '<span> {{ t(\'A\') }}</span><i></i> {{ t(\'B\') }}<input> {{ t(\'C\') }}<b>{{ t(\'D\') }} </b>';
  assert.deepEqual(edgeSpaceProblems(html, 'sample').map((problem) => problem.slice(0, 9)), ['sample:1:', 'sample:1:']);
  assert.deepEqual(edgeSpaceProblems(readFileSync(join(webRoot, 'index.html'), 'utf8'), 'index.html'), []);
});

test('every v-memo in the Web UI that shows UI text lists uiLocale.value', () => {
  const problems = sourceFiles().flatMap((path) => vMemoProblems(readFileSync(path, 'utf8'), relative(webRoot, path).split(sep).join('/')));
  assert.deepEqual(problems, []);
});

test('a raw double quote inside a double-quoted attribute with a translation call is refused', () => {
  const html = [
    '<i18n-text :text="t(\'Change &quot;{from}&quot; to\')"></i18n-text>',
    '<i18n-text :text="t(\'Change "{from}" to\')"></i18n-text>',
    '<span :title="t(\'Say " now\')">x</span>',
    '<span title="plain" :aria-label="t(\'Fine\')">{{ t(\'Say "hi"\') }}</span>',
    '<span title="Say "hi"">no call</span>'
  ].join('\n');
  assert.deepEqual(rawQuoteProblems(html, 'sample').map((problem) => problem.slice(0, 56)), [
    'sample:2: a double-quoted attribute holds a raw "; write',
    'sample:3: a double-quoted attribute holds a raw "; write'
  ]);
});

test('index.html has no raw double quote inside a double-quoted translated attribute', () => {
  assert.deepEqual(rawQuoteProblems(readFileSync(join(webRoot, 'index.html'), 'utf8'), 'index.html'), []);
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

// A module that marks keys with msg() imports it from utils/i18n.js, or, when
// a test loads it without utils/i18n.js (copied alone or with a fixed set of
// its imports), declares the same identity.
const LOCAL_MARKER = 'const msg = (text) => text;';
const MSG_IMPORT = /^import \{[^}]*\bmsg\b[^}]*\} from '(?:\.\.\/|\.\/)+utils\/i18n\.js';$/m;
const markerProblems = (source, file) => {
  if (!/(?<![\w$.])msg\(/.test(source)) return [];
  const local = source.split('\n').filter((line) => /^\s*(?:const|let|var|function)\s+msg\b/.test(line));
  if (MSG_IMPORT.test(source) && !local.length) return [];
  if (!MSG_IMPORT.test(source) && local.length === 1 && local[0].trim() === LOCAL_MARKER) return [];
  return [`${file}: msg() comes from utils/i18n.js, or, in a module a test loads without it, from \`${LOCAL_MARKER}\``];
};

test('msg() is the marker of utils/i18n.js or its identity where a test loads the module without it', () => {
  assert.deepEqual(markerProblems('import { msg } from \'../utils/i18n.js\';\nconst a = msg(\'A\');', 'a.js'), []);
  assert.deepEqual(markerProblems('const msg = (text) => text;\nconst a = msg(\'A\');', 'b.js'), []);
  assert.deepEqual(markerProblems('import x from \'./x.js\';\nconst msg = (text) => text;\nmsg(\'A\');', 'c.js'), []);
  assert.equal(markerProblems('import { msg } from \'../utils/i18n.js\';\nconst msg = (text) => text;\nmsg(\'A\');', 'e.js').length, 1);
  assert.equal(markerProblems('const msg = (text) => text.trim();\nmsg(\'A\');', 'd.js').length, 1);
  assert.deepEqual(sourceFiles().filter((path) => path.endsWith('.js'))
    .flatMap((path) => markerProblems(readFileSync(path, 'utf8'), relative(webRoot, path).split(sep).join('/'))), []);
});

test('every catalog matches the keys and placeholders of the sources', async () => {
  const english = new Map(sourceCalls().filter((call) => call.key).map((call) => [call.key, call.english]));
  // The error text keys are the values of the wording owner's tables (D-17), read without msg() markers.
  for (const key of await errorTextKeys()) english.set(key, key);
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
  assert.equal(msg('Canvas'), 'Canvas', 'msg() keeps the English, so comparisons keep working');
  assert.equal(tm(msg('Canvas')), 'Canvas');
  assert.equal(tm('Step {n}', { n: 2 }), 'Step 2');
  assert.equal(tm('a.gbk'), 'a.gbk', 'a value that is no key shows as written');
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
    const title = msg('Save Session');
    assert.equal(title, 'Save Session', 'a marked value stays English in every locale');
    assert.equal(tm(title), 'セッションを保存', 'tm() translates a marked value where it is shown');
    assert.equal(tm('Rows {from}–{to}', { from: 1, to: 2 }), '1–2 行');
    assert.equal(tm('a.gbk'), 'a.gbk');
    await assert.rejects(setUiLocale('../x'), /Unknown UI locale/);
    assert.equal(uiLocale.value, 'ja', 'a refused locale changes nothing');
  } finally {
    await setUiLocale(SOURCE_LOCALE);
  }
  assert.equal(t('Save Session'), 'Save Session');
});

test('a sentence with elements splits at the placeholders that have slots', () => {
  assert.deepEqual(textParts('Replaces your {type} default color {color}.', ['color'], { type: 'CDS' }), [
    { text: 'Replaces your CDS default color ' }, { slot: 'color' }, { text: '.' }
  ]);
  assert.deepEqual(textParts('{a} and {b}', ['a', 'b']), [{ slot: 'a' }, { text: ' and ' }, { slot: 'b' }]);
  assert.deepEqual(textParts('Line up with {reference}', ['reference']), [{ text: 'Line up with ' }, { slot: 'reference' }]);
  assert.deepEqual(textParts('Keep {missing} and {color}', ['color']), [{ text: 'Keep {missing} and ' }, { slot: 'color' }],
    'a placeholder without a slot or a value stays as written');
  assert.deepEqual(textParts('{type} is {color}', ['color'], { type: '{color}' }), [{ text: '{color} is ' }, { slot: 'color' }],
    'a value is filled after the split, so it never becomes a slot');
  assert.deepEqual(textParts('No placeholders', ['color']), [{ text: 'No placeholders' }]);
});

test('when setUiLocale calls overlap, only the latest one applies', async () => {
  /** A loader whose catalogs arrive when the test releases them. */
  const deferredLoader = () => {
    const pending = new Map();
    return {
      load: (locale) => new Promise((resolve) => pending.set(locale, () => resolve({ default: { 'Save Session': `${locale}: save` } }))),
      release: (locale) => pending.get(locale)()
    };
  };
  try {
    // The earlier call's catalog arrives last.
    let loader = deferredLoader();
    const first = setUiLocale('ja', loader.load);
    const second = setUiLocale('zh-Hans', loader.load);
    loader.release('zh-Hans');
    assert.equal(await second, true);
    loader.release('ja');
    assert.equal(await first, false);
    assert.equal(uiLocale.value, 'zh-Hans');
    assert.equal(t('Save Session'), 'zh-Hans: save');
    // The earlier call's catalog arrives first; English as the latest call needs no catalog.
    loader = deferredLoader();
    const third = setUiLocale('ja', loader.load);
    const fourth = setUiLocale(SOURCE_LOCALE, loader.load);
    assert.equal(await fourth, true);
    loader.release('ja');
    assert.equal(await third, false);
    assert.equal(uiLocale.value, SOURCE_LOCALE);
    assert.equal(t('Save Session'), 'Save Session');
  } finally {
    await setUiLocale(SOURCE_LOCALE);
  }
});

test('the hosted analytics notice of the About section is in the key list', () => {
  const keys = new Set(sourceCalls().map((call) => call.key));
  for (const key of [
    'Hosted Site Analytics',
    'The hosted {site} deployment uses Google Analytics 4 for aggregate page-usage metrics.',
    'Uploaded genome files and generated diagrams are still processed locally in your browser and are not sent to Google Analytics by gbdraw.'
  ]) assert.ok(keys.has(key), key);
});

test('the UI languages are English, Japanese, Simplified and Traditional Chinese, each named in itself', () => {
  assert.deepEqual(UI_LOCALES.map(({ tag, name }) => [tag, name]), [
    ['en', 'English'], ['ja', '日本語'], ['zh-Hans', '简体中文'], ['zh-Hant', '繁體中文']
  ]);
  assert.ok(Object.isFrozen(UI_LOCALES) && UI_LOCALES.every(Object.isFrozen));
  for (const { tag } of UI_LOCALES) assert.equal(uiLocaleOf(tag), tag);
  for (const value of ['', 'EN', 'ja-JP', 'zh', 'zh-CN', 'fr', null, undefined, 3, '../x']) {
    assert.equal(uiLocaleOf(value), null, JSON.stringify(value));
  }
});

test('a locale outside the UI languages is refused', async () => {
  await assert.rejects(setUiLocale('fr', async () => ({ default: {} })), /Unknown UI locale/);
  assert.equal(uiLocale.value, SOURCE_LOCALE);
});

test('two sentences are joined with a space in English and with nothing in Japanese and Chinese', async () => {
  assert.equal(sentenceSeparator(), ' ');
  assert.equal(joinSentences('Exported a.svg.', 'Find it in your browser\'s downloads.'), 'Exported a.svg. Find it in your browser\'s downloads.');
  try {
    for (const tag of ['ja', 'zh-Hans', 'zh-Hant']) {
      await setUiLocale(tag, async () => ({ default: {} }));
      assert.equal(sentenceSeparator(), '', tag);
      assert.equal(joinSentences('保存しました。', 'ダウンロードを確認してください。'), '保存しました。ダウンロードを確認してください。', tag);
    }
  } finally {
    await setUiLocale(SOURCE_LOCALE);
  }
});

test('a list keeps today\'s English and uses the locale\'s list pattern in Japanese and Chinese', async () => {
  assert.equal(listText([]), '');
  assert.equal(listText(['A']), 'A');
  assert.equal(listText(['A', 'B']), 'A and B');
  assert.equal(listText(['A', 'B', 'C']), 'A, B, and C');
  assert.equal(listText(['A', 'B', 'C', '2 more']), 'A, B, C, and 2 more');
  try {
    for (const tag of ['ja', 'zh-Hans', 'zh-Hant']) {
      await setUiLocale(tag, async () => ({ default: {} }));
      const expected = new Intl.ListFormat(tag, { type: 'conjunction' });
      assert.equal(listText(['A']), 'A', tag);
      assert.equal(listText(['A', 'B']), expected.format(['A', 'B']), tag);
      assert.equal(listText(['A', 'B', 'C']), expected.format(['A', 'B', 'C']), tag);
      assert.doesNotMatch(listText(['A', 'B', 'C']), /\band\b/, tag);
    }
  } finally {
    await setUiLocale(SOURCE_LOCALE);
  }
});

test('the shown locale sets the page language, and only a switch that applies does', async () => {
  const root = { lang: 'en', setAttribute(name, value) { assert.equal(name, 'lang'); this.lang = value; } };
  const previous = globalThis.document;
  globalThis.document = { documentElement: root };
  try {
    await setUiLocale('zh-Hant', async () => ({ default: {} }));
    assert.equal(root.lang, 'zh-Hant');
    await assert.rejects(setUiLocale('ja', async () => { throw new Error('offline'); }), /offline/);
    assert.equal(root.lang, 'zh-Hant', 'a catalog that fails to load changes nothing');
    assert.equal(uiLocale.value, 'zh-Hant');
    await setUiLocale(SOURCE_LOCALE);
    assert.equal(root.lang, 'en');
  } finally {
    await setUiLocale(SOURCE_LOCALE);
    globalThis.document = previous;
  }
});
