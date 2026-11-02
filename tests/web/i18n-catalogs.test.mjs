// The Web UI translation runtime (gbdraw/web/js/utils/i18n.js) and its
// catalogs (gbdraw/web/js/utils/i18n-catalogs/<locale>.js). The keys of the
// sources are read by tools/i18n/source-keys.mjs, which lists what a key is.
// This test fails on a call whose key is not a literal, on a `tm()` given a
// literal (that is a `t()` call), on a `tmx()` given a literal text (that is
// a `tx()` call), on a catalog key that no call produces, on a translation
// whose placeholders differ from the English, and on a value that is not a
// string, or, for a form of a `tn()` pair only, an object of plural forms
// (`PluralForms` in utils/i18n.js). It prints the coverage of each locale without failing: a new key
// ships in English until the batch translation before a release, and the
// release tier requires every key (`node tools/i18n/sync.mjs status
// --require-complete`, .claude/skills/translate-web-ui/SKILL.md).
import assert from 'node:assert/strict';
import { readFileSync, readdirSync } from 'node:fs';
import { join, sep } from 'node:path';
import test from 'node:test';
import { pathToFileURL } from 'node:url';

import {
  CONTEXT_SEPARATOR, SOURCE_LOCALE, UI_LOCALES, joinSentences, listText, msg, pluralTranslation, sentenceSeparator, setUiLocale, t, textParts, tm, tn, tx, uiLocale,
  uiLocaleOf
} from '../../gbdraw/web/js/utils/i18n.js';
import { COLLINEAR_COLOR_MODE_TEXT, memberValueLabelKeys } from '../../gbdraw/web/js/app/pairwise-match-popup.js';
import { COLLINEARITY_COLOR_MODES } from '../../tools/i18n/collinearity-words.generated.mjs';
import { detectedTopologyNote } from '../../gbdraw/web/js/services/record-display-model.js';
import { CATALOG_LOCALES, CATALOG_ROOT, catalogProblems, catalogStatus } from '../../tools/i18n/catalogs.mjs';
import {
  PLACEHOLDER, VOID_TAGS, WEB_ROOT as webRoot, compiledHtml, elementEnd, extractCalls, sourceCalls, sourceFiles, sourceKeys, webPath
} from '../../tools/i18n/source-keys.mjs';

const jsRoot = join(webRoot, 'js');
const catalogFiles = () => CATALOG_LOCALES.map((locale) => `${locale}.js`).sort();

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
    const body = source.slice(start.index, elementEnd(source, tag, start.index + start[0].length));
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

test('extraction reads the quoting forms of templates and modules', () => {
  const source = [
    '<button :title="t(\'Undo\')" aria-hidden="true">{{ t(\'Don\\\'t save\') }}</button>',
    '{{ t("Say \\"hi\\"") }} {{ t(`Plain`) }} {{ t(\'Rows {from}–{to}\', { from: a, to: b }) }}',
    'const label = tn(rows.filter((row) => row.on).length, \'{count} row\', \'{count} rows\');',
    'tx(\'Linear (size scaling)\', \'Linear\'); foo.t(notAKey); split(\'x\');'
  ].join('\n');
  assert.deepEqual(extractCalls(source, 'sample').map((call) => call.key ?? call.problem), [
    'Undo', 'Don\'t save', 'Say "hi"', 'Plain', 'Rows {from}–{to}', '{count} row', '{count} rows',
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

test('tmx() reads a derived text in a context and refuses a literal text', () => {
  const source = 'tmx(context, english); tmx(\'member role\', label); tmx(\'member role\', \'anchor\'); tmx("a\\"b", `x`);';
  assert.deepEqual(extractCalls(source, 'sample').map((call) => call.key ?? call.problem), [
    'sample:1: tmx() takes a derived text; write tx() for a literal',
    'sample:1: tmx() takes a derived text; write tx() for a literal'
  ]);
});

test('OV-1390: the member value labels are the generated Python values, each under its field context', () => {
  assert.deepEqual(memberValueLabelKeys().map(({ context, english }) => `${context}: ${english}`), [
    'member role: anchor', 'member role: coortholog', 'member role: inparalog', 'member role: low confidence',
    'member role: local paralog', 'member confidence: high', 'member confidence: medium', 'member confidence: low'
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
  const problems = sourceFiles().flatMap((path) => vMemoProblems(readFileSync(path, 'utf8'), webPath(path)));
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

// OV-1535: an aria-label, title, placeholder, alt or label is UI text that t()
// translates, so a selector comparing its value matches only in English. A
// stylesheet or a module styles and finds an element by a class or data
// attribute that its template owns. A presence test such as [title] is fine.
const TRANSLATED_ATTRIBUTE_SELECTOR = /\[\s*(?:aria-label|aria-description|aria-roledescription|aria-valuetext|title|placeholder|alt|label)\s*[~|^$*]?=/;
/** @param {string} source @param {string} file @returns {string[]} one `file:line: selector` problem per line */
const translatedSelectorProblems = (source, file) => source.split('\n').flatMap((line, index) => (
  TRANSLATED_ATTRIBUTE_SELECTOR.test(line) ? [`${file}:${index + 1}: ${line.trim()}`] : []
));
const webStylesheets = () => readdirSync(webRoot, { recursive: true }).map((path) => String(path).split(sep).join('/'))
  .filter((path) => path.endsWith('.css') && !path.startsWith('vendor/')).map((path) => join(webRoot, path));

test('a selector that compares a translated attribute is refused; a class, data attribute or presence test is not', () => {
  const css = [
    '.result-pane > [aria-label="Result Preview"] { flex: none; }',
    ".preview-controls [aria-label='Reset zoom'] { width: 3rem; }",
    'input[placeholder^="Search"], img[ alt *= "logo" ], [title|=x] {}',
    'abbr:where([title]) {} button[aria-pressed="true"] {} [data-result-preview] {} .preview-zoom-reset {}',
    "document.querySelector('[aria-label]').getAttribute('aria-label');"
  ].join('\n');
  assert.deepEqual(translatedSelectorProblems(css, 'sample').map((problem) => problem.split(':').slice(0, 2).join(':')), [
    'sample:1', 'sample:2', 'sample:3'
  ]);
});

test('no stylesheet, template or module of the Web UI selects by a translated attribute value', () => {
  const files = [...sourceFiles(), ...webStylesheets()];
  assert.ok(files.some((path) => path.endsWith(join('css', 'tailwind.input.css'))), 'the Tailwind input is read');
  const problems = files.flatMap((path) => translatedSelectorProblems(readFileSync(path, 'utf8'), webPath(path)));
  assert.deepEqual(problems, []);
});

// A static placeholder shows in every language as written, so it holds only a
// code value (slot_id, product,gene,locus_tag); an example sentence goes
// through t('e.g., {example}', …) like its siblings (OV-2329: "e.g. K-12").
test('index.html has no static placeholder that is an English example', () => {
  const html = readFileSync(join(webRoot, 'index.html'), 'utf8');
  const problems = html.split('\n').flatMap((line, index) => (
    /(?:^|\s)placeholder="(?:e\.g\.|for example)/i.test(line) ? [`index.html:${index + 1}: ${line.trim()}`] : []
  ));
  assert.deepEqual(problems, []);
});

test('every t(), tn() and tx() call in the Web UI has a literal key', () => {
  const problems = sourceCalls().filter((call) => call.problem).map((call) => call.problem);
  assert.deepEqual(problems, []);
});

/**
 * Source keys for the catalog check, each its own English: `[key]` for a
 * plain key, `[one, other]` for the two forms of a tn() pair.
 * @param {...string[]} groups
 */
const keysOf = (...groups) => new Map(groups.flatMap((group) => {
  const plural = group.length === 2 ? { one: group[0], other: group[1] } : null;
  return group.map((key) => [key, { key, english: key, plural }]);
}));

test('catalog check: placeholder sets and non-string or empty values; a stale key is no problem', () => {
  const english = keysOf(['Save Session'], ['{count} rows']);
  assert.deepEqual(catalogProblems(english, {
    'Save Session': 'セッションを保存',
    '{count} rows': '{n} 行',
    'Load Session': 'セッションを読み込む',
    'Undo': 1,
    'Empty': 1
  }), [
    '"{count} rows": placeholders {n} differ from the English {count}'
  ]);
  const withUndo = new Map([...english, ...keysOf(['Undo'], ['Blank'])]);
  assert.deepEqual(catalogProblems(withUndo, { Undo: 1, Blank: ' \n' }), [
    '"Undo": the value is not a string',
    '"Blank": the translation is empty'
  ]);
  assert.deepEqual(catalogStatus(english.keys(), { 'Save Session': 'x', 'Load Session': 'y' }).stale, ['Load Session']);
});

test('catalog check: plural forms only on the `other` form of a tn() pair, with `other`, each form checked like a string', () => {
  const english = keysOf(['{count} row', '{count} rows'], ['Save Session'], ['{count} file', '{count} files'], ['{count} item', '{count} items']);
  assert.deepEqual(catalogProblems(english, {
    '{count} row': '{count} запись',
    '{count} rows': { one: '{count} запись', few: '{count} записи', many: '{count} записей', other: '{count} записи' }
  }), [], 'an object of plural forms on the `other` form');
  assert.deepEqual(catalogProblems(english, {
    'Save Session': { one: 'Sitzung speichern', other: 'Sitzungen speichern' },
    '{count} row': { one: '{count} запись', other: '{count} записи' },
    '{count} rows': { few: '{count} записи', many: '{count} записей' },
    '{count} items': { other: '{count} штуки', dual: '{count}', 0: 'x' },
    '{count} files': { one: 7, few: '{n} файла', many: '  ', other: '{count} файла' }
  }), [
    '"Save Session": plural forms are allowed only for the other form of a tn() pair',
    '"{count} row": plural forms are allowed only for the other form of a tn() pair',
    '"{count} rows": the plural forms lack "other"',
    '"{count} items": 0, dual are not plural categories',
    '"{count} files" one: the value is not a string',
    '"{count} files" few: placeholders {n} differ from the English {count}',
    '"{count} files" many: the translation is empty'
  ]);
  assert.deepEqual(catalogProblems(english, { '{count} rows': ['{count} записи'], '{count} row': null }), [
    '"{count} rows": the value is not a string',
    '"{count} row": the value is not a string'
  ], 'an array or null is no object of plural forms');
});

// A module that marks keys with msg() imports it from utils/i18n.js, or, when
// a test loads it without utils/i18n.js (copied alone or with a fixed set of
// its imports), declares the same identity.
const LOCAL_MARKER = 'const msg = (text) => text;';
// A module in utils/ imports it as './i18n.js'.
const MSG_IMPORT = /^import \{[^}]*\bmsg\b[^}]*\} from '(?:(?:\.\.\/|\.\/)+utils\/|\.\/)i18n\.js';$/m;
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
    .flatMap((path) => markerProblems(readFileSync(path, 'utf8'), webPath(path))), []);
});

test('every catalog matches the keys and placeholders of the sources', async () => {
  const { keys: english } = await sourceKeys();
  assert.deepEqual(readdirSync(CATALOG_ROOT).sort(), catalogFiles(), 'one catalog per UI language except English');
  for (const file of catalogFiles()) {
    const { default: catalog } = await import(pathToFileURL(join(CATALOG_ROOT, file)).href);
    assert.ok(Object.isFrozen(catalog), `${file} exports a frozen object`);
    assert.deepEqual(catalogProblems(english, catalog), [], file);
    // Stale keys wait for the translation batch (sync.mjs apply drops them);
    // only the release tier fails on them.
    const { translated, stale } = catalogStatus(english.keys(), catalog);
    console.log(`${file}: ${translated.length} of ${english.size} keys translated, ${stale.length} stale`);
  }
  // The keys a module derives rather than writes as literals.
  assert.ok(english.has(`member role${CONTEXT_SEPARATOR}anchor`), 'the member value labels (OV-1390)');
  assert.ok(english.has(`rail${CONTEXT_SEPARATOR}Features`), 'the rail short names with a context (OV-1417)');
  console.log(`${english.size} keys in the sources; catalogs: ${catalogFiles().join(', ')}`);
});

// A locale's lint list (tests/web/i18n-lint/<locale>.tsv, docs/internal/I18N_STYLE.md):
// patterns of its style and glossary rules that a translation must not match.
const lintRoot = join(process.cwd(), 'tests', 'web', 'i18n-lint');
const LINT_COLUMNS = 'id\tregex\tflags\tmessage\tseverity';

/**
 * One row of a lint list: an `error` match fails the catalog test, a `warn` match is printed.
 * @typedef {object} LintRule
 * @property {string} id
 * @property {RegExp} pattern the `regex` and `flags` columns
 * @property {string} message
 * @property {'error' | 'warn'} severity
 */
/**
 * A value of a catalog that a lint rule matches.
 * @typedef {object} LintMatch
 * @property {LintRule} rule
 * @property {string} key
 * @property {string} value
 */

/**
 * @param {string} text the TSV: the header row, then one rule per row
 * @returns {LintRule[]}
 */
const parseLintRules = (text) => {
  const [header, ...rows] = text.split('\n').filter((line) => line !== '');
  assert.equal(header, LINT_COLUMNS);
  return rows.map((row) => {
    const cells = row.split('\t');
    assert.equal(cells.length, 5, row);
    const [id, regex, flags, message, severity] = cells;
    assert.ok(severity === 'error' || severity === 'warn', `${id}: severity ${severity}`);
    return { id, pattern: new RegExp(regex, flags), message, severity };
  });
};

/**
 * Every text of a catalog that a rule matches, each plural form on its own,
 * with each `{name}` placeholder masked as `{X}` first: a rule reads the
 * translation, never the placeholder names.
 * @param {readonly LintRule[]} rules
 * @param {import('../../gbdraw/web/js/utils/i18n.js').TranslationCatalog} catalog
 * @returns {LintMatch[]}
 */
const lintMatches = (rules, catalog) => Object.entries(catalog).flatMap(([key, entry]) => (
  (typeof entry === 'string' ? [entry] : Object.values(entry)).flatMap((value) => {
    const masked = value.replace(PLACEHOLDER, '{X}');
    return rules.filter((rule) => rule.pattern.test(masked)).map((rule) => ({ rule, key, value }));
  })
));

test('lint check: placeholder names are masked; a match reports its rule and key', () => {
  const rules = parseLintRules(`${LINT_COLUMNS}\nx-upload\t\\bupload\\b\tui\tWrite Unggah.\terror\nx-dots\t\\.\\.\\.\tu\tWrite one ….\twarn\n`);
  assert.deepEqual(lintMatches(rules, { 'Upload {upload}': 'Unggah {upload}', 'Upload': 'Upload', 'Wait...': 'Tunggu...' })
    .map(({ rule, key }) => [rule.id, rule.severity, key]), [['x-upload', 'error', 'Upload'], ['x-dots', 'warn', 'Wait...']]);
  assert.deepEqual(lintMatches(rules, { '{count} uploads': { few: '{count} upload', many: 'Tunggu...', other: 'x' } })
    .map(({ rule, value }) => [rule.id, value]), [['x-upload', '{count} upload'], ['x-dots', 'Tunggu...']], 'every plural form is read');
  assert.throws(() => parseLintRules('id\tregex\n'));
  assert.throws(() => parseLintRules(`${LINT_COLUMNS}\nx\ta\tu\tm\tinfo\n`));
});

test('every catalog has a lint list and matches none of its error rules', async () => {
  const lists = readdirSync(lintRoot).filter((name) => name.endsWith('.tsv')).sort();
  assert.deepEqual(lists, CATALOG_LOCALES.map((locale) => `${locale}.tsv`).sort());
  for (const list of lists) {
    const locale = list.slice(0, -'.tsv'.length);
    const rules = parseLintRules(readFileSync(join(lintRoot, list), 'utf8'));
    assert.equal(new Set(rules.map(({ id }) => id)).size, rules.length, `${list}: rule ids are unique`);
    const { default: catalog } = await import(pathToFileURL(join(CATALOG_ROOT, `${locale}.js`)).href);
    const matches = lintMatches(rules, catalog);
    const format = ({ rule, key, value }) => `${locale} ${rule.id} ${JSON.stringify(key)} -> ${JSON.stringify(value)}: ${rule.message}`;
    for (const match of matches.filter(({ rule }) => rule.severity === 'warn')) console.log(`warn: ${format(match)}`);
    assert.deepEqual(matches.filter(({ rule }) => rule.severity === 'error').map(format), [], list);
  }
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

// services/linear-comparisons.js marks each unplaced-draft issue with msg(); the list shows it with tm().
test('OV-1407: an unplaced comparison draft shows its issues in the UI language', () => {
  const html = readFileSync(join(webRoot, 'index.html'), 'utf8');
  assert.match(html, /<li v-for="issue in unplaced\.issues"[^>]*>\{\{ tm\(issue\.message\) \}\}<\/li>/);
});

test('OV-1387: the Record Size Mode option Linear is its own key; the diagram mode keeps Linear', () => {
  const html = readFileSync(join(webRoot, 'index.html'), 'utf8');
  const select = html.slice(html.indexOf('v-model="adv.multi_record_size_mode"'));
  assert.match(select.slice(0, select.indexOf('</select>')), /<option value="linear">\{\{ tx\('size scaling', 'Linear'\) \}\}<\/option>/);
  assert.match(html, /@click="selectDrawing\('linear'\)"[^\n]*\{\{ t\('Linear'\) \}\}/);
});

test('OV-1465: the Circular record selector has its own key; the topology checkbox keeps Circular record', () => {
  const html = readFileSync(join(webRoot, 'index.html'), 'utf8');
  const select = html.slice(html.indexOf(':value="form.circular_record_selector"'));
  assert.match(select.slice(0, select.indexOf('>')), /:aria-label="tx\('record selector', 'Circular record'\)"/);
  assert.match(html, /controller\.setTopology\(row, \$event\.target\.checked\)[^\n]*> \{\{ t\('Circular record'\) \}\}<\/label>/);
});

test('OV-1447: a stroke-width placeholder has its own key; the Circular track slot Width label keeps Width', () => {
  const html = readFileSync(join(webRoot, 'index.html'), 'utf8');
  const placeholders = [...html.matchAll(/:placeholder="([^"]*Width[^"]*)"/g)].map((match) => match[1]);
  assert.deepEqual(placeholders, ["tx('stroke width', 'Width')", "tx('stroke width', 'Width')"]);
  assert.match(readFileSync(join(jsRoot, 'app', 'circular-track-slots', 'measure-input.js'), 'utf8'), /label: msg\('Width'\)/);
});

// OV-1475: a sentence that counts goes through tn(), so each locale words one
// and many apart; a plain t() key never writes the plural as "(s)", and the
// counted sentences below have a `one` and an `other` form.
test('OV-1475: a counted sentence is a tn() call, never a plain key with (s)', () => {
  const plain = sourceCalls().filter((call) => call.key && call.name !== 'tn');
  assert.deepEqual(plain.filter(({ key }) => /\{count\}/.test(key) && /\((?:s|es)\)/.test(key)).map(({ where, key }) => `${where} ${key}`), []);
  const plainKeys = new Set(plain.map(({ key }) => key));
  for (const key of [
    'Remove {count} unmatched feature edits', '{count} similarity groups', '{drawn} of {total} records drawn',
    '{drawn} of {count} records (separate diagrams)', 'There are {count} other features with the legend name "{name}".'
  ]) assert.ok(!plainKeys.has(key), key);
  const keys = new Set(sourceCalls().map((call) => call.key));
  for (const key of [
    'Remove {count} unmatched feature edit', '{count} similarity group', '{drawn} of {count} record drawn',
    'There is {count} other feature with the legend name "{name}".', 'Loaded {count} row.', 'Applied to {count} label.'
  ]) assert.ok(keys.has(key), key);
  // A count written into a template literal bypasses the catalogs; only an
  // Error message (read by error-normalization.js) may still say "(s)".
  const counted = /`[^`]*\$\{[^}]+\}[^`]*?\w\((?:s|es)\)[^`]*`/g;
  const literals = sourceFiles().flatMap((path) => {
    const text = readFileSync(path, 'utf8');
    return [...text.matchAll(counted)]
      .filter((match) => !/Error\(\s*$/.test(text.slice(Math.max(0, match.index - 40), match.index)))
      .map((match) => `${webPath(path)} ${match[0].slice(0, 80)}`);
  });
  assert.deepEqual(literals, []);
});

// OV-1476: the English names only what the UI shows.
test('OV-1476: no key names a control, row or syntax the UI does not have', async () => {
  const keys = [...(await sourceKeys()).keys.keys()];
  const stale = [/All adjacent LOSAT/, /the GC (?:content|skew) row/, /\(\?i\)/, /style size/, /TLOSATX code/, /from file \{n\}/, /#index/, /^\{count\} needs? a choice$/];
  assert.deepEqual(keys.filter((key) => stale.some((pattern) => pattern.test(key))), []);
  for (const key of [
    'Selected LOSAT comparisons require LOSATP Pairwise mode. For Similarity groups or Collinear blocks, choose Run LOSAT for all adjacent pairs.',
    '{label} font size', '{count} record awaiting a choice', 'Use #<number> to distinguish records with the same ID.'
  ]) assert.ok(keys.includes(key), key);
});

// OV-1477: a capitalised label or a control name is not dropped into a
// sentence; each case is a whole sentence of its own.
test('OV-1477: track, mode, control and setting names are not placeholders inside sentences', () => {
  const keys = sourceCalls().map((call) => call.key).filter(Boolean);
  const inserted = [/\{track\}/, /\{control\}/, /No \{mode\} Result/, /^\{setting\} removes/, /\{mode\} track/, /\{topology\}/];
  assert.deepEqual(keys.filter((key) => inserted.some((pattern) => pattern.test(key))), []);
  for (const key of [
    'Hide GC content?', 'Keep GC skew', 'Hidden: Hide GC Content is on.', 'No Circular Result yet. Configure settings and click Generate.',
    'Choosing "{setting}" removes the lane for {count} feature you placed manually.', 'This change removes the lane for {count} features you placed manually.'
  ]) assert.ok(keys.includes(key), key);
});

// The cross-locale check (i18n-crossfix-20261022): the English names the
// control, option or value the UI shows and says what happens.
test('the English cites the controls the UI has and names what each text describes', async () => {
  const keys = [...(await sourceKeys()).keys.keys()];
  const stale = [
    /Export (?:outputs|uses) the current Result/, /^File order is unavailable/, /cannot use Overlay\.$/, /^Search field$/,
    /^Orient forward/, /direction changes(?: restored)?\.?$/, /comparison endpoints/i, /was not found in the current input/,
    /Zipped/, /Use Move outside or Move inside/, /'First' shows/, /name no record or feature here/, /per-record blanks/,
    /blank TSV cell/, /^Tick Font$/, /^Annotation feature value$/, /^No features selected$/, /^Color for \{file\}$/,
    /^Upload BLAST$/, /^Advanced similarity groups$/, /^Set on individual features$/, /member hits per protein/i, /^Depth \{n\}$/
  ];
  assert.deepEqual(keys.filter((key) => stale.some((pattern) => pattern.test(key))), []);
  for (const key of [
    'Field to search', 'Tick Font Size', 'Annotation qualifier value', 'No feature types selected', 'Comparison ring color for {file}',
    'Max hits per protein for groups', 'Color set on individual features',
    'Upload BLAST TSV', 'Advanced similarity group settings'
  ]) assert.ok(keys.includes(key), key);
});

// OV-1446: the legend-name field of a Depth series shows, as its placeholder,
// the default name the figure draws, which is English in every locale.
test('OV-1446: a Depth legend-name placeholder is the English default name of the figure', async () => {
  const html = readFileSync(join(webRoot, 'index.html'), 'utf8');
  assert.deepEqual([...html.matchAll(/:placeholder="[^"]*\bt\('Depth/g)].map((match) => match[0]), []);
  // The series cards' Legend title fields, then the Depth rows' legend-label fields (OV-1508).
  assert.deepEqual([...html.matchAll(/:placeholder="(depthSeriesDefaultName\([^"]*)"/g)].map((match) => match[1]), [
    'depthSeriesDefaultName(track.index)', 'depthSeriesDefaultName(track.index)',
    'depthSeriesDefaultName(normalizeDepthSlotTrackIndex(entry.slot))', 'depthSeriesDefaultName(normalizeDepthSlotTrackIndex(entry.slot))'
  ]);
  const { defaultDepthSeriesLabel } = await import('../../gbdraw/web/js/services/depth-track-state.js');
  try {
    await setUiLocale('ja', async () => ({ default: Object.freeze({ Depth: '深度', 'Depth {n}': '深度{n}' }) }));
    assert.deepEqual([defaultDepthSeriesLabel(0, 1), defaultDepthSeriesLabel(0, 2), defaultDepthSeriesLabel(1, 2)], ['Depth', 'Depth 1', 'Depth 2'],
      'one series: Depth; several: Depth N from 1, as gbdraw/analysis/depth_tracks.py names them');
  } finally {
    await setUiLocale(SOURCE_LOCALE);
  }
});

// OV-1506: a rail name too wide for its entry has a rail key of its own.
test('OV-1506: Layout, Labels and Features have rail keys of their own', async () => {
  const rail = [...(await sourceKeys()).keys.values()].filter(({ context }) => context === 'rail');
  assert.deepEqual(rail.map(({ context, english }) => `${context}:${english}`), ['rail:Layout', 'rail:Features', 'rail:Labels']);
});

test('the UI languages are the seventeen of the Language menu, Latin script first, then Cyrillic, then Thai, each named in itself', () => {
  assert.deepEqual(UI_LOCALES.map(({ tag, name }) => [tag, name]), [
    ['en', 'English'], ['id', 'Bahasa Indonesia'], ['de', 'Deutsch'], ['es', 'Español'], ['fr', 'Français'], ['it', 'Italiano'],
    ['pl', 'Polski'], ['pt', 'Português'], ['vi', 'Tiếng Việt'], ['tr', 'Türkçe'],
    ['ru', 'Русский'], ['uk', 'Українська'], ['th', 'ไทย'], ['ja', '日本語'], ['zh-Hans', '简体中文'], ['zh-Hant', '繁體中文'], ['ko', '한국어']
  ]);
  assert.ok(Object.isFrozen(UI_LOCALES) && UI_LOCALES.every(Object.isFrozen));
  for (const { tag } of UI_LOCALES) assert.equal(uiLocaleOf(tag), tag);
  for (const value of ['', 'EN', 'ja-JP', 'zh', 'zh-CN', 'pt-BR', 'pt-PT', 'es-ES', 'fr-FR', 'ko-KR', 'id-ID', 'in', 'vi-VN', 'th-TH', 'de-DE', 'DE', 'it-IT', 'pl-PL', 'tr-TR', 'ru-RU', 'uk-UA', 'ua', 'xx', null, undefined, 3, '../x']) {
    assert.equal(uiLocaleOf(value), null, JSON.stringify(value));
  }
});

test('a locale outside the UI languages is refused', async () => {
  await assert.rejects(setUiLocale('xx', async () => ({ default: {} })), /Unknown UI locale/);
  assert.equal(uiLocale.value, SOURCE_LOCALE);
});

test('two sentences are joined with a space in English, the other Latin-script languages, Thai and Korean and with nothing in Japanese and Chinese', async () => {
  assert.equal(sentenceSeparator(), ' ');
  assert.equal(joinSentences('Exported a.svg.', 'Find it in your browser\'s downloads.'), 'Exported a.svg. Find it in your browser\'s downloads.');
  try {
    for (const tag of ['id', 'de', 'es', 'fr', 'it', 'pl', 'pt', 'vi', 'tr', 'ru', 'uk', 'th', 'ko']) {
      await setUiLocale(tag, async () => ({ default: {} }));
      assert.equal(sentenceSeparator(), ' ', tag);
      assert.equal(joinSentences('A.', 'B.'), 'A. B.', tag);
    }
    for (const tag of ['ja', 'zh-Hans', 'zh-Hant']) {
      await setUiLocale(tag, async () => ({ default: {} }));
      assert.equal(sentenceSeparator(), '', tag);
      assert.equal(joinSentences('保存しました。', 'ダウンロードを確認してください。'), '保存しました。ダウンロードを確認してください。', tag);
    }
  } finally {
    await setUiLocale(SOURCE_LOCALE);
  }
});

test('a list keeps today\'s English and uses the locale\'s list pattern in every other locale', async () => {
  assert.equal(listText([]), '');
  assert.equal(listText(['A']), 'A');
  assert.equal(listText(['A', 'B']), 'A and B');
  assert.equal(listText(['A', 'B', 'C']), 'A, B, and C');
  assert.equal(listText(['A', 'B', 'C', '2 more']), 'A, B, C, and 2 more');
  try {
    for (const tag of ['id', 'de', 'es', 'fr', 'it', 'pl', 'pt', 'vi', 'tr', 'ru', 'uk', 'th', 'ja', 'zh-Hans', 'zh-Hant', 'ko']) {
      await setUiLocale(tag, async () => ({ default: {} }));
      const expected = new Intl.ListFormat(tag, { type: 'conjunction' });
      assert.equal(listText(['A']), 'A', tag);
      assert.equal(listText(['A', 'B']), expected.format(['A', 'B']), tag);
      assert.equal(listText(['A', 'B', 'C']), expected.format(['A', 'B', 'C']), tag);
      assert.doesNotMatch(listText(['A', 'B', 'C']), /\band\b/, tag);
    }
    // Thai separates the items with spaces and no comma, and writes และ against the last one (STYLE 11.3).
    for (const [tag, two, three] of [['id', 'A dan B', 'A, B, dan C'], ['vi', 'A và B', 'A, B và C'], ['th', 'AและB', 'A B และC']]) {
      await setUiLocale(tag, async () => ({ default: {} }));
      assert.deepEqual([listText(['A', 'B']), listText(['A', 'B', 'C'])], [two, three], tag);
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

test('the shown locale\'s plural rules pick the form, and the English text of that form is the key', async () => {
  const catalog = Object.freeze({ '{count} row': 'one:{count}', '{count} rows': 'other:{count}' });
  const forms = (counts) => counts.map((count) => tn(count, '{count} row', '{count} rows'));
  assert.deepEqual(forms([0, 1, 2]), ['0 rows', '1 row', '2 rows'], 'English is unchanged');
  try {
    for (const [tag, expected] of [
      // fr and pt take `one` for 0 too; es only for 1.
      ['fr', ['one:0', 'one:1', 'other:2', 'other:1000000']],
      ['pt', ['one:0', 'one:1', 'other:2', 'other:1000000']],
      // es and fr select `many` for 1000000, which takes `other`.
      ['es', ['other:0', 'one:1', 'other:2', 'other:1000000']],
      // de and tr take `one` for 1 only; it selects `many` for 1000000, pl, ru and uk `few` for 2 and `many` for 0
      // and 1000000, and a string value shows its `other` text for those.
      ['de', ['other:0', 'one:1', 'other:2', 'other:1000000']],
      ['it', ['other:0', 'one:1', 'other:2', 'other:1000000']],
      ['pl', ['other:0', 'one:1', 'other:2', 'other:1000000']],
      ['tr', ['other:0', 'one:1', 'other:2', 'other:1000000']],
      ['ru', ['other:0', 'one:1', 'other:2', 'other:1000000']],
      ['uk', ['other:0', 'one:1', 'other:2', 'other:1000000']],
      ['ja', ['other:0', 'other:1', 'other:2', 'other:1000000']],
      ['zh-Hans', ['other:0', 'other:1', 'other:2', 'other:1000000']],
      ['zh-Hant', ['other:0', 'other:1', 'other:2', 'other:1000000']],
      ['ko', ['other:0', 'other:1', 'other:2', 'other:1000000']],
      // Indonesian, Vietnamese and Thai have no plural forms; their `one` rows repeat the `other` text.
      ['id', ['other:0', 'other:1', 'other:2', 'other:1000000']],
      ['vi', ['other:0', 'other:1', 'other:2', 'other:1000000']],
      ['th', ['other:0', 'other:1', 'other:2', 'other:1000000']]
    ]) {
      await setUiLocale(tag, async () => ({ default: catalog }));
      assert.deepEqual(forms([0, 1, 2, 1000000]), expected, tag);
    }
    // A missing translation shows the form English picks.
    await setUiLocale('fr', async () => ({ default: { '{count} rows': 'other:{count}' } }));
    assert.deepEqual(forms([0, 1, 2]), ['0 rows', '1 row', 'other:2']);
    // With an object on the `other` form, the `one` form is the text for exactly 1 and every other count takes the
    // object's text of its category, else its `other`: fr selects `one` for 0 and `many` for 1000000.
    await setUiLocale('fr', async () => ({ default: { '{count} row': 'one:{count}', '{count} rows': { many: 'many:{count}', other: 'other:{count}' } } }));
    assert.deepEqual(forms([0, 1, 2, 1000000]), ['other:0', 'one:1', 'other:2', 'many:1000000']);
    await setUiLocale(SOURCE_LOCALE);
    assert.deepEqual(forms([0, 1]), ['0 rows', '1 row'], 'English again after a switch back');
  } finally {
    await setUiLocale(SOURCE_LOCALE);
  }
});

test('a language with few and many forms: the `one` form is the text for exactly 1, the object on the `other` form gives every other count', () => {
  // Russian (Intl.PluralRules('ru')): one for 1, 21 and 31, few for 2-4 and 22, many for 0, 5, 11 and 25. A `one`
  // text without {count} ("From this color rule:") must not show for 21, so 21 takes the object's `one`.
  const russian = new Intl.PluralRules('ru');
  const catalog = Object.freeze({
    '{count} row': 'exactly-one:{count}',
    '{count} rows': Object.freeze({ one: 'one:{count}', few: 'few:{count}', many: 'many:{count}', other: 'other:{count}' })
  });
  const shown = (entries, counts) => counts.map((count) => pluralTranslation(russian, entries, count, '{count} row', '{count} rows'));
  assert.deepEqual(shown(catalog, [1, 21, 31, 2, 3, 4, 22, 0, 5, 11, 25]), [
    'exactly-one:{count}', 'one:{count}', 'one:{count}', 'few:{count}', 'few:{count}', 'few:{count}', 'few:{count}',
    'many:{count}', 'many:{count}', 'many:{count}', 'many:{count}'
  ]);
  assert.equal(pluralTranslation(russian, catalog, 1.5, '{count} row', '{count} rows'), 'other:{count}', 'a fraction takes `other`');
  assert.deepEqual(shown({ '{count} row': 'exactly-one', '{count} rows': { few: 'few', many: 'many', other: 'other' } }, [1, 21, 2]),
    ['exactly-one', 'other', 'few'], 'a category the object leaves out takes `other`; the `one` form shows only for 1');
  assert.deepEqual(shown({ '{count} rows': { one: 'one', other: 'other' } }, [1, 21]), [null, 'one'], 'a missing `one` form leaves 1 to the English');
  // String values keep the category rule: the `one` form for every count the rules call `one`.
  assert.deepEqual(shown({ '{count} row': 'one:{count}', '{count} rows': 'rows:{count}' }, [1, 21, 2, 5]), ['one:{count}', 'one:{count}', 'rows:{count}', 'rows:{count}']);
  assert.deepEqual(shown({ '{count} rows': 'rows:{count}' }, [1, 2, 5]), [null, 'rows:{count}', 'rows:{count}'],
    'a string value shows for few and many; a missing key is null, and tn() shows the English');
});

test('every `one` text of a tn() call is a catalog key', () => {
  const keys = new Set(sourceCalls().map((call) => call.key));
  for (const key of ['From this color rule:', 'The chosen feature already points left.', '({count} record)']) assert.ok(keys.has(key), key);
});

test('OV-1415: the words a declined label table names as changed are a key, shown in the UI language', async () => {
  const { RULE_INPUT_CHANGE_WORDS } = await import('../../gbdraw/web/js/utils/error-normalization.js');
  assert.ok(new Set(sourceCalls().map((call) => call.key)).has(RULE_INPUT_CHANGE_WORDS));
  const source = readFileSync(join(jsRoot, 'app', 'feature-editor', 'label-actions.js'), 'utf8');
  assert.match(source, /\{ inputs: tm\(RULE_INPUT_CHANGE_WORDS\) \}/);
});

// OV-1690: the record list sentence ends inside its key, so a locale writes its
// own punctuation, and the next sentence follows the locale's separator.
test('OV-1690: the record list sentence is one key, joined to the next with the sentence separator', () => {
  const html = readFileSync(new URL('../../gbdraw/web/index.html', import.meta.url), 'utf8');
  const start = html.indexOf('<p id="record-list-description"');
  const block = html.slice(start, html.indexOf('</p>', start));
  assert.match(block, /<i18n-text :text="t\('\{file\}: \{records\}\.'\)">/);
  assert.match(block, /sentenceSeparator\(\) \+ settingApplicationSentence\('inputs'\)/);
  assert.doesNotMatch(block, /<\/span>\.|<\/strong>:/);
});

// A literal with words as a branch of a shown expression is English in every
// locale: `{{ on ? 'reverse' : 'forward' }}`. A literal inside a call
// (`t('…')`, `settingApplicationText(on ? 'live' : 'generate')`) and the right
// side of a comparison are not branches.
/**
 * The word literals an expression can show: its whole value, or a branch of its
 * outermost conditional.
 * @param {string} expression
 * @returns {string[]}
 */
const bareEnglishBranches = (expression) => {
  /** @type {string[]} */
  const found = [];
  let depth = 0;
  let afterBranch = true;
  for (let at = 0; at < expression.length; at += 1) {
    const char = expression[at];
    if (char === '\'' || char === '"' || char === '`') {
      let end = at + 1;
      while (end < expression.length && expression[end] !== char) end += expression[end] === '\\' ? 2 : 1;
      const literal = expression.slice(at, end + 1);
      if (depth === 0 && afterBranch && /[A-Za-z]/.test(literal.replace(/\$\{[^}]*\}/g, ''))) found.push(literal);
      afterBranch = false;
      at = end;
    } else if (/\s/.test(char)) {
      continue;
    } else {
      if ('([{'.includes(char)) depth += 1;
      if (')]}'.includes(char)) depth -= 1;
      afterBranch = depth === 0 && (char === ':' || (char === '?' && expression[at + 1] !== '.' && expression[at + 1] !== '?' && expression[at - 1] !== '?'));
    }
  }
  return found;
};
const SHOWN_BINDING = /\s:(aria-label|title|placeholder|label|alt|data-input-aria-label)="([^"]*)"/g;

test('the bare English branch check finds a literal branch and leaves calls and comparisons alone', () => {
  assert.deepEqual(bareEnglishBranches("on ? 'reverse' : 'forward'"), ["'reverse'", "'forward'"]);
  assert.deepEqual(bareEnglishBranches("a ? t('A') : b === 'ready' ? t('B') : 'Align…'"), ["'Align…'"]);
  assert.deepEqual(bareEnglishBranches('`Depth TSV for sequence ${index + 1}`'), ['`Depth TSV for sequence ${index + 1}`']);
  assert.deepEqual(bareEnglishBranches("label(on ? 'live' : 'generate')"), []);
  assert.deepEqual(bareEnglishBranches("a ? t('A') : t('B')"), []);
  assert.deepEqual(bareEnglishBranches("index ? ', ' : ''"), []);
  assert.deepEqual(bareEnglishBranches('unit ? ` ${unit}` : \'\''), []);
  assert.deepEqual(bareEnglishBranches("a?.b ?? c"), []);
});

test('no mustache or shown attribute of index.html has a bare English literal as a branch', () => {
  const html = readFileSync(new URL('../../gbdraw/web/index.html', import.meta.url), 'utf8');
  const lineOf = (/** @type {number} */ index) => html.slice(0, index).split('\n').length;
  const found = [
    ...[...html.matchAll(/\{\{([\s\S]*?)\}\}/g)].map((match) => [match.index, match[1]]),
    ...[...html.matchAll(SHOWN_BINDING)].map((match) => [match.index, match[2]])
  ].flatMap(([index, expression]) => bareEnglishBranches(String(expression)).map((literal) => `index.html:${lineOf(Number(index))}: ${literal}`));
  assert.deepEqual(found, []);
});

test('the record row names the topology its file declares through one key per value', () => {
  assert.deepEqual(['circular', 'linear', 'unknown', undefined].map(detectedTopologyNote),
    ['(file: circular)', '(file: linear)', '(file: topology unknown)', '(file: topology unknown)']);
});

// A "Choose" sentence names the options as their controls show them: a label of
// the page, or of the two option tables the page lists (comparison-ui, pairwise-match-popup).
test('each name in a Choose sentence of an error is a control label', () => {
  const read = (path) => readFileSync(new URL(`../../gbdraw/web/${path}`, import.meta.url), 'utf8');
  const labels = [read('index.html'), read('js/app/comparison-ui.js'), read('js/app/pairwise-match-popup.js')].join('\n');
  const sentences = [...read('js/utils/error-normalization.js').matchAll(/: '(Choose [A-Z][^']*)\.'/g)].map((match) => match[1])
    .filter((sentence) => !sentence.includes(' for '));
  assert.ok(sentences.length >= 11, String(sentences.length));
  const isLabel = (name) => labels.includes(`'${name}')`) || labels.includes(`>${name}</option>`);
  const missing = sentences.flatMap((sentence) => sentence.slice('Choose '.length).split(/, or |, | or /)
    .filter((name) => !isLabel(name)).map((name) => `${sentence}: ${name}`));
  assert.deepEqual(missing, []);
});

// Python's color modes reach the shared table through the generator, so a new mode fails here.
test('the Collinear color mode table lists the modes Python writes', () => {
  assert.deepEqual(Object.keys(COLLINEAR_COLOR_MODE_TEXT), [...COLLINEARITY_COLOR_MODES]);
  const page = readFileSync(new URL('../../gbdraw/web/index.html', import.meta.url), 'utf8');
  assert.match(page, /v-for="\[value, label\] in collinearColorModeOptions"/);
  assert.doesNotMatch(page, /<option value="orientation_identity">/);
});

// The copied summary uses the keys of the popup that shows it.
test('the copied Collinear block summary holds no English of its own', () => {
  const source = readFileSync(new URL('../../gbdraw/web/js/app/pairwise-match-popup.js', import.meta.url), 'utf8');
  assert.doesNotMatch(source, /summaryText: `Collinear block:/);
});

// The sentences of the Source recipe name no mode: the recipe is of one diagram.
test('the record topology a file declares and the option values of an error are whole keys', async () => {
  const keys = [...(await sourceKeys()).keys.keys()];
  for (const key of [
    '(file: circular)', '(file: linear)', '(file: topology unknown)',
    'Choose Error on overflow, Compress, or Clip.', 'Choose Foreground or Underlay.', 'Choose Auto (Default), Linear, or Equal.',
    'Choose Above, Features on axis, or Below.', 'Choose Auto or Above Feature.',
    'Choose Pairwise matches, Similarity groups, or Collinear blocks.', 'Choose Adjacent pairs or All records.',
    'Choose Auto, CDS, or Locus.', 'Choose All hits or Reciprocal best hits.', 'Choose Strand or Gene order.',
    'Choose Average identity, Orientation, or Orientation + identity.',
    'Depth TSV for file {file}, series {series}', 'Depth TSV for sequence {sequence}, series {series}',
    'BLAST TSV for #{from} to #{to}',
    'vs {record}', 'Not moved: {reason}', 'Not moved', 'Session file: {path} could not be read and is left out.',
    'record direction\u0004forward', 'record direction\u0004reverse', 'collinearity block kind\u0004syntenic'
  ]) assert.ok(keys.includes(key), key);
  const stale = [
    /^Choose error, compress/, /^Choose foreground or/, /^Choose auto,/, /^Choose above, middle/, /^Choose auto or above_feature/,
    /^Choose pairwise,/, /^Choose adjacent or all/, /^Choose all or rbh/, /^Choose strand or order/, /^Choose average_identity/,
    /^vs$/, /^Not moved \(/, /^skipped$/, /\{name\} \(\{path\}\)/
  ];
  assert.deepEqual(keys.filter((key) => stale.some((pattern) => pattern.test(key))), []);
});

// A sentence reads whole in each case: no fragment takes its place, and a count is a tn() pair.
test('sentences with a count or a qualifier are whole keys, and none is built from a translated fragment', async () => {
  const { keys } = await sourceKeys();
  const pairOf = (key) => keys.get(key)?.plural;
  for (const [one, other] of [
    ['The label reading "{text}"', 'All {count} labels reading "{text}"'],
    ['The "{name}" feature', 'All {count} "{name}" features'],
    ['one more', '{count} more'],
    ['and one more', 'and {count} more'],
    ['Running LOSAT: {completed}/{total} search complete', 'Running LOSAT: {completed}/{total} searches complete'],
    ['Running {program} comparisons for the circular diagram: {completed}/{total} job complete', 'Running {program} comparisons for the circular diagram: {completed}/{total} jobs complete'],
    ['Importing it removes the annotation set.', 'Importing it removes the {count} annotation sets.'],
    ['The chosen feature already points left.', 'The {count} chosen features already point left.'],
    ['The chosen feature already points right.', 'The {count} chosen features already point right.'],
    ['Reverses {records} so the chosen feature points left.', 'Reverses {records} so the {count} chosen features point left.'],
    ['Reverses {records} so the chosen feature points right.', 'Reverses {records} so the {count} chosen features point right.']
  ]) assert.deepEqual(pairOf(other), { one, other }, other);
  for (const key of [
    'The Source recipe references a browser-generated file: {files}. Use "Download reproducibility files" to save it next to the uploaded inputs before running the command in a terminal.',
    'The Source recipe references browser-generated files: {files}. Use "Download reproducibility files" to save them next to the uploaded inputs before running the command in a terminal.',
    'Exact replay references {files}, which is included in "Download reproducibility files".',
    'Exact replay references {files}, which are included in "Download reproducibility files".',
    'This tool runs entirely in your browser.', 'It uses WebAssembly (Pyodide).',
    '{feature} where a qualifier matches “{pattern}”', 'Any feature where a qualifier matches “{pattern}”',
    'All features that share:',
    'Room kept toward the center for the inside tracks: {px} px.',
    'Room kept toward the center for the inside tracks (among them GC content, GC skew, and the scale): {px} px.',
    'number of groups in a block\u0004Similarity groups'
  ]) assert.ok(keys.has(key), key);
  for (const key of [
    '{first} and {second}', 'the scale', 'a qualifier', 'entirely in your browser', 'This tool runs {where} using WebAssembly (Pyodide).', 'All features with the same'
  ]) assert.ok(!keys.has(key), key);
  assert.deepEqual([...keys.keys()].filter((key) => /file\(s\)/.test(key)), []);
});
