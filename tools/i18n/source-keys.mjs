// @ts-check
// The catalog keys of the Web UI sources: the one owner of "which English
// texts the catalogs translate". English is the key: the literal first
// argument of `t()` and of `msg()` (a text kept until it is shown, which
// `tm()` translates there), the `one` and the `other` argument of `tn()` (the
// shown locale's plural rules pick one of them), `context + '\u0004' + text`
// of `tx()`, the table values of the error wording owner
// (tests/web/helpers/error-text-keys.mjs), the labels of the similarity group
// member values that Python owns (`memberValueLabelKeys`, shown with `tmx()`),
// and the rail short names with a context key of their own
// (`settingsShortKeys`). tests/web/i18n-catalogs.test.mjs, tools/i18n/sync.mjs
// and tools/i18n/review-page-data.mjs read the keys from here.
import { readFileSync, readdirSync } from 'node:fs';
import { join, relative, sep } from 'node:path';
import { fileURLToPath } from 'node:url';

import { CONTEXT_SEPARATOR } from '../../gbdraw/web/js/utils/i18n.js';
import { memberValueLabelKeys } from '../../gbdraw/web/js/app/pairwise-match-popup.js';
import { SETTINGS_SECTIONS, settingsShortKeys } from '../../gbdraw/web/js/services/settings-sections.js';
import { errorTextKeys } from '../../tests/web/helpers/error-text-keys.mjs';

export const REPOSITORY_ROOT = fileURLToPath(new URL('../../', import.meta.url));
export const WEB_ROOT = join(REPOSITORY_ROOT, 'gbdraw', 'web');
const JS_ROOT = join(WEB_ROOT, 'js');
// The runtime's own source and the catalogs hold no UI text of their own.
const RUNTIME_FILES = new Set(['utils/i18n.js']);

/** @typedef {'t' | 'msg' | 'tn' | 'tx'} CallName */
/**
 * What a reader sees the text as, where the source shows it: an attribute or
 * element of index.html, or the property a module stores the text under.
 * @typedef {'button' | 'label' | 'option' | 'placeholder' | 'title' | 'aria-label' | 'heading' | 'help' | 'dialog' | 'notice' | 'other'} KeyRole
 */
/**
 * One place that produces a key.
 * @typedef {object} KeySite
 * @property {string} file path under gbdraw/web
 * @property {number} line 1-based
 * @property {KeyRole} role
 * @property {string} snippet the source line, trimmed to 220 characters
 * @property {string} area the part of the app: a settings section, `Page`, `Shared components`, or `Script: <file>`
 */
/**
 * A catalog key of the sources.
 * @typedef {object} SourceKey
 * @property {string} key the catalog key
 * @property {string} english the English text shown
 * @property {string | null} context the `tx()` context, or null
 * @property {{ one: string, other: string } | null} plural the `tn()` pair the key is one form of, or null
 * @property {KeySite[]} sites
 */
/**
 * A call read from a source: its key, or the problem that keeps the key from being read.
 * @typedef {{ where: string, name: CallName, key: string, english: string, index: number, plural: { one: string, other: string } | null } | { where: string, problem: string }} SourceCall
 */

// For each function, the argument positions that are literal keys, and the
// catalog keys they make, each with its English text.
/** @type {Record<CallName, { literals: number[], keys: (args: (string | null)[]) => [string, string][] }>} */
const CALLS = {
  t: { literals: [0], keys: ([text]) => [[String(text), String(text)]] },
  msg: { literals: [0], keys: ([text]) => [[String(text), String(text)]] },
  tn: { literals: [1, 2], keys: ([, one, other]) => [[String(one), String(one)], [String(other), String(other)]] },
  tx: { literals: [0, 1], keys: ([context, text]) => [[context + CONTEXT_SEPARATOR + text, String(text)]] }
};
const CALL_START = /(?<![\w$.])(msg|tn|tx|t)\(/g;
// `tm()` translates a value that `msg()` marked; a literal there is a `t()` call.
const TM_LITERAL = /(?<![\w$.])tm\(\s*['"`]/g;
// `tmx()` translates a derived text; a literal text there is a `tx()` call.
const TMX_LITERAL = /(?<![\w$.])tmx\(\s*(['"`])(?:\\.|(?!\1)[^\\])*\1\s*,\s*['"`]/g;
/** @type {Record<string, string>} */
const ESCAPES = { n: '\n', t: '\t', r: '\r', b: '\b', f: '\f', v: '\v', 0: '\0' };
export const PLACEHOLDER = /\{([A-Za-z_][A-Za-z0-9_]*)\}/g;

/**
 * The 1-based line of an index of `source`.
 * @param {string} source
 * @returns {(index: number) => number}
 */
const lineFinder = (source) => {
  const starts = [0];
  for (let at = source.indexOf('\n'); at >= 0; at = source.indexOf('\n', at + 1)) starts.push(at + 1);
  return (index) => {
    let [low, high] = [0, starts.length - 1];
    while (low < high) {
      const middle = (low + high + 1) >> 1;
      if (starts[middle] <= index) low = middle;
      else high = middle - 1;
    }
    return low + 1;
  };
};

/**
 * Reads the string literal that starts at `start` ('…', "…", or `…` without
 * `${`). Returns its value and the index after it, or null for a template
 * literal with a substitution.
 * @param {string} source
 * @param {number} start
 * @returns {{ value: string, end: number } | null}
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

/**
 * Skips one argument expression; returns the index of the `,` or `)` that ends it.
 * @param {string} source
 * @param {number} start
 */
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
 * The calls of `t`, `msg`, `tn`, `tx` in `source`: each with its key, or with
 * the problem that keeps the key from being read.
 * @param {string} source
 * @param {string} file the name that problems report
 * @returns {SourceCall[]}
 */
export const extractCalls = (source, file) => {
  /** @type {SourceCall[]} */
  const calls = [];
  const lineAt = lineFinder(source);
  for (const match of source.matchAll(CALL_START)) {
    const name = /** @type {CallName} */ (match[1]);
    const where = `${file}:${lineAt(match.index)}`;
    /** @type {(string | null)[]} */
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
    // An in-DOM template sees an entity decoded, an x-template does not: write the character.
    if (/&(?:[A-Za-z]+|#\d+|#x[0-9A-Fa-f]+);/.test(args.join(''))) {
      calls.push({ where, problem: `${where}: ${name}() key holds an HTML entity; write the character` });
      continue;
    }
    const plural = name === 'tn' ? { one: String(args[1]), other: String(args[2]) } : null;
    for (const [key, english] of CALLS[name].keys(args)) calls.push({ where, name, key, english, index: match.index, plural });
  }
  for (const match of source.matchAll(TM_LITERAL)) {
    const where = `${file}:${lineAt(match.index)}`;
    calls.push({ where, problem: `${where}: tm() takes a value that msg() marked; write t() for a literal` });
  }
  for (const match of source.matchAll(TMX_LITERAL)) {
    const where = `${file}:${lineAt(match.index)}`;
    calls.push({ where, problem: `${where}: tmx() takes a derived text; write tx() for a literal` });
  }
  return calls;
};

/** @type {Record<string, string>} */
const NAMED_ENTITIES = { quot: '"', amp: '&', lt: '<', gt: '>', apos: '\'', '#39': '\'', nbsp: ' ' };
const ENTITY = /&(#x[0-9A-Fa-f]+|#\d+|[A-Za-z][A-Za-z0-9]*);/g;

/**
 * index.html as Vue compiles it: Vue decodes entities itself, in attributes
 * and in {{ }}, in the in-DOM root and in an x-template alike. The entities of
 * the list above and numeric references are decoded the way the browser does;
 * any other named entity stays as written, so extraction can refuse it.
 * @param {string} html
 */
export const compiledHtml = (html) => html.replace(ENTITY, (whole, name) => {
  if (name[0] === '#' && name[1] !== 'x' && name[1] !== 'X') return String.fromCodePoint(Number(name.slice(1)));
  if (name[0] === '#' && /^#[xX]/.test(name)) return String.fromCodePoint(parseInt(name.slice(2), 16));
  return NAMED_ENTITIES[name] ?? whole;
});

export const VOID_TAGS = new Set(['area', 'base', 'br', 'col', 'embed', 'hr', 'img', 'input', 'link', 'meta', 'source', 'track', 'wbr']);

/**
 * The index of the end tag that closes the element whose start tag ends at
 * `from`, counting nested elements of the same name; the source length when
 * none closes it.
 * @param {string} source
 * @param {string} tag lower-case tag name
 * @param {number} from
 */
export const elementEnd = (source, tag, from) => {
  const tags = new RegExp(`<(/?)${tag}\\b(?:"[^"]*"|'[^']*'|[^'">])*?(/?)>`, 'gi');
  tags.lastIndex = from;
  let depth = 1;
  for (let next = tags.exec(source); next; next = tags.exec(source)) {
    if (next[1]) depth -= 1;
    else if (!next[2] && !VOID_TAGS.has(tag)) depth += 1;
    if (depth === 0) return next.index;
  }
  return source.length;
};

/** The Web UI sources that hold keys: index.html and every module except the runtime and the catalogs. */
export const sourceFiles = () => [
  join(WEB_ROOT, 'index.html'),
  ...readdirSync(JS_ROOT, { recursive: true })
    .map((path) => String(path))
    .filter((path) => path.endsWith('.js'))
    .filter((path) => !RUNTIME_FILES.has(path.split(sep).join('/')))
    .filter((path) => !path.split(sep).join('/').startsWith('utils/i18n-catalogs/'))
    .map((path) => join(JS_ROOT, path))
];

/** @param {string} path */
export const webPath = (path) => relative(WEB_ROOT, path).split(sep).join('/');

/** @param {string} path a source file; index.html is read as Vue compiles it */
export const sourceText = (path) => {
  const source = readFileSync(path, 'utf8');
  return path.endsWith('.html') ? compiledHtml(source) : source;
};

/** Every call of every source file. */
export const sourceCalls = () => sourceFiles().flatMap((path) => extractCalls(sourceText(path), webPath(path)));

/** @param {string} text @returns {string[]} the sorted placeholder names */
export const placeholdersOf = (text) => [...new Set([...text.matchAll(PLACEHOLDER)].map((match) => match[1]))].sort();

/** @type {[RegExp, KeyRole][]} */
const ATTRIBUTE_ROLES = [
  [/^placeholder$/, 'placeholder'], [/^title$/, 'title'], [/^aria-label$/, 'aria-label'], [/^label$/, 'label']
];
/** @type {[RegExp, KeyRole][]} */
const TAG_ROLES = [
  [/^button$/, 'button'], [/^label$/, 'label'], [/^(?:option|optgroup)$/, 'option'], [/^h[1-6]$/, 'heading'],
  [/^help-tip$/, 'help'], [/dialog$/, 'dialog']
];
/** @type {[RegExp, KeyRole][]} */
const PROPERTY_ROLES = [
  [/placeholder/i, 'placeholder'], [/aria-?label/i, 'aria-label'], [/label|short|name/i, 'label'],
  [/title|heading/i, 'title'], [/help|hint|description|detail|guidance|remedy|tip/i, 'help'],
  [/message|notice|toast|status|warning|error|reason|issue/i, 'notice'], [/button|action|confirm|cancel/i, 'button'],
  [/option/i, 'option']
];
/** @param {[RegExp, KeyRole][]} table @param {string} name @returns {KeyRole | null} */
const roleIn = (table, name) => table.find(([pattern]) => pattern.test(name))?.[1] ?? null;

const TAG = /<(\/?)([A-Za-z][\w-]*)\b((?:"[^"]*"|'[^']*'|[^'">])*?)(\/?)>/g;

/**
 * The role of an element of index.html, from its tag name and `role` attribute.
 * @param {string} tag
 * @param {string} attributes
 * @returns {KeyRole | null}
 */
const elementRole = (tag, attributes) => (/\srole="(?:alert)?dialog"/.test(attributes) ? 'dialog' : roleIn(TAG_ROLES, tag.toLowerCase()));

/**
 * What the text of the call at `index` is to a reader: in index.html the
 * attribute that holds the call, else the nearest enclosing element with a
 * role; in a module the property the text is stored under.
 * @param {string} source
 * @param {number} index
 * @param {boolean} html
 * @returns {KeyRole}
 */
export const roleAt = (source, index, html) => {
  const before = source.slice(Math.max(0, index - 4000), index);
  if (!html) {
    const property = /([A-Za-z_$][\w$]*)['"]?\s*[:=]\s*(?:[^,;{}()]*\?\s*)?$/.exec(before.slice(before.lastIndexOf('\n') + 1))?.[1] ?? '';
    return roleIn(PROPERTY_ROLES, property) ?? 'other';
  }
  /** @type {{ tag: string, attributes: string }[]} */
  const open = [];
  for (const [, end, tag, attributes, selfClosing] of before.matchAll(TAG)) {
    if (end) {
      const at = open.map((element) => element.tag).lastIndexOf(tag.toLowerCase());
      if (at >= 0) open.length = at;
    } else if (!selfClosing && !VOID_TAGS.has(tag.toLowerCase())) open.push({ tag: tag.toLowerCase(), attributes });
  }
  const lt = before.lastIndexOf('<');
  const startTag = /^<([A-Za-z][\w-]*)([\s\S]*)$/.exec(before.slice(lt));
  if (lt > before.lastIndexOf('>') && startTag) {
    // Inside a start tag: the attribute, then the element and its ancestors.
    const attribute = /(?:^|\s)(?:v-bind)?[:@]?([\w.-]+)\s*=\s*(?:"[^"]*|'[^']*)$/.exec(startTag[2])?.[1] ?? '';
    const role = roleIn(ATTRIBUTE_ROLES, attribute);
    if (role) return role;
    open.push({ tag: startTag[1].toLowerCase(), attributes: startTag[2] });
  }
  for (let at = open.length - 1; at >= 0; at -= 1) {
    const role = elementRole(open[at].tag, open[at].attributes);
    if (role) return role;
  }
  return 'other';
};

/**
 * The part of index.html each line belongs to: the settings section whose
 * element holds it (its title from SETTINGS_SECTIONS), a dialog, the shared
 * component templates, or the page.
 * @param {string} html
 * @returns {(line: number) => string}
 */
const htmlAreas = (html) => {
  /** @type {Map<string, string>} */
  const titles = new Map();
  for (const section of SETTINGS_SECTIONS) {
    titles.set(section.key, 'title' in section ? section.title : '');
  }
  for (const section of SETTINGS_SECTIONS) {
    if ('dock' in section) titles.set(section.key, titles.get(section.dock) ?? section.dock);
  }
  const lineAt = lineFinder(html);
  /** @type {{ from: number, to: number, area: string }[]} */
  const spans = [];
  const sections = /<([A-Za-z][\w-]*)\b(?:"[^"]*"|'[^']*'|[^'">])*?\sv-settings-section="'([a-z-]+)'"(?:"[^"]*"|'[^']*'|[^'">])*>/g;
  for (const start of html.matchAll(sections)) {
    const end = elementEnd(html, start[1].toLowerCase(), start.index + start[0].length);
    spans.push({ from: lineAt(start.index), to: lineAt(end), area: `Settings: ${titles.get(start[2]) ?? start[2]}` });
  }
  for (const start of html.matchAll(TAG)) {
    if (start[1] || start[4] || elementRole(start[2], start[3]) !== 'dialog') continue;
    const end = elementEnd(html, start[2].toLowerCase(), start.index + start[0].length);
    spans.push({ from: lineAt(start.index), to: lineAt(end), area: 'Dialogs' });
  }
  spans.sort((a, b) => a.from - b.from);
  const templates = html.search(/<script\b[^>]*type="text\/x-template"/);
  const templateLine = templates < 0 ? Infinity : lineAt(templates);
  return (line) => {
    if (line >= templateLine) return 'Shared components';
    const inner = spans.filter(({ from, to }) => from <= line && line <= to).pop();
    return inner ? inner.area : 'Page';
  };
};

/**
 * The line of a derived key's English text in `file`, or 1.
 * @param {string} file path under gbdraw/web
 * @param {string} english
 */
const derivedSite = (file, english) => {
  const source = readFileSync(join(WEB_ROOT, file), 'utf8');
  const at = [`'${english}'`, `"${english}"`].map((literal) => source.indexOf(literal)).find((index) => index >= 0) ?? -1;
  const line = at < 0 ? 1 : lineFinder(source)(at);
  return { file, line, snippet: source.split('\n')[line - 1].trim().slice(0, 220), area: `Script: ${file}` };
};

/**
 * Every catalog key of the sources, in source order (index.html, then the
 * modules, then the derived keys), with the places that produce it, and the
 * problems of the calls whose key cannot be read.
 * @returns {Promise<{ keys: Map<string, SourceKey>, problems: string[] }>}
 */
export const sourceKeys = async () => {
  /** @type {Map<string, SourceKey>} */
  const keys = new Map();
  /** @type {string[]} */
  const problems = [];
  /**
   * @param {string} key
   * @param {string} english
   * @param {string | null} context
   * @param {KeySite} site
   * @param {{ one: string, other: string } | null} plural
   */
  const add = (key, english, context, site, plural) => {
    const entry = keys.get(key) ?? { key, english, context, plural: null, sites: [] };
    entry.plural ??= plural;
    entry.sites.push(site);
    keys.set(key, entry);
  };
  for (const path of sourceFiles()) {
    const source = sourceText(path);
    const file = webPath(path);
    const html = file.endsWith('.html');
    const lines = source.split('\n');
    const lineAt = lineFinder(source);
    const area = html ? htmlAreas(source) : () => `Script: ${file}`;
    for (const call of extractCalls(source, file)) {
      if ('problem' in call) {
        problems.push(call.problem);
        continue;
      }
      const line = lineAt(call.index);
      const site = { file, line, role: roleAt(source, call.index, html), snippet: lines[line - 1].trim().slice(0, 220), area: area(line) };
      const context = call.name === 'tx' ? call.key.slice(0, call.key.indexOf(CONTEXT_SEPARATOR)) : null;
      add(call.key, call.english, context, site, call.plural);
    }
  }
  for (const key of await errorTextKeys()) add(key, key, null, { ...derivedSite('js/utils/error-normalization.js', key), role: 'notice' }, null);
  for (const { context, english } of memberValueLabelKeys()) {
    add(context + CONTEXT_SEPARATOR + english, english, context, { ...derivedSite('js/app/pairwise-match-popup.js', english), role: 'label' }, null);
  }
  for (const { context, english } of settingsShortKeys()) {
    add(context + CONTEXT_SEPARATOR + english, english, context, { ...derivedSite('js/services/settings-sections.js', english), role: 'label' }, null);
  }
  return { keys, problems };
};
