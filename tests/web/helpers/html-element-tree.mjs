// The element tree of an HTML text, for the node tests that read index.html
// statically (tests/web/setting-application.test.mjs and
// tests/web/live-random-walk-settings.test.mjs): one scanner, so the two read
// the same elements, attributes and parents. Quoted attribute values may hold
// `<` and `>` (Vue expressions); comments and the text of script and style
// elements hold no elements (a component template in a
// `<script type="text/x-template">` is read by its own test).

/**
 * An element: its lower-case tag, its attributes as written (a bound one keeps
 * its `:` or `@`), its parent (the root for a top-level element), its children
 * in document order, and the line of its opening tag.
 * @typedef {{ tag: string, attributes: Map<string, string>, parent: HtmlElement | null, children: HtmlElement[], line: number }} HtmlElement
 */

const VOID = new Set(['area', 'base', 'br', 'col', 'embed', 'hr', 'img', 'input', 'link', 'meta', 'param', 'source', 'track', 'wbr']);
const RAW_TEXT = new Set(['script', 'style']);

/** @param {string} text @returns {Map<string, string>} */
const parseAttributes = (text) => {
  const attributes = new Map();
  for (const match of text.matchAll(/([^\s=/"']+)(?:\s*=\s*("[^"]*"|'[^']*'|[^\s>"']+))?/g)) {
    const value = match[2] ?? '';
    attributes.set(match[1], value.replace(/^["']|["']$/g, ''));
  }
  return attributes;
};

/**
 * The tree of `html` under a root element `#root` (no attributes, line 0).
 * @param {string} html
 * @returns {HtmlElement}
 */
export const parseElementTree = (html) => {
  /** @param {number} offset */
  const lineOf = (offset) => html.slice(0, offset).split('\n').length;
  // The end of an opening tag, skipping '>' inside quoted attribute values.
  /** @param {number} from */
  const tagEnd = (from) => {
    /** @type {string | null} */
    let quote = null;
    for (let index = from; index < html.length; index += 1) {
      const char = html[index];
      if (quote) {
        if (char === quote) quote = null;
      } else if (char === '"' || char === "'") {
        quote = char;
      } else if (char === '>') {
        return index;
      }
    }
    throw new Error(`Unterminated tag at line ${lineOf(from)}`);
  };
  /** @type {HtmlElement} */
  const root = { tag: '#root', attributes: new Map(), parent: null, children: [], line: 0 };
  let current = root;
  let index = 0;
  while (index < html.length) {
    const open = html.indexOf('<', index);
    if (open < 0) break;
    if (html.startsWith('<!--', open)) {
      index = html.indexOf('-->', open) + 3;
      continue;
    }
    if (html.startsWith('</', open)) {
      const end = html.indexOf('>', open);
      const tag = html.slice(open + 2, end).trim().toLowerCase();
      /** @type {HtmlElement | null} */
      let node = current;
      while (node && node.tag !== tag) node = node.parent;
      if (node?.parent) current = node.parent;
      index = end + 1;
      continue;
    }
    if (!/[a-zA-Z]/.test(html[open + 1] || '')) {
      index = open + 1;
      continue;
    }
    const end = tagEnd(open);
    const raw = html.slice(open + 1, end);
    const tag = /** @type {RegExpMatchArray} */ (raw.match(/^[^\s/>]+/))[0].toLowerCase();
    /** @type {HtmlElement} */
    const element = { tag, attributes: parseAttributes(raw.slice(tag.length)), parent: current, children: [], line: lineOf(open) };
    current.children.push(element);
    index = end + 1;
    if (RAW_TEXT.has(tag)) {
      const close = html.indexOf(`</${tag}>`, index);
      index = close < 0 ? html.length : close + `</${tag}>`.length;
      continue;
    }
    if (!VOID.has(tag) && !raw.endsWith('/')) current = element;
  }
  return root;
};

/** Every element under `element`, in document order. @param {HtmlElement} element @returns {Generator<HtmlElement>} */
export const descendants = function* descendants(element) {
  for (const child of element.children) {
    yield child;
    yield* descendants(child);
  }
};

/** `element`, its parent, and so on up to the root. @param {HtmlElement} element @returns {Generator<HtmlElement>} */
export const ancestorsOrSelf = function* ancestorsOrSelf(element) {
  for (let node = /** @type {HtmlElement | null} */ (element); node; node = node.parent) yield node;
};
