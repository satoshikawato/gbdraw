// The Web app reads right to left when its UI language does (setUiLocale writes the
// page's dir). A side of the page is therefore named by the reading direction: start
// and end, not left and right. Tailwind: ms-/me- for ml-/mr-, ps-/pe- for pl-/pr-,
// start-/end- for left-/right-, text-start/text-end, border-s/border-e,
// rounded-s/-e/-ss/-se/-es/-ee. CSS: margin-inline-start/end, padding-inline-*,
// inset-inline-*, border-inline-*, border-start-start-radius and so on, and
// padding-block plus padding-inline for a four-value padding or margin whose left
// and right differ. A style binding: marginInlineEnd, insetInlineStart.
//
// This guard fails when the page gains a physical side outside the lists below:
// a declaration in the generated stylesheet (so a class used in index.html or a
// module counts) or in index.html <style>, a style binding of index.html, or a
// style write of a module. Besides one-sided properties it sees a four-value
// shorthand, a border radius, a box shadow with a horizontal offset, and a left or
// right position (object-position, background-position, transform-origin), and a
// horizontal translation (transform, translate). A pair that sets left and right
// alike (px-2, mx-4, rounded-t-md) is the same in both directions. A selector under
// :dir() is the direction's own rule; the other selectors of its list are read.
// The exported interactive SVG's own stylesheet
// (services/standalone-interactivity-assets.js) is not the page and is not read.
//
// Each list entry is a physical side that stays, with its reason: an element
// placed at viewport or figure coordinates, or centred on a point, is the same in
// both directions. An entry exempts the declarations it names, not the rest of its
// rule or binding. Adding an entry, or raising a count, widens what the guard
// allows; no design rule registers these lists yet
// (tools/web-design-rule-guards.json), so review such a change as one that
// widens a guard.
import assert from 'node:assert/strict';
import { readdirSync, readFileSync } from 'node:fs';
import { join, relative } from 'node:path';
import test from 'node:test';
import { fileURLToPath } from 'node:url';

const repository = fileURLToPath(new URL('../../', import.meta.url));
const web = join(repository, 'gbdraw/web');
const read = (path) => readFileSync(join(repository, path), 'utf8');

/** CSS rules (selector) whose listed physical declarations stay. */
const CSS_EXCEPTIONS = new Map([
  ['.settings-index-tip', [['transform: translateX(-50%)'], 'placed at viewport coordinates from its button; translateX(-50%) centres it on that point']],
  ['.settings-index-tip[data-side="left"]', [['transform: translate(-100%, -50%)'], 'beside a rail button in a right-to-left page: it ends at its point']],
  ['.left-1\\/2', [['left: 50%'], 'centres an element with -translate-x-1/2, the same in both directions (a class list has both: last test)']],
  ['.-translate-x-1\\/2', [['--tw-translate-x: -50%'], 'centres an element with left-1/2, the same in both directions']],
  ['.right-drawer[aria-hidden="true"]', [['transform: translateX(100%)'], 'slides the Editor out past the inline end; the :dir(rtl) rule beside it slides it the other way']],
  ['.drawer-toggle[aria-expanded="true"]', [['transform: translate( calc(-1 * var(--right-drawer-width)), -50% )'], 'moves the toggle with the open Editor; the :dir(rtl) rule beside it moves it the other way']],
  ['.history-list-row[data-place="current"]', [['box-shadow: inset 3px 0 0 #2563eb'], 'marks the current History step at the row\'s start; the :dir(rtl) rule beside it marks the other side']],
  ['.feature-popup-resize-handle', [['right: 0.25rem'], 'the resize grip of a popup at viewport coordinates: a drag widens the popup to the right from its fixed left edge (app/result-popups.js), so the grip stays bottom right']],
  ['.feature-popup-resize-handle::before', [['right: 0.35rem', 'border-right: 2px solid currentColor', 'box-shadow: 4px 4px 0 -2px currentColor'], 'the mark of that bottom-right grip']],
  ['#app-boot-splash .boot-wordmark', [['object-position: right'], 'chooses the part of the wordmark image that shows; an image does not mirror']]
]);

/** index.html style bindings: the physical side declarations that stay (the start of `key: value`). */
const BINDING_EXCEPTIONS = new Map([
  ['left: historyList.place.value.left', 'the History list, at viewport coordinates from its toggle (app/history-list.js)'],
  ['left: `${regionOverlay.readoutX', 'the region readout, beside the pointer'],
  ['left: tip.left', 'the text of a settings index button, at viewport coordinates from the button']
]);

/** Modules under gbdraw/web/js that write a physical side, with how many writes. */
const MODULE_EXCEPTIONS = new Map([
  ['app/result-popups.js', [1, 'a popup at the clicked point, kept inside the viewport (utils/popup-bounds.js)']],
  ['app/feature-selection.js', [1, 'the selection rectangle under the pointer']],
  ['app/feature-editor/svg-actions.js', [3, 'the alignment overlay rings and labels, at figure coordinates (the overlay reads left to right)']],
  ['app/feature-editor/ring-hover.js', [1, 'a ring popup beside the pointer']],
  ['services/export.js', [1, 'the PDF staging box, above and left of the page, outside its scroll area in either direction']],
  ['utils/tooltip-placement.js', [1, 'the house tooltip (a help tip, the whole hint of a cut field), centred on its anchor at viewport coordinates']],
  ['utils/clipboard.js', [1, 'the copy fallback text area, fixed off the left edge: a fixed box adds no scroll area in either direction']]
]);

const SIDE_PAIRS = [
  ['left', 'right'], ['margin-left', 'margin-right'], ['padding-left', 'padding-right'],
  ['scroll-margin-left', 'scroll-margin-right'], ['scroll-padding-left', 'scroll-padding-right'],
  ['border-top-left-radius', 'border-top-right-radius'], ['border-bottom-left-radius', 'border-bottom-right-radius']
];
const PHYSICAL_PROPERTY = /^(?:left|right|(?:margin|padding|scroll-margin|scroll-padding)-(?:left|right)|border-(?:left|right)(?:-[a-z]+)?|border-(?:top|bottom)-(?:left|right)-radius)$/;
const ZERO = /^-?0(?:px|%|rem)?$/;
const FOUR_SIDES = new Set(['margin', 'padding', 'inset', 'scroll-margin', 'scroll-padding', 'border-width', 'border-style', 'border-color']);
const LENGTH = /^-?[\d.]+(?:px|rem|em|%)?$/;

/** The space-separated parts of a value, keeping parentheses whole. */
const valueParts = (value) => {
  const parts = [];
  let depth = 0;
  let part = '';
  for (const character of value) {
    if (character === '(') depth += 1;
    if (character === ')') depth -= 1;
    if (/\s/.test(character) && depth === 0) {
      if (part) parts.push(part);
      part = '';
    } else part += character;
  }
  return part ? [...parts, part] : parts;
};

/** Whether a border radius rounds a left corner unlike the right one beside it. */
const unevenCorners = (value) => {
  const parts = valueParts(value.split('/')[0]);
  const [topLeft, topRight, bottomRight, bottomLeft] = parts.length === 2 ? [parts[0], parts[1], parts[0], parts[1]]
    : parts.length === 3 ? [parts[0], parts[1], parts[2], parts[1]] : parts;
  return parts.length >= 2 && (topLeft !== topRight || bottomLeft !== bottomRight);
};

/** Whether a shadow has a layer moved left or right. */
const sidewaysShadow = (value) => value.split(/,(?![^(]*\))/).some((layer) => {
  const [x] = valueParts(layer.trim()).filter((part) => part !== 'inset');
  return x !== undefined && LENGTH.test(x) && !ZERO.test(x);
});

/** The first argument of a translate(...) or translateX(...) call, at paren depth 0. */
const translateArguments = (value) => {
  const found = [];
  for (const match of value.matchAll(/translateX?\(/g)) {
    let depth = 0;
    let argument = '';
    for (const character of value.slice(match.index + match[0].length)) {
      if (depth === 0 && (character === ',' || character === ')')) break;
      if (character === '(') depth += 1;
      if (character === ')') depth -= 1;
      argument += character;
    }
    found.push(argument.trim());
  }
  return found;
};

/**
 * The physical sides a rule body sets: one-sided or unequal pairs, a horizontal
 * translation, a left or right text alignment, float or clear.
 * @param {string} body
 * @returns {string[]}
 */
const physicalDeclarations = (body) => {
  const declarations = new Map();
  for (const part of body.split(';')) {
    const colon = part.indexOf(':');
    if (colon < 0) continue;
    const property = part.slice(0, colon).trim().toLowerCase();
    const value = part.slice(colon + 1).replace(/!important/, '').trim().replace(/\s+/g, ' ');
    if (property.startsWith('@')) continue;
    declarations.set(property, value);
  }
  const found = [];
  const paired = new Set();
  for (const [a, b] of [...SIDE_PAIRS, ...['width', 'style', 'color'].map((x) => [`border-left-${x}`, `border-right-${x}`]), ['border-left', 'border-right']]) {
    if (declarations.has(a) && declarations.has(b) && declarations.get(a) === declarations.get(b)) paired.add(a).add(b);
  }
  for (const [property, value] of declarations) {
    if (PHYSICAL_PROPERTY.test(property) && !paired.has(property)) found.push(`${property}: ${value}`);
    if ((property === 'text-align' || property === 'float' || property === 'clear') && /\b(?:left|right)\b/.test(value)) found.push(`${property}: ${value}`);
    if (property === '--tw-translate-x' && !ZERO.test(value)) found.push(`${property}: ${value}`);
    if (FOUR_SIDES.has(property) && valueParts(value).length === 4 && valueParts(value)[1] !== valueParts(value)[3]) found.push(`${property}: ${value}`);
    if (property === 'border-radius' && unevenCorners(value)) found.push(`${property}: ${value}`);
    if ((property === 'box-shadow' || property === '--tw-shadow') && sidewaysShadow(value)) found.push(`${property}: ${value}`);
    if (['object-position', 'background-position', 'background-position-x', 'transform-origin'].includes(property) && /\b(?:left|right)\b/.test(value)) found.push(`${property}: ${value}`);
    if (property === 'transform' && translateArguments(value).some((x) => !ZERO.test(x) && x !== 'var(--tw-translate-x)')) found.push(`${property}: ${value}`);
    if (property === 'translate' && !ZERO.test(valueParts(value)[0] ?? '0') && value !== 'none') found.push(`${property}: ${value}`);
  }
  return found;
};

/**
 * Each innermost rule of a stylesheet, with its selector.
 * @param {string} css
 * @returns {{ selector: string, body: string }[]}
 */
const innermostRules = (css) => {
  const text = css.replace(/\/\*[\s\S]*?\*\//g, '');
  const rules = [];
  const opens = [];
  for (let i = 0; i < text.length; i += 1) {
    if (text[i] === '{') opens.push(i);
    if (text[i] !== '}' || opens.length === 0) continue;
    const open = opens.pop();
    const body = text.slice(open + 1, i);
    if (body.includes('{')) continue;
    let start = open - 1;
    while (start >= 0 && !'{};'.includes(text[start])) start -= 1;
    rules.push({ selector: text.slice(start + 1, open).trim().replace(/\s+/g, ' '), body });
  }
  return rules;
};

const html = read('gbdraw/web/index.html');
const styleBlocks = [...html.matchAll(/<style\b[^>]*>([\s\S]*?)<\/style>/g)].map((match) => match[1]);

/** The parts of a list separated by commas at parenthesis and brace depth 0. */
const topLevelParts = (text) => {
  const parts = [''];
  let depth = 0;
  for (const character of text) {
    if ('({['.includes(character)) depth += 1;
    if (')}]'.includes(character)) depth -= 1;
    if (character === ',' && depth === 0) parts.push('');
    else parts[parts.length - 1] += character;
  }
  return parts.map((part) => part.trim()).filter(Boolean);
};

/** A selector list without its selectors under :dir(), which are the direction's own rules. */
const selectorsOutsideDir = (selector) => topLevelParts(selector).filter((part) => !part.includes(':dir(')).join(', ');

/** @returns {{ where: string, selector: string, sides: string[] }[]} */
const cssFindings = () => [
  ['gbdraw/web/css/tailwind.css', read('gbdraw/web/css/tailwind.css')],
  ...styleBlocks.map((block) => ['gbdraw/web/index.html <style>', block])
].flatMap(([where, css]) => innermostRules(css)
  .map(({ selector, body }) => ({ where, selector: selectorsOutsideDir(selector), sides: physicalDeclarations(body) }))
  .filter(({ selector, sides }) => selector && sides.length > 0));

/** A finding's declarations that its CSS_EXCEPTIONS entry does not name. */
const unlistedSides = ({ selector, sides }) => sides.filter((side) => !(CSS_EXCEPTIONS.get(selector)?.[0] ?? []).includes(side));

const BINDING_SIDE = /^(?:left|right|marginLeft|marginRight|paddingLeft|paddingRight|margin-left|margin-right|padding-left|padding-right|['"](?:left|right|margin-left|margin-right|padding-left|padding-right)['"])\s*:/;
/** Each physical side declaration (`key: value`) of the style bindings of a page. */
const bindingSides = (page) => [...page.matchAll(/(?<![\w:-])(?::style|v-bind:style|style)="([^"]*)"/g)]
  .flatMap((match) => topLevelParts(match[1].trim().replace(/^\{([\s\S]*)\}$/, '$1')))
  .filter((declaration) => BINDING_SIDE.test(declaration));
const bindingFindings = () => bindingSides(html);

const MODULE_WRITES = [
  /\bstyle\.(?:left|right|marginLeft|marginRight|paddingLeft|paddingRight|borderLeft\w*|borderRight\w*|translate)\s*=(?!=)/g,
  /\bsetProperty\(\s*['"](?:left|right|translate|(?:margin|padding|border)-(?:left|right)[\w-]*)['"]/g,
  /\b(?:left|right|marginLeft|marginRight|paddingLeft|paddingRight|translate)\s*:\s*(?:`[^`\n]*`|'[^'\n]*'|"[^"\n]*"|[^,}\n]*\+\s*'px')/g,
  /\bstyle\.transform\s*=\s*[`'"][^`'"\n]*translate/g,
  /\bstyle\.cssText\s*=\s*[`'"][^`'"\n]*?(?:^|[;`'"\s])(?:left|right|translate|(?:margin|padding|border)-(?:left|right)[\w-]*)\s*:/g
];
const moduleFiles = (directory) => readdirSync(directory, { withFileTypes: true }).flatMap((entry) => {
  const path = join(directory, entry.name);
  if (entry.isDirectory()) return moduleFiles(path);
  return entry.name.endsWith('.js') ? [path] : [];
});
/** @returns {Map<string, number>} */
const moduleWrites = () => {
  const counts = new Map();
  for (const path of moduleFiles(join(web, 'js'))) {
    const name = relative(join(web, 'js'), path).split('\\').join('/');
    if (name === 'services/standalone-interactivity-assets.js') continue;
    const source = readFileSync(path, 'utf8');
    const count = MODULE_WRITES.reduce((sum, pattern) => sum + [...source.matchAll(pattern)].length, 0);
    if (count > 0) counts.set(name, count);
  }
  return counts;
};

const HINT = 'name the side by the reading direction (start/end; see the head of tests/web/physical-sides.test.mjs)';

test('the page styles name a side by the reading direction, except the listed declarations', () => {
  const unlisted = cssFindings().map((finding) => ({ ...finding, sides: unlistedSides(finding) })).filter(({ sides }) => sides.length > 0);
  assert.deepEqual(unlisted, [], `physical sides in CSS: ${HINT}`);
  const used = new Set(cssFindings().flatMap(({ selector, sides }) => sides.map((side) => `${selector} { ${side} }`)));
  assert.deepEqual([...CSS_EXCEPTIONS].flatMap(([selector, [sides]]) => sides.map((side) => `${selector} { ${side} }`)).filter((entry) => !used.has(entry)), [],
    'a listed declaration is no longer in its rule: remove it from CSS_EXCEPTIONS');
});

/** Whether a binding declaration is one BINDING_EXCEPTIONS names. */
const listedBinding = (declaration) => [...BINDING_EXCEPTIONS.keys()].some((key) => declaration.startsWith(key));

test('the style bindings of index.html name a side by the reading direction, except the listed ones', () => {
  const findings = bindingFindings();
  assert.deepEqual(findings.filter((declaration) => !listedBinding(declaration)), [], `physical sides in style bindings: ${HINT}`);
  assert.deepEqual([...BINDING_EXCEPTIONS.keys()].filter((key) => !findings.some((declaration) => declaration.startsWith(key))), [],
    'a listed binding no longer sets a physical side: remove it from BINDING_EXCEPTIONS');
});

test('a class list that centres with left-1/2 also has -translate-x-1/2, and the other way round', () => {
  const sources = [html, ...moduleFiles(join(web, 'js')).map((path) => readFileSync(path, 'utf8'))];
  const classLists = sources.flatMap((source) => [...source.matchAll(/\bclass="([^"]*)"/g)].map((match) => match[1].split(/\s+/)));
  const halves = classLists.filter((classes) => classes.includes('left-1/2') !== classes.includes('-translate-x-1/2'));
  assert.deepEqual(halves, [], 'left-1/2 without -translate-x-1/2 (or the other way round) sets a physical side');
});

test('the modules write a physical side only where the list says, as often as it says', () => {
  const observed = Object.fromEntries(moduleWrites());
  const allowed = Object.fromEntries([...MODULE_EXCEPTIONS].map(([name, [count]]) => [name, count]));
  assert.deepEqual(observed, allowed, `physical style writes in modules: ${HINT}; a lower count lowers its MODULE_EXCEPTIONS entry`);
});

test('the guard sees the physical sides it looks for', () => {
  assert.deepEqual(physicalDeclarations('margin-left: 1rem; color: red'), ['margin-left: 1rem']);
  assert.deepEqual(physicalDeclarations('padding-left: 1rem; padding-right: 1rem'), []);
  assert.deepEqual(physicalDeclarations('left: 0; right: 4px'), ['left: 0', 'right: 4px']);
  assert.deepEqual(physicalDeclarations('text-align: right'), ['text-align: right']);
  assert.deepEqual(physicalDeclarations('transform: translate(calc(-1 * var(--w)), -50%)'), ['transform: translate(calc(-1 * var(--w)), -50%)']);
  assert.deepEqual(physicalDeclarations('transform: translateX(0); transform: translate(var(--tw-translate-x), var(--tw-translate-y))'), []);
  assert.deepEqual(physicalDeclarations('margin-inline-start: 1rem; inset-inline-end: 0; text-align: start'), []);
  assert.deepEqual(physicalDeclarations('padding: 0 0.25rem 1.5rem 0; margin: 0 4px 3px 4px'), ['padding: 0 0.25rem 1.5rem 0']);
  assert.deepEqual(physicalDeclarations('padding: 6px 0 3px'), []);
  assert.deepEqual(physicalDeclarations('border-radius: 0 3px 3px 0'), ['border-radius: 0 3px 3px 0']);
  assert.deepEqual(physicalDeclarations('border-radius: 8px; border-radius: 4px 4px 0 0'), []);
  assert.deepEqual(physicalDeclarations('box-shadow: inset 3px 0 0 #2563eb'), ['box-shadow: inset 3px 0 0 #2563eb']);
  assert.deepEqual(physicalDeclarations('box-shadow: 0 4px 10px rgba(15, 23, 42, 0.25), inset 0 0 0 1px #fff'), []);
  assert.deepEqual(physicalDeclarations('object-position: right; transform-origin: center top'), ['object-position: right']);
  assert.deepEqual(physicalDeclarations('translate: 4px 0; scale: -1 1'), ['translate: 4px 0']);
  assert.deepEqual(physicalDeclarations('translate: 0 -50%'), []);
  assert.deepEqual(physicalDeclarations('border-top-left-radius: 4px; border-top-right-radius: 4px'), [], 'a symmetric rounded-t is the same in both directions');
  assert.deepEqual(physicalDeclarations('border-top-left-radius: 4px'), ['border-top-left-radius: 4px']);
  assert.deepEqual(selectorsOutsideDir('.a:dir(rtl), .b, :is(.c, .d):dir(rtl)'), '.b');
  assert.deepEqual(unlistedSides({ selector: '.settings-index-tip', sides: ['transform: translateX(-50%)', 'margin-left: 4px'] }), ['margin-left: 4px'],
    'an exception exempts only the declarations it names');
  assert.deepEqual(bindingSides(`<div :style="{ left: tip.left + 'px', marginRight: '4px', top: y }">`).filter((declaration) => !listedBinding(declaration)), ["marginRight: '4px'"],
    'a second side in a listed binding is seen');
  assert.deepEqual(innermostRules('@media (x) { .a { left: 1px } } /* .b { left: 2px } */ .c { right: 0 }')
    .map(({ selector }) => selector), ['.a', '.c']);
  const write = (source) => MODULE_WRITES.reduce((sum, pattern) => sum + [...source.matchAll(pattern)].length, 0);
  assert.equal(write("el.style.left = '4px'; el.style.insetInlineStart = '4px'; if (el.style.left == x) {}"), 1);
  assert.equal(write("Object.assign(el.style, { left: x + 'px', top: y + 'px' }); const box = { left: rect.left };"), 1);
  assert.equal(write("el.style.setProperty('margin-right', '2px')"), 1);
  assert.equal(write("Object.assign(el.style, { right: 'auto', left: `calc(${x}px - 1rem)` })"), 2);
  assert.equal(write("el.style.transform = 'translateX(4px)'; el.style.translate = '4px'; el.style.transform = transform;"), 2);
  assert.equal(write("el.style.cssText = 'position:absolute;left:0'; probe.style.cssText = 'height:0;width:4px'"), 1);
});

// PD-OI-106: identifiers, coordinates, sequences and colors read left to right in every language. A
// font-mono element shows one of them, except those that show a translated sentence.
const MONO_SENTENCES = ['data-generation-progress', '<ul'];
test('every font-mono element of index.html and the module templates reads left to right, except the translated sentences', () => {
  const sources = [['gbdraw/web/index.html', html], ...moduleFiles(join(web, 'js'))
    .filter((path) => !path.endsWith('standalone-interactivity-assets.js'))
    .map((path) => [relative(repository, path), readFileSync(path, 'utf8')])];
  const unpinned = sources.flatMap(([where, source]) => [...source.matchAll(/<[a-z][\w-]*\s[^<>]*?\bfont-mono\b[^<>]*>/g)]
    .map((match) => match[0])
    .filter((tag) => !/\s:?dir=/.test(tag) && !MONO_SENTENCES.some((mark) => tag.includes(mark)))
    .map((tag) => `${where}: ${tag.slice(0, 120)}`));
  assert.deepEqual(unpinned, [], 'give the element dir="ltr" (or :dir when it may show a translated word)');
});
