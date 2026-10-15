// OV-904: two Tailwind utilities that set the same property on one element (a
// static `bg-white` and a bound `bg-blue-50`) are decided by their order in
// css/tailwind.css, and that order follows the order in which the compiler met
// the classes, not the binding. No element of index.html may apply two such
// utilities at once: the binding gives each state one of them.
//
// The test reads the `class` and `:class` attributes of every element, treats
// the two branches of one ternary as exclusive and everything else as applied
// together, and finds pairs whose rules in tailwind.css set the same property
// in the same selector context. A pair is order-dependent when compiling the
// two classes in the other order swaps their rules; pairs whose order Tailwind
// fixes (different plugins or variants) never swap and pass.
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import test from 'node:test';

import {
  compileTailwindCss,
  TAILWIND_COMPILER_PATH,
  TAILWIND_INPUT_PATH,
  TAILWIND_OUTPUT_PATH
} from '../../tools/generate_web_tailwind_css.mjs';

const repositoryFile = (path) => readFileSync(new URL(`../../${path}`, import.meta.url), 'utf8');

// OV-905 pairs, kept until the Owner decides D-09; the value is the number of
// elements that carry the pair. Remove an entry when its elements change.
const PENDING_OWNER_DECISION = new Map([
  ['bg-amber-50 bg-white', 1], // Layout edit toggle, active
  ['border-amber-300 border-slate-200', 1], // Layout edit toggle, active
  ['text-amber-700 text-slate-600', 1], // Layout edit toggle, active
  ['text-blue-600 text-slate-600', 2] // Canvas padding and Search toggles, active
]);

// ---- class bindings ---------------------------------------------------------

const CLOSING = { '(': ')', '[': ']', '{': '}' };

// Visit each character of a JavaScript expression outside string literals with
// its bracket depth; stop when the visitor returns false.
const scan = (text, visit) => {
  let depth = 0;
  for (let i = 0; i < text.length; i += 1) {
    const char = text[i];
    if (char === "'" || char === '"' || char === '`') {
      for (i += 1; i < text.length && text[i] !== char; i += 1) if (text[i] === '\\') i += 1;
      continue;
    }
    if (CLOSING[char]) depth += 1;
    else if (char === ')' || char === ']' || char === '}') depth -= 1;
    if (visit(i, char, depth) === false) return;
  }
};
const splitTopLevel = (text, separator) => {
  const parts = [];
  let start = 0;
  scan(text, (i, char, depth) => {
    if (depth === 0 && char === separator) {
      parts.push(text.slice(start, i));
      start = i + 1;
    }
  });
  return [...parts, text.slice(start)];
};
const enclosedBy = (text, open) => {
  if (text[0] !== open) return false;
  let end = -1;
  scan(text, (i, char, depth) => {
    if (depth === 0) {
      end = i;
      return false;
    }
    return true;
  });
  return end === text.length - 1;
};
const isTernaryMark = (text, i) => text[i] === '?' && text[i + 1] !== '.' && text[i + 1] !== '?' && text[i - 1] !== '?';
const ternaryParts = (text) => {
  let question = -1;
  let colon = -1;
  let nested = 0;
  scan(text, (i, char, depth) => {
    if (depth !== 0) return true;
    if (question < 0) {
      if (isTernaryMark(text, i)) question = i;
      return true;
    }
    if (isTernaryMark(text, i)) nested += 1;
    else if (char === ':' && nested-- === 0) {
      colon = i;
      return false;
    }
    return true;
  });
  return colon < 0 ? null : [text.slice(question + 1, colon), text.slice(colon + 1)];
};
const stringTokens = (text) => [...text.matchAll(/'([^']*)'|"([^"]*)"|`([^`]*)`/g)]
  .flatMap((match) => (match[1] ?? match[2] ?? match[3]).replace(/\$\{[^}]*\}/g, ' ').split(/\s+/))
  .filter(Boolean);

// Each token of a `:class` value with the ternary branches it sits in.
const bindingTokens = (expression, branches, out, ids) => {
  const text = expression.trim();
  if (!text) return;
  if (enclosedBy(text, '(')) return bindingTokens(text.slice(1, -1), branches, out, ids);
  if (enclosedBy(text, '[')) {
    splitTopLevel(text.slice(1, -1), ',').forEach((item) => bindingTokens(item, branches, out, ids));
    return;
  }
  if (enclosedBy(text, '{')) {
    for (const entry of splitTopLevel(text.slice(1, -1), ',')) {
      const key = splitTopLevel(entry, ':')[0].trim();
      const names = /^['"`]/.test(key) ? stringTokens(key) : [key];
      names.filter(Boolean).forEach((token) => out.push({ token, branches }));
    }
    return;
  }
  const ternary = ternaryParts(text);
  if (ternary) {
    const id = ids.next += 1;
    ternary.forEach((branch, side) => bindingTokens(branch, [...branches, `${id}:${side}`], out, ids));
    return;
  }
  stringTokens(text).forEach((token) => out.push({ token, branches }));
};
const exclusive = (a, b) => a.branches.some((branch) => {
  const [id, side] = branch.split(':');
  return b.branches.includes(`${id}:${1 - Number(side)}`);
});

const TAG = /<!--[\s\S]*?-->|<([a-zA-Z][\w-]*)((?:\s+[^\s=>/"']+(?:\s*=\s*(?:"[^"]*"|'[^']*'|[^\s>"']+))?)*)\s*\/?>/g;
const ATTRIBUTE = /([^\s=>/"']+)(?:\s*=\s*(?:"([^"]*)"|'([^']*)'|([^\s>"']+)))?/g;

const classedElements = (html) => {
  const ids = { next: 0 };
  const elements = [];
  for (const match of html.matchAll(TAG)) {
    if (!match[1]) continue;
    const tokens = [];
    for (const attribute of match[2].matchAll(ATTRIBUTE)) {
      const value = (attribute[2] ?? attribute[3] ?? attribute[4] ?? '').replace(/&quot;/g, '"');
      if (attribute[1] === 'class') value.split(/\s+/).filter(Boolean).forEach((token) => tokens.push({ token, branches: [] }));
      else if (attribute[1] === ':class' || attribute[1] === 'v-bind:class') bindingTokens(value, [], tokens, ids);
    }
    if (tokens.length > 1) elements.push({ line: html.slice(0, match.index).split('\n').length, tokens });
  }
  return elements;
};

// ---- rules ------------------------------------------------------------------

const unescapeClass = (text) => text
  .replace(/\\([0-9a-fA-F]{1,6}) ?/g, (_, hex) => String.fromCodePoint(Number.parseInt(hex, 16)))
  .replace(/\\(.)/g, '$1');
const CLASS_SELECTOR = /\.((?:\\[0-9a-fA-F]{1,6} ?|\\.|[\w-])+)/g;

// token -> Map(context of selector, at-rules, property, and importance -> rule index)
const ruleIndex = (css) => {
  const index = new Map();
  let rule = 0;
  const walk = (text, atRules) => {
    let i = 0;
    while (i < text.length) {
      const open = text.indexOf('{', i);
      if (open < 0) return;
      const prelude = text.slice(i, open).trim();
      let depth = 1;
      let end = open + 1;
      for (; depth && end < text.length; end += 1) depth += text[end] === '{' ? 1 : text[end] === '}' ? -1 : 0;
      const body = text.slice(open + 1, end - 1);
      if (prelude.startsWith('@')) walk(body, `${atRules}${prelude} `);
      else {
        rule += 1;
        const declarations = body.split(';').map((d) => d.trim()).filter((d) => d && !d.startsWith('--'));
        for (const selector of prelude.split(/,(?![^(]*\))/)) {
          for (const match of selector.matchAll(CLASS_SELECTOR)) {
            const keys = index.get(unescapeClass(match[1])) ?? new Map();
            index.set(unescapeClass(match[1]), keys);
            const context = atRules + selector.replace(match[0], '&').trim();
            for (const declaration of declarations) {
              const property = declaration.slice(0, declaration.indexOf(':')).trim();
              const key = `${context} ${property}${/!important$/.test(declaration) ? ' !important' : ''}`;
              if (!keys.has(key)) keys.set(key, rule);
            }
          }
        }
      }
      i = end;
    }
  };
  walk(css.replace(/\/\*[\s\S]*?\*\//g, ''), '');
  return index;
};

// ---- conflicts --------------------------------------------------------------

const findOrderConflicts = async ({ html, css, components, compilerSource }) => {
  const rules = ruleIndex(css);
  const pairs = new Map();
  for (const { line, tokens } of classedElements(html)) {
    for (let x = 0; x < tokens.length; x += 1) {
      for (let y = x + 1; y < tokens.length; y += 1) {
        const [a, b] = [tokens[x], tokens[y]];
        if (a.token === b.token || exclusive(a, b) || components.has(a.token) || components.has(b.token)) continue;
        const keysA = rules.get(a.token);
        const keysB = rules.get(b.token);
        const shared = keysA && keysB ? [...keysA.keys()].filter((key) => keysB.has(key)) : [];
        if (!shared.length) continue;
        const name = [a.token, b.token].sort().join(' ');
        const pair = pairs.get(name) ?? { name, shared, lines: new Set() };
        pair.lines.add(line);
        pairs.set(name, pair);
      }
    }
  }
  if (!pairs.size) return [];
  const tokens = [...new Set([...pairs.keys()].flatMap((name) => name.split(' ')))].sort();
  const compile = (batches) => compileTailwindCss({ compilerSource, componentRules: '', batches }).then(ruleIndex);
  const forward = await compile([tokens]);
  const backward = await compile([...tokens].reverse().map((token) => [token]));
  const order = (index, a, b, key) => Math.sign(index.get(a).get(key) - index.get(b).get(key));
  return [...pairs.values()]
    .filter(({ name, shared }) => {
      const [a, b] = name.split(' ');
      return shared.some((key) => order(forward, a, b, key) !== order(backward, a, b, key));
    })
    .map(({ name, shared, lines }) => {
      const [a, b] = name.split(' ');
      return { name, lines: [...lines], winner: order(rules, a, b, shared[0]) > 0 ? a : b };
    });
};

const compilerSource = repositoryFile(TAILWIND_COMPILER_PATH);
const componentNames = (css) => new Set([...css.matchAll(/\.([a-zA-Z][\w-]*)/g)].map((match) => match[1]));

test('no element of index.html applies two order-dependent utilities for one property at once', async () => {
  const conflicts = await findOrderConflicts({
    html: repositoryFile('gbdraw/web/index.html'),
    css: repositoryFile(TAILWIND_OUTPUT_PATH),
    components: componentNames(repositoryFile(TAILWIND_INPUT_PATH)),
    compilerSource
  });
  const pending = new Map([...PENDING_OWNER_DECISION.keys()].map((name) => [name, 0]));
  const unexpected = [];
  for (const conflict of conflicts) {
    if (pending.has(conflict.name)) pending.set(conflict.name, conflict.lines.length);
    else unexpected.push(`${conflict.name} on index.html:${conflict.lines.join(',')} (${conflict.winner} wins today)`);
  }
  assert.deepEqual(unexpected, [], 'give each state one of the two classes in the binding');
  assert.deepEqual(pending, PENDING_OWNER_DECISION, 'update PENDING_OWNER_DECISION to the elements that still carry each pair');
});

test('the order guard tells bound conflicts from exclusive branches and fixed-order pairs', async () => {
  const html = [
    '<a class="bg-white" :class="on ? \'bg-blue-50\' : \'\'"></a>',
    '<b class="p-2" :class="on ? \'bg-blue-50\' : \'bg-white\'"></b>',
    '<i :class="[\'cursor-pointer\', { \'cursor-wait\': busy }]"></i>',
    '<s class="px-3 py-1 p-2 !grid-cols-1 grid-cols-2 hover:bg-transparent bg-slate-50"></s>'
  ].join('\n');
  const tokens = ['bg-white', 'bg-blue-50', 'p-2', 'px-3', 'py-1', 'cursor-pointer', 'cursor-wait', '!grid-cols-1', 'grid-cols-2', 'hover:bg-transparent', 'bg-slate-50'];
  const css = await compileTailwindCss({ compilerSource, componentRules: '', batches: [tokens] });
  const conflicts = await findOrderConflicts({ html, css, components: new Set(), compilerSource });
  assert.deepEqual(conflicts.map(({ name, lines }) => `${name} @${lines}`), ['bg-blue-50 bg-white @1', 'cursor-pointer cursor-wait @3']);
});
