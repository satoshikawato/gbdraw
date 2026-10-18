// The slice write scan (R15; D-15-44: the rows of a list slice are the
// slice's too): the writes of a drawing's editor members outside their slice,
// in a module's source and in the template. tests/web/drawing-slices.test.mjs
// runs it over every Web module and index.html; tests/web/mode-switch-isolation.test.mjs
// over the watchers. The type checker sees a write through a typed drawing,
// whose views are read-only. This scan covers the untyped holders: in each
// module outside services/drawing/, a write through `<drawing>.<member>` (an
// assignment to it, to its `.value`, an element, or a field, a mutating
// method, `++` or `--`, `delete`, or `Object.assign` into it) whose receiver
// is a drawing by name (an identifier that ends in `drawing` or `Drawing`,
// `activeDrawing()`, `drawings.<mode>`, or `drawings[<mode>]`), and any write
// through a computed key on such a receiver (`drawing[key]`,
// `drawings[<mode>][key]`, which may name a member); and in the template, in
// a double-quoted attribute, a `v-model` (or `v-model:<arg>`) on a member or
// on a row of a `v-for` over one (nested rows too), and a handler that writes
// through either. A name that is a member's (a helper's parameter or a
// destructured member) counts as the member. Blind spots: an alias under
// another name (a member or a row held in a local, a slice handle held in a
// local) and a drawing held under another name; a write by a helper the
// member is handed to (`Reflect.set`, `Object.defineProperty`,
// `Array.prototype.push.call`, a function or a component prop that mutates
// it) or by a destructuring assignment; a single-quoted template attribute.
// The type checker sees an alias where the holder is typed.
// The settings are not scanned: their controls write them (plan-r2 §13 item 3).
import { DRAWING_EDITOR_KEYS as MEMBERS } from './drawing-state.mjs';

const MUTATING_METHODS = ['push', 'pop', 'shift', 'unshift', 'splice', 'sort', 'reverse', 'fill', 'copyWithin', 'add', 'delete', 'clear', 'set'];
const ASSIGNMENT = /^\s*(?:(?:[-+*/%&|^]|\*\*|<<|>>>?|&&|\|\||\?\?)?=(?![=>])|\+\+|--)/;
const PREFIX_WRITE = /(?:\bdelete|\+\+|--)\s*$/;
const DRAWING = String.raw`(?:[\w$]*[dD]rawing|activeDrawing\s*\(\s*\)|drawings\s*(?:\??\.\s*[\w$]+|\[[^\]]*\]))`;
const DRAWING_RECEIVER = new RegExp(String.raw`${DRAWING}\s*\??$`);
const DRAWING_COMPUTED_KEY = new RegExp(String.raw`${DRAWING}\s*(?:\?\.)?\s*\[[^\]]*\]`, 'g');
// Comments and quoted strings blanked, so offsets stay and prose is no writer.
const codeOnly = (source) => source.replace(
  /\/\*[\s\S]*?\*\/|\/\/[^\n]*|'(?:\\.|[^'\\\n])*'|"(?:\\.|[^"\\\n])*"/g, (text) => text.replace(/[^\n]/g, ' ')
);
// Whether the code after a name writes through it.
const writesThrough = (after) => {
  let rest = after;
  let last = '';
  for (let step; (step = rest.match(/^\s*(?:\??\.\s*([\w$]+)|\[[^\]]*\])/)); rest = rest.slice(step[0].length)) last = step[1] || '';
  return (MUTATING_METHODS.includes(last) && /^\s*\(/.test(rest)) || ASSIGNMENT.test(rest);
};
const lineOf = (text, index) => text.slice(0, index).split('\n').length;

export const memberWriteSites = (source, members = MEMBERS) => {
  const code = codeOnly(source);
  const sites = [...code.matchAll(new RegExp(`\\??\\.\\s*(${members.join('|')})\\b`, 'g'))].flatMap((match) => {
    const before = code.slice(0, match.index);
    if (!DRAWING_RECEIVER.test(before)) return [];
    const receiver = before.match(/[\w$]+(?:\(\s*\))?(?:\??\.[\w$]+(?:\(\s*\))?|\[[^\]]*\])*\??$/)?.[0] ?? '';
    const start = match.index - receiver.length;
    const writes = writesThrough(code.slice(match.index + match[0].length))
      || PREFIX_WRITE.test(code.slice(0, start))
      || /\bObject\.assign\s*\(\s*$/.test(code.slice(0, start));
    return writes ? [`${match[1]}:${lineOf(code, match.index)}`] : [];
  });
  // A drawing's computed key may name any member: a write through it counts.
  const computedKey = [...code.matchAll(DRAWING_COMPUTED_KEY)].flatMap((match) => {
    const before = code.slice(0, match.index);
    const writes = writesThrough(code.slice(match.index + match[0].length))
      || PREFIX_WRITE.test(before) || /\bObject\.assign\s*\(\s*$/.test(before);
    return writes ? [`[computed]:${lineOf(code, match.index)}`] : [];
  });
  // A name that is a member's (a parameter, a destructured member, an alias
  // under its name): a write through it, `delete`, or `Object.assign` into it.
  // A rebinding of the name is no write of a member.
  const bare = [...code.matchAll(new RegExp(`(?<![\\w$.])(${members.join('|')})\\b`, 'g'))].flatMap((match) => {
    const before = code.slice(0, match.index);
    const after = code.slice(match.index + match[0].length);
    const step = /^\s*(?:\??\.|\[)/.test(after);
    const writes = (step && (writesThrough(after) || PREFIX_WRITE.test(before)))
      || /\bObject\.assign\s*\(\s*$/.test(before);
    return writes ? [`${match[1]}:${lineOf(code, match.index)}`] : [];
  });
  return [...sites, ...computedKey, ...bare].sort((left, right) => Number(left.split(':')[1]) - Number(right.split(':')[1]));
};

// The markup of the element whose start tag holds `index` (a void element: its tag).
const elementRange = (html, index) => {
  const start = html.lastIndexOf('<', index);
  const name = html.slice(start + 1).match(/^[\w-]+/)[0];
  if (/^(?:input|img|br|hr|meta|link)$/i.test(name)) return [start, html.indexOf('>', index) + 1];
  const tag = new RegExp(`<(/?)${name}\\b`, 'g');
  tag.lastIndex = start;
  for (let depth = 0, found; (found = tag.exec(html));) {
    if ((depth += found[1] ? -1 : 1) === 0) return [start, html.indexOf('>', found.index) + 1];
  }
  return [start, html.length];
};
const BINDINGS = /\s(v-model(?::[\w-]+)?(?:\.[\w.]+)?|@[\w.:-]+|v-on:[\w.:-]+)="([^"]*)"/g;

export const templateMemberWrites = (html, members = MEMBERS) => {
  /** Rows of a `v-for` over a member or over another such row, by the element range they name. */
  const rows = [];
  const named = (at) => new Set([...members, ...rows.filter(({ start, end }) => start <= at && at < end).map(({ name }) => name)]);
  for (const match of html.matchAll(/v-for="\(?\s*([\w$]+)[^"]*?\s+(?:in|of)\s+([\w$]+)[^"]*"/g)) {
    if (!named(match.index).has(match[2])) continue;
    const [start, end] = elementRange(html, match.index);
    rows.push({ name: match[1], start, end });
  }
  return [...html.matchAll(BINDINGS)].flatMap((match) => {
    const [text, attribute, expression] = match;
    const names = [...named(match.index)];
    const code = codeOnly(expression);
    const uses = [...code.matchAll(new RegExp(`(?<![\\w$.])(?:${names.join('|')})\\b`, 'g'))];
    const writes = attribute.startsWith('v-model') ? uses.some((use) => use.index === code.search(/\S/))
      : uses.some((use) => writesThrough(code.slice(use.index + use[0].length)));
    return writes ? [`${text.trim()}:${lineOf(html, match.index)}`] : [];
  });
};
