// Shared-value scope (R15): a value of one drawing lives as long as that
// drawing and travels with it; a value every drawing shares has a named
// reason to be shared. A per-drawing value kept where every drawing reads or
// writes it, or kept by a drawing's ID past the drawing's life, reaches another
// drawing: a second drawing reads or overwrites it, or a new drawing that takes
// a deleted drawing's ID inherits it.
//
// The guard asserts:
//   1. every generated-artifact `state` key is an artifact slot key (the slot
//      keeps it while another drawing is shown) or is registered below with the
//      reason it needs no slot;
//   2. every key of the shared classes (`PROJECT_KEYS`, `UI_KEYS` of
//      helpers/drawing-state.mjs, registered lists) names on its own line why
//      every drawing shares it;
//   3. no long-lived container (declared at module level or at the top level of
//      a `create*`/`setup*` factory) is keyed by a drawing's ID, except the
//      registered stores that own a drawing's lifetime. Key such a value by the
//      drawing object (`WeakMap<DrawingState, ...>`).
//
// Check 3 reads the source text, so it does not see: a container held as a
// property of a long-lived object (`cache.byId.set(drawing.id, v)`); a key
// through a local (`const id = drawing.id; m.set(id, v)`) or a drawing in a
// variable not named `drawing`, `arriving` or `departing`; `Set.add(drawing.id)`
// with no JSDoc type and no `.has`; `Object.create(null)` written with
// `[d.id]`; a factory that is not `create*`/`setup*` or is indented by other
// than two spaces; a typedef imported from another module; a composite key
// (`${drawing.id}:${uid}`). None of these holds a value per drawing on dev
// when this guard was added. A single value that names one drawing by its ID
// (not a container) is outside the check (R15 says how it is cleared).
import assert from 'node:assert/strict';
import { readdirSync, readFileSync, statSync } from 'node:fs';
import { dirname, join, relative } from 'node:path';
import test from 'node:test';
import { fileURLToPath } from 'node:url';
import { maskJavaScript } from '../../tools/web-change-source.mjs';
import { ARTIFACT_KEYS } from './helpers/drawing-state.mjs';
import { ARTIFACT_SLOT_KEYS } from '../../gbdraw/web/js/services/artifact-slot.js';

const ROOT = join(dirname(fileURLToPath(import.meta.url)), '..', '..');
const JS = join(ROOT, 'gbdraw', 'web', 'js');
const HELPER = join(ROOT, 'tests', 'web', 'helpers', 'drawing-state.mjs');

// 1. The generated-artifact keys that no slot keeps, each with the reason its
// value cannot reach another drawing.
// A registered design-rule allowlist (tools/web-design-rule-guards.json, R15,
// kind `set`): a new entry is an authority-only change.
export const ARTIFACT_KEYS_OUTSIDE_SLOT = Object.freeze([
  // The slot keeps both in `slot.legend`; every slot install writes them (`installShownSetLegend`).
  'originalLegendOrder', 'originalLegendColors',
  // Read again from the mounted Result's composition at every mount (`syncStateFromComposition`,
  // `setupDiagramDrag`); History records them per Result (`compositionUserDeltas`). Without a
  // Result nothing shows them (Reset Layout needs a shown figure).
  'legendCurrentOffset', 'diagramOffset', 'plotTitleUserOffset',
  // Every mount of a Result sets it: zero for a Generate, a Result selection or a drawing switch
  // (`syncLengthBarElement`), the restored step's own value for a History or Session restore. A drawing
  // that shows no Result may hold the departing drawing's value (also in its History steps), which no
  // reader shows: the bar and Reset Layout need a shown Result (OV-3893, not a product defect).
  'lengthBarUserOffset',
  // The slot's own `drawingId`: every slot install and owner-set install writes it from there.
  'generatedDrawingId'
]);

/**
 * The generated-artifact keys outside the slot that `allowed` does not
 * register, and the registered keys that are no longer artifact keys outside
 * the slot.
 * @param {readonly string[]} artifactKeys @param {readonly string[]} slotKeys @param {readonly string[]} allowed
 */
const artifactKeysOutsideSlot = (artifactKeys, slotKeys, allowed) => ({
  unregistered: artifactKeys.filter((key) => !slotKeys.includes(key) && !allowed.includes(key)),
  stale: allowed.filter((key) => slotKeys.includes(key) || !artifactKeys.includes(key))
});

test('a generated artifact travels with its drawing', () => {
  const { unregistered, stale } = artifactKeysOutsideSlot(ARTIFACT_KEYS, ARTIFACT_SLOT_KEYS, ARTIFACT_KEYS_OUTSIDE_SLOT);
  assert.deepEqual(unregistered, [], 'keep each in the artifact slot (ARTIFACT_SLOT_EMPTY, services/artifact-slot.js)');
  assert.deepEqual(stale, [], 'remove these from ARTIFACT_KEYS_OUTSIDE_SLOT');
});

test('the artifact key check finds a key outside the slot and a stale entry', () => {
  assert.deepEqual(artifactKeysOutsideSlot(['results', 'plan', 'offset'], ['results'], ['offset', 'results', 'gone']),
    { unregistered: ['plan'], stale: ['results', 'gone'] });
});

// 2. The shared classes: each key of a registered class literal sits on its
// own line with a `//` comment that says why every drawing shares it.
const SHARED_CLASSES = Object.freeze(['PROJECT_KEYS', 'UI_KEYS']);

/**
 * The entries of `symbol`'s array literal that do not stand alone on a line
 * with a reason after them, and every entry it holds.
 * @param {string} source @param {string} symbol
 */
const unexplainedEntries = (source, symbol) => {
  const literal = source.match(new RegExp(`export const ${symbol} = Object\\.freeze\\(\\[\\n([\\s\\S]*?)\\n\\]\\);`));
  assert.ok(literal, `${symbol} is a frozen array literal, one entry per line`);
  const entries = [];
  const unexplained = [];
  for (const line of literal[1].split('\n')) {
    const keys = [...line.replace(/\/\/.*$/, '').matchAll(/'([^']+)'/g)].map((match) => match[1]);
    entries.push(...keys);
    if (keys.length === 0) continue;
    if (keys.length > 1 || !/^\s*'[^']+',?\s*\/\/\s*\S/.test(line)) unexplained.push(...keys);
  }
  return { entries, unexplained };
};

test('a value every drawing shares names why', async () => {
  const source = readFileSync(HELPER, 'utf8');
  const helper = await import('./helpers/drawing-state.mjs');
  for (const symbol of SHARED_CLASSES) {
    const { entries, unexplained } = unexplainedEntries(source, symbol);
    assert.deepEqual(unexplained, [], `${symbol}: one key per line, followed by // the reason every drawing shares it`);
    assert.deepEqual(entries, [...helper[symbol]], `${symbol}: the scan reads the exported list`);
  }
});

test('the shared-key scan finds a key without a reason', () => {
  const source = "export const LIST = Object.freeze([\n  'a', // shared: why\n  'b',\n  'c', 'd', // two\n  // a comment line\n  'e' // last\n]);";
  assert.deepEqual(unexplainedEntries(source, 'LIST'), { entries: ['a', 'b', 'c', 'd', 'e'], unexplained: ['b', 'c', 'd'] });
});

// 3. Long-lived containers keyed by a drawing's ID. The stores that own a
// drawing's lifetime, by `module#name`, and how Delete, its Undo, and a new
// drawing that takes a freed ID reach each.
// A registered design-rule allowlist (tools/web-design-rule-guards.json, R15,
// kind `set`): a new entry is an authority-only change.
export const DRAWING_LIFETIME_STORES = Object.freeze([
  // The drawing store: `removeDrawing` drops the drawing; the Undo's `addDrawing` and a new drawing put a drawing in.
  'state.js#drawings',
  // The artifact slots: `removeDrawing` takes the slot out (`takeSlot`) into the step; its Undo puts it back
  // (`putSlot`); a new drawing starts with no slot, or with the one its own operation made (Duplicate's copy, the
  // slot a Redo holds).
  'services/artifact-slot.js#slots',
  // The History router's Histories: `removeMember` hands the History to the project step, its Undo gives it back;
  // `addMember` gives a new drawing an empty History and refuses an ID that has one.
  'services/history.js#managers'
]);

// A JSDoc type that a drawing's ID keys: `Map<DrawingId, ...>`, `ReadonlyMap<DrawingId, ...>`,
// `Set<DrawingId>`, `ReadonlySet<DrawingId>`, `Record<DrawingId, ...>`.
const DRAWING_ID_KEYED = /\b(?:Readonly)?(?:Map|Set)<\s*DrawingId\b|\bRecord<\s*DrawingId\b/;
// An expression that is a drawing's ID: `drawing.id`, `arriving.id`, `state.activeDrawing().id`,
// `drawingId`, `shownDrawingId`, `entry.drawingId`, `state.activeDrawingId.value`.
const DRAWING_ID = String.raw`[\w$.]*(?:[dD]rawingId(?:\.value)?|[dD]rawing(?:\(\))?\.id|(?:arriving|departing)\.id)`;

/** The `{...}` type that starts at `open` in `text`, braces balanced. @param {string} text @param {number} open */
const bracedType = (text, open) => {
  let depth = 0;
  for (let index = open; index < text.length; index += 1) {
    if (text[index] === '{') depth += 1;
    else if (text[index] === '}' && --depth === 0) return text.slice(open + 1, index);
  }
  return '';
};
/** Every type that follows `tag` in the JSDoc of `text`, with its name. @param {string} text @param {string} tag */
const taggedTypes = (text, tag) => [...text.matchAll(new RegExp(`@${tag}\\s*\\{`, 'g'))].map((match) => {
  const open = match.index + match[0].length - 1;
  const type = bracedType(text, open);
  const name = text.slice(open + type.length + 2).match(/^\s*\[?([\w$]+)/)?.[1] ?? '';
  return { type: type.replace(/\s*\n\s*\*\s*/g, ' '), name };
});

/**
 * The containers of a module that live as long as the module or its factory
 * (declared at module level, or at the top level of a `create*`/`setup*`
 * factory body) and that a drawing's ID keys, by JSDoc type (also through a
 * typedef of the module, or a factory option whose type returns one) or by use
 * (`.get`, `.set`, `.has`, `.delete`, or `[...]` with a drawing's ID).
 * @param {string} source
 * @returns {string[]}
 */
const drawingIdKeyedContainers = (source) => {
  const code = maskJavaScript(source);
  const lines = source.split('\n');
  const codeLines = code.split('\n');
  const keyedTypedefs = new Set();
  for (let size = -1; size !== keyedTypedefs.size;) {
    size = keyedTypedefs.size;
    for (const { type, name } of taggedTypes(source, 'typedef')) {
      if (DRAWING_ID_KEYED.test(type) || [...keyedTypedefs].some((keyed) => new RegExp(`\\b${keyed}\\b`).test(type))) {
        keyedTypedefs.add(name);
      }
    }
  }
  const keyed = (type) => DRAWING_ID_KEYED.test(type)
    || [...keyedTypedefs].some((name) => new RegExp(`\\b${name}\\b`).test(type));
  // An option or parameter whose type is a function that returns a keyed container (`reactive`).
  const keyedFactories = new Set([...taggedTypes(source, 'property'), ...taggedTypes(source, 'param')]
    .filter(({ type }) => type.includes('=>') && keyed(type.slice(type.lastIndexOf('=>'))))
    .map(({ name }) => name));
  // The names used with a drawing's ID as the key.
  const usedWithId = new Set([...code.matchAll(new RegExp(
    String.raw`([\w$]+)\s*(?:\.\s*(?:get|set|has|delete)\s*\(\s*${DRAWING_ID}\s*[,)]|\[\s*${DRAWING_ID}\s*\])`, 'g'))]
    .map((match) => match[1]));
  const found = new Set();
  let inFactory = false;
  codeLines.forEach((line, index) => {
    if (/^(?:export\s+)?(?:(?:const|let)\s+(?:create|setup)\w*\s*=|(?:async\s+)?function\s+(?:create|setup)\w*)/.test(line)) {
      inFactory = true;
    } else if (/^[}\])][^{]*$/.test(line)) inFactory = false;
    const declaration = line.match(inFactory ? /^ {2}(?:const|let|var)\s+([\w$]+)\s*=\s*(.*)$/
      : /^(?:export\s+)?(?:const|let|var)\s+([\w$]+)\s*=\s*(.*)$/);
    if (!declaration) return;
    const [, name, initializer] = declaration;
    let start = index;
    while (start > 0 && /^\s*(?:\/\*\*|\*|\/\/)/.test(lines[start - 1])) start -= 1;
    const jsdoc = [...lines.slice(start, index), lines[index]].join('\n');
    const typed = taggedTypes(jsdoc, 'type').some(({ type }) => keyed(type));
    const made = [...keyedFactories].some((factory) => new RegExp(`^${factory}\\s*\\(`).test(initializer));
    if (typed || made || usedWithId.has(name)) found.add(name);
  });
  return [...found];
};

const modules = (directory) => readdirSync(directory).flatMap((name) => {
  const path = join(directory, name);
  return statSync(path).isDirectory() ? modules(path) : path.endsWith('.js') ? [path] : [];
});

test('no long-lived container is keyed by a drawing\'s ID', () => {
  const found = modules(JS).flatMap((path) => drawingIdKeyedContainers(readFileSync(path, 'utf8'))
    .map((name) => `${relative(JS, path).split('\\').join('/')}#${name}`));
  assert.deepEqual(found.filter((entry) => !DRAWING_LIFETIME_STORES.includes(entry)), [],
    'key these by the drawing object (WeakMap<DrawingState, ...>) so they live exactly as long as the drawing');
  assert.deepEqual(DRAWING_LIFETIME_STORES.filter((entry) => !found.includes(entry)), [],
    'remove these from DRAWING_LIFETIME_STORES: they are gone or no longer keyed by a drawing\'s ID');
});

test('the drawing-ID container scan finds each form', () => {
  const scan = (text) => drawingIdKeyedContainers(text).sort();
  // By JSDoc type, at module level and at the top level of a factory.
  assert.deepEqual(scan([
    '/** @type {Map<DrawingId, Options>} */', 'const byId = new Map();',
    '// Each drawing\'s set.', '/** @type {ReadonlyMap<DrawingId, number>} */', 'export let counts = new Map();',
    'export const createThing = ({ reactive }) => {',
    '  /** @type {Set<DrawingId>} */', '  const seen = new Set();',
    '  const record = /** @type {Record<DrawingId, string>} */ ({});',
    '  return {};', '};'
  ].join('\n')), ['byId', 'counts', 'record', 'seen']);
  // Through a typedef, and through an option whose type returns a keyed container.
  assert.deepEqual(scan([
    '/** @typedef {Map<DrawingId, Rows>} RowsByDrawing */', '/** @typedef {{ rows: RowsByDrawing }} Holder */',
    '/**', ' * @typedef {object} Options', ' * @property {(initial: Map<DrawingId, Rows>) =>', ' *   Map<DrawingId, Rows>} reactive', ' */',
    '/** @type {Holder} */', 'const holder = { rows: new Map() };',
    'export function setupRows({ reactive }) {', '  const rows = reactive(new Map());', '  return rows;', '}'
  ].join('\n')), ['holder', 'rows']);
  // By use with a drawing's ID.
  assert.deepEqual(scan([
    'export const createSearch = () => {',
    '  const a = new Map();', '  const b = new Map();', '  const c = new Set();', '  const d = {};', '  const e = new Map();',
    '  const f = new Map();',
    '  a.get(drawing.id); b.set(drawingId, 1); c.has(entry.drawingId); d[drawingId] = 1;',
    '  e.delete(state.activeDrawing().id); f.get(state.activeDrawingId.value);',
    '  return {};', '};'
  ].join('\n')), ['a', 'b', 'c', 'd', 'e', 'f']);
  // Not long-lived (local to a call), keyed by the drawing object, keyed by something else, or named in a comment.
  assert.deepEqual(scan([
    'const load = (drawings) => {', '  /** @type {Map<DrawingId, Draft>} */', '  const drafts = new Map();', '  return drafts;', '};',
    'export const createSearch = () => {',
    '  const run = () => {', '    /** @type {Map<DrawingId, number>} */', '    const local = new Map();', '    local.get(drawing.id);', '  };',
    '  /** @type {WeakMap<DrawingState, Search>} */', '  const byDrawing = new WeakMap();', '  byDrawing.get(drawing);',
    '  const byUid = new Map();', '  byUid.get(row.uid);',
    '  const named = new Map(); // named.get(drawing.id)', '  const text = \'named.get(drawingId)\';',
    '  return { run };', '};'
  ].join('\n')), []);
});
