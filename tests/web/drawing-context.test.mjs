// The drawing context (plan OV-80, PR-0a to PR-1): a diagram mode's settings
// and edits live only in its drawing (`modeDrawing(state, mode)`,
// `state.activeDrawing()`). A service reads them from the drawing it is given;
// an owner resolves the drawing of each action. Each mode has its own drawing
// (PR-1).
//
// The guard asserts:
//   1. every `state` key and every drawing member is in exactly one class
//      (helpers/drawing-state.mjs);
//   2. each mode's drawing holds exactly its identity and the drawing keys,
//      the two drawings share no member object, `state` holds no drawing
//      key, and `activeDrawing()` returns the shown drawing
//      (tests/web/drawing-slices.test.mjs checks its one writer);
//   3. no module reads a drawing key from `state`;
//   4. a service never resolves the active drawing, and Generate resolves it
//      once, at its entry.
import assert from 'node:assert/strict';
import { readdirSync, readFileSync, statSync } from 'node:fs';
import { dirname, join, relative } from 'node:path';
import test from 'node:test';
import { fileURLToPath } from 'node:url';
import {
  ARTIFACT_KEYS,
  DRAWING_DERIVED_KEYS,
  DRAWING_DRAFT_KEYS,
  DRAWING_EDITOR_KEYS,
  DRAWING_HANDLE_KEYS,
  DRAWING_IDENTITY_KEYS,
  DRAWING_KEYS,
  PROJECT_KEYS,
  STORE_KEYS,
  TRANSIENT_KEYS,
  UI_KEYS, modeDrawing } from './helpers/drawing-state.mjs';

globalThis.window = { Vue: {
  ref: (value) => ({ value }), reactive: (value) => value,
  computed: (getter) => ({ get value() { return getter(); } }), nextTick: async () => {}
} };
const { state } = await import('../../gbdraw/web/js/state.js');
await (await import('./helpers/drawing-state.mjs')).addStartLinearDrawing();

const JS = join(dirname(fileURLToPath(import.meta.url)), '..', '..', 'gbdraw', 'web', 'js');

// The classes of the `state` keys and of the drawing members. A drawing's
// identity names a state key too (`mode`: the shown drawing's), so the two
// sets are classed apart; no drawing key is a state key (below).
const STATE_CLASSES = { ARTIFACT_KEYS, PROJECT_KEYS, UI_KEYS, TRANSIENT_KEYS, STORE_KEYS };
const DRAWING_CLASSES = {
  DRAWING_IDENTITY_KEYS, DRAWING_DRAFT_KEYS, DRAWING_EDITOR_KEYS, DRAWING_DERIVED_KEYS, DRAWING_HANDLE_KEYS
};

test('every state key and drawing member is in exactly one class', () => {
  for (const [classes, members] of [
    [STATE_CLASSES, Object.keys(state)],
    [DRAWING_CLASSES, Object.keys(modeDrawing(state, 'circular'))]
  ]) {
    const classOf = new Map();
    for (const [name, keys] of Object.entries(classes)) {
      for (const key of keys) {
        assert.ok(!classOf.has(key), `${key} is in ${classOf.get(key)} and ${name}`);
        classOf.set(key, name);
      }
    }
    const unclassified = members.filter((key) => !classOf.has(key));
    assert.deepEqual(unclassified, [], 'classify each new state key or drawing member in tests/web/helpers/drawing-state.mjs');
    assert.deepEqual([...classOf.keys()].filter((key) => !members.includes(key)), [], 'a classified key left state or the drawing');
  }
});

// R13: `state` is frozen, so no module adds or replaces a member at run time;
// tsc (R14, tests/web/typed-boundaries.test.mjs) refuses such a write in a
// module, since the frozen object's type is read-only.
test('state is frozen: no module adds or replaces a state member', () => {
  assert.ok(Object.isFrozen(state));
});

test('only the drawings hold the drawing keys, one drawing per mode', () => {
  const circular = modeDrawing(state, 'circular');
  const linear = modeDrawing(state, 'linear');
  assert.deepEqual(state.drawingList().map((drawing) => drawing.id), ['circular', 'linear']);
  assert.notEqual(circular, linear, 'each diagram mode has its own drawing');
  assert.ok(Object.isFrozen(circular) && Object.isFrozen(linear));
  assert.deepEqual(Object.keys(circular), [...DRAWING_IDENTITY_KEYS, ...DRAWING_KEYS]);
  assert.deepEqual(Object.keys(linear), [...DRAWING_IDENTITY_KEYS, ...DRAWING_KEYS]);
  assert.deepEqual([circular.id, circular.mode, linear.id, linear.mode], ['circular', 'circular', 'linear', 'linear']);
  // A member object shared by both drawings would carry an edit across modes.
  const shared = DRAWING_KEYS.filter((key) => typeof circular[key] === 'object' && circular[key] !== null
    && circular[key] === linear[key]);
  assert.deepEqual(shared, [], 'each drawing owns its members');
  assert.deepEqual(DRAWING_KEYS.filter((key) => Object.hasOwn(state, key)), [], 'a drawing key is on state');
  assert.equal(state.activeDrawing(), modeDrawing(state, state.activeDrawingId.value));
});

// `state.<key>`, `state?.<key>`, and `{ <key> } = state`, with comments removed.
const stateReads = (source, readKeys = DRAWING_KEYS) => {
  const code = source.replace(/\/\*[\s\S]*?\*\//g, '').replace(/(^|[^:'"`\\])\/\/.*$/gm, '$1');
  const keys = readKeys.join('|');
  const reads = [...code.matchAll(new RegExp(`\\bstate\\s*\\??\\.\\s*(${keys})\\b`, 'g'))].map((match) => match[1]);
  for (const [, names] of code.matchAll(/\{([^{}]*)\}\s*=\s*state\b/g)) {
    reads.push(...names.split(',').map((name) => name.trim().split(/[\s:=]/)[0]).filter((name) => readKeys.includes(name)));
  }
  return reads;
};

const modules = (directory) => readdirSync(directory).flatMap((name) => {
  const path = join(directory, name);
  return statSync(path).isDirectory() ? modules(path) : path.endsWith('.js') ? [path] : [];
});

test('no module reads a drawing key from state', () => {
  const found = {};
  for (const path of modules(JS)) {
    const module = relative(JS, path);
    if (module === 'state.js') continue;
    const reads = stateReads(readFileSync(path, 'utf8'));
    if (reads.length) found[module] = [...new Set(reads)].sort();
  }
  assert.deepEqual(found, {}, 'read these keys from the drawing of the action');
});

// The project input keys the inputs slice replaced (multi-drawing B, DESIGN
// §2): each is gone from `state`, and no module reads it from `state`.
// The five card files are each drawing's too (`drawing.files`, R16-4), so
// the project's `cardFiles` is gone.
const MOVED_INPUT_KEYS = Object.freeze(['linearSeqs', 'circularRecordList', 'circularRecordDiscovery', 'cInputType', 'lInputType', 'cardFiles']);
test('no module reads a project input key the drawing\'s inputs replaced', () => {
  assert.deepEqual(MOVED_INPUT_KEYS.filter((key) => Object.hasOwn(state, key)), [], 'a moved input key is on state');
  const found = {};
  for (const path of modules(JS)) {
    const source = readFileSync(path, 'utf8');
    const reads = stateReads(source, MOVED_INPUT_KEYS);
    if (reads.length) found[relative(JS, path)] = [...new Set(reads)].sort();
  }
  assert.deepEqual(found, {}, 'read the inputs from the drawing of the action (drawing.files, drawing.linearSeqs, ...)');
});

test('the drawing-key scan finds each read form', () => {
  assert.deepEqual(stateReads('state.form.x; state?.adv; const { legendEntries, mode } = state;'),
    ['form', 'adv', 'legendEntries']);
  assert.deepEqual(stateReads('drawing.form; sourceState.adv; // state.form\nstate.mode.value'), []);
});

// `activeDrawing`, or a lookup of a drawing by the shown ID
// (`drawingOf(state.activeDrawingId.value)`, `drawings.get(...)`), OV-1218.
const resolvesShownDrawing = (source) => /\bactiveDrawing\b/.test(source)
  || /\b(?:drawingOf\s*\(|drawings\s*(?:\.\s*get\s*\(|\[))[^()\]]*\bactiveDrawingId\b/.test(source);

test('a service never resolves the active drawing', () => {
  const resolving = modules(join(JS, 'services'))
    .filter((path) => resolvesShownDrawing(readFileSync(path, 'utf8')))
    .map((path) => relative(JS, path));
  assert.deepEqual(resolving, [], 'the owner of the action passes the drawing');
});

test('the active-drawing scan finds each resolving form', () => {
  for (const form of ['state.activeDrawing()', 'state.drawingOf(state.activeDrawingId.value)',
    'drawingOf( activeDrawingId.value )', 'state.drawings.get(state.activeDrawingId.value)',
    'drawings[state.activeDrawingId.value]']) {
    assert.ok(resolvesShownDrawing(form), form);
  }
  for (const form of ['drawing.id === state.activeDrawingId.value', 'shownDrawingId: state.activeDrawingId.value',
    'state.drawingOf(handle.drawingId)', 'managerOf(activeDrawingId())']) {
    assert.ok(!resolvesShownDrawing(form), form);
  }
});

test('Generate resolves its drawing once, at its entry', () => {
  // The run, its History checkpoint, and its decoration capture use the
  // drawing the click resolved; the other run-analysis entries draw the
  // committed request's mode or a fixed mode (ov80 pr0b picks).
  const source = readFileSync(join(JS, 'app', 'run-analysis.js'), 'utf8');
  const calls = [...source.matchAll(/\bactiveDrawing\s*\(/g)];
  assert.equal(calls.length, 1, 'run-analysis.js resolves the shown drawing only in runAnalysis');
  assert.match(source, /const runAnalysis = async \([^)]*\) => \{\s*const drawing = state\.activeDrawing\(\);/);
});
