// The drawing context (plan OV-80, PR-0a to PR-0c): a diagram mode's settings
// and edits live only in its drawing (`state.drawings[mode]`,
// `state.activeDrawing()`). A service reads them from the drawing it is given;
// an owner resolves the drawing of each action. Until the per-mode settings
// change, both modes share one drawing, so the conversion changes no behavior.
//
// The guard asserts:
//   1. every `state` key and every drawing member is in exactly one class
//      (helpers/drawing-state.mjs);
//   2. the drawing holds exactly the drawing keys, `state` holds none of them,
//      and `activeDrawing()` returns the drawing of the shown mode;
//   3. no module reads a drawing key from `state`;
//   4. a service never resolves the active drawing.
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
  DRAWING_KEYS,
  PROJECT_KEYS,
  STORE_KEYS,
  TRANSIENT_KEYS,
  UI_KEYS
} from './helpers/drawing-state.mjs';

globalThis.window = { Vue: {
  ref: (value) => ({ value }), reactive: (value) => value,
  computed: (getter) => ({ get value() { return getter(); } }), nextTick: async () => {}
} };
const { state } = await import('../../gbdraw/web/js/state.js');

const JS = join(dirname(fileURLToPath(import.meta.url)), '..', '..', 'gbdraw', 'web', 'js');

const CLASSES = {
  DRAWING_DRAFT_KEYS,
  DRAWING_EDITOR_KEYS,
  DRAWING_DERIVED_KEYS,
  ARTIFACT_KEYS,
  PROJECT_KEYS,
  UI_KEYS,
  TRANSIENT_KEYS,
  STORE_KEYS
};

test('every state key and drawing member is in exactly one class', () => {
  const classOf = new Map();
  for (const [name, keys] of Object.entries(CLASSES)) {
    for (const key of keys) {
      assert.ok(!classOf.has(key), `${key} is in ${classOf.get(key)} and ${name}`);
      classOf.set(key, name);
    }
  }
  const members = [...Object.keys(state), ...Object.keys(state.drawings.circular)];
  const unclassified = members.filter((key) => !classOf.has(key));
  assert.deepEqual(unclassified, [], 'classify each new state key or drawing member in tests/web/helpers/drawing-state.mjs');
  assert.deepEqual([...classOf.keys()].filter((key) => !members.includes(key)), [], 'a classified key left state and the drawing');
});

test('only the drawing holds the drawing keys, one drawing for both modes', () => {
  const { circular, linear } = state.drawings;
  assert.deepEqual(Object.keys(state.drawings), ['circular', 'linear']);
  assert.equal(circular, linear, 'both modes share one drawing until the per-mode settings change');
  assert.ok(Object.isFrozen(state.drawings) && Object.isFrozen(circular));
  assert.deepEqual(Object.keys(circular), [...DRAWING_KEYS]);
  assert.deepEqual(DRAWING_KEYS.filter((key) => Object.hasOwn(state, key)), [], 'a drawing key is on state');
  for (const mode of ['circular', 'linear']) {
    state.mode.value = mode;
    assert.equal(state.activeDrawing(), state.drawings[mode]);
  }
  state.mode.value = 'circular';
});

// `state.<key>`, `state?.<key>`, and `{ <key> } = state`, with comments removed.
const stateReads = (source) => {
  const code = source.replace(/\/\*[\s\S]*?\*\//g, '').replace(/(^|[^:'"`\\])\/\/.*$/gm, '$1');
  const keys = DRAWING_KEYS.join('|');
  const reads = [...code.matchAll(new RegExp(`\\bstate\\s*\\??\\.\\s*(${keys})\\b`, 'g'))].map((match) => match[1]);
  for (const [, names] of code.matchAll(/\{([^{}]*)\}\s*=\s*state\b/g)) {
    reads.push(...names.split(',').map((name) => name.trim().split(/[\s:=]/)[0]).filter((name) => DRAWING_KEYS.includes(name)));
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

test('the drawing-key scan finds each read form', () => {
  assert.deepEqual(stateReads('state.form.x; state?.adv; const { legendEntries, mode } = state;'),
    ['form', 'adv', 'legendEntries']);
  assert.deepEqual(stateReads('drawing.form; sourceState.adv; // state.form\nstate.mode.value'), []);
});

test('a service never resolves the active drawing', () => {
  const resolving = modules(join(JS, 'services'))
    .filter((path) => /\bactiveDrawing\b/.test(readFileSync(path, 'utf8')))
    .map((path) => relative(JS, path));
  assert.deepEqual(resolving, [], 'the owner of the action passes the drawing');
});
