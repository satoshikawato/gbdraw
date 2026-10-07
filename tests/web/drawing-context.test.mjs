// The drawing context (plan OV-80, PR-0a): a diagram mode's settings and edits
// live in its drawing (`state.drawings[mode]`, `state.activeDrawing()`), and a
// service reads them from the drawing it is given, not from `state`. Until the
// per-mode settings change, both modes share one drawing whose members are the
// `state` objects, so this conversion changes no behavior.
//
// The guard asserts:
//   1. every `state` key is in exactly one class (helpers/drawing-state.mjs);
//   2. the drawing holds exactly the drawing keys, as the `state` objects, and
//      `activeDrawing()` returns the drawing of the shown mode;
//   3. a converted module reads no drawing key from `state` (the list grows
//      with each conversion and becomes every module but state.js);
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

// Modules that read every drawing key from a drawing.
const CONVERTED_MODULES = [
  'app/candidate-render.js',
  'app/feature-editor/color-actions.js',
  'app/feature-editor/label-actions.js',
  'app/feature-editor/pattern-drafts.js',
  'app/feature-editor/placement-actions.js',
  'app/feature-editor/rule-actions.js',
  'app/feature-editor/visibility-actions.js',
  'app/feature-search/preview-actions.js',
  'app/legend-layout.js',
  'app/legend-layout/canvas-actions.js',
  'app/legend.js',
  'app/legend/entry-actions.js',
  'app/legend/sort-actions.js',
  'app/legend/stroke-actions.js',
  'app/legend/track-data-styles.js',
  'app/palettes.js',
  'app/results.js',
  'app/rule-matching.js',
  'app/run-analysis.js',
  'app/run-info.js',
  'app/svg-styles.js',
  'mode-profiles.js',
  'services/config.js',
  'services/conservation-series.js',
  'services/current-option-values.js',
  'services/feature-catalog.js',
  'services/feature-visibility.js',
  'services/gallery-session-migration.js',
  'services/gallery-session-publication.js',
  'services/history-snapshot.js',
  'services/layout-preferences.js',
  'services/losat.js',
  'services/reset.js',
  'services/session-active-config-contract.js',
  'services/session-request.js',
  'services/session-resources.js',
  'services/svg-result-ingestion.js',
  'services/svg-serialization.js'
];

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

test('every state key is in exactly one class', () => {
  const classOf = new Map();
  for (const [name, keys] of Object.entries(CLASSES)) {
    for (const key of keys) {
      assert.ok(!classOf.has(key), `${key} is in ${classOf.get(key)} and ${name}`);
      classOf.set(key, name);
    }
  }
  const unclassified = Object.keys(state).filter((key) => !classOf.has(key));
  assert.deepEqual(unclassified, [], 'classify each new state key in tests/web/helpers/drawing-state.mjs');
  assert.deepEqual([...classOf.keys()].filter((key) => !Object.hasOwn(state, key)), [], 'a classified key left state');
});

test('the drawing holds the drawing keys as the state objects, one drawing for both modes', () => {
  const { circular, linear } = state.drawings;
  assert.deepEqual(Object.keys(state.drawings), ['circular', 'linear']);
  assert.equal(circular, linear, 'both modes share one drawing until the per-mode settings change');
  assert.ok(Object.isFrozen(state.drawings) && Object.isFrozen(circular));
  assert.deepEqual(Object.keys(circular), [...DRAWING_KEYS]);
  for (const key of DRAWING_KEYS) assert.equal(circular[key], state[key], key);
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

test('a converted module reads no drawing key from state', () => {
  const found = {};
  for (const module of CONVERTED_MODULES) {
    const reads = stateReads(readFileSync(join(JS, module), 'utf8'));
    if (reads.length) found[module] = [...new Set(reads)].sort();
  }
  assert.deepEqual(found, {}, 'read these keys from the drawing the caller passes');
});

test('the drawing-key scan finds each read form', () => {
  assert.deepEqual(stateReads('state.form.x; state?.adv; const { legendEntries, mode } = state;'),
    ['form', 'adv', 'legendEntries']);
  assert.deepEqual(stateReads('drawing.form; sourceState.adv; // state.form\nstate.mode.value'), []);
});

const modules = (directory) => readdirSync(directory).flatMap((name) => {
  const path = join(directory, name);
  return statSync(path).isDirectory() ? modules(path) : path.endsWith('.js') ? [path] : [];
});

test('a service never resolves the active drawing', () => {
  const resolving = modules(join(JS, 'services'))
    .filter((path) => /\bactiveDrawing\b/.test(readFileSync(path, 'utf8')))
    .map((path) => relative(JS, path));
  assert.deepEqual(resolving, [], 'the owner of the action passes the drawing');
});
