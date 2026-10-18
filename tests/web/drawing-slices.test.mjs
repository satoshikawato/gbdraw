// Typed drawing slices (plan-r2 §3.2, R15): a drawing's editor edits belong to
// slices (gbdraw/web/js/services/drawing/), which alone write them. The guard
// asserts, without a member or slice literal:
//   1. each slice's view holds exactly its members, every drawing shows each
//      slice's view, and no member is in two slices;
//   2. a slice imports only leaves (R13 layer 0) and other slices;
//   3. Load, a History restore, and Reset write a slice through one install:
//      `install(capture(x))` gives back `x`, also after `reset`, on every
//      drawing of each Session 46 fixture.
import assert from 'node:assert/strict';
import { readdirSync, readFileSync } from 'node:fs';
import { dirname, join, relative } from 'node:path';
import test from 'node:test';
import { fileURLToPath } from 'node:url';
import { gunzipSync } from 'node:zlib';
import { File } from 'node:buffer';
import { webLayerOf } from '../../tools/web-owner-graph-detectors.mjs';
import { DRAWING_EDITOR_KEYS, DRAWING_SLICES } from './helpers/drawing-state.mjs';
import { installSessionImportWorker } from './helpers/session-import-node.mjs';
import { validateWebConfigOverrides } from './helpers/python-evaluator.mjs';
import { installFakeSvgDom } from './fake-svg-dom.mjs';

globalThis.window = {
  Vue: {
    ref: (value) => ({ value }), reactive: (value) => value,
    computed: (getter) => ({ get value() { return getter(); } }),
    nextTick: async () => {}
  },
  DOMPurify: { sanitize: (value) => value }
};
globalThis.document = { createElement: () => ({ click: () => {} }) };
installFakeSvgDom();
globalThis.File = File;
globalThis.alert = () => {};
installSessionImportWorker();

const { importSession, setUnmanagedConfigOverrideValidator } = await import('../../gbdraw/web/js/services/config.js');
setUnmanagedConfigOverrideValidator(async (payload) => ({ result: await validateWebConfigOverrides(payload) }));
const { transformLegacyResultSvg } = await import('../../gbdraw/web/js/app/app-setup.js');
const { state } = await import('../../gbdraw/web/js/state.js');

const ROOT = join(dirname(fileURLToPath(import.meta.url)), '..', '..');
const JS = join(ROOT, 'gbdraw', 'web', 'js');
const SLICE_DIRECTORY = join(JS, 'services', 'drawing');
const FIXTURES = join(ROOT, 'tests', 'fixtures', 'sessions');
const MODES = Object.keys(state.drawings);
const SLICES = Object.entries(DRAWING_SLICES);

test('each slice view holds its members, every drawing shows it, and no member is in two slices', async () => {
  const sliceOf = new Map();
  for (const [name, slice] of SLICES) {
    const { handle, view } = slice.create();
    assert.ok(Object.isFrozen(handle) && Object.isFrozen(view), `${name}: the handle and the view are frozen`);
    assert.deepEqual(Object.keys(view), [...slice.members], `${name}: the view holds the slice's members`);
    for (const member of slice.members) {
      assert.ok(!sliceOf.has(member), `${member} is in ${sliceOf.get(member)} and ${name}`);
      sliceOf.set(member, name);
      assert.ok(DRAWING_EDITOR_KEYS.includes(member), `${name}.${member} is an editor member`);
    }
  }
  const modules = readdirSync(SLICE_DIRECTORY).filter((file) => file.endsWith('.js')).sort();
  const exported = await Promise.all(modules.map((file) => import(join(SLICE_DIRECTORY, file))));
  const descriptors = new Set(exported.flatMap((module) => Object.values(module)
    .filter((value) => value && typeof value === 'object' && Array.isArray(value.members))));
  assert.deepEqual(new Set(SLICES.map(([, slice]) => slice)), descriptors, 'DRAWING_SLICES names each slice module once');
  for (const mode of MODES) {
    const drawing = state.drawings[mode];
    assert.deepEqual(Object.keys(drawing.slices), SLICES.map(([name]) => name), `${mode}: one handle per slice`);
    for (const [name, slice] of SLICES) {
      const captured = slice.capture(drawing.slices[name]);
      for (const member of slice.members) {
        assert.ok(Object.hasOwn(captured, member), `${mode} ${name}: capture holds ${member}`);
      }
    }
  }
});

// Static imports, with comments masked (a JSDoc `@import` is not an edge).
const importsOf = (source) => [...source.replace(/\/\*[\s\S]*?\*\//g, '').replace(/(^|[^:'"`\\])\/\/.*$/gm, '$1')
  .matchAll(/\bfrom\s+['"]([^'"]+)['"]/g)].map((match) => match[1]);

test('a slice imports only leaves and other slices', () => {
  const found = [];
  for (const file of readdirSync(SLICE_DIRECTORY).filter((name) => name.endsWith('.js'))) {
    const path = join(SLICE_DIRECTORY, file);
    for (const target of importsOf(readFileSync(path, 'utf8'))) {
      const module = relative(JS, join(dirname(path), target)).split('\\').join('/');
      if (webLayerOf(module) !== 0 && !module.startsWith('services/drawing/')) found.push(`${file} imports ${module}`);
    }
  }
  assert.deepEqual(found, []);
});

const fixtureNames = readdirSync(FIXTURES).filter((name) => /\.v46\.gbdraw-session\.json\.gz$/.test(name)).sort();
const load = (name) => importSession({
  target: { files: [new Blob([gunzipSync(readFileSync(join(FIXTURES, name)))], { type: 'application/json' })], value: 'selected' }
}, { transformLegacyResultSvg, confirmLoadRepairs: async () => true });

test('install(capture(x)) gives back x, also after reset, on the Session 46 fixtures', async () => {
  assert.ok(fixtureNames.length > 0, 'the Session 46 fixtures exist');
  for (const name of fixtureNames) {
    const loaded = await load(name);
    assert.equal(loaded?.status, 'ok', `${name}: ${loaded?.error?.message ?? loaded?.status}`);
    for (const mode of MODES) {
      for (const [sliceName, slice] of SLICES) {
        const handle = state.drawings[mode].slices[sliceName];
        const captured = slice.capture(handle);
        const json = JSON.stringify(captured);
        slice.install(handle, captured);
        assert.equal(JSON.stringify(slice.capture(handle)), json, `${name} ${mode} ${sliceName}: install(capture(x))`);
        slice.reset(handle);
        slice.install(handle, JSON.parse(json));
        assert.equal(JSON.stringify(slice.capture(handle)), json, `${name} ${mode} ${sliceName}: after reset`);
      }
    }
  }
});
