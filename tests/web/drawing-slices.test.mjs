// Typed drawing slices (plan-r2 §3.2, R15): a drawing's editor edits belong to
// slices (gbdraw/web/js/services/drawing/), which alone write them. The guard
// asserts, without a member or slice literal:
//   1. each slice's view holds exactly its members, every drawing shows each
//      slice's view, and no member is in two slices;
//   2. a slice imports only leaves (R13 layer 0, or a module that imports
//      nothing) and other slices;
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

const { exportSession, importSession, setUnmanagedConfigOverrideValidator } = await import('../../gbdraw/web/js/services/config.js');
const { normalizeAnnotationSets } = await import('../../gbdraw/web/js/services/annotation-state.js');
const { normalizeFeatureVisibilityRule } = await import('../../gbdraw/web/js/services/feature-visibility.js');
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
      const targetPath = join(dirname(path), target);
      const module = relative(JS, targetPath).split('\\').join('/');
      const leaf = webLayerOf(module) === 0 || importsOf(readFileSync(targetPath, 'utf8')).length === 0;
      if (!leaf && !module.startsWith('services/drawing/')) found.push(`${file} imports ${module}`);
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

// 4. Only a slice writes its members (R15; D-15-44: the rows of a list slice
// are the slice's too). The type checker sees a write through a typed drawing,
// whose views are read-only. This scan covers the untyped holders: in each
// module outside services/drawing/, a write through `<drawing>.<member>` (an
// assignment to it, to its `.value`, an element, or a field, a mutating
// method, `delete`, or `Object.assign` into it) whose receiver is a drawing by
// name (an identifier that ends in `drawing` or `Drawing`, `activeDrawing()`,
// `drawings.<mode>`, or `drawings[<mode>]`); and in the template, a `v-model`
// on a member or on a row of a `v-for` over one (nested rows too), and a
// handler that writes through either. Blind spots: an alias (a member or a row
// held in a local or destructured, a slice handle held in a local) and a
// drawing held under another name; the type checker sees those where the
// holder is typed.
const MEMBERS = SLICES.flatMap(([, slice]) => slice.members);
const MUTATING_METHODS = ['push', 'pop', 'shift', 'unshift', 'splice', 'sort', 'reverse', 'fill', 'copyWithin', 'add', 'delete', 'clear', 'set'];
const ASSIGNMENT = /^\s*(?:[-+*/%&|^]|\*\*|<<|>>>?|&&|\|\||\?\?)?=(?![=>])/;
const DRAWING_RECEIVER = /(?:[\w$]*[dD]rawing|activeDrawing\s*\(\s*\)|drawings\s*(?:\??\.\s*[\w$]+|\[[^\]]*\]))\s*\??$/;
// The template writes M3.13 replaces by writer ports; this list only shrinks.
const TEMPLATE_WRITES_BEFORE_M3_13 = Object.freeze([
  'v-model.number="canvasPadding.top"', 'v-model.number="canvasPadding.left"',
  'v-model.number="canvasPadding.right"', 'v-model.number="canvasPadding.bottom"'
]);

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
  return [...code.matchAll(new RegExp(`\\??\\.\\s*(${members.join('|')})\\b`, 'g'))].flatMap((match) => {
    const before = code.slice(0, match.index);
    if (!DRAWING_RECEIVER.test(before)) return [];
    const receiver = before.match(/[\w$]+(?:\(\s*\))?(?:\??\.[\w$]+(?:\(\s*\))?|\[[^\]]*\])*\??$/)?.[0] ?? '';
    const start = match.index - receiver.length;
    const writes = writesThrough(code.slice(match.index + match[0].length))
      || /\bdelete\s+$/.test(code.slice(0, start))
      || /\bObject\.assign\s*\(\s*$/.test(code.slice(0, start));
    return writes ? [`${match[1]}:${lineOf(code, match.index)}`] : [];
  });
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
const BINDINGS = /\s(v-model(?:\.[\w.]+)?|@[\w.:-]+|v-on:[\w.:-]+)="([^"]*)"/g;

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

const webModules = (url = new URL('./', `file://${JS}/`), prefix = '') => readdirSync(url, { withFileTypes: true }).flatMap((entry) => {
  const path = `${prefix}${entry.name}`;
  if (entry.isDirectory()) return webModules(new URL(`${entry.name}/`, url), `${path}/`);
  return entry.name.endsWith('.js') ? [path] : [];
});

test('only a slice writes its members: no module and no template binding writes one', () => {
  const found = webModules().filter((file) => !file.startsWith('services/drawing/'))
    .flatMap((file) => memberWriteSites(readFileSync(join(JS, file), 'utf8')).map((site) => `${file} ${site}`));
  assert.deepEqual(found, []);
  const html = readFileSync(join(ROOT, 'gbdraw', 'web', 'index.html'), 'utf8');
  const template = templateMemberWrites(html).map((site) => site.replace(/:\d+$/, ''));
  assert.deepEqual(template.filter((site) => !TEMPLATE_WRITES_BEFORE_M3_13.includes(site)), []);
  assert.deepEqual(TEMPLATE_WRITES_BEFORE_M3_13.filter((site) => !template.includes(site)), [], 'remove the bindings M3.13 replaced');
});

test('the slice write scan flags a planted write and leaves reads alone', () => {
  const [member, other] = MEMBERS;
  const writes = memberWriteSites([
    `const edit = (drawing, state, otherDrawing) => {`,
    `  drawing.${member}.push(row);`,
    `  state.activeDrawing().${member} = [];`,
    `  delete otherDrawing?.${member}[key];`,
    `  Object.assign(state.drawings.circular.${other}, next);`,
    `  state.drawings[mode].${other}.value.top = 1;`,
    `};`,
    `const read = (drawing, saved) => {`,
    `  // drawing.${member}.push(row) in prose is no writer`,
    `  saved.${member} = drawing.${member}.map((row) => row);`,
    `  return drawing.${member}.length === 0 || drawing.${member}.includes(row);`,
    `};`
  ].join('\n'), [member, other]);
  assert.deepEqual(writes, [`${member}:2`, `${member}:3`, `${member}:4`, `${other}:5`, `${other}:6`]);
  const rows = (inner) => `<div v-for="set in ${member}"><div v-for="(item, index) in set.rows">${inner}</div></div>`;
  assert.deepEqual(templateMemberWrites(rows(`<input :value="item.label" @input="setLabel(set, item, $event.target.value)">`), [member]), []);
  assert.equal(templateMemberWrites(rows('<input v-model="item.label">'), [member]).length, 1);
  assert.equal(templateMemberWrites(rows('<button @click="set.rows.splice(index, 1)"></button>'), [member]).length, 1);
  assert.equal(templateMemberWrites(`<input v-model.number="${member}.top">`, [member]).length, 1);
  assert.equal(templateMemberWrites('<div v-for="item in list"><input v-model="item.label"></div>', [member]).length, 0);
});

// 5. The handle's members are reachable only through slice-handle.js, which
// only the slices import.
test('only the slices import slice-handle.js, and its member key is written there only', () => {
  const handleModule = 'services/drawing/slice-handle.js';
  const key = readFileSync(join(JS, handleModule), 'utf8').match(/Symbol\.for\('([^']+)'\)/)[1];
  const importers = [];
  const keyHolders = [];
  for (const file of webModules()) {
    const source = readFileSync(join(JS, file), 'utf8');
    if (importsOf(source).some((target) => join(dirname(join(JS, file)), target) === join(JS, handleModule))) importers.push(file);
    if (file !== handleModule && source.includes(key)) keyHolders.push(file);
  }
  assert.deepEqual(importers.filter((file) => !file.startsWith('services/drawing/')), []);
  assert.deepEqual(keyHolders, []);
});

// 6. A slice's install, its view, and Save then Load carry an edit of its own
// mode only (M6 isolation, node half): each slice holds a non-default sample
// in the forms its writers and Load give, so a no-op install fails.
const SAMPLES = {
  visibilityRules: () => ({ featureVisibilityManualRules: [normalizeFeatureVisibilityRule({
    id: 'sample-rule', featureType: 'CDS', qualifier: 'locus_tag', value: '^A$', action: 'off'
  })] }),
  annotationSets: () => ({ annotationSets: normalizeAnnotationSets([{ id: 'sample-set', legendLabel: 'Sample', annotations: [
    { id: 'region_1', target: { kind: 'coordinateSpan', start: 10, end: 20 }, label: 'Sample region' }
  ] }]) }),
  // Load keeps a Linear group's name only for a group of the saved Result
  // (`orthogroupEditsOf`); the fixture holds none.
  groupNames: (mode) => ({
    orthogroupNameOverrides: mode === 'linear' ? {} : { 'sample-group': 'Sample name' },
    orthogroupDescriptionOverrides: mode === 'linear' ? {} : { 'sample-group': 'Sample description' },
    orthogroupDormantOverrides: { '["sample-member"]': { name: 'Dormant name', description: 'Dormant description' } }
  }),
  canvas: () => ({ canvasPadding: { top: 12, right: 3, bottom: 7, left: 5 } })
};
// A view member as plain data: a ref's value, a set's values.
const plain = (value) => {
  if (value instanceof Set) return [...value];
  if (value && typeof value === 'object' && !Array.isArray(value) && Object.keys(value).length === 1 && 'value' in value) return plain(value.value);
  return value;
};
const captures = (mode) => JSON.stringify(SLICES.map(([name, slice]) => slice.capture(state.drawings[mode].slices[name])));
// Each Save names its own file: Save asks before downloading a file name again.
let saves = 0;
const saveSession = async () => {
  saves += 1;
  const saved = await exportSession(`slices-${saves}`);
  assert.equal(saved.status, 'saved', JSON.stringify(saved));
  return new Blob([gunzipSync(Buffer.from(await saved.blob.arrayBuffer()))], { type: 'application/json' });
};

test('a slice edit shows in its drawing, survives Save and Load, and leaves the other drawing unchanged', async () => {
  assert.deepEqual(Object.keys(SAMPLES).sort(), SLICES.map(([name]) => name).sort(), 'every slice has a sample');
  assert.equal((await load(fixtureNames[0]))?.status, 'ok');
  for (const mode of MODES) {
    const other = MODES.find((candidate) => candidate !== mode);
    for (const [name, slice] of SLICES) {
      const drawing = state.drawings[mode];
      const sample = SAMPLES[name](mode);
      const json = JSON.stringify(sample);
      const untouched = captures(other);
      slice.install(drawing.slices[name], sample);
      assert.equal(JSON.stringify(slice.capture(drawing.slices[name])), json, `${mode} ${name}: install writes the sample`);
      assert.equal(JSON.stringify(Object.fromEntries(slice.members.map((member) => [member, plain(drawing[member])]))), json,
        `${mode} ${name}: the drawing shows it`);
      assert.equal(captures(other), untouched, `${mode} ${name}: the ${other} drawing is unchanged`);
      const session = await saveSession();
      slice.reset(drawing.slices[name]);
      assert.notEqual(JSON.stringify(slice.capture(drawing.slices[name])), json, `${mode} ${name}: the sample is no default`);
      const loaded = await importSession({ target: { files: [session], value: 'selected' } },
        { transformLegacyResultSvg, confirmLoadRepairs: async () => true });
      assert.equal(loaded?.status, 'ok', `${mode} ${name}: ${loaded?.error?.message ?? loaded?.status}`);
      assert.equal(JSON.stringify(slice.capture(state.drawings[mode].slices[name])), json, `${mode} ${name}: Save then Load`);
      assert.equal(captures(other), untouched, `${mode} ${name}: Load leaves the ${other} drawing as saved`);
      slice.reset(state.drawings[mode].slices[name]);
    }
  }
});
