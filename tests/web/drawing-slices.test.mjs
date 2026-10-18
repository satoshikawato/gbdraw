// Typed drawing slices (plan-r2 §3.2, R15): a drawing's editor edits belong to
// slices (gbdraw/web/js/services/drawing/), which alone write them; its
// settings belong to a slice too, which installs, captures, and resets them
// while their controls keep writing them (plan-r2 §13 item 3). The guard
// asserts, without a member or slice literal:
//   1. each slice's view holds exactly its members, every drawing shows each
//      slice's view, and no member is in two slices (the drawing's members
//      are the slices' by construction, helpers/drawing-state.mjs;
//      tests/web/drawing-context.test.mjs checks them against state.js);
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
import { memberWriteSites, refValueWriteSites, templateMemberWrites } from './helpers/drawing-writes.mjs';
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

const {
  exportSession, importSession, installShownDrawing, setUnmanagedConfigOverrideValidator
} = await import('../../gbdraw/web/js/services/config.js');
const { normalizeAnnotationSets } = await import('../../gbdraw/web/js/services/annotation-state.js');
const { normalizeFeatureVisibilityRule } = await import('../../gbdraw/web/js/services/feature-visibility.js');
const { normalizeSpecificRule } = await import('../../gbdraw/web/js/services/specific-color-rules.js');
const { normalizePaletteColors } = await import('../../gbdraw/web/js/utils/color-utils.js');
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

// 4. Only a slice writes its members: the slice write scan
// (helpers/drawing-writes.mjs) finds no write in any module outside
// services/drawing/ or in the template, and flags planted ones.
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
  assert.deepEqual(templateMemberWrites(html), []);
});

test('the slice write scan flags a planted write and leaves reads alone', () => {
  const [member, other] = DRAWING_EDITOR_KEYS;
  const writes = memberWriteSites([
    `const edit = (drawing, state, otherDrawing) => {`,
    `  drawing.${member}.push(row);`,
    `  state.activeDrawing().${member} = [];`,
    `  delete otherDrawing?.${member}[key];`,
    `  Object.assign(state.drawings.circular.${other}, next);`,
    `  state.drawings[mode].${other}.value.top = 1;`,
    `  drawing.${member}.count++;`,
    `  --state.activeDrawing().${other}.value.top;`,
    `  state.activeDrawing()[key].value = value;`,
    `  delete drawing[key];`,
    `  state.drawings[state.mode.value][key].value = 1;`,
    `  state.drawings.circular[side].top = 1;`,
    `};`,
    `const read = (drawing, saved) => {`,
    `  // drawing.${member}.push(row) in prose is no writer`,
    `  saved.${member} = drawing.${member}.map((row) => row);`,
    `  return drawing.${member}.length === 0 || drawing.${member}.includes(row) || drawing[key] === state.drawings[mode];`,
    `};`
  ].join('\n'), [member, other]);
  assert.deepEqual(writes, [`${member}:2`, `${member}:3`, `${member}:4`, `${other}:5`, `${other}:6`, `${member}:7`, `${other}:8`,
    '[computed]:9', '[computed]:10', '[computed]:11', '[computed]:12']);
  // A helper that holds a member under its name (the M3.9 blind spot).
  assert.deepEqual(memberWriteSites([
    `const prune = ({ ${member} = {}, ${other} }, key) => {`,
    `  delete ${member}[key];`,
    `  ${other}.value = [];`,
    `  Object.assign(${member}, next);`,
    `  ${member} = {};`,
    `  return ${member}[key] === ${other}.value;`,
    `};`
  ].join('\n'), [member, other]), [`${member}:2`, `${other}:3`, `${member}:4`]);
  const rows = (inner) => `<div v-for="set in ${member}"><div v-for="(item, index) in set.rows">${inner}</div></div>`;
  assert.deepEqual(templateMemberWrites(rows(`<input :value="item.label" @input="setLabel(set, item, $event.target.value)">`), [member]), []);
  assert.equal(templateMemberWrites(rows('<input v-model="item.label">'), [member]).length, 1);
  assert.equal(templateMemberWrites(rows('<button @click="set.rows.splice(index, 1)"></button>'), [member]).length, 1);
  assert.equal(templateMemberWrites(`<input v-model.number="${member}.top">`, [member]).length, 1);
  assert.equal(templateMemberWrites(`<color-field v-model:color="${member}.top"></color-field>`, [member]).length, 1);
  assert.equal(templateMemberWrites(`<button @click="${member}.count++"></button>`, [member]).length, 1);
  assert.equal(templateMemberWrites('<div v-for="item in list"><input v-model="item.label"></div>', [member]).length, 0);
});

// 4b. The shown drawing (plan-r2 A2): `state.activeDrawingId` names it, and
// `installShownDrawing` (services/config.js) is its one writer; the shown
// mode is the shown drawing's own (`drawing.mode`), so nothing writes `mode`.
// The write scan sees a `.value` write through the ref's name and the
// templates; it does not see a write through an alias held under another
// name (`const m = state.mode; m.value = ...`).
test('the shown drawing is drawings[activeDrawingId], mode follows it, and one writer installs it', () => {
  const shown = state.activeDrawingId.value;
  for (const id of MODES) {
    installShownDrawing(id);
    assert.equal(state.activeDrawingId.value, id);
    assert.equal(state.activeDrawing(), state.drawings[id]);
    assert.equal(state.drawings[id].id, id, 'a drawing is keyed by its ID');
    assert.equal(state.mode.value, state.drawings[id].mode, 'the shown mode is the shown drawing\'s');
  }
  installShownDrawing(shown);
  assert.throws(() => installShownDrawing('none'), /Unknown drawing/);
  const writer = 'services/config.js';
  const source = readFileSync(join(JS, writer), 'utf8');
  const start = source.indexOf('export const installShownDrawing = ');
  assert.ok(start >= 0, `${writer} exports installShownDrawing`);
  const [first, last] = [start, source.indexOf('\n};', start)].map((index) => source.slice(0, index).split('\n').length);
  const names = ['activeDrawingId', 'mode'];
  const found = webModules().flatMap((file) => refValueWriteSites(readFileSync(join(JS, file), 'utf8'), names)
    .filter((site) => !(file === writer && site.startsWith('activeDrawingId:')
      && Number(site.split(':')[1]) >= first && Number(site.split(':')[1]) <= last))
    .map((site) => `${file} ${site}`));
  assert.deepEqual(found, [], 'install the shown drawing with installShownDrawing');
  assert.ok(refValueWriteSites(source, names).length > 0, 'installShownDrawing writes the ID');
  const html = readFileSync(join(ROOT, 'gbdraw', 'web', 'index.html'), 'utf8');
  assert.deepEqual(templateMemberWrites(html, names), []);
});

test('the shown-drawing write scan finds each write form and leaves reads alone', () => {
  assert.deepEqual(refValueWriteSites([
    'state.mode.value = next;',
    'mode.value = next; // mode.value = prose',
    'if (state.mode.value === next || mode.value !== other) return;',
    'state.activeDrawingId.value = id;',
    'drawing.mode === next; linearComparisonPlan.mode = "none"; const mode = next;'
  ].join('\n'), ['activeDrawingId', 'mode']), ['mode:1', 'mode:2', 'activeDrawingId:4']);
  assert.equal(templateMemberWrites('<button @click="mode = \'linear\'"></button>', ['mode']).length, 1);
  assert.equal(templateMemberWrites('<settings-section-index :mode="mode"></settings-section-index>', ['mode']).length, 0);
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
  colorRules: () => ({
    manualSpecificRules: [normalizeSpecificRule({ feat: 'CDS', qual: 'product', val: 'kinase', color: '#123456', cap: 'Kinase', fromFile: true })]
  }),
  // Load gives a palette's colors the comparison colors (`normalizePaletteColors`).
  palette: () => ({
    selectedPalette: 'mint', currentColors: normalizePaletteColors({ CDS: '#123456', tRNA: '#654321' }),
    pendingPaletteName: 'queued', pendingPaletteColors: normalizePaletteColors({ CDS: '#abcdef' })
  }),
  labelTables: () => ({
    labelTextBulkOverrides: { 'hypothetical protein': 'HP' },
    canonicalLabelOverrideRows: [{ recordId: '*', featureType: 'CDS', qualifier: 'product', valueRegex: '^kinase$', labelText: 'Kinase' }],
    manualPriorityRules: [{ feat: 'CDS', order: 'gene,product' }],
    filterMode: 'Whitelist',
    manualBlacklist: 'putative',
    manualWhitelist: [{ feat: 'CDS', qual: 'product', key: 'kinase' }]
  }),
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
  canvas: () => ({ canvasPadding: { top: 12, right: 3, bottom: 7, left: 5 } }),
  // Load keeps the rows `canonicalFeatureOverrides` and the placement reader
  // of the drawing's mode accept.
  // Load reads the Legend rows as `normalizeEditorStateData` gives them.
  paint: () => ({
    featureColorOverrides: { '["record-1","CDS_1"]': { color: '#123456', caption: 'Kinase' } },
    featureStrokeOverrides: { 'record-1\u0000CDS_1': { strokeColor: '#654321', strokeWidth: 2 } },
    legendEntries: [{ caption: 'Renamed CDS', originalCaption: 'CDS', color: '#abcdef', featureIds: [] }],
    deletedLegendEntries: [{ caption: 'tRNA', originalCaption: 'tRNA', color: '#fedcba', featureIds: [] }],
    dormantLegendEntries: [{ caption: 'Coverage', originalCaption: 'depth', color: '#111111', featureIds: [] }],
    legendColorOverrides: { CDS: '#abcdef' },
    legendStrokeOverrides: { CDS: { strokeColor: '#222222', strokeWidth: 1.5 } },
    addedLegendCaptions: ['Kinase']
  }),
  // A drawing's settings as Load gives them, with values no default holds.
  settings: (mode) => {
    const settings = DRAWING_SLICES.settings.capture(state.drawings[mode].slices.settings);
    return {
      ...settings,
      form: { ...settings.form, plot_title: `Sample ${mode} title` },
      losatProgram: 'tblastx',
      linearRecordLayoutEnabled: !settings.linearRecordLayoutEnabled
    };
  },
  featureEdits: (mode) => ({
    featureOverrides: { '["record-1","CDS_1"]': {
      recordKey: 'record-1', biologicalFeatureId: 'CDS_1', featureVisibility: 'off', labelVisibility: null, labelText: 'Sample', labelSourceText: null
    } },
    featurePlacementOverrides: { '["record-1","CDS_2"]': {
      recordKey: 'record-1', biologicalFeatureId: 'CDS_2', placement: { kind: 'lane', side: mode === 'linear' ? 'above' : 'outward', level: 1 }
    } }
  })
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
      const loaded = JSON.stringify(slice.capture(drawing.slices[name]));
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
      const reloaded = await importSession({ target: { files: [session], value: 'selected' } },
        { transformLegacyResultSvg, confirmLoadRepairs: async () => true });
      assert.equal(reloaded?.status, 'ok', `${mode} ${name}: ${reloaded?.error?.message ?? reloaded?.status}`);
      assert.equal(JSON.stringify(slice.capture(state.drawings[mode].slices[name])), json, `${mode} ${name}: Save then Load`);
      assert.equal(captures(other), untouched, `${mode} ${name}: Load leaves the ${other} drawing as saved`);
      // Back to the loaded fixture, which a Save and Load keep as it is (Load
      // resolves the defaults Reset gives, such as the track placements).
      slice.install(state.drawings[mode].slices[name], JSON.parse(loaded));
    }
  }
});

// 7. Reset Settings returns the settings a Reset owns and keeps the record
// gap and rows, the OFF records, and the comparison edges that hold a file
// or a LOSAT file name, inactive (services/drawing/settings.js `reset`).
test('Reset Settings keeps the record layout, the OFF records, and the comparison files, inactive', async () => {
  const { resetSettings } = await import('../../gbdraw/web/js/services/reset.js');
  const settings = DRAWING_SLICES.settings;
  const handle = state.drawings.linear.slices.settings;
  const before = settings.capture(handle);
  const file = new File(['pair'], 'pair.tsv');
  const edge = (id, fields) => ({ id, queryUid: 'a', subjectUid: 'b', source: 'upload', included: true, fileActive: true,
    losatFilenameActive: true, file: null, losatFilename: '', ...fields });
  settings.install(handle, {
    form: { ...before.form, plot_title: 'Edited' }, losatProgram: 'blastp', linearRecordGap: 40,
    linearRecordRows: [{ uid: 'card-a', row: 2 }], recordsOff: ['card-b'],
    linearComparisonPlan: { mode: 'custom', defaultSource: 'upload', edges: [
      edge('with-file', { file }), edge('with-name', { losatFilename: 'pair.losat' }), edge('empty', { losatFilename: ' ' })
    ] }
  });
  resetSettings(state);
  const after = settings.capture(handle);
  assert.deepEqual([after.form.plot_title, after.losatProgram], ['', 'blastn']);
  assert.deepEqual([after.linearRecordGap, after.linearRecordRows, after.recordsOff], [40, [{ uid: 'card-a', row: 2 }], ['card-b']]);
  assert.deepEqual([after.linearComparisonPlan.mode, after.linearComparisonPlan.defaultSource], ['none', 'losat']);
  assert.deepEqual(after.linearComparisonPlan.edges.map(({ id, included, fileActive, losatFilenameActive }) => (
    [id, included, fileActive, losatFilenameActive])), [['with-file', false, false, false], ['with-name', false, false, false]]);
  assert.equal(after.linearComparisonPlan.edges[0].file, file, 'an edge keeps its file');
  settings.install(handle, before);
});
