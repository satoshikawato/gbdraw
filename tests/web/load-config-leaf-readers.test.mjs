// E1a (PR-E): a configuration value Load refused, or kept invalid, is read by
// the reader table with a dialog row before the dialog, in a Session 46 slice,
// a Session 27-44 draft and a configuration file; the loaded state saves.
import assert from 'node:assert/strict';
import { test } from 'node:test';
import { readFileSync } from 'node:fs';
import { gunzipSync } from 'node:zlib';
import { File } from 'node:buffer';
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
setUnmanagedConfigOverrideValidator(async (payload) => ({ result: await validateWebConfigOverrides(payload) }));
const { transformLegacyResultSvg } = await import('../../gbdraw/web/js/app/app-setup.js');
const { state } = await import('../../gbdraw/web/js/state.js');

const fixture = (name) => {
  const bytes = readFileSync(new URL(`../fixtures/sessions/${name}`, import.meta.url));
  return JSON.parse(name.endsWith('.gz') ? gunzipSync(bytes) : bytes);
};
const gallery = (name) => JSON.parse(readFileSync(new URL(`../../gbdraw/web/gallery/sessions/${name}`, import.meta.url), 'utf8'));
const VECTORS = JSON.parse(readFileSync(new URL('../fixtures/session_admission_vectors.json', import.meta.url), 'utf8'));
const input = (value) => ({
  target: { files: [new File([JSON.stringify(value)], 'session.json', { type: 'application/json' })], value: 'selected' }
});
const BASES = {
  v30: () => fixture('BGC0000708-BGC0000713.v30.gbdraw-session.json.gz'),
  v40: () => {
    const base = VECTORS.bases['web-40'];
    return base.document ? structuredClone(base.document) : fixture(base.fixture.replace('tests/fixtures/sessions/', ''));
  },
  v44: () => fixture('HmmtDNA_basic_circular.v44-schema8.gbdraw-session.json.gz'),
  v44L: () => fixture('lambda_basic_linear.v44-schema8.gbdraw-session.json.gz'),
  v46: () => fixture('legend-renamed-rows.v46.gbdraw-session.json.gz'),
  g47L: () => gallery('lambda_basic_linear.gbdraw-session.json')
};
const set = (value, path, leaf) => {
  const parts = path.split('.');
  const container = parts.slice(0, -1).reduce((current, part) => current[part], value);
  container[parts.at(-1)] = leaf;
  return value;
};
const load = async (file) => {
  const loaded = await importSession(input(file), { transformLegacyResultSvg, confirmLoadRepairs: async () => true });
  assert.ok(['ok', 'legacy'].includes(loaded.status), JSON.stringify(loaded.error));
  return loaded.repairs.map(({ code, path }) => `${code} ${path}`);
};
const saves = async () => {
  const result = await exportSession('E1a', { askChoice: async () => 'download' });
  assert.equal(result?.status, 'saved', JSON.stringify(result?.error ?? result));
};
// [base, the edited path, the stored value, the rows]: a base Session's
// value replaced, then loaded and saved (a Session 27-39 is saved only after
// Generate, so it is loaded only).
const CASES = [
  // OV-757, A17c: a default type's shape of no name keeps its default shape.
  ['v44', 'config.adv.feature_shapes.CDS', 'hexagon', ['setting-replaced config.adv.feature_shapes.CDS']],
  ['g47L', 'modes.linear.config.adv.feature_shapes.CDS', 'hexagon', ['setting-replaced modes.linear.config.adv.feature_shapes.CDS']],
  // A17d: an arrow ratio out of range is the default; its text reads as its number.
  ['g47L', 'modes.linear.config.adv.arrow_shaft_width_ratio', 2, ['setting-replaced modes.linear.config.adv.arrow_shaft_width_ratio']],
  ['v44', 'config.adv.arrow_head_length_ratio', -1, ['setting-replaced config.adv.arrow_head_length_ratio']],
  ['v44', 'config.adv.arrow_shaft_width_ratio', '0.5', []],
  // A17e: an option value of no choice is the default; an older name is read as its current name.
  ['v44L', 'config.form.linear_track_layout', 'zigzag', ['setting-replaced config.form.linear_track_layout']],
  ['v44L', 'config.form.linear_track_layout', 'spreadout', []],
  ['g47L', 'modes.linear.config.adv.label_placement', 'inside', ['setting-replaced modes.linear.config.adv.label_placement']],
  ['v44', 'config.adv.multi_record_size_mode', 'huge', ['setting-replaced config.adv.multi_record_size_mode']],
  ['g47L', 'modes.linear.config.losat.blastp.mode', 'bogus', ['setting-replaced modes.linear.config.losat.blastp.mode']],
  ['g47L', 'modes.linear.config.losat.blastp.hitLimitsByMode', null, []],
  ['g47L', 'modes.linear.config.losat.blastp.collinearInferOrthogroups', 'x',
    ['setting-replaced modes.linear.config.losat.blastp.collinearInferOrthogroups']],
  // A17h.
  ['v44', 'config.adv.feature_overlap_tolerance_bp', -1, ['setting-replaced config.adv.feature_overlap_tolerance_bp']],
  // A17f: a collection of another kind, and an entry that is no object.
  ['v44', 'config.rules', 'x', ['setting-replaced config.rules']],
  ['v44', 'config.rules', [5], ['setting-entries-dropped config.rules']],
  ['v44', 'config.losat', 'x', ['setting-replaced config.losat']],
  ['v44', 'config.unmanagedConfigOverrides', 'x', ['setting-replaced config.unmanagedConfigOverrides']],
  // A17b, TB19: a key no gbdraw writer wrote is dropped with a row.
  ['v44', 'config.adv.zzz', 1, ['setting-replaced config.adv.zzz']],
  ['g47L', 'modes.linear.config.form.zzz', 1, ['setting-replaced modes.linear.config.form.zzz']],
  ['g47L', 'modes.linear.config.zzz', 1, ['setting-replaced modes.linear.config.zzz']],
  ['v44', 'config.rules', [{ feat: 'CDS', qual: 'product', val: 'x', color: '#ff0000', cap: 'x', zzz: 1 }],
    ['setting-replaced config.rules[0].zzz']],
  // A17b': a value only the other mode's drawing draws.
  ['v46', 'modes.circular.config.form.linear_track_layout', 'zigzag', []],
  // A17g: a retired name reads as its current name at every version.
  ['v44', 'config.adv.depth_tick_interval', 5, []],
  ['g47L', 'modes.linear.config.adv.depth_tick_interval', 5, []],
  ['v44', 'config.adv.linear_show_accession', false, []],
  // OV-764: the drawing's lists and maps read per entry.
  ['g47L', 'modes.linear.config.recordsOff', [5], ['setting-entries-dropped modes.linear.config.recordsOff']],
  ['g47L', 'modes.linear.config.recordsOff', ['zzz', 'zzz'], []],
  ['v44L', 'config.recordsOff', 'x', ['setting-replaced config.recordsOff']],
  ['g47L', 'modes.linear.config.linearComparisonPlan', 'x', ['setting-replaced modes.linear.config.linearComparisonPlan']],
  ['g47L', 'modes.linear.config.linearComparisonPlan.mode', 'bogus', ['setting-replaced modes.linear.config.linearComparisonPlan.mode']],
  ['v44L', 'config.linearComparisonPlan.edges', [5], ['setting-replaced config.linearComparisonPlan.edges[0]']],
  ['g47L', 'modes.linear.config.losat.blastp.hitLimitsByMode', { pairwise: { candidateLimit: 'x' } },
    ['setting-replaced modes.linear.config.losat.blastp.hitLimitsByMode.pairwise.candidateLimit']],
  ['v44L', 'config.featurePlacementOverrides', { a: 5 }, ['setting-entries-dropped config.featurePlacementOverrides']],
  ['g47L', 'modes.linear.config.unmanagedConfigOverrides', { 'zzz.yyy': 1 },
    ['setting-replaced modes.linear.config.unmanagedConfigOverrides["zzz.yyy"]']],
  // A23': a drawing's feature edits read per row.
  ['g47L', 'modes.linear.features.featureOverrides', 'x', ['setting-replaced modes.linear.features.featureOverrides']],
  ['g47L', 'modes.linear.features.featureOverrides', { a: 5 }, ['setting-entries-dropped modes.linear.features.featureOverrides']],
  // OV-639, OV-611: a record display row read per row, its selector text first.
  ['v44', 'config.recordDisplayDrafts', {}, ['setting-replaced config.recordDisplayDrafts']],
  ['v40', 'config.recordDisplayDrafts', [{ scope: 'circular', sourceUid: 'circular', selector: ['#1'], recordId: 'same',
    topologyOverride: null, startCoordinate: 25, reverseComplementOverride: null, anchorIntent: null }],
  ['setting-entries-dropped config.recordDisplayDrafts']]
];

for (const [base, path, value, rows] of CASES) {
  test(`E1a: ${base} ${path} = ${JSON.stringify(value)} loads with ${rows.length ? rows.join(', ') : 'no row'}`, async () => {
    assert.deepEqual(await load(set(BASES[base](), path, value)), rows);
    if (!['v30', 'v40'].includes(base)) await saves();
  });
}

test('E1a: a retired name and an older option value are read as their current names', async () => {
  const session = BASES.v44L();
  delete session.config.adv.depth_large_tick_interval;
  session.config.adv.depth_tick_interval = 7;
  session.config.form.linear_track_layout = 'spreadout';
  assert.deepEqual(await load(session), []);
  assert.equal(state.drawings.linear.adv.depth_large_tick_interval, 7);
  assert.equal(state.drawings.linear.form.linear_track_layout, 'above');
});

// TB32: the keys of a Session 27-39 draft no gbdraw writer wrote.
test('E1a: a Session 30 lists the keys of its draft no gbdraw writer wrote (TB32)', async () => {
  const session = BASES.v30();
  session.config.bogusCfg = 1;
  session.config.form = { ...session.config.form, bogus_form: 1 };
  session.config.adv = { ...session.config.adv, bogus_adv: 1 };
  assert.deepEqual((await load(session)).sort(), [
    'setting-replaced config.adv.bogus_adv', 'setting-replaced config.bogusCfg', 'setting-replaced config.form.bogus_form'
  ]);
});

// A retired name stored with its current name: the current one wins, no row.
test('E1a: a retired name stored with its current name gives way to it', async () => {
  const session = BASES.v44L();
  session.config.adv.depth_large_tick_interval = 9;
  session.config.adv.depth_tick_interval = 7;
  assert.deepEqual(await load(session), []);
  assert.equal(state.drawings.linear.adv.depth_large_tick_interval, 9);
});

// The configuration-file path (E-OD04): the same rows at the file's own places.
test('E1a: a configuration file lists what its drawing cannot read, and saves', async () => {
  // The configuration applies to the shown drawing, which holds a Result.
  assert.deepEqual(await load(BASES.g47L()), []);
  const config = structuredClone(BASES.g47L().modes.linear.config);
  delete config.unmanagedConfigOverrides;
  Object.assign(config, {
    zzz: 1, recordDisplayDrafts: {}, recordsOff: [5], unmanagedConfigOverrides: { 'zzz.yyy': 1 },
    linearComparisonPlan: { mode: 'none', defaultSource: 'losat', edges: [5] }
  });
  config.adv = { ...config.adv, zzz: 1, arrow_shaft_width_ratio: 2 };
  state.mode.value = 'linear';
  assert.deepEqual((await load(config)).sort(), [
    'setting-entries-dropped recordsOff',
    'setting-replaced adv.arrow_shaft_width_ratio',
    'setting-replaced adv.zzz',
    'setting-replaced linearComparisonPlan.edges[0]',
    'setting-replaced recordDisplayDrafts',
    'setting-replaced unmanagedConfigOverrides["zzz.yyy"]',
    'setting-replaced zzz'
  ]);
  await saves();
});

// OV-788: a configuration file reads its track slot leaves as a Session does
// (OV-758): a value Load cannot use is listed, not refused.
test('a configuration file lists a track slot value it cannot use (OV-788)', async () => {
  const config = structuredClone(BASES.v46().modes.linear.config);
  config.adv.linear_track_slots[0].height = -7.5;
  state.mode.value = 'linear';
  assert.deepEqual(await load(config), ['setting-replaced adv.linear_track_slots[0].height']);
});
