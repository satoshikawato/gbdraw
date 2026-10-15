// E1b (PR-E): a value Load changed without a row, kept invalid, or refused is
// read before the dialog with a row (D1b section 1, E1a residuals); the loaded
// state saves (a Session 27-39 is saved only after Generate, so it is loaded only).
import assert from 'node:assert/strict';
import { test } from 'node:test';
import { readFileSync } from 'node:fs';
import { execFileSync } from 'node:child_process';
import { fileURLToPath } from 'node:url';
import { gunzipSync } from 'node:zlib';
import { File } from 'node:buffer';
import { installSessionImportWorker } from './helpers/session-import-node.mjs';
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

const REPO = new URL('../../', import.meta.url);
const { exportSession, importSession, setUnmanagedConfigOverrideValidator } = await import('../../gbdraw/web/js/services/config.js');
setUnmanagedConfigOverrideValidator((payload) => ({ result: JSON.parse(execFileSync('python', ['-c', `
import json, sys
from gbdraw.web_support.config_overrides import validate_web_config_overrides_json
p = json.load(sys.stdin)
print(validate_web_config_overrides_json(p['mode'], json.dumps(p['config']),
    json.dumps(p['configOverrides']), json.dumps(p['managedPaths']), p['source']))
`], { input: JSON.stringify(payload), encoding: 'utf8', cwd: fileURLToPath(REPO) })) }));
const { transformLegacyResultSvg } = await import('../../gbdraw/web/js/app/app-setup.js');
const { describeLoadRepairs } = await import('../../gbdraw/web/js/utils/load-repairs.js');
const { state } = await import('../../gbdraw/web/js/state.js');
const { orderedConservationSources } = await import('../../gbdraw/web/js/services/conservation-series.js');

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
  v44: () => fixture('HmmtDNA_basic_circular.v44-schema8.gbdraw-session.json.gz'),
  v44L: () => fixture('lambda_basic_linear.v44-schema8.gbdraw-session.json.gz'),
  so46: () => structuredClone(VECTORS.bases['settings-only-46'].document),
  g47L: () => gallery('lambda_basic_linear.gbdraw-session.json'),
  // tests/utils/two_mode_session.py: a Circular Session 46 with a Linear `otherModeResult`.
  twoMode: () => JSON.parse(execFileSync('python', ['-c',
    'import json; from tests.utils.two_mode_session import two_mode_session; print(json.dumps(two_mode_session()))'
  ], { encoding: 'utf8', cwd: fileURLToPath(REPO), maxBuffer: 1 << 30 }))
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
const saved = async () => {
  const result = await exportSession('E1b', { askChoice: async () => 'download' });
  assert.equal(result?.status, 'saved', JSON.stringify(result?.error ?? result));
  return JSON.parse(gunzipSync(Buffer.from(await result.blob.arrayBuffer())));
};
const LINEAR_SEQ = 'webFiles.bindings.linearSeqs.0';

// [base, the edited path, the stored value, the rows, a check of the loaded state]
const CASES = [
  // OV-568: a Linear card's crop and genetic code, read once by the Files reader.
  ['so46', `${LINEAR_SEQ}.region_start`, 'abc', ['setting-replaced webFiles.bindings.linearSeqs[0].region_start'],
    () => assert.equal(state.linearSeqs[0].region_start, null)],
  ['so46', `${LINEAR_SEQ}.region_end`, [1], ['setting-replaced webFiles.bindings.linearSeqs[0].region_end'],
    () => assert.equal(state.linearSeqs[0].region_end, null)],
  ['so46', `${LINEAR_SEQ}.losat_gencode`, 'x', ['setting-replaced webFiles.bindings.linearSeqs[0].losat_gencode'],
    () => assert.equal(state.linearSeqs[0].losat_gencode, 1)],
  // An orientation that is no boolean is not used (D1b: changed, U false), its text included.
  ['so46', `${LINEAR_SEQ}.region_reverse`, 'true', ['setting-replaced webFiles.bindings.linearSeqs[0].region_reverse'],
    () => assert.equal(state.linearSeqs[0].region_reverse, false)],
  ['so46', `${LINEAR_SEQ}.definition`, 5, ['setting-replaced webFiles.bindings.linearSeqs[0].definition'],
    () => assert.equal(state.linearSeqs[0].definition, '')],
  ['v30', 'files.linearSeqs.0.losat_gencode', 0, ['setting-replaced files.linearSeqs[0].losat_gencode']],
  // OV-595: the PNG DPI is one of the choices the DPI menu lists.
  ['g47L', 'ui.downloadDpi', 'abc', ['setting-replaced ui.downloadDpi'], () => assert.equal(Number(state.downloadDpi.value), 300)],
  ['g47L', 'ui.downloadDpi', '600', [], () => assert.equal(Number(state.downloadDpi.value), 600)],
  // OV-594: an input type of no kind; the request's own type wins without a row.
  ['g47L', 'ui.cInputType', 5, ['setting-replaced ui.cInputType'], () => assert.equal(state.cInputType.value, 'gb')],
  ['g47L', 'ui.lInputType', 5, [], () => assert.equal(state.lInputType.value, 'gb')],
  ['v44', 'ui.cInputType', 'fasta', [], () => assert.equal(state.cInputType.value, 'gb')],
  // K19: a settings-only Session's mode of no mode (the request's mode wins in a full one).
  ['so46', 'ui.mode', 'zigzag', ['setting-replaced ui.mode'], () => assert.equal(state.mode.value, 'circular')],
  ['g47L', 'ui.mode', 'zigzag', [], () => assert.equal(state.mode.value, 'linear')],
  // RL-D: a drawing's feature color edits and label rows read per entry.
  ['g47L', 'modes.linear.features.featureColorOverrides', { a: 5, b: 'notacolor', c: { color: '#ff0000', caption: 'x' } },
    ['setting-entries-dropped modes.linear.features.featureColorOverrides'],
    () => assert.deepEqual(state.drawings.linear.featureColorOverrides, { c: { color: '#ff0000', caption: 'x' } })],
  // OV-866: an edit map of another kind (a list) is dropped as its row says, not read by index.
  ['g47L', 'modes.linear.features.featureColorOverrides', ['#ff0000'],
    ['setting-replaced modes.linear.features.featureColorOverrides'],
    () => assert.deepEqual(state.drawings.linear.featureColorOverrides, {})],
  ['g47L', 'modes.linear.features.labelOverrideRows', [5, 'x', { recordId: '*', featureType: 'CDS', qualifier: 'product', valueRegex: 'x', labelText: 'y' }],
    ['setting-entries-dropped modes.linear.features.labelOverrideRows'],
    () => assert.equal(state.drawings.linear.canonicalLabelOverrideRows.value.length, 1)],
  // OV-578: the alignment plan and record translations of a settings-only file.
  ['so46', 'modes.linear.config.linearRecordLayout.similarityAlignment', 'junk',
    ['setting-replaced modes.linear.config.linearRecordLayout.similarityAlignment'], () => assert.equal(state.similarityAlignmentPlan.value, null)],
  ['so46', 'modes.linear.config.linearRecordLayout.recordTranslations', 'nope',
    ['setting-replaced modes.linear.config.linearRecordLayout.recordTranslations'], () => assert.deepEqual(state.linearRecordTranslations.value, [])],
  // E1a residuals: LOSAT execution settings are no slice field of a Session 46 (E-OD03).
  ['g47L', 'modes.linear.config.losat.executionMode', 'serial', ['setting-replaced modes.linear.config.losat.executionMode'],
    () => assert.equal(state.losatExecution.executionMode, 'auto')],
  // An older Linear visibility that is no boolean is dropped with a row on every path.
  ['v30', 'config.adv.linear_show_accession', 'yes', ['setting-replaced config.adv.linear_show_accession'],
    () => assert.equal(state.drawings.linear.adv.linear_accession_visibility, 'show')],
  ['v44', 'config.adv.linear_show_length', 'yes', ['setting-replaced config.adv.linear_show_length']],
  // Unknown fields of a qualifier priority row, an annotation set and a comparison edge (E-OD03, E-OD12).
  ['v44', 'config.qualifierPriorityRules', [{ feat: 'CDS', order: 'product', zzz: 1 }], ['setting-replaced config.qualifierPriorityRules[0].zzz']],
  ['v44', 'config.annotationSets', [{ id: 'a', annotations: [], zzz: 1 }], ['setting-replaced config.annotationSets[0].zzz']],
  // The two Linear font sizes are read as settings (a positive number or unset).
  ['g47L', 'modes.linear.config.adv.scale_font_size', 'big', ['setting-replaced modes.linear.config.adv.scale_font_size'],
    () => assert.equal(state.drawings.linear.adv.scale_font_size, null)],
  ['g47L', 'modes.linear.config.adv.ruler_label_font_size', '14', [], () => assert.equal(Number(state.drawings.linear.adv.ruler_label_font_size), 14)],
  // A7i at a Session 46: recordsOff of another kind, an entry of no key, an unbound uid.
  ['g47L', 'modes.linear.config.recordsOff', 'x', ['setting-replaced modes.linear.config.recordsOff']],
  ['so46', 'modes.linear.config.recordsOff', [5], ['setting-entries-dropped modes.linear.config.recordsOff']],
  ['so46', 'modes.linear.config.recordsOff', ['unbound'], [], () => assert.deepEqual([...state.drawings.linear.recordsOff], [])],
  // OV-822: a null over a settings group merges nothing (E-OD06), and the group keeps its defaults.
  ['g47L', 'modes.linear.config.losat.blastp', null, [], () => assert.equal(typeof state.drawings.linear.losat.blastp, 'object')],
  ['so46', 'modes.circular.config.losat.blastp', null, []]
];

for (const [base, path, value, rows, check] of CASES) {
  test(`E1b: ${base} ${path} = ${JSON.stringify(value)} loads with ${rows.length ? rows.join(', ') : 'no row'}`, async () => {
    assert.deepEqual(await load(set(BASES[base](), path, value)), rows);
    check?.();
    if (base !== 'v30') await saved();
  });
}

// K14 (D1c 4.3): record display drafts of no list are read before the dialog
// by the reader the apply uses, so Load applies what the dialog lists and
// nothing acts after the commit point. The dialog's wording is the owner's.
test('K14: record display drafts of no list are listed before the dialog and applied as listed', async () => {
  const configFile = () => {
    const config = structuredClone(BASES.g47L().modes.linear.config);
    state.mode.value = 'linear';
    return config;
  };
  const cases = [
    ['g47L', () => BASES.g47L(), 'modes.linear.config.recordDisplayDrafts', 'linear'],
    ['v44', () => BASES.v44(), 'config.recordDisplayDrafts', 'circular'],
    ['configuration file', configFile, 'recordDisplayDrafts', 'linear']
  ];
  await load(BASES.g47L());
  for (const [name, base, path, mode] of cases) {
    for (const value of [{}, 'x', 5, true]) {
      const drawing = state.drawings[mode];
      drawing.recordDisplayDrafts.splice(0, drawing.recordDisplayDrafts.length, {
        sourceUid: 'held', selector: '#1', recordId: 'held', topologyOverride: null, startCoordinate: null,
        reverseComplementOverride: null, anchorIntent: null
      });
      const dialogs = [];
      const loaded = await importSession(input(set(base(), path, value)), {
        transformLegacyResultSvg, confirmLoadRepairs: async (dialog) => { dialogs.push(dialog.items); return true; }
      });
      const label = `${name} ${JSON.stringify(value)}`;
      assert.ok(['ok', 'legacy'].includes(loaded.status), `${label}: ${JSON.stringify(loaded.error)}`);
      assert.deepEqual(loaded.repairs.map(({ code, path: at }) => `${code} ${at}`), [`setting-replaced ${path}`], label);
      assert.deepEqual(dialogs, [describeLoadRepairs(loaded.repairs)?.items], label);
      assert.ok(dialogs[0].some((line) => line.includes(`(${path})`)), `${label}: ${dialogs[0]}`);
      assert.deepEqual([...drawing.recordDisplayDrafts], [], label);
    }
  }
  await saved();
});

// OV-892 (E-OD26): a stored Label blacklist that is no text keeps today's
// text (`String(value)`, nothing for null), listed when it is an object, a list
// or a boolean; a number reads as its text with no row. One table row serves a
// Session slice, a Session 44 draft and a configuration file.
test('OV-892: a Label blacklist that is no text is listed and keeps its text', async () => {
  const configFile = () => {
    const config = structuredClone(BASES.g47L().modes.linear.config);
    state.mode.value = 'linear';
    return config;
  };
  const cases = [
    ['g47L', () => BASES.g47L(), 'modes.linear.config.blacklistText', 'linear'],
    ['v44', () => BASES.v44(), 'config.blacklistText', 'circular'],
    ['configuration file', configFile, 'blacklistText', 'linear']
  ];
  await load(BASES.g47L());
  for (const [name, base, path, mode] of cases) {
    for (const [value, text, listed] of [
      [{ a: 1 }, '[object Object]', true], [['x', 'y'], 'x,y', true], [true, 'true', true], [false, '', true],
      [5, '5', false], [null, '', false]
    ]) {
      const dialogs = [];
      const loaded = await importSession(input(set(base(), path, value)), {
        transformLegacyResultSvg, confirmLoadRepairs: async (dialog) => { dialogs.push(dialog.items); return true; }
      });
      const label = `${name} ${JSON.stringify(value)}`;
      assert.ok(['ok', 'legacy'].includes(loaded.status), `${label}: ${JSON.stringify(loaded.error)}`);
      assert.deepEqual(loaded.repairs.map(({ code, path: at }) => `${code} ${at}`), listed ? [`setting-replaced ${path}`] : [], label);
      assert.deepEqual(dialogs, listed ? [describeLoadRepairs(loaded.repairs)?.items] : [], label);
      assert.equal(state.drawings[mode].manualBlacklist.value, text, label);
    }
  }
});

test('E1b: an unknown field of a stored comparison edge is dropped with a row (E-OD12)', async () => {
  const session = BASES.g47L();
  session.modes.linear.config.linearComparisonPlan.edges = [{
    id: 'e1', queryUid: 'a', subjectUid: 'b', included: false, fileActive: false, losatFilenameActive: false,
    source: 'losat', losatFilename: '', zzz: 1
  }];
  assert.deepEqual(await load(session), ['setting-replaced modes.linear.config.linearComparisonPlan.edges[0].zzz']);
  assert.equal(Object.hasOwn(state.drawings.linear.linearComparisonPlan.edges[0], 'zzz'), false);
  await saved();
});

test('OV-655, TB23: a single file slot keeps the first file it reads and lists the rest and a malformed entry', async () => {
  let session = BASES.v30();
  const first = session.files.linearSeqs[0].gb;
  session.files.linearSeqs[0].gb = [first, { ...first, name: 'second.gbk' }];
  assert.deepEqual(await load(session), ['file-entry-dropped files.linearSeqs[0].gb']);
  assert.equal(state.linearSeqs[0].gb?.name, first.name);
  session = BASES.v30();
  session.files.linearSeqs[1].gb = { bogus: 1 };
  assert.deepEqual(await load(session), ['file-entry-dropped files.linearSeqs[1].gb']);
});

test('OV-633: the other Result set\'s ui and editor edits are read with rows', async () => {
  const session = BASES.twoMode();
  session.otherModeResult.ui.appliedPaletteColors = { CDS: 5, gene: '#00ff00' };
  session.otherModeResult.ui.generatedMultiRecordCanvas = 'yes';
  assert.deepEqual(await load(session), [
    'setting-replaced otherModeResult.ui.generatedMultiRecordCanvas',
    'setting-entries-dropped otherModeResult.ui.appliedPaletteColors'
  ]);
  await saved();
});

test('A20: a Result set whose Results Load cannot compose is dropped with its row before the dialog (E-OD16)', async () => {
  const malformed = '<svg xmlns="http://www.w3.org/2000/svg"><g></svg>';
  // The shown set: its Results and the members that name them; the settings load.
  let session = BASES.v44();
  session.results[0].content = malformed;
  assert.deepEqual(await load(session), ['session-part-dropped results']);
  assert.equal(state.results.value.length, 0);
  await saved();
});

test('OV-789: a Session 27-39 Result Load cannot admit drops the Results with their row before the dialog; the settings load (E-OD17)', async () => {
  const malformed = '<svg xmlns="http://www.w3.org/2000/svg"><g></svg>';
  for (const name of ['scale-interval-negative-linear-cli.v30.gbdraw-session.json.gz', 'BGC0000708-BGC0000713.v39.gbdraw-session.json.gz']) {
    // A Result that is no SVG, and one the composition transform cannot compose.
    for (const content of ['x', malformed]) {
      const session = fixture(name);
      session.results[0].content = content;
      session.ui.selectedResultIndex = 0;
      assert.deepEqual(await load(session), ['session-part-dropped results'], `${name} ${content}`);
      assert.equal(state.results.value.length, 0);
      assert.equal(state.biologicalFeatures.value.length, 0);
      assert.equal(state.mode.value, session.ui.mode);
    }
  }
});

test('OV-822: a configuration file\'s null LOSAT group loads with no row, and saves', async () => {
  assert.deepEqual(await load(BASES.g47L()), []);
  const config = structuredClone(BASES.g47L().modes.linear.config);
  config.losat = { blastp: null };
  state.mode.value = 'linear';
  assert.deepEqual(await load(config), []);
  assert.equal(typeof state.drawings.linear.losat.blastp, 'object');
  await saved();
});

test('A configuration file: a Legend position outside its mode is listed (receipt residual)', async () => {
  const session = BASES.g47L();
  await load(session);
  const config = structuredClone(session.modes.linear.config);
  config.form.legend = 'upper_right';
  assert.deepEqual(await load(config), ['legend-position-defaulted form.legend']);
  await saved();
});

// E4d (review round 1, E-OD29). R1-B1: Circular comparison series read the
// sequence lists by place (`sourceIndex`), so a dropped binding keeps its place
// empty, as a Depth list does (E-OD28), and no later file changes series.
test('R1-B1: a dropped comparison sequence binding keeps its place, so no later file changes series', async () => {
  const session = fixture('synthetic_conservation.gbdraw-session.json.gz');
  session.webFiles.bindings.c_conservation_fastas[1].resourceId = 'zz';
  session.webFiles.bindings.c_conservation_sequence_sources[1].resourceId = 'zz';
  assert.deepEqual(await load(session), [
    'file-entry-dropped webFiles.bindings.c_conservation_fastas[1]',
    'file-entry-dropped webFiles.bindings.c_conservation_sequence_sources[1]'
  ]);
  const names = (files) => files.map((file) => file?.name ?? null);
  assert.deepEqual(names(state.files.c_conservation_fastas), ['comparison-b.fasta', null, 'comparison-d.fasta']);
  assert.deepEqual(names(state.files.c_conservation_sequence_sources), ['comparison-b.fasta', null, 'comparison-d.fasta']);
});

// OV-893 (E-OD30): a comparison series reads its BLAST table by place too, and
// its readers take no empty place: a dropped table binding drops its series
// with its files, as Remove does, and the series after it keep their own.
test('OV-893: a dropped comparison table binding drops its series with its files; later series keep theirs', async () => {
  const session = fixture('synthetic_conservation.gbdraw-session.json.gz');
  session.webFiles.bindings.c_conservation_blasts[1].resourceId = 'zz';
  const dialogs = [];
  const result = await importSession(input(session), {
    transformLegacyResultSvg, confirmLoadRepairs: async (dialog) => { dialogs.push(...dialog.items); return true; }
  });
  assert.deepEqual(result.repairs.map(({ code, path }) => `${code} ${path}`), [
    'file-entry-dropped webFiles.bindings.c_conservation_blasts[1]',
    'comparison-series-dropped config.circularConservation.series[1]'
  ]);
  assert.ok(dialogs.some((line) => line.includes('"comparison-c"')), JSON.stringify(dialogs));
  const names = (files) => files.map((file) => file?.name ?? null);
  assert.deepEqual(names(state.files.c_conservation_blasts),
    ['comparison-b.circular_conservation.losatn.tsv', 'comparison-d.circular_conservation.losatn.tsv']);
  assert.deepEqual(names(state.files.c_conservation_fastas), ['comparison-b.fasta', 'comparison-d.fasta']);
  assert.deepEqual(names(state.files.c_conservation_sequence_sources), ['comparison-b.fasta', 'comparison-d.fasta']);
  const conservation = state.drawings.circular.circularConservation;
  assert.deepEqual(conservation.series.map(({ label, sourceIndex }) => [label, sourceIndex]), [['comparison-b', 0], ['comparison-d', 1]]);
  assert.deepEqual(orderedConservationSources(state.files.c_conservation_blasts, conservation)
    .map(({ label, file, sourceIndex }) => [label, file.name, state.files.c_conservation_sequence_sources[sourceIndex]?.name]), [
    ['comparison-b', 'comparison-b.circular_conservation.losatn.tsv', 'comparison-b.fasta'],
    ['comparison-d', 'comparison-d.circular_conservation.losatn.tsv', 'comparison-d.fasta']
  ]);
});

// An older Session (27-30) keeps its comparison files in `files` lists: an
// entry that holds a file Load cannot read is listed, as a binding is; an
// empty place is no row.
test('a comparison file entry of an older Session that cannot be read is listed', async () => {
  const session = fixture('precomputed-ring-unlabelled-cli.v30.gbdraw-session.json.gz');
  session.files.c_conservation_fastas = [5, null];
  session.files.c_conservation_sequence_sources = [{ name: 'x' }];
  assert.deepEqual((await load(session)).filter((row) => row.startsWith('file-entry-dropped')), [
    'file-entry-dropped files.c_conservation_fastas[0]',
    'file-entry-dropped files.c_conservation_sequence_sources[0]'
  ]);
  assert.deepEqual(state.files.c_conservation_fastas, [null, null]);
});

// OV-898 (E-OD35): a 27-30 draft without series names its rings by `labels`,
// by place; the labels stand in for the series, so a dropped place leaves out
// its label too, with a row when it had a label or an entry Load cannot read,
// and the rings after it keep their own labels.
test('OV-898: a dropped table of a draft without series leaves out its label; later rings keep theirs', async () => {
  const rings = (entry) => {
    const session = fixture('precomputed-ring-unlabelled-cli.v30.gbdraw-session.json.gz');
    const [table] = session.files.c_conservation_blasts;
    session.files.c_conservation_blasts = [entry, { ...table, name: 'second.tsv' }];
    session.config.circularConservation.labels = 'a,b';
    return session;
  };
  const drawn = () => orderedConservationSources(state.files.c_conservation_blasts, state.drawings.circular.circularConservation)
    .map(({ label, file }) => [label, file.name]);
  for (const entry of [5, null, {}]) {
    const result = await importSession(input(rings(entry)), { transformLegacyResultSvg, confirmLoadRepairs: async () => true });
    assert.deepEqual(result.repairs.filter(({ code }) => code === 'comparison-series-dropped').map(({ path, label }) => [path, label]),
      [['files.c_conservation_blasts[0]', 'a']], JSON.stringify(entry));
    assert.deepEqual(drawn(), [['b', 'second.tsv']], JSON.stringify(entry));
    assert.equal(state.drawings.circular.circularConservation.labels, 'b');
  }
  // An empty place with no label is no row; the labels are kept by place.
  const unlabelled = rings(null);
  unlabelled.files.c_conservation_blasts.push(null);
  assert.deepEqual((await load(unlabelled)).filter((row) => row.startsWith('comparison-series-dropped')),
    ['comparison-series-dropped files.c_conservation_blasts[0]']);
  assert.deepEqual(drawn(), [['b', 'second.tsv']]);
  // A place past the series of a draft with series keeps its label by place too.
  const past = rings(null);
  const [first] = fixture('precomputed-ring-unlabelled-cli.v30.gbdraw-session.json.gz').files.c_conservation_blasts;
  past.files.c_conservation_blasts = [first, null, past.files.c_conservation_blasts[1]];
  past.config.circularConservation.series = [{ label: 'A', fileName: first.name }];
  past.config.circularConservation.labels = 'A,b,c';
  const pastRows = (await importSession(input(past), { transformLegacyResultSvg, confirmLoadRepairs: async () => true })).repairs;
  assert.deepEqual(pastRows.filter(({ code }) => code === 'comparison-series-dropped').map(({ path, label }) => [path, label]),
    [['files.c_conservation_blasts[1]', 'b']]);
  assert.deepEqual(drawn(), [['A', first.name], ['c', 'second.tsv']]);
});

// OV-896: with no series left, the comparison is off, as Remove leaves it.
test('OV-896: a Load that leaves out every comparison series turns the comparison off', async () => {
  const session = fixture('synthetic_conservation.gbdraw-session.json.gz');
  session.webFiles.bindings.c_conservation_blasts.forEach((binding) => { binding.resourceId = 'zz'; });
  const result = await importSession(input(session), { transformLegacyResultSvg, confirmLoadRepairs: async () => true });
  assert.deepEqual(result.repairs.filter(({ code }) => code === 'comparison-series-dropped').map(({ path }) => path),
    [2, 1, 0].map((index) => `config.circularConservation.series[${index}]`));
  const conservation = state.drawings.circular.circularConservation;
  assert.deepEqual([conservation.series.length, state.files.c_conservation_blasts.length], [0, 0]);
  assert.equal(conservation.enabled, false);
});

// OV-899 (E-OD34): a dropped Linear comparison file binding leaves the
// comparison as the pair card's Remove leaves it (no file, none active), and
// the committed request names no resource the file lacks, so Save writes the
// loaded state and the saved file loads with no dialog.
test('OV-899: a dropped comparison file binding loads as a removed file; Save and reload succeed', async () => {
  const { setMainSessionComparisonFrameConverter } = await import('../../gbdraw/web/js/services/main-session-comparison-frame.js');
  // A Session 42's stored tables go through the Worker frame conversion; the table here is in the search frame.
  setMainSessionComparisonFrameConverter(async ({ tableText }) => ({ result: { tsv: tableText } }));
  try {
    const q42 = () => fixture('q-frame-main-web-upload.v42.gbdraw-session.json.gz');
    const unknownId = q42();
    unknownId.webFiles.bindings.linearComparisons[0].file.resourceId = 'zz';
    const removed = q42();
    delete removed.resources['comparison-nucleotide-1'];
    for (const [name, session, rows, comparisons] of [
      ['unknown resource id', unknownId, ['file-entry-dropped webFiles.bindings.linearComparisons[0].file'],
        ['comparison-nucleotide-1']],
      ['resource removed', removed, [
        'file-entry-dropped webFiles.bindings.linearComparisons[0].file',
        'request-reference-dropped renderRequest.comparisons[0]'
      ], []]
    ]) {
      assert.deepEqual(await load(session), rows, name);
      assert.deepEqual(state.drawings.linear.linearComparisonPlan.edges.map(({ file, fileActive, included, source }) => (
        { file, fileActive, included, source })), [{ file: null, fileActive: false, included: true, source: 'upload' }], name);
      const document = await saved();
      assert.deepEqual(document.webFiles.bindings.linearComparisons ?? [], [], name);
      assert.deepEqual(document.renderRequest.comparisons.map(({ resourceId }) => resourceId), comparisons, name);
      comparisons.forEach((id) => assert.ok(document.resources[id], `${name}: ${id}`));
      assert.deepEqual(await load(document), [], name);
    }
  } finally {
    setMainSessionComparisonFrameConverter(null);
  }
});

// R2-B1: a null request of a full Session reads as an absent one (E-OD29):
// the settings-only row names everything that is left out.
test('R2-B1: a full Session 44 or 47 whose request is null loads as settings only with the request row', async () => {
  for (const base of ['v44', 'g47L']) {
    const absent = BASES[base]();
    delete absent.renderRequest;
    const expected = await load(absent);
    assert.deepEqual(expected, ['request-dropped renderRequest'], base);
    assert.deepEqual(await load(set(BASES[base](), 'renderRequest', null)), expected, base);
    assert.equal(state.results.value.length, 0);
  }
  // A settings-only file that holds Results is read the same way.
  const results = BASES.g47L().results;
  assert.deepEqual(await load({ ...BASES.so46(), results }), ['request-dropped renderRequest']);
  assert.deepEqual(await load(BASES.so46()), []);
});

// R2-M5: an editor Legend of another kind loses the Result's original Legend
// order and colors; Load keeps today's outcome (E-OD08) and lists it.
test('R2-M5: a top-level editorState.legend that is no object is listed', async () => {
  const session = set(fixture('legend-renamed-rows.v46.gbdraw-session.json.gz'), 'editorState.legend', 'x');
  assert.deepEqual(await load(session), ['setting-replaced editorState.legend']);
  assert.equal((state.originalLegendOrder?.value ?? []).length, 0);
});

// E4f m2 (E-OD36): an empty value of an optional editor map ([] as {} or
// null) holds no edit, so it loses nothing and is no row; a value of another
// kind still is.
test('E4f m2: an empty editor map is no row; a map of another kind is listed', async () => {
  const session = (path, value) => set(fixture('legend-renamed-rows.v46.gbdraw-session.json.gz'), path, value);
  for (const empty of [[], {}, null]) {
    assert.deepEqual(await load(session('editorState.legend', empty)), [], JSON.stringify(empty));
  }
  const withLegend = (key, value) => {
    const file = fixture('legend-renamed-rows.v46.gbdraw-session.json.gz');
    file.editorState.legend = { ...(file.editorState.legend ?? {}), [key]: value };
    return file;
  };
  assert.deepEqual(await load(withLegend('originalColors', [])), []);
  assert.deepEqual(await load(withLegend('originalColors', 'x')), ['setting-replaced editorState.legend.originalColors']);
  for (const other of ['x', 5]) {
    assert.deepEqual(await load(session('editorState.legend', other)), ['setting-replaced editorState.legend'], String(other));
  }
});

// R2-m2: the verb agrees with each part's name (G1: the wording owner only).
test('R2-m2: a dropped Session part is named with the verb its name takes', () => {
  const text = JSON.stringify(describeLoadRepairs(['resources', 'editorState', 'modes', 'otherModeResult']
    .map((path) => ({ code: 'session-part-dropped', path, effect: 'dropped' }))));
  for (const plural of ['the stored input files (resources)', 'the saved diagram (editorState)', 'both drawings (modes)']) {
    assert.match(text, new RegExp(`${plural.replace(/[()]/g, '\\$&')} could not be read and are left out`), plural);
  }
  assert.match(text, /other mode \(otherModeResult\) could not be read and is left out/);
});
