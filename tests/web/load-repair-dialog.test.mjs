// BA-D25 (Load repair dialog): Load reads a file as far as it can and, when
// something cannot be fully restored, asks once before anything changes.
// `utils/load-repairs.js` owns the wording (G1); the Load port answers.
import assert from 'node:assert/strict';
import { test } from 'node:test';
import { readFileSync, readdirSync } from 'node:fs';
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
globalThis.document = {};
installFakeSvgDom();
globalThis.File = File;
globalThis.alert = () => {};
installSessionImportWorker();

const { importSession } = await import('../../gbdraw/web/js/services/config.js');
// A Session older than 40 has its Results transformed through the composition root's port.
const { transformLegacyResultSvg } = await import('../../gbdraw/web/js/app/app-setup.js');
const { state } = await import('../../gbdraw/web/js/state.js');
const { describeLoadRepairs, loadRepair } = await import('../../gbdraw/web/js/utils/load-repairs.js');

const fixture = (name) => JSON.parse(gunzipSync(readFileSync(new URL(`../fixtures/sessions/${name}`, import.meta.url))));
const input = (session) => ({
  target: { files: [new File([JSON.stringify(session)], 'session.json', { type: 'application/json' })], value: 'selected' }
});
// What a Load replaces: the settings and edits of both drawings, the shown
// mode and Results, the title, the input files.
const loadedState = () => JSON.stringify({
  mode: state.mode.value,
  title: state.sessionTitle.value,
  results: state.results.value.map((result) => result.content),
  files: Object.keys(state.files),
  drawings: Object.fromEntries(['circular', 'linear'].map((mode) => {
    const drawing = state.drawings[mode];
    return [mode, {
      form: drawing.form, adv: drawing.adv, featureOverrides: drawing.featureOverrides,
      legendEntries: drawing.legendEntries.value
    }];
  }))
});

// A Linear Session 44 with one Feature visibility edit that names no feature
// of its saved catalog: the feature-edit migration drops it.
const unmatchedEditSession = () => {
  const session = fixture('feature-edits-crop-rc.v44.gbdraw-session.json.gz');
  session.features.featureVisibilityOverrides.f00000000_record_2 = 'off';
  return session;
};

test('a lossy Load asks once before it changes anything, and Cancel leaves everything as it was (PLAN 12.7 (1), (2))', async () => {
  const before = await importSession(input(fixture('lambda_basic_linear.v44-schema8.gbdraw-session.json.gz')), {
    confirmLoadRepairs: async () => assert.fail('a file Load restores fully asks nothing')
  });
  assert.equal(before.status, 'ok', JSON.stringify(before.error));
  const snapshot = loadedState();

  const dialogs = [];
  const reached = [];
  const event = input(unmatchedEditSession());
  const canceled = await importSession(event, {
    confirmLoadRepairs: async (dialog) => {
      dialogs.push(dialog);
      // Nothing is applied while the question is open.
      assert.equal(loadedState(), snapshot);
      return false;
    },
    rollbackState: { capture: () => reached.push('snapshot'), restore: () => reached.push('restore') },
    afterImport: async () => reached.push('afterImport'),
    beforePreviewMount: () => { reached.push('preview'); return null; },
    installLoadedArtifactSlots: () => reached.push('slots')
  });
  assert.equal(canceled.status, 'canceled');
  assert.deepEqual(dialogs, [{
    title: 'Load this file?',
    description: 'Some parts of this file cannot be fully restored:',
    items: ['1 feature edit(s) from an older Session could not be matched to a feature of its saved diagram and were dropped.'],
    note: 'Save Session afterwards writes the current format. The file itself is not changed.'
  }]);
  assert.equal(loadedState(), snapshot);
  assert.deepEqual(reached, []);
  assert.equal(event.target.value, '');

  const loaded = await importSession(input(unmatchedEditSession()), { confirmLoadRepairs: async () => true });
  assert.equal(loaded.status, 'ok', JSON.stringify(loaded.error));
  assert.deepEqual(loaded.repairs.map(({ code, count }) => ({ code, count })), [{ code: 'feature-edit-unmatched', count: 1 }]);
  assert.deepEqual(loaded.notices, []);
  assert.notEqual(loadedState(), snapshot);
});

// OV-652 (BA-D11, BA-D25): a stored legend position outside its mode reads
// as the mode's default or as unset and is listed; case and outer spaces are
// no row.
test('a legend position outside its mode is listed at Load (OV-652)', async () => {
  const session = fixture('lambda_basic_linear.v44-schema8.gbdraw-session.json.gz');
  session.ui.generatedLegendPosition = 'upper_left';
  session.ui.layoutPreferences = {
    circular: { single: { legend: ' Upper_Left ' }, multi: { legend: null } }, linear: { legend: 'banana' }
  };
  const dialogs = [];
  const loaded = await importSession(input(session), { confirmLoadRepairs: async ({ items }) => { dialogs.push(items); return true; } });
  assert.equal(loaded.status, 'ok', JSON.stringify(loaded.error));
  assert.deepEqual(dialogs, [[
    'Legend: the position "upper_left" (ui.generatedLegendPosition) is not one of its mode\'s positions; bottom is used.',
    'Legend: the position "banana" (ui.layoutPreferences.linear.legend) is not one of its mode\'s positions; it is not used.'
  ]]);
  assert.equal(state.generatedLegendPosition.value, 'bottom');
  assert.equal(state.drawings.linear.form.legend, 'bottom');
  assert.equal(state.drawings.circular.layoutPreferences.circular.single.legend, 'upper_left');
});

test('the repair list is frozen at the dialog: no step after the commit adds a repair (G2)', async () => {
  const loaded = await importSession(input(unmatchedEditSession()), { confirmLoadRepairs: async () => true });
  assert.equal(loaded.status, 'ok', JSON.stringify(loaded.error));
  assert.ok(Array.isArray(loaded.repairs) && Object.isFrozen(loaded.repairs));
});

// BA-D15 (OV-569, OV-710): a legacy configuration file is read before it is
// applied; malformed record display settings are dropped and listed, and
// Cancel leaves the drawing as it was.
test('a legacy configuration with malformed record display settings loads with the key dropped and listed (PLAN 12.7 (3))', async () => {
  const legacy = { form: { legend: 'left' }, recordDisplayDrafts: [{ junk: 1 }, 5, 'x'] };
  const drawing = state.drawings[state.mode.value === 'linear' ? 'linear' : 'circular'];
  drawing.recordDisplayDrafts.splice(0, drawing.recordDisplayDrafts.length);
  const formBefore = JSON.stringify(drawing.form);
  const dialogs = [];
  const canceled = await importSession(input(legacy), { confirmLoadRepairs: async (dialog) => { dialogs.push(dialog); return false; } });
  assert.equal(canceled.status, 'canceled');
  assert.equal(JSON.stringify(drawing.form), formBefore);
  assert.deepEqual(dialogs.map(({ items }) => items), [
    ['Record display (recordDisplayDrafts): the stored record display settings could not be read and were dropped.']
  ]);
  const loaded = await importSession(input(legacy), { confirmLoadRepairs: async () => true });
  assert.equal(loaded.status, 'legacy', JSON.stringify(loaded.error));
  assert.deepEqual(loaded.repairs.map(({ code, path }) => [code, path]), [['record-display-drafts-dropped', 'recordDisplayDrafts']]);
  assert.deepEqual([...drawing.recordDisplayDrafts], []);
  assert.equal(drawing.form.legend, 'left');
});

test('the dialog words each code once, in area order, at most 12 lines (12.3)', () => {
  assert.equal(describeLoadRepairs([]), null);
  const tables = [1, 2, 3, 4, 5, 6, 7].map((line) => loadRepair({ code: 'table-row-dropped', path: 'label-whitelist', effect: 'dropped', line }));
  const dialog = describeLoadRepairs([
    loadRepair({ code: 'table-row-joined', path: 'specific-colors', effect: 'changed', line: 3 }),
    ...tables,
    loadRepair({ code: 'feature-edit-unmatched', path: 'features', effect: 'dropped', count: 2 }),
    loadRepair({ code: 'legend-color-dropped', path: 'editorState.legend.entries', effect: 'dropped', label: 'X' }),
    loadRepair({ code: 'legend-color-dropped', path: 'editorState.legend.deletedEntries', effect: 'dropped', label: 'W' }),
    loadRepair({ code: 'legend-color-dropped', path: 'modes.linear.editorState.legend.entries', effect: 'dropped', label: 'X' })
  ]);
  assert.deepEqual(dialog.items, [
    'Legend: entries "X", "W" (deleted) had a color the app does not accept and were dropped.',
    '2 feature edit(s) from an older Session could not be matched to a feature of its saved diagram and were dropped.',
    'Specific colors: line 3 had extra cells, joined into the last column with one space.',
    'Label whitelist: lines 1, 2, 3, 4, 5, and 2 more lacked a required column and were dropped.'
  ]);
  // D-43 (OV-300): each Legend entry dropped for its color is named once, a
  // row outside the shown list is marked, and the names are capped; a caption
  // is free text quoted as other user text is.
  const legend = (path, label) => loadRepair({ code: 'legend-color-dropped', path: `editorState.legend.${path}`, effect: 'dropped', label });
  assert.deepEqual(describeLoadRepairs([
    ...['A', 'B', 'C', 'D'].map((caption) => legend('entries', caption)),
    legend('deletedEntries', 'A'), legend('dormantEntries', 'E'), legend('entries', 'F'), legend('entries', 'G')
  ]).items, ['Legend: entries "A", "B", "C", "D", "A" (deleted), and 3 more had a color the app does not accept and were dropped.']);
  assert.deepEqual(describeLoadRepairs([legend('dormantEntries', 'E')]).items,
    ['Legend: entry "E" (not drawn) had a color the app does not accept and was dropped.']);
  assert.deepEqual(describeLoadRepairs([legend('entries', "5' UTR")]).items,
    ['Legend: entry "5\' UTR" had a color the app does not accept and was dropped.']);
  const many = describeLoadRepairs(['default-colors', 'specific-colors', 'label-whitelist', 'qualifier-priority',
    'feature-visibility', 'label-overrides'].flatMap((path) => ['table-row-joined', 'table-row-dropped', 'table-row-invalid']
    .map((code) => loadRepair({ code, path, effect: 'dropped', line: 1 }))));
  assert.equal(many.items.length, 12);
  assert.equal(many.items[11], 'and 7 more');
  const long = loadRepair({ code: 'legend-color-dropped', path: `p${'.x'.repeat(150)}`, effect: 'dropped', label: 'c'.repeat(80), stored: 's'.repeat(80) });
  assert.equal(long.path.length, 200);
  assert.equal(long.label.length, 60);
  assert.equal(long.stored.length, 60);
});

// G1: the lossy Load texts have one owner. No other module words a repair.
test('no module outside utils/load-repairs.js words what Load could not restore (G1)', () => {
  const root = new URL('../../gbdraw/web/js/', import.meta.url);
  const modules = (dir) => readdirSync(new URL(dir, root), { withFileTypes: true }).flatMap((entry) => (entry.isDirectory()
    ? modules(`${dir}${entry.name}/`) : entry.name.endsWith('.js') ? [`${dir}${entry.name}`] : []));
  const texts = [
    'could not be matched to a feature of its saved diagram',
    'Feature visibility edit(s) from an older Session',
    'could not be renamed to the hash of a feature',
    'had a color the app does not accept',
    'joined into the last column',
    'lacked a required column',
    'had a color the table does not accept',
    'cannot be fully restored'
  ];
  const owners = modules('').filter((path) => texts.some((text) => readFileSync(new URL(path, root), 'utf8').includes(text)));
  assert.deepEqual(owners, ['utils/load-repairs.js']);
});

/** Load `session`, answering [Load]; returns the result and the dialogs' items. */
const loadAnswering = async (session) => {
  const dialogs = [];
  const loaded = await importSession(input(session), { confirmLoadRepairs: async ({ items }) => { dialogs.push(items); return true; } });
  assert.equal(loaded.status, 'ok', JSON.stringify(loaded.error));
  return { loaded, items: dialogs.flat() };
};
const codes = (loaded) => loaded.repairs.map(({ code }) => code);

// PLAN 12.7 (4), inventory B13 (#1051): a Session 27-44 Legend edit that a
// later edit of the same row overwrites is listed; in a Session 46 the figure
// showed the later edit, so it is dropped without a row (Q1).
test('a Legend edit a later edit of its row overwrites is listed for a Session 27-44 and not for a Session 46', async () => {
  const older = fixture('lambda_basic_linear.v44-schema8.gbdraw-session.json.gz');
  older.editorState.legend = {
    entries: [{ caption: 'X', originalCaption: 'K', color: '#aaaaaa' }], colorOverrides: { X: '#111111', K: '#222222' }
  };
  const { loaded, items } = await loadAnswering(older);
  assert.deepEqual(items, ['Legend: the color edit stored under "X" (editorState.legend.colorOverrides) sets the row "K",'
    + ' which a later edit also sets; the later edit is used and this one was dropped.']);
  assert.deepEqual(state.drawings.linear.legendColorOverrides, { K: '#222222' });
  assert.deepEqual(codes(loaded), ['legend-edit-overwritten']);

  const session46 = fixture('legend-renamed-rows.v46.gbdraw-session.json.gz');
  // `GC skew (+)` is the caption of the row `GC content`: both edits set it.
  session46.modes.circular.editorState.legend.colorOverrides['GC content'] = '#222222';
  const quiet = await importSession(input(session46), { confirmLoadRepairs: async () => assert.fail('no row (Q1)') });
  assert.equal(quiet.status, 'ok', JSON.stringify(quiet.error));
  assert.deepEqual(codes(quiet), ['legend-edit-overwritten-session46']);
});

// Inventory B39: the re-key of a rule row whose saved caption is read without
// its color ID (F2) keeps an edit the new key already has; the moved one is listed.
test('a Legend edit the rule row re-key cannot move is listed (B39)', async () => {
  const session = JSON.parse(JSON.stringify(fixture('rule-captions-with-color-id.v46.gbdraw-session.json.gz'))
    .split('Complex IV [#1f77b4]').join('Complex III'));
  session.modes.circular.editorState.legend.strokeOverrides['Complex IV'] = { strokeColor: '#123456', strokeWidth: 1 };
  const { loaded, items } = await loadAnswering(session);
  assert.deepEqual(items, ['Legend: the stroke edit stored under "Complex IV [#ffe700]"'
    + ' (modes.circular.editorState.legend.strokeOverrides) moves to the row "Complex IV", which has an edit of its own;'
    + ' that edit is used and this one was dropped.']);
  assert.deepEqual(codes(loaded), ['legend-edit-rekey-dropped']);
  assert.deepEqual(state.drawings.circular.legendStrokeOverrides['Complex IV'], { strokeColor: '#123456', strokeWidth: 1 });
});

// OV-719 (inventory B48): a Session 46 keeps the first Legend entry of a
// caption, a Session 47 the first of a row key; the entry dropped is listed.
test('a Legend entry that repeats an earlier one is dropped and listed, by caption up to Session 46 and by row key since 47 (OV-719)', async () => {
  const session46 = fixture('legend-renamed-rows.v46.gbdraw-session.json.gz');
  session46.modes.circular.editorState.legend.entries.push(
    { caption: 'tRNA', originalCaption: 'tRNA 2', color: '#111111' },
    { caption: 'Q', originalCaption: 'CDS', color: '#222222' }
  );
  const read46 = await loadAnswering(session46);
  assert.deepEqual(read46.items, ['Legend: entry "tRNA" repeats an earlier entry of its list and was dropped.']);
  assert.deepEqual(state.drawings.circular.legendEntries.value.map(({ caption }) => caption).filter((caption) => ['tRNA', 'Q'].includes(caption)),
    ['tRNA', 'Q']);

  const session47 = fixture('legend-renamed-rows.v46.gbdraw-session.json.gz');
  session47.version = 47;
  session47.modes.circular.editorState.legend.entries.push(
    { caption: 'tRNA', originalCaption: 'tRNA 2', color: '#111111' },
    { caption: 'Q', originalCaption: 'CDS', color: '#222222' }
  );
  const read47 = await loadAnswering(session47);
  assert.deepEqual(read47.items, ['Legend: entry "Q" repeats an earlier entry of its list and was dropped.']);
  assert.deepEqual(state.drawings.circular.legendEntries.value.filter(({ caption }) => caption === 'tRNA').length, 2);
});

// Inventory B50: an older Session whose feature information cannot be
// restored (here: no embedded GenBank input) lists it; the Features card
// still says so after the Load.
test('an older Session whose feature information cannot be restored lists it (B50)', async () => {
  const session = fixture('hash-color-rule-linear-rc.v30.gbdraw-session.json.gz');
  session.features.extractedFeatures = [];
  Object.keys(session.files).forEach((key) => { session.files[key] = null; });
  const dialogs = [];
  const loaded = await importSession(input(session), {
    transformLegacyResultSvg, confirmLoadRepairs: async ({ items }) => { dialogs.push(items); return true; }
  });
  assert.equal(loaded.status, 'ok', JSON.stringify(loaded.error));
  assert.deepEqual(dialogs, [['Features: the feature information of the saved diagram could not be restored from this file;'
    + ' the SVG preview stays, and the next Generate restores it.']]);
  assert.match(state.featureExtractionError.value.summary, /does not include embedded GenBank inputs/);
});

// OV-712: the sequences of the saved comparison matches that cannot be read
// again (here: a malformed GenBank input) are listed, not only reported as
// `degradedRecovery`.
test('an older Session whose match sequences cannot be read again lists it (OV-712)', async () => {
  const session = fixture('scale-interval-negative-linear-cli.v30.gbdraw-session.json.gz');
  session.files.linearSeqs[0].gb.data = Buffer.from('LOCUS broken\nORIGIN\n 1 zzzz\n').toString('base64');
  const dialogs = [];
  const loaded = await importSession(input(session), {
    transformLegacyResultSvg, confirmLoadRepairs: async ({ items }) => { dialogs.push(items); return true; }
  });
  assert.equal(loaded.status, 'ok', JSON.stringify(loaded.error));
  assert.equal(loaded.degradedRecovery, true);
  assert.ok(dialogs.flat().includes('Comparisons: the sequences of the saved comparison matches could not be read again from the'
    + ' input files; the SVG preview stays, and the next Generate reads them.'), JSON.stringify(dialogs));
  assert.deepEqual(loaded.repairs.filter(({ code }) => code === 'match-sequences-unrestored').map(({ path }) => path), ['files']);
});
