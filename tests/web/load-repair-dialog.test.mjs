// BA-D25 (Load repair dialog): Load reads a file as far as it can and, when
// something cannot be fully restored, asks once before anything changes.
// `utils/load-repairs.js` owns the wording (G1); the Load port answers.
import assert from 'node:assert/strict';
import { test } from 'node:test';
import { mkdtempSync, readFileSync, readdirSync } from 'node:fs';
import { execFileSync, spawnSync } from 'node:child_process';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
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

const { importSession, setUnmanagedConfigOverrideValidator } = await import('../../gbdraw/web/js/services/config.js');
// The Worker's typed override check, run by Python (as session-cli-compatibility does).
setUnmanagedConfigOverrideValidator((payload) => ({ result: JSON.parse(execFileSync('python', ['-c', `
import json, sys
from gbdraw.web_support.config_overrides import validate_web_config_overrides_json
p = json.load(sys.stdin)
print(validate_web_config_overrides_json(p['mode'], json.dumps(p['config']),
    json.dumps(p['configOverrides']), json.dumps(p['managedPaths']), p['requireUnmanagedOnly']))
`], { input: JSON.stringify(payload), encoding: 'utf8', cwd: new URL('../../', import.meta.url).pathname })) }));
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

// Inventory B51 (Q2): a Session older than 27 stays refused, and the error
// says it is from an older gbdraw. A file that is no gbdraw Session is not
// called an old one.
test('a Session older than 27 is refused as an old version (B51)', async () => {
  for (const version of [8, 13, 26]) {
    const refused = await importSession(input({ format: 'gbdraw-session', version }), {
      confirmLoadRepairs: async () => assert.fail('no dialog for a refusal')
    });
    assert.equal(refused.status, 'error');
    assert.deepEqual([refused.error.code, refused.error.context.reason], ['INPUT_INVALID', 'SESSION_TOO_OLD']);
    assert.match(refused.error.summary, /gbdraw 0\.11\.0 or older/);
  }
  const other = await importSession(input({ format: 'other', version: 26 }), { confirmLoadRepairs: async () => true });
  assert.equal(other.status, 'error');
  assert.notEqual(other.error.context.reason, 'SESSION_TOO_OLD');
});

// PLAN 12.2: Gallery publication never loads through the dialog; it reads the
// tables strictly, so a Session with a row Load would repair is refused.
test('Gallery publication refuses a Session with a table row Load would repair', () => {
  const out = join(mkdtempSync(join(tmpdir(), 'gbdraw-publish-')), 'session.json');
  const run = spawnSync(process.execPath, [
    new URL('../../tools/publish_gallery_session.mjs', import.meta.url).pathname, 'prepare',
    new URL('../fixtures/sessions/whitelist-tab-keyword.v39.gbdraw-session.json.gz', import.meta.url).pathname, out
  ], { encoding: 'utf8' });
  assert.equal(run.status, 1);
  assert.match(run.stderr, /TABLE_INVALID/);
});

// B37, B43: the settings' readers of `applyConfigData` lose an out-of-set
// value; Load lists it before it applies the draft. Case and outer spaces,
// aliases, and a sort are no row.
test('a stored setting its reader cannot restore is listed before the Load (B37, B43)', async () => {
  const session = fixture('lambda_basic_linear.v44-schema8.gbdraw-session.json.gz');
  Object.assign(session.config.adv, {
    label_rendering: 'sideways', circular_label_placement: ' Radial ', multi_record_min_radius_ratio: 5,
    multi_record_positions: [{ selector: '#2', row: 2 }, { selector: '#1', row: 1 }, { selector: '#2', row: 3 }]
  });
  session.config.losat = { ...session.config.losat, totalThreadBudget: 'auto', executionMode: 'turbo' };
  const dialogs = [];
  const loaded = await importSession(input(session), { confirmLoadRepairs: async ({ items }) => { dialogs.push(items); return true; } });
  assert.equal(loaded.status, 'ok', JSON.stringify(loaded.error));
  assert.deepEqual(dialogs, [[
    'Settings: the value "sideways" (config.adv.label_rendering) cannot be used; auto is used.',
    'Settings: the value "5" (config.adv.multi_record_min_radius_ratio) cannot be used; 0.55 is used.',
    'Settings: the value "turbo" (config.losat.executionMode) cannot be used; auto is used.',
    'Settings: 1 entry of config.adv.multi_record_positions could not be read and was dropped.'
  ]]);
  assert.equal(state.drawings.linear.adv.label_rendering, 'auto');
  // The circular settings are the Circular drawing's (the split, plan 4.2).
  const { adv } = state.drawings.circular;
  assert.equal(adv.circular_label_placement, 'radial');
  assert.equal(adv.multi_record_min_radius_ratio, 0.55);
  assert.deepEqual(adv.multi_record_positions, [{ selector: '#1', row: 1 }, { selector: '#2', row: 2 }]);
  assert.equal(state.losatExecution.executionMode, 'auto');
  assert.equal(state.losatExecution.totalThreadBudget, 'safe');
});

// B36, B43, B44: the app settings of `ui`, and a pending palette that
// Instant Preview drops.
test('a Session 46 app setting or pending palette Load cannot restore is listed (B36, B43, B44)', async () => {
  const session = fixture('legend-renamed-rows.v46.gbdraw-session.json.gz');
  session.ui = {
    ...session.ui, featurePanelTab: 'banana', paletteInstantPreviewEnabled: true,
    losatExecution: { ...session.ui?.losatExecution, threadsPerJob: '-2' }
  };
  session.modes.linear.ui = { ...session.modes.linear.ui, pendingPaletteName: 'Pastel' };
  const dialogs = [];
  const loaded = await importSession(input(session), { confirmLoadRepairs: async ({ items }) => { dialogs.push(items); return true; } });
  assert.equal(loaded.status, 'ok', JSON.stringify(loaded.error));
  assert.deepEqual(dialogs, [[
    'Settings: the value "banana" (ui.featurePanelTab) cannot be used; colors is used.',
    'Settings: the value "-2" (ui.losatExecution.threadsPerJob) cannot be used; auto is used.',
    'Palette: the palette "Pastel" (modes.linear.ui.pendingPaletteName) waits to be applied, but Instant Preview is on,'
      + ' which applies a palette at once; it was dropped.'
  ]]);
  assert.equal(state.featurePanelTab.value, 'colors');
  assert.equal(state.losatExecution.threadsPerJob, 'auto');
  assert.equal(state.drawings.linear.pendingPaletteName.value, '');
});

test('a legacy configuration lists a setting its reader cannot restore (B48)', async () => {
  const loaded = await importSession(input({ form: { legend: 'left' }, adv: { gc_content_mode: 'ratio' } }), {
    confirmLoadRepairs: async () => true
  });
  assert.equal(loaded.status, 'legacy', JSON.stringify(loaded.error));
  assert.deepEqual(loaded.repairs.map(({ code, path, stored, used }) => [code, path, stored, used]), [
    ['setting-replaced', 'adv.gc_content_mode', 'ratio', 'deviation']
  ]);
});

// A History restore applies captured settings with the same readers and
// lists nothing (only Load asks): what it installs is unchanged (R11).
test('a History restore reads stored settings as before and asks nothing', async () => {
  const { applyConfigData } = await import('../../gbdraw/web/js/services/config.js');
  const drawing = state.drawings.circular;
  applyConfigData(drawing, {
    adv: {
      label_rendering: 'sideways', multi_record_column_gap_ratio: null, multi_record_row_gap_ratio: 'x',
      depth_min: 5, depth_max: 2, gc_content_min_percent: 60, gc_content_max_percent: 40, track_axis_gap: 'auto',
      depth_tracks: [{ label: 'D', height: -1 }]
    },
    losat: { blastp: { mode: 'x', collinearMaxConflictsInMergeGap: null, orthogroupMaxHits: 3, orthogroupHitPolicy: 'old' } },
    circularConservation: { source: 'UPLOAD', series: [null, { color: 'nonsense', losat_gencode: 0 }] },
    linearRecordLayout: { recordGap: null, rows: [{ uid: '', row: 1 }, { uid: 'a', row: 2 }] },
    losatProgram: 'BLASTP'
  }, { resolveTrackPlacements: false });
  const { adv, losat, circularConservation } = drawing;
  assert.deepEqual([adv.label_rendering, adv.multi_record_column_gap_ratio, adv.multi_record_row_gap_ratio, adv.depth_max,
    adv.gc_content_max_percent, adv.track_axis_gap], ['auto', 0, 0.05, null, 60, null]);
  assert.deepEqual(adv.depth_tracks.map(({ label, height }) => [label, height]), [['D', null]]);
  assert.deepEqual([losat.blastp.mode, losat.blastp.collinearMaxConflictsInMergeGap, losat.blastp.orthogroupMemberMaxHits,
    Object.hasOwn(losat.blastp, 'orthogroupMaxHits'), Object.hasOwn(losat.blastp, 'orthogroupHitPolicy')], ['orthogroup', 0, null, false, false]);
  assert.deepEqual([circularConservation.source, circularConservation.series.map(({ color, losat_gencode }) => [color, losat_gencode])],
    ['upload', [['#4e79a7', 1]]]);
  assert.deepEqual([drawing.linearRecordGap.value, [...drawing.linearRecordRows]], [0, [{ uid: 'a', row: 2 }]]);
  assert.equal(drawing.losatProgram.value, 'blastn');
});

// G3: every Gallery Session and every Session fixture loads without a repair;
// a repair there is a producer bug. Two fixtures are lossy on purpose. In
// node, a Session older than 40 lists its feature recovery, which needs the
// diagram Worker, and five fixtures need that Worker to load at all.
const LOSSY_FIXTURES = {
  'feature-edits-circular-copies.v44.gbdraw-session.json.gz': ['feature-visibility-narrowed'],
  'whitelist-tab-keyword.v39.gbdraw-session.json.gz': ['table-row-joined']
};
const NEEDS_DIAGRAM_WORKER = new Set([
  'feature-edits-circular.v33.gbdraw-session.json.gz', 'feature-edits-linear-crop-rc.v33.gbdraw-session.json.gz',
  'hash-color-rule-linear-rc.v30.gbdraw-session.json.gz', 'q-frame-main-web-losatn.v42.gbdraw-session.json.gz',
  'q-frame-main-web-upload.v42.gbdraw-session.json.gz'
]);
test('every Gallery Session and Session fixture loads without a repair (G3)', async () => {
  const sessions = (dir, keep) => readdirSync(new URL(dir, import.meta.url)).filter(keep).map((name) => [name, new URL(`${dir}${name}`, import.meta.url)]);
  const corpus = [
    ...sessions('../fixtures/sessions/', (name) => /\.json(\.gz)?$/.test(name) && !/provenance|expected|evidence|vectors/.test(name)),
    ...sessions('../../gbdraw/web/gallery/sessions/', (name) => /\.json(\.gz)?$/.test(name))
  ];
  assert.ok(corpus.length > 50, String(corpus.length));
  for (const [name, url] of corpus) {
    const bytes = readFileSync(url);
    const session = JSON.parse((bytes[0] === 0x1f && bytes[1] === 0x8b ? gunzipSync(bytes) : bytes).toString('utf8'));
    const loaded = await importSession(input(session), { transformLegacyResultSvg, confirmLoadRepairs: async () => true });
    if (NEEDS_DIAGRAM_WORKER.has(name)) {
      assert.equal(loaded.status, 'error', `${name} loads in node now: check it here`);
      continue;
    }
    assert.equal(loaded.status, 'ok', `${name}: ${JSON.stringify(loaded.error)}`);
    const codes = loaded.repairs.map(({ code }) => code)
      .filter((code) => !(session.version < 40 && code === 'feature-recovery-failed'));
    assert.deepEqual(codes, LOSSY_FIXTURES[name] ?? [], name);
  }
});
