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
