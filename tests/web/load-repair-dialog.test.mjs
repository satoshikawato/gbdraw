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
// Save Session downloads through a link (G2 saves the loaded Session).
globalThis.document = { createElement: () => ({ click: () => {} }) };
installFakeSvgDom();
globalThis.File = File;
globalThis.alert = () => {};
installSessionImportWorker();

const { exportSession, importSession, normalizeEditorStateData, setUnmanagedConfigOverrideValidator } = await import('../../gbdraw/web/js/services/config.js');
const { MODE_SCOPED_SETTINGS } = await import('../../gbdraw/web/js/mode-scoped-settings.generated.js');
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
// as the mode's default or as unset and is listed with the value used
// (BA-OD89); case and outer spaces are no row.
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
    'Legend: the position "banana" (ui.layoutPreferences.linear.legend) is not one of its mode\'s positions; bottom is used.'
  ]]);
  assert.equal(state.generatedLegendPosition.value, 'bottom');
  assert.equal(state.drawings.linear.form.legend, 'bottom');
  assert.equal(state.drawings.circular.layoutPreferences.circular.single.legend, 'upper_left');
});

// G2: what Load applies after the dialog is what the dialog listed. Each
// stored leaf of a Session 46 is set to a value of another kind (and, for
// text and numbers, to an unknown value of its kind; for an unset leaf, to a
// number and a boolean too); the file is loaded and saved again. A leaf saved
// differently must be named by a row, unless the apply ignores it for a
// reason listed below, and a file Load refuses is refused before anything
// changes, so before the dialog (BA-OD96). The Session is one the current Web writes, with an entry in
// each list and two Depth series (r3 S1).
const savedSession = async (title) => {
  const result = await exportSession(title, { askChoice: async () => 'download' });
  assert.equal(result?.status, 'saved', JSON.stringify(result));
  return JSON.parse(gunzipSync(Buffer.from(await result.blob.arrayBuffer())));
};
const leafPaths = (value, parts = []) => [
  ...(parts.length ? [parts] : []),
  ...(value && typeof value === 'object' ? Object.entries(value).flatMap(([key, child]) => (
    leafPaths(child, [...parts, Array.isArray(value) ? Number(key) : key]))) : [])
];
const dotted = (parts) => parts.map((part) => (typeof part === 'number' ? `[${part}]` : `.${part}`)).join('').slice(1);
const valueAt = (value, parts) => parts.reduce((current, part) => (current == null ? undefined : current[part]), value);
/** A value of another kind, an unknown value of the same kind, and a number and a boolean for an unset leaf. */
const MUTATIONS = {
  otherKind: (value) => (value === null ? 'x-mut' : Array.isArray(value) || typeof value === 'object' ? 5
    : typeof value === 'string' ? 5 : 'x-mut'),
  sameKind: (value) => (typeof value === 'string' ? 'x-mut' : typeof value === 'number' ? -7.5 : undefined),
  unsetNumber: (value) => (value === null ? -7.5 : undefined),
  unsetBoolean: (value) => (value === null ? true : undefined)
};
/** The stored values a Load sets aside without a row, with the reason. */
const UNLISTED_LEAVES = [
  // View state (BA-OD84 (3)).
  [/^ui\.(zoom|canvasPan|selectedResultIndex)$|^modes\.\w+\.ui\.selectedFeatureRecordIdx$/, 'view state'],
  // The committed request's mode is the shown mode.
  [/^ui\.mode$/, 'the request\'s mode wins'],
  // The pending palette's colors apply only with a pending palette.
  [/^modes\.\w+\.ui\.pendingPaletteColors$/, 'no pending palette'],
  // The committed mode's request sets a layout slot that is no object.
  [/^modes\.circular\.ui\.layoutPreferences$/, 'the request\'s layout wins'],
  // The track slot model's normalizers do not report what they replace yet
  // (residual: the slot readers return what they lose, PR-E); the stack's
  // axis index is read with its slots.
  [/^modes\.\w+\.config\.adv\.(circular|linear)_track_slots(\[\d+\]\.|_axis_index$)/, 'track slot leaves'],
  // A record row that names no record of the Session names nothing the
  // figure draws; it is dropped without a row (as B47, BA-OD88).
  [/^modes\.\w+\.config\.linearRecordLayout\.rows\[\d+\]\.uid$/, 'a row naming no record']
];
const MODE_ROWS = MODE_SCOPED_SETTINGS.rows;
/** A value the drawing of its slice never draws (the registry gives it to the other mode). */
const otherModeValue = (path) => {
  const match = path.match(/^modes\.(circular|linear)\.(config|ui)\.(\w+)(?:\.(\w+))?/);
  if (!match) return false;
  const [, mode, part, head, next] = match;
  const row = MODE_ROWS.find((entry) => entry.domain === `${part}.${head}` && entry.path === next)
    ?? MODE_ROWS.find((entry) => entry.domain === part && entry.path === head);
  return Boolean(row) && row.modes !== 'both' && row.modes !== mode;
};
// The editor state is read before the dialog by its own readers, which do not
// report every value they replace yet (PR-E): a leaf saved as that reading
// gives it is theirs; a leaf the apply changes further is not.
const readBeforeTheDialog = (session, parts, saved) => {
  const [, mode, part, ...rest] = parts;
  if (parts[0] !== 'modes' || part !== 'editorState') return false;
  const read = normalizeEditorStateData(session.modes[mode].editorState, { sessionVersion: session.version });
  return JSON.stringify(valueAt(read, rest)) === JSON.stringify(saved);
};
const named = (path, repairs) => repairs.some((repair) => repair.path === path || path.startsWith(`${repair.path}.`)
  || path.startsWith(`${repair.path}[`) || repair.path.startsWith(`${path}.`) || repair.path.startsWith(`${path}[`));
const G2_FIXTURE = 'legend-renamed-rows.v46.gbdraw-session.json.gz';
/** A Depth series' fields with nothing set. */
const DEFAULT_SERIES = { height: null, large_tick_interval: null, small_tick_interval: null, tick_font_size: null };
/** The G2 fixture as the current Web saves it with an entry in each list and two Depth series. */
const filledSession = async () => {
  const draft = fixture(G2_FIXTURE);
  for (const mode of ['circular', 'linear']) {
    const { config } = draft.modes[mode];
    config.rules = [{ feat: 'CDS', qual: 'product', val: 'hypothetical', color: '#778899', cap: 'Hypothetical' }];
    config.whitelist = [{ feat: 'CDS', qual: 'product', key: 'cytochrome c oxidase' }];
    config.qualifierPriorityRules = [{ feat: 'CDS', order: 'gene' }];
    config.annotationSets = [{
      id: 'a', annotations: [{ id: 'r', target: { kind: 'coordinateSpan', start: 1, end: 10 }, label: 'x', mark: 'bracket', lane: 0 }]
    }];
    config.circularConservation.series = [{ sourceIndex: 0, label: 'Danio' }];
    config.adv.depth_tracks = [config.adv.depth_tracks[0], { ...DEFAULT_SERIES, label: 'Coverage B', color: '#123456' }];
  }
  const loaded = await importSession(input(draft), { confirmLoadRepairs: async ({ items }) => assert.fail(items.join('\n')) });
  assert.equal(loaded.status, 'ok', JSON.stringify(loaded.error));
  for (const drawing of Object.values(state.drawings)) {
    drawing.linearRecordRows.splice(0, drawing.linearRecordRows.length, ...state.linearSeqs.map(({ uid }, index) => ({ uid, row: index + 1 })));
  }
  return savedSession('g2-filled');
};
/**
 * Loads `base` with each mutation of each leaf under `roots` and saves it
 * again. Returns the files refused after the commit and the leaves saved
 * differently without a row. `savedParts` maps a leaf of the file to its place
 * in the saved Session; a row names the leaf by its place in the file.
 */
const mutationFindings = async (base, roots, savedParts = (parts) => parts) => {
  const loadedBase = await importSession(input(base), { confirmLoadRepairs: async ({ items }) => assert.fail(items.join('\n')) });
  assert.ok(['ok', 'legacy'].includes(loadedBase.status), JSON.stringify(loadedBase.error));
  assert.ok(Object.isFrozen(loadedBase.repairs));
  const baseSaved = await savedSession('g2-base');
  const leaves = roots.flatMap((root) => leafPaths(valueAt(base, root)).map((rest) => [...root, ...rest]))
    // A leaf the Save writes differently without a change is not this check's.
    .filter((parts) => JSON.stringify(valueAt(baseSaved, savedParts(parts))) === JSON.stringify(valueAt(base, parts)));
  const silent = [];
  const refusedAfterTheCommit = [];
  let checked = 0;
  for (const [kind, mutate] of Object.entries(MUTATIONS)) {
    for (const parts of leaves) {
      const path = dotted(parts);
      const savedPath = dotted(savedParts(parts));
      const mutated = mutate(valueAt(base, parts));
      if (mutated === undefined || UNLISTED_LEAVES.some(([pattern]) => pattern.test(savedPath)) || otherModeValue(savedPath)) continue;
      const file = structuredClone(base);
      valueAt(file, parts.slice(0, -1))[parts.at(-1)] = mutated;
      let committed = false;
      const loaded = await importSession(input(file), {
        confirmLoadRepairs: async () => true,
        rollbackState: { capture: () => { committed = true; }, restore: () => {} }
      });
      if (!['ok', 'legacy'].includes(loaded.status)) {
        if (committed) refusedAfterTheCommit.push(`${kind} ${path}: ${JSON.stringify(mutated)}`);
        continue;
      }
      checked += 1;
      let savedFile;
      try {
        savedFile = await savedSession(`g2-${checked}`);
      } catch (error) {
        silent.push(`${kind} ${path}: ${JSON.stringify(mutated)} -> Save fails: ${error.message}`);
        continue;
      }
      const saved = valueAt(savedFile, savedParts(parts));
      const kept = JSON.stringify(saved) === JSON.stringify(mutated) || (typeof mutated === 'number' && saved === String(mutated));
      if (!kept && !named(path, loaded.repairs) && !readBeforeTheDialog(file, savedParts(parts), saved)) {
        silent.push(`${kind} ${path}: ${JSON.stringify(mutated)} -> ${JSON.stringify(saved)}`);
      }
    }
  }
  return { checked, refusedAfterTheCommit, silent };
};
test('no apply step after the dialog drops or changes a stored value the dialog did not list, apart from the track-slot leaves (G2)', async () => {
  const base = await filledSession();
  for (const mode of ['circular', 'linear']) {
    const { config } = base.modes[mode];
    for (const list of [config.rules, config.whitelist, config.annotationSets, config.qualifierPriorityRules,
      config.circularConservation.series, config.linearRecordLayout.rows]) assert.ok(list.length > 0, mode);
    assert.equal(config.adv.depth_tracks.length, 2, mode);
  }
  const roots = [['ui'], ...['circular', 'linear'].flatMap((mode) => ['config', 'ui', 'features', 'editorState'].map((part) => ['modes', mode, part]))];
  const { checked, refusedAfterTheCommit, silent } = await mutationFindings(base, roots);
  assert.ok(checked > 100, String(checked));
  assert.deepEqual({ refusedAfterTheCommit, silent }, { refusedAfterTheCommit: [], silent: [] });
});

// G2 for a configuration file (BA-OD99 b): the stored settings of the shown
// drawing of the same Session, as a standalone file, applied over that drawing.
test('no apply step after the dialog drops or changes a stored value of a configuration file the dialog did not list (G2)', async () => {
  const session = await filledSession();
  const loadedSession = await importSession(input(session), { confirmLoadRepairs: async () => assert.fail('the fixture loads as is') });
  assert.equal(loadedSession.status, 'ok', JSON.stringify(loadedSession.error));
  const mode = state.mode.value;
  const configFile = structuredClone(session.modes[mode].config);
  for (const list of [configFile.rules, configFile.whitelist]) assert.ok(list.length > 0, mode);
  const { checked, refusedAfterTheCommit, silent } = await mutationFindings(configFile, [[]], (parts) => ['modes', mode, 'config', ...parts]);
  assert.ok(checked > 100, String(checked));
  assert.deepEqual({ refusedAfterTheCommit, silent }, { refusedAfterTheCommit: [], silent: [] });
});

// r1's silent cases: a qualifier priority rule, a form text of another kind,
// and an annotation mark and lane are read after the dialog with the readers
// that list them before it (BA-OD87).
test('a priority rule, a form value of another kind and an annotation mark or lane are listed before the Load (G2, BA-OD87)', async () => {
  const session = fixture('legend-renamed-rows.v46.gbdraw-session.json.gz');
  session.modes.circular.config.qualifierPriorityRules = [
    { feat: 'CDS', order: '' }, { feat: 'CDS', order: 'product' }, { feat: 'CDS', order: 'gene' }
  ];
  session.modes.circular.config.form.plot_title = 5;
  session.modes.linear.config.annotationSets = [{
    id: 'a', annotations: [{ id: 'r', target: { kind: 'coordinateSpan', start: 1, end: 10 }, label: 'x', mark: 'circle', lane: 'top' }]
  }];
  const { items } = await loadAnswering(session);
  assert.deepEqual(items, [
    'Settings: the value of modes.circular.config.form.plot_title cannot be used; its default is used.',
    'Settings: the value "circle" (modes.linear.config.annotationSets[0].annotations[0].mark) cannot be used; bracket is used.',
    'Settings: the value "top" (modes.linear.config.annotationSets[0].annotations[0].lane) cannot be used; 0 is used.',
    'Settings: 2 entries of modes.circular.config.qualifierPriorityRules could not be read and were dropped.'
  ]);
  assert.deepEqual([...state.drawings.circular.manualPriorityRules], [{ feat: 'CDS', order: 'gene' }]);
  assert.deepEqual(state.drawings.linear.annotationSets[0].annotations.map(({ mark, lane }) => [mark, lane]), [['bracket', 0]]);
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

// BA-OD89: a row of a list names the leaf with that leaf's values; free text
// (a label, an ID) is named by its place only; a row names the value used,
// and a valid value from another source is no row. BA-OD91 (rule 2): a
// value read as another kind is listed unless it reads as its own text.
test('a row names the leaf and the value used, never free text, and a value another source replaces is no row (BA-OD89, BA-OD91)', async () => {
  const session = fixture('legend-renamed-rows.v46.gbdraw-session.json.gz');
  session.modes.circular.config.adv.depth_tracks = [{ label: 'patient_0042_tumor.depth.tsv', height: -1 }];
  session.modes.circular.config.adv.multi_record_row_gap_ratio = true;
  session.modes.linear.config.annotationSets = [{ id: 'a', annotations: [
    { id: 'r', target: { kind: 'coordinateSpan', start: 1, end: 10 } }, { id: 'r', target: { kind: 'coordinateSpan', start: 2, end: 9 } }
  ] }];
  session.modes.linear.config.linearRecordLayout = { enabled: 'true', recordGap: 24, rows: [
    { uid: 'a', row: 1 }, { uid: 'b', row: 3, canonicalRow: 4, canonicalColumn: 1 }, { uid: 'c', row: 2, canonicalRow: '2', canonicalColumn: 1 }
  ] };
  session.ui = { ...session.ui, generatedCircularPlotTitlePosition: 5, circularPlotTitlePosition: 'top' };
  const { items } = await loadAnswering(session);
  assert.deepEqual(items, [
    'Settings: the value "true" (modes.circular.config.adv.multi_record_row_gap_ratio) cannot be used; 0.05 is used.',
    'Settings: the value "-1" (modes.circular.config.adv.depth_tracks[0].height) cannot be used; its default is used.',
    'Settings: the value "4" (modes.linear.config.linearRecordLayout.rows[1].canonicalRow) cannot be used; it is left unset.',
    'Settings: the value "1" (modes.linear.config.linearRecordLayout.rows[1].canonicalColumn) cannot be used; it is left unset.',
    'Settings: the value of modes.linear.config.annotationSets[0].annotations[1].id cannot be used; it is changed.'
  ]);
  assert.equal(state.generatedCircularPlotTitlePosition.value, 'top');
  assert.equal(state.drawings.linear.linearRecordLayoutEnabled.value, true);
  // A list given as another kind is one row of its own.
  const other = fixture('legend-renamed-rows.v46.gbdraw-session.json.gz');
  other.modes.linear.config.linearRecordLayout = { enabled: true, recordGap: 24, rows: { 0: { uid: 'a', row: 1 } } };
  const listed = await loadAnswering(other);
  assert.deepEqual(listed.items, [
    'Settings: the value "{"0":{"uid":"a","row":1}}" (modes.linear.config.linearRecordLayout.rows) cannot be used; its default is used.'
  ]);
});

// B36, B43: the app settings of `ui` are listed. B44: a pending palette that
// Instant Preview drops had no effect on the figure: no row (BA-OD88).
test('a Session 46 app setting Load cannot restore is listed; a pending palette Instant Preview drops is not (B36, B43, B44)', async () => {
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
    'Settings: the value "-2" (ui.losatExecution.threadsPerJob) cannot be used; auto is used.'
  ]]);
  assert.deepEqual(loaded.repairs.filter(({ code }) => code === 'pending-palette-dropped').map(({ path }) => path),
    ['modes.linear.ui.pendingPaletteName']);
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

// BA-OD96, r3 S5: the configuration fields Load refuses (a retired field
// name, a value no current option names) are refused before the dialog, for
// a Session slice and for a legacy configuration that is lossy otherwise, and
// nothing changes. r3's case: a Session 39's configuration as a legacy file
// (its track stacks are current, so its retired names are not renamed).
const REFUSED_CONFIG_FIELDS = [
  { adv: { linear_show_accession: true } }, { adv: { depth_tick_interval: 5 } },
  { adv: { depth_tracks: [{ label: 'D', tick_interval: 5 }] } }, { losat: { blastp: { collinearMaxGeneGap: 2 } } },
  { form: { linear_track_layout: 'sideways' } }, { adv: { label_placement: 'inside' } }, { adv: { multi_record_size_mode: 'huge' } },
  { adv: { feature_shapes: { CDS: 'bogus' } } }, { adv: { feature_shapes: { gene: 'bogus' } } },
  { adv: { arrow_head_length_ratio: -1 } }, { adv: { arrow_shaft_width_ratio: 2 } },
  // r4 S1: the stored Linear label visibilities, which the merge copies.
  { adv: { linear_accession_visibility: 'bogus' } }, { adv: { linear_accession_visibility: null } },
  { adv: { linear_length_visibility: 5 } }
];
test('a configuration field Load refuses is refused before the dialog (BA-OD96)', async () => {
  const legacyOfV39 = fixture('BGC0000708-BGC0000713.v39.gbdraw-session.json.gz').config;
  delete legacyOfV39.unmanagedConfigOverrides;
  // Lossy otherwise: an unknown GC content mode is a dialog row.
  legacyOfV39.adv.gc_content_mode = 'ratio';
  const sessionWith = ({ form = {}, adv = {}, losat }) => {
    const session = fixture('legend-renamed-rows.v46.gbdraw-session.json.gz');
    const { config } = session.modes.circular;
    Object.assign(config.form, form);
    Object.assign(config.adv, adv);
    if (losat) config.losat = { ...config.losat, blastp: { ...config.losat?.blastp, ...losat.blastp } };
    return session;
  };
  // r4 S1: a configuration file with current track stacks, which no older
  // reading renames or checks first (a Session 46 drawing's settings).
  const currentLegacyWith = ({ form = {}, adv = {}, losat }) => {
    const config = structuredClone(fixture('legend-renamed-rows.v46.gbdraw-session.json.gz').modes.circular.config);
    delete config.unmanagedConfigOverrides;
    Object.assign(config.form, form);
    Object.assign(config.adv, { gc_content_mode: 'ratio' }, adv);
    if (losat) config.losat = { ...config.losat, blastp: { ...config.losat?.blastp, ...losat.blastp } };
    return config;
  };
  const files = [['r3 legacy v39', legacyOfV39],
    ...REFUSED_CONFIG_FIELDS.flatMap((fields) => [
      [`legacy ${JSON.stringify(fields)}`, { ...fields, adv: { gc_content_mode: 'ratio', ...fields.adv } }],
      [`current legacy ${JSON.stringify(fields)}`, currentLegacyWith(fields)],
      [`session ${JSON.stringify(fields)}`, sessionWith(fields)]
    ])];
  const outcomes = [];
  for (const [name, file] of files) {
    const snapshot = loadedState();
    const asked = [];
    let committed = false;
    const loaded = await importSession(input(file), {
      confirmLoadRepairs: async ({ items }) => { asked.push(items); return true; },
      rollbackState: { capture: () => { committed = true; }, restore: () => {} }
    });
    outcomes.push([name, loaded.status, asked.length, committed]);
    // A refusal after the commit is in `outcomes`; node's rollback is not exact.
    if (loaded.status === 'error' && !committed) assert.equal(loadedState(), snapshot, name);
  }
  // A legacy file whose track stacks are older than the current schema has
  // its retired field names renamed first, so it loads and asks for its GC
  // content mode.
  const retiredNames = REFUSED_CONFIG_FIELDS.slice(0, 4).map((fields) => JSON.stringify(fields));
  assert.deepEqual(outcomes, files.map(([name]) => (retiredNames.some((fields) => name === `legacy ${fields}`)
    ? [name, 'legacy', 1, true] : [name, 'error', 0, false])));
});

// r4 S1, BA-OD99 (b): in a configuration file, a color rule or a whitelist
// entry that is no object is dropped and listed before the dialog, and the
// apply reads the lists with the same readers, so nothing is refused after
// [Load]. A Session refuses such an entry before the dialog, as before.
test('a configuration file drops and lists a rule or whitelist entry that is no object; a Session refuses it (r4 S1)', async () => {
  const rule = { feat: 'CDS', qual: 'product', val: 'kinase', color: 'red', cap: 'Kinase' };
  const file = { form: { plot_title: 'Rules' }, rules: [null, rule, 5], whitelist: [null, { feat: 'CDS', qual: 'gene', key: 'dnaA' }] };
  let committed = false;
  const asked = [];
  const loaded = await importSession(input(file), {
    confirmLoadRepairs: async ({ items }) => { asked.push(items); return true; },
    rollbackState: { capture: () => { committed = true; }, restore: () => {} }
  });
  assert.equal(loaded.status, 'legacy', JSON.stringify(loaded.error));
  assert.ok(committed);
  assert.deepEqual(asked, [[
    'Settings: 2 entries of rules could not be read and were dropped.',
    'Settings: 1 entry of whitelist could not be read and was dropped.'
  ]]);
  const drawing = state.drawings[state.mode.value];
  assert.deepEqual(JSON.parse(JSON.stringify(drawing.manualSpecificRules)), [{ ...rule, color: '#FF0000' }]);
  assert.deepEqual(JSON.parse(JSON.stringify(drawing.manualWhitelist)), [{ feat: 'CDS', qual: 'gene', key: 'dnaA' }]);
  for (const list of ['rules', 'whitelist']) {
    const session = fixture('legend-renamed-rows.v46.gbdraw-session.json.gz');
    session.modes.circular.config[list] = [null];
    const snapshot = loadedState();
    const dialogs = [];
    const refused = await importSession(input(session), { confirmLoadRepairs: async ({ items }) => { dialogs.push(items); return true; } });
    assert.deepEqual([refused.status, dialogs.length, loadedState() === snapshot], ['error', 0, true], list);
  }
});

// The shape of a feature type the user added is a key the default
// shapes lack; Load kept only the default keys and dropped it without a row.
test('a feature shape of a feature type the user added survives Save and Load', async () => {
  const loaded = await importSession(input(fixture('lambda_basic_linear.v44-schema8.gbdraw-session.json.gz')), {
    confirmLoadRepairs: async () => assert.fail('a file Load restores fully asks nothing')
  });
  assert.equal(loaded.status, 'ok', JSON.stringify(loaded.error));
  const { adv } = state.drawings.linear;
  adv.features.push('gene');
  adv.feature_shapes.gene = 'arrow';
  const saved = await savedSession('shapes');
  const again = await importSession(input(saved), {
    confirmLoadRepairs: async () => assert.fail('a file Load restores fully asks nothing')
  });
  assert.equal(again.status, 'ok', JSON.stringify(again.error));
  assert.equal(state.drawings.linear.adv.feature_shapes.gene, 'arrow');
});

// r3 N1: a similarity-group name that is no text (an object, a boolean) is
// dropped and counted, not applied as its text; a number is its text. r3 N2:
// the older `priorityRules` is read only when `qualifierPriorityRules` is
// absent, as the apply reads it.
test('a group name that is no text is dropped and listed; priorityRules is read only without qualifierPriorityRules (r3 N1, N2)', async () => {
  const { configFieldRepairs, defaultSettingsContainers, readStringMap } = await import('../../gbdraw/web/js/services/config-field-readers.js');
  const itemsOf = (config) => describeLoadRepairs(configFieldRepairs(defaultSettingsContainers('linear'), config, 'config', { mode: 'linear' }))?.items ?? [];
  const names = { og1: { x: 1 }, og2: true, og3: 5, og4: 'Kept' };
  assert.deepEqual(readStringMap(names).value, { og3: '5', og4: 'Kept' });
  assert.deepEqual(itemsOf({ webEdits: { orthogroupNameOverrides: names } }),
    ['Settings: 2 entries of config.webEdits.orthogroupNameOverrides could not be read and were dropped.']);
  const unreadable = [{ feat: '', order: 'x' }];
  assert.deepEqual(itemsOf({ qualifierPriorityRules: [{ feat: 'CDS', order: 'gene' }], priorityRules: unreadable }), []);
  assert.deepEqual(itemsOf({ priorityRules: unreadable }),
    ['Settings: 1 entry of config.priorityRules could not be read and was dropped.']);
});

// B47: a group name edit of a group the Session does not hold names nothing
// the figure draws; B46: a cache entry Load does not keep is computed again
// by the next comparison. Both are dropped without a row (Q1, BA-OD88).
test('a similarity-group edit of an unknown group and a dropped cache entry are dropped without a row (B46, B47)', async () => {
  const session = fixture('BGC0000708-BGC0000713.v39.gbdraw-session.json.gz');
  session.orthogroupState.orthogroupNameOverrides = { [session.orthogroupState.groups[0].id]: 'Kept', ghost: 'Gone' };
  // One derived entry more than the cache keeps (16): the oldest is dropped.
  const [derived] = session.losatDerivedCache.entries;
  session.losatDerivedCache.entries.push(...Array.from({ length: 16 }, (_, index) => ({ ...derived, key: `${derived.key}-${index}` })));
  const dialogs = [];
  const loaded = await importSession(input(session), {
    transformLegacyResultSvg,
    confirmLoadRepairs: async (dialog) => {
      dialogs.push(dialog.items.filter((line) => !line.startsWith('Features: the feature information')));
      return true;
    }
  });
  assert.equal(loaded.status, 'ok', JSON.stringify(loaded.error));
  assert.deepEqual(dialogs.flat(), []);
  assert.deepEqual(loaded.repairs.filter(({ code }) => code === 'cache-entry-dropped' || code === 'orthogroup-edit-dropped')
    .map(({ path, count }) => [path, count]),
  [['orthogroupState.orthogroupNameOverrides', 1], ['losatDerivedCache.entries', 1]]);
  assert.deepEqual({ ...state.drawings.linear.orthogroupNameOverrides }, { [session.orthogroupState.groups[0].id]: 'Kept' });
});

// B45, BA-OD94: a drawing's Depth series and rows fit the Depth files Load
// restores by the rule the app keeps when a Depth file goes: a series stays,
// a Depth row of a series with no file stays as the user set it, and a
// managed row goes with its file. Only a manual row that names no series of
// its drawing is turned off and listed (its own stack only).
const depthEdits = (session, mode) => {
  const { adv } = session.modes[mode].config;
  adv.depth_tracks = [adv.depth_tracks[0], { ...DEFAULT_SERIES, label: 'Coverage B', color: '#123456' }];
  const row = (id, enabled, params) => ({
    id, renderer: 'depth', enabled, width: '12', radius: null, inner_gap_px: null, outer_gap_px: null, side: null, z: 0, params
  });
  adv.circular_track_slots = [...adv.circular_track_slots, row('my_depth', true, { track_index: 1 }),
    row('lost_on', true, { track_index: 4 }), row('lost_named', false, { track_index: 5 }), row('off_unnamed', false, {})];
  return session;
};
test('a Depth row that names no series of its drawing is listed before the Load; one without a Depth file is not (B45, BA-OD94)', async () => {
  const { loaded, items } = await loadAnswering(depthEdits(fixture('legend-renamed-rows.v46.gbdraw-session.json.gz'), 'circular'));
  assert.deepEqual(items, [
    ...['lost_on', 'lost_named'].map((id) => `Track stack: the Depth row "${id}" (modes.circular.config.adv.circular_track_slots)`
      + ' names a Depth series its drawing does not have; the row is turned off and no longer names that series.')
  ]);
  assert.deepEqual(codes(loaded), ['depth-row-disabled', 'depth-row-disabled']);
  const { adv } = state.drawings.circular;
  assert.deepEqual(adv.depth_tracks.map(({ label }) => label), ['Depth', 'Coverage B']);
  assert.deepEqual(adv.circular_track_slots.filter(({ renderer }) => renderer === 'depth')
    .map(({ id, enabled, params }) => [id, enabled, params.track_index ?? null]),
  [['my_depth', true, 1], ['lost_on', false, null], ['lost_named', false, null], ['off_unnamed', false, 0]]);
  // The Linear drawing's Circular stack is not drawn.
  const linear = await loadAnswering(depthEdits(fixture('legend-renamed-rows.v46.gbdraw-session.json.gz'), 'linear'));
  assert.deepEqual(linear.items, []);
  assert.deepEqual(state.drawings.linear.adv.depth_tracks.map(({ label }) => label), ['Depth', 'Coverage B']);
});

// r1 S1, r2 S1, truth check S1 (BA-OD94): in either mode, Add Depth track
// adds a series without a Depth file, and clearing a Depth file keeps a row
// with geometry on its series (`reconcileManagedDepthSlots`). The Session
// Save writes both, and its Load asks nothing and keeps them as saved.
const DEPTH_ROUND_TRIP = {
  circular: ['legend-renamed-rows.v46.gbdraw-session.json.gz', 'circular_track_slots',
    { id: 'my_depth', renderer: 'depth', enabled: true, width: '12', radius: null, inner_gap_px: null, outer_gap_px: null, side: 'inside', z: 0 }],
  linear: ['lambda_basic_linear.v44-schema8.gbdraw-session.json.gz', 'linear_track_slots',
    { id: 'my_depth', renderer: 'depth', enabled: true, side: 'below', height: '12px', spacing: '', z: 0 }]
};
for (const [mode, [name, stackKey, manualRow]] of Object.entries(DEPTH_ROUND_TRIP)) {
  test(`a ${mode} Session keeps a Depth series without a file and a row whose file was cleared through Save and Load (BA-OD94)`, async () => {
    const { normalizeDepthTrackConfig, padDepthFileSlots, reconcileManagedDepthSlots } = await import('../../gbdraw/web/js/services/depth-track-state.js');
    await loadAnswering(fixture(name));
    assert.equal(state.mode.value, mode);
    const { adv } = state.drawings[mode];
    assert.equal(adv.depth_tracks.length, 1);
    // Add Depth track (Linear pads each record's Depth slots to the new series).
    adv.depth_tracks.push({ ...normalizeDepthTrackConfig(null, 1), label: 'Coverage B', small_tick_interval: 5 });
    if (mode === 'linear') state.linearSeqs.forEach((seq) => { seq.depth = padDepthFileSlots(seq.depth, 2); });
    // The row's series 1 had a file; clearing it keeps the row with geometry.
    const stack = adv[stackKey];
    const { slots } = reconcileManagedDepthSlots({
      slots: [...stack, { ...manualRow, params: { track_index: 1, legend_label: 'Coverage B' } }], previousSourced: [1], sourced: []
    });
    stack.splice(0, stack.length, ...slots);
    const saved = await savedSession(`depth-${mode}`);
    assert.equal(saved.modes[mode].config.adv.depth_tracks.length, 2);
    const { loaded, items } = await loadAnswering(saved);
    assert.deepEqual(items, []);
    assert.deepEqual(codes(loaded), []);
    const again = await savedSession(`depth-${mode}-again`);
    assert.deepEqual(again.modes[mode].config.adv, saved.modes[mode].config.adv);
  });
}

// A History restore applies captured settings with the same readers and
// lists nothing (only Load asks): what it installs is unchanged (R11).
test('a History restore reads stored settings as before and asks nothing', async () => {
  const { applyConfigData } = await import('../../gbdraw/web/js/services/config.js');
  const drawing = state.drawings.circular;
  applyConfigData(drawing, {
    adv: {
      label_rendering: 'sideways', multi_record_column_gap_ratio: null, multi_record_row_gap_ratio: 'x',
      depth_min: 5, depth_max: 2, gc_content_min_percent: 60, gc_content_max_percent: 40, track_axis_gap: 'auto',
      depth_tracks: [{ label: 'D', height: -1 }], depth_normalize: 'x', gc_content_show_axis: 'x', pairwise_match_style: 'Curve',
      multi_record_positions: [{ selector: '#2', row: 2 }, { selector: '#1', row: 1 }, { selector: '#2', row: 3 }]
    },
    losat: {
      executionMode: 'turbo', totalThreadBudget: 'auto',
      blastp: {
        mode: 'x', collinearMaxConflictsInMergeGap: null, orthogroupMaxHits: 3, orthogroupHitPolicy: 'old',
        collinearColorMode: 'identity', collinearAnchorMode: 'Top-N', collinearSearchScope: 'adjacent-pairs'
      }
    },
    circularConservation: { source: 'UPLOAD', series: [null, [1], { color: 'nonsense', losat_gencode: 0 }] },
    linearRecordLayout: {
      enabled: 'yes', recordGap: null,
      rows: [{ uid: '', row: 1 }, { uid: 'a', row: 2 }, { uid: 'b', row: 3, canonicalRow: '3', canonicalColumn: 1 }]
    },
    importedComparisonResolution: { action: 'bogus' },
    losatProgram: 'BLASTP'
  }, { resolveTrackPlacements: false });
  const { adv, losat, circularConservation } = drawing;
  assert.deepEqual([adv.label_rendering, adv.multi_record_column_gap_ratio, adv.multi_record_row_gap_ratio, adv.depth_max,
    adv.gc_content_max_percent, adv.track_axis_gap, adv.depth_normalize, adv.gc_content_show_axis, adv.pairwise_match_style],
  ['auto', 0, 0.05, null, 60, null, false, true, 'curve']);
  assert.deepEqual(adv.multi_record_positions, [{ selector: '#1', row: 1 }, { selector: '#2', row: 2 }]);
  assert.deepEqual(adv.depth_tracks.map(({ label, height }) => [label, height]), [['D', null]]);
  assert.deepEqual([losat.blastp.mode, losat.blastp.collinearMaxConflictsInMergeGap, losat.blastp.orthogroupMemberMaxHits,
    Object.hasOwn(losat.blastp, 'orthogroupMaxHits'), Object.hasOwn(losat.blastp, 'orthogroupHitPolicy')], ['orthogroup', 0, null, false, false]);
  assert.deepEqual([losat.blastp.collinearColorMode, losat.blastp.collinearAnchorMode, losat.blastp.collinearSearchScope],
    ['average_identity', 'all', 'adjacent']);
  assert.deepEqual([state.losatExecution.executionMode, state.losatExecution.totalThreadBudget], ['auto', 'safe']);
  assert.deepEqual([circularConservation.source, circularConservation.series.map(({ color, losat_gencode, sourceIndex }) => [color, losat_gencode, sourceIndex])],
    ['upload', [['#4e79a7', 1, 0]]]);
  assert.deepEqual([drawing.linearRecordLayoutEnabled.value, drawing.linearRecordGap.value, [...drawing.linearRecordRows]],
    [true, 0, [{ uid: 'a', row: 2 }, { uid: 'b', row: 3, canonicalRow: 3, canonicalColumn: 1 }]]);
  assert.equal(drawing.importedComparisonIntent.action, null);
  assert.equal(drawing.losatProgram.value, 'blastn');
});

// G3: every Gallery Session and every Session fixture loads without a repair;
// a repair there is a producer bug. Two fixtures are lossy on purpose. Node
// supplies no diagram Worker: a Session older than 40 reads its features again
// through it, so eight fixtures list that recovery here (and only because the
// Worker is missing: each asks for it). The browser specs load them (OV-751):
// gui-audit-20260930-inputs and live-generate-parity assert no dialog for
// seven; session-loading-feedback answers the one table row of
// whitelist-tab-keyword.v39. Five fixtures need a Worker service to load at
// all; feature-identity-overrides and gui-audit-20260930-comparisons load them.
const LOSSY_FIXTURES = {
  'feature-edits-circular-copies.v44.gbdraw-session.json.gz': ['feature-visibility-narrowed'],
  'whitelist-tab-keyword.v39.gbdraw-session.json.gz': ['table-row-joined']
};
const FEATURES_NEED_DIAGRAM_WORKER = new Set([
  'BGC0000708-BGC0000713.schema-v2.gbdraw-session.json.gz', 'BGC0000708-BGC0000713.v39.gbdraw-session.json.gz',
  'cli-linear-protein.v30.gbdraw-session.json.gz', 'conservation-fasta.v39.gbdraw-session.json.gz',
  'precomputed-ring-unlabelled-cli.v30.gbdraw-session.json.gz', 'scale-interval-negative-linear-cli.v30.gbdraw-session.json.gz',
  'scale-interval-zero-circular-cli.v30.gbdraw-session.json.gz', 'whitelist-tab-keyword.v39.gbdraw-session.json.gz'
]);
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
  const ImportWorker = globalThis.Worker;
  let askedDiagramWorker = false;
  globalThis.Worker = class extends ImportWorker {
    /** @param {string | URL} url */
    constructor(url) { askedDiagramWorker ||= !String(url).endsWith('/session-import-worker.js'); super(url); }
  };
  for (const [name, url] of corpus) {
    const bytes = readFileSync(url);
    const session = JSON.parse((bytes[0] === 0x1f && bytes[1] === 0x8b ? gunzipSync(bytes) : bytes).toString('utf8'));
    askedDiagramWorker = false;
    const loaded = await importSession(input(session), { transformLegacyResultSvg, confirmLoadRepairs: async () => true });
    if (NEEDS_DIAGRAM_WORKER.has(name)) {
      assert.equal(loaded.status, 'error', `${name} loads in node now: check it here`);
      continue;
    }
    assert.equal(loaded.status, 'ok', `${name}: ${JSON.stringify(loaded.error)}`);
    const nodeOnly = FEATURES_NEED_DIAGRAM_WORKER.has(name) ? ['feature-recovery-failed'] : [];
    assert.deepEqual(loaded.repairs.map(({ code }) => code), [...(LOSSY_FIXTURES[name] ?? []), ...nodeOnly], name);
    assert.equal(askedDiagramWorker, nodeOnly.length > 0, `${name}: the Features row comes from the missing diagram Worker`);
  }
  globalThis.Worker = ImportWorker;
});
