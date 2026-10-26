// BA-D25 (Load repair dialog): Load reads a file as far as it can and, when
// something cannot be fully restored, asks once before anything changes.
// `utils/load-repairs.js` owns the wording (G1); the Load port answers.
import assert from 'node:assert/strict';
import { test } from 'node:test';
import { mkdtempSync, readFileSync, readdirSync } from 'node:fs';
import { spawnSync } from 'node:child_process';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { gunzipSync } from 'node:zlib';
import { File } from 'node:buffer';
import { installSessionImportWorker } from './helpers/session-import-node.mjs';
import { validateWebConfigOverrides } from './helpers/python-evaluator.mjs';
import { installFakeSvgDom } from './fake-svg-dom.mjs';
import { modeDrawing } from './helpers/drawing-state.mjs';
import { savedDrawingView, savedModesView, withSavedDrawingView } from './helpers/saved-drawing-view.mjs';
import { callerSession } from './helpers/session-caller.mjs';
import { sessionRollbackPort } from './helpers/session-rollback-port.mjs';

globalThis.window = {
  Vue: {
    ref: (value) => ({ value }), reactive: (value) => value,
    computed: (getter) => ({ get value() { return getter(); } }),
    nextTick: async () => {}
  },
  DOMPurify: { sanitize: (value) => value }
};
// Save Session downloads through a link (G2 saves the loaded Session); setUiLocale sets the page's lang.
globalThis.document = { createElement: () => ({ click: () => {} }), documentElement: { setAttribute: () => {} } };
installFakeSvgDom();
globalThis.File = File;
globalThis.alert = () => {};
installSessionImportWorker();

const sessionConfig = await import('../../gbdraw/web/js/services/config.js');
const { getCommittedCanonicalSession, setUnmanagedConfigOverrideValidator, installShownDrawing } = sessionConfig;
const { MODE_SCOPED_SETTINGS } = await import('../../gbdraw/web/js/mode-scoped-settings.generated.js');
const { readSessionDocument } = await import('../../gbdraw/web/js/services/session-authority.js');
// The Worker's typed override check, run by Python (as session-cli-compatibility does).
setUnmanagedConfigOverrideValidator(async (payload) => ({ result: await validateWebConfigOverrides(payload) }));
// A Session older than 40 has its Results transformed through the composition root's port.
const { transformLegacyResultSvg } = await import('../../gbdraw/web/js/app/app-setup.js');
// Node has no diagram Worker: Python's conversion of the protein results a
// Session up to 47 stored for a reversed record (OV-1145) is stubbed per test.
const { setStoredProteinResultFrameConverter } = await import('../../gbdraw/web/js/services/main-session-comparison-frame.js');
const { readSessionResourceText } = await import('../../gbdraw/web/js/services/session-resource-backing.js');
const { state } = await import('../../gbdraw/web/js/state.js');
const { bindDrawingInputFileTree } = await import('../../gbdraw/web/js/services/drawing/inputs.js');
const rollbackPort = await sessionRollbackPort(state);
const { exportSession, importSession } = callerSession(sessionConfig, state, rollbackPort);
// The real port, with a note when the Load captures or restores.
const spyOnPort = (onCapture, onRestore = () => {}) => ({
  capture: () => { onCapture(); return rollbackPort.capture(); },
  restore: (snapshot) => { onRestore(); return rollbackPort.restore(snapshot); }
});
const { describeLoadRepairs, loadRepair } = await import('../../gbdraw/web/js/utils/load-repairs.js');
const { createDefaultCircularTrackSlots } = await import('../../gbdraw/web/js/services/circular-track-slot-model.js');
const { createDefaultLinearTrackSlots } = await import('../../gbdraw/web/js/services/linear-track-slot-model.js');

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
  files: [...state.drawingList().flatMap((drawing) => Object.keys(drawing.files))],
  drawings: Object.fromEntries(['circular', 'linear'].map((mode) => {
    const drawing = modeDrawing(state, mode);
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
    rollbackState: spyOnPort(() => reached.push('snapshot'), () => reached.push('restore')),
    afterImport: async () => reached.push('afterImport'),
    beforePreviewMount: () => { reached.push('preview'); return null; },
    installLoadedArtifactSlots: () => reached.push('slots')
  });
  assert.equal(canceled.status, 'canceled');
  assert.deepEqual(dialogs, [{
    title: 'Load this file?',
    description: 'Some parts of this file cannot be fully restored:',
    items: ['1 feature edit from an older Session could not be matched to a feature of its saved diagram and was dropped.'],
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
  assert.equal(modeDrawing(state, 'linear').form.legend, 'bottom');
  assert.equal(modeDrawing(state, 'circular').layoutPreferences.circular.single.legend, 'upper_left');
});

// G2: what Load applies after the dialog is what the dialog listed. Each
// stored leaf of a Session 46 is set to a value of another kind (and, for
// text and numbers, to an unknown value of its kind; for an unset leaf, to a
// number and a boolean too); the file is loaded and saved again. A leaf saved
// differently must be named by a row, unless the apply ignores it for a
// reason listed below, and a file Load refuses is refused before anything
// changes, so before the dialog (BA-OD96). The Session is one the current Web writes, with an entry in
// each list and two Depth series.
const savedSession = async (title) => {
  const result = await exportSession(title, { askChoice: async () => 'download' });
  assert.equal(result?.status, 'saved', JSON.stringify(result));
  return JSON.parse(gunzipSync(Buffer.from(await result.blob.arrayBuffer())));
};
// The saved Session 49 read as one 46-48-shaped file of every drawing: each
// drawing's `modes.<mode>`, and the input bindings of every drawing in the
// one `webFiles.bindings` a Session 48 holds (each drawing saves its own
// mode's input slots; the shown drawing's card files win). Without the Linear
// drawing's record rows the file's layout rows would name no record.
const savedView = async (title) => {
  const document = await savedSession(title);
  const view = savedModesView(document);
  const bindingsOf = (entry) => savedDrawingView(document, entry.id).webFiles.bindings;
  const others = document.drawings.filter((entry) => entry.id !== document.ui.activeDrawingId);
  view.webFiles = { ...view.webFiles, bindings: Object.assign({}, ...others.map(bindingsOf), view.webFiles.bindings) };
  return view;
};
const leafPaths = (value, parts = []) => [
  ...(parts.length ? [parts] : []),
  ...(value && typeof value === 'object' ? Object.entries(value).flatMap(([key, child]) => (
    leafPaths(child, [...parts, Array.isArray(value) ? Number(key) : key]))) : [])
];
const dotted = (parts) => parts.map((part) => (typeof part === 'number' ? `[${part}]` : `.${part}`)).join('').slice(1);
const valueAt = (value, parts) => parts.reduce((current, part) => (current == null ? undefined : current[part]), value);
/**
 * A value of another kind, an unknown value of the same kind, a number and a
 * boolean for an unset leaf, and a valid number other than the stored one (a
 * fraction stays a fraction, an integer an integer).
 */
const MUTATIONS = {
  otherKind: (value) => (value === null ? 'x-mut' : Array.isArray(value) || typeof value === 'object' ? 5
    : typeof value === 'string' ? 5 : 'x-mut'),
  sameKind: (value) => (typeof value === 'string' ? 'x-mut' : typeof value === 'number' ? -7.5 : undefined),
  unsetNumber: (value) => (value === null ? -7.5 : undefined),
  unsetBoolean: (value) => (value === null ? true : undefined),
  otherNumber: (value) => (typeof value !== 'number' ? undefined : value > 0 && value < 1 ? value / 2 : value * 2 + 1)
};
/** The stored values a Load sets aside without a row, with the reason. */
const UNLISTED_LEAVES = [
  // View state (BA-OD84 (3)).
  [/^ui\.(zoom|canvasPan|selectedResultIndex)$|^modes\.\w+\.ui\.selectedFeatureRecordIdx$/, 'view state'],
  // The committed request's mode is the shown mode, and its input type wins
  // over a stored one of no type (D1b OV-594; the fixture's request is Circular).
  [/^ui\.(mode|cInputType)$/, 'the request\'s mode and input type win'],
  // The pending palette's colors apply only with a pending palette.
  [/^modes\.\w+\.ui\.pendingPaletteColors$/, 'no pending palette'],
  // The committed mode's request sets a layout slot that is no object.
  [/^modes\.circular\.ui\.layoutPreferences$/, 'the request\'s layout wins'],
  // The Stroke options disclosure an earlier Session saved with a Legend entry (OV-157).
  [/^modes\.\w+\.editorState\.legend\.\w+\[\d+\]\.showStroke$/, 'view state'],
  // A record row that names no record of the Session names nothing the
  // figure draws; it is dropped without a row (as B47, BA-OD88). So is a
  // record turned OFF that names no record of its drawing.
  [/^modes\.\w+\.config\.linearRecordLayout\.rows\[\d+\]\.uid$/, 'a row naming no record'],
  [/^(modes\.\w+\.config\.)?recordsOff\[\d+\]$/, 'a record key naming no record']
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
// A row names a key that holds a dot as `["key"]`.
const plainPath = (path) => path.replace(/\["([^"]*)"\]/g, '.$1');
const named = (path, repairs) => repairs.map((repair) => plainPath(repair.path)).some((at) => at === path || path.startsWith(`${at}.`)
  || path.startsWith(`${at}[`) || at.startsWith(`${path}.`) || at.startsWith(`${path}[`));
const G2_FIXTURE = 'legend-renamed-rows.v46.gbdraw-session.json.gz';
/** A Depth series' fields with nothing set. */
const DEFAULT_SERIES = { height: null, large_tick_interval: null, small_tick_interval: null, tick_font_size: null };
/**
 * The G2 fixture as the current Web saves it with an entry in each list and two Depth series. `overrides`: with an
 * unmanaged config override, whose check runs Python once per Load (about 2 s), so only its own leaves are mutated.
 */
const filledSession = async ({ overrides = false } = {}) => {
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
    const { features } = draft.modes[mode];
    features.featureColorOverrides = { 'feature-a': '#112233' };
    features.labelOverrideRows = [{ recordId: '*', featureType: 'CDS', qualifier: 'product', valueRegex: 'x', labelText: 'y' }];
  }
  // Each list the E1 readers read holds an entry, so G2 changes it. (A
  // comparison edge names two records; the fixture has one Linear record.)
  const circular = draft.modes.circular;
  const recordKey = 'circular:#1';
  circular.config.recordsOff = ['#1'];
  circular.config.featurePlacementOverrides = {
    [JSON.stringify([recordKey, 'f1234abcd'])]: { recordKey, biologicalFeatureId: 'f1234abcd', placement: { kind: 'main' } }
  };
  if (overrides) circular.config.unmanagedConfigOverrides = { 'objects.gc_content.percent_background_opacity': 0.42 };
  circular.config.recordDisplayDrafts = [{
    sourceUid: 'circular', selector: '#1', recordId: 'same', topologyOverride: null, startCoordinate: 25,
    reverseComplementOverride: null, anchorIntent: null
  }];
  draft.modes.linear.config.recordsOff = [draft.webFiles.bindings.linearSeqs[0].uid];
  draft.modes.linear.editorState.legend.strokeOverrides = { CDS: { strokeColor: '#123456', strokeWidth: 1 } };
  const loaded = await importSession(input(draft), { confirmLoadRepairs: async ({ items }) => assert.fail(items.join('\n')) });
  assert.equal(loaded.status, 'ok', JSON.stringify(loaded.error));
  for (const drawing of state.drawingList()) {
    drawing.linearRecordRows.splice(0, drawing.linearRecordRows.length, ...modeDrawing(state, 'linear').linearSeqs.map(({ uid }, index) => ({ uid, row: index + 1 })));
  }
  // The G2 file is the saved drawings as one 48-shaped file (every leaf has a place `modes.<mode>...`).
  return { ...(await savedView('g2-filled')), version: 48 };
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
  const baseSaved = await savedView('g2-base');
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
        rollbackState: spyOnPort(() => { committed = true; })
      });
      if (!['ok', 'legacy'].includes(loaded.status)) {
        if (committed) refusedAfterTheCommit.push(`${kind} ${path}: ${JSON.stringify(mutated)}`);
        continue;
      }
      checked += 1;
      let savedFile;
      try {
        savedFile = await savedView(`g2-${checked}`);
      } catch (error) {
        silent.push(`${kind} ${path}: ${JSON.stringify(mutated)} -> Save fails: ${error.message}`);
        continue;
      }
      const saved = valueAt(savedFile, savedParts(parts));
      // A number may be saved as its text, and a track slot measure as pixels.
      const kept = JSON.stringify(saved) === JSON.stringify(mutated)
        || (typeof mutated === 'number' && (saved === String(mutated) || saved === `${mutated}px`));
      if (!kept && !named(path, loaded.repairs)) {
        silent.push(`${kind} ${path}: ${JSON.stringify(mutated)} -> ${JSON.stringify(saved)}`);
      }
    }
  }
  return { checked, refusedAfterTheCommit, silent };
};
test('no apply step after the dialog drops or changes a stored value the dialog did not list (G2)', async () => {
  const base = await filledSession();
  for (const mode of ['circular', 'linear']) {
    const { config } = base.modes[mode];
    for (const list of [config.rules, config.whitelist, config.annotationSets, config.qualifierPriorityRules,
      config.circularConservation.series, config.linearRecordLayout.rows]) assert.ok(list.length > 0, mode);
    assert.equal(config.adv.depth_tracks.length, 2, mode);
  }
  const circular = base.modes.circular;
  for (const list of [circular.features.featureColorOverrides, circular.features.labelOverrideRows, circular.config.recordsOff,
    circular.config.featurePlacementOverrides, circular.config.recordDisplayDrafts,
    circular.config.adv.circular_track_slots, base.modes.linear.config.recordsOff, base.modes.linear.config.adv.linear_track_slots,
    base.modes.linear.editorState.legend.strokeOverrides]) assert.ok(Object.keys(list).length > 0);
  const roots = [['ui'], ...['circular', 'linear'].flatMap((mode) => ['config', 'ui', 'features', 'editorState'].map((part) => ['modes', mode, part]))];
  const { checked, refusedAfterTheCommit, silent } = await mutationFindings(base, roots);
  assert.ok(checked > 100, String(checked));
  assert.deepEqual({ refusedAfterTheCommit, silent }, { refusedAfterTheCommit: [], silent: [] });
  const withOverrides = await filledSession({ overrides: true });
  const overrides = await mutationFindings(withOverrides, [['modes', 'circular', 'config', 'unmanagedConfigOverrides']]);
  assert.ok(overrides.checked > 0, String(overrides.checked));
  assert.deepEqual({ refusedAfterTheCommit: overrides.refusedAfterTheCommit, silent: overrides.silent },
    { refusedAfterTheCommit: [], silent: [] });
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

// Values Load changed without a row: a qualifier priority rule, a form text of another kind,
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
  assert.deepEqual([...modeDrawing(state, 'circular').manualPriorityRules], [{ feat: 'CDS', order: 'gene' }]);
  assert.deepEqual(modeDrawing(state, 'linear').annotationSets[0].annotations.map(({ mark, lane }) => [mark, lane]), [['bracket', 0]]);
});

// BA-D15 (OV-569, OV-710): a legacy configuration file is read before it is
// applied; malformed record display settings are dropped and listed, and
// Cancel leaves the drawing as it was.
test('a legacy configuration with malformed record display settings loads with the rows dropped and listed (PLAN 12.7 (3), E-OD04)', async () => {
  const legacy = { form: { legend: 'left' }, recordDisplayDrafts: [{ junk: 1 }, 5, 'x'] };
  const drawing = modeDrawing(state, state.mode.value === 'linear' ? 'linear' : 'circular');
  drawing.recordDisplayDrafts.splice(0, drawing.recordDisplayDrafts.length);
  const formBefore = JSON.stringify(drawing.form);
  const dialogs = [];
  const canceled = await importSession(input(legacy), { confirmLoadRepairs: async (dialog) => { dialogs.push(dialog); return false; } });
  assert.equal(canceled.status, 'canceled');
  assert.equal(JSON.stringify(drawing.form), formBefore);
  assert.deepEqual(dialogs.map(({ items }) => items), [
    ['Settings: 3 entries of recordDisplayDrafts could not be read and were dropped.']
  ]);
  const loaded = await importSession(input(legacy), { confirmLoadRepairs: async () => true });
  assert.equal(loaded.status, 'legacy', JSON.stringify(loaded.error));
  assert.deepEqual(loaded.repairs.map(({ code, path }) => [code, path]), [['setting-entries-dropped', 'recordDisplayDrafts']]);
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
    'Legend: entries "X" and "W" (deleted) had a color the app does not accept and were dropped.',
    '2 feature edits from an older Session could not be matched to a feature of its saved diagram and were dropped.',
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
    'Feature visibility edits from an older Session',
    'could not be renamed to the hash of a feature',
    'had a color the app does not accept',
    'joined into the last column',
    'lacked a required column',
    'had a color the table does not accept',
    'cannot be fully restored'
  ];
  // The generated translation catalogs repeat the English as their keys; they word nothing.
  const owners = modules('').filter((path) => !path.startsWith('utils/i18n-catalogs/'))
    .filter((path) => texts.some((text) => readFileSync(new URL(path, root), 'utf8').includes(text)));
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
  assert.deepEqual(modeDrawing(state, 'linear').legendColorOverrides, { K: '#222222' });
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
  assert.deepEqual(modeDrawing(state, 'circular').legendStrokeOverrides['Complex IV'], { strokeColor: '#123456', strokeWidth: 1 });
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
  assert.deepEqual(modeDrawing(state, 'circular').legendEntries.value.map(({ caption }) => caption).filter((caption) => ['tRNA', 'Q'].includes(caption)),
    ['tRNA', 'Q']);

  const session47 = fixture('legend-renamed-rows.v46.gbdraw-session.json.gz');
  session47.version = 47;
  session47.modes.circular.editorState.legend.entries.push(
    { caption: 'tRNA', originalCaption: 'tRNA 2', color: '#111111' },
    { caption: 'Q', originalCaption: 'CDS', color: '#222222' }
  );
  const read47 = await loadAnswering(session47);
  assert.deepEqual(read47.items, ['Legend: entry "Q" repeats an earlier entry of its list and was dropped.']);
  assert.deepEqual(modeDrawing(state, 'circular').legendEntries.value.filter(({ caption }) => caption === 'tRNA').length, 2);
});

// OV-758 (A10, E-OD08): the track stacks' readers list what they replace
// before the dialog. A slot value Load refused is dropped and its default
// used; a stack with a slot that is no object, or one the slot check still
// refuses then, is read as the default stack, with its axis index (a place
// in the stack) unset. An axis index the stack clamps keeps today's value.
// A schema version Load does not know reads the slots as the current schema.
test('a track slot value Load cannot use is listed; a slot that cannot be read resets its stack (OV-758, A10)', async () => {
  const v46 = () => fixture('legend-renamed-rows.v46.gbdraw-session.json.gz');
  let session = v46();
  session.modes.circular.config.adv.circular_track_slots_axis_index = 9;
  let read = await loadAnswering(session);
  assert.deepEqual(read.items, ['Settings: the value "9" (modes.circular.config.adv.circular_track_slots_axis_index) cannot be used; 4 is used.']);
  assert.equal(modeDrawing(state, 'circular').adv.circular_track_slots_axis_index, 4);

  session = v46();
  const slots = session.modes.circular.config.adv.circular_track_slots;
  Object.assign(slots[1], { width: 'wide', z: 'top' });
  slots[2].side = 'zigzag';
  slots[3].params = 5;
  session.modes.linear.config.adv.linear_track_slots[0].height = -7.5;
  read = await loadAnswering(session);
  assert.deepEqual(read.items, [
    'Settings: the value "wide" (modes.circular.config.adv.circular_track_slots[1].width) cannot be used; it is left unset.',
    'Settings: the value "top" (modes.circular.config.adv.circular_track_slots[1].z) cannot be used; it is left unset.',
    'Settings: the value "zigzag" (modes.circular.config.adv.circular_track_slots[2].side) cannot be used; it is left unset.',
    'Settings: the value "5" (modes.circular.config.adv.circular_track_slots[3].params) cannot be used; its default is used.',
    'Settings: the value "-7.5" (modes.linear.config.adv.linear_track_slots[0].height) cannot be used; it is left unset.'
  ]);
  const circularSlots = modeDrawing(state, 'circular').adv.circular_track_slots;
  assert.deepEqual([circularSlots[1].width, circularSlots[1].z, circularSlots[3].params], [null, 0, {}]);
  assert.ok([undefined, null, ''].includes(modeDrawing(state, 'linear').adv.linear_track_slots[0].height));

  session = v46();
  session.modes.circular.config.adv.circular_track_slots[1] = 5;
  session.modes.circular.config.adv.circular_track_slots_axis_index = 2;
  read = await loadAnswering(session);
  assert.deepEqual(read.items, [
    'Settings: the value "2" (modes.circular.config.adv.circular_track_slots_axis_index) cannot be used; its default is used.',
    'Settings: the value of modes.circular.config.adv.circular_track_slots cannot be used; its default is used.'
  ]);
  const ids = (stack) => stack.map(({ id }) => id);
  assert.deepEqual(ids(modeDrawing(state, 'circular').adv.circular_track_slots), ids(createDefaultCircularTrackSlots()));

  for (const stack of ['x', [{ id: 'x', renderer: 'zzz' }]]) {
    session = fixture('lambda_basic_linear.v44-schema8.gbdraw-session.json.gz');
    session.config.adv.linear_track_slots = stack;
    read = await loadAnswering(session);
    assert.deepEqual(read.items, ['Settings: the value of config.adv.linear_track_slots cannot be used; its default is used.'], JSON.stringify(stack));
    assert.deepEqual(ids(modeDrawing(state, 'linear').adv.linear_track_slots), ids(createDefaultLinearTrackSlots()));
  }

  session = fixture('HmmtDNA_basic_circular.v44-schema8.gbdraw-session.json.gz');
  session.config.adv.circular_track_slots_schema_version = 1;
  const quiet = await importSession(input(session), { confirmLoadRepairs: async ({ items }) => assert.fail(items.join('\n')) });
  assert.equal(quiet.status, 'ok', JSON.stringify(quiet.error));
});

// TB8, TB9 (E-OD08): the editor state's readers list the Legend entries,
// color and stroke edits, orders and original stroke they drop. An edit with
// no value edits nothing, and `showStroke` is view state: no row.
test('a Legend entry, color or stroke edit, order or original stroke Load cannot read is listed (TB8, TB9)', async () => {
  const session = fixture('lambda_basic_linear.v44-schema8.gbdraw-session.json.gz');
  session.editorState = {
    ...session.editorState,
    legend: {
      entries: [{ caption: '', color: '#123456' }, 5, { caption: 'A', color: '#111111', showStroke: true }],
      colorOverrides: { X: 'notacolor', Y: '' },
      strokeOverrides: { Y: { strokeColor: 'bad', strokeWidth: -1 }, Z: { strokeColor: '#00aa00', strokeWidth: 'w' }, E: {} },
      originalOrder: [5, '', 'A', false]
    },
    featureStrokes: { overrides: { f1: { strokeColor: 'bad' } } },
    originalSvgStroke: { color: 'bad', width: -2 }
  };
  const { items } = await loadAnswering(session);
  assert.deepEqual(items, [
    'Settings: the value "w" (editorState.legend.strokeOverrides.Z.strokeWidth) cannot be used; it is left unset.',
    'Settings: the value "bad" (editorState.originalSvgStroke.color) cannot be used; its default is used.',
    'Settings: the value "-2" (editorState.originalSvgStroke.width) cannot be used; its default is used.',
    'Settings: 1 entry of editorState.legend.strokeOverrides could not be read and was dropped.',
    'Settings: 2 entries of editorState.legend.entries could not be read and were dropped.',
    'Settings: 1 entry of editorState.legend.originalOrder could not be read and was dropped.',
    'Settings: 1 entry of editorState.featureStrokes.overrides could not be read and was dropped.',
    'Legend: entry "X" had a color the app does not accept and was dropped.'
  ]);
  assert.deepEqual(modeDrawing(state, 'linear').legendStrokeOverrides, { Z: { strokeColor: '#00aa00' } });
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
  assert.match(state.featureEditorStatus.value.error, /does not include embedded GenBank inputs/);

  // OV-1764: the failure is the shown drawing's. The other drawing, with a Result, shows its own
  // count; switching back shows the failure again (the status switches with the artifact slot).
  assert.equal(state.featureEditorStatusText.value, 'Feature editor unavailable');
  const { createHistorySnapshotService } = await import('../../gbdraw/web/js/services/history-snapshot.js');
  const { createInputFileStore } = await import('../../gbdraw/web/js/services/history-files.js');
  const { createArtifactSlot } = await import('../../gbdraw/web/js/services/artifact-slot.js');
  const slots = createHistorySnapshotService({ state, fileStore: createInputFileStore(), currentResults: () => state.results.value });
  const failedSlot = slots.captureArtifactSlot();
  const circular = modeDrawing(state, 'circular');
  const rows = [{ svg_id: 'f1' }, { svg_id: 'f2' }, { svg_id: 'f3' }];
  slots.installArtifactSlot(createArtifactSlot({
    drawing: circular, values: { results: [{ name: 'c.svg', content: '<svg/>' }], extractedFeatures: rows }
  }), { drawing: circular });
  assert.equal(state.featureEditorStatusText.value, '3 features ready', 'the Circular drawing has its own Result');
  slots.installArtifactSlot(failedSlot, { drawing: modeDrawing(state, 'linear') });
  assert.equal(state.featureEditorStatusText.value, 'Feature editor unavailable');
  assert.match(state.featureEditorStatus.value.error, /does not include embedded GenBank inputs/);
});

// OV-1766: every Load writes the shown drawing's feature-editor status; a Load
// whose shown drawing has no features does not keep an earlier Load's failure.
test('a Load without features does not keep the previous Load\'s "Feature editor unavailable" (OV-1766)', async () => {
  const failed = fixture('hash-color-rule-linear-rc.v30.gbdraw-session.json.gz');
  failed.features.extractedFeatures = [];
  Object.keys(failed.files).forEach((key) => { failed.files[key] = null; });
  const confirmLoadRepairs = async () => true;
  assert.equal((await importSession(input(failed), { transformLegacyResultSvg, confirmLoadRepairs })).status, 'ok');
  assert.equal(state.featureEditorStatusText.value, 'Feature editor unavailable');
  const next = await importSession(input(fixture('file-card-imports-circular.v48.gbdraw-session.json.gz')), {
    transformLegacyResultSvg, confirmLoadRepairs
  });
  assert.equal(next.status, 'ok', JSON.stringify(next.error));
  assert.equal(state.extractedFeatures.value.length, 0);
  assert.deepEqual({ ...state.featureEditorStatus.value }, { status: 'idle', error: null });
  assert.equal(state.featureEditorStatusText.value, 'No features ready');
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

// OV-890: a Depth binding sits one list deeper (`c_depth[record][slot]`); one
// that names a resource the Session does not hold is dropped with the row the
// other bindings get, its slot left empty. With no resources at all the file
// still loads with that row: its settings and Results are restored (R2-m5).
test('a Depth binding that names no resource is dropped with its row (OV-890)', async () => {
  const depthId = 'depth-tracks-1-source-record-1';
  const session = fixture('two-mode-project.v44.gbdraw-session.json.gz');
  assert.deepEqual(session.webFiles.bindings.c_depth.map((row) => row.map((binding) => binding?.resourceId)), [[depthId]]);
  delete session.resources[depthId];
  const read = readSessionDocument(session);
  assert.ok(read.repairs.some(({ code, path }) => code === 'file-entry-dropped' && path === 'webFiles.bindings.c_depth[0][0]'),
    JSON.stringify(read.repairs));
  assert.deepEqual(read.value.webFiles.bindings.c_depth, [[null]]);
  const { loaded } = await loadAnswering(session);
  assert.ok(loaded.repairs.some(({ path }) => path === 'webFiles.bindings.c_depth[0][0]'), JSON.stringify(loaded.repairs));
  for (const resources of [null, [], undefined]) {
    const bare = fixture('two-mode-project.v44.gbdraw-session.json.gz');
    if (resources === undefined) delete bare.resources;
    else bare.resources = resources;
    const repairs = readSessionDocument(bare).repairs;
    assert.ok(repairs.some(({ path }) => path === 'webFiles.bindings.c_depth[0][0]'), JSON.stringify(repairs));
    const result = await importSession(input(bare), { confirmLoadRepairs: async () => true });
    assert.equal(result.status, 'ok', JSON.stringify(result.error));
    assert.ok(result.repairs.some(({ path }) => path === 'webFiles.bindings.c_depth[0][0]'), JSON.stringify(result.repairs));
  }
});

// OV-891 (E-OD24): a null definition lock reads as an omission, the
// historical OFF: the drawing holds false, and no row is listed (E-OD06).
test('a null definition lock loads as OFF with no row (OV-891)', async () => {
  const session = fixture('legend-renamed-rows.v46.gbdraw-session.json.gz');
  assert.equal(session.modes.linear.config.form.keep_definition_left_aligned, true);
  session.modes.linear.config.form.keep_definition_left_aligned = null;
  const loaded = await importSession(input(session), { confirmLoadRepairs: async () => assert.fail('no row (E-OD06)') });
  assert.equal(loaded.status, 'ok', JSON.stringify(loaded.error));
  assert.equal(modeDrawing(state, 'linear').form.keep_definition_left_aligned, false);
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
  assert.equal(modeDrawing(state, 'linear').adv.label_rendering, 'auto');
  // The circular settings are the Circular drawing's (the split, plan 4.2).
  const { adv } = modeDrawing(state, 'circular');
  assert.equal(adv.circular_label_placement, 'radial');
  assert.equal(adv.multi_record_min_radius_ratio, 0.55);
  assert.deepEqual(adv.multi_record_positions, [{ selector: '#1', row: 1 }, { selector: '#2', row: 2 }]);
  assert.equal(state.losatExecution.executionMode, 'auto');
  assert.equal(state.losatExecution.totalThreadBudget, 'safe');
});

// OV-765: a Session 27-44 keeps one draft that both drawings read, so a value
// both read is one place in the file: one record and one row. A Session 46
// keeps a draft per drawing: the same value in each is two places, two rows.
test('a value both drawings of a Session 27-44 read is one row; a Session 46 lists each drawing\'s place (OV-765)', async () => {
  const rowsOf = ({ repairs }) => repairs.map(({ code, path, stored, used }) => `${code} ${path} ${stored} -> ${used}`);
  const load = async (session, name) => {
    const dialogs = [];
    const loaded = await importSession(input(session), {
      transformLegacyResultSvg, confirmLoadRepairs: async ({ items }) => { dialogs.push(items); return true; }
    });
    assert.equal(loaded.status, 'ok', `${name}: ${JSON.stringify(loaded.error)}`);
    return { loaded, items: dialogs.flat() };
  };
  const lossy = (config) => {
    config.adv = { ...config.adv, label_rendering: 'sideways', feature_shapes: { gene: 'bogus' } };
  };
  for (const name of ['BGC0000708-BGC0000713.v30.gbdraw-session.json.gz', 'lambda_basic_linear.v44-schema8.gbdraw-session.json.gz',
    'HmmtDNA_basic_circular.v44-schema8.gbdraw-session.json.gz']) {
    const session = fixture(name);
    lossy(session.config);
    const { loaded, items } = await load(session, name);
    // The text '0.5' reads as its number (no row).
    assert.deepEqual(rowsOf(loaded), [
      'setting-replaced config.adv.label_rendering sideways -> auto',
      'setting-replaced config.adv.feature_shapes.gene bogus -> undefined'
    ], name);
    assert.equal(items.length, 2, name);
  }
  const session46 = fixture('legend-renamed-rows.v46.gbdraw-session.json.gz');
  const MODES = ['circular', 'linear'];
  MODES.forEach((mode) => lossy(session46.modes[mode].config));
  // A Linear-only setting makes the Linear part a drawing the move keeps; an
  // unused part would take nothing of the file (D-15-120).
  const linearConfig = session46.modes.linear.config;
  linearConfig.form = { ...linearConfig.form, linear_track_layout: 'above' };
  const { loaded, items } = await load(session46, 'Session 46');
  assert.deepEqual(rowsOf(loaded), MODES.flatMap((mode) => [
    `setting-replaced modes.${mode}.config.adv.label_rendering sideways -> auto`,
    `setting-replaced modes.${mode}.config.adv.feature_shapes.gene bogus -> undefined`
  ]));
  assert.equal(items.length, 4);
});

// A Session 27-33 drew its repeat regions as rectangles whatever its shape
// map held, so a shape map that is no object reads as the repeat-region
// shape with no row, as the Python twin reads it. A configuration file loads
// over the drawing: there the map is not used, and its row says the drawing
// keeps its shapes.
test('a Session 27-33 shape map that is no object reads as the repeat-region shape; a configuration file lists it', async () => {
  const session = fixture('BGC0000708-BGC0000713.v30.gbdraw-session.json.gz');
  const features = Array.isArray(session.config.adv?.features) ? session.config.adv.features : [];
  session.config.adv = { ...session.config.adv, features: [...new Set([...features, 'repeat_region'])], feature_shapes: ['arrow'] };
  const loaded = await importSession(input(session), { transformLegacyResultSvg, confirmLoadRepairs: async () => true });
  assert.equal(loaded.status, 'ok', JSON.stringify(loaded.error));
  assert.deepEqual(loaded.repairs.filter(({ path }) => path.includes('feature_shapes')), []);
  assert.equal(modeDrawing(state, state.mode.value).adv.feature_shapes.repeat_region, 'rectangle');
  const drawing = modeDrawing(state, state.mode.value);
  const shapes = JSON.stringify(drawing.adv.feature_shapes);
  const asked = [];
  const fromFile = await importSession(input({ adv: { feature_shapes: ['arrow'] } }), {
    confirmLoadRepairs: async ({ items }) => { asked.push(...items); return true; }
  });
  assert.equal(fromFile.status, 'legacy', JSON.stringify(fromFile.error));
  assert.deepEqual(asked, ['Settings: the value "[\"arrow\"]" (adv.feature_shapes) cannot be used; the drawing keeps its current value.']);
  assert.equal(JSON.stringify(modeDrawing(state, state.mode.value).adv.feature_shapes), shapes);
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
  assert.equal(modeDrawing(state, 'linear').linearRecordLayoutEnabled.value, true);
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
  assert.equal(modeDrawing(state, 'linear').pendingPaletteName.value, '');
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

// BA-OD96: the configuration fields Load refuses (a retired field
// name, a value no current option names) are refused before the dialog, for
// a Session slice and for a legacy configuration that is lossy otherwise, and
// nothing changes. One case is a Session 39's configuration as a legacy file
// (its track stacks are current, so its retired names are not renamed).
const READ_CONFIG_FIELDS = [
  { adv: { linear_show_accession: true } }, { adv: { depth_tick_interval: 5 } },
  { adv: { depth_tracks: [{ label: 'D', tick_interval: 5 }] } }, { losat: { blastp: { collinearMaxGeneGap: 2 } } },
  { form: { linear_track_layout: 'sideways' } }, { adv: { label_placement: 'inside' } }, { adv: { multi_record_size_mode: 'huge' } },
  // A default feature type's shape the merge installs (text or null) that names no shape.
  { adv: { feature_shapes: { CDS: 'bogus' } } }, { adv: { feature_shapes: { CDS: null } } },
  { adv: { arrow_head_length_ratio: -1 } }, { adv: { arrow_shaft_width_ratio: 2 } },
  // The stored Linear label visibilities the merge copies (text or null).
  { adv: { linear_accession_visibility: 'bogus' } }, { adv: { linear_accession_visibility: null } }
];
// A Session 46 with `fields` in its Circular drawing's settings.
const sessionWith = ({ form = {}, adv = {}, losat }) => {
  const session = fixture('legend-renamed-rows.v46.gbdraw-session.json.gz');
  const { config } = session.modes.circular;
  Object.assign(config.form, form);
  Object.assign(config.adv, adv);
  if (losat) config.losat = { ...config.losat, blastp: { ...config.losat?.blastp, ...losat.blastp } };
  return session;
};
// A configuration file with current track stacks, which no older
// reading renames or checks first (a Session 46 drawing's settings), made
// lossy by an unknown GC content mode.
const currentLegacyWith = ({ form = {}, adv = {}, losat }) => {
  const config = structuredClone(fixture('legend-renamed-rows.v46.gbdraw-session.json.gz').modes.circular.config);
  delete config.unmanagedConfigOverrides;
  Object.assign(config.form, form);
  Object.assign(config.adv, { gc_content_mode: 'ratio' }, adv);
  if (losat) config.losat = { ...config.losat, blastp: { ...config.losat?.blastp, ...losat.blastp } };
  return config;
};
// The same as a minimal configuration file, which older readings migrate first.
const legacyWith = (fields) => ({ ...fields, adv: { gc_content_mode: 'ratio', ...fields.adv } });
// E1a (OV-757, A17c-A17g): a configuration field Load refused (BA-OD96) is
// read in every file: a retired name as its current name, an older option
// value as its current name, and a value of no choice as the default with a
// row. A Session 46 lists only what its drawing draws (A17b'); the
// configuration files list their unknown GC content mode too.
test('a configuration field Load refused is read with its row before the dialog (E1a)', async () => {
  const legacyOfV39 = fixture('BGC0000708-BGC0000713.v39.gbdraw-session.json.gz').config;
  delete legacyOfV39.unmanagedConfigOverrides;
  // Lossy otherwise: an unknown GC content mode is a dialog row.
  legacyOfV39.adv.gc_content_mode = 'ratio';
  const files = [['r3 legacy v39', legacyOfV39],
    ...READ_CONFIG_FIELDS.flatMap((fields) => [
      [`legacy ${JSON.stringify(fields)}`, legacyWith(fields)],
      [`current legacy ${JSON.stringify(fields)}`, currentLegacyWith(fields)],
      [`session ${JSON.stringify(fields)}`, sessionWith(fields)]
    ])];
  const outcomes = [];
  for (const [name, file] of files) {
    const asked = [];
    const loaded = await importSession(input(file), { confirmLoadRepairs: async ({ items }) => { asked.push(items); return true; } });
    outcomes.push([name, loaded.status, asked.flat().length]);
    if (loaded.status === 'ok' || loaded.status === 'legacy') await savedSession(name);
  }
  // The Circular drawing of a Session 46 draws these; a retired name reads as its current name.
  const circularRows = [{ adv: { multi_record_size_mode: 'huge' } }, { adv: { feature_shapes: { CDS: 'bogus' } } },
    { adv: { arrow_head_length_ratio: -1 } }, { adv: { arrow_shaft_width_ratio: 2 } }].map((fields) => `session ${JSON.stringify(fields)}`);
  assert.deepEqual(outcomes, files.map(([name]) => (name.startsWith('session ')
    ? [name, 'ok', circularRows.includes(name) ? 1 : 0]
    // The shown (Circular) drawing reads a configuration file; it draws no Linear setting.
    : [name, 'legacy', /multi_record_size_mode|"CDS":"bogus"|arrow/.test(name) ? 2 : 1])));
});

// BA-OD100: Load refuses no file the Web loaded before the dialog existed. These files
// loaded on origin/dev 28c8277f (fix4 probe: one process per file, the value
// in a minimal and a current configuration file and in a Session 46); the
// merge did not use the value or the feature type's shape was not read.
const LOADED_BEFORE = [
  ...[{ CDS: 5 }, { CDS: true }, { gene: 'bogus' }, { gene: 5 }, { gene: null }, { gene: { a: 1 } }, { '': 'arrow' }, { ' ': 'bogus' }]
    .flatMap((feature_shapes) => {
      const fields = { adv: { feature_shapes } };
      return [[`legacy ${JSON.stringify(fields)}`, legacyWith(fields)], [`current legacy ${JSON.stringify(fields)}`, currentLegacyWith(fields)],
        [`session ${JSON.stringify(fields)}`, sessionWith(fields)]];
    }),
  // A Linear label visibility of another kind in a configuration file with current track stacks.
  ...[{ linear_accession_visibility: 5 }, { linear_accession_visibility: true }, { linear_accession_visibility: {} },
    { linear_length_visibility: 0 }, { linear_length_visibility: false }]
    .map((adv) => [`current legacy ${JSON.stringify(adv)}`, currentLegacyWith({ adv })])
];
test('no file the Web loaded before the dialog existed is refused (BA-OD100)', async () => {
  const outcomes = [];
  for (const [name, file] of LOADED_BEFORE) {
    const loaded = await importSession(input(file), { confirmLoadRepairs: async () => true });
    outcomes.push([name, loaded.status === 'ok' || loaded.status === 'legacy' ? 'loaded' : loaded.error?.message ?? loaded.status]);
    if (loaded.status === 'ok' || loaded.status === 'legacy') await savedSession(name);
  }
  assert.deepEqual(outcomes, LOADED_BEFORE.map(([name]) => [name, 'loaded']));
});

// BA-OD100: in a configuration file, a default feature type's shape that is
// no text and a Linear label visibility that is no text are not used (the
// drawing keeps its value, `kept`, named in the row), and the shape of a feature type
// the defaults lack that names no shape is dropped. The apply reads them with
// the readers the dialog lists them by.
test('a configuration file with a shape or Linear visibility that is no text, or an unknown shape of a new feature type, loads with its rows (BA-OD100)', async () => {
  // A Linear Session first: Save of a Session 49 needs the shown Linear drawing to hold a sequence.
  await loadAnswering(fixture('lambda_basic_linear.v44-schema8.gbdraw-session.json.gz'));
  const shownMode = state.mode.value;
  installShownDrawing('linear');
  try {
    const drawing = modeDrawing(state, 'linear');
    const before = { CDS: drawing.adv.feature_shapes.CDS, accession: drawing.adv.linear_accession_visibility,
      length: drawing.adv.linear_length_visibility };
    const file = currentLegacyWith({ adv: {} });
    Object.assign(file.adv, {
      feature_shapes: { ...file.adv.feature_shapes, CDS: 5, gene: 'bogus', misc_feature: ' Rectangle ' },
      linear_accession_visibility: true, linear_length_visibility: 0
    });
    const loaded = await importSession(input(file), { confirmLoadRepairs: async () => true });
    assert.equal(loaded.status, 'legacy', JSON.stringify(loaded.error));
    assert.deepEqual(loaded.repairs.filter(({ path }) => /feature_shapes|visibility/.test(path))
      .map(({ code, path, effect, stored, used }) => [code, path, effect, stored, used]), [
      ['setting-replaced', 'adv.feature_shapes.CDS', 'kept', '5', before.CDS],
      ['setting-replaced', 'adv.linear_accession_visibility', 'kept', 'true', before.accession],
      ['setting-replaced', 'adv.linear_length_visibility', 'kept', '0', before.length],
      ['setting-replaced', 'adv.feature_shapes.gene', 'dropped', 'bogus', undefined]
    ]);
    // The Load commits a new drawing state (D-15-112): read the shown one.
    const loadedDrawing = modeDrawing(state, 'linear');
    assert.deepEqual([loadedDrawing.adv.feature_shapes.CDS, Object.hasOwn(loadedDrawing.adv.feature_shapes, 'gene'),
      loadedDrawing.adv.feature_shapes.misc_feature, loadedDrawing.adv.linear_accession_visibility,
      loadedDrawing.adv.linear_length_visibility], [before.CDS, false, 'rectangle', before.accession, before.length]);
    await savedSession('BA-OD100');
  } finally {
    installShownDrawing(shownMode);
  }
});

// BA-OD99 (b): in a configuration file, a color rule or a whitelist
// entry that is no object is dropped and listed before the dialog, and the
// apply reads the lists with the same readers, so nothing is refused after
// [Load]. A Session drops and lists it too.
test('a configuration file or a Session drops and lists a rule or whitelist entry that is no object', async () => {
  const rule = { feat: 'CDS', qual: 'product', val: 'kinase', color: 'red', cap: 'Kinase' };
  const file = { form: { plot_title: 'Rules' }, rules: [null, rule, 5], whitelist: [null, { feat: 'CDS', qual: 'gene', key: 'dnaA' }] };
  let committed = false;
  const asked = [];
  const loaded = await importSession(input(file), {
    confirmLoadRepairs: async ({ items }) => { asked.push(items); return true; },
    rollbackState: spyOnPort(() => { committed = true; })
  });
  assert.equal(loaded.status, 'legacy', JSON.stringify(loaded.error));
  assert.ok(committed);
  assert.deepEqual(asked, [[
    'Settings: 2 entries of rules could not be read and were dropped.',
    'Settings: 1 entry of whitelist could not be read and was dropped.'
  ]]);
  const drawing = modeDrawing(state, state.mode.value);
  assert.deepEqual(JSON.parse(JSON.stringify(drawing.manualSpecificRules)), [{ ...rule, color: '#FF0000' }]);
  assert.deepEqual(JSON.parse(JSON.stringify(drawing.manualWhitelist)), [{ feat: 'CDS', qual: 'gene', key: 'dnaA' }]);
  for (const list of ['rules', 'whitelist']) {
    const session = fixture('legend-renamed-rows.v46.gbdraw-session.json.gz');
    session.modes.circular.config[list] = [null];
    const dialogs = [];
    const read = await importSession(input(session), { confirmLoadRepairs: async ({ items }) => { dialogs.push(items); return true; } });
    assert.deepEqual([read.status, dialogs], ['ok', [[`Settings: 1 entry of modes.circular.config.${list} could not be read and was dropped.`]]], list);
  }
});

// A configuration file is installed over the drawing: a settings group or a
// list it cannot use is not used, and the drawing keeps its current value, so
// its row says so, never "its default is used". A Session resets the drawing
// first, so there the default is used, and its row says that.
test('a configuration file group or list Load does not use keeps the drawing\'s value, and its row says so', async () => {
  const first = { form: { plot_title: 'Mine', legend: 'left' },
    adv: { label_font_size: 21, arrow_shaft_width_ratio: 0.5, arrow_head_length_ratio: 0.3 },
    losat: { blastp: { collinearInferOrthogroups: true, collinearMinAnchors: 4,
      hitLimitsByMode: { pairwise: { candidateLimit: 7, orthogroupMemberMaxHits: 3 } } } },
    circularConservation: { subject_gencode: 11, ring_width: 12 },
    rules: [{ feat: 'CDS', qual: 'product', val: 'x', color: '#ff0000', cap: 'Mine' }], whitelist: [{ feat: 'CDS', qual: 'gene', key: 'abc' }] };
  /** The values of the first file a kept group or list holds. */
  const held = (drawing) => JSON.stringify([
    drawing.form.plot_title, drawing.adv.label_font_size,
    drawing.losat.blastp.collinearInferOrthogroups, drawing.losat.blastp.collinearMinAnchors, drawing.losat.blastp.hitLimitsByMode.pairwise,
    drawing.circularConservation.subject_gencode, drawing.circularConservation.ring_width,
    drawing.manualSpecificRules.length, drawing.manualWhitelist.length
  ]);
  const firstHeld = JSON.stringify(['Mine', 21, true, 4, { candidateLimit: 7, orthogroupMemberMaxHits: 3 }, 11, 12, 1, 1]);
  // The arrow ratios: an `adv` that is not used keeps the drawing's; a file
  // whose `adv` does not hold them gives their defaults, the arrow an older file was drawn with.
  const ratios = (drawing) => [drawing.adv.arrow_shaft_width_ratio, drawing.adv.arrow_head_length_ratio];
  const ratiosAfter = (file) => (file.adv === undefined || (typeof file.adv === 'object' && !Array.isArray(file.adv)) ? [1, null] : [0.5, 0.3]);
  const kept = 'the drawing keeps its current value.';
  const cases = [
    [{ form: 5 }, `Settings: the value "5" (form) cannot be used; ${kept}`],
    [{ adv: 'x' }, `Settings: the value "x" (adv) cannot be used; ${kept}`],
    [{ adv: 5 }, `Settings: the value "5" (adv) cannot be used; ${kept}`],
    [{ adv: [] }, `Settings: the value "[]" (adv) cannot be used; ${kept}`],
    [{ losat: 5 }, `Settings: the value "5" (losat) cannot be used; ${kept}`],
    [{ losat: 'x' }, `Settings: the value "x" (losat) cannot be used; ${kept}`],
    [{ losat: { blastp: 5 } }, `Settings: the value "5" (losat.blastp) cannot be used; ${kept}`, 'linear'],
    [{ losat: { blastp: [] } }, `Settings: the value "[]" (losat.blastp) cannot be used; ${kept}`, 'linear'],
    [{ losat: { blastp: { collinearInferOrthogroups: true, hitLimitsByMode: 5 } } }, `Settings: the value "5" (losat.blastp.hitLimitsByMode) cannot be used; ${kept}`, 'linear'],
    [{ circularConservation: 5 }, `Settings: the value "5" (circularConservation) cannot be used; ${kept}`],
    [{ rules: 'x' }, `Settings: the value "x" (rules) cannot be used; ${kept}`],
    [{ whitelist: 5 }, `Settings: the value "5" (whitelist) cannot be used; ${kept}`],
    [{ adv: { feature_shapes: 5 } }, `Settings: the value "5" (adv.feature_shapes) cannot be used; ${kept}`],
    [{ form: { plot_title: {} } }, `Settings: the value of form.plot_title cannot be used; ${kept}`]
  ];
  const outcomes = [];
  // `losat.blastp` is a Linear setting and `circularConservation` a Circular one: each row is its mode's.
  const shownMode = state.mode.value;
  for (const [file, , mode = 'circular'] of cases) {
    installShownDrawing(mode);
    assert.equal((await importSession(input(first), { confirmLoadRepairs: async () => true })).status, 'legacy');
    const drawing = modeDrawing(state, state.mode.value);
    assert.deepEqual([held(drawing), ratios(drawing)], [firstHeld, [0.5, 0.3]]);
    const shapes = JSON.stringify(drawing.adv.feature_shapes);
    const asked = [];
    const loaded = await importSession(input(file), { confirmLoadRepairs: async ({ items }) => { asked.push(...items); return true; } });
    // The Load commits a new drawing state (D-15-112): read the shown one.
    const loadedDrawing = modeDrawing(state, state.mode.value);
    outcomes.push([JSON.stringify(file), loaded.status, asked, loaded.repairs.map(({ effect }) => effect),
      held(loadedDrawing), ratios(loadedDrawing), JSON.stringify(loadedDrawing.adv.feature_shapes) === shapes]);
  }
  installShownDrawing(shownMode);
  assert.deepEqual(outcomes, cases.map(([file, line]) => [JSON.stringify(file), 'legacy', [line], ['kept'], firstHeld, ratiosAfter(file), true]));
  // A Session's drawing is reset first: the same values are its default.
  for (const [field, value] of [['feature_shapes', 5], ['plot_title', {}]]) {
    const session = fixture('legend-renamed-rows.v46.gbdraw-session.json.gz');
    session.modes.circular.config[field === 'plot_title' ? 'form' : 'adv'][field] = value;
    const loaded = await importSession(input(session), { confirmLoadRepairs: async () => true });
    assert.deepEqual([loaded.status, loaded.repairs.filter(({ path }) => path.endsWith(field)).map(({ effect }) => effect)],
      ['ok', ['defaulted']], field);
  }
});

// An arrow ratio stored as its text reads as its number: no row, and no
// dialog, in a Session and in a configuration file. A ratio that is no number
// is refused before the dialog.
test('an arrow ratio stored as its text loads as its number with no row', async () => {
  const ratios = { arrow_shaft_width_ratio: '0.5', arrow_head_length_ratio: '0.3' };
  const noDialog = { confirmLoadRepairs: async () => assert.fail('a ratio read as its number asks nothing') };
  const session = fixture('legend-renamed-rows.v46.gbdraw-session.json.gz');
  Object.assign(session.modes.circular.config.adv, ratios);
  const fromSession = await importSession(input(session), noDialog);
  assert.deepEqual([fromSession.status, fromSession.repairs?.length], ['ok', 0], JSON.stringify(fromSession.error));
  const { adv } = modeDrawing(state, 'circular');
  assert.deepEqual([adv.arrow_shaft_width_ratio, adv.arrow_head_length_ratio], [0.5, 0.3]);
  const fromFile = await importSession(input({ form: { plot_title: 'Ratios' }, adv: { arrow_shaft_width_ratio: '0.25' } }), noDialog);
  assert.deepEqual([fromFile.status, fromFile.repairs?.length], ['legacy', 0], JSON.stringify(fromFile.error));
  assert.equal(modeDrawing(state, state.mode.value).adv.arrow_shaft_width_ratio, 0.25);
  const snapshot = loadedState();
  const refused = await importSession(input({ adv: { arrow_shaft_width_ratio: 'x' } }), noDialog);
  assert.deepEqual([refused.status, loadedState() === snapshot], ['error', true]);
});

// OV-824: a value the app's own writer stores loads with no row (BA-D25). A
// number input bound with `.number` stores '' when it is cleared, which reads
// as the value the drawing would have (a blank is the default, R7); the
// comparison filters hold a number or its text (the e-value input is text, and
// `resolveComparisonThresholds` reads both) and load as stored. The cleared
// fields are read from the template's bindings.
const TEMPLATE = readFileSync(new URL('../../gbdraw/web/index.html', import.meta.url), 'utf8');
// Save refuses a cleared Feature overlap tolerance, so no file holds it (OV-826).
const CLEARED_NUMBER_FIELDS = [...new Set([...TEMPLATE.matchAll(/<input\b[^>]*>/g)].map(([tag]) => (
  /\btype="number"/.test(tag) ? tag.match(/\bv-model\.number="((?:form|adv|losat\.blastp|circularConservation)\.\w+)"/)?.[1] : undefined
)).filter((at) => at && at !== 'adv.feature_overlap_tolerance_bp'))];
const WRITTEN_VALUES = [...new Map([
  ...CLEARED_NUMBER_FIELDS.map((at) => [at, '']),
  ['adv.evalue', 1], ['adv.min_bitscore', '50'], ['adv.identity', '70.5'], ['adv.alignment_length', '']
])];
test('a value the app writes into a setting of another kind loads with no row, in a Session and a configuration file (OV-824)', async () => {
  assert.ok(CLEARED_NUMBER_FIELDS.length >= 30, String(CLEARED_NUMBER_FIELDS.length));
  for (const at of ['adv.gc_content_tick_interval', 'losat.blastp.maxHits', 'circularConservation.subject_gencode']) {
    assert.ok(CLEARED_NUMBER_FIELDS.includes(at), at);
  }
  /** @param {Record<string, any>} config @param {string} at */
  const place = (config, at) => {
    const parts = at.split('.');
    return [parts.slice(0, -1).reduce((container, key) => container[key], config), parts.at(-1)];
  };
  /** Loads `file`; returns its status and the rows the dialog listed. */
  const load = async (file) => {
    const items = [];
    const loaded = await importSession(input(file), { confirmLoadRepairs: async (prompt) => { items.push(...prompt.items); return true; } });
    return [loaded.status, items];
  };
  // A blank reads as '' or as unset (null).
  const valuesOf = (mode) => WRITTEN_VALUES.map(([at]) => {
    const [container, key] = place(modeDrawing(state, mode), at);
    return container[key] === null ? '' : container[key];
  });
  // Every value loads as stored, but a cleared number field reads as the
  // field's default (a blank is the default, R7; OV-827), so a written blank
  // that reads as the default counts as blank.
  const expected = () => WRITTEN_VALUES.map(([, value]) => value);
  const { defaultSettingsContainers } = await import('../../gbdraw/web/js/services/config-field-readers.js');
  /** @param {'circular' | 'linear'} mode */
  const defaultsOf = (mode) => WRITTEN_VALUES.map(([at]) => {
    const [container, key] = place(defaultSettingsContainers(mode), at);
    return container[key] === null ? '' : container[key];
  });
  /** @param {unknown[]} values @param {unknown[]} defaults */
  const blanked = (values, defaults) => values.map((value, index) => (
    WRITTEN_VALUES[index][1] === '' && value === defaults[index] ? '' : value));
  const base = fixture('legend-renamed-rows.v46.gbdraw-session.json.gz');
  assert.deepEqual(await load(base), ['ok', []]);
  const before = { circular: defaultsOf('circular'), linear: defaultsOf('linear') };
  const session = structuredClone(base);
  for (const mode of ['circular', 'linear']) {
    for (const [at, value] of WRITTEN_VALUES) {
      const [container, key] = place(session.modes[mode].config, at);
      container[key] = value;
    }
  }
  assert.deepEqual(await load(session), ['ok', []]);
  /** The places whose value differs from the written one. @param {unknown[]} values */
  const differing = (values) => WRITTEN_VALUES.flatMap(([at, value], index) => (values[index] === value ? [] : [`${at}=${JSON.stringify(values[index])}`]));
  for (const mode of ['circular', 'linear']) {
    const values = blanked(valuesOf(mode), before[mode]);
    assert.deepEqual(values, expected(), `${mode}: ${differing(values).join(', ')}`);
  }
  // A configuration file is applied over the shown drawing.
  assert.deepEqual(await load(base), ['ok', []]);
  const shown = state.mode.value;
  const configFile = { form: {}, adv: {}, losat: { blastp: {} }, circularConservation: {} };
  for (const [at, value] of WRITTEN_VALUES) {
    const [container, key] = place(configFile, at);
    container[key] = value;
  }
  assert.deepEqual(await load(configFile), ['legacy', []]);
  assert.deepEqual(blanked(valuesOf(shown), before[shown]), expected());
  // A blank is no value of a check box or a list: it is listed and not used.
  for (const file of [{ form: { suppress_gc: '' }, adv: { gc_content_show_axis: '', depth_tracks: '' } }]) {
    const [status, rows] = await load(file);
    assert.deepEqual([status, rows.length], ['legacy', 3], rows.join('\n'));
  }
  const blankSession = structuredClone(base);
  blankSession.modes.circular.config.form.suppress_gc = '';
  const [blankStatus, blankRows] = await load(blankSession);
  assert.deepEqual([blankStatus, blankRows.length], ['ok', 1], blankRows.join('\n'));
  // A text that is no number is not of a number filter's kind, so its row
  // says the drawing keeps its value; the e-value input holds text, so any text loads.
  const drawing = modeDrawing(state, shown);
  const kept = drawing.adv.min_bitscore;
  const [status, rows] = await load({ adv: { min_bitscore: 'abc', evalue: 'abc' } });
  assert.deepEqual([status, rows.length, rows[0]?.startsWith('Settings: the value "abc" (adv.min_bitscore) cannot be used;')], ['legacy', 1, true], rows.join('\n'));
  // Each Load commits a new drawing state (D-15-112): read the shown one.
  assert.deepEqual([modeDrawing(state, shown).adv.min_bitscore, modeDrawing(state, shown).adv.evalue], [kept, 'abc']);
  // A number filter reads text as `resolveComparisonThresholds` does
  // (`DECIMAL_NUMBER_PATTERN`): a hexadecimal text, which Number() reads as 16,
  // is no number of the filter, so it is listed and not used (BA-OD113 N1).
  const [hexStatus, hexRows] = await load({ adv: { min_bitscore: '0x10' } });
  assert.deepEqual([hexStatus, hexRows.length, hexRows[0]?.startsWith('Settings: the value "0x10" (adv.min_bitscore) cannot be used;')],
    ['legacy', 1, true], hexRows.join('\n'));
  assert.equal(modeDrawing(state, shown).adv.min_bitscore, kept);
});

// OV-827 (BA-OD113): a blank over a number reads as the field's default
// (R7), with no row, in a configuration file installed over the shown drawing
// and in a Session; a blank is never 0 for a reader whose default is not 0.
test('a blank over a number reads as the field\'s default with no row, in a configuration file and a Session (OV-827)', async () => {
  await loadAnswering(fixture('lambda_basic_linear.v44-schema8.gbdraw-session.json.gz'));
  // Each Load commits a new drawing state (D-15-112): edit and read the shown one.
  const edited = () => {
    const drawing = modeDrawing(state, state.mode.value);
    Object.assign(drawing.adv, { gc_content_tick_interval: 50, multi_record_column_gap_ratio: 0.3, multi_record_row_gap_ratio: 0.2 });
    Object.assign(drawing.losat.blastp, { collinearMaxConflictsInMergeGap: 4, collinearMaxUnitGap: 3 });
    drawing.linearRecordGap.value = 40;
  };
  const read = (drawing = modeDrawing(state, state.mode.value)) => [drawing.adv.gc_content_tick_interval,
    drawing.adv.multi_record_column_gap_ratio, drawing.adv.multi_record_row_gap_ratio,
    drawing.losat.blastp.collinearMaxConflictsInMergeGap, drawing.losat.blastp.collinearMaxUnitGap, drawing.linearRecordGap.value];
  // '' (a cleared input) reads as the field's default; null is copied as
  // stored and read by the field's reader, whose unset value is its default
  // (the GC tick interval's is null, auto), and 0 only where the default is 0.
  // A blank Linear record gap is its default, 24 px (OV-828, BA-D30).
  const READS = { '': [20, 0.10, 0.05, 1, 0, 24], null: [null, 0.10, 0.05, 1, 0, 24] };
  for (const blank of ['', null]) {
    edited();
    const file = {
      adv: { gc_content_tick_interval: blank, multi_record_column_gap_ratio: blank, multi_record_row_gap_ratio: blank },
      losat: { blastp: { collinearMaxConflictsInMergeGap: blank, collinearMaxUnitGap: blank } },
      linearRecordLayout: { recordGap: blank }
    };
    const loaded = await importSession(input(file), { confirmLoadRepairs: async ({ items }) => assert.fail(items.join('\n')) });
    assert.equal(loaded.status, 'legacy', JSON.stringify(loaded.error));
    assert.deepEqual(read(), READS[String(blank)], `blank ${JSON.stringify(blank)}`);
  }
  // A Session's drawing is reset first: the same owner reads the same defaults.
  const session = fixture('legend-renamed-rows.v46.gbdraw-session.json.gz');
  for (const mode of ['circular', 'linear']) {
    Object.assign(session.modes[mode].config.adv, { gc_content_tick_interval: '', multi_record_column_gap_ratio: '' });
    session.modes[mode].config.losat.blastp.collinearMaxConflictsInMergeGap = '';
  }
  session.modes.linear.config.linearRecordLayout = { ...session.modes.linear.config.linearRecordLayout, recordGap: '' };
  const { items } = await loadAnswering(session);
  assert.deepEqual(items, []);
  for (const mode of ['circular', 'linear']) {
    const { adv, losat } = modeDrawing(state, mode);
    assert.deepEqual([adv.gc_content_tick_interval, adv.multi_record_column_gap_ratio, losat.blastp.collinearMaxConflictsInMergeGap], [20, 0.10, 1], mode);
  }
  assert.equal(modeDrawing(state, 'linear').linearRecordGap.value, 24);
});

// OV-823 (BA-OD110): a Session 44 whose GUI-unmanaged font override is now
// the Font setting, with an `adv` that is no object, adopts the font and lists
// the `adv` value it drops, in the words a dropped `adv` has elsewhere; Load
// refuses no file it loaded before.
test('a Session 44 font override adopted into an adv that is no object lists the dropped adv (OV-823)', async () => {
  const session = fixture('lambda_basic_linear.v44-schema8.gbdraw-session.json.gz');
  session.config.adv = 5;
  session.config.unmanagedConfigOverrides = { 'objects.text.font_family': "'Pres Font'" };
  const dialogs = [];
  const loaded = await importSession(input(session), {
    confirmLoadRepairs: async ({ items }) => { dialogs.push(items); return true; },
    confirmSessionFonts: async () => true
  });
  assert.equal(loaded.status, 'ok', JSON.stringify(loaded.error));
  assert.deepEqual(dialogs, [['Settings: the value "5" (config.adv) cannot be used; its default is used.']]);
  assert.deepEqual(['circular', 'linear'].map((mode) => modeDrawing(state, mode).adv.font_family), ["'Pres Font'", "'Pres Font'"]);
  // An object `adv` takes the font with no row.
  const kept = fixture('lambda_basic_linear.v44-schema8.gbdraw-session.json.gz');
  kept.config.unmanagedConfigOverrides = { 'objects.text.font_family': "'Pres Font'" };
  const plain = await importSession(input(kept), {
    confirmLoadRepairs: async ({ items }) => assert.fail(items.join('\n')), confirmSessionFonts: async () => true
  });
  assert.equal(plain.status, 'ok', JSON.stringify(plain.error));
});

// OV-826 (BA-OD112): a cleared Feature overlap tolerance (the number input
// stores '') is its default (R7): Save writes the default, the Session reads
// back with it and no dialog, and the request draws with it.
test('a cleared Feature overlap tolerance saves and reads as its default (OV-826)', async () => {
  for (const [mode, name] of [['linear', 'lambda_basic_linear.v44-schema8.gbdraw-session.json.gz'],
    ['circular', 'legend-renamed-rows.v46.gbdraw-session.json.gz']]) {
    await loadAnswering(fixture(name));
    assert.equal(state.mode.value, mode);
    modeDrawing(state, mode).adv.feature_overlap_tolerance_bp = '';
    const saved = await savedSession(`tolerance-${mode}`);
    assert.equal(savedDrawingView(saved, mode).modes[mode].config.adv.feature_overlap_tolerance_bp, 0, mode);
    const { items } = await loadAnswering(saved);
    assert.deepEqual(items, []);
    assert.equal(modeDrawing(state, mode).adv.feature_overlap_tolerance_bp, 0, mode);
  }
});

// OV-828 (BA-D30): a cleared (or null) Linear record gap is its default,
// 24 px: Save writes 24 and the Session reads back with it and no dialog.
test('a cleared Linear record gap saves and reads as its default, 24 px (OV-828)', async () => {
  for (const blank of ['', null]) {
    await loadAnswering(fixture('lambda_basic_linear.v44-schema8.gbdraw-session.json.gz'));
    modeDrawing(state, 'linear').linearRecordGap.value = blank;
    const saved = await savedSession(`record-gap-${JSON.stringify(blank)}`);
    assert.equal(savedDrawingView(saved, 'linear').modes.linear.config.linearRecordLayout.recordGap, 24, JSON.stringify(blank));
    const { items } = await loadAnswering(saved);
    assert.deepEqual(items, []);
    assert.equal(modeDrawing(state, 'linear').linearRecordGap.value, 24, JSON.stringify(blank));
  }
});

// The shape of a feature type the user added is a key the default
// shapes lack; Load kept only the default keys and dropped it without a row.
test('a feature shape of a feature type the user added survives Save and Load', async () => {
  const loaded = await importSession(input(fixture('lambda_basic_linear.v44-schema8.gbdraw-session.json.gz')), {
    confirmLoadRepairs: async () => assert.fail('a file Load restores fully asks nothing')
  });
  assert.equal(loaded.status, 'ok', JSON.stringify(loaded.error));
  const { adv } = modeDrawing(state, 'linear');
  adv.features.push('gene');
  adv.feature_shapes.gene = 'arrow';
  const saved = await savedSession('shapes');
  const again = await importSession(input(saved), {
    confirmLoadRepairs: async () => assert.fail('a file Load restores fully asks nothing')
  });
  assert.equal(again.status, 'ok', JSON.stringify(again.error));
  assert.equal(modeDrawing(state, 'linear').adv.feature_shapes.gene, 'arrow');
});

// A stored feature label style outside its list (LP-D102) reads as the
// drawing's mode's packaged style, which the generated mode profile gives
// (Linear an underline on a box, Circular neither), and is listed.
test('a stored label style outside its list reads as its mode\'s default and is listed', async () => {
  const { configFieldRepairs, defaultSettingsContainers } = await import('../../gbdraw/web/js/services/config-field-readers.js');
  const { labelStyleForMode } = await import('../../gbdraw/web/js/mode-profiles.js');
  for (const mode of /** @type {const} */ (['circular', 'linear'])) {
    const targets = defaultSettingsContainers(mode);
    const repairs = configFieldRepairs(targets, { adv: { label_line: 'over', label_protection: 'bogus' } }, 'config', { mode });
    const { label_line: line, label_protection: protection } = labelStyleForMode(mode);
    assert.deepEqual([targets.adv.label_line, targets.adv.label_protection], [line, protection], mode);
    assert.deepEqual(describeLoadRepairs(repairs)?.items, [
      `Settings: the value "over" (config.adv.label_line) cannot be used; ${line} is used.`,
      `Settings: the value "bogus" (config.adv.label_protection) cannot be used; ${protection} is used.`
    ], mode);
  }
});

// A similarity-group name that is no text (an object, a boolean) is
// dropped and counted, not applied as its text; a number is its text. The
// older `priorityRules` is read only when `qualifierPriorityRules` is
// absent, as the apply reads it.
test('a group name that is no text is dropped and listed; priorityRules is read only without qualifierPriorityRules', async () => {
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
  assert.deepEqual({ ...modeDrawing(state, 'linear').orthogroupNameOverrides }, { [session.orthogroupState.groups[0].id]: 'Kept' });
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
  const { adv } = modeDrawing(state, 'circular');
  assert.deepEqual(adv.depth_tracks.map(({ label }) => label), ['Depth', 'Coverage B']);
  assert.deepEqual(adv.circular_track_slots.filter(({ renderer }) => renderer === 'depth')
    .map(({ id, enabled, params }) => [id, enabled, params.track_index ?? null]),
  [['my_depth', true, 1], ['lost_on', false, null], ['lost_named', false, null], ['off_unnamed', false, 0]]);
  // The Linear drawing's Circular stack is not drawn.
  const linear = await loadAnswering(depthEdits(fixture('legend-renamed-rows.v46.gbdraw-session.json.gz'), 'linear'));
  assert.deepEqual(linear.items, []);
  assert.deepEqual(modeDrawing(state, 'linear').adv.depth_tracks.map(({ label }) => label), ['Depth', 'Coverage B']);
});

// BA-OD94: in either mode, Add Depth track
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
    const { adv } = modeDrawing(state, mode);
    assert.equal(adv.depth_tracks.length, 1);
    // Add Depth track (Linear pads each record's Depth slots to the new series).
    adv.depth_tracks.push({ ...normalizeDepthTrackConfig(null, 1), label: 'Coverage B', small_tick_interval: 5 });
    if (mode === 'linear') modeDrawing(state, 'linear').linearSeqs.forEach((seq) => { seq.depth = padDepthFileSlots(seq.depth, 2); });
    // The row's series 1 had a file; clearing it keeps the row with geometry.
    const stack = adv[stackKey];
    const { slots } = reconcileManagedDepthSlots({
      slots: [...stack, { ...manualRow, params: { track_index: 1, legend_label: 'Coverage B' } }], previousSourced: [1], sourced: []
    });
    stack.splice(0, stack.length, ...slots);
    const saved = await savedSession(`depth-${mode}`);
    assert.equal(savedDrawingView(saved, mode).modes[mode].config.adv.depth_tracks.length, 2);
    const { loaded, items } = await loadAnswering(saved);
    assert.deepEqual(items, []);
    assert.deepEqual(codes(loaded), []);
    const again = await savedSession(`depth-${mode}-again`);
    assert.deepEqual(savedDrawingView(again, mode).modes[mode].config.adv, savedDrawingView(saved, mode).modes[mode].config.adv);
  });
}

// OV-761 (BA-OD114): a Save of a mode not generated yet builds its request
// from the draft. An enabled Depth row whose series has no source (its file
// was cleared; the validator's `depth_source_missing`) draws nothing, so the
// request leaves it out; the Session's settings keep the row and Load restores
// it. Generate still stops on the row issue (PD-OI-083; pinned by "a Linear
// File clear leaves a manual Depth row whose row issue stops Generate" in
// gui-audit-20260930-options.playwright.spec.js).
const UNSOURCED_DEPTH = {
  linear: ['lambda_basic_linear.v44-schema8.gbdraw-session.json.gz', 'linear_track_slots',
    { id: 'my_depth', renderer: 'depth', enabled: true, side: 'below', height: '12px', spacing: '', z: 0, params: { track_index: 0 } }],
  circular: ['legend-renamed-rows.v46.gbdraw-session.json.gz', 'circular_track_slots',
    { id: 'my_depth', renderer: 'depth', enabled: true, side: 'inside', width: '12', z: 0, params: { track_index: 0 } }]
};
for (const [mode, [name, stackKey, manualRow]] of Object.entries(UNSOURCED_DEPTH)) {
  for (const cleared of [null, [null]]) {
    test(`a ${mode} Save before Generate leaves an enabled Depth row without a source out of the request and keeps it in the settings (OV-761, depth ${JSON.stringify(cleared)})`, async () => {
      const { canonicalRenderArtifactOwner } = await import('../../gbdraw/web/js/services/config.js');
      await loadAnswering(fixture(name));
      assert.equal(state.mode.value, mode);
      const { adv } = modeDrawing(state, mode);
      adv[`${stackKey}_enabled`] = true;
      adv[stackKey].push(structuredClone(manualRow));
      if (mode === 'linear') modeDrawing(state, 'linear').linearSeqs.forEach((seq) => { seq.depth = structuredClone(cleared); });
      else bindDrawingInputFileTree(state.inputFiles, modeDrawing(state, 'circular').slices.inputs, 'c_depth', structuredClone(cleared));
      // No Generate since: no committed request and no Result.
      canonicalRenderArtifactOwner.restore(null);
      state.results.value = [];
      const saved = await savedSession(`unsourced-depth-${mode}`);
      const view = savedDrawingView(saved, mode);
      const tracks = view.renderRequest.diagramOptions.tracks;
      const requested = (mode === 'linear' ? tracks.linearTrackSlots : tracks.circularTrackSlots) ?? [];
      assert.deepEqual(requested.filter(({ renderer }) => renderer === 'depth'), []);
      assert.ok(requested.length > 0);
      const stored = view.modes[mode].config.adv;
      assert.equal(stored[`${stackKey}_enabled`], true);
      assert.deepEqual(stored[stackKey].filter(({ id }) => id === 'my_depth').map(({ enabled }) => enabled), [true]);
      const { items } = await loadAnswering(saved);
      assert.deepEqual(items, []);
      assert.deepEqual(modeDrawing(state, mode).adv[stackKey].filter(({ id }) => id === 'my_depth').map(({ enabled }) => enabled), [true]);
    });
  }
}

// A History restore applies captured settings with the same readers and
// lists nothing (only Load asks): what it installs is unchanged (R11), but a
// blank (null) reads as the reader's default (OV-827), the Linear record
// gap's too (OV-828, BA-D30).
test('a History restore reads a blank as the reader\'s default, other values as before, and asks nothing', async () => {
  const { applyCapturedConfigData } = await import('../../gbdraw/web/js/services/config.js');
  const drawing = modeDrawing(state, 'circular');
  applyCapturedConfigData(drawing, {
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
  ['auto', 0.10, 0.05, null, 60, null, false, true, 'curve']);
  assert.deepEqual(adv.multi_record_positions, [{ selector: '#1', row: 1 }, { selector: '#2', row: 2 }]);
  assert.deepEqual(adv.depth_tracks.map(({ label, height }) => [label, height]), [['D', null]]);
  assert.deepEqual([losat.blastp.mode, losat.blastp.collinearMaxConflictsInMergeGap, losat.blastp.orthogroupMemberMaxHits,
    Object.hasOwn(losat.blastp, 'orthogroupMaxHits'), Object.hasOwn(losat.blastp, 'orthogroupHitPolicy')], ['orthogroup', 1, null, false, false]);
  assert.deepEqual([losat.blastp.collinearColorMode, losat.blastp.collinearAnchorMode, losat.blastp.collinearSearchScope],
    ['average_identity', 'all', 'adjacent']);
  assert.deepEqual([state.losatExecution.executionMode, state.losatExecution.totalThreadBudget], ['auto', 'safe']);
  assert.deepEqual([circularConservation.source, circularConservation.series.map(({ color, losat_gencode, sourceIndex }) => [color, losat_gencode, sourceIndex])],
    ['upload', [['#4e79a7', 1, 0]]]);
  assert.deepEqual([drawing.linearRecordLayoutEnabled.value, drawing.linearRecordGap.value, [...drawing.linearRecordRows]],
    [true, 24, [{ uid: 'a', row: 2 }, { uid: 'b', row: 3, canonicalRow: 3, canonicalColumn: 1 }]]);
  assert.equal(drawing.importedComparisonIntent.action, null);
  assert.equal(drawing.losatProgram.value, 'blastn');
});

// Request schemas 1 and 2 (PR-D): an invalid Linear comparison height of the
// stored full config is replaced by its default and listed (B17, F1), on the
// direct read of a schema-1 Session and on the Gallery promotion of a schema-2
// one; a flat override name loads as the decoder reads it (OV-799).
test('a schema 1/2 request lists an invalid comparison height of its full config and loads (B17)', async () => {
  for (const name of ['lambda_basic_linear.v31-schema1.gbdraw-session.json.gz', 'BGC0000708-BGC0000713.schema-v2.gbdraw-session.json.gz']) {
    const session = fixture(name);
    session.renderRequest.diagramOptions.config.canvas.linear.comparison_height = 0;
    const dialogs = [];
    const loaded = await importSession(input(session), {
      transformLegacyResultSvg, confirmLoadRepairs: async ({ items }) => { dialogs.push(items); return true; }
    });
    assert.equal(loaded.status, 'ok', `${name}: ${JSON.stringify(loaded.error)}`);
    assert.deepEqual(loaded.repairs.filter(({ code }) => code === 'setting-replaced').map(({ path, effect, stored }) => [path, effect, stored]),
      [['renderRequest.diagramOptions.config.canvas.linear.comparison_height', 'defaulted', '0']], name);
    assert.ok(dialogs.flat().some((line) => line.includes('comparison_height') && line.includes('its default is used')), name);
    assert.equal(modeDrawing(state, 'linear').adv.comparison_height, 60, name);
  }
});

test('a flat configuration override name of a schema-1 request loads as the decoder reads it (OV-799)', async () => {
  const session = fixture('HmmtDNA_basic_circular.v32-schema1.gbdraw-session.json.gz');
  delete session.renderRequest.diagramOptions.config;
  session.renderRequest.diagramOptions.configOverrides = { show_gc: false };
  const loaded = await importSession(input(session), { transformLegacyResultSvg, confirmLoadRepairs: async ({ items }) => assert.fail(items.join('\n')) });
  assert.equal(loaded.status, 'ok', JSON.stringify(loaded.error));
  assert.equal(modeDrawing(state, 'circular').form.suppress_gc, true);
});

// Only a stored request is read with the decoder's flat names (D-OD11): a
// draft's unmanaged overrides keep the base check. Load drops a path that
// check refuses, with a row (OV-764), and Save still refuses it.
test('a flat override name in a draft is not read: Load drops it with a row and Save refuses it (OV-799, OV-764)', async () => {
  for (const [overrides, reason] of [[{ gc_stroke_color: 'red' }, 'UNKNOWN_CONFIG_PATH'], [{ label_blacklist: 'foo,bar' }, 'UNKNOWN_CONFIG_PATH']]) {
    const [key] = Object.keys(overrides);
    const session = fixture('legend-renamed-rows.v46.gbdraw-session.json.gz');
    session.modes.circular.config.unmanagedConfigOverrides = overrides;
    const loaded = await importSession(input(session), { transformLegacyResultSvg, confirmLoadRepairs: async () => true });
    assert.equal(loaded.status, 'ok', JSON.stringify(loaded.error));
    assert.deepEqual(loaded.repairs.map(({ code, path, effect }) => [code, path, effect]),
      [['setting-replaced', `modes.circular.config.unmanagedConfigOverrides["${key}"]`, 'dropped']], JSON.stringify(overrides));
    assert.equal(Object.hasOwn(modeDrawing(state, 'circular').unmanagedConfigOverrides, key), false);
    modeDrawing(state, 'circular').unmanagedConfigOverrides[key] = overrides[key];
    await assert.rejects(exportSession('flat-override', { askChoice: async () => 'download' }), (error) => {
      assert.deepEqual([error.code, error.context?.field, error.context?.reason], ['INPUT_INVALID', 'configOverrides', reason], JSON.stringify(error));
      return true;
    }, JSON.stringify(overrides));
    delete modeDrawing(state, 'circular').unmanagedConfigOverrides[key];
  }
});

// G3: every Gallery Session and every Session fixture loads without a repair;
// a repair there is a producer bug. Two fixtures are lossy on purpose. Node
// supplies no diagram Worker: a Session older than 40 reads its features again
// through it, so the fixtures below list that recovery here (and only because
// the Worker is missing: each asks for it). The browser specs load them (OV-751):
// gui-audit-20260930-inputs and live-generate-parity assert no dialog for
// seven; session-loading-feedback answers the one table row of
// whitelist-tab-keyword.v39 and loads the two mode-switch v30 fixtures without
// a dialog: their recovery reads the inputs of the Result's mode (OV-677). Five fixtures need a Worker service to load at
// all; feature-identity-overrides and gui-audit-20260930-comparisons load them.
const LOSSY_FIXTURES = {
  'feature-edits-circular-copies.v44.gbdraw-session.json.gz': ['feature-visibility-narrowed'],
  'whitelist-tab-keyword.v39.gbdraw-session.json.gz': ['table-row-joined']
};
const FEATURES_NEED_DIAGRAM_WORKER = new Set([
  'BGC0000708-BGC0000713.schema-v2.gbdraw-session.json.gz', 'BGC0000708-BGC0000713.v29.gbdraw-session.json.gz',
  'BGC0000708-BGC0000713.v39.gbdraw-session.json.gz', 'cli-linear-protein.v30.gbdraw-session.json.gz', 'conservation-fasta.v39.gbdraw-session.json.gz',
  'lambda_basic_linear.v31-schema1.gbdraw-session.json.gz', 'precomputed-ring-unlabelled-cli.v30.gbdraw-session.json.gz', 'scale-interval-negative-linear-cli.v30.gbdraw-session.json.gz',
  'scale-interval-zero-circular-cli.v30.gbdraw-session.json.gz', 'whitelist-tab-keyword.v39.gbdraw-session.json.gz',
  'mode-switch-circular-result-linear-shown.v30.gbdraw-session.json.gz', 'mode-switch-linear-result-circular-shown.v30.gbdraw-session.json.gz',
  'legacy-font-size-zero-circular-cli.v30.gbdraw-session.json.gz', 'legacy-font-size-zero-linear-cli.v30.gbdraw-session.json.gz',
  // gui-audit-20260930-comparisons loads it in the browser (OV-1181).
  'BGC0000708-BGC0000709-rc.v39.gbdraw-session.json.gz'
]);
const NEEDS_DIAGRAM_WORKER = new Set([
  'feature-edits-circular.v33.gbdraw-session.json.gz', 'feature-edits-linear-crop-rc.v33.gbdraw-session.json.gz',
  'hash-color-rule-linear-rc.v30.gbdraw-session.json.gz', 'q-frame-main-web-losatn.v42.gbdraw-session.json.gz',
  'q-frame-main-web-upload.v42.gbdraw-session.json.gz'
]);
// G2b (E-OD02): the document Save writes after each of those Loads reads
// again with no repair of the document reader (`readSessionDocument`); the
// settings readers' rows are G3's, checked on the stored file above.
test('every Gallery Session and Session fixture loads without a repair (G3), and the document its Save writes reads again without a document repair (G2b)', async () => {
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
  // The protein result conversion answers "nothing changed"; it is asked for
  // exactly the Sessions that stored protein results of a reversed record.
  const askedProteinConversion = new Set();
  let current = '';
  setStoredProteinResultFrameConverter(async () => { askedProteinConversion.add(current); return { result: { resources: {} } }; });
  for (const [name, url] of corpus) {
    const bytes = readFileSync(url);
    const session = JSON.parse((bytes[0] === 0x1f && bytes[1] === 0x8b ? gunzipSync(bytes) : bytes).toString('utf8'));
    askedDiagramWorker = false;
    current = name;
    const loaded = await importSession(input(session), { transformLegacyResultSvg, confirmLoadRepairs: async () => true });
    if (NEEDS_DIAGRAM_WORKER.has(name)) {
      assert.equal(loaded.status, 'error', `${name} loads in node now: check it here`);
      continue;
    }
    assert.equal(loaded.status, 'ok', `${name}: ${JSON.stringify(loaded.error)}`);
    const nodeOnly = FEATURES_NEED_DIAGRAM_WORKER.has(name) ? ['feature-recovery-failed'] : [];
    assert.deepEqual(loaded.repairs.map(({ code }) => code), [...(LOSSY_FIXTURES[name] ?? []), ...nodeOnly], name);
    assert.equal(askedDiagramWorker, nodeOnly.length > 0, `${name}: the Features row comes from the missing diagram Worker`);
    // A Session 27-39 with a Result is saved only after Generate.
    const saved = await savedSession(`G2b ${name}`).catch((error) => {
      if (error?.code !== 'SESSION_SAVE_REQUIRES_GENERATE' || session.version >= 40) throw error;
      return null;
    });
    if (saved) assert.deepEqual(readSessionDocument(saved).repairs, [], `${name}: G2b`);
  }
  globalThis.Worker = ImportWorker;
  setStoredProteinResultFrameConverter(null);
  assert.deepEqual([...askedProteinConversion].sort(), [
    'BGC0000708-BGC0000709-rc.v39.gbdraw-session.json.gz', 'BGC0000708-BGC0000709-rc.v40-early.gbdraw-session.json.gz',
    'BGC0000708-BGC0000713.v47.gbdraw-session.json.gz', 'majanivirus_collinear-rc.v47.gbdraw-session.json.gz'
  ]);
});

// OV-1181: a Session 39 converts too. Its Load commits a detached copy of the
// stored document, which takes the returned bytes, and the Linear draft's
// comparisons are projected again from them.
test('a Session 39 with a protein table of a reversed record installs the converted bytes Python returns at Load', async () => {
  const session = fixture('BGC0000708-BGC0000709-rc.v39.gbdraw-session.json.gz');
  const tableId = 'comparison-resolved-1';
  const converted = `${resourceText(session.resources[tableId])}# converted\n`;
  const calls = [];
  setStoredProteinResultFrameConverter(async (payload) => {
    calls.push({ version: payload.version, resources: Object.keys(payload.resources).sort() });
    return { result: { resources: { [tableId]: {
      data: Buffer.from(converted, 'utf8').toString('base64'), size: Buffer.byteLength(converted)
    } } } };
  });
  try {
    const loaded = await importSession(input(session), { transformLegacyResultSvg, confirmLoadRepairs: async () => true });
    assert.equal(loaded.status, 'ok', JSON.stringify(loaded.error));
    assert.deepEqual(calls, [{ version: 39, resources: Object.keys(session.resources).sort() }]);
    const committed = getCommittedCanonicalSession();
    assert.equal(resourceText(committed.resources[tableId]), converted, 'the committed copy holds the converted bytes');
    assert.equal(committed.resources[tableId].size, Buffer.byteLength(converted));
    const draft = modeDrawing(state, 'linear').files.linearCanonicalComparisons.find((item) => item.kind === 'precomputedProteinComparison');
    assert.equal(await state.inputFiles.fileOf(draft.file).text(), converted, 'the Linear draft holds the converted bytes');
  } finally {
    setStoredProteinResultFrameConverter(null);
  }
});

const COLLINEAR_RC = 'majanivirus_collinear-rc.v47.gbdraw-session.json.gz';
const COLLINEARITY_ID = 'comparison-canonical-collinearity-1';
const resourceText = (resource) => Buffer.from(resource.data, 'base64').toString('utf8');

// OV-1145 B-LOAD: the Collinear Session 47 with a reversed record asks Python
// once, before History's base; the returned bytes replace the stored ones
// under their ID in what Load commits, the Linear draft and the next Save.
test('a Collinear Session 47 with a reversed record installs the converted bytes Python returns at Load', async () => {
  const session = fixture(COLLINEAR_RC);
  const converted = `${resourceText(session.resources[COLLINEARITY_ID])}\n`;
  const calls = [];
  setStoredProteinResultFrameConverter(async (payload) => {
    calls.push({ version: payload.version, resources: Object.keys(payload.resources).sort() });
    return { result: { resources: { [COLLINEARITY_ID]: {
      data: Buffer.from(converted, 'utf8').toString('base64'), size: Buffer.byteLength(converted)
    } } } };
  });
  try {
    const loaded = await importSession(input(session), { transformLegacyResultSvg, confirmLoadRepairs: async () => assert.fail('no row') });
    assert.equal(loaded.status, 'ok', JSON.stringify(loaded.error));
    assert.deepEqual(calls, [{ version: 47, resources: Object.keys(session.resources).sort() }]);
    const draft = modeDrawing(state, 'linear').files.linearCanonicalComparisons.find((item) => item.kind === 'collinearityResult');
    assert.equal(await readSessionResourceText(state.inputFiles.fileOf(draft.file)), converted, 'the Linear draft holds the converted bytes');
    const saved = await savedSession('collinear rc');
    assert.equal(resourceText(saved.resources[COLLINEARITY_ID]), converted, 'Save writes the converted bytes');
    assert.equal(saved.resources[COLLINEARITY_ID].size, Buffer.byteLength(converted));
  } finally {
    setStoredProteinResultFrameConverter(null);
  }
});

// PD-OI-046: a conversion that fails fails the Load; the Session, Results
// and input files committed before it stay. The shown artifact's run
// metadata is the app's rollback extension (createSessionImportRollbackState),
// which this harness does not compose.
test('a failed protein result conversion fails the Load and keeps what was committed', async () => {
  await loadAnswering(fixture('lambda_basic_linear.v44-schema8.gbdraw-session.json.gz'));
  const before = loadedState();
  const savedBefore = await savedSession('before');
  setStoredProteinResultFrameConverter(async () => { throw new Error('probe: conversion failed'); });
  try {
    const loaded = await importSession(input(fixture('BGC0000708-BGC0000713.v47.gbdraw-session.json.gz')), {
      transformLegacyResultSvg, confirmLoadRepairs: async () => true
    });
    assert.equal(loaded.status, 'error');
  } finally {
    setStoredProteinResultFrameConverter(null);
  }
  assert.equal(loadedState(), before);
  const savedAfter = await savedSession('before');
  assert.deepEqual({ ...savedAfter, createdAt: null, runMetadata: null }, { ...savedBefore, createdAt: null, runMetadata: null });
});

// A Session whose font this page cannot draw also asks the Font question
// (PD-OI-098) before it loads; it comes after the repair list, so [Cancel]
// there loads no font file either, and [Cancel] in the font question still
// loads nothing.
test('the font question of a lossy Session comes after the repair dialog', async () => {
  await loadAnswering(fixture('lambda_basic_linear.v44-schema8.gbdraw-session.json.gz'));
  const saved = await savedSession('fonts');
  const circularSaved = withSavedDrawingView(saved, 'circular', (view) => { view.modes.circular.config.adv.font_family = null; });
  const session = withSavedDrawingView(circularSaved, 'linear', (view) => {
    view.modes.linear.config.adv.font_family = "'Missing Probe Font'";
    view.modes.linear.config.adv.linear_label_spacing = -5;
  });
  /** @param {boolean} loadRepairs @param {boolean} loadFont */
  const answered = async (loadRepairs, loadFont) => {
    const asked = [];
    const loaded = await importSession(input(session), {
      confirmLoadRepairs: async () => { asked.push('repairs'); return loadRepairs; },
      confirmSessionFonts: async (fonts) => { asked.push(['fonts', fonts]); return loadFont; }
    });
    return { status: loaded.status, asked };
  };
  const fonts = ['fonts', [null, "'Missing Probe Font'"]];
  assert.deepEqual(await answered(false, true), { status: 'canceled', asked: ['repairs'] });
  assert.deepEqual(await answered(true, false), { status: 'canceled', asked: ['repairs', fonts] });
  assert.deepEqual(await answered(true, true), { status: 'ok', asked: ['repairs', fonts] });
  assert.equal(modeDrawing(state, 'linear').adv.font_family, "'Missing Probe Font'");
});

// Review L3 (D-17): a Session table's label in a repair line is a key of the error wording
// owner (`SESSION_TABLE_LABELS`) and is shown in the UI language, inside the translated line.
test('a table repair line shows the table label in the UI language', async () => {
  const { SOURCE_LOCALE, setUiLocale } = await import('../../gbdraw/web/js/utils/i18n.js');
  try {
    await setUiLocale('ja', async () => ({ default: { 'Specific colors': '«Specific colors»' } }));
    const dialog = describeLoadRepairs([loadRepair({ code: 'table-row-joined', path: 'specific-colors', effect: 'changed', line: 3 }),
      loadRepair({ code: 'table-row-dropped', path: 'unknown-table', effect: 'dropped', line: 2 })]);
    assert.deepEqual(dialog.items, ['«Specific colors»: line 3 had extra cells, joined into the last column with one space.',
      'unknown-table: line 2 lacked a required column and was dropped.']);
  } finally {
    await setUiLocale(SOURCE_LOCALE);
  }
});
