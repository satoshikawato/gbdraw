// Session 46 (plan OV-80 §4, §8 S1, S3, S4, S6, S7): each diagram mode's
// drawing is saved as its own slice (`modes.<mode>`) and loaded back into its
// own drawing; a slice holds registry fields only; Session 45 and a Session 46
// with a retired field are rejected.
import assert from 'node:assert/strict';
import { test } from 'node:test';
import { readFileSync, readdirSync } from 'node:fs';
import { gunzipSync } from 'node:zlib';
import { installSessionImportWorker } from './helpers/session-import-node.mjs';
import { withBootPalettes } from './helpers/boot-palettes.mjs';
import { modeDrawing, withModeDrawings } from './helpers/drawing-state.mjs';
import { callerSession } from './helpers/session-caller.mjs';
import { SESSION_DRAWINGS_LAYOUT } from '../../gbdraw/web/js/contracts/schema.generated.js';
import { savedDrawingView, withSavedDrawingView } from './helpers/saved-drawing-view.mjs';
import { sessionRollbackPort } from './helpers/session-rollback-port.mjs';
import { evaluatePython, validateWebConfigOverrides } from './helpers/python-evaluator.mjs';

globalThis.window = {
  Vue: {
    ref: (value) => ({ value }), reactive: (value) => value,
    computed: (getter) => ({ get value() { return getter(); } }),
    nextTick: async () => {}
  },
  DOMPurify: { sanitize: (value) => value }
};
globalThis.document = {
  body: { appendChild: () => {} },
  createElement: () => ({ addEventListener: () => {}, click: () => {}, parentNode: null })
};
const alerts = [];
globalThis.alert = (message) => alerts.push(String(message));
installSessionImportWorker();

const sessionConfig = await import('../../gbdraw/web/js/services/config.js');
const { SESSION_VERSION, assertDrawingInputs, installShownDrawing, getCommittedCanonicalSession } = sessionConfig;
// A Session 49 drawing's saved unmanaged overrides are checked as Load checks them (PD-OI-009).
sessionConfig.setUnmanagedConfigOverrideValidator(async (payload) => ({ result: await validateWebConfigOverrides(payload) }));
const { MODE_SLICE_CONTAINERS } = await import('../../gbdraw/web/js/services/mode-scoped-migration.js');
const { createDefaultAdv, createDefaultForm } = await import('../../gbdraw/web/js/services/session-active-config-contract.js');
const { state } = await import('../../gbdraw/web/js/state.js');
const { exportSession, importSession } = callerSession(sessionConfig, state, await sessionRollbackPort(state));

const MODES = ['circular', 'linear'];
const save = async (title) => {
  const result = await exportSession(title);
  assert.equal(result.status, 'saved', result.error?.message);
  return JSON.parse(gunzipSync(Buffer.from(await result.blob.arrayBuffer())));
};
// The Load dialogs asked (BA-D25); the port answers [Load].
const loadDialogs = [];
const load = (session) => importSession({
  target: { files: [new Blob([JSON.stringify(session)], { type: 'application/json' })], value: 'selected' }
}, { confirmLoadRepairs: async (dialog) => { loadDialogs.push(dialog); return true; } });
const resetDrawings = () => {
  for (const mode of MODES) {
    Object.assign(modeDrawing(state, mode).form, createDefaultForm());
    Object.assign(modeDrawing(state, mode).adv, createDefaultAdv(mode));
  }
};
// Every key of a slice, as `container.key` paths.
const sliceKeys = (value, container = '') => Object.entries(value).flatMap(([key, child]) => {
  const path = container ? `${container}.${key}` : key;
  return [[container, key], ...(Object.hasOwn(MODE_SLICE_CONTAINERS, path) && child && typeof child === 'object'
    ? sliceKeys(child, path) : [])];
});

test('Session 49 saves each drawing as its own entry and loads it back into that drawing (S1, S6)', async () => {
  assert.equal(SESSION_VERSION, SESSION_DRAWINGS_LAYOUT.version);
  resetDrawings();
  installShownDrawing('circular');
  modeDrawing(state, 'circular').adv.label_font_size = 13;
  modeDrawing(state, 'circular').form.separate_strands = false;
  modeDrawing(state, 'circular').adv.gc_content_mode = 'percent';
  modeDrawing(state, 'linear').adv.label_font_size = 7;
  modeDrawing(state, 'linear').adv.depth_min = 3;
  modeDrawing(state, 'linear').form.plot_title = 'Linear only';
  // OV-120: a renamed row the last Linear Generate did not draw.
  // A Legend row holds no `showStroke` since OV-157 (view state, #940).
  const waiting = { caption: 'Coverage', originalCaption: 'depth', color: '#7b2cbf', featureIds: [] };
  modeDrawing(state, 'linear').dormantLegendEntries.value = [waiting];
  const saved = await save('two drawings');

  assert.equal(saved.version, SESSION_DRAWINGS_LAYOUT.version);
  // Each drawing is read as its one-drawing view, found by its ID.
  const views = Object.fromEntries(MODES.map((mode) => [mode, savedDrawingView(saved, mode)]));
  assert.equal(views.circular.renderRequest ?? null, null);
  for (const retired of ['config', 'features']) assert.equal(Object.hasOwn(saved, retired), false, retired);
  for (const retired of ['layoutPreferences', 'canvasPadding', 'pendingPaletteName', 'pendingPaletteColors',
    'linearTypographyLinked']) {
    assert.equal(Object.hasOwn(saved.ui, retired), false, `ui.${retired}`);
  }
  assert.equal(Object.hasOwn(views.circular.editorState ?? {}, 'featureStrokes'), false);
  assert.deepEqual(saved.drawings.map(({ id }) => id).sort(), MODES);
  for (const mode of MODES) assert.deepEqual(Object.keys(views[mode].modes), [mode]);
  assert.equal(views.circular.modes.circular.config.adv.label_font_size, 13);
  assert.equal(views.circular.modes.circular.config.form.separate_strands, false);
  assert.equal(views.circular.modes.circular.config.adv.gc_content_mode, 'percent');
  assert.equal(views.linear.modes.linear.config.adv.label_font_size, 7);
  assert.equal(views.linear.modes.linear.config.adv.depth_min, 3);
  assert.equal(views.linear.modes.linear.config.form.plot_title, 'Linear only');
  assert.equal(views.linear.modes.linear.config.form.separate_strands, createDefaultForm().separate_strands);
  // The waiting row follows the shown rows of its slice, marked `dormant`.
  assert.deepEqual(views.linear.modes.linear.editorState.legend.entries, [{ ...waiting, dormant: true }]);
  assert.deepEqual(views.circular.modes.circular.editorState.legend.entries, []);
  assert.equal(typeof saved.ui.losatExecution, 'object');
  assert.equal(typeof saved.ui.richFeaturePopup, 'boolean');

  // S7: a slice holds registry fields only.
  for (const mode of MODES) {
    const unknown = sliceKeys(views[mode].modes[mode]).filter(([container, key]) => !MODE_SLICE_CONTAINERS[container]?.has(key));
    assert.deepEqual(unknown, [], `modes.${mode} holds only registry fields`);
  }

  resetDrawings();
  installShownDrawing('linear');
  const loaded = await load(saved);
  assert.equal(loaded.status, 'ok', JSON.stringify(loaded.error));
  assert.equal(state.mode.value, 'circular');
  assert.equal(modeDrawing(state, 'circular').adv.label_font_size, 13);
  assert.equal(modeDrawing(state, 'circular').form.separate_strands, false);
  assert.equal(modeDrawing(state, 'circular').adv.gc_content_mode, 'percent');
  assert.equal(modeDrawing(state, 'linear').adv.label_font_size, 7);
  assert.equal(modeDrawing(state, 'linear').adv.depth_min, 3);
  assert.equal(modeDrawing(state, 'linear').form.plot_title, 'Linear only');
  assert.notEqual(modeDrawing(state, 'circular').form.plot_title, 'Linear only');
  assert.deepEqual(modeDrawing(state, 'linear').dormantLegendEntries.value, [waiting]);
  assert.deepEqual(modeDrawing(state, 'linear').legendEntries.value, []);
  assert.deepEqual(modeDrawing(state, 'circular').dormantLegendEntries.value, []);

  // A slice may omit any field, and the slice of a mode not shown may be
  // absent: Load fills them with that mode's defaults. (A settings-only
  // Session keeps the slice of its shown mode.)
  const partial = withSavedDrawingView(saved, 'circular', (view) => {
    delete view.modes.circular.config.adv.label_font_size;
  });
  partial.drawings = partial.drawings.filter(({ id }) => id !== 'linear');
  resetDrawings();
  modeDrawing(state, 'linear').dormantLegendEntries.value = [];
  modeDrawing(state, 'circular').adv.label_font_size = 99;
  modeDrawing(state, 'linear').adv.depth_min = 9;
  const loadedPartial = await load(partial);
  assert.equal(loadedPartial.status, 'ok', JSON.stringify(loadedPartial.error));
  assert.equal(modeDrawing(state, 'circular').adv.label_font_size, createDefaultAdv('circular').label_font_size);
  assert.equal(modeDrawing(state, 'circular').form.separate_strands, false);
  assert.equal(modeDrawing(state, 'linear').adv.depth_min, createDefaultAdv('linear').depth_min);
  assert.equal(modeDrawing(state, 'linear').form.plot_title, '');
});

test('a file that is not a gbdraw Session is refused as such before its contents are read', async () => {
  resetDrawings();
  const reason = async (value) => {
    const outcome = await load(value);
    assert.equal(outcome.status, 'error', JSON.stringify(value));
    return outcome.error.context?.reason;
  };
  for (const value of [null, [], 'x', 5, {}, { format: 'gbdraw-session' }, { format: 'other', version: 46 }]) {
    assert.equal(await reason(value), 'SESSION_FORMAT', JSON.stringify(value));
  }
  // The version is read first: a version this app does not read names the fields.
  assert.equal(await reason({ format: 'other', version: 45 }), 'FIELDS');
  assert.equal(await reason({ format: 'gbdraw-session', version: 48 }), 'FIELDS');
  // A second fault does not hide the format (Python's reader checks the format first too).
  assert.equal(await reason({ format: 'other', version: 39, bogus: 1 }), 'SESSION_FORMAT');
  assert.equal(await reason({ format: 'other', version: 46, bogus: 1 }), 'SESSION_FORMAT');
});

// S3, S4 (E-OD03, A7d, A7f): Load refuses a version it does not read; it drops
// a Session 46's retired homes, whose place the slices fill, with no row, and a
// key no writer wrote with one.
test('Session 45 is refused; a Session 46 drops its retired homes with no row and an unknown key with one (S3, S4)', async () => {
  resetDrawings();
  // The one-drawing view of a saved drawing is the Session 46 shape this test loads.
  const current = await save('current');
  const viewOf = (id) => ({ ...structuredClone(savedDrawingView(current, id)), version: 46 });
  const saved = viewOf('circular');
  const rows = async (session, label) => {
    const outcome = await load(session);
    assert.equal(outcome.status, 'ok', `${label}: ${JSON.stringify(outcome.error)}`);
    return (outcome.repairs ?? []).map(({ code, path }) => `${code}@${path}`);
  };
  const refused = await load({ ...structuredClone(saved), version: 45 });
  assert.equal(refused.status, 'error', 'Session 45');
  assert.deepEqual(await rows({ ...structuredClone(saved), config: structuredClone(saved.modes.circular.config) }, 'top-level config'), []);
  assert.deepEqual(await rows({ ...structuredClone(saved), features: {} }, 'top-level features'), []);
  const retiredUi = structuredClone(saved);
  retiredUi.ui.canvasPadding = { top: 0, right: 0, bottom: 0, left: 0 };
  assert.deepEqual(await rows(retiredUi, 'ui.canvasPadding'), []);
  const retiredLegend = structuredClone(saved);
  retiredLegend.editorState.legend = { ...retiredLegend.editorState.legend, entries: [] };
  assert.deepEqual(await rows(retiredLegend, 'editorState.legend.entries'), []);
  // A field no slice holds (a flat draft's `modeProfiles`) is dropped with a row (E-OD03).
  const unknownSliceField = structuredClone(saved);
  unknownSliceField.modes.circular.config.modeProfiles = { activeMode: 'circular' };
  assert.deepEqual(await rows(unknownSliceField, 'modeProfiles'), ['setting-replaced@modes.circular.config.modeProfiles']);
  const thirdMode = structuredClone(saved);
  thirdMode.modes.radial = {};
  assert.deepEqual(await rows(thirdMode, 'modes.radial'), ['session-field-dropped@modes.radial']);
});

// OV-160: a named CSS color in a Legend color or a stroke override (also the
// stroke an SVG had before the edit) loads as its CSS hex value, the value the
// Python split writes; an unknown name loads as no color, as before.
test('named Legend and stroke colors load as their CSS hex values (OV-160)', async () => {
  resetDrawings();
  installShownDrawing('circular');
  const saved = await save('named colors');
  const named = withSavedDrawingView(saved, 'circular', (view) => {
    const legend = view.modes.circular.editorState.legend;
    legend.colorOverrides = { CDS: 'SeaShell' };
    legend.strokeOverrides = {
      CDS: { strokeColor: 'navy', strokeWidth: 2, originalStrokeColor: 'gray', originalStrokeWidth: 1 }
    };
    view.modes.circular.editorState.featureStrokes = {
      overrides: { feature_1: { strokeColor: 'rebeccapurple', strokeWidth: 3, originalStrokeColor: 'darkgrey' } }
    };
    view.editorState = { ...view.editorState, originalSvgStroke: { color: 'gray', width: 1 } };
  });
  const loaded = await load(named);
  assert.equal(loaded.status, 'ok', JSON.stringify(loaded.error));
  const drawing = modeDrawing(state, 'circular');
  assert.deepEqual({ ...drawing.legendColorOverrides }, { CDS: '#fff5ee' });
  assert.deepEqual({ ...drawing.legendStrokeOverrides.CDS },
    { strokeColor: '#000080', strokeWidth: 2, originalStrokeColor: '#808080', originalStrokeWidth: 1 });
  assert.deepEqual({ ...drawing.featureStrokeOverrides.feature_1 },
    { strokeColor: '#663399', strokeWidth: 3, originalStrokeColor: '#a9a9a9' });
  assert.equal(state.originalSvgStroke.value.color, '#808080');

  const unknown = withSavedDrawingView(saved, 'circular', (view) => {
    view.modes.circular.editorState.legend.strokeOverrides = {
      CDS: { strokeColor: 'notacolor', strokeWidth: 2, originalStrokeColor: 'notacolor' }
    };
  });
  const loadedUnknown = await load(unknown);
  assert.equal(loadedUnknown.status, 'ok', JSON.stringify(loadedUnknown.error));
  assert.deepEqual({ ...modeDrawing(state, 'circular').legendStrokeOverrides.CDS }, { strokeWidth: 2, originalStrokeColor: null });
});

// D-43 (OV-300): a stored Legend entry whose color lies outside the Default
// colors (-d) domain is dropped at Load and named in the Load dialog; a named
// color is stored as its table hex and a hex as written (D-30). It holds for a
// Session 46 mode slice and for the top-level editorState of an older Session.
test('a stored Legend entry color outside the Default colors domain is dropped and named at Load (OV-300)', async () => {
  const stored = [
    { caption: 'X', color: 'buttonface', featureIds: [] },
    { caption: 'Y', color: 'Red', featureIds: [] },
    { caption: 'Z', color: '#AABBCC', featureIds: [] }
  ];
  const shownRows = (drawing) => drawing.legendEntries.value.map(({ caption, color }) => ({ caption, color }));
  const kept = [{ caption: 'Y', color: '#FF0000' }, { caption: 'Z', color: '#AABBCC' }];

  resetDrawings();
  installShownDrawing('circular');
  const current = withSavedDrawingView(await save('legend colors'), 'circular', (view) => {
    const legend = view.modes.circular.editorState.legend;
    legend.entries = structuredClone(stored);
    legend.deletedEntries = [{ caption: 'W', color: 'currentColor', featureIds: [] }];
  });
  alerts.length = 0;
  loadDialogs.length = 0;
  const loaded = await load(current);
  assert.equal(loaded.status, 'ok', JSON.stringify(loaded.error));
  assert.deepEqual(shownRows(modeDrawing(state, 'circular')), kept);
  assert.deepEqual(modeDrawing(state, 'circular').deletedLegendEntries.value, []);
  assert.deepEqual(loadDialogs.map(({ items }) => items),
    [['Legend: entries "X" and "W" (deleted) had a color the app does not accept and were dropped.']]);
  assert.deepEqual(loaded.notices, []);
  assert.deepEqual(alerts, []);

  const older = JSON.parse(gunzipSync(readFileSync(new URL('../fixtures/sessions/settings-only.v42.json.gz', import.meta.url))));
  // A caption whose first row has a bad color but a later row a valid one is shown, not named.
  older.editorState.legend.entries = [{ caption: 'Z', color: 'notacolor', featureIds: [] }, ...structuredClone(stored)];
  resetDrawings();
  alerts.length = 0;
  loadDialogs.length = 0;
  const loadedOlder = await load(older);
  assert.equal(loadedOlder.status, 'ok', JSON.stringify(loadedOlder.error));
  assert.deepEqual(shownRows(state.activeDrawing()), kept);
  assert.deepEqual(loadDialogs.map(({ items }) => items),
    [['Legend: entry "X" had a color the app does not accept and was dropped.']]);
  assert.deepEqual(alerts, []);
});

// PD-OI-061 revision 5: Session 47 keys a slice's Legend color and stroke
// edits by the row key, so two rows of one name keep their own edits; a
// Session 46 stored them under the shown caption and is read by row key
// (tests/fixtures/sessions/legend-row-key-vectors.json).
test('Session 47 keeps the edits of two rows of one name; a Session 46 is read by row key', async () => {
  resetDrawings();
  installShownDrawing('circular');
  const drawing = modeDrawing(state, 'circular');
  const rows = [
    { caption: 'GC', originalCaption: 'GC content', color: '#a0a0a0', featureIds: [] },
    { caption: 'GC', originalCaption: 'GC skew', color: '#5050f0', featureIds: [] }
  ];
  drawing.legendEntries.value = structuredClone(rows);
  Object.assign(drawing.legendColorOverrides, { 'GC content': '#ff0000', 'GC skew': '#0000ff' });
  Object.assign(drawing.legendStrokeOverrides, { 'GC skew': { strokeColor: '#00aa00', strokeWidth: 2 } });
  const saved = await save('two rows of one name');
  const legend = savedDrawingView(saved, 'circular').modes.circular.editorState.legend;
  assert.deepEqual(legend.entries.map(({ caption, originalCaption }) => [caption, originalCaption]),
    [['GC', 'GC content'], ['GC', 'GC skew']]);
  assert.deepEqual(legend.colorOverrides, { 'GC content': '#ff0000', 'GC skew': '#0000ff' });

  resetDrawings();
  drawing.legendEntries.value = [];
  for (const store of [drawing.legendColorOverrides, drawing.legendStrokeOverrides]) {
    Object.keys(store).forEach((key) => delete store[key]);
  }
  const loaded = await load(saved);
  assert.equal(loaded.status, 'ok', JSON.stringify(loaded.error));
  // A Load makes new drawing states (B2d-2): read the loaded drawing again.
  let loadedDrawing = modeDrawing(state, 'circular');
  assert.notEqual(loadedDrawing, drawing);
  assert.deepEqual(loadedDrawing.legendEntries.value.map(({ caption, originalCaption }) => [caption, originalCaption]),
    [['GC', 'GC content'], ['GC', 'GC skew']]);
  assert.deepEqual({ ...loadedDrawing.legendColorOverrides }, { 'GC content': '#ff0000', 'GC skew': '#0000ff' });
  assert.deepEqual({ ...loadedDrawing.legendStrokeOverrides['GC skew'] }, { strokeColor: '#00aa00', strokeWidth: 2 });

  // DESIGN 1a: GC skew renamed Z, then GC content renamed GC skew and colored red.
  const older = { ...structuredClone(savedDrawingView(saved, 'circular')), version: 46 };
  older.modes.circular.editorState.legend = {
    ...older.modes.circular.editorState.legend,
    entries: [
      { caption: 'GC skew', originalCaption: 'GC content', color: '#a0a0a0', featureIds: [] },
      { caption: 'Z', originalCaption: 'GC skew', color: '#5050f0', featureIds: [] }
    ],
    colorOverrides: { 'GC skew': '#ff0000' },
    strokeOverrides: { Z: { strokeColor: '#00aa00', strokeWidth: 2 } }
  };
  const loadedOlder = await load(older);
  assert.equal(loadedOlder.status, 'ok', JSON.stringify(loadedOlder.error));
  loadedDrawing = modeDrawing(state, 'circular');
  assert.deepEqual({ ...loadedDrawing.legendColorOverrides }, { 'GC content': '#ff0000' });
  assert.deepEqual(Object.keys(loadedDrawing.legendStrokeOverrides), ['GC skew']);
});

// A Session does not save the label fit report, so a loaded Result has no
// "Leader lines cross" notice (PD-OI-097).
test('Session Load clears the label fit notice of the Result shown before', async () => {
  resetDrawings();
  installShownDrawing('circular');
  const gallery = JSON.parse(readFileSync(new URL('../../gbdraw/web/gallery/sessions/HmmtDNA_basic_circular.gbdraw-session.json', import.meta.url), 'utf8'));
  state.labelFitReports.value = [{ resultIndex: 0, resultName: 'fit.svg', rowToken: 'circular:1:0',
    report: { drawn: true, notice: ['Leader lines cross (5 pairs).'] } }];
  assert.equal((await load(gallery)).status, 'ok');
  assert.ok(state.results.value.length > 0, 'the Session has a Result');
  assert.deepEqual(state.labelFitReports.value, []);
});

// Record selection D-05: each slice keeps the keys of its OFF records; an
// omitted list draws every record. Load drops a key that names no record of
// its drawing without a row (BA-OD88, E-OD13), reads a repeat once, and lists
// a list of the wrong shape.
test('a slice keeps its OFF records, and omission draws every record (record-selection D-05)', async () => {
  resetDrawings();
  installShownDrawing('circular');
  const gallery = JSON.parse(readFileSync(new URL('../../gbdraw/web/gallery/sessions/HmmtDNA_basic_circular.gbdraw-session.json', import.meta.url), 'utf8'));
  assert.equal((await load(gallery)).status, 'ok');
  const uid = modeDrawing(state, 'linear').linearSeqs[0].uid;
  modeDrawing(state, 'linear').recordsOff.splice(0, Infinity, uid);
  modeDrawing(state, 'circular').recordsOff.splice(0, Infinity, '#2');
  const saved = await save('records off');
  const configOf = (document, mode) => savedDrawingView(document, mode).modes[mode].config;
  assert.deepEqual(configOf(saved, 'linear').recordsOff, [uid]);
  assert.deepEqual(configOf(saved, 'circular').recordsOff, ['#2']);

  modeDrawing(state, 'linear').recordsOff.splice(0);
  modeDrawing(state, 'circular').recordsOff.splice(0);
  assert.equal((await load(saved)).status, 'ok');
  assert.deepEqual([...modeDrawing(state, 'linear').recordsOff], [uid]);
  assert.deepEqual([...modeDrawing(state, 'circular').recordsOff], ['#2']);
  assert.deepEqual(modeDrawing(state, 'linear').drawnLinearSeqs.value, []);

  const omitted = withSavedDrawingView(withSavedDrawingView(saved, 'linear', (view) => {
    delete view.modes.linear.config.recordsOff;
  }), 'circular', (view) => {
    delete view.modes.circular.config.recordsOff;
  });
  assert.equal((await load(omitted)).status, 'ok');
  assert.deepEqual([...modeDrawing(state, 'linear').recordsOff], []);
  assert.deepEqual([...modeDrawing(state, 'circular').recordsOff], []);
  // Review F6: with every record ON a Save writes no list, so the Session keeps its bytes.
  const allOn = await save('records on');
  assert.equal('recordsOff' in configOf(allOn, 'linear'), false);
  assert.equal('recordsOff' in configOf(allOn, 'circular'), false);

  const loaded = async (mode, edit, label, { linear, circular, rows = [] }) => {
    const session = withSavedDrawingView(saved, mode, (view) => edit(view.modes[mode].config));
    const result = await load(session);
    assert.equal(result.status, 'ok', label);
    assert.deepEqual(result.repairs.map(({ code, path }) => `${code}@${path}`), rows, label);
    assert.deepEqual([...modeDrawing(state, 'linear').recordsOff], linear, label);
    assert.deepEqual([...modeDrawing(state, 'circular').recordsOff], circular, label);
  };
  await loaded('linear', (config) => { config.recordsOff = ['not-a-bound-card']; }, 'unbound Linear uid',
    { linear: [], circular: ['#2'] });
  await loaded('linear', (config) => { config.recordsOff = [uid, uid]; }, 'duplicate key', { linear: [uid], circular: ['#2'] });
  await loaded('circular', (config) => { config.recordsOff = ['contig_2']; }, 'Circular key is not #N',
    { linear: [uid], circular: [] });
  // Review F8: a list of the wrong shape is listed; every record is drawn.
  await loaded('circular', (config) => { config.recordsOff = '#2'; }, 'not a list',
    { linear: [uid], circular: [], rows: [`setting-replaced@drawings[${saved.drawings.findIndex(({ id }) => id === 'circular')}].config.recordsOff`] });
});

// Review #992 finding 1: a Session without a slice of a mode (a Gallery
// Session of one mode) draws every record of that mode too, and saves again.
test('a Session without a mode slice turns every record of that mode ON (record selection)', async () => {
  resetDrawings();
  const gallery = (name) => JSON.parse(readFileSync(new URL(`../../gbdraw/web/gallery/sessions/${name}.gbdraw-session.json`, import.meta.url), 'utf8'));
  assert.equal((await load(gallery('lambda_basic_linear'))).status, 'ok');
  modeDrawing(state, 'linear').recordsOff.splice(0, Infinity, modeDrawing(state, 'linear').linearSeqs[0].uid);
  modeDrawing(state, 'circular').recordsOff.splice(0, Infinity, '#2');
  assert.equal((await load(gallery('HmmtDNA_basic_circular'))).status, 'ok');
  assert.deepEqual([...modeDrawing(state, 'linear').recordsOff], []);
  assert.deepEqual([...modeDrawing(state, 'circular').recordsOff], []);
  await save('circular only, saved again');
});

// Review #992 finding 6: Reset Settings keeps which records are ON in both
// drawings (PD-OI-070, PD-OI-091).
test('Reset Settings keeps the OFF records of both drawings (record selection)', async () => {
  const { resetSettings } = await import('../../gbdraw/web/js/services/reset.js');
  modeDrawing(state, 'linear').recordsOff.splice(0, Infinity, 'card-b');
  modeDrawing(state, 'circular').recordsOff.splice(0, Infinity, '#3');
  state.drawingList().forEach((drawing) => resetSettings(state, drawing));
  assert.deepEqual([...modeDrawing(state, 'linear').recordsOff], ['card-b']);
  assert.deepEqual([...modeDrawing(state, 'circular').recordsOff], ['#3']);
  modeDrawing(state, 'linear').recordsOff.splice(0);
  modeDrawing(state, 'circular').recordsOff.splice(0);
});

// R15-21 (OV-1107, OV-1108): Reset Settings removes a loaded Session's saved
// label table (also when the bulk label edits are empty) and the dormant
// Similarity group names, in both drawings.
test('Reset Settings removes the saved label table and the dormant Similarity group names of both drawings', async () => {
  const { resetSettings } = await import('../../gbdraw/web/js/services/reset.js');
  const row = { recordId: '*', featureType: '*', qualifier: 'product', valueRegex: '^a$', labelText: 'A' };
  MODES.forEach((mode) => {
    const drawing = modeDrawing(state, mode);
    drawing.canonicalLabelOverrideRows.value = [row];
    drawing.orthogroupDormantOverrides.OG_1 = { name: 'Saved name' };
  });
  state.drawingList().forEach((drawing) => resetSettings(state, drawing));
  MODES.forEach((mode) => {
    const drawing = modeDrawing(state, mode);
    assert.deepEqual(drawing.canonicalLabelOverrideRows.value, [], `${mode} label table`);
    assert.deepEqual({ ...drawing.orthogroupDormantOverrides }, {}, `${mode} dormant group names`);
  });
});

// Record selection 2.5: a Linear drawing with every card OFF (only a
// hand-edited Session reaches it) fails Generate and Save with NONE_DRAWN.
test('a Linear drawing with no ON record names the fix', () => {
  const sourceState = withModeDrawings({
    files: {}, lInputType: { value: 'gb' }
  }, { circular: {}, linear: { recordsOff: ['a', 'b'], linearSeqs: [{ uid: 'a', gb: {} }, { uid: 'b', gb: {} }] } });
  assert.throws(() => assertDrawingInputs(modeDrawing(sourceState, 'linear'), sourceState),
    (error) => error.code === 'RECORD_SELECTION' && error.context.reason === 'NONE_DRAWN');
  modeDrawing(sourceState, 'linear').recordsOff = ['a'];
  assert.doesNotThrow(() => assertDrawingInputs(modeDrawing(sourceState, 'linear'), sourceState));
});

// OV-1280 (D-19, D-ORIGIN-2): the Web writer's Session 48 after the five file
// cards keeps each card's rows and what its file set through Load and Save,
// and the file-imports.js helpers each card's Remove calls take back only that
// file's effect (the Remove path itself: auxiliary-file-history Playwright cases).
test('Session 48 keeps what each file card loaded through Load and Save, and the card-mark helpers take back only that', async () => {
  const fileImports = await import('../../gbdraw/web/js/services/file-imports.js');
  const fixture = JSON.parse(gunzipSync(readFileSync(new URL('../fixtures/sessions/file-card-imports-circular.v48.gbdraw-session.json.gz', import.meta.url))));
  resetDrawings();
  const loaded = await load(fixture);
  assert.equal(loaded.status, 'ok', JSON.stringify(loaded.error));
  const drawing = modeDrawing(state, 'circular');
  const plain = (value) => JSON.parse(JSON.stringify(value));
  const read = () => ({
    rules: plain(drawing.manualSpecificRules), priority: plain(drawing.manualPriorityRules),
    whitelist: plain(drawing.manualWhitelist), blacklist: drawing.manualBlacklist.value, fileImports: plain(drawing.fileImports)
  });
  const stored = fixture.modes.circular.config;
  assert.deepEqual(read(), {
    rules: stored.rules, priority: stored.qualifierPriorityRules, whitelist: stored.whitelist,
    blacklist: stored.blacklistText, fileImports: stored.fileImports
  });
  assert.deepEqual(['d_color', 't_color', 'qualifier_priority', 'whitelist', 'blacklist'].map((key) => Boolean(state.cardFiles[key])),
    [true, true, true, true, true]);
  // Save writes them back as they were loaded.
  const saved = await save('file cards');
  for (const field of ['rules', 'qualifierPriorityRules', 'whitelist', 'blacklistText', 'fileImports']) {
    assert.deepEqual(savedDrawingView(saved, 'circular').modes.circular.config[field], stored[field], field);
  }
  // Remove on each card, as the card's owner takes it back.
  assert.deepEqual(fileImports.withoutImportedRows(drawing.manualSpecificRules, 't_color').map((rule) => [rule.val, rule.color]),
    [['hand', '#445566'], ['beta', '#999999']]);
  const priority = plain(drawing.manualPriorityRules);
  fileImports.removeImportedPriorityRows(priority);
  assert.deepEqual(priority, [{ feat: 'CDS', order: 'gene' }, { feat: 'tRNA', order: 'product' }, { feat: 'rRNA', order: 'locus_tag' }]);
  assert.deepEqual(fileImports.withoutImportedRows(drawing.manualWhitelist, 'whitelist').map((rule) => rule.key), ['mine', 'betas']);
  assert.equal(fileImports.withoutImportedBlacklistWords(drawing.manualBlacklist.value, drawing.fileImports.blacklist), 'hypothetical, mine');
  const colors = fileImports.withoutImportedDefaultColors(drawing.currentColors.value, drawing.fileImports.d_color, { CDS: '#aaaaaa', tRNA: '#bbbbbb' });
  assert.deepEqual([colors.CDS, colors.tRNA], ['#111111', '#abcdef']);
});

// B2d-2: Load builds one drawing per drawing of the Session (its 49 document)
// and commits them at once; the fast path reads 46-48 without Python.
const fixtureSession = (name) => JSON.parse(gunzipSync(readFileSync(new URL(`../fixtures/sessions/${name}.gbdraw-session.json.gz`, import.meta.url))));
const drawingsVector = (name) => JSON.parse(readFileSync(new URL('../fixtures/sessions/drawings-vectors.json', import.meta.url)))
  .upgrades.find((entry) => entry.name === name);
// Every message the session-import Worker receives: a file to read, or an ack.
const workerMessages = () => {
  const messages = [];
  const { postMessage } = globalThis.Worker.prototype;
  globalThis.Worker.prototype.postMessage = function record(message) {
    messages.push(message);
    return postMessage.call(this, message);
  };
  return { messages, restore: () => { globalThis.Worker.prototype.postMessage = postMessage; } };
};

test('Load of a two-mode Session 48 makes its drawings in the order and with the shown drawing the fast path names', async () => {
  const vector = drawingsVector('two-mode-linear-shown.v48');
  resetDrawings();
  installShownDrawing('circular');
  const before = state.drawingList();
  const worker = workerMessages();
  let loaded;
  try {
    loaded = await load(fixtureSession('two-mode-linear-shown.v48'));
  } finally {
    worker.restore();
  }
  assert.equal(loaded.status, 'ok', JSON.stringify(loaded.error));
  assert.deepEqual(state.drawingList().map((drawing) => drawing.id), vector.drawings.map((drawing) => drawing.id));
  assert.deepEqual(state.drawingList().map((drawing) => drawing.mode), vector.drawings.map((drawing) => drawing.mode));
  assert.equal(state.activeDrawing().id, vector.activeDrawingId);
  // Each drawing is a new drawing state, not the one shown before the Load.
  for (const drawing of state.drawingList()) assert.equal(before.includes(drawing), false, drawing.id);
  // The Worker read the file and ran nothing else (no Python helper operation).
  assert.ok(worker.messages.length > 0);
  for (const message of worker.messages) {
    assert.deepEqual(Object.keys(message).filter((key) => !['operationId', 'file', 'ack'].includes(key)), [], JSON.stringify(Object.keys(message)));
  }
});

test('a Load that fails in its commit leaves the drawings, their order and the shown drawing as they were', async () => {
  resetDrawings();
  installShownDrawing('linear');
  modeDrawing(state, 'linear').adv.label_font_size = 11;
  const before = state.drawingList();
  const loaded = await importSession({
    target: { files: [new Blob([JSON.stringify(fixtureSession('two-mode-linear-shown.v48'))], { type: 'application/json' })], value: 'selected' }
  }, {
    confirmLoadRepairs: async () => true,
    afterLoad: async () => { throw new Error('install failed'); }
  });
  assert.equal(loaded.status, 'error');
  assert.equal(state.drawingList().length, before.length);
  state.drawingList().forEach((drawing, index) => assert.equal(drawing, before[index], drawing.id));
  assert.equal(state.activeDrawing().id, 'linear');
  assert.equal(modeDrawing(state, 'linear').adv.label_font_size, 11);
});

test('a Load of a Session with a drawing of one mode keeps a new drawing of the other mode for the mode buttons (T-1)', async () => {
  const vector = drawingsVector('file-card-imports-circular.v48');
  assert.deepEqual(vector.drawings.map((drawing) => drawing.mode), ['circular']);
  resetDrawings();
  modeDrawing(state, 'linear').adv.label_font_size = 9;
  const before = state.drawingList();
  const loaded = await load(fixtureSession('file-card-imports-circular.v48'));
  assert.equal(loaded.status, 'ok', JSON.stringify(loaded.error));
  assert.deepEqual(state.drawingList().map((drawing) => drawing.id), ['circular', 'linear']);
  const linear = modeDrawing(state, 'linear');
  assert.equal(before.includes(linear), false);
  // A new drawing at its mode's defaults (D-15-120): the pre-Load edit is gone.
  assert.deepEqual({ ...linear.form }, createDefaultForm());
  assert.equal(linear.adv.label_font_size, createDefaultAdv('linear').label_font_size);
  installShownDrawing(linear.id);
  assert.equal(state.activeDrawing(), linear);
  installShownDrawing('circular');
});

test('Load reads each drawing\'s saved Result set from its own entry of the 49 document, not from `otherModeResult`', () => {
  // The Load body after preflight reads the upgraded document's drawings; the
  // file's location of a set beside the committed one is only a repair path,
  // which the upgrade gives by drawing ID (`setSources`).
  const source = readFileSync(new URL('../../gbdraw/web/js/services/config.js', import.meta.url), 'utf8');
  const body = source.slice(source.indexOf('const importSessionDocument = async'), source.indexOf('export const importSession ='));
  assert.ok(body.length > 1000);
  const reads = body.split('\n').filter((line) => /\botherModeResult\b/.test(line) && !/^\s*\/\/|^\s*\*/.test(line));
  assert.deepEqual(reads.map((line) => line.trim()), []);
  for (const name of ['otherSetMode', 'otherArtifactSlot', 'opensOtherSet', 'ofMode']) {
    assert.equal(new RegExp(`\\b${name}\\b`).test(body), false, name);
  }
});

test('Load reads each drawing\'s draft from its own part of the file by drawing ID, not from mode-keyed slices', () => {
  // A Session 46-48 drawing's settings, per-feature edits, editor edits and
  // `ui` come from its entry of the 49 document; a Session 27-44 split is
  // keyed by its drawings' IDs. The file's one Linear part (comparisons,
  // LOSAT caches, group names) is the drawing the file names for it.
  const source = readFileSync(new URL('../../gbdraw/web/js/services/config.js', import.meta.url), 'utf8');
  const body = source.slice(source.indexOf('const importSessionDocument = async'), source.indexOf('export const importSession ='));
  assert.ok(body.length > 1000);
  for (const name of ['modeSlices', 'modeConfigs', 'modeEditorStates', 'linearDraftConfig', 'circularDraftConfig']) {
    assert.equal(new RegExp(`\\b${name}\\b`).test(body), false, name);
  }
  assert.equal(/draftConfigOf\((?:'linear'|'circular'|mode|generatedMode)\)/.test(body), false, 'draftConfigOf by mode');
  assert.equal(/legendRekeysOf\((?:mode|generatedMode)\)/.test(body), false, 'legendRekeysOf by mode');
  assert.equal(/committed\.find\(\(drawing\) => drawing\.mode === 'linear'\)/.test(body), false, 'the Linear drawing by mode');
});

test('Load reads each Result set\'s request, reset receipt and catalog where the file holds it; each drawing holds its own, by drawing ID', async () => {
  // tests/utils/two_mode_session.py: a Circular Session 48 with a Linear `otherModeResult`.
  const session = await evaluatePython('from tests.utils.two_mode_session import two_mode_session\nresult = two_mode_session()');
  assert.deepEqual([session.version, session.renderRequest.mode, session.otherModeResult.renderRequest.mode], [48, 'circular', 'linear']);
  // A receipt that binds no plan, in each set: each reads as none, at its own path.
  session.editorState.alignmentResetReceipt = { schema: 1 };
  session.otherModeResult.editorState.alignmentResetReceipt = { schema: 1 };
  // An entry the request reader drops, in the Linear set only.
  Object.assign(session.otherModeResult.renderRequest.records[0], { bogus: true });
  // The Load hands its caller the artifact of each drawing it does not show.
  const stashed = [];
  const loaded = await importSession({
    target: { files: [new Blob([JSON.stringify(session)], { type: 'application/json' })], value: 'selected' }
  }, {
    confirmLoadRepairs: async () => true,
    installLoadedArtifactSlots: (slots) => { stashed.push(...slots.stashed); }
  });
  assert.equal(loaded.status, 'ok', JSON.stringify(loaded.error));
  assert.deepEqual(loaded.repairs.map(({ code, path }) => `${code} ${path}`), [
    'setting-entries-dropped otherModeResult.renderRequest.records[0]',
    'session-field-dropped editorState.alignmentResetReceipt',
    'session-field-dropped otherModeResult.editorState.alignmentResetReceipt'
  ]);
  assert.equal(state.activeDrawing().id, 'circular');
  const linear = stashed.find((slot) => slot.drawingId === 'linear');
  assert.ok(linear, JSON.stringify(stashed.map((slot) => slot.drawingId)));
  // [drawing ID, its request, its features, its receipt]: the shown drawing's, then the kept one's.
  const held = [
    ['circular', getCommittedCanonicalSession()?.renderRequest, state.extractedFeatures.value, state.similarityAlignmentResetReceipt.value],
    ['linear', linear.runtimeState.canonical.committedCanonicalSession.renderRequest, linear.values.extractedFeatures,
      linear.values.similarityAlignmentResetReceipt]
  ];
  for (const [id, request, features, receipt] of held) {
    const set = id === 'circular' ? session : session.otherModeResult;
    assert.equal(request.mode, set.renderRequest.mode, id);
    assert.equal(Object.hasOwn(request.records[0], 'bogus'), false, id);
    assert.equal(features.length, set.editorState.featureCatalog.items[0].features.length, id);
    assert.equal(receipt, null, id);
  }
  // The readers read one set each; only `storedResultSets` and its writer
  // name where the file holds a set, and no catalog is passed by position.
  const source = readFileSync(new URL('../../gbdraw/web/js/services/config.js', import.meta.url), 'utf8');
  const code = (name) => {
    const start = source.indexOf(`const ${name} = `);
    assert.ok(start >= 0, name);
    return source.slice(start, source.indexOf('\n};\n', start)).split('\n').filter((line) => !/^\s*(\/\/|\*|\/\*\*)/.test(line)).join('\n');
  };
  for (const name of ['readResultSetRequest', 'readResultSetResetReceipt', 'readResultSetFeatureCatalog', 'readSessionRequests']) {
    assert.equal(/\botherModeResult\b/.test(code(name)), false, name);
  }
  for (const name of ['otherModeCatalog', 'otherCatalog', 'readSessionResetReceipts', 'readSessionFeatureCatalogs']) {
    assert.equal(new RegExp(`\\b${name}\\b`).test(source), false, name);
  }
});

// Session 49 (D-15-109): the Web reads the drawings document Python writes.
const v49Session = (name) => JSON.parse(gunzipSync(readFileSync(new URL(`../fixtures/sessions/${name}.v49.gbdraw-session.json.gz`, import.meta.url))));
const loadKeepingSlots = async (session) => {
  const stashed = [];
  let opening = null;
  const loaded = await importSession({
    target: { files: [new Blob([JSON.stringify(session)], { type: 'application/json' })], value: 'selected' }
  }, {
    confirmLoadRepairs: async () => true,
    installLoadedArtifactSlots: (slots) => { opening = slots.opening; stashed.push(...slots.stashed); }
  });
  return { loaded, stashed, opening };
};
// Each drawing's held request and features, by drawing ID: the shown one's from state, any other from its slot.
const heldSetOf = (id, { stashed, opening }) => {
  if (id === state.activeDrawing().id && !opening) {
    return { request: getCommittedCanonicalSession()?.renderRequest ?? null, features: state.extractedFeatures.value };
  }
  const slot = [opening, ...stashed].find((entry) => entry?.drawingId === id);
  return slot
    ? { request: slot.runtimeState.canonical.committedCanonicalSession.renderRequest, features: slot.values.extractedFeatures }
    : null;
};
const assertLoadedDrawings = (session, outcome) => {
  assert.equal(outcome.loaded.status, 'ok', JSON.stringify(outcome.loaded.error));
  const ids = session.drawings.map((entry) => entry.id);
  assert.deepEqual(state.drawingList().map((drawing) => drawing.id).slice(0, ids.length), ids);
  assert.deepEqual(state.drawingList().map((drawing) => drawing.mode).slice(0, ids.length), session.drawings.map((entry) => entry.mode));
  assert.equal(state.activeDrawing().id, session.ui.activeDrawingId);
  for (const entry of session.drawings) {
    const drawing = state.drawingOf(entry.id);
    if (entry.config?.form?.plot_title !== undefined) assert.equal(drawing.form.plot_title, entry.config.form.plot_title, entry.id);
    if (!entry.renderRequest || entry.results.length === 0) continue;
    const held = heldSetOf(entry.id, outcome);
    assert.ok(held, `${entry.id} holds its Result set`);
    assert.equal(held.request.mode, entry.mode, entry.id);
    assert.deepEqual(held.request.records.map((record) => record.recordKey), entry.renderRequest.records.map((record) => record.recordKey), entry.id);
    assert.equal(held.features.length, entry.editorState.featureCatalog.items[0].features.length, entry.id);
  }
};

test('a Python-written Session 49 loads each drawing\'s request, Results and draft by drawing ID', async () => {
  for (const name of ['drawings-circular-only', 'drawings-circular-linear', 'drawings-linear-pair-losat', 'drawings-per-drawing-cards']) {
    const session = v49Session(name);
    // A draft value of its own in each drawing that holds a draft.
    session.drawings.forEach((entry) => { if (entry.config?.form) entry.config.form.plot_title = `title of ${entry.id}`; });
    resetDrawings();
    assertLoadedDrawings(session, await loadKeepingSlots(session));
  }
});

test('a Session 49 with two Circular drawings and a Linear one loads all three in tab order and keeps each Result in its own slot', async () => {
  const session = v49Session('drawings-circular-linear');
  const circular = session.drawings.find((entry) => entry.id === 'circular');
  session.drawings.push({ ...structuredClone(circular), id: 'circular-2', name: 'Circular 2' });
  session.drawings.forEach((entry) => { if (entry.config?.form) entry.config.form.plot_title = `title of ${entry.id}`; });
  session.ui.activeDrawingId = 'circular-2';
  resetDrawings();
  const outcome = await loadKeepingSlots(session);
  assertLoadedDrawings(session, outcome);
  // Each drawing that is not shown keeps its Result in its own slot.
  assert.deepEqual([outcome.opening, ...outcome.stashed].filter(Boolean).map((slot) => slot.drawingId).sort(), ['circular', 'linear']);
});

test('a Session 49 without its drawings or its resources is refused as no Session (OV-1768)', async () => {
  for (const field of ['drawings', 'resources']) {
    const session = v49Session('drawings-circular-only');
    delete session[field];
    resetDrawings();
    const outcome = await load(session);
    assert.equal(outcome.status, 'error', field);
    assert.equal(outcome.error.context?.reason, 'SESSION_FORMAT', `${field}: ${JSON.stringify(outcome.error)}`);
  }
});

// B2d-3a-2: each drawing of a Session 49 takes its own files, bindings, LOSAT
// entries and comparison intent by its drawing ID, opened or not.
// A Circular drawing and two Linear ones, each with files of its own; the
// Circular drawing's resources are renamed so the two fixtures share none.
const threeDrawingSession = () => {
  const session = v49Session('drawings-linear-pair-losat');
  const circularOnly = v49Session('drawings-circular-only');
  const rename = (value) => JSON.parse(JSON.stringify(value).replaceAll('"record-1-genbank"', '"circular-record-1"'));
  session.resources['circular-record-1'] = rename(circularOnly.resources['record-1-genbank']);
  session.drawings.unshift(rename(circularOnly.drawings[0]));
  const second = session.drawings.find((entry) => entry.id === 'linear-2');
  second.webFiles.bindings.linearSeqs[0].gb.name = 'R2c-second.gb';
  return session;
};
const filesOfDrawings = () => Object.fromEntries(state.drawingList().map((drawing) => [drawing.id, {
  c_gb: drawing.files.c_gb?.name ?? null,
  linearSeqs: drawing.linearSeqs.map((seq) => seq.gb?.name ?? null)
}]));

test('a Session 49 gives each drawing its own files by drawing ID, the ones it does not open too', async () => {
  for (const opened of ['linear', 'linear-2', 'circular']) {
    const session = threeDrawingSession();
    session.ui.activeDrawingId = opened;
    resetDrawings();
    const outcome = await loadKeepingSlots(session);
    assert.equal(outcome.loaded.status, 'ok', JSON.stringify(outcome.loaded.error));
    const files = filesOfDrawings();
    assert.equal(files.circular.c_gb, 'cli-web-mito.gb', opened);
    assert.deepEqual(files.linear.linearSeqs, ['R2c.gb', 'R3c.gb'], opened);
    assert.deepEqual(files['linear-2'].linearSeqs, ['R2c-second.gb', 'R3c.gb'], opened);
  }
});

test('a Session 49 resource that only a drawing it does not open binds keeps its row when it cannot be read', async () => {
  const session = v49Session('drawings-per-drawing-cards');
  session.ui.activeDrawingId = 'circular-2';
  session.resources['resource-0001'] = 5;
  resetDrawings();
  const loaded = await load(session);
  assert.equal(loaded.status, 'ok', JSON.stringify(loaded.error));
  assert.ok(loaded.repairs.some((row) => row.code === 'file-entry-dropped' && row.path === 'resources.resource-0001'),
    JSON.stringify(loaded.repairs));
});

test('a Session 49 whose Linear drawings hold different LOSAT entries keeps each drawing\'s entries', async () => {
  const session = v49Session('drawings-linear-pair-losat');
  const [first, second] = session.drawings;
  const entry = structuredClone(first.losatCache.entries[0]);
  second.losatCache.entries = [{ ...entry, key: `${entry.key}-second`, filename: 'second.tsv' }];
  for (const opened of [first.id, second.id]) {
    session.ui.activeDrawingId = opened;
    resetDrawings();
    const loaded = await load(session);
    assert.equal(loaded.status, 'ok', JSON.stringify(loaded.error));
    assert.deepEqual([...state.losatCache.value.keys()].sort(), [entry.key, `${entry.key}-second`].sort(), opened);
  }
});

test('a Session 49 without a Linear drawing gives the Linear drawing the app adds the Linear defaults', async () => {
  const session = v49Session('drawings-circular-only');
  resetDrawings();
  modeDrawing(state, 'linear').form.plot_title = 'left over';
  modeDrawing(state, 'linear').adv.linear_track_layout = 'below';
  const loaded = await load(session);
  assert.equal(loaded.status, 'ok', JSON.stringify(loaded.error));
  const linear = state.drawingList().find((drawing) => drawing.mode === 'linear');
  assert.ok(linear, 'the app shows a Linear drawing');
  assert.equal(linear.form.plot_title, createDefaultForm().plot_title);
  assert.equal(linear.adv.linear_track_layout, createDefaultAdv('linear').linear_track_layout);
});

test('a Session 49 drawing\'s Web-only artifact fields are dropped as a Session 40-48\'s top-level ones are', async () => {
  const session = v49Session('drawings-circular-linear');
  const index = session.drawings.findIndex((entry) => entry.id !== session.ui.activeDrawingId);
  session.drawings[index].features = { ...session.drawings[index].features, extractedFeatures: [] };
  session.drawings[index].orthogroupState = { ...session.drawings[index].orthogroupState, groups: [{ id: 'g1' }] };
  resetDrawings();
  const loaded = await load(session);
  assert.equal(loaded.status, 'ok', JSON.stringify(loaded.error));
  const paths = loaded.repairs.map((row) => row.path);
  assert.ok(paths.includes(`drawings[${index}].features.extractedFeatures`), JSON.stringify(paths));
  assert.equal(paths.some((path) => path.includes('orthogroupState')), false, JSON.stringify(paths));
});

// B2d-3a-3a: one per-drawing projection. Every Session 49 drawing, opened
// or not, is projected by one function from its own one-drawing view.
test('one function builds the request projection input of every Session 49 drawing, opened or not', () => {
  const source = readFileSync(new URL('../../gbdraw/web/js/services/config.js', import.meta.url), 'utf8');
  // The projection input names the input slot its drawing binds; one builder writes it.
  assert.equal(source.match(/\bbindingInputSlot:/g)?.length ?? 0, 1, 'one projection-input builder');
  assert.match(source, /const drawingProjectionInputOf = /);
  assert.ok((source.match(/\bdrawingProjectionInputOf\(/g)?.length ?? 0) >= 2, 'the opened drawing and every other drawing call it');
});

test('a Session 49 drawing Load does not open keeps its dormant Similarity group names', async () => {
  const session = v49Session('drawings-linear-pair-losat');
  const [first, second] = session.drawings;
  second.config.webEdits.orthogroupDormantOverrides = { '["a","b"]': { name: 'Dormant of linear-2' } };
  first.config.webEdits.orthogroupDormantOverrides = { '["c","d"]': { name: 'Dormant of linear' } };
  for (const opened of [first.id, second.id]) {
    session.ui.activeDrawingId = opened;
    resetDrawings();
    const loaded = await load(session);
    assert.equal(loaded.status, 'ok', JSON.stringify(loaded.error));
    assert.deepEqual({ ...state.drawingOf(second.id).orthogroupDormantOverrides }, second.config.webEdits.orthogroupDormantOverrides, opened);
    assert.deepEqual({ ...state.drawingOf(first.id).orthogroupDormantOverrides }, first.config.webEdits.orthogroupDormantOverrides, opened);
  }
});

// B2d-3a-3a Part 2: session-authority lists a file's Result sets only
// through `storedResultSets` (the one place that names where a Session
// 27-49 file holds a set); its readers key the stored request schemas by
// set path.
test('session-authority reads `otherModeResult` only in the Result set enumeration owner', () => {
  const source = readFileSync(new URL('../../gbdraw/web/js/services/session-authority.js', import.meta.url), 'utf8');
  const bodyOf = (name) => {
    const start = source.indexOf(`export const ${name} = `);
    assert.ok(start >= 0, name);
    return source.slice(start, source.indexOf('\n};\n', start));
  };
  const owners = ['storedResultSets', 'withStoredResultSet'].map(bodyOf);
  const outside = owners.reduce((text, body) => text.replace(body, ''), source)
    .split('\n').filter((line) => !/^\s*(\/\/|\*|\/\*\*)/.test(line))
    // The field registries name the key; they read no set.
    .filter((line) => !/^\s*otherModeResult: \{ authority:|positionalArrays:/.test(line));
  assert.deepEqual(outside.filter((line) => /\botherModeResult\b/.test(line)), []);
});

// OV-1786: Load restores every Session 49 drawing through its own request
// projection, so a drawing gets the same settings whichever drawing opens
// (DESIGN `drawingViewOf`).
test('a Session 49 drawing gets the same settings whichever drawing Load opens (OV-1786)', async () => {
  const session = v49Session('drawings-circular-linear');
  // The Linear drawing has no feature-edit draft: it takes its request's
  // per-feature edits (OV-220), opened or not.
  const linear = session.drawings.find((entry) => entry.mode === 'linear');
  assert.equal(linear.features ?? null, null);
  linear.renderRequest.diagramOptions.featureOverrides = [{ recordKey: 'record-1', biologicalFeatureId: 'f1234~1',
    featureVisibility: null, labelVisibility: null, labelText: 'Request label' }];
  // The drawing's settings (its draft's form and adv), its other settings a
  // request can carry, and its per-feature edits.
  const SETTINGS = ['form', 'adv', 'rules', 'unmanagedConfigOverrides',
    'featureOverrides', 'featureVisibilityManualRules', 'canonicalLabelOverrideRows'];
  const settingsOf = (id) => {
    const drawing = state.drawingOf(id);
    return JSON.parse(JSON.stringify(Object.fromEntries(SETTINGS.map((name) => [name, drawing[name] ?? null]))));
  };
  const loads = [];
  for (const opened of session.drawings.map((entry) => entry.id)) {
    session.ui.activeDrawingId = opened;
    resetDrawings();
    const loaded = await load(session);
    assert.equal(loaded.status, 'ok', JSON.stringify(loaded.error));
    loads.push(Object.fromEntries(session.drawings.map((entry) => [entry.id, settingsOf(entry.id)])));
  }
  const [first, ...rest] = loads;
  const differs = rest.flatMap((other) => session.drawings.flatMap(({ id }) => SETTINGS.flatMap((name) => {
    const [a, b] = [first[id][name] ?? {}, other[id][name] ?? {}];
    return [...new Set([...Object.keys(a), ...Object.keys(b)])]
      .filter((key) => JSON.stringify(a[key]) !== JSON.stringify(b[key])).map((key) => `${id}.${name}.${key}`);
  })));
  assert.deepEqual(differs, []);
});

// DESIGN §6 (B2-D6): after the Session 46-48 upgrade, a 46-48 file and a 49
// file are both a Session 49, so no Web reader of a Session compares its
// version with a threshold of 46 or later (the writer's SESSION_VERSION, a
// band constant, LAYOUT.version, a literal 46-49), lists such versions, or
// holds a hand-written Session version of 46 or later. The band of a file is
// read through the session-authority band predicates, and only the upgrade
// carries a source version (`sourceFacts`). The sites below read a file's
// band before any upgrade runs, each for its reason; every other site fails.
// Thresholds of 44 and below serve the Session 27-44 path until Python
// converts it (B2d-4) and are not counted.
const BAND_READERS = {
  'services/session-authority.js versionHoldsDrawings': 'band predicate: a Session 49 is a list of drawings',
  'services/session-authority.js versionKeepsModeSlices': 'band predicate: a Session 46 or later keeps mode slices',
  'services/session-authority.js holdsCaptionLegendRows': 'band predicate: a Session 46 or older keys Legend rows by caption',
  'services/session-authority.js SETTINGS_ONLY_FORM_VERSIONS': 'band list: the versions with a settings-only form',
  'services/session-authority.js BINDING_SCHEMA_2_VERSIONS': 'band list: the versions whose bindings may be schema 2',
  'services/session-drawings-upgrade.js sessionLegendEditsByRowKey': 'before the upgrade: rekeys a Session 46\'s Legend edits',
  'services/session-drawings-upgrade.js upgradeRowKeyedModeScopedSession': 'the upgrade: admits 46-48, writes sourceFacts',
  'services/svg-result-ingestion.js createSavedResultPlan':
    'a Session 27-44 saved no paint records with its Results (B2d-4: Python supplies them, else one sourceFacts fact)',
  'services/main-session-comparison-frame.js storesDrawnFrameProteinResults':
    'Gallery publication, before any conversion: the stored file\'s protein frame (flip: Gallery reads 49)',
  'services/gallery-session-publication.js validateEnvelope': 'Gallery publication writes the writer\'s version (flip: 49)',
  'services/gallery-session-publication.js admit': 'Gallery publication admits a stored file by its band (flip: 46-48 via the upgrade)'
};
test('no Web reader of a Session compares its version with a Session 46 or later threshold outside the band owner', () => {
  const root = new URL('../../gbdraw/web/js/', import.meta.url);
  const sources = readdirSync(root, { recursive: true }).map(String)
    .filter((name) => name.endsWith('.js') && !name.startsWith('contracts/') && !name.endsWith('.generated.js'))
    .map((name) => ({ name, source: readFileSync(new URL(name, root), 'utf8') }));
  const NAMES = String.raw`MODE_SCOPED_SESSION_VERSION|ROW_KEYED_LEGEND_SESSION_VERSION|FILE_IMPORT_ORIGIN_SESSION_VERSION`
    + String.raw`|DRAWN_FRAME_PROTEIN_RESULT_SESSION_VERSION|SESSION_VERSION|CURRENT_VERSION|(?:SESSION_DRAWINGS_)?LAYOUT\.version`;
  const THRESHOLD = String.raw`(?:${NAMES}|4[6-9])\b`;
  const compare = new RegExp(String.raw`[Vv]ersion\)?\s*(?:[<>]=?|[!=]==)\s*${THRESHOLD}|\b${THRESHOLD}\s*(?:[<>]=?|[!=]==)\s*\S*[Vv]ersion`, 'g');
  // A version list: an array that names a band constant, or an `includes` list of a literal 46-49.
  const list = new RegExp(String.raw`\[[^\][]*\b(?:${NAMES})\b[^\][]*\]|\[[^\][]*\b4[6-9]\b[^\][]*\]\s*\.includes\(`, 'g');
  const handVersion = /\bconst [A-Z_]*SESSION[A-Z_]*VERSION\s*=\s*(4[6-9]|[5-9]\d)\b/g;
  const declaration = /^(\s*)(?:export\s+)?(?:const|let|function|async\s+function)\s+([A-Za-z_$][\w$]*)/;
  /** @param {string} line */
  const indentOf = (line) => line.match(/^\s*/)[0].length;
  // The declaration a line belongs to: its own at the top level, else the
  // nearest one whose block (the lines indented deeper than it) holds it.
  const ownerOf = (lines, index) => {
    const own = lines[index].match(declaration);
    if (own && own[1].length === 0) return own[2];
    for (let at = index - 1; at >= 0; at -= 1) {
      const found = lines[at].match(declaration);
      if (found && lines.slice(at + 1, index + 1).every((line) => !line.trim() || indentOf(line) > found[1].length)) return found[2];
    }
    return '(module)';
  };
  const sites = new Map();
  for (const { name, source } of sources) {
    const lines = source.split('\n');
    const lineOf = (offset) => source.slice(0, offset).split('\n').length - 1;
    for (const pattern of [compare, list, handVersion]) {
      for (const match of source.matchAll(pattern)) {
        const index = lineOf(match.index);
        if (/^\s*(\/\/|\*|\/\*\*)/.test(lines[index])) continue;
        const site = `${name} ${ownerOf(lines, index)}`;
        sites.set(site, [...(sites.get(site) ?? []), `${name}:${index + 1} ${lines[index].trim().slice(0, 120)}`]);
      }
    }
  }
  const unlisted = [...sites].filter(([site]) => !Object.hasOwn(BAND_READERS, site)).flatMap(([, at]) => at);
  assert.deepEqual(unlisted, [], 'a version compare outside the band owner: read a band predicate or a sourceFacts fact');
  assert.deepEqual(Object.keys(BAND_READERS).filter((site) => !sites.has(site)), [], 'a listed site compares no version: remove it');
});

// OV-1791: a Session 49 Load lists one row per fault in a drawing's stored
// unmanaged overrides, at the path the file holds (`drawings[i]`), for the
// opened drawing as for the others.
test('a Session 49 Load lists each drawing\'s stored override fault once, at its drawings[i] path', async () => {
  const session = v49Session('drawings-circular-linear');
  session.ui.activeDrawingId = 'circular';
  const faulted = session.drawings.flatMap((entry, index) => (entry.config ? [{ entry, index }] : []));
  assert.ok(faulted.some(({ entry }) => entry.id === 'circular'));
  faulted.forEach(({ entry }) => { entry.config.unmanagedConfigOverrides = { [`no_such.${entry.id}.key`]: 1 }; });
  const loaded = await importSession({
    target: { files: [new Blob([JSON.stringify(session)], { type: 'application/json' })], value: 'selected' }
  }, { confirmLoadRepairs: async () => true });
  assert.equal(loaded.status, 'ok', JSON.stringify(loaded.error));
  assert.deepEqual(loaded.repairs.filter((repair) => repair.path.includes('unmanagedConfigOverrides'))
    .map((repair) => `${repair.code} ${repair.path}`).sort(), faulted
    .map(({ entry, index }) => `setting-replaced drawings[${index}].config.unmanagedConfigOverrides["no_such.${entry.id}.key"]`).sort());
});

// Save writes Session 49 (B2d-3, DESIGN §4.4): every used drawing in tab order and the shown one.
test('Save writes each drawing of a Session 49 as its own entry; a fresh Load restores the drawings, each Result and the shown one', async () => {
  const session = v49Session('drawings-linear-pair-losat');
  resetDrawings();
  const first = await loadKeepingSlots(session);
  assertLoadedDrawings(session, first);
  const kept = [first.opening, ...first.stashed].filter(Boolean);
  const saved = await (async () => {
    const result = await exportSession(null, { readKeptArtifacts: () => kept });
    assert.equal(result.status, 'saved', result.error?.message);
    return JSON.parse(gunzipSync(Buffer.from(await result.blob.arrayBuffer())));
  })();
  assert.equal(saved.version, 49);
  for (const retired of ['modes', 'otherModeResult', 'renderRequest', 'results', 'config']) assert.equal(Object.hasOwn(saved, retired), false, retired);
  assert.equal(Object.hasOwn(saved.ui, 'mode'), false);
  assert.deepEqual(saved.drawings.map((entry) => entry.id).slice(0, 2), ['linear', 'linear-2']);
  assert.equal(saved.ui.activeDrawingId, 'linear-2');
  saved.drawings.slice(0, 2).forEach((entry) => assert.equal(entry.results.length, 1, entry.id));
  resetDrawings();
  assertLoadedDrawings(saved, await loadKeepingSlots(saved));
});

test('Save writes each Session 49 drawing\'s own CLI options, found by drawing ID, the drawing Load does not open too', async () => {
  const session = v49Session('drawings-linear-pair-losat');
  const optionsOf = (id) => ({ rawArgs: ['--gbk', `${id}.gb`, '-o', id] });
  session.drawings.forEach((entry) => { entry.cliOptions = optionsOf(entry.id); });
  resetDrawings();
  const first = await loadKeepingSlots(session);
  assertLoadedDrawings(session, first);
  const result = await exportSession('cli-options', { readKeptArtifacts: () => [first.opening, ...first.stashed].filter(Boolean) });
  assert.equal(result.status, 'saved', result.error?.message);
  const saved = JSON.parse(gunzipSync(Buffer.from(await result.blob.arrayBuffer())));
  for (const id of ['linear', 'linear-2']) {
    assert.deepEqual(saved.drawings.find((entry) => entry.id === id).cliOptions, optionsOf(id), id);
  }
});

test('Save checks each Session 49 drawing it writes with Save\'s checks, not Load\'s: an input type of no type is refused', async () => {
  const session = v49Session('drawings-linear-pair-losat');
  resetDrawings();
  const first = await loadKeepingSlots(session);
  assertLoadedDrawings(session, first);
  // Load reads such a type as GenBank with a row (OV-594); Save refuses to write one.
  state.activeDrawing().inputType.value = 'embl';
  await assert.rejects(
    exportSession('input-type', { readKeptArtifacts: () => [first.opening, ...first.stashed].filter(Boolean) }),
    /INPUT_INVALID/
  );
});

test('Save writes into each Session 49 drawing the LOSAT entries of its own Linear rows (K-3)', async () => {
  const session = v49Session('drawings-linear-pair-losat');
  const first = session.drawings.find((entry) => entry.id === 'linear');
  const entry = first.losatCache.entries[0];
  // linear-2 holds other rows: each uid of linear's renamed throughout it, with its own entry.
  const renamed = (text) => text.replaceAll(entry.queryUid, `${entry.queryUid}-2`).replaceAll(entry.subjectUid, `${entry.subjectUid}-2`);
  const index = session.drawings.findIndex((candidate) => candidate.id === 'linear-2');
  session.drawings[index] = JSON.parse(renamed(JSON.stringify(session.drawings[index])));
  session.drawings[index].losatCache.entries = [JSON.parse(renamed(JSON.stringify({ ...entry, key: `${entry.key}-second`, filename: 'second.tsv' })))];
  resetDrawings();
  const loaded = await loadKeepingSlots(session);
  assertLoadedDrawings(session, loaded);
  const result = await exportSession('losat-entries', { readKeptArtifacts: () => [loaded.opening, ...loaded.stashed].filter(Boolean) });
  assert.equal(result.status, 'saved', result.error?.message);
  const saved = JSON.parse(gunzipSync(Buffer.from(await result.blob.arrayBuffer())));
  const keysOf = (id) => saved.drawings.find((candidate) => candidate.id === id).losatCache.entries.map(({ key }) => key);
  assert.deepEqual(keysOf('linear'), [entry.key]);
  assert.deepEqual(keysOf('linear-2'), [`${entry.key}-second`]);
});

// OV-337: Load leaves a Linear row without a Depth file as the Session binds it
// (`depth: null`), so Save writes it back unchanged. Load pads Depth rows to
// the Depth track count only when a row holds a Depth file, as for Circular.
test('Load and Save keep a Linear row without a Depth file as depth null (OV-337)', async () => {
  const session = v49Session('drawings-linear-pair-losat');
  const saveKept = async (title, outcome) => {
    assert.equal(outcome.loaded.status, 'ok', JSON.stringify(outcome.loaded.error));
    const result = await exportSession(title, { readKeptArtifacts: () => [outcome.opening, ...outcome.stashed].filter(Boolean) });
    assert.equal(result.status, 'saved', result.error?.message);
    return JSON.parse(gunzipSync(Buffer.from(await result.blob.arrayBuffer())));
  };
  resetDrawings();
  const first = await saveKept('round trip', await loadKeepingSlots(session));
  const linearRows = (document) => savedDrawingView(document, 'linear').webFiles.bindings.linearSeqs;
  const rows = linearRows(first);
  assert.deepEqual(rows.map((row) => row.depth), [null, null]);
  assert.equal(savedDrawingView(first, 'linear').modes.linear.config.adv.depth_tracks.length, 1, 'the Linear drawing has a Depth series');
  const reloaded = await loadKeepingSlots(first);
  assert.deepEqual(state.drawingOf('linear').linearSeqs.map((seq) => seq.depth), [null, null]);
  const second = await saveKept('round trip again', reloaded);
  assert.deepEqual(linearRows(second).map((row) => row.depth), rows.map((row) => row.depth));
});

// OV-1811: the drawing Load adds for a mode the Session lacks (T-1) is a new drawing at its mode's
// defaults (D-15-120), so Save does not count it as used.
test('Save does not write the unused Circular drawing Load adds for the mode buttons (R16-8, OV-1811)', () => withBootPalettes(state, async () => {
  resetDrawings();
  const first = await loadKeepingSlots(v49Session('drawings-linear-pair-losat'));
  const result = await exportSession('ov-1811', { readKeptArtifacts: () => [first.opening, ...first.stashed].filter(Boolean) });
  assert.equal(result.status, 'saved', result.error?.message);
  const saved = JSON.parse(gunzipSync(Buffer.from(await result.blob.arrayBuffer())));
  assert.deepEqual(saved.drawings.map((entry) => entry.id), ['linear', 'linear-2']);
}));

// OV-1811, a Session 48 the move gives one drawing (its Linear part unused): the Linear drawing Load
// adds takes nothing of that part, so Save writes the Circular drawing alone.
test('Save does not write the unused Linear drawing Load adds after a Session 48 Load (R16-8, OV-1811)', () => withBootPalettes(state, async () => {
  assert.deepEqual(drawingsVector('file-card-imports-circular.v48').drawings.map((drawing) => drawing.mode), ['circular']);
  resetDrawings();
  const loaded = await load(fixtureSession('file-card-imports-circular.v48'));
  assert.equal(loaded.status, 'ok', JSON.stringify(loaded.error));
  assert.deepEqual(state.drawingList().map((drawing) => drawing.id), ['circular', 'linear']);
  const saved = await save('ov-1811 48');
  assert.deepEqual(saved.drawings.map((entry) => entry.id), ['circular']);
}));
