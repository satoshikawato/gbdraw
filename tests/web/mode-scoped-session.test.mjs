// Session 46 (plan OV-80 §4, §8 S1, S3, S4, S6, S7): each diagram mode's
// drawing is saved as its own slice (`modes.<mode>`) and loaded back into its
// own drawing; a slice holds registry fields only; Session 45 and a Session 46
// with a retired field are rejected.
import assert from 'node:assert/strict';
import { test } from 'node:test';
import { readFileSync, readdirSync } from 'node:fs';
import { gunzipSync } from 'node:zlib';
import { gallerySessionDocument, readSavedSession, sessionLoadEvent } from './helpers/session-documents.mjs';
import { installSessionImportWorker } from './helpers/session-import-node.mjs';
import { withBootPalettes } from './helpers/boot-palettes.mjs';
import { buildPaletteColorOverrideRows } from '../../gbdraw/web/js/utils/color-utils.js';
import { modeDrawing, withModeDrawings } from './helpers/drawing-state.mjs';
import { callerSession } from './helpers/session-caller.mjs';
import { SESSION_DRAWINGS_LAYOUT } from '../../gbdraw/web/js/contracts/schema.generated.js';
import { drawingViewAt, savedDrawingView, withSavedDrawingView } from './helpers/saved-drawing-view.mjs';
import { sessionRollbackPort } from './helpers/session-rollback-port.mjs';
import { evaluatePython, validateWebConfigOverrides, pythonLegacySessionUpgrader } from './helpers/python-evaluator.mjs';
import { createHash } from 'node:crypto';
import {
  emptyLosatEvidence, losatRowsOnePerKey, validateProteinRawEntryReferences
} from '../../gbdraw/web/js/services/losat-cache.js';
import { withArtifactSlotValues } from '../../gbdraw/web/js/services/artifact-slot.js';

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
// A Session 27-44 moves to a Session 51 in Python, as the composition root's
// Worker port does it (DESIGN 4.3).
const { setLegacySessionUpgrader } = await import('../../gbdraw/web/js/services/legacy-session-upgrade.js');
setLegacySessionUpgrader(pythonLegacySessionUpgrader());

const sessionConfig = await import('../../gbdraw/web/js/services/config.js');
const { SESSION_VERSION, assertDrawingInputs, installShownDrawing, getCommittedCanonicalSession } = sessionConfig;
// A Session 51 drawing's saved unmanaged overrides are checked as Load checks them (PD-OI-009).
sessionConfig.setUnmanagedConfigOverrideValidator(async (payload) => ({ result: await validateWebConfigOverrides(payload) }));
const { MODE_SLICE_CONTAINERS } = await import('../../gbdraw/web/js/services/mode-scoped-migration.js');
const { createDefaultAdv, createDefaultForm } = await import('../../gbdraw/web/js/services/session-active-config-contract.js');
const { state } = await import('../../gbdraw/web/js/state.js');
const { getSvgDimensions } = await import('../../gbdraw/web/js/services/svg-dimensions.js');
const { requestFigureSizeLeaves, shownFigureSize } = await import('../../gbdraw/web/js/services/figure-size-form.js');
await (await import('./helpers/drawing-state.mjs')).addStartLinearDrawing();
const { exportSession, importSession } = callerSession(sessionConfig, state, await sessionRollbackPort(state));

const MODES = ['circular', 'linear'];
const save = async (title) => {
  const result = await exportSession(title);
  assert.equal(result.status, 'saved', result.error?.message);
  return await readSavedSession(result.blob);
};
// The Load dialogs asked (BA-D25); the port answers [Load].
const loadDialogs = [];
const load = async (session) => importSession(await sessionLoadEvent(await session), { confirmLoadRepairs: async (dialog) => { loadDialogs.push(dialog); return true; } });
// The Circular and Linear pair a case works on (a one-mode Load makes no drawing of the other mode).
const addPairDrawings = await (await import('./helpers/drawing-state.mjs')).modeDrawingsAdder();
const resetDrawings = () => {
  addPairDrawings();
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

test('Session 51 saves each drawing as its own entry and loads it back into that drawing (S1, S6)', async () => {
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

  // A slice may omit any field, which Load fills with its mode's defaults, and
  // a Session may hold no drawing of a mode: Load makes none (D-15-120). (A
  // settings-only Session keeps the slice of its shown mode.)
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
  assert.deepEqual(state.drawingList().map((drawing) => drawing.id), ['circular']);
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

// A saved drawing's one-drawing view (a Session 51: its Circular input in `circularSources`) as the Session 46-50
// a writer of that version saved: only its fields (`drawingViewAt`), and the input in the `c_gb` family.
const olderView = (view, version) => {
  const session = drawingViewAt(view, version);
  const bindings = session.webFiles?.bindings;
  if (bindings && Object.hasOwn(bindings, 'circularSources')) {
    const [card] = bindings.circularSources;
    delete bindings.circularSources;
    Object.assign(bindings, { c_gb: card?.gb ?? null, c_gff: card?.gff ?? null, c_fasta: card?.fasta ?? null });
  }
  return session;
};

// S3, S4 (E-OD03, A7d, A7f): Load refuses a version it does not read; it drops
// a Session 46's retired homes, whose place the slices fill, with no row, and a
// key no writer wrote with one.
test('Session 45 is refused; a Session 46 drops its retired homes with no row and an unknown key with one (S3, S4)', async () => {
  resetDrawings();
  // The one-drawing view of a saved drawing is the Session 46 shape this test loads.
  const current = await save('current');
  const viewOf = (id) => olderView(savedDrawingView(current, id), 46);
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

// A Session reader admits only what its format had (CLAUDE.md "Persisted-format compatibility"): the Figure size
// and the window across the origin are fields Session 51 adds to a drawing's form, so in a Session 46-50 slice
// they are fields no writer of that version wrote, dropped with a row each like any other.
test('a Session 50 slice drops the form fields only Session 51 writes, with a row each', async () => {
  resetDrawings();
  const current = await save('session 50 fields');
  const saved = olderView(savedDrawingView(current, 'circular'), 50);
  const form = saved.modes.circular.config.form;
  const rows = async (session) => {
    const outcome = await load(session);
    assert.equal(outcome.status, 'ok', JSON.stringify(outcome.error));
    return (outcome.repairs ?? []).map(({ code, path }) => `${code}@${path}`).sort();
  };
  assert.deepEqual(await rows(structuredClone(saved)), []);
  Object.assign(form, { figure_width: '180mm', figure_height: null,
    circular_windows: [{ start: 90, end: 10, reverseComplement: false, column: null }] });
  assert.deepEqual(await rows(structuredClone(saved)), ['circular_windows', 'figure_height', 'figure_width']
    .map((field) => `setting-replaced@modes.circular.config.form.${field}`));
  const loaded = modeDrawing(state, 'circular').form;
  assert.deepEqual([loaded.figure_width, loaded.figure_height, loaded.circular_windows],
    [createDefaultForm().figure_width, createDefaultForm().figure_height, []]);
  // A Session 51 drawing holds them.
  assert.deepEqual(await rows(withSavedDrawingView(current, 'circular', (view) => {
    view.modes.circular.config.form.figure_width = '180mm';
  })), []);
  assert.equal(modeDrawing(state, 'circular').form.figure_width, '180mm');
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
  const older = olderView(savedDrawingView(saved, 'circular'), 46);
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
  const gallery = await gallerySessionDocument('HmmtDNA_basic_circular');
  state.labelFitReports.value = [{ resultIndex: 0, resultName: 'fit.svg', rowToken: 'circular:1:0',
    report: { drawn: true, notice: ['Leader lines cross (5 pairs).'] } }];
  state.inputSummaries.value = [{ resultIndex: 0, resultName: 'fit.svg', lines: ['Record r (a.gb), record 1 of 1, 8 bp: drew 0 features.'],
    warnings: [{ code: 'no_drawable_features', message: 'Record r (a.gb) has no feature to draw: it has no features.' }] }];
  assert.equal((await load(gallery)).status, 'ok');
  assert.ok(state.results.value.length > 0, 'the Session has a Result');
  assert.deepEqual(state.labelFitReports.value, []);
  assert.deepEqual(state.inputSummaries.value, [], 'a loaded Result has no input summary');
});

// Record selection D-05: each slice keeps the keys of its OFF records; an
// omitted list draws every record. Load drops a key that names no record of
// its drawing without a row (BA-OD88, E-OD13), reads a repeat once, and lists
// a list of the wrong shape.
test('a slice keeps its OFF records, and omission draws every record (record-selection D-05)', async () => {
  resetDrawings();
  installShownDrawing('circular');
  const gallery = await gallerySessionDocument('HmmtDNA_basic_circular');
  assert.equal((await load(gallery)).status, 'ok');
  addPairDrawings();
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
  // A record key has one shape in both modes (a card uid or a record identity), so a Circular key is no longer
  // required to be `#N`; an entry that is no key (padded text) is dropped.
  await loaded('circular', (config) => { config.recordsOff = ['contig_2']; }, 'Circular key of the key shape',
    { linear: [uid], circular: ['contig_2'] });
  await loaded('circular', (config) => { config.recordsOff = [' contig_2 ']; }, 'Circular entry that is no key',
    { linear: [uid], circular: [], rows: [`setting-entries-dropped@drawings[${saved.drawings.findIndex(({ id }) => id === 'circular')}].config.recordsOff`] });
  // Review F8: a list of the wrong shape is listed; every record is drawn.
  await loaded('circular', (config) => { config.recordsOff = '#2'; }, 'not a list',
    { linear: [uid], circular: [], rows: [`setting-replaced@drawings[${saved.drawings.findIndex(({ id }) => id === 'circular')}].config.recordsOff`] });
});

// Review #992 finding 1: a Session without a slice of a mode (a Gallery
// Session of one mode) draws every record of that mode too, and saves again.
test('a Session without a mode slice turns every record of that mode ON (record selection)', async () => {
  resetDrawings();
  const gallery = (name) => gallerySessionDocument(name);
  assert.equal((await load(gallery('lambda_basic_linear'))).status, 'ok');
  modeDrawing(state, 'linear').recordsOff.splice(0, Infinity, modeDrawing(state, 'linear').linearSeqs[0].uid);
  assert.equal((await load(gallery('HmmtDNA_basic_circular'))).status, 'ok');
  assert.deepEqual(state.drawingList().map((drawing) => drawing.mode), ['circular']);
  assert.deepEqual([...modeDrawing(state, 'circular').recordsOff], []);
  await save('circular only, saved again');
});

// Review #992 finding 6: Reset Settings keeps which records are ON in both
// drawings (PD-OI-070, PD-OI-091).
test('Reset Settings keeps the OFF records of both drawings (record selection)', async () => {
  const { resetSettings } = await import('../../gbdraw/web/js/services/reset.js');
  addPairDrawings();
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
  addPairDrawings();
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
  assert.deepEqual(['d_color', 't_color', 'qualifier_priority', 'whitelist', 'blacklist'].map((key) => Boolean(state.drawingOf('circular').files[key])),
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

// B2d-2: Load builds one drawing per drawing of the Session (its 51 document)
// and commits them at once; the fast path reads 46-50 without Python.
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
  const loaded = await importSession(await sessionLoadEvent(fixtureSession('two-mode-linear-shown.v48')), {
    confirmLoadRepairs: async () => true,
    afterLoad: async () => { throw new Error('install failed'); }
  });
  assert.equal(loaded.status, 'error');
  assert.equal(state.drawingList().length, before.length);
  state.drawingList().forEach((drawing, index) => assert.equal(drawing, before[index], drawing.id));
  assert.equal(state.activeDrawing().id, 'linear');
  assert.equal(modeDrawing(state, 'linear').adv.label_font_size, 11);
});

test('Load reads each drawing\'s saved Result set from its own entry of the 51 document, not from `otherModeResult`', () => {
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
  // A Session 46-50 drawing's settings, per-feature edits, editor edits and
  // `ui` come from its entry of the 51 document; a Session 27-44 split is
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
  const loaded = await importSession(await sessionLoadEvent(session), {
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

// Session 51 (D-15-109): the Web reads the drawings document Python writes.
// The version-51 fixtures are Session file containers, read once in their JSON form.
const V51_FIXTURES = Object.fromEntries(await Promise.all(readdirSync(new URL('../fixtures/sessions/', import.meta.url))
  .filter((name) => /\.v51[^/]*\.gbdraw$/.test(name))
  .map(async (name) => [name, await readSavedSession(new Blob([readFileSync(new URL(`../fixtures/sessions/${name}`, import.meta.url))]))])));
const v51File = (file) => structuredClone(V51_FIXTURES[file]);
const v49Session = (name) => v51File(`${name}.v51.gbdraw`);
const loadKeepingSlots = async (session) => {
  const stashed = [];
  let opening = null;
  const loaded = await importSession(await sessionLoadEvent(session), {
    confirmLoadRepairs: async () => true,
    installLoadedArtifactSlots: (slots) => { opening = slots.opening; stashed.push(...slots.stashed); }
  });
  return { loaded, stashed, opening };
};
// Each drawing's LOSAT evidence after a Load (D-221): the shown one's from state, any other's from its slot.
const loadedLosatEvidenceOf = (id, { stashed, opening }) => {
  if (id === state.activeDrawing().id && !opening) {
    return { losatCacheInfo: state.losatCacheInfo.value, losatCache: state.losatCache.value };
  }
  return [opening, ...stashed].find((entry) => entry?.drawingId === id)?.values ?? emptyLosatEvidence();
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

test('a Python-written Session 51 loads each drawing\'s request, Results and draft by drawing ID', async () => {
  for (const name of ['drawings-circular-only', 'drawings-circular-linear', 'drawings-linear-pair-losat', 'drawings-per-drawing-cards']) {
    const session = v49Session(name);
    // A draft value of its own in each drawing that holds a draft.
    session.drawings.forEach((entry) => { if (entry.config?.form) entry.config.form.plot_title = `title of ${entry.id}`; });
    resetDrawings();
    assertLoadedDrawings(session, await loadKeepingSlots(session));
  }
});

test('a Session 51 with two Circular drawings and a Linear one loads all three in tab order and keeps each Result in its own slot', async () => {
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

test('a Session 51 without its drawings or its resources is refused as no Session (OV-1768)', async () => {
  for (const field of ['drawings', 'resources']) {
    const session = v49Session('drawings-circular-only');
    delete session[field];
    resetDrawings();
    const outcome = await load(session);
    assert.equal(outcome.status, 'error', field);
    assert.equal(outcome.error.context?.reason, 'SESSION_FORMAT', `${field}: ${JSON.stringify(outcome.error)}`);
  }
});

// B2d-3a-2: each drawing of a Session 51 takes its own files, bindings, LOSAT
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
  circularSources: drawing.circularSources.map((source) => source.gb?.name ?? null),
  linearSeqs: drawing.linearSeqs.map((seq) => seq.gb?.name ?? null)
}]));

test('a Session 51 gives each drawing its own files by drawing ID, the ones it does not open too', async () => {
  for (const opened of ['linear', 'linear-2', 'circular']) {
    const session = threeDrawingSession();
    session.ui.activeDrawingId = opened;
    resetDrawings();
    const outcome = await loadKeepingSlots(session);
    assert.equal(outcome.loaded.status, 'ok', JSON.stringify(outcome.loaded.error));
    const files = filesOfDrawings();
    assert.deepEqual(files.circular.circularSources, ['cli-web-mito.gb'], opened);
    assert.deepEqual(files.linear.linearSeqs, ['R2c.gb', 'R3c.gb'], opened);
    assert.deepEqual(files['linear-2'].linearSeqs, ['R2c-second.gb', 'R3c.gb'], opened);
  }
});

// A container holds only descriptors (the reader refuses anything else as damaged), so a resource that
// cannot be read reaches the Load only in the JSON file of a Session 27-50.
test('a Session 50 resource that its binding names keeps its row when it cannot be read', async () => {
  const session = JSON.parse(gunzipSync(readFileSync(new URL('../fixtures/sessions/cli-linear-wrapped.v50.gbdraw-session.json.gz', import.meta.url))));
  // Only a binding names it: the request draws from the other resource.
  session.resources['colors-0001'] = 5;
  session.webFiles.bindings.d_color = { resourceId: 'colors-0001', name: 'colors.tsv', type: 'text/tab-separated-values', lastModified: 0 };
  resetDrawings();
  const loaded = await load(session);
  assert.equal(loaded.status, 'ok', JSON.stringify(loaded.error));
  assert.ok(loaded.repairs.some((row) => row.code === 'file-entry-dropped' && row.path === 'resources.colors-0001'),
    JSON.stringify(loaded.repairs));
});

test('a Session 51 whose Linear drawings hold different LOSAT entries keeps each drawing\'s entries', async () => {
  const session = v49Session('drawings-linear-pair-losat');
  const [first, second] = session.drawings;
  const entry = structuredClone(first.losatCache.entries[0]);
  second.losatCache.entries = [{ ...entry, key: `${entry.key}-second`, filename: 'second.tsv' }];
  for (const opened of [first.id, second.id]) {
    session.ui.activeDrawingId = opened;
    resetDrawings();
    const outcome = await loadKeepingSlots(session);
    assert.equal(outcome.loaded.status, 'ok', JSON.stringify(outcome.loaded.error));
    // Each drawing holds its own entries, whichever is shown (D-221).
    assert.deepEqual([first.id, second.id].map((id) => [...loadedLosatEvidenceOf(id, outcome).losatCache.keys()]),
      [[entry.key], [`${entry.key}-second`]], opened);
  }
});

// #1777 round 14 (OV-3210): two Linear drawings that hold one pair's entry get a row each, and each drawing's
// per-pair readers read its own row only: the pair once, with that drawing's row.
test('two Linear drawings that hold one LOSAT pair each read the pair once, from their own row (OV-3210)', async () => {
  const { createLinearLayoutReaders } = await import('../../gbdraw/web/js/app/linear-comparison-panel.js');
  const session = v49Session('drawings-linear-pair-losat');
  const [first, second] = session.drawings;
  second.losatCache.entries = structuredClone(first.losatCache.entries);
  resetDrawings();
  const outcome = await loadKeepingSlots(session);
  assert.equal(outcome.loaded.status, 'ok', JSON.stringify(outcome.loaded.error));
  const [{ edgeKey }] = first.losatCache.entries;
  const rowsOf = [first, second].map((drawing) => loadedLosatEvidenceOf(drawing.id, outcome).losatCacheInfo);
  assert.notEqual(rowsOf[0], rowsOf[1], 'each drawing has its own rows (D-221)');
  for (const [index, drawing] of [first, second].entries()) {
    const { linearLosatCacheInfoByEdgeKey } = createLinearLayoutReaders({
      state, shownDrawingOfMode: () => state.drawingOf(drawing.id), losatCacheInfo: { value: rowsOf[index] }
    });
    assert.deepEqual(Object.values(linearLosatCacheInfoByEdgeKey.value).map((row) => row.edgeKey), [edgeKey], drawing.id);
  }
});

// D-15-120: a Session 51 without a Linear drawing loads without one; the Linear drawing the app
// held before the Load goes with the drawings the Load replaces.
test('a Session 51 without a Linear drawing loads no Linear drawing', async () => {
  const session = v49Session('drawings-circular-only');
  resetDrawings();
  modeDrawing(state, 'linear').form.plot_title = 'left over';
  const loaded = await load(session);
  assert.equal(loaded.status, 'ok', JSON.stringify(loaded.error));
  assert.deepEqual(state.drawingList().map((drawing) => drawing.mode), session.drawings.map((entry) => entry.mode));
  assert.equal(state.drawingList().some((drawing) => drawing.mode === 'linear'), false);
});

test('a Session 51 drawing\'s Web-only artifact fields are dropped as a Session 40-50\'s top-level ones are', async () => {
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

// B2d-3a-3a: one per-drawing projection. Every Session 51 drawing, opened
// or not, is projected by one function from its own one-drawing view.
test('one function builds the request projection input of every Session 51 drawing, opened or not', () => {
  const source = readFileSync(new URL('../../gbdraw/web/js/services/config.js', import.meta.url), 'utf8');
  // The projection input names the input slot its drawing binds; one builder writes it.
  assert.equal(source.match(/\bbindingInputSlot:/g)?.length ?? 0, 1, 'one projection-input builder');
  assert.match(source, /const drawingProjectionInputOf = /);
  assert.ok((source.match(/\bdrawingProjectionInputOf\(/g)?.length ?? 0) >= 2, 'the opened drawing and every other drawing call it');
});

test('a Session 51 drawing Load does not open keeps its dormant Similarity group names', async () => {
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
// 27-50 file holds a set); its readers key the stored request schemas by
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

// OV-1786: Load restores every Session 51 drawing through its own request
// projection, so a drawing gets the same settings whichever drawing opens
// (DESIGN `drawingViewOf`).
test('a Session 51 drawing gets the same settings whichever drawing Load opens (OV-1786)', async () => {
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

// DESIGN §6 (B2-D6): after the Session 46-50 upgrade, a 46-50 file and a 51
// file are both a Session 51, so no Web reader of a Session compares its
// version with a threshold of 46 or later (the writer's SESSION_VERSION, a
// band constant, LAYOUT.version, a literal 46-50), lists such versions, or
// holds a hand-written Session version of 46 or later. The band of a file is
// read through the session-authority band predicates, and only the upgrade
// carries a source version (`sourceFacts`). The sites below read a file's
// band before any upgrade runs, each for its reason; every other site fails.
// Thresholds of 44 and below serve the Session 27-44 path until Python
// converts it (B2d-4) and are not counted.
const BAND_READERS = {
  'services/session-authority.js versionHoldsDrawings': 'band predicate: a Session 51 is a list of drawings',
  'services/session-authority.js versionKeepsModeSlices': 'band predicate: a Session 46 or later keeps mode slices',
  'services/session-authority.js holdsCaptionLegendRows': 'band predicate: a Session 46 or older keys Legend rows by caption',
  'services/session-drawings-upgrade.js sessionLegendEditsByRowKey': 'before the upgrade: rekeys a Session 46\'s Legend edits',
  'services/session-drawings-upgrade.js upgradeRowKeyedModeScopedSession': 'the upgrade: admits 46-50, writes sourceFacts',
  'services/session-drawings-upgrade.js sessionSourceFactsOf':
    'the source facts by version (one owner, shared by the 46-50 upgrade and the 27-44 Load)',
  'services/main-session-comparison-frame.js storesDrawnFrameProteinResults':
    'Gallery publication, before any conversion: the stored file\'s protein frame (flip: Gallery reads 50)',
  'services/gallery-session-publication.js validateEnvelope': 'Gallery publication writes the writer\'s version (flip: 50)',
  'services/gallery-session-publication.js admit': 'Gallery publication admits a stored file by its band (flip: 46-50 via the upgrade)'
};
test('no Web reader of a Session compares its version with a Session 46 or later threshold outside the band owner', () => {
  const root = new URL('../../gbdraw/web/js/', import.meta.url);
  const sources = readdirSync(root, { recursive: true }).map(String)
    .filter((name) => name.endsWith('.js') && !name.startsWith('contracts/') && !name.endsWith('.generated.js'))
    .map((name) => ({ name, source: readFileSync(new URL(name, root), 'utf8') }));
  const NAMES = String.raw`MODE_SCOPED_SESSION_VERSION|ROW_KEYED_LEGEND_SESSION_VERSION|FILE_IMPORT_ORIGIN_SESSION_VERSION`
    + String.raw`|DRAWN_FRAME_PROTEIN_RESULT_SESSION_VERSION|SESSION_VERSION|CURRENT_VERSION|(?:SESSION_DRAWINGS_)?LAYOUT\.version`;
  const THRESHOLD = String.raw`(?:${NAMES}|4[6-9]|50)\b`;
  const compare = new RegExp(String.raw`[Vv]ersion\)?\s*(?:[<>]=?|[!=]==)\s*${THRESHOLD}|\b${THRESHOLD}\s*(?:[<>]=?|[!=]==)\s*\S*[Vv]ersion`, 'g');
  // A version list: an array that names a band constant, or an `includes` list of a literal 46-50.
  const list = new RegExp(String.raw`\[[^\][]*\b(?:${NAMES})\b[^\][]*\]|\[[^\][]*\b(?:4[6-9]|50)\b[^\][]*\]\s*\.includes\(`, 'g');
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

// OV-1791: a Session 51 Load lists one row per fault in a drawing's stored
// unmanaged overrides, at the path the file holds (`drawings[i]`), for the
// opened drawing as for the others.
test('a Session 51 Load lists each drawing\'s stored override fault once, at its drawings[i] path', async () => {
  const session = v49Session('drawings-circular-linear');
  session.ui.activeDrawingId = 'circular';
  const faulted = session.drawings.flatMap((entry, index) => (entry.config ? [{ entry, index }] : []));
  assert.ok(faulted.some(({ entry }) => entry.id === 'circular'));
  faulted.forEach(({ entry }) => { entry.config.unmanagedConfigOverrides = { [`no_such.${entry.id}.key`]: 1 }; });
  const loaded = await importSession(await sessionLoadEvent(session), { confirmLoadRepairs: async () => true });
  assert.equal(loaded.status, 'ok', JSON.stringify(loaded.error));
  assert.deepEqual(loaded.repairs.filter((repair) => repair.path.includes('unmanagedConfigOverrides'))
    .map((repair) => `${repair.code} ${repair.path}`).sort(), faulted
    .map(({ entry, index }) => `setting-replaced drawings[${index}].config.unmanagedConfigOverrides["no_such.${entry.id}.key"]`).sort());
});

// Save writes Session 51 (B2d-3, DESIGN §4.4): every used drawing in tab order and the shown one.
test('Save writes each drawing of a Session 51 as its own entry; a fresh Load restores the drawings, each Result and the shown one', async () => {
  const session = v49Session('drawings-linear-pair-losat');
  resetDrawings();
  const first = await loadKeepingSlots(session);
  assertLoadedDrawings(session, first);
  const kept = [first.opening, ...first.stashed].filter(Boolean);
  const saved = await (async () => {
    const result = await exportSession(null, { readKeptArtifacts: () => kept });
    assert.equal(result.status, 'saved', result.error?.message);
    return await readSavedSession(result.blob);
  })();
  assert.equal(saved.version, SESSION_VERSION);
  for (const retired of ['modes', 'otherModeResult', 'renderRequest', 'results', 'config']) assert.equal(Object.hasOwn(saved, retired), false, retired);
  assert.equal(Object.hasOwn(saved.ui, 'mode'), false);
  assert.deepEqual(saved.drawings.map((entry) => entry.id).slice(0, 2), ['linear', 'linear-2']);
  assert.equal(saved.ui.activeDrawingId, 'linear-2');
  saved.drawings.slice(0, 2).forEach((entry) => assert.equal(entry.results.length, 1, entry.id));
  resetDrawings();
  assertLoadedDrawings(saved, await loadKeepingSlots(saved));
});

test('Save writes each Session 51 drawing\'s own CLI options, found by drawing ID, the drawing Load does not open too', async () => {
  const session = v49Session('drawings-linear-pair-losat');
  const optionsOf = (id) => ({ rawArgs: ['--gbk', `${id}.gb`, '-o', id] });
  session.drawings.forEach((entry) => { entry.cliOptions = optionsOf(entry.id); });
  resetDrawings();
  const first = await loadKeepingSlots(session);
  assertLoadedDrawings(session, first);
  const result = await exportSession('cli-options', { readKeptArtifacts: () => [first.opening, ...first.stashed].filter(Boolean) });
  assert.equal(result.status, 'saved', result.error?.message);
  const saved = await readSavedSession(result.blob);
  for (const id of ['linear', 'linear-2']) {
    assert.deepEqual(saved.drawings.find((entry) => entry.id === id).cliOptions, optionsOf(id), id);
  }
});

// A wrap length is one per Linear drawing (Session 51): Wrap Rows, the
// drawing's `form.linear_wrap_length` (PD-OI-113). Load reads each drawing's
// own value, the unopened one's too, found by drawing ID; a length a 0.14.5
// Web Save kept as a preserved override moves into the form of the drawing
// that holds it. Save writes it in that drawing's form, and no preserved
// override holds it.
test('a Session 51 keeps each Linear drawing\'s own wrap length override through Load and Save', async () => {
  const WRAP = 'canvas.linear.wrap_length';
  const session = v49Session('drawings-linear-pair-losat');
  const entryOf = (document, id) => document.drawings.find((entry) => entry.id === id);
  entryOf(session, 'linear-2').config.unmanagedConfigOverrides = { [WRAP]: 30000 };
  entryOf(session, 'linear').config.form.linear_wrap_length = 'auto';
  resetDrawings();
  const first = await loadKeepingSlots(session);
  assertLoadedDrawings(session, first);
  assert.equal(state.drawingOf('linear-2').form.linear_wrap_length, 30000);
  assert.equal(state.drawingOf('linear').form.linear_wrap_length, 'auto');
  for (const id of ['linear', 'linear-2']) assert.deepEqual({ ...state.drawingOf(id).unmanagedConfigOverrides }, {}, id);
  const result = await exportSession('wrap-length', { readKeptArtifacts: () => [first.opening, ...first.stashed].filter(Boolean) });
  assert.equal(result.status, 'saved', result.error?.message);
  const saved = await readSavedSession(result.blob);
  assert.equal(entryOf(saved, 'linear-2').config.form.linear_wrap_length, 30000);
  assert.equal(entryOf(saved, 'linear').config.form.linear_wrap_length, 'auto');
  for (const id of ['linear', 'linear-2']) assert.deepEqual(entryOf(saved, id).config.unmanagedConfigOverrides ?? {}, {}, id);
  // The saved form loads again into each drawing.
  resetDrawings();
  const again = await loadKeepingSlots(saved);
  assertLoadedDrawings(saved, again);
  assert.equal(state.drawingOf('linear-2').form.linear_wrap_length, 30000);
  assert.equal(state.drawingOf('linear').form.linear_wrap_length, 'auto');
});

// D-36: a stored wrap length the form cannot hold as stored is Python's to read (the config model's validator, at
// Load's configuration check): a length typed with a unit, another case, spaces. A value Python refuses is dropped
// with one Load repair row; the drawing's Wrap Rows is off.
test('a preserved wrap length takes Python\'s reading, and one Python refuses is dropped with a Load repair row', async () => {
  const WRAP = 'canvas.linear.wrap_length';
  for (const [stored, read] of [['50kb', 50000], ['Auto', 'auto'], [' 20000 ', 20000], ['None', null], ['abc', null]]) {
    const session = v49Session('drawings-linear-pair-losat');
    session.drawings.find((entry) => entry.id === 'linear-2').config.unmanagedConfigOverrides = { [WRAP]: stored };
    resetDrawings();
    loadDialogs.length = 0;
    const loaded = await load(session);
    assert.equal(loaded.status, 'ok', JSON.stringify(loaded.error));
    assert.equal(state.drawingOf('linear-2').form.linear_wrap_length, read, stored);
    assert.deepEqual({ ...state.drawingOf('linear-2').unmanagedConfigOverrides }, {}, stored);
    const items = loadDialogs.flatMap((dialog) => dialog.items).filter((item) => item.includes(WRAP));
    assert.equal(items.length, stored === 'abc' ? 1 : 0, `${stored}: ${JSON.stringify(items)}`);
    if (stored === 'abc') assert.ok(items[0].includes('abc'), JSON.stringify(items));
  }
});

// D-36, F3: a v0.14.5 Web Save of a Python API Session keeps the length as typed in the Linear draft and in the
// committed request (tests/fixtures/sessions/HmmtDNA-wrap-5kb.v50-web0145.provenance.json); v0.14.5 drew it wrapped.
// Load reads both as Python does, with no repair row, and the saved Session loads the same.
test('a v0.14.5 Session that stores "5kb" in its draft and request loads wrapped at 5,000 bp, with no repair row', async () => {
  const WRAP = 'canvas.linear.wrap_length';
  const fixture = JSON.parse(gunzipSync(readFileSync(new URL('../fixtures/sessions/HmmtDNA-wrap-5kb.v50-web0145.gbdraw-session.json.gz', import.meta.url))));
  assert.equal(fixture.renderRequest.diagramOptions.configOverrides[WRAP], '5kb');
  assert.equal(fixture.modes.linear.config.unmanagedConfigOverrides[WRAP], '5kb');
  for (const session of [fixture, null]) {
    const document = session ?? await save('wrap-5kb');
    resetDrawings();
    loadDialogs.length = 0;
    const loaded = await load(document);
    assert.equal(loaded.status, 'ok', JSON.stringify(loaded.error));
    const drawing = state.activeDrawing();
    assert.equal(drawing.form.linear_wrap_length, 5000);
    assert.deepEqual({ ...drawing.unmanagedConfigOverrides }, {});
    const items = loadDialogs.flatMap((dialog) => dialog.items).filter((item) => item.includes(WRAP));
    assert.deepEqual(items, [], session ? 'the v0.14.5 Session' : 'its Save');
  }
});

// Issue #1427: a Web Save of a Circular drawing sized 180 x 180 mm
// (tests/fixtures/sessions/HmmtDNA-figure-size-180mm.v51-web.provenance.json) keeps both sides as typed in the
// form and in the committed request; Load reads them back with no repair row, and its Save keeps them.
test('a Session 51 with a 180 x 180 mm figure size loads both sides into the form, and its Save keeps them', async () => {
  const fixture = v51File('HmmtDNA-figure-size-180mm.v51-web.gbdraw');
  const sides = { figure_width: '180mm', figure_height: '180mm' };
  const stored = fixture.drawings.find((entry) => entry.id === 'circular');
  assert.deepEqual({ figure_width: stored.config.form.figure_width, figure_height: stored.config.form.figure_height }, sides);
  for (const session of [fixture, null]) {
    const document = session ?? await save('figure-size-180mm');
    resetDrawings();
    loadDialogs.length = 0;
    const loaded = await load(document);
    assert.equal(loaded.status, 'ok', JSON.stringify(loaded.error));
    const drawing = state.drawingOf('circular');
    assert.deepEqual({ figure_width: drawing.form.figure_width, figure_height: drawing.form.figure_height }, sides);
    const overrides = (session ?? document).drawings.find((entry) => entry.id === 'circular').renderRequest.diagramOptions.configOverrides;
    assert.deepEqual({ width: overrides['canvas.figure_width'], height: overrides['canvas.figure_height'] }, { width: '180mm', height: '180mm' });
    assert.deepEqual(loadDialogs.flatMap((dialog) => dialog.items), [], session ? 'the fixture' : 'its Save');
  }
});

// OV-3259: a Web Session of gbdraw 0.14.5 saved its palette's rows with the user's colors (`config.colors`); its
// committed request holds only the user's rows. After OV-3114 changed misc_feature and default in 53 palettes, Load
// reads a saved row that still holds its palette's 0.14.5 color as the palette's color now, as the CLI replay of the
// committed request draws it; a color the user set stays.
const v0145Session = (edit = () => {}) => {
  const session = v51File('HmmtDNA-figure-size-180mm.v51-web.gbdraw');
  edit(session.drawings.find(({ id }) => id === 'circular'), session);
  return session;
};
const loadedColorsOf = async (session) => {
  resetDrawings();
  loadDialogs.length = 0;
  const loaded = await load(session);
  assert.equal(loaded.status, 'ok', JSON.stringify(loaded.error));
  return { ...state.drawingOf('circular').currentColors.value };
};
const currentPalette = (name) => state.paletteDefinitions.value[name];

test('a 0.14.5 Session reads its palette\'s saved misc_feature and default rows as the palette\'s colors now (OV-3259)', () => withBootPalettes(state, async () => {
  const session = v0145Session();
  const stored = session.drawings[0].config.colors;
  assert.deepEqual([stored.misc_feature, stored.default], ['#d3d3d3', '#d3d3d3']);
  const colors = await loadedColorsOf(session);
  assert.deepEqual([colors.misc_feature, colors.default], [currentPalette('default').misc_feature, currentPalette('default').default]);
  // Every other row is as saved; no row differs from the palette, so a Generate writes no -d row.
  for (const key of Object.keys(stored).filter((key) => !['misc_feature', 'default'].includes(key))) assert.equal(colors[key], stored[key], key);
  assert.deepEqual(buildPaletteColorOverrideRows({ colors, paletteColors: currentPalette('default') }), []);
  // Without the colors its Result shows (a Session no Web Generate saved), a row is read as 0.14.5 read it.
  const unapplied = await loadedColorsOf(v0145Session((drawing) => { delete drawing.ui.appliedPaletteColors; }));
  assert.equal(unapplied.misc_feature, currentPalette('default').misc_feature);
}));

test('Load keeps a saved color the user set: a live edit, a committed row, and a palette whose rows did not change (OV-3259)', () => withBootPalettes(state, async () => {
  // A live Instant Preview edit after the last Generate: the draft and the shown colors hold it, the request does not.
  const live = await loadedColorsOf(v0145Session((drawing) => {
    Object.assign(drawing.config.colors, { CDS: '#123456', misc_feature: '#123456' });
    Object.assign(drawing.ui.appliedPaletteColors, { CDS: '#123456', misc_feature: '#123456' });
  }));
  assert.deepEqual([live.CDS, live.misc_feature, live.default], ['#123456', '#123456', currentPalette('default').default]);
  // A misc_feature #d3d3d3 row the committed Default colors table sets: a user color the last Generate drew.
  const table = 'feature_type\tcolor\nmisc_feature\t#d3d3d3\n';
  const committed = await loadedColorsOf(v0145Session((drawing, session) => {
    session.resources['colors-default-colors'] = {
      kind: 'canonical-tsv', name: 'colors-default-colors-colors-default-colors.tsv', type: 'text/tab-separated-values',
      size: Buffer.byteLength(table), lastModified: 0, encoding: 'base64', data: Buffer.from(table).toString('base64')
    };
    drawing.renderRequest.diagramOptions.colors.defaultColors = { resourceId: 'colors-default-colors', representation: 'canonicalTsv' };
  }));
  assert.deepEqual([committed.misc_feature, committed.default], ['#d3d3d3', currentPalette('default').default]);
  // A palette OV-3114 did not change (earthtones): every saved row stays.
  const earthtones = { ...currentPalette('earthtones') };
  const kept = await loadedColorsOf(v0145Session((drawing) => {
    drawing.config.palette = 'earthtones';
    drawing.config.colors = { ...earthtones };
    drawing.ui.appliedPaletteName = 'earthtones';
    drawing.ui.appliedPaletteColors = { ...earthtones };
    drawing.renderRequest.diagramOptions.colors.defaultColorsPalette = 'earthtones';
  }));
  assert.deepEqual(kept, earthtones);
}));

// OV-3257: Save wrote a drawing's whole color table, its palette's rows included, so a later palette change never
// reached the Session. A drawing that changed none of its palette's colors saves none; one with a user color keeps it.
test('Save stores no palette colors for a drawing whose colors are its palette\'s, so a later palette reaches it (OV-3257)', () => withBootPalettes(state, async () => {
  const { editDefaultColor } = await import('../../gbdraw/web/js/services/drawing/palette.js');
  const slicesOf = (session) => session.drawings.find(({ id }) => id === 'circular').config;
  assert.deepEqual(buildPaletteColorOverrideRows({ colors: await loadedColorsOf(v0145Session()), paletteColors: currentPalette('default') }), []);
  const untouched = await save('palette colors');
  assert.equal(Object.hasOwn(slicesOf(untouched), 'colors'), false);
  assert.deepEqual(await loadedColorsOf(untouched), currentPalette('default'));
  // A later palette: the Session takes its colors.
  const palettes = state.paletteDefinitions.value;
  state.paletteDefinitions.value = { ...palettes, default: { ...palettes.default, CDS: '#010203' } };
  try {
    assert.equal((await loadedColorsOf(untouched)).CDS, '#010203');
  } finally {
    state.paletteDefinitions.value = palettes;
  }
  // A user color is saved and loads as saved.
  await loadedColorsOf(v0145Session());
  editDefaultColor(state.drawingOf('circular').slices.palette, 'CDS', '#123456');
  const edited = await save('a user color');
  assert.equal(slicesOf(edited).colors.CDS, '#123456');
  assert.equal((await loadedColorsOf(edited)).CDS, '#123456');
}));

// The request Generate builds from a drawing's draft. It reads the Legend rows of the loaded Result's SVG; Node has no
// SVG parser, and no row matters to the feature density cases.
const draftRequest = async (drawing) => {
  const parser = globalThis.DOMParser;
  globalThis.DOMParser = class { parseFromString() { return { documentElement: { localName: 'none' } }; } };
  try {
    return (await sessionConfig.projectDraftRequest(state, drawing, {
      linearRecordCatalog: null, recordDisplayRows: [], omittedAnnotationRecordKeys: []
    })).envelope.renderRequest;
  } finally {
    globalThis.DOMParser = parser;
  }
};
const densityFixture = () => v51File('HmmtDNA-feature-density.v51-cli.gbdraw');

// OV-3174, the Circular size shortcut: with Use custom stack off, a feature width sends the stack the Web builds, and
// that stack keeps a row per feature density track after the features row, as Python's implicit stack does.
test('a Circular size shortcut keeps the feature density rows in the stack it sends (OV-3174)', async () => {
  resetDrawings();
  loadDialogs.length = 0;
  assert.equal((await load(densityFixture())).status, 'ok');
  const drawing = state.activeDrawing();
  assert.equal(drawing.adv.circular_track_slots_enabled, false);
  drawing.adv.feature_width_circular = 30;
  const slots = (await draftRequest(drawing)).diagramOptions.tracks.circularTrackSlots;
  const ids = slots.map(({ id }) => id);
  assert.deepEqual(ids.slice(ids.indexOf('features') + 1, ids.indexOf('features') + 3), ['density_genes', 'density_trna'], JSON.stringify(ids));
  assert.deepEqual(slots.filter(({ renderer }) => renderer === 'feature_density').map(({ params }) => params.track), ['genes', 'trna']);
});

// OV-3173: a CLI Session whose stack names feature density rows loads them as density rows with their tracks, and the
// request Generate builds sends them back unchanged.
test('a CLI Session whose stack names feature density rows loads them and sends them back (OV-3173)', async () => {
  const fixture = densityFixture();
  const rows = ['features:features', 'density_genes:feature_density@track=genes', 'gc_content:dinucleotide_content'];
  const request = fixture.drawings[0].renderRequest.diagramOptions.tracks;
  request.circularTrackSlots = rows.map((token) => {
    const [head, options] = token.split('@');
    const [id, renderer] = head.split(':');
    return { kind: 'circularTrackSlot', id, renderer, enabled: true, side: null, radius: null, width: null, z: 0,
      params: options ? Object.fromEntries([options.split('=')]) : {}, innerGapPx: null, outerGapPx: null };
  });
  request.circularTrackAxisIndex = 0;
  resetDrawings();
  loadDialogs.length = 0;
  assert.equal((await load(fixture)).status, 'ok');
  const drawing = state.activeDrawing();
  assert.deepEqual(drawing.adv.circular_track_slots.map(({ id, renderer, params }) => [id, renderer, params.track ?? null]),
    [['features', 'features', null], ['density_genes', 'feature_density', 'genes'], ['gc_content', 'dinucleotide_content', null]]);
  const sent = (await draftRequest(drawing)).diagramOptions.tracks.circularTrackSlots;
  assert.deepEqual(sent.map(({ id, renderer }) => `${id}:${renderer}`), ['features:features', 'density_genes:feature_density', 'gc_content:dinucleotide_content']);
  assert.deepEqual(loadDialogs.flatMap((dialog) => dialog.items), []);
});

// OV-3031: a Session the CLI wrote with two feature density tracks
// (tests/fixtures/sessions/HmmtDNA-feature-density.v51-cli.provenance.json) loads the tracks into the drawing's draft
// (`adv.feature_density_tracks`); the request Generate builds from the draft carries them, Save writes them into the
// drawing's config and its request, and Load of that Save reads them back, with no repair row.
test('a CLI Session with feature density tracks keeps them through Load, the Generate request and Save', async () => {
  const fixture = densityFixture();
  const tracks = fixture.drawings[0].renderRequest.diagramOptions.featureDensityTracks;
  assert.deepEqual(tracks.map(({ id }) => id), ['genes', 'trna']);
  resetDrawings();
  loadDialogs.length = 0;
  assert.equal((await load(fixture)).status, 'ok');
  const drawing = state.activeDrawing();
  assert.deepEqual(drawing.adv.feature_density_tracks, tracks);
  assert.deepEqual((await draftRequest(drawing)).diagramOptions.featureDensityTracks, tracks);
  const saved = await save('feature-density');
  const entry = saved.drawings.find(({ id }) => id === drawing.id);
  assert.deepEqual(entry.config.adv.feature_density_tracks, tracks);
  assert.deepEqual(entry.renderRequest.diagramOptions.featureDensityTracks, tracks);
  resetDrawings();
  assert.equal((await load(saved)).status, 'ok');
  assert.deepEqual(state.activeDrawing().adv.feature_density_tracks, tracks);
  assert.deepEqual((await draftRequest(state.activeDrawing())).diagramOptions.featureDensityTracks, tracks);
  assert.deepEqual(loadDialogs.flatMap((dialog) => dialog.items), []);
  // A Session without tracks: its draft holds none, and its requests have no such key, as before the tracks existed.
  const plain = v51File('HmmtDNA-figure-size-180mm.v51-web.gbdraw');
  resetDrawings();
  assert.equal((await load(plain)).status, 'ok');
  assert.deepEqual(state.activeDrawing().adv.feature_density_tracks, []);
  assert.equal(Object.hasOwn((await draftRequest(state.activeDrawing())).diagramOptions, 'featureDensityTracks'), false);
  const resaved = (await save('no-density')).drawings.find(({ id }) => id === 'circular');
  assert.deepEqual(resaved.renderRequest, plain.drawings.find(({ id }) => id === 'circular').renderRequest);
});

// Issue #1427, review finding 2: a Session the CLI wrote with --figure_width 180mm
// (tests/fixtures/sessions/HmmtDNA-figure-width-180mm.v51-cli.provenance.json) holds Python's full config, with the
// height "auto", and a Result whose root states its size in mm. Load reads the width into the field and the height as
// Auto, with no repair row, and the size line reads the mm root through its px viewBox as 180 mm.
test('a CLI Session sized 180 mm wide loads the width into the field, and the size line reads its mm root', async () => {
  const fixture = v51File('HmmtDNA-figure-width-180mm.v51-cli.gbdraw');
  const { canvas } = fixture.drawings[0].renderRequest.diagramOptions.config;
  assert.deepEqual([canvas.figure_width, canvas.figure_height], ['180mm', 'auto']);
  resetDrawings();
  loadDialogs.length = 0;
  const loaded = await load(fixture);
  assert.equal(loaded.status, 'ok', JSON.stringify(loaded.error));
  const drawing = state.activeDrawing();
  assert.deepEqual([drawing.form.figure_width, drawing.form.figure_height], ['180mm', null]);
  assert.deepEqual(loadDialogs.flatMap((dialog) => dialog.items), []);
  const [rootTag] = /<svg[^>]*>/.exec(state.results.value[0].content);
  const attribute = (name) => new RegExp(`\\s${name}="([^"]*)"`).exec(rootTag)?.[1] ?? null;
  assert.equal(attribute('width'), '180mm');
  const px = getSvgDimensions({ getAttribute: attribute });
  const requested = requestFigureSizeLeaves(sessionConfig.getCommittedCanonicalRenderRequest());
  assert.deepEqual(requested, { width: '180mm', height: null });
  assert.deepEqual(shownFigureSize({ widthPx: px.width, heightPx: px.height, requested }),
    { first: '180 × 133.1 mm', second: '680 × 503 px' });
});

// Italic gene names are one setting per drawing (`config.adv.italic_gene_labels`):
// of two Linear drawings, the one that holds it on keeps it through Load and Save.
test('a Session 51 keeps each Linear drawing\'s own Italic gene names through Load and Save', async () => {
  const session = v49Session('drawings-linear-pair-losat');
  session.drawings.find((entry) => entry.id === 'linear-2').config.adv.italic_gene_labels = true;
  resetDrawings();
  const first = await loadKeepingSlots(session);
  assertLoadedDrawings(session, first);
  assert.equal(state.drawingOf('linear-2').adv.italic_gene_labels, true);
  assert.equal(state.drawingOf('linear').adv.italic_gene_labels, false);
  const result = await exportSession('italic-gene-labels', { readKeptArtifacts: () => [first.opening, ...first.stashed].filter(Boolean) });
  assert.equal(result.status, 'saved', result.error?.message);
  const saved = await readSavedSession(result.blob);
  const italicOf = (id) => saved.drawings.find((entry) => entry.id === id).config.adv.italic_gene_labels;
  assert.equal(italicOf('linear-2'), true);
  assert.equal(italicOf('linear'), false);
});

// The identity colour scale is one setting per drawing (`config.adv.identity_scale`):
// of two Linear drawings, the one that holds a scale keeps it through Load and Save, and
// the other keeps no scale.
test('a Session 51 keeps each Linear drawing\'s own identity colour scale through Load and Save', async () => {
  const session = v49Session('drawings-linear-pair-losat');
  session.drawings.find((entry) => entry.id === 'linear-2').config.adv.identity_scale = [40, 90];
  resetDrawings();
  const first = await loadKeepingSlots(session);
  assertLoadedDrawings(session, first);
  assert.deepEqual(state.drawingOf('linear-2').adv.identity_scale, [40, 90]);
  assert.equal(state.drawingOf('linear').adv.identity_scale, null);
  const result = await exportSession('identity-scale', { readKeptArtifacts: () => [first.opening, ...first.stashed].filter(Boolean) });
  assert.equal(result.status, 'saved', result.error?.message);
  const saved = await readSavedSession(result.blob);
  const scaleOf = (id) => saved.drawings.find((entry) => entry.id === id).config.adv.identity_scale;
  assert.deepEqual(scaleOf('linear-2'), [40, 90]);
  assert.equal(scaleOf('linear'), null);
});

test('Save checks each Session 51 drawing it writes with Save\'s checks, not Load\'s: an input type of no type is refused', async () => {
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

test('Save writes into each Session 51 drawing the LOSAT entries of its own Linear rows (K-3)', async () => {
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
  const saved = await readSavedSession(result.blob);
  const keysOf = (id) => saved.drawings.find((candidate) => candidate.id === id).losatCache.entries.map(({ key }) => key);
  assert.deepEqual(keysOf('linear'), [entry.key]);
  assert.deepEqual(keysOf('linear-2'), [`${entry.key}-second`]);
});

// #1777 round 16 (D-197, OV-3213): the protein identity manifest is each drawing's. Load installs each drawing's own
// field (the shown one in the state, another in its slot) and Save writes each drawing's back; a Session whose writer
// gave every drawing one shared manifest (Session 46 to 50, upgraded by copying it into each drawing) loads it into
// each drawing alike. The Session's shape and version do not change.
test('Load and Save keep each Session 51 drawing\'s own protein identity manifest (D-197)', async () => {
  const feature = (letter) => `f_${letter.repeat(64)}`;
  /** A valid schema-2 manifest of one record instance per uid, named by `tag`. */
  const manifestOf = (tag, uids) => ({
    schema: 2,
    proteinSets: Object.fromEntries(uids.map((uid, index) => [`sha256:set-${tag}-${uid}`,
      { schema: 1, proteins: [{ featureAnalysisId: feature('ab'[index]), aaSha256: 'ab'[index].repeat(64) }] }])),
    recordAnalyses: Object.fromEntries(uids.map((uid) => [`sha256:analysis-${tag}-${uid}`,
      { schema: 1, recordSourceId: uid, proteinSetHash: `sha256:set-${tag}-${uid}` }])),
    recordInstances: Object.fromEntries(uids.map((uid, index) => [uid, {
      schema: 2, recordAnalysisId: `sha256:analysis-${tag}-${uid}`, runtimeBindingHash: `sha256:runtime-${tag}-${uid}`,
      displayBindingHash: `sha256:display-${tag}-${uid}`, runtimeIds: { [feature('ab'[index])]: `h_${'ab'[index].repeat(26)}` },
      featureMetadata: { [feature('ab'[index])]: { displayAlias: `protein-${tag}-${uid}`, exportOrdinal: null } }
    }]))
  });
  for (const shared of [false, true]) {
    const session = v49Session('drawings-linear-pair-losat');
    const uids = session.drawings.find((entry) => entry.id === 'linear').losatCache.entries.flatMap(({ queryUid, subjectUid }) => [queryUid, subjectUid]);
    const manifests = { linear: manifestOf('first', uids), 'linear-2': manifestOf(shared ? 'first' : 'second', uids) };
    session.drawings.forEach((entry) => { entry.proteinIdentityManifest = manifests[entry.id]; });
    session.ui.activeDrawingId = 'linear';
    resetDrawings();
    const loaded = await loadKeepingSlots(session);
    assertLoadedDrawings(session, loaded);
    const slotOf = (id) => [loaded.opening, ...loaded.stashed].find((slot) => slot?.drawingId === id);
    assert.deepEqual(state.proteinIdentityManifest.value, manifests.linear, `shared ${shared}: the shown drawing's`);
    assert.deepEqual(slotOf('linear-2').values.proteinIdentityManifest, manifests['linear-2'], `shared ${shared}: linear-2's slot`);
    assert.notEqual(slotOf('linear-2').values.proteinIdentityManifest, state.proteinIdentityManifest.value, 'two objects');
    const result = await exportSession(`protein-manifests-${shared}`, { readKeptArtifacts: () => [loaded.opening, ...loaded.stashed].filter(Boolean) });
    assert.equal(result.status, 'saved', result.error?.message);
    const saved = await readSavedSession(result.blob);
    assert.equal(saved.version, session.version);
    for (const id of ['linear', 'linear-2']) {
      assert.deepEqual(saved.drawings.find((entry) => entry.id === id).proteinIdentityManifest, manifests[id], `shared ${shared}: ${id} saved`);
    }
  }
});

// #1777 round 18 (OV-3273, D-201): a copy of a Linear drawing keeps its card uids, so after its region changes one
// record instance key names other proteins in each drawing's manifest. A Load validates each drawing's saved protein
// entries with that drawing's own manifest, and Save reuses that check only for the manifest it ran with: Load, Save,
// Load, Save writes into each drawing only the protein entries its own manifest resolves, which the CLI requires.
test('Load and Save keep each Linear drawing\'s protein LOSAT entries to those its own manifest resolves (OV-3273)', async () => {
  const fixture = v51File('BGC0000708-BGC0000709-copy-crop.v51-web.gbdraw');
  const proteinKeysOf = (session) => Object.fromEntries(session.drawings.map((entry) => [entry.id,
    entry.losatCache.entries.filter((cached) => cached.queryRecordInstanceKey !== undefined).map(({ key }) => key).sort()]));
  assert.deepEqual(Object.values(proteinKeysOf(fixture)).map((keys) => keys.length), [4, 4]);
  /** @param {Blob} file A saved Session file. */
  const pythonAdmits = async (file) => JSON.parse(String(await evaluatePython(`
import base64, json, os, tempfile
from gbdraw.session_io import load_session
with tempfile.NamedTemporaryFile(suffix='.gbdraw', delete=False) as handle:
    handle.write(base64.b64decode(payload))
try:
    load_session(handle.name)
    result = json.dumps('admitted')
except Exception as error:
    result = json.dumps(str(error))
finally:
    os.unlink(handle.name)
`, Buffer.from(await file.arrayBuffer()).toString('base64'))));
  let session = fixture;
  for (const round of [1, 2]) {
    resetDrawings();
    const loaded = await loadKeepingSlots(session);
    assert.equal(loaded.loaded.status, 'ok', JSON.stringify(loaded.loaded.error));
    let reused = 0;
    globalThis.__GBDRAW_TEST_HOOKS__ = {
      onStructuralMetric: ({ name }) => { if (name === 'sessionSaveProteinRawTextValidationReuseCount') reused += 1; }
    };
    let result;
    try {
      result = await exportSession(`copy-crop-${round}`, { readKeptArtifacts: () => [loaded.opening, ...loaded.stashed].filter(Boolean) });
    } finally {
      delete globalThis.__GBDRAW_TEST_HOOKS__;
    }
    assert.equal(result.status, 'saved', result.error?.message);
    // Save checks again only what Load did not check with the drawing's manifest: each drawing's own 4 are reused.
    assert.equal(reused, 8, `save ${round}: validations Save reused`);
    session = await readSavedSession(result.blob);
    for (const entry of session.drawings) {
      for (const cached of entry.losatCache.entries.filter((value) => value.queryRecordInstanceKey !== undefined)) {
        assert.equal(validateProteinRawEntryReferences(cached, entry.proteinIdentityManifest), true,
          `save ${round}: ${entry.id} writes ${cached.key}, which its manifest resolves`);
      }
    }
    assert.deepEqual(proteinKeysOf(session), proteinKeysOf(fixture), `save ${round}: each drawing keeps its own entries`);
    assert.equal(await pythonAdmits(result.blob), 'admitted', `save ${round}: the CLI admits it`);
  }
  // Each Save keeps each drawing's manifest as the Web wrote it.
  for (const entry of session.drawings) {
    assert.deepEqual(entry.proteinIdentityManifest, fixture.drawings.find(({ id }) => id === entry.id).proteinIdentityManifest, entry.id);
  }
});

// #1777 round 19 (OV-3277, D-201): a Session 40-50 (what 0.14.5 writes) holds one LOSAT list and one manifest, the
// file's; Load validates the list with it and marks each adopted entry with it, so the first Save checks none again.
// A manifest Load did not check (a new object, as a Generate installs) is never trusted: Save checks every entry.
// The stub converter returns no resource: these Sessions' reversed protein Results need Python's frame conversion,
// which the LOSAT cache entries do not pass through.
test('a Session 47 Load then Save reuses Load\'s check of each protein LOSAT entry, and only with the manifest Load checked (OV-3277)', async () => {
  const { setStoredProteinResultFrameConverter } = await import('../../gbdraw/web/js/services/main-session-comparison-frame.js');
  setStoredProteinResultFrameConverter(async () => ({ result: { resources: {} } }));
  try {
    for (const [name, proteins] of [['BGC0000708-BGC0000713.v47', 25], ['majanivirus_collinear-rc.v47', 4]]) {
      resetDrawings();
      const loaded = await loadKeepingSlots(fixtureSession(name));
      assert.equal(loaded.loaded.status, 'ok', JSON.stringify(loaded.loaded.error));
      const saveCounting = async (title) => {
        const counts = { reused: 0, validated: 0 };
        globalThis.__GBDRAW_TEST_HOOKS__ = {
          onStructuralMetric: ({ name: metric }) => {
            if (metric === 'sessionSaveProteinRawTextValidationReuseCount') counts.reused += 1;
            if (metric === 'sessionSaveProteinRawTextValidationCount') counts.validated += 1;
          }
        };
        let result;
        try {
          result = await exportSession(title, { readKeptArtifacts: () => [loaded.opening, ...loaded.stashed].filter(Boolean) });
        } finally {
          delete globalThis.__GBDRAW_TEST_HOOKS__;
        }
        assert.equal(result.status, 'saved', result.error?.message);
        const saved = await readSavedSession(result.blob);
        const savedProteins = saved.drawings.flatMap((entry) => entry.losatCache?.entries ?? [])
          .filter((entry) => entry.queryRecordInstanceKey !== undefined);
        return [savedProteins.length, counts.reused, counts.validated];
      };
      assert.deepEqual(await saveCounting(`reuse-${name}`), [proteins, proteins, 0], `${name}: Load's manifest`);
      state.proteinIdentityManifest.value = JSON.parse(JSON.stringify(state.proteinIdentityManifest.value));
      assert.deepEqual(await saveCounting(`recheck-${name}`), [proteins, 0, proteins], `${name}: a manifest Load did not check`);
    }
  } finally {
    setStoredProteinResultFrameConverter(null);
  }
});

// #1540 (OV-2379): rows the gbdraw 0.14.4 Web engine searched (`runtime.version: null`) for a search the
// current build gives alike are reused. Load gives each Linear drawing its entries with that record (D-221), the
// reuse rule takes each of them, and Save writes each record back into its own drawing.
test('a Session 51 keeps each Linear drawing\'s 0.14.4 Web LOSAT entries reusable through Load and Save', async () => {
  const { isReusableLosatEntry } = await import('../../gbdraw/web/js/services/losat-cache.js');
  const session = v49Session('drawings-linear-pair-losat');
  const first = session.drawings.find((entry) => entry.id === 'linear');
  const entry = first.losatCache.entries[0];
  entry.runtime = { kind: 'losat', source: 'wasm', version: null, program: entry.program };
  const renamed = (text) => text.replaceAll(entry.queryUid, `${entry.queryUid}-2`).replaceAll(entry.subjectUid, `${entry.subjectUid}-2`);
  const index = session.drawings.findIndex((candidate) => candidate.id === 'linear-2');
  session.drawings[index] = JSON.parse(renamed(JSON.stringify(session.drawings[index])));
  session.drawings[index].losatCache.entries = [JSON.parse(renamed(JSON.stringify({ ...entry, key: `${entry.key}-second`, filename: 'second.tsv' })))];
  for (const opened of ['linear-2', 'linear']) {
    session.ui.activeDrawingId = opened;
    resetDrawings();
    const loaded = await loadKeepingSlots(session);
    assertLoadedDrawings(session, loaded);
    for (const [id, key] of [['linear', entry.key], ['linear-2', `${entry.key}-second`]]) {
      const cached = loadedLosatEvidenceOf(id, loaded).losatCache.get(key);
      assert.deepEqual(cached?.runtime, entry.runtime, `${opened}: ${key}`);
      // OV-2377: the 0.14.4 Wasm's rows are reused for LOSATP only. The fixture's entries are LOSATN, so they
      // are searched again; the LOSATP spelling of the same record is reusable.
      assert.equal(isReusableLosatEntry(cached), false, `${opened}: ${key}`);
      assert.equal(isReusableLosatEntry({ ...cached, program: 'blastp' }), true, `${opened}: ${key}`);
    }
    const result = await exportSession(`losat-runtime-${opened}`, { readKeptArtifacts: () => [loaded.opening, ...loaded.stashed].filter(Boolean) });
    assert.equal(result.status, 'saved', JSON.stringify(result));
    const saved = await readSavedSession(result.blob);
    const runtimesOf = (id) => saved.drawings.find((candidate) => candidate.id === id).losatCache.entries.map(({ runtime }) => runtime);
    assert.deepEqual(runtimesOf('linear'), [entry.runtime], opened);
    assert.deepEqual(runtimesOf('linear-2'), [entry.runtime], opened);
  }
});

// DESIGN §1, B2-D9: a drawing writes the LOSAT entries its request's inputs
// name. linear-2 binds linear's rows but has no request (no Generate), so it
// names none; the shown drawing linear keeps its entry.
test('Save writes no LOSAT entry into a Session 51 drawing without a request, though it binds the entry\'s rows', async () => {
  const session = v49Session('drawings-linear-pair-losat');
  const entry = session.drawings.find((candidate) => candidate.id === 'linear').losatCache.entries[0];
  const second = session.drawings.find((candidate) => candidate.id === 'linear-2');
  Object.assign(second, { renderRequest: null, results: [], editorState: { featureCatalog: null } });
  delete second.runMetadata;
  session.ui.activeDrawingId = 'linear';
  resetDrawings();
  const loaded = await loadKeepingSlots(session);
  assert.equal(loaded.loaded.status, 'ok', JSON.stringify(loaded.loaded.error));
  assert.deepEqual(state.drawingOf('linear-2').linearSeqs.map(({ uid }) => uid), [entry.queryUid, entry.subjectUid]);
  const result = await exportSession('losat-no-request', { readKeptArtifacts: () => [loaded.opening, ...loaded.stashed].filter(Boolean) });
  assert.equal(result.status, 'saved', result.error?.message);
  const saved = await readSavedSession(result.blob);
  const keysOf = (id) => saved.drawings.find((candidate) => candidate.id === id).losatCache.entries.map(({ key }) => key);
  assert.equal(saved.drawings.find((candidate) => candidate.id === 'linear-2').renderRequest, null);
  assert.deepEqual(keysOf('linear'), [entry.key]);
  assert.deepEqual(keysOf('linear-2'), []);
});

// K-3: a Circular conservation entry names no Linear row; only a Circular
// request whose conservation rings read sequences names it. Neither drawing
// here has such rings, so neither writes it.
test('Save writes a Circular conservation LOSAT entry into no drawing whose request has no conservation sequences (K-3)', async () => {
  const session = v49Session('drawings-circular-linear');
  const conservation = {
    schema: 2, kind: 'raw-losat', identityKind: 'nucleotide', program: 'blastn', flow: 'circular-conservation',
    outfmt: '6', args: ['--task', 'megablast'], queryCanonicalHash: 'a'.repeat(64), subjectCanonicalHash: 'b'.repeat(64),
    text: 'q\ts\t100.000\t10\t0\t0\t1\t10\t1\t10\t0.0\t20\n', key: 'c'.repeat(64), filename: 'ring.circular_conservation.losatn.tsv', display: true
  };
  session.drawings.find((candidate) => candidate.id === 'circular').losatCache.entries = [conservation];
  resetDrawings();
  const loaded = await loadKeepingSlots(session);
  assert.equal(loaded.loaded.status, 'ok', JSON.stringify(loaded.loaded.error));
  assert.deepEqual([...loadedLosatEvidenceOf('circular', loaded).losatCache.keys()], [conservation.key]);
  const result = await exportSession('losat-conservation', { readKeptArtifacts: () => [loaded.opening, ...loaded.stashed].filter(Boolean) });
  assert.equal(result.status, 'saved', result.error?.message);
  const saved = await readSavedSession(result.blob);
  assert.deepEqual(Object.fromEntries(saved.drawings.map(({ id, losatCache }) => [id, losatCache.entries.map(({ key }) => key)])),
    { circular: [], linear: [] });
});

// OV-3201 (#1777 review 8), D-221: each drawing's rows are its own evidence's, so a ring edit clears only its
// drawing's rows. A ring entry is each Circular drawing's that holds it, one row each, never a Linear drawing's; a
// Session 27-50's one list gives its ring entries to its Circular part.
const ringCacheEntry = () => ({
  schema: 2, kind: 'raw-losat', identityKind: 'nucleotide', program: 'blastn', flow: 'circular-conservation',
  outfmt: '6', args: ['--task', 'megablast'], queryCanonicalHash: 'a'.repeat(64), subjectCanonicalHash: 'b'.repeat(64),
  text: 'q\ts\t100.000\t10\t0\t0\t1\t10\t1\t10\t0.0\t20\n', key: 'c'.repeat(64), filename: 'ring.circular_conservation.losatn.tsv', display: true
});
/** @param {string} id @param {Awaited<ReturnType<typeof loadKeepingSlots>>} outcome */
const ringRowKeysOf = (id, outcome) => loadedLosatEvidenceOf(id, outcome).losatCacheInfo
  .map((row) => row.key.slice(0, 4)).filter((key) => key === 'cccc');
test('a Load gives a ring entry one cache-info row for each Circular drawing that holds it, none for a Linear one (OV-3201)', async () => {
  for (const [name, holders] of /** @type {const} */ ([
    ['drawings-per-drawing-cards', ['circular', 'circular-2']],
    ['drawings-circular-linear', ['circular', 'linear']]
  ])) {
    const session = v49Session(name);
    session.drawings.filter((candidate) => holders.includes(candidate.id)).forEach((candidate) => { candidate.losatCache.entries = [ringCacheEntry()]; });
    resetDrawings();
    const outcome = await loadKeepingSlots(session);
    assert.equal(outcome.loaded.status, 'ok', JSON.stringify(outcome.loaded.error));
    // Each panel reads its drawing's results, one per cache key (#1777 review 9, OV-3203, D-223): 1 row in each
    // Circular drawing's panel, none in the Linear one's.
    assert.deepEqual(holders.map((id) => ringRowKeysOf(id, outcome)),
      holders.map((id) => (id.startsWith('circular') ? ['cccc'] : [])), name);
    assert.deepEqual(holders.map((id) => losatRowsOnePerKey(loadedLosatEvidenceOf(id, outcome).losatCacheInfo).length),
      holders.map((id) => (id.startsWith('circular') ? 1 : 0)), name);
  }
  const legacy = fixtureSession('two-mode-linear-shown.v48');
  legacy.losatCache.entries.push(ringCacheEntry());
  resetDrawings();
  const outcome = await loadKeepingSlots(legacy);
  assert.equal(outcome.loaded.status, 'ok', JSON.stringify(outcome.loaded.error));
  assert.deepEqual(state.drawingList().map(({ id }) => [id, ringRowKeysOf(id, outcome)]),
    state.drawingList().map(({ id, mode }) => [id, mode === 'circular' ? ['cccc'] : []]),
    'a Session 27-50\'s ring entry names its Circular part');
});

// K-3: a ring's query hash is the SHA-256 of the FASTA the Python ring reader
// (the Worker's) returns; node has no Worker, so the fixture's other ring is
// hashed here through Python.
const readRingInPython = async (_operation, { files: [{ bytes }] }) => ({
  result: JSON.parse(String(await evaluatePython(`
import base64, os, tempfile
from gbdraw.web_support.comparison_sequences import read_comparison_sequence_json
with tempfile.NamedTemporaryFile(delete=False) as handle:
    handle.write(base64.b64decode(payload))
try:
    result = read_comparison_sequence_json(handle.name)
finally:
    os.unlink(handle.name)
`, Buffer.from(new Uint8Array(bytes)).toString('base64'))))
});
// The Worker helper requests (ring reads) started while `callback` runs.
const helperRequestsDuring = async (callback) => {
  const requests = [];
  globalThis.__GBDRAW_TEST_HOOKS__ = {
    onStructuralMetric: ({ name, operation }) => { if (name === 'diagramHelperRequest') requests.push(operation); }
  };
  try {
    return { value: await callback(), requests };
  } finally {
    delete globalThis.__GBDRAW_TEST_HOOKS__;
  }
};
const ringQueryHash = async (text) => {
  const { result } = await readRingInPython('', { files: [{ bytes: Buffer.from(text) }] });
  return createHash('sha256').update(result.fasta, 'utf8').digest('hex');
};
// A Session 51 of `synthetic_conservation` (three LOSAT rings and their entries)
// with a second Circular drawing `circular-2` whose ring 3 reads other bytes
// (one base changed) and whose entries are rings 1 and 2's and that ring 3's.
let conservationSaves = 0;
const twoConservationDrawings = async () => {
  const fixture = JSON.parse(gunzipSync(readFileSync(new URL('../fixtures/sessions/synthetic_conservation.gbdraw-session.json.gz', import.meta.url))));
  resetDrawings();
  const first = await loadKeepingSlots(fixture);
  assert.equal(first.loaded.status, 'ok', JSON.stringify(first.loaded.error));
  // Each call saves under its own title, so no repeat-download question cancels it.
  const result = await exportSession(`conservation-${++conservationSaves}`, { readKeptArtifacts: () => [first.opening, ...first.stashed].filter(Boolean) });
  assert.equal(result.status, 'saved', result.error?.message);
  const session = await readSavedSession(result.blob);
  const circular = session.drawings.find((entry) => entry.id === 'circular');
  const entries = circular.losatCache.entries;
  assert.deepEqual(entries.map(({ key }) => key), fixture.losatCache.entries.map(({ key }) => key));
  // A Session 51 binds each ring's genome as a list of files.
  const ringId = circular.webFiles.bindings.c_conservation_fastas[2][0].resourceId;
  const ringText = Buffer.from(session.resources[ringId].data, 'base64').toString();
  const otherText = ringText.replace('\nACGT', '\nTCGT');
  assert.notEqual(otherText, ringText);
  session.resources['ring-other'] = { ...session.resources[ringId], name: 'ring-other-comparison-d.fasta', data: Buffer.from(otherText).toString('base64') };
  const other = { ...entries[2], key: 'e'.repeat(64), filename: 'ring-other.circular_conservation.losatn.tsv', queryCanonicalHash: await ringQueryHash(otherText) };
  assert.notEqual(other.queryCanonicalHash, entries[2].queryCanonicalHash);
  const second = JSON.parse(JSON.stringify(circular).replaceAll(`"${ringId}"`, '"ring-other"'));
  Object.assign(second, { id: 'circular-2', name: 'Circular 2' });
  second.losatCache.entries = [entries[0], entries[1], other];
  session.drawings.splice(session.drawings.indexOf(circular) + 1, 0, second);
  return { session, keys: { ring1: entries[0].key, ring2: entries[1].key, ring3: entries[2].key, other: other.key } };
};
const savedEntryKeys = async (title, loaded) => {
  const result = await exportSession(title, { readKeptArtifacts: () => [loaded.opening, ...loaded.stashed].filter(Boolean) });
  assert.equal(result.status, 'saved', result.error?.message);
  const saved = await readSavedSession(result.blob);
  return Object.fromEntries(saved.drawings.map(({ id, losatCache }) => [id, losatCache.entries.map(({ key }) => key)]));
};

// K-3 (D-15-124): Save keeps a conservation entry only in the drawing whose
// request names it, by the entry's subject and query hashes: after a Load, the
// pairs of the entries the file wrote into that drawing. Save reads no ring
// file (OV-1908: the first Save after a Load booted the Worker to hash them).
test('Save after a Load reads no ring file and writes each Circular drawing its own conservation entries (K-3)', async () => {
  const { session, keys } = await twoConservationDrawings();
  resetDrawings();
  const loaded = await loadKeepingSlots(session);
  assert.equal(loaded.loaded.status, 'ok', JSON.stringify(loaded.loaded.error));
  const { value: saved, requests } = await helperRequestsDuring(() => savedEntryKeys('two rings', loaded));
  assert.deepEqual(requests, []);
  assert.deepEqual(saved.circular, [keys.ring1, keys.ring2, keys.ring3]);
  assert.deepEqual(saved['circular-2'], [keys.ring1, keys.ring2, keys.other]);
});

// OV-3475 (F1), D-221: Save writes each drawing's own evidence. After Clear in circular-2 (its LOSAT evidence
// emptied in its slot), it writes no entry, not even the ring searches circular still holds.
test('a drawing whose shared ring searches were cleared writes none of them (OV-3475)', async () => {
  const { session, keys } = await twoConservationDrawings();
  resetDrawings();
  const loaded = await loadKeepingSlots(session);
  assert.equal(loaded.loaded.status, 'ok', JSON.stringify(loaded.loaded.error));
  const cleared = (slot) => (slot?.drawingId === 'circular-2'
    ? withArtifactSlotValues(slot, state.drawingOf('circular-2'), emptyLosatEvidence())
    : slot);
  const saved = await savedEntryKeys('cleared circular-2', {
    opening: cleared(loaded.opening), stashed: loaded.stashed.map(cleared)
  });
  assert.deepEqual(saved.circular, [keys.ring1, keys.ring2, keys.ring3]);
  assert.deepEqual(saved['circular-2'], [], 'circular-2 lists no result after its Clear');
});

// K-3 (D-15-124): a drawing's entries are its request's. A ring file bound in
// place of the request's after its Generate (no new Generate) changes nothing
// Save writes: the request's entry stays, the bound file's is not added.
test('Save keeps the request\'s conservation entry for a ring whose bound file is another (K-3)', async () => {
  const { session, keys } = await twoConservationDrawings();
  const second = session.drawings.find((entry) => entry.id === 'circular-2');
  // circular-2's request still names `ring-other`; its ring 3 binds circular's ring 3 file.
  const ringId = session.drawings.find((entry) => entry.id === 'circular').webFiles.bindings.c_conservation_fastas[2][0].resourceId;
  for (const field of ['c_conservation_fastas', 'c_conservation_sequence_sources']) {
    second.webFiles.bindings[field][2] = [{ ...second.webFiles.bindings[field][2][0], resourceId: ringId }];
  }
  resetDrawings();
  const loaded = await loadKeepingSlots(session);
  assert.equal(loaded.loaded.status, 'ok', JSON.stringify(loaded.loaded.error));
  const saved = await savedEntryKeys('replaced ring', loaded);
  assert.deepEqual(saved.circular, [keys.ring1, keys.ring2, keys.ring3]);
  assert.deepEqual(saved['circular-2'], [keys.ring1, keys.ring2, keys.other]);
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
    return await readSavedSession(result.blob);
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

// OV-1811: Save writes the drawings a one-mode Session 51 Load made, and no drawing of the other mode.
test('Save after a one-mode Session 51 Load writes no drawing of the other mode (R16-8, OV-1811)', () => withBootPalettes(state, async () => {
  resetDrawings();
  const first = await loadKeepingSlots(v49Session('drawings-linear-pair-losat'));
  const result = await exportSession('ov-1811', { readKeptArtifacts: () => [first.opening, ...first.stashed].filter(Boolean) });
  assert.equal(result.status, 'saved', result.error?.message);
  const saved = await readSavedSession(result.blob);
  assert.deepEqual(saved.drawings.map((entry) => entry.id), ['linear', 'linear-2']);
}));

// OV-1811, a Session 48 the move gives one drawing (its Linear part unused): Load makes no Linear
// drawing (D-15-120), so Save writes the Circular drawing alone.
test('a Session 48 Load whose Linear part is unused makes no Linear drawing; Save writes the Circular one (R16-8, OV-1811)', () => withBootPalettes(state, async () => {
  assert.deepEqual(drawingsVector('file-card-imports-circular.v48').drawings.map((drawing) => drawing.mode), ['circular']);
  resetDrawings();
  const loaded = await load(fixtureSession('file-card-imports-circular.v48'));
  assert.equal(loaded.status, 'ok', JSON.stringify(loaded.error));
  assert.deepEqual(state.drawingList().map((drawing) => drawing.id), ['circular']);
  const saved = await save('ov-1811 48');
  assert.deepEqual(saved.drawings.map((entry) => entry.id), ['circular']);
}));

// D-15-120: a Load makes only the drawings the Session holds; the bar's offer tab, not the Load,
// makes a drawing of the other mode.
test('a one-mode Session 51 Load lists only its own drawings (D-15-120)', () => withBootPalettes(state, async () => {
  resetDrawings();
  const loaded = await loadKeepingSlots(v49Session('drawings-linear-pair-losat'));
  assert.equal(loaded.loaded.status, 'ok', JSON.stringify(loaded.loaded.error));
  assert.deepEqual(state.drawingList().map((drawing) => drawing.id), ['linear', 'linear-2']);
  assert.deepEqual(state.drawingTabs.value.map((tab) => tab.id), ['linear', 'linear-2']);
}));

// PD-OI-115: a Web Save of a generated Circular HmmtDNA drawing
// (tests/fixtures/sessions/HmmtDNA-region-axis.v51-web.provenance.json) keeps each record's recordKey and axis in its
// track-slot geometry. Load reads them back with no repair row, so Select region reads the axis without Generate,
// and its Save keeps them.
test('a Session 51 whose track-slot geometry has the axis loads it back for Select region, and its Save keeps it', async () => {
  const { readRegionAxisRecords } = await import('../../gbdraw/web/js/services/region-axis.js');
  const fixture = v51File('HmmtDNA-region-axis.v51-web.gbdraw');
  const [stored] = fixture.drawings[0].runMetadata.trackSlotGeometry.records;
  for (const session of [fixture, null]) {
    const document = session ?? await save('region-axis');
    resetDrawings();
    loadDialogs.length = 0;
    const loaded = await load(document);
    assert.equal(loaded.status, 'ok', JSON.stringify(loaded.error));
    const records = readRegionAxisRecords(state.trackSlotResolvedGeometry.value, 0);
    assert.deepEqual(records.map(({ recordKey, recordId, axis }) => ({ recordKey, recordId, axis })),
      [{ recordKey: 'record-1', recordId: 'NC_012920.1', axis: stored.axis }], session ? 'the fixture' : 'its Save');
    assert.equal(records[0].axis.source.drawnLength, 16_569);
    assert.deepEqual(loadDialogs.flatMap((dialog) => dialog.items), []);
  }
});
