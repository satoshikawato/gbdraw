// Session 46 (plan OV-80 §4, §8 S1, S3, S4, S6, S7): each diagram mode's
// drawing is saved as its own slice (`modes.<mode>`) and loaded back into its
// own drawing; a slice holds registry fields only; Session 45 and a Session 46
// with a retired field are rejected.
import assert from 'node:assert/strict';
import { test } from 'node:test';
import { readFileSync } from 'node:fs';
import { gunzipSync } from 'node:zlib';
import { installSessionImportWorker } from './helpers/session-import-node.mjs';

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

const { SESSION_VERSION, exportSession, importSession } = await import('../../gbdraw/web/js/services/config.js');
const { MODE_SLICE_CONTAINERS } = await import('../../gbdraw/web/js/services/mode-scoped-migration.js');
const { createDefaultAdv, createDefaultForm } = await import('../../gbdraw/web/js/services/session-active-config-contract.js');
const { state } = await import('../../gbdraw/web/js/state.js');

const MODES = ['circular', 'linear'];
const save = async (title) => {
  const result = await exportSession(title);
  assert.equal(result.status, 'saved', result.error?.message);
  return JSON.parse(gunzipSync(Buffer.from(await result.blob.arrayBuffer())));
};
const load = (session) => importSession({
  target: { files: [new Blob([JSON.stringify(session)], { type: 'application/json' })], value: 'selected' }
});
const resetDrawings = () => {
  for (const mode of MODES) {
    Object.assign(state.drawings[mode].form, createDefaultForm());
    Object.assign(state.drawings[mode].adv, createDefaultAdv(mode));
  }
};
// Every key of a slice, as `container.key` paths.
const sliceKeys = (value, container = '') => Object.entries(value).flatMap(([key, child]) => {
  const path = container ? `${container}.${key}` : key;
  return [[container, key], ...(Object.hasOwn(MODE_SLICE_CONTAINERS, path) && child && typeof child === 'object'
    ? sliceKeys(child, path) : [])];
});

test('Session 46 saves each drawing as its mode slice and loads it back into that drawing (S1, S6)', async () => {
  assert.equal(SESSION_VERSION, 46);
  resetDrawings();
  state.mode.value = 'circular';
  state.drawings.circular.adv.label_font_size = 13;
  state.drawings.circular.form.separate_strands = false;
  state.drawings.circular.adv.gc_content_mode = 'percent';
  state.drawings.linear.adv.label_font_size = 7;
  state.drawings.linear.adv.depth_min = 3;
  state.drawings.linear.form.plot_title = 'Linear only';
  // OV-120: a renamed row the last Linear Generate did not draw.
  // A Legend row holds no `showStroke` since OV-157 (view state, #940).
  const waiting = { caption: 'Coverage', originalCaption: 'depth', color: '#7b2cbf', featureIds: [] };
  state.drawings.linear.dormantLegendEntries.value = [waiting];
  const saved = await save('two drawings');

  assert.equal(saved.version, 46);
  assert.equal(saved.renderRequest, null);
  for (const retired of ['config', 'features']) assert.equal(Object.hasOwn(saved, retired), false, retired);
  for (const retired of ['layoutPreferences', 'canvasPadding', 'pendingPaletteName', 'pendingPaletteColors',
    'linearTypographyLinked']) {
    assert.equal(Object.hasOwn(saved.ui, retired), false, `ui.${retired}`);
  }
  assert.equal(Object.hasOwn(saved.editorState, 'featureStrokes'), false);
  assert.deepEqual(Object.keys(saved.modes).sort(), MODES);
  assert.equal(saved.modes.circular.config.adv.label_font_size, 13);
  assert.equal(saved.modes.circular.config.form.separate_strands, false);
  assert.equal(saved.modes.circular.config.adv.gc_content_mode, 'percent');
  assert.equal(saved.modes.linear.config.adv.label_font_size, 7);
  assert.equal(saved.modes.linear.config.adv.depth_min, 3);
  assert.equal(saved.modes.linear.config.form.plot_title, 'Linear only');
  assert.equal(saved.modes.linear.config.form.separate_strands, createDefaultForm().separate_strands);
  // The waiting row follows the shown rows of its slice, marked `dormant`.
  assert.deepEqual(saved.modes.linear.editorState.legend.entries, [{ ...waiting, dormant: true }]);
  assert.deepEqual(saved.modes.circular.editorState.legend.entries, []);
  assert.equal(typeof saved.ui.losatExecution, 'object');
  assert.equal(typeof saved.ui.richFeaturePopup, 'boolean');

  // S7: a slice holds registry fields only.
  for (const mode of MODES) {
    const unknown = sliceKeys(saved.modes[mode]).filter(([container, key]) => !MODE_SLICE_CONTAINERS[container]?.has(key));
    assert.deepEqual(unknown, [], `modes.${mode} holds only registry fields`);
  }

  resetDrawings();
  state.mode.value = 'linear';
  const loaded = await load(saved);
  assert.equal(loaded.status, 'ok', JSON.stringify(loaded.error));
  assert.equal(state.mode.value, 'circular');
  assert.equal(state.drawings.circular.adv.label_font_size, 13);
  assert.equal(state.drawings.circular.form.separate_strands, false);
  assert.equal(state.drawings.circular.adv.gc_content_mode, 'percent');
  assert.equal(state.drawings.linear.adv.label_font_size, 7);
  assert.equal(state.drawings.linear.adv.depth_min, 3);
  assert.equal(state.drawings.linear.form.plot_title, 'Linear only');
  assert.notEqual(state.drawings.circular.form.plot_title, 'Linear only');
  assert.deepEqual(state.drawings.linear.dormantLegendEntries.value, [waiting]);
  assert.deepEqual(state.drawings.linear.legendEntries.value, []);
  assert.deepEqual(state.drawings.circular.dormantLegendEntries.value, []);

  // A slice may omit any field, and the slice of a mode not shown may be
  // absent: Load fills them with that mode's defaults. (A settings-only
  // Session keeps the slice of its shown mode.)
  const partial = structuredClone(saved);
  delete partial.modes.circular.config.adv.label_font_size;
  delete partial.modes.linear;
  resetDrawings();
  state.drawings.linear.dormantLegendEntries.value = [];
  state.drawings.circular.adv.label_font_size = 99;
  state.drawings.linear.adv.depth_min = 9;
  const loadedPartial = await load(partial);
  assert.equal(loadedPartial.status, 'ok', JSON.stringify(loadedPartial.error));
  assert.equal(state.drawings.circular.adv.label_font_size, createDefaultAdv('circular').label_font_size);
  assert.equal(state.drawings.circular.form.separate_strands, false);
  assert.equal(state.drawings.linear.adv.depth_min, createDefaultAdv('linear').depth_min);
  assert.equal(state.drawings.linear.form.plot_title, '');
});

test('Session 45 and a Session 46 with a retired field or a non-registry slice field are rejected (S3, S4)', async () => {
  resetDrawings();
  const saved = await save('current');
  const rejected = async (session, label) => {
    const outcome = await load(session);
    assert.equal(outcome.status, 'error', label);
  };
  await rejected({ ...structuredClone(saved), version: 45 }, 'Session 45');
  await rejected({ ...structuredClone(saved), config: structuredClone(saved.modes.circular.config) }, 'top-level config');
  await rejected({ ...structuredClone(saved), features: {} }, 'top-level features');
  const retiredUi = structuredClone(saved);
  retiredUi.ui.canvasPadding = { top: 0, right: 0, bottom: 0, left: 0 };
  await rejected(retiredUi, 'ui.canvasPadding');
  const retiredLegend = structuredClone(saved);
  retiredLegend.editorState.legend = { ...retiredLegend.editorState.legend, entries: [] };
  await rejected(retiredLegend, 'editorState.legend.entries');
  const unknownSliceField = structuredClone(saved);
  unknownSliceField.modes.linear.config.modeProfiles = {};
  await rejected(unknownSliceField, 'modes.linear.config.modeProfiles');
  const thirdMode = structuredClone(saved);
  thirdMode.modes.radial = {};
  await rejected(thirdMode, 'modes.radial');
});

// OV-160: a named CSS color in a Legend color or a stroke override (also the
// stroke an SVG had before the edit) loads as its CSS hex value, the value the
// Python split writes; an unknown name loads as no color, as before.
test('named Legend and stroke colors load as their CSS hex values (OV-160)', async () => {
  resetDrawings();
  state.mode.value = 'circular';
  const saved = await save('named colors');
  const named = structuredClone(saved);
  const legend = named.modes.circular.editorState.legend;
  legend.colorOverrides = { CDS: 'SeaShell' };
  legend.strokeOverrides = {
    CDS: { strokeColor: 'navy', strokeWidth: 2, originalStrokeColor: 'gray', originalStrokeWidth: 1 }
  };
  named.modes.circular.editorState.featureStrokes = {
    overrides: { feature_1: { strokeColor: 'rebeccapurple', strokeWidth: 3, originalStrokeColor: 'darkgrey' } }
  };
  named.editorState.originalSvgStroke = { color: 'gray', width: 1 };
  const loaded = await load(named);
  assert.equal(loaded.status, 'ok', JSON.stringify(loaded.error));
  const drawing = state.drawings.circular;
  assert.deepEqual({ ...drawing.legendColorOverrides }, { CDS: '#fff5ee' });
  assert.deepEqual({ ...drawing.legendStrokeOverrides.CDS },
    { strokeColor: '#000080', strokeWidth: 2, originalStrokeColor: '#808080', originalStrokeWidth: 1 });
  assert.deepEqual({ ...drawing.featureStrokeOverrides.feature_1 },
    { strokeColor: '#663399', strokeWidth: 3, originalStrokeColor: '#a9a9a9' });
  assert.equal(state.originalSvgStroke.value.color, '#808080');

  const unknown = structuredClone(saved);
  unknown.modes.circular.editorState.legend.strokeOverrides = {
    CDS: { strokeColor: 'notacolor', strokeWidth: 2, originalStrokeColor: 'notacolor' }
  };
  const loadedUnknown = await load(unknown);
  assert.equal(loadedUnknown.status, 'ok', JSON.stringify(loadedUnknown.error));
  assert.deepEqual({ ...state.drawings.circular.legendStrokeOverrides.CDS }, { strokeWidth: 2, originalStrokeColor: null });
});

// D-43 (OV-300): a stored Legend entry whose color lies outside the Default
// colors (-d) domain is dropped at Load and named in the Load notice; a named
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
  state.mode.value = 'circular';
  const current = structuredClone(await save('legend colors'));
  current.modes.circular.editorState.legend.entries = structuredClone(stored);
  current.modes.circular.editorState.legend.deletedEntries = [{ caption: 'W', color: 'currentColor', featureIds: [] }];
  alerts.length = 0;
  const loaded = await load(current);
  assert.equal(loaded.status, 'ok', JSON.stringify(loaded.error));
  assert.deepEqual(shownRows(state.drawings.circular), kept);
  assert.deepEqual(state.drawings.circular.deletedLegendEntries.value, []);
  assert.match(alerts.at(-1),
    / Legend: entries "X", "W" \(deleted\) had a color the app does not accept and were dropped\.$/);

  const older = JSON.parse(gunzipSync(readFileSync(new URL('../fixtures/sessions/settings-only.v42.json.gz', import.meta.url))));
  // A caption whose first row has a bad color but a later row a valid one is shown, not named.
  older.editorState.legend.entries = [{ caption: 'Z', color: 'notacolor', featureIds: [] }, ...structuredClone(stored)];
  resetDrawings();
  alerts.length = 0;
  const loadedOlder = await load(older);
  assert.equal(loadedOlder.status, 'ok', JSON.stringify(loadedOlder.error));
  assert.deepEqual(shownRows(state.activeDrawing()), kept);
  assert.match(alerts.at(-1), / Legend: entry "X" had a color the app does not accept and was dropped\.$/);
  assert.equal(alerts.at(-1).match(/Legend:/g).length, 1);
});

// OV-337: Load leaves a Linear row without a Depth file as the Session binds it
// (`depth: null`), so Save writes it back unchanged. Load pads Depth rows to
// the Depth track count only when a row holds a Depth file, as for Circular.
test('Load and Save keep a Linear row without a Depth file as depth null (OV-337)', async () => {
  resetDrawings();
  state.mode.value = 'circular';
  const gallery = JSON.parse(readFileSync(new URL('../../gbdraw/web/gallery/sessions/HmmtDNA_basic_circular.gbdraw-session.json', import.meta.url), 'utf8'));
  assert.equal((await load(gallery)).status, 'ok');
  const first = await save('gallery round trip');
  const rows = first.webFiles.bindings.linearSeqs;
  assert.deepEqual(rows.map((row) => row.depth), [null]);
  assert.equal(first.modes.linear.config.adv.depth_tracks.length, 1, 'the Linear drawing has a Depth series');
  assert.equal((await load(first)).status, 'ok');
  assert.deepEqual(state.linearSeqs.map((seq) => seq.depth), [null]);
  const second = await save('gallery round trip again');
  assert.deepEqual(second.webFiles.bindings.linearSeqs.map((row) => row.depth), rows.map((row) => row.depth));
});
