// Session 46 (plan OV-80 §4, §8 S1, S3, S4, S6, S7): each diagram mode's
// drawing is saved as its own slice (`modes.<mode>`) and loaded back into its
// own drawing; a slice holds registry fields only; Session 45 and a Session 46
// with a retired field are rejected.
import assert from 'node:assert/strict';
import { test } from 'node:test';
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
globalThis.alert = () => {};
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
  const waiting = { caption: 'Coverage', originalCaption: 'depth', color: '#7b2cbf', showStroke: false, featureIds: [] };
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
