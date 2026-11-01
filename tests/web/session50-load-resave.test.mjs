// A Session 50 Load keeps every drawing's request metadata (`webFiles` beside
// its bindings) from its own entry, never from the drawing the Load opened, so
// a re-Save writes it back (OV-1872), and a re-Save of drawings without a
// Result writes what the first Save wrote (OV-1902).
import assert from 'node:assert/strict';
import { test } from 'node:test';
import { readFileSync } from 'node:fs';
import { gunzipSync } from 'node:zlib';
import { File } from 'node:buffer';
import { installSessionImportWorker } from './helpers/session-import-node.mjs';
import { validateWebConfigOverrides } from './helpers/python-evaluator.mjs';
import { installFakeSvgDom } from './fake-svg-dom.mjs';
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
globalThis.document = { body: { appendChild: () => {} }, createElement: () => ({ click: () => {}, addEventListener: () => {}, parentNode: null }) };
installFakeSvgDom();
globalThis.File = File;
globalThis.alert = () => {};
installSessionImportWorker();

const sessionConfig = await import('../../gbdraw/web/js/services/config.js');
sessionConfig.setUnmanagedConfigOverrideValidator(async (payload) => ({ result: await validateWebConfigOverrides(payload) }));
const { state } = await import('../../gbdraw/web/js/state.js');
await (await import('./helpers/drawing-state.mjs')).addStartLinearDrawing();
const { createArtifactSlot } = await import('../../gbdraw/web/js/services/artifact-slot.js');
const { createPaletteLoader } = await import('../../gbdraw/web/js/app/palettes.js');
const { exportSession, importSession } = callerSession(sessionConfig, state, await sessionRollbackPort(state));

const v49Session = () => JSON.parse(gunzipSync(readFileSync(new URL(
  '../fixtures/sessions/drawings-circular-linear.v50.gbdraw-session.json.gz', import.meta.url
))));
// The Load hands its caller the artifact slot of each drawing it does not
// show (app-setup keeps them; Save writes each drawing's request metadata from
// its slot's committed request).
const loadKeepingSlots = async (session) => {
  const stashed = [];
  const loaded = await importSession({
    target: { files: [new File([JSON.stringify(session)], 'session.json', { type: 'application/json' })], value: 'selected' }
  }, { confirmLoadRepairs: async () => true, installLoadedArtifactSlots: (slots) => { stashed.push(...slots.stashed); } });
  assert.equal(loaded.status, 'ok', JSON.stringify(loaded.error));
  return stashed;
};
/** `webFiles` without its bindings (Save writes the bindings from the drawing's inputs). */
const requestMetadata = (webFiles) => {
  const { bindings: _bindings, resourceOriginalNames: _names, ...metadata } = webFiles ?? {};
  return metadata;
};

const savedDocument = async (title, readKeptArtifacts) => {
  // A second Save of the same title downloads again.
  const result = await exportSession(title, { readKeptArtifacts, askChoice: async () => 'download' });
  assert.equal(result?.status, 'saved', JSON.stringify(result));
  return JSON.parse(gunzipSync(Buffer.from(await result.blob.arrayBuffer())));
};

// OV-1902: a drawing switch keeps the drawing left in its artifact slot, empty
// when it has no Result; a Load keeps none for a drawing without a request.
// What Save writes of a drawing without a Result does not depend on that slot,
// so a Save of the loaded drawings writes them again unchanged. Runs first, on
// the new app's drawings.
test('Save, Load and Save of drawings without a Result write the same drawings (OV-1902)', async () => {
  // The boot palette install, as app.js runs it after setup.
  const palettes = JSON.parse(readFileSync(new URL('../../gbdraw/web/gallery/palettes/palettes.json', import.meta.url), 'utf8'));
  assert.ok(createPaletteLoader({ state }).applyPalettes(palettes.palettes));
  // Each drawing holds one Depth series as set, so the case does not depend
  // on how setup or Load fill an empty list.
  for (const drawing of state.drawingList()) {
    drawing.adv.depth_tracks = [{
      color: '#4A90E2', height: null, label: '', large_tick_interval: null, small_tick_interval: null, tick_font_size: null
    }];
  }
  state.drawingOf('circular').form.labels_mode = 'both';
  state.drawingOf('linear').form.show_labels_linear = 'all';
  assert.equal(state.activeDrawing().id, 'circular');
  const switchKept = [createArtifactSlot({ drawing: state.drawingOf('linear'), values: {} })];
  const first = await savedDocument('ov-1902', () => switchKept);
  assert.deepEqual(first.drawings.map((entry) => entry.id), ['circular', 'linear']);
  const stashed = await loadKeepingSlots(first);
  const second = await savedDocument('ov-1902', () => stashed);
  assert.deepEqual(second.drawings, first.drawings);
});

// OV-1904: Save writes the shown drawing's inputs without a Generate as its
// draft request, which a Load keeps as the drawing's request (a slot per
// drawing with a request), so Save writes that drawing's Result editor and
// view fields as it does for a committed request, and a second Save writes
// them again. Runs on the drawings the case above loaded.
test('Save, Load and Save of a shown drawing saved with a draft request write the same Result fields (OV-1904)', async () => {
  assert.equal(state.activeDrawing().id, 'circular');
  state.activeDrawing().files.c_gb = state.inputFiles.registerFile(new File(
    [readFileSync(new URL('../fixtures/sessions/cli-web-mito.gb', import.meta.url))], 'cli-web-mito.gb', { type: 'text/plain' }
  ));
  const resultFields = (document) => {
    const { editorState, runMetadata, ui } = document.drawings.find((entry) => entry.id === 'circular');
    return { editorState, runMetadata, ui };
  };
  const first = await savedDocument('ov-1904', () => []);
  assert.equal(first.drawings.find((entry) => entry.id === 'circular').renderRequest?.mode, 'circular');
  const stashed = await loadKeepingSlots(first);
  const second = await savedDocument('ov-1904', () => stashed);
  assert.deepEqual(resultFields(second), resultFields(first));
});

for (const opened of ['linear', 'circular']) {
  test(`a Session 50 Load keeps each drawing's own request metadata in the drawing it does not open (${opened} opened)`, async () => {
    const session = v49Session();
    session.ui.activeDrawingId = opened;
    const stashed = await loadKeepingSlots(session);
    assert.equal(state.activeDrawing().id, opened);
    const other = session.drawings.find((entry) => entry.id !== opened);
    const slot = stashed.find((entry) => entry.drawingId === other.id);
    assert.ok(slot, JSON.stringify(stashed.map((entry) => entry.drawingId)));
    assert.deepEqual(requestMetadata(slot.runtimeState.canonical.committedCanonicalSession.webFiles), requestMetadata(other.webFiles));
  });
}

// OV-1871: what a Session 50 Load cannot restore of one drawing is listed
// before anything changes (BA-D25), and the drawing stays without its set
// (`withStoredResultSet`); only a newer version refuses (〔拒否〕).
const loadRows = async (session) => {
  let asked = false;
  const loaded = await importSession({
    target: { files: [new File([JSON.stringify(session)], 'session.json', { type: 'application/json' })], value: 'selected' }
  }, { confirmLoadRepairs: async () => { asked = true; return true; } });
  assert.equal(loaded.status, 'ok', JSON.stringify(loaded.error ?? loaded.status));
  assert.equal(asked, true, 'the Load lists its rows before it changes anything');
  return loaded.repairs.map(({ code, path }) => `${code} ${path}`);
};
for (const [name, mutate] of [
  // A drawing without a request holds no Results (Session format): a null
  // request beside Results reads as one Load cannot read (E-OD29).
  ['a null request beside Results', (entry) => { entry.renderRequest = null; }],
  ['a missing request beside Results', (entry) => { delete entry.renderRequest; }]
]) {
  test(`a Session 50 Load lists the opened drawing's request of ${name} and keeps the drawing (OV-1871)`, async () => {
    const session = v49Session();
    session.ui.activeDrawingId = 'circular';
    const index = session.drawings.findIndex((entry) => entry.id === 'circular');
    mutate(session.drawings[index]);
    assert.deepEqual(await loadRows(session), [`request-dropped drawings[${index}].renderRequest`]);
  });
}
test('a Session 50 Load lists the request of a drawing it does not open, without a request beside Results (OV-1871)', async () => {
  const session = v49Session();
  session.ui.activeDrawingId = 'linear';
  const index = session.drawings.findIndex((entry) => entry.id === 'circular');
  session.drawings[index].renderRequest = null;
  assert.deepEqual(await loadRows(session), [`other-result-set-dropped drawings[${index}].renderRequest`]);
});
test('a Session 50 whose resources are gone lists the opened drawing\'s request and its bindings, not a refusal (OV-1871)', async () => {
  const session = v49Session();
  session.ui.activeDrawingId = 'circular';
  session.resources = {};
  const rows = await loadRows(session);
  assert.ok(rows.includes('request-dropped drawings[0].renderRequest'), JSON.stringify(rows));
});

// V50a (D-15-175): dev's Session 49 is the mode-scoped layout that adds the
// Circular slice's `adv.show_definition` (PD-OI-109). Load reads it like a
// Session 46-48 into drawings, so its Circular drawing keeps Show Definition
// off, and Save writes the drawings layout (Session 50) with it off.
test('a mode-scoped Session 49 with Show Definition off loads into a Circular drawing with it off and saves as 50 with it off', async () => {
  const session = JSON.parse(gunzipSync(readFileSync(new URL(
    '../fixtures/sessions/hidden-definition-circular.v49.gbdraw-session.json.gz', import.meta.url
  ))));
  assert.equal(session.version, 49);
  assert.equal(session.modes.circular.config.adv.show_definition, false);
  const stashed = await loadKeepingSlots(session);
  assert.equal(state.drawingOf('circular').adv.show_definition, false);
  const saved = await savedDocument('show-definition-off', () => stashed);
  assert.equal(saved.version, 50);
  const circular = saved.drawings.find((entry) => entry.id === 'circular');
  assert.equal(circular?.config.adv.show_definition, false);
});
