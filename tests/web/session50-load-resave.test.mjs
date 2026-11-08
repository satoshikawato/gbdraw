// A Session 51 Load keeps every drawing's request metadata (`webFiles` beside
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
  '../fixtures/sessions/drawings-circular-linear.v51.gbdraw-session.json.gz', import.meta.url
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
  const { bindCircularSourceFile, createCircularSource, replaceCircularSources } = await import('../../gbdraw/web/js/services/drawing/inputs.js');
  replaceCircularSources(state.activeDrawing().slices.inputs, [createCircularSource({ uid: 'record-1' })]);
  bindCircularSourceFile(state.inputFiles, state.activeDrawing().slices.inputs, 'record-1', 'gb', new File(
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
  test(`a Session 51 Load keeps each drawing's own request metadata in the drawing it does not open (${opened} opened)`, async () => {
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

// OV-1871: what a Session 51 Load cannot restore of one drawing is listed
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
  test(`a Session 51 Load lists the opened drawing's request of ${name} and keeps the drawing (OV-1871)`, async () => {
    const session = v49Session();
    session.ui.activeDrawingId = 'circular';
    const index = session.drawings.findIndex((entry) => entry.id === 'circular');
    mutate(session.drawings[index]);
    assert.deepEqual(await loadRows(session), [`request-dropped drawings[${index}].renderRequest`]);
  });
}
test('a Session 51 Load lists the request of a drawing it does not open, without a request beside Results (OV-1871)', async () => {
  const session = v49Session();
  session.ui.activeDrawingId = 'linear';
  const index = session.drawings.findIndex((entry) => entry.id === 'circular');
  session.drawings[index].renderRequest = null;
  assert.deepEqual(await loadRows(session), [`other-result-set-dropped drawings[${index}].renderRequest`]);
});
test('a Session 51 whose resources are gone lists the opened drawing\'s request and its bindings, not a refusal (OV-1871)', async () => {
  const session = v49Session();
  session.ui.activeDrawingId = 'circular';
  session.resources = {};
  const rows = await loadRows(session);
  assert.ok(rows.includes('request-dropped drawings[0].renderRequest'), JSON.stringify(rows));
});

// V50a (D-15-175): dev's Session 49 is the mode-scoped layout that adds the
// Circular slice's `adv.show_definition` (PD-OI-109). Load reads it like a
// Session 46-48 into drawings, so its Circular drawing keeps Show Definition
// off, and Save writes the drawings layout (Session 51) with it off.
test('a mode-scoped Session 49 with Show Definition off loads into a Circular drawing with it off and saves as 51 with it off', async () => {
  const session = JSON.parse(gunzipSync(readFileSync(new URL(
    '../fixtures/sessions/hidden-definition-circular.v49.gbdraw-session.json.gz', import.meta.url
  ))));
  assert.equal(session.version, 49);
  assert.equal(session.modes.circular.config.adv.show_definition, false);
  const stashed = await loadKeepingSlots(session);
  assert.equal(state.drawingOf('circular').adv.show_definition, false);
  const saved = await savedDocument('show-definition-off', () => stashed);
  assert.equal(saved.version, sessionConfig.SESSION_VERSION);
  const circular = saved.drawings.find((entry) => entry.id === 'circular');
  assert.equal(circular?.config.adv.show_definition, false);
});

// D-15-190: tag v0.14.5 writes a flat Session 50 (the `modes` layout). Load
// reads it like Sessions 46-49 into drawings; the Linear drawing keeps the
// wrap length, italic gene names and the identity scale, and Save writes them
// in that drawing of a Session 51.
const releasedFlat50 = () => JSON.parse(gunzipSync(readFileSync(new URL(
  '../fixtures/sessions/cli-linear-wrapped.v50.gbdraw-session.json.gz', import.meta.url
))));
/**
 * A config leaf of a drawing's stored request: its full `config` (a CLI
 * Session's) or its `configOverrides` (a Web Save's).
 * @param {{ renderRequest: { diagramOptions: { config?: Record<string, any>, configOverrides?: Record<string, unknown> } } }} drawing
 * @param {string} path
 */
const requestLeaf = (drawing, path) => {
  const { config, configOverrides } = drawing.renderRequest.diagramOptions;
  return config ? path.split('.').reduce((node, key) => node[key], config) : configOverrides[path];
};

test('a released flat Session 50 loads into a Linear drawing with its wrap length, italic names and identity scale, and saves as 51', async () => {
  const session = releasedFlat50();
  assert.equal(session.version, 50);
  const stashed = await loadKeepingSlots(session);
  const linear = state.drawingOf('linear');
  assert.equal(linear.adv.italic_gene_labels, true);
  assert.deepEqual(linear.adv.identity_scale, [60, 95]);
  // The CLI's leaf is the drawing's Wrap Rows value (PD-OI-113).
  assert.equal(linear.form.linear_wrap_length, 5000);
  const saved = await savedDocument('released-flat-50', () => stashed);
  assert.equal(saved.version, 51);
  assert.deepEqual(saved.drawings.map((entry) => entry.mode), ['linear']);
  const [drawing] = saved.drawings;
  assert.equal(drawing.config.adv.italic_gene_labels, true);
  assert.deepEqual(drawing.config.adv.identity_scale, [60, 95]);
  assert.equal(requestLeaf(drawing, 'labels.linear.italic_gene_labels'), true);
  assert.deepEqual(requestLeaf(drawing, 'objects.blast_match.identity_scale'), [60, 95]);
  assert.equal(requestLeaf(drawing, 'canvas.linear.wrap_length'), 5000);
  assert.equal(drawing.config.form.linear_wrap_length, 5000);
});

// PD-OI-113: tag v0.14.5's Web had no Wrap Rows, so its Save kept a wrap
// length as a preserved override of the Linear draft. Load moves it into the
// drawing's Wrap Rows value; Save writes it there and in the request, and no
// preserved override holds it.
test('a wrap length a v0.14.5 Web Save kept as a preserved override loads into Wrap Rows', async () => {
  const session = JSON.parse(gunzipSync(readFileSync(new URL(
    '../fixtures/sessions/web-linear-wrap-preserved.v50.gbdraw-session.json.gz', import.meta.url
  ))));
  assert.equal(session.version, 50);
  assert.deepEqual(session.modes.linear.config.unmanagedConfigOverrides, { 'canvas.linear.wrap_length': 5000 });
  assert.equal(Object.hasOwn(session.modes.linear.config.form, 'linear_wrap_length'), false);
  const stashed = await loadKeepingSlots(session);
  const linear = state.drawingOf('linear');
  assert.equal(linear.form.linear_wrap_length, 5000);
  assert.deepEqual({ ...linear.unmanagedConfigOverrides }, {});
  const saved = await savedDocument('web-wrap-preserved', () => stashed);
  const drawing = saved.drawings.find((entry) => entry.id === 'linear');
  assert.equal(drawing.config.form.linear_wrap_length, 5000);
  assert.deepEqual(drawing.config.unmanagedConfigOverrides ?? {}, {});
  assert.equal(requestLeaf(drawing, 'canvas.linear.wrap_length'), 5000);
});

test('a Session 50 that lists drawings (never released) is refused with a message', async () => {
  const saved = await savedDocument('combined-50', () => []);
  const loaded = await importSession({
    target: { files: [new File([JSON.stringify({ ...saved, version: 50 })], 'session.json', { type: 'application/json' })], value: 'selected' }
  }, { confirmLoadRepairs: async () => true, installLoadedArtifactSlots: () => {} });
  assert.notEqual(loaded.status, 'ok');
  console.log('old combined 50:', JSON.stringify(loaded.error));
  assert.equal(loaded.error?.context?.reason, 'FIELDS');
});

// OV-2493: a drawing a Session 51 Load does not open takes its request's
// legend and title sides, as the drawing it opens does. OV-2491 (Owner
// decision R17-38): a Linear request that states no `output` reads with the
// sides Python draws it with (`ABSENT_LINEAR_OUTPUT`, generated from
// `LinearOutputOptions()`), opened or not; a Python-saved Linear Session
// states none.
test('each drawing of a Session 51 Load takes its request\'s layout, a Linear one without output Python\'s', async () => {
  const { ABSENT_LINEAR_OUTPUT } = await import('../../gbdraw/web/js/contracts/schema.generated.js');
  const linearLayout = async (session) => {
    await loadKeepingSlots(session);
    return { ...state.drawingOf('linear').layoutPreferences.linear };
  };
  for (const opened of ['circular', 'linear']) {
    const stated = v49Session();
    stated.ui.activeDrawingId = opened;
    assert.deepEqual(await linearLayout(stated), { legend: 'left', plotTitlePosition: 'bottom' }, `opened ${opened}`);
    const absent = v49Session();
    absent.ui.activeDrawingId = opened;
    delete absent.drawings.find((entry) => entry.id === 'linear').renderRequest.diagramOptions.output;
    assert.deepEqual(await linearLayout(absent), ABSENT_LINEAR_OUTPUT, `opened ${opened}, no output`);
  }
  const python = JSON.parse(gunzipSync(readFileSync(new URL(
    '../fixtures/sessions/python-linear-absent-output.v51.gbdraw-session.json.gz', import.meta.url
  ))));
  assert.equal(Object.hasOwn(python.drawings[0].renderRequest.diagramOptions, 'output'), false);
  assert.deepEqual(await linearLayout(python), ABSENT_LINEAR_OUTPUT);
});

// Owner Q1 (R17-39, A): a Linear legend the mode cannot use (`upper_left`,
// written by hand) stays `bottom`, and the pre-Load dialog lists the repair,
// even when the drawing's request states no `output` (which reads `right`).
test('a stored Linear legend outside its mode reads bottom, not the absent output\'s right, and is listed', async () => {
  for (const opened of ['circular', 'linear']) {
    const session = v49Session();
    session.ui.activeDrawingId = opened;
    const entry = session.drawings.find((drawing) => drawing.id === 'linear');
    delete entry.renderRequest.diagramOptions.output;
    entry.ui.layoutPreferences = { legend: 'upper_left', plotTitlePosition: 'bottom' };
    let dialog = null;
    const loaded = await importSession({
      target: { files: [new File([JSON.stringify(session)], 'session.json', { type: 'application/json' })], value: 'selected' }
    }, { confirmLoadRepairs: async (asked) => { dialog = asked; return true; }, installLoadedArtifactSlots: () => {} });
    assert.equal(loaded.status, 'ok', JSON.stringify(loaded.error));
    assert.deepEqual({ ...state.drawingOf('linear').layoutPreferences.linear }, { legend: 'bottom', plotTitlePosition: 'bottom' }, `opened ${opened}`);
    assert.deepEqual(dialog?.items, [
      'Legend: the position "upper_left" (drawings[1].ui.layoutPreferences.legend) is not one of its mode\'s positions; bottom is used.'
    ], `opened ${opened}`);
  }
});
