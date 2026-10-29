// PD-OI-104 scenario revision 5, PD-OI-103 rev 2, PD-OI-102: the GC content,
// GC skew and Depth toggles, the slot switches that hide or show one of those
// tracks, and "Use custom stack" go live (the GT set). A live apply writes the
// committed request with only the GT leaves replaced, so that it equals a
// Generate of the same draft (R3, R1 (c)); a show re-inserts the parts the
// last Generate built (the kept track parts), and every pending value of a
// row (a width, a colour, a track type) and every other pending setting
// stays out (PD-OI-103 rev 2).
import assert from 'node:assert/strict';
import { test } from 'node:test';
import { modeDrawing } from './helpers/drawing-state.mjs';

globalThis.window = {
  Vue: {
    ref: (value) => ({ value }), reactive: (value) => value,
    computed: (getter) => ({ get value() { return getter(); } }),
    nextTick: async () => {}
  },
  DOMPurify: { sanitize: (value) => value }
};

const { state } = await import('../../gbdraw/web/js/state.js');
const config = await import('../../gbdraw/web/js/services/config.js');
const request = await import('../../gbdraw/web/js/services/session-request.js');
const owner = await import('../../gbdraw/web/js/services/setting-application.js');
const legendRows = await import('../../gbdraw/web/js/services/track-legend-rows.js');
const { getResourcePayloadOwner, resourceIdOwnedBy, setResourcePayloadOwner } = await import('../../gbdraw/web/js/services/resource-payload-owner.js');
const { resolveLinearComparisonPlan } = await import('../../gbdraw/web/js/services/linear-comparisons.js');
const { createDefaultLinearTrackSlots } = await import('../../gbdraw/web/js/services/linear-track-slot-model.js');
const { applyCircularSuppressControlsToSlots, createDefaultCircularTrackSlots } = await import('../../gbdraw/web/js/services/circular-track-slot-model.js');
const { createAnnotationSet } = await import('../../gbdraw/web/js/services/annotation-state.js');
const validation = await import('../../gbdraw/web/js/services/track-slot-validation.js');
// The request projections a Generate runs (`canonicalRequestProjectionCount`).
let projections = 0;
globalThis.__GBDRAW_TEST_HOOKS__ = { onStructuralMetric: ({ name }) => { if (name === 'canonicalRequestProjectionCount') projections += 1; } };

const K = Object.freeze({
  GC: 'diagramOptions.configOverrides.canvas.show_gc',
  SKEW: 'diagramOptions.configOverrides.canvas.show_skew',
  DEPTH: 'diagramOptions.depthTracks',
  CSLOTS: 'diagramOptions.tracks.circularTrackSlots',
  LSLOTS: 'diagramOptions.tracks.linearTrackSlots'
});
const GT = Object.values(K);

const record = (id) => `LOCUS       ${id.padEnd(16)}          40 bp    DNA     linear   UNK 01-JAN-1980
DEFINITION  GT fixture ${id}.
ACCESSION   ${id}
VERSION     ${id}
FEATURES             Location/Qualifiers
     CDS             1..30
                     /product="test protein"
ORIGIN
        1 atgcatgcat gcatgcatgc atgcatgcat gcatgcatgc
//
`;
const textFile = (name, text) => ({
  name, type: 'text/plain', size: new TextEncoder().encode(text).byteLength, lastModified: 0, data: btoa(text)
});
const genbank = (ids) => textFile('gt.gb', ids.map(record).join(''));
const depthText = (id, seed) => Array.from({ length: 40 }, (_, index) => `${id}\t${index + 1}\t${(index * seed) % 17}`).join('\n');
const depthFile = (id, seed) => textFile(`depth-${seed}.tsv`, depthText(id, seed));
const seq = (uid, id) => ({
  uid, gb: genbank([id]), depth: null, losat_gencode: 1, region_record_id: '', region_start: null, region_end: null,
  region_reverse: false, definition: '', record_subtitle: '', file_definition: `File ${id}`, file_subtitle: ''
});

// The draft inputs: the Linear cards and the Circular Depth files.
const inputs = { linearSeqs: [seq('record-1', 'GT1'), seq('record-2', 'GT2')], circularDepth: null };
const circularGenbank = genbank(['GT1']);
const DEPTH_A = depthFile('GT1', 3);
const DEPTH_B = depthFile('GT1', 5);

const restore = [];
const set = (object, key, value) => {
  const had = Object.hasOwn(object, key);
  const old = object[key];
  restore.push(() => { if (had) object[key] = old; else delete object[key]; });
  object[key] = value;
};
const undoAll = () => { while (restore.length) restore.pop()(); };

const draftDepth = (mode) => (mode === 'circular'
  ? inputs.circularDepth
  : inputs.linearSeqs.map((card) => card.depth));
// The files a Generate of `mode` sends now (serializeActiveRenderFiles): the
// Depth files when the draft draws Depth, else the hidden ones (F1).
const filesOf = (mode, { hiddenDepth = true } = {}) => {
  const drawing = modeDrawing(state, mode);
  const draws = legendRows.draftDrawsDepth(mode, drawing);
  const hidden = !draws && hiddenDepth ? { hiddenDepth: mode === 'circular' ? { circular: inputs.circularDepth } : { linear: draftDepth(mode) } } : {};
  return mode === 'circular'
    ? { c_gb: circularGenbank, c_depth: draws ? inputs.circularDepth : undefined, linearSeqs: [], ...hidden }
    : { linearSeqs: inputs.linearSeqs.map((card) => ({ ...card, depth: draws ? card.depth : null })), linearComparisons: [], ...hidden };
};
// The committed Session a Generate of `mode` builds from that mode's drawing now.
const generate = (mode, { keepTrackParts = true, hiddenDepth = true } = {}) => {
  const drawing = modeDrawing(state, mode);
  config.installShownDrawing(mode);
  state.cInputType.value = 'gb';
  state.lInputType.value = 'gb';
  const filesData = filesOf(mode, { hiddenDepth });
  const comparisonPlanSnapshot = mode === 'linear'
    ? resolveLinearComparisonPlan({
        plan: drawing.linearComparisonPlan, sequences: filesData.linearSeqs, layout: [],
        losatProgram: drawing.losatProgram.value, blastpMode: drawing.losat.blastp.mode
      })
    : null;
  const built = request.buildCanonicalRenderRequest({ state, drawing, filesData, comparisonPlanSnapshot, keepTrackParts });
  return { renderRequest: built.renderRequest, resources: built.resources, webFiles: built.webFiles || {} };
};
const project = (committed, mode, requestPaths = GT) => request.projectLiveSettingWrites({
  committed, drawing: modeDrawing(state, mode), requestPaths, inputs
});
// The committed Session after an apply, as the commit writes it.
const applied = (committed, writes) => {
  const renderRequest = structuredClone(committed.renderRequest);
  request.applyLiveSettingWrites(renderRequest, writes);
  const next = {
    ...committed, renderRequest,
    resources: request.applyLiveResourceWrites(committed.resources, writes),
    webFiles: request.applyLiveWebFileWrites(committed.webFiles, writes)
  };
  request.carryKeptTrackParts(committed.renderRequest, next.renderRequest, writes);
  return next;
};
const live = (committed, mode, requestPaths = GT) => applied(committed, project(committed, mode, requestPaths));

// Draft values that apply on Generate and stay pending through every case:
// a Depth colour and tick interval (G3), the track type, a row width.
const pendingValues = (mode) => {
  const drawing = modeDrawing(state, mode);
  set(drawing.adv, 'label_font_size', 29);
  set(drawing.adv, 'depth_tracks', (drawing.adv.depth_tracks || []).map((track) => ({ ...track })));
  if (drawing.adv.depth_tracks[0]) {
    set(drawing.adv.depth_tracks[0], 'color', '#123456');
    set(drawing.adv.depth_tracks[0], 'large_tick_interval', 7);
  }
  if (mode === 'circular') set(drawing.form, 'track_type', drawing.form.track_type === 'middle' ? 'spreadout' : 'middle');
  else set(drawing.form, 'linear_track_layout', drawing.form.linear_track_layout === 'above' ? 'middle' : 'above');
};
// Live equals Generate of the same draft without its pending values: the
// request, the resource ids, and every Depth resource by payload identity.
const assertLiveEqualsGenerate = (next, mode, revertPending, message) => {
  revertPending();
  const expected = generate(mode, { keepTrackParts: false });
  assert.deepEqual(next.renderRequest, expected.renderRequest, message);
  assert.deepEqual(Object.keys(next.resources).sort(), Object.keys(expected.resources).sort(), `${message}: resources`);
  for (const [id, descriptor] of Object.entries(expected.resources)) {
    if (descriptor.kind !== 'depth-track-file') continue;
    assert.deepEqual({ ...next.resources[id] }, { ...descriptor }, `${message}: ${id}`);
  }
  assert.deepEqual(next.webFiles.resourceOriginalNames || {}, expected.webFiles.resourceOriginalNames || {}, `${message}: names`);
};
// Runs `body` with the pending values set, then checks live = Generate.
const withPending = (mode, body) => {
  const mark = restore.length;
  pendingValues(mode);
  const next = body();
  const revert = () => { while (restore.length > mark) restore.pop()(); };
  return { next, revert };
};

const circularPlain = () => {
  const drawing = modeDrawing(state, 'circular');
  set(drawing.form, 'multi_record_canvas', false);
  set(drawing.adv, 'circular_track_slots_enabled', false);
  set(drawing.form, 'suppress_gc', false);
  set(drawing.form, 'suppress_skew', false);
  set(drawing.form, 'show_depth', false);
  set(inputs, 'circularDepth', null);
  return drawing;
};
const circularStack = ({ depthSeries = 2, extra = [] } = {}) => {
  const drawing = circularPlain();
  set(inputs, 'circularDepth', [DEPTH_A, DEPTH_B].slice(0, depthSeries));
  set(drawing.form, 'show_depth', true);
  const slots = createDefaultCircularTrackSlots({ showDepth: true, depthTrackCount: depthSeries });
  slots.filter((slot) => slot.renderer === 'depth').forEach((slot, index) => { slot.params = { ...slot.params, track_index: index }; });
  set(drawing.adv, 'circular_track_slots', [...slots, ...extra]);
  set(drawing.adv, 'circular_track_slots_enabled', true);
  set(drawing.adv, 'depth_tracks', Array.from({ length: depthSeries }, (_, index) => ({ label: `Depth ${index + 1}` })));
  return drawing;
};
const linearSimple = () => {
  const drawing = modeDrawing(state, 'linear');
  set(drawing.adv, 'linear_track_slots_enabled', false);
  set(drawing.form, 'show_gc', true);
  set(drawing.form, 'show_skew', true);
  set(drawing.form, 'show_depth', false);
  inputs.linearSeqs.forEach((card) => set(card, 'depth', null));
  return drawing;
};
const linearStack = ({ depthSeries = 2 } = {}) => {
  const drawing = linearSimple();
  inputs.linearSeqs.forEach((card, index) => set(card, 'depth', [DEPTH_A, DEPTH_B].slice(0, depthSeries).map(
    (file, series) => depthFile(`GT${index + 1}`, 3 + 2 * series)
  )));
  set(drawing.form, 'show_depth', true);
  const slots = createDefaultLinearTrackSlots({ showGc: true, showSkew: true, showDepth: true, depthTrackCount: depthSeries });
  slots.filter((slot) => slot.renderer === 'depth').forEach((slot, index) => { slot.params = { ...slot.params, track_index: index }; });
  set(drawing.adv, 'linear_track_slots', slots);
  set(drawing.adv, 'linear_track_slots_enabled', true);
  set(drawing.adv, 'depth_tracks', Array.from({ length: depthSeries }, (_, index) => ({ label: `Depth ${index + 1}` })));
  return drawing;
};
// Hide GC Content / Hide GC Skew as their control sets them: the form, and in
// a Circular stack the rows it marks (`applyCircularSuppressControlsToSlots`).
const setSuppress = (drawing, key, value) => {
  set(drawing.form, key, value);
  if (drawing.adv.circular_track_slots_enabled) {
    set(drawing.adv, 'circular_track_slots', applyCircularSuppressControlsToSlots(drawing.adv.circular_track_slots, drawing.form, drawing.adv.nt));
  }
};
const rowOf = (drawing, mode, predicate) => (mode === 'circular' ? drawing.adv.circular_track_slots : drawing.adv.linear_track_slots)
  .find(predicate);

test('N16: track-visibility is a live row whose draft paths request the five GT keys', () => {
  const row = owner.SETTING_APPLICATION.find((candidate) => candidate.id === 'track-visibility');
  assert.equal(row.application, 'live');
  assert.equal(row.note, undefined);
  assert.ok(row.partialNote);
  const step = (...path) => owner.liveSettingRequestPaths([{ path }]).sort();
  for (const path of [['config', 'form', 'show_gc'], ['config', 'form', 'suppress_skew'], ['config', 'form', 'show_depth'],
    ['files', 'c_depth'], ['files', 'linearSeqs', 0, 'depth'], ['config', 'adv', 'circular_track_slots', 2, 'enabled'],
    ['config', 'adv', 'linear_track_slots_enabled']]) {
    assert.deepEqual(step(...path).filter((key) => GT.includes(key)), [...GT].sort(), path.join('.'));
  }
  for (const key of GT) {
    assert.ok(owner.allLiveSettingRequestPaths().includes(key), key);
    assert.ok(request.LIVE_SETTING_REQUEST_PATHS.includes(key), key);
  }
  // A G2 Height still requests its own key (and the GT keys by prefix).
  assert.ok(step('config', 'adv', 'linear_track_slots', 1, 'height').includes(K.LSLOTS));
  // depth-sources and track-slots stay on Generate with their notes.
  for (const id of ['depth-sources', 'track-slots']) {
    const other = owner.SETTING_APPLICATION.find((candidate) => candidate.id === id);
    assert.equal(other.application, 'generate', id);
    assert.ok(other.note, id);
  }
});

test('N1: the flags hide and show live in Circular plain and Linear simple, equal to Generate', () => {
  for (const [mode, setup, flag] of [['circular', circularPlain, 'suppress_gc'], ['circular', circularPlain, 'suppress_skew'],
    ['linear', linearSimple, 'show_gc'], ['linear', linearSimple, 'show_skew']]) {
    const drawing = setup();
    const committed = generate(mode);
    set(drawing.form, flag, !drawing.form[flag]);
    const hidden = withPending(mode, () => live(committed, mode));
    assertLiveEqualsGenerate(hidden.next, mode, hidden.revert, `${mode} ${flag} hide`);
    set(drawing.form, flag, !drawing.form[flag]);
    const shown = withPending(mode, () => live(hidden.next, mode));
    assertLiveEqualsGenerate(shown.next, mode, shown.revert, `${mode} ${flag} show`);
    undoAll();
  }
});

test('M2: right after a Generate the GT keys write nothing; also after a Load (no kept parts)', () => {
  for (const [mode, setup] of [['circular', circularStack], ['linear', linearStack], ['circular', circularPlain]]) {
    setup();
    const committed = generate(mode);
    assert.deepEqual(project(committed, mode), [], `${mode} after Generate`);
    const loaded = structuredClone(committed);
    assert.equal(request.keptTrackPartsOf(loaded.renderRequest), null, 'a cloned (loaded) request has no kept parts');
    assert.deepEqual(project({ ...committed, renderRequest: loaded.renderRequest }, mode), [], `${mode} after Load`);
    undoAll();
  }
});

test('M4: keepTrackParts changes nothing Generate sends, also with hidden Depth files', () => {
  for (const [mode, setup] of [['circular', circularStack], ['linear', linearStack], ['circular', circularPlain]]) {
    setup();
    if (mode === 'circular') set(inputs, 'circularDepth', [DEPTH_A]);
    const drawing = modeDrawing(state, mode);
    set(drawing.form, 'show_depth', false);
    for (const row of (mode === 'circular' ? drawing.adv.circular_track_slots : drawing.adv.linear_track_slots) || []) {
      if (row.renderer === 'depth') set(row, 'enabled', false);
    }
    const metrics = () => projections;
    const beforeOff = metrics();
    const off = generate(mode, { keepTrackParts: false, hiddenDepth: false });
    const countOff = metrics() - beforeOff;
    const beforeOn = metrics();
    const on = generate(mode, { keepTrackParts: true, hiddenDepth: true });
    assert.equal(metrics() - beforeOn, countOff, `${mode}: one request projection`);
    assert.deepEqual(on.renderRequest, off.renderRequest, `${mode}: request`);
    assert.deepEqual(on.resources, off.resources, `${mode}: resources`);
    assert.deepEqual(on.webFiles, off.webFiles, `${mode}: web files`);
    assert.ok(request.keptTrackPartsOf(on.renderRequest), `${mode}: kept parts set`);
    assert.equal(request.keptTrackPartsOf(off.renderRequest), null, `${mode}: none without the option`);
    undoAll();
  }
});

test('N4/N5: Show Depth off removes every Depth part, and on re-inserts them from the kept parts', () => {
  for (const [mode, setup] of [['circular', circularStack], ['linear', linearStack]]) {
    const drawing = setup();
    const committed = generate(mode);
    set(drawing.form, 'show_depth', false);
    (mode === 'circular' ? drawing.adv.circular_track_slots : drawing.adv.linear_track_slots)
      .filter((row) => row.renderer === 'depth').forEach((row) => set(row, 'enabled', false));
    const hidden = withPending(mode, () => live(committed, mode));
    assert.ok(!Object.values(hidden.next.resources).some((descriptor) => descriptor.kind === 'depth-track-file'), mode);
    assertLiveEqualsGenerate(hidden.next, mode, hidden.revert, `${mode} Depth hide`);
    assert.ok(request.keptTrackPartsOf(hidden.next.renderRequest)?.depth, `${mode}: kept Depth parts follow the commit`);
    undoAll();
    setup();
    set(modeDrawing(state, mode).form, 'show_depth', false);
    (mode === 'circular' ? modeDrawing(state, mode).adv.circular_track_slots : modeDrawing(state, mode).adv.linear_track_slots)
      .filter((row) => row.renderer === 'depth').forEach((row) => set(row, 'enabled', false));
    const off = generate(mode);
    set(modeDrawing(state, mode).form, 'show_depth', true);
    (mode === 'circular' ? modeDrawing(state, mode).adv.circular_track_slots : modeDrawing(state, mode).adv.linear_track_slots)
      .filter((row) => row.renderer === 'depth').forEach((row) => set(row, 'enabled', true));
    const shown = withPending(mode, () => live(off, mode));
    assertLiveEqualsGenerate(shown.next, mode, shown.revert, `${mode} Depth show after a Generate with it off`);
    undoAll();
  }
});

test('N6: a Depth show waits for a new File or one more series, and shows the same File objects', () => {
  const drawing = circularStack({ depthSeries: 1 });
  const committed = generate('circular');
  set(drawing.form, 'show_depth', false);
  rowOf(drawing, 'circular', (row) => row.renderer === 'depth').enabled = false;
  const hidden = live(committed, 'circular');
  rowOf(drawing, 'circular', (row) => row.renderer === 'depth').enabled = true;
  set(drawing.form, 'show_depth', true);
  set(inputs, 'circularDepth', [depthFile('GT1', 3)]);
  const waits = live(hidden, 'circular');
  assert.ok(!Object.values(waits.resources).some((descriptor) => descriptor.kind === 'depth-track-file'), 'a new File waits');
  assert.equal(owner.settingApplicationNote('track-visibility', (path) => request.liveSettingNoteApplies(
    waits.renderRequest, path, drawing, inputs, waits.resources)).text,
  owner.SETTING_APPLICATION.find((row) => row.id === 'track-visibility').partialNote);
  set(inputs, 'circularDepth', [DEPTH_A]);
  const shows = live(hidden, 'circular');
  assert.equal(resourceIdOwnedBy(shows.resources, DEPTH_A, 'depth-track-file'), 'depth-tracks-1-source');
  undoAll();
});

test('N8: one of two Depth slots leaves and returns between its neighbours by id; a Depth row off at Generate turns on live', () => {
  for (const [mode, setup] of [['circular', circularStack], ['linear', linearStack]]) {
    const drawing = setup();
    const committed = generate(mode);
    const second = rowOf(drawing, mode, (row) => row.renderer === 'depth' && row.params.track_index === 1);
    set(second, 'enabled', false);
    const off = withPending(mode, () => live(committed, mode));
    assertLiveEqualsGenerate(off.next, mode, off.revert, `${mode} second Depth slot off`);
    set(second, 'enabled', true);
    const on = withPending(mode, () => live(off.next, mode));
    assertLiveEqualsGenerate(on.next, mode, on.revert, `${mode} second Depth slot on`);
    undoAll();
    const again = setup();
    set(rowOf(again, mode, (row) => row.renderer === 'depth' && row.params.track_index === 1), 'enabled', false);
    const generatedOff = generate(mode);
    set(rowOf(again, mode, (row) => row.renderer === 'depth' && row.params.track_index === 1), 'enabled', true);
    const turnedOn = withPending(mode, () => live(generatedOff, mode));
    assertLiveEqualsGenerate(turnedOn.next, mode, turnedOn.revert, `${mode} Depth row off at Generate turns on`);
    undoAll();
  }
});

test('N9/N10: a GC row\'s own switch hides and shows live in a Circular and a Linear stack', () => {
  for (const [mode, setup] of [['circular', circularStack], ['linear', linearStack]]) {
    const drawing = setup();
    const committed = generate(mode);
    const gc = rowOf(drawing, mode, (row) => row.renderer === 'dinucleotide_content');
    set(gc, 'enabled', false);
    const off = withPending(mode, () => live(committed, mode));
    assertLiveEqualsGenerate(off.next, mode, off.revert, `${mode} GC row off`);
    set(gc, 'enabled', true);
    const on = withPending(mode, () => live(off.next, mode));
    assertLiveEqualsGenerate(on.next, mode, on.revert, `${mode} GC row on`);
    undoAll();
  }
});

test('N3 (PQ-4 A): in a Circular stack Hide GC off after a live hide waits; the flag still writes', () => {
  const drawing = circularStack();
  const committed = generate('circular');
  setSuppress(drawing, 'suppress_gc', true);
  const hidden = live(committed, 'circular');
  assert.ok(!hidden.renderRequest.diagramOptions.tracks.circularTrackSlots.some((slot) => slot.renderer === 'dinucleotide_content'));
  setSuppress(drawing, 'suppress_gc', false);
  const waits = live(hidden, 'circular');
  assert.equal(waits.renderRequest.diagramOptions.configOverrides['canvas.show_gc'], true);
  assert.ok(!waits.renderRequest.diagramOptions.tracks.circularTrackSlots.some((slot) => slot.renderer === 'dinucleotide_content'));
  assert.equal(request.liveSettingNoteApplies(waits.renderRequest, K.CSLOTS, drawing, inputs, waits.resources), false);
  undoAll();
});

test('N11/N12: Use custom stack switches from the kept parts; a pending width and an added row stay out', () => {
  const drawing = linearStack();
  set(drawing.adv, 'linear_track_slots_enabled', false);
  const committed = generate('linear');
  set(drawing.adv, 'linear_track_slots_enabled', true);
  const on = withPending('linear', () => {
    const gc = rowOf(drawing, 'linear', (row) => row.renderer === 'dinucleotide_content');
    set(gc, 'height', '33');
    set(drawing.adv, 'linear_track_slots', [...drawing.adv.linear_track_slots]);
    return live(committed, 'linear');
  });
  assert.equal(request.liveSettingNoteApplies(on.next.renderRequest, K.LSLOTS, drawing, inputs, on.next.resources), true);
  on.revert();
  // The stack as placed at the last Generate, with the G2 Height pass.
  assert.ok(Array.isArray(on.next.renderRequest.diagramOptions.tracks.linearTrackSlots));
  set(drawing.adv, 'linear_track_slots_enabled', false);
  set(drawing.form, 'show_gc', false);
  const off = withPending('linear', () => live(on.next, 'linear'));
  assertLiveEqualsGenerate(off.next, 'linear', off.revert, 'Linear stack off with Show GC Content off');
  undoAll();
});

test('N12b: Use custom stack off then on after a Generate with it on returns the whole stack, its features row too', () => {
  for (const [mode, setup, key] of [['linear', linearStack, 'linear_track_slots_enabled'], ['circular', circularStack, 'circular_track_slots_enabled']]) {
    const drawing = setup();
    const committed = generate(mode);
    set(drawing.adv, key, false);
    const off = withPending(mode, () => live(committed, mode));
    assertLiveEqualsGenerate(off.next, mode, off.revert, `${mode} stack off`);
    set(drawing.adv, key, true);
    const on = withPending(mode, () => live(off.next, mode));
    assertLiveEqualsGenerate(on.next, mode, on.revert, `${mode} stack on again`);
    undoAll();
  }
});

test('N13: Depth resources go by reference through two hide/show cycles and the commit rebuilds the resource table', () => {
  const drawing = circularStack({ depthSeries: 1 });
  const committed = generate('circular');
  config.canonicalRenderArtifactOwner.restore({ committedCanonicalSession: committed, activeSessionResourceTable: null });
  for (let cycle = 0; cycle < 2; cycle += 1) {
    for (const showDepth of [false, true]) {
      set(drawing.form, 'show_depth', showDepth);
      rowOf(drawing, 'circular', (row) => row.renderer === 'depth').enabled = showDepth;
      const current = config.getCommittedCanonicalSession();
      const writes = project(current, 'circular');
      const candidate = request.projectCommittedEditorIntent
        ? request.applyLiveResourceWrites(current.resources, writes) : null;
      if (showDepth) assert.equal(resourceIdOwnedBy(candidate, DEPTH_A, 'depth-track-file'), 'depth-tracks-1-source');
      config.commitLiveSettingValues(writes);
      const next = config.getCommittedCanonicalSession();
      assert.equal(resourceIdOwnedBy(next.resources, DEPTH_A, 'depth-track-file'),
        showDepth ? 'depth-tracks-1-source' : null, `cycle ${cycle} ${showDepth}`);
    }
  }
  config.canonicalRenderArtifactOwner.restore(null);
  undoAll();
});

test('N14: kept parts follow the request object: the commit carries them, a clone (Save, Load) has none', () => {
  circularStack();
  const committed = generate('circular');
  const kept = request.keptTrackPartsOf(committed.renderRequest);
  assert.ok(kept);
  const saved = JSON.parse(JSON.stringify(committed));
  assert.equal(request.keptTrackPartsOf(saved.renderRequest), null);
  assert.ok(!JSON.stringify(saved).includes('keptTrackParts'));
  undoAll();
  // A Save after a Generate with Show Depth off holds no hidden Depth resource.
  const drawing = circularStack({ depthSeries: 1 });
  set(drawing.form, 'show_depth', false);
  rowOf(drawing, 'circular', (row) => row.renderer === 'depth').enabled = false;
  const off = generate('circular');
  assert.ok(!Object.values(off.resources).some((descriptor) => descriptor.kind === 'depth-track-file'));
  assert.ok(request.keptTrackPartsOf(off.renderRequest)?.depth);
  undoAll();
});

test('N15: one writer per leaf over the GT keys and G2\'s Linear heights', () => {
  const drawing = linearStack();
  const committed = generate('linear');
  set(drawing.form, 'show_depth', false);
  rowOf(drawing, 'linear', (row) => row.renderer === 'depth' && row.params.track_index === 0).enabled = false;
  set(rowOf(drawing, 'linear', (row) => row.renderer === 'dinucleotide_content'), 'height', '41');
  const writes = project(committed, 'linear', [...request.LIVE_SETTING_REQUEST_PATHS]);
  const leaves = writes.map(([path]) => path.join('.'));
  assert.deepEqual(leaves, [...new Set(leaves)], 'each leaf once');
  undoAll();
});

test('N2 (one hide owner): Hide GC Skew removes the GC-nt skew row and keeps an AT skew row in the Result and the Legend list', () => {
  const at = { id: 'at_skew', renderer: 'dinucleotide_skew', enabled: true, params: { nt: 'AT' } };
  const drawing = circularStack({ extra: [at] });
  const committed = generate('circular');
  setSuppress(drawing, 'suppress_skew', true);
  const hidden = withPending('circular', () => live(committed, 'circular'));
  const slots = hidden.next.renderRequest.diagramOptions.tracks.circularTrackSlots;
  assert.ok(slots.some((slot) => slot.id === 'at_skew'), 'the AT row stays');
  assert.ok(!slots.some((slot) => slot.renderer === 'dinucleotide_skew' && slot.id !== 'at_skew'));
  assertLiveEqualsGenerate(hidden.next, 'circular', hidden.revert, 'Hide GC Skew with an AT row');
  const rows = new Map([
    ['GC skew (+)', { key: 'GC skew (+)', slotId: 'gc_skew', kind: 'skew', nt: 'GC', paletteKey: 'skew_high', ownColor: null }],
    ['AT skew (+)', { key: 'AT skew (+)', slotId: 'at_skew', kind: 'skew', nt: 'AT', paletteKey: 'skew_high', ownColor: null }]
  ]);
  assert.deepEqual(legendRows.hiddenTrackLegendKeys(rows, committed.renderRequest, drawing), ['GC skew (+)']);
  undoAll();
});

test('N18 (request-less Results): nothing says Live edit and the Legend list keeps every row', () => {
  const drawing = circularPlain();
  set(drawing.form, 'suppress_gc', true);
  const rows = new Map([['GC content', { key: 'GC content', slotId: null, kind: 'content', nt: 'GC', paletteKey: 'gc_content', ownColor: null }]]);
  assert.deepEqual(legendRows.hiddenTrackLegendKeys(rows, null, drawing), []);
  config.canonicalRenderArtifactOwner.restore(null);
  assert.equal(config.committedSettingNotePredicate()(K.GC), true, 'no Result: every path can');
  undoAll();
});

test('M1: a Linear Height item note keeps G2\'s answer while a GT projection is incomplete', () => {
  const drawing = linearStack();
  const committed = generate('linear');
  // A row added since the last Generate makes the GT projection incomplete.
  set(drawing.adv, 'linear_track_slots', [...drawing.adv.linear_track_slots,
    { id: 'gc_new', renderer: 'dinucleotide_content', enabled: true, params: {} }]);
  const gcRow = committed.renderRequest.diagramOptions.tracks.linearTrackSlots.find((slot) => slot.renderer === 'dinucleotide_content');
  assert.equal(request.liveSettingNoteApplies(committed.renderRequest, K.LSLOTS, drawing, inputs, committed.resources, { slotId: gcRow.id }), true);
  assert.equal(request.liveSettingNoteApplies(committed.renderRequest, K.LSLOTS, drawing, inputs, committed.resources), false);
  undoAll();
});

test('N18 (request-less Results): a Result shown without a committed request says Applies on Generate', () => {
  config.canonicalRenderArtifactOwner.restore(null);
  const shown = state.results.value;
  state.results.value = [{ name: 'loaded' }];
  try {
    for (const key of [...GT, 'diagramOptions.output.legend', 'diagramOptions.plotTitle']) {
      assert.equal(config.committedSettingNotePredicate()(key), false, key);
    }
    assert.notEqual(owner.settingApplicationNote('track-visibility', config.committedSettingNotePredicate()).text,
      owner.LIVE_NOTE?.text ?? 'Live edit');
  } finally {
    state.results.value = shown;
  }
});

// D-199 (ii), D-214 (e): after a Load, the entry a live hide removed is kept
// and its show (or the hide's Undo) re-inserts it by slot id, live, also in
// a Circular request with slots, which does not show which Use custom stack
// switch drew it.
test('N7: after a Load a live Depth hide completes the kept parts, and its Undo shows live', () => {
  for (const [mode, setup] of [['circular', circularStack], ['linear', linearStack]]) {
    const drawing = setup({ depthSeries: 1 });
    const generated = generate(mode);
    const loaded = structuredClone(generated);
    // A Load gives the draft the files its resources hold: each loaded
    // descriptor stands for the draft File.
    Object.keys(loaded.resources).forEach((id) => setResourcePayloadOwner(loaded.resources[id], getResourcePayloadOwner(generated.resources[id])));
    assert.equal(request.keptTrackPartsOf(loaded.renderRequest), null);
    set(drawing.form, 'show_depth', false);
    rowOf(drawing, mode, (row) => row.renderer === 'depth').enabled = false;
    const hidden = live(loaded, mode);
    assert.ok(request.keptTrackPartsOf(hidden.renderRequest)?.depth, `${mode}: the hide completes the kept parts`);
    rowOf(drawing, mode, (row) => row.renderer === 'depth').enabled = true;
    set(drawing.form, 'show_depth', true);
    assert.equal(request.liveSettingNoteApplies(hidden.renderRequest, K.DEPTH, drawing, inputs, hidden.resources), true);
    const undone = withPending(mode, () => live(hidden, mode));
    assertLiveEqualsGenerate(undone.next, mode, undone.revert, `${mode}: Undo of the hide after a Load`);
    undoAll();
  }
});

test('N17: the shortcut plan with Depth forced, then filtered, equals Generate with Show Depth off', () => {
  for (const trackType of ['tuckin', 'middle', 'spreadout']) {
    for (const series of [1, 2]) {
      const drawing = circularPlain();
      set(drawing.form, 'track_type', trackType);
      set(inputs, 'circularDepth', [DEPTH_A, DEPTH_B].slice(0, series));
      set(drawing.adv, 'depth_tracks', Array.from({ length: series }, (_, index) => ({ label: `Depth ${index + 1}` })));
      const off = generate('circular');
      assert.deepEqual(project(off, 'circular'), [], `${trackType}/${series}: the kept plan filtered is the Generate`);
      set(drawing.form, 'show_depth', true);
      const shown = live(off, 'circular');
      assertLiveEqualsGenerate(shown, 'circular', () => {}, `${trackType}/${series}: Depth shown`);
      set(drawing.form, 'show_depth', false);
      assertLiveEqualsGenerate(live(shown, 'circular'), 'circular', () => {}, `${trackType}/${series}: hidden again`);
      undoAll();
    }
  }
});

test('N19: the placement owners write what the Generate placement wrote', async () => {
  const { placeCircularStackForGenerate } = await import('../../gbdraw/web/js/services/circular-track-slot-model.js');
  const { placeLinearStackForGenerate } = await import('../../gbdraw/web/js/services/linear-track-slot-model.js');
  const circular = placeCircularStackForGenerate({
    slots: [
      { id: 'features', renderer: 'features', enabled: true, params: {} },
      { id: 'gc_content', renderer: 'dinucleotide_content', enabled: true, params: {} },
      { id: 'depth_off', renderer: 'depth', enabled: false, params: {} },
      { id: 'depth_a', renderer: 'depth', enabled: true, params: {} },
      { id: 'depth_b', renderer: 'depth', enabled: true, params: {} }
    ],
    nt: 'GC', trackType: 'tuckin', axisIndex: 99, form: { suppress_gc: true, suppress_skew: false }
  });
  const byId = new Map(circular.slots.map((slot) => [slot.id, slot]));
  assert.equal(byId.get('depth_off').params.track_index, undefined, 'a Depth row off takes no ordinal');
  assert.deepEqual([byId.get('depth_a').params.track_index, byId.get('depth_b').params.track_index], [0, 1]);
  assert.equal(byId.get('gc_content').enabled, false, 'Hide GC Content marks the GC row');
  assert.ok(circular.axisIndex >= 0 && circular.axisIndex <= circular.slots.length, 'the Axis boundary is clamped');
  const linear = placeLinearStackForGenerate({
    slots: [
      { id: 'features', renderer: 'features', enabled: true, params: {} },
      { id: 'depth_a', renderer: 'depth', enabled: true, params: {} },
      { id: 'depth_err', renderer: 'depth', enabled: false, depth_binding_error: 'x', params: { track_index: 4 } },
      { id: 'depth_b', renderer: 'depth', enabled: true, params: { track_index: 1 } }
    ],
    nt: 'GC', trackLayout: 'middle', axisIndex: 99,
    depthTracks: [{ label: 'A', height: '27' }, { label: 'B' }], depthSeriesCount: 2
  });
  const lin = new Map(linear.slots.map((slot) => [slot.id, slot]));
  assert.deepEqual([lin.get('depth_a').params.track_index, lin.get('depth_b').params.track_index], [0, 1]);
  assert.equal(lin.get('depth_err').params.track_index, undefined, 'a row with a binding error loses its index');
  assert.equal(lin.get('depth_a').height, '27px', 'the Depth height comes from its series');
  assert.ok(linear.axisIndex >= 0 && linear.axisIndex <= linear.slots.length);
});

test('N20 (F1): Generate serializes the hidden Depth files only when the draft draws no Depth, each File read once', async () => {
  let reads = 0;
  const bytes = new TextEncoder().encode(depthText('GT1', 3));
  const fileLike = { name: 'depth.tsv', type: 'text/plain', lastModified: 0, arrayBuffer: async () => { reads += 1; return bytes.buffer.slice(0); } };
  const drawing = circularPlain();
  const source = { files: { c_depth: [fileLike] }, linearSeqs: [] };
  set(drawing.form, 'show_depth', true);
  const on = await config.serializeActiveRenderFiles('circular', source, drawing, { comparisonPlan: null, hiddenDepth: true });
  assert.equal(on.hiddenDepth, undefined, 'a draft that draws Depth sends it as c_depth');
  assert.ok(on.c_depth);
  set(drawing.form, 'show_depth', false);
  const off = await config.serializeActiveRenderFiles('circular', source, drawing, { comparisonPlan: null, hiddenDepth: true });
  assert.equal(off.c_depth, null);
  assert.deepEqual(off.hiddenDepth, { circular: on.c_depth });
  const plain = await config.serializeActiveRenderFiles('circular', source, drawing, { comparisonPlan: null });
  assert.equal(plain.hiddenDepth, undefined, 'only Generate asks for them');
  assert.equal(reads, 1, 'the File is read once');
  undoAll();
});

test('N13b: a commit that changes the resources rebuilds the adopted resource table', async () => {
  const { adoptCurrentSessionResources, adoptedSessionResourceDescriptor } = await import('../../gbdraw/web/js/services/session-resource-backing.js');
  const drawing = circularStack({ depthSeries: 1 });
  const committed = generate('circular');
  const table = adoptCurrentSessionResources(committed.resources);
  config.canonicalRenderArtifactOwner.restore({ committedCanonicalSession: committed, activeSessionResourceTable: table });
  set(drawing.form, 'show_depth', false);
  rowOf(drawing, 'circular', (row) => row.renderer === 'depth').enabled = false;
  config.commitLiveSettingValues(project(committed, 'circular'));
  const afterHide = config.canonicalRenderArtifactOwner.capture().activeSessionResourceTable;
  assert.notEqual(afterHide, table);
  assert.throws(() => adoptedSessionResourceDescriptor(afterHide, 'depth-tracks-1-source'), /Missing canonical resource/);
  set(drawing.form, 'show_depth', true);
  rowOf(drawing, 'circular', (row) => row.renderer === 'depth').enabled = true;
  config.commitLiveSettingValues(project(config.getCommittedCanonicalSession(), 'circular'));
  const afterShow = config.canonicalRenderArtifactOwner.capture().activeSessionResourceTable;
  assert.ok(adoptedSessionResourceDescriptor(afterShow, 'depth-tracks-1-source'));
  config.canonicalRenderArtifactOwner.restore(null);
  undoAll();
});

// D-212: a hide the switch brings (Linear stack off with Show GC Content
// off) still applies from the committed request; the rest of the switch waits.
test('N11b (Load): with no kept parts, Use custom stack off applies whole on a Linear stack (C\'); the other switches wait', () => {
  for (const [mode, setup, key, on] of [
    ['linear', linearStack, 'linear_track_slots_enabled', false],
    ['linear', linearSimple, 'linear_track_slots_enabled', true],
    ['circular', circularStack, 'circular_track_slots_enabled', false],
    ['circular', circularPlain, 'circular_track_slots_enabled', true]
  ]) {
    for (const showGc of [false, true]) {
      const label = `${mode} stack ${on ? 'on' : 'off'}, GC ${showGc ? 'shown' : 'hidden'}`;
      const drawing = setup();
      // The rows the switch turns on (a plain draft keeps its stack rows while it is off).
      if (on && mode === 'linear') set(drawing.adv, 'linear_track_slots', createDefaultLinearTrackSlots({ showGc: true, showSkew: true }));
      if (on && mode === 'circular') set(drawing.adv, 'circular_track_slots', createDefaultCircularTrackSlots({ showDepth: false }));
      if (mode === 'linear') set(drawing.form, 'show_gc', showGc);
      else setSuppress(drawing, 'suppress_gc', !showGc);
      const committed = generate(mode);
      const loaded = { ...committed, renderRequest: structuredClone(committed.renderRequest) };
      assert.equal(request.keptTrackPartsOf(loaded.renderRequest), null);
      set(drawing.adv, key, on);
      const applies = mode === 'linear' && !on;
      if (applies) {
        const off = withPending(mode, () => live(loaded, mode));
        assertLiveEqualsGenerate(off.next, mode, off.revert, `${label}: the switch applies whole`);
      } else {
        assert.deepEqual(project(loaded, mode), [], `${label}: writes nothing`);
      }
      for (const path of [K.GC, mode === 'circular' ? K.CSLOTS : K.LSLOTS]) {
        assert.equal(request.liveSettingNoteApplies(loaded.renderRequest, path, drawing, inputs, loaded.resources), applies,
          `${label}: ${path}`);
      }
      undoAll();
    }
  }
  // The switch the request was drawn with: a hide after the Load stays live.
  const drawing = linearStack();
  const committed = generate('linear');
  const loaded = { ...committed, renderRequest: structuredClone(committed.renderRequest) };
  set(rowOf(drawing, 'linear', (row) => row.renderer === 'dinucleotide_content'), 'enabled', false);
  const hidden = withPending('linear', () => live(loaded, 'linear'));
  assertLiveEqualsGenerate(hidden.next, 'linear', hidden.revert, 'Linear GC row off after a Load');
  undoAll();
});

test('H2: committed slots stated as spec strings are read through the slot parser and stay strings', () => {
  for (const [mode, setup, key, slotsKey, axisKey, hide] of [
    ['circular', circularPlain, 'circular_track_slots_enabled', 'circularTrackSlots', 'circularTrackAxisIndex',
      (drawing) => setSuppress(drawing, 'suppress_gc', true)],
    ['linear', linearSimple, 'linear_track_slots_enabled', 'linearTrackSlots', 'linearTrackAxisIndex',
      (drawing) => set(rowOf(drawing, 'linear', (row) => row.renderer === 'dinucleotide_content'), 'enabled', false)]
  ]) {
    const drawing = setup();
    set(drawing.adv, mode === 'circular' ? 'circular_track_slots' : 'linear_track_slots', mode === 'circular'
      ? createDefaultCircularTrackSlots({ showDepth: false }) : createDefaultLinearTrackSlots({ showGc: true, showSkew: true }));
    set(drawing.adv, key, true);
    const committed = generate(mode);
    const renderRequest = structuredClone(committed.renderRequest);
    const specs = renderRequest.diagramOptions.tracks[slotsKey].map((slot) => `${slot.id}:${slot.renderer}`);
    renderRequest.diagramOptions.tracks[slotsKey] = specs;
    const loaded = { ...committed, renderRequest };
    assert.deepEqual(project(loaded, mode), [], `${mode}: a no-op writes nothing`);
    hide(drawing);
    const next = applied(loaded, project(loaded, mode));
    const expected = generate(mode, { keepTrackParts: false }).renderRequest.diagramOptions.tracks;
    const written = next.renderRequest.diagramOptions.tracks;
    assert.deepEqual(written[slotsKey], specs.filter((spec) => spec.split(':')[1] !== 'dinucleotide_content'),
      `${mode}: the kept rows stay the committed strings`);
    assert.deepEqual(written[slotsKey].map((spec) => spec.split(':')[0]), expected[slotsKey].map((slot) => slot.id), `${mode}: ids`);
    assert.equal(written[axisKey], expected[axisKey], `${mode}: axis`);
    // The hide completed the kept parts from the loaded plan; nothing waits (D-214 (e)).
    assert.equal(request.liveSettingNoteApplies(next.renderRequest, mode === 'circular' ? K.CSLOTS : K.LSLOTS, drawing, inputs, next.resources),
      true);
    undoAll();
  }
});

test('A committed null Axis index (a Python-written request: each slot its own side) stays null through a hide', () => {
  const drawing = circularStack({ depthSeries: 1 });
  const committed = generate('circular');
  const renderRequest = structuredClone(committed.renderRequest);
  renderRequest.diagramOptions.tracks.circularTrackAxisIndex = null;
  const loaded = { ...committed, renderRequest };
  setSuppress(drawing, 'suppress_gc', true);
  const next = applied(loaded, project(loaded, 'circular'));
  const tracks = next.renderRequest.diagramOptions.tracks;
  assert.ok(!tracks.circularTrackSlots.some((slot) => slot.renderer === 'dinucleotide_content'), 'the hide applies');
  assert.equal(tracks.circularTrackAxisIndex, null);
  undoAll();
});

test('L6: one note render runs the track projection once for the five GT keys', () => {
  circularStack();
  const committed = generate('circular');
  config.canonicalRenderArtifactOwner.restore({ committedCanonicalSession: committed, activeSessionResourceTable: null });
  const hooks = globalThis.__GBDRAW_TEST_HOOKS__;
  const previous = hooks.onStructuralMetric;
  let runs = 0;
  hooks.onStructuralMetric = (event) => { if (event.name === 'liveTrackPlanProjectionCount') runs += 1; previous(event); };
  try {
    const applies = config.committedSettingNotePredicate
      ? config.committedSettingNotePredicate() : (path) => config.committedLiveSettingApplies(path);
    assert.equal(owner.settingApplicationNote('track-visibility', applies).text, owner.settingApplicationNote('track-visibility').text);
    assert.equal(runs, 1);
    // The next render reads the draft again.
    owner.settingApplicationNote('track-visibility', config.committedSettingNotePredicate());
    assert.equal(runs, 2);
  } finally {
    hooks.onStructuralMetric = previous;
    config.canonicalRenderArtifactOwner.restore(null);
    undoAll();
  }
});

// Review 2 (D-212): a hide always applies from the committed request; the
// Legend list hides a track's rows only when the apply removes the track.
// Whether a request draws a track: a slot track by its slot id, a default
// content or skew track by its flag, a default Depth track by its series.
const resultDraws = (renderRequest, track) => {
  if (track.slotId !== null) return legendRows.drawnTracks(renderRequest).some((drawn) => drawn.slotId === track.slotId);
  if (track.kind === 'depth') return Array.isArray(renderRequest.diagramOptions.depthTracks);
  return renderRequest.diagramOptions.configOverrides?.[track.kind === 'content' ? 'canvas.show_gc' : 'canvas.show_skew'] !== false;
};
// The Legend list (`hiddenTrackLegendKeys` on the committed request, with the
// drawing the apply judges hides with, as the app reads it) hides exactly the
// rows of the tracks the live apply removes; returns the apply.
const assertListFollowsApply = (committed, mode, label) => {
  const drawing = modeDrawing(state, mode);
  const tracks = legendRows.drawnTracks(committed.renderRequest);
  const rows = new Map(tracks.map((track, index) => [`row ${index}`, { key: `row ${index}`, ...track }]));
  const judged = request.liveTrackHideDrawing(committed.renderRequest, drawing, inputs, committed.resources);
  const listedHidden = new Set(legendRows.hiddenTrackLegendKeys(rows, committed.renderRequest, judged));
  const next = live(committed, mode);
  for (const [key, track] of rows) {
    assert.equal(listedHidden.has(key), !resultDraws(next.renderRequest, track), `${label}: the list follows the apply for ${track.kind} ${track.slotId}`);
  }
  return next;
};
const loadOf = (committed) => ({ ...committed, renderRequest: structuredClone(committed.renderRequest) });
const drawsKind = (renderRequest, renderer) => (renderRequest.diagramOptions.tracks.circularTrackSlots
  ?? renderRequest.diagramOptions.tracks.linearTrackSlots ?? []).some((slot) => slot.renderer === renderer);

test('H1 (review 2 P1): after a Load of a Circular Session drawn with track widths, Hide GC Content applies; also without the Load', () => {
  for (const loaded of [true, false]) {
    const label = loaded ? 'after a Load' : 'after a Generate';
    const drawing = circularPlain();
    set(drawing.adv, 'gc_content_width_circular', 40);
    const generated = generate('circular');
    const committed = loaded ? loadOf(generated) : generated;
    assert.ok(drawsKind(committed.renderRequest, 'dinucleotide_content'), `${label}: the widths drew slots`);
    setSuppress(drawing, 'suppress_gc', true);
    const next = assertListFollowsApply(committed, 'circular', label);
    assert.ok(!drawsKind(next.renderRequest, 'dinucleotide_content'), `${label}: the GC content slot leaves`);
    const hidden = withPending('circular', () => live(committed, 'circular'));
    assertLiveEqualsGenerate(hidden.next, 'circular', hidden.revert, `${label}: Hide GC Content`);
    undoAll();
  }
});

test('H1 (review 2): a hide after a Load applies with a pending width, and after a pending Use custom stack switch', () => {
  // A Circular stack Session; the draft gets a Circular track width (pending) after the Load.
  {
    const drawing = circularStack({ depthSeries: 1 });
    const committed = loadOf(generate('circular'));
    set(drawing.adv, 'gc_content_width_circular', 40);
    setSuppress(drawing, 'suppress_gc', true);
    const next = assertListFollowsApply(committed, 'circular', 'Circular stack, pending width');
    assert.ok(!drawsKind(next.renderRequest, 'dinucleotide_content'), 'the GC content slot leaves');
    assert.ok(drawsKind(next.renderRequest, 'dinucleotide_skew') && drawsKind(next.renderRequest, 'depth'), 'the other tracks stay');
    undoAll();
  }
  // A Linear stack Session; Use custom stack turned on in a Linear simple
  // Session after the Load waits whole (no plan of the stack), and Show GC
  // Content off then hides the default GC track (review 3 moved the stack-off
  // case to C', where the switch applies).
  {
    const drawing = linearSimple();
    set(drawing.adv, 'linear_track_slots', createDefaultLinearTrackSlots({ showGc: true, showSkew: true }));
    const committed = loadOf(generate('linear'));
    set(drawing.adv, 'linear_track_slots_enabled', true);
    assert.deepEqual(project(committed, 'linear'), [], 'the switch alone writes nothing');
    set(drawing.form, 'show_gc', false);
    const next = assertListFollowsApply(committed, 'linear', 'Linear stack on (waits), then Show GC Content off');
    assert.equal(next.renderRequest.diagramOptions.configOverrides['canvas.show_gc'], false, 'the GC content track leaves');
    assert.equal(next.renderRequest.diagramOptions.tracks.linearTrackSlots, null, 'no stack');
    for (const path of [K.GC, K.LSLOTS]) {
      assert.equal(request.liveSettingNoteApplies(next.renderRequest, path, drawing, inputs, next.resources), false,
        `${path}: the switch still says Applies on Generate`);
    }
    undoAll();
  }
});

test('M1 (review 2 P2): the drawn switch is never read from the draft widths; a pending width cleared and the stack turned on wait', () => {
  const drawing = circularPlain();
  set(drawing.adv, 'gc_content_width_circular', 40);
  set(drawing.adv, 'circular_track_slots', createDefaultCircularTrackSlots({ showDepth: false }));
  const committed = loadOf(generate('circular'));
  set(drawing.adv, 'gc_content_width_circular', null);
  set(drawing.adv, 'circular_track_slots_enabled', true);
  assert.deepEqual(project(committed, 'circular'), [], 'nothing is written for the switch');
  for (const path of [K.GC, K.CSLOTS]) {
    assert.equal(request.liveSettingNoteApplies(committed.renderRequest, path, drawing, inputs, committed.resources), false,
      `${path} says Applies on Generate`);
  }
  undoAll();
});

test('M2 (review 2): while a show waits the Legend list lists no row of the track and hides none', () => {
  // PQ-4 A: Hide GC Content off after a live hide in a Circular stack waits.
  {
    const drawing = circularStack();
    const committed = generate('circular');
    setSuppress(drawing, 'suppress_gc', true);
    const hidden = assertListFollowsApply(committed, 'circular', 'Hide GC Content');
    setSuppress(drawing, 'suppress_gc', false);
    const waits = assertListFollowsApply(hidden, 'circular', 'Hide GC Content off (waits)');
    assert.ok(!legendRows.drawnTracks(waits.renderRequest).some((track) => track.kind === 'content'), 'no GC content row is listed');
    assert.equal(request.liveSettingNoteApplies(waits.renderRequest, K.CSLOTS, drawing, inputs, waits.resources), false);
    undoAll();
  }
  // After a Load of a Circular stack (the drawn switch unknown): a Depth hide applies, and its show (D-214 (e)).
  {
    const drawing = circularStack({ depthSeries: 1 });
    const committed = loadOf(generate('circular'));
    set(drawing.form, 'show_depth', false);
    set(rowOf(drawing, 'circular', (row) => row.renderer === 'depth'), 'enabled', false);
    const hidden = assertListFollowsApply(committed, 'circular', 'Depth hide after a Load');
    assert.ok(!drawsKind(hidden.renderRequest, 'depth'));
    set(drawing.form, 'show_depth', true);
    set(rowOf(drawing, 'circular', (row) => row.renderer === 'depth'), 'enabled', true);
    const shown = assertListFollowsApply(hidden, 'circular', 'Depth show after a Load');
    assert.ok(drawsKind(shown.renderRequest, 'depth'), 'the show applies');
    assert.equal(request.liveSettingNoteApplies(shown.renderRequest, K.DEPTH, drawing, inputs, shown.resources), true);
    undoAll();
  }
});

// The error a Generate of `mode`'s draft throws now (none: the case fails).
const generateError = (mode) => {
  try {
    generate(mode, { keepTrackParts: false });
  } catch (error) {
    return error;
  }
  assert.fail('Generate refuses the draft');
};
const sameError = (expected) => (error) => {
  assert.equal(error.constructor, expected.constructor);
  assert.equal(error.message, expected.message);
  assert.deepEqual(error.issues, expected.issues);
  return true;
};

// CI round 3 (tobacco): a committed row that names a hidden track slot
// (an Overlay tint row's Extend to slot or anchor) leaves with it when the
// draft does not draw it either; while the draft draws it, the apply is
// refused as Generate refuses it (TRACK_INVALID).
const tintRow = () => ({ id: 'tint', renderer: 'annotations', enabled: true, side: 'overlay', z: -1,
  params: { set_id: 'tint_set', anchor_slot: 'features', layer: 'underlay', to_slot: 'gc_content', marks: ['highlight'] } });
// The annotation sets are the slice's array: the tint set replaces its items.
const withTintSet = (drawing) => {
  const sets = drawing.annotationSets;
  const before = sets.splice(0, sets.length);
  restore.push(() => sets.splice(0, sets.length, ...before));
  sets.push(createAnnotationSet({ id: 'tint_set', annotations: [{ id: 'span', target: { kind: 'coordinateSpan', start: 1, end: 20 }, mark: 'highlight' }] }));
};

test('A tint row that extends to the GC content slot leaves with it when its own switch is off too', () => {
  for (const loaded of [false, true]) {
    const label = loaded ? 'after a Load' : 'after a Generate';
    const drawing = circularStack({ depthSeries: 1, extra: [tintRow()] });
    withTintSet(drawing);
    const generated = generate('circular');
    const committed = loaded ? loadOf(generated) : generated;
    assert.ok(committed.renderRequest.diagramOptions.tracks.circularTrackSlots.some((slot) => slot.id === 'tint'), `${label}: the tint is drawn`);
    set(rowOf(drawing, 'circular', (row) => row.id === 'tint'), 'enabled', false);
    set(rowOf(drawing, 'circular', (row) => row.id === 'gc_content'), 'enabled', false);
    const hidden = withPending('circular', () => live(committed, 'circular'));
    const ids = hidden.next.renderRequest.diagramOptions.tracks.circularTrackSlots.map((slot) => slot.id);
    assert.ok(!ids.includes('gc_content') && !ids.includes('tint'), `${label}: ${ids}`);
    assertLiveEqualsGenerate(hidden.next, 'circular', hidden.revert, `${label}: GC row and tint off`);
    undoAll();
    // The tint still on: Generate refuses the stack, and so does the apply.
    const again = circularStack({ depthSeries: 1, extra: [tintRow()] });
    withTintSet(again);
    const base = generate('circular');
    set(rowOf(again, 'circular', (row) => row.id === 'gc_content'), 'enabled', false);
    const refusal = generateError('circular');
    assert.throws(() => project(loaded ? loadOf(base) : base, 'circular'), sameError(refusal), `${label}: the apply is refused as Generate`);
    undoAll();
  }
});

// Review 3 M1, M2: a refused apply is data. The setting note and the Legend
// list read the projection without an error (the note says the hide does not
// apply, the list keeps the track Generate would refuse to remove); only the
// writer throws, with the error Generate gives the same draft.
test('M1/M2 (review 3): Hide GC Content with a tint row extending to GC: the writer throws Generate\'s error; the note and the Legend list follow the refusal', () => {
  for (const loaded of [false, true]) {
    const label = loaded ? 'after a Load' : 'after a Generate';
    const drawing = circularStack({ depthSeries: 1, extra: [tintRow()] });
    withTintSet(drawing);
    const generated = generate('circular');
    const committed = loaded ? loadOf(generated) : generated;
    setSuppress(drawing, 'suppress_gc', true);
    const applies = (path) => request.liveSettingNoteApplies(committed.renderRequest, path, drawing, inputs, committed.resources);
    assert.equal(applies(K.GC), false, `${label}: the note says the hide does not apply`);
    const note = owner.settingApplicationNote('track-visibility', applies);
    assert.notEqual(note.text, owner.settingApplicationNote('track-visibility').text, `${label}: the note is not Live edit`);
    const refusal = generateError('circular');
    assert.throws(() => project(committed, 'circular'), sameError(refusal), `${label}: the writer throws Generate's error`);
    const judged = request.liveTrackHideDrawing(committed.renderRequest, drawing, inputs, committed.resources);
    assert.equal(judged, null, `${label}: nothing applies`);
    const tracks = legendRows.drawnTracks(committed.renderRequest);
    const rows = new Map(tracks.map((track, index) => [`row ${index}`, { key: `row ${index}`, ...track }]));
    assert.ok(legendRows.hiddenTrackLegendKeys(rows, committed.renderRequest, drawing).length > 0, `${label}: the hide answer hides GC`);
    assert.deepEqual(legendRows.hiddenTrackLegendKeys(rows, committed.renderRequest, judged), [], `${label}: the list keeps every row`);
    config.canonicalRenderArtifactOwner.restore({ committedCanonicalSession: committed, activeSessionResourceTable: null });
    try {
      assert.equal(config.committedTrackHideDrawing(drawing), null, `${label}: the app reads the refusal`);
    } finally {
      config.canonicalRenderArtifactOwner.restore(null);
    }
    undoAll();
  }
});

test('M2 (review 3): one reader of an Overlay row\'s references; a falsy value names no slot', () => {
  assert.equal(validation.overlaySlotReference({ anchor_slot: ' features ', to_slot: 0 }, 'anchor_slot'), 'features');
  assert.equal(validation.overlaySlotReference({ anchor_slot: 'features', to_slot: 0 }, 'to_slot'), '');
  // Generate's plan check and the projection's check read the same issue.
  const rows = [{ id: 'features', renderer: 'features', side: 'inside', z: 0, params: {} },
    { id: 'tint', renderer: 'annotations', side: 'overlay', z: -1,
      params: { set_id: 'tint_set', anchor_slot: 'features', to_slot: 'gc_content', layer: 'underlay' } }];
  const plan = validation.validateCustomTrackPlan({ mode: 'circular', slots: rows, axisIndex: 1, trackType: 'tuckin', annotationSetIds: ['tint_set'] });
  const issues = validation.overlayReferenceIssues(rows, 'circular', (row) => rows.indexOf(row)).map(({ issue }) => issue);
  assert.deepEqual(issues, plan.rowIssues.get(1));
});

// Review 3 R3-3, D-214 (e): after a Load of a Circular stack, a GC content
// row's own switch off then on shows the row again by slot id, equal to
// Generate; PQ-4 A stays: Hide GC Content turned off after a live hide after
// the Load waits, in a Circular stack and with Circular track widths.
test('R3-3 (review 3): after a Load a Circular row switched off then on shows live; Hide GC Content off still waits', () => {
  {
    const drawing = circularStack({ depthSeries: 1 });
    const committed = loadOf(generate('circular'));
    const row = rowOf(drawing, 'circular', (each) => each.renderer === 'dinucleotide_content');
    set(row, 'enabled', false);
    const hidden = live(committed, 'circular');
    assert.ok(!drawsKind(hidden.renderRequest, 'dinucleotide_content'), 'the row hides');
    row.enabled = true;
    assert.equal(request.liveSettingNoteApplies(hidden.renderRequest, K.CSLOTS, drawing, inputs, hidden.resources), true);
    const shown = withPending('circular', () => live(hidden, 'circular'));
    assertLiveEqualsGenerate(shown.next, 'circular', shown.revert, 'the row shows again');
    undoAll();
  }
  for (const [label, setup] of [
    ['Circular stack', () => circularStack({ depthSeries: 1 })],
    ['Circular track widths', () => { const drawing = circularPlain(); set(drawing.adv, 'gc_content_width_circular', 40); return drawing; }]
  ]) {
    const drawing = setup();
    const committed = loadOf(generate('circular'));
    setSuppress(drawing, 'suppress_gc', true);
    const hidden = live(committed, 'circular');
    assert.ok(!drawsKind(hidden.renderRequest, 'dinucleotide_content'), `${label}: Hide GC Content applies`);
    setSuppress(drawing, 'suppress_gc', false);
    assert.equal(request.liveSettingNoteApplies(hidden.renderRequest, K.CSLOTS, drawing, inputs, hidden.resources), false,
      `${label}: Hide GC Content off says Applies on Generate`);
    const after = live(hidden, 'circular');
    assert.ok(!drawsKind(after.renderRequest, 'dinucleotide_content'), `${label}: the GC content track waits (PQ-4 A)`);
    undoAll();
  }
});

// Review 3 H1, D-214 (C'): after a Load of a Linear stack Session, Use custom
// stack turned off applies whole from the committed request (no slots, the
// draft's flags, the committed Depth parts), equal to Generate, also with
// the mode's default Show GC Content and Show GC Skew (off in Linear).
test('H1 (review 3): Linear stack after a Load, Use custom stack off alone applies whole, equal to Generate', () => {
  for (const showDepth of [true, false]) {
    const label = `Show Depth ${showDepth ? 'on' : 'off'}`;
    const drawing = linearStack({ depthSeries: 1 });
    set(drawing.form, 'show_gc', false);
    set(drawing.form, 'show_skew', false);
    set(drawing.form, 'show_depth', showDepth);
    const committed = loadOf(generate('linear'));
    set(drawing.adv, 'linear_track_slots_enabled', false);
    assert.equal(request.liveSettingNoteApplies(committed.renderRequest, K.LSLOTS, drawing, inputs, committed.resources), true, `${label}: note`);
    const off = withPending('linear', () => live(committed, 'linear'));
    assert.equal(off.next.renderRequest.diagramOptions.tracks.linearTrackSlots, null, `${label}: no stack`);
    assertLiveEqualsGenerate(off.next, 'linear', off.revert, `${label}: the switch off`);
    // And back on: the loaded stack returns from the parts the switch completed.
    set(drawing.adv, 'linear_track_slots_enabled', true);
    set(drawing.form, 'show_depth', true);
    const on = withPending('linear', () => live(off.next, 'linear'));
    assertLiveEqualsGenerate(on.next, 'linear', on.revert, `${label}: the switch on again`);
    undoAll();
  }
});

// Review 3 (f), D-214: the switch cannot apply alone while the draft's
// Region Annotation sets are not the committed ones (Reset removes them; the
// removal waits for Generate): it waits whole with the hides it brings, and
// a hide of the drawn stack still applies.
test('(f) (review 3): Use custom stack waits whole while the draft annotation sets differ; a row hide still applies', () => {
  const drawing = circularStack({ depthSeries: 1, extra: [tintRow()] });
  withTintSet(drawing);
  const committed = generate('circular');
  drawing.annotationSets.splice(0, drawing.annotationSets.length);
  set(drawing.adv, 'circular_track_slots_enabled', false);
  set(drawing.form, 'show_depth', false);
  assert.deepEqual(project(committed, 'circular'), [], 'the switch and the Depth hide it brings write nothing');
  assert.equal(request.liveSettingNoteApplies(committed.renderRequest, K.CSLOTS, drawing, inputs, committed.resources), false);
  set(rowOf(drawing, 'circular', (row) => row.renderer === 'dinucleotide_skew'), 'enabled', false);
  const next = assertListFollowsApply(committed, 'circular', 'a skew row hide while the switch waits');
  const ids = next.renderRequest.diagramOptions.tracks.circularTrackSlots.map((slot) => slot.id);
  assert.ok(!drawsKind(next.renderRequest, 'dinucleotide_skew'), 'the skew row hide applies');
  assert.ok(ids.includes('tint') && drawsKind(next.renderRequest, 'depth'), `the stack stays otherwise: ${ids}`);
  undoAll();
});
