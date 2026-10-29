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
  const hidden = !draws && hiddenDepth ? { hiddenDepth: mode === 'circular' ? [inputs.circularDepth] : draftDepth(mode) } : {};
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
  assert.equal(config.committedLiveSettingApplies(K.GC), true, 'no Result: every path can');
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
      assert.equal(config.committedLiveSettingApplies(key), false, key);
    }
    assert.notEqual(owner.settingApplicationNote('track-visibility', config.committedLiveSettingApplies).text,
      owner.LIVE_NOTE?.text ?? 'Live edit');
  } finally {
    state.results.value = shown;
  }
});

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
  assert.equal(off.hiddenDepth.length, 1);
  assert.deepEqual(off.hiddenDepth[0], on.c_depth);
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
