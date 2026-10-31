// Custom stack row edits through the Circular and Linear editors (TK-04,
// TK-07, TK-09, TK-10, GX-01): a renderer change keeps only the params the new
// renderer accepts, an invalid Depth track index is not replaced silently,
// Hide GC Skew reaches only the GC rows, Reset follows the loaded Depth
// series, and Duplicate is unavailable while a Session operation runs.
import assert from 'node:assert/strict';
import { cp, mkdtemp, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';
import test from 'node:test';
import { withDrawings } from './helpers/drawing-state.mjs';
import { loadPageVue } from './helpers/vue-runtime.mjs';

const repoRoot = process.cwd();
const tempRoot = await mkdtemp(join(tmpdir(), 'gbdraw-track-slot-row-edits-'));
for (const directory of ['app', 'utils', 'contracts']) {
  await cp(join(repoRoot, 'gbdraw', 'web', 'js', directory), join(tempRoot, directory), { recursive: true });
}
// The slice editors and svg-serialization resolve a drawing's ID by its mode.
await cp(join(repoRoot, 'gbdraw', 'web', 'js', 'services', 'drawing-view.js'), join(tempRoot, 'services', 'drawing-view.js'));
await writeFile(join(tempRoot, 'package.json'), '{"type":"module"}', 'utf8');
// track-slot-display.js reads the generated size-class boundaries.
await cp(join(repoRoot, 'gbdraw', 'web', 'js', 'mode-profiles.generated.js'), join(tempRoot, 'mode-profiles.generated.js'));
// The services the slot editors import.
for (const leaf of ['conservation-series.js', 'current-option-values.js', 'depth-track-state.js', 'linear-label-visibility.js', 'losat-normalization.js',
  'track-slot-display.js', 'track-slot-validation.js', 'circular-track-slot-model.js', 'linear-track-slot-model.js']) {
  await cp(join(repoRoot, 'gbdraw', 'web', 'js', 'services', leaf), join(tempRoot, 'services', leaf));
}

const { createCircularTrackSlotEditor } = await import(pathToFileURL(join(tempRoot, 'app', 'circular-track-slots.js')));
const { createLinearTrackSlotEditor } = await import(pathToFileURL(join(tempRoot, 'app', 'linear-track-slots.js')));
const { paramsKeptOnRendererChange } = await import(pathToFileURL(join(tempRoot, 'services', 'track-slot-validation.js')));
const { askingTrackRowSwitch } = await import(pathToFileURL(join(tempRoot, 'app', 'track-slot-edits.js')));

const depthFile = { name: 'a.depth.tsv', size: 1, lastModified: 1 };

const circularState = (slots, { showDepth = false, depthFiles = [] } = {}) => withDrawings({
  mode: 'circular',
  form: { track_type: 'tuckin', show_depth: showDepth, suppress_gc: false, suppress_skew: false, show_scale: true, separate_strands: true },
  adv: {
    circular_track_slots_enabled: true,
    circular_track_slots_axis_index: 1,
    circular_track_slots: structuredClone(slots),
    nt: 'GC',
    features: ['CDS'],
    feature_shapes: { CDS: 'arrow' },
    depth_tracks: [],
    feature_width_circular: null,
    depth_width_circular: null,
    gc_content_width_circular: null,
    gc_content_radius_circular: null,
    gc_skew_width_circular: null,
    gc_skew_radius_circular: null
  },
  files: { c_depth: depthFiles, c_conservation_blasts: [], c_conservation_fastas: [] },
  circularConservation: { enabled: false, source: 'upload', labels: '', series: [] },
  annotationSets: [],
  circularRecordList: { value: [] }
});

const linearState = (slots, { showDepth = false, depthFiles = [] } = {}) => withDrawings({
  mode: 'linear', form: { linear_track_layout: 'middle', show_depth: showDepth, show_gc: true, show_skew: true },
  adv: {
    nt: 'GC',
    linear_track_slots_enabled: true,
    linear_track_slots_axis_index: 0,
    linear_track_slots: structuredClone(slots),
    depth_tracks: [],
    depth_height: null,
    gc_height: null
  },
  linearSeqs: [{ depth: depthFiles }],
  annotationSets: []
});

const circularRow = (id, renderer, params = {}, side = 'inside') => ({ id, renderer, enabled: true, side, z: 0, params });
const linearRow = (id, renderer, params = {}, side = 'below') => ({ id, renderer, enabled: true, side, height: '', spacing: '', z: 0, params });

const MODES = {
  circular: (slots, options, askChoice = async () => 'hide') => {
    const state = circularState(slots, options);
    const editor = createCircularTrackSlotEditor({ state, askChoice });
    return {
      state,
      editor,
      slots: () => state.adv.circular_track_slots,
      setRenderer: (slot, renderer) => editor.updateCircularTrackSlotRenderer(slot, renderer),
      issue: (index) => editor.circularTrackSlotIssue(state.adv.circular_track_slots[index], index),
      reset: () => editor.resetCircularTrackSlotsFromSimpleControls(),
      canDuplicate: (slot) => editor.canDuplicateCircularTrackSlot(slot),
      duplicate: (index) => editor.duplicateCircularTrackSlot(index),
      setPlacement: (slot, placement) => editor.updateCircularTrackSlotPlacement(slot, placement),
      setOverlayParam: (slot, key, value) => editor.setCircularOverlayAnnotationParam(slot, key, value),
      // The composition root's switch: it asks before turning off a row other rows reach.
      setEnabled: askingTrackRowSwitch({ askChoice, runUndoable: (_label, fn) => fn(), busy: () => null,
        rows: () => state.adv.circular_track_slots, mode: 'circular' }, editor.setCircularTrackSlotEnabled),
      offAxis: 'outside',
      referenceOptions: (slot, key) => editor.circularAnnotationAnchorOptions(slot, key).map((option) => option.id),
      referenceKnown: (slot, key) => editor.circularAnnotationAnchorIsKnown(slot, key)
    };
  },
  linear: (slots, options, askChoice = async () => 'hide') => {
    const state = linearState(slots, options);
    const editor = createLinearTrackSlotEditor({ state });
    return {
      state,
      editor,
      slots: () => state.adv.linear_track_slots,
      setRenderer: (slot, renderer) => editor.updateLinearTrackSlotRenderer(slot, renderer),
      issue: (index) => editor.linearTrackSlotIssue(state.adv.linear_track_slots[index], index),
      reset: () => editor.resetLinearTrackSlotsFromSimpleControls(),
      canDuplicate: (slot) => editor.canDuplicateLinearTrackSlot(slot),
      duplicate: (index) => editor.duplicateLinearTrackSlot(index),
      setPlacement: (slot, placement) => editor.updateLinearTrackSlotPlacement(slot, placement),
      setOverlayParam: (slot, key, value) => editor.setLinearOverlayAnnotationParam(slot, key, value),
      setEnabled: askingTrackRowSwitch({ askChoice, runUndoable: (_label, fn) => fn(), busy: () => null,
        rows: () => state.adv.linear_track_slots, mode: 'linear' }, editor.setLinearTrackSlotEnabled),
      offAxis: 'above',
      referenceOptions: (slot, key) => editor.linearAnnotationAnchorOptions(slot, key).map((option) => option.id),
      referenceKnown: (slot, key) => editor.linearAnnotationAnchorIsKnown(slot, key)
    };
  }
};

const ROWS = {
  circular: () => [
    circularRow('features', 'features', { lane_direction: 'outside' }, 'outside'),
    circularRow('ticks', 'ticks', { tick_label_layout: 'label_out_tick_in' }),
    circularRow('gc_skew', 'dinucleotide_skew', { nt: 'GC', positive_color: '#ff0000' })
  ],
  linear: () => [
    linearRow('features', 'features', {}, 'overlay'),
    linearRow('gc_skew', 'dinucleotide_skew', { nt: 'GC', positive_color: '#ff0000' })
  ]
};

for (const [mode, create] of Object.entries(MODES)) {
  // TK-04: the params of the previous renderer made the row invalid, also
  // after switching back.
  test(`${mode}: a renderer change keeps only the params the new renderer accepts`, () => {
    const editor = create(ROWS[mode](), { depthFiles: [depthFile] });
    const index = editor.slots().findIndex((slot) => slot.id === 'gc_skew');
    const slot = editor.slots()[index];

    editor.setRenderer(slot, 'depth');
    // Circular leaves an unset index (series 0); Linear writes 0.
    assert.equal(editor.slots()[index].params.track_index ?? 0, 0);
    assert.deepEqual(Object.keys(editor.slots()[index].params).filter((key) => key !== 'track_index'), []);
    assert.equal(editor.issue(index), '');

    // A Depth row's legend label is its series title (syncDepthSlotLabels).
    editor.slots()[index].params.legend_label = 'hmmt_depth_a';
    editor.setRenderer(editor.slots()[index], 'dinucleotide_content');
    assert.deepEqual({ ...editor.slots()[index].params }, { nt: 'GC' });
    assert.equal(editor.issue(index), '');

    editor.slots()[index].params.nt = 'AT';
    editor.setRenderer(editor.slots()[index], 'dinucleotide_skew');
    assert.deepEqual({ ...editor.slots()[index].params }, { nt: 'AT' }, 'nt belongs to both dinucleotide renderers');
    assert.equal(editor.issue(index), '');
  });

  // OV-1186 (LP-D70 b): a comparison-ring row's label, color and fileName
  // describe its ring, as a legend_label names its track; a renderer change
  // drops them with the ring's other keys.
  test(`${mode}: a renderer change drops the ring keys of a comparison row`, () => {
    const ring = {
      managed: 'circular_conservation', series_key: 'k1', source_index: '0', track_index: '1',
      label: 'Query A', color: '#ff0000', fileName: 'a.fa', legend_label: 'Ring A'
    };
    const renderers = mode === 'circular'
      ? ['spacer', 'features', 'depth', 'dinucleotide_content', 'dinucleotide_skew', 'ticks', 'annotations']
      : ['spacer', 'features', 'depth', 'dinucleotide_content', 'dinucleotide_skew', 'annotations'];
    for (const renderer of renderers) {
      const kept = paramsKeptOnRendererChange(ring, renderer);
      assert.deepEqual(['label', 'color', 'fileName', 'legend_label'].filter((key) => Object.hasOwn(kept, key)), [], renderer);
    }
  });

  // PR 2 of the underlay span (LP-D17, LP-D65): Extend to offers the slots the
  // anchor select offers except the anchor itself (SAME-REJECT, LP-D21), and
  // to_slot leaves the row with Overlay and with the annotations renderer.
  test(`${mode}: Extend to offers the anchor's slots without the anchor and is cleared off Overlay`, () => {
    const overlayRow = { ...(mode === 'circular' ? circularRow : linearRow)('regions', 'annotations',
      { set_id: 'review', anchor_slot: 'features', to_slot: 'gc_skew', layer: 'underlay' }, 'overlay'), z: -1 };
    const editor = create([...ROWS[mode](), overlayRow]);
    const regions = () => editor.slots().find((slot) => slot.id === 'regions');
    assert.deepEqual(editor.referenceOptions(regions()), ['features', 'gc_skew']);
    assert.deepEqual(editor.referenceOptions(regions(), 'to_slot'), ['gc_skew']);
    assert.equal(editor.referenceKnown(regions(), 'to_slot'), true);
    regions().params.to_slot = 'features';
    assert.equal(editor.referenceKnown(regions(), 'to_slot'), false, 'the anchor is not a valid Extend to slot');
    regions().params.to_slot = 'gc_skew';
    editor.setPlacement(regions(), editor.offAxis);
    assert.equal(Object.hasOwn(regions().params, 'to_slot'), false);
    assert.equal(Object.hasOwn(regions().params, 'anchor_slot'), false);
    editor.setPlacement(regions(), 'overlay');
    regions().params.to_slot = 'gc_skew';
    editor.setRenderer(regions(), 'spacer');
    assert.equal(Object.hasOwn(regions().params, 'to_slot'), false);
  });
  // OV-1297 (LP-D91): the GUI has no z control, so Placement Overlay, Layer,
  // the anchor and Extend to give an Overlay annotation row the z its layer
  // needs: one below the lowest (Underlay) or one above the highest
  // (Foreground) z of the slots it names; a z that already meets it is kept.
  test(`${mode}: Overlay annotation edits give the row the z its layer needs`, () => {
    const make = mode === 'circular' ? circularRow : linearRow;
    const rows = ROWS[mode]();
    rows.find((slot) => slot.id === 'gc_skew').z = -2;
    const editor = create([...rows, make('regions', 'annotations', { set_id: 'review', layer: 'foreground' }, mode === 'circular' ? 'outside' : 'above')]);
    editor.state.annotationSets.push({ id: 'review', name: 'review', annotations: [] });
    const index = () => editor.slots().findIndex((slot) => slot.id === 'regions');
    const regions = () => editor.slots()[index()];

    editor.setPlacement(regions(), 'overlay');
    editor.setOverlayParam(regions(), 'anchor_slot', 'features');
    assert.deepEqual([regions().params.layer, regions().z], ['foreground', 1]);
    assert.equal(editor.issue(index()), '', 'Overlay / Foreground draws above its anchor');

    editor.setOverlayParam(regions(), 'layer', 'underlay');
    assert.equal(regions().z, -1);
    assert.equal(editor.issue(index()), '', 'Overlay / Underlay draws below its anchor');

    editor.setOverlayParam(regions(), 'to_slot', 'gc_skew');
    assert.equal(regions().z, -3, 'below the lowest of the anchor and the Extend to slot');
    assert.equal(editor.issue(index()), '');

    editor.setOverlayParam(regions(), 'layer', 'foreground');
    assert.equal(regions().z, 1, 'above the highest of the anchor and the Extend to slot');
    assert.equal(editor.issue(index()), '');

    regions().z = 5;
    editor.setOverlayParam(regions(), 'to_slot', '');
    assert.equal(Object.hasOwn(regions().params, 'to_slot'), false, 'None leaves to_slot unset');
    assert.equal(regions().z, 5, 'a z that meets the layer rule is kept');
  });

  // OV-1297: Placement Overlay alone binds the only anchor offered and gives
  // the row its z, so the row is valid before any other edit.
  test(`${mode}: Placement Overlay with one offered anchor leaves a valid row`, () => {
    const make = mode === 'circular' ? circularRow : linearRow;
    const rows = ROWS[mode]().filter((slot) => slot.id !== 'gc_skew');
    const editor = create([...rows, make('regions', 'annotations', { set_id: 'review', layer: 'foreground' }, mode === 'circular' ? 'outside' : 'above')]);
    editor.state.annotationSets.push({ id: 'review', name: 'review', annotations: [] });
    const index = () => editor.slots().findIndex((slot) => slot.id === 'regions');
    const regions = () => editor.slots()[index()];

    editor.setPlacement(regions(), 'overlay');
    assert.deepEqual([regions().params.anchor_slot, regions().z], ['features', 1]);
    assert.equal(editor.issue(index()), '');
  });

  // OV-1297 review F1: a disabled anchor still counts, so a Layer change made
  // while it is off leaves a row that is valid once it is on again.
  test(`${mode}: a Layer change while the anchor is disabled keeps the row valid when it is enabled again`, () => {
    const make = mode === 'circular' ? circularRow : linearRow;
    const editor = create([...ROWS[mode](), make('regions', 'annotations', { set_id: 'review', layer: 'foreground' }, mode === 'circular' ? 'outside' : 'above')]);
    editor.state.annotationSets.push({ id: 'review', name: 'review', annotations: [] });
    const index = () => editor.slots().findIndex((slot) => slot.id === 'regions');
    const regions = () => editor.slots()[index()];
    const gcSkew = () => editor.slots().find((slot) => slot.id === 'gc_skew');

    editor.setPlacement(regions(), 'overlay');
    editor.setOverlayParam(regions(), 'anchor_slot', 'gc_skew');
    assert.equal(regions().z, 1);
    editor.setEnabled(gcSkew(), false);
    editor.setOverlayParam(regions(), 'layer', 'underlay');
    assert.equal(regions().z, -1);
    editor.setEnabled(gcSkew(), true);
    assert.equal(gcSkew().enabled, true);
    assert.equal(editor.issue(index()), '');
  });

  // GX-01: Duplicate uses the availability term of the other row controls.
  test(`${mode}: Duplicate is unavailable while a Session operation runs`, () => {
    const editor = create(ROWS[mode]());
    const slot = editor.slots().find((candidate) => candidate.id === 'gc_skew');
    assert.equal(editor.canDuplicate(slot), true);
    editor.state.sessionOperationAvailability = () => ({ code: 'SESSION_BUSY' });
    assert.equal(editor.canDuplicate(slot), false);
  });

  // OV-1099 (Owner decision 2026-10-18): both modes name a copy from the source
  // row's own ID, without one trailing number, and the first unused number from
  // 2, so a copy of a copy is numbered on, not suffixed again.
  test(`${mode}: Duplicate names the copy from the source row's ID with the next unused number`, () => {
    const editor = create(ROWS[mode]());
    const ids = () => editor.slots().map((slot) => slot.id);
    editor.duplicate(ids().indexOf('gc_skew'));
    assert.deepEqual(ids().slice(-2), ['gc_skew', 'gc_skew_2']);
    editor.duplicate(ids().indexOf('gc_skew_2'));
    assert.deepEqual(ids().slice(-3), ['gc_skew', 'gc_skew_2', 'gc_skew_3']);
    const copy = editor.slots().at(-1);
    assert.equal(copy.params.positive_color, '#ff0000', 'the copy keeps the source params');
  });

  test(`${mode}: Duplicate of a row with its own ID keeps that ID as the base`, () => {
    const rows = ROWS[mode]();
    rows.at(-1).id = 'my_skew';
    const editor = create(rows);
    const ids = () => editor.slots().map((slot) => slot.id);
    editor.duplicate(ids().indexOf('my_skew'));
    assert.deepEqual(ids().slice(-2), ['my_skew', 'my_skew_2']);
    editor.duplicate(ids().indexOf('my_skew'));
    assert.deepEqual(ids().slice(-3), ['my_skew', 'my_skew_3', 'my_skew_2'], 'a second copy of the source takes the next free number');
    editor.duplicate(ids().indexOf('my_skew_2'));
    assert.deepEqual(ids().slice(-4), ['my_skew', 'my_skew_3', 'my_skew_2', 'my_skew_4'], 'a copy of a copy is numbered on');
  });

  test(`${mode}: a copy never takes the bare base ID, even when it is free`, () => {
    const rows = ROWS[mode]();
    rows.at(-1).id = 'my_skew_2';
    const editor = create(rows);
    const ids = () => editor.slots().map((slot) => slot.id);
    editor.duplicate(ids().indexOf('my_skew_2'));
    assert.deepEqual(ids().slice(-2), ['my_skew_2', 'my_skew_3']);
  });

  // TK-10: Reset rebuilds from the loaded Depth series; Show Depth is not an
  // input (docs/REFERENCE/web-app.md, Custom Track Slots).
  test(`${mode}: Reset keeps one Depth row per loaded series while Show Depth is off`, () => {
    for (const showDepth of [false, true]) {
      const editor = create(ROWS[mode](), { showDepth, depthFiles: [depthFile] });
      editor.reset();
      const depthRows = editor.slots().filter((slot) => slot.renderer === 'depth');
      assert.deepEqual(depthRows.map((slot) => slot.params.track_index ?? 0), [0], `Show Depth ${showDepth}`);
    }
    const withoutDepth = create(ROWS[mode](), { showDepth: false, depthFiles: [] });
    withoutDepth.reset();
    assert.equal(withoutDepth.slots().some((slot) => slot.renderer === 'depth'), false);
  });
}

test('circular: Reset to preset keeps the Depth row while Show Depth is off', () => {
  const editor = MODES.circular(ROWS.circular(), { showDepth: false, depthFiles: [depthFile] });
  editor.editor.resetCircularTrackSlotsToPreset('middle');
  assert.equal(editor.slots().filter((slot) => slot.renderer === 'depth').length, 1);
});

test('circular: a renderer change back to Ticks leaves a valid Ticks row', () => {
  const editor = MODES.circular(ROWS.circular());
  const index = editor.slots().findIndex((slot) => slot.id === 'ticks');
  editor.setRenderer(editor.slots()[index], 'dinucleotide_content');
  assert.deepEqual({ ...editor.slots()[index].params }, { nt: 'GC' });
  assert.equal(editor.issue(index), '');
  editor.setRenderer(editor.slots()[index], 'ticks');
  assert.equal(editor.slots()[index].params.nt, undefined);
  assert.equal(editor.issue(index), '');
});

// TK-09: Hide GC Skew draws no GC skew in the simple path; in the custom stack
// it disables only the rows of the drawing's dinucleotide.
test('circular: Hide GC Skew disables only the GC skew rows and names them in the right number', async () => {
  // The app's choice dialog asks; the stub records the question and hides.
  const questions = [];
  const editor = MODES.circular([
    ...ROWS.circular(),
    circularRow('gc_skew_2', 'dinucleotide_skew', { nt: 'AT' })
  ], undefined, async (request) => { questions.push(request); return 'hide'; });
  await editor.editor.setCircularSkewSuppressed(true);
  const rows = Object.fromEntries(editor.slots().map((slot) => [slot.id, slot.enabled]));
  assert.equal(rows.gc_skew, false);
  assert.equal(rows.gc_skew_2, true, 'the AT skew row stays enabled');
  assert.equal(editor.editor.circularTrackSlotSuppressMessage(editor.slots().find((slot) => slot.id === 'gc_skew_2')), '');
  assert.equal(questions.length, 1);
  assert.equal(questions[0].title, 'Hide GC skew?');
  assert.match(questions[0].description, /has 1 GC skew track\. Hiding GC skew turns that track off\./);
  assert.deepEqual(questions[0].actions.map((action) => action.label), ['Hide GC skew']);

  await editor.editor.setCircularSkewSuppressed(false);
  assert.equal(editor.slots().every((slot) => slot.enabled), true);

  editor.editor.duplicateCircularTrackSlot(editor.slots().findIndex((slot) => slot.id === 'gc_skew'));
  await editor.editor.setCircularSkewSuppressed(true);
  assert.match(questions[1].description, /has 2 GC skew tracks\. Hiding GC skew turns those tracks off\./);

  // Keep leaves the setting and the rows unchanged.
  const kept = MODES.circular(ROWS.circular(), undefined, async () => 'cancel');
  assert.equal(await kept.editor.setCircularSkewSuppressed(true), false);
  assert.equal(kept.state.form.suppress_skew, false);
  assert.equal(kept.slots().every((slot) => slot.enabled), true);
});

// The question is asked before any History step (app-setup.js gives the editor
// History's `runUndoable` and the app's one choice prompt): the chosen hide is
// one step recorded at the answer, and Keep records none. A hide chosen while
// another action is still running applies at once as its own step, and what
// that action writes afterwards is a step under its name (OV-460).
test('circular: a hide that asks is one History step recorded at the answer, also while another action runs', async () => {
  const { createChoicePrompt } = await import(pathToFileURL(join(tempRoot, 'app', 'app-messages.js')));
  const { createHistoryManager } = await import('../../gbdraw/web/js/services/history.js');
  const state = Object.assign(circularState(ROWS.circular()), { loaded: 0 });
  const { ask, resolve } = createChoicePrompt({ reactive: (value) => value });
  const history = createHistoryManager({
    buildIntent: () => structuredClone({ hide: state.form.suppress_skew, slots: state.adv.circular_track_slots, loaded: state.loaded }),
    applyIntent: (intent) => {
      state.form.suppress_skew = intent.hide;
      state.adv.circular_track_slots.splice(0, Infinity, ...structuredClone(intent.slots));
      state.loaded = intent.loaded;
    },
    buildCheckpoint: () => ({}), applyCheckpoint: () => {}
  });
  const editor = createCircularTrackSlotEditor({ state, askChoice: ask, runUndoable: history.runUndoable });
  const read = () => [state.form.suppress_skew, state.adv.circular_track_slots.find((slot) => slot.id === 'gc_skew').enabled, state.loaded];
  const steps = () => [history.getUndoCount(), history.undoLabel()];
  const tick = () => new Promise((settle) => setImmediate(settle));
  await history.initializeIntentBaseline();
  /** @param {string} choice */
  const hide = async (choice) => {
    const asked = editor.setCircularSkewSuppressed(true);
    await tick();
    assert.deepEqual(read().slice(0, 2), [false, true], 'nothing is written before the answer');
    resolve(choice);
    return asked;
  };

  const asking = editor.setCircularSkewSuppressed(true);
  await tick();
  assert.equal(history.mutationPending(), false, 'no step is open while the question is asked');
  resolve('cancel');
  assert.equal(await asking, false);
  assert.deepEqual([read(), history.getUndoCount()], [[false, true, 0], 0]);

  await hide('hide');
  assert.deepEqual([read(), steps()], [[true, false, 0], [1, 'Hide GC Skew']]);
  await history.undo();
  assert.deepEqual(read(), [false, true, 0]);
  await history.redo();
  assert.deepEqual(read(), [true, false, 0]);
  await history.undo();

  let release;
  const held = new Promise((settle) => { release = settle; });
  const load = history.runUndoable('Load label edits', async () => { await held; state.loaded = 1; });
  await tick();
  await hide('hide');
  assert.deepEqual([read(), steps(), history.mutationPending()], [[true, false, 0], [1, 'Hide GC Skew'], false]);
  release();
  await load;
  assert.deepEqual([read(), steps()], [[true, false, 1], [2, 'Load label edits']]);
  await history.undo();
  await history.undo();
  assert.deepEqual(read(), [false, true, 0]);
  await history.redo();
  await history.redo();
  assert.deepEqual(read(), [true, false, 1]);
});

// TK-07: `1.5` or `-1` became track index 1 without a message.
test('circular: an invalid Depth track index stays visible with a field error and the row keeps its index', async () => {
  const Vue = loadPageVue();
  globalThis.window = { Vue };
  const { reactive, nextTick } = Vue;
  const { DepthTrackIndexInput } = await import(pathToFileURL(join(tempRoot, 'app', 'circular-track-slots', 'track-index-input.js')));
  const props = reactive({ modelValue: 1, slotId: 'depth_2', controlId: 'slot-track-index', helpId: 'help-index', disabled: false });
  const emissions = [];
  const field = DepthTrackIndexInput.setup(props, { emit: (event, value) => {
    emissions.push(value);
    props.modelValue = value;
  } });

  for (const entry of ['1.5', '-1', '']) {
    field.commit(entry);
    assert.equal(field.text.value, entry);
    assert.equal(field.error.value, 'Enter a whole number of 0 or more. The row keeps track index 1.');
    assert.equal(field.describedBy.value, 'help-index slot-track-index-error');
  }
  assert.deepEqual(emissions, [], 'the row keeps its last valid index');

  field.commit('2');
  assert.deepEqual(emissions, [2]);
  assert.equal(field.text.value, '2');
  assert.equal(field.error.value, '');
  assert.equal(field.describedBy.value, 'help-index');

  field.commit('-1');
  props.modelValue = 0; // Undo restores another index and drops the rejected entry.
  await nextTick();
  assert.equal(field.text.value, '0');
  assert.equal(field.error.value, '');
});

// TK-09 review: the Dinucleotide setting and a row's dinucleotide are inputs of
// the Hide GC projection, so their edits run it; the editor state equals what
// Generate projects (run-analysis applies the same projection to the request).
test('circular: Dinucleotide and row dinucleotide edits re-run the Hide GC Skew projection', async () => {
  const { applyCircularSuppressControlsToSlots, applyCircularTrackOrderPlacements } = await import(
    pathToFileURL(join(tempRoot, 'services', 'circular-track-slot-model.js'))
  );
  {
    const editor = MODES.circular(ROWS.circular());
    const { state } = editor;
    const requestSlots = () => applyCircularSuppressControlsToSlots(
      applyCircularTrackOrderPlacements(state.adv.circular_track_slots, state.adv.nt, state.form.track_type, state.adv.circular_track_slots_axis_index),
      state.form,
      state.adv.nt
    );
    const enabled = () => editor.slots().find((slot) => slot.id === 'gc_skew').enabled;
    await editor.editor.setCircularSkewSuppressed(true);
    assert.equal(enabled(), false);

    editor.editor.setCircularDinucleotide('AT');
    assert.equal(state.adv.nt, 'AT');
    assert.equal(enabled(), true, 'the GC row no longer uses the Dinucleotide setting');
    assert.deepEqual(JSON.parse(JSON.stringify(editor.slots())), JSON.parse(JSON.stringify(requestSlots())));

    const row = editor.slots().find((slot) => slot.id === 'gc_skew');
    editor.editor.setCircularTrackSlotNt(row, ' AT ');
    assert.equal(row.params.nt, 'AT');
    assert.equal(enabled(), false, 'the row now uses the Dinucleotide setting');
    assert.deepEqual(JSON.parse(JSON.stringify(editor.slots())), JSON.parse(JSON.stringify(requestSlots())));
  }
});

// PD-OI-104 scenario revision 6 (D-217, D-219): hiding a drawn track row that
// another row starts at or extends to (an Overlay row's anchor or Extend to
// slot) asks; the question names those rows as a refused Generate does
// ("Track row N", with the slot id) and offers to turn them off with the hide
// as one History step; Cancel changes nothing.
const REACH_ROWS = {
  circular: () => [
    circularRow('features', 'features', { lane_direction: 'outside' }, 'outside'),
    circularRow('gc_content', 'dinucleotide_content', { nt: 'GC' }),
    circularRow('depth', 'depth', { track_index: 0 }),
    { ...circularRow('tint', 'annotations', { set_id: 'tint_set', anchor_slot: 'features', to_slot: 'gc_content', layer: 'underlay' }, 'overlay'), z: -1 },
    { ...circularRow('mark', 'annotations', { set_id: 'tint_set', anchor_slot: 'depth', layer: 'foreground' }, 'overlay'), z: 1 }
  ],
  linear: () => [
    linearRow('features', 'features', {}, 'overlay'),
    linearRow('gc_content', 'dinucleotide_content', { nt: 'GC' }),
    linearRow('depth', 'depth', { track_index: 0 }),
    { ...linearRow('tint', 'annotations', { set_id: 'tint_set', anchor_slot: 'features', to_slot: 'gc_content', layer: 'underlay' }, 'overlay'), z: -1 },
    { ...linearRow('mark', 'annotations', { set_id: 'tint_set', anchor_slot: 'depth', layer: 'foreground' }, 'overlay'), z: 1 }
  ]
};
const enabledOf = (editor) => Object.fromEntries(editor.slots().map((slot) => [slot.id, slot.enabled]));

for (const mode of ['circular', 'linear']) {
  test(`${mode}: turning off a GC or Depth row that another row reaches asks, names that row, and turns it off too`, async () => {
    for (const [rowId, reachId, number] of [['gc_content', 'tint', 4], ['depth', 'mark', 5]]) {
      const questions = [];
      const editor = MODES[mode](REACH_ROWS[mode](), { showDepth: true, depthFiles: [depthFile] }, async (request) => { questions.push(request); return 'hide'; });
      const row = editor.slots().find((slot) => slot.id === rowId);
      const box = { checked: false };
      const asked = editor.setEnabled(row, false, { target: box });
      assert.ok(asked instanceof Promise, `${rowId}: it asks`);
      assert.equal(box.checked, true, `${rowId}: the checkbox shows the row on until the answer`);
      await asked;
      const rowNumber = editor.slots().indexOf(row) + 1;
      assert.equal(questions.length, 1);
      assert.equal(questions[0].title, `Turn off Track row ${rowNumber}?`);
      assert.deepEqual(questions[0].items, [`Track row ${number} (${reachId})`]);
      assert.match(questions[0].description, /This row starts at or extends to a track being turned off\./);
      assert.deepEqual(questions[0].actions.map((action) => action.label), [`Turn off Track row ${rowNumber} and Track row ${number}`]);
      const enabled = enabledOf(editor);
      assert.equal(enabled[rowId], false, `${rowId}: the row is off`);
      assert.equal(enabled[reachId], false, `${rowId}: the row that reaches it is off`);
      assert.equal(Object.entries(enabled).filter(([, on]) => !on).length, 2, `${rowId}: no other row changes`);
    }
    // Cancel changes nothing; a row no other row reaches turns off without a question.
    const kept = MODES[mode](REACH_ROWS[mode](), { showDepth: true, depthFiles: [depthFile] }, async () => 'cancel');
    assert.equal(await kept.setEnabled(kept.slots().find((slot) => slot.id === 'gc_content'), false), false);
    assert.ok(Object.values(enabledOf(kept)).every(Boolean), 'Cancel keeps every row on');
    const quiet = MODES[mode](REACH_ROWS[mode](), { showDepth: true, depthFiles: [depthFile] }, async () => assert.fail('no question'));
    quiet.setEnabled(quiet.slots().find((slot) => slot.id === 'tint'), false);
    assert.equal(enabledOf(quiet).tint, false);
  });

  test(`${mode}: the chosen turning off is one History step whose Undo turns both rows on again`, async () => {
    const { createChoicePrompt } = await import(pathToFileURL(join(tempRoot, 'app', 'app-messages.js')));
    const { createHistoryManager } = await import('../../gbdraw/web/js/services/history.js');
    const state = mode === 'circular' ? circularState(REACH_ROWS.circular()) : linearState(REACH_ROWS.linear());
    const slots = () => (mode === 'circular' ? state.adv.circular_track_slots : state.adv.linear_track_slots);
    const { ask, resolve } = createChoicePrompt({ reactive: (value) => value });
    const history = createHistoryManager({
      buildIntent: () => structuredClone({ slots: slots() }),
      applyIntent: (intent) => { slots().splice(0, Infinity, ...structuredClone(intent.slots)); },
      buildCheckpoint: () => ({}), applyCheckpoint: () => {}
    });
    const editor = mode === 'circular' ? createCircularTrackSlotEditor({ state }) : createLinearTrackSlotEditor({ state });
    const setEnabled = askingTrackRowSwitch({ askChoice: ask, runUndoable: history.runUndoable, busy: () => null, rows: slots, mode },
      mode === 'circular' ? editor.setCircularTrackSlotEnabled : editor.setLinearTrackSlotEnabled);
    const read = () => slots().filter((slot) => ['gc_content', 'tint'].includes(slot.id)).map((slot) => slot.enabled);
    await history.initializeIntentBaseline();
    const asking = setEnabled(slots().find((slot) => slot.id === 'gc_content'), false);
    await new Promise((settle) => setImmediate(settle));
    assert.deepEqual(read(), [true, true], 'nothing is written before the answer');
    resolve('hide');
    await asking;
    assert.deepEqual([read(), history.getUndoCount(), history.undoLabel()], [[false, false], 1, 'Turn off track rows']);
    await history.undo();
    assert.deepEqual(read(), [true, true]);
  });
}

test('circular: Hide GC Content names the rows that reach the GC content track and turns them off with it', async () => {
  const questions = [];
  const editor = MODES.circular(REACH_ROWS.circular(), { showDepth: true, depthFiles: [depthFile] }, async (request) => { questions.push(request); return 'hide'; });
  await editor.editor.setCircularGcSuppressed(true);
  assert.equal(questions[0].title, 'Hide GC content?');
  assert.match(questions[0].description, /has 1 GC content track\. Hiding GC content turns that track off\. This row starts at or extends to a track being turned off\./);
  assert.deepEqual(questions[0].items, ['Track row 4 (tint)']);
  assert.deepEqual(questions[0].actions.map((action) => action.label), ['Hide GC content and turn off Track row 4']);
  assert.equal(editor.state.form.suppress_gc, true);
  const enabled = enabledOf(editor);
  assert.deepEqual([enabled.gc_content, enabled.tint, enabled.mark, enabled.depth], [false, false, true, true]);
  // Keep changes nothing.
  const kept = MODES.circular(REACH_ROWS.circular(), undefined, async () => 'cancel');
  assert.equal(await kept.editor.setCircularGcSuppressed(true), false);
  assert.equal(kept.state.form.suppress_gc, false);
  assert.ok(Object.values(enabledOf(kept)).every(Boolean));
});

// The hide button is one whole sentence per hide (a translation never joins a
// control name into a sentence).
test('circular: the Hide GC Content button with the rows it turns off is one translated sentence', async () => {
  const { setUiLocale } = await import(pathToFileURL(join(tempRoot, 'utils', 'i18n.js')));
  await setUiLocale('ja', async () => ({ default: {
    'Hide GC content and turn off {rows}': '{rows} をオフにして GC content を隠す',
    'Track row {row}': 'トラック行 {row}'
  } }));
  try {
    const questions = [];
    const editor = MODES.circular(REACH_ROWS.circular(), { showDepth: true, depthFiles: [depthFile] }, async (request) => { questions.push(request); return 'cancel'; });
    await editor.editor.setCircularGcSuppressed(true);
    assert.deepEqual(questions[0].actions.map((action) => action.label), ['トラック行 4 をオフにして GC content を隠す']);
  } finally {
    await setUiLocale('en');
  }
});
