import assert from 'node:assert/strict';
import test from 'node:test';
import { readFileSync } from 'node:fs';
import {
  buildLegendStyleRetirement,
  createDrawingLegendRetirement,
  trackDataLegendCaptions
} from '../../gbdraw/web/js/app/legend/track-data-styles.js';
import { PAINT_SLICE } from '../../gbdraw/web/js/services/drawing/paint.js';

// A drawing whose paint is these members (the fixture's own objects).
const paintDrawing = (members) => ({ slices: { paint: PAINT_SLICE.create(members).handle } });
import {
  depthSeriesLegendCaptions, ensureDepthTrackConfigCount, getDepthTrackFileBaseName, isDepthTrackAutoLabel
} from '../../gbdraw/web/js/services/depth-track-state.js';

test('track data names the legend labels of annotation sets and the sourced Depth series (OV-65)', () => {
  const captions = trackDataLegendCaptions({
    annotationSets: [
      { legendLabel: ' Set ', annotations: [{ legendLabel: 'Region X' }, { legendLabel: null }] },
      { legendLabel: '', annotations: [] }
    ],
    depthTracks: [{ label: 'sample' }, { label: 'other' }],
    depthSlots: [
      { renderer: 'depth', params: { track_index: 0, legend_label: 'Coverage' } },
      { renderer: 'depth', params: { track_index: 1, legend_label: 'Unsourced' } },
      { renderer: 'features', params: { legend_label: 'Not Depth' } }
    ],
    sourcedDepthTrackIndexes: [0]
  });
  assert.deepEqual([...captions].sort(), ['Coverage', 'Region X', 'Set', 'sample']);
});

// The Depth rows Python draws for series without a label: Depth for one,
// Depth 1, Depth 2 ... for several (`_default_depth_track_labels`; the shared
// vector's unnamed Circular cases). A series the Web adds keeps no label, so
// its name follows the series count as the CLI draws it (OV-1508).
const UNNAMED_DEPTH_CASES = JSON.parse(readFileSync(new URL('../fixtures/legend_track_slot_rows.json', import.meta.url), 'utf8'))
  .cases.filter(({ mode, slots, depthTracks }) => mode === 'circular' && slots === null
    && Array.isArray(depthTracks) && depthTracks.every(({ label }) => label === ''));

test('a Depth series the Web adds without a name keeps no label and is named as Python names it (OV-1508)', () => {
  assert.deepEqual(UNNAMED_DEPTH_CASES.map(({ depthTracks }) => depthTracks.length), [1, 3]);
  for (const { name, depthTracks, rows } of UNNAMED_DEPTH_CASES) {
    const pythonNames = rows.filter(([, slot]) => slot === 'depth').map(([key]) => key);
    const configs = ensureDepthTrackConfigCount([], depthTracks.length);
    assert.deepEqual(configs.map(({ label }) => label), depthTracks.map(() => ''), name);
    const indexes = depthTracks.map((_, index) => index);
    assert.deepEqual([...trackDataLegendCaptions({ depthTracks: configs, sourcedDepthTrackIndexes: indexes })], pythonNames, name);
    assert.deepEqual([...trackDataLegendCaptions({ sourcedDepthTrackIndexes: indexes })], pythonNames, name);
    assert.deepEqual(depthSeriesLegendCaptions({ depthTracks: configs, slots: [] }), pythonNames, name);
  }
});

test('a default name of its series, an empty label, or its file name is a label the next file may replace', () => {
  assert.equal(getDepthTrackFileBaseName({ name: 'dir/reads.depth.tsv' }), 'reads.depth');
  assert.equal(getDepthTrackFileBaseName(null), '');
  // A Session saved before OV-1508 holds Depth for the first of several series.
  for (const [label, index] of [['', 0], ['Depth', 0], ['Depth 1', 0], ['Depth 2', 1], ['Depth 3', 2]]) {
    assert.equal(isDepthTrackAutoLabel(label, index), true, `${label} at ${index}`);
  }
  for (const [label, index] of [['Depth', 1], ['Depth 3', 1], ['Reads', 0]]) {
    assert.equal(isDepthTrackAutoLabel(label, index), false, `${label} at ${index}`);
  }
  assert.equal(isDepthTrackAutoLabel('reads.depth', 0, { name: 'reads.depth.tsv' }), true);
});

test('styles follow the caption: a data change retires only the captions the data stops naming', () => {
  const sets = [{ legendLabel: '', annotations: [{ legendLabel: 'Keep' }, { legendLabel: 'Drop' }] }];
  const legendColorOverrides = { Keep: '#111111', Drop: '#222222', Feature: '#333333' };
  const legendStrokeOverrides = { Drop: { strokeColor: '#000000' }, Feature: { strokeWidth: 2 } };
  const retire = buildLegendStyleRetirement({
    drawing: paintDrawing({ legendColorOverrides, legendStrokeOverrides, legendEntries: { value: [] } }),
    projectLegendEntries: () => assert.fail('no rename to retire'),
    namedCaptions: () => trackDataLegendCaptions({ annotationSets: sets })
  });
  const result = retire(() => {
    sets[0].annotations.pop();
    return 'changed';
  });
  assert.equal(result, 'changed');
  assert.deepEqual(legendColorOverrides, { Keep: '#111111', Feature: '#333333' });
  assert.deepEqual(legendStrokeOverrides, { Feature: { strokeWidth: 2 } });
});

test('a change that keeps every caption retires nothing', () => {
  const legendColorOverrides = { Keep: '#111111' };
  const legendEntries = { value: [{ caption: 'Kept name', originalCaption: 'Keep' }] };
  const retire = buildLegendStyleRetirement({
    drawing: paintDrawing({ legendColorOverrides, legendStrokeOverrides: {}, legendEntries }),
    projectLegendEntries: () => assert.fail('no rename to retire'),
    namedCaptions: () => new Set(['Keep'])
  });
  retire(() => undefined);
  assert.deepEqual(legendColorOverrides, { Keep: '#111111' });
  assert.deepEqual(legendEntries.value, [{ caption: 'Kept name', originalCaption: 'Keep' }]);
});

// OV-87: names follow the caption as styles do. The rename of a row whose caption
// the data stops naming is retired with the styles stored under its key (D-25).
test('a data change retires the Legend rename of a caption it stops naming, and the styles of its row', () => {
  const named = new Set(['depth', 'Region X']);
  const legendColorOverrides = { depth: '#7b2cbf', 'Region X': '#111111', CDS: '#222222' };
  const legendStrokeOverrides = { depth: { strokeWidth: 2 } };
  const entries = [
    { caption: 'CDS', originalCaption: 'CDS', color: '#222222' },
    { caption: 'Coverage', originalCaption: 'depth', color: '#7b2cbf' },
    { caption: 'Region Y', originalCaption: 'Region X', color: '#111111' }
  ];
  const legendEntries = { value: entries };
  let projections = 0;
  const retire = buildLegendStyleRetirement({
    drawing: paintDrawing({ legendColorOverrides, legendStrokeOverrides, legendEntries }),
    projectLegendEntries: () => { projections += 1; },
    namedCaptions: () => new Set(named)
  });
  retire(() => named.delete('depth'));
  assert.deepEqual(legendEntries.value.map(({ originalCaption, caption }) => `${originalCaption}=>${caption}`),
    ['CDS=>CDS', 'depth=>depth', 'Region X=>Region Y']);
  assert.equal(legendEntries.value[1].color, '#7b2cbf', 'the entry keeps what the displayed row shows');
  assert.deepEqual(legendColorOverrides, { 'Region X': '#111111', CDS: '#222222' });
  assert.deepEqual(legendStrokeOverrides, {});
  assert.equal(projections, 1, 'the displayed Result shows the retired rename');
  assert.equal(entries[1].caption, 'Coverage', 'the History step keeps the entries it captured');
});

test('a new name that the new data names keeps its styles when the rename is retired', () => {
  const named = new Set(['depth']);
  const legendColorOverrides = { Coverage: '#7b2cbf' };
  const legendEntries = { value: [{ caption: 'Coverage', originalCaption: 'depth' }] };
  const retire = buildLegendStyleRetirement({
    drawing: paintDrawing({ legendColorOverrides, legendStrokeOverrides: {}, legendEntries }),
    projectLegendEntries: () => {},
    namedCaptions: () => new Set(named)
  });
  retire(() => { named.delete('depth'); named.add('Coverage'); });
  assert.deepEqual(legendEntries.value, [{ caption: 'depth', originalCaption: 'depth' }]);
  assert.deepEqual(legendColorOverrides, { Coverage: '#7b2cbf' });
});

// OV-120: a renamed row a Generate hid (Show Depth off) waits in the drawing;
// removing its data retires the waiting rename and its styles like a shown one.
test('a data change retires a waiting Legend rename of a caption it stops naming (OV-120)', () => {
  const named = new Set(['depth', 'Region X']);
  const legendColorOverrides = { depth: '#7b2cbf', 'Region X': '#111111' };
  const legendStrokeOverrides = { depth: { strokeWidth: 2 } };
  const dormantLegendEntries = {
    value: [
      { caption: 'Coverage', originalCaption: 'depth', color: '#7b2cbf' },
      { caption: 'Region Y', originalCaption: 'Region X', color: '#111111' }
    ]
  };
  const retire = buildLegendStyleRetirement({
    drawing: paintDrawing({ legendColorOverrides, legendStrokeOverrides, legendEntries: { value: [] }, dormantLegendEntries }),
    projectLegendEntries: () => assert.fail('no shown rename to retire'),
    namedCaptions: () => new Set(named)
  });
  retire(() => named.delete('depth'));
  assert.deepEqual(dormantLegendEntries.value, [{ caption: 'Region Y', originalCaption: 'Region X', color: '#111111' }]);
  assert.deepEqual(legendColorOverrides, { 'Region X': '#111111' });
  assert.deepEqual(legendStrokeOverrides, {});
});

// OV-1198, OV-1215: the rows of a drawn track the draft hides leave the list
// with the track and return with it. The port asks the Legend list owner to
// list them (here a stand-in that, like `listLegendRows`, lists the Result's
// rows except the hidden ones, named by the list's rename, else a dormant one).
const trackPort = ({ entries, dormant = [], hidden, result = ['CDS', 'depth', 'GC content'] }) => {
  const legendEntries = { value: entries };
  const dormantLegendEntries = { value: dormant };
  const state = { hidden };
  let projections = 0;
  const keyOf = (entry) => entry.originalCaption || entry.caption;
  const renameOf = (key, list) => list.find((entry) => keyOf(entry) === key && entry.caption !== key)?.caption;
  const run = buildLegendStyleRetirement({
    drawing: paintDrawing({ legendColorOverrides: {}, legendStrokeOverrides: {}, legendEntries, dormantLegendEntries }),
    namedCaptions: () => new Set(),
    hiddenRows: () => (state.hidden ? ['depth'] : []),
    relistLegendEntries: () => {
      const listed = new Set(legendEntries.value.map(keyOf));
      legendEntries.value = result.filter((key) => !(state.hidden && key === 'depth')).map((key) => ({
        caption: (listed.has(key) ? renameOf(key, legendEntries.value) : renameOf(key, dormantLegendEntries.value)) || key,
        originalCaption: key
      }));
    },
    projectLegendEntries: () => { projections += 1; }
  });
  const captions = () => legendEntries.value.map((entry) => entry.caption);
  return { legendEntries, dormantLegendEntries, state, run, captions, projections: () => projections };
};

test('a hidden track\'s row leaves the list, and a renamed one waits as after a Generate (OV-120)', () => {
  const port = trackPort({
    entries: [
      { caption: 'CDS', originalCaption: 'CDS' },
      { caption: 'Coverage', originalCaption: 'depth' },
      { caption: 'GC content', originalCaption: 'GC content' }
    ],
    hidden: false
  });
  port.run(() => { port.state.hidden = true; });
  assert.deepEqual(port.captions(), ['CDS', 'GC content']);
  assert.deepEqual(port.dormantLegendEntries.value, [{ caption: 'Coverage', originalCaption: 'depth' }]);
  assert.equal(port.projections(), 1);
  // Nothing is left to move on a later port run that hides nothing more.
  port.run(() => undefined);
  assert.equal(port.projections(), 1);
});

test('a track shown again returns its renamed row and leaves no dormant copy (review H1)', () => {
  const port = trackPort({
    entries: [
      { caption: 'CDS', originalCaption: 'CDS' },
      { caption: 'Coverage', originalCaption: 'depth' },
      { caption: 'GC content', originalCaption: 'GC content' }
    ],
    hidden: false
  });
  port.run(() => { port.state.hidden = true; });
  port.run(() => { port.state.hidden = false; });
  assert.deepEqual(port.captions(), ['CDS', 'Coverage', 'GC content']);
  assert.deepEqual(port.dormantLegendEntries.value, [], 'a Session saved now holds the row once');
  assert.equal(port.projections(), 2);
  // The user's later rename back to the key is not undone by a stale dormant entry.
  port.legendEntries.value = port.legendEntries.value.map((entry) => ({ ...entry, caption: entry.originalCaption }));
  port.run(() => { port.state.hidden = true; });
  port.run(() => { port.state.hidden = false; });
  assert.deepEqual(port.captions(), ['CDS', 'depth', 'GC content']);
});

test('a change that hides no drawn track leaves the list as it was and projects nothing', () => {
  const port = trackPort({ entries: [{ caption: 'depth', originalCaption: 'depth' }], hidden: false });
  port.run(() => undefined);
  assert.deepEqual(port.legendEntries.value, [{ caption: 'depth', originalCaption: 'depth' }]);
  assert.equal(port.projections(), 0);
});

// Review S2: drawings C (shown) and "C copy" share the Depth file `a.tsv`, and
// "C copy" has a Legend color on the auto-named series "a". Replacing the file
// from C with "Every drawing" renames the series in both; each drawing's port
// measures its own captions, so "C copy"'s paint under "a" is retired like C's.
// Only the shown drawing relists and projects (they read the displayed Result).
test('the Legend port retires the captions of the drawing it is given, also when that drawing is not shown (review S2)', () => {
  const depthDrawing = (id, color, entries = []) => {
    const members = { legendColorOverrides: { a: color }, legendStrokeOverrides: {}, legendEntries: { value: entries } };
    return { id, label: 'a', members, slices: { paint: PAINT_SLICE.create(members).handle } };
  };
  const shown = depthDrawing('circular', '#111111');
  const copy = depthDrawing('circular-2', '#222222', [{ caption: 'Coverage', originalCaption: 'a' }]);
  const calls = [];
  const retire = createDrawingLegendRetirement({
    shownDrawingId: () => shown.id,
    namedCaptions: (drawing) => new Set([drawing.label]),
    hiddenRows: () => [],
    relistLegendEntries: () => calls.push('relist'),
    projectLegendEntries: () => calls.push('project')
  });
  retire(copy, () => { copy.label = 'b'; });
  assert.deepEqual(copy.members.legendColorOverrides, {});
  assert.deepEqual(copy.members.legendEntries.value, [{ caption: 'a', originalCaption: 'a' }]);
  assert.deepEqual(shown.members.legendColorOverrides, { a: '#111111' }, 'the shown drawing keeps its paint');
  assert.deepEqual(calls, [], 'a drawing that is not shown neither relists nor projects the displayed Result');

  shown.members.legendEntries.value = [{ caption: 'Coverage', originalCaption: 'a' }];
  retire(shown, () => { shown.label = 'b'; });
  assert.deepEqual(shown.members.legendColorOverrides, {});
  assert.deepEqual(calls, ['project'], 'the shown drawing projects its retired rename');
});
