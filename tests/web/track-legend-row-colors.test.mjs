// OV-484: a track's Legend row shows the color of the track slot that draws
// it. The Legend panel, the live compile and the palette's track painter read
// it from the tracks of the request that drew the Result
// (services/track-legend-rows.js), never from the shape of the row's caption.
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { test } from 'node:test';

import { compileDirectEditorMutationPlan, draftLegendPanelColors, draftLegendRowColors } from '../../gbdraw/web/js/app/candidate-render.js';
import { buildLinearTrackSlotPayload, normalizeLinearTrackSlots } from '../../gbdraw/web/js/services/linear-track-slot-model.js';
import {
  draftDrawsDepth, drawnTrackHidden, drawnTracks, hiddenTrackLegendKeys, loadedTracks, skewTrackColor, trackLegendRows, trackRowColor
} from '../../gbdraw/web/js/services/track-legend-rows.js';
import { fixture as paintFixture } from './helpers/svg-style-fixture.mjs';
import { modeDrawing } from './helpers/drawing-state.mjs';

const readJson = (path) => JSON.parse(readFileSync(new URL(path, import.meta.url), 'utf8'));
const HMMT_AT_SKEW = readJson('../../gbdraw/web/gallery/sessions/HmmtDNA_ATskew.gbdraw-session.json');
// Python's rows of a Result by their keys, where a case reads no drawn fill.
const drawnRows = (keys) => new Map(keys.map((key) => [key, { key, color: null }]));
const requestOf = ({ mode, dinucleotide, slots, depthTracks }) => ({
  mode,
  diagramOptions: {
    dinucleotide,
    tracks: mode === 'linear' ? { linearTrackSlots: slots } : { circularTrackSlots: slots },
    ...(depthTracks ? { depthTracks } : {})
  }
});

// The slots of the HmmtDNA_ATskew Gallery Session (its AT skew slot has its
// own colors) and a second GC skew slot that repeats the label, with the rows
// Python draws for them.
const SLOTS = [
  'gc_content:dinucleotide_content', 'gc_skew:dinucleotide_skew',
  'a_skew_2:dinucleotide_skew@nt=AT,positive_color=#deaf6e,negative_color=#7294e3,legend_label=AT skew',
  'gc_skew_3:dinucleotide_skew@negative_color=#ff00ff'
];
const TRACKS = drawnTracks(requestOf({ mode: 'circular', dinucleotide: 'GC', slots: SLOTS }));
const DRAWN = [
  ['CDS', '#84b9ec'], ['GC content', '#8daf38'], ['GC skew (+)', '#f6a8f1'], ['GC skew (-)', '#bd96df'],
  ['AT skew (+)', '#deaf6e'], ['AT skew (-)', '#7294e3'], ['GC skew (+) (2)', '#f6a8f1'], ['GC skew (-) (2)', '#ff00ff']
];
const PALETTE = { CDS: '#84b9ec', gc_content: '#8daf38', skew_high: '#f6a8f1', skew_low: '#bd96df' };
const legendState = () => ({
  legendEntries: DRAWN.map(([caption, color]) => ({ caption, originalCaption: caption, color })),
  originalLegendOrder: DRAWN.map(([key]) => key),
  pythonRows: new Map(DRAWN.map(([key, color]) => [key, { key, color }]))
});
const trackColors = (pairs) => pairs.filter(([key]) => key !== 'CDS');

test('the Legend panel shows a skew slot\'s own colors on its rows (OV-484)', () => {
  const { legendEntries, originalLegendOrder, pythonRows } = legendState();
  const panel = draftLegendPanelColors({
    legendEntries, originalLegendOrder, legendColorOverrides: {}, rules: [], pythonRows, paletteColors: PALETTE, drawnTracks: TRACKS
  });
  assert.deepEqual(trackColors([...panel.listed]), trackColors(DRAWN), 'each track row shows the color Python drew it in');
});

const liveLegendFills = (paletteColors, tracks = TRACKS) => {
  const { legendEntries, originalLegendOrder, pythonRows } = legendState();
  const operations = compileDirectEditorMutationPlan({
    catalogAdmission: {
      resultNames: ['diagram.svg'], renderedTargetsByOverrideKey: new Map(), resultIndexesByRenderedId: new Map(),
      renderedFeaturesByResult: [new Map()]
    },
    legendEntries, originalLegendOrder,
    livePreview: { domains: ['legendFills'], paletteColors, drawnContext: null, pythonRows, features: [], drawnTracks: tracks }
  }).operationsByResult[0];
  return operations.legendFills.map(({ caption, color }) => [caption, color]);
};

test('a live palette edit leaves the rows of a skew slot with its own colors, and recolors every other skew row (OV-484)', () => {
  assert.deepEqual(trackColors(liveLegendFills({ ...PALETTE, skew_high: '#112233', skew_low: '#445566' })), [
    ['GC content', '#8daf38'], ['GC skew (+)', '#112233'], ['GC skew (-)', '#445566'],
    ['AT skew (+)', '#deaf6e'], ['AT skew (-)', '#7294e3'],
    // The second GC skew slot: the palette's positive color, its own negative color.
    ['GC skew (+) (2)', '#112233'], ['GC skew (-) (2)', '#ff00ff']
  ]);
});

test('a row the request draws no track for takes no palette color from its caption', () => {
  assert.deepEqual(trackColors(liveLegendFills({ ...PALETTE, skew_high: '#112233' }, [])), []);
});

test('the tracks of the HmmtDNA_ATskew Gallery Session name the slot and the own colors of each track', () => {
  assert.deepEqual(drawnTracks(HMMT_AT_SKEW.renderRequest), [
    { slotId: 'gc_content', kind: 'content', label: 'GC content', positiveColor: null, negativeColor: null },
    { slotId: 'gc_skew', kind: 'skew', label: 'GC skew', positiveColor: null, negativeColor: null },
    { slotId: 'a_skew_2', kind: 'skew', label: 'AT skew', positiveColor: '#deaf6e', negativeColor: '#7294e3' }
  ]);
});

// The rows Python draws for a request's slot entries (tests/test_legend_row_facts.py
// runs the same vectors through Python's slot parsers and Legend table). The
// Web app reads the whole table of the Result: each key's palette key, slot
// and fill, and no track for a row whose palette key the vector leaves null.
// A Depth row (`depth`) names its series' slot; its series fills it (OV-1198).
const readRows = (tracks, expected, palette, shared = []) => {
  const rows = trackLegendRows(tracks, drawnRows(expected.map(([key]) => key)), new Set(shared));
  assert.equal(rows.size, expected.filter(([, paletteKey]) => paletteKey).length, 'no key beside the rows of the Result');
  return expected.map(([key, , , fill]) => {
    const row = rows.get(key);
    if (!row) return [key, null, null, fill];
    return row.kind === 'depth' ? [key, 'depth', row.slotId, fill] : [key, row.paletteKey, row.slotId, trackRowColor(row, palette)];
  });
};

for (const { name, palette, rows: expected, shared, ...request } of readJson('../fixtures/legend_track_slot_rows.json').cases) {
  test(`track rows of slot entries are Python's: ${name}`, () => {
    assert.deepEqual(readRows(drawnTracks(requestOf(request)), expected, palette, shared), expected);
  });
}

// The rows of a request without explicit slots (Python's `_generated_legend_rows`).
test('the track rows of a request without slots are Python\'s default rows', () => {
  const trackKeys = ['gc_content', 'skew_high', 'skew_low'];
  for (const { name, dinucleotide, rows } of readJson('../fixtures/legend_generated_rows.json').cases) {
    const expected = rows.map(([key, paletteKey]) => [key, trackKeys.includes(paletteKey) ? paletteKey : null, null, '']);
    for (const mode of ['circular', 'linear']) {
      assert.deepEqual(readRows(drawnTracks(requestOf({ mode, dinucleotide, slots: null })), expected, {}), expected, `${name} (${mode})`);
    }
  }
});

// A request without explicit slots draws one skew track, whose group Python
// tags with its preset slot id (`gc_skew`, whatever the dinucleotide): the
// default track answers for it, in both modes.
test('a palette edit paints the skew track of a request without slots', () => {
  for (const mode of ['circular', 'linear']) {
    for (const dinucleotide of ['GC', 'AT']) {
      const tracks = drawnTracks(requestOf({ mode, dinucleotide, slots: null }));
      assert.deepEqual(['positive', 'negative'].map((side) => skewTrackColor(tracks, 'gc_skew', side, PALETTE)),
        [PALETTE.skew_high, PALETTE.skew_low], `${mode} ${dinucleotide}`);
    }
  }
  const { actions, skewFills } = paintFixture({ default: '#d3d3d3', skew_high: '#abcdef', skew_low: '#fedcba' }, []);
  actions.applyPaletteToSvg();
  assert.deepEqual(skewFills, ['#abcdef', '#fedcba']);
});

// A Result loaded without a request (a Session 27 to 30, which the Web shows
// from its saved settings) is described by the track settings it was loaded
// with. Its SVG is the old writer's: the skew group carries no slot id.
const asSaved = (track) => ({ ...track, saved: true });
test('the tracks of a Result loaded without a request are those of its loaded settings', () => {
  for (const mode of ['circular', 'linear']) {
    for (const dinucleotide of ['GC', 'AT']) {
      const tracks = loadedTracks({ mode, dinucleotide, slots: null });
      assert.deepEqual(tracks, drawnTracks(requestOf({ mode, dinucleotide, slots: null })).map(asSaved), `${mode} ${dinucleotide}`);
      assert.deepEqual(['positive', 'negative'].map((side) => skewTrackColor(tracks, null, side, PALETTE)),
        [PALETTE.skew_high, PALETTE.skew_low], `${mode} ${dinucleotide}`);
    }
  }
  assert.deepEqual(loadedTracks(null), [], 'no loaded settings name no track');
  // Linear slots as the form holds them after Load: the tracks are those of
  // the request Generate writes from the same slots.
  const slots = normalizeLinearTrackSlots([
    { id: 'features', renderer: 'features', enabled: true, side: 'overlay', params: {} },
    { id: 'share', renderer: 'dinucleotide_content', enabled: true, side: 'below', params: { legend_label: 'Share' } },
    { id: 'own', renderer: 'dinucleotide_skew', enabled: true, side: 'below', params: { positive_color: '#ff0000', negative_color: '#0000ff' } },
    { id: 'named', renderer: 'dinucleotide_skew', enabled: true, side: 'below', params: { nt: 'GC', legend_label: 'Skew' } },
    { id: 'off', renderer: 'dinucleotide_skew', enabled: false, side: 'below', params: {} }
  ], 'AT', 'middle');
  const tracks = loadedTracks({ mode: 'linear', dinucleotide: 'AT', slots });
  assert.deepEqual(tracks, drawnTracks(requestOf({
    mode: 'linear', dinucleotide: 'AT', slots: slots.filter((slot) => slot.enabled).map(buildLinearTrackSlotPayload)
  })).map(asSaved));
  assert.deepEqual(tracks.map(({ slotId, label, positiveColor }) => [slotId, label, positiveColor]),
    [['share', 'Share', null], ['own', 'AT skew', '#ff0000'], ['named', 'Skew', null]]);
});

// The loaded settings are the Session's saved draft: a slot added after the
// last Generate did not draw the Result. Such a track claims a row in its own
// color only where Python drew the row in that color; a request's track drew
// its rows, so it claims them as before.
test('a saved slot claims a row in its own color only where the Result drew the row in it', () => {
  const slots = normalizeLinearTrackSlots([
    { id: 'features', renderer: 'features', enabled: true, side: 'overlay', params: {} },
    { id: 'own', renderer: 'dinucleotide_skew', enabled: true, side: 'below', params: { positive_color: '#ff0000', negative_color: '#0000ff' } }
  ], 'GC', 'middle');
  const saved = loadedTracks({ mode: 'linear', dinucleotide: 'GC', slots });
  const requested = drawnTracks(requestOf({ mode: 'linear', dinucleotide: 'GC', slots: slots.map(buildLinearTrackSlotPayload) }));
  const KEYS = ['GC skew (+)', 'GC skew (-)'];
  const rowColors = (tracks, drawn) => {
    const draft = draftLegendRowColors({
      rules: [], originalLegendOrder: KEYS, paletteColors: PALETTE, drawnTracks: tracks,
      pythonRows: new Map(KEYS.map((key, index) => [key, { key, color: drawn[index] }]))
    });
    return KEYS.map((key) => draft.colorOf(key));
  };
  const palette = [PALETTE.skew_high, PALETTE.skew_low];
  assert.deepEqual(rowColors(saved, palette), rowColors([], palette), 'drawn in the palette colors: the rows stay unclaimed');
  assert.deepEqual(rowColors(saved, ['#FF0000', 'blue']), ['#ff0000', '#0000ff'], 'drawn in the slot\'s colors');
  assert.deepEqual(rowColors(requested, palette), ['#ff0000', '#0000ff'], 'a request\'s track is not affected');
});

// What the tracks read of a request is checked: a request of an older schema
// may lack `tracks`, and a value that is no request names no track.
test('a request without tracks draws the default tracks, and a value that is no request names none', () => {
  const defaults = [
    { slotId: null, kind: 'content', label: 'GC content', positiveColor: null, negativeColor: null },
    { slotId: null, kind: 'skew', label: 'GC skew', positiveColor: null, negativeColor: null }
  ];
  assert.deepEqual(drawnTracks({ mode: 'circular', diagramOptions: {} }), defaults);
  assert.deepEqual(drawnTracks({ mode: 'linear', diagramOptions: { tracks: { linearTrackSlots: null } } }), defaults);
  for (const [reason, value] of [
    ['no request (a Session 27 to 30)', null], ['no diagramOptions', { mode: 'circular' }], ['an unknown mode', { mode: 'grid', diagramOptions: {} }],
    ['tracks that is no object', { mode: 'circular', diagramOptions: { tracks: [] } }],
    ['slots that are no list', { mode: 'circular', diagramOptions: { tracks: { circularTrackSlots: 'gc_skew:dinucleotide_skew' } } }],
    ['a slot entry that is no spec', { mode: 'circular', diagramOptions: { tracks: { circularTrackSlots: [7] } } }],
    ['a slot the parser rejects', { mode: 'linear', diagramOptions: { tracks: { linearTrackSlots: ['gc_skew:unknown_renderer'] } } }]
  ]) {
    assert.deepEqual(drawnTracks(value), [], reason);
  }
});

// A user may delete a Legend row: the rows that are left keep their tracks.
test('the rows left after a Legend row is deleted keep their slots', () => {
  const rows = trackLegendRows(TRACKS, drawnRows(['GC content', 'GC skew (-)', 'AT skew (+)', 'AT skew (-)', 'GC skew (-) (2)']));
  assert.deepEqual([...rows.values()].map(({ key, slotId }) => [key, slotId]), [
    ['GC content', 'gc_content'], ['GC skew (-)', 'gc_skew'], ['GC skew (-) (2)', 'gc_skew_3'], ['AT skew (+)', 'a_skew_2'], ['AT skew (-)', 'a_skew_2']
  ]);
});

// OV-488: track slots apply on Generate (R1), so the palette paints a skew
// track with the colors of the slot that drew it, not with the draft's.
test('a palette edit paints a skew track with the colors of the slot that drew it, not the draft\'s (OV-488)', () => {
  const { actions, state, skewFills } = paintFixture({ default: '#d3d3d3', skew_high: '#abcdef', skew_low: '#fedcba' }, [], {
    slots: ['gc_skew:dinucleotide_skew@positive_color=#112233']
  });
  modeDrawing(state, 'circular').adv.circular_track_slots_enabled = true;
  modeDrawing(state, 'circular').adv.circular_track_slots = [
    { id: 'gc_skew', renderer: 'dinucleotide_skew', params: { positive_color: '#00aa00', negative_color: '#445566' } }
  ];
  actions.applyPaletteToSvg();
  assert.deepEqual(skewFills, ['#112233', '#fedcba']);
});

// No color is better than another slot's: a skew track the request names no
// slot for keeps the colors it was drawn in.
test('a palette edit leaves a skew track the request names no slot for', () => {
  assert.equal(skewTrackColor(TRACKS, 'unknown_slot', 'positive', PALETTE), '');
  const { actions, skewFills } = paintFixture({ default: '#d3d3d3', skew_high: '#abcdef', skew_low: '#fedcba' }, [], { slots: ['other:dinucleotide_skew'] });
  actions.applyPaletteToSvg();
  assert.deepEqual(skewFills, ['#6dded3', '#ad72e3']);
});

// OV-1198, OV-1215: the rows that leave with a track the draft hides, by the
// predicate the track hide reads (Circular Hide GC Content / Hide GC Skew,
// Linear Show GC Content / Show GC Skew, Show Depth in both modes).
test('the rows of a hidden drawn track are the rows of its kind', () => {
  const vector = readJson('../fixtures/legend_track_slot_rows.json').cases;
  const linear = vector.find((each) => each.mode === 'linear' && each.depthTracks);
  const rows = trackLegendRows(drawnTracks(requestOf(linear)), drawnRows(linear.rows.map(([key]) => key)));
  const shown = { show_gc: true, show_skew: true, show_depth: true };
  /** @param {Record<string, unknown>} form @param {Record<string, unknown>} [adv] */
  const draft = (form, adv = {}) => ({ form, adv });
  assert.deepEqual(hiddenTrackLegendKeys(rows, 'linear', draft(shown)), []);
  assert.deepEqual(hiddenTrackLegendKeys(rows, 'linear', draft({ ...shown, show_depth: false })), ['Depth 1', 'Reads B']);
  assert.deepEqual(hiddenTrackLegendKeys(rows, 'linear', draft({ ...shown, show_gc: false })), ['GC content']);
  assert.deepEqual(hiddenTrackLegendKeys(rows, 'linear', draft({ ...shown, show_skew: false })), ['GC skew (+)', 'GC skew (-)']);
  const circular = trackLegendRows(drawnTracks(requestOf({ mode: 'circular', dinucleotide: 'GC', slots: null })),
    drawnRows(['CDS', 'GC content', 'GC skew (+)', 'GC skew (-)']));
  const circularShown = { suppress_gc: false, suppress_skew: false, show_depth: true };
  assert.deepEqual(hiddenTrackLegendKeys(circular, 'circular', draft({ ...circularShown, suppress_gc: true })), ['GC content']);
  assert.deepEqual(hiddenTrackLegendKeys(circular, 'circular', draft({ ...circularShown, suppress_skew: true })), ['GC skew (+)', 'GC skew (-)']);
  assert.equal(drawnTrackHidden('content', 'circular', draft({ show_gc: false })), false, 'Circular reads Hide GC Content');
});

// R2-H1: with a custom stack on, its enabled slots name the tracks the request
// draws; the toggles the stack disables keep a stale value that hides nothing.
test('a custom track stack, not its disabled toggles, says which drawn track is hidden', () => {
  const hiddenForm = { show_gc: false, show_skew: false, show_depth: false, suppress_gc: false, suppress_skew: false };
  const depthSlot = { id: 'depth_1', renderer: 'depth', enabled: true };
  for (const mode of /** @type {const} */ (['circular', 'linear'])) {
    const on = mode === 'linear' ? 'linear_track_slots_enabled' : 'circular_track_slots_enabled';
    const slots = mode === 'linear' ? 'linear_track_slots' : 'circular_track_slots';
    const stack = (/** @type {unknown[]} */ list) => ({ form: hiddenForm, adv: { [on]: true, [slots]: list } });
    assert.equal(draftDrawsDepth(mode, stack([depthSlot])), true, mode);
    assert.equal(drawnTrackHidden('depth', mode, stack([depthSlot])), false, mode);
    assert.equal(drawnTrackHidden('depth', mode, stack([{ ...depthSlot, enabled: false }])), true, mode);
    assert.equal(drawnTrackHidden('depth', mode, { form: { show_depth: true }, adv: { [on]: true, [slots]: [] } }), true, mode);
    assert.equal(drawnTrackHidden('depth', mode, { form: hiddenForm, adv: { [slots]: [depthSlot] } }), true, `${mode} stack off`);
  }
  const linearStack = { form: hiddenForm, adv: { linear_track_slots_enabled: true } };
  assert.equal(drawnTrackHidden('content', 'linear', linearStack), false);
  assert.equal(drawnTrackHidden('skew', 'linear', linearStack), false);
  const circularStack = { form: { suppress_gc: true, suppress_skew: true }, adv: { circular_track_slots_enabled: true } };
  assert.equal(drawnTrackHidden('content', 'circular', circularStack), true, 'Circular keeps Hide GC Content');
  assert.equal(drawnTrackHidden('skew', 'circular', circularStack), true, 'Circular keeps Hide GC Skew');
});
