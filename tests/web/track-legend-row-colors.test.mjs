// OV-484: a track's Legend row shows the color of the track slot that draws
// it. The Legend panel, the live compile and the palette's track painter read
// it from the tracks of the request that drew the Result
// (services/track-legend-rows.js), never from the shape of the row's caption.
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { test } from 'node:test';

import { compileDirectEditorMutationPlan, draftLegendPanelColors } from '../../gbdraw/web/js/app/candidate-render.js';
import { drawnTracks, skewTrackColor, trackLegendRows, trackRowColor } from '../../gbdraw/web/js/services/track-legend-rows.js';
import { fixture as paintFixture } from './helpers/svg-style-fixture.mjs';

const readJson = (path) => JSON.parse(readFileSync(new URL(path, import.meta.url), 'utf8'));
const HMMT_AT_SKEW = readJson('../../gbdraw/web/gallery/sessions/HmmtDNA_ATskew.gbdraw-session.json');
const requestOf = ({ mode, dinucleotide, slots }) => ({
  mode, diagramOptions: { dinucleotide, tracks: mode === 'linear' ? { linearTrackSlots: slots } : { circularTrackSlots: slots } }
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
const readRows = (tracks, expected, palette) => {
  const rows = trackLegendRows(tracks, expected.map(([key]) => key));
  assert.equal(rows.size, expected.filter(([, paletteKey]) => paletteKey).length, 'no key beside the rows of the Result');
  return expected.map(([key, , , fill]) => {
    const row = rows.get(key);
    return row ? [key, row.paletteKey, row.slotId, trackRowColor(row, palette)] : [key, null, null, fill];
  });
};

for (const { name, palette, rows: expected, ...request } of readJson('../fixtures/legend_track_slot_rows.json').cases) {
  test(`track rows of slot entries are Python's: ${name}`, () => {
    assert.deepEqual(readRows(drawnTracks(requestOf(request)), expected, palette), expected);
  });
}

// The rows of a request without explicit slots (Python's `_generated_legend_fills`).
test('the track rows of a request without slots are Python\'s default rows', () => {
  const trackKeys = ['gc_content', 'skew_high', 'skew_low'];
  for (const { name, dinucleotide, rows } of readJson('../fixtures/legend_generated_rows.json').cases) {
    const expected = rows.map(([key, paletteKey]) => [key, trackKeys.includes(paletteKey) ? paletteKey : null, null, '']);
    for (const mode of ['circular', 'linear']) {
      assert.deepEqual(readRows(drawnTracks(requestOf({ mode, dinucleotide, slots: null })), expected, {}), expected, `${name} (${mode})`);
    }
  }
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
    ['no request', null], ['no diagramOptions', { mode: 'circular' }], ['an unknown mode', { mode: 'grid', diagramOptions: {} }],
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
  const rows = trackLegendRows(TRACKS, ['GC content', 'GC skew (-)', 'AT skew (+)', 'AT skew (-)', 'GC skew (-) (2)']);
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
  state.drawings.circular.adv.circular_track_slots_enabled = true;
  state.drawings.circular.adv.circular_track_slots = [
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
