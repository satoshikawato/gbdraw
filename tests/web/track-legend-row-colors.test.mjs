// OV-484: a track's Legend row shows the color of the track slot that draws
// it. The Legend panel, the live compile and the palette's track painter read
// it from the track rows of the request that drew the Result
// (services/track-legend-rows.js), never from the shape of the row's caption.
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { test } from 'node:test';

import { compileDirectEditorMutationPlan, draftLegendPanelColors } from '../../gbdraw/web/js/app/candidate-render.js';
import { fixture as paintFixture } from './helpers/svg-style-fixture.mjs';

const readJson = (path) => JSON.parse(readFileSync(new URL(path, import.meta.url), 'utf8'));
const HMMT_AT_SKEW = readJson('../../gbdraw/web/gallery/sessions/HmmtDNA_ATskew.gbdraw-session.json');
const trackRowsOf = async (request) => (await import('../../gbdraw/web/js/services/track-legend-rows.js')).trackLegendRows(request);

// The track rows of the HmmtDNA_ATskew Gallery Session, written out: the AT
// skew slot has its own colors, and a second GC skew slot repeats the label.
const trackRow = (key, slotId, paletteKey, ownColor = null) => [key, { key, slotId, paletteKey, ownColor }];
const TRACK_ROWS = new Map([
  trackRow('GC content', 'gc_content', 'gc_content'),
  trackRow('GC skew (+)', 'gc_skew', 'skew_high'), trackRow('GC skew (-)', 'gc_skew', 'skew_low'),
  trackRow('AT skew (+)', 'a_skew_2', 'skew_high', '#deaf6e'), trackRow('AT skew (-)', 'a_skew_2', 'skew_low', '#7294e3'),
  trackRow('GC skew (+) (2)', 'gc_skew_3', 'skew_high'), trackRow('GC skew (-) (2)', 'gc_skew_3', 'skew_low', '#ff00ff')
]);
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
    legendEntries, originalLegendOrder, legendColorOverrides: {}, rules: [], pythonRows, paletteColors: PALETTE, trackRows: TRACK_ROWS
  });
  assert.deepEqual(trackColors([...panel.listed]), trackColors(DRAWN), 'each track row shows the color Python drew it in');
});

const liveLegendFills = (paletteColors, trackRows = TRACK_ROWS) => {
  const { legendEntries, originalLegendOrder, pythonRows } = legendState();
  const operations = compileDirectEditorMutationPlan({
    catalogAdmission: {
      resultNames: ['diagram.svg'], renderedTargetsByOverrideKey: new Map(), resultIndexesByRenderedId: new Map(),
      renderedFeaturesByResult: [new Map()]
    },
    legendEntries, originalLegendOrder,
    livePreview: { domains: ['legendFills'], paletteColors, drawnContext: null, pythonRows, features: [], trackRows }
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
  assert.deepEqual(trackColors(liveLegendFills({ ...PALETTE, skew_high: '#112233' }, new Map())), []);
});

test('the track rows of the HmmtDNA_ATskew Gallery Session name the slot and the own colors of each row', async () => {
  const rows = await trackRowsOf(HMMT_AT_SKEW.renderRequest);
  const read = (key) => { const row = rows.get(key); return row && [row.slotId, row.paletteKey, row.ownColor]; };
  assert.deepEqual(read('GC content'), ['gc_content', 'gc_content', null]);
  assert.deepEqual(read('GC skew (+)'), ['gc_skew', 'skew_high', null]);
  assert.deepEqual(read('AT skew (+)'), ['a_skew_2', 'skew_high', '#deaf6e']);
  assert.deepEqual(read('AT skew (-)'), ['a_skew_2', 'skew_low', '#7294e3']);
});

// The slot rows Python draws (tests/test_legend_row_facts.py runs the same
// vectors through Python's Legend table).
test('the track rows read from a request\'s slots are Python\'s, with Python\'s fills', async () => {
  const { trackRowColor } = await import('../../gbdraw/web/js/services/track-legend-rows.js');
  for (const { name, mode, dinucleotide, palette, slots, rows: expected } of readJson('../fixtures/legend_track_slot_rows.json').cases) {
    const tracks = mode === 'linear' ? { linearTrackSlots: slots } : { circularTrackSlots: slots };
    const rows = await trackRowsOf({ mode, diagramOptions: { dinucleotide, tracks } });
    for (const [key, paletteKey, slotId, fill] of expected) {
      const row = rows.get(key);
      assert.ok(row, `${name}: ${key}`);
      assert.deepEqual([row.paletteKey, row.slotId, trackRowColor(row, palette)], [paletteKey, slotId, fill], `${name}: ${key}`);
    }
  }
});

// The rows of a request without explicit slots (Python's `_generated_legend_fills`).
test('the track rows of a request without slots are Python\'s default rows', async () => {
  for (const { name, dinucleotide, rows: expected } of readJson('../fixtures/legend_generated_rows.json').cases) {
    for (const mode of ['circular', 'linear']) {
      const rows = await trackRowsOf({ mode, diagramOptions: { dinucleotide, tracks: { circularTrackSlots: null, linearTrackSlots: null } } });
      for (const [key, paletteKey] of expected.filter(([, palette]) => ['gc_content', 'skew_high', 'skew_low'].includes(palette))) {
        assert.deepEqual([rows.get(key)?.paletteKey, rows.get(key)?.slotId], [paletteKey, null], `${name} (${mode}): ${key}`);
      }
      assert.equal(rows.has('Depth'), false, 'no palette key colors Depth');
    }
  }
  assert.equal((await trackRowsOf(null)).size, 0, 'no request, no track rows');
  const olderSchema = await trackRowsOf({ mode: 'circular', diagramOptions: {} });
  assert.deepEqual([...olderSchema.keys()], ['GC content', 'GC content (+)', 'GC content (-)', 'GC skew', 'GC skew (+)', 'GC skew (-)']);
});

// OV-488: track slots apply on Generate (R1), so the palette paints a skew
// track with the colors of the slot that drew it, not with the draft's.
test('a palette edit paints a skew track with the colors of the slot that drew it, not the draft\'s (OV-488)', () => {
  const drawn = new Map([...TRACK_ROWS].map(([key, row]) => [key, row.slotId === 'gc_skew' && row.paletteKey === 'skew_high' ? { ...row, ownColor: '#112233' } : row]));
  const { actions, state, skewFills } = paintFixture({ default: '#d3d3d3', skew_high: '#abcdef', skew_low: '#fedcba' }, [], { trackRows: drawn });
  state.drawings.circular.adv.circular_track_slots_enabled = true;
  state.drawings.circular.adv.circular_track_slots = [
    { id: 'gc_skew', renderer: 'dinucleotide_skew', params: { positive_color: '#00aa00', negative_color: '#445566' } }
  ];
  actions.applyPaletteToSvg();
  assert.deepEqual(skewFills, ['#112233', '#fedcba']);
});
