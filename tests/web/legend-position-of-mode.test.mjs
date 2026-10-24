// OV-652 (BA-D07): a legend position is of the set of its diagram mode. Python
// names the sets (gbdraw/api/options.py), `contracts/schema.generated.js`
// carries them, and `legendPositionOfMode` in services/layout-preferences.js
// is the one Web reader: a value outside the mode's set reads like an absent
// one, so every site keeps the fallback it has for an empty value.
import assert from 'node:assert/strict';
import { readFile, readdir } from 'node:fs/promises';
import { join } from 'node:path';
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

const webRoot = join(process.cwd(), 'gbdraw', 'web');
const { state, createDefaultForm } = await import('../../gbdraw/web/js/state.js');
const layout = await import('../../gbdraw/web/js/services/layout-preferences.js');
const { SETTINGS_SLICE } = await import('../../gbdraw/web/js/services/drawing/settings.js');
const { LEGEND_POSITIONS } = await import('../../gbdraw/web/js/contracts/schema.generated.js');
const { applyUiStateData, installShownDrawing } = await import('../../gbdraw/web/js/services/config.js');
const { buildCanonicalRenderRequest } = await import('../../gbdraw/web/js/services/session-request.js');
const { resolveLinearComparisonPlan } = await import('../../gbdraw/web/js/services/linear-comparisons.js');

const storedLinear = (legend) => layout.resolveActiveLayoutPreference(
  layout.normalizeLayoutPreferences({ linear: { legend, plotTitlePosition: 'bottom' } }), 'linear'
).legend;
const storedCircular = (single, multi) => layout.normalizeLayoutPreferences({
  circular: { single: { legend: single }, multi: { legend: multi } }
}).circular;

test('the generated sets are the positions Python accepts for each mode', () => {
  assert.deepEqual(LEGEND_POSITIONS, {
    circular: ['left', 'right', 'top', 'bottom', 'upper_left', 'upper_right', 'lower_left', 'lower_right', 'none'],
    linear: ['left', 'right', 'top', 'bottom', 'none']
  });
});

test('the owner returns a position of the mode, or null', () => {
  for (const mode of ['circular', 'linear']) {
    for (const position of LEGEND_POSITIONS[mode]) {
      assert.equal(layout.legendPositionOfMode(position, mode), position);
      assert.equal(layout.legendPositionOfMode(` ${position.toUpperCase()} `, mode), position);
    }
    for (const outside of ['banana', '', '   ', null, undefined, 5, true, ['top'], { legend: 'top' }, 'constructor']) {
      assert.equal(layout.legendPositionOfMode(outside, mode), null, `${mode} ${JSON.stringify(outside)}`);
    }
  }
  assert.equal(layout.legendPositionOfMode('upper_left', 'linear'), null);
  assert.equal(layout.legendPositionOfMode('upper_left', 'circular'), 'upper_left');
});

test('a stored Linear position outside the Linear set reads as the Linear default (the OV-652 probe)', () => {
  assert.equal(storedLinear('Upper_Left'), 'bottom');
  assert.equal(storedLinear(' banana '), 'bottom');
  assert.equal(storedLinear('top'), 'top');
  assert.equal(storedLinear(' Top '), 'top');
  assert.equal(storedLinear(''), 'bottom');
});

test('a stored Circular corner stays, and a value outside the Circular set reads as the Circular default', () => {
  assert.equal(storedCircular('Upper_Left', null).single.legend, 'upper_left');
  assert.equal(storedCircular(' lower_right ', null).single.legend, 'lower_right');
  assert.equal(storedCircular('banana', null).single.legend, 'left');
});

test('a Circular multi-record position outside the set follows the single-record position', () => {
  const circular = storedCircular('upper_right', ' banana ');
  assert.equal(circular.multi.legend, null);
  assert.equal(layout.resolveCircularLayoutPreference({ circular }, true).legend, 'upper_right');
  assert.equal(storedCircular('upper_right', ' Bottom ').multi.legend, 'bottom');
  assert.equal(storedCircular('upper_right', '').multi.legend, null);
});

test('the legacy migration takes the next source when a field holds a value outside the set', () => {
  // Without a stored single-record position, the single-record slot takes the active position.
  assert.deepEqual(layout.migrateLegacyLayoutPreferences(
    { circularLegendPosition: 'banana', legend: 'Upper_Right' }
  ).circular, {
    single: { legend: 'left', plotTitlePosition: 'none' },
    multi: { legend: 'upper_right', plotTitlePosition: 'none' }
  });
  assert.equal(layout.migrateLegacyLayoutPreferences(
    { circularLegendPosition: 'banana', legend: 'banana' }, { activeLegend: 'top' }
  ).circular.single.legend, 'top');
  assert.equal(layout.migrateLegacyLayoutPreferences(
    { circularSingleRecordLegendPosition: 'right', circularMultiRecordLegendPosition: 'banana', circularLegendPosition: 'top' }
  ).circular.multi.legend, 'top');
  assert.equal(layout.migrateLegacyLayoutPreferences(
    { linearLegendPosition: 'upper_left' }, { mode: 'linear', activeLegend: 'top' }
  ).linear.legend, 'top');
  assert.equal(layout.migrateLegacyLayoutPreferences(
    { linearLegendPosition: 'upper_left' }, { mode: 'circular', activeLegend: 'upper_left' }
  ).linear.legend, 'bottom');
  // A Session that release 0.13.0 saved in Linear after a Circular configuration was loaded.
  const released = layout.restoredLayoutPreferences(
    {
      mode: 'linear', legend: 'upper_left', linearLegendPosition: 'upper_left', circularLegendPosition: 'left',
      circularSingleRecordLegendPosition: 'left', circularMultiRecordLegendPosition: 'left'
    },
    { mode: 'linear', multiRecord: false, active: { legend: 'bottom', plotTitlePosition: 'bottom' } }
  );
  assert.deepEqual(
    [released.linear.legend, released.circular.single.legend, released.circular.multi.legend],
    ['bottom', 'left', 'left']
  );
  // The Linear fallback to `ui.legend` still applies when `linearLegendPosition` holds no Linear position.
  assert.equal(layout.restoredLayoutPreferences(
    { mode: 'linear', legend: 'top', linearLegendPosition: 'upper_left' },
    { mode: 'linear', multiRecord: false, active: { legend: 'bottom', plotTitlePosition: 'bottom' } }
  ).linear.legend, 'top');
});

test('a legacy legend field outside its mode does not count as a stored layout of a canonical Session', () => {
  const restore = (ui) => layout.restoredLayoutPreferences(ui, {
    mode: 'linear', multiRecord: false,
    projected: { linear: { legend: 'top', plotTitlePosition: 'top' } },
    active: { legend: 'bottom', plotTitlePosition: 'bottom' }
  });
  const unset = restore({});
  for (const ui of [
    { legend: 'banana' }, { legend: '' }, { linearLegendPosition: 'upper_left' },
    { circularLegendPosition: 'banana' }, { circularMultiRecordLegendPosition: 'centre' }
  ]) {
    assert.deepEqual(restore(ui), unset, JSON.stringify(ui));
  }
  // A position of the field's mode is still a stored layout: it migrates.
  assert.equal(restore({ legend: ' Left ' }).linear.legend, 'left');
  // `legend` is also a Circular source: a corner keeps the Circular slot, and the Linear slot reads it as absent.
  const corner = restore({ legend: 'Upper_Left' });
  assert.deepEqual([corner.circular.multi.legend, corner.linear.legend], ['upper_left', 'top']);
  assert.equal(restore({ circularSingleRecordLegendPosition: 'upper_right' }).circular.single.legend, 'upper_right');
  // The plot title fields keep their presence check.
  assert.notDeepEqual(restore({ linearPlotTitlePosition: 'top' }), unset);
});

test('a drawing never holds a position outside its mode: the Position control and a loaded slot', () => {
  const linear = modeDrawing(state, 'linear');
  linear.form.legend = 'top';
  linear.form.legend = 'Upper_Left';
  assert.equal(linear.form.legend, 'bottom');
  assert.equal(linear.layoutPreferences.linear.legend, 'bottom');
  SETTINGS_SLICE.install(linear.slices.settings, { layoutPreferences: layout.normalizeLayoutPreferences({ linear: { legend: 'Upper_Left' } }) });
  assert.equal(linear.form.legend, 'bottom');
  const circular = modeDrawing(state, 'circular');
  circular.form.multi_record_canvas = false;
  circular.form.legend = ' Lower_Left ';
  assert.equal(circular.form.legend, 'lower_left');
  circular.form.legend = 'banana';
  assert.equal(circular.form.legend, 'left');
});

test('a stored generated position outside its mode reads as the side an unset one had', () => {
  const read = (generatedMode, generatedLegendPosition) => {
    state.generatedLegendPosition.value = 'right';
    applyUiStateData(modeDrawing(state, generatedMode), { generatedMode, generatedLegendPosition });
    return state.generatedLegendPosition.value;
  };
  assert.equal(read('linear', 'Upper_Left'), 'bottom');
  assert.equal(read('linear', ' Top '), 'top');
  assert.equal(read('circular', 'Upper_Left'), 'upper_left');
  assert.equal(read('circular', 'banana'), 'left');
  assert.equal(read('circular', ''), 'right');
});

const genbankText = `LOCUS       LEGENDPOS                  40 bp    DNA     linear   UNK 01-JAN-1980
DEFINITION  Legend position fixture.
ACCESSION   LEGENDPOS
VERSION     LEGENDPOS
KEYWORDS    .
SOURCE      .
  ORGANISM  .
            .
FEATURES             Location/Qualifiers
     CDS             1..30
                     /product="test protein"
ORIGIN
        1 atgcatgcat gcatgcatgc atgcatgcat gcatgcatgc
//
`;
const genbank = () => ({
  name: 'legend-position.gb', type: 'text/plain', size: new TextEncoder().encode(genbankText).byteLength,
  lastModified: 0, data: btoa(genbankText)
});
// The request of `mode` from `drawing` (default: that mode's drawing), as Generate builds it.
const requestLegend = (mode, drawing = modeDrawing(state, mode)) => {
  installShownDrawing(mode);
  modeDrawing(state, 'circular').inputType.value = 'gb';
  modeDrawing(state, 'linear').inputType.value = 'gb';
  drawing.circularRecordList.value = [];
  const filesData = mode === 'circular'
    ? { c_gb: genbank(), linearSeqs: [] }
    : {
        linearSeqs: [{
          uid: 'record-1', gb: genbank(), depth: null, losat_gencode: 1, region_record_id: '',
          region_start: null, region_end: null, region_reverse: false
        }],
        linearComparisons: []
      };
  const comparisonPlanSnapshot = mode === 'linear'
    ? resolveLinearComparisonPlan({
        plan: drawing.linearComparisonPlan, sequences: filesData.linearSeqs, layout: [],
        losatProgram: drawing.losatProgram.value, blastpMode: drawing.losat.blastp.mode
      })
    : null;
  return buildCanonicalRenderRequest({ state, drawing, filesData, comparisonPlanSnapshot })
    .renderRequest.diagramOptions.output.legend;
};
// A drawing whose `form.legend` is a plain value, as a test or a tool may pass one.
const plainDrawing = (mode, legend) => ({
  ...modeDrawing(state, mode), form: { ...createDefaultForm(), ...modeDrawing(state, mode).form, legend }
});

test('the request writer writes a position of the request\'s mode for every drawing', () => {
  for (const [mode, fallback] of [['linear', 'bottom'], ['circular', 'left']]) {
    modeDrawing(state, mode).form.multi_record_canvas = false;
    // A slot written without the owner, as a Session Load wrote it before OV-652.
    const slot = mode === 'linear'
      ? modeDrawing(state, 'linear').layoutPreferences.linear
      : modeDrawing(state, 'circular').layoutPreferences.circular.single;
    slot.legend = ' banana ';
    assert.equal(requestLegend(mode), fallback, `${mode}: a slot outside the set`);
    slot.legend = 'top';
    assert.equal(requestLegend(mode), 'top');
    for (const legend of [' banana ', '', undefined, null, 5]) {
      assert.equal(requestLegend(mode, plainDrawing(mode, legend)), fallback, `${mode} ${JSON.stringify(legend)}`);
    }
    assert.equal(requestLegend(mode, plainDrawing(mode, ' Right ')), 'right');
  }
  assert.equal(requestLegend('linear', plainDrawing('linear', 'Upper_Left')), 'bottom');
  assert.equal(requestLegend('circular', plainDrawing('circular', 'Upper_Left')), 'upper_left');
});

test('the Position control lists the generated set of each mode', async () => {
  const html = await readFile(join(webRoot, 'index.html'), 'utf8');
  const select = html.match(/<select id="legend-position"[\s\S]*?<\/select>/);
  assert.ok(select, 'index.html has the Legend Position control');
  const options = [...select[0].matchAll(/<option value="([^"]*)"([^>]*)>/g)]
    .map(([, value, attributes]) => ({ value, circularOnly: attributes.includes(`v-if="mode === 'circular'"`) }));
  assert.deepEqual(options.map(({ value }) => value).sort(), [...LEGEND_POSITIONS.circular].sort());
  assert.deepEqual(
    options.filter(({ circularOnly }) => !circularOnly).map(({ value }) => value).sort(),
    [...LEGEND_POSITIONS.linear].sort()
  );
});

test('no other Web module normalizes a legend position', async () => {
  const modules = (await readdir(join(webRoot, 'js'), { recursive: true }))
    .filter((path) => path.endsWith('.js') && !path.startsWith('workers'));
  const twins = [];
  for (const path of modules) {
    const source = await readFile(join(webRoot, 'js', path), 'utf8');
    if (/\bnormalize\w*LegendPosition\b/.test(source)) twins.push(path);
  }
  assert.deepEqual(twins, []);
  const writer = await readFile(join(webRoot, 'js', 'services', 'session-request.js'), 'utf8');
  assert.doesNotMatch(writer, /legend: String\(drawing\.form\.legend/);
});
