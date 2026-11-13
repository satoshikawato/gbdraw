// @ts-check
// The Web's move of a Session 46-48 to drawings (`upgradeModeScopedSession`,
// services/session-drawings-upgrade.js) against the `upgrades` cases of
// tests/fixtures/sessions/drawings-vectors.json, which Python's
// `drawings_document_from_mode_scoped` also reads (tests/test_session_drawings_vectors.py).
import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { gunzipSync } from 'node:zlib';

const { upgradeModeScopedSession } = await import('../../gbdraw/web/js/services/session-drawings-upgrade.js');
const root = new URL('../../', import.meta.url);
const vectors = JSON.parse(readFileSync(new URL('tests/fixtures/sessions/drawings-vectors.json', root), 'utf8'));

/**
 * What a vector's `circular` pins of each Circular drawing (`_circular_cards`,
 * tests/test_session_drawings_vectors.py): its File cards, files by resource
 * ID, and the record selectors of its draft.
 * @param {Record<string, any>} document
 */
const circularCards = (document) => document.drawings
  .filter((/** @type {Record<string, any>} */ drawing) => drawing.mode === 'circular')
  .map((/** @type {Record<string, any>} */ drawing) => {
    const rows = drawing.webFiles?.bindings?.circularSources;
    const config = drawing.config || {};
    return {
      id: drawing.id,
      circularSources: rows === undefined ? null : rows.map((/** @type {Record<string, any>} */ row) => ({
        uid: row.uid, gb: row.gb?.resourceId ?? null, gff: row.gff?.resourceId ?? null, fasta: row.fasta?.resourceId ?? null
      })),
      recordsOff: config.recordsOff ?? null,
      circular_record_selector: config.form?.circular_record_selector ?? null,
      multi_record_positions: config.adv?.multi_record_positions ?? null
    };
  });

for (const vector of vectors.upgrades) {
  test(`a Session 46-50 moves to the vector drawings: ${vector.name}`, () => {
    const input = typeof vector.input === 'string'
      ? JSON.parse(gunzipSync(readFileSync(new URL(vector.input, root))).toString('utf8'))
      : structuredClone(vector.input);
    for (const { path: [...parents], value } of vector.edits || []) {
      const key = parents.pop();
      parents.reduce((holder, part) => holder[part], input)[key] = value;
    }
    const { document, sourceFacts } = upgradeModeScopedSession(input, input.version);
    // The Circular input moves to File cards (circular-multi-file D-01, D-08).
    assert.deepEqual(circularCards(document), vector.circular || []);
    assert.deepEqual(
      document.drawings.map((/** @type {Record<string, unknown>} */ drawing) => (
        { id: drawing.id, mode: drawing.mode, keys: Object.keys(drawing).sort() }
      )),
      vector.drawings
    );
    assert.equal(document.ui.activeDrawingId, vector.activeDrawingId);
    // The vectors hold the facts Python's twin writes; every Session 46-48
    // also saved records of Python's paint with its Results, stored no rings
    // setting, and stored a Depth series tick never set as null (the Web's own
    // facts).
    assert.deepEqual(sourceFacts,
      { ...vector.sourceFacts, resultsRecordPythonPaint: true, ringSettingStored: false, depthTicksUnsetAsNull: true });
  });
}

// The rule cases of tests/test_session_drawings_vectors.py, on the same fixtures.
/** @param {string} name @returns {Record<string, any>} */
const fixture = (name) => JSON.parse(gunzipSync(readFileSync(new URL(`tests/fixtures/sessions/${name}`, root))).toString('utf8'));
const settingsOnly = () => fixture('file-card-imports-circular.v48.gbdraw-session.json.gz');
/** @param {Record<string, any>} data */
const idsOf = (data) => upgradeModeScopedSession(structuredClone(data), data.version).document.drawings
  .map((/** @type {{ id: string }} */ drawing) => drawing.id);

test('other-mode settings count only when they differ from its default and the shown mode', () => {
  const data = settingsOnly();
  assert.deepEqual(idsOf(data), ['circular']);
  data.modes.linear.config.colors = data.modes.circular.config.colors;
  assert.deepEqual(idsOf(data), ['circular']);
  data.modes.linear.config.form.show_labels_linear = 'all';
  assert.deepEqual(idsOf(data), ['circular', 'linear']);
});

test('an other mode that holds the default Depth series is not used when the shown mode has two', () => {
  // OV-1901 (D-15-122): 0.14.0 wrote one default Depth series in every mode,
  // so rule d sees it as the mode's default; Python's rule d agrees
  // (tests/test_session_drawings_vectors.py, the same mutation).
  const data = settingsOnly();
  const series = { label: '', color: '#4A90E2', height: null, large_tick_interval: null, small_tick_interval: null, tick_font_size: null };
  data.modes.circular.config.adv.depth_tracks = [series, { ...series, color: '#E45756' }];
  data.modes.linear.config.adv.depth_tracks = [series];
  assert.deepEqual(idsOf(data), ['circular']);
  data.modes.linear.config.adv.depth_tracks = [];
  assert.deepEqual(idsOf(data), ['circular', 'linear']);
});

test('an other mode that holds the Depth series Sessions 46 and 47 named is not used', () => {
  // main wrote this default before unnamed series kept an empty label (5f3aa020d),
  // in `tests/fixtures/sessions/legend-renamed-rows.v46...`; Python's rule d agrees
  // (tests/test_session_drawings_vectors.py, the same mutation).
  const data = settingsOnly();
  const series = { label: 'Depth', color: '#4A90E2', height: null, large_tick_interval: null, small_tick_interval: null, tick_font_size: null };
  data.modes.circular.config.adv.depth_tracks = [{ ...series, color: '#E45756' }];
  data.modes.linear.config.adv.depth_tracks = [series];
  assert.deepEqual(idsOf(data), ['circular']);
  data.modes.linear.config.adv.depth_tracks = [{ ...series, label: 'Reads' }];
  assert.deepEqual(idsOf(data), ['circular', 'linear']);
});

// OV-3114: 0.14.5 and earlier wrote repeat_region, misc_feature and default
// as one gray in the default palette. A Session 46-50's other mode that holds
// those colors is not used (Python's rule d agrees, the same mutation in
// tests/test_session_drawings_vectors.py); a drawing of the current app that
// holds them has colors that differ from its default, so Save keeps it.
const EARLIER_DEFAULT_COLORS = Object.freeze({
  CDS: '#54bcf8', rRNA: '#71ee7d', tRNA: '#e8b441', tmRNA: '#ded44e', ncRNA: '#c4fac3',
  repeat_region: '#d3d3d3', misc_feature: '#d3d3d3', default: '#d3d3d3',
  skew_high: '#6dded3', skew_low: '#ad72e3', gc_content: '#a1a1a1',
  pairwise_match: '#d3d3d3', pairwise_match_min: '#FFE7E7', pairwise_match_max: '#FF7272',
  collinear_block_plus_min: '#f0f1f5', collinear_block_plus: '#8b9cc1',
  collinear_block_minus_min: '#FFE7E7', collinear_block_minus: '#E15759'
});

test('an other mode that holds the default colors of 0.14.5 is not used', () => {
  const data = settingsOnly();
  data.modes.linear.config.colors = { ...EARLIER_DEFAULT_COLORS };
  assert.deepEqual(idsOf(data), ['circular']);
  data.modes.linear.config.colors = { ...EARLIER_DEFAULT_COLORS, default: '#123456' };
  assert.deepEqual(idsOf(data), ['circular', 'linear']);
});

test('Save keeps a drawing whose colors are the default colors of 0.14.5', async () => {
  const { sliceDiffersFromDefault } = await import('../../gbdraw/web/js/services/session-drawings-upgrade.js');
  const { MODE_SLICE_DEFAULTS } = await import('../../gbdraw/web/js/mode-slice-defaults.generated.js');
  const current = MODE_SLICE_DEFAULTS.defaults.linear.find((row) => row.domain === 'config' && row.path === 'colors');
  assert.equal(sliceDiffersFromDefault({ config: { colors: { ...EARLIER_DEFAULT_COLORS } } }, 'linear'), true);
  assert.equal(sliceDiffersFromDefault({ config: { colors: current?.value } }, 'linear'), false);
});

test('a Linear record layout counts only when it arranges records', () => {
  const data = settingsOnly();
  data.modes.linear.config.linearRecordLayout.rows = [{ row: 1, uid: 'linear-seq-a' }, { row: 2, uid: 'linear-seq-b' }];
  assert.deepEqual(idsOf(data), ['circular']);
  data.modes.linear.config.linearRecordLayout.rows[1].row = 1;
  assert.deepEqual(idsOf(data), ['circular', 'linear']);
});

test('other-mode keyed rows count entry by entry', () => {
  const data = settingsOnly();
  const [[key, entry]] = Object.entries(fixture('rule-captions-with-color-id.v46.gbdraw-session.json.gz')
    .modes.circular.features.featureColorOverrides);
  data.modes.circular.features.featureColorOverrides = { [key]: entry };
  data.modes.linear.features.featureColorOverrides = { [key]: entry };
  assert.deepEqual(idsOf(data), ['circular']);
  data.modes.linear.features.featureColorOverrides = { [key]: entry, [`${key}-2`]: entry };
  assert.deepEqual(idsOf(data), ['circular', 'linear']);
});

test('an other-mode input binding counts and an empty sequence row does not', () => {
  const data = settingsOnly();
  assert.deepEqual(idsOf(data), ['circular']);
  data.webFiles.bindings.linearSeqs[0].depth = { ...data.webFiles.bindings.d_color };
  assert.deepEqual(idsOf(data), ['circular', 'linear']);
});

test('a shown mode without a Result opens the other mode\'s Result first', () => {
  // The shown-at-Save mode keeps its drawing; the drawing with
  // the stored Result opens and comes first.
  const data = fixture('legend-renamed-rows.v46.gbdraw-session.json.gz');
  data.ui.mode = 'linear';
  const { document } = upgradeModeScopedSession(data, data.version);
  assert.deepEqual(document.drawings.map((drawing) => [drawing.id, drawing.renderRequest === null]), [
    ['circular', false], ['linear', true]
  ]);
  assert.equal(document.ui.activeDrawingId, 'circular');
});

test('an integral float equals its integer default', () => {
  // OV-1668: a Python writer's `2.0` parses as `2`; Python's rule d agrees
  // (tests/test_session_drawings_vectors.py, the same mutation).
  const data = settingsOnly();
  data.modes.linear.config.adv.linear_track_slots_schema_version = 2;
  const text = JSON.stringify(data).replace('"linear_track_slots_schema_version":2', '"linear_track_slots_schema_version":2.0');
  assert.ok(text.includes('2.0'));
  assert.deepEqual(idsOf(JSON.parse(text)), ['circular']);
});

test('source facts name a Session 48 without a Web draft', () => {
  const data = settingsOnly();
  delete data.modes;
  assert.deepEqual(upgradeModeScopedSession(data, data.version).sourceFacts,
    { cardRowsAreFileRows: true, proteinResultsInDrawnFrame: false, resultsRecordPythonPaint: true, ringSettingStored: false,
      depthTicksUnsetAsNull: true });
});

for (const version of [44, 51]) {
  test(`only Sessions 46-50 move to drawings (${version})`, () => {
    const data = settingsOnly();
    data.version = version;
    assert.throws(() => upgradeModeScopedSession(data, version), (/** @type {{ code?: string }} */ error) => error.code === 'INPUT_INVALID');
  });
}
