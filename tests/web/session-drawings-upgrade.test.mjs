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

for (const vector of vectors.upgrades) {
  test(`a Session 46-48 moves to the vector drawings: ${vector.name}`, () => {
    const input = typeof vector.input === 'string'
      ? JSON.parse(gunzipSync(readFileSync(new URL(vector.input, root))).toString('utf8'))
      : structuredClone(vector.input);
    const { document, sourceFacts } = upgradeModeScopedSession(input, input.version);
    assert.deepEqual(
      document.drawings.map((/** @type {Record<string, unknown>} */ drawing) => (
        { id: drawing.id, mode: drawing.mode, keys: Object.keys(drawing).sort() }
      )),
      vector.drawings
    );
    assert.equal(document.ui.activeDrawingId, vector.activeDrawingId);
    // The vectors hold the facts Python's twin writes; every Session 46-48
    // also saved records of Python's paint with its Results (the Web's own fact).
    assert.deepEqual(sourceFacts, { ...vector.sourceFacts, resultsRecordPythonPaint: true });
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
    { cardRowsAreFileRows: true, proteinResultsInDrawnFrame: false, resultsRecordPythonPaint: true });
});

for (const version of [44, 50]) {
  test(`only Sessions 46-49 move to drawings (${version})`, () => {
    const data = settingsOnly();
    data.version = version;
    assert.throws(() => upgradeModeScopedSession(data, version), (/** @type {{ code?: string }} */ error) => error.code === 'INPUT_INVALID');
  });
}
