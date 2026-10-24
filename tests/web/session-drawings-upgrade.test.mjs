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
    assert.deepEqual(sourceFacts, vector.sourceFacts);
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

test('the shown mode comes first and is active', () => {
  const data = fixture('legend-renamed-rows.v46.gbdraw-session.json.gz');
  data.ui.mode = 'linear';
  const { document } = upgradeModeScopedSession(data, data.version);
  assert.deepEqual(document.drawings.map((drawing) => [drawing.id, drawing.renderRequest === null]), [
    ['linear', true], ['circular', false]
  ]);
  assert.equal(document.ui.activeDrawingId, 'linear');
});

test('source facts name a Session 48 without a Web draft', () => {
  const data = settingsOnly();
  delete data.modes;
  assert.deepEqual(upgradeModeScopedSession(data, data.version).sourceFacts,
    { cardRowsAreFileRows: true, proteinResultsInDrawnFrame: false });
});

for (const version of [44, 49]) {
  test(`only Sessions 46-48 move to drawings (${version})`, () => {
    const data = settingsOnly();
    data.version = version;
    assert.throws(() => upgradeModeScopedSession(data, version), (/** @type {{ code?: string }} */ error) => error.code === 'INPUT_INVALID');
  });
}
