// @ts-check
// The Web reading of a Session 46 Legend's color and stroke edits by row key
// (Session 47), against the vectors it shares with the Python twin
// (tests/fixtures/sessions/legend-row-key-vectors.json, tests/test_session_io.py).
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import test from 'node:test';

import { legendEditsByRowKey, sessionLegendEditsByRowKey } from '../../gbdraw/web/js/services/session-drawings-upgrade.js';
import { MODE_SCOPED_SESSION_VERSION, ROW_KEYED_LEGEND_SESSION_VERSION, SESSION_VERSION } from '../../gbdraw/web/js/contracts/schema.generated.js';

/** @type {{ readsSessionVersion: number, cases: { name: string, legend: Record<string, unknown>, expected: Record<string, unknown>, outcomes: unknown[] }[] }} */
const vectors = JSON.parse(readFileSync(new URL('../fixtures/sessions/legend-row-key-vectors.json', import.meta.url), 'utf8'));

test('the vectors read Session 46, the version before the row-keyed Legend', () => {
  assert.equal(vectors.readsSessionVersion, MODE_SCOPED_SESSION_VERSION);
  assert.equal(ROW_KEYED_LEGEND_SESSION_VERSION, vectors.readsSessionVersion + 1);
});

for (const { name, legend, expected, outcomes } of vectors.cases) {
  test(`legend row key vector: ${name}`, () => {
    const source = structuredClone(legend);
    assert.deepEqual(legendEditsByRowKey(source), { legend: expected, outcomes });
    assert.deepEqual(source, legend);
  });
}

test('a Session 46 reads each mode slice by row key and keeps the rest of the document', () => {
  const [chain] = vectors.cases;
  const session = {
    version: 46,
    editorState: { legend: { originalOrder: ['GC content', 'GC skew'] } },
    modes: {
      circular: { config: { form: {} }, editorState: { legend: chain.legend, featureStrokes: { overrides: {} } } },
      linear: { config: { form: {} } }
    }
  };
  const { session: read, outcomes } = sessionLegendEditsByRowKey(session, 46);
  assert.deepEqual(outcomes, chain.outcomes.map((outcome) => ({ ...outcome, slice: 'circular' })));
  assert.deepEqual(read.modes.circular.editorState.legend, chain.expected);
  assert.deepEqual(read.modes.circular.editorState.featureStrokes, { overrides: {} });
  assert.deepEqual(read.modes.linear, session.modes.linear);
  assert.deepEqual(read.editorState, session.editorState);
  assert.equal(read.version, 46);
});

test('a Session 27-44 reads its flat draft by row key; a Session 47 is read as it is', () => {
  const [chain] = vectors.cases;
  const flat = { version: 44, editorState: { legend: chain.legend } };
  const older = sessionLegendEditsByRowKey(flat, 44);
  assert.deepEqual(older.session.editorState.legend, chain.expected);
  assert.deepEqual(older.outcomes, chain.outcomes.map((outcome) => ({ ...outcome, slice: 'draft' })));
  const current = { version: SESSION_VERSION, modes: { circular: { editorState: { legend: chain.legend } } } };
  assert.deepEqual(sessionLegendEditsByRowKey(current, SESSION_VERSION), { session: current, outcomes: [] });
  assert.equal(sessionLegendEditsByRowKey(current, SESSION_VERSION).session, current);
});
