// Legacy protein artifacts a written Session keeps, shared with Python (OV-1132, OV-1143):
// tests/fixtures/legacy_protein_artifact_vectors.json, read by tests/test_legacy_protein_artifact_vectors.py too.
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { test } from 'node:test';

import {
  serializableLegacyProteinCandidateEnvelope,
  writtenLegacyProteinArtifacts
} from '../../gbdraw/web/js/services/losat-cache.js';

const { cases } = JSON.parse(readFileSync(new URL('../fixtures/legacy_protein_artifact_vectors.json', import.meta.url), 'utf8'));

/** @param {string} state @param {number} index */
const candidate = (state, index) => ({
  schema: 1,
  state,
  rejectionReason: state === 'rejected' ? 'no match' : null,
  originalEntry: { schema: 2, kind: 'raw-losat', key: `legacy-${index}`, program: 'blastp', text: '' }
});

for (const item of cases) {
  test(`legacy protein artifacts: ${item.id}`, () => {
    const evidence = Array.from({ length: item.evidence }, (_unused, index) => ({
      schema: 1, kind: 'derived-losatp-payload', key: `evidence-${index}`, payload: {}
    }));
    const written = writtenLegacyProteinArtifacts({
      candidates: { schema: 1, entries: item.candidates.map(candidate) },
      evidence: { schema: 1, entries: evidence },
      promoted: item.promoted,
      hasCurrentDerived: item.currentDerived
    });
    const kept = serializableLegacyProteinCandidateEnvelope(written.candidates).entries;
    assert.deepEqual(kept.map((entry) => entry.state), item.written.candidates);
    assert.deepEqual(kept.map((entry) => entry.originalEntry.key),
      item.candidates.map((_state, index) => `legacy-${index}`).filter((_key, index) => !item.promoted.includes(index)));
    assert.deepEqual(written.evidence.entries, evidence.slice(0, item.written.evidence));
  });
}
