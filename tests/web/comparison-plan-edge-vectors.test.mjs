// @ts-check
// The Web Save check of a Linear comparison-plan edge decides every case of the vector that Python's admission
// (tests/test_comparison_plan_edge_vectors.py) decides.
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import test from 'node:test';

import { validateCurrentComparisonAuthority } from '../../gbdraw/web/js/services/session-authority.js';

const VECTORS = JSON.parse(readFileSync(new URL('../fixtures/linear_comparison_plan_edge_cases.json', import.meta.url), 'utf8'));
// The first Session whose edges may hold windowLosatFilenames (WINDOW_LOSAT_FILENAMES_SESSION_VERSION).
const SESSION_VERSION = 51;

for (const { name, edges, accepted } of /** @type {{ name: string, edges: unknown[], accepted: boolean }[]} */ (VECTORS.cases)) {
  test(`Save admits a comparison-plan edge as Python does: ${name}`, () => {
    const session = /** @type {any} */ ({
      version: SESSION_VERSION,
      modes: { linear: { config: { linearComparisonPlan: { mode: 'selected', defaultSource: 'losat', edges } } } }
    });
    if (accepted) assert.doesNotThrow(() => validateCurrentComparisonAuthority(session));
    else assert.throws(() => validateCurrentComparisonAuthority(session));
  });
}
