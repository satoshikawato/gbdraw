// The shared vector of request schemas 1 and 2 (PR-D): the corpus of stored
// requests and the arm cases, each case a corpus request with its edits
// applied. The node tests and the cross-runtime pytest
// (tests/test_legacy_render_request_lift.py, through
// `node tests/web/helpers/legacy-request-vectors.mjs <vector file>`) read it
// here, so one applier writes each case's input.
import { readFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';

/** @param {string | URL} path the vector file */
export const readVectors = (path) => JSON.parse(readFileSync(path, 'utf8'));

/**
 * @param {unknown} root
 * @param {(string | number)[]} path
 */
export const valueAt = (root, path) => path.reduce(
  (/** @type {any} */ value, key) => (value !== null && typeof value === 'object' && Object.hasOwn(value, key) ? value[key] : undefined),
  root
);

/**
 * The stored request of one case: its corpus request with the edits applied.
 * @param {{ requests: { id: string, renderRequest: unknown }[] }} vectors
 * @param {{ base: string, edits: { path: (string | number)[], value?: unknown, delete?: boolean }[] }} vectorCase
 */
export const caseRequest = (vectors, { base, edits }) => {
  const stored = vectors.requests.find(({ id }) => id === base);
  if (!stored) throw new Error(`No corpus request ${base}.`);
  const request = structuredClone(stored.renderRequest);
  for (const edit of edits) {
    const parent = valueAt(request, edit.path.slice(0, -1));
    const key = edit.path.at(-1);
    if (parent === null || typeof parent !== 'object') throw new Error(`No parent for ${edit.path.join('.')}.`);
    if (edit.delete) delete parent[key];
    else parent[key] = structuredClone(edit.value);
  }
  return request;
};

// As a script: every corpus request and case input with its lift, or the
// code and context of the lift's refusal, as JSON.
if (process.argv[1] === fileURLToPath(import.meta.url)) {
  const vectors = readVectors(process.argv[2]);
  const { liftLegacyRenderRequest } = await import('../../../gbdraw/web/js/services/legacy-render-request.js');
  /** @param {string} id @param {unknown} renderRequest */
  const lift = (id, renderRequest) => {
    try {
      return { id, input: renderRequest, ...liftLegacyRenderRequest(renderRequest) };
    } catch (error) {
      const { code, context } = /** @type {{ code?: unknown, context?: unknown }} */ (error);
      if (typeof code !== 'string') throw error;
      return { id, input: renderRequest, refused: { code, context } };
    }
  };
  process.stdout.write(JSON.stringify({
    requests: vectors.requests.map(({ id, renderRequest }) => lift(id, renderRequest)),
    cases: vectors.cases.map((vectorCase) => lift(vectorCase.id, caseRequest(vectors, vectorCase)))
  }));
}
