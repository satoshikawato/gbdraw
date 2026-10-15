// @ts-check
import { CANONICAL_REQUEST_SCHEMA } from '../contracts/schema.generated.js';

/**
 * Whether a committed or projected request is of the current schema. Load
 * reads every stored request to it except schema 1, which it keeps as stored
 * until that schema's reader lifts it; the Gallery migration projects a stored
 * schema 1, 2 or 5 before it rebuilds the request. Each of those answers
 * false, as the schema tests this predicate replaced did (`>= 6` to `>= 9`).
 * It goes with the schema-1 arms.
 * @param {unknown} request
 * @returns {boolean}
 */
export const isCurrentRequestSchema = (request) => (
  request !== null && typeof request === 'object' && !Array.isArray(request)
  && /** @type {{ schema?: unknown }} */ (request).schema === CANONICAL_REQUEST_SCHEMA
);
