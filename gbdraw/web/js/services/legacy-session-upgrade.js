// @ts-check
// Load of a Session 27-44 (DESIGN §4.3): Python moves it to a Session 51 in
// one Worker call (`legacy_session_to_drawings_json`,
// gbdraw/web_support/legacy_session.py); the Web has no converter of its own.
// Load of a Session 46-50 whose Linear request holds a LOSATP search: Python
// gives the request the results of the LOSAT rows the file saved, as it does
// in a Session 27-44's move (`saved_protein_results_json`, OV-3545).
// This module holds the ports the composition root sets and reads the replies.
import { diagnosticError } from '../utils/error-normalization.js';
import { readSentLoadRepair } from '../utils/load-repairs.js';
import { isModeScopedBand, storedResultSets, withStoredResultSet } from './session-authority.js';

/** @import { LoadRepair } from '../utils/load-repairs.js' */
/** @import { SessionChecked, SessionEnvelope } from '../contracts/session.generated.js' */
/** @import { UnvalidatedFields } from './session-authority.js' */

/**
 * Python's reply: the Session 51, the stored Results a Session 51 cannot
 * hold (Sessions 27-39), by drawing ID, where the file held the Result set
 * each drawing took (`resultSetSources`, by drawing ID: the file's one set is
 * at its top level, `''`), the Load repairs of the draft migrations, and how
 * many `hash=` annotation targets moved to their source feature, which the
 * notice after Load names (`LegacySessionDrawings`,
 * gbdraw/web_support/legacy_session.py). The document is read by Load's
 * band-51 reader, as a file's would be.
 * @typedef {object} LegacySessionUpgrade
 * @property {Record<string, unknown>} document
 * @property {ReadonlyMap<string, readonly unknown[]>} storedResults
 * @property {ReadonlyMap<string, ''>} resultSetSources
 * @property {readonly LoadRepair[]} warnings
 * @property {number} movedAnnotationTargets
 */

/**
 * What Load sends Python: the file's document as its readers read it, and
 * the schema of the request the file stored (the reader sends it read), or
 * null for the document's own.
 * @typedef {{ document: unknown, storedRequestSchema: number | null }} LegacySessionUpgradeRequest
 */
/** @typedef {(request: LegacySessionUpgradeRequest) => Promise<unknown>} LegacySessionUpgrader */

/** @type {LegacySessionUpgrader | null} */
let legacySessionUpgrader = null;
/** @param {LegacySessionUpgrader | null} upgrader */
export const setLegacySessionUpgrader = (upgrader) => {
  legacySessionUpgrader = typeof upgrader === 'function' ? upgrader : null;
};

const LEGACY_SESSION_UPGRADE = /** @type {const} */ ({ stage: 'helper', operation: 'upgradeLegacySession' });

/** @param {unknown} value @returns {value is Record<string, unknown>} */
const isRecord = (value) => Boolean(value) && typeof value === 'object' && !Array.isArray(value);

/**
 * The Session 51 of a Session 27-44, from Python, once. A Worker failure
 * throws its own diagnostic; a reply of another shape throws RESULT_INVALID
 * of this operation; without the port, HELPER_PROTOCOL.
 * @param {unknown} document The file's document as the Session input reader admitted it.
 * @param {{ storedRequestSchema?: unknown, upgrade?: LegacySessionUpgrader | null }} [options]
 *   `storedRequestSchema`: the schema of the request the file stored (OV-2052).
 * @returns {Promise<LegacySessionUpgrade>}
 */
export const upgradeLegacySession = async (document, { storedRequestSchema = null, upgrade = legacySessionUpgrader } = {}) => {
  if (!upgrade) throw diagnosticError('HELPER_PROTOCOL', {}, LEGACY_SESSION_UPGRADE);
  const response = await upgrade({
    document,
    storedRequestSchema: Number.isSafeInteger(storedRequestSchema) ? Number(storedRequestSchema) : null
  });
  const reply = isRecord(response) ? response.result : null;
  const invalid = () => diagnosticError('RESULT_INVALID', {}, LEGACY_SESSION_UPGRADE);
  if (!isRecord(reply) || !isRecord(reply.document) || !isRecord(reply.storedResults)
    || !isRecord(reply.resultSetSources) || !Array.isArray(reply.warnings)
    || !Number.isSafeInteger(reply.movedAnnotationTargets) || Number(reply.movedAnnotationTargets) < 0) {
    throw invalid();
  }
  /** @type {Map<string, readonly unknown[]>} */
  const storedResults = new Map();
  for (const [drawingId, results] of Object.entries(reply.storedResults)) {
    if (!Array.isArray(results)) throw invalid();
    storedResults.set(drawingId, results);
  }
  /** @type {Map<string, ''>} */
  const resultSetSources = new Map();
  for (const [drawingId, at] of Object.entries(reply.resultSetSources)) {
    if (at !== '') throw invalid();
    resultSetSources.set(drawingId, at);
  }
  const warnings = reply.warnings.map(readSentLoadRepair);
  if (warnings.some((warning) => warning === null)) throw invalid();
  return {
    document: reply.document,
    storedResults,
    resultSetSources,
    warnings: /** @type {LoadRepair[]} */ (warnings),
    movedAnnotationTargets: Number(reply.movedAnnotationTargets)
  };
};

/**
 * Python's reply for a Session 46-50 (`SavedProteinResultRequests`,
 * gbdraw/web_support/legacy_session.py): each committed request that took the
 * results of the file's saved LOSAT rows, by drawing ID (a Session 46-50
 * names a drawing by its mode), and the resources they add (the result
 * tables), by ID, each with its bytes in base64 as a JSON Session holds them.
 * @typedef {object} SavedProteinResultRequests
 * @property {Readonly<Record<string, Record<string, unknown>>>} requests
 * @property {Readonly<Record<string, Record<string, unknown>>>} resources
 */
/** What Load sends Python: the file's document as the file holds it. @typedef {{ document: SessionEnvelope }} SavedProteinResultsRequest */
/** @typedef {(request: SavedProteinResultsRequest) => Promise<unknown>} SavedProteinResultReader */

/** @type {SavedProteinResultReader | null} */
let savedProteinResultReader = null;
/** @param {SavedProteinResultReader | null} reader */
export const setSavedProteinResultReader = (reader) => {
  savedProteinResultReader = typeof reader === 'function' ? reader : null;
};

const SAVED_PROTEIN_RESULTS = /** @type {const} */ ({ stage: 'helper', operation: 'storeSavedProteinResults' });

/**
 * Whether a stored request is a Linear one whose LOSATP comparison is a
 * search; Python decides whether the saved rows answer it.
 * @param {unknown} stored
 */
const searchesProteins = (stored) => {
  const request = /** @type {UnvalidatedFields<'mode' | 'comparisons'>} */ (stored);
  return request?.mode === 'linear' && Array.isArray(request.comparisons) && request.comparisons.some((item) => {
    const comparison = /** @type {UnvalidatedFields<'kind' | 'mode'>} */ (item);
    return comparison?.kind === 'generatedProteinComparison' && comparison.mode !== 'none';
  });
};

/**
 * A Session 46-50 as Load reads it: a Linear drawing whose committed LOSATP
 * search the file's saved LOSAT rows answer holds their results, as a
 * Generate commits them, so the browser's render redraws it before the first
 * Generate (OV-3545). One Worker call, only for a file whose Linear request
 * holds a search and that saved LOSAT rows; any other document is returned as
 * given. A Worker failure throws its own diagnostic; a reply of another shape
 * throws RESULT_INVALID of this operation; without the port, HELPER_PROTOCOL.
 * @param {SessionEnvelope} document The file's document as the file holds it.
 * @param {{ read?: SavedProteinResultReader | null }} [options]
 * @returns {Promise<SessionEnvelope>}
 */
export const withSavedProteinResults = async (document, { read = savedProteinResultReader } = {}) => {
  if (!isModeScopedBand(/** @type {SessionChecked} */ (/** @type {unknown} */ (document)))) return document;
  const searches = storedResultSets(document).filter(({ set }) => searchesProteins(set.renderRequest));
  // A Session 46-50 keeps its LOSAT caches at the top level, which every drawing reads.
  const { losatCache } = /** @type {{ losatCache?: unknown }} */ (document);
  const rows = /** @type {UnvalidatedFields<'entries'>} */ (losatCache)?.entries;
  if (searches.length === 0 || !Array.isArray(rows) || rows.length === 0) return document;
  if (!read) throw diagnosticError('HELPER_PROTOCOL', {}, SAVED_PROTEIN_RESULTS);
  const response = await read({ document });
  const reply = isRecord(response) ? response.result : null;
  const invalid = () => diagnosticError('RESULT_INVALID', {}, SAVED_PROTEIN_RESULTS);
  const held = isRecord(document.resources) ? document.resources : {};
  if (!isRecord(reply) || !isRecord(reply.requests) || !isRecord(reply.resources)
    || !Object.values(reply.requests).every(isRecord)
    || !Object.entries(reply.resources).every(([id, resource]) => isRecord(resource) && !Object.hasOwn(held, id))) {
    throw invalid();
  }
  const { requests, resources } = /** @type {SavedProteinResultRequests} */ (reply);
  if (Object.keys(requests).length === 0) return document;
  // The projected set comes first (`storedResultSets`), so replacing it keeps the other mode's set.
  let next = document;
  const replaced = new Set();
  for (const { at, set } of searches) {
    const drawingId = String(/** @type {UnvalidatedFields<'mode'>} */ (set.renderRequest)?.mode);
    if (!Object.hasOwn(requests, drawingId)) continue;
    next = withStoredResultSet(next, at, { ...set, renderRequest: requests[drawingId] });
    replaced.add(drawingId);
  }
  // Python returns a request only for a drawing Load sent with a search.
  if (replaced.size !== Object.keys(requests).length) throw invalid();
  return { ...next, resources: { ...held, ...resources } };
};
