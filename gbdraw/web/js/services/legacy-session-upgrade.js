// @ts-check
// Load of a Session 27-44 (DESIGN §4.3): Python moves it to a Session 51 in
// one Worker call (`legacy_session_to_drawings_json`,
// gbdraw/web_support/legacy_session.py); the Web has no converter of its own.
// Load of a Session 46-50 whose Linear request holds a LOSATP search: Python
// gives the request the results of the LOSAT rows the file saved, as it does
// in a Session 27-44's move (`saved_protein_results_json`, OV-3545).
// This module holds the port the composition root sets, which runs the
// operation each step names, and reads the replies.
import { diagnosticError } from '../utils/error-normalization.js';
import { readSentLoadRepair } from '../utils/load-repairs.js';

/** @import { LoadRepair } from '../utils/load-repairs.js' */
/** @import { SessionEnvelope } from '../contracts/session.generated.js' */

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
/** What Load sends Python for a Session 46-50: the file's document as the file holds it. @typedef {{ document: SessionEnvelope }} SavedProteinResultsRequest */
/**
 * The Worker port of Load's Python steps for an older Session: it runs the
 * helper operation it is given with that operation's request.
 * @typedef {{
 *   (operation: 'upgradeLegacySession', request: LegacySessionUpgradeRequest): Promise<unknown>,
 *   (operation: 'storeSavedProteinResults', request: SavedProteinResultsRequest): Promise<unknown>
 * }} LegacySessionUpgrader
 */

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
  const response = await upgrade(LEGACY_SESSION_UPGRADE.operation, {
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

const SAVED_PROTEIN_RESULTS = /** @type {const} */ ({ stage: 'helper', operation: 'storeSavedProteinResults' });
/** Codes of a failed exchange with the Worker, which stop the Load. */
const PROTOCOL_FAILURES = new Set(['RESULT_INVALID', 'HELPER_PROTOCOL']);

/**
 * Python refused the file (its reader is stricter than the Web's): an error
 * of this operation that is neither a protocol failure nor a Cancel.
 * @param {unknown} error
 */
const refusedByPython = (error) => isRecord(error) && error.operation === SAVED_PROTEIN_RESULTS.operation
  && !PROTOCOL_FAILURES.has(String(error.code)) && error.canceled !== true;

/**
 * The requests of a Session 46-50 that take the results of its saved LOSAT
 * rows, from Python, once (OV-3545); Load replaces them in their drawings.
 * null when Python refuses the file: Load then reads it as before, with its
 * search. Another Worker failure throws its own diagnostic; a reply of
 * another shape, or one that adds a resource the file holds, throws
 * RESULT_INVALID of this operation; without the port, HELPER_PROTOCOL.
 * @param {SessionEnvelope} document The file's document as the file holds it.
 * @param {{ upgrade?: LegacySessionUpgrader | null }} [options]
 * @returns {Promise<SavedProteinResultRequests | null>}
 */
export const readSavedProteinResults = async (document, { upgrade = legacySessionUpgrader } = {}) => {
  if (!upgrade) throw diagnosticError('HELPER_PROTOCOL', {}, SAVED_PROTEIN_RESULTS);
  /** @type {unknown} */
  let response;
  try {
    response = await upgrade(SAVED_PROTEIN_RESULTS.operation, { document });
  } catch (error) {
    if (refusedByPython(error)) return null;
    throw error;
  }
  const reply = isRecord(response) ? response.result : null;
  const held = isRecord(document.resources) ? document.resources : {};
  if (!isRecord(reply) || !isRecord(reply.requests) || !isRecord(reply.resources)
    || !Object.values(reply.requests).every(isRecord)
    || !Object.entries(reply.resources).every(([id, resource]) => isRecord(resource) && !Object.hasOwn(held, id))) {
    throw diagnosticError('RESULT_INVALID', {}, SAVED_PROTEIN_RESULTS);
  }
  return /** @type {SavedProteinResultRequests} */ (reply);
};
