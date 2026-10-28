// @ts-check
// Load of a Session 27-44 (DESIGN §4.3): Python moves it to a Session 50 in
// one Worker call (`legacy_session_to_drawings_json`,
// gbdraw/web_support/legacy_session.py); the Web has no converter of its own.
// This module holds the port the composition root sets and reads the reply.
import { diagnosticError } from '../utils/error-normalization.js';
import { readSentLoadRepair } from '../utils/load-repairs.js';

/** @import { LoadRepair } from '../utils/load-repairs.js' */

/**
 * Python's reply: the Session 50, the stored Results a Session 50 cannot
 * hold (Sessions 27-39), by drawing ID, where the file held the Result set
 * each drawing took (`resultSetSources`, by drawing ID: the file's one set is
 * at its top level, `''`), the Load repairs of the draft migrations, and how
 * many `hash=` annotation targets moved to their source feature, which the
 * notice after Load names (`LegacySessionDrawings`,
 * gbdraw/web_support/legacy_session.py). The document is read by Load's
 * band-50 reader, as a file's would be.
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
 * The Session 50 of a Session 27-44, from Python, once. A Worker failure
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
