// The catalog keys of the user-facing error text (D-17). The wording owner,
// gbdraw/web/js/utils/error-normalization.js, keeps its English in tables that
// other tests read as text (`CODE: '…'`), so its values carry no msg()
// marker. It has no imports, so a copy loaded as a data URL can return the
// values of its private tables; every text of an error text part comes from
// one of them (tests/web/error-text.test.mjs checks that).
import { readFileSync } from 'node:fs';

const SOURCE = new URL('../../../gbdraw/web/js/utils/error-normalization.js', import.meta.url);

const PROBE = `
export const errorTextKeyTables = {
  DEFINITIONS: Object.values(DEFINITIONS).flatMap(([message, , guidance]) => (guidance ? [message, guidance] : [message])),
  REASONS: Object.values(REASONS),
  FIELD_LABELS: Object.values(FIELD_LABELS),
  ORDINAL_LABELS: Object.values(ORDINAL_LABELS),
  SESSION_TABLE_LABELS: Object.values(SESSION_TABLE_LABELS),
  SWITCHED_TRACKS: Object.values(SWITCHED_TRACKS).map(({ remedy }) => remedy),
  KEPT_INSIDE_AMONG: Object.values(KEPT_INSIDE_AMONG),
  PHRASES: Object.values(PHRASES),
  OPERATION_TITLES: Object.values(OPERATION_TITLES),
  RECOVERY_GUIDANCE: Object.values(RECOVERY_GUIDANCE)
};
export const errorVocabulary = {
  codes: Object.keys(DEFINITIONS), reasons: Object.keys(REASONS), labelledFields: Object.keys(FIELD_LABELS),
  sessionTables: Object.keys(SESSION_TABLE_LABELS), switchedTracks: Object.keys(SWITCHED_TRACKS),
  operations: Object.keys(OPERATION_TITLES), recoveries: Object.keys(RECOVERY_GUIDANCE)
};`;

const probe = () => import(`data:text/javascript,${encodeURIComponent(`${readFileSync(SOURCE, 'utf8')}\n${PROBE}`)}`);

/**
 * The keys of each table of the wording owner.
 * @returns {Promise<Record<string, string[]>>}
 */
export const errorTextKeyTables = async () => (await probe()).errorTextKeyTables;

/**
 * The codes, reasons, labelled fields, Session tables, switched tracks,
 * operations and recoveries the wording owner words.
 * @returns {Promise<Record<string, string[]>>}
 */
export const errorVocabulary = async () => (await probe()).errorVocabulary;

/**
 * Every catalog key of the error text, once.
 * @returns {Promise<Set<string>>}
 */
export const errorTextKeys = async () => new Set(Object.values(await errorTextKeyTables()).flat());
