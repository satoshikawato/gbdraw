// @ts-check
// The queries of the Sequence search panel (PD-OI-117): the query text, read by Python into queries (name, length,
// colour), their removal and colours, and the queries Use as search query adds. They are the panel's own, transient
// values: not in the Session and no History step.
import { normalizeUserFacingError } from '../../utils/error-normalization.js';
import { panelQueries, readSearchQueriesReply } from '../../services/sequence-search.js';

/** @import { SearchQuery } from '../../services/sequence-search.js' */

/** @typedef {ReturnType<typeof normalizeUserFacingError>} SearchError */

/**
 * The fields of the panel's state that the queries owner writes.
 * @typedef {object} SearchQueryState
 * @property {string} queryText The query box as typed.
 * @property {SearchQuery[]} queries The queries of the latest read.
 * @property {SearchError | null} queryError Why the latest read refused the text.
 * @property {boolean} reading A read runs.
 */

/**
 * Runs a Worker helper operation (`runDiagramHelperOperation`).
 * @typedef {(operation: string, payload: Record<string, unknown>) => Promise<{ result?: unknown }>} RunHelper
 */

/**
 * @param {SearchQueryState} panel
 * @param {RunHelper} runHelper
 * @param {string} operation The Worker helper that reads the queries (`readSearchQueries`).
 */
export const createSearchQueries = (panel, runHelper, operation) => {
  let readSequence = 0;
  // The text of the latest read and that read: a click on Search blurs the query box, whose change event asks for the
  // same text again, so a text is read once and Search waits for its read.
  let readFor = '';
  /** @type {Promise<void>} */
  let pendingRead = Promise.resolve();
  /** @param {string} text */
  const readText = async (text) => {
    const sequence = ++readSequence;
    panel.queryError = null;
    if (text.trim() === '') {
      panel.queries = [];
      return;
    }
    panel.reading = true;
    try {
      const bytes = new TextEncoder().encode(text).buffer;
      const { result } = await runHelper(operation, { files: [{ role: 'source', bytes }] });
      if (sequence !== readSequence) return;
      panel.queries = panelQueries(readSearchQueriesReply(result), panel.queries);
    } catch (error) {
      if (sequence !== readSequence) return;
      panel.queries = [];
      panel.queryError = normalizeUserFacingError(error, { stage: 'helper', operation });
    } finally {
      if (sequence === readSequence) panel.reading = false;
    }
  };
  /** @param {string} text */
  const readQueries = (text) => {
    panel.queryText = text;
    if (text !== readFor) {
      readFor = text;
      pendingRead = readText(text);
    }
    return pendingRead;
  };
  /** A text not read yet, or whose read still runs: Search reads it first (OV-2579). */
  const unread = () => panel.queryText !== readFor || panel.reading;
  /** The query box as typed, before its change event reads it. @param {string} text */
  const editQueryText = (text) => { panel.queryText = text; };
  /** @param {{ target?: HTMLInputElement | null } | null | undefined} event */
  const chooseQueryFile = async (event) => {
    const file = event?.target?.files?.[0];
    if (!file) return;
    try {
      await readQueries(await file.text());
    } finally {
      if (event?.target) event.target.value = '';
    }
  };
  /**
   * Adds a query, named by its FASTA header, to the query text and reads it again. A query of that name is already
   * listed: nothing is added twice.
   * @param {{ name: string, fasta: string }} query
   */
  const addQuery = async ({ name, fasta }) => {
    if (panel.queries.some((query) => query.name === name)) return;
    const text = panel.queryText.trim() === '' ? '' : panel.queryText.replace(/\s*$/, '\n');
    await readQueries(`${text}${fasta}`);
  };
  /**
   * Removes a query from the list and from the query text, which then holds the FASTA of the queries left (as read), so
   * a later read or Use as search query does not bring it back.
   * @param {string} name
   */
  const removeQuery = (name) => {
    panel.queries = panel.queries.filter((query) => query.name !== name);
    panel.queryText = panel.queries.map((query) => query.fasta).join('');
    // A read still running for an earlier text is dropped.
    readSequence += 1;
    panel.reading = false;
    readFor = panel.queryText;
    pendingRead = Promise.resolve();
  };
  /** Forgets the reads, and drops one still running (a Session Load empties the panel). */
  const clearReads = () => {
    readSequence += 1;
    readFor = '';
    pendingRead = Promise.resolve();
  };
  /** @param {string} name @param {string} color */
  const setQueryColor = (name, color) => {
    panel.queries = panel.queries.map((query) => (query.name === name ? { ...query, color } : query));
  };
  return { readQueries, unread, editQueryText, chooseQueryFile, addQuery, removeQuery, setQueryColor, clearReads };
};
