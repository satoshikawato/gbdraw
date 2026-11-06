// @ts-check
// The Sequence search section (PD-OI-116): the panel's queries and settings,
// the search, and the hit table of the latest search. The queries and settings
// are the panel's own, transient values: not in the Session, no History step,
// and never written to Pairwise Comparisons. A search runs LOSAT through the
// Pairwise Comparisons runtime (`runSearchLosat`), Python turns its rows into
// an annotation set (`services/sequence-search.js`), and adding that set is one
// History step of the shown drawing that the automatic rerender draws at once
// (`searchSetIdsOfStep`, the live projection in services/session-request.js).
import { computed, reactive } from '../utils/reactivity.js';
import { msg, t, tm, tn } from '../utils/i18n.js';
import { normalizeUserFacingError } from '../utils/error-normalization.js';
import { comparisonStateForMode, resolveComparisonThresholds } from '../mode-profiles.js';
import { draftAnnotationSetsOfRequest } from '../services/annotation-state.js';
import { downloadTextFile } from '../services/text-download.js';
import { SEARCH_GRADIENT_END_COLOR, SEARCH_GRADIENT_START_COLOR } from '../utils/sequence-search-colors.generated.js';
import {
  nextSearchSetId, outlinedHitIds, panelQueries, readSearchQueriesReply, readSequenceSearchReply, searchSetIdsOfStep,
  searchSpans, searchSubjects
} from '../services/sequence-search.js';

/** @import { SearchColorMode, SearchHit, SearchQuery, SearchStepChange, SequenceSearchInput } from '../services/sequence-search.js' */
/** @import { RegionAxisRecord } from '../services/region-axis.js' */
/** @import { RecordSpanInput, RecordSpanSequencesReply } from '../services/diagram-worker-protocol.js' */
/** @import { AnnotationSet } from '../services/annotation-state.js' */
/** @import { DrawingId } from '../services/drawing-view.js' */

/** @typedef {ReturnType<typeof normalizeUserFacingError>} SearchError */
/** @typedef {{ min_bitscore: string | number, evalue: string | number, identity: string | number, alignment_length: string | number }} SearchThresholdDraft */

/**
 * What the latest search found: the set it added (null when no hit passed the thresholds), its hits and queries.
 * @typedef {object} SearchOutcome
 * @property {DrawingId} drawingId
 * @property {string | null} setId
 * @property {SearchHit[]} hits
 * @property {{ name: string, length: number, color: string }[]} queries
 * @property {number} rowsLeftOut Rows LOSAT wrote that the thresholds removed.
 * @property {Record<string, string>} recordIds The record ID of each searched record, by its request key (as Python read
 *   the spans), which the hit table shows.
 * @property {string} rawText The rows as LOSAT wrote them (Save Raw LOSAT TSV).
 */

/**
 * The ports the composition root gives the panel.
 * @typedef {object} SequenceSearchPorts
 * @property {{ activeDrawing: () => { id: DrawingId, mode: 'circular' | 'linear', annotationSets: readonly { id: string }[] } }} state
 * @property {() => RegionAxisRecord[]} readAxisRecords Every drawn record of the shown drawing's Results.
 * @property {(spans: readonly RecordSpanInput[]) => Promise<RecordSpanSequencesReply>} readRecordSpans
 * @property {(search: { queryFasta: string, subjectFasta: string, task: string }) => Promise<{ status: 'ok', text: string } | { status: 'canceled' } | { status: 'busy', reason: string }>} runSearchLosat
 * @property {(operation: string, payload: Record<string, unknown>) => Promise<{ result?: unknown }>} runHelper
 * @property {{ readQueries: string, buildHitSet: string }} helperOperations The Worker helper operations (`readSearchQueries`,
 *   `searchHitAnnotationSet` in services/diagram-worker-protocol.js) the panel runs.
 * @property {(label: string, write: () => unknown) => unknown} addSetStep Runs `write` as one History step of the shown drawing.
 * @property {(drawingId: DrawingId, set: AnnotationSet) => void} addSet The annotation-sets slice writer of that drawing.
 */

export const SEARCH_TASKS = Object.freeze(['megablast', 'blastn', 'dc-megablast']);

/** @param {SequenceSearchPorts} ports */
export const createSequenceSearch = ({
  state, readAxisRecords, readRecordSpans, runSearchLosat, runHelper, helperOperations, addSetStep, addSet
}) => {
  const panel = reactive({
    queryText: '',
    /** @type {SearchQuery[]} */
    queries: [],
    /** @type {SearchError | null} */
    queryError: null,
    reading: false,
    task: 'megablast',
    /** @type {SearchColorMode} */
    colorMode: 'query',
    /** The thresholds of each mode, from that mode's Pairwise Comparisons defaults. @type {Record<'circular' | 'linear', SearchThresholdDraft>} */
    thresholds: { circular: { ...comparisonStateForMode('circular') }, linear: { ...comparisonStateForMode('linear') } },
    running: false,
    /** @type {SearchError | null} */
    searchError: null,
    /** Why the latest Search did not run (a Generate or a Session load is in progress); '' otherwise. */
    busy: '',
    /** @type {SearchOutcome | null} */
    outcome: null,
    folded: false,
    /** The query whose hits the table shows ('' for all), and the hit outlined on the drawing. */
    chip: '',
    /** @type {string} */
    selectedHit: ''
  });
  // The sets each drawing's searches added, by id, with the set as written: a History step that adds one of them
  // (its own, its Redo) or whose Undo removes it draws at once.
  /** @type {Map<DrawingId, Map<string, string>>} */
  const searched = new Map();
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
      const { result } = await runHelper(helperOperations.readQueries, { files: [{ role: 'source', bytes }] });
      if (sequence !== readSequence) return;
      panel.queries = panelQueries(readSearchQueriesReply(result), panel.queries);
    } catch (error) {
      if (sequence !== readSequence) return;
      panel.queries = [];
      panel.queryError = normalizeUserFacingError(error, { stage: 'helper', operation: helperOperations.readQueries });
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
   * Use as search query (Region popup): adds a span's bases, read by Python on the source strand and named by their
   * FASTA header (`<record>:<start>-<end>`), to the query text and reads it again; the panel unfolds to show it. A
   * query of that name is already listed: nothing is added twice.
   * @param {{ name: string, fasta: string }} query
   */
  const addQuery = async ({ name, fasta }) => {
    panel.folded = false;
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
  /** @param {string} name @param {string} color */
  const setQueryColor = (name, color) => {
    panel.queries = panel.queries.map((query) => (query.name === name ? { ...query, color } : query));
  };

  // Why Search is off, in the order the user meets them; '' when it can run.
  const searchDisabledReason = computed(() => {
    // Only a Result names the drawn bases (its record axis, as for Select region); the axis is reactive, the committed
    // Session that every Result comes from is not, so the axis alone decides.
    if (readAxisRecords().length === 0) return t('Generate to search');
    // A text not read yet is read when Search runs.
    const unread = panel.queryText !== readFor;
    if (panel.queryError && !unread) return t('Correct the query to search');
    if (panel.queries.length === 0 && (!unread || panel.queryText.trim() === '')) return t('Paste or choose a query FASTA');
    return '';
  });

  /** @param {Parameters<typeof normalizeUserFacingError>[0]} error */
  const failed = (error) => {
    panel.searchError = normalizeUserFacingError(error, { stage: 'losat', operation: helperOperations.buildHitSet });
  };

  const search = async () => {
    if (searchDisabledReason.value || panel.running) return;
    panel.searchError = null;
    panel.busy = '';
    panel.running = true;
    try {
      await readQueries(panel.queryText);
      // The read may refuse the text or find no query; its reason then shows.
      if (searchDisabledReason.value) return;
      const drawing = state.activeDrawing();
      const queries = panel.queries.map(({ name, length, color }) => ({ name, length, color }));
      const queryFasta = panel.queries.map((query) => query.fasta).join('');
      const thresholds = resolveComparisonThresholds(panel.thresholds[drawing.mode], drawing.mode);
      const records = readAxisRecords();
      const read = await readRecordSpans(searchSpans(records));
      const run = await runSearchLosat({ queryFasta, subjectFasta: read.fasta, task: panel.task });
      // A cancelled search, or one a Generate or a Session load kept from running, adds nothing and keeps the Result.
      if (run.status === 'busy') panel.busy = tm(run.reason);
      if (run.status !== 'ok' || state.activeDrawing().id !== drawing.id) return;
      const setId = nextSearchSetId(drawing.annotationSets);
      /** @type {SequenceSearchInput} */
      const input = {
        table: run.text,
        setId,
        subjects: searchSubjects(read, records),
        thresholds: {
          evalue: Number(thresholds.evalue), bitscore: Number(thresholds.bitscore),
          identity: Number(thresholds.identity), alignmentLength: Number(thresholds.alignmentLength)
        },
        program: 'blastn',
        queries: queries.map(({ name, length }) => ({ name, length })),
        colorMode: panel.colorMode,
        queryColors: Object.fromEntries(queries.map(({ name, color }) => [name, color]))
      };
      const reply = readSequenceSearchReply((await runHelper(helperOperations.buildHitSet, { search: input })).result);
      if (state.activeDrawing().id !== drawing.id) return;
      const [set] = reply.annotationSet ? draftAnnotationSetsOfRequest([/** @type {Partial<AnnotationSet>} */ (reply.annotationSet)]) : [];
      panel.outcome = {
        drawingId: drawing.id, setId: set ? set.id : null, hits: reply.hits, queries,
        rowsLeftOut: reply.rowsRead - reply.rowsKept, rawText: run.text,
        recordIds: Object.fromEntries(read.spans.map((span) => [span.recordKey, span.recordId]))
      };
      panel.chip = '';
      panel.selectedHit = '';
      panel.folded = true;
      if (!set) return;
      const sets = searched.get(drawing.id) ?? new Map();
      searched.set(drawing.id, sets.set(set.id, JSON.stringify(set)));
      addSetStep(msg('Sequence search'), () => addSet(drawing.id, set));
    } catch (error) {
      failed(error);
    } finally {
      panel.running = false;
    }
  };

  /**
   * The sets a History step of the shown drawing added that its searches wrote (`searchSetIdsOfStep`).
   * @param {readonly SearchStepChange[] | null | undefined} changes
   */
  const searchSetIdsOf = (changes) => searchSetIdsOfStep(changes, searched.get(state.activeDrawing().id) ?? new Map());

  // The latest search of the shown drawing, while its set exists (a search without hits has no set).
  const shownOutcome = computed(() => {
    const outcome = panel.outcome;
    const drawing = state.activeDrawing();
    if (!outcome || outcome.drawingId !== drawing.id) return null;
    return outcome.setId === null || drawing.annotationSets.some((set) => set.id === outcome.setId) ? outcome : null;
  });
  const hitRows = computed(() => (shownOutcome.value?.hits ?? []).filter((hit) => !panel.chip || hit.query === panel.chip));
  const chips = computed(() => (shownOutcome.value?.queries ?? []).map(({ name }) => ({
    name, count: (shownOutcome.value?.hits ?? []).filter((hit) => hit.query === name).length
  })));
  // "No hits for <query> (N bp)", and how many alignments the thresholds removed.
  const noHitText = computed(() => {
    const outcome = shownOutcome.value;
    if (!outcome || outcome.setId !== null) return '';
    const names = outcome.queries.map(({ name, length }) => t('No hits for {query} ({length} bp)', { query: name, length }));
    return [...names, tn(outcome.rowsLeftOut, 'The thresholds removed {count} alignment.', 'The thresholds removed {count} alignments.', { count: outcome.rowsLeftOut })].join(' ');
  });
  // The hits the drawing outlines, by the shown set's id and the hits' annotation ids (app/search-hit-outline.js).
  const outline = computed(() => {
    const outcome = shownOutcome.value;
    const annotationIds = outcome?.setId ? outlinedHitIds(outcome.hits, panel.chip, panel.selectedHit) : [];
    return outcome?.setId && annotationIds.length ? { setId: outcome.setId, annotationIds } : null;
  });
  // The folded line names the program once: "blastn" for the blastn task, "blastn, task megablast" for another task.
  const programText = computed(() => (panel.task === 'blastn' ? 'blastn' : t('{program}, task {task}', { program: 'blastn', task: panel.task })));
  /** The record a hit is on, by its ID (the request key when Python read no ID). @param {SearchHit} hit */
  const recordOf = (hit) => shownOutcome.value?.recordIds[hit.recordKey] || hit.recordKey;
  /** @param {string} name */
  const chooseChip = (name) => { panel.chip = name; panel.selectedHit = ''; };
  /** @param {SearchHit} hit */
  const selectHit = (hit) => { panel.selectedHit = panel.selectedHit === hit.annotationId ? '' : hit.annotationId; };
  const saveRawTsv = () => {
    const outcome = shownOutcome.value;
    if (outcome) downloadTextFile(`${(outcome.setId ?? 'sequence_search').replace(/\s+/g, '_')}.losatn.tsv`, outcome.rawText);
  };

  // Reactive, so the template reads the computed values unwrapped (`sequenceSearch.hitRows`).
  return reactive({
    panel, readQueries, editQueryText, chooseQueryFile, addQuery, removeQuery, setQueryColor, searchDisabledReason, search,
    searchSetIdsOf, shownOutcome, hitRows, chips, noHitText, outline, programText, recordOf, chooseChip, selectHit, saveRawTsv,
    /** The hits to outline, for the outline's port (app/search-hit-outline.js). */
    outlineTarget: () => outline.value,
    tasks: SEARCH_TASKS,
    // The panel's bar of the gradient by query position, the colours of the Legend's "Query position" bar.
    gradientBar: `linear-gradient(to right, ${SEARCH_GRADIENT_START_COLOR}, ${SEARCH_GRADIENT_END_COLOR})`,
    /** @param {boolean} folded */
    setFolded: (folded) => { panel.folded = folded; }
  });
};
