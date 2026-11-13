// @ts-check
// The Sequence search section (PD-OI-117): the panel's queries and settings,
// the search, and the hit table of the latest search. The queries and settings
// are the panel's own, transient values: not in the Session, no History step,
// and never written to Pairwise Comparisons. A search runs LOSAT through the
// Pairwise Comparisons runtime (`runLosatJobs`), Python turns its rows into
// an annotation set (`services/sequence-search.js`), and adding that set is one
// History step of the shown drawing that the automatic rerender draws at once
// (`searchSetIdsOfStep`, the live projection in services/session-request.js).
// The queries (`sequence-search/queries.js`) and the hit table with its outline
// (`sequence-search/hit-table.js`) are owners of their own.
import { computed, reactive } from '../utils/reactivity.js';
import { t, tm } from '../utils/i18n.js';
import { normalizeUserFacingError } from '../utils/error-normalization.js';
import { comparisonStateForMode, resolveComparisonThresholds } from '../mode-profiles.js';
import { draftAnnotationSetsOfRequest } from '../services/annotation-state.js';
import { boundAnnotation } from './annotations/record-selector.js';
import { SEARCH_GRADIENT_END_COLOR, SEARCH_GRADIENT_START_COLOR } from '../utils/sequence-search-colors.generated.js';
import {
  nextSearchSetId, readSequenceSearchReply, resultRecordBinding, searchLosatPlan, searchPrograms, searchRawText, searchSetIdsOfStep,
  searchSpans, searchSubjects
} from '../services/sequence-search.js';
import { createSearchQueries } from './sequence-search/queries.js';
import { createSearchHitTable } from './sequence-search/hit-table.js';

/** @import { ResultRecordBindings, SearchColorMode, SearchHit, SearchQuery, SequenceSearchInput } from '../services/sequence-search.js' */
/** @import { HistoryIntentChange } from '../services/history.js' */
/** @import { RegionAxisRecord } from '../services/region-axis.js' */
/** @import { LosatJob, LosatJobResult } from '../services/losat.js' */
/** @import { RecordSpanRead, RecordSpanSequencesReply } from '../services/diagram-worker-protocol.js' */
/** @import { AnnotationSet } from '../services/annotation-state.js' */
/** @import { DrawingId } from '../services/drawing-view.js' */
/** @import { RunHelper, SearchError } from './sequence-search/queries.js' */
/** @import { SearchOutcome } from './sequence-search/hit-table.js' */

/**
 * The thresholds the panel's fields hold: a number field holds its number, or null when it is blank (`v-blank-null`),
 * which takes the mode default (`resolveComparisonThresholds`); the E-value field holds its text.
 * @typedef {{ min_bitscore: number | null, evalue: string, identity: number | null, alignment_length: number | null }} SearchThresholdDraft
 */

/**
 * The ports the composition root gives the panel.
 * @typedef {object} SequenceSearchPorts
 * @property {{ activeDrawing: () => { id: DrawingId, mode: 'circular' | 'linear', annotationSets: readonly { id: string }[] } }} state
 * @property {() => RegionAxisRecord[]} readAxisRecords Every drawn record of the shown drawing's Results.
 * @property {(spans: readonly RecordSpanRead[]) => Promise<RecordSpanSequencesReply>} readRecordSpans
 * @property {() => readonly { id: string }[]} committedSets The annotation sets of the committed request the shown
 *   Results were made from.
 * @property {() => ResultRecordBindings} readRecordBindings How a row binds each record of the shown Result
 *   (`resultRecordBindings`): the binding key Add to annotation set writes, and its selector in the committed request.
 * @property {(plan: { jobs: LosatJob[], sequences: Map<string, string> }) => Promise<{ status: 'ok', results: LosatJobResult[] } | { status: 'canceled' } | { status: 'busy', reason: string }>} runLosatJobs
 *   LOSAT jobs under the full-window overlay with its progress and Cancel (`runLosatJobsInOverlay`).
 * @property {RunHelper} runHelper
 * @property {{ readQueries: string, buildHitSet: string }} helperOperations The Worker helper operations (`readSearchQueries`,
 *   `searchHitAnnotationSet` in services/diagram-worker-protocol.js) the panel runs.
 * @property {(drawingId: DrawingId, set: AnnotationSet) => unknown} addSearchSet Adds the set to that drawing (the
 *   annotation-sets slice writer) as one History step of the shown drawing.
 * @property {() => SVGSVGElement | null} shownRoot The mounted root of the shown Result, which the outline follows.
 * @property {() => void} closeRegionPopup
 * @property {() => Promise<unknown>} revealSection Opens the panel's section.
 */

export const SEARCH_TASKS = Object.freeze(['megablast', 'blastn', 'dc-megablast']);

/**
 * The set of a search with each hit bound to its Result record as the record select binds a row: by the binding key of
 * the record's settings owner (an OFF record's too, so its rows wait with it at Generate), its target naming the record
 * in the committed request, which the live draw resolves again by that key (services/session-request.js).
 * @param {AnnotationSet} set
 * @param {readonly SearchHit[]} hits
 * @param {ResultRecordBindings} bindings
 * @returns {AnnotationSet}
 */
const boundSearchSet = (set, hits, bindings) => {
  const recordKeys = new Map(hits.map((hit) => [hit.annotationId, hit.recordKey]));
  return { ...set, annotations: set.annotations.map((row) => {
    return boundAnnotation(row, resultRecordBinding(bindings, recordKeys.get(row.id) ?? ''));
  }) };
};

// The panel as it opens, and after a Session Load.
const initialPanel = () => ({
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

/** @param {SequenceSearchPorts} ports */
export const createSequenceSearch = ({
  state, readAxisRecords, readRecordSpans, committedSets, readRecordBindings, runLosatJobs, runHelper, helperOperations,
  addSearchSet, shownRoot, closeRegionPopup, revealSection
}) => {
  const panel = reactive(initialPanel());
  // A Session Load counts here from its start (`loadStarts`), before it installs the loaded document: a search it
  // overtook adds nothing to that document and writes nothing in the panel it empties.
  let loads = 0;
  // The sets each drawing's searches added, by id, with the set as written: a History step that adds one of them
  // (its own, its Redo) or whose Undo removes it draws at once.
  /** @type {Map<DrawingId, Map<string, string>>} */
  const searched = new Map();
  const queries = createSearchQueries(panel, runHelper, helperOperations.readQueries);
  const table = createSearchHitTable(panel, { activeDrawing: state.activeDrawing, shownRoot });

  // Why Search is off, in the order the user meets them; '' when it can run.
  const searchDisabledReason = computed(() => {
    // Only a Result names the drawn bases (its record axis, as for Select region); the axis is reactive, the committed
    // Session that every Result comes from is not, so the axis alone decides.
    if (readAxisRecords().length === 0) return t('Generate to search');
    // A text not read yet, or whose read still runs, is read when Search runs: Search waits for it (OV-2579).
    const unread = queries.unread();
    if (panel.queryError && !unread) return t('Correct the query to search');
    if (panel.queries.length === 0 && (!unread || panel.queryText.trim() === '')) return t('Paste or choose a query FASTA');
    return '';
  });

  const search = async () => {
    if (searchDisabledReason.value || panel.running) return;
    panel.searchError = null;
    panel.busy = '';
    panel.running = true;
    const load = loads;
    try {
      await queries.readQueries(panel.queryText);
      // The read may refuse the text or find no query; its reason then shows.
      if (searchDisabledReason.value || load !== loads) return;
      const drawing = state.activeDrawing();
      // The records of the shown Result, as last generated (receipt): never the draft's (review r2 B1).
      const bindings = readRecordBindings();
      const listed = panel.queries;
      const searchQueries = listed.map(({ name, length, program, color }) => ({ name, length, program, color }));
      const thresholds = resolveComparisonThresholds(panel.thresholds[drawing.mode], drawing.mode);
      const records = readAxisRecords();
      const read = await readRecordSpans(searchSpans(records, listed));
      // A nucleotide query is searched with blastn, a protein query with tblastn: one job per program (PD-OI-117).
      const run = await runLosatJobs(searchLosatPlan({ queries: listed, subjectFasta: read.fasta, task: panel.task }));
      // A cancelled search, or one a Generate or a Session load kept from running, adds nothing and keeps the Result.
      if (load !== loads) return;
      if (run.status === 'busy') panel.busy = tm(run.reason);
      if (run.status !== 'ok' || state.activeDrawing().id !== drawing.id) return;
      const text = searchRawText(run.results);
      /** @type {SequenceSearchInput} */
      const input = {
        table: text,
        // Unique against the draft, the committed Session, and every set a search added (B2): a freed id that the
        // committed Session still draws, or that History can bring back, is not reused.
        setId: nextSearchSetId([...drawing.annotationSets, ...committedSets()].map(({ id }) => id)
          .concat([...(searched.get(drawing.id)?.keys() ?? [])])),
        subjects: searchSubjects(read, bindings),
        thresholds: {
          evalue: Number(thresholds.evalue), bitscore: Number(thresholds.bitscore),
          identity: Number(thresholds.identity), alignmentLength: Number(thresholds.alignmentLength)
        },
        queries: searchQueries.map(({ name, length, program }) => ({ name, length, program })),
        colorMode: panel.colorMode,
        queryColors: Object.fromEntries(searchQueries.map(({ name, color }) => [name, color]))
      };
      const reply = readSequenceSearchReply((await runHelper(helperOperations.buildHitSet, { search: input })).result);
      if (state.activeDrawing().id !== drawing.id || load !== loads) return;
      const [found] = reply.annotationSet ? draftAnnotationSetsOfRequest([/** @type {Partial<AnnotationSet>} */ (reply.annotationSet)]) : [];
      const set = found && boundSearchSet(found, reply.hits, bindings);
      panel.outcome = {
        drawingId: drawing.id, setId: set ? set.id : null, hits: reply.hits, queries: searchQueries,
        rowsLeftOut: reply.rowsRead - reply.rowsKept, rawText: text,
        recordIds: Object.fromEntries(read.spans.map((span) => [span.recordKey, span.recordId]))
      };
      panel.chip = '';
      panel.selectedHit = '';
      panel.folded = true;
      if (!set) return;
      const sets = searched.get(drawing.id) ?? new Map();
      searched.set(drawing.id, sets.set(set.id, JSON.stringify(set)));
      addSearchSet(drawing.id, set);
    } catch (error) {
      if (load === loads) panel.searchError = normalizeUserFacingError(error, { stage: 'losat', operation: helperOperations.buildHitSet });
    } finally {
      panel.running = false;
    }
  };

  /** Adds a query (`queries.addQuery`); the panel unfolds to show it. @param {{ name: string, fasta: string }} query */
  const addQuery = (query) => {
    panel.folded = false;
    return queries.addQuery(query);
  };
  /**
   * Use as search query (Region popup): Python reads the span on the source strand, named by its FASTA header
   * `<record>:<start>-<end>`; the panel adds it as a query, the popup closes, and the panel's section opens.
   * @param {RecordSpanRead} span
   */
  const useSpanAsQuery = async (span) => {
    const read = await readRecordSpans([span]);
    await addQuery({ name: read.spans[0].header, fasta: read.fasta });
    // The popup closes first, so the focus the section takes stays there.
    closeRegionPopup();
    await revealSection();
  };

  /**
   * The sets a History step of the shown drawing added that its searches wrote (`searchSetIdsOfStep`).
   * @param {readonly HistoryIntentChange[] | null | undefined} changes
   */
  const searchSetIdsOf = (changes) => searchSetIdsOfStep(changes, searched.get(state.activeDrawing().id) ?? new Map());

  // The programs of the listed queries (blastn before any query is listed). The Task applies to blastn only.
  const programs = computed(() => (panel.queries.length ? searchPrograms(panel.queries) : ['blastn']));
  const taskApplies = computed(() => programs.value.includes('blastn'));
  // The folded line names each program once: "blastn" for the blastn task, "blastn, task megablast" for another task.
  const programText = computed(() => programs.value
    .map((program) => (program !== 'blastn' || panel.task === 'blastn' ? program : t('{program}, task {task}', { program, task: panel.task })))
    .join(' · '));

  // A Session Load starts: a search still running is overtaken.
  const loadStarts = () => { loads += 1; };
  // The Session Load replaced the document: the panel empties (receipt), and its searches' sets are no search's any more.
  const reset = () => {
    loadStarts();
    Object.assign(panel, initialPanel(), { running: panel.running });
    searched.clear();
    queries.clearReads();
  };

  // Reactive, so the template reads the computed values unwrapped (`sequenceSearch.hitRows`).
  return reactive({
    ...queries, ...table, panel, addQuery, useSpanAsQuery, searchDisabledReason, search, searchSetIdsOf, programs, taskApplies, programText,
    loadStarts, reset,
    tasks: SEARCH_TASKS,
    // The panel's bar of the gradient by query position, the colours of the Legend's "Query position" bar.
    gradientBar: `linear-gradient(to right, ${SEARCH_GRADIENT_START_COLOR}, ${SEARCH_GRADIENT_END_COLOR})`,
    /** @param {boolean} folded */
    setFolded: (folded) => { panel.folded = folded; }
  });
};
