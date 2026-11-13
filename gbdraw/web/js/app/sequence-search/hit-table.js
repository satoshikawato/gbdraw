// @ts-check
// The hit table of the Sequence search panel (PD-OI-117): the rows of the latest search of the shown drawing while
// its set exists, the chips that filter them by query, the hits a chip or a row click outlines on the drawing
// (`hit-outline.js`; no History step, not saved, not exported, and the figure keeps every hit), the no-hit message, and
// Save Raw LOSAT TSV.
import { computed } from '../../utils/reactivity.js';
import { numberText, t, tn } from '../../utils/i18n.js';
import { downloadTextFile } from '../../services/text-download.js';
import { outlinedHitIds, searchPrograms } from '../../services/sequence-search.js';
import { createSearchHitOutline } from './hit-outline.js';

/** @import { SearchHit, SearchProgram } from '../../services/sequence-search.js' */
/** @import { DrawingId } from '../../services/drawing-view.js' */

/**
 * What the latest search found: the set it added (null when no hit passed the thresholds), its hits and queries.
 * @typedef {object} SearchOutcome
 * @property {DrawingId} drawingId
 * @property {string | null} setId
 * @property {SearchHit[]} hits
 * @property {{ name: string, length: number, program: SearchProgram, color: string }[]} queries
 * @property {number} rowsLeftOut Rows LOSAT wrote that the thresholds removed.
 * @property {Record<string, string>} recordIds The record ID of each searched record, by its request key (as Python read
 *   the spans), which the hit table shows.
 * @property {string} rawText The rows as LOSAT wrote them (Save Raw LOSAT TSV).
 */

/**
 * The fields of the panel's state that the table reads; it writes `chip` and `selectedHit`.
 * @typedef {object} SearchTableState
 * @property {SearchOutcome | null} outcome The latest search's outcome (the search writes it).
 * @property {string} chip The query whose hits the table shows ('' for all).
 * @property {string} selectedHit The annotation id of the hit a row click outlines ('' for none).
 */

/**
 * @typedef {object} SearchTablePorts
 * @property {() => { id: DrawingId, annotationSets: readonly { id: string }[] }} activeDrawing The shown drawing.
 * @property {() => SVGSVGElement | null} shownRoot The mounted root of the shown Result, which the outline follows.
 */

/** @param {SearchTableState} panel @param {SearchTablePorts} ports */
export const createSearchHitTable = (panel, { activeDrawing, shownRoot }) => {
  // The latest search of the shown drawing, while its set exists (a search without hits has no set).
  const shownOutcome = computed(() => {
    const outcome = panel.outcome;
    const drawing = activeDrawing();
    if (!outcome || outcome.drawingId !== drawing.id) return null;
    return outcome.setId === null || drawing.annotationSets.some((set) => set.id === outcome.setId) ? outcome : null;
  });
  const hitRows = computed(() => (shownOutcome.value?.hits ?? []).filter((hit) => !panel.chip || hit.query === panel.chip));
  // The aligned length of a hit as its row shows it: bases for a nucleotide query, residues (aa) for a protein query.
  const programs = computed(() => new Map((shownOutcome.value?.queries ?? []).map(({ name, program }) => [name, program])));
  /** @param {SearchHit} hit */
  const hitLengthText = (hit) => (programs.value.get(hit.query) === 'tblastn'
    ? t('{length} aa', { length: String(hit.alignmentLength) })
    : String(hit.alignmentLength));
  const chips = computed(() => (shownOutcome.value?.queries ?? []).map(({ name }) => ({
    name, count: (shownOutcome.value?.hits ?? []).filter((hit) => hit.query === name).length
  })));
  // "No hits for <query> (N bp)" ("N aa" for a protein query), and how many alignments the thresholds removed.
  const noHitText = computed(() => {
    const outcome = shownOutcome.value;
    if (!outcome || outcome.setId !== null) return '';
    const names = outcome.queries.map(({ name, length, program }) => (program === 'tblastn'
      ? t('No hits for {query} ({length} aa)', { query: name, length: numberText(length) })
      : t('No hits for {query} ({length} bp)', { query: name, length: numberText(length) })));
    return [...names, tn(outcome.rowsLeftOut, 'The thresholds removed {count} alignment.', 'The thresholds removed {count} alignments.', { count: numberText(outcome.rowsLeftOut) })].join(' ');
  });
  // The hits the drawing outlines, by the shown set's id and the hits' annotation ids.
  const outline = computed(() => {
    const outcome = shownOutcome.value;
    const annotationIds = outcome?.setId ? outlinedHitIds(outcome.hits, panel.chip, panel.selectedHit) : [];
    return outcome?.setId && annotationIds.length ? { setId: outcome.setId, annotationIds } : null;
  });
  /** The record a hit is on, by its ID (the request key when Python read no ID). @param {SearchHit} hit */
  const recordOf = (hit) => shownOutcome.value?.recordIds[hit.recordKey] || hit.recordKey;
  /** @param {string} name */
  const chooseChip = (name) => { panel.chip = name; panel.selectedHit = ''; };
  /** @param {SearchHit} hit */
  const selectHit = (hit) => { panel.selectedHit = panel.selectedHit === hit.annotationId ? '' : hit.annotationId; };
  // Named by its programs as Pairwise Comparisons names its tables: `.losatn.tsv` (blastn), `.tlosatn.tsv` (tblastn).
  const saveRawTsv = () => {
    const outcome = shownOutcome.value;
    if (!outcome) return;
    const programs = searchPrograms(outcome.queries).map((program) => (program === 'tblastn' ? 'tlosatn' : 'losatn')).join('-');
    downloadTextFile(`${(outcome.setId ?? 'sequence_search').replace(/\s+/g, '_')}.${programs}.tsv`, outcome.rawText);
  };
  const hitOutline = createSearchHitOutline({ shownRoot, target: () => outline.value });
  return { shownOutcome, hitRows, hitLengthText, chips, noHitText, outline, recordOf, chooseChip, selectHit, saveRawTsv, hitOutline };
};
