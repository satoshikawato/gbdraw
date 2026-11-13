// @ts-check
// Find in other files (PD-OI-123): what the dialog sends to Python and reads back. Python reads the Files of the
// shown drawing's request and owns every rule (the reference genes, the search plan, the windows, the drawing); this
// module types the Worker payloads, checks each reply against the description Python generates
// (`contracts/locate.generated.js`), turns the search plan into LOSAT jobs, and reads a reference row the user typed.
import { LOCATE_DESCRIPTION } from '../contracts/locate.generated.js';
import { shapeViolations } from '../utils/contract-shape.js';
import { diagnosticError } from '../utils/error-normalization.js';
import { numberText } from '../utils/i18n.js';
import { windowSpanOf } from './linear-card-windows.js';

/** @import { LosatJob, LosatJobResult } from './losat.js' */
/** @import { CanonicalRecord, CanonicalRenderRequest } from '../contracts/request.generated.js' */
/** @import { GenerationPayload } from './diagram-generation.js' */
/**
 * @import { LocateDrawingInput, LocateDrawingReady, LocateEditRefused, LocateFilesReply, LocateOutcomeInput,
 *   LocatePlanInput, LocatePlanReply, LocateReferenceInput, LocateReferenceReady, LocateSearchInput,
 *   LocateSearchPlanReply, LocateWindowInput, LocateWindowRefused, WindowGeneJson } from '../contracts/locate.generated.js'
 */

/** The Worker operation of every step. */
export const LOCATE_OPERATION = 'findInOtherFiles';

/**
 * One step of the Worker operation (`LocateStep` in gbdraw/web_support/locate.py).
 * @typedef {{ step: 'files' } | { step: 'reference', reference: LocateReferenceInput }
 *   | { step: 'search', search: LocateSearchInput } | { step: 'plan', plan: LocatePlanInput }
 *   | { step: 'drawing', drawing: LocateDrawingInput }} LocateStep
 */

/**
 * The Worker payload of a step: the shown drawing's draft request as the next Generate would send it, whose resources
 * the page stages, and the step.
 * @typedef {object} LocateHelperPayload
 * @property {CanonicalRenderRequest} canonicalRequest
 * @property {GenerationPayload['resources']} resources
 * @property {LocateStep} locate
 */

/** @typedef {LocateReferenceReady | LocateWindowRefused} LocateReferenceReply */
/** @typedef {LocateDrawingReady | LocateEditRefused} LocateDrawingReply */

const replyInvalid = () => diagnosticError('RESULT_INVALID', {}, { stage: 'helper', operation: LOCATE_OPERATION });

/** @param {unknown} value @returns {value is Record<string, unknown>} */
const isRecord = (value) => Boolean(value) && typeof value === 'object' && !Array.isArray(value);

/**
 * A reply checked against the root `root` of Python's description; a reply of another shape is refused.
 * @template T
 * @param {string} root
 * @param {unknown} value
 * @returns {T}
 */
const checked = (root, value) => {
  const name = LOCATE_DESCRIPTION.roots[root];
  if (!name || shapeViolations(value, { kind: 'ref', name }, LOCATE_DESCRIPTION.types, 'reply').length) throw replyInvalid();
  return /** @type {T} */ (value);
};

/** @param {unknown} value @returns {LocateFilesReply} */
export const readFilesReply = (value) => checked('LocateFilesReply', value);
/** @param {unknown} value @returns {LocateReferenceReply} */
export const readReferenceReply = (value) => (isRecord(value) && value.status === 'refused'
  ? checked('LocateWindowRefused', value) : checked('LocateReferenceReady', value));
/** @param {unknown} value @returns {LocateSearchPlanReply} */
export const readSearchPlanReply = (value) => checked('LocateSearchPlanReply', value);
/** @param {unknown} value @returns {LocatePlanReply} */
export const readPlanReply = (value) => checked('LocatePlanReply', value);
/** @param {unknown} value @returns {LocateDrawingReply} */
export const readDrawingReply = (value) => (isRecord(value) && value.status === 'refused'
  ? checked('LocateEditRefused', value) : checked('LocateDrawingReady', value));

/**
 * The LOSAT jobs of a search plan: one per File, in the plan's order, each searching the reference genes (one query
 * sequence for all) in every record of its File with the plan's program and options.
 * @param {LocateSearchPlanReply} plan
 * @returns {{ jobs: (LosatJob & { fileKey: string })[], sequences: Map<string, string> }}
 */
export const locateLosatJobs = (plan) => {
  const sequences = new Map();
  const jobs = plan.jobs.map((job, pairIndex) => {
    const querySequenceKey = 'locate-query';
    const subjectSequenceKey = `locate-subject-${pairIndex}`;
    sequences.set(querySequenceKey, job.queryFasta);
    sequences.set(subjectSequenceKey, job.subjectFasta);
    return {
      pairIndex, cacheKey: job.rawKey, program: plan.search, querySequenceKey, subjectSequenceKey,
      extraArgs: [...plan.cacheArgs], fileKey: job.fileKey
    };
  });
  return { jobs, sequences };
};

/**
 * Each job's rows, or why it failed, by its File: a result is matched to its job by `pairIndex`.
 * @param {readonly (LosatJob & { fileKey: string })[]} jobs
 * @param {readonly (LosatJobResult & { error?: string })[]} results
 * @returns {LocateOutcomeInput[]}
 */
export const locateOutcomes = (jobs, results) => jobs.map((job) => {
  const result = results.find((item) => item.pairIndex === job.pairIndex);
  if (!result) throw replyInvalid();
  return { fileKey: job.fileKey, text: result.error ? '' : String(result.text ?? ''), error: result.error ?? null };
});

/**
 * One reference row as the user typed it: Start and End, and "Across the origin" for a circular record.
 * @typedef {object} ReferenceRowDraft
 * @property {string} start
 * @property {string} end
 * @property {boolean} acrossOrigin
 */

/** A whole number the user typed (thousands separators allowed), or null. @param {string} text @returns {number | null} */
export const wholeNumber = (text) => (/^\s*[1-9][\d,]*\s*$/.test(text) ? Number(text.replace(/[\s,]/g, '')) : null);

/**
 * The window a reference row gives, read as the record card reads its Region fields (`windowSpanOf`, PD-OI-124): with
 * "Across the origin" a Start greater than the End crosses the origin; without it, an End before the Start reverses
 * the span.
 * Null with the reason when the row is not a window yet.
 * @param {ReferenceRowDraft} row
 * @returns {{ window: LocateWindowInput } | { reason: 'EMPTY' | 'NOT_WHOLE' | 'ORIGIN_NEEDS_START_AFTER_END' }}
 */
export const referenceWindowOf = (row) => {
  if (row.start.trim() === '' || row.end.trim() === '') return { reason: 'EMPTY' };
  const start = wholeNumber(row.start);
  const end = wholeNumber(row.end);
  if (start === null || end === null) return { reason: 'NOT_WHOLE' };
  const window = windowSpanOf(start, end, row.acrossOrigin);
  return window ? { window } : { reason: 'ORIGIN_NEEDS_START_AFTER_END' };
};

/**
 * A window's span as the dialog shows it: thousands separators, and "(−)" when it is reverse-complemented.
 * @param {{ start: number, end: number, reverse_complement: boolean }} window
 */
export const windowSpanText = (window) => `${numberText(window.start)}–${numberText(window.end)}${
  window.reverse_complement ? ' (−)' : ''}`;

/** A window's number in the drawing as a circled number (①..⑳, then "(21)"). @param {number} number */
export const circledNumber = (number) => (number >= 1 && number <= 20 ? String.fromCodePoint(0x2460 + number - 1) : `(${number})`);

/**
 * The first and last gene a window holds in the order of the reference, by the names Python gives them (as the
 * report of `gbdraw locate` names them).
 * @param {readonly WindowGeneJson[]} genes
 * @param {Readonly<Record<string, string>>} names
 * @returns {{ count: number, first: string, last: string }}
 */
export const windowGeneSpan = (genes, names) => {
  const ordered = [...genes].sort((a, b) => a.order - b.order);
  const nameOf = (/** @type {WindowGeneJson | undefined} */ gene) => (gene ? names[gene.protein_id] ?? gene.protein_id : '');
  return { count: genes.length, first: nameOf(ordered[0]), last: nameOf(ordered[ordered.length - 1]) };
};

/**
 * The name of the Web file each record of a drawing's request reads (the GenBank file, else the GFF3 file, else the
 * FASTA file), by its record key. A File's key is the record key of the first request record of its source
 * (`d:1` for a card the request writes as one record per record), so every File key is named.
 * @param {readonly CanonicalRecord[]} records The request's records.
 * @param {ReadonlyMap<string, { name: string }>} filesByResource The Web file of each resource.
 * @returns {Record<string, string>}
 */
export const fileNamesOf = (records, filesByResource) => Object.fromEntries(records.map(({ recordKey, source }) => {
  const resourceId = source.kind === 'genbank' ? source.resourceId
    : source.kind === 'gffFasta' ? source.gffResourceId : source.fastaResourceId;
  return [recordKey, filesByResource.get(resourceId)?.name ?? recordKey];
}));
