// @ts-check
// A drawing's annotation record catalog: the records its annotations can name
// (a Circular drawing's record set, or a Linear drawing's drawn cards), why a
// catalog is not ready, and the catalog Generate and Add to annotation set read
// once every Linear card's records are read (`prepareLinearRecordCatalog`).
// The composition root passes the state, its reader of one input's records,
// and the Linear record selector's members as ports. It keeps the rebind of the
// shown drawing's annotations to the catalog (`rebindAnnotationRecords`, R10).
import { circularDiscoveryForInput } from '../services/record-discovery.js';
import { circularInputSlots, inputSourceFiles } from '../services/input-types.js';
import { circularDrawingRecordSet } from '../services/record-options.js';
import { diagnosticError, normalizeUserFacingError } from '../utils/error-normalization.js';
import { annotationSourceKey, buildAnnotationRecordCatalog } from './annotations/record-catalog.js';
/** @import { DrawingState, LinearSeq, SessionBusy, WebState } from '../state.js' */
/** @import { InputFileRef } from '../services/history-files.js' */
/** @import { SequenceRecord } from '../services/record-discovery.js' */
/** @import { UserFacingError } from '../utils/error-normalization.js' */
/** @import { AnnotationCatalogSource, AnnotationRecordCatalog } from './annotations/record-catalog.js' */
/** @import { RecordSelectorStatus, RefreshOutcome } from './linear-record-selector.js' */

/**
 * What the record selector holds for one Linear card: the status, records and
 * error model its record select shows.
 * @typedef {{ status: RecordSelectorStatus, records: SequenceRecord[], error: UserFacingError | '' }} LinearCardCatalogState
 */

/**
 * @typedef {object} DrawingRecordCatalogPorts
 * @property {WebState} state
 * @property {(inputType: string, primaryFile: InputFileRef | null, pairedFile: InputFileRef | null) => Promise<SequenceRecord[]>} readInputRecords
 *   The records of one input, read from the file store.
 * @property {(drawing: DrawingState, seq: LinearSeq) => LinearCardCatalogState} linearRecordCatalogStateFor
 * @property {(drawing: DrawingState, seq: LinearSeq) => UserFacingError | ''} linearRecordErrorModelFor
 *   The read error a Linear card's record select stores, or '' when it has none.
 * @property {(drawing: DrawingState) => SessionBusy | Promise<RefreshOutcome>} refreshLinearRecords
 *   Reads the records of the drawing's Linear cards.
 */

/**
 * Registers nothing. The root creates it once the record selector exists,
 * before Generate's runner, which takes `prepareLinearRecordCatalog`.
 * @param {DrawingRecordCatalogPorts} ports
 */
export const createDrawingRecordCatalog = ({
  state, readInputRecords, linearRecordCatalogStateFor, linearRecordErrorModelFor, refreshLinearRecords
}) => {
  /** @param {DrawingState} drawing */
  const getCircularRecordDiscoveryState = (drawing) => circularDiscoveryForInput(drawing);
  /** @param {DrawingState} drawing */
  const circularAnnotationSourceKey = (drawing) => annotationSourceKey({
    scope: 'circular',
    inputType: drawing.inputType.value,
    ...inputSourceFiles(drawing.inputType.value, circularInputSlots(drawing.files))
  });
  /**
   * A Linear record card as the annotation record catalog reads it.
   * @param {DrawingState} drawing the drawing that holds the card
   * @param {LinearSeq} seq
   * @param {boolean} drawn
   * @returns {AnnotationCatalogSource}
   */
  const linearAnnotationSource = (drawing, seq, drawn) => {
    const inputType = drawing.inputType.value;
    const { primaryFile, pairedFile, complete } = inputSourceFiles(inputType, seq);
    return {
      sourceKey: annotationSourceKey({
        scope: 'linear',
        uid: seq.uid,
        inputType,
        primaryFile,
        pairedFile
      }),
      selector: seq.region_record_id,
      drawn,
      hasInput: complete,
      ...linearRecordCatalogStateFor(drawing, seq)
    };
  };
  /**
   * The records one drawing's annotations can name: a Circular drawing's
   * record set, or a Linear drawing's drawn cards (`linearSourcesOverride`
   * replaces the cards a Linear drawing reads).
   * @param {DrawingState} drawing
   * @param {AnnotationCatalogSource[] | null} [linearSourcesOverride]
   */
  const getAnnotationRecordCatalog = (drawing, linearSourcesOverride = null) => {
    const circularSource = () => {
      const circularDiscovery = getCircularRecordDiscoveryState(drawing);
      const circularRecordSet = circularDrawingRecordSet(circularDiscovery.records, drawing);
      return {
        sourceKey: circularAnnotationSourceKey(drawing),
        hasInput: circularDiscovery.hasInput,
        status: circularDiscovery.status,
        error: circularDiscovery.error,
        records: circularRecordSet.records,
        omittedRecords: circularRecordSet.omittedRecords
      };
    };
    const linearSources = () => {
      const recordsOff = new Set(drawing.recordsOff);
      return linearSourcesOverride || drawing.linearSeqs.map((seq) => linearAnnotationSource(drawing, seq, !recordsOff.has(seq.uid)));
    };
    return buildAnnotationRecordCatalog(drawing.mode === 'linear'
      ? { mode: 'linear', linearSources: linearSources() }
      : { mode: 'circular', circularSource: circularSource() });
  };

  /** @param {Pick<AnnotationRecordCatalog, 'issues'>} catalog */
  const catalogIssueError = (catalog) => {
    const issue = catalog.issues[0];
    return issue ? diagnosticError(issue.code, issue.context) : diagnosticError('INPUT_UNREADABLE');
  };
  /**
   * Why a drawing's record catalog is not ready, as its record list says it: for Linear, the first drawn card's own
   * read error (its record select's), else the catalog's first issue. Generate and Add to annotation set read it here.
   * @param {DrawingState} drawing
   * @param {Pick<AnnotationRecordCatalog, 'issues'>} catalog
   */
  const recordCatalogError = (drawing, catalog) => (drawing.mode === 'linear'
    // A stored '' (no error) has no `code`: the read is undefined, as for null.
    && drawing.drawnLinearSeqs.value.map(seq => linearRecordErrorModelFor(drawing, seq)).find(error => /** @type {UserFacingError} */ (error)?.code))
    || catalogIssueError(catalog);
  async function prepareLinearRecordCatalog({ privateCandidate = false } = {}) {
    const drawing = state.activeDrawing();
    if (drawing.mode !== 'linear') return { catalog: null, error: '' };
    const hasAutomaticSequence = drawing.linearSeqs.some((sequence) => {
      if (String(sequence?.region_record_id || '').trim()) return false;
      return inputSourceFiles(drawing.inputType.value, sequence).complete;
    });
    const hasRegionAnnotations = drawing.annotationSets.some((/** @type {{ annotations?: unknown } | null | undefined} */ set) => (
      Array.isArray(set?.annotations) && set.annotations.length > 0
    ));
    if (!hasAutomaticSequence && !hasRegionAnnotations) {
      return { catalog: null, error: '' };
    }
    if (privateCandidate) {
      const linearRecordsOff = new Set(drawing.recordsOff);
      const sources = await Promise.all(drawing.linearSeqs.map(async (seq) => {
        const drawn = !linearRecordsOff.has(seq.uid);
        const inputType = drawing.inputType.value;
        const { primaryFile, pairedFile } = inputSourceFiles(inputType, seq);
        const sourceKey = annotationSourceKey({ scope: 'linear', uid: seq.uid, inputType, primaryFile, pairedFile });
        try {
          const records = await readInputRecords(inputType, primaryFile, pairedFile);
          return { sourceKey, selector: seq.region_record_id, drawn, hasInput: Boolean(primaryFile), status: 'ready', records };
        } catch (error) {
          /** @type {UserFacingError | ''} The model the card's record select would store for this failure. */
          const model = normalizeUserFacingError(error, { operation: 'listSequenceRecords', stage: 'helper' }) || '';
          return { sourceKey, selector: seq.region_record_id, drawn, hasInput: Boolean(primaryFile), status: 'error', error: model, records: [] };
        }
      }));
      const catalog = getAnnotationRecordCatalog(drawing, sources);
      return catalog.status === 'ready' ? { catalog, error: '' }
        : { catalog: null, error: catalogIssueError(catalog) };
    }
    let catalog = getAnnotationRecordCatalog(drawing);
    if (catalog.status !== 'ready') {
      try {
        await refreshLinearRecords(drawing);
      } catch (error) {
        return { catalog: null, error: normalizeUserFacingError(error, { operation: 'listSequenceRecords', stage: 'helper' }) };
      }
      catalog = getAnnotationRecordCatalog(drawing);
    }
    return catalog.status === 'ready' ? { catalog, error: '' } : { catalog: null, error: recordCatalogError(drawing, catalog) };
  }

  return {
    getCircularRecordDiscoveryState, circularAnnotationSourceKey, linearAnnotationSource, getAnnotationRecordCatalog,
    recordCatalogError, prepareLinearRecordCatalog
  };
};
