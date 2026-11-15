// @ts-check
/** @import { DrawingState, WebState } from '../state.js' */
/** @import { CommittedRecord, CommittedRecordRequest } from '../services/feature-identity.js' */
// Source-bound editable rotation intent. The request service owns serialization.
import { circularSourceInputs, inputSourceFiles, storedRecordSourceResourceIds } from '../services/input-types.js';
import {
  bindRecordDisplayDrafts,
  buildRecordDisplayRows,
  canonicalRecordReverseComplement,
  captureLinearCardOrientation,
  cropFormReverseComplement,
  linearCardCropForm,
  parseRecordDisplayStart,
  recordDisplayCrop,
  recordDisplayKey,
  recordDisplaySurface,
  requestedRecordTransform,
  requireReverseComplementOverride,
  validateAnchorIntent
} from '../services/record-display-model.js';
import { orientLinearCard, restoreLinearCardOrientation } from '../services/linear-card-windows.js';
import { matchesSessionResourceDescriptor } from '../services/session-resource-backing.js';
import { cloneJsonData } from '../services/json-clone.js';
import { featureIdentityKey, requestFeaturePlacements } from '../services/feature-placement.js';
import { committedRecordForKey, committedRecordUsesSource } from '../services/feature-identity.js';
import { circularDiscoveryForInput } from '../services/record-discovery.js';
import { circularDrawingAloneRecord, circularDrawingCroppedRecords } from '../services/record-options.js';
/** @import { linearDrawingWindowsOf } from '../services/session-request.js' */
/** @import { DrawingId } from '../services/drawing-view.js' */
import { RECORD_READ_ERROR_LABEL } from './linear-record-selector.js';
import { ownEntry } from '../utils/own-entry.js';
import { draftFeatureOverlapTolerance } from '../services/session-active-config-contract.js';
import { msg, t, tm } from '../utils/i18n.js';
import { diagnosticError } from '../utils/error-normalization.js';
import { errorText } from '../utils/error-text.js';
/** @import { LinearSeq } from '../state.js' */
/** @import { LinearSeqWindow } from '../services/drawing/inputs.js' */
/** @import { LinearCardOrientation } from '../services/record-display-model.js' */
/** @import { RecordSelectorStatus } from './linear-record-selector.js' */
/** @import { SequenceRecord } from '../services/record-discovery.js' */
/** @import { InputFileRef } from '../services/history-files.js' */
/** @import { RecordOrientationDecision } from '../services/diagram-worker-protocol.js' */
/** @import { RecordOrientationWrite } from './record-auto-orientation.js' */
/** @import { AnchorIntent, RecordDisplayCrop, RecordDisplayDraft, RecordDisplayDraftIdentity, RecordDisplayRow, RecordDisplayScope } from '../services/record-display-model.js' */
/** @import { RecordEntryInput } from '../services/record-options.js' */
/** @import { RecordPlacementStart } from '../services/session-request.js' */

// Session Load reads no record bytes (776a2f93), so a loaded Result can have a
// current source whose records are not read yet. That popup target failure is
// not stale: its kind lets an explicit record action read them.
export const RECORD_TARGET_NOT_DISCOVERED = 'not-discovered';
/** @param {string} kind @param {string} message */
const recordTargetError = (kind, message) => Object.assign(new Error(message), { kind });

/** @import { CircularListRecord } from '../services/record-discovery.js' */

/**
 * A record source of one drawing (`sourcesOf`): a File card and its records,
 * named by the drawing's ID.
 * @typedef {object} DrawingRecordSource
 * @property {RecordDisplayScope} scope
 * @property {DrawingId} drawingId
 * @property {string} sourceUid
 * @property {unknown} source
 * @property {unknown} paired
 * @property {RecordEntryInput[]} records
 * @property {unknown} selector
 * @property {(row: DisplayRowOfSource) => RowCrop} cropOf Whether a record row of the source draws a crop, so it takes
 *   no display start, and the crop it draws.
 * @property {boolean} reverse
 * @property {boolean} shown Whether the drawing shows its records' rows (`circularDisplaySources`).
 */

/**
 * Whether a record row draws a crop, and the crop as far as it can be read (`RecordDisplayCrop`): null while it draws
 * none, or while a crop Generate refuses or a window of no stated end cannot be read.
 * @typedef {{ cropped: boolean, crop: RecordDisplayCrop | null }} RowCrop
 */
/**
 * A record row a source's records build (`buildRecordDisplayRows`): a Circular record keeps its identity.
 * @typedef {RecordDisplayRow & { identity?: unknown }} DisplayRowOfSource
 */
/** @type {RowCrop} */
const UNCROPPED = Object.freeze({ cropped: false, crop: null });

/**
 * The record sources of a Circular drawing: its source cards, each with the
 * drawing's crop and orientation. The rows are the records the drawing's
 * record set offers (`circularDrawingAloneRecord`, OV-3340): while it offers
 * one record alone, that record's card shows it alone and the other cards
 * show none, as the request draws it alone (OV-NEW-1).
 * @param {DrawingState} drawing @returns {DrawingRecordSource[]}
 */
const circularDisplaySources = (drawing) => {
  /** @type {CircularListRecord[]} */
  const list = drawing.circularRecordList.value;
  const alone = circularDrawingAloneRecord(list, drawing);
  const ranged = circularDrawingCroppedRecords(list, drawing);
  /** @param {DisplayRowOfSource} row @returns {RowCrop} */
  const cropOf = (row) => {
    const identity = String(row.identity ?? '');
    if (alone !== null) {
      return alone.crop !== null && alone.record.value === identity
        ? { cropped: true, crop: recordDisplayCrop(alone.crop) } : UNCROPPED;
    }
    const window = ranged.get(identity);
    return window === undefined ? UNCROPPED : { cropped: true, crop: window === null ? null : recordDisplayCrop(window) };
  };
  return circularSourceInputs(drawing.inputType.value, drawing.circularSources).map(({ uid, primaryFile, pairedFile }) => ({
    scope: 'circular', drawingId: drawing.id, sourceUid: uid, source: primaryFile, paired: pairedFile,
    records: list.filter((record) => record.uid === uid)
      .map((record) => ({ ...record, recordId: record.record_id, recordLength: record.record_length })),
    // A record offered alone is named in its card's file by its `#k`.
    selector: alone !== null && alone.record.uid === uid ? alone.record.selector : '',
    shown: alone === null || alone.record.uid === uid,
    cropOf,
    reverse: cropFormReverseComplement({ windows: drawing.form.circular_windows, reverse: drawing.form.circular_reverse })
  }));
};

/**
 * The windows a Linear drawing's request writes for a card, with the record the card draws (`linearDrawingWindowsOf`
 * of services/session-request.js, a port: the request owns them).
 * @typedef {typeof linearDrawingWindowsOf<LinearSeq>} LinearWindowsOf
 */

/**
 * Whether a record row of a Linear File card draws a crop (its Start/End fields or its windows), and the first window
 * it draws, as the request writes the card's windows (`linearWindowsOf`, with the row's record); a crop Generate
 * refuses is shown by the card's own error and gives no crop here.
 * @param {LinearWindowsOf} linearWindowsOf @param {DrawingState} drawing @param {LinearSeq} seq
 * @param {DisplayRowOfSource} row
 * @returns {RowCrop}
 */
const linearCardRowCrop = (linearWindowsOf, drawing, seq, row) => {
  if (!(seq.region_start || seq.region_end || seq.windows.length)) return UNCROPPED;
  const record = { length: row.recordLength ?? null, id: String(row.recordId ?? ''), count: null, fileCount: null };
  try {
    const [first] = linearWindowsOf(drawing, () => record)(seq);
    return { cropped: true, crop: first === undefined ? null : recordDisplayCrop(first) };
  } catch {
    return { cropped: true, crop: null };
  }
};

/**
 * The record sources of a Linear drawing: its cards, each read from that drawing's inputs, with its crop and
 * orientation.
 * @param {DrawingState} drawing
 * @param {(drawing: DrawingState, sequence: LinearSeq) => SequenceRecord[]} linearRecordsFor
 * @param {LinearWindowsOf} linearWindowsOf
 * @returns {DrawingRecordSource[]}
 */
const linearDisplaySources = (drawing, linearRecordsFor, linearWindowsOf) => drawing.linearSeqs.map((seq) => {
  const { primaryFile, pairedFile } = inputSourceFiles(drawing.inputType.value, seq);
  return { scope: 'linear', drawingId: drawing.id, sourceUid: seq.uid, source: primaryFile, paired: pairedFile,
    records: linearRecordsFor(drawing, seq), selector: seq.region_record_id,
    cropOf: (row) => linearCardRowCrop(linearWindowsOf, drawing, seq, row),
    reverse: cropFormReverseComplement(linearCardCropForm(seq)), shown: true };
});

/**
 * A row of the committed request: the draft row with what the request carries for the record.
 * @typedef {DrawnRecordRow & {
 *   committedDisplay: { isCircular: boolean | null, startCoordinate: number | null },
 *   committedReverseComplement: boolean, committedCropped: boolean,
 *   canonicalRecordKey: string, canonicalSource: Record<string, unknown>, recordKey: string
 * }} CommittedDisplayRow
 */

/**
 * A record's orientation in an alignment.
 * @typedef {{ recordKey: string, reverseComplement: boolean }} AlignmentOrientation
 */

/**
 * A Linear File card's orientation fields as a checkpoint keeps them (`captureLinearCardOrientation`): the windows it lists,
 * copied, and its Reverse complement.
 * @typedef {LinearCardOrientation<LinearSeqWindow>} CardOrientation
 */

/**
 * What `captureAlignmentOrientationIntent` keeps of an orientation: the File card's orientation fields and the record's
 * drafts.
 * @typedef {{ recordKey: string, card: CardOrientation, drafts: TargetDraftCheckpoint[] }} AlignmentOrientationCheckpoint
 */

/**
 * A resolved record transform to write on a row: the start, the orientation and
 * the provenance of an eligible anchor resolution (`resolvedTransform`,
 * record-display/feature-record-rotation.js, its one writer).
 * `writeResolvedTransform` checks each field again before a draft takes it.
 * @typedef {{ startCoordinate: number, reverseComplement: boolean, anchorIntent: AnchorIntent }} ResolvedRecordTransform
 */

/**
 * A display start a record placement writes, with the committed row whose draft takes it (`recordStartTargets`).
 * @typedef {{ row: CommittedDisplayRow, startCoordinate: number }} RecordStartTarget
 */

/**
 * The draft `captureTargetDraft` keeps: a record's identity and the mode whose drawing holds its draft.
 * A row is one; a draft names no mode, so its caller adds it.
 * @typedef {RecordDisplayDraftIdentity & { scope: string, drawingId: DrawingId }} RecordDraftTarget
 */

/**
 * A record row with the drawing that holds its drafts (the source it was built from carries it).
 * @typedef {RecordDisplayRow & { drawingId: DrawingId }} DrawnRecordRow
 */

/**
 * What `captureTargetDraft` keeps of a row's draft: restoring puts it back.
 * @typedef {object} TargetDraftCheckpoint
 * @property {string} scope the kind of row (`circular` source or `linear` sequence)
 * @property {DrawingId} drawingId the drawing whose draft it is, found by identity
 * @property {string} key
 * @property {number} index
 * @property {RecordDisplayDraft | null} draft
 * @property {string} [ownerUid]
 * @property {CardOrientation} [ownerOrientation] the orientation fields of the File card `ownerUid`
 */

/**
 * @typedef {object} RecordDisplayControlsOptions
 * @property {Pick<WebState, 'activeDrawing' | 'drawingList' | 'drawingOf' | 'featureCatalog' | 'inputFiles'
 *   | 'processing' | 'semanticFileWatchersSuppressed' | 'sessionImportPending' | 'sessionOperationAvailability'>} state
 *   the members of the Web state (`state.js`) the controls read
 * @property {<T>(getter: () => T) => { value: T }} computed Vue `computed`
 * @property {(source: any, callback: (...args: any[]) => void, options?: Record<string, any>) => any} watch Vue `watch`
 * @property {(drawing: DrawingState, sequence: LinearSeq) => SequenceRecord[]} linearRecordsFor a card's records, read in the drawing that holds it
 * @property {LinearWindowsOf} linearWindowsOf the windows a Linear drawing's request writes for a card
 * @property {(drawing: DrawingState, sequence: LinearSeq | undefined) => RecordSelectorStatus} linearRecordStatusFor
 * @property {(drawing: DrawingState, sequence: LinearSeq | undefined) => string} linearRecordErrorFor the card's error text, '' when none
 * @property {(label: string, apply: () => unknown) => Promise<unknown>} runUndoable History's one undoable step:
 *   the result of `apply`, or History's busy value (`{ status: 'busy', reason }`) when it refuses the step without running it
 * @property {() => Record<string, any> | null} getCommittedRequest
 * @property {() => Record<string, any> | null} getCommittedSession
 */

// R13: the Linear record selector's three readers and History's undoable step
// arrive as ports; this owner holds neither owner.
/**
 * A source's record rows, each with its drawing, crop, orientation and pair; the selector '' builds a row for every
 * record of the source's files.
 * @param {DrawingRecordSource} source @param {unknown} [selector]
 */
const sourceRows = (source, selector = source.selector) => buildRecordDisplayRows({ ...source, selector }).map((row) => (
  { ...row, drawingId: source.drawingId, ...source.cropOf(row), reverse: source.reverse, paired: source.paired }));
/**
 * The record rows a drawing shows, from its sources (`sourcesOf`). A source whose selector does not read builds none:
 * its selector control shows the error.
 * @param {readonly DrawingRecordSource[]} sources
 */
const shownRows = (sources) => sources.filter((source) => source.shown).flatMap((source) => {
  try { return source.source ? sourceRows(source) : []; } catch { return []; }
});
/** The rows a drawing shows (`shownRows`). @typedef {ReturnType<typeof shownRows>} ShownRows */
/** @typedef {(drawing: DrawingState) => ShownRows} ShownRowsOf */
/**
 * The rows each drawing shows, built once per drawing for one use: a capture of many records reads them once.
 * @param {(drawing: DrawingState) => DrawingRecordSource[]} sourcesOf
 * @returns {ShownRowsOf}
 */
const shownRowsOnce = (sourcesOf) => {
  /** @type {Map<DrawingState, ShownRows>} */
  const built = new Map();
  return (drawing) => {
    if (!built.has(drawing)) built.set(drawing, shownRows(sourcesOf(drawing)));
    return /** @type {ShownRows} */ (built.get(drawing));
  };
};

/**
 * @param {RecordDisplayControlsOptions} options
 */
export const createRecordDisplayControls = ({
  state, computed, watch, linearRecordsFor, linearWindowsOf, linearRecordStatusFor, linearRecordErrorFor, runUndoable,
  getCommittedRequest, getCommittedSession
}) => {
  /** @type {CommittedDisplayRow[]} */
  let committedRows = [];
  /** @type {Record<string, any> | null} */
  let boundRequest = null;
  /** @type {Pick<DrawingRecordSource, 'drawingId' | 'scope' | 'sourceUid' | 'source' | 'paired'>[]} */
  let boundSources = [];
  // The record sources of one drawing, each named by that drawing's ID: a
  // Circular drawing's source cards with its crop and orientation, a Linear
  // drawing's cards: each read from that drawing's inputs, by card uid.
  /** @param {DrawingState} drawing @returns {DrawingRecordSource[]} */
  const sourcesOf = (drawing) => (drawing.mode === 'circular'
    ? circularDisplaySources(drawing) : linearDisplaySources(drawing, linearRecordsFor, linearWindowsOf));
  // The sources of every drawing, in the drawings' order.
  const sources = computed(() => state.drawingList().flatMap(sourcesOf));
  const allRows = computed(() => sources.value.flatMap((source) => (source.source ? sourceRows(source, '') : [])));
  const rows = computed(() => shownRows(sourcesOf(state.activeDrawing())));
  // A drawing's drafts on the rows of their sources (OV-278). An edit of a
  // draft first writes them there, within its own step.
  /** @param {DrawingState} drawing */
  const boundDrafts = (drawing) => bindRecordDisplayDrafts(drawing.recordDisplayDrafts, allRows.value, drawing.id);
  /** @param {DrawingState} drawing */
  const settleDrafts = (drawing) => {
    const bound = boundDrafts(drawing);
    if (bound !== drawing.recordDisplayDrafts) drawing.recordDisplayDrafts.splice(0, bound.length, ...bound);
  };
  /** @param {DrawnRecordRow} row */
  const draftFor = (row) => boundDrafts(state.drawingOf(row.drawingId)).find((/** @type {RecordDisplayDraft} */ draft) => recordDisplayKey(draft) === row.key) || {};
  // Linear orientation has one owner, the File card: the windows it lists, else
  // its region_reverse (`orientLinearCard`, as Generate reads them; CO-03,
  // N-09). A row that is the only record its card draws writes it there; the
  // row override remains only for Circular and for a Linear card that draws
  // several records.
  // The rows are counted in the row's own drawing: another drawing can hold a card of the same uid.
  /**
   * @param {DrawingState} drawing the row's drawing @param {{ scope: string, sourceUid: string } | null | undefined} row
   * @param {ShownRowsOf} [rowsOf] The rows of a drawing; a capture of many records passes one that builds them once.
   */
  const linearOrientationOwner = (drawing, row, rowsOf = shownRowsOnce(sourcesOf)) => {
    if (row?.scope !== 'linear') return null;
    const sequence = drawing.linearSeqs.find((/** @type {LinearSeq} */ seq) => seq.uid === row.sourceUid);
    return sequence && rowsOf(drawing).filter((entry) => entry.scope === 'linear'
      && entry.sourceUid === row.sourceUid).length === 1 ? sequence : null;
  };
  /** @param {DrawingState} drawing @param {RecordDisplayRow} row @param {RecordDisplayDraft} draft @param {boolean | null} reverseComplement */
  const writeReverseComplement = (drawing, row, draft, reverseComplement) => {
    const owner = reverseComplement === null ? null : linearOrientationOwner(drawing, row);
    if (owner && reverseComplement !== null) orientLinearCard(owner, reverseComplement);
    draft.reverseComplementOverride = owner ? null : reverseComplement;
  };
  // A row's draft in its drawing, added when it has none yet.
  /** @param {DrawingState} drawing @param {RecordDisplayRow} row @returns {RecordDisplayDraft} */
  const draftOf = (drawing, row) => {
    const found = drawing.recordDisplayDrafts.find((/** @type {RecordDisplayDraft} */ entry) => recordDisplayKey(entry) === row.key);
    if (found) return found;
    drawing.recordDisplayDrafts.push({ sourceUid: row.sourceUid, selector: row.selector,
      recordId: row.recordId, topologyOverride: null, startCoordinate: null,
      reverseComplementOverride: null, anchorIntent: null });
    return drawing.recordDisplayDrafts[drawing.recordDisplayDrafts.length - 1];
  };
  /** @param {DrawnRecordRow} row */
  const surfaceFor = (row) => {
    const surface = recordDisplaySurface(row, draftFor(row), row);
    const busy = state.sessionOperationAvailability();
    return busy ? { ...surface, startEnabled: false, disabledReason: busy.reason } : surface;
  };
  /** @param {DrawingState} drawing @param {RecordDisplayRow} row @param {Partial<RecordDisplayDraft>} patch @param {string} label */
  const edit = (drawing, row, patch, label) => runUndoable(label, () => {
    settleDrafts(drawing);
    const draft = draftOf(drawing, row);
    const { reverseComplementOverride, ...rest } = patch;
    Object.assign(draft, rest);
    if (Object.hasOwn(patch, 'reverseComplementOverride')) {
      // Every patch that names the key carries a boolean or null.
      writeReverseComplement(drawing, row, draft, /** @type {boolean | null} */ (reverseComplementOverride));
    }
  });
  /** @param {DrawnRecordRow} row @param {ResolvedRecordTransform} transform */
  const writeResolvedTransform = (row, { startCoordinate, reverseComplement, anchorIntent }) => {
    const drawing = state.drawingOf(row.drawingId);
    const busy = state.sessionOperationAvailability();
    if (busy) return busy;
    if (typeof reverseComplement !== 'boolean') {
      throw new Error('Resolved record orientation must be boolean.');
    }
    const resolvedStart = parseRecordDisplayStart(startCoordinate);
    const resolvedIntent = validateAnchorIntent(anchorIntent);
    settleDrafts(drawing);
    const draft = draftOf(drawing, row);
    Object.assign(draft, { startCoordinate: resolvedStart, anchorIntent: resolvedIntent });
    writeReverseComplement(drawing, row, draft, reverseComplement);
  };
  /**
   * @template {{ recordKey: string }} O
   * @param {DrawingState} drawing
   * @param {O[]} orientations
   */
  const alignmentRows = (drawing, orientations) => {
    refreshCommittedRows();
    const request = getCommittedRequest();
    return orientations.map((orientation) => {
      const records = request?.mode === 'linear'
        ? request.records.filter((/** @type {CommittedRecord} */ record) => record.recordKey === orientation.recordKey) : [];
      const inputs = sourcesOf(drawing).filter(source => source.scope === 'linear'
        && source.sourceUid === orientation.recordKey);
      // The card of the one input, in the drawing's rows.
      const sequence = inputs.length === 1
        ? drawing.linearSeqs.find((card) => card.uid === inputs[0].sourceUid) : undefined;
      if (records.length !== 1 || !sequence) {
        throw new Error('Alignment record transform target is stale.');
      }
      const record = records[0];
      const input = inputs[0];
      // The resource of each file of the record's source, in source order: the File card's, then a GFF3 file's FASTA.
      const [primaryId, pairedId] = storedRecordSourceResourceIds(record.source);
      if (!boundSources.some(source => source.drawingId === input.drawingId && source.scope === input.scope
          && source.sourceUid === input.sourceUid && source.source === input.source && source.paired === input.paired)
        || !input.source || !matchesSavedSource(input.source, primaryId)
        || (pairedId !== undefined && (!input.paired || !matchesSavedSource(input.paired, pairedId)))) {
        throw new Error('Alignment record source binding changed.');
      }
      const selector = record.selector;
      const drafts = drawing.recordDisplayDrafts.filter((/** @type {RecordDisplayDraft} */ draft) => draft.sourceUid === input.sourceUid
        && (selector?.kind === 'recordIndex'
          ? draft.selector === `#${selector.index + 1}`
          : selector?.kind === 'recordId' ? draft.recordId === selector.value : draft.selector === '#1'));
      return { recordKey: orientation.recordKey, orientation, drafts, sequence };
    });
  };
  /** @param {DrawingState} drawing the Linear drawing the alignment run pinned @param {AlignmentOrientation[]} orientations @param {ShownRowsOf} [rowsOf] The rows of the drawing, built once for the capture. @returns {AlignmentOrientationCheckpoint[]} */
  const captureAlignmentOrientationIntent = (drawing, orientations, rowsOf = shownRowsOnce(sourcesOf)) => {
    // A draft names no mode: these are the Linear drawing's (PD-OI-090, OV-579).
    // On their rows first (OV-278): a draft that the capture's own settle moved would not be found.
    settleDrafts(drawing);
    return alignmentRows(drawing, orientations)
      .map(({recordKey, sequence, drafts}) => ({ recordKey, card: captureLinearCardOrientation(sequence),
        drafts: drafts.map((/** @type {RecordDisplayDraft} */ draft) => captureTargetDraft({ ...draft, scope: 'linear', drawingId: drawing.id }, rowsOf)) }));
  };
  /** @param {DrawingState} drawing @param {AlignmentOrientationCheckpoint[]} checkpoints */
  const restoreAlignmentOrientationIntent = (drawing, checkpoints) => {
    alignmentRows(drawing, checkpoints).forEach(({sequence}, index) => {
      restoreLinearCardOrientation(sequence, checkpoints[index].card);
      checkpoints[index].drafts.forEach(restoreTargetDraft);
    });
  };
  /** @param {DrawingState} drawing @param {AlignmentOrientation[]} orientations */
  const commitAlignmentOrientations = (drawing, orientations) => {
    alignmentRows(drawing, orientations).forEach(({sequence, orientation, drafts}) => {
      orientLinearCard(sequence, orientation.reverseComplement);
      drafts.forEach((/** @type {RecordDisplayDraft} */ draft) => { draft.reverseComplementOverride = null; });
    });
  };
  /** @param {unknown} file the source file of a File card (`isCurrentFeature` passes it on untyped) @param {string | undefined} resourceId */
  const matchesSavedSource = (file, resourceId) => {
    // A resource id is Session data: only an own key of the table names a saved resource.
    const expected = ownEntry(getCommittedSession()?.resources, String(resourceId));
    // A File card's source is its binding's ref (or null); its bytes are the store's.
    return matchesSessionResourceDescriptor(state.inputFiles.fileOf(/** @type {InputFileRef | null} */ (file)), expected);
  };
  const refreshCommittedRows = () => {
    if (state.semanticFileWatchersSuppressed?.value || state.sessionImportPending?.value) return;
    const request = getCommittedRequest();
    if (!request) { committedRows = []; boundRequest = null; return; }
    if (request !== boundRequest) {
      boundRequest = request;
      // The committed request is the shown drawing's: its cards, by drawing and uid (another drawing, one that Find in
      // other files created, can hold cards of the same uids; OV-3375).
      boundSources = sourcesOf(state.activeDrawing())
        .map(({ drawingId, scope, sourceUid, source, paired }) => ({ drawingId, scope, sourceUid, source, paired }));
    }
    committedRows = allRows.value.filter((row) => row.scope === request.mode
      && boundSources.some((source) => source.drawingId === row.drawingId && source.scope === row.scope
        && source.sourceUid === row.sourceUid && source.source === row.source && source.paired === row.paired)).flatMap((row) => {
      const selected = request.records.find((/** @type {CommittedRecord} */ record) => {
        const selector = record.selector;
        const matchesSelector = !selector || (selector.kind === 'recordIndex'
          ? row.selector === `#${selector.index + 1}` : row.recordId === selector.value);
        // A record of a card is keyed by the card's uid or its identity `<uid>:<k>`, in both modes.
        return matchesSelector
          && (record.recordKey === row.sourceUid || record.recordKey === `${row.sourceUid}:${Number(row.selector.slice(1))}`);
      });
      if (!selected) return [];
      const source = selected.source;
      const [primaryId, pairedId] = storedRecordSourceResourceIds(source);
      if (!matchesSavedSource(row.source, primaryId) || (row.paired && !matchesSavedSource(row.paired, pairedId))) return [];
      return [{ ...row,
        committedDisplay: selected.display || { isCircular: null, startCoordinate: null },
        committedReverseComplement: canonicalRecordReverseComplement(selected),
        committedCropped: selected.windows.length > 0,
        canonicalRecordKey: selected.recordKey,
        canonicalSource: selected.source,
        recordKey: selected.cardinality === 'all'
        && allRows.value.filter((entry) => entry.drawingId === row.drawingId && entry.sourceUid === row.sourceUid).length > 1
        ? `${selected.recordKey}:${Number(row.selector.slice(1))}` : selected.recordKey }];
    });
  };
  watch([() => state.featureCatalog.value, () => state.processing?.value, allRows], refreshCommittedRows, { flush: 'post' });
  // The transition of every write of a drawing's row source (a file, an input
  // type, the Linear cards): `mutate` writes it, then the drafts of each source
  // of `drawing` it replaced or removed go from `drawing` (R10). Load, its
  // rollback, and a History restore install the drafts they read or captured.
  /** @template T @param {DrawingState} drawing @param {() => T} mutate @returns {T} */
  const changeRecordSources = (drawing, mutate) => {
    const previous = sourcesOf(drawing);
    const result = mutate();
    const current = sourcesOf(drawing);
    const drafts = drawing.recordDisplayDrafts;
    previous.filter((old) => {
      const next = current.find((entry) => entry.sourceUid === old.sourceUid);
      return old.source && (!next || next.source !== old.source || next.paired !== old.paired);
    }).forEach((replaced) => {
      for (let index = drafts.length - 1; index >= 0; index -= 1) {
        if (drafts[index].sourceUid === replaced.sourceUid) drafts.splice(index, 1);
      }
    });
    return result;
  };
  const hasPendingChanges = computed(() => {
    const drawing = state.activeDrawing();
    void state.processing?.value;
    refreshCommittedRows();
    const request = getCommittedRequest();
    if (!request) return false;
    const options = request.diagramOptions || {};
    const tolerance = options.configOverrides?.['canvas.feature_overlap_tolerance_bp']
      ?? options.config?.canvas?.feature_overlap_tolerance_bp ?? 0;
    const placements = options.featurePlacements || [];
    // The draft keeps the rows of other records (R2); compare the rows this request carries.
    const draftPlacementCount = () => {
      try { return requestFeaturePlacements(drawing.featurePlacementOverrides, request.mode, request.records || []).length; }
      catch { return -1; }
    };
    if (request.mode !== drawing.mode || tolerance !== draftFeatureOverlapTolerance(drawing.adv)
      || placements.length !== draftPlacementCount()
      || placements.some((/** @type {{ recordKey: string, biologicalFeatureId: string, placement: Record<string, unknown> }} */ row) => {
        /** @type {Readonly<Record<string, unknown>> | undefined} */
        const target = drawing.featurePlacementOverrides[
          featureIdentityKey(row.recordKey, row.biologicalFeatureId)]?.placement;
        return !target || ['kind', 'side', 'level'].some((field) => target[field] !== row.placement[field]);
      })) return true;
    return committedRows.some((row) => {
      try {
        return JSON.stringify(requestedRecordTransform(row, draftFor(row), row))
          !== JSON.stringify({
            display: row.committedDisplay,
            reverseComplement: row.committedReverseComplement
          });
      }
      catch { return true; }
    });
  });
  /** @param {DrawingRecordSource} source */
  const discoveryFor = (source) => {
    if (source.scope === 'circular') return circularDiscoveryForInput(state.drawingOf(source.drawingId));
    const drawing = state.drawingOf(source.drawingId);
    const sequence = drawing.linearSeqs.find((/** @type {LinearSeq} */ seq) => seq.uid === source.sourceUid);
    return { status: linearRecordStatusFor(drawing, sequence), error: linearRecordErrorFor(drawing, sequence) };
  };
  // The committed record's current sources have no read records: say whether
  // the read failed or has not happened, instead of calling the target stale.
  /** @param {string} recordKey */
  const unreadRecordsError = (recordKey) => {
    // The committed request is the adopted canonical request: its `records` is always set.
    const request = /** @type {CommittedRecordRequest | null} */ (getCommittedRequest());
    const record = committedRecordForKey(request, recordKey);
    // The committed request is the shown drawing's: its sources, not another drawing's of the same uids.
    const discoveries = record ? sourcesOf(state.activeDrawing())
      .filter((source) => source.source && committedRecordUsesSource(request, record, source)).map(discoveryFor) : [];
    if (!discoveries.length || discoveries.some(({ status }) => status === 'ready')) return null;
    const failed = discoveries.find(({ status }) => status === 'error');
    return failed
      ? recordTargetError('discovery-failed',
        (typeof failed.error === 'string' ? failed.error : failed.error && errorText(failed.error.summaryParts)) || tm(RECORD_READ_ERROR_LABEL))
      : recordTargetError(RECORD_TARGET_NOT_DISCOVERED, t('The records of this feature have not been read yet.'));
  };
  // The display start and orientation that Generate would apply from this
  // committed row's draft, when either differs from the committed request.
  /** @param {CommittedDisplayRow} row */
  const pendingTransformFor = (row) => {
    let requested;
    try { requested = requestedRecordTransform(row, draftFor(row), row); } catch { return null; }
    const startCoordinate = requested.display.startCoordinate;
    return startCoordinate === (row.committedDisplay?.startCoordinate ?? null)
      && requested.reverseComplement === row.committedReverseComplement
      ? null
      : { startCoordinate, reverseComplement: requested.reverseComplement };
  };
  /** @param {Record<string, unknown> | null | undefined} feature */
  const targetForFeature = (feature) => {
    refreshCommittedRows();
    const recordKey = String(feature?.record_key || '');
    const matches = committedRows.filter((row) => row.recordKey === recordKey);
    if (matches.length !== 1) {
      const unread = !matches.length && unreadRecordsError(recordKey);
      if (unread) throw unread;
      throw new Error('The popup feature target is stale or ambiguous.');
    }
    const row = matches[0];
    const members = committedRows
      .filter((member) => member.scope === row.scope
        && member.sourceUid === row.sourceUid
        && member.canonicalRecordKey === row.canonicalRecordKey)
      .map((member) => ({
        canonicalRecordKey: member.canonicalRecordKey,
        recordKey: member.recordKey,
        selector: member.selector,
        recordId: member.recordId,
        recordLength: member.recordLength,
        committedDisplay: { ...member.committedDisplay },
        committedReverseComplement: member.committedReverseComplement
      }));
    return {
      row,
      target: {
        scope: row.scope,
        sourceUid: row.sourceUid,
        selector: row.selector,
        recordId: row.recordId,
        recordLength: row.recordLength,
        recordKey: row.recordKey,
        canonicalRecordKey: row.canonicalRecordKey,
        source: { ...row.canonicalSource },
        effectiveCircular: row.committedDisplay?.isCircular
          ?? row.detectedTopology === 'circular',
        cropped: row.committedCropped,
        committedReverseComplement: row.committedReverseComplement,
        pendingTransform: pendingTransformFor(row),
        members
      }
    };
  };
  // Orientations decided for the committed request's records (by record key),
  // written as one History step the way the record popup's Reverse complement
  // writes one: the File card's checkbox, or the row of a card that draws
  // several records. They were decided from the comparisons the last Generate
  // drew, so a row whose orientation changed since then is stale.
  /** @param {readonly RecordOrientationDecision[]} decisions @param {string} label @returns {Promise<RecordOrientationWrite>} */
  const writeDecidedOrientations = async (decisions, label) => {
    // The button acts on the shown drawing, whose last Generate it read.
    const drawing = state.activeDrawing();
    if (drawing.mode !== 'linear') return { status: 'notLinear' };
    const busy = state.sessionOperationAvailability();
    if (busy) return { status: 'busy', reason: busy.reason };
    refreshCommittedRows();
    const rowsByKey = new Map(committedRows.filter((row) => row.scope === 'linear').map((row) => [row.recordKey, row]));
    const drawnChanged = (/** @type {CommittedDisplayRow} */ row) => {
      try { return requestedRecordTransform(row, draftFor(row), row).reverseComplement !== row.committedReverseComplement; }
      catch { return true; }
    };
    if (decisions.length !== rowsByKey.size || decisions.some(({ recordKey }) => !rowsByKey.has(recordKey))
      || [...rowsByKey.values()].some(drawnChanged)) return { status: 'stale' };
    const changes = decisions.flatMap(({ recordKey, reverseComplement }) => {
      const row = /** @type {CommittedDisplayRow} */ (rowsByKey.get(recordKey));
      return reverseComplement === row.committedReverseComplement ? [] : [{ row, reverseComplement }];
    });
    if (!changes.length) return { status: 'unchanged' };
    const WRITTEN = Symbol('written');
    // History refuses the step, without running it, while a Session operation blocks edits.
    const outcome = await runUndoable(label, () => {
      settleDrafts(drawing);
      changes.forEach(({ row, reverseComplement }) => {
        const draft = draftOf(drawing, row);
        draft.anchorIntent = null;
        writeReverseComplement(drawing, row, draft, reverseComplement);
      });
      return WRITTEN;
    });
    if (outcome !== WRITTEN) {
      // History refuses a step only with its busy value, whose reason it always words.
      const { reason } = /** @type {{ reason?: unknown }} */ (outcome && typeof outcome === 'object' ? outcome : {});
      if (typeof reason !== 'string' || !reason) {
        throw diagnosticError('RESULT_INVALID', {}, { stage: 'helper', operation: 'decideAutoOrientation' });
      }
      return { status: 'busy', reason };
    }
    return { status: 'written', recordIds: changes.map(({ row }) => String(row.recordId || row.recordKey)) };
  };
  /**
   * @param {RecordDraftTarget} row @param {ShownRowsOf} [rowsOf] The rows of a drawing, for the orientation owner.
   * @returns {TargetDraftCheckpoint}
   */
  const captureTargetDraft = (row, rowsOf) => {
    const drawing = state.drawingOf(row.drawingId);
    settleDrafts(drawing);
    const key = recordDisplayKey(row);
    const index = drawing.recordDisplayDrafts.findIndex(
      (/** @type {RecordDisplayDraft} */ draft) => recordDisplayKey(draft) === key
    );
    const owner = linearOrientationOwner(drawing, row, rowsOf);
    return {
      scope: row.scope,
      drawingId: row.drawingId,
      key,
      index,
      draft: index >= 0 ? cloneJsonData(drawing.recordDisplayDrafts[index]) : null,
      ...(owner ? { ownerUid: owner.uid, ownerOrientation: captureLinearCardOrientation(owner) } : {})
    };
  };
  /** @param {TargetDraftCheckpoint} checkpoint */
  const restoreTargetDraft = (checkpoint) => {
    const drawing = state.drawingOf(checkpoint.drawingId);
    if (checkpoint.ownerUid && checkpoint.ownerOrientation) {
      const owner = drawing.linearSeqs.find((/** @type {LinearSeq} */ seq) => seq.uid === checkpoint.ownerUid);
      if (owner) restoreLinearCardOrientation(owner, checkpoint.ownerOrientation);
    }
    const index = drawing.recordDisplayDrafts.findIndex(
      (/** @type {RecordDisplayDraft} */ draft) => recordDisplayKey(draft) === checkpoint.key
    );
    if (index >= 0) drawing.recordDisplayDrafts.splice(index, 1);
    if (checkpoint.draft) {
      const insertAt = Math.max(0, Math.min(
        checkpoint.index,
        drawing.recordDisplayDrafts.length
      ));
      drawing.recordDisplayDrafts.splice(insertAt, 0, cloneJsonData(checkpoint.draft));
    }
  };
  // The display starts a record placement writes (Align at feature's rotated records), each with the committed row of
  // its record key in the Linear drawing, so Generate draws them from the drafts (the caller reads the records
  // first, so a loaded Session's records have their rows). A record without one such row has no File card of its
  // own: the start is then refused, not written to a draft that the next Generate would not read.
  /** @param {DrawingState} drawing @param {readonly RecordPlacementStart[]} starts @returns {RecordStartTarget[]} */
  const recordStartTargets = (drawing, starts) => {
    if (starts.length) refreshCommittedRows();
    return starts.map(({ recordKey, startCoordinate }) => {
      const matches = committedRows.filter((row) => row.drawingId === drawing.id && row.recordKey === recordKey);
      if (matches.length !== 1) {
        throw diagnosticError('INPUT_INVALID', { field: 'anchor_feature', reason: 'ANCHOR_RECORD_CARDS' });
      }
      return { row: matches[0], startCoordinate: /** @type {number} */ (parseRecordDisplayStart(startCoordinate)) };
    });
  };
  /**
   * @param {readonly RecordStartTarget[]} targets
   * @param {ShownRowsOf} [rowsOf] The rows of each drawing, built once for the whole capture.
   * @returns {TargetDraftCheckpoint[]}
   */
  const captureRecordStarts = (targets, rowsOf = shownRowsOnce(sourcesOf)) => targets.map(({ row }) => captureTargetDraft(row, rowsOf));
  /** @param {readonly TargetDraftCheckpoint[]} checkpoints */
  const restoreRecordStarts = (checkpoints) => checkpoints.forEach(restoreTargetDraft);
  // The start is the record's own now: like an edit of its Display start, it replaces any earlier anchor intent.
  /** @param {DrawingState} drawing @param {readonly RecordStartTarget[]} targets */
  const commitRecordStarts = (drawing, targets) => {
    if (!targets.length) return;
    settleDrafts(drawing);
    targets.forEach(({ row, startCoordinate }) => {
      Object.assign(draftOf(drawing, row), { startCoordinate, anchorIntent: null });
    });
  };
  return { rows, allRows, draftFor, surfaceFor, hasPendingChanges, changeRecordSources,
    // The data of isCurrentFeature (services/feature-identity.js), refreshed
    // for the committed request before it is read (R13).
    sourceBinding: () => {
      refreshCommittedRows();
      // The committed request is the adopted canonical request: its `records` is always set.
      return { request: /** @type {CommittedRecordRequest | null} */ (getCommittedRequest()),
        sources: sourcesOf(state.activeDrawing()), boundSources, matchesSavedSource };
    },
    // The drawn records of the committed request, each with its File card
    // (`sourceUid`) and its request record (`canonicalRecordKey`).
    committedRecordRows: () => {
      refreshCommittedRows();
      return committedRows;
    },
    rowsFor: (/** @type {string} */ uid) => rows.value.filter((row) => row.sourceUid === uid),
    setTopology: (/** @type {RecordDisplayRow} */ row, /** @type {boolean | null} */ value) => {
      const drawing = state.activeDrawing();
      return edit(drawing, row, { topologyOverride: value }, msg('Change record topology'));
    },
    setStart: (/** @type {RecordDisplayRow} */ row, /** @type {unknown} */ value) => {
      const drawing = state.activeDrawing();
      return state.sessionOperationAvailability() || edit(drawing, row, {
        startCoordinate: parseRecordDisplayStart(value), anchorIntent: null
      }, msg('Change record display start'));
    },
    resetStart: (/** @type {RecordDisplayRow} */ row) => {
      const drawing = state.activeDrawing();
      return edit(drawing, row, { startCoordinate: null, anchorIntent: null }, msg('Reset record display start'));
    },
    setReverseComplement: (/** @type {RecordDisplayRow} */ row, /** @type {unknown} */ value) => {
      const drawing = state.activeDrawing();
      return state.sessionOperationAvailability() || edit(drawing, row, {
        reverseComplementOverride: requireReverseComplementOverride(value), anchorIntent: null
      }, msg('Change record orientation'));
    },
    // Apply on Generate: one undoable draft step, like the sidebar's edits.
    setResolvedTransform: (/** @type {DrawnRecordRow} */ row, /** @type {ResolvedRecordTransform} */ transform) => runUndoable(
      msg('Rotate record to feature on Generate'),
      () => writeResolvedTransform(row, transform)
    ),
    commitResolvedTransform: writeResolvedTransform,
    captureAlignmentOrientationIntent,
    restoreAlignmentOrientationIntent,
    commitAlignmentOrientations,
    writeDecidedOrientations,
    captureTargetDraft,
    restoreTargetDraft,
    recordStartTargets,
    captureRecordStarts,
    restoreRecordStarts,
    commitRecordStarts,
    targetForFeature };
};
