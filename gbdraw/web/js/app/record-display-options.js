// @ts-check
/** @import { DrawingState } from '../state.js' */
/** @import { CommittedRecord, CommittedRecordRequest } from '../services/feature-identity.js' */
// Source-bound editable rotation intent. The request service owns serialization.
import {
  bindRecordDisplayDrafts,
  buildRecordDisplayRows,
  parseRecordDisplayStart,
  recordDisplayKey,
  recordDisplaySurface,
  requestedRecordTransform,
  requireReverseComplementOverride,
  validateAnchorIntent
} from '../services/record-display-model.js';
import { matchesSessionResourceDescriptor } from '../services/session-resource-backing.js';
import { cloneJsonData } from '../services/json-clone.js';
import { featureIdentityKey, requestFeaturePlacements } from '../services/feature-placement.js';
import { committedRecordForKey, committedRecordUsesSource } from '../services/feature-identity.js';
import { circularDiscoveryForInput } from '../services/record-discovery.js';
import { circularDrawingCropped } from '../services/record-options.js';
/** @import { DrawingId } from '../services/drawing-view.js' */
import { RECORD_READ_ERROR_LABEL } from './linear-record-selector.js';
import { ownEntry } from '../utils/own-entry.js';
import { draftFeatureOverlapTolerance } from '../services/session-active-config-contract.js';
import { msg, t, tm } from '../utils/i18n.js';
import { errorText } from '../utils/error-text.js';
/** @import { RecordDisplaySource } from '../services/feature-identity.js' */
/** @import { LinearSeq } from '../state.js' */
/** @import { InputFileRef } from '../services/history-files.js' */
/** @import { AnchorIntent, RecordDisplayDraft, RecordDisplayDraftIdentity, RecordDisplayRow, RecordDisplayScope } from '../services/record-display-model.js' */
/** @import { RecordEntryInput } from '../services/record-options.js' */

// Session Load reads no record bytes (776a2f93), so a loaded Result can have a
// current source whose records are not read yet. That popup target failure is
// not stale: its kind lets an explicit record action read them.
export const RECORD_TARGET_NOT_DISCOVERED = 'not-discovered';
/** @param {string} kind @param {string} message */
const recordTargetError = (kind, message) => Object.assign(new Error(message), { kind });

/**
 * An entry of the Circular record list, as the record discovery stores it (`run-analysis.js`).
 * @typedef {object} CircularListRecord
 * @property {string} [selector]
 * @property {string} [record_id]
 * @property {string} [recordId]
 * @property {number | null} [record_length]
 * @property {number | null} [recordLength]
 * @property {unknown} [detectedTopology]
 * @property {string} [recordKey]
 */

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
 * @property {boolean} cropped
 * @property {boolean} reverse
 */

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
 * What `captureAlignmentOrientationIntent` keeps of an orientation: the File card's reverse and the record's drafts.
 * @typedef {AlignmentOrientation & { drafts: TargetDraftCheckpoint[] }} AlignmentOrientationCheckpoint
 */

/**
 * A resolved record transform to write on a row: the start, the orientation and
 * the provenance of an eligible anchor resolution (`resolvedTransform`,
 * record-display/feature-record-rotation.js, its one writer).
 * `writeResolvedTransform` checks each field again before a draft takes it.
 * @typedef {{ startCoordinate: number, reverseComplement: boolean, anchorIntent: AnchorIntent }} ResolvedRecordTransform
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
 * @property {boolean} [ownerReverse]
 */

/**
 * @typedef {object} RecordDisplayControlsOptions
 * @property {Record<string, any>} state the Web state; its shape belongs to `state.js`
 * @property {<T>(getter: () => T) => { value: T }} computed Vue `computed`
 * @property {(source: any, callback: (...args: any[]) => void, options?: Record<string, any>) => any} watch Vue `watch`
 * @property {(sequence: Record<string, any>) => Record<string, any>[]} linearRecordsFor
 * @property {(sequence: Record<string, any> | undefined) => string} linearRecordStatusFor
 * @property {(sequence: Record<string, any> | undefined) => any} linearRecordErrorFor
 * @property {(label: string, apply: () => any) => any} runUndoable History's one undoable step
 * @property {() => Record<string, any> | null} getCommittedRequest
 * @property {() => Record<string, any> | null} getCommittedSession
 */

// R13: the Linear record selector's three readers and History's undoable step
// arrive as ports; this owner holds neither owner.
/**
 * @param {RecordDisplayControlsOptions} options
 */
export const createRecordDisplayControls = ({
  state, computed, watch, linearRecordsFor, linearRecordStatusFor, linearRecordErrorFor, runUndoable,
  getCommittedRequest, getCommittedSession
}) => {
  /** @type {CommittedDisplayRow[]} */
  let committedRows = [];
  /** @type {Record<string, any> | null} */
  let boundRequest = null;
  /** @type {RecordDisplaySource[]} */
  let boundSources = [];
  // The record sources of one drawing, each named by that drawing's ID: a
  // Circular drawing's file card with its crop and orientation, a Linear
  // drawing's cards: each read from that drawing's inputs.
  /** @param {DrawingState} drawing @returns {DrawingRecordSource[]} */
  const sourcesOf = (drawing) => (drawing.mode === 'circular'
    ? [{ scope: 'circular', drawingId: drawing.id, sourceUid: 'circular',
        source: state.cInputType.value === 'gff' ? drawing.files.c_gff : drawing.files.c_gb,
        paired: state.cInputType.value === 'gff' ? drawing.files.c_fasta : null,
        records: drawing.circularRecordList.value.map((/** @type {CircularListRecord} */ record) => ({ ...record,
          recordId: record.record_id ?? record.recordId, recordLength: record.record_length ?? record.recordLength })), selector: drawing.form.circular_record_selector,
        cropped: circularDrawingCropped(drawing.circularRecordList.value, drawing),
        reverse: drawing.form.circular_reverse }]
    : drawing.linearSeqs.map((/** @type {LinearSeq} */ seq) => ({ scope: 'linear', drawingId: drawing.id, sourceUid: seq.uid,
        source: state.lInputType.value === 'gff' ? seq.gff : seq.gb,
        paired: state.lInputType.value === 'gff' ? seq.fasta : null,
        records: linearRecordsFor(seq), selector: seq.region_record_id,
        cropped: Boolean(seq.region_start || seq.region_end), reverse: seq.region_reverse })));
  // The sources of every drawing, in the drawings' order.
  const sources = computed(() => {
    /** @type {DrawingState[]} */
    const drawings = state.drawingList();
    return drawings.flatMap(sourcesOf);
  });
  const allRows = computed(() => sources.value.flatMap((source) => source.source
    ? buildRecordDisplayRows({ ...source, selector: '' }).map((row) => ({ ...row, drawingId: source.drawingId,
        cropped: source.cropped, reverse: source.reverse, paired: source.paired })) : []));
  const rows = computed(() => sourcesOf(state.activeDrawing())
    .flatMap((source) => {
      try {
        return source.source ? buildRecordDisplayRows(source).map((row) => ({ ...row, drawingId: source.drawingId,
          cropped: source.cropped, reverse: source.reverse, paired: source.paired })) : [];
      } catch { return []; } // Existing selector control owns its visible validation error.
    }));
  // A drawing's drafts on the rows of their sources (OV-278). An edit of a
  // draft first writes them there, within its own step.
  /** @param {DrawingState} drawing */
  const boundDrafts = (drawing) => bindRecordDisplayDrafts(drawing.recordDisplayDrafts, allRows.value);
  /** @param {DrawingState} drawing */
  const settleDrafts = (drawing) => {
    const bound = boundDrafts(drawing);
    if (bound !== drawing.recordDisplayDrafts) drawing.recordDisplayDrafts.splice(0, bound.length, ...bound);
  };
  /** @param {DrawnRecordRow} row */
  const draftFor = (row) => boundDrafts(state.drawingOf(row.drawingId)).find((/** @type {RecordDisplayDraft} */ draft) => recordDisplayKey(draft) === row.key) || {};
  // Linear orientation has one owner, the File card's region_reverse (CO-03,
  // N-09). A row that is the only record its card draws writes it there; the
  // row override remains only for Circular and for a Linear card that draws
  // several records.
  /** @param {DrawingState} drawing the row's drawing @param {{ scope: string, sourceUid: string } | null | undefined} row */
  const linearOrientationOwner = (drawing, row) => {
    if (row?.scope !== 'linear') return null;
    const sequence = drawing.linearSeqs.find((/** @type {LinearSeq} */ seq) => seq.uid === row.sourceUid);
    return sequence && rows.value.filter((entry) => entry.scope === 'linear'
      && entry.sourceUid === row.sourceUid).length === 1 ? sequence : null;
  };
  /** @param {DrawingState} drawing @param {RecordDisplayRow} row @param {RecordDisplayDraft} draft @param {boolean | null} reverseComplement */
  const writeReverseComplement = (drawing, row, draft, reverseComplement) => {
    const owner = reverseComplement === null ? null : linearOrientationOwner(drawing, row);
    if (owner && reverseComplement !== null) owner.region_reverse = reverseComplement;
    draft.reverseComplementOverride = owner ? null : reverseComplement;
  };
  /** @param {DrawnRecordRow} row */
  const surfaceFor = (row) => {
    const surface = recordDisplaySurface(row, draftFor(row), row);
    const busy = state.sessionOperationAvailability?.();
    return busy ? { ...surface, startEnabled: false, disabledReason: busy.reason } : surface;
  };
  /** @param {DrawingState} drawing @param {RecordDisplayRow} row @param {Partial<RecordDisplayDraft>} patch @param {string} label */
  const edit = (drawing, row, patch, label) => runUndoable(label, () => {
    settleDrafts(drawing);
    let draft = drawing.recordDisplayDrafts.find((/** @type {RecordDisplayDraft} */ entry) => recordDisplayKey(entry) === row.key);
    if (!draft) {
      draft = { sourceUid: row.sourceUid, selector: row.selector,
        recordId: row.recordId, topologyOverride: null, startCoordinate: null,
        reverseComplementOverride: null, anchorIntent: null };
      drawing.recordDisplayDrafts.push(draft);
      draft = drawing.recordDisplayDrafts[drawing.recordDisplayDrafts.length - 1];
    }
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
    const busy = state.sessionOperationAvailability?.();
    if (busy) return busy;
    if (typeof reverseComplement !== 'boolean') {
      throw new Error('Resolved record orientation must be boolean.');
    }
    const resolvedStart = parseRecordDisplayStart(startCoordinate);
    const resolvedIntent = validateAnchorIntent(anchorIntent);
    settleDrafts(drawing);
    let draft = drawing.recordDisplayDrafts.find(
      (/** @type {RecordDisplayDraft} */ entry) => recordDisplayKey(entry) === row.key
    );
    if (!draft) {
      draft = {
        sourceUid: row.sourceUid,
        selector: row.selector,
        recordId: row.recordId,
        topologyOverride: null,
        startCoordinate: null,
        reverseComplementOverride: null,
        anchorIntent: null
      };
      drawing.recordDisplayDrafts.push(draft);
      draft = drawing.recordDisplayDrafts[drawing.recordDisplayDrafts.length - 1];
    }
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
      const inputs = sources.value.filter(source => source.scope === 'linear'
        && source.sourceUid === orientation.recordKey);
      // The card of the one input, in the drawing's rows.
      const sequence = inputs.length === 1
        ? drawing.linearSeqs.find((card) => card.uid === inputs[0].sourceUid) : undefined;
      if (records.length !== 1 || !sequence) {
        throw new Error('Alignment record transform target is stale.');
      }
      const record = records[0];
      const input = inputs[0];
      if (!boundSources.some(source => source.scope === input.scope && source.sourceUid === input.sourceUid
          && source.source === input.source && source.paired === input.paired)
        || !input.source || !matchesSavedSource(input.source, record.source.resourceId || record.source.gffResourceId)
        || (record.source.kind === 'gffFasta'
          && (!input.paired || !matchesSavedSource(input.paired, record.source.fastaResourceId)))) {
        throw new Error('Alignment record source binding changed.');
      }
      const selector = record.region?.selector || record.selector;
      const drafts = drawing.recordDisplayDrafts.filter((/** @type {RecordDisplayDraft} */ draft) => draft.sourceUid === input.sourceUid
        && (selector?.kind === 'recordIndex'
          ? draft.selector === `#${selector.index + 1}`
          : selector?.kind === 'recordId' ? draft.recordId === selector.value : draft.selector === '#1'));
      return { recordKey: orientation.recordKey, orientation, drafts, sequence };
    });
  };
  /** @param {DrawingState} drawing the Linear drawing the alignment run pinned @param {AlignmentOrientation[]} orientations @returns {AlignmentOrientationCheckpoint[]} */
  const captureAlignmentOrientationIntent = (drawing, orientations) => {
    // A draft names no mode (PD-OI-086): these are the Linear drawing's (OV-579).
    // On their rows first (OV-278): a draft that the capture's own settle moved would not be found.
    settleDrafts(drawing);
    return alignmentRows(drawing, orientations)
      .map(({recordKey, sequence, drafts}) => ({ recordKey, reverseComplement: Boolean(sequence.region_reverse),
        drafts: drafts.map((/** @type {RecordDisplayDraft} */ draft) => captureTargetDraft({ ...draft, scope: 'linear', drawingId: drawing.id })) }));
  };
  /** @param {DrawingState} drawing @param {AlignmentOrientationCheckpoint[]} checkpoints */
  const restoreAlignmentOrientationIntent = (drawing, checkpoints) => {
    alignmentRows(drawing, checkpoints).forEach(({sequence}, index) => {
      sequence.region_reverse = checkpoints[index].reverseComplement;
      checkpoints[index].drafts.forEach(restoreTargetDraft);
    });
  };
  /** @param {DrawingState} drawing @param {AlignmentOrientation[]} orientations */
  const commitAlignmentOrientations = (drawing, orientations) => {
    alignmentRows(drawing, orientations).forEach(({sequence, orientation, drafts}) => {
      sequence.region_reverse = orientation.reverseComplement;
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
      boundSources = sources.value.map(({ scope, sourceUid, source, paired }) => ({ scope, sourceUid, source, paired }));
    }
    committedRows = allRows.value.filter((row) => row.scope === request.mode
      && boundSources.some((source) => source.scope === row.scope && source.sourceUid === row.sourceUid
        && source.source === row.source && source.paired === row.paired)).flatMap((row) => {
      const selected = request.records.find((/** @type {CommittedRecord} */ record, /** @type {number} */ index) => {
        const selector = record.region?.selector || record.selector;
        const matchesSelector = !selector || (selector.kind === 'recordIndex'
          ? row.selector === `#${selector.index + 1}` : row.recordId === selector.value);
        return matchesSelector && (row.scope === 'circular'
          ? request.records.length === 1 || index === Number(row.selector.slice(1)) - 1
          : record.recordKey === row.sourceUid || record.recordKey === `${row.sourceUid}:${Number(row.selector.slice(1))}`);
      });
      if (!selected) return [];
      const source = selected.source;
      if (!matchesSavedSource(row.source, source.resourceId || source.gffResourceId)
        || (row.paired && !matchesSavedSource(row.paired, source.fastaResourceId))) return [];
      return [{ ...row,
        committedDisplay: selected.display || { isCircular: null, startCoordinate: null },
        committedReverseComplement: selected.region
          ? Boolean(selected.region.reverseComplement)
          : Boolean(selected.presentation?.reverseComplement),
        committedCropped: Boolean(selected.region),
        canonicalRecordKey: selected.recordKey,
        canonicalSource: selected.source,
        recordKey: selected.cardinality === 'all'
        && allRows.value.filter((entry) => entry.sourceUid === row.sourceUid).length > 1
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
    if (source.scope === 'circular') return circularDiscoveryForInput(state.cInputType.value, state.drawingOf(source.drawingId));
    const sequence = state.drawingOf(source.drawingId).linearSeqs.find((/** @type {LinearSeq} */ seq) => seq.uid === source.sourceUid);
    return { status: linearRecordStatusFor(sequence), error: linearRecordErrorFor(sequence) };
  };
  // The committed record's current sources have no read records: say whether
  // the read failed or has not happened, instead of calling the target stale.
  /** @param {string} recordKey */
  const unreadRecordsError = (recordKey) => {
    // The committed request is the adopted canonical request: its `records` is always set.
    const request = /** @type {CommittedRecordRequest | null} */ (getCommittedRequest());
    const record = committedRecordForKey(request, recordKey);
    const discoveries = record ? sources.value
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
  /** @param {RecordDraftTarget} row @returns {TargetDraftCheckpoint} */
  const captureTargetDraft = (row) => {
    const drawing = state.drawingOf(row.drawingId);
    settleDrafts(drawing);
    const key = recordDisplayKey(row);
    const index = drawing.recordDisplayDrafts.findIndex(
      (/** @type {RecordDisplayDraft} */ draft) => recordDisplayKey(draft) === key
    );
    const owner = linearOrientationOwner(drawing, row);
    return {
      scope: row.scope,
      drawingId: row.drawingId,
      key,
      index,
      draft: index >= 0 ? cloneJsonData(drawing.recordDisplayDrafts[index]) : null,
      ...(owner ? { ownerUid: owner.uid, ownerReverse: Boolean(owner.region_reverse) } : {})
    };
  };
  /** @param {TargetDraftCheckpoint} checkpoint */
  const restoreTargetDraft = (checkpoint) => {
    const drawing = state.drawingOf(checkpoint.drawingId);
    if (checkpoint.ownerUid) {
      const owner = drawing.linearSeqs.find((/** @type {LinearSeq} */ seq) => seq.uid === checkpoint.ownerUid);
      if (owner) owner.region_reverse = checkpoint.ownerReverse;
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
  return { rows, allRows, draftFor, surfaceFor, hasPendingChanges, changeRecordSources,
    // The data of isCurrentFeature (services/feature-identity.js), refreshed
    // for the committed request before it is read (R13).
    sourceBinding: () => {
      refreshCommittedRows();
      // The committed request is the adopted canonical request: its `records` is always set.
      return { request: /** @type {CommittedRecordRequest | null} */ (getCommittedRequest()),
        sources: sources.value, boundSources, matchesSavedSource };
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
      return state.sessionOperationAvailability?.() || edit(drawing, row, {
        startCoordinate: parseRecordDisplayStart(value), anchorIntent: null
      }, msg('Change record display start'));
    },
    resetStart: (/** @type {RecordDisplayRow} */ row) => {
      const drawing = state.activeDrawing();
      return edit(drawing, row, { startCoordinate: null, anchorIntent: null }, msg('Reset record display start'));
    },
    setReverseComplement: (/** @type {RecordDisplayRow} */ row, /** @type {unknown} */ value) => {
      const drawing = state.activeDrawing();
      return state.sessionOperationAvailability?.() || edit(drawing, row, {
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
    captureTargetDraft,
    restoreTargetDraft,
    targetForFeature };
};
