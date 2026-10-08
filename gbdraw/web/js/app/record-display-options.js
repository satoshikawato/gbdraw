// @ts-check
/** @import { DrawingState } from '../state.js' */
// Source-bound editable rotation intent. The request service owns serialization.
import {
  buildRecordDisplayRows,
  parseRecordDisplayStart,
  recordDisplayKey,
  recordDisplaySurface,
  requestedRecordTransform,
  requireReverseComplementOverride,
  validateAnchorIntent
} from '../services/record-display-model.js';
import { resolveFeatureAnchor } from './record-display/feature-anchor.js';
import { matchesSessionResourceDescriptor } from '../services/session-resource-backing.js';
import { cloneJsonData } from '../services/json-clone.js';
import { featureIdentityKey, requestFeaturePlacements } from '../services/feature-placement.js';
import { committedRecordForKey, committedRecordUsesSource } from '../services/feature-identity.js';
import { circularDiscoveryForInput } from '../services/record-discovery.js';
import { RECORD_READ_ERROR_LABEL } from './linear-record-selector.js';

// Session Load reads no record bytes (776a2f93), so a loaded Result can have a
// current source whose records are not read yet. That popup target failure is
// not stale: its kind lets an explicit record action read them.
export const RECORD_TARGET_NOT_DISCOVERED = 'not-discovered';
const recordTargetError = (kind, message) => Object.assign(new Error(message), { kind });

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
  let committedRows = [];
  /** @type {Record<string, any> | null} */
  let boundRequest = null;
  let boundSources = [];
  const sources = computed(() => {
    const drawing = state.activeDrawing();
    return [
      { scope: 'circular', sourceUid: 'circular',
        source: state.cInputType.value === 'gff' ? state.files.c_gff : state.files.c_gb,
        paired: state.cInputType.value === 'gff' ? state.files.c_fasta : null,
        records: state.circularRecordList.value.map((record) => ({ ...record,
          recordId: record.record_id ?? record.recordId, recordLength: record.record_length ?? record.recordLength })), selector: drawing.form.circular_record_selector,
        cropped: drawing.form.circular_region_start != null || drawing.form.circular_region_end != null,
        reverse: drawing.form.circular_reverse },
      ...state.linearSeqs.map((seq) => ({ scope: 'linear', sourceUid: seq.uid,
        source: state.lInputType.value === 'gff' ? seq.gff : seq.gb,
        paired: state.lInputType.value === 'gff' ? seq.fasta : null,
        records: linearRecordsFor(seq), selector: seq.region_record_id,
        cropped: Boolean(seq.region_start || seq.region_end), reverse: seq.region_reverse }))
    ];
  });
  const allRows = computed(() => sources.value.flatMap((source) => source.source
    ? buildRecordDisplayRows({ ...source, selector: '' }).map((row) => ({ ...row,
        cropped: source.cropped, reverse: source.reverse, paired: source.paired })) : []));
  const rows = computed(() => sources.value.filter((source) => source.scope === state.mode.value)
    .flatMap((source) => {
      try {
        return source.source ? buildRecordDisplayRows(source).map((row) => ({ ...row,
          cropped: source.cropped, reverse: source.reverse, paired: source.paired })) : [];
      } catch { return []; } // Existing selector control owns its visible validation error.
    }));
  // A row's drafts are those of its mode's drawing.
  /** @param {Record<string, any>} row @returns {DrawingState} */
  const drawingOfRow = (row) => (row?.scope === 'linear' ? state.drawings.linear : state.drawings.circular);
  const draftFor = (row) => drawingOfRow(row).recordDisplayDrafts.find((draft) => recordDisplayKey(draft) === row.key) || {};
  // Linear orientation has one owner, the File card's region_reverse (CO-03,
  // N-09). A row that is the only record its card draws writes it there; the
  // row override remains only for Circular and for a Linear card that draws
  // several records.
  const linearOrientationOwner = (row) => {
    if (row?.scope !== 'linear') return null;
    const sequence = state.linearSeqs.find((seq) => seq.uid === row.sourceUid);
    return sequence && rows.value.filter((entry) => entry.scope === 'linear'
      && entry.sourceUid === row.sourceUid).length === 1 ? sequence : null;
  };
  const writeReverseComplement = (row, draft, reverseComplement) => {
    const owner = reverseComplement === null ? null : linearOrientationOwner(row);
    if (owner) owner.region_reverse = reverseComplement;
    draft.reverseComplementOverride = owner ? null : reverseComplement;
  };
  const surfaceFor = (row) => {
    const surface = recordDisplaySurface(row, draftFor(row), row);
    const busy = state.sessionOperationAvailability?.();
    return busy ? { ...surface, startEnabled: false, disabledReason: busy.reason } : surface;
  };
  /** @param {DrawingState} drawing */
  const edit = (drawing, row, patch, label) => runUndoable(label, () => {
    let draft = drawing.recordDisplayDrafts.find((entry) => recordDisplayKey(entry) === row.key);
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
      writeReverseComplement(row, draft, reverseComplementOverride);
    }
  });
  const writeResolvedTransform = (row, { startCoordinate, reverseComplement, anchorIntent }) => {
    const drawing = drawingOfRow(row);
    const busy = state.sessionOperationAvailability?.();
    if (busy) return busy;
    if (typeof reverseComplement !== 'boolean') {
      throw new Error('Resolved record orientation must be boolean.');
    }
    const resolvedStart = parseRecordDisplayStart(startCoordinate);
    const resolvedIntent = validateAnchorIntent(anchorIntent);
    let draft = drawing.recordDisplayDrafts.find(
      (entry) => recordDisplayKey(entry) === row.key
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
    writeReverseComplement(row, draft, reverseComplement);
  };
  /** @param {DrawingState} drawing */
  const alignmentRows = (drawing, orientations) => {
    refreshCommittedRows();
    const request = getCommittedRequest();
    return orientations.map(orientation => {
      const records = request?.mode === 'linear'
        ? request.records.filter(record => record.recordKey === orientation.recordKey) : [];
      const inputs = sources.value.filter(source => source.scope === 'linear'
        && source.sourceUid === orientation.recordKey);
      if (records.length !== 1 || inputs.length !== 1) {
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
      const drafts = drawing.recordDisplayDrafts.filter(draft => draft.sourceUid === input.sourceUid
        && (selector?.kind === 'recordIndex'
          ? draft.selector === `#${selector.index + 1}`
          : selector?.kind === 'recordId' ? draft.recordId === selector.value : draft.selector === '#1'));
      return { recordKey: orientation.recordKey, orientation, drafts,
        sequence: state.linearSeqs.find(sequence => sequence.uid === input.sourceUid) };
    });
  };
  const captureAlignmentOrientationIntent = (orientations) => {
    const drawing = state.drawings.linear;
    return alignmentRows(drawing, orientations)
      .map(({recordKey, sequence, drafts}) => ({ recordKey,
        reverseComplement: Boolean(sequence.region_reverse), drafts: drafts.map(captureTargetDraft) }));
  };
  const restoreAlignmentOrientationIntent = (checkpoints) => {
    const drawing = state.drawings.linear;
    alignmentRows(drawing, checkpoints).forEach(({sequence}, index) => {
      sequence.region_reverse = checkpoints[index].reverseComplement;
      checkpoints[index].drafts.forEach(restoreTargetDraft);
    });
  };
  const commitAlignmentOrientations = (orientations) => {
    const drawing = state.drawings.linear;
    alignmentRows(drawing, orientations).forEach(({sequence, orientation, drafts}) => {
      sequence.region_reverse = orientation.reverseComplement;
      drafts.forEach(draft => { draft.reverseComplementOverride = null; });
    });
  };
  const matchesSavedSource = (file, resourceId) => {
    const expected = getCommittedSession()?.resources?.[resourceId];
    return matchesSessionResourceDescriptor(file, expected);
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
      const selected = request.records.find((record, index) => {
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
  const shortcutState = (row, shortcut) => {
    refreshCommittedRows();
    if (!surfaceFor(row).startEnabled) return { enabled: false, reason: surfaceFor(row).disabledReason };
    try {
      const selectedFeatures = state.selectedFeatures.value;
      if (selectedFeatures.length !== 1) {
        throw new Error('Select one feature bound to the current source record.');
      }
      const start = selectedFeatureDisplayStart({ row,
        committedRow: committedRows.find((entry) => entry.key === row.key && entry.paired === row.paired),
        feature: selectedFeatures[0], shortcut });
      return { enabled: true, start, reason: '' };
    } catch (error) { return { enabled: false, reason: error.message }; }
  };
  watch(() => sources.value.map(({ scope, sourceUid, source, paired }) => ({ scope, sourceUid, source, paired })),
    (current, previous = []) => {
      if (state.semanticFileWatchersSuppressed?.value || state.sessionImportRollbackInProgress?.value) return;
      const replaced = previous.filter((old) => {
        const next = current.find((entry) => entry.scope === old.scope && entry.sourceUid === old.sourceUid);
        return old.source && (!next || next.source !== old.source || next.paired !== old.paired);
      });
      // A replaced source's drafts go from its mode's drawing.
      replaced.forEach((source) => {
        const drafts = drawingOfRow(source).recordDisplayDrafts;
        for (let index = drafts.length - 1; index >= 0; index -= 1) {
          if (drafts[index].sourceUid === source.sourceUid) drafts.splice(index, 1);
        }
      });
    }, { flush: 'sync' });
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
    if (request.mode !== state.mode.value || tolerance !== drawing.adv.feature_overlap_tolerance_bp
      || placements.length !== draftPlacementCount()
      || placements.some((row) => {
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
  const discoveryFor = (source) => {
    if (source.scope === 'circular') return circularDiscoveryForInput(state);
    const sequence = state.linearSeqs.find((seq) => seq.uid === source.sourceUid);
    return { status: linearRecordStatusFor(sequence), error: linearRecordErrorFor(sequence) };
  };
  // The committed record's current sources have no read records: say whether
  // the read failed or has not happened, instead of calling the target stale.
  const unreadRecordsError = (recordKey) => {
    const request = getCommittedRequest();
    const record = committedRecordForKey(request, recordKey);
    const discoveries = record ? sources.value
      .filter((source) => source.source && committedRecordUsesSource(request, record, source)).map(discoveryFor) : [];
    if (!discoveries.length || discoveries.some(({ status }) => status === 'ready')) return null;
    const failed = discoveries.find(({ status }) => status === 'error');
    return failed
      ? recordTargetError('discovery-failed',
        (typeof failed.error === 'string' ? failed.error : failed.error?.summary) || RECORD_READ_ERROR_LABEL)
      : recordTargetError(RECORD_TARGET_NOT_DISCOVERED, 'The records of this feature have not been read yet.');
  };
  // The display start and orientation that Generate would apply from this
  // committed row's draft, when either differs from the committed request.
  const pendingTransformFor = (row) => {
    let requested;
    try { requested = requestedRecordTransform(row, draftFor(row), row); } catch { return null; }
    const startCoordinate = requested.display.startCoordinate;
    return startCoordinate === (row.committedDisplay?.startCoordinate ?? null)
      && requested.reverseComplement === row.committedReverseComplement
      ? null
      : { startCoordinate, reverseComplement: requested.reverseComplement };
  };
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
  const captureTargetDraft = (row) => {
    const drawing = drawingOfRow(row);
    const key = recordDisplayKey(row);
    const index = drawing.recordDisplayDrafts.findIndex(
      (draft) => recordDisplayKey(draft) === key
    );
    const owner = linearOrientationOwner(row);
    return {
      scope: row.scope,
      key,
      index,
      draft: index >= 0 ? cloneJsonData(drawing.recordDisplayDrafts[index]) : null,
      ...(owner ? { ownerUid: owner.uid, ownerReverse: Boolean(owner.region_reverse) } : {})
    };
  };
  const restoreTargetDraft = (checkpoint) => {
    const drawing = drawingOfRow(checkpoint);
    if (checkpoint.ownerUid) {
      const owner = state.linearSeqs.find((seq) => seq.uid === checkpoint.ownerUid);
      if (owner) owner.region_reverse = checkpoint.ownerReverse;
    }
    const index = drawing.recordDisplayDrafts.findIndex(
      (draft) => recordDisplayKey(draft) === checkpoint.key
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
  return { rows, allRows, draftFor, surfaceFor, shortcutState, hasPendingChanges,
    applyShortcut: (row, shortcut) => {
      const drawing = state.activeDrawing();
      const busy = state.sessionOperationAvailability?.();
      if (busy) return busy;
      const result = shortcutState(row, shortcut);
      if (!result.enabled) throw new Error(result.reason);
      return edit(drawing, row, { startCoordinate: result.start }, 'Use selected feature for record display start');
    },
    // The data of isCurrentFeature (services/feature-identity.js), refreshed
    // for the committed request before it is read (R13).
    sourceBinding: () => {
      refreshCommittedRows();
      return { request: getCommittedRequest(), sources: sources.value, boundSources, matchesSavedSource };
    },
    rowsFor: (uid) => rows.value.filter((row) => row.sourceUid === uid),
    setTopology: (row, value) => {
      const drawing = state.activeDrawing();
      return edit(drawing, row, { topologyOverride: value }, 'Change record topology');
    },
    setStart: (row, value) => {
      const drawing = state.activeDrawing();
      return state.sessionOperationAvailability?.() || edit(drawing, row, {
        startCoordinate: parseRecordDisplayStart(value), anchorIntent: null
      }, 'Change record display start');
    },
    resetStart: (row) => {
      const drawing = state.activeDrawing();
      return edit(drawing, row, { startCoordinate: null, anchorIntent: null }, 'Reset record display start');
    },
    setReverseComplement: (row, value) => {
      const drawing = state.activeDrawing();
      return state.sessionOperationAvailability?.() || edit(drawing, row, {
        reverseComplementOverride: requireReverseComplementOverride(value), anchorIntent: null
      }, 'Change record orientation');
    },
    // Apply on Generate: one undoable draft step, like the sidebar's edits.
    setResolvedTransform: (row, transform) => runUndoable(
      'Rotate record to feature on Generate',
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

export const selectedFeatureDisplayStart = ({ row, committedRow, feature, shortcut }) => {
  if (!committedRow || recordDisplayKey(row) !== recordDisplayKey(committedRow)
    || !row.source || row.source !== committedRow.source
    || row.recordLength !== committedRow.recordLength || !committedRow.recordKey) {
    throw new Error('Select one feature bound to the current source record.');
  }
  if (feature?.record_key !== committedRow.recordKey) throw new Error('Selected feature belongs to another record.');
  const parts = feature.location_parts;
  if (!['five-prime', 'midpoint'].includes(shortcut)) throw new Error('Unknown feature shortcut.');
  const result = resolveFeatureAnchor({
    recordLength: row.recordLength,
    effectiveCircular: true,
    cropped: false,
    currentReverseComplement: false,
    identity: {
      recordKey: committedRow.recordKey,
      biologicalFeatureId: feature.biological_feature_id
    },
    parts,
    profile: feature.anchorProfile,
    intent: { placement: 'anchor', anchor: shortcut, offsetBp: 0, orientForward: false }
  });
  if (!result.eligibility.enabled) throw new Error(result.eligibility.message);
  return result.startCoordinate;
};
