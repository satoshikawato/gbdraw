// Source-bound editable rotation intent. The request service owns serialization.
import { resolveDisambiguatedRecordSelection } from './record-options.js';
import { matchesSessionResourceDescriptor } from '../services/session-resource-backing.js';

export const recordDisplayKey = ({ scope, sourceUid, selector }) => {
  if (!['circular', 'linear'].includes(scope) || !sourceUid || !/^#[1-9]\d*$/.test(selector)) {
    throw new Error('Record display requires a scope, source UID, and exact record selector.');
  }
  return JSON.stringify([scope, sourceUid, selector]);
};

export const buildRecordDisplayRows = ({ scope, sourceUid, source, records, selector = '' }) => {
  const selection = resolveDisambiguatedRecordSelection(records, selector);
  if (!['unspecified', 'resolved'].includes(selection.status)) {
    throw new Error(`Record selector is ${selection.status}: ${selector}.`);
  }
  return (selection.record ? [selection.record] : selection.entries).map((record) => {
    const row = { ...record, scope, sourceUid, source };
    return { ...row, key: recordDisplayKey(row) };
  });
};

export const reconcileRecordDisplayDrafts = (drafts, discoveredRows, replacedSourceUids = []) => {
  // Pass all discovered rows, including inactive selectors. A source replacement
  // purges its drafts even when the existing source card keeps its UID.
  const keys = new Set(discoveredRows.map(recordDisplayKey));
  return drafts.filter((draft) => !replacedSourceUids.includes(draft.sourceUid)
    && keys.has(recordDisplayKey(draft)));
};

export const parseRecordDisplayStart = (value) => {
  if (value === null || (typeof value === 'string' && !value.trim())) return null;
  if (!['string', 'number'].includes(typeof value)) throw new Error('Display start must be an integer.');
  const number = Number(value);
  if (!Number.isSafeInteger(number) || number < 1) {
    throw new Error('Display start must be a positive integer.');
  }
  return number;
};

export const validateRecordDisplayDrafts = (drafts) => {
  if (!Array.isArray(drafts)) throw new Error('Record display drafts must be an array.');
  const identities = new Set();
  for (const draft of drafts) {
    if (!draft || Object.keys(draft).sort().join(',') !== 'recordId,scope,selector,sourceUid,startCoordinate,topologyOverride'
      || typeof draft.recordId !== 'string' || draft.recordId.includes('\0')
      || typeof draft.sourceUid !== 'string' || draft.sourceUid.includes('\0')
      || (draft.topologyOverride !== null && typeof draft.topologyOverride !== 'boolean')
      || (draft.startCoordinate !== null && (!Number.isSafeInteger(draft.startCoordinate) || draft.startCoordinate < 1))) {
      throw new Error('Invalid record display draft; only source-bound requested intent is supported.');
    }
    const key = recordDisplayKey(draft);
    if (identities.has(key)) throw new Error('Duplicate record display draft identity.');
    identities.add(key);
  }
  return drafts;
};

export const recordDisplaySurface = (row, draft = {}, { cropped = false, reverse = false } = {}) => {
  const override = draft.topologyOverride ?? null;
  if (override !== null && typeof override !== 'boolean') throw new Error('Topology override must be boolean or null.');
  const effectiveCircular = override ?? row.detectedTopology === 'circular';
  const lengthKnown = Number.isSafeInteger(row.recordLength) && row.recordLength > 0;
  return {
    effectiveCircular,
    startEnabled: effectiveCircular && lengthKnown && !cropped,
    disabledReason: cropped ? 'Display start is unavailable for a cropped record.'
      : !lengthKnown ? 'Record length is unavailable.'
        : !effectiveCircular ? 'Display start requires a circular record.' : '',
    currentStart: reverse ? (lengthKnown ? row.recordLength : null) : 1
  };
};

export const requestedRecordDisplay = (row, draft = {}, context = {}) => {
  const surface = recordDisplaySurface(row, draft, context);
  const start = parseRecordDisplayStart(draft.startCoordinate ?? null);
  if (surface.startEnabled && start !== null && start > row.recordLength) {
    throw new Error(`Display start for ${row.recordId} must be between 1 and ${row.recordLength}.`);
  }
  return { isCircular: draft.topologyOverride ?? null,
    startCoordinate: surface.startEnabled ? start : null };
};

export const createRecordDisplayControls = ({ state, computed, watch, linearRecordSelector, history, getCommittedRequest, getCommittedSession }) => {
  let committedRows = [];
  let boundRequest = null;
  let boundSources = [];
  const sources = computed(() => [
    { scope: 'circular', sourceUid: 'circular',
      source: state.cInputType.value === 'gff' ? state.files.c_gff : state.files.c_gb,
      paired: state.cInputType.value === 'gff' ? state.files.c_fasta : null,
      records: state.circularRecordList.value.map((record) => ({ ...record,
        recordId: record.record_id ?? record.recordId, recordLength: record.record_length ?? record.recordLength })), selector: state.form.circular_record_selector,
      cropped: state.form.circular_region_start != null || state.form.circular_region_end != null,
      reverse: state.form.circular_reverse },
    ...state.linearSeqs.map((seq) => ({ scope: 'linear', sourceUid: seq.uid,
      source: state.lInputType.value === 'gff' ? seq.gff : seq.gb,
      paired: state.lInputType.value === 'gff' ? seq.fasta : null,
      records: linearRecordSelector.recordsFor(seq), selector: seq.region_record_id,
      cropped: Boolean(seq.region_start || seq.region_end), reverse: seq.region_reverse }))
  ]);
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
  const draftFor = (row) => state.recordDisplayDrafts.find((draft) => recordDisplayKey(draft) === row.key) || {};
  const surfaceFor = (row) => recordDisplaySurface(row, draftFor(row), row);
  const edit = (row, patch, label) => history.runUndoable(label, () => {
    let draft = state.recordDisplayDrafts.find((entry) => recordDisplayKey(entry) === row.key);
    if (!draft) {
      draft = { scope: row.scope, sourceUid: row.sourceUid, selector: row.selector,
        recordId: row.recordId, topologyOverride: null, startCoordinate: null };
      state.recordDisplayDrafts.push(draft);
      draft = state.recordDisplayDrafts[state.recordDisplayDrafts.length - 1];
    }
    Object.assign(draft, patch);
  });
  const matchesSavedSource = (file, resourceId) => {
    const expected = getCommittedSession()?.resources?.[resourceId];
    return matchesSessionResourceDescriptor(file, expected);
  };
  const recordForKey = (key) => getCommittedRequest()?.records.find((record) => record.recordKey === key
    || (record.cardinality === 'all' && key.startsWith(`${record.recordKey}:`)
      && /^[1-9]\d*$/.test(key.slice(record.recordKey.length + 1))));
  const recordUsesSource = (record, source) => source.scope === getCommittedRequest()?.mode
    && (source.scope === 'circular' || record.recordKey === source.sourceUid
      || (record.selector?.kind === 'recordIndex'
        && record.recordKey === `${source.sourceUid}:${record.selector.index + 1}`));
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
      return [{ ...row, committedDisplay: selected.display || { isCircular: null, startCoordinate: null }, recordKey: selected.cardinality === 'all'
        && allRows.value.filter((entry) => entry.sourceUid === row.sourceUid).length > 1
        ? `${selected.recordKey}:${Number(row.selector.slice(1))}` : selected.recordKey }];
    });
  };
  watch([() => state.featureCatalog.value, () => state.processing?.value, allRows], refreshCommittedRows, { flush: 'post' });
  const shortcutState = (row, shortcut) => {
    refreshCommittedRows();
    if (!surfaceFor(row).startEnabled) return { enabled: false, reason: surfaceFor(row).disabledReason };
    try {
      const start = selectedFeatureDisplayStart({ row,
        committedRow: committedRows.find((entry) => entry.key === row.key && entry.paired === row.paired),
        selectedFeatures: state.selectedFeatures.value, shortcut });
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
      Object.entries(state.featurePlacementOverrides).forEach(([key, row]) => {
        const record = recordForKey(row.recordKey);
        if (record && replaced.some((source) => recordUsesSource(record, source))) {
          delete state.featurePlacementOverrides[key];
        }
      });
      for (let index = state.recordDisplayDrafts.length - 1; index >= 0; index -= 1) {
        if (replaced.some((source) => source.scope === state.recordDisplayDrafts[index].scope
          && source.sourceUid === state.recordDisplayDrafts[index].sourceUid)) state.recordDisplayDrafts.splice(index, 1);
      }
    }, { flush: 'sync' });
  const hasPendingChanges = computed(() => {
    void state.processing?.value;
    refreshCommittedRows();
    const request = getCommittedRequest();
    if (!request) return false;
    const options = request.diagramOptions || {};
    const tolerance = options.configOverrides?.['canvas.feature_overlap_tolerance_bp']
      ?? options.config?.canvas?.feature_overlap_tolerance_bp ?? 0;
    const placements = options.featurePlacements || [];
    if (request.mode !== state.mode.value || tolerance !== state.adv.feature_overlap_tolerance_bp
      || placements.length !== Object.keys(state.featurePlacementOverrides).length
      || placements.some((row) => {
        const target = state.featurePlacementOverrides[JSON.stringify([row.recordKey, row.biologicalFeatureId])]?.placement;
        return !target || ['kind', 'side', 'level'].some((field) => target[field] !== row.placement[field]);
      })) return true;
    return committedRows.some((row) => {
      try { return JSON.stringify(requestedRecordDisplay(row, draftFor(row), row)) !== JSON.stringify(row.committedDisplay); }
      catch { return true; }
    });
  });
  return { rows, allRows, draftFor, surfaceFor, shortcutState, hasPendingChanges,
    applyShortcut: (row, shortcut) => {
      const result = shortcutState(row, shortcut);
      if (!result.enabled) throw new Error(result.reason);
      return edit(row, { startCoordinate: result.start }, 'Use selected feature for record display start');
    },
    isCurrentFeature: (feature) => {
      refreshCommittedRows();
      const record = recordForKey(feature.record_key);
      if (!record) return false;
      return sources.value.some((source) => source.source && recordUsesSource(record, source)
        && boundSources.some((bound) => bound.scope === source.scope && bound.sourceUid === source.sourceUid
          && bound.source === source.source && bound.paired === source.paired)
        && matchesSavedSource(source.source, record.source.resourceId || record.source.gffResourceId)
        && (!source.paired || matchesSavedSource(source.paired, record.source.fastaResourceId)));
    },
    rowsFor: (uid) => rows.value.filter((row) => row.sourceUid === uid),
    setTopology: (row, value) => edit(row, { topologyOverride: value }, 'Change record topology'),
    setStart: (row, value) => edit(row, { startCoordinate: parseRecordDisplayStart(value) }, 'Change record display start'),
    resetStart: (row) => edit(row, { startCoordinate: null }, 'Reset record display start') };
};

export const selectedFeatureDisplayStart = ({ row, committedRow, selectedFeatures, shortcut }) => {
  if (!committedRow || recordDisplayKey(row) !== recordDisplayKey(committedRow)
    || !row.source || row.source !== committedRow.source
    || row.recordLength !== committedRow.recordLength || !committedRow.recordKey
    || selectedFeatures.length !== 1) {
    throw new Error('Select one feature bound to the current source record.');
  }
  const feature = selectedFeatures[0];
  if (feature?.record_key !== committedRow.recordKey) throw new Error('Selected feature belongs to another record.');
  const parts = feature.location_parts;
  const strand = feature.strand;
  if (!['+', '-'].includes(strand) || !Array.isArray(parts) || !parts.length
    || !Number.isSafeInteger(row.recordLength) || row.recordLength <= 0
    || parts.some((part) => part.strand !== strand
      || !Number.isSafeInteger(part.start) || !Number.isSafeInteger(part.end)
      || part.start < 0 || part.end <= part.start || part.end > row.recordLength)) {
    throw new Error('Feature shortcut requires nonempty source parts with one known strand.');
  }
  if (!['five-prime', 'midpoint'].includes(shortcut)) throw new Error('Unknown feature shortcut.');
  const length = parts.reduce((sum, part) => sum + part.end - part.start, 0);
  let offset = shortcut === 'five-prime' ? 0 : Math.floor((length - 1) / 2);
  for (const part of parts) {
    if (offset < part.end - part.start) {
      return strand === '+' ? part.start + 1 + offset : part.end - offset;
    }
    offset -= part.end - part.start;
  }
};
