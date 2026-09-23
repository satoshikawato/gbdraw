import {
  featureIdentity,
  orthogroupIdStatus
} from '../services/feature-identity.js';

const { computed, ref } = window.Vue;

const MODES = new Set(['position', 'position_and_orientation']);
const STATUSES = new Set(['reference', 'aligned', 'skipped']);
const RATIONALES = new Set([
  'reference',
  'user_selected',
  'only_usable_candidate',
  'unique_direct_rbh',
  'skipped_by_user',
  'skipped_no_candidate',
  'skipped_unmappable'
]);

const isObject = (value) => (
  value !== null && typeof value === 'object' && !Array.isArray(value)
);

const exactObject = (value, keys, path) => {
  if (!isObject(value) || Object.keys(value).sort().join('\0') !== [...keys].sort().join('\0')) {
    throw new Error(`${path} has invalid fields.`);
  }
  return value;
};

const text = (value, path) => {
  if (typeof value !== 'string' || !value.trim() || value.includes('\0')) {
    throw new Error(`${path} must be non-empty text without NUL.`);
  }
  return value.trim();
};

const nullableInteger = (value, path) => {
  if (value === null) return null;
  if (!Number.isSafeInteger(value) || value < 0) {
    throw new Error(`${path} must be a non-negative integer or null.`);
  }
  return value;
};

const sameJson = (left, right) => JSON.stringify(left) === JSON.stringify(right);

const anchorKey = (anchor) => JSON.stringify(anchor);

const strandLabel = (strand) => (
  strand === 1 ? '+' : strand === -1 ? '-' : 'unknown'
);

const rationaleLabels = Object.freeze({
  reference: 'Exact reference',
  user_selected: 'Selected by user',
  only_usable_candidate: 'Only usable candidate',
  unique_direct_rbh: 'Unique direct RBH',
  skipped_by_user: 'Skipped by user',
  skipped_no_candidate: 'No usable candidate',
  skipped_unmappable: 'Candidate center outside the displayed crop'
});

const cloneJson = (value) => {
  if (typeof structuredClone === 'function') return structuredClone(value);
  return JSON.parse(JSON.stringify(value));
};

const deepFreeze = (value) => {
  if (!value || typeof value !== 'object' || Object.isFrozen(value)) return value;
  Object.values(value).forEach(deepFreeze);
  return Object.freeze(value);
};

const candidateView = (candidate, request) => {
  const key = anchorKey(candidate.anchor);
  const member = request.members.find(({ anchor }) => anchorKey(anchor) === key);
  if (!member) throw new Error('Alignment candidate has no current member facts.');
  const evidence = request.directEdges
    .filter(({ query, subject }) => (
      sameJson(query, request.reference) && sameJson(subject, candidate.anchor)
      || sameJson(subject, request.reference) && sameJson(query, candidate.anchor)
    ))
    .map(({ edgeKind }) => edgeKind.toUpperCase());
  return {
    key,
    anchor: candidate.anchor,
    featureIdentifier: candidate.anchor.biologicalFeatureId,
    coordinates: `${Number(member.sourceStart) + 1}..${Number(member.sourceEnd)}`,
    displayedStrand: strandLabel(candidate.displayedStrand),
    representative: candidate.representative,
    role: candidate.role || 'member',
    directEvidence: evidence.length ? evidence : ['None']
  };
};

const ambiguityViews = (response, request, previous = []) => {
  const byRecord = new Map(previous.map((entry) => [entry.recordKey, entry]));
  response.records
    .filter(({ kind }) => kind === 'ambiguous')
    .forEach((ambiguity) => {
      byRecord.set(ambiguity.recordKey, {
        recordKey: ambiguity.recordKey,
        candidates: ambiguity.candidates.map((candidate) => candidateView(candidate, request))
      });
    });
  const choices = new Map(request.choices.map((choice) => [choice.recordKey, choice]));
  return request.records
    .filter(({ recordKey }) => byRecord.has(recordKey))
    .map(({ recordKey }) => {
      const view = byRecord.get(recordKey);
      const choice = choices.get(recordKey) || null;
      return {
        ...view,
        choice: choice && choice.kind === 'select'
          ? { kind: 'select', candidateKey: anchorKey(choice.anchor) }
          : choice && choice.kind === 'skip'
            ? { kind: 'skip', candidateKey: null }
            : null
      };
    });
};

const successfulSummary = (plan) => {
  const aligned = plan.records.filter(({ status }) => status === 'aligned').length;
  const explicitlySkipped = plan.records.filter(
    ({ rationale }) => rationale === 'skipped_by_user'
  ).length;
  const noCandidate = plan.records.filter(({ rationale }) => (
    rationale === 'skipped_no_candidate' || rationale === 'skipped_unmappable'
  )).length;
  const unchanged = plan.records.filter(({ status }) => status === 'reference').length;
  const reversed = plan.records.filter(
    ({ effectiveReverseComplement }) => effectiveReverseComplement !== null
  ).length;
  return deepFreeze({
    aligned,
    unchanged,
    explicitlySkipped,
    noCandidate,
    reversed,
    text: `Alignment applied: ${aligned} aligned, ${unchanged} unchanged, `
      + `${explicitlySkipped} explicitly skipped, ${noCandidate} no-candidate, `
      + `${reversed} reversed.`
  });
};

const baseReverseComplement = (record) => Boolean(
  record?.region
    ? record.region.reverseComplement
    : record?.presentation?.reverseComplement
);

const inspectActivePlan = (plan, request) => {
  if (!plan || plan.schema !== 1 || !Array.isArray(plan.records)) return null;
  const records = new Map(
    (Array.isArray(request?.records) ? request.records : [])
      .map((record) => [String(record?.recordKey || ''), record])
  );
  return deepFreeze({
    groupId: plan.groupId,
    mode: plan.mode,
    modeLabel: plan.mode === 'position_and_orientation' ? 'Align & orient' : 'Align',
    reference: {
      ...plan.reference,
      label: plan.reference.biologicalFeatureId
    },
    records: plan.records.map((decision) => {
      const baseReverse = baseReverseComplement(records.get(decision.recordKey));
      const effectiveReverse = decision.effectiveReverseComplement === null
        ? baseReverse
        : decision.effectiveReverseComplement;
      return {
        recordKey: decision.recordKey,
        anchorLabel: decision.anchor?.biologicalFeatureId || 'Skip',
        status: decision.status,
        rationale: decision.rationale,
        rationaleLabel: rationaleLabels[decision.rationale] || decision.rationale,
        reversedFromSource: Boolean(effectiveReverse)
      };
    })
  });
};

const anchorFromSource = (source, path) => {
  const identity = featureIdentity(source);
  if (
    !identity.usable
    || !identity.recordKey.supplied
    || !identity.biologicalId.supplied
  ) {
    throw new Error(`${path} does not have one canonical feature identity.`);
  }
  return {
    recordKey: identity.recordKey.value,
    biologicalFeatureId: identity.biologicalId.value,
    sourceFeatureIndex: identity.sourceIndex.supplied ? identity.sourceIndex.value : null,
    stableFeatureSvgId: identity.stableId.supplied ? identity.stableId.value : null
  };
};

const validateAnchor = (value, path) => {
  const raw = exactObject(value, [
    'recordKey',
    'biologicalFeatureId',
    'sourceFeatureIndex',
    'stableFeatureSvgId'
  ], path);
  const stableId = raw.stableFeatureSvgId === null
    ? null
    : text(raw.stableFeatureSvgId, `${path}.stableFeatureSvgId`);
  return {
    recordKey: text(raw.recordKey, `${path}.recordKey`),
    biologicalFeatureId: text(raw.biologicalFeatureId, `${path}.biologicalFeatureId`),
    sourceFeatureIndex: nullableInteger(raw.sourceFeatureIndex, `${path}.sourceFeatureIndex`),
    stableFeatureSvgId: stableId
  };
};

const validateDecision = (value, path) => {
  const raw = exactObject(value, [
    'kind',
    'recordKey',
    'status',
    'rationale',
    'anchor',
    'effectiveReverseComplement'
  ], path);
  if (raw.kind !== 'decision' || !STATUSES.has(raw.status) || !RATIONALES.has(raw.rationale)) {
    throw new Error(`${path} contains an unknown decision value.`);
  }
  const recordKey = text(raw.recordKey, `${path}.recordKey`);
  const anchor = raw.anchor === null ? null : validateAnchor(raw.anchor, `${path}.anchor`);
  if (anchor && anchor.recordKey !== recordKey) {
    throw new Error(`${path}.anchor belongs to another record.`);
  }
  const orientation = raw.effectiveReverseComplement;
  if (orientation !== null && typeof orientation !== 'boolean') {
    throw new Error(`${path}.effectiveReverseComplement must be boolean or null.`);
  }
  const alignedRationales = new Set([
    'user_selected', 'only_usable_candidate', 'unique_direct_rbh'
  ]);
  const skippedRationales = new Set([
    'skipped_by_user', 'skipped_no_candidate', 'skipped_unmappable'
  ]);
  const valid = raw.status === 'reference'
    ? anchor !== null && raw.rationale === 'reference' && orientation === null
    : raw.status === 'aligned'
      ? anchor !== null && alignedRationales.has(raw.rationale)
      : anchor === null && skippedRationales.has(raw.rationale) && orientation === null;
  if (!valid) throw new Error(`${path} contains an invalid decision combination.`);
  return {
    kind: 'decision',
    recordKey,
    status: raw.status,
    rationale: raw.rationale,
    anchor,
    effectiveReverseComplement: orientation
  };
};

const validateCandidate = (value, path, recordKey) => {
  const raw = exactObject(value, [
    'anchor', 'displayedStrand', 'hidden', 'representative', 'role'
  ], path);
  const anchor = validateAnchor(raw.anchor, `${path}.anchor`);
  if (
    anchor.recordKey !== recordKey
    || ![null, -1, 1].includes(raw.displayedStrand)
    || typeof raw.hidden !== 'boolean'
    || typeof raw.representative !== 'boolean'
    || typeof raw.role !== 'string'
    || raw.role.includes('\0')
  ) throw new Error(`${path} contains invalid candidate facts.`);
  return {
    anchor,
    displayedStrand: raw.displayedStrand,
    hidden: raw.hidden,
    representative: raw.representative,
    role: raw.role.trim()
  };
};

const validateAmbiguity = (value, path) => {
  const raw = exactObject(value, [
    'kind', 'recordKey', 'candidates', 'directRbhCandidates'
  ], path);
  if (raw.kind !== 'ambiguous' || !Array.isArray(raw.candidates) || raw.candidates.length < 2) {
    throw new Error(`${path} is not a valid ambiguity.`);
  }
  const recordKey = text(raw.recordKey, `${path}.recordKey`);
  const candidates = raw.candidates.map((candidate, index) => (
    validateCandidate(candidate, `${path}.candidates[${index}]`, recordKey)
  ));
  const candidateKeys = new Set(candidates.map(({ anchor }) => JSON.stringify(anchor)));
  if (candidateKeys.size !== candidates.length || !Array.isArray(raw.directRbhCandidates)) {
    throw new Error(`${path} contains duplicate candidates or invalid direct evidence.`);
  }
  const directRbhCandidates = raw.directRbhCandidates.map((anchor, index) => (
    validateAnchor(anchor, `${path}.directRbhCandidates[${index}]`)
  ));
  if (directRbhCandidates.some((anchor) => !candidateKeys.has(JSON.stringify(anchor)))) {
    throw new Error(`${path} direct evidence does not identify a candidate.`);
  }
  return { kind: 'ambiguous', recordKey, candidates, directRbhCandidates };
};

const validatePlan = (value, response, path) => {
  const raw = exactObject(value, ['schema', 'mode', 'groupId', 'reference', 'records'], path);
  if (raw.schema !== 1 || raw.mode !== response.mode || raw.groupId !== response.groupId) {
    throw new Error(`${path} does not match its resolution.`);
  }
  const reference = validateAnchor(raw.reference, `${path}.reference`);
  if (!sameJson(reference, response.reference) || !Array.isArray(raw.records)) {
    throw new Error(`${path}.reference or records are invalid.`);
  }
  const records = raw.records.map((decision, index) => {
    const normalized = validateDecision(
      { ...decision, kind: 'decision' },
      `${path}.records[${index}]`
    );
    const { kind: _kind, ...planDecision } = normalized;
    return planDecision;
  });
  return { schema: 1, mode: raw.mode, groupId: raw.groupId, reference, records };
};

const validateResolution = (value, request) => {
  const raw = exactObject(value, [
    'schema', 'status', 'mode', 'groupId', 'reference', 'records', 'plan'
  ], 'alignment helper response');
  if (
    raw.schema !== 1
    || !['resolved', 'ambiguous'].includes(raw.status)
    || raw.mode !== request.mode
    || raw.groupId !== request.groupId
    || !Array.isArray(raw.records)
  ) throw new Error('Alignment helper response does not match its request.');
  const reference = validateAnchor(raw.reference, 'alignment helper response.reference');
  if (!sameJson(reference, request.reference)) {
    throw new Error('Alignment helper response changed the exact reference.');
  }
  const records = raw.records.map((record, index) => (
    record?.kind === 'decision'
      ? validateDecision(record, `alignment helper response.records[${index}]`)
      : record?.kind === 'ambiguous'
        ? validateAmbiguity(record, `alignment helper response.records[${index}]`)
        : (() => { throw new Error('Alignment helper response contains an unknown record kind.'); })()
  ));
  const expectedKeys = request.records.map(({ recordKey }) => recordKey);
  if (
    records.length !== expectedKeys.length
    || records.some((record, index) => record.recordKey !== expectedKeys[index])
  ) throw new Error('Alignment helper response changed displayed-record coverage.');
  const knownAnchors = new Set(request.members.map(({ anchor }) => JSON.stringify(anchor)));
  if (records.some((record) => (
    record.kind === 'decision'
      ? record.anchor !== null && !knownAnchors.has(JSON.stringify(record.anchor))
      : record.candidates.some(({ anchor }) => !knownAnchors.has(JSON.stringify(anchor)))
  ))) throw new Error('Alignment helper response contains an unknown group member.');
  const ambiguities = records.filter(({ kind }) => kind === 'ambiguous');
  if ((raw.status === 'ambiguous') !== (ambiguities.length > 0)) {
    throw new Error('Alignment helper response status does not match its records.');
  }
  const response = {
    schema: 1,
    status: raw.status,
    mode: raw.mode,
    groupId: raw.groupId,
    reference,
    records,
    plan: null
  };
  if (raw.status === 'resolved') {
    if (raw.plan === null) throw new Error('Resolved alignment helper response has no plan.');
    response.plan = validatePlan(raw.plan, response, 'alignment helper response.plan');
    const planRecords = response.plan.records;
    const decisions = records.map(({ kind: _kind, ...decision }) => decision);
    if (!sameJson(planRecords, decisions)) {
      throw new Error('Resolved plan differs from the helper decisions.');
    }
  } else if (raw.plan !== null) {
    throw new Error('Ambiguous alignment helper response must not contain a plan.');
  }
  return deepFreeze(response);
};

const strandValue = (value, path) => {
  if (value === 1 || value === '+' || String(value).trim() === '1') return 1;
  if (value === -1 || value === '-' || String(value).trim() === '-1') return -1;
  if (value === null || value === undefined || String(value).trim() === '') return null;
  throw new Error(`${path} has an invalid source strand.`);
};

const recordLengths = (catalog) => {
  const lengths = new Map();
  (Array.isArray(catalog?.items) ? catalog.items : []).forEach((item) => {
    const keys = Array.isArray(item?.recordKeys) ? item.recordKeys : [];
    const sources = Array.isArray(item?.sequenceSources) ? item.sequenceSources : [];
    if (keys.length !== sources.length) return;
    keys.forEach((recordKey, index) => {
      const sequence = sources[index]?.sequence;
      if (typeof sequence === 'string' && sequence.length > 0) {
        lengths.set(String(recordKey), sequence.length);
      }
    });
  });
  return lengths;
};

const buildHelperRequest = ({ group, members, reference, mode, request, catalog, choices }) => {
  const groupStatus = orthogroupIdStatus(group);
  if (!groupStatus.valid || !groupStatus.supplied) {
    throw new Error('The selected Similarity Group identity is invalid.');
  }
  if (!MODES.has(mode)) throw new Error('The requested alignment mode is unsupported.');
  if (!request || request.mode !== 'linear' || !Array.isArray(request.records)) {
    throw new Error('Generate a Linear diagram before aligning a Similarity Group.');
  }
  const lengths = recordLengths(catalog);
  const records = request.records.map((record, index) => {
    const recordKey = text(record?.recordKey, `renderRequest.records[${index}].recordKey`);
    const region = record?.region === null || record?.region === undefined
      ? null
      : {
          start: Number(record.region.start),
          end: Number(record.region.end),
          reverseComplement: Boolean(record.region.reverseComplement)
        };
    if (region && (!Number.isSafeInteger(region.start) || !Number.isSafeInteger(region.end))) {
      throw new Error(`renderRequest.records[${index}].region is invalid.`);
    }
    return {
      recordKey,
      recordLength: lengths.get(recordKey) || null,
      region,
      presentation: {
        reverseComplement: Boolean(record?.presentation?.reverseComplement)
      }
    };
  });
  const displayedRecords = new Set(records.map(({ recordKey }) => recordKey));
  const normalizedMembers = members.map((member, index) => {
    const anchor = anchorFromSource(member, `Similarity Group member ${index + 1}`);
    const start = Number(member?.start);
    const end = Number(member?.end);
    if (
      !displayedRecords.has(anchor.recordKey)
      || !Number.isSafeInteger(start)
      || start < 0
      || !Number.isSafeInteger(end)
      || end < start
    ) throw new Error(`Similarity Group member ${index + 1} has invalid display facts.`);
    return {
      source: member,
      canonicalKey: `${anchor.recordKey}\0${anchor.biologicalFeatureId}`,
      payload: {
        groupId: groupStatus.value,
        anchor,
        sourceStart: start,
        sourceEnd: end,
        sourceStrand: strandValue(member?.strand, `Similarity Group member ${index + 1}`),
        identityIsUnique: true,
        hidden: Boolean(member?.hidden),
        representative: Boolean(member?.representative),
        role: String(member?.role || '')
      }
    };
  });
  const identityCounts = new Map();
  normalizedMembers.forEach(({ canonicalKey }) => {
    identityCounts.set(canonicalKey, (identityCounts.get(canonicalKey) || 0) + 1);
  });
  normalizedMembers.forEach((member) => {
    member.payload.identityIsUnique = identityCounts.get(member.canonicalKey) === 1;
  });
  const exactReference = anchorFromSource(reference, 'Alignment reference');
  const referenceKey = `${exactReference.recordKey}\0${exactReference.biologicalFeatureId}`;
  if (identityCounts.get(referenceKey) !== 1) {
    throw new Error('The clicked feature is not a unique member of this Similarity Group.');
  }
  const exactMember = normalizedMembers.find(({ canonicalKey }) => canonicalKey === referenceKey);
  if (!sameJson(exactMember.payload.anchor, exactReference)) {
    throw new Error('The reference identity conflicts with current group metadata.');
  }

  const endpointIndex = new Map();
  normalizedMembers.forEach(({ source, payload }) => {
    const recordIndex = Number(source?.recordIndex);
    [source?.proteinId, source?.sourceProteinId].forEach((proteinId) => {
      const id = String(proteinId || '').trim();
      if (!Number.isSafeInteger(recordIndex) || recordIndex < 0 || !id) return;
      const key = `${recordIndex}\0${id}`;
      const entries = endpointIndex.get(key) || [];
      if (!entries.some((anchor) => sameJson(anchor, payload.anchor))) entries.push(payload.anchor);
      endpointIndex.set(key, entries);
    });
  });
  const directEdges = (Array.isArray(group?.orthologEdges) ? group.orthologEdges : [])
    .map((edge, index) => {
      const edgeGroup = text(edge?.orthogroupId, `orthologEdges[${index}].orthogroupId`);
      if (edgeGroup !== groupStatus.value) {
        throw new Error(`orthologEdges[${index}] belongs to another Similarity Group.`);
      }
      const endpoint = (side) => {
        const recordIndex = Number(edge?.[`${side}RecordIndex`]);
        const proteinId = text(edge?.[`${side}ProteinId`], `orthologEdges[${index}].${side}ProteinId`);
        const matches = endpointIndex.get(`${recordIndex}\0${proteinId}`) || [];
        if (matches.length !== 1) {
          throw new Error(`orthologEdges[${index}] does not resolve to one current member.`);
        }
        return matches[0];
      };
      return {
        groupId: edgeGroup,
        query: endpoint('query'),
        subject: endpoint('subject'),
        edgeKind: text(edge?.edgeKind, `orthologEdges[${index}].edgeKind`)
      };
    });
  return deepFreeze({
    schema: 1,
    mode,
    groupId: groupStatus.value,
    records,
    reference: exactReference,
    members: normalizedMembers.map(({ payload }) => payload),
    directEdges,
    choices: cloneJson(choices)
  });
};

const replacementTranslations = (state, recordKeys) => {
  const current = Array.isArray(state.linearRecordTranslations?.value)
    ? state.linearRecordTranslations.value
    : [];
  const byRecord = new Map(current.map((entry) => [String(entry?.recordKey || ''), entry]));
  if (
    byRecord.size === recordKeys.length
    && recordKeys.every((recordKey) => {
      const entry = byRecord.get(recordKey);
      return entry && Number.isFinite(entry.x) && Number.isFinite(entry.y);
    })
  ) {
    return recordKeys.map((recordKey) => {
      const entry = byRecord.get(recordKey);
      return { recordKey, x: Number(entry.x), y: Number(entry.y) };
    });
  }
  return recordKeys.map((recordKey) => ({ recordKey, x: 0, y: 0 }));
};

export const createSimilarityAlignmentActions = ({
  state,
  getOrthogroupById,
  getEnrichedOrthogroupMembers,
  getCommittedRequest,
  runAnalysis,
  cancelRunAnalysis = null,
  runHelperOperation,
  resolveOperation,
  previewCandidate = null,
  clearCandidatePreview = null,
  onError = null
}) => {
  if (
    !state
    || typeof getOrthogroupById !== 'function'
    || typeof getEnrichedOrthogroupMembers !== 'function'
    || typeof getCommittedRequest !== 'function'
    || typeof runAnalysis !== 'function'
    || typeof runHelperOperation !== 'function'
    || typeof resolveOperation !== 'string'
    || !resolveOperation
  ) throw new Error('Similarity alignment actions require their current artifact owners.');

  const draft = ref(null);
  const status = ref('idle');
  const error = ref(null);
  const summary = ref(null);
  const drawerReferenceKey = ref('');
  let actionId = 0;
  let activeRequest = null;
  let activeApply = null;
  let reportedAmbiguities = [];

  const publishError = (value) => {
    const normalized = value instanceof Error ? value : new Error(String(value || 'Alignment failed.'));
    error.value = normalized;
    if (typeof onError === 'function') onError(normalized);
    return normalized;
  };

  const clearDraft = () => {
    draft.value = null;
    activeRequest = null;
    reportedAmbiguities = [];
    if (typeof clearCandidatePreview === 'function') clearCandidatePreview();
  };

  const applyPlan = async (plan, request, expectedActionId) => {
    if (expectedActionId !== actionId) return { status: 'stale' };
    status.value = 'applying';
    const recordKeys = request.records.map(({ recordKey }) => recordKey);
    const promise = runAnalysis({
      canonicalStateOverride: {
        similarityAlignmentPlan: cloneJson(plan),
        linearRecordTranslations: replacementTranslations(state, recordKeys)
      }
    });
    activeApply = promise;
    let outcome;
    try {
      outcome = await promise;
    } catch (cause) {
      outcome = { status: 'error', error: cause };
    } finally {
      if (activeApply === promise) activeApply = null;
    }
    if (expectedActionId !== actionId) return { status: 'stale' };
    if (outcome?.status === 'ok') {
      summary.value = successfulSummary(plan);
      clearDraft();
      error.value = null;
      status.value = 'idle';
      return { status: 'ok' };
    }
    clearDraft();
    status.value = 'idle';
    if (outcome?.error) publishError(outcome.error);
    return { status: outcome?.status || 'error' };
  };

  const resolveRequest = async (request, expectedActionId, { autoApply }) => {
    status.value = 'resolving';
    let response;
    try {
      const helper = await runHelperOperation(
        resolveOperation,
        { request }
      );
      if (expectedActionId !== actionId) return { status: 'stale' };
      response = validateResolution(helper?.result, request);
    } catch (cause) {
      if (expectedActionId !== actionId) return { status: 'stale' };
      clearDraft();
      status.value = 'idle';
      publishError(cause);
      return { status: 'error' };
    }
    activeRequest = request;
    error.value = null;
    reportedAmbiguities = ambiguityViews(response, request, reportedAmbiguities);
    if (response.status === 'resolved') {
      draft.value = deepFreeze({
        response,
        choices: request.choices,
        ambiguities: reportedAmbiguities
      });
      status.value = 'ready';
      return autoApply
        ? applyPlan(response.plan, request, expectedActionId)
        : { status: 'ready' };
    }
    draft.value = deepFreeze({
      response,
      choices: request.choices,
      ambiguities: reportedAmbiguities
    });
    status.value = 'ambiguous';
    return { status: 'ambiguous' };
  };

  const start = async ({ groupId, reference, mode = 'position', source }) => {
    const id = String(groupId || '').trim();
    if (!reference) {
      publishError(new Error(
        source === 'drawer'
          ? 'Select an exact reference feature before aligning this Similarity Group.'
          : 'The clicked feature does not resolve to one exact Similarity Group member.'
      ));
      return { status: 'rejected' };
    }
    const expectedActionId = ++actionId;
    if (activeApply) {
      if (typeof cancelRunAnalysis === 'function') cancelRunAnalysis();
      await activeApply;
    }
    clearDraft();
    error.value = null;
    try {
      const group = getOrthogroupById(id);
      if (!group) throw new Error('The selected Similarity Group is unavailable.');
      const request = buildHelperRequest({
        group,
        members: getEnrichedOrthogroupMembers(group),
        reference,
        mode,
        request: getCommittedRequest(),
        catalog: state.featureCatalog?.value,
        choices: []
      });
      return resolveRequest(request, expectedActionId, { autoApply: true });
    } catch (cause) {
      if (expectedActionId !== actionId) return { status: 'stale' };
      status.value = 'idle';
      publishError(cause);
      return { status: 'rejected' };
    }
  };

  const answer = async (recordKey, kind, anchor = null) => {
    if (!activeRequest || !draft.value || !['ambiguous', 'ready'].includes(status.value)) {
      return { status: 'rejected' };
    }
    const ambiguity = reportedAmbiguities.find((record) => record.recordKey === recordKey);
    if (!ambiguity) return { status: 'rejected' };
    let selectedAnchor = null;
    if (kind === 'select') {
      selectedAnchor = validateAnchor(anchor, 'Selected alignment candidate');
      if (!ambiguity.candidates.some((candidate) => sameJson(candidate.anchor, selectedAnchor))) {
        return { status: 'rejected' };
      }
    } else if (kind !== 'skip') {
      return { status: 'rejected' };
    }
    const expectedActionId = ++actionId;
    const choices = activeRequest.choices.filter((choice) => choice.recordKey !== recordKey);
    choices.push({ recordKey, kind, anchor: selectedAnchor });
    const request = deepFreeze({ ...cloneJson(activeRequest), choices });
    return resolveRequest(request, expectedActionId, { autoApply: false });
  };

  const cancel = () => {
    actionId += 1;
    if (activeApply && typeof cancelRunAnalysis === 'function') cancelRunAnalysis();
    clearDraft();
    error.value = null;
    status.value = 'idle';
    return { status: 'canceled' };
  };

  const drawerReferenceOptions = (groupId) => {
    const id = String(groupId || '').trim();
    const group = getOrthogroupById(id);
    if (!group) return [];
    const candidates = getEnrichedOrthogroupMembers(group).map((member) => {
      try {
        const anchor = anchorFromSource(member, 'Drawer reference');
        return {
          anchor,
          canonicalKey: `${anchor.recordKey}\0${anchor.biologicalFeatureId}`,
          key: anchorKey(anchor),
          label: `${anchor.recordKey} · ${anchor.biologicalFeatureId} · `
            + `${Number(member?.start) + 1}..${Number(member?.end)} `
            + `(${strandLabel(strandValue(member?.strand, 'Drawer reference strand'))})`
        };
      } catch (_error) {
        return null;
      }
    }).filter(Boolean);
    const counts = new Map();
    candidates.forEach(({ canonicalKey }) => {
      counts.set(canonicalKey, (counts.get(canonicalKey) || 0) + 1);
    });
    return candidates
      .filter(({ canonicalKey }) => counts.get(canonicalKey) === 1)
      .map(({ canonicalKey: _canonicalKey, ...candidate }) => candidate);
  };

  const selectedDrawerReference = (groupId) => (
    drawerReferenceOptions(groupId)
      .find(({ key }) => key === drawerReferenceKey.value)?.anchor || null
  );

  const drawerDisabledReason = (groupId) => {
    if (['resolving', 'applying'].includes(status.value)) {
      return 'Wait for the current alignment operation to finish.';
    }
    const options = drawerReferenceOptions(groupId);
    if (!options.length) return 'No exact reference features are available for this group.';
    if (!selectedDrawerReference(groupId)) {
      return 'Select an exact reference record and feature before aligning.';
    }
    return '';
  };

  const activePlanInspector = computed(() => {
    try {
      return inspectActivePlan(state.similarityAlignmentPlan?.value, getCommittedRequest());
    } catch (_error) {
      return null;
    }
  });

  const unresolvedCount = computed(() => (
    Array.isArray(draft.value?.ambiguities)
      ? draft.value.ambiguities.filter(({ choice }) => choice === null).length
      : 0
  ));

  const applyDisabledReason = computed(() => {
    if (status.value === 'applying') return 'Applying the alignment plan.';
    if (status.value === 'resolving') return 'Checking the current Select and Skip choices.';
    if (unresolvedCount.value > 0) {
      return `${unresolvedCount.value} ambiguous record${unresolvedCount.value === 1 ? '' : 's'} still require Select or Skip.`;
    }
    if (status.value !== 'ready' || !draft.value?.response?.plan) {
      return 'Resolve every ambiguous record before applying.';
    }
    return '';
  });

  return {
    draft,
    status,
    error,
    summary,
    drawerReferenceKey,
    activePlanInspector,
    dialogOpen: computed(() => Boolean(
      draft.value?.ambiguities?.length && status.value !== 'idle'
    )),
    unresolvedCount,
    applyDisabledReason,
    canApply: computed(() => status.value === 'ready' && Boolean(draft.value?.response?.plan)),
    startFromPopup: (options) => start({ ...options, source: 'popup' }),
    startFromDrawer: (options) => start({
      ...options,
      reference: selectedDrawerReference(options?.groupId),
      source: 'drawer'
    }),
    drawerReferenceOptions,
    setDrawerReference: (groupId, key) => {
      const option = drawerReferenceOptions(groupId).find((entry) => entry.key === key);
      drawerReferenceKey.value = option?.key || '';
      return Boolean(option);
    },
    drawerDisabledReason,
    selectCandidate: (recordKey, anchor) => answer(recordKey, 'select', anchor),
    skipRecord: (recordKey) => answer(recordKey, 'skip'),
    applyDraft: () => (
      status.value === 'ready' && draft.value?.response?.plan && activeRequest
        ? applyPlan(draft.value.response.plan, activeRequest, actionId)
        : Promise.resolve({ status: 'rejected' })
    ),
    cancel,
    previewCandidate: (anchor) => {
      if (typeof previewCandidate === 'function') previewCandidate(anchor);
    },
    clearCandidatePreview: () => {
      if (typeof clearCandidatePreview === 'function') clearCandidatePreview();
    }
  };
};
