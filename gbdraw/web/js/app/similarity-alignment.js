import { canonicalRecordReverseComplement } from './record-display-options.js';
import { validateSimilarityAlignmentResetReceipt } from '../services/session-active-config-contract.js';
import {
  featureIdentity,
  orthogroupIdStatus
} from '../services/feature-identity.js';
import { materializeRecordTranslations } from './legend-layout/composition-actions.js';
import { isInternalProteinDisplayId } from './feature-utils.js';

const { computed, ref } = window.Vue;

const REVIEW_REASONS = new Set([
  'only_usable_candidate', 'unique_direct_rbh',
  'unique_representative', 'deterministic_candidate_1'
]);
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

const sameKeyedPlan = (left, right) => {
  if (left?.schema !== 2 || right?.schema !== 2 ||
      !Array.isArray(left.records) || !Array.isArray(right.records)) return false;
  const keyed = (plan) => ({ ...plan, records: [...plan.records].sort(
    (a, b) => String(a?.recordKey || '').localeCompare(String(b?.recordKey || ''))
  ) });
  return sameJson(keyed(left), keyed(right));
};

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

const reviewReasonLabels = Object.freeze({
  only_usable_candidate: 'The only usable candidate in this record.',
  unique_direct_rbh: 'The only candidate with a direct reciprocal best hit.',
  unique_representative: 'The only representative candidate.',
  deterministic_candidate_1: 'The first candidate in Python’s stable order.'
});

const cloneJson = (value) => {
  if (value === undefined) return undefined;
  return JSON.parse(JSON.stringify(value));
};

const deepFreeze = (value) => {
  if (!value || typeof value !== 'object' || Object.isFrozen(value)) return value;
  Object.values(value).forEach(deepFreeze);
  return Object.freeze(value);
};

const displayText = (...values) => values
  .flatMap((value) => Array.isArray(value) ? value : [value])
  .map((value) => String(value ?? '').trim())
  .find((value) => value && !isInternalProteinDisplayId(value)) || '';

const featureLabel = (anchor, displayFacts, fallback = '') => {
  const facts = displayFacts.get(anchorKey(anchor)) || {};
  const feature = facts.sequenceFeature || {};
  const gene = displayText(feature.gene, facts.gene, feature.qualifiers?.gene);
  const locusTag = displayText(feature.locus_tag, feature.locusTag,
    facts.locus_tag, facts.locusTag, feature.qualifiers?.locus_tag);
  const product = displayText(feature.product, facts.product, feature.qualifiers?.product);
  return [gene, locusTag].filter(Boolean).join(' · ')
    || product || displayText(fallback) || anchor.biologicalFeatureId;
};

const candidateView = (candidate, request, displayFacts = new Map()) => {
  const key = anchorKey(candidate.anchor);
  const member = request.members.find(({ anchor }) => anchorKey(anchor) === key);
  if (!member) throw new Error('Alignment candidate has no current member facts.');
  const coordinates = `${(candidate.sourceStart + 1).toLocaleString('en-US')}..${candidate.sourceEnd.toLocaleString('en-US')} bp`;
  return {
    key,
    anchor: candidate.anchor,
    label: featureLabel(candidate.anchor, displayFacts, candidate.displayName),
    coordinates,
    displayedStrand: strandLabel(candidate.displayedStrand),
    featureIdentifier: candidate.anchor.biologicalFeatureId,
    representative: candidate.representative,
    role: candidate.role || 'member',
    directEvidence: candidate.directEvidence.length ? candidate.directEvidence : ['None'],
    strandRelation: candidate.strandRelation,
    usable: candidate.usable
  };
};

const referenceView = (response, request, displayFacts, recordLabels) => {
  const member = request.members.find(({ anchor }) => sameJson(anchor, response.reference));
  if (!member) throw new Error('Alignment reference has no current member facts.');
  return {
    label: featureLabel(response.reference, displayFacts),
    featureIdentifier: response.reference.biologicalFeatureId,
    recordLabel: recordLabels.get(response.reference.recordKey) || response.reference.recordKey,
    recordKey: response.reference.recordKey,
    coordinates: (member.sourceStart + 1).toLocaleString('en-US') + '..'
      + member.sourceEnd.toLocaleString('en-US') + ' bp',
    displayedStrand: strandLabel(response.referenceDisplayedStrand)
  };
};

const reviewRows = (response, request, displayFacts, recordLabels) => response.records
  .filter(({ recordKey }) => recordKey !== response.reference.recordKey)
  .map((record) => {
    const selectedAnchor = record.kind === 'ambiguous'
      ? record.recommendedAnchor
      : record.status === 'aligned' ? record.anchor : null;
    const reason = record.kind === 'ambiguous'
      ? record.recommendationReason
      : record.status === 'skipped' ? record.rationale : record.reviewReason;
    const recommendedKey = selectedAnchor && REVIEW_REASONS.has(reason)
      ? anchorKey(selectedAnchor) : null;
    return {
      recordKey: record.recordKey,
      recordLabel: recordLabels.get(record.recordKey) || `Record ${request.records.findIndex(
        ({ recordKey }) => recordKey === record.recordKey
      ) + 1}`,
      candidates: record.candidates.filter(({ usable }) => usable)
        .map((candidate) => candidateView(candidate, request, displayFacts)),
      reason,
      reasonLabel: reason
        ? reviewReasonLabels[reason] || rationaleLabels[reason] : 'Unchanged',
      recommendedKey,
      recommendationReasonLabel: recommendedKey ? reviewReasonLabels[reason] : null,
      choice: record.kind === 'ambiguous' ? null : selectedAnchor
        ? { kind: 'select', candidateKey: anchorKey(selectedAnchor) }
        : { kind: 'skip', candidateKey: null },
      unchanged: !selectedAnchor,
      repairRequired: false
    };
  });

const successfulSummary = (plan, reversed) => {
  const aligned = plan.records.filter(({ status }) => status === 'aligned').length;
  const explicitlySkipped = plan.records.filter(
    ({ rationale }) => rationale === 'skipped_by_user'
  ).length;
  const noCandidate = plan.records.filter(({ rationale }) => (
    rationale === 'skipped_no_candidate' || rationale === 'skipped_unmappable'
  )).length;
  const unchanged = plan.records.filter(({ status }) => status === 'reference').length;
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

const baseReverseComplement = canonicalRecordReverseComplement;

const inspectActivePlan = (plan, linearSeqs) => {
  if (!plan || plan.schema !== 2 || !Array.isArray(plan.records)) return null;
  const records = new Map(
    (Array.isArray(linearSeqs) ? linearSeqs : [])
      .map((record) => [String(record?.uid || ''), Boolean(record.region_reverse)])
  );
  return deepFreeze({
    groupId: plan.groupId,
    modeLabel: 'Align',
    reference: {
      ...plan.reference,
      label: plan.reference.biologicalFeatureId
    },
    records: plan.records.map((decision) => {
      return {
        recordKey: decision.recordKey,
        anchorLabel: decision.anchor?.biologicalFeatureId || 'Skip',
        status: decision.status,
        rationale: decision.rationale,
        rationaleLabel: rationaleLabels[decision.rationale] || decision.rationale,
        reversedFromSource: records.get(decision.recordKey) || false
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

const validateCandidate = (value, path, recordKey) => {
  const raw = exactObject(value, [
    'anchor', 'displayName', 'sourceStart', 'sourceEnd', 'displayCenter',
    'displayedStrand', 'hidden', 'representative', 'role', 'usable',
    'directEvidence', 'strandRelation'
  ], path);
  const anchor = validateAnchor(raw.anchor, `${path}.anchor`);
  if (anchor.recordKey !== recordKey || typeof raw.displayName !== 'string'
    || !Number.isSafeInteger(raw.sourceStart) || !Number.isSafeInteger(raw.sourceEnd)
    || raw.sourceStart < 0 || raw.sourceEnd < raw.sourceStart
    || (raw.displayCenter !== null && !Number.isFinite(raw.displayCenter))
    || ![null, -1, 1].includes(raw.displayedStrand)
    || typeof raw.hidden !== 'boolean' || typeof raw.representative !== 'boolean'
    || typeof raw.usable !== 'boolean' || typeof raw.role !== 'string'
    || !Array.isArray(raw.directEvidence)
    || raw.directEvidence.some((evidence) => typeof evidence !== 'string')) {
    throw new Error(`${path} contains invalid candidate facts.`);
  }
  if (!['same', 'opposite', 'unknown'].includes(raw.strandRelation)) {
    throw new Error(`${path}.strandRelation is invalid.`);
  }
  return { ...raw, anchor };
};

const validateDecision = (value, path) => {
  const raw = exactObject(value, [
    'kind', 'recordKey', 'status', 'rationale', 'reviewReason', 'anchor',
    'candidates'
  ], path);
  if (raw.kind !== 'decision' || !STATUSES.has(raw.status)
    || !RATIONALES.has(raw.rationale)
    || (raw.reviewReason !== null && !REVIEW_REASONS.has(raw.reviewReason))) {
    throw new Error(`${path} contains an unknown decision value.`);
  }
  const recordKey = text(raw.recordKey, `${path}.recordKey`);
  const anchor = raw.anchor === null ? null : validateAnchor(raw.anchor, `${path}.anchor`);
  const candidates = raw.candidates.map((candidate, index) => (
    validateCandidate(candidate, `${path}.candidates[${index}]`, recordKey)
  ));
  if (anchor && anchor.recordKey !== recordKey) throw new Error(`${path}.anchor belongs to another record.`);
  const alignedRationales = new Set(['user_selected', 'only_usable_candidate', 'unique_direct_rbh']);
  const skippedRationales = new Set(['skipped_by_user', 'skipped_no_candidate', 'skipped_unmappable']);
  const valid = raw.status === 'reference'
    ? anchor !== null && raw.rationale === 'reference' && raw.reviewReason === null
    : raw.status === 'aligned'
      ? anchor !== null && alignedRationales.has(raw.rationale)
        && candidates.some((candidate) => candidate.usable && sameJson(candidate.anchor, anchor))
        && raw.reviewReason === (raw.rationale === 'user_selected' ? null : raw.rationale)
      : anchor === null && skippedRationales.has(raw.rationale)
        && raw.reviewReason === null;
  if (!valid) throw new Error(`${path} contains an invalid decision combination.`);
  return { ...raw, anchor, candidates };
};

const validateAmbiguity = (value, path) => {
  const raw = exactObject(value, [
    'kind', 'recordKey', 'candidates', 'directRbhCandidates',
    'recommendedAnchor', 'recommendationReason'
  ], path);
  if (raw.kind !== 'ambiguous' || !Array.isArray(raw.candidates)
    || raw.candidates.length < 2 || !['unique_representative', 'deterministic_candidate_1']
      .includes(raw.recommendationReason)) throw new Error(`${path} is not a valid ambiguity.`);
  const recordKey = text(raw.recordKey, `${path}.recordKey`);
  const candidates = raw.candidates.map((candidate, index) => (
    validateCandidate(candidate, `${path}.candidates[${index}]`, recordKey)
  ));
  const keys = new Set(candidates.map(({ anchor }) => anchorKey(anchor)));
  const recommendedAnchor = validateAnchor(raw.recommendedAnchor, `${path}.recommendedAnchor`);
  if (keys.size !== candidates.length || !candidates.some((candidate) => (
    candidate.usable && sameJson(candidate.anchor, recommendedAnchor)
  )) || !Array.isArray(raw.directRbhCandidates)) {
    throw new Error(`${path} contains invalid candidates or recommendation.`);
  }
  const directRbhCandidates = raw.directRbhCandidates.map((anchor, index) => (
    validateAnchor(anchor, `${path}.directRbhCandidates[${index}]`)
  ));
  if (directRbhCandidates.some((anchor) => !keys.has(anchorKey(anchor)))) {
    throw new Error(`${path} contains unknown direct evidence.`);
  }
  return { ...raw, recordKey, candidates, directRbhCandidates, recommendedAnchor };
};

const validatePlan = (value, response, path) => {
  const raw = exactObject(value, ['schema', 'groupId', 'reference', 'records'], path);
  if (raw.schema !== 2 || raw.groupId !== response.groupId || !Array.isArray(raw.records)
    || !sameJson(validateAnchor(raw.reference, `${path}.reference`), response.reference)) {
    throw new Error(`${path} does not match its resolution.`);
  }
  const records = raw.records.map((record, index) => {
    const decision = response.records[index];
    const fields = exactObject(record, [
      'recordKey', 'status', 'rationale', 'anchor'
    ], `${path}.records[${index}]`);
    if (!decision || decision.kind !== 'decision'
      || !sameJson(fields, {
        recordKey: decision.recordKey, status: decision.status,
        rationale: decision.rationale, anchor: decision.anchor
      })) throw new Error(`${path} differs from helper decisions.`);
    return fields;
  });
  return { schema: 2, groupId: raw.groupId, reference: raw.reference, records };
};

const validateProjection = (value, request) => {
  if (value === null) return null;
  const raw = exactObject(value, ['binding', 'geometryOrientations', 'records'], 'alignment projection');
  const keys = request.records.map(({ recordKey }) => recordKey);
  if (!/^[a-f0-9]{64}$/.test(raw.binding) || !Array.isArray(raw.records)
    || raw.records.length !== keys.length) throw new Error('Invalid alignment projection binding or coverage.');
  exactObject(raw.geometryOrientations, keys, 'alignment projection.geometryOrientations');
  const known = new Set(request.members.map(({ anchor }) => anchorKey(anchor)));
  const finite = (number) => { if (!Number.isFinite(number)) throw new Error('Nonfinite alignment placement.'); };
  raw.records.forEach((record, index) => {
    exactObject(record, ['recordKey', 'beforeReverseComplement', 'base', 'beforeAxisY',
      'beforeAnchors', 'variants'], 'alignment projection record');
    if (record.recordKey !== keys[index] || typeof raw.geometryOrientations[record.recordKey] !== 'boolean'
      || record.beforeReverseComplement !== baseReverseComplement(request.records[index])) {
      throw new Error('Alignment projection changed record orientation binding.');
    }
    exactObject(record.base, ['x', 'y'], 'alignment projection base');
    [record.base.x, record.base.y, record.beforeAxisY].forEach(finite);
    const validateAnchors = (anchors, variant) => {
      if (!Array.isArray(anchors)) throw new Error('Alignment projection anchors must be an array.');
      const identities = new Set();
      anchors.forEach((entry) => {
        exactObject(entry, variant ? ['anchor', 'displayedStrand', 'displayCenter', 'centerX']
          : ['anchor', 'centerX'], 'alignment projection anchor');
        const anchor = validateAnchor(entry.anchor, 'alignment projection anchor identity');
        const key = anchorKey(anchor);
        if (anchor.recordKey !== record.recordKey || !known.has(key) || identities.has(key)) {
          throw new Error('Alignment projection changed anchor coverage.');
        }
        identities.add(key);
        if (variant) {
          if (![null, -1, 1].includes(entry.displayedStrand)
            || (entry.displayCenter === null) !== (entry.centerX === null)) {
            throw new Error('Alignment projection contains invalid strand or center facts.');
          }
          if (entry.displayCenter !== null) finite(entry.displayCenter);
        }
        if (!variant || entry.centerX !== null) finite(entry.centerX);
      });
      return identities;
    };
    const beforeAnchors = validateAnchors(record.beforeAnchors, false);
    if (!Array.isArray(record.variants) || record.variants.length !== 2) {
      throw new Error('Alignment projection requires both absolute directions.');
    }
    record.variants.forEach((variant, direction) => {
      exactObject(variant, ['reverseComplement', 'axisY', 'anchors'], 'alignment projection variant');
      if (variant.reverseComplement !== Boolean(direction)) throw new Error('Invalid alignment direction variant.');
      finite(variant.axisY);
      const anchors = validateAnchors(variant.anchors, true);
      const expected = request.members.filter(({ anchor }) => anchor.recordKey === record.recordKey);
      if (anchors.size !== expected.length || expected.some(({ anchor }) => !anchors.has(anchorKey(anchor)))) {
        throw new Error('Alignment projection omitted a source anchor.');
      }
      if (variant.reverseComplement === record.beforeReverseComplement
        && variant.anchors.some(({ anchor, centerX }) => (centerX !== null) !== beforeAnchors.has(anchorKey(anchor)))) {
        throw new Error('Alignment projection changed current anchor mappability.');
      }
    });
  });
  return raw;
};

// One transient resolver is shared by local previews and final helper admission.
export const projectSimilarityAlignmentDirections = ({ resolution, intent, expectedBinding, choices = [] }) => {
  const facts = resolution?.projection;
  if (!facts) throw new Error('Source-bound alignment projection facts are required.');
  if (typeof expectedBinding !== 'string') throw new Error('An explicit alignment binding is required.');
  if (facts.binding !== expectedBinding) return deepFreeze({ status: 'stale', reason: 'binding_changed', binding: facts.binding });
  const mode = intent?.mode;
  if (!['keep', 'right', 'left', 'custom'].includes(mode)) throw new Error('Invalid alignment direction mode.');
  exactObject(intent, mode === 'custom' ? ['mode', 'byRecordKey'] : ['mode'], 'alignment direction intent');
  if (mode === 'custom') {
    if (!isObject(intent.byRecordKey) || Object.entries(intent.byRecordKey).some(([key, value]) => (
      !facts.records.some(({ recordKey }) => recordKey === key) || !['keep', 'right', 'left'].includes(value)
    ))) throw new Error('Invalid Custom alignment choices.');
  }
  if (!Array.isArray(choices)) throw new Error('Alignment choices must be an array.');
  const selectedChoices = new Map();
  choices.forEach((choice) => {
    exactObject(choice, ['recordKey', 'kind', 'anchor'], 'local alignment choice');
    const record = resolution.records.find(({ recordKey }) => recordKey === choice.recordKey);
    if (!record || choice.recordKey === resolution.reference.recordKey || selectedChoices.has(choice.recordKey)
      || !['select', 'skip'].includes(choice.kind)) throw new Error('Invalid local alignment choice.');
    if (choice.kind === 'skip' ? choice.anchor !== null : !record.candidates.some((candidate) => (
      candidate.usable && sameJson(candidate.anchor, validateAnchor(choice.anchor, 'local alignment anchor'))
    ))) throw new Error('Local alignment selection is not a usable Python candidate.');
    selectedChoices.set(choice.recordKey, choice);
  });
  const records = facts.records.map((fact, index) => {
    const row = resolution.records[index];
    const choice = selectedChoices.get(fact.recordKey);
    const decision = choice
      ? { status: choice.kind === 'select' ? 'aligned' : 'skipped', anchor: choice.anchor, rationale: 'skipped_by_user' }
      : row.kind === 'ambiguous' ? { status: 'ambiguous', anchor: null, rationale: 'selection_required' } : row;
    const current = fact.variants[Number(fact.beforeReverseComplement)];
    const selected = decision.anchor === null ? null : current.anchors.find(({ anchor }) => sameJson(anchor, decision.anchor));
    let exclusion = ['skipped', 'ambiguous'].includes(decision.status) ? decision.rationale : null;
    if (!exclusion && (!selected || selected.centerX === null)) exclusion = 'unusable_anchor';
    if (!exclusion && selected.displayedStrand === null) exclusion = 'unknown_strand';
    const requested = mode === 'custom' ? intent.byRecordKey[fact.recordKey] || 'keep' : mode;
    const flip = !exclusion && requested !== 'keep'
      && selected.displayedStrand !== (requested === 'right' ? 1 : -1);
    const afterReverseComplement = flip ? !fact.beforeReverseComplement : fact.beforeReverseComplement;
    const variant = fact.variants[Number(afterReverseComplement)];
    const afterAnchor = selected ? variant.anchors.find(({ anchor }) => sameJson(anchor, selected.anchor)) : null;
    if (!exclusion && afterAnchor?.centerX === null) throw new Error('Requested direction makes the selected anchor unusable.');
    return { recordKey: fact.recordKey, anchor: decision.anchor, status: decision.status,
      beforeReverseComplement: fact.beforeReverseComplement, afterReverseComplement,
      beforeArrow: selected?.displayedStrand ?? null, afterArrow: afterAnchor?.displayedStrand ?? null,
      exclusion, beforeCenterX: fact.beforeAnchors.find(({ anchor }) => sameJson(anchor, decision.anchor))?.centerX ?? null,
      afterCenterX: afterAnchor?.centerX ?? null,
      translation: { x: fact.base.x, y: fact.base.y + fact.beforeAxisY - variant.axisY } };
  });
  const reference = records.find(({ recordKey }) => recordKey === resolution.reference.recordKey);
  const fact = facts.records.find(({ recordKey }) => recordKey === reference.recordKey);
  if (!Number.isFinite(reference.beforeCenterX) || !Number.isFinite(reference.afterCenterX)) {
    throw new Error('The exact reference center is unavailable.');
  }
  const beforeX = reference.beforeCenterX + fact.base.x;
  reference.translation.x = beforeX - reference.afterCenterX;
  const deltaX = reference.translation.x - fact.base.x;
  if (!Number.isFinite(deltaX) || records.some(({ translation }) => !Number.isFinite(translation.y))) {
    throw new Error('Alignment correction is nonfinite.');
  }
  const signature = JSON.stringify({
    records: records.map(({ recordKey, anchor, status, afterReverseComplement }) => ({ recordKey, anchor, status, afterReverseComplement })),
    referenceBeforeX: beforeX, referenceDeltaX: deltaX
  });
  return deepFreeze({ status: 'projected', binding: facts.binding, records,
    selectionComplete: records.every(({ status }) => status !== 'ambiguous'),
    reference: { anchor: resolution.reference, beforeX, afterX: reference.afterCenterX + reference.translation.x, deltaX },
    geometryValidated: records.every(({ recordKey, afterReverseComplement }) => (
      facts.geometryOrientations[recordKey] === afterReverseComplement
    )), signature });
};

export const validateSimilarityAlignmentResolution = (value, request) => {
  const raw = exactObject(value, [
    'schema', 'status', 'groupId', 'reference', 'referenceDisplayedStrand',
    'referenceDisplayCenter', 'records', 'plan', 'projection'
  ], 'alignment helper response');
  if (raw.schema !== 2 || !['resolved', 'ambiguous'].includes(raw.status)
    || raw.groupId !== request.groupId || !Array.isArray(raw.records)
    || ![null, -1, 1].includes(raw.referenceDisplayedStrand)
    || !Number.isFinite(raw.referenceDisplayCenter)) {
    throw new Error('Alignment helper response does not match its request.');
  }
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
  const referenceDecision = records.find(({ recordKey }) => recordKey === reference.recordKey);
  if (referenceDecision?.kind !== 'decision' || referenceDecision.status !== 'reference'
    || !sameJson(referenceDecision.anchor, reference)
    || records.filter(({ status }) => status === 'reference').length !== 1) {
    throw new Error('Alignment helper response changed the exact reference decision.');
  }
  const expectedKeys = request.records.map(({ recordKey }) => recordKey);
  if (records.length !== expectedKeys.length
    || records.some((record, index) => record.recordKey !== expectedKeys[index])) {
    throw new Error('Alignment helper response changed displayed-record coverage.');
  }
  const known = new Set(request.members.map(({ anchor }) => anchorKey(anchor)));
  if (records.some((record) => (
    record.kind === 'decision'
      ? record.anchor !== null && !known.has(anchorKey(record.anchor))
      : record.candidates.some(({ anchor }) => !known.has(anchorKey(anchor)))
  ))) throw new Error('Alignment helper response contains an unknown group member.');
  const ambiguous = records.some(({ kind }) => kind === 'ambiguous');
  if ((raw.status === 'ambiguous') !== ambiguous) {
    throw new Error('Alignment helper response status does not match its records.');
  }
  const response = {
    schema: 2, status: raw.status, groupId: raw.groupId, reference,
    referenceDisplayedStrand: raw.referenceDisplayedStrand,
    referenceDisplayCenter: raw.referenceDisplayCenter,
    records, plan: null, projection: validateProjection(raw.projection, request)
  };
  response.projection?.records.forEach((fact, index) => {
    const current = fact.variants[Number(fact.beforeReverseComplement)];
    const candidates = records[index].candidates;
    candidates.forEach((candidate) => {
      const projected = current.anchors.find(({ anchor }) => sameJson(anchor, candidate.anchor));
      if (!projected || projected.displayedStrand !== candidate.displayedStrand
        || projected.displayCenter !== candidate.displayCenter) {
        throw new Error('Alignment projection disagrees with Python candidate facts.');
      }
    });
    if (fact.recordKey === reference.recordKey) {
      const selected = current.anchors.find(({ anchor }) => sameJson(anchor, reference));
      if (!selected || selected.displayedStrand !== raw.referenceDisplayedStrand
        || selected.displayCenter !== raw.referenceDisplayCenter) {
        throw new Error('Alignment projection changed the exact reference facts.');
      }
    }
  });
  if (raw.status === 'resolved') {
    if (raw.plan === null) throw new Error('Resolved alignment helper response has no plan.');
    response.plan = validatePlan(raw.plan, response, 'alignment helper response.plan');
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

const recordLengths = (catalog, recordCatalog, linearSeqs) => {
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
  if (recordCatalog?.status === 'ready') {
    const discovered = new Map();
    (recordCatalog.records || []).forEach((record) => {
      const index = record?.sourceIndex;
      if (!Number.isSafeInteger(index) || index < 0) return;
      discovered.set(index, discovered.has(index) ? null : record);
    });
    discovered.forEach((record, index) => {
      const key = String(linearSeqs?.[index]?.uid || '');
      if (key && record && Number.isSafeInteger(record.recordLength)
          && record.recordLength > 0 && !lengths.has(key)) {
        lengths.set(key, record.recordLength);
      }
    });
  }
  return lengths;
};

const buildHelperRequest = ({ group, members, reference, request, catalog,
  recordCatalog, linearSeqs, choices }) => {
  const groupStatus = orthogroupIdStatus(group);
  if (!groupStatus.valid || !groupStatus.supplied) {
    throw new Error('The selected Similarity Group identity is invalid.');
  }

  if (!request || request.mode !== 'linear' || !Array.isArray(request.records)) {
    throw new Error('Generate a Linear diagram before aligning a Similarity Group.');
  }
  const lengths = recordLengths(catalog, recordCatalog, linearSeqs);
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
  if (!anchorsAgree(exactMember.payload.anchor, exactReference)) {
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
    schema: 2,
    groupId: groupStatus.value,
    records,
    reference: exactMember.payload.anchor,
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

const anchorsAgree = (left, right) => {
  if (!left || !right) return false;
  if (
    left.recordKey !== right.recordKey
    || left.biologicalFeatureId !== right.biologicalFeatureId
  ) return false;
  return ['sourceFeatureIndex', 'stableFeatureSvgId'].every((field) => (
    left[field] === null
    || left[field] === undefined
    || right[field] === null
    || right[field] === undefined
    || left[field] === right[field]
  ));
};

export const createSimilarityAlignmentActions = ({
  state,
  getOrthogroupById,
  getEnrichedOrthogroupMembers,
  getCommittedRequest,
  getRecordCatalog = null,
  getCommittedSession,
  projectCommittedAlignment,
  runCommittedCanonicalCandidate,
  recordDisplayControls,
  cancelRunAnalysis = null,
  runHelperOperation,
  resolveOperation,
  getCurrentSvg = null,
  previewCandidate = null,
  clearCandidatePreview = null,
  onError = null
}) => {
  if (
    !state
    || typeof getOrthogroupById !== 'function'
    || typeof getEnrichedOrthogroupMembers !== 'function'
    || typeof getCommittedRequest !== 'function'
    || typeof getCommittedSession !== 'function'
    || typeof projectCommittedAlignment !== 'function'
    || typeof runCommittedCanonicalCandidate !== 'function'
    || typeof runHelperOperation !== 'function'
    || typeof resolveOperation !== 'string'
    || !resolveOperation
  ) throw new Error('Similarity alignment actions require their current artifact owners.');

  const draft = ref(null);
  const status = ref('idle');
  const busy = computed(() => status.value === 'resolving' || status.value === 'applying');
  const error = ref(null);
  const summary = ref(null);
  const notice = ref('');
  const repair = ref(null);
  const drawerReferenceKey = ref('');
  const resetDialogOpen = ref(false);
  const automaticApply = ref(false);
  const resetScope = ref('positions');
  let actionId = 0;
  let activeRequest = null;
  let activeApply = null;
  let activeBaseline = null;
  let displayFacts = new Map();
  let recordLabels = new Map();
  let artifactStamp = null;
  let materializedRecordDrag = false;
  let pendingRecordDragBaseline = null;

  const publishError = (value, notify = true) => {
    const normalized = value instanceof Error ? value : new Error(String(value || 'Alignment failed.'));
    error.value = normalized;
    if (notify && typeof onError === 'function') onError(normalized);
    return normalized;
  };

  const clearDraft = () => {
    draft.value = null;
    activeRequest = null;
    displayFacts = new Map();
    recordLabels = new Map();
    artifactStamp = null;
    if (typeof clearCandidatePreview === 'function') clearCandidatePreview();
  };

  const currentRequest = () => getCommittedRequest();
  const captureArtifact = (groupId) => ({
    records: JSON.stringify(currentRequest()?.records || []),
    group: JSON.stringify(getOrthogroupById(groupId)),
    selectedGroup: state.selectedOrthogroupId?.value,
    catalog: state.featureCatalog?.value,
    result: state.results?.value?.[state.selectedResultIndex?.value ?? 0],
    translations: JSON.stringify(state.linearRecordTranslations?.value || [])
  });
  const artifactIsCurrent = () => {
    if (!artifactStamp || !activeRequest) return false;
    const current = captureArtifact(activeRequest.groupId);
    return current.records === artifactStamp.records
      && current.group === artifactStamp.group
      && current.selectedGroup === artifactStamp.selectedGroup
      && current.catalog === artifactStamp.catalog
      && current.result === artifactStamp.result
      && current.translations === artifactStamp.translations;
  };
  const rejectStaleDraft = () => {
    clearDraft();
    activeBaseline = null;
    status.value = 'idle';
    publishError(new Error('The diagram, source, crop, or Similarity Group changed. Start alignment again.'));
    return { status: 'stale' };
  };

  const currentRecordKeys = () => (
    (Array.isArray(currentRequest()?.records) ? currentRequest().records : [])
      .map(({ recordKey }) => String(recordKey || ''))
      .filter(Boolean)
  );

  const currentSvg = () => (
    typeof getCurrentSvg === 'function' ? getCurrentSvg() : null
  );

  const materializedTranslations = (recordKeys, plan = null) => {
    const base = replacementTranslations(state, recordKeys);
    const svg = currentSvg();
    return svg ? materializeRecordTranslations(svg, base, recordKeys, plan) : base;
  };

  const baseline = ({ materializePlan = false } = {}) => {
    const canonical = getCommittedSession();
    const request = cloneJson(canonical?.renderRequest);
    if (!request) throw new Error('Generate a Linear diagram before aligning.');
    const recordKeys = request.records.map(({ recordKey }) => recordKey);
    const translations = materializePlan
      ? materializedTranslations(recordKeys, request.layout?.similarityAlignment)
      : cloneJson(request.layout?.recordTranslations || replacementTranslations(state, recordKeys));
    const orientations = request.records.map(record => ({ recordKey: record.recordKey,
      reverseComplement: baseReverseComplement(record) }));
    const materialized = projectCommittedAlignment({ committed: canonical,
      plan: materializePlan ? null : request.layout?.similarityAlignment,
      translations, orientations });
    return { canonical: materialized, request: materialized.renderRequest, translations };
  };

  const installBaseState = (value) => {
    if (!value) return;
    state.linearRecordTranslations.value = cloneJson(value.translations);
  };

  const publishNotice = (message) => {
    notice.value = String(message || '');
  };

  const clearCommittedPlan = (reason) => {
    const hadPlan = Boolean(state.similarityAlignmentPlan?.value);
    if (activeApply && typeof cancelRunAnalysis === 'function') cancelRunAnalysis();
    actionId += 1;
    clearDraft();
    status.value = 'idle';
    activeBaseline = null;
    state.similarityAlignmentResetReceipt.value = null;
    resetDialogOpen.value = false;
    if (!hadPlan) return false;
    state.similarityAlignmentPlan.value = null;
    summary.value = null;
    repair.value = null;
    publishNotice(`Alignment cleared: ${reason}`);
    return true;
  };

  const rowChoices = () => (draft.value?.rows || []).filter(row => row.choice).map(row => ({
    recordKey: row.recordKey, kind: row.choice.kind,
    anchor: row.choice.kind === 'select'
      ? row.candidates.find(({ key }) => key === row.choice.candidateKey)?.anchor : null
  }));
  const projectDraft = (response = draft.value?.response, choices = rowChoices()) => (
    projectSimilarityAlignmentDirections({ resolution: response,
      expectedBinding: draft.value?.binding || response.projection.binding,
      intent: draft.value?.intent || { mode: 'keep' }, choices })
  );
  const installReviewDraft = (response, request) => {
    if (!response.projection) throw new Error('Source-bound alignment facts are required.');
    draft.value = deepFreeze({ response, binding: response.projection.binding,
      intent: { mode: 'keep' },
      reference: referenceView(response, request, displayFacts, recordLabels),
      rows: reviewRows(response, request, displayFacts, recordLabels) });
    status.value = 'reviewing';
  };

  const runAlignmentCandidate = async ({ canonical, label, alignmentResetBefore = null,
    alignmentResetReceipt = undefined }) => {
    const beforeRecords = new Map(currentRequest().records.map(record => [record.recordKey, record]));
    const changed = canonical.renderRequest.records.filter(record => (
      baseReverseComplement(record) !== baseReverseComplement(beforeRecords.get(record.recordKey))
    )).map(record => ({recordKey: record.recordKey, reverseComplement: baseReverseComplement(record)}));
    return runCommittedCanonicalCandidate({ canonical, label, alignmentResetBefore, alignmentResetReceipt,
      captureIntentCheckpoint: () => recordDisplayControls.captureAlignmentOrientationIntent(changed),
      restoreIntentCheckpoint: checkpoint => recordDisplayControls.restoreAlignmentOrientationIntent(checkpoint),
      commitIntent: () => recordDisplayControls.commitAlignmentOrientations(changed) });
  };

  const applyPlan = async (response, request, expectedActionId, projected) => {
    if (expectedActionId !== actionId) return { status: 'stale' };
    if (!artifactIsCurrent()) return rejectStaleDraft();
    if (!projected.selectionComplete || !projected.geometryValidated) {
      throw new Error('Alignment selection or final geometry is not validated.');
    }
    const orientations = projected.records.map(({ recordKey, afterReverseComplement }) => ({
      recordKey, reverseComplement: afterReverseComplement }));
    const canonical = projectCommittedAlignment({ committed: activeBaseline.canonical,
      plan: response.plan, orientations,
      translations: projected.records.map(({ recordKey, translation }) => ({recordKey, ...translation})) });
    const promise = runAlignmentCandidate({canonical, label: 'Align Similarity Group',
      alignmentResetBefore: activeBaseline.canonical});
    activeApply = promise;
    let outcome;
    try { outcome = await promise; }
    catch (cause) { outcome = {status:'error', error:cause}; }
    finally { if (activeApply === promise) activeApply = null; }
    if (expectedActionId !== actionId) return {status:'stale'};
    if (outcome?.status === 'ok') {
      summary.value = successfulSummary(response.plan, projected.records.filter(record => (
        record.beforeReverseComplement !== record.afterReverseComplement)).length);
      repair.value = null; notice.value = ''; activeBaseline = null;
      clearDraft(); error.value = null; status.value = 'idle';
      return {status:'ok'};
    }
    status.value = 'reviewing';
    const cause = outcome?.error || state.errorLog?.value;
    publishError(new Error(cause?.summary || cause?.message || 'Alignment generation failed. Review the draft and retry Apply.'), !state.errorLog?.value?.summary);
    return {status: outcome?.status || 'error'};
  };

  const resolveRequest = async (request, expectedActionId, mode) => {
    let response;
    try {
      const helper = await runHelperOperation(resolveOperation, { request,
        projection: { canonicalRequest: activeBaseline.request, orientations: Object.fromEntries(
          request.records.map(record => [record.recordKey, baseReverseComplement(record)])) },
        resources: activeBaseline.canonical.resources });
      if (expectedActionId !== actionId) return { status: 'stale' };
      response = validateSimilarityAlignmentResolution(helper?.result, request);
      if (!artifactIsCurrent()) return rejectStaleDraft();
    } catch (cause) {
      if (expectedActionId !== actionId) return { status: 'stale' };
      clearDraft();
      activeBaseline = null;
      status.value = 'idle';
      publishError(cause);
      return { status: 'error' };
    }
    activeRequest = request;
    error.value = null;
    installReviewDraft(response, request);
    if (mode === 'align' && response.status === 'resolved') {
      status.value = 'applying';
      automaticApply.value = true;
      try { return await applyPlan(response, request, expectedActionId, projectDraft(response, [])); }
      finally { automaticApply.value = false; }
    }
    return { status: 'reviewing' };
  };

  const start = async ({ groupId, reference, source, mode = 'align' }) => {
    if (busy.value) return { status: 'busy' };
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
    clearDraft();
    error.value = null;
    status.value = 'resolving';
    try {
      const group = getOrthogroupById(id);
      if (!group) throw new Error('The selected Similarity Group is unavailable.');
      activeBaseline = baseline({ materializePlan: true });
      const members = getEnrichedOrthogroupMembers(group);
      const request = buildHelperRequest({
        group, members, reference,
        request: activeBaseline.request,
        catalog: state.featureCatalog?.value,
        recordCatalog: getRecordCatalog?.(),
        linearSeqs: state.linearSeqs,
        choices: []
      });
      displayFacts = new Map(members.map((member) => [
        anchorKey(anchorFromSource(member, 'Alignment member')), member
      ]));
      recordLabels = new Map(request.records.map(({ recordKey }, index) => {
        const sequence = (state.linearSeqs || []).find((entry) => String(entry?.uid) === recordKey);
        return [recordKey, displayText(sequence?.definition, sequence?.accession,
          sequence?.gb?.name, sequence?.gff?.name) || `Record ${index + 1}`];
      }));
      activeRequest = request;
      artifactStamp = captureArtifact(request.groupId);
      return resolveRequest(request, expectedActionId, mode);
    } catch (cause) {
      if (expectedActionId !== actionId) return { status: 'stale' };
      clearDraft();
      status.value = 'idle';
      activeBaseline = null;
      publishError(cause);
      return { status: 'rejected' };
    }
  };

  const editRow = (recordKey, update) => {
    if (!activeRequest || !draft.value || status.value !== 'reviewing') {
      return { status: 'rejected' };
    }
    const row = draft.value.rows.find((entry) => entry.recordKey === recordKey);
    if (!row) return { status: 'rejected' };
    let next = row;
    if (update.kind === 'select') {
      const selectedAnchor = validateAnchor(update.anchor, 'Selected alignment candidate');
      const candidate = row.candidates.find(({ anchor }) => sameJson(anchor, selectedAnchor));
      if (!candidate) return { status: 'rejected' };
      next = { ...row, choice: { kind: 'select', candidateKey: candidate.key },
        reason: 'user_selected', reasonLabel: 'Selected by user',
        unchanged: false, repairRequired: false };
    } else if (update.kind === 'skip') {
      if (!row.candidates.length) return { status: 'rejected' };
      next = { ...row, choice: { kind: 'skip', candidateKey: null },
        reason: 'skipped_by_user', reasonLabel: 'Skipped by user',
        unchanged: true, repairRequired: false };
    } else return { status: 'rejected' };
    draft.value = deepFreeze({
      ...draft.value,
      rows: draft.value.rows.map((entry) => entry === row ? next : entry)
    });
    error.value = null;
    return { status: 'selected' };
  };

  const applyDraft = async () => {
    if (!activeRequest || status.value !== 'reviewing' || unresolvedCount.value > 0) {
      return { status: 'rejected' };
    }
    if (!artifactIsCurrent()) return rejectStaleDraft();
    const expectedActionId = actionId;
    const choices = rowChoices();
    const request = deepFreeze({ ...cloneJson(activeRequest), choices });
    const preview = projectDraft();
    const orientations = Object.fromEntries(preview.records.map(record => [record.recordKey, record.afterReverseComplement]));
    status.value = 'applying';
    try {
      const helper = await runHelperOperation(resolveOperation, { request,
        projection: { canonicalRequest: activeBaseline.request, orientations },
        resources: activeBaseline.canonical.resources });
      if (expectedActionId !== actionId) return { status: 'stale' };
      const response = validateSimilarityAlignmentResolution(helper?.result, request);
      if (!artifactIsCurrent()) return rejectStaleDraft();
      const final = projectDraft(response, choices);
      if (final.status === 'stale') return rejectStaleDraft();
      if (!final.selectionComplete || response.status !== 'resolved' || !response.plan) {
        throw new Error('The resolver did not resolve every record. Review the choices and retry.');
      }
      if (final.signature !== preview.signature || !final.geometryValidated) {
        draft.value = deepFreeze({ ...draft.value, response });
        status.value = 'reviewing';
        publishError(new Error('Validated directions or reference placement changed. Review the updated preview and Apply again.'));
        return {status:'reviewing'};
      }
      return applyPlan(response, request, expectedActionId, final);
    } catch (cause) {
      if (expectedActionId !== actionId) return { status: 'stale' };
      status.value = 'reviewing'; publishError(cause); return { status: 'error' };
    }
  };

  const cancel = () => {
    actionId += 1;
    if (activeApply && typeof cancelRunAnalysis === 'function') cancelRunAnalysis();
    clearDraft();
    error.value = null;
    status.value = 'idle';
    activeBaseline = null;
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
      return inspectActivePlan(state.similarityAlignmentPlan?.value, state.linearSeqs);
    } catch (_error) {
      return null;
    }
  });

  const resetPreview = computed(() => {
    const receipt = state.similarityAlignmentResetReceipt.value;
    const request = currentRequest();
    const targets = (receipt?.directions || []).map(delta => {
      const record = request?.records.find(({recordKey}) => recordKey === delta.recordKey);
      const sequence = state.linearSeqs.find(({uid}) => uid === delta.recordKey);
      return {recordKey: delta.recordKey, label: displayText(record?.presentation?.label, sequence?.definition, sequence?.file_definition, sequence?.accession) || delta.recordKey,
        current: baseReverseComplement(record), restored: delta.before,
        laterManualEdit: baseReverseComplement(record) !== delta.after};
    });
    return { targets, disabledReason: !receipt
      ? 'This Session has no historical alignment direction evidence. Reset positions is available.'
      : !targets.length ? 'The latest Align made no direction changes. Reset positions is available.' : '' };
  });
  const resetAlignment = async (scope = 'positions') => {
    if (!state.similarityAlignmentPlan?.value) return {status:'noop'};
    if (!['positions', 'positions-and-directions'].includes(scope)) return {status:'rejected'};
    if (scope === 'positions-and-directions' && resetPreview.value.disabledReason) return {status:'rejected'};
    if (busy.value) return {status:'busy'};
    const expectedActionId = ++actionId;
    status.value = 'applying';
    let promise;
    try {
      const current = baseline();
      const receipt = await validateSimilarityAlignmentResetReceipt(state.similarityAlignmentResetReceipt.value,
        current.canonical);
      const translations = cloneJson(current.translations);
      if (receipt?.referenceDeltaX) {
        const delta = receipt.referenceDeltaX;
        const reference = translations.find(({recordKey}) => recordKey === delta.recordKey);
        if (!reference) throw new Error('Alignment reset reference placement is unavailable.');
        reference.x -= delta.deltaX;
      }
      const restore = new Map(scope === 'positions-and-directions'
        ? receipt.directions.map(delta => [delta.recordKey, delta.before]) : []);
      const orientations = current.request.records.map(record => ({recordKey: record.recordKey,
        reverseComplement: restore.has(record.recordKey) ? restore.get(record.recordKey) : baseReverseComplement(record)}));
      if (orientations.some((direction, index) => direction.reverseComplement
        !== baseReverseComplement(current.request.records[index]))) {
        const plan = current.request.layout.similarityAlignment;
        const catalogFeatures = (state.featureCatalog?.value?.items || [])
          .flatMap(item => item.biologicalFeatures || []);
        const members = plan.records.filter(record => record.anchor).map(({anchor}) => {
          const matches = catalogFeatures.filter(feature => feature.recordKey === anchor.recordKey
            && feature.biologicalFeatureId === anchor.biologicalFeatureId
            && (anchor.sourceFeatureIndex === null || feature.sourceFeatureIndex === anchor.sourceFeatureIndex));
          if (matches.length !== 1) throw new Error('Alignment Reset source anchor is unavailable or ambiguous.');
          return {...matches[0], ...anchor};
        });
        const request = buildHelperRequest({group:{id:plan.groupId,orthologEdges:[]}, members,
          reference:plan.reference, request:current.request, catalog:state.featureCatalog?.value,
          recordCatalog:getRecordCatalog?.(), linearSeqs:state.linearSeqs, choices:planChoices(plan)});
        const vector = Object.fromEntries(orientations.map(direction => [direction.recordKey,direction.reverseComplement]));
        const helper = await runHelperOperation(resolveOperation, {request,
          projection:{canonicalRequest:current.request,orientations:vector},resources:current.canonical.resources});
        if (expectedActionId !== actionId) return {status:'stale'};
        const facts = validateSimilarityAlignmentResolution(helper?.result,request).projection;
        if (!facts || orientations.some(direction => facts.geometryOrientations[direction.recordKey]
          !== direction.reverseComplement)) throw new Error('Alignment Reset geometry is not validated.');
        translations.forEach(translation => {
          const fact = facts.records.find(record => record.recordKey === translation.recordKey);
          translation.y += fact.beforeAxisY - fact.variants[Number(vector[translation.recordKey])].axisY;
        });
      }
      const canonical = projectCommittedAlignment({committed: current.canonical, plan:null, translations, orientations});
      promise = runAlignmentCandidate({canonical, label: scope === 'positions'
        ? 'Reset alignment positions' : 'Reset alignment positions and direction changes', alignmentResetReceipt:null});
      activeApply = promise;
      const outcome = await promise;
      if (expectedActionId !== actionId) return {status:'stale'};
      if (outcome?.status === 'ok') {
        clearDraft(); repair.value = null; summary.value = null; resetDialogOpen.value = false;
        error.value = null;
        publishNotice(scope === 'positions' ? 'Alignment reset: record positions restored; record directions unchanged.'
          : 'Alignment reset: positions and the latest Align direction changes restored.');
      } else publishError(new Error(outcome?.error?.message || state.errorLog?.value?.summary || 'Alignment Reset failed. Retry the same scope.'));
      return outcome;
    } catch (cause) { publishError(cause); return {status:'error', error:cause}; }
    finally { if (activeApply === promise) activeApply = null; if (expectedActionId === actionId) status.value = 'idle'; }
  };

  const planChoices = (plan, staleRecordKeys = new Set()) => (
    plan.records
      .filter(({ status, recordKey }) => status !== 'reference' && !staleRecordKeys.has(recordKey))
      .map((decision) => ({
        recordKey: decision.recordKey,
        kind: decision.status === 'aligned' ? 'select' : 'skip',
        anchor: decision.status === 'aligned' ? decision.anchor : null
      }))
  );

  const markStaleReference = (reason) => {
    clearDraft();
    status.value = 'idle';
    repair.value = deepFreeze({
      kind: 'reference',
      groupId: state.similarityAlignmentPlan?.value?.groupId || '',
      reason: String(reason || 'The exact reference is no longer available.')
    });
    publishNotice(`Alignment needs repair: ${repair.value.reason}`);
    return { status: 'blocked', reason: 'stale-reference' };
  };

  const prepareStaleTargetRepair = async (
    plan,
    request,
    staleRecordKeys,
    expectedActionId
  ) => {
    const repairRequest = deepFreeze({
      ...cloneJson(request),
      choices: planChoices(plan, staleRecordKeys)
    });
    let response;
    activeBaseline = baseline({materializePlan:true});
    try {
      const helper = await runHelperOperation(resolveOperation, { request: repairRequest,
        projection:{canonicalRequest:activeBaseline.request, orientations:null}, resources:activeBaseline.canonical.resources });
      if (expectedActionId !== actionId) return { status: 'stale' };
      response = validateSimilarityAlignmentResolution(helper?.result, repairRequest);
    } catch (cause) {
      if (expectedActionId !== actionId) return { status: 'stale' };
      return markStaleReference(cause?.message || cause);
    }
    activeRequest = repairRequest;
    artifactStamp = captureArtifact(repairRequest.groupId);
    const rows = reviewRows(response, repairRequest, displayFacts, recordLabels).map((row) => (
      staleRecordKeys.has(row.recordKey)
        ? { ...row, choice: null, repairRequired: true }
        : row
    ));
    draft.value = deepFreeze({
      response, binding: response.projection.binding, intent:{mode:'keep'},
      reference: referenceView(response, repairRequest, displayFacts, recordLabels),
      rows, repair: true
    });
    repair.value = deepFreeze({
      kind: 'targets',
      groupId: plan.groupId,
      reason: 'One or more saved target anchors are no longer usable. Select or Skip each target.'
    });
    status.value = 'reviewing';
    publishNotice(`Alignment needs repair: ${repair.value.reason}`);
    return { status: 'blocked', reason: 'stale-target' };
  };

  const validateBeforeGenerate = async () => {
    const expectedActionId = ++actionId;
    const plan = state.similarityAlignmentPlan?.value;
    if (!plan) {
      repair.value = null;
      return { status: 'ok' };
    }
    const group = getOrthogroupById(plan.groupId);
    if (!group) {
      // A saved Result can retain its exact plan after comparison groups have
      // left the editor catalog. The typed render still validates every anchor.
      if (sameKeyedPlan(plan, currentRequest()?.layout?.similarityAlignment)) {
        repair.value = null;
        return { status: 'ok' };
      }
      return markStaleReference('The saved Similarity Group is no longer available.');
    }
    const members = getEnrichedOrthogroupMembers(group);
    let request;
    try {
      request = buildHelperRequest({
        group,
        members,
        reference: plan.reference,
        request: currentRequest(),
        catalog: state.featureCatalog?.value,
        recordCatalog: getRecordCatalog?.(),
        linearSeqs: state.linearSeqs,
        choices: planChoices(plan)
      });
    } catch (cause) {
      return markStaleReference(cause?.message || cause);
    }
    const currentAnchors = members.map((member) => {
      try { return anchorFromSource(member, 'Current Similarity Group member'); } catch (_error) { return null; }
    }).filter(Boolean);
    const staleRecordKeys = new Set(
      plan.records
        .filter(({ status, anchor }) => (
          status === 'aligned'
          && !currentAnchors.some((candidate) => anchorsAgree(candidate, anchor))
        ))
        .map(({ recordKey }) => recordKey)
    );
    if (staleRecordKeys.size > 0) {
      return prepareStaleTargetRepair(plan, request, staleRecordKeys, expectedActionId);
    }
    try {
      const helper = await runHelperOperation(resolveOperation, { request });
      if (expectedActionId !== actionId) return { status: 'stale' };
      validateSimilarityAlignmentResolution(helper?.result, request);
    } catch (_cause) {
      if (expectedActionId !== actionId) return { status: 'stale' };
      const alignedKeys = new Set(
        plan.records.filter(({ status }) => status === 'aligned').map(({ recordKey }) => recordKey)
      );
      if (alignedKeys.size > 0) {
        return prepareStaleTargetRepair(plan, request, alignedKeys, expectedActionId);
      }
      return markStaleReference('The exact reference is no longer usable with current facts.');
    }
    repair.value = null;
    return { status: 'ok' };
  };

  const retainForStableReorder = (recordKeys) => {
    const orderedKeys = (Array.isArray(recordKeys) ? recordKeys : []).map(String);
    const plan = state.similarityAlignmentPlan?.value;
    if (plan) {
      const decisions = new Map(plan.records.map((entry) => [entry.recordKey, entry]));
      if (orderedKeys.every((recordKey) => decisions.has(recordKey))) {
        state.similarityAlignmentPlan.value = {
          ...cloneJson(plan),
          records: orderedKeys.map((recordKey) => cloneJson(decisions.get(recordKey)))
        };
      }
    }
    const translations = new Map(
      (state.linearRecordTranslations?.value || []).map((entry) => [entry.recordKey, entry])
    );
    if (orderedKeys.every((recordKey) => translations.has(recordKey))) {
      state.linearRecordTranslations.value = orderedKeys.map(
        (recordKey) => cloneJson(translations.get(recordKey))
      );
    }
  };

  const beforeRecordDrag = () => {
    pendingRecordDragBaseline = state.similarityAlignmentPlan?.value
      ? baseline({ materializePlan: true })
      : null;
    return Boolean(pendingRecordDragBaseline);
  };

  const afterRecordDrag = ({ moved = true } = {}) => {
    if (!pendingRecordDragBaseline || !moved) {
      pendingRecordDragBaseline = null;
      return false;
    }
    const materializedPlan = state.similarityAlignmentPlan?.value;
    installBaseState(pendingRecordDragBaseline);
    pendingRecordDragBaseline = null;
    materializedRecordDrag = clearCommittedPlan('record moved manually.');
    if (!materializedRecordDrag) return false;
    const recordKeys = currentRecordKeys();
    state.linearRecordTranslations.value = materializedTranslations(
      recordKeys,
      materializedPlan
    );
    materializedRecordDrag = false;
    return true;
  };

  const directionPreview = computed(() => draft.value ? projectDraft() : null);
  const setDirectionIntent = (intent) => {
    if (!draft.value || status.value !== 'reviewing') return {status:'rejected'};
    projectSimilarityAlignmentDirections({resolution:draft.value.response,
      expectedBinding:draft.value.binding, choices:rowChoices(), intent});
    draft.value = deepFreeze({...draft.value, intent:cloneJson(intent)});
    error.value = null;
    return {status:'selected'};
  };
  const setDirectionMode = (mode) => setDirectionIntent(mode === 'custom'
    ? {mode, byRecordKey:{}} : {mode});
  const setCustomDirection = (recordKey, direction) => setDirectionIntent({mode:'custom',
    byRecordKey:{...draft.value?.intent?.byRecordKey, [recordKey]:direction}});

  const unresolvedCount = computed(() => (
    draft.value?.rows?.filter(({ choice }) => choice === null).length || 0
  ));

  const applyDisabledReason = computed(() => {
    if (status.value === 'applying') return 'Applying the alignment plan.';
    if (status.value === 'resolving') return 'Resolving the current alignment.';
    if (unresolvedCount.value > 0) {
      return `${unresolvedCount.value} record${unresolvedCount.value === 1 ? '' : 's'} still require Select or Skip.`;
    }
    if (status.value !== 'reviewing') return 'Open an alignment review before applying.';
    return '';
  });

  return {
    draft,
    status,
    busy,
    error,
    summary,
    notice,
    repair,
    drawerReferenceKey,
    activePlanInspector,
    dialogOpen: computed(() => !automaticApply.value && Boolean(draft.value && status.value !== 'idle')),
    isDraftArtifactCurrent: artifactIsCurrent,
    unresolvedCount,
    applyDisabledReason,
    canApply: computed(() => status.value === 'reviewing' && unresolvedCount.value === 0),
    validateBeforeGenerate,
    resetAlignment,
    clearForMutation: clearCommittedPlan,
    retainForStableReorder,
    beforeRecordDrag,
    afterRecordDrag,
    directionPreview,
    setDirectionMode,
    setCustomDirection,
    resetPreview,
    resetDialogOpen,
    resetScope,
    openReset: () => { resetScope.value = 'positions'; resetDialogOpen.value = true; },
    cancelReset: () => { if (!busy.value) resetDialogOpen.value = false; },
    applyReset: () => resetAlignment(resetScope.value),
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
    selectCandidate: (recordKey, anchor) => editRow(recordKey, { kind: 'select', anchor }),
    skipRecord: (recordKey) => editRow(recordKey, { kind: 'skip' }),
    applyDraft,
    cancel,
    previewCandidate: (anchor) => {
      if (typeof previewCandidate === 'function') previewCandidate(anchor);
    },
    clearCandidatePreview: () => {
      if (typeof clearCandidatePreview === 'function') clearCandidatePreview();
    }
  };
};
