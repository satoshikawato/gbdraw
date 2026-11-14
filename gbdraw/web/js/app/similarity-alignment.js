// @ts-check
/** @import { AlignmentResetReceipt } from '../services/session-active-config-contract.js' */
/**
 * @import { CanonicalAlignmentAnchorIdentity, CanonicalAlignmentRecordDecision, CanonicalRecordWindow,
 *   CanonicalSimilarityAlignmentPlan } from '../contracts/request.generated.js'
 */
import { normalizeUserFacingError } from '../utils/error-normalization.js';
/** @import { UserFacingError } from '../utils/error-normalization.js' */
/** @import { DrawnFeatureIdsOf } from './legend-layout/composition-actions.js' */
import { msg, numberText, t, tm, tn } from '../utils/i18n.js';
import { plainTextLinearRecordLabel } from '../services/linear-comparisons.js';
import { linearRecordNames } from '../services/linear-sources.js';
import { recordTranslationsFor } from '../services/linear-record-layout.js';
import { canonicalRecordReverseComplement, cropFormReverseComplement, linearCardCropForm } from '../services/record-display-model.js';
/** @import { LinearSeq, WebState } from '../state.js' */
import { validateSimilarityAlignmentResetReceipt } from '../services/session-active-config-contract.js';
import {
  featureIdentity,
  orthogroupIdStatus
} from '../services/feature-identity.js';
import { CompositionMetadataError, materializeRecordTranslations } from './legend-layout/composition-actions.js';
import { isInternalProteinDisplayId, locationPartBound } from '../services/feature-utils.js';
import { describeDirectionPreview, recordsNeedingChoice, resetPreviewView, shortRecordName } from '../services/alignment-review-view.js';
import {
  anchorsAgree, classifyPlanAgainstCatalog, resetSourceMembers, restoredAxisTranslations, restoredPlacement
} from '../services/alignment-plan-facts.js';
import { cloneJsonData as cloneJson } from '../services/json-clone.js';
import { SIMILARITY_ALIGNMENT_PLAN_SCHEMA } from '../contracts/schema.generated.js';
import { ownEntry } from '../utils/own-entry.js';
import { catalogNameIndex } from '../services/feature-catalog.js';

/**
 * @import { FeatureCatalog } from '../services/feature-catalog.js'
 * @import { AliasSource, CommittedRecord } from '../services/feature-identity.js'
 * @import { OrthogroupGroup, OrthogroupMember } from '../services/orthogroup-feature-metadata.js'
 * @import { CommittedCanonicalSession } from '../services/session-request.js'
 */

const { computed, ref, watch } = window.Vue;

/**
 * What the name helpers read of a catalog feature or an orthogroup member; any other field is ignored.
 * @typedef {{
 *   gene?: unknown, locus_tag?: unknown, locusTag?: unknown, product?: unknown, protein_id?: unknown,
 *   proteinId?: unknown, displayProteinId?: unknown, sourceProteinId?: unknown, start?: unknown,
 *   end?: unknown, strand?: unknown, type?: unknown, qualifiers?: Record<string, unknown> | null,
 *   sequenceFeature?: FeatureFacts | null
 * }} FeatureFacts
 * @typedef {{ label: string, name: string, detail: string }} RecordNameEntry A record's name for the Review panel and the Reset preview.
 * @typedef {import('../services/alignment-review-view.js').DisplayedStrand} DisplayedStrand
 * @typedef {{ groupId?: unknown, reference?: AliasSource | null, source?: string, mode?: string }} AlignmentStartOptions
 * @typedef {ReturnType<typeof reviewRows>[number]} ReviewRow One record of the Review panel.
 * @typedef {{
 *   mode?: string,
 *   records?: {
 *     recordKey?: unknown,
 *     windows: CanonicalRecordWindow[],
 *     presentation: { reverseComplement: boolean }
 *   }[]
 * } | null | undefined} CommittedAlignmentRequest The committed render request, as the helper request reads it.
 * @typedef {AliasSource & {
 *   start?: unknown, end?: unknown, strand?: unknown, hidden?: unknown, representative?: unknown, role?: unknown
 * }} OrthogroupMemberFacts A Similarity Group member.
 * @typedef {{
 *   recordKey: string, anchor: AnchorIdentity | null, status: string,
 *   beforeReverseComplement: boolean, afterReverseComplement: boolean,
 *   beforeArrow: DisplayedStrand, afterArrow: DisplayedStrand, exclusion: string | null,
 *   beforeCenterX: number | null, afterCenterX: number | null, translation: { x: number, y: number }
 * }} ProjectedRecord
 * @typedef {{
 *   status: 'projected', binding: string, records: ProjectedRecord[], selectionComplete: boolean,
 *   reference: { anchor: AnchorIdentity, beforeX: number, afterX: number, deltaX: number },
 *   geometryValidated: boolean, signature: string
 * }} ProjectedDirections
 * @typedef {{ status: 'stale', reason: string, binding: string }} StaleDirections
 * @typedef {{ mode: string, byRecordKey?: Record<string, string> }} DirectionIntent
 *   The Review panel's direction choice; `projectSimilarityAlignmentDirections` validates it.
 * @typedef {CanonicalAlignmentAnchorIdentity} AnchorIdentity The exact feature identity that a plan or a candidate names.
 * @typedef {{
 *   anchor: AnchorIdentity, displayName: string, sourceStart: number, sourceEnd: number,
 *   displayCenter: number | null, displayedStrand: DisplayedStrand, hidden: boolean,
 *   representative: boolean, role: string, usable: boolean, directEvidence: string[],
 *   strandRelation: 'same' | 'opposite' | 'unknown'
 * }} AlignmentCandidate A candidate of the helper response, as `validateCandidate` checks it.
 * @typedef {{
 *   kind: 'decision', recordKey: string, status: 'reference' | 'aligned' | 'skipped', rationale: string,
 *   reviewReason: string | null, anchor: AnchorIdentity | null, candidates: AlignmentCandidate[]
 * }} AlignmentDecision
 * @typedef {{
 *   kind: 'ambiguous', status?: undefined, recordKey: string, candidates: AlignmentCandidate[],
 *   directRbhCandidates: AnchorIdentity[], recommendedAnchor: AnchorIdentity, recommendationReason: string
 * }} AlignmentAmbiguity
 * @typedef {CanonicalAlignmentRecordDecision} AlignmentPlanRecord
 * @typedef {CanonicalSimilarityAlignmentPlan} AlignmentPlan
 * @typedef {{ anchor: AnchorIdentity, centerX: number }} ProjectionBeforeAnchor
 * @typedef {{
 *   anchor: AnchorIdentity, displayedStrand: DisplayedStrand, displayCenter: number | null, centerX: number | null
 * }} ProjectionVariantAnchor
 * @typedef {{ reverseComplement: boolean, axisY: number, anchors: ProjectionVariantAnchor[] }} ProjectionVariant
 * @typedef {{
 *   recordKey: string, beforeReverseComplement: boolean, base: { x: number, y: number }, beforeAxisY: number,
 *   beforeAnchors: ProjectionBeforeAnchor[], variants: ProjectionVariant[]
 * }} ProjectionRecord
 * @typedef {{
 *   binding: string, geometryOrientations: Record<string, boolean>, records: ProjectionRecord[]
 * }} AlignmentProjection Python's source-bound facts for both directions of every record.
 * @typedef {{
 *   schema: 2, status: 'resolved' | 'ambiguous', groupId: string, reference: AnchorIdentity,
 *   referenceDisplayedStrand: DisplayedStrand, referenceDisplayCenter: number,
 *   referenceSourceStart: number, referenceSourceEnd: number,
 *   records: (AlignmentDecision | AlignmentAmbiguity)[], plan: AlignmentPlan | null,
 *   projection: AlignmentProjection
 * }} AlignmentResolution The validated helper response: Python places every candidate on the committed request
 *   that the call names, and returns its projection (D5-13).
 * @typedef {{ recordKey: string, kind: 'select' | 'skip', anchor: AnchorIdentity | null }} AlignmentChoice
 * @typedef {{
 *   groupId: string, anchor: AnchorIdentity,
 *   identityIsUnique: boolean, hidden: boolean, representative: boolean, role: string
 * }} AlignmentMemberFacts A group member as a helper request names it: by its identity. Python reads its source
 *   coordinates and strand from the committed request (Python's `_member_candidate`); the page sends no copy.
 * @typedef {{
 *   recordKey: string, windows: CanonicalRecordWindow[], presentation: { reverseComplement: boolean }
 * }} AlignmentRequestRecord A displayed record of a helper request: its windows and orientation as the committed
 *   request stores them, which the page reads to check the projection. Python's `AlignmentHelperRecord` reads only
 *   the key: the committed request owns each record's windows and source length (D5).
 * @typedef {{
 *   schema: 2, groupId: string, records: AlignmentRequestRecord[], reference: AnchorIdentity,
 *   members: AlignmentMemberFacts[], directEdges: never[], choices: AlignmentChoice[]
 * }} AlignmentHelperRequest What the Worker helper resolves.
 * @typedef {{
 *   canonicalRequest: Record<string, unknown>, orientations: Record<string, boolean> | null, geometry?: boolean
 * }} AlignmentHelperProjection The committed request a helper call places every candidate on (Python's
 *   `AlignmentHelperProjection`): `orientations` are the directions the projection geometry is drawn after (null:
 *   every record flipped); `geometry: false` asks for the decisions only, and the reply's projection is null.
 */

/** @type {Set<string | null>} */
const REVIEW_REASONS = new Set([
  'only_usable_candidate', 'unique_direct_rbh',
  'unique_representative', 'deterministic_candidate_1'
]);
/** @type {Set<AlignmentDecision['status']>} */
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

/**
 * @param {unknown} value
 * @returns {value is Record<string, unknown>}
 */
const isObject = (value) => (
  value !== null && typeof value === 'object' && !Array.isArray(value)
);

/**
 * The value itself, known to be an object with exactly `keys`. Its fields stay unchecked:
 * a validator reads them from the result, so each check narrows the field it reads.
 * @template {string} K
 * @param {unknown} value
 * @param {readonly K[]} keys
 * @param {string} path
 * @returns {Record<K, unknown>}
 */
const exactObject = (value, keys, path) => {
  if (!isObject(value) || Object.keys(value).sort().join('\0') !== [...keys].sort().join('\0')) {
    throw new Error(`${path} has invalid fields.`);
  }
  return value;
};

// The checks of a reply field as type guards: the field has the type once its check has passed.
/**
 * @param {unknown} value
 * @returns {value is number}
 */
const isSafeInteger = (value) => Number.isSafeInteger(value);

/**
 * @param {unknown} value
 * @returns {value is number}
 */
const isFiniteNumber = (value) => Number.isFinite(value);

/**
 * @template {string | number | null} T
 * @param {readonly T[]} values
 * @param {unknown} value
 * @returns {value is T}
 */
const oneOf = (values, value) => /** @type {readonly unknown[]} */ (values).includes(value);

/**
 * @template T
 * @param {ReadonlySet<T>} values
 * @param {unknown} value
 * @returns {value is T}
 */
const inSet = (values, value) => /** @type {ReadonlySet<unknown>} */ (values).has(value);

/**
 * @param {unknown} value
 * @param {string} path
 */
const text = (value, path) => {
  if (typeof value !== 'string' || !value.trim() || value.includes('\0')) {
    throw new Error(`${path} must be non-empty text without NUL.`);
  }
  return value.trim();
};

/**
 * @param {unknown} value
 * @param {string} path
 * @returns {number | null}
 */
const nullableInteger = (value, path) => {
  if (value === null) return null;
  if (!isSafeInteger(value) || value < 0) {
    throw new Error(`${path} must be a non-negative integer or null.`);
  }
  return value;
};

/**
 * @param {unknown} left
 * @param {unknown} right
 */
const sameJson = (left, right) => JSON.stringify(left) === JSON.stringify(right);

/** @param {unknown} anchor */
const anchorKey = (anchor) => JSON.stringify(anchor);

/** @param {unknown} strand */
const strandLabel = (strand) => (
  strand === 1 ? '+' : strand === -1 ? '-' : t('unknown')
);

/** @type {Readonly<Record<string, string>>} */
const rationaleLabels = Object.freeze({
  reference: msg('Exact reference'),
  user_selected: msg('Selected by user'),
  only_usable_candidate: msg('Only usable candidate'),
  unique_direct_rbh: msg('Unique direct RBH'),
  skipped_by_user: msg('Skipped by user'),
  skipped_no_candidate: msg('No usable candidate'),
  skipped_unmappable: msg('Candidate center outside the displayed crop')
});

// Direct evidence kinds between a candidate and the exact reference (Python edge kinds).
/** @type {Readonly<Record<string, string>>} */
const evidenceLabels = Object.freeze({
  rbh: msg('reciprocal best hit'),
  coortholog: msg('near-reciprocal hit'),
  record_local_paralog: msg('within-record paralog hit')
});

// Why a candidate is recommended, as a clause after its name ("Recommended: livU, …").
/** @type {Readonly<Record<string, string>>} */
const reviewReasonLabels = Object.freeze({
  only_usable_candidate: msg('the only usable candidate in this record'),
  unique_direct_rbh: msg('the only candidate with a reciprocal best hit'),
  unique_representative: msg('the only group representative among the candidates'),
  deterministic_candidate_1: msg('the first candidate in a fixed order; none stands out')
});

/**
 * @template T
 * @param {T} value
 * @returns {T}
 */
const deepFreeze = (value) => {
  if (!value || typeof value !== 'object' || Object.isFrozen(value)) return value;
  Object.values(value).forEach(deepFreeze);
  return Object.freeze(value);
};

/** @param {unknown[]} values */
const displayText = (...values) => values
  .flatMap((value) => Array.isArray(value) ? value : [value])
  .map((value) => String(value ?? '').trim())
  .find((value) => value && !isInternalProteinDisplayId(value)) || '';

/**
 * The facts of a feature's source: an enriched member's catalog feature (`sequenceFeature`, which the orthogroup owner
 * finds by identity), else the facts given. A member's own copy may be stale: a Session an older gbdraw wrote saved
 * its strand in the drawn direction (OV-3442).
 * @param {FeatureFacts} facts
 * @returns {FeatureFacts}
 */
const sourceFeature = (facts) => facts.sequenceFeature || facts;

// Human-readable names shared by Review, the direction and Reset previews, and
// the active-plan inspector. Record keys and feature hashes stay binding-only.
/**
 * @param {{ biologicalFeatureId?: string }} anchor
 * @param {FeatureFacts} [facts]
 * @param {string} [fallback]
 */
export const featureNameParts = (anchor, facts = {}, fallback = '') => {
  const feature = sourceFeature(facts);
  const qualifiers = feature.qualifiers || {};
  const gene = displayText(feature.gene, facts.gene, qualifiers.gene);
  const locusTag = displayText(feature.locus_tag, feature.locusTag,
    facts.locus_tag, facts.locusTag, qualifiers.locus_tag);
  const proteinId = displayText(facts.displayProteinId, feature.protein_id, feature.proteinId,
    facts.proteinId, facts.sourceProteinId, qualifiers.protein_id);
  const product = displayText(feature.product, facts.product, qualifiers.product);
  const name = gene || product || locusTag || proteinId
    || (fallback === anchor.biologicalFeatureId ? '' : displayText(fallback));
  if (name) {
    const ids = [locusTag, proteinId].filter((value) => value && value !== name);
    return { name, ids, label: [name, ...ids].join(' · ') };
  }
  const start = locationPartBound(feature.start == null ? facts : feature, 'start');
  const end = locationPartBound(feature.end == null ? facts : feature, 'end');
  const coordinates = Number.isSafeInteger(start) && Number.isSafeInteger(end)
    ? `${numberText(start + 1)}..${numberText(end)} bp` : '';
  const typeName = [displayText(feature.type, facts.type) || t('Feature'), coordinates].filter(Boolean).join(' ');
  return { name: typeName, ids: [], label: typeName };
};

// The joined label names a feature in accessible names; the parts lay it out.
/**
 * @param {{ biologicalFeatureId?: string }} anchor
 * @param {Map<string, FeatureFacts>} displayFacts
 * @param {string} [fallback]
 */
const featureView = (anchor, displayFacts, fallback = '') => {
  const parts = featureNameParts(anchor, displayFacts.get(anchorKey(anchor)) || {}, fallback);
  return { label: parts.label, name: parts.name, ids: parts.ids.join(' · ') };
};

// Presentation label or definition plus accession; the display position is
// appended only to tell otherwise identical records apart.
/**
 * @param {{
 *   recordKeys: string[],
 *   request?: Record<string, any> | null,
 *   linearSeqs?: readonly LinearSeq[],
 *   catalog?: FeatureCatalog | null
 * }} options
 */
const recordNames = ({ recordKeys, request = null, linearSeqs = [], catalog = null }) => {
  const { rowsByRecord } = catalogNameIndex(catalog);
  // A record's accession: the first of its catalog features that names one.
  /** @param {string} recordKey */
  const catalogAccession = (recordKey) => {
    for (const feature of rowsByRecord.get(recordKey) || []) {
      const id = displayText(feature?.record_id, feature?.recordId);
      if (id) return id;
    }
    return '';
  };
  // A record without its own Definition takes its File's (OV-436). The accession and
  // the label stay this panel's own: the accession is always added and is the short name.
  const fileNames = linearRecordNames(Array.isArray(linearSeqs) ? linearSeqs : []);
  const parts = recordKeys.map((recordKey, index) => {
    const record = request?.records?.find((/** @type {{ recordKey?: unknown, presentation?: { label?: unknown } | null } | null} */ entry) => entry?.recordKey === recordKey);
    const sequence = (Array.isArray(linearSeqs) ? linearSeqs : [])
      .find((entry) => String(entry?.uid || '') === recordKey);
    const definition = displayText(record?.presentation?.label, fileNames.get(recordKey)?.definition);
    const name = definition ? plainTextLinearRecordLabel(definition) : '';
    // The catalog is read only for a record whose File names no accession.
    const accession = displayText(sequence?.accession) || catalogAccession(recordKey);
    const fileName = displayText(sequence?.gb?.name, sequence?.gff?.name);
    const label = [...new Set([name, accession].filter(Boolean))].join(' · ')
      || (fileName ? plainTextLinearRecordLabel(fileName) : '')
      || t('Record {n}', { n: index + 1 });
    return { label, name: name || accession || label, accession: name && accession !== name ? accession : '' };
  });
  /** @type {Map<string, number>} */
  const counts = new Map();
  parts.forEach(({ label }) => counts.set(label, (counts.get(label) || 0) + 1));
  // Each entry keeps the joined label for accessible names and its parts for layout.
  return new Map(recordKeys.map((recordKey, index) => {
    const { label, name, accession } = parts[index];
    // Every label was counted above.
    const duplicate = /** @type {number} */ (counts.get(label)) > 1 ? t('Record {n}', { n: index + 1 }) : '';
    return [recordKey, {
      label: duplicate ? `${label} · ${duplicate}` : label,
      name,
      detail: [accession, duplicate].filter(Boolean).join(' · ')
    }];
  }));
};

/**
 * @param {AlignmentCandidate} candidate
 * @param {AlignmentHelperRequest} request
 * @param {Map<string, FeatureFacts>} [displayFacts]
 */
const candidateView = (candidate, request, displayFacts = new Map()) => {
  const key = anchorKey(candidate.anchor);
  const member = request.members.find(({ anchor }) => anchorKey(anchor) === key);
  if (!member) throw new Error('Alignment candidate has no current member facts.');
  const coordinates = `${numberText(candidate.sourceStart + 1)}..${numberText(candidate.sourceEnd)} bp`;
  const feature = featureView(candidate.anchor, displayFacts, candidate.displayName);
  return {
    key,
    anchor: candidate.anchor,
    label: feature.label,
    name: feature.name,
    ids: feature.ids,
    // A kind is helper text: one the table does not hold is shown as written (OV-573).
    evidence: candidate.directEvidence.map((kind) => Object.hasOwn(evidenceLabels, kind) ? tm(evidenceLabels[kind]) : kind),
    coordinates,
    displayedStrand: strandLabel(candidate.displayedStrand),
    representative: candidate.representative,
    role: candidate.role || 'member',
    directEvidence: candidate.directEvidence.length ? candidate.directEvidence : ['None'],
    strandRelation: candidate.strandRelation,
    usable: candidate.usable
  };
};

/**
 * @param {AlignmentResolution} response
 * @param {AlignmentHelperRequest} request
 * @param {Map<string, FeatureFacts>} displayFacts
 * @param {Map<string, RecordNameEntry>} recordLabels
 */
const referenceView = (response, request, displayFacts, recordLabels) => {
  if (!request.members.some(({ anchor }) => sameJson(anchor, response.reference))) {
    throw new Error('Alignment reference has no current member facts.');
  }
  const feature = featureView(response.reference, displayFacts);
  const record = recordLabels.get(response.reference.recordKey);
  return {
    label: feature.label,
    name: feature.name,
    ids: feature.ids,
    recordLabel: record?.label,
    recordName: record?.name,
    recordDetail: record?.detail,
    recordKey: response.reference.recordKey,
    // Python's source coordinates of the reference, read from the committed request.
    coordinates: numberText(response.referenceSourceStart + 1) + '..'
      + numberText(response.referenceSourceEnd) + ' bp',
    displayedStrand: strandLabel(response.referenceDisplayedStrand)
  };
};

/**
 * @param {AlignmentResolution} response
 * @param {AlignmentHelperRequest} request
 * @param {Map<string, FeatureFacts>} displayFacts
 * @param {Map<string, RecordNameEntry>} recordLabels
 */
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
    const names = recordLabels.get(record.recordKey);
    return {
      recordKey: record.recordKey,
      recordLabel: names?.label,
      recordName: names?.name,
      recordDetail: names?.detail,
      candidates: record.candidates.filter(({ usable }) => usable)
        .map((candidate) => candidateView(candidate, request, displayFacts)),
      reason,
      reasonLabel: reason
        ? tm(reviewReasonLabels[reason] || rationaleLabels[reason]) : t('Unchanged'),
      recommendedKey,
      recommendationReasonLabel: recommendedKey
        // `recommendedKey` is set only when `REVIEW_REASONS` has the reason, so it is text.
        ? tm(reviewReasonLabels[/** @type {string} */ (reason)]) : null,
      choice: record.kind === 'ambiguous' ? null : selectedAnchor
        ? { kind: 'select', candidateKey: anchorKey(selectedAnchor) }
        : { kind: 'skip', candidateKey: null },
      unchanged: !selectedAnchor,
      repairRequired: false
    };
  });

/**
 * @param {AlignmentPlan} plan
 * @param {number} reversed
 */
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
    text: t('Alignment applied: {aligned} aligned, {unchanged} unchanged, {explicitlySkipped} explicitly skipped, {noCandidate} without a candidate, {reversed} reversed.',
      { aligned, unchanged, explicitlySkipped, noCandidate, reversed })
  });
};

const baseReverseComplement = canonicalRecordReverseComplement;

/**
 * @param {AlignmentPlan | null} plan
 * @param {{
 *   linearSeqs?: LinearSeq[],
 *   catalog?: FeatureCatalog | null,
 *   request?: Record<string, any> | null
 * }} [options]
 */
const inspectActivePlan = (plan, { linearSeqs = [], catalog = null, request = null } = {}) => {
  if (!plan || plan.schema !== SIMILARITY_ALIGNMENT_PLAN_SCHEMA || !Array.isArray(plan.records)) return null;
  const records = new Map(
    (Array.isArray(linearSeqs) ? linearSeqs : [])
      .map((record) => [record.uid, cropFormReverseComplement(linearCardCropForm(record))])
  );
  const factsOf = catalogNameIndex(catalog).featureOf;
  const nameOf = (/** @type {{ recordKey?: string, biologicalFeatureId?: string }} */ anchor) => featureNameParts(anchor, factsOf(anchor) || {});
  const names = recordNames({ recordKeys: plan.records.map(({ recordKey }) => recordKey),
    request, linearSeqs, catalog });
  const reference = nameOf(plan.reference);
  const referenceRecord = names.get(plan.reference.recordKey);
  return deepFreeze({
    groupId: plan.groupId,
    // The collapsed inspector line: what the alignment holds fixed.
    summary: t('Aligned to {name} in {record}', { name: reference.name, record: shortRecordName(referenceRecord) }),
    reference: {
      ...plan.reference,
      label: reference.label,
      recordLabel: referenceRecord?.label
    },
    records: plan.records.map((decision) => {
      const recordName = names.get(decision.recordKey);
      return {
        recordKey: decision.recordKey,
        recordLabel: recordName?.label,
        recordName: recordName?.name,
        recordDetail: recordName?.detail,
        anchorLabel: decision.anchor ? nameOf(decision.anchor).label : t('Skip'),
        status: decision.status,
        rationale: decision.rationale,
        rationaleLabel: rationaleLabels[decision.rationale] ? tm(rationaleLabels[decision.rationale]) : decision.rationale,
        reversedFromSource: records.get(decision.recordKey) || false
      };
    })
  });
};

/**
 * @param {AliasSource | null | undefined} source
 * @param {string} path
 * @returns {AnchorIdentity}
 */
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

/**
 * @param {unknown} value
 * @param {string} path
 * @returns {AnchorIdentity}
 */
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

/**
 * @param {unknown} value
 * @param {string} path
 * @param {string} recordKey
 * @returns {AlignmentCandidate}
 */
const validateCandidate = (value, path, recordKey) => {
  const raw = exactObject(value, [
    'anchor', 'displayName', 'sourceStart', 'sourceEnd', 'displayCenter',
    'displayedStrand', 'hidden', 'representative', 'role', 'usable',
    'directEvidence', 'strandRelation'
  ], path);
  const anchor = validateAnchor(raw.anchor, `${path}.anchor`);
  if (anchor.recordKey !== recordKey || typeof raw.displayName !== 'string'
    || !isSafeInteger(raw.sourceStart) || !isSafeInteger(raw.sourceEnd)
    || raw.sourceStart < 0 || raw.sourceEnd < raw.sourceStart
    || (raw.displayCenter !== null && !isFiniteNumber(raw.displayCenter))
    || !oneOf([null, -1, 1], raw.displayedStrand)
    || typeof raw.hidden !== 'boolean' || typeof raw.representative !== 'boolean'
    || typeof raw.usable !== 'boolean' || typeof raw.role !== 'string'
    || !Array.isArray(raw.directEvidence)
    || raw.directEvidence.some((/** @type {unknown} */ evidence) => typeof evidence !== 'string')) {
    throw new Error(`${path} contains invalid candidate facts.`);
  }
  if (!oneOf(['same', 'opposite', 'unknown'], raw.strandRelation)) {
    throw new Error(`${path}.strandRelation is invalid.`);
  }
  // The checked fields, in the helper's order.
  return {
    anchor, displayName: raw.displayName, sourceStart: raw.sourceStart, sourceEnd: raw.sourceEnd,
    displayCenter: raw.displayCenter, displayedStrand: raw.displayedStrand, hidden: raw.hidden,
    representative: raw.representative, role: raw.role, usable: raw.usable,
    directEvidence: raw.directEvidence, strandRelation: raw.strandRelation
  };
};

/**
 * @param {unknown} value
 * @param {string} path
 * @returns {AlignmentDecision}
 */
const validateDecision = (value, path) => {
  const raw = exactObject(value, [
    'kind', 'recordKey', 'status', 'rationale', 'reviewReason', 'anchor',
    'candidates'
  ], path);
  if (raw.kind !== 'decision' || !inSet(STATUSES, raw.status)
    || !inSet(RATIONALES, raw.rationale)
    || (raw.reviewReason !== null && !inSet(REVIEW_REASONS, raw.reviewReason))
    || !Array.isArray(raw.candidates)) {
    throw new Error(`${path} contains an unknown decision value.`);
  }
  const recordKey = text(raw.recordKey, `${path}.recordKey`);
  const anchor = raw.anchor === null ? null : validateAnchor(raw.anchor, `${path}.anchor`);
  /** @type {AlignmentCandidate[]} */
  const candidates = raw.candidates.map((/** @type {unknown} */ candidate, /** @type {number} */ index) => (
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
  return {
    kind: raw.kind,
    // The reply's own spelling, which `text` accepted above.
    recordKey: /** @type {string} */ (raw.recordKey),
    status: raw.status, rationale: raw.rationale, reviewReason: raw.reviewReason, anchor, candidates
  };
};

/**
 * @param {unknown} value
 * @param {string} path
 * @returns {AlignmentAmbiguity}
 */
const validateAmbiguity = (value, path) => {
  const raw = exactObject(value, [
    'kind', 'recordKey', 'candidates', 'directRbhCandidates',
    'recommendedAnchor', 'recommendationReason'
  ], path);
  if (raw.kind !== 'ambiguous' || !Array.isArray(raw.candidates)
    || raw.candidates.length < 2 || !oneOf(['unique_representative', 'deterministic_candidate_1'],
      raw.recommendationReason)) throw new Error(`${path} is not a valid ambiguity.`);
  const recordKey = text(raw.recordKey, `${path}.recordKey`);
  /** @type {AlignmentCandidate[]} */
  const candidates = raw.candidates.map((/** @type {unknown} */ candidate, /** @type {number} */ index) => (
    validateCandidate(candidate, `${path}.candidates[${index}]`, recordKey)
  ));
  const keys = new Set(candidates.map(({ anchor }) => anchorKey(anchor)));
  const recommendedAnchor = validateAnchor(raw.recommendedAnchor, `${path}.recommendedAnchor`);
  if (keys.size !== candidates.length || !candidates.some((candidate) => (
    candidate.usable && sameJson(candidate.anchor, recommendedAnchor)
  )) || !Array.isArray(raw.directRbhCandidates)) {
    throw new Error(`${path} contains invalid candidates or recommendation.`);
  }
  /** @type {AnchorIdentity[]} */
  const directRbhCandidates = raw.directRbhCandidates.map((/** @type {unknown} */ anchor, /** @type {number} */ index) => (
    validateAnchor(anchor, `${path}.directRbhCandidates[${index}]`)
  ));
  if (directRbhCandidates.some((anchor) => !keys.has(anchorKey(anchor)))) {
    throw new Error(`${path} contains unknown direct evidence.`);
  }
  return { kind: raw.kind, recordKey, candidates, directRbhCandidates, recommendedAnchor,
    recommendationReason: raw.recommendationReason };
};

/**
 * @param {unknown} value
 * @param {Pick<AlignmentResolution, 'groupId' | 'reference' | 'records'>} response
 * @param {string} path
 * @returns {AlignmentPlan}
 */
const validatePlan = (value, response, path) => {
  const raw = exactObject(value, ['schema', 'groupId', 'reference', 'records'], path);
  if (raw.schema !== SIMILARITY_ALIGNMENT_PLAN_SCHEMA || raw.groupId !== response.groupId || !Array.isArray(raw.records)
    || !sameJson(validateAnchor(raw.reference, `${path}.reference`), response.reference)) {
    throw new Error(`${path} does not match its resolution.`);
  }
  /** @type {AlignmentPlanRecord[]} */
  const records = raw.records.map((/** @type {unknown} */ record, /** @type {number} */ index) => {
    const decision = response.records[index];
    const fields = exactObject(record, [
      'recordKey', 'status', 'rationale', 'anchor'
    ], `${path}.records[${index}]`);
    if (!decision || decision.kind !== 'decision'
      || !sameJson(fields, {
        recordKey: decision.recordKey, status: decision.status,
        rationale: decision.rationale, anchor: decision.anchor
      })) throw new Error(`${path} differs from helper decisions.`);
    // A reply is JSON, and this record equals the checked decision's four fields as JSON.
    return /** @type {AlignmentPlanRecord} */ (fields);
  });
  return { schema: SIMILARITY_ALIGNMENT_PLAN_SCHEMA, groupId: raw.groupId,
    // `validateAnchor` accepted this anchor above; the plan keeps the reply's own object.
    reference: /** @type {AnchorIdentity} */ (raw.reference), records };
};

/**
 * @param {unknown} value
 * @param {AlignmentHelperRequest} request
 * @returns {AlignmentProjection}
 */
const validateProjection = (value, request) => {
  const raw = exactObject(value, ['binding', 'geometryOrientations', 'records'], 'alignment projection');
  const keys = request.records.map(({ recordKey }) => recordKey);
  if (typeof raw.binding !== 'string' || !/^[a-f0-9]{64}$/.test(raw.binding) || !Array.isArray(raw.records)
    || raw.records.length !== keys.length) throw new Error('Invalid alignment projection binding or coverage.');
  const orientations = exactObject(raw.geometryOrientations, keys, 'alignment projection.geometryOrientations');
  const known = new Set(request.members.map(({ anchor }) => anchorKey(anchor)));
  const finite = (/** @type {unknown} */ number) => { if (!Number.isFinite(number)) throw new Error('Nonfinite alignment placement.'); };
  raw.records.forEach((/** @type {unknown} */ recordValue, /** @type {number} */ index) => {
    const record = exactObject(recordValue, ['recordKey', 'beforeReverseComplement', 'base', 'beforeAxisY',
      'beforeAnchors', 'variants'], 'alignment projection record');
    if (record.recordKey !== keys[index] || typeof orientations[record.recordKey] !== 'boolean'
      || record.beforeReverseComplement !== baseReverseComplement(request.records[index])) {
      throw new Error('Alignment projection changed record orientation binding.');
    }
    const base = exactObject(record.base, ['x', 'y'], 'alignment projection base');
    [base.x, base.y, record.beforeAxisY].forEach(finite);
    /**
     * @param {unknown} anchors
     * @param {boolean} variant
     */
    const validateAnchors = (anchors, variant) => {
      if (!Array.isArray(anchors)) throw new Error('Alignment projection anchors must be an array.');
      const identities = new Set();
      anchors.forEach((/** @type {unknown} */ entryValue) => {
        const entry = exactObject(entryValue, variant ? ['anchor', 'displayedStrand', 'displayCenter', 'centerX']
          : ['anchor', 'centerX'], 'alignment projection anchor');
        const anchor = validateAnchor(entry.anchor, 'alignment projection anchor identity');
        const key = anchorKey(anchor);
        if (anchor.recordKey !== record.recordKey || !known.has(key) || identities.has(key)) {
          throw new Error('Alignment projection changed anchor coverage.');
        }
        identities.add(key);
        if (variant) {
          if (!oneOf([null, -1, 1], entry.displayedStrand)
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
    record.variants.forEach((/** @type {unknown} */ variantValue, /** @type {number} */ direction) => {
      const variant = exactObject(variantValue, ['reverseComplement', 'axisY', 'anchors'], 'alignment projection variant');
      if (variant.reverseComplement !== Boolean(direction)) throw new Error('Invalid alignment direction variant.');
      finite(variant.axisY);
      const anchors = validateAnchors(variant.anchors, true);
      const expected = request.members.filter(({ anchor }) => anchor.recordKey === record.recordKey);
      if (anchors.size !== expected.length || expected.some(({ anchor }) => !anchors.has(anchorKey(anchor)))) {
        throw new Error('Alignment projection omitted a source anchor.');
      }
      if (variant.reverseComplement === record.beforeReverseComplement
        // `validateAnchors(variant.anchors, true)` above checked the list and each entry.
        && /** @type {ProjectionVariantAnchor[]} */ (variant.anchors).some(({ anchor, centerX }) => (centerX !== null) !== beforeAnchors.has(anchorKey(anchor)))) {
        throw new Error('Alignment projection changed current anchor mappability.');
      }
    });
  });
  // The reply's own object: the checks above cover every field of `AlignmentProjection`.
  return /** @type {AlignmentProjection} */ (raw);
};

// One transient resolver is shared by local previews and final helper admission.
/**
 * @param {{
 *   resolution: AlignmentResolution, intent: DirectionIntent, expectedBinding: unknown,
 *   choices?: AlignmentChoice[]
 * }} options
 * @returns {ProjectedDirections | StaleDirections}
 */
export const projectSimilarityAlignmentDirections = ({ resolution, intent, expectedBinding, choices = [] }) => {
  const facts = resolution?.projection;
  if (!facts) throw new Error('Source-bound alignment projection facts are required.');
  if (typeof expectedBinding !== 'string') throw new Error('An explicit alignment binding is required.');
  if (facts.binding !== expectedBinding) return deepFreeze({ status: 'stale', reason: 'binding_changed', binding: facts.binding });
  const mode = intent?.mode;
  if (!['keep', 'right', 'left', 'custom'].includes(mode)) throw new Error('Invalid alignment direction mode.');
  exactObject(intent, mode === 'custom' ? ['mode', 'byRecordKey'] : ['mode'], 'alignment direction intent');
  // Read only in Custom mode, where it is checked as an object and every entry below.
  const byRecordKey = /** @type {Record<string, string>} */ (intent.byRecordKey);
  if (mode === 'custom') {
    if (!isObject(intent.byRecordKey) || Object.entries(byRecordKey).some(([key, value]) => (
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
    // Without an exclusion, `selected` is a usable anchor (checked above).
    if (!exclusion && /** @type {ProjectionVariantAnchor} */ (selected).displayedStrand === null) exclusion = 'unknown_strand';
    // A record key is free text from a Session (`constructor` must not read `Object.prototype`).
    const requested = mode === 'custom'
      ? (Object.hasOwn(byRecordKey, fact.recordKey) && byRecordKey[fact.recordKey]) || 'keep' : mode;
    const flip = !exclusion && requested !== 'keep'
      && /** @type {ProjectionVariantAnchor} */ (selected).displayedStrand !== (requested === 'right' ? 1 : -1);
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
  // The reference record is among the records, and its centers are checked finite next.
  const reference = /** @type {ProjectedRecord & { beforeCenterX: number, afterCenterX: number }} */ (
    records.find(({ recordKey }) => recordKey === resolution.reference.recordKey));
  const fact = /** @type {ProjectionRecord} */ (facts.records.find(({ recordKey }) => recordKey === reference.recordKey));
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

/**
 * The projection of the check before a Generate with a plan: the committed request its decisions are placed on, and
 * no geometry (the reply's projection is null).
 * @param {Record<string, unknown>} canonicalRequest
 * @returns {AlignmentHelperProjection}
 */
const decisionsOnlyProjection = (canonicalRequest) => ({ canonicalRequest, orientations: null, geometry: false });

/**
 * The fields of a helper reply.
 * @param {unknown} value
 */
const helperReplyFields = (value) => exactObject(value, [
  'schema', 'status', 'groupId', 'reference', 'referenceDisplayedStrand',
  'referenceDisplayCenter', 'referenceSourceStart', 'referenceSourceEnd', 'records', 'plan', 'projection'
], 'alignment helper response');

/**
 * The decisions of a helper reply, checked against its request: every validator of a reply reads them.
 * @param {ReturnType<typeof helperReplyFields>} raw
 * @param {AlignmentHelperRequest} request
 * @returns {Omit<AlignmentResolution, 'projection'>}
 */
const validateDecisions = (raw, request) => {
  if (raw.schema !== SIMILARITY_ALIGNMENT_PLAN_SCHEMA || !oneOf(['resolved', 'ambiguous'], raw.status)
    || raw.groupId !== request.groupId || !Array.isArray(raw.records)
    || !oneOf([null, -1, 1], raw.referenceDisplayedStrand)
    || !isFiniteNumber(raw.referenceDisplayCenter)
    || !isSafeInteger(raw.referenceSourceStart) || !isSafeInteger(raw.referenceSourceEnd)
    || raw.referenceSourceStart < 0 || raw.referenceSourceEnd < raw.referenceSourceStart) {
    throw new Error('Alignment helper response does not match its request.');
  }
  const reference = validateAnchor(raw.reference, 'alignment helper response.reference');
  if (!sameJson(reference, request.reference)) {
    throw new Error('Alignment helper response changed the exact reference.');
  }
  // Any value can be asked for its `kind`; only a record of a known kind goes on to its validator.
  /** @type {(AlignmentDecision | AlignmentAmbiguity)[]} */
  const records = raw.records.map((/** @type {{ kind?: unknown } | null | undefined} */ record, /** @type {number} */ index) => (
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
  /** @type {Omit<AlignmentResolution, 'projection'>} */
  const response = {
    schema: SIMILARITY_ALIGNMENT_PLAN_SCHEMA, status: raw.status, groupId: raw.groupId, reference,
    referenceDisplayedStrand: raw.referenceDisplayedStrand,
    referenceDisplayCenter: raw.referenceDisplayCenter,
    referenceSourceStart: raw.referenceSourceStart, referenceSourceEnd: raw.referenceSourceEnd,
    records, plan: /** @type {AlignmentPlan | null} */ (null)
  };
  if (raw.status === 'resolved') {
    if (raw.plan === null) throw new Error('Resolved alignment helper response has no plan.');
    response.plan = validatePlan(raw.plan, response, 'alignment helper response.plan');
  } else if (raw.plan !== null) {
    throw new Error('Ambiguous alignment helper response must not contain a plan.');
  }
  return response;
};

/**
 * A review's reply: the decisions, and the projection whose current placements agree with them.
 * @param {unknown} value
 * @param {AlignmentHelperRequest} request
 * @returns {AlignmentResolution}
 */
export const validateSimilarityAlignmentResolution = (value, request) => {
  const raw = helperReplyFields(value);
  const decided = validateDecisions(raw, request);
  const { records, reference } = decided;
  /** @type {AlignmentResolution} */
  const response = { ...decided, projection: validateProjection(raw.projection, request) };
  response.projection.records.forEach((fact, index) => {
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
  return deepFreeze(response);
};

/**
 * The reply of the check before a Generate with a plan, which asks for no projection geometry (`geometry: false`):
 * the decisions only.
 * @param {unknown} value
 * @param {AlignmentHelperRequest} request
 */
export const validateSimilarityAlignmentCheck = (value, request) => {
  const raw = helperReplyFields(value);
  if (raw.projection !== null) throw new Error('The alignment check reply must not contain a projection.');
  validateDecisions(raw, request);
};

/**
 * @param {unknown} value
 * @param {string} path
 * @returns {DisplayedStrand}
 */
const strandValue = (value, path) => {
  if (value === 1 || value === '+' || String(value).trim() === '1') return 1;
  if (value === -1 || value === '-' || String(value).trim() === '-1') return -1;
  if (value === null || value === undefined
    || ['', 'undefined', 'unstranded', 'mixed'].includes(String(value).trim())) return null;
  throw new Error(`${path} has an invalid source strand.`);
};

/**
 * A feature's source location as an option names it, `start..end (strand)`, from the facts its name reads.
 * @param {FeatureFacts} facts
 */
const sourceLocationText = (facts) => {
  const feature = sourceFeature(facts);
  const start = locationPartBound(feature, 'start');
  const end = locationPartBound(feature, 'end');
  return `${Number.isSafeInteger(start) ? numberText(start + 1) : ''}..`
    + `${Number.isSafeInteger(end) ? numberText(end) : ''} `
    + `(${strandLabel(strandValue(feature.strand, 'Drawer reference strand'))})`;
};

/**
 * The baseline a repair starts from: the records where the shown Result draws them. A Result an older gbdraw drew has
 * no record keys to read a stale plan's placements by; the repair then starts from the committed placements, which
 * the next Generate draws (OV-3440; PD-OI-028: a stale plan is repaired, the Result kept).
 * @template T
 * @param {(options: { materializePlan: boolean, fromResult?: boolean }) => T} baseline
 * @returns {T}
 */
const repairBaseline = (baseline) => {
  try {
    return baseline({ materializePlan: true });
  } catch (cause) {
    if (!(cause instanceof CompositionMetadataError)) throw cause;
    return baseline({ materializePlan: true, fromResult: false });
  }
};

/**
 * The materialized Session for an operation that replaces every record placement, such as Reset record positions: a
 * Result an older gbdraw drew gives the committed placements, as a repair does (`repairBaseline`, OV-3538).
 * @param {(options: { materializePlan: boolean, fromResult?: boolean }) => { canonical: unknown }} baseline
 * @returns {() => CommittedCanonicalSession}
 */
const repairBaseSessionReader = (baseline) => () => /** @type {CommittedCanonicalSession} */ (
  repairBaseline(baseline).canonical);

/**
 * @param {{
 *   group: AliasSource | null | undefined, members: OrthogroupMemberFacts[],
 *   reference: AliasSource | null | undefined, request: CommittedAlignmentRequest, choices: AlignmentChoice[]
 * }} options
 * @returns {AlignmentHelperRequest}
 */
const buildHelperRequest = ({ group, members, reference, request, choices }) => {
  const groupStatus = orthogroupIdStatus(group);
  if (!groupStatus.valid || !groupStatus.supplied) {
    throw new Error('The selected Similarity Group identity is invalid.');
  }

  if (!request || request.mode !== 'linear' || !Array.isArray(request.records)) {
    throw new Error('Generate a Linear diagram before aligning a Similarity Group.');
  }
  const records = request.records.map((record, index) => {
    const recordKey = text(record?.recordKey, `renderRequest.records[${index}].recordKey`);
    // The windows go as the request stores them; the page checks the projection against them.
    return {
      recordKey,
      windows: record.windows.map(({ start, end, reverseComplement, column }) => ({ start, end, reverseComplement, column })),
      presentation: { reverseComplement: record.presentation.reverseComplement }
    };
  });
  const displayedRecords = new Set(records.map(({ recordKey }) => recordKey));
  const normalizedMembers = members.map((member, index) => {
    const anchor = anchorFromSource(member, `Similarity Group member ${index + 1}`);
    if (!displayedRecords.has(anchor.recordKey)) {
      throw new Error(`Similarity Group member ${index + 1} has invalid display facts.`);
    }
    return {
      source: member,
      canonicalKey: `${anchor.recordKey}\0${anchor.biologicalFeatureId}`,
      payload: {
        groupId: groupStatus.value,
        anchor,
        identityIsUnique: true,
        hidden: Boolean(member?.hidden),
        representative: Boolean(member?.representative),
        role: String(member?.role || '')
      }
    };
  });
  /** @type {Map<string, number>} */
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
  // `identityCounts` holds `referenceKey` exactly once, so the member is found.
  const exactMember = /** @type {(typeof normalizedMembers)[number]} */ (
    normalizedMembers.find(({ canonicalKey }) => canonicalKey === referenceKey));
  if (!anchorsAgree(exactMember.payload.anchor, exactReference)) {
    throw new Error('The reference identity conflicts with current group metadata.');
  }

  return deepFreeze({
    schema: SIMILARITY_ALIGNMENT_PLAN_SCHEMA,
    groupId: groupStatus.value,
    records,
    reference: exactMember.payload.anchor,
    members: normalizedMembers.map(({ payload }) => payload),
    // Direct evidence comes from the committed orthogroup resource in the Worker.
    directEdges: [],
    choices: cloneJson(choices)
  });
};

/**
 * @param {SimilarityAlignmentActionsOptions['state']} state
 * @param {string[]} recordKeys
 */
const replacementTranslations = (state, recordKeys) => recordTranslationsFor(
  state.linearRecordTranslations.value,
  recordKeys
);


/**
 * One alignment run that the composition root executes as a candidate (R13):
 * the root commits `orientations` to the record display intent.
 * @typedef {object} AlignmentCandidateRun
 * @property {Record<string, any>} canonical canonical Session projected by the alignment
 * @property {string} operation
 * @property {string} label
 * @property {Record<string, any> | null} alignmentResetBefore
 * @property {AlignmentResetReceipt | null | undefined} alignmentResetReceipt
 * @property {{ recordKey: string, reverseComplement: boolean }[]} orientations
 */

/**
 * The ports and state that the alignment owner receives. Ports are function
 * types; the shapes of requests and helper payloads belong to their owners
 * (the committed request to Python, R7) and stay `Record<string, any>`.
 * @typedef {object} SimilarityAlignmentActionsOptions
 * @property {Pick<WebState, 'activeDrawing' | 'featureCatalog' | 'linearRecordTranslations' | 'mode'
 *   | 'selectedOrthogroupId' | 'similarityAlignmentPlan' | 'similarityAlignmentResetReceipt'>} state
 *   the Web state members the alignment actions read and write
 * @property {(groupId: string) => OrthogroupGroup | null} getOrthogroupById
 * @property {(group: OrthogroupGroup | string | null) => OrthogroupMember[]} getEnrichedOrthogroupMembers
 * @property {DrawnFeatureIdsOf} drawnFeatureIdsOf the drawn element IDs of an anchor's feature (the orthogroup owner's
 *   catalog row lookup), which bind a Result an older gbdraw drew to the plan's records
 * @property {{ readonly value: SelectedGroupDetail | null }} selectedOrthogroupDetail
 *   the orthogroup owner's detail of the drawer's selected group (OV-920)
 * @property {() => Record<string, any> | null} getCommittedRequest
 * @property {() => Record<string, any> | null} getCommittedSession
 * @property {(args: Record<string, any>) => Record<string, any>} projectCommittedAlignment
 * @property {(run: AlignmentCandidateRun) => Promise<Record<string, any>>} runRecordAlignment
 * @property {(() => void) | null} [cancelRunAnalysis]
 * @property {(operation: string, payload: Record<string, any>) => Promise<Record<string, any>>} runHelperOperation
 * @property {string} resolveOperation
 * @property {(() => SVGSVGElement | null) | null} [getCurrentSvg]
 * @property {((anchor: Record<string, any>) => any) | null} [previewCandidate]
 * @property {(() => void) | null} [clearCandidatePreview]
 * @property {((error: UserFacingError) => void) | null} [onError]
 */

/**
 * The orthogroup owner's detail of the drawer's selected group (OV-920): the
 * group, the id its controls pass, and its enriched members.
 * @typedef {object} SelectedGroupDetail
 * @property {NonNullable<ReturnType<SimilarityAlignmentActionsOptions['getOrthogroupById']>>} group
 * @property {unknown} id
 * @property {ReturnType<SimilarityAlignmentActionsOptions['getEnrichedOrthogroupMembers']>} members
 */

/**
 * One committed alignment baseline: the canonical Session that `projectCommittedAlignment`
 * returned, its render request, and the record translations it carries.
 * @typedef {object} AlignmentBaseline
 * @property {Record<string, any>} canonical
 * @property {any} request
 * @property {any} translations
 */

/**
 * What `captureArtifact` records to tell whether the sources of a draft changed: the records,
 * the Similarity Group, the committed resources, the translations. Not the Result or the feature
 * catalog: a live setting change draws a new Result, and the automatic rerender a new catalog,
 * from the same sources (PD-OI-102, OV-1120). The resources object is the identity of those
 * sources: the rerender keeps it (`commitLiveSettingValues`), and a Generate, a Session load or
 * an Undo across one installs another, also from the same files (PD-OI-039). Every path that
 * changes the catalog's features is one of these.
 * @typedef {object} AlignmentArtifactStamp
 * @property {string} records
 * @property {string} group
 * @property {string} selectedGroup
 * @property {unknown} resources
 * @property {string} translations
 */

/**
 * @param {SimilarityAlignmentActionsOptions} options
 */
export const createSimilarityAlignmentActions = ({
  state,
  getOrthogroupById,
  getEnrichedOrthogroupMembers,
  drawnFeatureIdsOf,
  selectedOrthogroupDetail,
  getCommittedRequest,
  getCommittedSession,
  projectCommittedAlignment,
  runRecordAlignment,
  cancelRunAnalysis = null,
  runHelperOperation,
  resolveOperation,
  getCurrentSvg = null,
  previewCandidate = null,
  clearCandidatePreview = null,
  onError = null
}) => {
  // The cards of the shown Linear drawing, and those it draws; none when
  // another drawing is shown. Resolved at each read.
  const shownLinearCards = () => {
    const drawing = state.activeDrawing();
    return drawing.mode === 'linear' ? drawing.linearSeqs : [];
  };
  if (
    !state
    || typeof getOrthogroupById !== 'function'
    || typeof getEnrichedOrthogroupMembers !== 'function'
    || typeof drawnFeatureIdsOf !== 'function'
    || !selectedOrthogroupDetail
    || typeof getCommittedRequest !== 'function'
    || typeof getCommittedSession !== 'function'
    || typeof projectCommittedAlignment !== 'function'
    || typeof runRecordAlignment !== 'function'
    || typeof runHelperOperation !== 'function'
    || typeof resolveOperation !== 'string'
    || !resolveOperation
  ) throw new Error('Similarity alignment actions require their current artifact owners.');

  const draft = ref(null);
  const status = ref('idle');
  const busy = computed(() => status.value === 'resolving' || status.value === 'applying');
  const error = ref(null);
  const summary = ref(null);
  // A "cleared" notice holds only while no plan is active, so an Undo that
  // restores the plan also retires it (OV-382).
  const noticeState = ref({ text: '', whileNoPlan: false });
  const notice = computed(() => (
    noticeState.value.whileNoPlan && state.similarityAlignmentPlan.value ? '' : noticeState.value.text
  ));
  const repair = ref(null);
  // The user's drawer reference pick and the committed plan it was made
  // against. A new plan (Apply, Session load, Undo) drops the pick (OV-381).
  const drawerReferencePick = ref(/** @type {{ key: string, plan: unknown } | null} */ (null));
  const resetDialogOpen = ref(false);
  const automaticApply = ref(false);
  const resetScope = ref('positions');
  let actionId = 0;
  /** @type {AlignmentHelperRequest | null} */
  let activeRequest = null;
  /** @type {Promise<Record<string, any>> | null} */
  let activeApply = null;
  /** @type {AlignmentBaseline | null} */
  let activeBaseline = null;
  /** @type {Map<string, FeatureFacts>} */
  let displayFacts = new Map();
  /** @type {Map<string, RecordNameEntry>} */
  let recordLabels = new Map();
  /** @type {AlignmentArtifactStamp | null} */
  let artifactStamp = null;
  let materializedRecordDrag = false;
  /** @type {AlignmentBaseline | null} */
  let pendingRecordDragBaseline = null;

  /**
   * @param {unknown} value
   * @param {boolean} [notify]
   */
  const publishError = (value, notify = true) => {
    // The value is always truthy (`value || { code: 'UNKNOWN' }`), and only a falsy value normalizes to null.
    const normalized = /** @type {UserFacingError} */ (
      normalizeUserFacingError(value || { code: 'UNKNOWN' }, { operation: 'align', stage: 'helper' })
    );
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
  /** @param {string} groupId */
  const captureArtifact = (groupId) => ({
    records: JSON.stringify(currentRequest()?.records || []),
    group: JSON.stringify(getOrthogroupById(groupId)),
    selectedGroup: state.selectedOrthogroupId.value,
    resources: getCommittedSession()?.resources,
    translations: JSON.stringify(state.linearRecordTranslations.value)
  });
  const artifactIsCurrent = () => {
    if (!artifactStamp || !activeRequest) return false;
    const current = captureArtifact(activeRequest.groupId);
    return current.records === artifactStamp.records
      && current.group === artifactStamp.group
      && current.selectedGroup === artifactStamp.selectedGroup
      && current.resources === artifactStamp.resources
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
    // The same committed request: its `records` was just checked to be an array.
    (Array.isArray(currentRequest()?.records) ? /** @type {Record<string, any>} */ (currentRequest()).records : [])
      .map((/** @type {{ recordKey?: unknown }} */ { recordKey }) => String(recordKey || ''))
      .filter(Boolean)
  );

  const currentSvg = () => (
    typeof getCurrentSvg === 'function' ? getCurrentSvg() : null
  );

  /**
   * @param {string[]} recordKeys
   * @param {{ records?: unknown } | null} [plan]
   */
  const materializedTranslations = (recordKeys, plan = null) => {
    const base = replacementTranslations(state, recordKeys);
    const svg = currentSvg();
    return svg ? materializeRecordTranslations(svg, base, recordKeys, plan, drawnFeatureIdsOf) : base;
  };

  // `drafted`: the translations a draft's baseline materialized, kept as they are. `fromResult`: read the translations
  // off the shown Result; false keeps the committed ones (a reply without geometry reads none).
  /** @param {{ materializePlan?: boolean, drafted?: unknown, fromResult?: boolean }} [options] */
  const baseline = ({ materializePlan = false, drafted = null, fromResult = materializePlan } = {}) => {
    const canonical = getCommittedSession();
    const request = cloneJson(canonical?.renderRequest);
    if (!request) throw new Error('Generate a Linear diagram before aligning.');
    const recordKeys = request.records.map((/** @type {CommittedRecord} */ { recordKey }) => recordKey);
    const translations = drafted ? cloneJson(drafted)
      : fromResult
        ? materializedTranslations(recordKeys, request.layout?.similarityAlignment)
        : cloneJson(request.layout?.recordTranslations || replacementTranslations(state, recordKeys));
    const orientations = request.records.map((/** @type {CommittedRecord} */ record) => ({ recordKey: record.recordKey,
      reverseComplement: baseReverseComplement(record) }));
    const materialized = projectCommittedAlignment({ committed: canonical,
      plan: materializePlan ? null : request.layout?.similarityAlignment,
      translations, orientations, starts: [] });
    return { canonical: materialized, request: materialized.renderRequest, translations };
  };

  // The draft's baseline on the committed Session as it is now. The stamp keeps its records and
  // translations; a live setting that the automatic rerender committed during the review is
  // in the rest of that Session, so the helper, Apply and the Reset receipt read it (PD-OI-102,
  // OV-1120). Only an open draft calls it, which has its baseline (see `applyPlan`).
  const rebaseDraftBaseline = () => {
    activeBaseline = baseline({ materializePlan: true,
      drafted: /** @type {AlignmentBaseline} */ (activeBaseline).translations });
  };

  /** @param {AlignmentBaseline | null | undefined} value */
  const installBaseState = (value) => {
    if (!value) return;
    state.linearRecordTranslations.value = cloneJson(value.translations);
  };

  /**
   * @param {unknown} message
   * @param {boolean} [whileNoPlan]
   */
  const publishNotice = (message, whileNoPlan = false) => {
    noticeState.value = { text: String(message || ''), whileNoPlan };
  };

  /** @param {string} reason A `msg`-marked English clause; the caller marks it. */
  const announceClearedPlan = (reason) => {
    summary.value = null;
    repair.value = null;
    publishNotice(t('Alignment cleared: {reason}', { reason: tm(reason) }), true);
  };

  /** @param {string} reason A `msg`-marked English clause; the caller marks it. */
  const clearCommittedPlan = (reason) => {
    const hadPlan = Boolean(state.similarityAlignmentPlan.value);
    if (activeApply && typeof cancelRunAnalysis === 'function') cancelRunAnalysis();
    actionId += 1;
    clearDraft();
    status.value = 'idle';
    activeBaseline = null;
    state.similarityAlignmentResetReceipt.value = null;
    resetDialogOpen.value = false;
    if (!hadPlan) return false;
    state.similarityAlignmentPlan.value = null;
    announceClearedPlan(reason);
    return true;
  };

  /**
   * Another record placement committed a step that holds no plan (Align at
   * feature, PD-OI-028 manual record movement): forget the draft and say why.
   * @param {string} reason A `msg`-marked English clause; the caller marks it.
   */
  const notePlanReplaced = (reason) => {
    clearCommittedPlan(reason);
    announceClearedPlan(reason);
  };

  /** @returns {AlignmentChoice[]} */
  const rowChoices = () => (draft.value?.rows || []).filter((/** @type {ReviewRow} */ row) => row.choice
    && (row.choice.kind === 'select' || row.candidates.length > 0)).map((/** @type {ReviewRow & { choice: NonNullable<ReviewRow['choice']> }} */ row) => ({
    recordKey: row.recordKey, kind: row.choice.kind,
    anchor: row.choice.kind === 'select'
      ? row.candidates.find(({ key }) => key === row.choice.candidateKey)?.anchor : null
  }));
  /**
   * @param {AlignmentResolution} [response]
   * @param {AlignmentChoice[]} [choices]
   */
  const projectDraft = (response = draft.value?.response, choices = rowChoices()) => (
    projectSimilarityAlignmentDirections({ resolution: response,
      // The fallback is read only without a draft, where the response has its projection (`installReviewDraft`).
      expectedBinding: draft.value?.binding || response.projection.binding,
      intent: draft.value?.intent || { mode: 'keep' }, choices })
  );
  // The draft as the review opens with it. The records without a choice then are a fact of
  // the draft: they keep their section while the user chooses.
  /**
   * @param {AlignmentResolution} response
   * @param {AlignmentHelperRequest} request
   * @param {ReviewRow[]} [rows]
   * @param {{ repair?: boolean }} [extra]
   */
  const installReviewDraft = (response, request,
    rows = reviewRows(response, request, displayFacts, recordLabels), extra = {}) => {
    draft.value = deepFreeze({ response, binding: response.projection.binding,
      intent: { mode: 'keep' },
      reference: referenceView(response, request, displayFacts, recordLabels),
      rows, openedWithoutChoice: recordsNeedingChoice(rows), ...extra });
    status.value = 'reviewing';
  };

  /**
   * @param {{
   *   canonical: Record<string, any>,
   *   label: string,
   *   alignmentResetBefore?: Record<string, any> | null,
   *   alignmentResetReceipt?: AlignmentResetReceipt | null
   * }} options
   */
  const runAlignmentCandidate = async ({ canonical, label, alignmentResetBefore = null,
    alignmentResetReceipt = undefined }) => {
    // Both callers ran `baseline()` first, which throws when no render request is committed.
    const beforeRecords = new Map(/** @type {Record<string, any>} */ (currentRequest()).records.map((/** @type {CommittedRecord} */ record) => [record.recordKey, record]));
    const changed = canonical.renderRequest.records.filter((/** @type {CommittedRecord} */ record) => (
      baseReverseComplement(record) !== baseReverseComplement(beforeRecords.get(record.recordKey))
    )).map((/** @type {CommittedRecord} */ record) => ({recordKey: record.recordKey, reverseComplement: baseReverseComplement(record)}));
    // R13 port: the root commits `orientations` to the record display intent.
    return runRecordAlignment({ canonical, operation: 'align', label, alignmentResetBefore, alignmentResetReceipt,
      orientations: changed });
  };

  /**
   * @param {AlignmentResolution} response
   * @param {number} expectedActionId
   * @param {ProjectedDirections} projected
   */
  const applyPlan = async (response, expectedActionId, projected) => {
    if (expectedActionId !== actionId) return { status: 'stale' };
    if (!artifactIsCurrent()) return rejectStaleDraft();
    if (!projected.selectionComplete || !projected.geometryValidated) {
      throw new Error('Alignment selection or final geometry is not validated.');
    }
    rebaseDraftBaseline();
    const orientations = projected.records.map(({ recordKey, afterReverseComplement }) => ({
      recordKey, reverseComplement: afterReverseComplement }));
    // An applying draft has its baseline: `start` and `prepareStaleTargetRepair` set it with
    // `activeRequest`, and every clear of one is followed by a stale `actionId`, checked above.
    const canonical = projectCommittedAlignment({
      committed: /** @type {AlignmentBaseline} */ (activeBaseline).canonical,
      plan: response.plan, orientations, starts: [],
      translations: projected.records.map(({ recordKey, translation }) => ({recordKey, ...translation})) });
    const promise = runAlignmentCandidate({canonical, label: msg('Align Similarity Group'),
      alignmentResetBefore: /** @type {AlignmentBaseline} */ (activeBaseline).canonical});
    activeApply = promise;
    let outcome;
    try { outcome = await promise; }
    catch (cause) { outcome = cause?.canceled === true ? {status:'canceled'} : {status:'error', error:cause}; }
    finally { if (activeApply === promise) activeApply = null; }
    if (expectedActionId !== actionId) return {status:'stale'};
    if (outcome?.status === 'ok') {
      // A resolved response always has its plan (`validateSimilarityAlignmentResolution`).
      summary.value = successfulSummary(/** @type {AlignmentPlan} */ (response.plan), projected.records.filter(record => (
        record.beforeReverseComplement !== record.afterReverseComplement)).length);
      repair.value = null; publishNotice(''); activeBaseline = null;
      clearDraft(); error.value = null; status.value = 'idle';
      return {status:'ok'};
    }
    status.value = 'reviewing';
    // A run refused while a Generate runs (`busy`) is no failure: it publishes no error.
    if (['busy', 'canceled', 'stale', 'superseded'].includes(outcome?.status)) return outcome;
    const cause = publishError(outcome?.error);
    return { status: 'error', error: cause, recovery: outcome?.recovery };
  };

  /**
   * @param {AlignmentHelperRequest} request
   * @param {number} expectedActionId
   * @param {string} mode
   */
  const resolveRequest = async (request, expectedActionId, mode) => {
    let response;
    try {
      // `start` sets the baseline right before it resolves; nothing clears it in between.
      const helper = await runHelperOperation(resolveOperation, { request,
        projection: { canonicalRequest: /** @type {AlignmentBaseline} */ (activeBaseline).request,
          orientations: Object.fromEntries(
            request.records.map(record => [record.recordKey, baseReverseComplement(record)])) },
        resources: /** @type {AlignmentBaseline} */ (activeBaseline).canonical.resources });
      if (expectedActionId !== actionId) return { status: 'stale' };
      response = validateSimilarityAlignmentResolution(helper?.result, request);
      if (!artifactIsCurrent()) return rejectStaleDraft();
    } catch (cause) {
      if (expectedActionId !== actionId) return { status: 'stale' };
      clearDraft();
      activeBaseline = null;
      status.value = 'idle';
      if (cause?.canceled === true) return { status: 'canceled' };
      return { status: 'error', error: publishError(cause) };
    }
    activeRequest = request;
    error.value = null;
    // A resolved one-click Align applies without a review: none is shown unless it fails.
    automaticApply.value = mode === 'align' && response.status === 'resolved';
    installReviewDraft(response, request);
    if (automaticApply.value) {
      status.value = 'applying';
      // The draft was just installed from this response, so the binding matches and the projection is never stale.
      try { return await applyPlan(response, expectedActionId, /** @type {ProjectedDirections} */ (projectDraft(response, []))); }
      finally { automaticApply.value = false; }
    }
    return { status: 'reviewing' };
  };

  /** @param {AlignmentStartOptions} options */
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
        choices: []
      });
      setReviewNames(request, members, activeBaseline.request);
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

  /**
   * @param {string} recordKey
   * @param {{ kind: string, anchor?: unknown }} update
   */
  const editRow = (recordKey, update) => {
    if (!activeRequest || !draft.value || status.value !== 'reviewing') {
      return { status: 'rejected' };
    }
    /** @type {ReviewRow | undefined} */
    const row = draft.value.rows.find((/** @type {ReviewRow} */ entry) => entry.recordKey === recordKey);
    if (!row) return { status: 'rejected' };
    let next = row;
    if (update.kind === 'select') {
      const selectedAnchor = validateAnchor(update.anchor, 'Selected alignment candidate');
      const candidate = row.candidates.find(({ anchor }) => sameJson(anchor, selectedAnchor));
      if (!candidate) return { status: 'rejected' };
      next = { ...row, choice: { kind: 'select', candidateKey: candidate.key },
        reason: 'user_selected', reasonLabel: t('Selected by user'),
        unchanged: false, repairRequired: false };
    } else if (update.kind === 'skip') {
      if (!row.candidates.length) return { status: 'rejected' };
      next = { ...row, choice: { kind: 'skip', candidateKey: null },
        reason: 'skipped_by_user', reasonLabel: t('Skipped by user'),
        unchanged: true, repairRequired: false };
    } else return { status: 'rejected' };
    draft.value = deepFreeze({
      ...draft.value,
      rows: draft.value.rows.map((/** @type {ReviewRow} */ entry) => entry === row ? next : entry)
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
    // The draft's own response and binding, so the projection is never stale.
    const preview = /** @type {ProjectedDirections} */ (projectDraft());
    const orientations = Object.fromEntries(preview.records.map(record => [record.recordKey, record.afterReverseComplement]));
    status.value = 'applying';
    try {
      rebaseDraftBaseline();
      // An open draft (`activeRequest`, checked above) always has its baseline (see `applyPlan`).
      const helper = await runHelperOperation(resolveOperation, { request,
        projection: { canonicalRequest: /** @type {AlignmentBaseline} */ (activeBaseline).request, orientations },
        resources: /** @type {AlignmentBaseline} */ (activeBaseline).canonical.resources });
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
      return applyPlan(response, expectedActionId, final);
    } catch (cause) {
      if (expectedActionId !== actionId) return { status: 'stale' };
      status.value = 'reviewing';
      if (cause?.canceled === true) return { status: 'canceled' };
      return { status: 'error', error: publishError(cause) };
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

  // Record and feature names of the committed request, shared by the drawer
  // reference options and the Reset preview; record keys and hashes stay binding-only.
  const committedNames = computed(() => {
    const request = currentRequest();
    const catalog = state.featureCatalog.value;
    return {
      records: recordNames({ recordKeys: (request?.records || []).map((/** @type {CommittedRecord} */ { recordKey }) => recordKey),
        request, linearSeqs: shownLinearCards(), catalog }),
      factsOf: catalogNameIndex(catalog).featureOf
    };
  });

  // The exact reference options of a group's enriched members.
  /** @param {SelectedGroupDetail['members']} members */
  const referenceOptionsOf = (members) => {
    const { records, factsOf } = committedNames.value;
    // `filter(Boolean)` drops the members whose anchor does not resolve (the `null` results).
    const candidates = /** @type {{ anchor: ReturnType<typeof anchorFromSource>, canonicalKey: string, key: string, label: string, location: string }[]} */ (members.map((member) => {
      try {
        const anchor = anchorFromSource(member, 'Drawer reference');
        const record = shortRecordName(records.get(anchor.recordKey)) || t('Unavailable record');
        // The option's name and location read one value: the catalog row of the anchor, else the member's own.
        const facts = factsOf(anchor) || member;
        return {
          anchor,
          canonicalKey: `${anchor.recordKey}\0${anchor.biologicalFeatureId}`,
          key: anchorKey(anchor),
          label: `${record} · ${featureNameParts(anchor, facts).name}`,
          location: sourceLocationText(facts)
        };
      } catch (_error) {
        return null;
      }
    }).filter(Boolean));
    const counts = new Map();
    candidates.forEach(({ canonicalKey }) => {
      counts.set(canonicalKey, (counts.get(canonicalKey) || 0) + 1);
    });
    const unique = candidates.filter(({ canonicalKey }) => counts.get(canonicalKey) === 1);
    // Two features with the same record and name (inparalogs) also show their location.
    const labels = new Map();
    unique.forEach(({ label }) => labels.set(label, (labels.get(label) || 0) + 1));
    return unique.map(({ canonicalKey: _canonicalKey, location, ...candidate }) => (
      labels.get(candidate.label) > 1 ? { ...candidate, label: `${candidate.label} · ${location}` } : candidate
    ));
  };

  // OV-920: the options of the drawer's selected group, from the members the
  // orthogroup owner enriched for its detail; every render reads them.
  /** @type {{ readonly value: { groupId: string, options: ReturnType<typeof referenceOptionsOf> } }} */
  const selectedGroupOptions = computed(() => {
    const detail = selectedOrthogroupDetail.value;
    const groupId = String(detail?.id || '').trim();
    const group = getOrthogroupById(groupId);
    if (!detail || !group) return { groupId, options: [] };
    // The drawer's controls name the group by `id`, as `drawerReferenceOptions` reads it.
    return {
      groupId,
      options: referenceOptionsOf(group === detail.group ? detail.members : getEnrichedOrthogroupMembers(group))
    };
  });

  /** @param {unknown} groupId @returns {ReturnType<typeof referenceOptionsOf>} */
  const drawerReferenceOptions = (groupId) => {
    const id = String(groupId || '').trim();
    const selected = selectedGroupOptions.value;
    if (selected.groupId && selected.groupId === id) return selected.options;
    const group = getOrthogroupById(id);
    if (!group) return [];
    return referenceOptionsOf(getEnrichedOrthogroupMembers(group));
  };

  // OV-381: a new plan (Apply, Session load, Undo) selects its own group in the
  // drawer, whichever group was selected or saved before, so the drawer can
  // offer the plan's reference below.
  watch(() => state.similarityAlignmentPlan.value, (/** @type {{ groupId?: unknown } | null | undefined} */ plan) => {
    const groupId = String(plan?.groupId || '');
    if (groupId && getOrthogroupById(groupId)) state.selectedOrthogroupId.value = groupId;
  });

  // OV-381: until the user picks, the drawer offers the committed plan's exact
  // reference (from a popup Align or a loaded Session), never a stale one.
  const drawerReferenceKey = computed(() => {
    const plan = state.similarityAlignmentPlan.value;
    const pick = drawerReferencePick.value;
    if (pick && pick.plan === plan) return pick.key;
    return plan?.reference && !repair.value ? anchorKey(plan.reference) : '';
  });

  /** @param {ReturnType<typeof referenceOptionsOf>} options */
  const pickedReference = (options) => (
    options.find(({ key }) => key === drawerReferenceKey.value)?.anchor || null
  );
  /** @param {unknown} groupId */
  const selectedDrawerReference = (groupId) => pickedReference(drawerReferenceOptions(groupId));

  /** @param {ReturnType<typeof referenceOptionsOf>} options */
  const disabledReasonOf = (options) => {
    if (['resolving', 'applying'].includes(status.value)) {
      return t('Wait for the current alignment operation to finish.');
    }
    if (!options.length) return t('No exact reference features are available for this group.');
    if (!pickedReference(options)) {
      return t('Select an exact reference record and feature before aligning.');
    }
    return '';
  };
  /** @param {unknown} groupId */
  const drawerDisabledReason = (groupId) => disabledReasonOf(drawerReferenceOptions(groupId));

  // What the drawer shows for its selected group: one value per change of the
  // group, its members, the picked reference, or the operation status.
  const drawer = computed(() => {
    const { groupId, options } = selectedGroupOptions.value;
    return { groupId, options, disabledReason: disabledReasonOf(options) };
  });

  const activePlanInspector = computed(() => {
    try {
      return inspectActivePlan(state.similarityAlignmentPlan.value, { linearSeqs: shownLinearCards(),
        catalog: state.featureCatalog.value, request: currentRequest() });
    } catch (_error) {
      return null;
    }
  });

  const resetPreview = computed(() => resetPreviewView({ receipt: state.similarityAlignmentResetReceipt.value,
    // The committed render request (Python's, R7).
    request: /** @type {Parameters<typeof resetPreviewView>[0]['request']} */ (currentRequest()),
    names: committedNames.value.records }));
  /** @param {string} [scope] */
  const resetAlignment = async (scope = 'positions') => {
    if (!state.similarityAlignmentPlan.value) return {status:'noop'};
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
      // `resetPreview` rejects the directions scope above unless the receipt is set.
      const placement = restoredPlacement({ translations: current.translations, records: current.request.records,
        receipt, scope });
      const { orientations } = placement;
      let { translations } = placement;
      if (placement.directionsChange) {
        /** @type {AlignmentPlan} */
        const plan = current.request.layout.similarityAlignment;
        const members = resetSourceMembers(plan, catalogNameIndex(state.featureCatalog.value).rows);
        const request = buildHelperRequest({group:{id:plan.groupId}, members,
          reference:plan.reference, request:current.request, choices:planChoices(plan)});
        const vector = Object.fromEntries(orientations.map(direction => [direction.recordKey,direction.reverseComplement]));
        const helper = await runHelperOperation(resolveOperation, {request,
          projection:{canonicalRequest:current.request,orientations:vector},resources:current.canonical.resources});
        if (expectedActionId !== actionId) return {status:'stale'};
        translations = restoredAxisTranslations({ translations, vector,
          projection: validateSimilarityAlignmentResolution(helper?.result,request).projection });
      }
      const canonical = projectCommittedAlignment({committed: current.canonical, plan:null, translations, orientations, starts: []});
      promise = runAlignmentCandidate({canonical, label: scope === 'positions'
        ? msg('Reset alignment positions') : msg('Reset alignment positions and directions'), alignmentResetReceipt:null});
      activeApply = promise;
      const outcome = await promise;
      if (expectedActionId !== actionId) return {status:'stale'};
      if (outcome?.status === 'ok') {
        clearDraft(); repair.value = null; summary.value = null; resetDialogOpen.value = false;
        error.value = null;
        publishNotice(scope === 'positions' ? t('Alignment reset: record positions restored; record directions unchanged.')
          : t('Alignment reset: record positions and directions restored to before the latest Align.'));
      } else if (!['busy', 'canceled', 'stale', 'superseded'].includes(outcome?.status)) publishError(outcome?.error);
      return outcome;
    } catch (cause) {
      if (cause?.canceled === true) return { status: 'canceled' };
      if (expectedActionId !== actionId) return { status: 'stale' };
      return { status: 'error', error: publishError(cause) };
    }
    finally { if (activeApply === promise) activeApply = null; if (expectedActionId === actionId) status.value = 'idle'; }
  };

  /**
   * @param {AlignmentPlan} plan
   * @param {Set<string>} [staleRecordKeys]
   * @returns {AlignmentChoice[]}
   */
  const planChoices = (plan, staleRecordKeys = new Set()) => (
    plan.records
      .filter(({ status, recordKey }) => status !== 'reference' && !staleRecordKeys.has(recordKey))
      .map((decision) => ({
        recordKey: decision.recordKey,
        kind: decision.status === 'aligned' ? 'select' : 'skip',
        anchor: decision.status === 'aligned' ? decision.anchor : null
      }))
  );

  /**
   * @param {AlignmentHelperRequest} request
   * @param {OrthogroupMemberFacts[]} members
   * @param {CommittedAlignmentRequest} canonicalRequest
   */
  const setReviewNames = (request, members, canonicalRequest) => {
    displayFacts = new Map(members.flatMap((member) => {
      try { return [[anchorKey(anchorFromSource(member, 'Alignment member')), member]]; } catch (_error) { return []; }
    }));
    recordLabels = recordNames({ recordKeys: request.records.map(({ recordKey }) => recordKey),
      request: canonicalRequest, linearSeqs: shownLinearCards(), catalog: state.featureCatalog.value });
  };

  /** @param {unknown} reason */
  const markStaleReference = (reason) => {
    clearDraft();
    status.value = 'idle';
    repair.value = deepFreeze({
      kind: 'reference',
      groupId: state.similarityAlignmentPlan.value?.groupId || '',
      reason: String(reason || t('The exact reference is no longer available.'))
    });
    publishNotice(t('Alignment needs repair: {reason}', { reason: repair.value.reason }));
    return { status: 'blocked', reason: 'stale-reference' };
  };

  /**
   * @param {AlignmentPlan} plan
   * @param {AlignmentHelperRequest} request
   * @param {Set<string>} staleRecordKeys
   * @param {number} expectedActionId
   */
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
    activeBaseline = repairBaseline(baseline);
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
    installReviewDraft(response, repairRequest, rows, { repair: true });
    repair.value = deepFreeze({
      kind: 'targets',
      groupId: plan.groupId,
      reason: t('One or more saved target anchors are no longer usable. Select or Skip each target.')
    });
    status.value = 'reviewing';
    publishNotice(t('Alignment needs repair: {reason}', { reason: repair.value.reason }));
    return { status: 'blocked', reason: 'stale-target' };
  };

  const validateBeforeGenerate = async () => {
    const expectedActionId = ++actionId;
    const plan = state.similarityAlignmentPlan.value;
    if (!plan) {
      repair.value = null;
      return { status: 'ok' };
    }
    // The plan is Linear draft state: only a Linear request draws it.
    if (state.mode.value !== 'linear') return { status: 'ok' };
    const group = getOrthogroupById(plan.groupId);
    const members = group ? getEnrichedOrthogroupMembers(group) : [];
    const check = classifyPlanAgainstCatalog({ plan, catalogFeatures: catalogNameIndex(state.featureCatalog.value).rows,
      groupAnchors: members.flatMap((member) => {
        try { return [anchorFromSource(member, 'Current Similarity Group member')]; } catch (_error) { return []; }
      }) });
    if (check.kind === 'stale-reference') return markStaleReference(check.reason);
    if (check.kind === 'current') {
      repair.value = null;
      return { status: 'ok' };
    }
    const { goneRecordKeys } = check;
    let request;
    try {
      request = buildHelperRequest({
        group,
        members,
        reference: plan.reference,
        request: currentRequest(),
        choices: planChoices(plan)
      });
    } catch (cause) {
      return markStaleReference(cause?.message || cause);
    }
    setReviewNames(request, members, currentRequest());
    if (goneRecordKeys.size > 0) {
      return prepareStaleTargetRepair(plan, request, goneRecordKeys, expectedActionId);
    }
    try {
      // Python places the plan's anchors on the committed request (D5-13); the check plans and builds no layout, so it
      // reads no placement off the shown Result (one an older gbdraw drew has no record keys to read it by).
      const committed = baseline({ materializePlan: true, fromResult: false });
      const helper = await runHelperOperation(resolveOperation, { request,
        projection: decisionsOnlyProjection(committed.request), resources: committed.canonical.resources });
      if (expectedActionId !== actionId) return { status: 'stale' };
      validateSimilarityAlignmentCheck(helper?.result, request);
    } catch (_cause) {
      if (expectedActionId !== actionId) return { status: 'stale' };
      const alignedKeys = new Set(
        plan.records.filter((/** @type {AlignmentPlanRecord} */ { status }) => status === 'aligned').map((/** @type {AlignmentPlanRecord} */ { recordKey }) => recordKey)
      );
      if (alignedKeys.size > 0) {
        return prepareStaleTargetRepair(plan, request, alignedKeys, expectedActionId);
      }
      return markStaleReference(t('The exact reference is no longer usable with current facts.'));
    }
    repair.value = null;
    return { status: 'ok' };
  };

  /** @param {unknown} recordKeys */
  const retainForStableReorder = (recordKeys) => {
    const orderedKeys = (Array.isArray(recordKeys) ? recordKeys : []).map(String);
    const plan = state.similarityAlignmentPlan.value;
    if (plan) {
      const decisions = new Map(plan.records.map((entry) => [entry.recordKey, entry]));
      if (orderedKeys.every((recordKey) => decisions.has(recordKey))) { // So each `get` below finds its decision.
        state.similarityAlignmentPlan.value = {
          ...cloneJson(plan),
          records: orderedKeys.map((recordKey) => cloneJson(/** @type {AlignmentPlanRecord} */ (decisions.get(recordKey))))
        };
      }
    }
    const translations = new Map(
      state.linearRecordTranslations.value.map((entry) => [entry?.recordKey, entry])
    );
    if (orderedKeys.every((recordKey) => translations.has(recordKey))) {
      state.linearRecordTranslations.value = orderedKeys.map(
        (recordKey) => cloneJson(translations.get(recordKey))
      );
    }
  };

  const beforeRecordDrag = () => {
    pendingRecordDragBaseline = state.similarityAlignmentPlan.value
      ? baseline({ materializePlan: true })
      : null;
    return Boolean(pendingRecordDragBaseline);
  };

  const afterRecordDrag = ({ moved = true } = {}) => {
    if (!pendingRecordDragBaseline || !moved) {
      pendingRecordDragBaseline = null;
      return false;
    }
    const materializedPlan = state.similarityAlignmentPlan.value;
    installBaseState(pendingRecordDragBaseline);
    pendingRecordDragBaseline = null;
    materializedRecordDrag = clearCommittedPlan(msg('record moved manually.'));
    if (!materializedRecordDrag) return false;
    const recordKeys = currentRecordKeys();
    state.linearRecordTranslations.value = materializedTranslations(
      recordKeys,
      materializedPlan
    );
    materializedRecordDrag = false;
    return true;
  };

  const directionPreview = computed(() => {
    if (!draft.value) return null;
    // The draft's own response and binding, so the projection is never stale.
    const projected = /** @type {ProjectedDirections} */ (projectDraft());
    const { reference, rows, intent } = draft.value;
    return describeDirectionPreview({ projected, reference, rows, intent });
  });
  /** @param {DirectionIntent} intent */
  const setDirectionIntent = (intent) => {
    if (!draft.value || status.value !== 'reviewing') return {status:'rejected'};
    projectSimilarityAlignmentDirections({resolution:draft.value.response,
      expectedBinding:draft.value.binding, choices:rowChoices(), intent});
    draft.value = deepFreeze({...draft.value, intent:cloneJson(intent)});
    error.value = null;
    return {status:'selected'};
  };
  /** @param {string} mode */
  const setDirectionMode = (mode) => setDirectionIntent(mode === 'custom'
    ? {mode, byRecordKey:{}} : {mode});
  /**
   * @param {string} recordKey
   * @param {string} direction
   */
  const setCustomDirection = (recordKey, direction) => setDirectionIntent({mode:'custom',
    byRecordKey:{...draft.value?.intent?.byRecordKey, [recordKey]:direction}});
  // The Custom choice of a record, for its select: a record key is data, so only an own key
  // names a choice (OV-572); a record without one keeps its direction, as the projection reads it.
  /** @param {string} recordKey */
  const customDirection = (recordKey) => ownEntry(draft.value?.intent?.byRecordKey, recordKey) || 'keep';

  const unresolvedCount = computed(() => recordsNeedingChoice(draft.value?.rows || []).length);

  const applyDisabledReason = computed(() => {
    if (status.value === 'applying') return t('Applying the alignment plan.');
    if (status.value === 'resolving') return t('Resolving the current alignment.');
    if (unresolvedCount.value > 0) {
      return tn(unresolvedCount.value, 'Choose a candidate or Skip for {count} record.', 'Choose a candidate or Skip for {count} records.');
    }
    if (status.value !== 'reviewing') return t('Open an alignment review before applying.');
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
    notePlanReplaced,
    /**
     * The committed Session with the active plan materialized into the record
     * translations the Result draws, and no plan: the base a record placement starts from.
     * @returns {CommittedCanonicalSession}
     */
    materializedSession: () => /** @type {CommittedCanonicalSession} */ (baseline({ materializePlan: true }).canonical),
    repairBaseSession: repairBaseSessionReader(baseline),
    retainForStableReorder,
    beforeRecordDrag,
    afterRecordDrag,
    directionPreview,
    setDirectionMode,
    setCustomDirection,
    customDirection,
    resetPreview,
    resetDialogOpen,
    resetScope,
    openReset: () => {
      if (!state.similarityAlignmentPlan.value || busy.value) return;
      error.value = null; resetScope.value = 'positions'; resetDialogOpen.value = true;
    },
    cancelReset: () => { if (!busy.value) resetDialogOpen.value = false; },
    applyReset: () => resetAlignment(resetScope.value),
    /** @param {AlignmentStartOptions} options */
    startFromPopup: (options) => start({ ...options, source: 'popup' }),
    /** @param {{ groupId?: unknown, mode?: string } | undefined} options */
    startFromDrawer: (options) => start({
      ...options,
      reference: selectedDrawerReference(options?.groupId),
      source: 'drawer'
    }),
    drawerReferenceOptions,
    drawer,
    /**
     * @param {unknown} groupId
     * @param {string} key
     */
    setDrawerReference: (groupId, key) => {
      const option = drawerReferenceOptions(groupId).find((entry) => entry.key === key);
      drawerReferencePick.value = { key: option?.key || '', plan: state.similarityAlignmentPlan.value };
      return Boolean(option);
    },
    drawerDisabledReason,
    /**
     * @param {string} recordKey
     * @param {unknown} anchor `editRow` validates it
     */
    selectCandidate: (recordKey, anchor) => editRow(recordKey, { kind: 'select', anchor }),
    /** @param {string} recordKey */
    skipRecord: (recordKey) => editRow(recordKey, { kind: 'skip' }),
    applyDraft,
    cancel,
    /** @param {AliasSource} anchor */
    previewCandidate: (anchor) => {
      if (typeof previewCandidate === 'function') previewCandidate(anchor);
    },
    clearCandidatePreview: () => {
      if (typeof clearCandidatePreview === 'function') clearCandidatePreview();
    }
  };
};
