// @ts-check
import { featureIdentity } from './feature-identity.js';
import { LEGACY_FEATURE_CATALOG_SCHEMA, SIMILARITY_ALIGNMENT_PLAN_SCHEMA } from '../contracts/schema.generated.js';

// The released Sessions this module reads were not validated: every field is
// `unknown`, and an entry of a list may be missing.
/**
 * A member of a saved orthogroup (a catalog member or a legacy state member).
 * @typedef {object} LegacyMember
 * @property {unknown} [recordKey]
 * @property {unknown} [recordIndex]
 * @property {unknown} [representative]
 * @property {unknown} [biologicalFeatureId]
 * @property {unknown} [stableFeatureSvgId]
 * @property {unknown} [stable_feature_svg_id]
 * @property {unknown} [featureSvgId]
 * @property {unknown} [featureIndex]
 * @property {unknown} [sourceProteinId]
 * @property {unknown} [proteinId]
 * @property {unknown} [label]
 */
/** @typedef {{ id?: unknown, members?: unknown }} LegacyGroup */
/**
 * @typedef {object} LegacyBiologicalFeature
 * @property {unknown} [recordKey]
 * @property {unknown} [biologicalFeatureId]
 * @property {unknown} [protein_id]
 * @property {unknown} [source_protein_id]
 * @property {unknown} [sourceFeatureIndex]
 * @property {unknown} [stableFeatureId]
 */
/** @typedef {{ recordKey?: unknown, biologicalFeatureId?: unknown, svgId?: unknown }} LegacyRenderedFeature */
/** @typedef {{ biologicalFeatures?: unknown, features?: unknown, orthogroups?: unknown }} LegacyCatalogItem */
/** @typedef {{ schema?: unknown, items?: unknown }} LegacyFeatureCatalog */
/** @typedef {{ recordKey?: unknown }} LegacyRecord */

/**
 * @param {unknown} value
 * @param {string} path
 * @returns {string}
 */
const requiredText = (value, path) => {
  if (typeof value !== 'string' || !value.trim() || value.includes('\0')) {
    throw new Error(`${path} must be non-empty text without NUL.`);
  }
  return value.trim();
};

/**
 * @param {LegacyFeatureCatalog | null | undefined} featureCatalog
 * @returns {(LegacyCatalogItem | null | undefined)[]}
 */
const catalogItems = (featureCatalog) => {
  if (!featureCatalog || featureCatalog.schema !== LEGACY_FEATURE_CATALOG_SCHEMA || !Array.isArray(featureCatalog.items)) {
    return [];
  }
  return featureCatalog.items;
};

/**
 * @param {{ groups?: unknown } | null | undefined} orthogroupState
 * @returns {(LegacyGroup | null | undefined)[]}
 */
const legacyGroups = (orthogroupState) => (
  Array.isArray(orthogroupState?.groups) ? orthogroupState.groups : []
);

/**
 * Materialize one released string-valued alignment selection at the legacy boundary.
 * @param {{
 *   target: unknown,
 *   records: readonly LegacyRecord[],
 *   featureCatalog: LegacyFeatureCatalog | null | undefined,
 *   legacyOrthogroupState?: { groups?: unknown } | null
 * }} input
 */
export const materializeLegacySimilarityAlignment = ({
  target,
  records,
  featureCatalog,
  legacyOrthogroupState = null
}) => {
  const legacyTarget = requiredText(target, 'legacy similarity alignment target');
  const recordKeys = records.map((record, index) => requiredText(
    record.recordKey,
    `renderRequest.records[${index}].recordKey`
  ));
  const items = catalogItems(featureCatalog);
  /** @type {(LegacyBiologicalFeature | null | undefined)[]} */
  const biologicalFeatures = items.flatMap((item) => (
    Array.isArray(item?.biologicalFeatures) ? item.biologicalFeatures : []
  ));
  /** @type {(LegacyRenderedFeature | null | undefined)[]} */
  const renderedFeatures = items.flatMap((item) => (
    Array.isArray(item?.features) ? item.features : []
  ));
  /** @type {(LegacyGroup | null | undefined)[]} */
  const catalogGroups = items.flatMap((item) => (
    Array.isArray(item?.orthogroups) ? item.orthogroups : []
  ));
  const useLegacyGroups = catalogGroups.length === 0;
  const groups = useLegacyGroups
    ? legacyGroups(legacyOrthogroupState)
    : catalogGroups;
  /** @param {LegacyMember | null | undefined} member */
  const memberRecordKey = (member) => {
    if (!useLegacyGroups) return member?.recordKey;
    // `isSafeInteger` holds in the later operands: `member` exists and `recordIndex` is a number.
    if (!Number.isSafeInteger(member?.recordIndex) ||
        /** @type {{ recordIndex: number }} */ (member).recordIndex < 0 ||
        /** @type {{ recordIndex: number }} */ (member).recordIndex >= recordKeys.length) {
      throw new Error('Legacy similarity alignment member has no stable record mapping.');
    }
    return recordKeys[/** @type {{ recordIndex: number }} */ (member).recordIndex];
  };
  const groupsById = groups.filter((group) => String(group?.id || '') === legacyTarget);
  /** @type {LegacyGroup | null | undefined} */
  let group = null;
  /** @type {LegacyMember | null | undefined} */
  let selectedMember = null;
  if (groupsById.length === 1) {
    // The length check above makes this the one matching group.
    group = /** @type {LegacyGroup} */ (groupsById[0]);
    /** @type {(LegacyMember | null | undefined)[]} */
    const members = Array.isArray(group.members) ? group.members : [];
    selectedMember = members.find((member) => member?.representative === true) || members[0] || null;
  } else if (groupsById.length > 1) {
    throw new Error('Legacy similarity alignment group is ambiguous in the saved feature catalog.');
  } else {
    /** @type {{ group: LegacyGroup | null | undefined, member: LegacyMember | null | undefined }[]} */
    const matches = [];
    groups.forEach((candidateGroup) => {
      (Array.isArray(candidateGroup?.members) ? candidateGroup.members : []).forEach((member) => {
        const recordKey = memberRecordKey(member);
        const biological = biologicalFeatures.find((feature) => (
          feature?.recordKey === recordKey &&
          feature?.biologicalFeatureId === member?.biologicalFeatureId
        ));
        const rendered = renderedFeatures.find((feature) => (
          feature?.recordKey === recordKey &&
          feature?.biologicalFeatureId === member?.biologicalFeatureId
        ));
        const identities = new Set([
          member?.biologicalFeatureId,
          member?.stableFeatureSvgId,
          member?.stable_feature_svg_id,
          member?.featureSvgId,
          member?.sourceProteinId,
          member?.proteinId,
          member?.label,
          biological?.biologicalFeatureId,
          biological?.protein_id,
          biological?.source_protein_id,
          rendered?.svgId
        ].map((value) => String(value || '').trim()).filter(Boolean));
        if (identities.has(legacyTarget)) matches.push({ group: candidateGroup, member });
      });
    });
    if (matches.length !== 1) {
      throw new Error(
        matches.length === 0
          ? 'Legacy similarity alignment target is not represented by saved stable metadata.'
          : 'Legacy similarity alignment target matches multiple saved group members.'
      );
    }
    ({ group, member: selectedMember } = matches[0]);
  }
  if (!group || !selectedMember) {
    throw new Error('Legacy similarity alignment has no materializable reference member.');
  }
  /** @param {LegacyMember} member */
  const anchorFor = (member) => {
    if (useLegacyGroups) {
      const biologicalFeatureId = requiredText(
        member?.stableFeatureSvgId ||
          member?.stable_feature_svg_id ||
          member?.featureSvgId,
        'legacy member stable feature identity'
      );
      if (!Number.isSafeInteger(member?.featureIndex) || /** @type {number} */ (member.featureIndex) < 0) { // number after `isSafeInteger`
        throw new Error('Legacy similarity alignment member lacks a source feature index.');
      }
      return {
        recordKey: requiredText(memberRecordKey(member), 'legacy member recordKey'),
        biologicalFeatureId,
        // Legacy featureIndex was an orthogroup subset ordinal, not a source feature index.
        sourceFeatureIndex: null,
        stableFeatureSvgId: biologicalFeatureId
      };
    }
    const biologicalMatches = biologicalFeatures.filter((feature) => (
      feature?.recordKey === member?.recordKey &&
      feature?.biologicalFeatureId === member?.biologicalFeatureId
    ));
    if (biologicalMatches.length !== 1) {
      throw new Error(
        'Legacy similarity alignment member lacks unique saved biological identity metadata.'
      );
    }
    // The unique match has a string `recordKey` equal to the member's (`requiredText` below), so it is not a missing entry.
    const [biological] = /** @type {[LegacyBiologicalFeature]} */ (biologicalMatches);
    return {
      recordKey: requiredText(member.recordKey, 'legacy member recordKey'),
      biologicalFeatureId: requiredText(
        member.biologicalFeatureId,
        'legacy member biologicalFeatureId'
      ),
      // A number after `isSafeInteger`.
      sourceFeatureIndex: Number.isSafeInteger(biological.sourceFeatureIndex) &&
        /** @type {number} */ (biological.sourceFeatureIndex) >= 0
        ? biological.sourceFeatureIndex
        : null,
      stableFeatureSvgId: requiredText(
        biological.stableFeatureId || member.biologicalFeatureId,
        'legacy member stable feature identity'
      )
    };
  };
  const reference = anchorFor(selectedMember);
  const members = Array.isArray(group.members) ? group.members : [];
  const decisions = recordKeys.map((recordKey) => {
    if (recordKey === reference.recordKey) {
      return {
        recordKey,
        status: 'reference',
        rationale: 'reference',
        anchor: reference
      };
    }
    const candidates = members.filter((member) => memberRecordKey(member) === recordKey);
    const representatives = candidates.filter(
      (candidate) => candidate?.representative === true
    );
    if (representatives.length > 1 || (candidates.length > 1 && representatives.length !== 1)) {
      throw new Error(
        'Legacy similarity alignment has ambiguous members for one displayed record.'
      );
    }
    const member = representatives[0] || candidates[0] || null;
    if (!member) return {
      recordKey,
      status: 'skipped',
      rationale: 'skipped_no_candidate',
      anchor: null
    };
    return {
      recordKey,
      status: 'aligned',
      rationale: 'only_usable_candidate',
      anchor: anchorFor(member)
    };
  });
  return {
    schema: SIMILARITY_ALIGNMENT_PLAN_SCHEMA,
    groupId: requiredText(group.id, 'legacy orthogroup id'),
    reference,
    records: decisions
  };
};

/**
 * Supply canonical identities only where a released legacy member has one stable mapping.
 * @param {readonly (LegacyGroup | null | undefined)[] | null | undefined} groups
 * @param {readonly (LegacyRecord | null | undefined)[] | null | undefined} records
 */
export const migrateLegacyOrthogroupMembers = (groups, records) => {
  const recordKeys = (Array.isArray(records) ? records : []).map((record) => record?.recordKey);
  return (Array.isArray(groups) ? groups : []).map((group) => ({
    ...group,
    members: (/** @type {(LegacyMember | null | undefined)[]} */ (
      Array.isArray(group?.members) ? group.members : []
    )).map((member) => {
      const identity = featureIdentity(member);
      // `recordIndex.supplied` is checked below, and a supplied valid index has a number value.
      if (!identity.valid || identity.recordKey.supplied || identity.biologicalId.supplied ||
          !identity.recordIndex.supplied || !identity.sourceIndex.supplied ||
          !identity.stableId.supplied || typeof recordKeys[/** @type {number} */ (identity.recordIndex.value)] !== 'string' ||
          !recordKeys[/** @type {number} */ (identity.recordIndex.value)]) return member;
      // Every operand above held, so `featureIdentity` read supplied fields of an object member.
      const { featureIndex: _legacyIndex, feature_index: _legacySnakeIndex,
        sourceFeatureIndex: _legacySourceIndex,
        source_feature_index: _legacySnakeSourceIndex, ...stableMember } = /** @type {Record<string, unknown>} */ (member);
      return {
        ...stableMember,
        recordKey: recordKeys[/** @type {number} */ (identity.recordIndex.value)],
        biologicalFeatureId: identity.stableId.value
      };
    })
  }));
};
