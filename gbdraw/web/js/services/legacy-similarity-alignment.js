import { featureIdentity } from './feature-identity.js';

const requiredText = (value, path) => {
  if (typeof value !== 'string' || !value.trim() || value.includes('\0')) {
    throw new Error(`${path} must be non-empty text without NUL.`);
  }
  return value.trim();
};

const catalogItems = (featureCatalog) => {
  if (!featureCatalog || featureCatalog.schema !== 3 || !Array.isArray(featureCatalog.items)) {
    return [];
  }
  return featureCatalog.items;
};

const legacyGroups = (orthogroupState) => (
  Array.isArray(orthogroupState?.groups) ? orthogroupState.groups : []
);

/** Materialize one released string-valued alignment selection at the legacy boundary. */
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
  const biologicalFeatures = items.flatMap((item) => (
    Array.isArray(item?.biologicalFeatures) ? item.biologicalFeatures : []
  ));
  const renderedFeatures = items.flatMap((item) => (
    Array.isArray(item?.features) ? item.features : []
  ));
  const catalogGroups = items.flatMap((item) => (
    Array.isArray(item?.orthogroups) ? item.orthogroups : []
  ));
  const useLegacyGroups = catalogGroups.length === 0;
  const groups = useLegacyGroups
    ? legacyGroups(legacyOrthogroupState)
    : catalogGroups;
  const memberRecordKey = (member) => {
    if (!useLegacyGroups) return member?.recordKey;
    if (!Number.isSafeInteger(member?.recordIndex) ||
        member.recordIndex < 0 || member.recordIndex >= recordKeys.length) {
      throw new Error('Legacy similarity alignment member has no stable record mapping.');
    }
    return recordKeys[member.recordIndex];
  };
  const groupsById = groups.filter((group) => String(group?.id || '') === legacyTarget);
  let group = null;
  let selectedMember = null;
  if (groupsById.length === 1) {
    group = groupsById[0];
    const members = Array.isArray(group.members) ? group.members : [];
    selectedMember = members.find((member) => member?.representative === true) || members[0] || null;
  } else if (groupsById.length > 1) {
    throw new Error('Legacy similarity alignment group is ambiguous in the saved feature catalog.');
  } else {
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
  const anchorFor = (member) => {
    if (useLegacyGroups) {
      const biologicalFeatureId = requiredText(
        member?.stableFeatureSvgId ||
          member?.stable_feature_svg_id ||
          member?.featureSvgId,
        'legacy member stable feature identity'
      );
      if (!Number.isSafeInteger(member?.featureIndex) || member.featureIndex < 0) {
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
    const [biological] = biologicalMatches;
    return {
      recordKey: requiredText(member.recordKey, 'legacy member recordKey'),
      biologicalFeatureId: requiredText(
        member.biologicalFeatureId,
        'legacy member biologicalFeatureId'
      ),
      sourceFeatureIndex: Number.isSafeInteger(biological.sourceFeatureIndex) &&
        biological.sourceFeatureIndex >= 0
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
        anchor: reference,
        orientationPolicy: 'preserve',
        effectiveReverseComplement: null
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
      anchor: null,
      orientationPolicy: 'preserve',
      effectiveReverseComplement: null
    };
    const record = records.find((item) => item.recordKey === recordKey);
    return {
      recordKey,
      status: 'aligned',
      rationale: 'only_usable_candidate',
      anchor: anchorFor(member),
      orientationPolicy: 'preserve',
      effectiveReverseComplement: Boolean(record?.presentation?.reverseComplement)
        !== Boolean(record?.region?.reverseComplement)
    };
  });
  return {
    schema: 2,
    groupId: requiredText(group.id, 'legacy orthogroup id'),
    reference,
    records: decisions
  };
};

/** Supply canonical identities only where a released legacy member has one stable mapping. */
export const migrateLegacyOrthogroupMembers = (groups, records) => {
  const recordKeys = (Array.isArray(records) ? records : []).map((record) => record?.recordKey);
  return (Array.isArray(groups) ? groups : []).map((group) => ({
    ...group,
    members: (Array.isArray(group?.members) ? group.members : []).map((member) => {
      const identity = featureIdentity(member);
      if (!identity.valid || identity.recordKey.supplied || identity.biologicalId.supplied ||
          !identity.recordIndex.supplied || !identity.sourceIndex.supplied ||
          !identity.stableId.supplied || typeof recordKeys[identity.recordIndex.value] !== 'string' ||
          !recordKeys[identity.recordIndex.value]) return member;
      const { featureIndex: _legacyIndex, feature_index: _legacySnakeIndex,
        sourceFeatureIndex: _legacySourceIndex,
        source_feature_index: _legacySnakeSourceIndex, ...stableMember } = member;
      return {
        ...stableMember,
        recordKey: recordKeys[identity.recordIndex.value],
        biologicalFeatureId: identity.stableId.value
      };
    })
  }));
};
