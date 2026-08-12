import {
  FEATURE_SELECTOR,
  getFeatureElements,
  getFeatureIdentity
} from './feature-editor/svg-actions.js';
import { buildFeatureSequenceFastas } from './feature-sequence-fasta.js';
import {
  groupMetadataScopeLabel,
  normalizeGroupMetadataScope
} from './losat-normalization.js';
import {
  isInternalProteinDisplayId,
  resolveDisplayProteinId
} from './feature-utils.js';
import { downloadTextFile } from '../services/text-download.js';
import { copyTextToClipboard } from '../utils/clipboard.js';
import {
  RECORD_INDEX_KEYS,
  STABLE_FEATURE_ID_KEYS,
  featureIdentity,
  identityMatches,
  nonnegativeIntegerAliasStatus,
  orthogroupIdStatus,
  renderedFeatureIdentity,
  textAliasStatus,
  uniqueOrthogroupEntries
} from '../services/feature-identity.js';

export { resolveUniqueOrthogroupMemberForFeature } from '../services/feature-identity.js';

const { computed } = window.Vue;

const normalizeText = (value) => String(value ?? '').trim();
const normalizeLower = (value) => normalizeText(value).toLowerCase();

const memberStableFeatureId = (member) => {
  const status = textAliasStatus(member, STABLE_FEATURE_ID_KEYS);
  return status.valid ? status.value : '';
};

const normalizeSequence = (value) => String(value ?? '').replace(/\s+/g, '').toUpperCase();

const firstSequenceText = (...values) => {
  for (const value of values) {
    const normalized = normalizeSequence(value);
    if (normalized) return normalized;
  }
  return '';
};

const makeSafeFilename = (value, fallback = 'orthogroup') => {
  const cleaned = normalizeText(value).replace(/[^\w.-]+/g, '_').replace(/^_+|_+$/g, '');
  return cleaned || fallback;
};

const sequenceKindLabel = (sequenceKind) => (sequenceKind === 'aa' ? 'aa' : 'nt');
const sequenceExtension = (sequenceKind) => (sequenceKindLabel(sequenceKind) === 'aa' ? 'faa' : 'fna');

const numericOrthogroupId = (id) => {
  const match = String(id || '').match(/^og_(\d+)$/i);
  return match ? Number(match[1]) : Number.POSITIVE_INFINITY;
};

const orthogroupIdValue = (group) => {
  const status = orthogroupIdStatus(group);
  return status.valid && status.supplied ? status.value : '';
};

const compareOrthogroupId = (left, right) => {
  const leftId = orthogroupIdValue(left);
  const rightId = orthogroupIdValue(right);
  const leftNumber = numericOrthogroupId(leftId);
  const rightNumber = numericOrthogroupId(rightId);
  if (leftNumber !== rightNumber) return leftNumber - rightNumber;
  return leftId.localeCompare(rightId);
};

const getGroupMembers = (group) => (Array.isArray(group?.members) ? group.members : []);

const getMemberSequence = (member, sequenceKind) => (
  sequenceKindLabel(sequenceKind) === 'aa'
    ? firstSequenceText(member?.aminoAcidSequence, member?.amino_acid_sequence, member?.proteinSequence, member?.sequence)
    : firstSequenceText(member?.nucleotideSequence, member?.nucleotide_sequence)
);

const memberLocationText = (member) => {
  const start = Number(member?.start);
  const end = Number(member?.end);
  if (!Number.isFinite(start) || !Number.isFinite(end)) return '';
  const normalizedStrand = normalizeMemberStrand(member?.strand);
  const strand = normalizedStrand ? `(${normalizedStrand})` : '';
  return `${start + 1}-${end}${strand}`;
};

const normalizeMemberStrand = (strand) => {
  if (strand === -1 || String(strand).trim() === '-1') return '-';
  if (strand === 1 || String(strand).trim() === '1') return '+';
  return normalizeText(strand);
};

const buildMemberFeaturePayload = (member) => {
  const sourceFeature = member?.sequenceFeature && typeof member.sequenceFeature === 'object'
    ? member.sequenceFeature
    : {};
  const displayProteinId = resolveDisplayProteinId(
    sourceFeature,
    member,
    memberLocationText(member)
  );
  return {
    ...sourceFeature,
    record_id: sourceFeature.record_id || sourceFeature.recordId || member?.recordId || member?.record_id,
    start: sourceFeature.start ?? member?.start,
    end: sourceFeature.end ?? member?.end,
    strand: sourceFeature.strand || normalizeMemberStrand(member?.strand),
    source_protein_id: displayProteinId,
    protein_id: displayProteinId,
    product: sourceFeature.product || member?.product,
    note: sourceFeature.note || member?.note,
    gene: sourceFeature.gene || member?.gene,
    organism: sourceFeature.organism || member?.organism,
    nucleotide_sequence: getMemberSequence(member, 'nt'),
    amino_acid_sequence: getMemberSequence(member, 'aa')
  };
};

const memberFastaText = (member, sequenceKind, orthogroupId = '') => {
  const feature = buildMemberFeaturePayload(member);
  const fastas = buildFeatureSequenceFastas(feature, {
    nucleotideSequence: getMemberSequence(member, 'nt'),
    aminoAcidSequence: getMemberSequence(member, 'aa')
  });
  const text = sequenceKindLabel(sequenceKind) === 'aa' ? fastas.aminoAcidFasta : fastas.nucleotideFasta;
  return text ? `${text}\n` : '';
};

const getMemberSearchText = (member) => [
  member?.displayProteinId,
  member?.proteinId,
  member?.sourceProteinId,
  member?.gene,
  member?.product,
  member?.note,
  member?.recordId,
  member?.label,
  memberStableFeatureId(member)
]
  .filter((value) => !isInternalProteinDisplayId(value))
  .map(normalizeLower)
  .join(' ');

const getCandidateSearchText = (candidate) => [
  candidate?.text,
  candidate?.source
].map(normalizeLower).join(' ');

const setOriginalStroke = (el) => {
  if (!el.hasAttribute('data-og-original-stroke')) {
    el.setAttribute('data-og-original-stroke', el.getAttribute('stroke') || '');
    el.setAttribute('data-og-original-stroke-width', el.getAttribute('stroke-width') || '');
  }
};

const domIdentityAgrees = (element, candidate) => {
  const stableId = textAliasStatus({
    stableId: element?.getAttribute?.('data-gbdraw-stable-feature-id')
  }, ['stableId']);
  const recordIndex = nonnegativeIntegerAliasStatus({
    recordIndex: element?.getAttribute?.('data-gbdraw-record-index')
  }, ['recordIndex']);
  if (!stableId.valid || !recordIndex.valid) return false;
  if (
    stableId.supplied &&
    (!candidate.stableId.supplied || stableId.value !== candidate.stableId.value)
  ) return false;
  if (
    recordIndex.supplied &&
    (!candidate.recordIndex.supplied || recordIndex.value !== candidate.recordIndex.value)
  ) return false;
  return true;
};

const buildRenderedFeatureIndex = (svg, features = []) => {
  const metadataByRenderedId = new Map();
  (Array.isArray(features) ? features : []).forEach((feature) => {
    const candidate = renderedFeatureIdentity(feature);
    if (!candidate.usable || getFeatureElements(svg, candidate.renderedId.value).length === 0) return;
    const values = metadataByRenderedId.get(candidate.renderedId.value) || [];
    values.push(candidate);
    metadataByRenderedId.set(candidate.renderedId.value, values);
  });

  const candidates = [];
  metadataByRenderedId.forEach((values, renderedId) => {
    if (values.length !== 1) return;
    const elements = getFeatureElements(svg, renderedId);
    if (elements.every((element) => domIdentityAgrees(element, values[0]))) {
      candidates.push(values[0]);
    }
  });

  const domByRenderedId = new Map();
  Array.from(svg?.querySelectorAll?.(FEATURE_SELECTOR) || []).forEach((element) => {
    const renderedId = normalizeText(getFeatureIdentity(element));
    if (!renderedId || metadataByRenderedId.has(renderedId)) return;
    const source = {
      renderedFeatureSvgId: renderedId,
      stableFeatureSvgId: element.getAttribute?.('data-gbdraw-stable-feature-id'),
      recordIndex: element.getAttribute?.('data-gbdraw-record-index')
    };
    const candidate = featureIdentity(source);
    if (!candidate.usable || !candidate.renderedId.supplied) return;
    const existing = domByRenderedId.get(renderedId);
    if (!existing) {
      domByRenderedId.set(renderedId, candidate);
      return;
    }
    if (
      !identityMatches(existing, candidate, { includeRendered: true }) ||
      !identityMatches(candidate, existing, { includeRendered: true })
    ) {
      domByRenderedId.set(renderedId, null);
    }
  });
  domByRenderedId.forEach((candidate) => {
    if (candidate) candidates.push(candidate);
  });
  return candidates;
};

const renderedFeatureIdForMember = (member, renderedIndex) => {
  const identity = featureIdentity(member);
  if (!identity.usable) return '';
  const matches = (Array.isArray(renderedIndex) ? renderedIndex : [])
    .filter((candidate) => identityMatches(identity, candidate, { includeRendered: true }));
  return matches.length === 1 ? matches[0].renderedId.value : '';
};

export const createOrthogroupEditor = ({ state, runAnalysis }) => {
  const {
    orthogroups,
    orthogroupNameOverrides,
    orthogroupDescriptionOverrides,
    selectedOrthogroupId,
    orthogroupSearch,
    orthogroupSortMode,
    selectedOrthogroupAlignmentFeature,
    svgContainer,
    linearSeqs,
    extractedFeatures,
    biologicalFeatures
  } = state;

  const resolvableOrthogroupEntries = () => uniqueOrthogroupEntries(orthogroups.value);

  const getOrthogroupById = (orthogroupId) => {
    const id = normalizeText(orthogroupId);
    if (!id) return null;
    const matches = resolvableOrthogroupEntries()
      .filter((entry) => entry.id === id);
    return matches.length === 1 ? matches[0].group : null;
  };

  const resolveOrthogroupName = (groupOrId) => {
    const group = typeof groupOrId === 'string' ? getOrthogroupById(groupOrId) : groupOrId;
    const id = group ? orthogroupIdValue(group) : normalizeText(groupOrId);
    if (!id) return '';
    return normalizeText(orthogroupNameOverrides[id]) || normalizeText(group?.name) || id;
  };

  const resolveOrthogroupDescription = (groupOrId) => {
    const group = typeof groupOrId === 'string' ? getOrthogroupById(groupOrId) : groupOrId;
    const id = group ? orthogroupIdValue(group) : normalizeText(groupOrId);
    if (!id) return '';
    return normalizeText(orthogroupDescriptionOverrides[id]) || normalizeText(group?.description);
  };

  const orthogroupScope = (groupOrId) => {
    const group = typeof groupOrId === 'string' ? getOrthogroupById(groupOrId) : groupOrId;
    return normalizeGroupMetadataScope(group?.scope);
  };

  const orthogroupScopeLabel = (groupOrId) => groupMetadataScopeLabel(orthogroupScope(groupOrId));

  const isOrthogroupRenamed = (groupOrId) => {
    const id = typeof groupOrId === 'string'
      ? normalizeText(groupOrId)
      : orthogroupIdValue(groupOrId);
    return Boolean(
      id &&
      (
        Object.prototype.hasOwnProperty.call(orthogroupNameOverrides, id) ||
        Object.prototype.hasOwnProperty.call(orthogroupDescriptionOverrides, id)
      )
    );
  };

  const orthogroupCount = computed(() => resolvableOrthogroupEntries().length);

  const selectedAlignmentTargetLabel = computed(() => {
    const target = normalizeText(selectedOrthogroupAlignmentFeature.value);
    if (!target) return '';
    const group = getOrthogroupById(target);
    return group ? `${resolveOrthogroupName(group)} (${target})` : target;
  });

  const filteredOrthogroups = computed(() => {
    const query = normalizeLower(orthogroupSearch.value);
    const groups = resolvableOrthogroupEntries().map(({ group }) => group);
    const filtered = query
      ? groups.filter((group) => {
          const candidates = Array.isArray(group?.nameCandidates) ? group.nameCandidates : [];
          const members = getGroupMembers(group);
          const haystack = [
            orthogroupIdValue(group),
            group?.name,
            group?.description,
            resolveOrthogroupName(group),
            resolveOrthogroupDescription(group),
            ...candidates.map(getCandidateSearchText),
            ...members.map(getMemberSearchText)
          ].map(normalizeLower).join(' ');
          return haystack.includes(query);
        })
      : groups;

    const modeRaw = normalizeText(orthogroupSortMode.value) || 'id';
    const mode = ['id', 'name', 'member_count', 'record_coverage'].includes(modeRaw) ? modeRaw : 'id';
    filtered.sort((left, right) => {
      if (mode === 'name') {
        return resolveOrthogroupName(left).localeCompare(resolveOrthogroupName(right)) || compareOrthogroupId(left, right);
      }
      if (mode === 'member_count') {
        const leftCount = Number(left?.member_count || getGroupMembers(left).length || 0);
        const rightCount = Number(right?.member_count || getGroupMembers(right).length || 0);
        return rightCount - leftCount || compareOrthogroupId(left, right);
      }
      if (mode === 'record_coverage') {
        const leftCoverage = Number(left?.record_coverage_count || 0);
        const rightCoverage = Number(right?.record_coverage_count || 0);
        return rightCoverage - leftCoverage || compareOrthogroupId(left, right);
      }
      return compareOrthogroupId(left, right);
    });
    return filtered;
  });

  const selectedOrthogroup = computed(() => {
    const selected = getOrthogroupById(selectedOrthogroupId.value);
    if (selected) return selected;
    return filteredOrthogroups.value[0] || null;
  });

  const featureSequenceLookup = computed(() => {
    const hasBiologicalFeatures = Array.isArray(biologicalFeatures?.value) && biologicalFeatures.value.length > 0;
    const features = hasBiologicalFeatures
      ? biologicalFeatures.value
      : (Array.isArray(extractedFeatures?.value) ? extractedFeatures.value : []);
    const renderedIdentities = (Array.isArray(extractedFeatures?.value) ? extractedFeatures.value : [])
      .map(renderedFeatureIdentity)
      .filter((identity) => identity.usable);
    return features.map((feature) => {
      const identity = hasBiologicalFeatures
        ? featureIdentity(feature, { allowLegacySvgStable: true })
        : renderedFeatureIdentity(feature);
      const renderedMatches = identity.usable
        ? renderedIdentities.filter((rendered) => (
          identityMatches(identity, rendered, { includeRendered: true })
        ))
        : [];
      const renderedIdentityInvalid = renderedMatches.length > 1 || (
        identity.renderedId.supplied && renderedMatches.length !== 1
      );
      const resolvedIdentity = renderedIdentityInvalid
        ? { ...identity, usable: false }
        : (
            renderedMatches.length === 1 && !identity.renderedId.supplied
              ? { ...identity, renderedId: renderedMatches[0].renderedId }
              : identity
          );
      return {
        identity: resolvedIdentity,
        entry: {
          nucleotideSequence: firstSequenceText(feature?.nucleotideSequence, feature?.nucleotide_sequence),
          aminoAcidSequence: firstSequenceText(feature?.aminoAcidSequence, feature?.amino_acid_sequence),
          sequenceFeature: feature,
          sequenceWarnings: Array.isArray(feature?.sequence_warnings)
            ? feature.sequence_warnings
            : (Array.isArray(feature?.sequenceWarnings) ? feature.sequenceWarnings : [])
        }
      };
    });
  });

  const enrichOrthogroupMember = (member) => {
    const memberSequenceFeature = member?.sequenceFeature || null;
    const memberIdentity = featureIdentity(member);
    const matches = memberIdentity.usable
      ? (Array.isArray(featureSequenceLookup.value) ? featureSequenceLookup.value : [])
        .filter((candidate) => (
          identityMatches(memberIdentity, candidate.identity, { includeRendered: true })
        ))
      : [];
    const sequenceEntry = matches.length === 1 ? matches[0].entry : null;
    const sequenceFeature = sequenceEntry?.sequenceFeature || memberSequenceFeature;
    const resolvedDisplayProteinId = resolveDisplayProteinId(
      sequenceFeature,
      member,
      memberLocationText(member)
    );
    if (!sequenceEntry) return { ...member, displayProteinId: resolvedDisplayProteinId };
    return {
      ...member,
      displayProteinId: resolvedDisplayProteinId,
      nucleotideSequence: firstSequenceText(member?.nucleotideSequence, member?.nucleotide_sequence, sequenceEntry.nucleotideSequence),
      aminoAcidSequence: firstSequenceText(member?.aminoAcidSequence, member?.amino_acid_sequence, sequenceEntry.aminoAcidSequence),
      sequenceFeature: member?.sequenceFeature || sequenceEntry.sequenceFeature,
      sequenceWarnings: Array.isArray(member?.sequenceWarnings)
        ? member.sequenceWarnings
        : (Array.isArray(member?.sequence_warnings) ? member.sequence_warnings : sequenceEntry.sequenceWarnings)
    };
  };

  const getEnrichedOrthogroupMembers = (groupOrId) => {
    const group = typeof groupOrId === 'string' ? getOrthogroupById(groupOrId) : groupOrId;
    return getGroupMembers(group).map(enrichOrthogroupMember);
  };

  const groupOrthogroupMembersByRecord = (members) => {
    const byRecord = new Map();
    (Array.isArray(members) ? members : []).forEach((member) => {
      const recordIndex = nonnegativeIntegerAliasStatus(member, RECORD_INDEX_KEYS);
      const key = recordIndex.valid && recordIndex.supplied ? recordIndex.value : -1;
      if (!byRecord.has(key)) byRecord.set(key, []);
      byRecord.get(key).push(member);
    });
    return Array.from(byRecord.entries())
      .sort((left, right) => left[0] - right[0])
      .map(([recordIndex, members]) => ({
        recordIndex,
        recordLabel: recordIndex >= 0
          ? (
              linearSeqs[recordIndex]?.definition ||
              linearSeqs[recordIndex]?.gb?.name ||
              linearSeqs[recordIndex]?.gff?.name ||
              `Record ${recordIndex + 1}`
            )
          : 'Record',
        members
      }));
  };

  const selectedOrthogroupMembersByRecord = computed(() => {
    const group = selectedOrthogroup.value;
    return groupOrthogroupMembersByRecord(getEnrichedOrthogroupMembers(group));
  });

  const getOrthogroupSequenceCount = (groupOrId, sequenceKind) =>
    getEnrichedOrthogroupMembers(groupOrId).filter((member) => getMemberSequence(member, sequenceKind)).length;

  const hasOrthogroupSequence = (groupOrId, sequenceKind) => getOrthogroupSequenceCount(groupOrId, sequenceKind) > 0;

  const hasOrthogroupMemberSequence = (member, sequenceKind) =>
    Boolean(getMemberSequence(enrichOrthogroupMember(member), sequenceKind));

  const buildOrthogroupFasta = (groupOrId, sequenceKind) => {
    const group = typeof groupOrId === 'string' ? getOrthogroupById(groupOrId) : groupOrId;
    const orthogroupId = group ? orthogroupIdValue(group) : normalizeText(groupOrId);
    return getEnrichedOrthogroupMembers(group)
      .map((member) => memberFastaText(member, sequenceKind, orthogroupId))
      .filter(Boolean)
      .join('');
  };

  const buildOrthogroupMemberFasta = (member, sequenceKind, groupOrId = selectedOrthogroup.value) => {
    const group = typeof groupOrId === 'string' ? getOrthogroupById(groupOrId) : groupOrId;
    const orthogroupId = group ? orthogroupIdValue(group) : normalizeText(groupOrId);
    return memberFastaText(enrichOrthogroupMember(member), sequenceKind, orthogroupId);
  };

  const orthogroupSequenceFilename = (groupOrId, sequenceKind) => {
    const group = typeof groupOrId === 'string' ? getOrthogroupById(groupOrId) : groupOrId;
    const id = group ? orthogroupIdValue(group) : normalizeText(groupOrId);
    const name = makeSafeFilename(resolveOrthogroupName(group) || id, id || 'orthogroup');
    const stem = makeSafeFilename(`${id || 'orthogroup'}_${name}_${sequenceKindLabel(sequenceKind)}`);
    return `${stem}.${sequenceExtension(sequenceKind)}`;
  };

  const orthogroupMemberSequenceFilename = (member, sequenceKind, groupOrId = selectedOrthogroup.value) => {
    const group = typeof groupOrId === 'string' ? getOrthogroupById(groupOrId) : groupOrId;
    const id = (group ? orthogroupIdValue(group) : normalizeText(groupOrId)) || 'orthogroup';
    const stableId = memberStableFeatureId(member);
    const memberId = normalizeText(
      resolveDisplayProteinId(
        member?.sequenceFeature,
        member,
        memberLocationText(member)
      ) ||
      (isInternalProteinDisplayId(stableId) ? '' : stableId) ||
      'member'
    );
    const stem = makeSafeFilename(`${id}_${memberId}_${sequenceKindLabel(sequenceKind)}`);
    return `${stem}.${sequenceExtension(sequenceKind)}`;
  };

  const copyOrthogroupSequences = async (groupOrId = selectedOrthogroup.value, sequenceKind = 'nt') => {
    const text = buildOrthogroupFasta(groupOrId, sequenceKind);
    if (!text) return;
    await copyTextToClipboard(text);
  };

  const downloadOrthogroupSequences = (groupOrId = selectedOrthogroup.value, sequenceKind = 'nt') => {
    const text = buildOrthogroupFasta(groupOrId, sequenceKind);
    if (!text) return;
    downloadTextFile(
      orthogroupSequenceFilename(groupOrId, sequenceKind),
      text,
      'text/plain;charset=utf-8'
    );
  };

  const copyOrthogroupMemberSequence = async (member, sequenceKind = 'nt', groupOrId = selectedOrthogroup.value) => {
    const text = buildOrthogroupMemberFasta(member, sequenceKind, groupOrId);
    if (!text) return;
    await copyTextToClipboard(text);
  };

  const downloadOrthogroupMemberSequence = (member, sequenceKind = 'nt', groupOrId = selectedOrthogroup.value) => {
    const text = buildOrthogroupMemberFasta(member, sequenceKind, groupOrId);
    if (!text) return;
    downloadTextFile(
      orthogroupMemberSequenceFilename(member, sequenceKind, groupOrId),
      text,
      'text/plain;charset=utf-8'
    );
  };

  const selectOrthogroup = (orthogroupId) => {
    const id = normalizeText(orthogroupId);
    if (!id || !getOrthogroupById(id)) return false;
    selectedOrthogroupId.value = id;
    return true;
  };

  const setOrthogroupNameOverride = (orthogroupId, value) => {
    const id = normalizeText(orthogroupId);
    if (!id) return;
    const group = getOrthogroupById(id);
    if (!group) return;
    const normalized = normalizeText(value);
    const base = normalizeText(group?.name);
    if (!normalized || normalized === base) {
      delete orthogroupNameOverrides[id];
      return;
    }
    orthogroupNameOverrides[id] = normalized;
  };

  const setOrthogroupDescriptionOverride = (orthogroupId, value) => {
    const id = normalizeText(orthogroupId);
    if (!id) return;
    const group = getOrthogroupById(id);
    if (!group) return;
    const normalized = normalizeText(value);
    const base = normalizeText(group?.description);
    if (!normalized || normalized === base) {
      delete orthogroupDescriptionOverrides[id];
      return;
    }
    orthogroupDescriptionOverrides[id] = normalized;
  };

  const resetOrthogroupRename = (orthogroupId = selectedOrthogroupId.value) => {
    const id = normalizeText(orthogroupId);
    if (!id) return;
    delete orthogroupNameOverrides[id];
    delete orthogroupDescriptionOverrides[id];
  };

  const clearOrthogroupHighlight = () => {
    if (!svgContainer.value) return;
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return;
    svg.querySelectorAll('[data-og-original-stroke], [data-og-original-stroke-width]').forEach((el) => {
      const stroke = el.getAttribute('data-og-original-stroke');
      const strokeWidth = el.getAttribute('data-og-original-stroke-width');
      if (stroke) el.setAttribute('stroke', stroke);
      else el.removeAttribute('stroke');
      if (strokeWidth) el.setAttribute('stroke-width', strokeWidth);
      else el.removeAttribute('stroke-width');
      el.removeAttribute('data-og-original-stroke');
      el.removeAttribute('data-og-original-stroke-width');
    });
  };

  const highlightOrthogroupById = (orthogroupId = selectedOrthogroupId.value) => {
    const group = getOrthogroupById(orthogroupId);
    if (!group || !svgContainer.value) return;
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return;
    clearOrthogroupHighlight();
    const renderedIndex = buildRenderedFeatureIndex(svg, extractedFeatures?.value);
    const featureIds = new Set(
      getGroupMembers(group)
        .map((member) => renderedFeatureIdForMember(member, renderedIndex))
        .filter(Boolean)
    );
    featureIds.forEach((featureId) => {
      getFeatureElements(svg, featureId).forEach((el) => {
        setOriginalStroke(el);
        el.setAttribute('stroke', '#2563eb');
        el.setAttribute('stroke-width', '2.4');
      });
    });
  };

  const alignOrthogroupById = async (orthogroupId = selectedOrthogroupId.value) => {
    const id = normalizeText(orthogroupId);
    if (!id || !getOrthogroupById(id)) return;
    selectedOrthogroupAlignmentFeature.value = id;
    if (typeof runAnalysis === 'function') await runAnalysis();
  };

  const resetOrthogroupAlignment = async () => {
    if (!selectedOrthogroupAlignmentFeature.value) return;
    selectedOrthogroupAlignmentFeature.value = '';
    if (typeof runAnalysis === 'function') await runAnalysis();
  };

  return {
    orthogroupCount,
    selectedAlignmentTargetLabel,
    filteredOrthogroups,
    selectedOrthogroup,
    selectedOrthogroupMembersByRecord,
    getEnrichedOrthogroupMembers,
    groupOrthogroupMembersByRecord,
    getOrthogroupById,
    resolveOrthogroupName,
    resolveOrthogroupDescription,
    orthogroupScope,
    orthogroupScopeLabel,
    isOrthogroupRenamed,
    getOrthogroupSequenceCount,
    hasOrthogroupSequence,
    hasOrthogroupMemberSequence,
    copyOrthogroupSequences,
    downloadOrthogroupSequences,
    copyOrthogroupMemberSequence,
    downloadOrthogroupMemberSequence,
    selectOrthogroup,
    setOrthogroupNameOverride,
    setOrthogroupDescriptionOverride,
    resetOrthogroupRename,
    highlightOrthogroupById,
    clearOrthogroupHighlight,
    alignOrthogroupById,
    resetOrthogroupAlignment
  };
};
