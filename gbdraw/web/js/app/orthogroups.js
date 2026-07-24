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
import { downloadTextFile } from '../services/text-download.js';

const { computed } = window.Vue;

const normalizeText = (value) => String(value ?? '').trim();
const normalizeLower = (value) => normalizeText(value).toLowerCase();

const normalizeRecordIndex = (value) => {
  if (value === null || value === undefined || value === '') return null;
  const normalized = Number(value);
  return Number.isInteger(normalized) && normalized >= 0 ? normalized : null;
};

const memberStableFeatureId = (member) => normalizeText(
  member?.stableFeatureSvgId ||
  member?.stable_feature_svg_id ||
  member?.stableFeatureId ||
  member?.stable_feature_id ||
  member?.featureSvgId ||
  member?.feature_svg_id
);

const featureStableFeatureId = (feature) => normalizeText(
  feature?.stable_feature_id ||
  feature?.stableFeatureId ||
  feature?.stable_svg_id ||
  feature?.stableFeatureSvgId ||
  feature?.svg_id ||
  feature?.svgId
);

const featureRecordIndex = (feature) => {
  for (const value of [feature?.fileIdx, feature?.recordIndex, feature?.record_index, feature?.record_idx]) {
    const recordIndex = normalizeRecordIndex(value);
    if (recordIndex !== null) return recordIndex;
  }
  return null;
};

const memberRecordIndex = (member) => normalizeRecordIndex(member?.recordIndex ?? member?.record_index);
const stableRecordKey = (recordIndex, stableId) => (
  recordIndex !== null && stableId ? `${recordIndex}:${stableId}` : ''
);

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

const copyTextToClipboard = async (text) => {
  const value = String(text ?? '');
  if (navigator.clipboard?.writeText) {
    await navigator.clipboard.writeText(value);
    return;
  }
  const textarea = document.createElement('textarea');
  textarea.value = value;
  textarea.setAttribute('readonly', '');
  textarea.style.position = 'fixed';
  textarea.style.left = '-9999px';
  textarea.style.top = '0';
  document.body.appendChild(textarea);
  textarea.select();
  try {
    document.execCommand('copy');
  } finally {
    document.body.removeChild(textarea);
  }
};

const numericOrthogroupId = (id) => {
  const match = String(id || '').match(/^og_(\d+)$/i);
  return match ? Number(match[1]) : Number.POSITIVE_INFINITY;
};

const compareOrthogroupId = (left, right) => {
  const leftId = normalizeText(left?.id);
  const rightId = normalizeText(right?.id);
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
  return {
    ...sourceFeature,
    record_id: sourceFeature.record_id || sourceFeature.recordId || member?.recordId || member?.record_id,
    start: sourceFeature.start ?? member?.start,
    end: sourceFeature.end ?? member?.end,
    strand: sourceFeature.strand || normalizeMemberStrand(member?.strand),
    source_protein_id: sourceFeature.source_protein_id || sourceFeature.sourceProteinId || member?.sourceProteinId || member?.source_protein_id,
    protein_id: sourceFeature.protein_id || sourceFeature.proteinId || member?.proteinId || member?.protein_id,
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
  member?.proteinId,
  member?.sourceProteinId,
  member?.gene,
  member?.product,
  member?.note,
  member?.recordId,
  member?.label,
  memberStableFeatureId(member)
].map(normalizeLower).join(' ');

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

const addRenderedId = (index, key, renderedId) => {
  if (!key || !renderedId) return;
  const values = index.get(key) || new Set();
  values.add(renderedId);
  index.set(key, values);
};

const buildRenderedFeatureIndex = (svg, features = []) => {
  const byStableRecord = new Map();
  const byStable = new Map();
  Array.from(svg?.querySelectorAll?.(FEATURE_SELECTOR) || []).forEach((element) => {
    const renderedId = normalizeText(getFeatureIdentity(element));
    if (!renderedId) return;
    const stableId = normalizeText(
      element.getAttribute?.('data-gbdraw-stable-feature-id') || renderedId
    );
    const recordIndex = normalizeRecordIndex(
      element.getAttribute?.('data-gbdraw-record-index')
    );
    addRenderedId(byStableRecord, stableRecordKey(recordIndex, stableId), renderedId);
    addRenderedId(byStable, stableId, renderedId);
  });
  (Array.isArray(features) ? features : []).forEach((feature) => {
    const renderedId = normalizeText(
      feature?.rendered_svg_id || feature?.renderedSvgId || feature?.svg_id
    );
    const stableId = featureStableFeatureId(feature);
    const recordIndex = featureRecordIndex(feature);
    if (!renderedId || !stableId || getFeatureElements(svg, renderedId).length === 0) return;
    addRenderedId(byStableRecord, stableRecordKey(recordIndex, stableId), renderedId);
    addRenderedId(byStable, stableId, renderedId);
  });
  return { byStableRecord, byStable };
};

const renderedFeatureIdForMember = (member, renderedIndex) => {
  const stableId = memberStableFeatureId(member);
  if (!stableId) return '';
  const recordIndex = memberRecordIndex(member);
  const recordMatches = recordIndex === null
    ? null
    : renderedIndex.byStableRecord.get(stableRecordKey(recordIndex, stableId));
  if (recordIndex !== null) {
    return recordMatches?.size === 1 ? recordMatches.values().next().value : '';
  }
  const stableMatches = renderedIndex.byStable.get(stableId);
  return stableMatches?.size === 1 ? stableMatches.values().next().value : '';
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
    showRightDrawer,
    rightDrawerTab,
    showFeaturePanel,
    showLegendPanel,
    linearSeqs,
    extractedFeatures,
    biologicalFeatures
  } = state;

  const getOrthogroupById = (orthogroupId) => {
    const id = normalizeText(orthogroupId);
    if (!id) return null;
    return (Array.isArray(orthogroups.value) ? orthogroups.value : [])
      .find((group) => normalizeText(group?.id) === id) || null;
  };

  const resolveOrthogroupName = (groupOrId) => {
    const group = typeof groupOrId === 'string' ? getOrthogroupById(groupOrId) : groupOrId;
    const id = normalizeText(group?.id || groupOrId);
    if (!id) return '';
    return normalizeText(orthogroupNameOverrides[id]) || normalizeText(group?.name) || id;
  };

  const resolveOrthogroupDescription = (groupOrId) => {
    const group = typeof groupOrId === 'string' ? getOrthogroupById(groupOrId) : groupOrId;
    const id = normalizeText(group?.id || groupOrId);
    if (!id) return '';
    return normalizeText(orthogroupDescriptionOverrides[id]) || normalizeText(group?.description);
  };

  const orthogroupScope = (groupOrId) => {
    const group = typeof groupOrId === 'string' ? getOrthogroupById(groupOrId) : groupOrId;
    return normalizeGroupMetadataScope(group?.scope);
  };

  const orthogroupScopeLabel = (groupOrId) => groupMetadataScopeLabel(orthogroupScope(groupOrId));

  const isOrthogroupRenamed = (groupOrId) => {
    const id = normalizeText(typeof groupOrId === 'string' ? groupOrId : groupOrId?.id);
    return Boolean(
      id &&
      (
        Object.prototype.hasOwnProperty.call(orthogroupNameOverrides, id) ||
        Object.prototype.hasOwnProperty.call(orthogroupDescriptionOverrides, id)
      )
    );
  };

  const orthogroupCount = computed(() => (Array.isArray(orthogroups.value) ? orthogroups.value.length : 0));

  const selectedAlignmentTargetLabel = computed(() => {
    const target = normalizeText(selectedOrthogroupAlignmentFeature.value);
    if (!target) return '';
    const group = getOrthogroupById(target);
    return group ? `${resolveOrthogroupName(group)} (${target})` : target;
  });

  const filteredOrthogroups = computed(() => {
    const query = normalizeLower(orthogroupSearch.value);
    const groups = Array.isArray(orthogroups.value) ? [...orthogroups.value] : [];
    const filtered = query
      ? groups.filter((group) => {
          const candidates = Array.isArray(group?.nameCandidates) ? group.nameCandidates : [];
          const members = getGroupMembers(group);
          const haystack = [
            group?.id,
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
    return filteredOrthogroups.value[0] || (Array.isArray(orthogroups.value) ? orthogroups.value[0] : null) || null;
  });

  const featureSequenceLookup = computed(() => {
    const lookup = new Map();
    const features = Array.isArray(biologicalFeatures?.value) && biologicalFeatures.value.length > 0
      ? biologicalFeatures.value
      : (Array.isArray(extractedFeatures?.value) ? extractedFeatures.value : []);
    const stableCounts = new Map();
    features.forEach((feature) => {
      const stableId = featureStableFeatureId(feature);
      if (stableId) stableCounts.set(stableId, (stableCounts.get(stableId) || 0) + 1);
    });
    features.forEach((feature) => {
      const stableId = featureStableFeatureId(feature);
      if (!stableId) return;
      const recordIndex = featureRecordIndex(feature);
      const entry = {
        nucleotideSequence: firstSequenceText(feature?.nucleotideSequence, feature?.nucleotide_sequence),
        aminoAcidSequence: firstSequenceText(feature?.aminoAcidSequence, feature?.amino_acid_sequence),
        sequenceFeature: feature,
        sequenceWarnings: Array.isArray(feature?.sequence_warnings)
          ? feature.sequence_warnings
          : (Array.isArray(feature?.sequenceWarnings) ? feature.sequenceWarnings : [])
      };
      if (!entry.nucleotideSequence && !entry.aminoAcidSequence) return;
      const recordKey = stableRecordKey(recordIndex, stableId);
      if (recordKey) lookup.set(recordKey, entry);
      if ((stableCounts.get(stableId) || 0) === 1) lookup.set(stableId, entry);
    });
    return lookup;
  });

  const enrichOrthogroupMember = (member) => {
    const stableId = memberStableFeatureId(member);
    if (!stableId) return member;
    const recordIndex = memberRecordIndex(member);
    const lookup = featureSequenceLookup.value instanceof Map ? featureSequenceLookup.value : new Map();
    const sequenceEntry = (
      recordIndex !== null ? lookup.get(stableRecordKey(recordIndex, stableId)) : null
    ) || lookup.get(stableId) || null;
    if (!sequenceEntry) return member;
    return {
      ...member,
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
      const recordIndex = Number(member?.recordIndex);
      const key = Number.isInteger(recordIndex) ? recordIndex : -1;
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
    const orthogroupId = normalizeText(group?.id || groupOrId);
    return getEnrichedOrthogroupMembers(group)
      .map((member) => memberFastaText(member, sequenceKind, orthogroupId))
      .filter(Boolean)
      .join('');
  };

  const buildOrthogroupMemberFasta = (member, sequenceKind, groupOrId = selectedOrthogroup.value) => {
    const group = typeof groupOrId === 'string' ? getOrthogroupById(groupOrId) : groupOrId;
    return memberFastaText(enrichOrthogroupMember(member), sequenceKind, normalizeText(group?.id || groupOrId));
  };

  const orthogroupSequenceFilename = (groupOrId, sequenceKind) => {
    const group = typeof groupOrId === 'string' ? getOrthogroupById(groupOrId) : groupOrId;
    const id = normalizeText(group?.id || groupOrId);
    const name = makeSafeFilename(resolveOrthogroupName(group) || id, id || 'orthogroup');
    const stem = makeSafeFilename(`${id || 'orthogroup'}_${name}_${sequenceKindLabel(sequenceKind)}`);
    return `${stem}.${sequenceExtension(sequenceKind)}`;
  };

  const orthogroupMemberSequenceFilename = (member, sequenceKind, groupOrId = selectedOrthogroup.value) => {
    const group = typeof groupOrId === 'string' ? getOrthogroupById(groupOrId) : groupOrId;
    const id = normalizeText(group?.id || groupOrId) || 'orthogroup';
    const memberId = normalizeText(member?.sourceProteinId || member?.proteinId || memberStableFeatureId(member) || 'member');
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
    if (!id || !getOrthogroupById(id)) return;
    selectedOrthogroupId.value = id;
  };

  const setOrthogroupNameOverride = (orthogroupId, value) => {
    const id = normalizeText(orthogroupId);
    if (!id) return;
    const group = getOrthogroupById(id);
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
    if (!id) return;
    selectedOrthogroupAlignmentFeature.value = id;
    if (typeof runAnalysis === 'function') await runAnalysis();
  };

  const resetOrthogroupAlignment = async () => {
    if (!selectedOrthogroupAlignmentFeature.value) return;
    selectedOrthogroupAlignmentFeature.value = '';
    if (typeof runAnalysis === 'function') await runAnalysis();
  };

  const openRightDrawerTab = (tab) => {
    const normalized = ['legend', 'features', 'orthogroups'].includes(tab) ? tab : 'features';
    if (normalized === 'orthogroups' && orthogroupCount.value === 0) return;
    rightDrawerTab.value = normalized;
    showRightDrawer.value = true;
    showFeaturePanel.value = normalized === 'features';
    showLegendPanel.value = normalized === 'legend';
  };

  const closeRightDrawer = () => {
    showRightDrawer.value = false;
    showFeaturePanel.value = false;
    showLegendPanel.value = false;
  };

  const openOrthogroupInDrawer = (orthogroupId) => {
    const id = normalizeText(orthogroupId);
    if (id && getOrthogroupById(id)) selectedOrthogroupId.value = id;
    openRightDrawerTab('orthogroups');
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
    resetOrthogroupAlignment,
    openRightDrawerTab,
    closeRightDrawer,
    openOrthogroupInDrawer
  };
};
