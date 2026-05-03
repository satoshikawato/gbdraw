const { computed } = window.Vue;

const normalizeText = (value) => String(value ?? '').trim();
const normalizeLower = (value) => normalizeText(value).toLowerCase();

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

const getMemberSearchText = (member) => [
  member?.proteinId,
  member?.sourceProteinId,
  member?.gene,
  member?.product,
  member?.note,
  member?.recordId,
  member?.label,
  member?.featureSvgId
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
    linearSeqs
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

  const selectedOrthogroupMembersByRecord = computed(() => {
    const group = selectedOrthogroup.value;
    const byRecord = new Map();
    getGroupMembers(group).forEach((member) => {
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
  });

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
    const featureIds = new Set(
      getGroupMembers(group)
        .map((member) => normalizeText(member?.featureSvgId))
        .filter(Boolean)
    );
    featureIds.forEach((featureId) => {
      svg.querySelectorAll(`#${CSS.escape(featureId)}`).forEach((el) => {
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
    getOrthogroupById,
    resolveOrthogroupName,
    resolveOrthogroupDescription,
    isOrthogroupRenamed,
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
