// @ts-check
import {
  FEATURE_SELECTOR,
  getFeatureElements,
  getFeatureIdentity
} from './feature-editor/svg-actions.js';
import { buildFeatureSequenceFastas } from '../services/feature-sequence-fasta.js';
import {
  groupMetadataScopeLabel,
  normalizeGroupMetadataScope
} from '../services/losat-normalization.js';
import {
  isInternalProteinDisplayId,
  locationPartBound,
  locationRangeText,
  resolveDisplayProteinId
} from '../services/feature-utils.js';
import { downloadTextFile } from '../services/text-download.js';
import { linearRecordNames } from '../services/linear-sources.js';
import { copyTextToClipboard } from '../utils/clipboard.js';
import { ownEntry } from '../utils/own-entry.js';
import { GROUP_NAMES_SLICE, clearGroupRename, setGroupDescription, setGroupName } from '../services/drawing/group-names.js';
import {
  ORTHOGROUP_ID_KEYS,
  RECORD_INDEX_KEYS,
  STABLE_FEATURE_ID_KEYS,
  featureIdentity,
  identityMatches,
  indexFeatureIdentities,
  nonnegativeIntegerAliasStatus,
  orthogroupIdStatus,
  renderedFeatureIdentity,
  textAliasStatus,
  uniqueOrthogroupEntries
} from '../services/feature-identity.js';

/** @import { FeatureIdentity } from '../services/feature-identity.js' */
/** @import { FeatureSvgRoot } from '../services/feature-dom.js' */

/**
 * One member of a Similarity group as the analysis rows hold it; Python owns the
 * fields (R7), so they are read as `unknown`.
 * @typedef {Record<string, unknown>} OrthogroupMember
 */
/**
 * One Similarity group row. Its id is any of the alias keys `orthogroupIdStatus` reads.
 * @typedef {Record<string, unknown> & {
 *   members?: OrthogroupMember[],
 *   nameCandidates?: { text?: unknown, source?: unknown }[]
 * }} OrthogroupGroup
 */
/** @typedef {OrthogroupGroup | string | null | undefined} GroupOrId A group row, or the id of one. */

export { resolveUniqueOrthogroupMemberForFeature } from '../services/feature-identity.js';

const { computed, nextTick, reactive, watch, toDisplayString } = window.Vue;

const COPY_FEEDBACK_DURATION_MS = 1500;
// Keyed by `success` and `failure`, which are the only statuses this module
// stores; a missing status reads undefined.
/** @type {Readonly<Record<string, string | undefined>>} */
const COPY_FEEDBACK_LABELS = Object.freeze({
  success: 'Copied!',
  failure: 'Copy failed'
});

/** @param {unknown} value */
const normalizeText = (value) => String(value ?? '').trim();
/** @param {unknown} value */
const normalizeLower = (value) => normalizeText(value).toLowerCase();

/** @param {OrthogroupMember | null | undefined} member */
const memberStableFeatureId = (member) => {
  const status = textAliasStatus(member, STABLE_FEATURE_ID_KEYS);
  return status.valid ? status.value : '';
};

/** @param {OrthogroupMember | null | undefined} member */
const memberCopyFeedbackIdentity = (member) => {
  const identity = featureIdentity(member);
  if (identity.usable) {
    return JSON.stringify([
      identity.recordIndex.supplied ? identity.recordIndex.value : null,
      identity.recordKey.supplied ? identity.recordKey.value : '',
      identity.biologicalId.supplied ? identity.biologicalId.value : '',
      identity.sourceIndex.supplied ? identity.sourceIndex.value : null,
      identity.stableId.supplied ? identity.stableId.value : ''
    ]);
  }
  return JSON.stringify([
    normalizeText(member?.recordId ?? member?.record_id),
    memberStableFeatureId(member),
    normalizeText(member?.start),
    normalizeText(member?.end),
    normalizeText(member?.strand),
    normalizeText(resolveDisplayProteinId(/** @type {Record<string, unknown> | null | undefined} */ (member?.sequenceFeature), member, ''))
  ]);
};

/** @param {unknown} value */
const normalizeSequence = (value) => String(value ?? '').replace(/\s+/g, '').toUpperCase();

/** @param {unknown[]} values */
const firstSequenceText = (...values) => {
  for (const value of values) {
    const normalized = normalizeSequence(value);
    if (normalized) return normalized;
  }
  return '';
};

/** @param {unknown} value @param {string} [fallback] */
const makeSafeFilename = (value, fallback = 'orthogroup') => {
  const cleaned = normalizeText(value).replace(/[^\w.-]+/g, '_').replace(/^_+|_+$/g, '');
  return cleaned || fallback;
};

/** @param {unknown} sequenceKind */
const sequenceKindLabel = (sequenceKind) => (sequenceKind === 'aa' ? 'aa' : 'nt');
/** @param {unknown} sequenceKind */
const sequenceExtension = (sequenceKind) => (sequenceKindLabel(sequenceKind) === 'aa' ? 'faa' : 'fna');

/** @param {unknown} id */
const numericOrthogroupId = (id) => {
  const match = String(id || '').match(/^og_(\d+)$/i);
  return match ? Number(match[1]) : Number.POSITIVE_INFINITY;
};

/** @param {OrthogroupGroup | null | undefined} group */
const orthogroupIdValue = (group) => {
  const status = orthogroupIdStatus(group);
  return status.valid && status.supplied ? status.value : '';
};

/** @param {OrthogroupGroup} left @param {OrthogroupGroup} right */
const compareOrthogroupId = (left, right) => {
  const leftId = orthogroupIdValue(left);
  const rightId = orthogroupIdValue(right);
  const leftNumber = numericOrthogroupId(leftId);
  const rightNumber = numericOrthogroupId(rightId);
  if (leftNumber !== rightNumber) return leftNumber - rightNumber;
  return leftId.localeCompare(rightId);
};

/** @param {OrthogroupGroup | null | undefined} group @returns {OrthogroupMember[]} */
const getGroupMembers = (group) => (Array.isArray(group?.members) ? group.members : []);

/** @param {OrthogroupMember | null | undefined} member @param {unknown} sequenceKind */
const getMemberSequence = (member, sequenceKind) => (
  sequenceKindLabel(sequenceKind) === 'aa'
    ? firstSequenceText(member?.aminoAcidSequence, member?.amino_acid_sequence, member?.proteinSequence, member?.sequence)
    : firstSequenceText(member?.nucleotideSequence, member?.nucleotide_sequence)
);

/** @param {OrthogroupMember | null | undefined} member */
const memberLocationText = (member) => {
  const start = locationPartBound(member, 'start');
  const end = locationPartBound(member, 'end');
  if (!Number.isFinite(start) || !Number.isFinite(end)) return '';
  const normalizedStrand = normalizeMemberStrand(member?.strand);
  const strand = normalizedStrand ? `(${normalizedStrand})` : '';
  return `${start + 1}-${end}${strand}`;
};

/** @param {unknown} strand */
const normalizeMemberStrand = (strand) => {
  if (strand === -1 || String(strand).trim() === '-1') return '-';
  if (strand === 1 || String(strand).trim() === '1') return '+';
  return normalizeText(strand);
};

/** @param {OrthogroupMember | null | undefined} member */
const buildMemberFeaturePayload = (member) => {
  // A non-null object here; its fields are read as `unknown`.
  const sourceFeature = /** @type {Record<string, unknown>} */ (member?.sequenceFeature && typeof member.sequenceFeature === 'object'
    ? member.sequenceFeature
    : {});
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

/** @param {OrthogroupMember | null | undefined} member @param {unknown} sequenceKind */
const memberFastaText = (member, sequenceKind) => {
  const feature = buildMemberFeaturePayload(member);
  const fastas = buildFeatureSequenceFastas(feature, {
    nucleotideSequence: getMemberSequence(member, 'nt'),
    aminoAcidSequence: getMemberSequence(member, 'aa')
  });
  const text = sequenceKindLabel(sequenceKind) === 'aa' ? fastas.aminoAcidFasta : fastas.nucleotideFasta;
  return text ? `${text}\n` : '';
};

/** @param {OrthogroupMember | null | undefined} member */
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

/** @param {{ text?: unknown, source?: unknown } | null | undefined} candidate */
const getCandidateSearchText = (candidate) => [
  candidate?.text,
  candidate?.source
].map(normalizeLower).join(' ');

/** @param {Element} el */
const setOriginalStroke = (el) => {
  if (!el.hasAttribute('data-og-original-stroke')) {
    el.setAttribute('data-og-original-stroke', el.getAttribute('stroke') || '');
    el.setAttribute('data-og-original-stroke-width', el.getAttribute('stroke-width') || '');
  }
};

/** @param {Element | null | undefined} element @param {FeatureIdentity} candidate */
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

/** @param {FeatureSvgRoot | null | undefined} svg @param {(Record<string, unknown> | null | undefined)[] | null} [features] */
const buildRenderedFeatureIndex = (svg, features = []) => {
  /** @type {Map<string, FeatureIdentity[]>} */
  const metadataByRenderedId = new Map();
  (Array.isArray(features) ? features : []).forEach((feature) => {
    const candidate = renderedFeatureIdentity(feature);
    if (!candidate.usable || getFeatureElements(svg, candidate.renderedId.value).length === 0) return;
    const values = metadataByRenderedId.get(candidate.renderedId.value) || [];
    values.push(candidate);
    metadataByRenderedId.set(candidate.renderedId.value, values);
  });

  /** @type {FeatureIdentity[]} */
  const candidates = [];
  metadataByRenderedId.forEach((values, renderedId) => {
    if (values.length !== 1) return;
    const elements = getFeatureElements(svg, renderedId);
    if (elements.every((element) => domIdentityAgrees(element, values[0]))) {
      candidates.push(values[0]);
    }
  });

  /** @type {Map<string, FeatureIdentity | null>} */
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

/** @param {OrthogroupMember} member @param {FeatureIdentity[] | null | undefined} renderedIndex */
const renderedFeatureIdForMember = (member, renderedIndex) => {
  const identity = featureIdentity(member);
  if (!identity.usable) return '';
  const matches = (Array.isArray(renderedIndex) ? renderedIndex : [])
    .filter((candidate) => identityMatches(identity, candidate, { includeRendered: true }));
  return matches.length === 1 ? matches[0].renderedId.value : '';
};

/**
 * The reactive state fields the similarity-group editor reads and writes.
 * A `{ value }` field is a Vue ref; the override maps and `linearSeqs` are
 * reactive objects.
 * @typedef {object} OrthogroupEditorState
 * @property {{ value: OrthogroupGroup[] }} orthogroups
 * @property {{ linear: OrthogroupEditorDrawing }} drawings Similarity groups are the Linear drawing's (registry: `config.webEdits`, linear only).
 * @property {{ value: string }} selectedOrthogroupId
 * @property {{ value: string }} orthogroupSearch
 * @property {{ value: string }} orthogroupSortMode
 * @property {{ value: Record<string, any> | null }} clickedFeature
 * @property {{ value: boolean }} showRightDrawer
 * @property {{ value: string }} rightDrawerTab
 * @property {{ value: HTMLElement | null }} svgContainer
 * @property {Record<string, any>[]} linearSeqs
 * @property {{ value: Record<string, any>[] }} extractedFeatures
 * @property {{ value: Record<string, any>[] }} biologicalFeatures
 * @property {() => { status: string, reason: string } | null} [sessionOperationAvailability]
 */

/** @import { GroupNamesSlice, GroupNamesView } from '../services/drawing/group-names.js' */
/**
 * The group names of a drawing (`DrawingState` of state.js), and the handle their writers take.
 * @typedef {GroupNamesView & { readonly slices: { readonly groupNames: GroupNamesSlice } }} OrthogroupEditorDrawing
 */

/**
 * @typedef {() => { records?: { recordKey?: unknown }[] } | null} CommittedRequestPort
 *   The request of the shown Result, whose `records[i].recordKey` is the card uid of record index i.
 */

/**
 * @param {{ state: OrthogroupEditorState, getCommittedRequest?: CommittedRequestPort }} options
 */
export const createOrthogroupEditor = ({ state, getCommittedRequest = () => null }) => {
  const {
    orthogroups,
    selectedOrthogroupId,
    orthogroupSearch,
    orthogroupSortMode,
    clickedFeature,
    showRightDrawer,
    rightDrawerTab,
    svgContainer,
    linearSeqs,
    extractedFeatures,
    biologicalFeatures
  } = state;

  const copyFeedbackByKey = reactive({});
  /** @type {Map<string, number>} */
  const copyFeedbackTimers = new Map();
  /** @type {Map<string, number>} */
  const copyOperationIds = new Map();
  let nextCopyOperationId = 0;

  /**
   * @param {GroupOrId} groupOrId
   * @param {unknown} sequenceKind
   * @param {OrthogroupMember | null} [member]
   * @param {string} [context]
   */
  const orthogroupCopyFeedbackKey = (
    groupOrId,
    sequenceKind,
    member = null,
    context = 'drawer'
  ) => {
    const group = typeof groupOrId === 'string' ? getOrthogroupById(groupOrId) : groupOrId;
    const groupId = group ? orthogroupIdValue(group) : normalizeText(groupOrId);
    return JSON.stringify([
      normalizeText(context),
      member ? 'member' : 'group',
      groupId,
      member ? memberCopyFeedbackIdentity(member) : '',
      sequenceKindLabel(sequenceKind)
    ]);
  };

  /** @param {string} key @param {number} operationId @param {string} status */
  const updateCopyFeedback = (key, operationId, status) => {
    if (copyOperationIds.get(key) !== operationId) return;
    const previousTimer = copyFeedbackTimers.get(key);
    if (previousTimer !== undefined) window.clearTimeout(previousTimer);
    copyFeedbackByKey[key] = { operationId, status };
    const timer = window.setTimeout(() => {
      if (
        copyOperationIds.get(key) !== operationId ||
        copyFeedbackByKey[key]?.operationId !== operationId
      ) return;
      delete copyFeedbackByKey[key];
      copyFeedbackTimers.delete(key);
      copyOperationIds.delete(key);
    }, COPY_FEEDBACK_DURATION_MS);
    copyFeedbackTimers.set(key, timer);
  };

  /** @param {string} key */
  const beginCopyFeedback = (key) => {
    const previousTimer = copyFeedbackTimers.get(key);
    if (previousTimer !== undefined) window.clearTimeout(previousTimer);
    copyFeedbackTimers.delete(key);
    delete copyFeedbackByKey[key];
    nextCopyOperationId += 1;
    copyOperationIds.set(key, nextCopyOperationId);
    return nextCopyOperationId;
  };

  const clearCopyFeedback = () => {
    copyFeedbackTimers.forEach((timer) => window.clearTimeout(timer));
    copyFeedbackTimers.clear();
    copyOperationIds.clear();
    Object.keys(copyFeedbackByKey).forEach((key) => delete copyFeedbackByKey[key]);
  };

  /** @param {string} text @param {string} key */
  const copyWithFeedback = async (text, key) => {
    if (!text) return false;
    const operationId = beginCopyFeedback(key);
    try {
      await copyTextToClipboard(text);
      updateCopyFeedback(key, operationId, 'success');
      return true;
    } catch {
      updateCopyFeedback(key, operationId, 'failure');
      return false;
    }
  };

  /**
   * @param {GroupOrId} groupOrId
   * @param {unknown} sequenceKind
   * @param {unknown} originalLabel
   * @param {OrthogroupMember | null} [member]
   * @param {string} [context]
   */
  const orthogroupCopyFeedbackLabel = (
    groupOrId,
    sequenceKind,
    originalLabel,
    member = null,
    context = 'drawer'
  ) => {
    const key = orthogroupCopyFeedbackKey(groupOrId, sequenceKind, member, context);
    const status = copyFeedbackByKey[key]?.status;
    return (status && COPY_FEEDBACK_LABELS[status]) || String(originalLabel || 'Copy');
  };

  if (typeof watch === 'function') {
    watch(
      [orthogroups, selectedOrthogroupId, clickedFeature, showRightDrawer, rightDrawerTab].filter(Boolean),
      clearCopyFeedback,
      { flush: 'sync' }
    );
    // The open Similarity groups list shows the selected group's row, whoever
    // selected it: a new alignment plan (OV-381), Show in editor, or the user.
    // The tab panel's v-show applies after post watchers, hence nextTick.
    watch(
      () => [selectedOrthogroupId.value, showRightDrawer?.value && rightDrawerTab?.value === 'orthogroups'],
      async (/** @type {[unknown, unknown]} */ [id, listShown]) => {
        if (!id || !listShown || typeof document === 'undefined') return;
        await nextTick();
        document.querySelector(`[data-similarity-group-row="${CSS.escape(String(id))}"]`)
          ?.scrollIntoView({ block: 'nearest' });
      }
    );
  }

  /** @type {{ readonly value: { id: string, group: OrthogroupGroup }[] }} */
  const resolvableOrthogroupEntries = computed(() => {
    const groups = Array.isArray(orthogroups.value) ? orthogroups.value : [];
    // Track absent aliases too; adding a conflicting ID must invalidate the list.
    for (const group of groups) {
      if (group && typeof group === 'object') {
        for (const key of ORTHOGROUP_ID_KEYS) void (key in group);
      }
    }
    return uniqueOrthogroupEntries(groups);
  });

  /** @param {unknown} orthogroupId @returns {OrthogroupGroup | null} */
  const getOrthogroupById = (orthogroupId) => {
    const id = normalizeText(orthogroupId);
    if (!id) return null;
    return resolvableOrthogroupEntries.value
      .find((/** @type {{ id: string, group: OrthogroupGroup }} */ entry) => entry.id === id)?.group || null;
  };

  /** @param {GroupOrId} groupOrId */
  const resolveOrthogroupName = (groupOrId) => {
    const drawing = state.drawings.linear;
    const group = typeof groupOrId === 'string' ? getOrthogroupById(groupOrId) : groupOrId;
    const id = group ? orthogroupIdValue(group) : normalizeText(groupOrId);
    if (!id) return '';
    // A Similarity group id is data, so the rename tables are read by own key.
    return normalizeText(ownEntry(drawing.orthogroupNameOverrides, id)) || normalizeText(group?.name) || id;
  };

  /** @param {GroupOrId} groupOrId */
  const resolveOrthogroupDescription = (groupOrId) => {
    const drawing = state.drawings.linear;
    const group = typeof groupOrId === 'string' ? getOrthogroupById(groupOrId) : groupOrId;
    const id = group ? orthogroupIdValue(group) : normalizeText(groupOrId);
    if (!id) return '';
    return normalizeText(ownEntry(drawing.orthogroupDescriptionOverrides, id)) || normalizeText(group?.description);
  };

  // The id of a group as the rename tables key it.
  /** @param {GroupOrId} groupOrId */
  const renameKey = (groupOrId) => (
    typeof groupOrId === 'string' ? normalizeText(groupOrId) : orthogroupIdValue(groupOrId)
  );

  // The text of the editor's Description field: the user's own edit of the
  // group, and nothing for a group without one (Python's description is not
  // put into the field).
  /** @param {GroupOrId} groupOrId */
  const orthogroupDescriptionFieldText = (groupOrId) => (
    ownEntry(state.drawings.linear.orthogroupDescriptionOverrides, renameKey(groupOrId)) || ''
  );

  /** @param {GroupOrId} groupOrId */
  const orthogroupScope = (groupOrId) => {
    const group = typeof groupOrId === 'string' ? getOrthogroupById(groupOrId) : groupOrId;
    return normalizeGroupMetadataScope(group?.scope);
  };

  /** @param {GroupOrId} groupOrId */
  const orthogroupScopeLabel = (groupOrId) => groupMetadataScopeLabel(orthogroupScope(groupOrId));

  // The ids the rename tables hold. `Object.keys` is the Vue dependency: it
  // follows a key that is added or removed, and not a change of its text, so
  // typing a description does not rebuild the list rows.
  /** @type {{ readonly value: Set<string> }} */
  const renamedOrthogroupIds = computed(() => {
    const drawing = state.drawings.linear;
    return new Set([
      ...Object.keys(drawing.orthogroupNameOverrides), ...Object.keys(drawing.orthogroupDescriptionOverrides)
    ]);
  });

  /** @param {GroupOrId} groupOrId */
  const isOrthogroupRenamed = (groupOrId) => {
    const id = renameKey(groupOrId);
    return Boolean(id) && renamedOrthogroupIds.value.has(id);
  };

  /** @type {{ readonly value: number }} */
  const orthogroupCount = computed(() => resolvableOrthogroupEntries.value.length);

  /** @type {{ readonly value: OrthogroupGroup[] }} */
  const filteredOrthogroups = computed(() => {
    const query = normalizeLower(orthogroupSearch.value);
    const groups = resolvableOrthogroupEntries.value.map(({ group }) => group);
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

  /** @type {{ readonly value: OrthogroupGroup | null }} */
  const selectedOrthogroup = computed(() => {
    const selected = getOrthogroupById(selectedOrthogroupId.value);
    if (selected) return selected;
    return filteredOrthogroups.value[0] || null;
  });

  const orthogroupRows = computed(() => {
    /** @type {(count: string, noun: string) => string} */
    const counted = (count, noun) => `${count} ${count === '1' ? noun : `${noun}s`}`;
    return filteredOrthogroups.value.map((group) => ({
      group, id: group.id,
      name: resolveOrthogroupName(group),
      // Only a scope other than the usual cross-record group is worth a word on the row.
      scopeLabel: orthogroupScope(group) === 'global' ? '' : orthogroupScopeLabel(group),
      countsLabel: [
        counted(toDisplayString(group.member_count || (group.members ? group.members.length : 0)), 'gene'),
        counted(toDisplayString(group.record_coverage_count || 0), 'record')
      ].join(' · '),
      renamed: isOrthogroupRenamed(group)
    }));
  });

  const featureSequenceLookup = computed(() => {
    const hasBiologicalFeatures = Array.isArray(biologicalFeatures?.value) && biologicalFeatures.value.length > 0;
    const features = hasBiologicalFeatures
      ? biologicalFeatures.value
      : (Array.isArray(extractedFeatures?.value) ? extractedFeatures.value : []);
    const renderedIdentities = (Array.isArray(extractedFeatures?.value) ? extractedFeatures.value : [])
      .map(renderedFeatureIdentity)
      .filter((identity) => identity.usable);
    const findRendered = indexFeatureIdentities(renderedIdentities, (identity) => identity);
    const entries = features.map((feature) => {
      const identity = hasBiologicalFeatures
        ? featureIdentity(feature, { allowLegacySvgStable: true })
        : renderedFeatureIdentity(feature);
      const renderedMatches = findRendered(identity);
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
        get entry() {
          return {
            nucleotideSequence: firstSequenceText(feature?.nucleotideSequence, feature?.nucleotide_sequence),
            aminoAcidSequence: firstSequenceText(feature?.aminoAcidSequence, feature?.amino_acid_sequence),
            sequenceFeature: feature,
            sequenceWarnings: Array.isArray(feature?.sequence_warnings)
              ? feature.sequence_warnings
              : (Array.isArray(feature?.sequenceWarnings) ? feature.sequenceWarnings : [])
          };
        }
      };
    });
    return indexFeatureIdentities(entries, (entry) => entry.identity);
  });

  /** @param {OrthogroupMember} member */
  const enrichOrthogroupMember = (member) => {
    const memberSequenceFeature = member?.sequenceFeature || null;
    const memberIdentity = featureIdentity(member);
    const matches = featureSequenceLookup.value(memberIdentity);
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

  /** @param {GroupOrId} groupOrId */
  const getEnrichedOrthogroupMembers = (groupOrId) => {
    const group = typeof groupOrId === 'string' ? getOrthogroupById(groupOrId) : groupOrId;
    return getGroupMembers(group).map(enrichOrthogroupMember);
  };

  // A member's record index is a position in the records of the Result's
  // request, which names each record's card (record selection: a card turned
  // OFF after Generate does not shift the labels).
  /** @param {unknown} members */
  const groupOrthogroupMembersByRecord = (members) => {
    const requestRecords = getCommittedRequest()?.records || [];
    /** @param {number} recordIndex */
    const cardOf = (recordIndex) => {
      const uid = String(requestRecords[recordIndex]?.recordKey ?? '');
      return uid ? linearSeqs.find((sequence) => sequence.uid === uid) : undefined;
    };
    const names = linearRecordNames(linearSeqs);
    /** @type {Map<number, OrthogroupMember[]>} */
    const byRecord = new Map();
    (Array.isArray(members) ? members : []).forEach((member) => {
      const recordIndex = nonnegativeIntegerAliasStatus(member, RECORD_INDEX_KEYS);
      // `value` is a number when the alias is valid and supplied.
      const key = recordIndex.valid && recordIndex.supplied ? /** @type {number} */ (recordIndex.value) : -1;
      if (!byRecord.has(key)) byRecord.set(key, []);
      // Set on the line above when absent.
      /** @type {OrthogroupMember[]} */ (byRecord.get(key)).push(member);
    });
    return Array.from(byRecord.entries())
      .sort((left, right) => left[0] - right[0])
      .map(([recordIndex, members]) => ({
        recordIndex,
        // The record as every panel names it (OV-369, OV-436).
        recordLabel: recordIndex >= 0
          ? names.get(String(cardOf(recordIndex)?.uid ?? ''))?.label || `Record ${recordIndex + 1}`
          : 'Record',
        members
      }));
  };

  /** @param {GroupOrId} groupOrId @param {unknown} sequenceKind */
  const getOrthogroupSequenceCount = (groupOrId, sequenceKind) =>
    getEnrichedOrthogroupMembers(groupOrId).filter((member) => getMemberSequence(member, sequenceKind)).length;

  /** @param {OrthogroupMember} member @param {unknown} sequenceKind */
  const hasOrthogroupMemberSequence = (member, sequenceKind) =>
    Boolean(getMemberSequence(enrichOrthogroupMember(member), sequenceKind));

  // OV-920: what the drawer shows of the selected group, built when the group
  // or a sequence source changes. Every render of the app reads it, with the
  // drawer closed too, so the template asks no helper that enriches members.
  /**
   * @typedef {object} SelectedOrthogroupDetail
   * @property {OrthogroupGroup} group
   * @property {unknown} id the id the drawer's controls pass (`group.id`)
   * @property {OrthogroupMember[]} members the enriched members
   * @property {ReturnType<typeof groupOrthogroupMembersByRecord>} membersByRecord
   * @property {{ nt: boolean, aa: boolean }} hasSequence whether a member has an nt or aa sequence
   * @property {(member: OrthogroupMember, sequenceKind: 'nt' | 'aa') => boolean} hasMemberSequence
   */
  /** @type {{ readonly value: SelectedOrthogroupDetail | null }} */
  const selectedOrthogroupDetail = computed(() => {
    const group = selectedOrthogroup.value;
    if (!group) return null;
    const members = getEnrichedOrthogroupMembers(group);
    /** @type {Map<OrthogroupMember, { nt: boolean, aa: boolean }>} */
    const memberSequences = new Map(members.map((member) => [member, {
      nt: hasOrthogroupMemberSequence(member, 'nt'),
      aa: hasOrthogroupMemberSequence(member, 'aa')
    }]));
    return {
      group,
      id: group.id,
      members,
      membersByRecord: groupOrthogroupMembersByRecord(members),
      hasSequence: {
        nt: members.some((member) => getMemberSequence(member, 'nt')),
        aa: members.some((member) => getMemberSequence(member, 'aa'))
      },
      hasMemberSequence: (/** @type {OrthogroupMember} */ member, /** @type {'nt' | 'aa'} */ sequenceKind) => (
        Boolean(memberSequences.get(member)?.[sequenceKind])
      )
    };
  });

  /** @param {GroupOrId} groupOrId @param {unknown} sequenceKind */
  const buildOrthogroupFasta = (groupOrId, sequenceKind) => {
    const group = typeof groupOrId === 'string' ? getOrthogroupById(groupOrId) : groupOrId;
    return getEnrichedOrthogroupMembers(group)
      .map((member) => memberFastaText(member, sequenceKind))
      .filter(Boolean)
      .join('');
  };

  /** @param {OrthogroupMember} member @param {unknown} sequenceKind */
  const buildOrthogroupMemberFasta = (member, sequenceKind) =>
    memberFastaText(enrichOrthogroupMember(member), sequenceKind);

  /** @param {GroupOrId} groupOrId @param {unknown} sequenceKind */
  const orthogroupSequenceFilename = (groupOrId, sequenceKind) => {
    const group = typeof groupOrId === 'string' ? getOrthogroupById(groupOrId) : groupOrId;
    const id = group ? orthogroupIdValue(group) : normalizeText(groupOrId);
    const name = makeSafeFilename(resolveOrthogroupName(group) || id, id || 'orthogroup');
    const stem = makeSafeFilename(`${id || 'orthogroup'}_${name}_${sequenceKindLabel(sequenceKind)}`);
    return `${stem}.${sequenceExtension(sequenceKind)}`;
  };

  /** @param {OrthogroupMember} member @param {unknown} sequenceKind @param {GroupOrId} [groupOrId] */
  const orthogroupMemberSequenceFilename = (member, sequenceKind, groupOrId = selectedOrthogroup.value) => {
    const group = typeof groupOrId === 'string' ? getOrthogroupById(groupOrId) : groupOrId;
    const id = (group ? orthogroupIdValue(group) : normalizeText(groupOrId)) || 'orthogroup';
    const stableId = memberStableFeatureId(member);
    const memberId = normalizeText(
      resolveDisplayProteinId(
        // A feature row object when present; the resolver reads it by key.
        /** @type {Record<string, unknown> | null | undefined} */ (member?.sequenceFeature),
        member,
        memberLocationText(member)
      ) ||
      (isInternalProteinDisplayId(stableId) ? '' : stableId) ||
      'member'
    );
    const stem = makeSafeFilename(`${id}_${memberId}_${sequenceKindLabel(sequenceKind)}`);
    return `${stem}.${sequenceExtension(sequenceKind)}`;
  };

  /** @param {GroupOrId} [groupOrId] @param {unknown} [sequenceKind] @param {string} [context] */
  const copyOrthogroupSequences = async (
    groupOrId = selectedOrthogroup.value,
    sequenceKind = 'nt',
    context = 'drawer'
  ) => {
    const text = buildOrthogroupFasta(groupOrId, sequenceKind);
    const key = orthogroupCopyFeedbackKey(groupOrId, sequenceKind, null, context);
    return copyWithFeedback(text, key);
  };

  /** @param {GroupOrId} [groupOrId] @param {unknown} [sequenceKind] */
  const downloadOrthogroupSequences = (groupOrId = selectedOrthogroup.value, sequenceKind = 'nt') => {
    const text = buildOrthogroupFasta(groupOrId, sequenceKind);
    if (!text) return;
    downloadTextFile(
      orthogroupSequenceFilename(groupOrId, sequenceKind),
      text,
      'text/plain;charset=utf-8'
    );
  };

  /** @param {OrthogroupMember} member @param {unknown} [sequenceKind] @param {GroupOrId} [groupOrId] @param {string} [context] */
  const copyOrthogroupMemberSequence = async (
    member,
    sequenceKind = 'nt',
    groupOrId = selectedOrthogroup.value,
    context = 'drawer'
  ) => {
    const text = buildOrthogroupMemberFasta(member, sequenceKind);
    const key = orthogroupCopyFeedbackKey(groupOrId, sequenceKind, member, context);
    return copyWithFeedback(text, key);
  };

  /** @param {OrthogroupMember} member @param {unknown} [sequenceKind] @param {GroupOrId} [groupOrId] */
  const downloadOrthogroupMemberSequence = (member, sequenceKind = 'nt', groupOrId = selectedOrthogroup.value) => {
    const text = buildOrthogroupMemberFasta(member, sequenceKind);
    if (!text) return;
    downloadTextFile(
      orthogroupMemberSequenceFilename(member, sequenceKind, groupOrId),
      text,
      'text/plain;charset=utf-8'
    );
  };

  /** @param {unknown} orthogroupId */
  const selectOrthogroup = (orthogroupId) => {
    const id = normalizeText(orthogroupId);
    if (!id || !getOrthogroupById(id)) return false;
    selectedOrthogroupId.value = id;
    return true;
  };

  /** @param {unknown} orthogroupId @param {unknown} value */
  const setOrthogroupNameOverride = (orthogroupId, value) => {
    const drawing = state.drawings.linear;
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const id = normalizeText(orthogroupId);
    if (!id) return;
    const group = getOrthogroupById(id);
    if (!group) return;
    const normalized = normalizeText(value);
    setGroupName(drawing.slices.groupNames, id, normalized === normalizeText(group?.name) ? '' : normalized);
  };

  /** @param {unknown} orthogroupId @param {unknown} value */
  const setOrthogroupDescriptionOverride = (orthogroupId, value) => {
    const drawing = state.drawings.linear;
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const id = normalizeText(orthogroupId);
    if (!id) return;
    const group = getOrthogroupById(id);
    if (!group) return;
    const normalized = normalizeText(value);
    setGroupDescription(drawing.slices.groupNames, id, normalized === normalizeText(group?.description) ? '' : normalized);
  };

  /** @param {unknown} [orthogroupId] */
  const resetOrthogroupRename = (orthogroupId = selectedOrthogroupId.value) => {
    const drawing = state.drawings.linear;
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const id = normalizeText(orthogroupId);
    if (!id) return;
    clearGroupRename(drawing.slices.groupNames, id);
  };

  // Dormant names are listed and cleared here (D-21).
  const orthogroupDormantNames = computed(() => Object.values(state.drawings.linear.orthogroupDormantOverrides || {})
    .map((entry) => normalizeText(entry?.name) || normalizeText(entry?.description)).filter(Boolean));
  const clearOrthogroupDormantOverrides = () => {
    const drawing = state.drawings.linear;
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    GROUP_NAMES_SLICE.install(drawing.slices.groupNames, { orthogroupDormantOverrides: {} });
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

  return {
    orthogroupCount,
    filteredOrthogroups,
    selectedOrthogroup,
    selectedOrthogroupDetail,
    getEnrichedOrthogroupMembers,
    groupOrthogroupMembersByRecord,
    // The location text of a member row of `groupOrthogroupMembersByRecord`.
    memberRangeText: locationRangeText,
    getOrthogroupById,
    resolveOrthogroupName,
    resolveOrthogroupDescription,
    orthogroupDescriptionFieldText,
    orthogroupScope,
    orthogroupScopeLabel,
    orthogroupRows,
    isOrthogroupRenamed,
    getOrthogroupSequenceCount,
    hasOrthogroupMemberSequence,
    orthogroupCopyFeedbackLabel,
    clearCopyFeedback,
    copyOrthogroupSequences,
    downloadOrthogroupSequences,
    copyOrthogroupMemberSequence,
    downloadOrthogroupMemberSequence,
    selectOrthogroup,
    setOrthogroupNameOverride,
    setOrthogroupDescriptionOverride,
    resetOrthogroupRename,
    orthogroupDormantNames,
    clearOrthogroupDormantOverrides,
    highlightOrthogroupById,
    clearOrthogroupHighlight,
  };
};
