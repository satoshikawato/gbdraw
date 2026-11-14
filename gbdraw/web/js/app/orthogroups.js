// @ts-check
import {
  FEATURE_SELECTOR,
  getFeatureElements,
  getFeatureIdentity
} from './feature-editor/svg-actions.js';
import {
  groupMetadataScopeLabel,
  normalizeGroupMetadataScope
} from '../services/losat-normalization.js';
import {
  isInternalProteinDisplayId,
  locationPartBound,
  locationRangeText,
  normalizeMemberStrand,
  resolveDisplayProteinId
} from '../services/feature-utils.js';
import { linearRecordNames } from '../services/linear-sources.js';
import { featureSequenceReads, unreadSequence } from '../services/catalog-sequences.js';
import { memberFastaFeature } from '../services/feature-sequence-fasta.js';
import { ownEntry } from '../utils/own-entry.js';
import { toRaw } from '../utils/vue-raw.js';
import { clearDormantGroupNames, clearGroupRename, setGroupDescription, setGroupName } from '../services/drawing/group-names.js';
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
  uniqueOrthogroupEntries,
  withLegacyCanonicalIdentity
} from '../services/feature-identity.js';
import { msg, t, tm } from '../utils/i18n.js';

/** @import { AliasSource, FeatureIdentity } from '../services/feature-identity.js' */
/** @import { FeatureSvgRoot } from '../services/feature-dom.js' */
/** @import { OrthogroupGroup, OrthogroupMember } from '../services/orthogroup-feature-metadata.js' */
/** @import { FeatureLike } from '../services/feature-utils.js' */
/** @import { FeatureSequenceRead, SequenceReadState } from '../services/catalog-sequences.js' */
/** @import { MemberAction } from './sequence-reads.js' */
/** @typedef {OrthogroupGroup | string | null | undefined} GroupOrId A group row, or the id of one. */

export { resolveUniqueOrthogroupMemberForFeature } from '../services/feature-identity.js';

const { computed, nextTick, reactive, watch, toDisplayString } = window.Vue;

const COPY_FEEDBACK_DURATION_MS = 1500;
/** What a Copy button says for a moment after its text was copied. */
const COPIED_LABEL = msg('Copied!');
/**
 * Where a Similarity group's sequence actions show their read: the feature popup and the Similarity groups drawer.
 * @typedef {'feature-popup' | 'drawer'} MemberReadContext
 */

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

/** @param {OrthogroupMember | null | undefined} member */
const memberLocationText = (member) => {
  const start = locationPartBound(member, 'start');
  const end = locationPartBound(member, 'end');
  if (!Number.isFinite(start) || !Number.isFinite(end)) return '';
  const normalizedStrand = normalizeMemberStrand(member?.strand);
  const strand = normalizedStrand ? `(${normalizedStrand})` : '';
  return `${start + 1}-${end}${strand}`;
};

/**
 * The read of an enriched member's sequences: the identity and bases of the catalog row it resolved to, or null
 * when it resolved to none or the row has no identity. Such a member is the one a sequence action can read, so it
 * is also what "has a sequence" means before the read (D5): whether it holds one is known only by reading it.
 * @param {OrthogroupMember} member
 * @returns {FeatureSequenceRead | null}
 */
const memberSequenceRead = (member) => {
  const row = member.sequenceFeature;
  return row && typeof row === 'object' ? featureSequenceReads([/** @type {FeatureLike} */ (row)])[0] ?? null : null;
};

/**
 * The "Copied!" each Copy button shows for a moment, by feedback key. An action supersedes the earlier one on its
 * key, whose end then shows nothing. The read itself is shown by the context's read state (`orthogroupMemberRead`).
 */
const makeCopiedFeedback = () => {
  /** @type {Record<string, true>} */
  const copiedKeys = reactive({});
  /** @type {Map<string, number>} */
  const timers = new Map();
  /** @type {Map<string, number>} */
  const operationIds = new Map();
  let nextOperationId = 0;
  /** @param {string} key */
  const end = (key) => {
    const timer = timers.get(key);
    if (timer !== undefined) window.clearTimeout(timer);
    timers.delete(key);
    delete copiedKeys[key];
  };
  /** @param {string} key */
  const begin = (key) => {
    end(key);
    nextOperationId += 1;
    operationIds.set(key, nextOperationId);
    return nextOperationId;
  };
  /** @param {string} key @param {number} operationId */
  const copied = (key, operationId) => {
    if (operationIds.get(key) !== operationId) return;
    end(key);
    copiedKeys[key] = true;
    timers.set(key, window.setTimeout(() => {
      if (operationIds.get(key) === operationId) end(key);
    }, COPY_FEEDBACK_DURATION_MS));
  };
  const clear = () => {
    [...timers.keys()].forEach(end);
    operationIds.clear();
    Object.keys(copiedKeys).forEach((key) => delete copiedKeys[key]);
  };
  /** @param {string} key */
  const isCopied = (key) => copiedKeys[key] === true;
  return { begin, copied, clear, isCopied };
};

/**
 * The sequence actions of a Similarity group's members: each reads through `memberAction` (the one owner of a
 * popup's member read, app/sequence-reads.js) with the read state of its context, and a Copy shows "Copied!" on its
 * button. Clearing gives each context a new state, so a read still running writes only to the state it started
 * with, which no popup shows any more.
 * @param {MemberAction} memberAction
 */
const makeMemberSequenceActions = (memberAction) => {
  const feedback = makeCopiedFeedback();
  /** @type {Record<MemberReadContext, SequenceReadState>} */
  const reads = reactive({ 'feature-popup': unreadSequence(), drawer: unreadSequence() });
  /**
   * @param {{ key: string, context: MemberReadContext, action: 'copy' | 'download', sequenceKind: unknown,
   *   members: readonly OrthogroupMember[], filename: string }} request
   */
  const run = async ({ key, context, action, sequenceKind, members, filename }) => {
    const rows = members.map((member) => ({
      read: memberSequenceRead(member),
      fastaFeature: memberFastaFeature(member, /** @type {FeatureLike | null} */ (member.sequenceFeature ?? null), memberLocationText(member))
    }));
    const operationId = feedback.begin(key);
    const done = await memberAction(reads[context], rows, sequenceKindLabel(sequenceKind), action, filename);
    if (done && action === 'copy') feedback.copied(key, operationId);
    return done;
  };
  return {
    run,
    /** @param {MemberReadContext} context */
    readOf: (context) => reads[context],
    isCopied: feedback.isCopied,
    clear: () => {
      feedback.clear();
      reads['feature-popup'] = unreadSequence();
      reads.drawer = unreadSequence();
    }
  };
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

/**
 * The lookup of the catalog rows a feature identity names: the biological rows, else the drawn ones, each under its
 * identity with the drawn ID its one drawn row gives it. A row Load read again for a Session before 40 keeps that
 * Session's names (record index and stable ID) and its own edit keys; a member Load gave canonical names finds it
 * under the same names, read for the row from the committed request by the one rule (`withLegacyCanonicalIdentity`,
 * OV-3539).
 * @param {readonly Record<string, unknown>[] | null | undefined} biological
 * @param {readonly Record<string, unknown>[]} extracted
 * @param {readonly unknown[]} recordKeys The committed request's record keys, in record order.
 */
const catalogRowLookup = (biological, extracted, recordKeys) => {
  const hasBiologicalFeatures = Array.isArray(biological) && biological.length > 0;
  const features = hasBiologicalFeatures ? biological : extracted;
  const renderedIdentities = extracted
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
    return { identity: withLegacyCanonicalIdentity(resolvedIdentity, recordKeys), feature };
  });
  return indexFeatureIdentities(entries, (entry) => entry.identity);
};

/**
 * The catalog row lookup (`catalogRowLookup`) of the shown rows, rebuilt only when the rows or the committed request's
 * record keys change. It depends on the two array refs, which every writer replaces (services/config.js,
 * services/history-snapshot.js), and reads the rows through `toRaw`, as state.js does, so Vue records no dependency on
 * each alias key of every row that `featureIdentity` reads. The record keys stay the same list while a committed change
 * keeps them (a placement, a setting), so such a change does not read every row again.
 * @param {Pick<WebState, 'extractedFeatures' | 'biologicalFeatures'>} rows
 * @param {CommittedRequestPort} getCommittedRequest
 */
const catalogRowLookupOf = ({ extractedFeatures, biologicalFeatures }, getCommittedRequest) => {
  const recordKeys = computed((/** @type {unknown[] | undefined} */ previous) => {
    const keys = (getCommittedRequest()?.records || []).map((record) => record.recordKey);
    return previous?.length === keys.length && keys.every((key, index) => key === previous[index]) ? previous : keys;
  });
  return computed(() => catalogRowLookup(
    toRaw(biologicalFeatures.value),
    Array.isArray(extractedFeatures.value) ? toRaw(extractedFeatures.value) : [],
    recordKeys.value
  ));
};

/**
 * The drawn element ID of an alignment anchor's feature, from the one catalog row the anchor names (`catalogRowLookup`).
 * The alignment binds a Result an older gbdraw drew through it: a record drawn reverse complemented there draws its
 * features under their hash in the drawn record, not under the stable ID the anchor names (OV-3538).
 * @param {{ readonly value: ReturnType<typeof catalogRowLookup> }} lookup
 * @returns {(anchor: AliasSource) => string[]}
 */
const anchorDrawnFeatureIds = (lookup) => (anchor) => {
  const matches = lookup.value(featureIdentity(anchor));
  const renderedId = matches.length === 1 ? matches[0].identity.renderedId : null;
  return renderedId !== null && renderedId.supplied ? [renderedId.value] : [];
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
 * The app state the similarity-group editor reads and writes. `activeDrawing` is the shown drawing; Similarity groups
 * are the Linear drawing's (registry: `config.webEdits`, linear only), so every reader here uses it only when `mode` is
 * `'linear'`.
 * @typedef {Pick<WebState, 'orthogroups' | 'activeDrawing' | 'selectedOrthogroupId' | 'orthogroupSearch'
 *   | 'orthogroupSortMode' | 'clickedFeature' | 'showRightDrawer' | 'rightDrawerTab' | 'svgContainer'
 *   | 'extractedFeatures' | 'biologicalFeatures' | 'sessionOperationAvailability'>} OrthogroupEditorState
 */

/** @import { GroupNamesSlice, GroupNamesView } from '../services/drawing/group-names.js' */
/**
 * The group names of a drawing (`DrawingState` of state.js), and the handle their writers take.
 * @typedef {GroupNamesView & { readonly slices: { readonly groupNames: GroupNamesSlice } }} OrthogroupEditorDrawing
 */
/** @import { LinearSeq, WebState } from '../state.js' */
/**
 * The shown drawing: its group names, its mode, and its record rows.
 * @typedef {OrthogroupEditorDrawing & { readonly mode: string, readonly linearSeqs: readonly LinearSeq[] }} OrthogroupEditorShownDrawing
 */

/**
 * @typedef {() => { records?: { recordKey?: unknown }[] } | null} CommittedRequestPort
 *   The request of the shown Result, whose `records[i].recordKey` is the card uid of record index i.
 */

/**
 * @param {{ state: OrthogroupEditorState, getCommittedRequest?: CommittedRequestPort,
 *   memberAction: MemberAction }} options
 *   `memberAction` reads members' sequences from Python for the shown Result and copies or downloads their FASTA,
 *   with the read state a popup shows: the one owner the Similarity group popup of a match uses too.
 */
export const createOrthogroupEditor = ({ state, getCommittedRequest = () => null, memberAction }) => {
  const {
    orthogroups,
    selectedOrthogroupId,
    orthogroupSearch,
    orthogroupSortMode,
    clickedFeature,
    showRightDrawer,
    rightDrawerTab,
    svgContainer,
    extractedFeatures
  } = state;

  const sequenceActions = makeMemberSequenceActions(memberAction);
  const { readOf: orthogroupMemberRead, clear: clearCopyFeedback } = sequenceActions;

  /**
   * The feedback key of one sequence action's button.
   * @param {GroupOrId} groupOrId
   * @param {unknown} sequenceKind
   * @param {OrthogroupMember | null} [member]
   * @param {string} [context]
   * @param {'copy' | 'download'} [action]
   */
  const orthogroupCopyFeedbackKey = (groupOrId, sequenceKind, member = null, context = 'drawer', action = 'copy') => {
    const group = typeof groupOrId === 'string' ? getOrthogroupById(groupOrId) : groupOrId;
    const groupId = group ? orthogroupIdValue(group) : normalizeText(groupOrId);
    return JSON.stringify([
      normalizeText(context),
      member ? 'member' : 'group',
      groupId,
      member ? memberCopyFeedbackIdentity(member) : '',
      sequenceKindLabel(sequenceKind),
      action
    ]);
  };

  /**
   * @param {GroupOrId} groupOrId
   * @param {unknown} sequenceKind
   * @param {unknown} originalLabel
   * @param {OrthogroupMember | null} [member]
   * @param {string} [context]
   * @param {'copy' | 'download'} [action]
   */
  const orthogroupCopyFeedbackLabel = (
    groupOrId,
    sequenceKind,
    originalLabel,
    member = null,
    context = 'drawer',
    action = 'copy'
  ) => {
    return sequenceActions.isCopied(orthogroupCopyFeedbackKey(groupOrId, sequenceKind, member, context, action))
      ? tm(COPIED_LABEL)
      : String(originalLabel || t('Copy'));
  };

  if (typeof watch === 'function') {
    watch(
      [orthogroups, selectedOrthogroupId, clickedFeature, showRightDrawer, rightDrawerTab],
      clearCopyFeedback,
      { flush: 'sync' }
    );
    // The open Similarity groups list shows the selected group's row, whoever
    // selected it: a new alignment plan (OV-381), Show in editor, or the user.
    // The tab panel's v-show applies after post watchers, hence nextTick.
    watch(
      () => [selectedOrthogroupId.value, showRightDrawer.value && rightDrawerTab.value === 'orthogroups'],
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

  // The shown drawing when it is the Linear one, else null (the group-name
  // tables belong to the Linear drawing). Resolved at each read, never cached.
  /** @returns {OrthogroupEditorShownDrawing | null} */
  const shownLinearDrawing = () => {
    const drawing = state.activeDrawing();
    return drawing.mode === 'linear' ? drawing : null;
  };

  // The group-name readers take the drawing their caller resolved: the shown
  // Linear drawing, or null for any other.
  /**
   * @param {GroupOrId} groupOrId
   * @param {OrthogroupEditorShownDrawing | null} drawing
   */
  const resolveOrthogroupName = (groupOrId, drawing) => {
    const group = typeof groupOrId === 'string' ? getOrthogroupById(groupOrId) : groupOrId;
    const id = group ? orthogroupIdValue(group) : normalizeText(groupOrId);
    if (!id) return '';
    // A Similarity group id is data, so the rename tables are read by own key.
    return normalizeText(drawing && ownEntry(drawing.orthogroupNameOverrides, id)) || normalizeText(group?.name) || id;
  };

  /**
   * @param {GroupOrId} groupOrId
   * @param {OrthogroupEditorShownDrawing | null} drawing
   */
  const resolveOrthogroupDescription = (groupOrId, drawing) => {
    const group = typeof groupOrId === 'string' ? getOrthogroupById(groupOrId) : groupOrId;
    const id = group ? orthogroupIdValue(group) : normalizeText(groupOrId);
    if (!id) return '';
    return normalizeText(drawing && ownEntry(drawing.orthogroupDescriptionOverrides, id)) || normalizeText(group?.description);
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
  const orthogroupDescriptionFieldText = (groupOrId) => {
    const drawing = shownLinearDrawing();
    return drawing === null ? '' : ownEntry(drawing.orthogroupDescriptionOverrides, renameKey(groupOrId)) || '';
  };

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
    const drawing = shownLinearDrawing();
    if (drawing === null) return new Set();
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
    const drawing = shownLinearDrawing();
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
            resolveOrthogroupName(group, drawing),
            resolveOrthogroupDescription(group, drawing),
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
        return resolveOrthogroupName(left, drawing).localeCompare(resolveOrthogroupName(right, drawing)) || compareOrthogroupId(left, right);
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
    const drawing = shownLinearDrawing();
    return filteredOrthogroups.value.map((group) => ({
      group, id: group.id,
      name: resolveOrthogroupName(group, drawing),
      // Only a scope other than the usual cross-record group is worth a word on the row.
      scopeLabel: orthogroupScope(group) === 'global' ? '' : orthogroupScopeLabel(group),
      countsLabel: [
        counted(toDisplayString(group.member_count || getGroupMembers(group).length), 'gene'),
        counted(toDisplayString(group.record_coverage_count || 0), 'record')
      ].join(' · '),
      renamed: isOrthogroupRenamed(group)
    }));
  });

  // The catalog row of a member, by its identity.
  const featureCatalogLookup = catalogRowLookupOf(state, getCommittedRequest);

  // A member with the catalog row it resolved to (`sequenceFeature`); it carries no sequence (D5).
  /** @param {OrthogroupMember} member @returns {OrthogroupMember} */
  const enrichOrthogroupMember = (member) => {
    const matches = featureCatalogLookup.value(featureIdentity(member));
    const sequenceFeature = (matches.length === 1 ? matches[0].feature : null) || member?.sequenceFeature || null;
    const displayProteinId = resolveDisplayProteinId(sequenceFeature, member, memberLocationText(member));
    return sequenceFeature ? { ...member, displayProteinId, sequenceFeature } : { ...member, displayProteinId };
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
    // The cards of the shown Linear drawing, whose Result this is.
    const drawing = shownLinearDrawing();
    const linearSeqs = drawing === null ? [] : drawing.linearSeqs;
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
          ? names.get(String(cardOf(recordIndex)?.uid ?? ''))?.label || t('Record {n}', { n: recordIndex + 1 })
          : t('Record'),
        members
      }));
  };

  // A member's sequences are known only by reading them, so a member with a catalog identity counts as having
  // them, for both kinds (`memberSequenceRead`).
  /** @param {GroupOrId} groupOrId */
  const getOrthogroupSequenceCount = (groupOrId) =>
    getEnrichedOrthogroupMembers(groupOrId).filter((member) => memberSequenceRead(member)).length;

  /** @param {OrthogroupMember} member */
  const hasOrthogroupMemberSequence = (member) => memberSequenceRead(enrichOrthogroupMember(member)) !== null;

  // OV-920: what the drawer shows of the selected group, built when the group
  // or a sequence source changes. Every render of the app reads it, with the
  // drawer closed too, so the template asks no helper that enriches members.
  /**
   * @typedef {object} SelectedOrthogroupDetail
   * @property {OrthogroupGroup} group
   * @property {unknown} id the id the drawer's controls pass (`group.id`)
   * @property {OrthogroupMember[]} members the enriched members
   * @property {ReturnType<typeof groupOrthogroupMembersByRecord>} membersByRecord
   * @property {boolean} hasSequence whether a member has a sequence to read
   * @property {(member: OrthogroupMember) => boolean} hasMemberSequence
   */
  /** @type {{ readonly value: SelectedOrthogroupDetail | null }} */
  const selectedOrthogroupDetail = computed(() => {
    const group = selectedOrthogroup.value;
    if (!group) return null;
    const members = getEnrichedOrthogroupMembers(group);
    const readable = new Set(members.filter((member) => memberSequenceRead(member)));
    return {
      group,
      id: group.id,
      members,
      membersByRecord: groupOrthogroupMembersByRecord(members),
      hasSequence: readable.size > 0,
      hasMemberSequence: (/** @type {OrthogroupMember} */ member) => readable.has(member)
    };
  });

  /** @param {GroupOrId} groupOrId @param {unknown} sequenceKind */
  const orthogroupSequenceFilename = (groupOrId, sequenceKind) => {
    const group = typeof groupOrId === 'string' ? getOrthogroupById(groupOrId) : groupOrId;
    const id = group ? orthogroupIdValue(group) : normalizeText(groupOrId);
    const name = makeSafeFilename(resolveOrthogroupName(group, shownLinearDrawing()) || id, id || 'orthogroup');
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

  /**
   * Copies or downloads the FASTA of a group's members, or of one member, read from Python when pressed, with the
   * read state of `context` (Reading, No sequence available, or the failure with Details and Retry).
   * @param {'copy' | 'download'} action
   * @param {GroupOrId} groupOrId
   * @param {unknown} sequenceKind
   * @param {OrthogroupMember | null} member
   * @param {MemberReadContext} context
   */
  const runOrthogroupSequenceAction = (action, groupOrId, sequenceKind, member, context) => {
    const shown = member ? enrichOrthogroupMember(member) : null;
    return sequenceActions.run({
      key: orthogroupCopyFeedbackKey(groupOrId, sequenceKind, member, context, action), context, action, sequenceKind,
      members: shown ? [shown] : getEnrichedOrthogroupMembers(groupOrId),
      filename: shown ? orthogroupMemberSequenceFilename(shown, sequenceKind, groupOrId) : orthogroupSequenceFilename(groupOrId, sequenceKind)
    });
  };
  /** @param {GroupOrId} [groupOrId] @param {unknown} [sequenceKind] @param {MemberReadContext} [context] */
  const copyOrthogroupSequences = (groupOrId = selectedOrthogroup.value, sequenceKind = 'nt', context = 'drawer') =>
    runOrthogroupSequenceAction('copy', groupOrId, sequenceKind, null, context);
  /** @param {GroupOrId} [groupOrId] @param {unknown} [sequenceKind] @param {MemberReadContext} [context] */
  const downloadOrthogroupSequences = (groupOrId = selectedOrthogroup.value, sequenceKind = 'nt', context = 'drawer') =>
    runOrthogroupSequenceAction('download', groupOrId, sequenceKind, null, context);
  /** @param {OrthogroupMember} member @param {unknown} [sequenceKind] @param {GroupOrId} [groupOrId] @param {MemberReadContext} [context] */
  const copyOrthogroupMemberSequence = (member, sequenceKind = 'nt', groupOrId = selectedOrthogroup.value, context = 'drawer') =>
    runOrthogroupSequenceAction('copy', groupOrId, sequenceKind, member, context);
  /** @param {OrthogroupMember} member @param {unknown} [sequenceKind] @param {GroupOrId} [groupOrId] @param {MemberReadContext} [context] */
  const downloadOrthogroupMemberSequence = (member, sequenceKind = 'nt', groupOrId = selectedOrthogroup.value, context = 'drawer') =>
    runOrthogroupSequenceAction('download', groupOrId, sequenceKind, member, context);

  /** @param {unknown} orthogroupId */
  const selectOrthogroup = (orthogroupId) => {
    const id = normalizeText(orthogroupId);
    if (!id || !getOrthogroupById(id)) return false;
    selectedOrthogroupId.value = id;
    return true;
  };

  /** @param {unknown} orthogroupId @param {unknown} value */
  const setOrthogroupNameOverride = (orthogroupId, value) => {
    const drawing = shownLinearDrawing();
    if (drawing === null) return;
    const sessionBusy = state.sessionOperationAvailability();
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
    const drawing = shownLinearDrawing();
    if (drawing === null) return;
    const sessionBusy = state.sessionOperationAvailability();
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
    const drawing = shownLinearDrawing();
    if (drawing === null) return;
    const sessionBusy = state.sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    const id = normalizeText(orthogroupId);
    if (!id) return;
    clearGroupRename(drawing.slices.groupNames, id);
  };

  // Dormant names are listed and cleared here (D-21).
  const orthogroupDormantNames = computed(() => Object.values(shownLinearDrawing()?.orthogroupDormantOverrides || {})
    .map((entry) => normalizeText(entry?.name) || normalizeText(entry?.description)).filter(Boolean));
  const clearOrthogroupDormantOverrides = () => {
    const drawing = shownLinearDrawing();
    if (drawing === null) return;
    const sessionBusy = state.sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    clearDormantGroupNames(drawing.slices.groupNames);
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
    const renderedIndex = buildRenderedFeatureIndex(svg, extractedFeatures.value);
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
    drawnFeatureIdsOf: anchorDrawnFeatureIds(featureCatalogLookup),
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
    orthogroupMemberRead,
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
