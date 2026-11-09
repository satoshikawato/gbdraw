// @ts-check
// The Similarity group actions of the clicked feature's popup in a Linear drawing: whether they
// apply, the group's detail the popup shows, Align by it, highlight it, open it in the Editor's
// Similarity groups tab, and choose the alignment's reference there again. The orthogroup owner
// (app/orthogroups.js) and the alignment owner (app/similarity-alignment.js) keep the groups and
// the plan; the composition root passes their members, the Editor drawer's tab, and the review's
// focus helpers (app/similarity-alignment-review.js) as ports.
/** @import { DrawingState, WebState } from '../state.js' */
/** @import { GroupOrId } from './orthogroups.js' */
/** @import { OrthogroupGroup, OrthogroupMember } from '../services/orthogroup-feature-metadata.js' */
/** @import { AlignmentStartOptions } from './similarity-alignment.js' */
import { resolveUniqueOrthogroupMemberForFeature } from './orthogroups.js';
import { groupCountText } from '../services/pairwise-match-popup.js';

const { computed } = window.Vue;

/**
 * @typedef {object} OrthogroupPopupPorts
 * @property {WebState} state
 * @property {() => DrawingState | null} displayedLinearDrawing The drawing of the displayed Result when it is Linear.
 * @property {(orthogroupId: unknown) => OrthogroupGroup | null} getOrthogroupById
 * @property {(groupOrId: GroupOrId) => OrthogroupMember[]} getEnrichedOrthogroupMembers
 * @property {(members: unknown) => { recordIndex: number, recordLabel: string, members: OrthogroupMember[] }[]} groupOrthogroupMembersByRecord
 * @property {(groupOrId: GroupOrId, drawing: DrawingState | null) => string} resolveOrthogroupName
 * @property {(groupOrId: GroupOrId) => number} getOrthogroupSequenceCount The members whose sequences can be read.
 * @property {(groupOrId: GroupOrId) => string} orthogroupScope
 * @property {(groupOrId: GroupOrId) => string} orthogroupScopeLabel
 * @property {(orthogroupId: string) => void} highlightOrthogroupById
 * @property {() => void} removeOrthogroupHighlight The orthogroup owner's clear of the highlight.
 * @property {(orthogroupId: unknown) => boolean} selectOrthogroup
 * @property {(options: AlignmentStartOptions) => Promise<unknown>} startAlignmentFromPopup
 * @property {(groupId: unknown, key: string) => boolean} setAlignmentDrawerReference
 * @property {{ readonly value: { groupId?: string } | null }} alignmentRepair The plan the alignment owner asks to repair.
 * @property {(tab: string) => unknown} openRightDrawerTab
 * @property {(event: Event | null | undefined) => void} rememberSimilarityAlignmentInvoker
 * @property {<T>(outcome: T) => Promise<T>} finishSimilarityAlignmentStart
 */

/**
 * Registers the popup's two computeds; the root calls it where they stood.
 * @param {OrthogroupPopupPorts} ports
 */
export const createOrthogroupPopup = ({
  state, displayedLinearDrawing, getOrthogroupById, getEnrichedOrthogroupMembers, groupOrthogroupMembersByRecord,
  resolveOrthogroupName, getOrthogroupSequenceCount, orthogroupScope, orthogroupScopeLabel, highlightOrthogroupById,
  removeOrthogroupHighlight, selectOrthogroup, startAlignmentFromPopup, setAlignmentDrawerReference, alignmentRepair,
  openRightDrawerTab, rememberSimilarityAlignmentInvoker, finishSimilarityAlignmentStart
}) => {
  const { clickedFeature, mode } = state;
  const canUseClickedOrthogroupActions = computed(() => {
    const drawing = state.activeDrawing();
    const cf = clickedFeature.value;
    return Boolean(
      cf &&
      mode.value === 'linear' &&
      drawing.hasActiveLinearLosatIntent.value &&
      drawing.linearComparisonResolution.value.valid &&
      drawing.losatProgram.value === 'blastp' &&
      drawing.losat.blastp?.mode === 'orthogroup' &&
      drawing.inputType.value === 'gb' &&
      cf.feat.type === 'CDS' &&
      cf.orthogroupId
    );
  });

  const clickedOrthogroupDetail = computed(() => {
    const cf = clickedFeature.value;
    const orthogroupId = String(cf?.orthogroupId || '').trim();
    if (!orthogroupId) return null;
    const group = getOrthogroupById(orthogroupId);
    if (!group) return null;
    const members = getEnrichedOrthogroupMembers(group);
    const currentMember = resolveUniqueOrthogroupMemberForFeature(cf?.feat, members);
    const membersByRecord = groupOrthogroupMembersByRecord(members);
    const memberCount = Number(group.member_count || members.length || 0);
    const recordCoverage = Number(group.record_coverage_count || membersByRecord.length || 0);
    return {
      id: orthogroupId,
      displayName: resolveOrthogroupName(group, displayedLinearDrawing()),
      memberCount,
      recordCoverage,
      countText: groupCountText(memberCount, recordCoverage),
      genesText: groupCountText(memberCount, 0),
      sequenceCount: getOrthogroupSequenceCount(group),
      currentMember,
      membersByRecord,
      scopeLabel: orthogroupScope(group) === 'global'
        ? ''
        : orthogroupScopeLabel(group)
    };
  });

  const alignByClickedOrthogroup = async (event = null, mode = 'align') => {
    const detail = clickedOrthogroupDetail.value;
    if (!detail?.id) return { status: 'rejected' };
    rememberSimilarityAlignmentInvoker(event);
    const outcome = await startAlignmentFromPopup({
      groupId: detail.id,
      reference: detail.currentMember,
      mode
    });
    return finishSimilarityAlignmentStart(outcome);
  };

  const highlightClickedOrthogroup = () => {
    const cf = clickedFeature.value;
    const orthogroupId = String(cf?.orthogroupId || '').trim();
    if (!orthogroupId) return;
    highlightOrthogroupById(orthogroupId);
  };

  const clearOrthogroupHighlight = () => {
    removeOrthogroupHighlight();
  };

  /** @param {string} orthogroupId */
  const openOrthogroupInDrawer = (orthogroupId) => {
    if (!selectOrthogroup(orthogroupId)) return false;
    openRightDrawerTab('orthogroups');
    return true;
  };

  const reselectSimilarityAlignmentReference = () => {
    const groupId = alignmentRepair.value?.groupId
      || state.similarityAlignmentPlan.value?.groupId
      || '';
    setAlignmentDrawerReference(groupId, '');
    return openOrthogroupInDrawer(groupId);
  };

  const openClickedOrthogroupInEditor = () => {
    const orthogroupId = String(clickedFeature.value?.orthogroupId || '').trim();
    if (!openOrthogroupInDrawer(orthogroupId)) return;
    clickedFeature.value = null;
  };
  return {
    canUseClickedOrthogroupActions, clickedOrthogroupDetail, alignByClickedOrthogroup, highlightClickedOrthogroup,
    clearOrthogroupHighlight, openOrthogroupInDrawer, reselectSimilarityAlignmentReference, openClickedOrthogroupInEditor
  };
};
