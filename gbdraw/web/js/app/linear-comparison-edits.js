// @ts-check
// The comparison plan edits of a Linear drawing: the one plan replacement
// (which clears the alignment candidate and the generated comparisons), the
// record row edits that keep the plan in step with the rows, the global
// action, the LOSAT program and mode, and the edits of a pair (add, omit,
// endpoints, source, BLAST TSV, LOSAT filename). The pure plan helpers are
// module functions; the composition root passes the shown-drawing readers and
// admission, History's undoable step, the alignment owner's clear, the record
// discovery, and the Comparison card's readers and focus paths as ports. The
// root creates these edits once the alignment owner exists.
import { SETTINGS_SLICE } from '../services/drawing/settings.js';
import {
  LINEAR_COMPARISON_MODES,
  LINEAR_COMPARISON_SOURCES,
  adjacentRowPairs,
  createLinearComparisonEdge,
  linearComparisonEdgeKey,
  materializeResolvedEdgesAsSelectedPlan,
  normalizeLinearComparisonPlan,
  reconcileLinearComparisonPlan,
  resolveLinearComparisonPlan,
  setComparisonEdgeFile,
  swapComparisonEdgeFile
} from '../services/linear-comparisons.js';
import {
  moveLinearRecordInRow,
  reconcileLinearRecordLayout,
  setLinearRecordRow as updateLinearRecordRow
} from '../services/linear-record-layout.js';
import { createDefaultLosatpHitLimits } from '../services/session-active-config-contract.js';
import { msg } from '../utils/i18n.js';
import { projectLinearComparisonLosatModeSelection, projectLinearComparisonLosatpModeSelection } from './comparison-ui.js';
import { effectiveLinearComparisonLayout } from './linear-comparison-panel.js';
/** @import { DiagramMode, DrawingState, SessionBusy, WebState } from '../state.js' */
/** @import { InputFileCandidate, InputFileRef } from '../services/history-files.js' */
/** @import { LinearComparisonEdge, LinearComparisonPair, LinearComparisonPlan, LinearComparisonResolvedEdge, LinearComparisonSource } from '../services/linear-comparisons.js' */
/** @import { LinearRecordRow } from '../services/linear-record-layout.js' */
/** @import { LinearComparisonUiIssue } from './comparison-ui.js' */

/**
 * A mode panel's mutating handler: the session's busy answer, nothing when the
 * shown drawing is of another mode, else `action` on the shown drawing.
 * @typedef {<T>(drawingMode: DiagramMode, action: (drawing: DrawingState) => T) => T | SessionBusy | undefined} ShownAction
 */
/** @typedef {(drawingMode: DiagramMode) => DrawingState | null} ShownDrawingOfMode The shown drawing when it is of that mode. */
/** @typedef {(label: string, fn: () => unknown) => Promise<unknown>} RunUndoable History's undoable step. */
/** @typedef {(reason: string) => unknown} ClearAlignmentForMutation The alignment owner clears its plan; `reason` is a `msg` clause. */
/**
 * The invalidation of a drawing's generated comparisons; `preserveLosatCacheInfo` keeps its LOSAT cache rows,
 * reindexed to the drawn records.
 * @typedef {(drawing: DrawingState, options?: { preserveLosatCacheInfo?: boolean }) => void} InvalidateLinearComparisonArtifacts
 */
/**
 * The replacement of a drawing's comparison plan with `nextPlan`, an unvalidated plan.
 * @typedef {(drawing: DrawingState, nextPlan: unknown, options?: { invalidate?: boolean }) => void} ReplaceLinearComparisonPlan
 */
/** @typedef {(options?: { preserveLosatCacheInfo?: boolean }) => boolean} SyncLinearRecordLayout */
/** @typedef {(edgeKey: string) => LinearComparisonPair | null} LinearComparisonPairForEdgeKey */

/**
 * @param {LinearComparisonEdge | undefined} left
 * @param {LinearComparisonEdge | undefined} right
 */
const sameLinearComparisonEdge = (left, right) => (
  left?.id === right?.id &&
  left?.queryUid === right?.queryUid &&
  left?.subjectUid === right?.subjectUid &&
  left?.included === right?.included &&
  left?.fileActive === right?.fileActive &&
  left?.losatFilenameActive === right?.losatFilenameActive &&
  left?.source === right?.source &&
  left?.file === right?.file &&
  left?.losatFilename === right?.losatFilename
);

// OV-380: a pair set to a source that has no data yet (Upload BLAST TSV
// without a file) compares nothing new, so it keeps the alignment. Every
// other pair must keep its current comparison; the pending pairs count as
// their current one.
/** @typedef {import('../services/linear-comparisons.js').LinearComparisonResolvedEdge} ResolvedEdge */
/** @param {DrawingState} drawing @param {import('../services/linear-comparisons.js').LinearComparisonPlan} normalized */
const onlyPendingSourcesChanged = (drawing, normalized) => {
  /** @param {ResolvedEdge} edge */
  const pending = (edge) => edge.source === LINEAR_COMPARISON_SOURCES.UPLOAD && !(edge.fileActive && edge.file);
  /** @type {readonly ResolvedEdge[]} */
  const before = drawing.linearComparisonResolution.value.edges;
  /** @type {readonly ResolvedEdge[]} */
  const after = resolveLinearComparisonPlan({
    plan: normalized, sequences: drawing.linearSeqs, layout: effectiveLinearComparisonLayout(drawing),
    losatProgram: drawing.losatProgram.value, blastpMode: drawing.losat.blastp?.mode
  }).edges;
  if (!after.some(pending)) return false;
  const beforeByKey = new Map(before.map((edge) => [edge.edgeKey, edge]));
  const effective = after.map((edge) => (pending(edge) ? beforeByKey.get(edge.edgeKey) : edge));
  return effective.length === before.length && effective.every((edge, index) => (
    edge?.edgeKey === before[index].edgeKey && edge.source === before[index].source
    && edge.file === before[index].file && edge.losatFilename === before[index].losatFilename
  ));
};

/** @param {DrawingState} drawing */
const selectedPlanForEdit = (drawing) => {
  if (drawing.linearComparisonPlan.mode === LINEAR_COMPARISON_MODES.SELECTED) {
    return normalizeLinearComparisonPlan(drawing.linearComparisonPlan);
  }
  return materializeResolvedEdgesAsSelectedPlan(
    drawing.linearComparisonPlan,
    drawing.linearComparisonResolution.value
  );
};

/**
 * @param {readonly LinearComparisonEdge[]} edges
 * @param {string} id
 */
const findEdgeIndex = (edges, id) => edges.findIndex((edge) => edge.id === id);
/**
 * @param {readonly LinearComparisonEdge[]} edges
 * @param {Pick<LinearComparisonPair, 'edgeKey' | 'edgeId'> & Partial<Pick<LinearComparisonPair, 'draft' | 'resolved'>>} pair
 */
const findEdgeIndexForPair = (edges, pair) => {
  const ownerId = pair?.draft?.id || pair?.resolved?.id || pair?.edgeId || '';
  const ownerIndex = ownerId ? findEdgeIndex(edges, ownerId) : -1;
  if (ownerIndex >= 0) return ownerIndex;
  const matchingIndexes = edges
    .map((edge, index) => (
      linearComparisonEdgeKey(edge.queryUid, edge.subjectUid) === pair?.edgeKey ? index : -1
    ))
    .filter((index) => index >= 0);
  return matchingIndexes.length === 1 ? matchingIndexes[0] : -1;
};

/**
 * @param {LinearComparisonPlan} next
 * @param {{ id?: string, queryUid: string, subjectUid: string, source?: LinearComparisonSource }} edge
 */
const upsertSelectedComparison = (next, {
  id = '',
  queryUid,
  subjectUid,
  source = next.defaultSource
}) => {
  const edgeKey = linearComparisonEdgeKey(queryUid, subjectUid);
  let index = id ? findEdgeIndex(next.edges, id) : -1;
  if (index < 0) {
    const matchingIndexes = next.edges
      .map((edge, edgeIndex) => (
        linearComparisonEdgeKey(edge.queryUid, edge.subjectUid) === edgeKey ? edgeIndex : -1
      ))
      .filter((edgeIndex) => edgeIndex >= 0);
    if (!id || matchingIndexes.length === 1) index = matchingIndexes[0] ?? -1;
  }
  if (index < 0) {
    next.edges.push(createLinearComparisonEdge({
      queryUid,
      subjectUid,
      included: true,
      source
    }));
    return next.edges[next.edges.length - 1];
  }
  const edge = next.edges[index];
  edge.queryUid = String(queryUid || '');
  edge.subjectUid = String(subjectUid || '');
  edge.included = true;
  edge.source = source;
  return edge;
};
/**
 * An omitted edge stays in the plan only to retain its file or LOSAT filename.
 * @param {LinearComparisonPlan} next
 * @param {number} index
 */
const omitPlanEdge = (next, index) => {
  if (index < 0) return;
  const edge = next.edges[index];
  edge.included = false;
  if (!edge.file && !String(edge.losatFilename || '').trim()) next.edges.splice(index, 1);
};

/**
 * @typedef {object} LinearComparisonPlanWriterPorts
 * @property {ShownAction} shownAction
 * @property {ShownDrawingOfMode} shownDrawingOfMode
 * @property {RunUndoable} runUndoable
 * @property {ClearAlignmentForMutation} clearAlignmentForMutation
 * @property {InvalidateLinearComparisonArtifacts} invalidateLinearComparisonArtifacts
 * @property {(drawing: DrawingState) => Promise<void>} materializeAutomaticLinearRecords
 *   Reads the records of the drawing's rows that wait for their discovery.
 */

/**
 * The plan replacement and the record row edits that keep the plan in step with the rows.
 * @param {LinearComparisonPlanWriterPorts} ports
 */
export const createLinearComparisonPlanWriter = ({
  shownAction, shownDrawingOfMode, runUndoable, clearAlignmentForMutation, invalidateLinearComparisonArtifacts,
  materializeAutomaticLinearRecords
}) => {
  /**
   * @param {DrawingState} drawing
   * @param {unknown} nextPlan An unvalidated plan; `normalizeLinearComparisonPlan` validates it.
   */
  const replaceLinearComparisonPlan = (drawing, nextPlan, { invalidate = true } = {}) => {
    const normalized = normalizeLinearComparisonPlan(nextPlan);
    const changed = normalized.mode !== drawing.linearComparisonPlan.mode
      || normalized.defaultSource !== drawing.linearComparisonPlan.defaultSource
      || normalized.edges.length !== drawing.linearComparisonPlan.edges.length
      || normalized.edges.some((edge, index) => (
        !sameLinearComparisonEdge(edge, drawing.linearComparisonPlan.edges[index])
      ));
    if (invalidate && changed && !onlyPendingSourcesChanged(drawing, normalized)) {
      clearAlignmentForMutation(msg('comparison configuration changed.'));
    }
    SETTINGS_SLICE.install(drawing.slices.settings, { linearComparisonPlan: normalized });
    if (invalidate) invalidateLinearComparisonArtifacts(drawing);
  };

  /**
   * @param {DrawingState} drawing
   * @param {(next: LinearComparisonPlan) => void} mutator
   */
  const mutateLinearComparisonPlan = (drawing, mutator) => runUndoable(msg('Change comparisons'), () => {
    const next = normalizeLinearComparisonPlan(drawing.linearComparisonPlan);
    mutator(next);
    replaceLinearComparisonPlan(drawing, next);
  });

  /** @param {DrawingState} drawing */
  const syncLinearComparisonRecords = (drawing, { invalidate = true } = {}) => {
    const next = reconcileLinearComparisonPlan(drawing.linearComparisonPlan, drawing.linearSeqs);
    const currentEdges = drawing.linearComparisonPlan.edges;
    const unchanged = (
      next.mode === drawing.linearComparisonPlan.mode &&
      next.defaultSource === drawing.linearComparisonPlan.defaultSource &&
      next.edges.length === currentEdges.length &&
      next.edges.every((edge, index) => sameLinearComparisonEdge(edge, currentEdges[index]))
    );
    if (!unchanged) replaceLinearComparisonPlan(drawing, next, { invalidate });
    return !unchanged;
  };

  const syncLinearRecordLayout = ({ preserveLosatCacheInfo = false } = {}) => {
    const drawing = shownDrawingOfMode('linear');
    if (!drawing) return false;
    const next = reconcileLinearRecordLayout(drawing.linearSeqs, drawing.linearRecordRows);
    const rowsUnchanged = next.length === drawing.linearRecordRows.length && next.every((entry, index) => (
      entry.uid === drawing.linearRecordRows[index]?.uid && entry.row === drawing.linearRecordRows[index]?.row
    ));
    if (!rowsUnchanged) drawing.linearRecordRows.splice(0, drawing.linearRecordRows.length, ...next);
    const comparisonsChanged = syncLinearComparisonRecords(drawing, { invalidate: false });
    if (!rowsUnchanged || comparisonsChanged) {
      invalidateLinearComparisonArtifacts(drawing, { preserveLosatCacheInfo });
    }
    return !rowsUnchanged || comparisonsChanged;
  };
  /** @param {string} uid @param {unknown} row */
  const setLinearRecordRow = (uid, row) => shownAction('linear', (drawing) => {
    syncLinearRecordLayout({ preserveLosatCacheInfo: true });
    const previous = drawing.linearRecordRows.find((/** @type {LinearRecordRow} */ entry) => entry.uid === uid)?.row;
    updateLinearRecordRow(drawing.linearRecordRows, uid, row);
    if (drawing.linearRecordRows.find((/** @type {LinearRecordRow} */ entry) => entry.uid === uid)?.row !== previous) {
      invalidateLinearComparisonArtifacts(drawing, { preserveLosatCacheInfo: true });
    }
  });
  /** @param {boolean} enabled */
  const setLinearRecordLayoutEnabled = (enabled) => shownAction('linear', (drawing) => {
    const nextEnabled = Boolean(enabled);
    if (drawing.linearRecordLayoutEnabled.value === nextEnabled) return;
    drawing.linearRecordLayoutEnabled.value = nextEnabled;
    syncLinearRecordLayout({ preserveLosatCacheInfo: true });
    invalidateLinearComparisonArtifacts(drawing, { preserveLosatCacheInfo: true });
    if (nextEnabled) return materializeAutomaticLinearRecords(drawing);
  });
  /** @param {string} uid @param {number} direction */
  const moveLinearRecordWithinRow = (uid, direction) => shownAction('linear', (drawing) => {
    const next = moveLinearRecordInRow(drawing.linearSeqs, drawing.linearRecordRows, uid, direction);
    drawing.linearRecordRows.splice(0, drawing.linearRecordRows.length, ...next);
    syncLinearComparisonRecords(drawing, { invalidate: false });
    invalidateLinearComparisonArtifacts(drawing, { preserveLosatCacheInfo: true });
  });
  return {
    replaceLinearComparisonPlan, mutateLinearComparisonPlan, syncLinearRecordLayout, setLinearRecordRow,
    setLinearRecordLayoutEnabled, moveLinearRecordWithinRow
  };
};

/**
 * @typedef {object} LinearComparisonModesPorts
 * @property {ShownAction} shownAction
 * @property {RunUndoable} runUndoable
 * @property {(drawing: DrawingState, mutator: (next: LinearComparisonPlan) => void) => Promise<unknown>} mutateLinearComparisonPlan
 *   A plan edit as one History step.
 * @property {(disclosureKey: string) => Promise<HTMLDetailsElement | null>} openLinearComparisonDisclosure
 * @property {(edgeKey: string, options?: { focus?: boolean }) => boolean} chooseLinearComparisonPairFile
 * @property {(target?: LinearComparisonUiIssue) => Promise<boolean>} focusLinearComparisonIssue
 * @property {ClearAlignmentForMutation} clearAlignmentForMutation
 * @property {InvalidateLinearComparisonArtifacts} invalidateLinearComparisonArtifacts
 */

/**
 * The global action, the LOSAT program and mode of the shown drawing, and the actions of the issue list.
 * @param {LinearComparisonModesPorts} ports
 */
export const createLinearComparisonModes = ({
  shownAction, runUndoable, mutateLinearComparisonPlan, openLinearComparisonDisclosure, chooseLinearComparisonPairFile,
  focusLinearComparisonIssue, clearAlignmentForMutation, invalidateLinearComparisonArtifacts
}) => {
  // An action of an entry in the Comparison card's issue list; a
  // setting it changes is one Undo step.
  /** @param {LinearComparisonUiIssue} issue @param {string} actionKey */
  const runLinearComparisonIssueAction = (issue, actionKey) => {
    if (actionKey === 'all-losat') return setLinearComparisonGlobalAction('losat');
    if (actionKey === 'pairwise') {
      return runUndoable(msg('Change setting'), () => setLinearComparisonLosatpMode('pairwise'));
    }
    if (actionKey === 'choose-file' && chooseLinearComparisonPairFile(issue.edgeKey, { focus: true })) return true;
    return focusLinearComparisonIssue(issue);
  };

  /** @param {unknown} action */
  const setLinearComparisonGlobalAction = async (action) => shownAction('linear', async (drawing) => {
    const normalized = String(action || '').trim().toLowerCase();
    const result = await mutateLinearComparisonPlan(drawing, (next) => {
      if (normalized === 'none') {
        next.mode = LINEAR_COMPARISON_MODES.NONE;
        return;
      }
      next.mode = LINEAR_COMPARISON_MODES.ADJACENT;
      next.defaultSource = normalized === LINEAR_COMPARISON_SOURCES.UPLOAD
        ? LINEAR_COMPARISON_SOURCES.UPLOAD
        : LINEAR_COMPARISON_SOURCES.LOSAT;
    });
    if (normalized === LINEAR_COMPARISON_SOURCES.LOSAT) {
      await openLinearComparisonDisclosure('settings');
    }
    return result;
  });

  /** @param {string} modeKey */
  const setLinearComparisonLosatMode = (modeKey) => shownAction('linear', (drawing) => {
    const selection = projectLinearComparisonLosatModeSelection({ modeKey });
    if (!selection.selectable || !selection.patch) return false;
    const nextProgram = selection.patch.losatProgram;
    if (drawing.losatProgram.value === nextProgram) return true;
    clearAlignmentForMutation(msg('comparison program changed.'));
    drawing.losatProgram.value = nextProgram;
    invalidateLinearComparisonArtifacts(drawing);
    return true;
  });

  /** @param {string} modeKey */
  const setLinearComparisonLosatpMode = (modeKey) => shownAction('linear', (drawing) => {
    const selection = projectLinearComparisonLosatpModeSelection({
      plan: drawing.linearComparisonPlan,
      modeKey
    });
    if (!selection.selectable || !selection.patch) return false;
    const nextBlastpMode = selection.patch.blastpMode;
    if (drawing.losat.blastp?.mode === nextBlastpMode) return true;
    clearAlignmentForMutation(msg('comparison mode changed.'));
    const hitLimits = drawing.losat.blastp.hitLimitsByMode ||= createDefaultLosatpHitLimits();
    hitLimits[drawing.losat.blastp.mode] = {
      candidateLimit: drawing.losat.blastp.candidateLimit,
      orthogroupMemberMaxHits: drawing.losat.blastp.orthogroupMemberMaxHits
    };
    Object.assign(drawing.losat.blastp, hitLimits[nextBlastpMode]);
    drawing.losat.blastp.mode = nextBlastpMode;
    invalidateLinearComparisonArtifacts(drawing);
    return true;
  });
  return {
    runLinearComparisonIssueAction, setLinearComparisonGlobalAction, setLinearComparisonLosatMode,
    setLinearComparisonLosatpMode
  };
};

/**
 * @typedef {object} LinearComparisonPairEditsPorts
 * @property {WebState} state
 * @property {ShownAction} shownAction
 * @property {ReplaceLinearComparisonPlan} replaceLinearComparisonPlan
 * @property {SyncLinearRecordLayout} syncLinearRecordLayout
 * @property {(edgeKey: string) => Promise<void>} focusLinearComparisonPair Brings a pair's row into view.
 * @property {LinearComparisonPairForEdgeKey} linearComparisonPairForEdgeKey
 */

/**
 * The selected pairs of the shown drawing: add, omit, clear, a pair's endpoints and source, and its BLAST TSV.
 * @param {LinearComparisonPairEditsPorts} ports
 */
export const createLinearComparisonPairEdits = ({
  state, shownAction, replaceLinearComparisonPlan, syncLinearRecordLayout, focusLinearComparisonPair,
  linearComparisonPairForEdgeKey
}) => {
  const addLinearComparison = async () => shownAction('linear', async (drawing) => {
    const drawn = drawing.drawnLinearSeqs.value;
    if (drawn.length < 2) return;
    syncLinearRecordLayout();
    const [firstPair] = adjacentRowPairs(drawn, effectiveLinearComparisonLayout(drawing));
    const [queryUid, subjectUid] = firstPair || [drawn[0].uid, drawn[1].uid];
    const next = selectedPlanForEdit(drawing);
    upsertSelectedComparison(next, { queryUid, subjectUid });
    replaceLinearComparisonPlan(drawing, next);
    await focusLinearComparisonPair(linearComparisonEdgeKey(queryUid, subjectUid));
  });
  /**
   * Sets the pairs to No comparison as one plan replacement: a pair's radio,
   * and the pairs of the BLAST TSV missing dialog (D-54).
   * @param {DrawingState} drawing
   * @param {string[]} edgeKeys
   */
  const omitLinearComparisonPairs = (drawing, edgeKeys) => {
    const next = selectedPlanForEdit(drawing);
    edgeKeys.forEach((edgeKey) => {
      const pair = linearComparisonPairForEdgeKey(edgeKey);
      if (pair) omitPlanEdge(next, findEdgeIndexForPair(next.edges, pair));
    });
    replaceLinearComparisonPlan(drawing, next);
  };
  /** @param {string} id */
  const omitLinearComparison = (id) => shownAction('linear', (drawing) => {
    const next = selectedPlanForEdit(drawing);
    const index = findEdgeIndex(next.edges, id);
    if (index < 0) return;
    omitPlanEdge(next, index);
    replaceLinearComparisonPlan(drawing, next);
  });
  const clearSelectedLinearComparisons = () => shownAction('linear', (drawing) => {
    const next = selectedPlanForEdit(drawing);
    next.edges = next.edges
      .filter((edge) => edge.file || String(edge.losatFilename || '').trim())
      .map((edge) => ({ ...edge, included: false }));
    replaceLinearComparisonPlan(drawing, next);
  });
  /** @param {string} id @param {'queryUid' | 'subjectUid'} endpoint @param {unknown} uid */
  const setLinearComparisonEndpoint = (id, endpoint, uid) => shownAction('linear', (drawing) => {
    if (!['queryUid', 'subjectUid'].includes(endpoint)) return;
    const next = selectedPlanForEdit(drawing);
    const edge = next.edges.find((entry) => entry.id === id);
    if (!edge) return;
    edge[endpoint] = String(uid || '');
    edge.included = true;
    replaceLinearComparisonPlan(drawing, next);
  });
  /** @param {string} id @param {unknown} source */
  const setLinearComparisonSource = (id, source) => shownAction('linear', (drawing) => {
    const normalized = source === LINEAR_COMPARISON_SOURCES.LOSAT
      ? LINEAR_COMPARISON_SOURCES.LOSAT
      : LINEAR_COMPARISON_SOURCES.UPLOAD;
    const next = selectedPlanForEdit(drawing);
    const edge = next.edges.find((entry) => entry.id === id);
    if (!edge) return;
    edge.source = normalized;
    edge.included = true;
    replaceLinearComparisonPlan(drawing, next);
  });
  /**
   * The shown drawing's upload into a comparison's file input, found by the
   * comparison's id: the file becomes the pair's active upload.
   * @param {string} id @param {InputFileCandidate} file
   */
  const setLinearComparisonFile = (id, file) => shownAction('linear', (drawing) => {
    const next = selectedPlanForEdit(drawing);
    const edge = next.edges.find((entry) => entry.id === id);
    if (!edge) return;
    const ref = state.inputFiles.registerFile(file);
    setComparisonEdgeFile(edge, ref);
    edge.source = LINEAR_COMPARISON_SOURCES.UPLOAD;
    edge.included = Boolean(ref) || edge.included;
    replaceLinearComparisonPlan(drawing, next);
  });
  /**
   * A replacement of the file a comparison binds, in any drawing (PD-OI-090
   * 〔入力ファイル〕): only the file changes (R17-28).
   * @param {DrawingState} drawing @param {string} id @param {InputFileCandidate | InputFileRef} file
   */
  const swapLinearComparisonFile = (drawing, id, file) => {
    const next = selectedPlanForEdit(drawing);
    const edge = next.edges.find((entry) => entry.id === id);
    if (!edge) return;
    const ref = state.inputFiles.registerFile(file);
    if (!ref) return;
    swapComparisonEdgeFile(edge, ref);
    replaceLinearComparisonPlan(drawing, next);
  };
  /** @param {string} id */
  const reuseLinearComparisonFile = (id) => shownAction('linear', (drawing) => {
    const next = selectedPlanForEdit(drawing);
    const edge = next.edges.find((entry) => entry.id === id);
    if (!edge?.file) return;
    edge.fileActive = true;
    edge.source = LINEAR_COMPARISON_SOURCES.UPLOAD;
    edge.included = true;
    replaceLinearComparisonPlan(drawing, next);
  });
  /** @param {string} id */
  const deactivateLinearComparisonFile = (id) => shownAction('linear', (drawing) => {
    const next = selectedPlanForEdit(drawing);
    const edge = next.edges.find((entry) => entry.id === id);
    if (!edge?.file) return;
    edge.fileActive = false;
    if (edge.source === LINEAR_COMPARISON_SOURCES.UPLOAD) edge.included = false;
    replaceLinearComparisonPlan(drawing, next);
  });
  return {
    addLinearComparison, omitLinearComparisonPairs, omitLinearComparison, clearSelectedLinearComparisons,
    setLinearComparisonEndpoint, setLinearComparisonSource, setLinearComparisonFile, swapLinearComparisonFile,
    reuseLinearComparisonFile, deactivateLinearComparisonFile
  };
};

/**
 * @typedef {object} LinearComparisonGapEditsPorts
 * @property {() => SessionBusy | null} sessionOperationAvailability
 * @property {ShownAction} shownAction
 * @property {ReplaceLinearComparisonPlan} replaceLinearComparisonPlan
 * @property {SyncLinearRecordLayout} syncLinearRecordLayout
 * @property {LinearComparisonPairForEdgeKey} linearComparisonPairForEdgeKey
 * @property {(drawing: DrawingState, edgeKeys: string[]) => void} omitLinearComparisonPairs
 * @property {(id: string, file: InputFileCandidate) => unknown} setLinearComparisonFile
 */

/**
 * The LOSAT filename of a pair (by its id, or by the edge key of a resolved pair), Add all pairs, the action of a
 * pair between two rows, and the BLAST TSV chosen on such a pair.
 * @param {LinearComparisonGapEditsPorts} ports
 */
export const createLinearComparisonGapEdits = ({
  sessionOperationAvailability, shownAction, replaceLinearComparisonPlan, syncLinearRecordLayout,
  linearComparisonPairForEdgeKey, omitLinearComparisonPairs, setLinearComparisonFile
}) => {
  /** @param {string} id @param {unknown} value */
  const setLinearComparisonLosatFilename = (id, value) => shownAction('linear', (drawing) => {
    const next = selectedPlanForEdit(drawing);
    const edge = next.edges.find((entry) => entry.id === id);
    if (!edge) return;
    edge.losatFilename = String(value || '');
    edge.losatFilenameActive = Boolean(edge.losatFilename.trim());
    replaceLinearComparisonPlan(drawing, next, { invalidate: false });
  });
  /** @param {string} id */
  const reuseLinearComparisonLosatFilename = (id) => shownAction('linear', (drawing) => {
    const next = selectedPlanForEdit(drawing);
    const edge = next.edges.find((entry) => entry.id === id);
    if (!edge || !String(edge.losatFilename || '').trim()) return;
    edge.losatFilenameActive = true;
    edge.source = LINEAR_COMPARISON_SOURCES.LOSAT;
    edge.included = true;
    replaceLinearComparisonPlan(drawing, next);
  });
  /** @param {string} id */
  const deactivateLinearComparisonLosatFilename = (id) => shownAction('linear', (drawing) => {
    const next = selectedPlanForEdit(drawing);
    const edge = next.edges.find((entry) => entry.id === id);
    if (!edge) return;
    edge.losatFilenameActive = false;
    replaceLinearComparisonPlan(drawing, next, { invalidate: false });
  });
  /**
   * @param {DrawingState} drawing
   * @param {string} edgeKey
   * @param {(edge: LinearComparisonEdge) => void} updater
   */
  const updateResolvedLosatFilenameDraft = (drawing, edgeKey, updater) => {
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    const resolved = drawing.linearComparisonResolution.value.edges.find((/** @type {LinearComparisonResolvedEdge} */ edge) => edge.edgeKey === edgeKey);
    if (!resolved) return;
    const next = normalizeLinearComparisonPlan(drawing.linearComparisonPlan);
    const pair = linearComparisonPairForEdgeKey(edgeKey) || {
      edgeKey,
      edgeId: resolved.id,
      resolved
    };
    let index = findEdgeIndexForPair(next.edges, pair);
    if (index < 0) {
      next.edges.push(createLinearComparisonEdge({
        queryUid: resolved.queryUid,
        subjectUid: resolved.subjectUid,
        included: false,
        source: /** @type {LinearComparisonSource} */ (LINEAR_COMPARISON_SOURCES.LOSAT)
      }));
      index = next.edges.length - 1;
    }
    updater(next.edges[index]);
    replaceLinearComparisonPlan(drawing, next, { invalidate: false });
  };
  /** @param {string} edgeKey @param {unknown} value */
  const setResolvedLinearComparisonLosatFilename = (edgeKey, value) => shownAction('linear', (drawing) => {
    updateResolvedLosatFilenameDraft(drawing, edgeKey, (edge) => {
      edge.losatFilename = String(value || '');
      edge.losatFilenameActive = Boolean(edge.losatFilename.trim());
    });
  });
  /** @param {string} edgeKey */
  const reuseResolvedLinearComparisonLosatFilename = (edgeKey) => shownAction('linear', (drawing) => {
    updateResolvedLosatFilenameDraft(drawing, edgeKey, (edge) => {
      if (String(edge.losatFilename || '').trim()) edge.losatFilenameActive = true;
    });
  });
  /** @param {string} edgeKey */
  const deactivateResolvedLinearComparisonLosatFilename = (edgeKey) => shownAction('linear', (drawing) => {
    updateResolvedLosatFilenameDraft(drawing, edgeKey, (edge) => {
      edge.losatFilenameActive = false;
    });
  });
  const addLinearComparisonBatch = (allPairs = false) => shownAction('linear', (drawing) => {
    syncLinearRecordLayout();
    const next = selectedPlanForEdit(drawing);
    adjacentRowPairs(drawing.drawnLinearSeqs.value, effectiveLinearComparisonLayout(drawing), allPairs).forEach(([queryUid, subjectUid]) => {
      upsertSelectedComparison(next, { queryUid, subjectUid });
    });
    replaceLinearComparisonPlan(drawing, next);
  });
  /** @param {string} edgeKey @param {string} action */
  const setLinearComparisonGapAction = (edgeKey, action) => shownAction('linear', (drawing) => {
    const pair = linearComparisonPairForEdgeKey(edgeKey);
    if (!pair) return;
    if (action === 'none') {
      omitLinearComparisonPairs(drawing, [edgeKey]);
      return;
    }
    const next = selectedPlanForEdit(drawing);
    upsertSelectedComparison(next, {
      id: pair.edgeId,
      queryUid: pair.queryUid,
      subjectUid: pair.subjectUid,
      source: action === LINEAR_COMPARISON_SOURCES.UPLOAD
        ? LINEAR_COMPARISON_SOURCES.UPLOAD
        : LINEAR_COMPARISON_SOURCES.LOSAT
    });
    replaceLinearComparisonPlan(drawing, next);
  });
  /** @param {string} edgeKey @param {File | null} file */
  const setLinearComparisonCardFile = (edgeKey, file) => {
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    let pair = linearComparisonPairForEdgeKey(edgeKey);
    let draft = pair?.draft || null;
    if (!draft) {
      setLinearComparisonGapAction(edgeKey, LINEAR_COMPARISON_SOURCES.UPLOAD);
      pair = linearComparisonPairForEdgeKey(edgeKey);
      draft = pair?.draft || null;
    }
    if (draft) setLinearComparisonFile(draft.id, file);
  };
  return {
    setLinearComparisonLosatFilename, reuseLinearComparisonLosatFilename, deactivateLinearComparisonLosatFilename,
    setResolvedLinearComparisonLosatFilename, reuseResolvedLinearComparisonLosatFilename,
    deactivateResolvedLinearComparisonLosatFilename, addLinearComparisonBatch, setLinearComparisonGapAction,
    setLinearComparisonCardFile
  };
};
