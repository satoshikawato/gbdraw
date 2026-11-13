// @ts-check
// The Comparison card of the shown Linear drawing: the comparison height
// check, the card's projection (its intent, pairs, issues), the pair timeline
// and the record names its rows show, the LOSAT cache rows of the pairs, and
// the paths that open a disclosure of the card and bring a pair, an issue, or
// a pair's BLAST TSV uploader into view. The File card row input and the
// record row tokens are here too. The composition root passes the state, the
// shown-drawing readers, the config service's LOSAT cache-info write, and the
// settings pane's reveal as ports. The comparison plan's edits are in
// app/linear-comparison-edits.js.
import { bindDrawingCanonicalComparisons } from '../services/drawing/inputs.js';
import {
  adjacentRowPairs,
  buildLinearComparisonTimeline,
  linearComparisonEdgeKey
} from '../services/linear-comparisons.js';
import { linearRecordPositionTokens } from '../services/linear-record-layout.js';
import { linearRecordNames } from '../services/linear-sources.js';
import { catalogUnitPlaces, isDrawnWindowPairRow, linearWindowPairsByCardEdge } from '../services/linear-units.js';
import { inputSourceFiles } from '../services/input-types.js';
import {
  comparisonCauseRepeats, diagnosticError, losatFailureRepeats, normalizeUserFacingError
} from '../utils/error-normalization.js';
import { errorText } from '../utils/error-text.js';
import { t } from '../utils/i18n.js';
import { classifyOptionalPositiveNumber } from '../utils/optional-positive-number.js';
import { projectLinearComparisonUi } from './comparison-ui.js';
/** @import { DiagramMode, DrawingState, LinearSeq, WebState } from '../state.js' */
/** @import { LinearComparisonPair, LinearComparisonResolvedEdge } from '../services/linear-comparisons.js' */
/** @import { LinearRecordRow } from '../services/linear-record-layout.js' */
/** @import { LinearWindowPair } from '../services/linear-units.js' */
/** @import { EditLosatRows, LosatCacheInfoEntry } from '../services/losat-cache.js' */
/** @import { UserFacingError } from '../utils/error-normalization.js' */
/** @import { RevealSettingsSection } from './settings-pane.js' */

const { computed, nextTick, watch } = window.Vue;

/** @typedef {(drawingMode: DiagramMode) => DrawingState | null} ShownDrawingOfMode The shown drawing when it is of that mode. */
/** @typedef {{ readonly value: unknown }} LosatCacheInfoRows The shown drawing's LOSAT result rows. */
/**
 * The control of the Comparison card that an issue names, and the pair it is about.
 * @typedef {{ disclosureKey: string, focusTargetKey: string, edgeKey?: string, edgeId?: string }} LinearComparisonFocusTarget
 */
/** @typedef {(target: HTMLElement, options?: { align?: 'start', focus?: HTMLElement }) => Promise<void>} RevealSetting */

/** @param {DrawingState} drawing */
export const effectiveLinearComparisonLayout = (drawing) => (
  drawing.linearRecordLayoutEnabled.value ? drawing.linearRecordRows : []
);
/** @param {DrawingState} drawing @param {string} uid @param {number} fallback */
export const linearRecordRowFor = (drawing, uid, fallback) => (
  drawing.linearRecordRows.find((/** @type {LinearRecordRow} */ entry) => entry.uid === uid)?.row || fallback
);

/**
 * @typedef {object} LinearComparisonViewPorts
 * @property {WebState} state
 * @property {{ readonly value: DiagramMode }} mode The shown drawing's mode.
 * @property {EditLosatRows} editLosatRows
 * @property {ShownDrawingOfMode} shownDrawingOfMode
 */

/**
 * The Comparison card of the shown drawing: the comparison height check, the card's projection and its pressed
 * global action, whether LOSAT can run, and the invalidation of a drawing's generated comparisons and LOSAT cache rows
 * after its plan or record rows changed.
 * @param {LinearComparisonViewPorts} ports
 */
export const createLinearComparisonView = ({ state, mode, editLosatRows, shownDrawingOfMode }) => {
  const comparisonHeightValidationError = computed(() => {
    const drawing = state.activeDrawing();
    if (
      mode.value !== 'linear' ||
      drawing.linearComparisonResolution.value.hasComparisonIntent !== true
    ) return '';
    return classifyOptionalPositiveNumber(drawing.adv.comparison_height).status === 'invalid'
      // A diagnostic error object is truthy, so it always normalizes to a model.
      ? errorText(/** @type {UserFacingError} */ (normalizeUserFacingError(diagnosticError('INPUT_INVALID', { field: 'match_height', reason: 'POSITIVE_OR_AUTO' }))).summaryParts)
      : '';
  });

  /** @param {DrawingState} drawing */
  const reindexLinearLosatCacheInfo = (drawing) => {
    // Resolution indexes are positions among the drawn records.
    const indexByUid = new Map(
      drawing.drawnLinearSeqs.value.map((sequence, index) => [String(sequence?.uid || ''), index])
    );
    const resolvedByEdgeKey = new Map(
      drawing.linearComparisonResolution.value.edges.map((/** @type {LinearComparisonResolvedEdge} */ edge) => [edge.edgeKey, edge])
    );
    // A window pair's row names the drawn units of the last Generate, whose indexes only that request counts: it
    // stays while the records of its units are drawn (the drawn Result's units, OV-2742).
    const places = catalogUnitPlaces(state.featureCatalog.value);
    editLosatRows(drawing.id, (rows) => rows.flatMap((entry) => {
      const edgeKey = String(entry?.edgeKey || '');
      if (!edgeKey) return [entry];
      if (isDrawnWindowPairRow(entry, places, (uid) => indexByUid.has(uid))) return [entry];
      const resolved = resolvedByEdgeKey.get(edgeKey);
      const queryUid = String(resolved?.queryUid || entry?.queryUid || '');
      const subjectUid = String(resolved?.subjectUid || entry?.subjectUid || '');
      const queryIndex = indexByUid.get(queryUid);
      const subjectIndex = indexByUid.get(subjectUid);
      if (!Number.isInteger(queryIndex) || !Number.isInteger(subjectIndex)) return [];
      return [{
        ...entry,
        edgeKey,
        queryUid,
        subjectUid,
        queryIndex,
        subjectIndex,
        ordinal: Number.isInteger(Number(resolved?.ordinal))
          ? Number(resolved?.ordinal)
          : entry.ordinal
      }];
    }));
  };

  /** @param {DrawingState} drawing */
  const invalidateLinearComparisonArtifacts = (drawing, { preserveLosatCacheInfo = false } = {}) => {
    bindDrawingCanonicalComparisons(state.inputFiles, drawing.slices.inputs, []);
    if (preserveLosatCacheInfo) reindexLinearLosatCacheInfo(drawing);
    else editLosatRows(drawing.id, (rows) => rows.filter((entry) => !entry?.edgeKey));
  };

  const linearComparisonUi = computed(() => {
    const drawing = state.activeDrawing();
    return projectLinearComparisonUi({
      plan: drawing.linearComparisonPlan,
      resolution: drawing.linearComparisonResolution.value,
      adjacentEdgeKeys: adjacentRowPairs(drawing.drawnLinearSeqs.value, effectiveLinearComparisonLayout(drawing), true)
        .map(([queryUid, subjectUid]) => linearComparisonEdgeKey(queryUid, subjectUid)),
      losatProgram: drawing.losatProgram.value,
      blastpMode: drawing.losat.blastp.mode,
      inferOrthogroups: drawing.losat.blastp.collinearInferOrthogroups === true,
      anchorMode: drawing.losat.blastp.collinearAnchorMode,
      filters: {
        min_bitscore: drawing.adv.min_bitscore,
        evalue: drawing.adv.evalue,
        identity: drawing.adv.identity,
        alignment_length: drawing.adv.alignment_length
      }
    });
  });
  // The pressed global action follows the projected intent (UJ-05).
  const linearComparisonGlobalAction = computed(() => {
    const { intentKey } = linearComparisonUi.value;
    return intentKey === 'custom' ? 'selected' : intentKey;
  });
  const canRunLinearLosat = computed(() => {
    const drawing = shownDrawingOfMode('linear');
    return drawing !== null && drawing.drawnLinearSeqs.value.filter((sequence) => (
      inputSourceFiles(drawing.inputType.value, sequence).complete
    )).length >= 2;
  });
  return {
    comparisonHeightValidationError, invalidateLinearComparisonArtifacts, linearComparisonUi,
    linearComparisonGlobalAction, canRunLinearLosat
  };
};

/**
 * @typedef {object} LinearComparisonTimelinePorts
 * @property {WebState} state
 * @property {() => readonly LinearSeq[]} shownLinearSeqs The record rows of the shown Linear drawing.
 */

/**
 * The pair timeline of the shown drawing, a pair by its edge key, and the record names its rows show.
 * @param {LinearComparisonTimelinePorts} ports
 */
export const createLinearComparisonTimeline = ({ state, shownLinearSeqs }) => {
  const linearComparisonTimeline = computed(() => {
    const drawing = state.activeDrawing();
    return buildLinearComparisonTimeline({
      sequences: drawing.linearSeqs,
      recordsOff: drawing.recordsOff,
      layout: effectiveLinearComparisonLayout(drawing),
      plan: drawing.linearComparisonPlan,
      resolution: drawing.linearComparisonResolution.value
    });
  });
  /** @param {string} edgeKey */
  const linearComparisonPairForEdgeKey = (edgeKey) => {
    for (const row of linearComparisonTimeline.value.rows) {
      const pair = row.boundaryAfter?.pairs.find((/** @type {LinearComparisonPair} */ entry) => entry.edgeKey === edgeKey);
      if (pair) return pair;
    }
    return null;
  };
  // Every panel names a record the same way (OV-436).
  const linearRecordNameMap = computed(() => linearRecordNames(shownLinearSeqs()));
  /** @param {string} uid */
  const linearComparisonRecordLabel = (uid) => linearRecordNameMap.value.get(uid)?.label || t('Record');
  // The parts <record-name> lays out: the name shrinks, the accession stays whole.
  /** @param {string} uid */
  const linearComparisonRecordName = (uid) => {
    const entry = linearRecordNameMap.value.get(uid);
    return { name: entry?.name || t('Record'), detail: entry?.accession || '' };
  };
  return {
    linearComparisonTimeline, linearComparisonPairForEdgeKey, linearComparisonRecordLabel, linearComparisonRecordName
  };
};

/**
 * @typedef {object} LinearComparisonFocusPorts
 * @property {WebState} state
 * @property {{ readonly value: DiagramMode }} mode The shown drawing's mode.
 * @property {{ readonly value: UserFacingError | null }} errorLog The banner's error.
 * @property {{ readonly value: { readonly issues: readonly LinearComparisonFocusTarget[] } }} linearComparisonUi
 *   The Comparison card's projection: its issues.
 * @property {RevealSetting} revealSetting The settings pane shows and focuses a control.
 * @property {RevealSettingsSection} revealSettingsSection The settings pane shows a section and focuses its heading.
 */

/**
 * The paths into the Comparison card: a disclosure that opens itself while LOSAT compares the shown drawing, a pair's
 * row, an issue's control, a pair's BLAST TSV chooser, and the cause of a Generate error.
 * @param {LinearComparisonFocusPorts} ports
 */
export const createLinearComparisonFocus = ({ state, mode, errorLog, linearComparisonUi, revealSetting, revealSettingsSection }) => {
  /** @param {string} disclosureKey */
  const openLinearComparisonDisclosure = async (disclosureKey) => {
    if (!disclosureKey) return null;
    const details = /** @type {HTMLDetailsElement[]} */ ([...document.querySelectorAll('[data-linear-comparison-disclosure]')])
      .find((element) => element.dataset.linearComparisonDisclosure === disclosureKey);
    if (!details) return null;
    details.open = true;
    await nextTick();
    return details;
  };
  watch([mode, () => state.activeDrawing().hasActiveLinearLosatIntent.value], (/** @type {[DiagramMode, boolean]} */ [activeMode, activeLosat]) => {
    if (activeMode === 'linear' && activeLosat) openLinearComparisonDisclosure('settings');
  }, { flush: 'post' });

  /** @param {string} edgeKey */
  const focusLinearComparisonPair = async (edgeKey) => {
    await openLinearComparisonDisclosure('selected-pairs');
    const container = /** @type {HTMLElement[]} */ ([...document.querySelectorAll('[data-edge-key]')])
      .find((element) => element.dataset.edgeKey === edgeKey);
    const control = /** @type {HTMLElement | null | undefined} */ (container?.querySelector('input, select, button'));
    if (control) await revealSetting(control);
  };

  /** @param {{ disclosureKey: string, focusTargetKey: string, edgeKey?: string, edgeId?: string } | null | undefined} target */
  const focusLinearComparisonIssue = async (target = linearComparisonUi.value.issues[0]) => {
    if (!target) return false;
    const details = await openLinearComparisonDisclosure(target.disclosureKey);
    /** @type {HTMLElement | null} */
    let container = details;
    if (target.edgeKey) {
      container = /** @type {HTMLElement[]} */ ([...document.querySelectorAll('[data-edge-key]')])
        .find((element) => element.dataset.edgeKey === target.edgeKey) || container;
    }
    if (target.edgeId && target.focusTargetKey === 'pair-row') {
      container = /** @type {HTMLElement[]} */ ([...document.querySelectorAll('[data-linear-unplaced-draft]')])
        .find((element) => element.dataset.linearUnplacedDraft === target.edgeId) || container;
    }
    const marked = container?.querySelector(
      `[data-comparison-focus="${target.focusTargetKey}"]`
    );
    const focusable = marked?.matches?.('input, select, button, [tabindex]')
      ? marked
      : marked?.querySelector?.('input, select, button, [tabindex]');
    const control = /** @type {HTMLElement | null | undefined} */ (focusable || container?.querySelector?.('input, select, button, [tabindex]'));
    if (control) await revealSetting(control);
    return true;
  };

  // The one way to choose a pair's BLAST TSV outside its row:
  // the pair's own uploader opens its file chooser inside the click, so the
  // file takes the uploader's availability check and History step.
  /** @param {string} edgeKey @param {{ focus?: boolean }} [options] */
  const chooseLinearComparisonPairFile = (edgeKey, { focus = false } = {}) => {
    void openLinearComparisonDisclosure('selected-pairs');
    const uploader = /** @type {HTMLElement[]} */ ([...document.querySelectorAll('[data-linear-comparison-pair-upload]')])
      .find((element) => element.dataset.linearComparisonPairUpload === edgeKey);
    const button = /** @type {HTMLElement | null | undefined} */ (uploader?.querySelector('[role="button"]'));
    if (focus && button) void revealSetting(button);
    /** @type {HTMLInputElement | null | undefined} */ (uploader?.querySelector('input[type="file"]'))?.click();
    return Boolean(uploader);
  };

  // Whether the Generate error offers "Show in Comparison": in Linear for
  // every error that leads to the comparison settings; in Circular, whose
  // comparison rings have no issue list, for a comparison setting the same
  // Generate stops on again (a LOSAT failure that repeats, OV-3317, OV-3327,
  // or a genetic code LOSAT does not read, OV-3326), whose options and
  // inputs are in Pairwise Comparisons.
  const comparisonCauseShown = computed(() => {
    const error = normalizeUserFacingError(errorLog.value);
    return Boolean(error?.actions.includes('edit-comparison')) && (mode.value === 'linear' || (error !== null && comparisonCauseRepeats(error)));
  });
  // The Generate error's "Show in Comparison": in Circular, the Pairwise
  // Comparisons section; in Linear, the refused card's LOSAT Gencode, found
  // by the card's identity (OV-3326), the first issue's control, the LOSAT
  // mode when the same search would stop the same way, the Execution setting
  // for another LOSAT runtime failure, or the Comparison card when the cause
  // has no control.
  const showComparisonCause = async () => {
    if (mode.value !== 'linear') {
      await revealSettingsSection('pairwise-comparisons');
      return;
    }
    const error = normalizeUserFacingError(errorLog.value);
    const code = error?.code;
    const cardUid = error?.context.reason === 'GENETIC_CODE' ? error.context.cardUid : undefined;
    const gencodeInput = cardUid === undefined ? undefined
      : /** @type {HTMLElement[]} */ ([...document.querySelectorAll('[data-linear-gencode-input]')])
        .find((element) => element.dataset.linearGencodeInput === cardUid);
    if (gencodeInput) {
      await revealSetting(gencodeInput);
      return;
    }
    const runtimeCause = error && losatFailureRepeats(error) ? { disclosureKey: 'settings', focusTargetKey: 'losat-options' }
      : code === 'LOSAT_THREADING_UNAVAILABLE' || code === 'LOSAT_RUNTIME'
        ? { disclosureKey: 'settings', focusTargetKey: 'losat-execution' } : null;
    if (await focusLinearComparisonIssue(linearComparisonUi.value.issues[0] || runtimeCause)) return;
    const card = /** @type {HTMLElement | null} */ (document.querySelector('[data-linear-comparison-card]'));
    if (!card) return;
    await revealSetting(card, {
      align: 'start', focus: /** @type {HTMLElement | null} */ (card.querySelector('button:not([disabled])')) ?? undefined
    });
  };
  return {
    openLinearComparisonDisclosure, focusLinearComparisonPair, focusLinearComparisonIssue, chooseLinearComparisonPairFile,
    comparisonCauseShown, showComparisonCause
  };
};

/**
 * @typedef {object} LinearLayoutReadersPorts
 * @property {WebState} state
 * @property {ShownDrawingOfMode} shownDrawingOfMode
 * @property {LosatCacheInfoRows} losatCacheInfo
 */

/**
 * The File card's row input of the shown Linear drawing, its record row tokens, the LOSAT cache rows by pair, and the
 * window pairs of each card pair.
 * @param {LinearLayoutReadersPorts} ports
 */
export const createLinearLayoutReaders = ({ state, shownDrawingOfMode, losatCacheInfo }) => {
  /**
   * The File card's row input: the row of `uid` in the shown Linear drawing,
   * or `fallback` (its position) when the drawing gives it none.
   * @param {string} uid @param {number} fallback
   */
  const shownLinearRecordRowFor = (uid, fallback) => {
    const drawing = shownDrawingOfMode('linear');
    return drawing === null ? fallback : linearRecordRowFor(drawing, uid, fallback);
  };
  const linearLayoutTokens = computed(() => {
    const drawing = state.activeDrawing();
    return (
      drawing.linearRecordLayoutEnabled.value
        ? linearRecordPositionTokens(drawing.linearSeqs, drawing.linearRecordRows)
        : []
    );
  });
  // The shown Linear drawing's rows, its own (D-221, OV-3210).
  /** @returns {LosatCacheInfoEntry[]} */
  const cacheInfoRows = () => (shownDrawingOfMode('linear') === null || !Array.isArray(losatCacheInfo.value) ? [] : losatCacheInfo.value);
  const linearLosatCacheInfoByEdgeKey = computed(() => Object.fromEntries(
    cacheInfoRows()
      .filter((entry) => String(entry?.edgeKey || ''))
      .map((entry) => [String(entry.edgeKey), entry])
  ));
  // The window pairs of each card pair, by its edge key: one row each, with its own Save Raw LOSAT TSV and file name
  // (OV-2742). The units come from the drawn Result's catalog, which the Generate that wrote the rows committed.
  /** @type {{ readonly value: Record<string, LinearWindowPair[]> }} */
  const linearLosatWindowPairs = computed(() => Object.fromEntries(
    linearWindowPairsByCardEdge(cacheInfoRows(), catalogUnitPlaces(state.featureCatalog.value))
  ));
  /**
   * A side of a window pair as its row names it after `#`: the card's number, and the window's when the card draws
   * several ("1 window 2").
   * @param {number} cardIndex @param {number | null} window
   */
  const linearWindowPairSide = (cardIndex, window) => (
    window === null ? String(cardIndex + 1) : t('{n} window {window}', { n: cardIndex + 1, window })
  );
  return {
    shownLinearRecordRowFor, linearLayoutTokens, linearLosatCacheInfoByEdgeKey, linearLosatWindowPairs,
    linearWindowPairSide
  };
};
