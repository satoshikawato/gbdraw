// @ts-check
// The actions of the Result and of the Generate choices that come with it:
// copy the run info commands, choose what a loaded Session's comparison
// becomes, answer the missing BLAST TSV dialog, copy the clicked feature, and
// export the shown Result. The composition root passes the state and the
// owners' functions they use as ports.
import { copyTextToClipboard } from '../utils/clipboard.js';
import { normalizeUserFacingError } from '../utils/error-normalization.js';
import { SOURCE_LOCALE, msg, numberText, t, tm, uiLocale } from '../utils/i18n.js';
import { downloadTextFile } from '../services/text-download.js';
import { exportLinearWrap } from '../services/linear-wrap-form.js';
import { requestedFigureSize } from '../services/figure-size-form.js';
import { missingUploadPairsToResolve } from '../services/linear-comparisons.js';
import {
  IMPORTED_COMPARISON_ACTIONS,
  IMPORTED_COMPARISON_DISPOSITIONS,
  resolveImportedComparisonAction
} from '../services/imported-comparison-intent.js';
import { formatElapsedMs, reproducibilityLabel } from './run-info.js';
/** @import { DrawingState, WebState } from '../state.js' */
/** @import { DrawingId } from '../services/drawing-view.js' */
/** @import { LinearComparisonIssue, LinearComparisonPair } from '../services/linear-comparisons.js' */
/** @import { LinearDraftRecord } from '../services/linear-wrap-form.js' */
/** @import { InputSummary } from '../services/input-summary.js' */
/** @import { CanonicalRenderRequest } from '../contracts/request.generated.js' */
/** @import { UserFacingError } from '../utils/error-normalization.js' */
/** @import { ShowNotice } from './app-messages.js' */
/** @import { ExportLimitReply, ExportSnapshot, FontEmbedding, SvgRootReply } from '../services/export.js' */
/** @import { PdfFontsReply } from '../services/pdf-fonts.js' */
/** @import { CatalogSequences } from '../services/catalog-sequences.js' */

const { computed, reactive, ref } = window.Vue;

/** @typedef {{ value: UserFacingError | null }} ErrorLogRef The banner's error; null when there is none. */

/**
 * The run info of the shown Result whose commands the run info panel copies.
 * @typedef {object} CopiedRunInfo
 * @property {unknown} [elapsedMs]
 * @property {{ level?: unknown } | null} [reproducibility]
 * @property {unknown[]} [helperFiles]
 * @property {unknown} [command]
 * @property {unknown} [sessionCommand]
 * @property {{ command?: unknown } | null} [sourceRecipe]
 * @property {{ command?: unknown } | null} [exactReplay]
 */

/**
 * The run info panel's Copy buttons and texts.
 * @param {{ lastRunInfo: { value: CopiedRunInfo | null } }} options The run info of the shown Result.
 */
export const createRunInfoCopy = ({ lastRunInfo }) => {
  const runInfoCopyStatus = ref('');
  const exactReplayCopyStatus = ref('');

  /** @param {{ elapsedMs?: unknown } | null | undefined} info */
  const runInfoElapsedText = (info) => formatElapsedMs(info?.elapsedMs);
  /** @param {{ reproducibility?: { level?: unknown } | null } | null | undefined} info */
  const runInfoReproducibilityText = (info) => reproducibilityLabel(info?.reproducibility?.level);
  const runInfoHasCliHelperFiles = computed(() =>
    Array.isArray(lastRunInfo.value?.helperFiles) && lastRunInfo.value.helperFiles.length > 0
  );
  /** @param {unknown} commandValue @param {{ value: string }} status */
  const copyRunInfoCommand = async (commandValue, status) => {
    const command = String(commandValue || '');
    if (!command) return;
    try {
      await copyTextToClipboard(command);
      const copied = t('Copied!');
      status.value = copied;
      setTimeout(() => {
        if (status.value === copied) status.value = '';
      }, 1600);
    } catch (error) {
      console.warn('Failed to copy the requested command.', normalizeUserFacingError(error));
      const failed = t('Copy failed');
      status.value = failed;
      setTimeout(() => {
        if (status.value === failed) status.value = '';
      }, 2200);
    }
  };
  const copyRunCommand = () => copyRunInfoCommand(
    lastRunInfo.value?.sourceRecipe?.command || lastRunInfo.value?.command,
    runInfoCopyStatus
  );
  const copyExactReplayCommand = () => copyRunInfoCommand(
    lastRunInfo.value?.exactReplay?.command || lastRunInfo.value?.sessionCommand,
    exactReplayCopyStatus
  );
  return {
    runInfoCopyStatus,
    exactReplayCopyStatus,
    runInfoElapsedText,
    runInfoReproducibilityText,
    runInfoHasCliHelperFiles,
    copyRunCommand,
    copyExactReplayCommand
  };
};

/**
 * @typedef {object} ImportedComparisonChoicePorts
 * @property {WebState} state
 * @property {(label: string, fn: () => unknown) => Promise<unknown>} runUndoable History's undoable step.
 * @property {ErrorLogRef} errorLog
 * @property {(drawing: DrawingState, nextPlan: unknown) => unknown} replaceLinearComparisonPlan
 *   The comparison owner's replacement of a drawing's plan.
 */

/**
 * The choice of what a loaded Session's saved comparison becomes.
 * @param {ImportedComparisonChoicePorts} ports
 */
export const createImportedComparisonChoice = ({ state, runUndoable, errorLog, replaceLinearComparisonPlan }) => {
  // The buttons that choose what happens to a saved comparison, the History
  // step of the choice, and the line that names the chosen one. English names
  // it as the Session stores it (INHERIT, REPLACE, CLEAR; PD-OI-106 keeps the
  // English text); another language names it by its button.
  /** @type {Readonly<Record<string, string>>} */
  const importedComparisonActionLabels = Object.freeze({
    [IMPORTED_COMPARISON_ACTIONS.INHERIT]: msg('Inherit saved comparison'),
    [IMPORTED_COMPARISON_ACTIONS.REPLACE]: msg('Replace with current controls'),
    [IMPORTED_COMPARISON_ACTIONS.CLEAR]: msg('Clear comparison')
  });
  /** @param {DrawingState} drawing @param {string} action */
  const chooseImportedComparisonAction = (drawing, action) => runUndoable(
    importedComparisonActionLabels[action],
    () => {
      const outcome = resolveImportedComparisonAction({
        intent: drawing.importedComparisonIntent,
        action,
        draftResolution: drawing.linearComparisonResolution.value
      });
      if (!outcome.ok) {
        errorLog.value = normalizeUserFacingError(outcome.message, { operation: 'generate', stage: 'request-validation' });
        return false;
      }
      if (outcome.action === IMPORTED_COMPARISON_ACTIONS.CLEAR) {
        replaceLinearComparisonPlan(drawing, { mode: 'none', defaultSource: 'losat', edges: [] });
      }
      errorLog.value = null;
      return true;
    }
  );
  const inheritImportedComparison = () => {
    const drawing = state.activeDrawing();
    return chooseImportedComparisonAction(
      drawing,
      IMPORTED_COMPARISON_ACTIONS.INHERIT
    );
  };
  const replaceImportedComparison = () => {
    const drawing = state.activeDrawing();
    return chooseImportedComparisonAction(
      drawing,
      IMPORTED_COMPARISON_ACTIONS.REPLACE
    );
  };
  const clearImportedComparison = () => {
    const drawing = state.activeDrawing();
    return chooseImportedComparisonAction(
      drawing,
      IMPORTED_COMPARISON_ACTIONS.CLEAR
    );
  };
  const importedComparisonNeedsResolution = computed(() => (
    state.activeDrawing().importedComparisonIntent.disposition !== IMPORTED_COMPARISON_DISPOSITIONS.EDITABLE
  ));
  const importedComparisonActionText = computed(() => {
    /** @type {string | null} */
    const action = state.activeDrawing().importedComparisonIntent.action;
    if (!action || uiLocale.value === SOURCE_LOCALE || !Object.hasOwn(importedComparisonActionLabels, action)) return action ?? '';
    return tm(importedComparisonActionLabels[action]);
  });
  const importedComparisonCanInherit = computed(() => (
    state.activeDrawing().importedComparisonIntent.disposition
      === IMPORTED_COMPARISON_DISPOSITIONS.PRESERVED_READ_ONLY
  ));
  const importedComparisonCanReplace = computed(() => {
    const drawing = state.activeDrawing();
    return (
      drawing.linearComparisonResolution.value.valid
      && drawing.linearComparisonResolution.value.hasComparisonIntent
    );
  });
  return {
    importedComparisonActionLabels,
    inheritImportedComparison,
    replaceImportedComparison,
    clearImportedComparison,
    importedComparisonNeedsResolution,
    importedComparisonActionText,
    importedComparisonCanInherit,
    importedComparisonCanReplace
  };
};

/**
 * @typedef {object} MissingBlastTsvChoicePorts
 * @property {WebState} state
 * @property {(label: () => string, handler: (choice: string) => unknown, cancel: () => unknown) => (choice: string) => unknown} withHistory
 *   The dialog choice's History step.
 * @property {(close: () => unknown) => void} closeAfterChoice The dialog choice's close once its step ends.
 * @property {{ value: boolean }} choicePending Whether a dialog choice's step is running.
 * @property {(drawing: DrawingState, edgeKeys: string[]) => void} omitLinearComparisonPairs
 * @property {(edgeKey: string) => LinearComparisonPair | null} linearComparisonPairForEdgeKey
 * @property {() => unknown} runAnalysis Generate.
 */

/**
 * The dialog Generate opens for Upload pairs without a BLAST TSV file.
 * @param {MissingBlastTsvChoicePorts} ports
 */
export const createMissingBlastTsvChoice = ({
  state, withHistory, closeAfterChoice, choicePending, omitLinearComparisonPairs, linearComparisonPairForEdgeKey, runAnalysis
}) => {
  // Generate asks about Upload pairs without a file: each pair's row
  // takes its own file; Skip sets the pairs still without one to No comparison
  // in one History step, then Generates.
  /** @type {{ show: boolean, drawingId: DrawingId, pairs: ReturnType<typeof missingUploadPairsToResolve> }} */
  const missingBlastTsvDialog = reactive({ show: false, drawingId: state.activeDrawingId.value, pairs: [] });
  const closeMissingBlastTsvDialog = () => { missingBlastTsvDialog.show = false; };
  const setMissingBlastTsvPairsToNone = withHistory(
    () => msg('Set pairs to No comparison'),
    () => {
      // The pairs as they are now: Undo may have changed them while the dialog was open.
      const drawing = state.drawingOf(missingBlastTsvDialog.drawingId);
      const pairs = missingUploadPairsToResolve(drawing.linearComparisonResolution.value);
      if (pairs.length) omitLinearComparisonPairs(drawing, pairs.map((entry) => entry.edgeKey));
      closeAfterChoice(closeMissingBlastTsvDialog);
      return true;
    },
    closeMissingBlastTsvDialog
  );
  // The dialog's rows: the pairs it opened with, each with its file as it is now.
  const missingBlastTsvRows = computed(() => {
    const missing = new Set(/** @type {readonly LinearComparisonIssue[]} */ (state.drawingOf(missingBlastTsvDialog.drawingId).linearComparisonResolution.value?.errors || [])
      .filter((issue) => issue.code === 'missing-upload').map((issue) => issue.edgeKey));
    return missingBlastTsvDialog.pairs.map((pair) => {
      const found = linearComparisonPairForEdgeKey(pair.edgeKey);
      return { ...pair, queryUid: found?.queryUid || '', subjectUid: found?.subjectUid || '',
        fileName: found?.draft?.fileActive ? String(found.draft.file?.name || '') : '', missing: missing.has(pair.edgeKey) };
    });
  });
  /** @param {'generate' | 'none' | 'cancel'} choice */
  const resolveMissingBlastTsv = async (choice) => {
    if (choice === 'generate') {
      if (choicePending.value) return undefined;
      closeMissingBlastTsvDialog();
      return runAnalysis();
    }
    return await setMissingBlastTsvPairsToNone(choice) === true ? runAnalysis() : undefined;
  };
  return { missingBlastTsvDialog, missingBlastTsvRows, resolveMissingBlastTsv };
};

/**
 * The feature popup's header and its Copy buttons.
 * @param {Pick<WebState, 'clickedFeature'>} options The open feature popup's payload.
 */
export const createClickedFeatureCopy = ({ clickedFeature }) => {
  // The popup header names the gene (else the feature label), then
  // its product, then "<record ID> · <location>" with grouped digits (parts
  // joined by "; ", since a grouped number holds a comma), all from payload
  // fields only.
  const clickedFeatureHeading = computed(() => {
    const cf = clickedFeature.value;
    const title = String(cf?.gene || cf?.label || '');
    return { title, product: cf?.product && cf.product !== title ? cf.product : '' };
  });
  const clickedFeatureSummary = computed(() => {
    const cf = clickedFeature.value;
    const location = String(cf?.location || '').replaceAll(', ', '; ')
      .replace(/\d{4,}/g, (digits) => numberText(Number(digits)));
    return [cf?.recordId, location].filter(Boolean).join(' · ');
  });

  // One Copy per tab: the header lines and the rows, or every qualifier.
  const copyClickedFeatureDetails = () => {
    const cf = clickedFeature.value;
    if (!cf) return;
    const { title, product } = clickedFeatureHeading.value;
    copyTextToClipboard([title, product, clickedFeatureSummary.value, ...cf.detailRows.map((row) => `${row.label}: ${row.value}`)]
      .filter(Boolean).join('\n'));
  };
  const copyClickedFeatureQualifiers = () => {
    const rows = clickedFeature.value?.qualifierRows || [];
    copyTextToClipboard(rows.flatMap((/** @type {{ key: string, values: string[] }} */ row) => row.values.map((value) => `/${row.key}="${value}"`)).join('\n'));
  };
  return { clickedFeatureHeading, clickedFeatureSummary, copyClickedFeatureDetails, copyClickedFeatureQualifiers };
};

/**
 * @typedef {(state: WebState, drawing: DrawingState, options: { interactive: boolean }) => ExportSnapshot} CaptureSvgExport
 *   The capture of the shown Result, taken before the export modules load.
 * @typedef {(operation: string, payload: Record<string, unknown>) => Promise<{ result: unknown }>} RunHelperOperation
 *   One operation of the Worker's helpers; `result` is the operation's reply.
 * @typedef {{ EMBED_FONTS: string, PDF_FONTS: string, CHECK_EXPORT_LIMIT: string }} ExportHelperOperations
 *   The names of the helper operations an export runs.
 */

/**
 * @typedef {object} ResultExportPorts
 * @property {WebState} state
 * @property {ErrorLogRef} errorLog
 * @property {ShowNotice} showNotice
 * @property {CaptureSvgExport} captureSvgExport
 * @property {() => Promise<typeof import('../services/export.js')>} loadExportService
 *   The composition root's lazy load of the export modules.
 * @property {RunHelperOperation} runDiagramHelperOperation
 * @property {ExportHelperOperations} DIAGRAM_HELPER_OPERATIONS
 * @property {() => DrawingState} displayedArtifactDrawing The drawing the shown Result belongs to.
 * @property {{ value: boolean }} embedInExports Embed font in exported SVG.
 * @property {{ value: InputSummary | null }} shownInputSummary Python's input summary of the shown Result.
 * @property {(drawing: DrawingState) => unknown} shownLinearWrapLeaf
 * @property {(drawing: DrawingState) => CanonicalRenderRequest | null} shownCommittedRequest
 * @property {(drawing: DrawingState) => Promise<void>} awaitLinearRecordDiscovery
 * @property {(drawing: DrawingState) => LinearDraftRecord[]} committedLinearRecordsOf
 * @property {CatalogSequences['readEmbedJson']} readEmbedJson Python's encoding of the shown Result's catalog item
 *   with every sequence an Interactive SVG embeds (D5); null when another Result is shown before it arrives.
 */

/**
 * The exports of the shown Result, and the text download.
 * @param {ResultExportPorts} ports
 */
export const createResultExport = ({
  state, errorLog, showNotice, captureSvgExport, loadExportService, runDiagramHelperOperation, DIAGRAM_HELPER_OPERATIONS,
  displayedArtifactDrawing, embedInExports, shownInputSummary, shownLinearWrapLeaf, shownCommittedRequest,
  awaitLinearRecordDiscovery, committedLinearRecordsOf, readEmbedJson
}) => {
  /** @param {unknown} filename @param {unknown} text @param {string} [type] */
  const downloadText = (filename, text, type = 'text/plain;charset=utf-8') => {
    const value = String(text ?? '');
    if (!value) return;
    downloadTextFile(String(filename || 'gbdraw.txt'), value, type);
  };

  let latestExportOperation = 0;
  const failedInteractiveSvgExport = ref(null);
  const canRetryInteractiveSvgExport = computed(() => Boolean(failedInteractiveSvgExport.value
    && errorLog.value === failedInteractiveSvgExport.value));
  /** @param {'downloadSVG' | 'downloadInteractiveSVG' | 'downloadPNG' | 'downloadPDF'} methodName @param {string} operation */
  const runExportAction = async (methodName, operation) => {
    const operationId = ++latestExportOperation;
    const previousError = errorLog.value;
    try {
      const captured = captureSvgExport(state, displayedArtifactDrawing(), { interactive: methodName === 'downloadInteractiveSVG' });
      // An Interactive SVG embeds every sequence, which Python reads from the request the shown Result was made
      // from (D5) before the file is built: a failed read writes no file and shows why (Retry); a Result shown
      // before the read ends writes nothing and says nothing.
      const { interactivity } = captured;
      const embeddedSequences = interactivity?.catalogItem ? await readEmbedJson(interactivity.catalogItem) : null;
      if (interactivity?.catalogItem && !embeddedSequences) return { status: 'stale' };
      const snapshot = interactivity && embeddedSequences ? {
        ...captured, interactivity: { ...interactivity, embeddedSequences }
      } : captured;
      const exportService = await loadExportService();
      const exportMethod = exportService?.[methodName];
      if (typeof exportMethod !== 'function') {
        throw new Error('The export service did not provide the requested action.');
      }
      // Custom fonts in exports (D-05, D-13): the Worker holds the fonts of
      // this page, so only a page with a loaded font embeds. SVG and
      // Interactive SVG follow Embed font in exported SVG; a PNG always draws
      // with the fonts, and a PDF always registers the faces its text uses.
      // The size check of a PNG or PDF reads the exported Result's wrap
      // (D-52, D-55): Python reads its stored leaf, and decides from it
      // whether a record is cut, and the row length it advises. A loaded
      // Result has no input summary: its records are discovered first, and
      // Python crops each by the region its request drew (D-70, D-73).
      const exportedDrawing = displayedArtifactDrawing();
      const exportedSummary = shownInputSummary.value ?? undefined;
      const exportedWrapLeaf = shownLinearWrapLeaf(exportedDrawing);
      const exportedRequest = shownCommittedRequest(exportedDrawing);
      const exportedFigureSize = exportedRequest ? requestedFigureSize(exportedRequest) : null;
      const linearWrap = () => exportLinearWrap({
        leaf: exportedWrapLeaf,
        summary: exportedSummary,
        discoverRecords: async () => {
          await awaitLinearRecordDiscovery(exportedDrawing);
          return committedLinearRecordsOf(exportedDrawing);
        }
      });
      /** @type {import('../services/export.js').ExportPorts} */
      const exportPorts = {
        embedFonts: state.loadedFonts.value.length > 0
          ? async (svg, kind) => /** @type {FontEmbedding} */ (
            (await runDiagramHelperOperation(DIAGRAM_HELPER_OPERATIONS.EMBED_FONTS, { svg, kind })).result)
          : null,
        embedInSvg: embedInExports.value,
        resolvePdfFonts: async (runs, known) => /** @type {PdfFontsReply} */ ((
          await runDiagramHelperOperation(DIAGRAM_HELPER_OPERATIONS.PDF_FONTS, { runs, known })
        ).result),
        checkExportLimit: async (request) => /** @type {ExportLimitReply} */ ((
          await runDiagramHelperOperation(DIAGRAM_HELPER_OPERATIONS.CHECK_EXPORT_LIMIT, {
            ...request, linearWrap: await linearWrap()
          })
        ).result),
        // A Result drawn to a requested size states it in its SVG as the CLI writes it (Python's
        // `export_root_size`, through the export check); a Result drawn at Auto exports its root as shown.
        sizedSvgRoot: exportedFigureSize === null ? null : async ({ width, height }) => /** @type {SvgRootReply} */ ((
          await runDiagramHelperOperation(DIAGRAM_HELPER_OPERATIONS.CHECK_EXPORT_LIMIT, {
            format: 'svg', width, height, dpi: null, linearWrap: null, figureSize: exportedFigureSize
          })
        ).result).root
      };
      /** @type {import('../services/export.js').ExportOutcome} */
      const result = await exportMethod(snapshot, exportPorts);
      if (operationId === latestExportOperation && errorLog.value === previousError && previousError?.operation?.startsWith('export-')) errorLog.value = null;
      if (result.fontWarnings.length > 0) {
        showNotice({ text: t('Exported {file}.', { file: result.filename }), items: result.fontWarnings });
      }
      return result;
    } catch (error) {
      const normalized = normalizeUserFacingError(error, { operation, stage: 'export-capture' });
      if (operationId !== latestExportOperation || errorLog.value !== previousError) return { status: 'stale' };
      errorLog.value = normalized;
      failedInteractiveSvgExport.value = methodName === 'downloadInteractiveSVG' ? normalized : null;
      return { status: 'error', error: normalized };
    }
  };

  const downloadSVG = () => runExportAction('downloadSVG', 'export-svg');
  const downloadInteractiveSVG = () => (
    runExportAction('downloadInteractiveSVG', 'export-svg')
  );
  const downloadPNG = () => runExportAction('downloadPNG', 'export-png');
  const downloadPDF = () => runExportAction('downloadPDF', 'export-pdf');
  return { downloadText, canRetryInteractiveSvgExport, downloadSVG, downloadInteractiveSVG, downloadPNG, downloadPDF };
};
