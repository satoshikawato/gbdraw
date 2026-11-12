// @ts-check
// Load and Export Label TSV (the bulk label table): the file read, the labels
// a row selects in every Result of a batch, the declines of a load whose
// inputs changed while it waited, the failure and its Retry, and the export.
// The label owner (`label-actions.js`) writes the label edits a load makes
// and keeps the notice of the last load; this owner reaches it through the
// ports the composition root (`feature-editor.js`) passes (R13).
/** @import { NoticeHandle, ShowNotice } from '../app-messages.js' */
/** @import { DrawingState, WebState } from '../../state.js' */
/** @import { RulePreparation, RulePreparationOptions } from '../rule-matching.js' */
/** @import { UserFacingError } from '../../utils/error-normalization.js' */
/** @import { LabelImportEdit, LabelImportTarget } from '../../services/label-import-plan.js' */
/** @import { LabelFeature } from './label-actions.js' */
import { normalizeUserFacingError } from '../../utils/error-normalization.js';
import { DRAWN_SELECTOR_QUALIFIERS, drawnSelectorUnknown, labelRulePayloads } from '../rule-matching.js';
import { parseLabelOverrideTsv } from '../../services/label-override-table.js';
import { labelImportNotice, labelIntentSignature, labelTableExportLines, planLabelImport } from '../../services/label-import-plan.js';
import { featureIdentityKeyOf } from '../../services/feature-placement.js';
import { resultRenderedFeatures } from '../../services/feature-catalog.js';
import { downloadTextFile } from '../../services/text-download.js';
import { createRetryCopyHold, holdFileCopy, readFileText } from '../../services/file-content-cache.js';
import { LABEL_FEATURE_ID_ATTRIBUTE } from '../../services/drawn-groups.js';
import { t } from '../../utils/i18n.js';

/**
 * The input a Label TSV change reads: a file input, or the object the retry
 * passes in its place. Edit opens a file input's picker again.
 * @typedef {{ files: ArrayLike<File> | null, value: string, click?: () => void }} LabelTableInput
 */
/**
 * A failed load that offers Retry while its error shows and the inputs it
 * read are current.
 * @typedef {object} LabelImportFailure
 * @property {UserFacingError} error The error it shows.
 * @property {ReturnType<RulePreparation['snapshot']>} snapshot The rule inputs the load read.
 * @property {string} intent The label intent the load read (`labelIntentSignature`).
 * @property {number} revision The load's revision; a newer load ends the offer.
 * @property {LabelTableInput} input The input Edit opens again.
 * @property {() => Promise<unknown>} retry
 */

/**
 * @typedef {object} LabelTableActionsOptions
 * @property {Pick<WebState, 'activeDrawing' | 'errorLog' | 'featureCatalog' | 'mode' | 'results' | 'selectedResultIndex'
 *   | 'sessionOperationAvailability' | 'svgContainer'>} state
 *   The members of the app state (state.js) this owner and the helpers it passes `state` to read.
 * @property {<T>(value: T) => { value: T }} ref Vue `ref`
 * @property {<T>(getter: () => T) => { value: T }} computed Vue `computed`
 * @property {RulePreparation['isCurrent']} isRuleInputCurrent Whether a snapshot of the rule inputs is still the current inputs.
 * @property {RulePreparation['snapshot']} snapshotRuleInputs The inputs the rule matches depend on.
 * @property {RulePreparationOptions['evaluate']} evaluateLabelRules
 *   The root's Python evaluation of Label TSV rows against the displayed labels (R7).
 * @property {ShowNotice} showNotice The root's notice banner: a table import or export reports there.
 * @property {(drawing: DrawingState, svg: SVGSVGElement) => LabelImportTarget[]} labelTargetsOf
 *   The label owner's labels of the displayed Result: each editable label, its feature, and its source text.
 * @property {(drawing: DrawingState, element: Element, feature: LabelFeature | null) => string} labelSourceText
 *   The label owner's source text of a label of another Result of the batch.
 * @property {(drawing: DrawingState, applicable: readonly LabelImportEdit[]) => void} applyLabelTable
 *   The label owner's replacement of the drawing's label edits by those a loaded table makes.
 * @property {(drawing: DrawingState, notice: NoticeHandle | void) => void} rememberLabelImportNotice
 *   Hands the label owner the notice of the load that applied, which Reset label edits ends.
 */

/** @param {unknown} name @param {string} [fallback] */
const makeSafeFilename = (name, fallback = 'gbdraw') => {
  const cleaned = String(name || '')
    .replace(/[^\w.-]+/g, '_')
    .replace(/^_+|_+$/g, '');
  return cleaned || fallback;
};

// A `location` or `record_location` row cannot be decided for a label whose
// drawn values the Result does not record (feature catalog 3 or 4 on a
// cropped, reverse-complemented, or rotated record). The import is declined
// whole, so the label edits it would replace stay (R4).
/**
 * @param {readonly LabelImportTarget[]} labels
 * @param {readonly import('../../services/label-override-table.js').LabelOverrideTsvRow[]} rows
 */
const undecidedLabelCount = (labels, rows) => labels.filter(({ feature }) => (
  // `drawnSelectorUnknown` holds only for a feature that has a `drawnSelector` key, so `feature` is set below.
  drawnSelectorUnknown(feature) && rows.some((row) => (
    DRAWN_SELECTOR_QUALIFIERS.has(row.qualifier.toLowerCase())
    && (row.recordId === '*' || row.recordId === String(/** @type {LabelFeature} */ (feature).record_id ?? ''))
    && (row.featureType === '*' || row.featureType === String(/** @type {LabelFeature} */ (feature).type ?? ''))
  ))
)).length;

// A Label TSV row selects labels in every Result of a batch. The label owner
// binds the displayed Result's labels as the editor does; the committed
// content of another Result carries the renderer's label bindings (B6).
/**
 * The readers the labels a row selects come from: the app state and the label owner's ports.
 * @typedef {Pick<LabelTableActionsOptions, 'state' | 'labelTargetsOf' | 'labelSourceText'>} LabelTargetReaders
 */
/**
 * @param {LabelTargetReaders} readers
 * @param {DrawingState} drawing @param {SVGSVGElement} svg
 * @returns {LabelImportTarget[]}
 */
const labelImportTargets = ({ state, labelTargetsOf, labelSourceText }, drawing, svg) => {
  const targets = labelTargetsOf(drawing, svg);
  const seen = new Set(targets.map(({ feature }) => featureIdentityKeyOf(feature)).filter(Boolean));
  state.results.value.forEach((/** @type {{ content?: unknown } | null | undefined} */ result, /** @type {number} */ index) => {
    if (index === state.selectedResultIndex.value || typeof result?.content !== 'string') return;
    const features = resultRenderedFeatures(state, drawing.mode, index);
    const root = new DOMParser().parseFromString(result.content, 'image/svg+xml').documentElement;
    root.querySelectorAll(`text[${LABEL_FEATURE_ID_ATTRIBUTE}]`).forEach((element) => {
      const featureId = String(element.getAttribute(LABEL_FEATURE_ID_ATTRIBUTE) || '').trim();
      const feature = featureId ? features?.get(featureId) || null : null;
      const key = featureIdentityKeyOf(feature);
      if (!key || seen.has(key)) return;
      seen.add(key);
      targets.push({ featureId, feature, sourceText: labelSourceText(drawing, element, feature) });
    });
  });
  return targets;
};

/** @param {LabelTableActionsOptions} options */
export const createLabelTableActions = ({
  state, ref, computed, isRuleInputCurrent, snapshotRuleInputs, evaluateLabelRules, showNotice,
  labelTargetsOf, labelSourceText, applyLabelTable, rememberLabelImportNotice
}) => {
  const shownSvg = () => state.svgContainer.value?.querySelector('svg') || null;
  const targetReaders = { state, labelTargetsOf, labelSourceText };

  let labelImportRevision = 0;
  const labelImportFailure = ref(/** @type {LabelImportFailure | null} */ (null));
  const keepRetryCopy = createRetryCopyHold();
  const canRetryLabelImportFailure = computed(() => {
    const drawing = state.activeDrawing();
    const failure = labelImportFailure.value;
    return Boolean(failure
      && state.errorLog.value === failure.error
      && failure.revision === labelImportRevision
      && isRuleInputCurrent(failure.snapshot)
      && failure.intent === labelIntentSignature(drawing));
  });
  const retryLabelImportFailure = () => (canRetryLabelImportFailure.value && labelImportFailure.value
    ? labelImportFailure.value.retry() : false);
  const editLabelImportFailure = () => {
    if (canRetryLabelImportFailure.value) labelImportFailure.value?.input.click?.();
  };
  /** @param {{ target: LabelTableInput | null, sourceInput?: LabelTableInput }} event A file input's change event, or the retry's stand-in for it. */
  const loadLabelOverrideTable = async (event) => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    const input = event.target;
    const file = input?.files?.[0];
    if (!input || !file) return;
    const sourceSvg = shownSvg();
    const before = snapshotRuleInputs();
    const labelIntent = labelIntentSignature(drawing);
    const revision = ++labelImportRevision;
    const previousError = state.errorLog.value;
    const endRead = holdFileCopy(file);

    try {
      const text = await readFileText(file);
      if (input.files?.[0] !== file || shownSvg() !== sourceSvg) return;
      if (state.sessionOperationAvailability()) return state.sessionOperationAvailability();
      const rows = parseLabelOverrideTsv(text);
      const labels = sourceSvg ? labelImportTargets(targetReaders, drawing, sourceSvg) : [];
      const evaluation = await evaluateLabelRules({ kind: 'label', rules: rows, features: labelRulePayloads(rows, labels) });
      if (state.sessionOperationAvailability()) return state.sessionOperationAvailability();
      if (revision !== labelImportRevision || input.files?.[0] !== file) return;
      const declined = labelIntent !== labelIntentSignature(drawing) ? 'label-intent'
        : shownSvg() !== sourceSvg ? 'diagram'
          : !isRuleInputCurrent(before) ? 'rule-inputs' : !sourceSvg ? 'no-diagram' : null;
      const undecided = declined ? 0 : undecidedLabelCount(labels, rows);
      const plan = declined || undecided > 0 ? null : planLabelImport(labels, rows, evaluation.winners);
      if (!plan || plan.applicable.length === 0) {
        showNotice(labelImportNotice(rows.length, declined ? { declined } : plan ? { plan } : { undecided }));
        return;
      }
      applyLabelTable(drawing, plan.applicable);
      if (state.errorLog.value === labelImportFailure.value?.error) state.errorLog.value = null;
      labelImportFailure.value = null;
      keepRetryCopy(null);
      rememberLabelImportNotice(drawing, showNotice(labelImportNotice(rows.length, { plan })));
    } catch (error) {
      if (revision !== labelImportRevision || input.files?.[0] !== file || !isRuleInputCurrent(before)
        || state.errorLog.value !== previousError) return;
      // The read, the parse, and the evaluation throw an Error, so the model is set.
      const model = /** @type {UserFacingError} */ (normalizeUserFacingError(error, { operation: 'evaluateRules', stage: 'resource-staging' }));
      state.errorLog.value = model;
      const retryInput = event.sourceInput || input;
      labelImportFailure.value = { error: model, snapshot: before, intent: labelIntent, revision, input: retryInput,
        retry: () => loadLabelOverrideTable({ target: { files: [file], value: '' }, sourceInput: retryInput }) };
      keepRetryCopy(file, { error: model, drawing });
    } finally {
      endRead();
      if (input.files?.[0] === file) input.value = '';
    }
  };

  const downloadLabelOverrideTable = () => {
    const drawing = state.activeDrawing();
    const rows = labelTableExportLines(drawing.canonicalLabelOverrideRows.value, drawing.labelTextBulkOverrides);
    if (rows.length === 0) {
      showNotice({ text: t('There are no label rules to export.') });
      return;
    }
    const selectedIdx = state.selectedResultIndex.value;
    const resultName = selectedIdx >= 0 && selectedIdx < state.results.value.length
      ? String(state.results.value[selectedIdx]?.name || '') : '';
    downloadTextFile(`${makeSafeFilename(resultName, 'gbdraw')}.label_table.tsv`, `${rows.join('\n')}\n`, 'text/tab-separated-values');
  };

  return { loadLabelOverrideTable, downloadLabelOverrideTable, canRetryLabelImportFailure, retryLabelImportFailure, editLabelImportFailure };
};
