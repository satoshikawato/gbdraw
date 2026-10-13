// @ts-check
// The "Labels do not fit" dialog (LP-D43, LP-D45 SHOW-ALL, LP-D48
// WARN-NOTICE): the rows of Python's label fit report in its order, each with
// the settings it changes and the canvas size, or disabled with its reason. A
// row that leaves crossing leaders says so and draws with the notice (LP-D57).
// A Generate that drew nothing opens it; crossing leaders or labels over the
// definition text open it from the notice bar over the drawn figure. Rows a
// bounded search did not reach are computed here, one after another, by the
// Worker (OD-F14, the font scan bound). Apply writes the row's settings as one History step
// (PD-OI-088) and runs Generate; Cancel writes nothing.
import { labelFitRowOffered, normalizeLabelFitRow } from '../services/label-fit.js';
import { diagnosticError } from '../utils/error-normalization.js';

/** @import { LabelFitReport, LabelFitRow, LabelFitRowKind, LabelFitRowToken, ConfigPath, LabelFitSetting, LabelFitValue } from '../services/label-fit.js' */

/**
 * @typedef {object} LabelFitDialogRow What a row of the dialog shows.
 * @property {LabelFitRowKind} kind
 * @property {string} title
 * @property {'found' | 'found_with_notice' | 'not_found' | 'not_computed' | 'computing' | 'failed'} status
 * @property {boolean} offered Apply can choose it (`found`, or `found_with_notice`).
 * @property {readonly LabelFitSetting[]} settings
 * @property {{ width: number, height: number } | null} canvasSize
 * @property {string | null} reason
 * @property {string | null} note What the row leaves (crossing leaders, labels over the definition text).
 */
/**
 * @typedef {object} LabelFitDialogState
 * @property {boolean} show
 * @property {boolean} drawn True over a drawn figure (crossing leaders or labels over the definition text).
 * @property {{ recordId: string, labelKindTitle: string, labelCount: number, fontSize: number,
 *   counts: string }[]} records Python's counts text per record (OV-699: a refused record says why).
 * @property {LabelFitDialogRow[]} rows
 * @property {LabelFitRowKind | ''} selectedKind
 */
/** @typedef {'circular' | 'linear'} LabelFitMode */
/**
 * @typedef {object} LabelFitDialogPorts
 * @property {(value: LabelFitDialogState) => LabelFitDialogState} reactive Vue `reactive`.
 * @property {(rowToken: LabelFitRowToken, kind: LabelFitRowKind) => Promise<unknown>} computeRow
 *   The Worker's `computeLabelFitRow` reply for one row.
 * @property {(values: Readonly<Record<ConfigPath, LabelFitValue>>, mode: LabelFitMode) => void} writeRowSettings
 *   Writes a row's settings into the drawing of `mode` (the request mapping's draft fields).
 * @property {(label: () => string, handler: (choice: string) => unknown, cancel: () => unknown) =>
 *   (choice: string) => unknown} withHistory The modal choice's one History step (`createDialogChoice`).
 * @property {(close: () => unknown) => void} closeAfterChoice
 * @property {() => boolean} choicePending
 * @property {() => Promise<unknown>} runAnalysis Generate.
 * @property {(error: unknown) => string} errorSummary The user-facing summary of a failure.
 */

/**
 * @param {LabelFitRow} row
 * @returns {LabelFitDialogRow}
 */
const displayRow = (row) => ({
  kind: row.kind,
  title: row.title,
  status: row.status,
  offered: labelFitRowOffered(row.status),
  settings: row.settings,
  canvasSize: row.canvasSize,
  reason: row.reason,
  note: row.note
});

/** @param {LabelFitDialogPorts} ports */
export const createLabelFitDialog = ({
  reactive, computeRow, writeRowSettings, withHistory, closeAfterChoice, choicePending, runAnalysis, errorSummary
}) => {
  const dialog = reactive({
    show: false, drawn: false, records: [], rows: [], selectedKind: ''
  });
  // The rows' values stay with the report, out of the reactive display state.
  /** @type {Map<LabelFitRowKind, LabelFitRow>} */
  let rowsByKind = new Map();
  /** @type {LabelFitMode} */
  let dialogMode = 'circular';
  let openCount = 0;

  const close = () => { dialog.show = false; openCount += 1; };
  /** @param {LabelFitRowKind} kind @param {Partial<LabelFitDialogRow>} update */
  const updateRow = (kind, update) => {
    const index = dialog.rows.findIndex((row) => row.kind === kind);
    if (index >= 0) dialog.rows[index] = { ...dialog.rows[index], ...update };
  };
  // Python orders the rows that fit before those that draw with the notice.
  const selectFirstOffered = () => {
    const selected = rowsByKind.get(/** @type {LabelFitRowKind} */ (dialog.selectedKind));
    if (selected && labelFitRowOffered(selected.status)) return;
    dialog.selectedKind = dialog.rows.find((row) => row.offered)?.kind || '';
  };

  /**
   * @param {number} opened
   * @param {LabelFitRowToken | null} rowToken
   */
  const computeMissingRows = async (opened, rowToken) => {
    for (const row of [...rowsByKind.values()]) {
      if (row.status !== 'not_computed') continue;
      if (opened !== openCount) return;
      if (!rowToken) continue;
      updateRow(row.kind, { status: 'computing' });
      try {
        const reply = /** @type {{ row?: unknown } | null} */ (await computeRow(rowToken, row.kind));
        const computed = normalizeLabelFitRow(reply?.row);
        if (opened !== openCount) return;
        if (computed.kind !== row.kind) throw diagnosticError('RESULT_INVALID', {}, { stage: 'result-admission' });
        rowsByKind.set(row.kind, computed);
        updateRow(row.kind, displayRow(computed));
        selectFirstOffered();
      } catch (error) {
        if (opened !== openCount) return;
        updateRow(row.kind, { status: 'failed', reason: errorSummary(error) });
      }
    }
  };

  /**
   * Opens the dialog on a report; with a row token, computes the rows the
   * report did not search, one after another.
   * @param {LabelFitReport} report
   * @param {{ rowToken?: LabelFitRowToken | null, mode: LabelFitMode }} options
   */
  const openLabelFitDialog = (report, { rowToken = null, mode }) => {
    openCount += 1;
    const opened = openCount;
    dialogMode = mode;
    rowsByKind = new Map(report.rows.map((row) => [row.kind, row]));
    dialog.drawn = report.drawn;
    dialog.records = report.records.map((record) => ({
      recordId: record.recordId,
      labelKindTitle: record.labelKindTitle,
      labelCount: record.labelCount,
      fontSize: record.fontSize,
      counts: record.counts
    }));
    dialog.rows = report.rows.map(displayRow);
    dialog.selectedKind = '';
    selectFirstOffered();
    dialog.show = true;
    return computeMissingRows(opened, rowToken);
  };

  const applyRow = withHistory(
    () => `Fit labels: ${rowsByKind.get(/** @type {LabelFitRowKind} */ (dialog.selectedKind))?.title || ''}`,
    () => {
      const row = rowsByKind.get(/** @type {LabelFitRowKind} */ (dialog.selectedKind));
      if (!row || !labelFitRowOffered(row.status) || !row.values) return false;
      writeRowSettings(row.values, dialogMode);
      closeAfterChoice(close);
      return true;
    },
    close
  );

  /** @param {'apply' | 'cancel'} choice */
  const resolveLabelFitChoice = async (choice) => {
    if (choice === 'cancel') {
      if (choicePending()) return undefined;
      return applyRow('cancel');
    }
    return await applyRow('apply') === true ? runAnalysis() : undefined;
  };

  return { labelFitDialog: dialog, openLabelFitDialog, resolveLabelFitChoice };
};
