// @ts-check
// The "Labels do not fit" dialog (LP-D43, LP-D45 SHOW-ALL, LP-D48
// WARN-NOTICE): the rows of Python's label fit report in its order, each with
// the settings it changes and the canvas size, or disabled with its reason. A
// row that leaves crossing leaders says so and draws with the notice (LP-D57).
// A Generate that drew nothing opens it; crossing leaders or labels over the
// definition text open it from the notice bar over the drawn figure. Rows a
// bounded search did not reach are computed here, one after another, by the
// Worker (OD-F14, the font scan bound), each for at most Python's row bound
// (LP-D60). Apply and Cancel first let the row being computed finish, which
// starts no further row. Apply writes the row's settings as one History step
// (PD-OI-088) and runs Generate; Cancel ("Keep as drawn" over a drawn figure)
// writes nothing. A Circular centre definition over a track of a drawn figure
// also offers "Edit Species / Strain…" (PD-OI-097): it closes the dialog as
// Cancel does and shows the field that writes that record's name lines, found
// by the record's drawn index; it applies nothing. Labels that need much more
// room than the map has were not placed (CAPACITY-STOP): the dialog has no rows
// and lists the ways to make room Python names, with Close only; it words the
// room they need from Python's values (`labelCapacityText`, OV-2179).
import { LABEL_CAPACITY_REMEDY_WORDS, labelFitRowOffered, normalizeLabelFitRow } from '../services/label-fit.js';
import { diagnosticError } from '../utils/error-normalization.js';

/**
 * @import {
 *   LabelFitReport, LabelFitRow, LabelFitRowKind, LabelFitRowToken, ConfigPath, LabelFitSetting, LabelFitTrackOptions,
 *   LabelFitValue
 * } from '../services/label-fit.js'
 */

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
 * @property {boolean} capacity True when a record was not placed for lack of room (CAPACITY-STOP).
 * @property {string[]} remedies The ways to make room, as the Web words them (i18n keys), in Python's order.
 * @property {{ recordId: string, recordIndex: number | null, labelKind: string, labelKindTitle: string,
 *   labelCount: number, fontSize: number, counts: string, capacityRatio: number | null }[]} records Python's counts
 *   text per record (OV-699: a refused record says why); a record not placed for lack of room (`capacityRatio`) is
 *   worded by `labelCapacityText` instead.
 * @property {LabelFitDialogRow[]} rows
 * @property {LabelFitRowKind | ''} selectedKind
 * @property {boolean} finishing Apply or Cancel waits for the row being computed.
 * @property {LabelFitDefinitionRecord[]} definitionRecords The drawn records whose centre definition text overlaps
 *   a track; none when nothing was drawn.
 */
/**
 * @typedef {object} LabelFitDefinitionRecord A drawn record whose name lines "Edit Species / Strain…" shows.
 * @property {number} recordIndex The index of its drawn definition group (`data-gbdraw-record-index`), its identity.
 * @property {string} recordId Shown on the button when several records have one.
 */
/**
 * @typedef {object} LabelFitRowSettings What Apply writes of an offered row.
 * @property {Readonly<Record<ConfigPath, LabelFitValue>>} values The config overrides, by config path.
 * @property {Readonly<LabelFitTrackOptions>} trackOptions The request track options, by request name.
 */
/** @typedef {'circular' | 'linear'} LabelFitMode */
/**
 * @typedef {object} LabelFitDialogPorts
 * @property {(value: LabelFitDialogState) => LabelFitDialogState} reactive Vue `reactive`.
 * @property {(rowToken: LabelFitRowToken, kind: LabelFitRowKind) => Promise<unknown>} computeRow
 *   The Worker's `computeLabelFitRow` reply for one row.
 * @property {(settings: LabelFitRowSettings, mode: LabelFitMode) => void} writeRowSettings
 *   Writes a row's settings into the drawing of `mode` (the request mapping's draft fields).
 * @property {(label: () => string, handler: (choice: string) => unknown, cancel: () => unknown) =>
 *   (choice: string) => unknown} withHistory The modal choice's one History step (`createDialogChoice`).
 * @property {(close: () => unknown) => void} closeAfterChoice
 * @property {() => boolean} choicePending
 * @property {() => Promise<unknown>} runAnalysis Generate.
 * @property {(error: unknown) => string} errorSummary The user-facing summary of a failure.
 * @property {(recordIndex: number) => boolean} showRecordText Shows and focuses the field that writes the name
 *   lines of the drawn record with this index (the record-definition popup's owner); false when no field does.
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
  reactive, computeRow, writeRowSettings, withHistory, closeAfterChoice, choicePending, runAnalysis, errorSummary,
  showRecordText
}) => {
  const dialog = reactive({
    show: false, drawn: false, capacity: false, remedies: [], records: [], rows: [], selectedKind: '', finishing: false,
    definitionRecords: []
  });
  // The rows' values stay with the report, out of the reactive display state.
  /** @type {Map<LabelFitRowKind, LabelFitRow>} */
  let rowsByKind = new Map();
  /** @type {LabelFitMode} */
  let dialogMode = 'circular';
  /** @type {LabelFitRowToken | null} */
  let dialogRowToken = null;
  let openCount = 0;
  /** @type {Promise<unknown> | null} The Worker's reply for the row being computed. */
  let runningRow = null;

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
      if (opened !== openCount || dialog.finishing) return;
      if (!rowToken) continue;
      updateRow(row.kind, { status: 'computing' });
      const computing = runningRow = computeRow(rowToken, row.kind);
      try {
        const reply = /** @type {{ row?: unknown } | null} */ (await computing);
        const computed = normalizeLabelFitRow(reply?.row);
        if (opened !== openCount) return;
        if (computed.kind !== row.kind) throw diagnosticError('RESULT_INVALID', {}, { stage: 'result-admission' });
        rowsByKind.set(row.kind, computed);
        updateRow(row.kind, displayRow(computed));
        selectFirstOffered();
      } catch (error) {
        if (opened !== openCount) return;
        updateRow(row.kind, { status: 'failed', reason: errorSummary(error) });
      } finally {
        if (runningRow === computing) runningRow = null;
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
    dialogRowToken = rowToken;
    rowsByKind = new Map(report.rows.map((row) => [row.kind, row]));
    dialog.drawn = report.drawn;
    dialog.capacity = report.records.some((record) => record.capacityRatio !== null);
    dialog.remedies = report.remedies.map((remedy) => LABEL_CAPACITY_REMEDY_WORDS[remedy]);
    dialog.records = report.records.map((record) => ({
      recordId: record.recordId,
      recordIndex: record.recordIndex,
      labelKind: record.labelKind,
      labelKindTitle: record.labelKindTitle,
      labelCount: record.labelCount,
      fontSize: record.fontSize,
      counts: record.counts,
      capacityRatio: record.capacityRatio
    }));
    // Only a drawn figure has the definition to edit: after a Generate that drew nothing the preview may show an
    // older figure. One button per drawn record, by its index (two records may share an ID).
    /** @type {Map<number, LabelFitDefinitionRecord>} */
    const definitionRecords = new Map();
    if (report.drawn) {
      for (const { labelKind, recordIndex, recordId } of report.records) {
        if (labelKind === 'circular_definition' && recordIndex !== null && !definitionRecords.has(recordIndex)) {
          definitionRecords.set(recordIndex, { recordIndex, recordId });
        }
      }
    }
    dialog.definitionRecords = [...definitionRecords.values()];
    dialog.rows = report.rows.map(displayRow);
    dialog.selectedKind = '';
    dialog.finishing = false;
    selectFirstOffered();
    dialog.show = true;
    return computeMissingRows(opened, rowToken);
  };

  const applyRow = withHistory(
    () => `Fit labels: ${rowsByKind.get(/** @type {LabelFitRowKind} */ (dialog.selectedKind))?.title || ''}`,
    () => {
      const row = rowsByKind.get(/** @type {LabelFitRowKind} */ (dialog.selectedKind));
      if (!row || !labelFitRowOffered(row.status) || !row.values) return false;
      writeRowSettings({ values: row.values, trackOptions: row.trackOptions }, dialogMode);
      closeAfterChoice(close);
      return true;
    },
    close
  );

  /**
   * Waits for the row being computed (the Worker runs one call at a time; at
   * most Python's row bound), then applies the choice.
   * @param {'apply' | 'cancel'} choice
   */
  const resolveLabelFitChoice = async (choice) => {
    if (dialog.finishing || choicePending()) return undefined;
    const opened = openCount;
    const running = runningRow;
    if (running) {
      dialog.finishing = true;
      await running.then(() => undefined, () => undefined);
      dialog.finishing = false;
      if (opened !== openCount) return undefined;
    }
    const chosen = await applyRow(choice);
    // A choice History refused leaves the dialog open: the rows the wait stopped go on.
    if (running && opened === openCount) void computeMissingRows(opened, dialogRowToken);
    return choice === 'apply' && chosen === true ? runAnalysis() : undefined;
  };

  /**
   * "Edit Species / Strain…": closes the dialog as Cancel does, then shows the
   * record's name field (no Apply, no Generate).
   * @param {number} recordIndex The record's drawn index (`LabelFitDefinitionRecord`).
   */
  const editLabelFitRecordText = async (recordIndex) => {
    await resolveLabelFitChoice('cancel');
    return dialog.show ? false : showRecordText(recordIndex);
  };

  return { labelFitDialog: dialog, openLabelFitDialog, resolveLabelFitChoice, editLabelFitRecordText };
};
