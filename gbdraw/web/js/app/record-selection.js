// @ts-check
// The owner of which records a drawing draws (record selection, D-01..D-10):
// the Draw this record checkbox of a Linear record card, the record list of a
// file (Choose records…), the question before the last drawn record of a file
// goes OFF (D-06), Delete settings of OFF records (D-08), and the list that
// opens after the upload of a file with many records (D-04). It writes each
// drawing's `recordsOff`; services/record-draw-selection.js answers "is this
// record drawn". The dialogs hold display values only (R3).
/**
 * @import { RecordDrawKey, RecordListEntry, RecordListRow, RecordListSort } from '../services/record-draw-selection.js'
 */
import {
  RECORD_LIST_SORTS, nextRecordsOff, offChangeLeavesSourceEmpty, recordListRows
} from '../services/record-draw-selection.js';

/**
 * @typedef {'linear' | 'circular'} RecordSelectionMode
 * @typedef {object} RecordSource One input file of a mode and its records, in file order.
 * @property {string} key Linear: the File's uid (`linearSourceGroups`); Circular: 'circular'.
 * @property {string} name The file name(s) the dialogs show.
 * @property {RecordListEntry[]} records
 * @typedef {(mode: RecordSelectionMode) => RecordDrawKey[]} RecordsOffPort
 *   The reactive OFF list of the mode's drawing, which this owner rewrites.
 * @typedef {(mode: RecordSelectionMode, sourceKey: string) => RecordSource | null} RecordSourcePort
 *   One source of the mode with its records, or null when it is gone.
 * @typedef {(label: string, change: () => unknown) => Promise<unknown>} RunUndoablePort History's one step.
 * @typedef {(label: () => string, handler: (choice: string) => unknown, cancel: () => unknown)
 *   => (choice: string) => unknown} DialogChoicePort A choice dialog's one History step (`createDialogChoice`).
 * @typedef {(close: () => void) => void} CloseAfterDialogChoicePort Closes once the choice's step ends.
 * @typedef {(mode: RecordSelectionMode) => void} AfterRecordSetChangePort
 *   Invalidates what the drawn record set derives (Linear comparison results, the alignment plan).
 * @typedef {(mode: RecordSelectionMode, sourceKey: string) => unknown} RemoveSourcePort
 *   Removes a whole file inside the open History step (D-06 Remove File).
 * @typedef {(mode: RecordSelectionMode, keys: RecordDrawKey[]) => unknown} DeleteRecordSettingsPort
 *   Deletes the settings and edits of records inside the open History step (D-08).
 * @typedef {object} RecordSelectionOptions
 * @property {<T extends object>(value: T) => T} reactive Vue `reactive`.
 * @property {<T>(getter: () => T) => { readonly value: T }} computed Vue `computed`.
 * @property {RecordsOffPort} recordsOff
 * @property {RecordSourcePort} source
 * @property {RunUndoablePort} runUndoable
 * @property {DialogChoicePort} withDialogChoice
 * @property {CloseAfterDialogChoicePort} closeAfterDialogChoice
 * @property {AfterRecordSetChangePort} afterRecordSetChange
 * @property {RemoveSourcePort} removeSource
 * @property {DeleteRecordSettingsPort} deleteRecordSettings
 * @property {RecordListRequest[]} autoOpenRequests The files whose upload opens their list (D-04), in
 *   upload order; the composition root adds one when a discovery of an upload finds more records than
 *   `WEB_UX_PROFILE.recordList.autoOpenAbove`.
 */

/**
 * @typedef {{ mode: RecordSelectionMode, sourceKey: string }} RecordListRequest
 * @typedef {{ open: boolean, mode: RecordSelectionMode, sourceKey: string, query: string, sort: RecordListSort,
 *   returnFocus: HTMLElement | null }} RecordListDialog
 * @typedef {{ open: boolean, mode: RecordSelectionMode, sourceKey: string, name: string }} RemoveFileDialog
 */

/** @param {RecordSelectionOptions} options */
export const createRecordSelection = ({
  reactive, computed, recordsOff, source, runUndoable, withDialogChoice, closeAfterDialogChoice,
  afterRecordSetChange, removeSource, deleteRecordSettings, autoOpenRequests
}) => {
  const offSets = computed(() => ({
    linear: new Set(recordsOff('linear')),
    circular: new Set(recordsOff('circular'))
  }));
  /** @param {RecordSelectionMode} mode @param {RecordDrawKey} key */
  const isDrawn = (mode, key) => !offSets.value[mode].has(key);
  /** @param {RecordSelectionMode} mode @param {readonly RecordDrawKey[]} keys */
  const drawnCount = (mode, keys) => keys.filter((key) => isDrawn(mode, key)).length;

  /** @type {RecordListDialog} */
  const recordList = reactive({
    open: false, mode: /** @type {RecordSelectionMode} */ ('linear'), sourceKey: '', query: '',
    sort: /** @type {RecordListSort} */ ('file'), returnFocus: /** @type {HTMLElement | null} */ (null)
  });
  /** @type {RemoveFileDialog} */
  const removeFileDialog = reactive({
    open: false, mode: /** @type {RecordSelectionMode} */ ('linear'), sourceKey: '', name: ''
  });
  // The list shown: the one opened with Choose records…, else the first
  // upload that asks for its list (D-04: one list at a time, in upload order).
  const shownList = computed(() => {
    if (recordList.open) return { mode: recordList.mode, sourceKey: recordList.sourceKey, requested: false };
    const request = autoOpenRequests.find((entry) => source(entry.mode, entry.sourceKey));
    return request ? { ...request, requested: true } : null;
  });
  const listSource = computed(() => (shownList.value ? source(shownList.value.mode, shownList.value.sourceKey) : null));
  const listMode = () => shownList.value?.mode || recordList.mode;
  /** The rows the list shows: filtered by record ID and sorted for display (D-10). */
  const recordListView = computed(() => {
    const shown = listSource.value;
    const mode = listMode();
    /** @type {RecordListRow[]} */
    const rows = shown ? recordListRows({
      records: shown.records, recordsOff: recordsOff(mode), query: recordList.query, sort: recordList.sort
    }) : [];
    return {
      open: Boolean(shown),
      mode,
      sourceKey: shownList.value?.sourceKey || '',
      name: shown?.name || '',
      total: shown?.records.length || 0,
      drawn: shown ? drawnCount(mode, shown.records.map((record) => record.key)) : 0,
      rows,
      offWithSettings: rows.filter((row) => !row.drawn && row.hasSettings).map((row) => row.key)
    };
  });

  /**
   * Turns records ON or OFF as one History step. Turning off the last drawn
   * record of a file writes nothing and asks whether to remove the file (D-06).
   * @param {RecordSelectionMode} mode
   * @param {string} sourceKey
   * @param {RecordDrawKey[]} keys
   * @param {boolean} drawn
   * @returns {Promise<unknown>}
   */
  const setRecordsDrawn = async (mode, sourceKey, keys, drawn) => {
    const shown = source(mode, sourceKey);
    if (!shown || keys.length === 0) return false;
    const current = recordsOff(mode);
    if (offChangeLeavesSourceEmpty({
      sourceKeys: shown.records.map((record) => record.key), recordsOff: current, keys, drawn
    })) {
      Object.assign(removeFileDialog, { open: true, mode, sourceKey, name: shown.name });
      return { status: 'question' };
    }
    const next = nextRecordsOff(current, keys, drawn);
    if (next.length === current.length && next.every((key, index) => key === current[index])) return false;
    return runUndoable(drawn ? 'Draw records' : 'Leave out records', () => {
      recordsOff(mode).splice(0, Infinity, ...next);
      afterRecordSetChange(mode);
      return true;
    });
  };
  /**
   * The Draw this record checkbox and a list row checkbox. The control shows
   * the state again after the step (unchanged when the step was refused).
   * @param {RecordSelectionMode} mode
   * @param {string} sourceKey
   * @param {RecordDrawKey} key
   * @param {Event} event
   */
  const toggleRecord = async (mode, sourceKey, key, event) => {
    const control = /** @type {HTMLInputElement | null} */ (event?.target || null);
    const outcome = await setRecordsDrawn(mode, sourceKey, [key], Boolean(control?.checked));
    if (control) control.checked = isDrawn(mode, key);
    return outcome;
  };

  /** @param {boolean} drawn Select all (true) or Select none (false) on the rows shown. */
  const selectShown = (drawn) => setRecordsDrawn(
    listMode(), recordListView.value.sourceKey, recordListView.value.rows.map((row) => row.key), drawn
  );

  /**
   * Delete settings (D-08) of OFF records: their card settings, rotation,
   * feature edits and placements, annotations, and comparison pairs, as one
   * History step. An ON record is never touched.
   * @param {RecordSelectionMode} mode
   * @param {RecordDrawKey[]} keys
   */
  const deleteSettings = (mode, keys) => {
    const off = keys.filter((key) => !isDrawn(mode, key));
    if (off.length === 0) return Promise.resolve(false);
    return runUndoable('Delete record settings', () => deleteRecordSettings(mode, off));
  };
  const deleteShownOffSettings = () => deleteSettings(listMode(), recordListView.value.offWithSettings);

  /**
   * @param {RecordSelectionMode} mode
   * @param {string} sourceKey
   * @param {HTMLElement | null} [returnFocus]
   */
  const openRecordList = (mode, sourceKey, returnFocus = null) => {
    if (!source(mode, sourceKey)) return false;
    Object.assign(recordList, { open: true, mode, sourceKey, query: '', sort: 'file', returnFocus });
    return true;
  };
  const closeRecordList = () => {
    const shown = shownList.value;
    const target = recordList.returnFocus;
    Object.assign(recordList, { open: false, sourceKey: '', query: '', sort: 'file', returnFocus: null });
    if (shown?.requested) {
      // The requests before it named files that are gone.
      const index = autoOpenRequests.findIndex((entry) => entry.mode === shown.mode && entry.sourceKey === shown.sourceKey);
      autoOpenRequests.splice(0, index + 1);
    }
    if (target?.isConnected) target.focus({ preventScroll: true });
  };
  /** @param {string} value */
  const setRecordListQuery = (value) => { recordList.query = String(value ?? ''); };
  /** @param {string} value */
  const setRecordListSort = (value) => {
    if (RECORD_LIST_SORTS.includes(/** @type {RecordListSort} */ (value))) recordList.sort = /** @type {RecordListSort} */ (value);
  };

  const closeRemoveFileDialog = () => { removeFileDialog.open = false; };
  /** D-06: Remove File removes the whole file as one step; Cancel changes nothing. */
  const resolveRemoveFile = withDialogChoice(() => 'Remove File', () => {
    const { mode, sourceKey } = removeFileDialog;
    const removed = removeSource(mode, sourceKey);
    closeAfterDialogChoice(() => {
      closeRemoveFileDialog();
      const shown = shownList.value;
      if (shown && shown.mode === mode && shown.sourceKey === sourceKey) closeRecordList();
    });
    return removed;
  }, closeRemoveFileDialog);

  return {
    isDrawn, drawnCount, toggleRecord, setRecordsDrawn, selectShown, deleteSettings, deleteShownOffSettings,
    recordList, recordListView, openRecordList, closeRecordList, setRecordListQuery, setRecordListSort,
    removeFileDialog, resolveRemoveFile
  };
};
