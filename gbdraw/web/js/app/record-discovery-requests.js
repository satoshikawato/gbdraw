// @ts-check
// What a drawing's next record discoveries do, and the ends of discovery the
// composition root does not keep: which Linear rows the next discovery may
// expand, which uploads infer record definitions and ask for the record list,
// the drawing a discovery ran for, the reader of one input's records that both
// discovery owners use, the end of a Circular drawing's discovery, and the read
// of a Linear drawing's records before Generate or the record layout uses them.
// The root keeps the expansion of discovered Linear records and the Linear
// record selector's discovery handler: they read the Linear File edit
// (`applyLinearSeqMutation`) and the annotation rebind, which the root creates
// after the selector that takes the handler. The root passes the state, the
// annotation rebind, and the Linear record selector's refresh as ports.
import { inputTypeSpec } from '../services/input-types.js';
import { discoverInputRecords } from '../services/record-discovery.js';
import { WEB_UX_PROFILE } from '../web-ux-profile.js';
/** @import { DrawingState, SessionBusy, WebState } from '../state.js' */
/** @import { DrawingId } from '../services/drawing-view.js' */
/** @import { InputFileRef } from '../services/history-files.js' */
/** @import { CircularRecordsDiscovered } from './run-analysis.js' */
/** @import { RefreshOutcome } from './linear-record-selector.js' */
/** @import { RecordListRequest, RecordSelectionMode } from './record-selection.js' */

const { reactive } = window.Vue;

/**
 * @typedef {object} RecordDiscoveryRequestPorts
 * @property {WebState} state
 */

/**
 * Registers nothing and takes only the state. The root creates it before the
 * Linear record selector, whose discovery handler and record reader read what
 * it returns.
 * @param {RecordDiscoveryRequestPorts} ports
 */
export const createRecordDiscoveryRequests = ({ state }) => {
  /**
   * What a drawing's next record discoveries do, for the rows and input of
   * that drawing only (OV-2148: a Duplicate keeps its row uids, so a uid alone
   * names a row of several drawings).
   * - `expansions`: the Linear rows whose next discovery may expand them, each
   *   with whether the row's own names were entered for the file being read.
   *   They were not when the file replaced another: those names belong to it.
   * - `uploads`, `circularUpload`: the Linear Files and the Circular input just
   *   uploaded: their first discovery infers record definitions (Linear
   *   GenBank) and asks for the record list of a file with many records (D-04).
   *   A Session load, a rollback, and a History step never add one.
   * - `listRequests`: those asks, in upload order; the drawing's panel opens them.
   * @typedef {{ expansions: Map<string, boolean>, uploads: Set<string>, circularUpload: InputFileRef | null,
   *   listRequests: RecordListRequest[] }} PendingRecordDiscovery
   */
  /** @type {WeakMap<DrawingState, PendingRecordDiscovery>} */
  const pendingDiscoveries = new WeakMap();
  /** @param {DrawingState} drawing @returns {PendingRecordDiscovery} */
  const pendingDiscoveryOf = (drawing) => {
    const known = pendingDiscoveries.get(drawing);
    if (known) return known;
    /** @type {PendingRecordDiscovery} */
    const pending = { expansions: new Map(), uploads: new Set(), circularUpload: null, listRequests: reactive([]) };
    pendingDiscoveries.set(drawing, pending);
    return pending;
  };
  /**
   * @param {DrawingState} drawing
   * @param {RecordSelectionMode} listMode @param {string} sourceKey @param {number} recordCount
   */
  const requestRecordListAfterUpload = (drawing, listMode, sourceKey, recordCount) => {
    if (recordCount > WEB_UX_PROFILE.recordList.autoOpenAbove) {
      pendingDiscoveryOf(drawing).listRequests.push({ mode: listMode, sourceKey });
    }
  };
  /**
   * The drawing a discovery ran for, or null when it was deleted while the
   * discovery ran: such a discovery changes nothing (OV-2148).
   * @param {DrawingId} drawingId @returns {DrawingState | null}
   */
  const listedDrawingOf = (drawingId) => (
    state.drawingTabs.value.some((tab) => tab.id === drawingId) ? state.drawingOf(drawingId) : null
  );
  /**
   * The records of one input: a GenBank file, a GFF3 file with its FASTA, or a FASTA alone.
   * The bindings name the files by ref; their bytes come from the store.
   * @param {string} inputType
   * @param {InputFileRef | null} primaryFile
   * @param {InputFileRef | null} pairedFile
   */
  const readInputRecords = (inputType, primaryFile, pairedFile) => discoverInputRecords(
    inputTypeSpec(inputType), state.inputFiles.fileOf(primaryFile), state.inputFiles.fileOf(pairedFile)
  );
  return { pendingDiscoveryOf, requestRecordListAfterUpload, listedDrawingOf, readInputRecords };
};

/**
 * @typedef {object} RecordDiscoveryRunPorts
 * @property {(drawing: DrawingState) => { expansions: Map<string, boolean>, circularUpload: InputFileRef | null }} pendingDiscoveryOf
 *   What the drawing's next record discoveries do (`createRecordDiscoveryRequests`).
 * @property {(drawing: DrawingState, listMode: RecordSelectionMode, sourceKey: string, recordCount: number) => void} requestRecordListAfterUpload
 *   Asks for the record list of an uploaded file with many records.
 * @property {(drawingId: DrawingId) => DrawingState | null} listedDrawingOf
 *   The drawing a discovery ran for, or null when it was deleted.
 * @property {() => void} rebindAnnotationRecords Rebinds the shown drawing's annotations to its record catalog (R10).
 * @property {(drawing: DrawingState) => SessionBusy | Promise<RefreshOutcome>} refreshLinearRecords
 *   Reads the records of the drawing's Linear cards (the Linear record selector's refresh).
 */

/**
 * Registers nothing. The root creates it once the annotation rebind and the
 * Linear record selector exist, before Generate's runner and the Linear
 * comparison plan writer take what it returns.
 * @param {RecordDiscoveryRunPorts} ports
 */
export const createRecordDiscoveryRuns = ({
  pendingDiscoveryOf, requestRecordListAfterUpload, listedDrawingOf, rebindAnnotationRecords, refreshLinearRecords
}) => {
  /**
   * A Circular drawing's discovery settled (success or error): the first after
   * its input's upload asks for the record list.
   * @param {CircularRecordsDiscovered} discovered
   */
  const handleCircularRecordsDiscovered = ({ drawingId, primaryFile, pairedFile, recordCount }) => {
    rebindAnnotationRecords();
    const drawing = listedDrawingOf(drawingId);
    if (!drawing) return;
    const pending = pendingDiscoveryOf(drawing);
    const upload = pending.circularUpload;
    if (!upload || (upload !== primaryFile && upload !== pairedFile)) return;
    pending.circularUpload = null;
    requestRecordListAfterUpload(drawing, 'circular', 'circular', recordCount);
  };
  /** @param {DrawingState} drawing */
  const materializeAutomaticLinearRecords = async (drawing) => {
    if (drawing.mode !== 'linear') return;
    // A row that already waits for its discovery keeps what its writer said
    // about its names (a replacing file does not own the row's record name).
    const { expansions } = pendingDiscoveryOf(drawing);
    drawing.linearSeqs.forEach((seq) => {
      if (!seq.region_record_id && !expansions.has(seq.uid)) expansions.set(seq.uid, true);
    });
    await refreshLinearRecords(drawing);
  };
  return { handleCircularRecordsDiscovered, materializeAutomaticLinearRecords };
};
