// @ts-check
// The comparison rings of a Circular drawing (similarity rings, kept as
// `circularConservation`): their source files, the series each file draws,
// Show Comparison Rings (`enabled`, a stored setting that follows the ring
// files as Show Depth follows the Depth files), the managed stack rows that
// follow them, and the ring controls of the Circular panel. The composition root passes the state,
// the shown-drawing admission, the stack editor's ring-row functions, the
// LOSAT cache-info write, and History's undoable step as ports.
import {
  bindDrawingInputFileList, bindDrawingInputFileTree, ringGenomeFileLists, setConservationBlastsSource
} from '../services/drawing/inputs.js';
import { readComparisonSequence } from '../services/record-discovery.js';
import {
  applyComparisonSequenceRecordLabel,
  circularRingSources,
  conservationSourceDescriptors,
  defaultConservationSeriesLabel,
  moveConservationSeriesEntry,
  normalizeFileList,
  orderedConservationSources,
  orderedOptionalConservationFiles,
  parseConservationLabelText,
  reconcileConservationSeries,
  withoutConservationSeries
} from '../services/conservation-series.js';
import { shownAfterSourceChange } from '../services/depth-track-state.js';
import { msg, t } from '../utils/i18n.js';
import { askTrackSourceRemoval } from './track-slot-edits.js';
/** @import { DiagramMode, DrawingState, SessionBusy, WebState } from '../state.js' */
/** @import { AskingStackPorts, SwitchedTrackRow } from './track-slot-edits.js' */
/** @import { ConservationSlotEntry } from './circular-track-slots.js' */
/** @import { FileLike } from '../services/session-resource-backing.js' */
/** @import { InputFileRef } from '../services/history-files.js' */
/** @import { ConservationFile, ConservationSourcePlace, RawConservationSeriesEntry } from '../services/conservation-series.js' */
/** @import { LosatCacheInfoEntry } from '../services/config.js' */
/** @import { DepthSourceChange } from '../services/depth-track-state.js' */

const { computed } = window.Vue;

/**
 * A mode panel's mutating handler: the session's busy answer, nothing when the
 * shown drawing is of another mode, else `action` on the shown drawing.
 * @typedef {<T>(drawingMode: DiagramMode, action: (drawing: DrawingState) => T) => T | SessionBusy | undefined} ShownAction
 */
/** @typedef {(drawingMode: DiagramMode) => DrawingState | null} ShownDrawingOfMode The shown drawing when it is of that mode. */
/** @typedef {(next: LosatCacheInfoEntry[]) => void} WriteLosatCacheInfo The write of the project's LOSAT cache-info rows. */
/** @typedef {(drawing: DrawingState) => void} DrawingRingStep A ring transition of one drawing. */
/**
 * The ask before a change of the comparison sources takes out managed ring rows other rows reach.
 * @typedef {(
 *   drawing: DrawingState, sourceFiles: readonly ConservationSourcePlace<ConservationFile>[],
 *   series: readonly RawConservationSeriesEntry[], label: string,
 *   apply: (turnOffReaching: () => void) => unknown
 * ) => unknown} AskRingSourceRemoval
 */
/**
 * The ring source files of a drawing, set by a change of the kind `change` names: the BLAST tables of uploaded
 * rings, else each ring's genome files (`circularRingSources`).
 * @typedef {(drawing: DrawingState, nextFiles: readonly (InputFileRef | null | readonly InputFileRef[])[],
 *   change: DepthSourceChange) => void} SetRingSourceFiles
 */

/** @param {DrawingState} drawing */
const isCircularConservationUploadSource = (drawing) => (
  String(drawing.circularConservation.source || '').trim().toLowerCase() === 'upload'
);
/** @param {DrawingState} drawing */
const isDerivedCircularConservationReplay = (drawing) => (
  !isCircularConservationUploadSource(drawing) &&
  drawing.files.c_conservation_blasts_source === 'losat-cache' &&
  normalizeFileList(drawing.files.c_conservation_blasts).length > 0
);
/**
 * The comparison genome files of each ring row, in series order: the row's
 * LOSAT query files, or for an uploaded BLAST table its optional comparison
 * sequence files (the span export reads them). A row lists them and its file
 * controls change them (`setRingRowFiles`).
 * @param {DrawingState} drawing
 * @returns {InputFileRef[][]}
 */
const ringRowFiles = (drawing) => {
  const { series } = drawing.circularConservation;
  if (isCircularConservationUploadSource(drawing)) {
    /** @type {InputFileRef[][]} */
    const companions = ringGenomeFileLists(drawing.files.c_conservation_sequence_sources);
    return series.map((/** @type {RawConservationSeriesEntry} */ entry, /** @type {number} */ index) => (
      [...(companions[seriesSourceIndex(entry, index)] || [])]
    ));
  }
  /** @type {InputFileRef[][]} */
  const rings = ringGenomeFileLists(drawing.files.c_conservation_fastas);
  // A LOSAT-cache replay keeps each ring's place beside its table, a ring
  // without its files included; otherwise every ring has files.
  const places = isDerivedCircularConservationReplay(drawing) ? rings : rings.filter((files) => files.length > 0);
  return orderedOptionalConservationFiles(places, drawing.circularConservation)
    .slice(0, series.length).map((files) => [...(files || [])]);
};
/**
 * The source index a series entry names (its place among the sources), else its place in the series.
 * @param {RawConservationSeriesEntry | null | undefined} entry @param {number} index
 */
const seriesSourceIndex = (entry, index) => (Number.isInteger(Number(entry?.sourceIndex)) ? Number(entry?.sourceIndex) : index);
/** @param {DrawingState} drawing */
const hasRingSourceFiles = (drawing) => circularRingSources(drawing).length > 0;
/** A ring file added as a new series (Add Seq, an upload). @type {DepthSourceChange} */
const RING_FILE_ADDED = Object.freeze({ setsFile: true, series: 'added' });
/** A series removed with its ring file. @type {DepthSourceChange} */
const RING_FILE_REMOVED = Object.freeze({ setsFile: false, series: 'removed' });
/** A change of the ring source or of its files that adds or removes no series. @type {DepthSourceChange} */
const RING_SOURCE_CHANGED = Object.freeze({ setsFile: false, series: 'same' });
/** A ring file replaced in its place, or added to a ring's genome. @type {DepthSourceChange} */
const RING_FILE_REPLACED = Object.freeze({ setsFile: true, series: 'same' });
/**
 * Whether a ring row's file control may remove one of its `count` files: a
 * LOSAT ring keeps its last (Remove series removes the ring).
 * @param {DrawingState} drawing @param {number} count
 */
const canRemoveRingRowFile = (drawing, count) => (isCircularConservationUploadSource(drawing) ? count > 0 : count > 1);

/**
 * A change of the ring source files (`write` binds them) and the setting that
 * follows it: Show Comparison Rings turns on with the first ring file or a
 * ring file added or set in the shown drawing, and off with the last, by the
 * rule of Show Depth (`shownAfterSourceChange`, OV-2818); the caller then
 * reconciles the rings.
 * @param {WebState} state @param {DrawingState} drawing @param {DepthSourceChange} change @param {() => void} write
 */
const changeRingSourceFiles = (state, drawing, change, write) => {
  const wasSourced = hasRingSourceFiles(drawing);
  write();
  const shown = shownAfterSourceChange({
    wasSourced,
    isSourced: hasRingSourceFiles(drawing),
    inShownDrawing: drawing.id === state.activeDrawing().id,
    change
  });
  if (shown !== null) drawing.circularConservation.enabled = shown;
};

/**
 * The series a change of the ring files leaves (`reconcileCircularConservation`).
 * @param {DrawingState} drawing @param {readonly ConservationSourcePlace<ConservationFile>[]} sourceFiles
 */
const seriesOfRingFiles = (drawing, sourceFiles) => reconcileConservationSeries({
  sourceFiles,
  previousSeries: drawing.circularConservation.series,
  legacyLabels: parseConservationLabelText(drawing.circularConservation.labels)
});

/**
 * @typedef {object} CircularConservationSourcesPorts
 * @property {WebState} state
 * @property {WriteLosatCacheInfo} writeLosatCacheInfo
 * @property {DrawingRingStep} syncCircularConservationSlots The stack editor's managed ring rows follow the series.
 * @property {(drawing: DrawingState, entries: readonly ConservationSlotEntry[]) => Set<string>} circularConservationRowsLeaving
 *   The slot ids of the managed ring rows a change of the sources takes out of the stack.
 * @property {(slotsKey: 'circular_track_slots') => AskingStackPorts<SwitchedTrackRow>} askingStackPorts
 */

/**
 * The rings' source transitions and their reconcile, for each drawing.
 * @param {CircularConservationSourcesPorts} ports
 */
export const createCircularConservationSources = ({
  state, writeLosatCacheInfo, syncCircularConservationSlots, circularConservationRowsLeaving, askingStackPorts
}) => {
  /** @param {DrawingState} drawing */
  const clearDerivedCircularConservationBlasts = (drawing) => {
    if (drawing.files.c_conservation_blasts_source !== 'losat-cache') return;
    bindDrawingInputFileList(state.inputFiles, drawing.slices.inputs, 'c_conservation_blasts', []);
    setConservationBlastsSource(drawing.slices.inputs, null);
  };
  /** @type {SetRingSourceFiles} */
  const setCircularConservationSourceFiles = (drawing, nextFiles, change) => {
    changeRingSourceFiles(state, drawing, change, () => {
      if (isCircularConservationUploadSource(drawing)) {
        bindDrawingInputFileList(state.inputFiles, drawing.slices.inputs, 'c_conservation_blasts',
          /** @type {readonly (InputFileRef | null)[]} */ (nextFiles));
        setConservationBlastsSource(drawing.slices.inputs, null);
      } else {
        clearDerivedCircularConservationBlasts(drawing);
        bindDrawingInputFileTree(state.inputFiles, drawing.slices.inputs, 'c_conservation_fastas', nextFiles);
      }
    });
    writeLosatCacheInfo([]);
    reconcileCircularConservation(drawing);
  };
  /**
   * A Circular drawing's similarity-ring BLAST files, the files kept in their
   * places: the uploaded rings through their source edit, the rows of a LOSAT
   * cache replay in place, each followed by the rings' reconcile.
   * @param {DrawingState} drawing @param {readonly (InputFileRef | null)[]} nextFiles
   */
  const replaceCircularConservationBlasts = (drawing, nextFiles) => {
    if (isCircularConservationUploadSource(drawing)) {
      setCircularConservationSourceFiles(drawing, nextFiles, RING_FILE_REPLACED);
      return;
    }
    changeRingSourceFiles(state, drawing, RING_FILE_REPLACED, () => {
      bindDrawingInputFileList(state.inputFiles, drawing.slices.inputs, 'c_conservation_blasts', nextFiles);
    });
    reconcileCircularConservation(drawing);
  };
  /**
   * A change of the comparison sources asks when a managed ring row another row
   * reaches leaves with its source (`askTrackSourceRemoval`).
   * @param {DrawingState} drawing @param {readonly ConservationSourcePlace<ConservationFile>[]} sourceFiles The ring files it leaves.
   * @param {readonly RawConservationSeriesEntry[]} series The series it leaves. @param {string} label
   * @param {(turnOffReaching: () => void) => unknown} apply
   */
  const askRingSourceRemoval = (drawing, sourceFiles, series, label, apply) => askTrackSourceRemoval(askingStackPorts('circular_track_slots'), {
    kind: 'comparison',
    loss: {
      leaving: circularConservationRowsLeaving(
        drawing, orderedConservationSources(sourceFiles, { ...drawing.circularConservation, series })
      ),
      stayingOn: new Set(),
      emptySeries: () => []
    },
    sources: () => [drawing.files.c_conservation_blasts, drawing.files.c_conservation_fastas, ...drawing.circularConservation.series],
    label,
    apply
  });
  // The comparison rings of a Circular drawing follow their source files: the
  // series and the managed stack rows. Each transition that changes the
  // source, its files, or the legacy series names calls it (R10): the ring
  // controls, a Generate that adopts a LOSAT result, Reset Settings, and
  // Session Load (`afterImport`). It leaves Show Comparison Rings as it is: a
  // change of the files sets it first (`changeRingSourceFiles`). A History
  // restore and the Load rollback install what they captured.
  /** @param {DrawingState} drawing */
  const reconcileCircularConservation = (drawing) => {
    const sourceFiles = circularRingSources(drawing);
    if (isDerivedCircularConservationReplay(drawing)) {
      refreshCircularConservationSlots(drawing);
      return;
    }
    const legacyLabels = parseConservationLabelText(drawing.circularConservation.labels);
    const nextSeries = reconcileConservationSeries({
      sourceFiles,
      previousSeries: drawing.circularConservation.series,
      legacyLabels
    });
    drawing.circularConservation.series.splice(0, drawing.circularConservation.series.length, ...nextSeries);
    refreshCircularConservationSlots(drawing);
  };
  // Reset Settings: the kept comparison files keep a ring row each, with the
  // default settings, and Show Comparison Rings stays off as Reset left it, as
  // Show Depth does for kept Depth files (OV-2818); the LOSAT result the old
  // settings searched goes with those settings (OV-2758).
  /** @param {DrawingState} drawing */
  const rederiveCircularConservation = (drawing) => {
    clearDerivedCircularConservationBlasts(drawing);
    reconcileCircularConservation(drawing);
  };
  // Session Load (`afterImport`): a Session whose Show Comparison Rings is not
  // a setting (`ringSettingStored` false: before Session 51, the drawings
  // layout; their writers turned the rings on by the ring files) shows the
  // rings when it has ring files; a Session 51 keeps the setting it saved.
  /** @param {DrawingState} drawing @param {boolean} ringSettingStored */
  const reconcileLoadedCircularConservation = (drawing, ringSettingStored) => {
    if (!ringSettingStored) drawing.circularConservation.enabled = hasRingSourceFiles(drawing);
    reconcileCircularConservation(drawing);
  };
  // A ring row's name or color, or a record label read for it, shows in its
  // managed stack row.
  /** @param {DrawingState} drawing */
  const refreshCircularConservationSlots = (drawing) => {
    if (drawing.adv.circular_track_slots_enabled === true) {
      syncCircularConservationSlots(drawing);
    }
  };
  return {
    clearDerivedCircularConservationBlasts, setCircularConservationSourceFiles, replaceCircularConservationBlasts,
    askRingSourceRemoval, reconcileCircularConservation, rederiveCircularConservation, reconcileLoadedCircularConservation,
    refreshCircularConservationSlots
  };
};

/**
 * @typedef {object} CircularConservationSeriesPanelPorts
 * @property {WebState} state
 * @property {ShownAction} shownAction
 * @property {ShownDrawingOfMode} shownDrawingOfMode
 * @property {WriteLosatCacheInfo} writeLosatCacheInfo
 * @property {{ value: HTMLInputElement | null }} circularConservationFastaInput The Add Seq file input.
 * @property {DrawingRingStep} reconcileCircularConservation
 * @property {DrawingRingStep} refreshCircularConservationSlots
 * @property {AskRingSourceRemoval} askRingSourceRemoval
 */

/**
 * The ring series panel of the shown Circular drawing: the uploaded ring files, the series rows and their order, and
 * the ring controls' transitions.
 * @param {CircularConservationSeriesPanelPorts} ports
 */
export const createCircularConservationSeriesPanel = ({
  state, shownAction, shownDrawingOfMode, writeLosatCacheInfo, circularConservationFastaInput,
  reconcileCircularConservation, refreshCircularConservationSlots, askRingSourceRemoval
}) => {
  /**
   * The shown Circular drawing's uploaded ring files; it asks first when a
   * managed ring row another row reaches leaves with its file.
   * @param {FileLike | FileLike[] | null | undefined} nextFiles
   */
  const setCircularConservationUploadFiles = (nextFiles) => shownAction('circular', (drawing) => {
    const sourceFiles = normalizeFileList(nextFiles);
    /** @param {() => void} [turnOffReaching] */
    const apply = (turnOffReaching) => {
      changeRingSourceFiles(state, drawing, sourceFiles.length > 0 ? RING_FILE_ADDED : RING_FILE_REMOVED, () => {
        bindDrawingInputFileList(state.inputFiles, drawing.slices.inputs, 'c_conservation_blasts', sourceFiles);
        setConservationBlastsSource(drawing.slices.inputs, null);
      });
      bindDrawingInputFileTree(state.inputFiles, drawing.slices.inputs, 'c_conservation_sequence_sources', []);
      writeLosatCacheInfo([]);
      reconcileCircularConservation(drawing);
      turnOffReaching?.();
    };
    return askRingSourceRemoval(drawing, sourceFiles, seriesOfRingFiles(drawing, sourceFiles), msg('Change uploaded file'), apply) ?? apply();
  });
  const circularConservationSeriesRows = computed(() => {
    const drawing = state.activeDrawing();
    const rowFiles = ringRowFiles(drawing);
    return (Array.isArray(drawing.circularConservation.series) ? drawing.circularConservation.series : []).map((/** @type {RawConservationSeriesEntry} */ entry, /** @type {number} */ index) => ({
      index,
      filename: String(entry?.fileName || `source_${Number(index) + 1}`).trim(),
      sourceLabel: isCircularConservationUploadSource(drawing)
        ? t('BLAST {n}', { n: Number(index) + 1 })
        : t('Comparison {n}', { n: Number(index) + 1 }),
      sourceIndex: seriesSourceIndex(entry, index),
      // The row's genome files (`ringRowFiles`): names, and whether a file may move or go.
      files: (rowFiles[index] || []).map((file) => String(file?.name || '')),
      canRemoveFile: canRemoveRingRowFile(drawing, (rowFiles[index] || []).length),
      defaultLabel: defaultConservationSeriesLabel(
        { name: entry?.fileName },
        Number.isInteger(Number(entry?.sourceIndex)) ? Number(entry.sourceIndex) : index
      )
    }));
  });
  /** @param {number} index @param {number} direction */
  const canMoveCircularConservationSeries = (index, direction) => {
    const drawing = shownDrawingOfMode('circular');
    if (!drawing) return false;
    const idx = Number(index);
    const target = idx + Math.sign(Number(direction));
    return (
      Array.isArray(drawing.circularConservation.series) &&
      Number.isInteger(idx) &&
      idx >= 0 &&
      idx < drawing.circularConservation.series.length &&
      target >= 0 &&
      target < drawing.circularConservation.series.length
    );
  };
  /** @param {number} index @param {number} direction */
  const moveCircularConservationSeries = (index, direction) => shownAction('circular', (drawing) => {
    if (moveConservationSeriesEntry(drawing.circularConservation.series, index, direction)) {
      refreshCircularConservationSlots(drawing);
    }
  });
  const openCircularConservationComparisonFilePicker = () => {
    circularConservationFastaInput.value?.click();
  };
  // The ring controls' transitions: the source radio (it changes the ring
  // files the drawing uses, so Show Comparison Rings follows them), Show
  // Comparison Rings, and a ring row's name or color.
  /** @param {'losat' | 'upload'} source */
  const changeCircularConservationSource = (source) => shownAction('circular', (drawing) => {
    changeRingSourceFiles(state, drawing, RING_SOURCE_CHANGED, () => { drawing.circularConservation.source = source; });
    reconcileCircularConservation(drawing);
  });
  /** @param {boolean} shown */
  const setCircularConservationShown = (shown) => shownAction('circular', (drawing) => {
    drawing.circularConservation.enabled = shown;
    refreshCircularConservationSlots(drawing);
  });
  const refreshShownCircularConservationSlots = () => {
    const drawing = shownDrawingOfMode('circular');
    if (drawing) refreshCircularConservationSlots(drawing);
  };
  return {
    setCircularConservationUploadFiles, circularConservationSeriesRows, canMoveCircularConservationSeries,
    moveCircularConservationSeries, openCircularConservationComparisonFilePicker, changeCircularConservationSource,
    setCircularConservationShown, refreshShownCircularConservationSlots
  };
};

/**
 * @typedef {object} CircularConservationComparisonFilesPorts
 * @property {WebState} state
 * @property {ShownAction} shownAction
 * @property {WriteLosatCacheInfo} writeLosatCacheInfo
 * @property {(label: string, fn: () => unknown) => Promise<unknown>} runUndoable History's undoable step.
 * @property {DrawingRingStep} clearDerivedCircularConservationBlasts
 * @property {SetRingSourceFiles} setCircularConservationSourceFiles
 * @property {DrawingRingStep} reconcileCircularConservation
 * @property {DrawingRingStep} refreshCircularConservationSlots
 * @property {AskRingSourceRemoval} askRingSourceRemoval
 */

/**
 * The file controls of a ring row (D-02): Add file, the arrows, and the remove
 * button change the genome files of that ring, one History step each.
 * @param {Pick<CircularConservationComparisonFilesPorts, 'state' | 'shownAction' | 'writeLosatCacheInfo'
 *   | 'clearDerivedCircularConservationBlasts' | 'reconcileCircularConservation'>} ports
 */
const createRingRowFileEdits = ({
  state, shownAction, writeLosatCacheInfo, clearDerivedCircularConservationBlasts, reconcileCircularConservation
}) => {
  /**
   * The genome files of one ring row (`ringRowFiles`), in their new order. An
   * uploaded table's row writes its comparison sequence files at its source's
   * place. A LOSAT ring writes the rings in series order, each series entry
   * then keyed by its ring's files (so a row keeps its name and color while
   * its files change), and a LOSAT-cache replay's tables go: its rings search
   * again at the next Generate.
   * @param {DrawingState} drawing @param {number} rowIndex @param {InputFileRef[]} nextFiles
   * @param {DepthSourceChange} change
   */
  const setRingRowFiles = (drawing, rowIndex, nextFiles, change) => {
    const rows = ringRowFiles(drawing);
    if (!Number.isInteger(rowIndex) || rowIndex < 0 || rowIndex >= rows.length) return;
    if (isCircularConservationUploadSource(drawing)) {
      /** @type {InputFileRef[][]} */
      const companions = ringGenomeFileLists(drawing.files.c_conservation_sequence_sources);
      const place = seriesSourceIndex(drawing.circularConservation.series[rowIndex], rowIndex);
      while (companions.length <= place) companions.push([]);
      companions[place] = nextFiles;
      bindDrawingInputFileTree(state.inputFiles, drawing.slices.inputs, 'c_conservation_sequence_sources', companions);
      return;
    }
    rows[rowIndex] = nextFiles;
    changeRingSourceFiles(state, drawing, change, () => {
      clearDerivedCircularConservationBlasts(drawing);
      bindDrawingInputFileTree(state.inputFiles, drawing.slices.inputs, 'c_conservation_fastas', rows);
    });
    const descriptors = conservationSourceDescriptors(rows);
    let sourceIndex = 0;
    drawing.circularConservation.series.forEach((/** @type {RawConservationSeriesEntry} */ entry, /** @type {number} */ index) => {
      if (!rows[index]?.length) return;
      const { sourceKey, fileName } = descriptors[sourceIndex];
      Object.assign(entry, { sourceKey, fileName, sourceIndex });
      sourceIndex += 1;
    });
    writeLosatCacheInfo([]);
    reconcileCircularConservation(drawing);
  };
  /**
   * Add file on a ring row: the picked files join the end of the row's genome.
   * @param {number} rowIndex @param {Event & { target: HTMLInputElement }} event
   */
  const addRingRowFiles = (rowIndex, event) => shownAction('circular', (drawing) => {
    const picked = Array.from(event?.target?.files || []).filter(Boolean);
    if (event?.target) event.target.value = '';
    if (picked.length === 0) return;
    const added = picked.map((file) => /** @type {InputFileRef} */ (state.inputFiles.registerFile(file)));
    setRingRowFiles(drawing, Number(rowIndex), [...(ringRowFiles(drawing)[Number(rowIndex)] || []), ...added], RING_FILE_REPLACED);
  });
  /**
   * Move a file of a ring row up (-1) or down (1) in the row's genome.
   * @param {number} rowIndex @param {number} fileIndex @param {number} direction
   */
  const moveRingRowFile = (rowIndex, fileIndex, direction) => shownAction('circular', (drawing) => {
    const files = ringRowFiles(drawing)[Number(rowIndex)] || [];
    const target = Number(fileIndex) + Math.sign(Number(direction));
    if (!Number.isInteger(Number(fileIndex)) || target < 0 || target >= files.length || Number(fileIndex) < 0) return;
    const next = [...files];
    [next[Number(fileIndex)], next[target]] = [next[target], next[Number(fileIndex)]];
    setRingRowFiles(drawing, Number(rowIndex), next, RING_SOURCE_CHANGED);
  });
  /**
   * Remove a file of a ring row. A LOSAT ring keeps at least one file (Remove
   * series removes the ring); an uploaded table's comparison sequence is
   * optional.
   * @param {number} rowIndex @param {number} fileIndex
   */
  const removeRingRowFile = (rowIndex, fileIndex) => shownAction('circular', (drawing) => {
    const files = ringRowFiles(drawing)[Number(rowIndex)] || [];
    if (!canRemoveRingRowFile(drawing, files.length) || !(Number(fileIndex) in files)) return;
    setRingRowFiles(drawing, Number(rowIndex), files.filter((_, index) => index !== Number(fileIndex)), RING_SOURCE_CHANGED);
  });
  return { addRingRowFiles, moveRingRowFile, removeRingRowFile };
};

/**
 * The comparison sequences of the shown Circular drawing: Add Seq (with its record label read), the file controls
 * of a ring row (`createRingRowFileEdits`), and the removal of a series.
 * @param {CircularConservationComparisonFilesPorts} ports
 */
export const createCircularConservationComparisonFiles = ({
  state, shownAction, writeLosatCacheInfo, runUndoable, clearDerivedCircularConservationBlasts,
  setCircularConservationSourceFiles, reconcileCircularConservation, refreshCircularConservationSlots, askRingSourceRemoval
}) => {
  const pendingComparisonRecordLabels = new Set();
  // Also waits for a read that starts meanwhile (another Add Seq).
  const settleComparisonRecordLabels = async () => {
    while (pendingComparisonRecordLabels.size > 0) {
      await Promise.allSettled([...pendingComparisonRecordLabels]);
    }
  };
  /** @param {Event & { target: HTMLInputElement }} event */
  const addCircularConservationComparisonFile = (event) => shownAction('circular', (drawing) => {
    const target = event?.target || null;
    const selectedFile = Array.from(target?.files || []).filter(Boolean)[0] || null;
    if (!selectedFile) return;
    const selectedRef = /** @type {InputFileRef} */ (state.inputFiles.registerFile(selectedFile));
    // D12: the Python ring reader names a GenBank or DDBJ row after it is added
    // (a read failure keeps the file-name default; Generate reports it).
    // Generate waits for these reads before it reads the row labels.
    const pending = readComparisonSequence({ files: [selectedFile] })
      .then(({ recordLabel }) => {
        const applied = applyComparisonSequenceRecordLabel({
          series: drawing.circularConservation.series,
          sourceFiles: circularRingSources(drawing),
          file: selectedRef,
          recordLabel
        });
        refreshCircularConservationSlots(drawing);
        return applied;
      })
      .catch(() => false)
      .finally(() => pendingComparisonRecordLabels.delete(pending));
    pendingComparisonRecordLabels.add(pending);
    // B22 (R11): the input is data-history-managed. This one Add Seq step
    // captures its before now, the row appears at once, and the step commits
    // after the label read, so Undo and Redo restore the label.
    const step = runUndoable(msg('Change uploaded file'), settleComparisonRecordLabels);
    changeRingSourceFiles(state, drawing, RING_FILE_ADDED, () => {
      clearDerivedCircularConservationBlasts(drawing);
      // A new ring of the picked file; Add file on its row adds the genome's other files.
      bindDrawingInputFileTree(state.inputFiles, drawing.slices.inputs, 'c_conservation_fastas',
        [...ringGenomeFileLists(drawing.files.c_conservation_fastas).filter((files) => files.length > 0), [selectedRef]]);
    });
    writeLosatCacheInfo([]);
    reconcileCircularConservation(drawing);
    if (target) target.value = '';
    return step;
  });
  /** @param {number} index */
  const removeCircularConservationSource = (index) => shownAction('circular', (drawing) => {
    const idx = Number(index);
    if (!Number.isInteger(idx) || idx < 0 || idx >= drawing.circularConservation.series.length) return;
    if (isDerivedCircularConservationReplay(drawing)) {
      const orderedBlasts = orderedConservationSources(
        drawing.files.c_conservation_blasts,
        drawing.circularConservation
      ).map((entry) => entry.file);
      const orderedFastas = orderedOptionalConservationFiles(
        drawing.files.c_conservation_fastas,
        drawing.circularConservation
      );
      const orderedSequenceSources = orderedOptionalConservationFiles(
        drawing.files.c_conservation_sequence_sources,
        drawing.circularConservation
      );
      const kept = withoutConservationSeries({
        blasts: orderedBlasts,
        fastas: orderedFastas,
        sequenceSources: orderedSequenceSources,
        series: drawing.circularConservation.series
      }, idx);
      /** @param {() => void} [turnOffReaching] */
      const applyKept = (turnOffReaching) => {
        drawing.circularConservation.series.splice(0, drawing.circularConservation.series.length, ...kept.series);
        drawing.circularConservation.labels = kept.labels;
        changeRingSourceFiles(state, drawing, RING_FILE_REMOVED, () => {
          bindDrawingInputFileList(state.inputFiles, drawing.slices.inputs, 'c_conservation_blasts', kept.blasts);
          bindDrawingInputFileTree(state.inputFiles, drawing.slices.inputs, 'c_conservation_fastas', kept.fastas.map((files) => files || []));
          bindDrawingInputFileTree(
            state.inputFiles, drawing.slices.inputs, 'c_conservation_sequence_sources', kept.sequenceSources.map((files) => files || [])
          );
          setConservationBlastsSource(drawing.slices.inputs, kept.blasts.length > 0 ? 'losat-cache' : null);
        });
        writeLosatCacheInfo([]);
        refreshCircularConservationSlots(drawing);
        turnOffReaching?.();
      };
      return askRingSourceRemoval(drawing, normalizeFileList(kept.blasts), kept.series, msg('Remove series'), applyKept) ?? applyKept();
    }
    const entry = drawing.circularConservation.series[idx];
    const sourceFiles = circularRingSources(drawing);
    const descriptors = conservationSourceDescriptors(sourceFiles);
    let sourceIndex = descriptors.findIndex((descriptor) => descriptor.sourceKey === String(entry?.sourceKey || ''));
    if (sourceIndex < 0) {
      const fileName = String(entry?.fileName || '').trim();
      sourceIndex = descriptors.findIndex((descriptor) => descriptor.fileName === fileName);
    }
    if (sourceIndex < 0 && idx < sourceFiles.length) sourceIndex = idx;
    if (sourceIndex < 0 || sourceIndex >= sourceFiles.length) return;
    const keptFiles = sourceFiles.filter((_, fileIndex) => fileIndex !== sourceIndex);
    /** @param {() => void} [turnOffReaching] */
    const apply = (turnOffReaching) => {
      if (isCircularConservationUploadSource(drawing)) {
        bindDrawingInputFileTree(state.inputFiles, drawing.slices.inputs, 'c_conservation_sequence_sources',
          drawing.files.c_conservation_sequence_sources.filter((_, fileIndex) => fileIndex !== sourceIndex));
      }
      setCircularConservationSourceFiles(drawing, keptFiles, RING_FILE_REMOVED);
      turnOffReaching?.();
    };
    return askRingSourceRemoval(drawing, keptFiles, seriesOfRingFiles(drawing, keptFiles), msg('Remove series'), apply) ?? apply();
  });
  return {
    settleComparisonRecordLabels, addCircularConservationComparisonFile, removeCircularConservationSource,
    ...createRingRowFileEdits({ state, shownAction, writeLosatCacheInfo, clearDerivedCircularConservationBlasts, reconcileCircularConservation })
  };
};
