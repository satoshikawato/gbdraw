// @ts-check
// The comparison rings of a Circular drawing (similarity rings, kept as
// `circularConservation`): their source files, the series each file draws,
// Show Comparison Rings (`enabled`, a stored setting that follows the ring
// files as Show Depth follows the Depth files), the managed stack rows that
// follow them, and the ring controls of the Circular panel. The composition root passes the state,
// the shown-drawing admission, the stack editor's ring-row functions, the
// LOSAT cache-info write, and History's undoable step as ports.
import { bindDrawingInputFileList, bindDrawingInputFileTree, setConservationBlastsSource } from '../services/drawing/inputs.js';
import { readComparisonSequence } from '../services/record-discovery.js';
import {
  applyComparisonSequenceRecordLabel,
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
/** @import { InputFileCandidate, InputFileRef } from '../services/history-files.js' */
/** @import { ConservationFile, RawConservationSeriesEntry } from '../services/conservation-series.js' */
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
 *   drawing: DrawingState, sourceFiles: ConservationFile[], series: readonly RawConservationSeriesEntry[], label: string,
 *   apply: (turnOffReaching: () => void) => unknown
 * ) => unknown} AskRingSourceRemoval
 */
/**
 * The ring source files of a drawing, set by a change of the kind `change` names.
 * @typedef {(drawing: DrawingState, nextFiles: InputFileCandidate | readonly InputFileCandidate[], change: DepthSourceChange) => void} SetRingSourceFiles
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
/** @param {DrawingState} drawing */
const getCircularConservationSourceFiles = (drawing) => (
  isCircularConservationUploadSource(drawing) || isDerivedCircularConservationReplay(drawing)
    ? normalizeFileList(drawing.files.c_conservation_blasts)
    : normalizeFileList(drawing.files.c_conservation_fastas)
);
/** @param {DrawingState} drawing */
const hasRingSourceFiles = (drawing) => getCircularConservationSourceFiles(drawing).length > 0;
/** A ring file added as a new series (Add Seq, an upload). @type {DepthSourceChange} */
const RING_FILE_ADDED = Object.freeze({ setsFile: true, series: 'added' });
/** A series removed with its ring file. @type {DepthSourceChange} */
const RING_FILE_REMOVED = Object.freeze({ setsFile: false, series: 'removed' });
/** A change of the ring source or of its files that adds or removes no series. @type {DepthSourceChange} */
const RING_SOURCE_CHANGED = Object.freeze({ setsFile: false, series: 'same' });
/** A ring file replaced in its place. @type {DepthSourceChange} */
const RING_FILE_REPLACED = Object.freeze({ setsFile: true, series: 'same' });

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
 * @param {DrawingState} drawing @param {ConservationFile[]} sourceFiles
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
        bindDrawingInputFileList(state.inputFiles, drawing.slices.inputs, 'c_conservation_blasts', nextFiles);
        setConservationBlastsSource(drawing.slices.inputs, null);
      } else {
        clearDerivedCircularConservationBlasts(drawing);
        bindDrawingInputFileList(state.inputFiles, drawing.slices.inputs, 'c_conservation_fastas', nextFiles);
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
   * @param {DrawingState} drawing @param {ConservationFile[]} sourceFiles The ring files it leaves.
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
    const sourceFiles = getCircularConservationSourceFiles(drawing);
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
    return (Array.isArray(drawing.circularConservation.series) ? drawing.circularConservation.series : []).map((/** @type {RawConservationSeriesEntry} */ entry, /** @type {number} */ index) => ({
      index,
      filename: String(entry?.fileName || `source_${Number(index) + 1}`).trim(),
      sourceLabel: isCircularConservationUploadSource(drawing)
        ? t('BLAST {n}', { n: Number(index) + 1 })
        : t('Comparison {n}', { n: Number(index) + 1 }),
      sourceIndex: Number.isInteger(Number(entry?.sourceIndex)) ? Number(entry.sourceIndex) : index,
      comparisonSequenceFilename: String(
        drawing.files.c_conservation_sequence_sources[
          Number.isInteger(Number(entry?.sourceIndex)) ? Number(entry.sourceIndex) : index
        ]?.name || ''
      ),
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
 * The comparison sequences of the shown Circular drawing: Add Seq (with its record label read), a companion file,
 * and the removal of a series.
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
    const pending = readComparisonSequence({ file: selectedFile })
      .then(({ recordLabel }) => {
        const applied = applyComparisonSequenceRecordLabel({
          series: drawing.circularConservation.series,
          sourceFiles: getCircularConservationSourceFiles(drawing),
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
      bindDrawingInputFileTree(state.inputFiles, drawing.slices.inputs, 'c_conservation_fastas',
        [...normalizeFileList(drawing.files.c_conservation_fastas), selectedRef]);
    });
    writeLosatCacheInfo([]);
    reconcileCircularConservation(drawing);
    if (target) target.value = '';
    return step;
  });
  /** @param {number} sourceIndex @param {Event & { target: HTMLInputElement }} event */
  const setCircularConservationCompanionFile = (sourceIndex, event) => shownAction('circular', (drawing) => {
    const index = Number(sourceIndex);
    if (!Number.isInteger(index) || index < 0) return;
    const selectedFile = Array.from(event?.target?.files || []).filter(Boolean)[0] || null;
    const next = [...drawing.files.c_conservation_sequence_sources];
    while (next.length <= index) next.push(null);
    next[index] = state.inputFiles.registerFile(selectedFile);
    bindDrawingInputFileTree(state.inputFiles, drawing.slices.inputs, 'c_conservation_sequence_sources', next);
    if (event?.target) event.target.value = '';
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
          bindDrawingInputFileTree(state.inputFiles, drawing.slices.inputs, 'c_conservation_fastas', kept.fastas);
          bindDrawingInputFileTree(state.inputFiles, drawing.slices.inputs, 'c_conservation_sequence_sources', kept.sequenceSources);
          setConservationBlastsSource(drawing.slices.inputs, kept.blasts.length > 0 ? 'losat-cache' : null);
        });
        writeLosatCacheInfo([]);
        refreshCircularConservationSlots(drawing);
        turnOffReaching?.();
      };
      return askRingSourceRemoval(drawing, normalizeFileList(kept.blasts), kept.series, msg('Remove series'), applyKept) ?? applyKept();
    }
    const entry = drawing.circularConservation.series[idx];
    const sourceFiles = getCircularConservationSourceFiles(drawing);
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
    settleComparisonRecordLabels, addCircularConservationComparisonFile, setCircularConservationCompanionFile,
    removeCircularConservationSource
  };
};
