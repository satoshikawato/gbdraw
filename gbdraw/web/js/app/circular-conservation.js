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
import { showChosenFileRefusal, showReadRefusal } from '../services/file-content-cache.js';
import {
  applyComparisonSequenceRecordLabel,
  circularRingSources,
  conservationSourceDescriptors,
  conservationSourcesOfSeries,
  defaultConservationSeriesLabel,
  NO_RING_KEY_CHANGES,
  ringGenomePlaces,
  ringKeyChangesOf,
  moveConservationSeriesEntry,
  normalizeFileList,
  orderedConservationSources,
  parseConservationLabelText,
  readRingColorScheme,
  reassignDrawnRingColors,
  reassignRingColorsOnReorder,
  seriesEntriesOfUnclaimedSources,
  reconcileConservationSeries,
  effectiveRingColor,
  ringOwnColor,
  ringSwatchColor,
  seriesNamingSources,
  thenRingKeyChanges,
  withoutConservationSeries
} from '../services/conservation-series.js';
import { normalizeOptionalHexColor } from '../utils/color-utils.js';
import { shownAfterSourceChange } from '../services/depth-track-state.js';
import { msg, t } from '../utils/i18n.js';
import { askTrackSourceRemoval } from './track-slot-edits.js';
/** @import { DiagramMode, DrawingState, SessionBusy, WebState } from '../state.js' */
/** @import { AskingStackPorts, SwitchedTrackRow } from './track-slot-edits.js' */
/** @import { ConservationSlotEntry } from './circular-track-slots.js' */
/** @import { FileLike } from '../services/session-resource-backing.js' */
/** @import { InputFileRef } from '../services/history-files.js' */
/**
 * @import {
 *   ConservationFile, ConservationSeriesEntry, ConservationSourceDescriptor, ConservationSourcePlace, RawConservationSeriesEntry,
 *   RingKeyChanges
 * }
 *   from '../services/conservation-series.js'
 */
/** @import { EditLosatRows } from '../services/losat-cache.js' */
/** @import { SearchedRingAdoption } from '../services/circular-conservation-keys.js' */
/** @import { DepthSourceChange } from '../services/depth-track-state.js' */

const { computed } = window.Vue;

/**
 * A mode panel's mutating handler: the session's busy answer, nothing when the
 * shown drawing is of another mode, else `action` on the shown drawing.
 * @typedef {<T>(drawingMode: DiagramMode, action: (drawing: DrawingState) => T) => T | SessionBusy | undefined} ShownAction
 */
/** @typedef {(drawingMode: DiagramMode) => DrawingState | null} ShownDrawingOfMode The shown drawing when it is of that mode. */
/** @typedef {(drawing: DrawingState) => void} DrawingRingStep A ring transition of one drawing. */
/**
 * The series a change of the comparison sources leaves and how it moves the rings' keys (`RingKeyChanges`).
 * @typedef {{ series: readonly RawConservationSeriesEntry[], keyChanges: RingKeyChanges }} RingSeriesChange
 */
/**
 * The ask before a change of the comparison sources takes out managed ring rows other rows reach.
 * @typedef {(
 *   drawing: DrawingState, sourceFiles: readonly ConservationSourcePlace<ConservationFile>[],
 *   change: RingSeriesChange, label: string,
 *   apply: (turnOffReaching: () => void) => unknown
 * ) => unknown} AskRingSourceRemoval
 */
/**
 * The ring source files of a drawing, set by a change of the kind `change` names: the BLAST tables of uploaded
 * rings, else each ring's genome files (`circularRingSources`); `keyChanges` are the key changes the transition made
 * before it (`reconcileCircularConservation`).
 * @typedef {(drawing: DrawingState, nextFiles: readonly (InputFileRef | null | readonly InputFileRef[])[],
 *   change: DepthSourceChange, keyChanges?: RingKeyChanges) => void} SetRingSourceFiles
 */
/**
 * A ring transition of one drawing that ends in the rings' reconcile, with the key changes it made before it.
 * @typedef {(drawing: DrawingState, keyChanges?: RingKeyChanges) => void} ReconcileRings
 */
/**
 * The LOSAT rings' genome files set to `lists`, one list per place, each ring kept at its place (Add Seq adds a place
 * after the others, a ring row's file edit changes the files at one).
 * @typedef {(drawing: DrawingState, lists: readonly (readonly InputFileRef[])[], change: DepthSourceChange) => void} SetRingGenomeLists
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
 * The comparison genome files of each ring row, in series order, at the
 * places Generate matches the rows to (`ringGenomePlaces`). A row lists them
 * and its file controls change them (`setRingRowFiles`).
 * @param {DrawingState} drawing
 * @returns {InputFileRef[][]}
 */
const ringRowFiles = (drawing) => {
  const { lists, places } = ringGenomePlaces(drawing);
  return places.map((place) => (place === null ? [] : [...(lists[place] || [])]));
};
/**
 * The source index a series entry names (its place among the sources), else its place in the series.
 * @param {RawConservationSeriesEntry | null | undefined} entry @param {number} index
 */
const seriesSourceIndex = (entry, index) => (Number.isInteger(Number(entry?.sourceIndex)) ? Number(entry?.sourceIndex) : index);
/** @param {DrawingState} drawing */
const hasRingSourceFiles = (drawing) => conservationSourceDescriptors(circularRingSources(drawing)).length > 0;
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
 * After a Generate's adoption of its LOSAT result, Generate's own rings (`ringAdoption`, OV-3195): adopted ring k is the
 * ring it searched at k, so it takes that ring's genetic code, and its own colour or none and its scheme index
 * (PD-OI-110), which the request does not carry (it holds the colour each ring drew); and, while the draft's
 * rings are still those it searched (`builtFromDraft`), the stack row of the draft's ring at that ring's place (a ring's
 * genome files and, in a LOSAT-cache replay, its table sit at its place, `ringGenomePlaces`) moves to adopted ring k's
 * key. Every other ring row of the draft goes, so no row lands on another ring.
 * @param {DrawingState} drawing
 * @param {readonly ConservationSourceDescriptor[]} draftSources The draft's ring sources before the adoption.
 * @param {SearchedRingAdoption} ringAdoption
 * @returns {RingKeyChanges}
 */
const adoptedRingKeyChanges = (drawing, draftSources, ringAdoption) => {
  const { searched, builtFromDraft } = ringAdoption;
  drawing.circularConservation.series.forEach((/** @type {RawConservationSeriesEntry} */ entry, /** @type {number} */ index) => {
    const ring = searched[index];
    if (!ring) return;
    entry.losat_gencode = ring.losat_gencode;
    entry.color = ring.ownColor;
    entry.schemeIndex = ring.schemeIndex;
  });
  const { matched } = conservationSourcesOfSeries(circularRingSources(drawing), drawing.circularConservation);
  return new Map(draftSources.map((source) => {
    const index = builtFromDraft ? searched.findIndex((ring) => ring.placeIndex === source.placeIndex) : -1;
    return /** @type {[string, string]} */ ([source.sourceKey, index >= 0 ? String(matched[index]?.sourceKey || '') : '']);
  }));
};
/**
 * A change of a drawing's ring files takes its LOSAT result rows away; its searches stay for its next Generate, and
 * other drawings' rows are theirs (OV-3196, D-221).
 * @param {EditLosatRows} editLosatRows @param {DrawingState} drawing
 */
const clearRingCacheInfo = (editLosatRows, drawing) => editLosatRows(drawing.id, () => []);

/**
 * A change of the ring sources (`write`) that keeps each ring at its place: each series entry takes the identity of
 * the source at the place the one matcher gave it before the change, and the change's key changes are each row's
 * key before (the matcher's, which the row carries) to that source's key. A LOSAT-cache replay's table sits at its
 * ring's place (`ringGenomePlaces`), so leaving the replay keeps each ring on its genome files.
 * @param {DrawingState} drawing @param {() => void} write
 * @returns {RingKeyChanges}
 */
const keepRingsAtPlaces = (drawing, write) => {
  const { matched: before } = conservationSourcesOfSeries(circularRingSources(drawing), drawing.circularConservation);
  write();
  const after = conservationSourceDescriptors(circularRingSources(drawing));
  const series = seriesNamingSources(drawing.circularConservation.series,
    before.map((source) => (source ? after.find((descriptor) => descriptor.placeIndex === source.placeIndex) : null)));
  drawing.circularConservation.series.splice(0, drawing.circularConservation.series.length, ...series);
  return ringKeyChangesOf(before, series);
};
/**
 * The series a change of the ring files leaves and how the change moves the rings' keys, from the key changes the
 * transition made before (`reconcileCircularConservation`).
 * @param {DrawingState} drawing @param {readonly ConservationSourcePlace<ConservationFile>[]} sourceFiles
 * @param {RawConservationSeriesEntry[]} [series] The series the change starts from: the drawing's, or without a removed row.
 * @param {RingKeyChanges} [earlier]
 * @returns {{ series: ConservationSeriesEntry[], keyChanges: RingKeyChanges }}
 */
const seriesOfRingFiles = (drawing, sourceFiles, series = drawing.circularConservation.series, earlier = NO_RING_KEY_CHANGES) => {
  const reconciled = reconcileConservationSeries({
    sourceFiles,
    previousSeries: series,
    legacyLabels: parseConservationLabelText(drawing.circularConservation.labels)
  });
  return { series: reconciled.series, keyChanges: thenRingKeyChanges(earlier, reconciled.keyChanges) };
};
/**
 * The series a LOSAT-cache replay's reconcile leaves, whose tables stay at their places: each entry names the table the
 * matcher gives it (`seriesNamingSources`), and each table no entry claims gains an entry, labelled by the stored label
 * list where an older Session named it there (`seriesEntriesOfUnclaimedSources`). The rings' keys move only as
 * `keyChanges` says.
 * @param {DrawingState} drawing @param {RingKeyChanges} keyChanges
 * @returns {RingSeriesChange}
 */
const seriesOfReplayTables = (drawing, keyChanges) => {
  const { series, labels } = drawing.circularConservation;
  const { matched, unclaimed } = conservationSourcesOfSeries(circularRingSources(drawing), drawing.circularConservation);
  const named = seriesNamingSources(series, matched);
  const gained = seriesEntriesOfUnclaimedSources(unclaimed, parseConservationLabelText(labels), series.length, named);
  return { series: [...named, ...gained], keyChanges };
};

/**
 * @typedef {object} CircularConservationSourcesPorts
 * @property {WebState} state
 * @property {EditLosatRows} editLosatRows
 * @property {DrawingRingStep} syncCircularConservationSlots The stack editor's managed ring rows follow the series.
 * @property {(drawing: DrawingState, keyChanges: RingKeyChanges) => void} rekeyCircularConservationRows
 *   The stack editor's managed ring rows move to their rings' new keys.
 * @property {(drawing: DrawingState, entries: readonly ConservationSlotEntry[], keyChanges: RingKeyChanges) => Set<string>}
 *   circularConservationRowsLeaving The slot ids of the managed ring rows a change of the sources takes out of the stack.
 * @property {(slotsKey: 'circular_track_slots') => AskingStackPorts<SwitchedTrackRow>} askingStackPorts
 */

/**
 * The rings' source transitions and their reconcile, for each drawing. Each transition keeps the rings' identity
 * invariant (`seriesNamingSources`): each series entry names one source of the drawing's ring sources by its key, file
 * name and place, and each managed stack row carries the key the one matcher gives its entry.
 * @param {CircularConservationSourcesPorts} ports
 */
export const createCircularConservationSources = ({
  state, editLosatRows, syncCircularConservationSlots, rekeyCircularConservationRows, circularConservationRowsLeaving,
  askingStackPorts
}) => {
  /** @param {DrawingState} drawing */
  const clearDerivedCircularConservationBlasts = (drawing) => {
    if (drawing.files.c_conservation_blasts_source !== 'losat-cache') return;
    bindDrawingInputFileList(state.inputFiles, drawing.slices.inputs, 'c_conservation_blasts', []);
    setConservationBlastsSource(drawing.slices.inputs, null);
  };
  /** @type {SetRingGenomeLists} */
  const setRingGenomeLists = (drawing, lists, change) => {
    const keyChanges = keepRingsAtPlaces(drawing, () => changeRingSourceFiles(state, drawing, change, () => {
      clearDerivedCircularConservationBlasts(drawing);
      bindDrawingInputFileTree(state.inputFiles, drawing.slices.inputs, 'c_conservation_fastas', lists);
    }));
    clearRingCacheInfo(editLosatRows, drawing);
    reconcileCircularConservation(drawing, keyChanges);
  };
  /** @type {SetRingSourceFiles} */
  const setCircularConservationSourceFiles = (drawing, nextFiles, change, keyChanges = NO_RING_KEY_CHANGES) => {
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
    clearRingCacheInfo(editLosatRows, drawing);
    reconcileCircularConservation(drawing, keyChanges);
  };
  /**
   * A Circular drawing's similarity-ring BLAST files, the files kept in their
   * places: the uploaded rings through their source edit, the rows of a LOSAT
   * cache replay in place, each followed by the rings' reconcile.
   * @param {DrawingState} drawing @param {readonly (InputFileRef | null)[]} nextFiles
   */
  const replaceCircularConservationBlasts = (drawing, nextFiles) => {
    const uploaded = isCircularConservationUploadSource(drawing);
    const keyChanges = keepRingsAtPlaces(drawing, () => changeRingSourceFiles(state, drawing, RING_FILE_REPLACED, () => {
      bindDrawingInputFileList(state.inputFiles, drawing.slices.inputs, 'c_conservation_blasts', nextFiles);
      if (uploaded) setConservationBlastsSource(drawing.slices.inputs, null);
    }));
    if (uploaded) clearRingCacheInfo(editLosatRows, drawing);
    reconcileCircularConservation(drawing, keyChanges);
  };
  /**
   * A change of the comparison sources asks when a managed ring row another row
   * reaches leaves with its source (`askTrackSourceRemoval`).
   * @param {DrawingState} drawing @param {readonly ConservationSourcePlace<ConservationFile>[]} sourceFiles The ring files it leaves.
   * @param {RingSeriesChange} change The series it leaves and how it moves the rings' keys. @param {string} label
   * @param {(turnOffReaching: () => void) => unknown} apply
   */
  const askRingSourceRemoval = (drawing, sourceFiles, { series, keyChanges }, label, apply) => askTrackSourceRemoval(askingStackPorts('circular_track_slots'), {
    kind: 'comparison',
    loss: {
      leaving: circularConservationRowsLeaving(
        drawing, orderedConservationSources(sourceFiles, { ...drawing.circularConservation, series: [...series] }), keyChanges
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
  // The managed ring rows move to their rings' new keys here, once per
  // transition (`rekeyCircularConservationRows`): `keyChanges` are those the
  // transition made before (a removed ring's row, a ring row's file edit),
  // followed by those of the series the one matcher leaves.
  // It ends each transition with the rings' identity invariant
  // (`seriesNamingSources`). A LOSAT-cache replay keeps its series (its tables
  // stand for the rings Generate searched): each entry takes the identity of
  // the table the matcher gives it, so after a Generate's adoption an entry
  // names its table as its row does, and the rows move only by `keyChanges`.
  /** @type {ReconcileRings} */
  const reconcileCircularConservation = (drawing, keyChanges = NO_RING_KEY_CHANGES) => {
    const reconciled = isDerivedCircularConservationReplay(drawing)
      ? seriesOfReplayTables(drawing, keyChanges)
      : seriesOfRingFiles(drawing, circularRingSources(drawing), drawing.circularConservation.series, keyChanges);
    drawing.circularConservation.series.splice(0, drawing.circularConservation.series.length, ...reconciled.series);
    // The series owns the ring labels. The stored label list of an older Session named its rings past its series: both
    // reconciles gave each of them an entry, so the list is read, and a list left would name the next new ring (OV-3166,
    // OV-3200 in a replay).
    drawing.circularConservation.labels = '';
    rekeyCircularConservationRows(drawing, reconciled.keyChanges);
    refreshCircularConservationSlots(drawing);
  };
  // Reset Settings: the kept comparison files keep a ring row each, with the
  // default settings, and Show Comparison Rings stays off as Reset left it, as
  // Show Depth does for kept Depth files (OV-2818); the LOSAT result the old
  // settings searched goes with those settings (OV-2758).
  /** @param {DrawingState} drawing */
  const rederiveCircularConservation = (drawing) => {
    reconcileCircularConservation(drawing, keepRingsAtPlaces(drawing, () => clearDerivedCircularConservationBlasts(drawing)));
  };
  /**
   * The one gate of the rings' adoption. Only a Generate that searched rings of this drawing (`ringAdoption` with
   * searched rings) asks `adopt` to install them (`adopt` says whether a LOSAT result did), then reconciles the rings
   * with Generate's own rings (`adoptedRingKeyChanges`): each takes its searched ring's genetic code and stack row. Any
   * other adoption (a Linear Generate, a Generate that searched no ring, a candidate of the committed request) leaves
   * the draft's rings unread and as they are.
   * @param {DrawingState} drawing @param {(installRings: boolean) => boolean} adopt
   * @param {SearchedRingAdoption | null} ringAdoption
   */
  const adoptCircularRingResult = (drawing, adopt, ringAdoption) => {
    if (!ringAdoption || ringAdoption.searched.length === 0) return void adopt(false);
    const draftSources = conservationSourceDescriptors(circularRingSources(drawing));
    if (adopt(true)) reconcileCircularConservation(drawing, adoptedRingKeyChanges(drawing, draftSources, ringAdoption));
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
    setRingGenomeLists, setCircularConservationSourceFiles, replaceCircularConservationBlasts, adoptCircularRingResult,
    askRingSourceRemoval, reconcileCircularConservation, rederiveCircularConservation, reconcileLoadedCircularConservation,
    refreshCircularConservationSlots
  };
};

/**
 * @typedef {object} CircularConservationSeriesPanelPorts
 * @property {WebState} state
 * @property {ShownAction} shownAction
 * @property {ShownDrawingOfMode} shownDrawingOfMode
 * @property {EditLosatRows} editLosatRows
 * @property {{ value: HTMLInputElement | null }} circularConservationFastaInput The Add Seq file input.
 * @property {DrawingRingStep} reconcileCircularConservation
 * @property {DrawingRingStep} refreshCircularConservationSlots
 * @property {AskRingSourceRemoval} askRingSourceRemoval
 */

/**
 * The ring series panel of the shown Circular drawing: the uploaded ring files, the series rows and their order, and
 * the ring controls' transitions. It returns the panel's members under the names the template uses.
 * @param {CircularConservationSeriesPanelPorts} ports
 */
export const createCircularConservationSeriesPanel = ({
  state, shownAction, shownDrawingOfMode, editLosatRows, circularConservationFastaInput,
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
      clearRingCacheInfo(editLosatRows, drawing);
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
      ),
      // The colour the ring draws (`effectiveRingColor`; the configured colour under Same color), and whether it is
      // the ring's own (「方式の色に戻す」 shows on those rows only, D-03).
      color: ringSwatchColor(effectiveRingColor(entry, drawing.circularConservation.colorScheme), drawing),
      handSet: ringOwnColor(entry) !== null
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
      if (reassignRingColorsOnReorder(drawing)) refreshCircularConservationSlots(drawing);
    }
  });
  // Ring colors (PD-OI-110): the scheme, the re-assignment on reorder, the one-time re-assignment, a ring's own
  // colour, and 「方式の色に戻す」 (D-03). Each is one History step of its control; the card applies on Generate.
  /** @param {unknown} scheme */
  const setCircularConservationColorScheme = (scheme) => shownAction('circular', (drawing) => {
    drawing.circularConservation.colorScheme = readRingColorScheme(scheme);
    refreshCircularConservationSlots(drawing);
  });
  /** @param {boolean} on */
  const setCircularConservationReassignOnReorder = (on) => shownAction('circular', (drawing) => {
    drawing.circularConservation.reassignOnReorder = on === true;
  });
  const reassignCircularConservationColors = () => shownAction('circular', (drawing) => {
    if (reassignDrawnRingColors(drawing)) refreshCircularConservationSlots(drawing);
  });
  /** @param {number} index @param {unknown} value A colour picked for the ring: its own colour from now on. */
  const setCircularConservationRingColor = (index, value) => shownAction('circular', (drawing) => {
    const entry = drawing.circularConservation.series[Number(index)];
    const color = normalizeOptionalHexColor(value);
    if (!entry || !color) return;
    entry.color = color;
    refreshCircularConservationSlots(drawing);
  });
  /** @param {number} index */
  const useCircularConservationSchemeColor = (index) => shownAction('circular', (drawing) => {
    const entry = drawing.circularConservation.series[Number(index)];
    if (!entry || ringOwnColor(entry) === null) return;
    entry.color = null;
    refreshCircularConservationSlots(drawing);
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
    setCircularConservationShown, refreshCircularConservationSlots: refreshShownCircularConservationSlots,
    setCircularConservationColorScheme, setCircularConservationReassignOnReorder, reassignCircularConservationColors,
    setCircularConservationRingColor, useCircularConservationSchemeColor
  };
};

/**
 * The end of a comparison sequence's failed record label read (D12): its row keeps the file-name default, and
 * Generate reports the failure; a gzip file past the expanded cap (PD-OI-116) is refused now, in the error panel, as
 * the other file fields refuse it. This read is the one owner of that refusal for Add Seq, with or without a copy
 * store.
 * @param {{ value: unknown }} errorLog The page's error panel (`state.errorLog`).
 * @returns {(error: unknown) => false}
 */
const keepRecordLabelDefault = (errorLog) => (error) => {
  showReadRefusal(errorLog, error);
  return false;
};

/**
 * A comparison file chosen with a ring row's Add file, held by the input store. Its read refusal, a gzip file past the
 * expanded cap (PD-OI-116), shows once the store's copy counts its expanded size; with no copy store, its first read
 * refuses it.
 * @param {WebState} state
 * @param {File} file
 */
const registerChosenComparisonFile = (state, file) => {
  const ref = state.inputFiles.registerFile(file);
  void showChosenFileRefusal(state.errorLog, file);
  return ref;
};

/**
 * @typedef {object} CircularConservationComparisonFilesPorts
 * @property {WebState} state
 * @property {ShownAction} shownAction
 * @property {EditLosatRows} editLosatRows
 * @property {(label: string, fn: () => unknown) => Promise<unknown>} runUndoable History's undoable step.
 * @property {SetRingGenomeLists} setRingGenomeLists
 * @property {SetRingSourceFiles} setCircularConservationSourceFiles
 * @property {ReconcileRings} reconcileCircularConservation
 * @property {DrawingRingStep} refreshCircularConservationSlots
 * @property {AskRingSourceRemoval} askRingSourceRemoval
 */

/**
 * The file controls of a ring row (D-02): Add file and the remove button change
 * the genome files of that ring, one History step each. The files keep the
 * order they were added in (D-106): their order changes no hit drawn.
 * @param {Pick<CircularConservationComparisonFilesPorts, 'state' | 'shownAction' | 'setRingGenomeLists'>} ports
 */
const createRingRowFileEdits = ({ state, shownAction, setRingGenomeLists }) => {
  /**
   * The genome files of one ring row (`ringGenomePlaces`) after an edit, written
   * at the row's place. An uploaded table's row writes its comparison sequence
   * files. A LOSAT ring keeps its place (`setRingGenomeLists`), so a row keeps its
   * name, color and stack row while its files change, and a LOSAT-cache replay's
   * tables go: its rings search again at the next Generate.
   * @param {DrawingState} drawing @param {number} rowIndex @param {InputFileRef[]} nextFiles
   * @param {DepthSourceChange} change
   */
  const setRingRowFiles = (drawing, rowIndex, nextFiles, change) => {
    const { key, lists, places } = ringGenomePlaces(drawing);
    const place = places[rowIndex];
    if (place === undefined || place === null) return;
    while (lists.length <= place) lists.push([]);
    lists[place] = nextFiles;
    if (key === 'c_conservation_sequence_sources') {
      bindDrawingInputFileTree(state.inputFiles, drawing.slices.inputs, key, lists);
      return;
    }
    setRingGenomeLists(drawing, lists, change);
  };
  /**
   * Add file on a ring row: the picked files join the end of the row's genome.
   * @param {number} rowIndex @param {Event & { target: HTMLInputElement }} event
   */
  const addRingRowFiles = (rowIndex, event) => shownAction('circular', (drawing) => {
    const picked = Array.from(event?.target?.files || []).filter(Boolean);
    if (event?.target) event.target.value = '';
    if (picked.length === 0) return;
    const added = picked.map((file) => /** @type {InputFileRef} */ (registerChosenComparisonFile(state, file)));
    setRingRowFiles(drawing, Number(rowIndex), [...(ringRowFiles(drawing)[Number(rowIndex)] || []), ...added], RING_FILE_REPLACED);
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
  return { addRingRowFiles, removeRingRowFile };
};

/**
 * The comparison sequences of the shown Circular drawing: Add Seq (with its record label read), the file controls
 * of a ring row (`createRingRowFileEdits`), and the removal of a series.
 * @param {CircularConservationComparisonFilesPorts} ports
 */
export const createCircularConservationComparisonFiles = ({
  state, shownAction, editLosatRows, runUndoable, setRingGenomeLists,
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
    // (`keepRecordLabelDefault` ends a failed read). Generate waits for these reads before it reads the row labels.
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
      .catch(keepRecordLabelDefault(state.errorLog))
      .finally(() => pendingComparisonRecordLabels.delete(pending));
    pendingComparisonRecordLabels.add(pending);
    // B22 (R11): the input is data-history-managed. This one Add Seq step
    // captures its before now, the row appears at once, and the step commits
    // after the label read, so Undo and Redo restore the label.
    const step = runUndoable(msg('Change uploaded file'), settleComparisonRecordLabels);
    // A new ring of the picked file, after the others at their places (an empty one kept, `circularRingSources`);
    // Add file on its row adds the genome's other files.
    setRingGenomeLists(drawing, [...ringGenomeFileLists(drawing.files.c_conservation_fastas), [selectedRef]], RING_FILE_ADDED);
    if (target) target.value = '';
    return step;
  });
  /** @param {number} index */
  const removeCircularConservationSource = (index) => shownAction('circular', (drawing) => {
    const idx = Number(index);
    if (!Number.isInteger(idx) || idx < 0 || idx >= drawing.circularConservation.series.length) return;
    if (isDerivedCircularConservationReplay(drawing)) {
      const replayed = conservationSourcesOfSeries(drawing.files.c_conservation_blasts, drawing.circularConservation);
      const orderedBlasts = orderedConservationSources(
        drawing.files.c_conservation_blasts,
        drawing.circularConservation
      ).map((entry) => entry.file);
      // Each table's genome files sit at its place (`ringGenomePlaces`).
      const { places } = ringGenomePlaces(drawing);
      /** @param {readonly (readonly InputFileRef[])[]} list */
      const inSeriesOrder = (list) => {
        const lists = ringGenomeFileLists(list);
        return places.map((place) => (place === null ? null : lists[place] || null));
      };
      const orderedFastas = inSeriesOrder(drawing.files.c_conservation_fastas);
      const orderedSequenceSources = inSeriesOrder(drawing.files.c_conservation_sequence_sources);
      const kept = withoutConservationSeries({
        blasts: orderedBlasts,
        fastas: orderedFastas,
        sequenceSources: orderedSequenceSources,
        series: drawing.circularConservation.series
      }, idx);
      // Each kept ring's row moves to its table's key among the kept tables; the removed ring's row goes.
      const keptChange = {
        series: kept.series,
        keyChanges: ringKeyChangesOf(replayed.matched, conservationSourcesOfSeries(normalizeFileList(kept.blasts), kept).matched, idx)
      };
      /** @param {() => void} [turnOffReaching] */
      const applyKept = (turnOffReaching) => {
        drawing.circularConservation.series.splice(0, drawing.circularConservation.series.length, ...kept.series);
        // An older Session's label list names rings past its series by their places, which the removal shifts.
        drawing.circularConservation.labels = '';
        changeRingSourceFiles(state, drawing, RING_FILE_REMOVED, () => {
          bindDrawingInputFileList(state.inputFiles, drawing.slices.inputs, 'c_conservation_blasts', kept.blasts);
          bindDrawingInputFileTree(state.inputFiles, drawing.slices.inputs, 'c_conservation_fastas', kept.fastas.map((files) => files || []));
          bindDrawingInputFileTree(
            state.inputFiles, drawing.slices.inputs, 'c_conservation_sequence_sources', kept.sequenceSources.map((files) => files || [])
          );
          setConservationBlastsSource(drawing.slices.inputs, kept.blasts.length > 0 ? 'losat-cache' : null);
        });
        clearRingCacheInfo(editLosatRows, drawing);
        reconcileCircularConservation(drawing, keptChange.keyChanges);
        turnOffReaching?.();
      };
      return askRingSourceRemoval(drawing, normalizeFileList(kept.blasts), keptChange, msg('Remove series'), applyKept) ?? applyKept();
    }
    const sourceFiles = circularRingSources(drawing);
    // The row's source as Generate matches it (`conservationSourcesOfSeries`).
    const source = conservationSourcesOfSeries(sourceFiles, drawing.circularConservation).matched[idx];
    if (!source) return;
    const { placeIndex } = source;
    const keptFiles = sourceFiles.filter((_, fileIndex) => fileIndex !== placeIndex);
    // The removed row's entry leaves before the reconcile, so it takes no other ring's files; the rings after its
    // place move up one place with their files, so their stored places still name them.
    const keptSeries = drawing.circularConservation.series.filter(
      (/** @type {RawConservationSeriesEntry} */ _, /** @type {number} */ entryIndex) => entryIndex !== idx
    ).map((/** @type {RawConservationSeriesEntry} */ entry) => (
      Number(entry.sourceIndex) > placeIndex ? { ...entry, sourceIndex: Number(entry.sourceIndex) - 1 } : entry
    ));
    // The removed ring's stack row carries its source's key, which a kept ring may take: it goes first.
    const removedRow = new Map([[source.sourceKey, '']]);
    /** @param {() => void} [turnOffReaching] */
    const apply = (turnOffReaching) => {
      drawing.circularConservation.series.splice(0, drawing.circularConservation.series.length, ...keptSeries);
      if (isCircularConservationUploadSource(drawing)) {
        bindDrawingInputFileTree(state.inputFiles, drawing.slices.inputs, 'c_conservation_sequence_sources',
          drawing.files.c_conservation_sequence_sources.filter((_, fileIndex) => fileIndex !== placeIndex));
      }
      setCircularConservationSourceFiles(drawing, keptFiles, RING_FILE_REMOVED, removedRow);
      turnOffReaching?.();
    };
    return askRingSourceRemoval(
      drawing, keptFiles, seriesOfRingFiles(drawing, keptFiles, keptSeries, removedRow), msg('Remove series'), apply
    ) ?? apply();
  });
  return {
    settleComparisonRecordLabels, addCircularConservationComparisonFile, removeCircularConservationSource,
    ...createRingRowFileEdits({ state, shownAction, setRingGenomeLists })
  };
};
