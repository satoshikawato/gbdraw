// @ts-check
// The Depth tracks of a drawing: its Depth file rows and series, the Depth
// panels of the shown drawing (rows, coverage, labels, colors), and the Depth
// file and series transitions. The pure readers and writers of a drawing are
// module functions; the composition root passes the state, the shown-drawing
// readers and admission, the Legend retirement, the stack editors' Depth
// source functions, the stack question, and History's undoable step as ports.
// A series removal (`removeCircularDepthTrack`, `removeLinearDepthTrack`)
// stays in app-setup.js: it writes the track slots as a listed reconcile (R10).
import { bindDrawingInputFileTree } from '../services/drawing/inputs.js';
import { linearSourceDepthStatus } from '../services/linear-sources.js';
import {
  activeDepthTrackIndices,
  clearDepthTrackSourceAt,
  compactDepthFileSlots,
  depthFileSlotsFromValue,
  depthTrackCoverageCount,
  depthTrackMatrixWidth,
  ensureDepthTrackConfigCount as ensureDepthTrackConfigCountEntries,
  depthTrackConfigForEdit,
  depthTrackEditDefaults,
  depthSeriesDefaultColor,
  isRecordMajorDepthFileMatrix,
  normalizeDepthTrackConfig as normalizeDepthTrackConfigEntry,
  normalizeRecordMajorDepthFileRows,
  padDepthFileSlots,
  representativeDepthFiles,
  defaultDepthSeriesLabel,
  getDepthTrackFileBaseName,
  isDepthSeriesAutoLabel,
  syncDepthSlotLabels,
  uploadedDepthFileCount
} from '../services/depth-track-state.js';
import { msg, t, tn } from '../utils/i18n.js';
import { askTrackSourceRemoval } from './track-slot-edits.js';
/** @import { DiagramMode, DrawingState, LinearSeq, SessionBusy, WebState } from '../state.js' */
/** @import { AskingStackPorts, SwitchedTrackRow } from './track-slot-edits.js' */
/** @import { FileLike } from '../services/session-resource-backing.js' */
/** @import { InputFileCandidate, InputFileRef } from '../services/history-files.js' */
/** @import { LinearSourceGroup } from '../services/linear-sources.js' */
/** @import { DepthFileSlots, DepthSourceChange, DepthSourceLoss, RawDepthTrackConfig, SyncDepthSlotLabelsOptions } from '../services/depth-track-state.js' */
/** @import { CircularTrackSlot } from '../services/circular-track-slot-model.js' */
/** @import { LinearTrackSlot } from '../services/linear-track-slot-model.js' */

const { computed } = window.Vue;

/**
 * A mode panel's mutating handler: the session's busy answer, nothing when the
 * shown drawing is of another mode, else `action` on the shown drawing.
 * @typedef {<T>(drawingMode: DiagramMode, action: (drawing: DrawingState) => T) => T | SessionBusy | undefined} ShownAction
 */
/** @typedef {(drawingMode: DiagramMode) => DrawingState | null} ShownDrawingOfMode The shown drawing when it is of that mode. */
/**
 * The Legend retirement a Depth transition runs in: the Legend styles and names follow the captions the change leaves.
 * @typedef {<T>(drawing: DrawingState, change: () => T) => T} RetireLegendStyles
 */
/**
 * A stack editor's Depth source change: `mutate` writes the Depth files, and the managed Depth rows follow.
 * @typedef {(drawing: DrawingState, mutate: () => unknown, change?: DepthSourceChange) => unknown} ChangeDepthSources
 */

/** @param {DepthFileSlots | DepthFileSlots[]} value */
const depthFileCount = (value) => (
  isRecordMajorDepthFileMatrix(value)
    ? representativeDepthFiles(value).filter(Boolean).length
    : uploadedDepthFileCount(value)
);
/** @param {DepthFileSlots | DepthFileSlots[]} value */
export const depthTrackCountLabel = (value) => {
  const count = depthFileCount(value);
  return tn(count, '{count} TSV', '{count} TSVs');
};
/** @param {DrawingState} drawing */
const circularDrawingHasDepthFiles = (drawing) => depthFileCount(drawing.files.c_depth) > 0;
/** @param {DrawingState} drawing */
const linearDrawingHasDepthFiles = (drawing) => drawing.linearSeqs.some((seq) => depthFileCount(seq.depth) > 0);
/** @param {LinearSeq | null | undefined} seq */
export const hasLinearDepthFiles = (seq) => depthFileCount(seq?.depth) > 0;
/** @param {DrawingState} drawing */
const circularDepthRecordCount = (drawing) => {
  const discoveredCount = drawing.circularRecordList.value.length;
  if (discoveredCount > 0) return discoveredCount;
  const depth = drawing.files.c_depth;
  return isRecordMajorDepthFileMatrix(depth) && Array.isArray(depth)
    ? Math.max(1, depth.length)
    : 1;
};
/** @param {DrawingState} drawing */
export const circularDepthRows = (drawing) => normalizeRecordMajorDepthFileRows(
  drawing.files.c_depth,
  circularDepthRecordCount(drawing)
);
/** @param {DrawingState} drawing */
export const circularDepthRepresentatives = (drawing) => representativeDepthFiles(circularDepthRows(drawing));
/** @param {DepthFileSlots | DepthFileSlots[]} slots @param {number} [uiCount] */
const sourceDepthTrackCount = (slots, uiCount = 1) => Math.max(
  1,
  Number(uiCount) || 1,
  representativeDepthFiles(slots).length
);
// A Circular drawing's Depth panel count: its stored series (`adv.depth_tracks`,
// which Add Depth track and a Load grow) or its Depth files, whichever is more.
// Each drawing has its own, and Undo and Reset restore it with the series (OV-2200).
/** @param {DrawingState} drawing */
const circularDepthTrackCount = (drawing) => sourceDepthTrackCount(drawing.files.c_depth, drawing.adv.depth_tracks.length);
// The Depth panels read the drawing's series and never write them (OV-109):
// the Depth transitions and the Session reconcile grow, trim and label the
// list, and a series it does not hold yet shows its defaults.
/** @param {DrawingState} drawing @param {number} count */
const rowsForDepthTrackCount = (drawing, count) => {
  const normalizedCount = Math.max(1, Number(count) || 1);
  return Array.from({ length: normalizedCount }, (_, index) => ({
    index,
    key: `depth-track-${index}`,
    config: drawing.adv.depth_tracks[index] || normalizeDepthTrackConfig(drawing, null, index)
  }));
};
// A Linear drawing's Depth files, row by row (each record row holds its own).
/** @param {DrawingState} drawing */
export const linearDepthRows = (drawing) => drawing.linearSeqs.map((seq) => depthFileSlotsFromValue(seq.depth));
// The Depth file rows a removal reads, by identity: a Depth change writes new ones.
/** @param {DrawingState} drawing */
export const circularDepthSources = (drawing) => [drawing.files.c_depth];
/** @param {DrawingState} drawing */
export const linearDepthSources = (drawing) => drawing.linearSeqs.map((seq) => seq.depth);
/** @param {DrawingState} drawing */
export const linearDepthLogicalWidth = (drawing) => depthTrackMatrixWidth(linearDepthRows(drawing));
/** @param {DrawingState} drawing */
const linearDepthTrackUiCount = (drawing) => Math.max(1, linearDepthLogicalWidth(drawing));
/** @param {DrawingState} drawing @param {number} width */
const padLinearDepthRows = (drawing, width) => {
  const targetWidth = Math.max(0, Number(width) || 0);
  drawing.linearSeqs.forEach((seq) => {
    seq.depth = padDepthFileSlots(seq.depth, targetWidth);
  });
};
/** @param {DrawingState} drawing @param {RawDepthTrackConfig | null | undefined} entry @param {number} index */
const normalizeDepthTrackConfig = (drawing, entry, index) => (
  normalizeDepthTrackConfigEntry(entry, index, depthTrackEditDefaults(drawing.adv))
);
// Each drawing holds a series for each of its own mode's Depth sources. A
// new drawing holds its mode default's one series (`createDefaultAdv`,
// OV-1901); a drawing's Reset, which keeps the Depth files, fits the reset
// series to them (`resetDrawing`). Each Depth control fits its drawing
// itself; Load and History restore install the saved or captured series
// (OV-1164).
/** @param {DrawingState} drawing */
export const activeDepthTrackCount = (drawing) => {
  if (drawing.mode === 'linear') {
    return linearDepthTrackUiCount(drawing);
  }
  return circularDepthTrackCount(drawing);
};
/** @param {DrawingState} drawing */
export const ensureDepthTrackConfigCount = (drawing, count = activeDepthTrackCount(drawing)) => {
  const targetCount = Math.max(1, Number(count) || 1);
  const normalized = ensureDepthTrackConfigCountEntries(
    drawing.adv.depth_tracks,
    targetCount,
    depthTrackEditDefaults(drawing.adv)
  );
  drawing.adv.depth_tracks.splice(0, drawing.adv.depth_tracks.length, ...normalized);
};
/** @param {DrawingState} drawing @param {number} trackIndex */
const linearDepthTrackCoverageLabelOf = (drawing, trackIndex) => {
  const covered = depthTrackCoverageCount(linearDepthRows(drawing), trackIndex);
  const total = drawing.linearSeqs.length;
  return tn(total, '{covered}/{count} record', '{covered}/{count} records', { covered });
};
/** @param {CircularTrackSlot | LinearTrackSlot | null | undefined} slot */
export const normalizeDepthSlotTrackIndex = (slot) => {
  const rawTrackIndex = Number(slot?.params?.track_index);
  return Number.isInteger(rawTrackIndex) && rawTrackIndex >= 0 ? rawTrackIndex : 0;
};
// A series' settings for an edit of the shown drawing, or null when the index
// names no series of its mode: the edit is dropped and `depth_tracks` does not
// grow (TK-03).
/** @param {DrawingState} drawing @param {number} index */
const depthTrackConfigForIndex = (drawing, index) => {
  if (!Array.isArray(drawing.adv.depth_tracks)) drawing.adv.depth_tracks = [];
  return depthTrackConfigForEdit(
    drawing.adv.depth_tracks, index, activeDepthTrackCount(drawing), depthTrackEditDefaults(drawing.adv)
  );
};
// A series' saved settings, or its defaults while the drawing holds none (read only, OV-109).
/** @param {DrawingState} drawing @param {number} index */
const readDepthTrackConfig = (drawing, index) => (
  drawing.adv.depth_tracks[index] || normalizeDepthTrackConfig(drawing, null, index)
);
/** @param {DrawingState} drawing */
const depthTrackSlotCollections = (drawing) => [
  Array.isArray(drawing.adv.circular_track_slots) ? drawing.adv.circular_track_slots : [],
  Array.isArray(drawing.adv.linear_track_slots) ? drawing.adv.linear_track_slots : []
];
/** @param {DrawingState} drawing @param {number} index */
export const syncDepthTrackSlotLabelsForTrack = (drawing, index) => {
  depthTrackSlotCollections(drawing).forEach((slots) => {
    syncDepthSlotLabels(/** @type {SyncDepthSlotLabelsOptions} */ ({
      slots,
      depthTracks: drawing.adv.depth_tracks,
      activeCount: drawing.adv.depth_tracks.length
    }));
  });
  void index;
};
// A series label is auto while it is a default name or the name of a file the
// series had in its drawing; the drawing's files decide it (OV-2202).
/** @param {DrawingState} drawing @param {DepthFileSlots} previousFiles @param {DepthFileSlots} nextFiles @param {number} removedIndex */
export const refreshDepthTrackLabelsAfterRemoval = (drawing, previousFiles, nextFiles, removedIndex) => {
  depthFileSlotsFromValue(nextFiles).forEach((file, newIndex) => {
    const oldIndex = newIndex >= removedIndex ? newIndex + 1 : newIndex;
    const config = drawing.adv.depth_tracks[newIndex];
    if (!config) return;
    const oldFile = depthFileSlotsFromValue(previousFiles)[oldIndex] || null;
    const currentLabel = String(config.label ?? '').trim();
    if (isDepthSeriesAutoLabel(currentLabel, oldIndex, [oldFile])) {
      // A file names its series; a series without one keeps no label (OV-1508).
      config.label = getDepthTrackFileBaseName(file);
    }
  });
};
/**
 * @param {DrawingState} drawing
 * @param {number} index
 * @param {InputFileRef | null} file
 * @param {readonly (InputFileRef | null)[]} previousFiles The series' files in the drawing before the change
 */
const updateDepthTrackLabelFromFile = (drawing, index, file, previousFiles) => {
  if (!file) return;
  const config = drawing.adv.depth_tracks[index];
  if (!config) return;
  if (isDepthSeriesAutoLabel(String(config.label ?? '').trim(), index, previousFiles)) {
    config.label = getDepthTrackFileBaseName(file);
    syncDepthTrackSlotLabelsForTrack(drawing, index);
  }
};
/** @param {LinearSeq | null | undefined} seq @param {number} index */
export const getLinearDepthFile = (seq, index) => depthFileSlotsFromValue(seq?.depth)[Number(index)] || null;

/**
 * @typedef {object} DepthTrackPanelsPorts
 * @property {WebState} state
 * @property {ShownDrawingOfMode} shownDrawingOfMode
 */

/**
 * The Depth panels of the shown drawing: whether it has Depth files, its series rows and their coverage, and each
 * series' label, name, and color. They read the drawing's series and never write them (OV-109).
 * @param {DepthTrackPanelsPorts} ports
 */
export const createDepthTrackPanels = ({ state, shownDrawingOfMode }) => {
  const hasCircularDepthFiles = computed(() => circularDrawingHasDepthFiles(state.activeDrawing()));
  const canShowDepthTrack = computed(() => {
    const drawing = state.activeDrawing();
    return drawing.mode === 'linear' ? linearDrawingHasDepthFiles(drawing) : circularDrawingHasDepthFiles(drawing);
  });
  const circularDepthTrackRows = computed(() => {
    const drawing = shownDrawingOfMode('circular');
    if (!drawing) return [];
    return rowsForDepthTrackCount(drawing, circularDepthTrackCount(drawing));
  });
  const linearDepthTrackRows = () => {
    const drawing = shownDrawingOfMode('linear');
    if (!drawing) return [];
    return rowsForDepthTrackCount(drawing, linearDepthTrackUiCount(drawing));
  };
  /** @param {LinearSourceGroup} source */
  const linearSourceDepthRows = (source) => linearDepthTrackRows().map((track) => ({
    ...track,
    status: linearSourceDepthStatus(source, track.index)
  }));
  /** @param {LinearSourceGroup} source */
  const linearSourceDepthSummary = (source) => {
    const tracks = linearSourceDepthRows(source);
    const attachedCount = tracks.filter((track) => track.status.selectedCount > 0).length;
    if (attachedCount === 0) return t('No depth track attached');
    return tn(attachedCount, '{count} depth track attached', '{count} depth tracks attached');
  };
  const depthTrackRows = computed(() => {
    const drawing = state.activeDrawing();
    return rowsForDepthTrackCount(drawing, activeDepthTrackCount(drawing));
  });
  // The shown Linear drawing's Depth series, for its panel.
  /** @param {number} trackIndex */
  const linearDepthTrackCoverageLabel = (trackIndex) => {
    const drawing = shownDrawingOfMode('linear');
    return drawing ? linearDepthTrackCoverageLabelOf(drawing, trackIndex) : '';
  };
  const linearDepthTrackIndexOptions = () => {
    const drawing = shownDrawingOfMode('linear');
    if (!drawing) return [];
    const rows = linearDepthRows(drawing);
    const active = new Set(activeDepthTrackIndices(rows));
    return Array.from({ length: linearDepthLogicalWidth(drawing) }, (_, trackIndex) => ({
      trackIndex,
      label: depthSeriesName(trackIndex),
      coverage: linearDepthTrackCoverageLabelOf(drawing, trackIndex),
      disabled: !active.has(trackIndex)
    }));
  };
  /** @param {number} index */
  const getDepthTrackLabel = (index) => {
    const config = readDepthTrackConfig(state.activeDrawing(), Math.max(0, Number(index) || 0));
    return String(config?.label ?? '');
  };
  // The name the figure draws for Depth series `index` of the shown drawing
  // without a label, among the series the request draws (OV-1508). Figure
  // text: English in every UI language.
  /** @param {number} index */
  const depthSeriesDefaultName = (index) => defaultDepthSeriesLabel(index, activeDepthTrackCount(state.activeDrawing()));
  // The name of Depth series `index`: its label, else the default name.
  /** @param {number} index */
  const depthSeriesName = (index) => getDepthTrackLabel(index) || depthSeriesDefaultName(index);
  /** @param {number} index */
  const getDepthTrackColor = (index) => {
    const idx = Math.max(0, Number(index) || 0);
    const config = readDepthTrackConfig(state.activeDrawing(), idx);
    return String(config?.color || depthSeriesDefaultColor(idx, state.activeDrawing().adv));
  };
  /** @param {CircularTrackSlot | LinearTrackSlot | null | undefined} slot */
  const getDepthTrackLegendLabelForSlot = (slot) => (
    getDepthTrackLabel(normalizeDepthSlotTrackIndex(slot))
  );
  /** @param {number} index */
  const getCircularDepthFile = (index) => circularDepthRepresentatives(state.activeDrawing())[Number(index)] || null;
  return {
    hasCircularDepthFiles, canShowDepthTrack, circularDepthTrackRows, linearDepthTrackRows, linearSourceDepthRows,
    linearSourceDepthSummary, depthTrackRows, linearDepthTrackCoverageLabel, linearDepthTrackIndexOptions,
    getDepthTrackLabel, depthSeriesDefaultName, depthSeriesName, getDepthTrackColor, getDepthTrackLegendLabelForSlot,
    getCircularDepthFile
  };
};

/**
 * @typedef {object} DepthSeriesEditsPorts
 * @property {WebState} state
 * @property {{ readonly value: DiagramMode }} mode The shown drawing's mode.
 * @property {() => SessionBusy | null} sessionOperationAvailability
 * @property {RetireLegendStyles} retireLegendStylesOfUnnamedCaptions
 */

/**
 * The series edits of the shown drawing: a series' color and label, a Depth row's Legend label, and a Depth row's
 * label after its series changed.
 * @param {DepthSeriesEditsPorts} ports
 */
export const createDepthSeriesEdits = ({
  state, mode, sessionOperationAvailability, retireLegendStylesOfUnnamedCaptions
}) => {
  /** @param {number} index @param {unknown} value */
  const setDepthTrackColor = (index, value) => {
    const drawing = state.activeDrawing();
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    const idx = Math.max(0, Number(index) || 0);
    const color = String(value ?? '').trim();
    const config = depthTrackConfigForIndex(drawing, idx);
    if (!config) return;
    config.color = color || depthSeriesDefaultColor(idx, drawing.adv);
  };
  /** @param {number} index @param {unknown} value */
  const setDepthTrackLabel = (index, value) => {
    const drawing = state.activeDrawing();
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    const idx = Math.max(0, Number(index) || 0);
    const config = depthTrackConfigForIndex(drawing, idx);
    if (!config) return;
    retireLegendStylesOfUnnamedCaptions(drawing, () => {
      config.label = String(value ?? '');
      syncDepthTrackSlotLabelsForTrack(drawing, idx);
    });
  };
  /** @param {CircularTrackSlot | LinearTrackSlot | null | undefined} slot @param {unknown} value */
  const setDepthTrackLegendLabelForSlot = (slot, value) => {
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    if (!slot) return;
    const idx = normalizeDepthSlotTrackIndex(slot);
    const drawing = state.activeDrawing();
    // A slot whose index names no series has no series label to edit (TK-03).
    if (!depthTrackConfigForIndex(drawing, idx)) return;
    retireLegendStylesOfUnnamedCaptions(drawing, () => {
      slot.params = slot.params && typeof slot.params === 'object' ? { ...slot.params } : {};
      const label = String(value ?? '');
      if (label.trim()) {
        slot.params.legend_label = label;
      } else {
        delete slot.params.legend_label;
      }
      setDepthTrackLabel(idx, label);
    });
  };
  /** @param {LinearTrackSlot | null | undefined} slot */
  const syncDepthTrackSlotLabel = (slot) => {
    const drawing = state.activeDrawing();
    if (!slot || slot.renderer !== 'depth') return;
    const trackIndex = normalizeDepthSlotTrackIndex(slot);
    const hasSource = mode.value === 'linear'
      ? activeDepthTrackIndices(linearDepthRows(drawing)).includes(trackIndex)
      : Boolean(circularDepthRepresentatives(drawing)[trackIndex]);
    if (hasSource) delete slot.depth_binding_error;
    syncDepthTrackSlotLabelsForTrack(drawing, trackIndex);
  };
  return { setDepthTrackColor, setDepthTrackLabel, setDepthTrackLegendLabelForSlot, syncDepthTrackSlotLabel };
};

/**
 * @typedef {object} DepthFilesPorts
 * @property {WebState} state
 * @property {ShownAction} shownAction
 * @property {RetireLegendStyles} retireLegendStylesOfUnnamedCaptions
 * @property {(slotsKey: 'circular_track_slots' | 'linear_track_slots') => AskingStackPorts<SwitchedTrackRow>} askingStackPorts
 *   The ports of the question a removal asks when it leaves track rows that Generate refuses.
 * @property {ChangeDepthSources} changeCircularDepthSources
 * @property {(drawing: DrawingState, nextDepth: DrawingState['files']['c_depth']) => DepthSourceLoss} circularDepthSourceLoss
 *   What a Circular Depth file removal does to the stack.
 * @property {ChangeDepthSources} changeLinearDepthSources
 * @property {(drawing: DrawingState, nextDepth: (DepthFileSlots | undefined)[]) => DepthSourceLoss} linearDepthSourceLoss
 *   What a Linear Depth file removal does to the stack.
 * @property {(label: string, fn: () => unknown) => Promise<unknown>} runUndoable History's undoable step.
 */

/**
 * The Depth file and series transitions of a drawing: a Depth file's replacement in a Circular drawing or in Linear
 * record rows (the shown drawing's controls and the replacement for every drawing write through them), and Add
 * Depth track.
 * @param {DepthFilesPorts} ports
 */
export const createDepthFiles = ({
  state, shownAction, retireLegendStylesOfUnnamedCaptions, askingStackPorts, changeCircularDepthSources,
  circularDepthSourceLoss, changeLinearDepthSources, linearDepthSourceLoss, runUndoable
}) => {
  /**
   * A Circular drawing's Depth file at series `index`: the shown drawing's
   * control and the replacement for every drawing write through it.
   * @param {DrawingState} drawing @param {number} index @param {InputFileCandidate | InputFileRef} file
   */
  const replaceCircularDepthFile = (drawing, index, file) => {
    const idx = Math.max(0, Number(index) || 0);
    const ref = state.inputFiles.registerFile(file);
    const rows = circularDepthRows(drawing);
    const previousFiles = rows.map((row) => row[idx] || null);
    const nextDepth = rows.map((row) => {
      const next = [...row];
      next[idx] = ref;
      return compactDepthFileSlots(next);
    });
    /** @param {() => void} [turnOffReaching] */
    const apply = (turnOffReaching) => {
      ensureDepthTrackConfigCount(drawing, idx + 1);
      // A new file can rename an auto-labeled series, so the label follows in the same retirement.
      retireLegendStylesOfUnnamedCaptions(drawing, () => {
        changeCircularDepthSources(drawing, () => {
          bindDrawingInputFileTree(state.inputFiles, drawing.slices.inputs, 'c_depth', nextDepth);
          turnOffReaching?.();
        }, { setsFile: ref !== null, series: 'same' });
        if (ref) updateDepthTrackLabelFromFile(drawing, idx, ref, previousFiles);
      });
    };
    const asked = ref ? null : askTrackSourceRemoval(askingStackPorts('circular_track_slots'), {
      kind: 'depth-file',
      loss: circularDepthSourceLoss(drawing, nextDepth),
      sources: () => circularDepthSources(drawing),
      label: msg('Change uploaded file'),
      apply
    });
    return asked ?? apply();
  };
  /** @param {number} index @param {FileLike | null} file */
  const setCircularDepthFile = (index, file) => shownAction('circular', (drawing) => replaceCircularDepthFile(drawing, index, file));
  /**
   * The Depth file at series `index` of Linear record rows: the shown drawing's
   * controls and the replacement for every drawing write through it. A removal
   * asks first when it leaves track rows that Generate refuses (`askTrackSourceRemoval`).
   * @param {DrawingState} drawing @param {LinearSeq[]} sequences @param {number} index
   * @param {InputFileCandidate | InputFileRef | null} file
   * @param {string} [label] The History step of the removal, when it asks.
   * @param {(apply: () => unknown) => unknown} [run] Runs the removal when nothing asks.
   */
  const setLinearDepthFiles = (drawing, sequences, index, file, label = msg('Change uploaded file'), run = (apply) => apply()) => {
    const targets = Array.from(sequences || []).filter(Boolean);
    if (!targets.length) return;
    const idx = Math.max(0, Number(index) || 0);
    const ref = state.inputFiles.registerFile(file);
    const logicalWidth = Math.max(linearDepthLogicalWidth(drawing), idx + 1);
    const previousFiles = linearDepthRows(drawing).map((row) => row[idx] || null);
    /** @param {() => void} [turnOffReaching] */
    const apply = (turnOffReaching) => {
      padLinearDepthRows(drawing, logicalWidth);
      ensureDepthTrackConfigCount(drawing, logicalWidth);
      retireLegendStylesOfUnnamedCaptions(drawing, () => {
        changeLinearDepthSources(drawing, () => {
          targets.forEach((seq) => {
            const slots = depthFileSlotsFromValue(seq.depth);
            if (ref) {
              slots[idx] = ref;
              seq.depth = slots;
            } else {
              seq.depth = clearDepthTrackSourceAt(slots, idx, logicalWidth);
            }
          });
          turnOffReaching?.();
        }, { setsFile: ref !== null, series: 'same' });
        if (ref) updateDepthTrackLabelFromFile(drawing, idx, ref, previousFiles);
      });
    };
    if (ref) return apply();
    const cleared = new Set(targets.map((seq) => seq.uid));
    const asked = askTrackSourceRemoval(askingStackPorts('linear_track_slots'), {
      kind: 'depth-file',
      loss: linearDepthSourceLoss(drawing, drawing.linearSeqs.map((seq) => (
        cleared.has(seq.uid) ? clearDepthTrackSourceAt(seq.depth, idx, logicalWidth) : depthFileSlotsFromValue(seq.depth)
      ))),
      sources: () => linearDepthSources(drawing),
      label,
      apply
    });
    return asked ?? run(apply);
  };
  /** @param {LinearSeq} seq @param {number} index @param {FileLike | null} file */
  const setLinearDepthFile = (seq, index, file) => shownAction('linear', (drawing) => {
    return setLinearDepthFiles(drawing, [seq], index, file);
  });
  /** @param {LinearSourceGroup} source @param {number} index @param {FileLike | null} file @param {string} [label] @param {(apply: () => unknown) => unknown} [run] */
  const setLinearSourceDepthFile = (source, index, file, label, run) => shownAction('linear', (drawing) => setLinearDepthFiles(
    drawing,
    (source?.records || []).map(({ sequence }) => sequence),
    index,
    file,
    label,
    run
  ));
  /** @param {LinearSourceGroup} source @param {number} index */
  const clearLinearSourceDepthFile = (source, index) => setLinearSourceDepthFile(
    source,
    index,
    null,
    msg('Clear File Depth TSV'),
    (apply) => runUndoable(msg('Clear File Depth TSV'), apply)
  );
  const addCircularDepthTrack = () => shownAction('circular', (drawing) => {
    ensureDepthTrackConfigCount(drawing, circularDepthTrackCount(drawing) + 1);
    retireLegendStylesOfUnnamedCaptions(drawing, () => changeCircularDepthSources(drawing, () => {}, { setsFile: false, series: 'added' }));
  });
  const addLinearDepthTrack = () => shownAction('linear', (drawing) => {
    const nextCount = linearDepthTrackUiCount(drawing) + 1;
    padLinearDepthRows(drawing, nextCount);
    ensureDepthTrackConfigCount(drawing, nextCount);
    retireLegendStylesOfUnnamedCaptions(drawing, () => changeLinearDepthSources(drawing, () => {}, { setsFile: false, series: 'added' }));
  });
  return {
    replaceCircularDepthFile, setCircularDepthFile, setLinearDepthFiles, setLinearDepthFile, setLinearSourceDepthFile,
    clearLinearSourceDepthFile, addCircularDepthTrack, addLinearDepthTrack
  };
};
