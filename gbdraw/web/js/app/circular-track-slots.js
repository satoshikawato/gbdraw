// @ts-check
/** @import { DrawingState } from '../state.js' */
/** @import { ChangeTrackLayout, ChangeTrackVisibility } from './track-slot-edits.js' */
/** @import { AskChoice } from './app-messages.js' */
import {
  CONSERVATION_SLOT_MANAGER,
  isManagedConservationSlot,
  normalizeFileList,
  orderedConservationSources,
  safeConservationSlotId
} from '../services/conservation-series.js';
import {
  isDefaultManagedDepthSlot,
  depthSourceLoss,
  reconcileManagedDepthSlots,
  representativeDepthFiles
} from '../services/depth-track-state.js';
import { MODE_PROFILE_DATA } from '../mode-profiles.generated.js';
import { resolveTrackSlotSkewColorValue } from './track-slot-colors.js';
import {
  findTrackSlotGeometry,
  findTrackSlotGeometryRecord,
  formatPxAuto,
  formatRadiusFactorAuto,
  isManualSlotValue,
  normalizeOptionalText,
  overlayFollowsNote,
  sizeClassForRecordLengths,
  tickAnchorRadiusFactor
} from '../services/track-slot-display.js';
import { askTrackHide, featureSlotEdits, trackLayoutWithVisibility, trackRowsText } from './track-slot-edits.js';
import { draftDrawsDepth } from '../services/track-legend-rows.js';
import { ownEntry } from '../utils/own-entry.js';
import {
  GLOBAL_SUPPRESS_PARAM,
  SUPPRESS_FORM_KEY_BY_RENDERER,
  SUPPRESS_KEY_BY_RENDERER,
  SUPPRESS_RENDERER_BY_KEY,
  applyCircularSuppressControlsToSlots,
  circularTrackSlotHiddenBySuppressForm,
  slotUsesDrawingDinucleotide
} from '../services/circular-track-slot-model.js';
import { listText, msg, t, tm, tn } from '../utils/i18n.js';
import { errorText } from '../utils/error-text.js';
import {
  OVERLAY_ANNOTATION_PARAMS,
  nextCopyId,
  nextNumberedId,
  overlayAnnotationZ,
  paramsKeptOnRendererChange,
  parseOptionalCircularScalar,
  parseOptionalPixel,
  setOverlayAnnotationParam,
  trackAddBlockReason,
  validateCustomTrackPlan
} from '../services/track-slot-validation.js';
import {
  applyCircularGeometryShortcuts,
  applyCircularTrackOrderPlacements,
  clampCircularTrackAxisIndex,
  cloneParams,
  createDefaultCircularTrackSlots,
  DEFAULT_SLOT_IDS,
  DEFAULT_TICK_LABEL_LAYOUT,
  effectiveSlotPlacement,
  enforceSingleOnAxisSlot,
  inferLegacyAxisIndexFromFeature,
  laneDirectionForPreset,
  laneDirectionForSide,
  makeSlot,
  normalizeCircularTrackPreset,
  normalizeCircularTrackSlot,
  normalizeCircularTrackSlots,
  normalizeColorParam,
  normalizeLaneDirection,
  normalizeNt,
  normalizeOptionalPlacement,
  normalizePlacement,
  normalizeSlotSide,
  normalizeTickLabelLayout,
  normalizeTrackIndex,
  resolveCircularTrackFeaturePlacement,
  sideForLaneDirection,
  CIRCULAR_TRACK_RENDERERS as SUPPORTED_RENDERERS,
  syncSlotPlacementFromSide,
  syncSlotsFromAxisIndex,
  tickLabelLayoutFromSides,
  normalizeGapText
} from '../services/circular-track-slot-model.js';
/** @import { CircularTrackSlot } from '../services/circular-track-slot-model.js' */
/** @import { OverlaySlotReference, TrackSlotRow } from '../services/track-slot-validation.js' */
/** @import { FileLike } from '../services/session-resource-backing.js' */
/** @import { DepthFileSlots } from '../services/depth-track-state.js' */

/**
 * What the slot helpers read of the Web state (`state.js`): the Depth and
 * Conservation inputs and the shown mode.
 * @typedef {object} TrackSlotStateView
 * @property {{ value?: string }} [mode]
 * @property {{
 *   c_depth?: DepthFileSlots | DepthFileSlots[],
 *   c_conservation_blasts?: FileLike | FileLike[] | null,
 *   c_conservation_blasts_source?: unknown,
 *   c_conservation_fastas?: FileLike | FileLike[] | null
 * }} [files]
 */

const UI_RENDERERS = SUPPORTED_RENDERERS.filter(
  (renderer) => renderer !== 'sequence_conservation'
);

/** @type {Record<string, string>} */
const RENDERER_LABELS = {
  features: msg('Features'),
  ticks: msg('Ticks'),
  dinucleotide_content: msg('Dinucleotide content'),
  dinucleotide_skew: msg('Dinucleotide skew'),
  depth: msg('Depth'),
  sequence_conservation: msg('Pairwise comparison'),
  annotations: msg('Annotations'),
  spacer: msg('Spacer')
};

const STACK_ENTRY_AXIS = 'axis';
const STACK_ENTRY_SLOT = 'slot';
/**
 * The shown texts of a Hide GC Content or Hide GC Skew control. Each is a
 * whole sentence, so no locale fits a track or control name into a sentence
 * written around another word (OV-1477).
 * @typedef {object} SuppressTexts
 * @property {() => string} title the title of the dialog that asks before hiding
 * @property {(count: number) => string} description the dialog text: the custom stack rows the hide turns off
 * @property {() => string} hide the dialog's hide button
 * @property {(rows: string) => string} hideWith the dialog's hide button when it also turns off `rows`, the stack rows that reach the track
 * @property {() => string} keep the dialog's cancel button
 * @property {() => string} hiddenBy the note on a stack row that the control hides
 */
/** @type {Record<string, SuppressTexts>} */
const SUPPRESS_TEXT_BY_RENDERER = {
  dinucleotide_content: {
    title: () => t('Hide GC content?'),
    description: (count) => tn(count,
      'Your custom track stack has {count} GC content track. Hiding GC content turns that track off.',
      'Your custom track stack has {count} GC content tracks. Hiding GC content turns those tracks off.'),
    hide: () => t('Hide GC content'),
    hideWith: (rows) => t('Hide GC content and turn off {rows}', { rows }),
    keep: () => t('Keep GC content'),
    hiddenBy: () => t('Hidden: Hide GC Content is on.')
  },
  dinucleotide_skew: {
    title: () => t('Hide GC skew?'),
    description: (count) => tn(count,
      'Your custom track stack has {count} GC skew track. Hiding GC skew turns that track off.',
      'Your custom track stack has {count} GC skew tracks. Hiding GC skew turns those tracks off.'),
    hide: () => t('Hide GC skew'),
    hideWith: (rows) => t('Hide GC skew and turn off {rows}', { rows }),
    keep: () => t('Keep GC skew'),
    hiddenBy: () => t('Hidden: Hide GC Skew is on.')
  }
};
/** @type {Record<string, string>} */
const SUPPRESS_CONTROL_LABEL_BY_RENDERER = {
  dinucleotide_content: msg('Hide GC Content'),
  dinucleotide_skew: msg('Hide GC Skew')
};

/** @type {Record<string, string>} */
const PRESET_LABELS = {
  tuckin: msg('Tuckin'),
  middle: msg('Middle'),
  spreadout: msg('Spreadout')
};
// The Auto sizes a row asks for before Generate draws it (circular_auto_slot_sizes,
// gbdraw/diagrams/circular/presets.py): its width and gaps. Its radius is the
// layout's own result, so a row shows none until it is drawn (OV-1712).
const AUTO_SLOT_SIZES = MODE_PROFILE_DATA.autoSizes.circularSlots;
const ANNOTATION_MARK_OPTIONS = Object.freeze([
  'line',
  'bracket',
  'band',
  'highlight'
]);

/** @param {unknown} preset */
const formatPresetName = (preset) => tm(PRESET_LABELS[normalizeCircularTrackPreset(preset)] || PRESET_LABELS.tuckin);

/** @param {unknown} laneDirection */
const laneDirectionLabel = (laneDirection) => {
  const lane = normalizeLaneDirection(laneDirection);
  if (lane === 'outside') return t('Features outside axis');
  if (lane === 'split') return t('Features centered on axis');
  return t('Features inside axis');
};

/**
 * The Auto width of a row of `renderer`: by size class, and for a ticks row
 * also by preset and Separate Strands (never 0 px, TK-15); `null` for a
 * renderer without one.
 * @param {string} renderer @param {string} lengthParam @param {DrawingState} drawing
 * @returns {number | null}
 */
const autoWidthPx = (renderer, lengthParam, drawing) => {
  const size = lengthParam === 'short' ? 'short' : 'long';
  if (renderer === 'ticks') {
    const byStrands = ownEntry(AUTO_SLOT_SIZES.ticksWidthPx[size], normalizeCircularTrackPreset(drawing.form.track_type));
    return byStrands ? byStrands[drawing.form.separate_strands ? 'separate' : 'combined'] : null;
  }
  return ownEntry(AUTO_SLOT_SIZES.widthPx[size], renderer) ?? null;
};

/**
 * The size class of the estimates: that of the records they describe, and long
 * while no length is known.
 * @param {DrawingState} drawing the Circular drawing
 * @param {number[]} recordLengths the lengths of the records the estimates describe
 */
const getPreviewLengthParam = (drawing, recordLengths) => (
  sizeClassForRecordLengths('circular', recordLengths, drawing) || 'long'
);

/** @param {CircularTrackSlot | null | undefined} slot */
const slotHasManualGeometry = (slot) => (
  normalizeOptionalText(slot?.width) !== null ||
  normalizeOptionalText(slot?.radius) !== null ||
  normalizeOptionalText(slot?.inner_gap_px) !== null ||
  normalizeOptionalText(slot?.outer_gap_px) !== null
);

/** @param {TrackSlotStateView | null | undefined} state */
const circularAvailableDepthTrackCountForState = (state) => {
  const files = representativeDepthFiles(state?.files?.c_depth);
  return files.some(Boolean) ? files.length : 0;
};

/** @param {TrackSlotStateView | null | undefined} state */
const circularSourcedDepthTrackIndexesForState = (state) => (
  representativeDepthFiles(state?.files?.c_depth)
    .flatMap((file, trackIndex) => (file ? [trackIndex] : []))
);

/** @param {CircularTrackSlot | null | undefined} slot @param {unknown} [placement] */
const applyPlacementDefaults = (slot, placement = 'inside') => {
  if (!slot) return;
  const requestedSide = normalizePlacement(placement);
  slot.side = requestedSide;
  slot.params = cloneParams(slot.params);
  if (slot.renderer === 'features') {
    slot.params.lane_direction = laneDirectionForSide(requestedSide);
  }
};

/**
 * @param {unknown} slots
 * @param {unknown} [axisIndex]
 * @param {string} [preset]
 */
export const circularTrackAxisIndexForEnabledSlots = (
  slots,
  axisIndex = null,
  preset = 'tuckin'
) => {
  const source = Array.isArray(slots) ? slots : [];
  const resolvedAxis = (
    clampCircularTrackAxisIndex(axisIndex, source.length) ??
    inferLegacyAxisIndexFromFeature(source, preset)
  );
  return source
    .slice(0, resolvedAxis)
    .filter((slot) => slot?.enabled !== false)
    .length;
};

/**
 * @param {string} renderer
 * @param {({ id?: unknown, renderer?: unknown } | null)[]} [existingSlots]
 * @param {string} [nt]
 * @param {unknown} [placement] A side name; null leaves the side unset.
 * @returns {CircularTrackSlot}
 */
export const createCircularTrackSlotForRenderer = (renderer, existingSlots = [], nt = 'GC', placement = null) => {
  const normalizedRenderer = SUPPORTED_RENDERERS.includes(renderer) ? renderer : 'dinucleotide_skew';
  const id = nextNumberedId(
    DEFAULT_SLOT_IDS[normalizedRenderer] || normalizedRenderer,
    Array.isArray(existingSlots) ? existingSlots : []
  );

  /** @type {Record<string, unknown>} */
  const params = {};
  const side = normalizeSlotSide(placement);
  if (normalizedRenderer === 'ticks') {
    params.tick_label_layout = DEFAULT_TICK_LABEL_LAYOUT;
  } else if (normalizedRenderer === 'features' && side !== null) {
    params.lane_direction = laneDirectionForSide(side);
  } else if (normalizedRenderer === 'depth') {
    params.track_index = Math.max(
      0,
      (Array.isArray(existingSlots) ? existingSlots : []).filter((slot) => slot?.renderer === 'depth').length
    );
  } else if (normalizedRenderer === 'annotations') {
    params.set_id = '';
    params.overflow = 'error';
    params.show_labels = true;
    params.layer = 'foreground';
  }
  void nt;

  const slot = makeSlot({
    id,
    renderer: normalizedRenderer,
    side: normalizedRenderer === 'annotations' && side === null ? 'outside' : side,
    params
  });
  if (side !== null) applyPlacementDefaults(slot, side);
  return slot;
};

/** @param {string[]} options @param {string} key @param {unknown} value */
const appendOption = (options, key, value) => {
  const text = normalizeOptionalText(value);
  if (text === null) return;
  options.push(`${key}=${text}`);
};

/**
 * @param {TrackSlotRow} slot
 * @param {string} [defaultNt]
 * @param {string} [preset]
 * @param {{ includeSide?: boolean }} [optionsOverride]
 */
export const buildCircularTrackSlotSpec = (slot, defaultNt = 'GC', preset = 'tuckin', optionsOverride = {}) => {
  const normalized = normalizeCircularTrackSlot(slot, 0, defaultNt, preset);
  const options = [];
  const params = normalized.params || {};
  const normalizedPreset = normalizeCircularTrackPreset(preset);
  const includeSide = optionsOverride?.includeSide !== false;

  if (!normalized.enabled) options.push('enabled=false');
  appendOption(options, 'w', normalized.width);
  appendOption(options, 'r', normalized.radius);
  appendOption(options, 'inner_gap_px', parseOptionalPixel(
    normalized.inner_gap_px, `Circular track '${normalized.id}' inner_gap_px`, { allowZero: true }
  ));
  appendOption(options, 'outer_gap_px', parseOptionalPixel(
    normalized.outer_gap_px, `Circular track '${normalized.id}' outer_gap_px`, { allowZero: true }
  ));
  if (includeSide || normalizePlacement(normalized.side) === 'overlay') {
    appendOption(options, 'side', normalized.side);
  }
  if (Number.isFinite(Number(normalized.z)) && Number(normalized.z) !== 0) {
    options.push(`z=${Number(normalized.z)}`);
  }

  if (normalized.renderer === 'ticks') {
    appendOption(options, 'tick_label_layout', params.tick_label_layout);
    if (normalizeOptionalText(params.preset) !== null && normalizeCircularTrackPreset(params.preset) !== normalizedPreset) {
      appendOption(options, 'preset', params.preset);
    }
  } else if (normalized.renderer === 'features') {
    const placement = resolveCircularTrackFeaturePlacement(normalized, normalizedPreset);
    appendOption(options, 'lane_direction', placement.laneDirection);
  } else if (normalized.renderer === 'dinucleotide_content' || normalized.renderer === 'dinucleotide_skew') {
    const nt = normalizeOptionalText(params.nt);
    if (nt !== null && normalizeNt(nt) !== normalizeNt(defaultNt)) {
      options.push(`nt=${normalizeNt(nt)}`);
    }
    if (normalized.renderer === 'dinucleotide_skew') {
      appendOption(options, 'positive_color', params.positive_color);
      appendOption(options, 'negative_color', params.negative_color);
    }
  } else if (normalized.renderer === 'depth') {
    const trackIndex = normalizeTrackIndex(params.track_index);
    if (trackIndex !== null && trackIndex !== 0) {
      options.push(`track_index=${trackIndex}`);
    }
  } else if (normalized.renderer === 'sequence_conservation') {
    appendOption(options, 'track_index', params.track_index);
    appendOption(options, 'source_index', params.source_index);
  } else if (normalized.renderer === 'annotations') {
    appendOption(options, 'set_id', params.set_id);
    if (Array.isArray(params.marks) && params.marks.length > 0) {
      appendOption(options, 'marks', params.marks.join('|'));
    }
    appendOption(options, 'lane_gap_px', params.lane_gap_px);
    appendOption(options, 'padding_px', params.padding_px);
    appendOption(options, 'overflow', params.overflow);
    appendOption(options, 'show_labels', params.show_labels === false ? 'false' : 'true');
    appendOption(options, 'anchor_slot', params.anchor_slot);
    appendOption(options, 'to_slot', params.to_slot);
    appendOption(options, 'layer', params.layer);
    if (params.cover_anchor === true) {
      appendOption(options, 'cover_anchor', 'true');
    }
  }
  appendOption(options, 'legend_label', params.legend_label);

  return `${normalized.id}:${normalized.renderer}${options.length ? `@${options.join(',')}` : ''}`;
};

/** @param {TrackSlotStateView | null | undefined} state @param {DrawingState} drawing */
const conservationSourceFilesForState = (state, drawing) => {
  const blasts = normalizeFileList(state?.files?.c_conservation_blasts);
  if (
    String(drawing?.circularConservation?.source || '').trim().toLowerCase() === 'upload' ||
    (
      state?.files?.c_conservation_blasts_source === 'losat-cache' &&
      blasts.length > 0
    )
  ) {
    return blasts;
  }
  return normalizeFileList(state?.files?.c_conservation_fastas);
};

/** @param {TrackSlotStateView | null | undefined} state @param {DrawingState} drawing */
const conservationEntriesForState = (state, drawing) => {
  if (drawing?.circularConservation?.enabled !== true) return [];
  return orderedConservationSources(
    conservationSourceFilesForState(state, drawing),
    drawing.circularConservation
  );
};

/**
 * The Conservation series entry a managed row is drawn from; only the fields read here.
 * @typedef {{ sourceKey?: unknown, orderIndex?: number, label?: unknown, defaultLabel?: unknown, color?: unknown, fileName?: unknown }} ConservationSlotEntry
 */

/** @param {{ params?: { series_key?: unknown } | null } | null | undefined} slot */
const managedConservationSlotKey = (slot) => String(slot?.params?.series_key || '').trim();

/**
 * The managed ring rows no source draws: a row whose source has gone, and a
 * second row of one source.
 * @param {readonly CircularTrackSlot[]} slots
 * @param {readonly ConservationSlotEntry[]} entries
 * @returns {Set<CircularTrackSlot>}
 */
const managedRowsWithoutSource = (slots, entries) => {
  const desiredKeys = new Set(entries.map((entry) => String(entry.sourceKey || '')));
  /** @type {Set<string>} */
  const presentKeys = new Set();
  return new Set(slots.filter((slot) => {
    if (!isManagedConservationSlot(slot)) return false;
    const key = managedConservationSlotKey(slot);
    if (!key || !desiredKeys.has(key) || presentKeys.has(key)) return true;
    presentKeys.add(key);
    return false;
  }));
};

/**
 * @template {CircularTrackSlot | null | undefined} T
 * @param {T} slot
 * @param {ConservationSlotEntry | null | undefined} entry
 * @param {number} orderIndex
 * @returns {T}
 */
const refreshManagedConservationSlot = (slot, entry, orderIndex) => {
  if (!slot || !entry) return slot;
  slot.renderer = 'sequence_conservation';
  slot.side = normalizePlacement(slot.side, 'inside') === 'overlay' ? 'inside' : normalizePlacement(slot.side, 'inside');
  slot.params = cloneParams(slot.params);
  slot.params.managed = CONSERVATION_SLOT_MANAGER;
  slot.params.series_key = String(entry.sourceKey || '');
  slot.params.source_index = String(Number(entry.orderIndex ?? orderIndex));
  slot.params.track_index = String(Number(entry.orderIndex ?? orderIndex) + 1);
  slot.params.label = String(entry.label || entry.defaultLabel || entry.fileName || `Comparison ${Number(orderIndex) + 1}`);
  slot.params.color = String(entry.color || '');
  slot.params.fileName = String(entry.fileName || '');
  return slot;
};

/** @param {ConservationSlotEntry} entry @param {number} orderIndex @param {Set<string>} existingIds */
const makeManagedConservationSlot = (entry, orderIndex, existingIds) => refreshManagedConservationSlot(
  makeSlot({
    id: safeConservationSlotId(entry, orderIndex, existingIds),
    renderer: 'sequence_conservation',
    side: 'inside',
    params: {
      managed: CONSERVATION_SLOT_MANAGER
    }
  }),
  entry,
  orderIndex
);

/** @param {unknown} target @param {unknown} source */
const replaceObjectContents = (target, source) => {
  if (
    !target || typeof target !== 'object' || Array.isArray(target) ||
    !source || typeof source !== 'object' || Array.isArray(source)
  ) {
    return source;
  }
  Object.keys(target).forEach((key) => {
    // `target` is a plain object past the guard above.
    if (!Object.prototype.hasOwnProperty.call(source, key)) delete /** @type {Record<string, unknown>} */ (target)[key];
  });
  Object.assign(target, source);
  return target;
};

/** @param {DrawingState} drawing */
const circularGeometryShortcutsForState = (drawing) => ({
  featureWidth: drawing?.adv?.feature_width_circular,
  depthWidth: drawing?.adv?.depth_width_circular,
  gcContentWidth: drawing?.adv?.gc_content_width_circular,
  gcContentRadius: drawing?.adv?.gc_content_radius_circular,
  gcSkewWidth: drawing?.adv?.gc_skew_width_circular,
  gcSkewRadius: drawing?.adv?.gc_skew_radius_circular
});

/**
 * @typedef {object} CircularTrackSlotEditorOptions
 * @property {Record<string, any>} state the Web state; its shape belongs to `state.js`
 * @property {ChangeTrackLayout} [changeTrackLayout] R10 port (default: apply directly)
 * @property {AskChoice} [askChoice] The root's choice dialog: hiding a GC track the custom stack draws asks first.
 * @property {(label: string, fn: () => unknown) => unknown} [runUndoable] History's undoable step, which the chosen hide runs in.
 * @property {() => number[]} [estimateRecordLengths] The lengths of the records the estimates describe: those of the
 *   diagram of the next Generate that draws the shown Result's record (default: none known).
 * @property {ChangeTrackVisibility} [changeTrackVisibility] The Legend port a Hide GC Content or Hide GC Skew
 *   change and every stack edit run in: the drawn track's Legend rows leave and return with it (OV-1198, OV-1215;
 *   default: apply directly).
 */

// `changeTrackLayout` is the feature placement owner's transition, injected as
// a port (R10, Q3, R13): every stack edit that can change the feature slot
// runs through it.
/**
 * @param {CircularTrackSlotEditorOptions} options
 */
export const createCircularTrackSlotEditor = ({
  state, changeTrackLayout = (apply) => apply(), askChoice = async () => 'cancel', runUndoable = (_label, fn) => fn(),
  estimateRecordLengths = () => [], changeTrackVisibility = (change) => change()
}) => {
  /** @type {WeakMap<CircularTrackSlot, string>} */
  const editorKeys = new WeakMap();
  // A slot's Gaps disclosure stays open while its row is rebuilt (Generate and
  // Session load replace the slot objects, so the row is keyed anew).
  /** @type {Set<string>} */
  const openSlotOptions = new Set();
  /** @param {{ id?: unknown } | null | undefined} slot */
  const circularTrackSlotOptionsOpen = (slot) => openSlotOptions.has(String(slot?.id));
  /** @param {{ id?: unknown } | null | undefined} slot @param {Event} event */
  const setCircularTrackSlotOptionsOpen = (slot, event) => {
    const details = /** @type {HTMLDetailsElement} */ (event.currentTarget);
    if (details.open) openSlotOptions.add(String(slot?.id));
    else openSlotOptions.delete(String(slot?.id));
  };
  // The dinucleotide of the Circular drawing: what the default tracks draw, and
  // a content or skew slot that sets none (OV-487). The slot summary, the
  // preset summary, and the placeholders of a slot's fields name it from here.
  /** @returns {string} */
  const circularDrawingDinucleotide = () => normalizeNt(state.drawingOfMode('circular').adv.nt);
  // What the closed disclosure holds; a GC or AT track names its dinucleotide.
  /** @param {{ renderer?: unknown, params?: { nt?: unknown } } | null | undefined} slot */
  const circularTrackSlotOptionsSummary = (slot) => {
    if (slot?.renderer !== 'dinucleotide_content' && slot?.renderer !== 'dinucleotide_skew') return t('Gaps');
    const nt = normalizeNt(slot?.params?.nt, circularDrawingDinucleotide());
    return slot.renderer === 'dinucleotide_skew'
      ? t('{nt} · legend label, gaps, colors', { nt })
      : t('{nt} · legend label, gaps', { nt });
  };
  let nextEditorKey = 1;
  /** @param {CircularTrackSlot | null | undefined} slot */
  const circularTrackSlotEditorKey = (slot) => {
    if (!slot || typeof slot !== 'object') return 'circular-slot-invalid';
    let key = editorKeys.get(slot);
    if (!key) {
      key = `circular-editor-slot-${nextEditorKey}`;
      nextEditorKey += 1;
      editorKeys.set(slot, key);
    }
    return key;
  };

  /** @param {DrawingState} drawing @param {CircularTrackSlot[]} slots */
  const axisIndexForCurrentSlots = (drawing, slots) => {
    const current = clampCircularTrackAxisIndex(drawing.adv.circular_track_slots_axis_index, slots.length);
    if (current !== null) {
      drawing.adv.circular_track_slots_axis_index = current;
      return current;
    }
    const inferred = inferLegacyAxisIndexFromFeature(slots, drawing.form.track_type);
    drawing.adv.circular_track_slots_axis_index = inferred;
    return inferred;
  };

  /** @param {DrawingState} drawing */
  const normalizedSlotsForCurrentState = (drawing) => applyCircularTrackOrderPlacements(
    drawing.adv.circular_track_slots,
    drawing.adv.nt,
    drawing.form.track_type,
    drawing.adv.circular_track_slots_axis_index
  );

  /** @param {DrawingState} drawing */
  const annotationSetIds = (drawing) => (
    (Array.isArray(drawing.annotationSets) ? drawing.annotationSets : [])
      .map((set) => String(set?.id || '').trim())
      .filter(Boolean)
  );

  /** @param {DrawingState} drawing */
  const circularTrackValidationPlan = (drawing) => validateCustomTrackPlan({
    mode: 'circular',
    slots: drawing.adv.circular_track_slots,
    axisIndex: drawing.adv.circular_track_slots_axis_index,
    trackType: drawing.form.track_type,
    depthTrackCount: circularAvailableDepthTrackCountForState(state),
    depthSourcedTrackIndexes: circularSourcedDepthTrackIndexesForState(state),
    annotationSetIds: annotationSetIds(drawing),
    conservationSeries: conservationEntriesForState(state, drawing)
  });

  /** @param {CircularTrackSlot | null | undefined} slot @param {unknown} [index] */
  const circularTrackSlotIssue = (slot, index = null) => {
    const drawing = state.drawingOfMode('circular');
    const resolvedIndex = Number.isInteger(Number(index))
      ? Number(index)
      : drawing.adv.circular_track_slots.findIndex((/** @type {CircularTrackSlot} */ candidate) => candidate === slot);
    if (resolvedIndex < 0) return '';
    return errorText((circularTrackValidationPlan(drawing).rowIssues.get(resolvedIndex) || [])
      // The Width and Radius fields show their own invalid value (TK-12).
      .filter((issue) => !(issue.code === 'geometry_invalid' && ['width', 'radius'].includes(issue.field)))
      .flatMap((issue) => issue.parts));
  };

  /** The plan's issues, each shown as `errorText(issue.parts)` and keyed by its code (codes do not repeat). */
  const circularTrackGlobalIssues = () => circularTrackValidationPlan(state.drawingOfMode('circular')).globalIssues;

  /**
   * The slots an Overlay annotation row may name: for its anchor, every
   * enabled slot with a band; for Extend to, the same without its anchor.
   * @param {CircularTrackSlot | null} [slot]
   * @param {OverlaySlotReference} [key]
   * @returns {{ id: string, label: string }[]}
   */
  const circularAnnotationAnchorOptions = (slot = null, key = 'anchor_slot') => (
    state.drawingOfMode('circular').adv.circular_track_slots
      .filter((/** @type {CircularTrackSlot} */ candidate) => (
        candidate &&
        candidate !== slot &&
        candidate.enabled !== false &&
        !['annotations', 'ticks', 'spacer'].includes(String(candidate.renderer || '').trim()) &&
        String(candidate.id || '').trim() &&
        (key !== 'to_slot' || String(candidate.id).trim() !== String(slot?.params?.anchor_slot || '').trim())
      ))
      .map((/** @type {CircularTrackSlot} */ candidate) => ({
        id: String(candidate.id).trim(),
        label: t('{id} · {renderer}', { id: String(candidate.id).trim(), renderer: circularTrackRendererLabel(candidate.renderer) })
      }))
  );

  /** @param {CircularTrackSlot | null | undefined} slot @param {OverlaySlotReference} [key] */
  const circularAnnotationAnchorIsKnown = (slot, key = 'anchor_slot') => {
    const reference = String(slot?.params?.[key] || '').trim();
    return !reference || circularAnnotationAnchorOptions(slot, key).some((option) => option.id === reference);
  };

  /**
   * A row placed on Overlay keeps a valid anchor or takes the only one offered,
   * and the z its layer needs (OV-1297).
   * @param {DrawingState} drawing The drawing whose slot it is.
   * @param {CircularTrackSlot | null | undefined} slot
   */
  const bindOnlyCircularAnnotationAnchor = (drawing, slot) => {
    if (!slot || slot.renderer !== 'annotations') return;
    slot.params = cloneParams(slot.params);
    const options = circularAnnotationAnchorOptions(slot);
    const current = String(slot.params.anchor_slot || '').trim();
    if (!options.some((option) => option.id === current)) {
      if (options.length === 1) slot.params.anchor_slot = options[0].id;
      else delete slot.params.anchor_slot;
    }
    slot.z = overlayAnnotationZ(slot, drawing.adv.circular_track_slots);
  };

  /**
   * The Layer, anchor and Extend to selects of an Overlay annotation row (OV-1297).
   * @param {CircularTrackSlot | null | undefined} slot
   * @param {OverlaySlotReference | 'layer'} key
   * @param {unknown} value
   */
  const setCircularOverlayAnnotationParam = (slot, key, value) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (slot?.renderer !== 'annotations') return;
    setOverlayAnnotationParam(slot, state.drawingOfMode('circular').adv.circular_track_slots, key, value);
  };

  /** Why Add track cannot add this renderer, or '' (OV-1254). @param {unknown} renderer */
  const circularTrackAddBlockReason = (renderer) => {
    const normalizedRenderer = String(renderer || '').trim();
    return trackAddBlockReason(normalizedRenderer, {
      supported: UI_RENDERERS.includes(normalizedRenderer),
      annotationSetCount: annotationSetIds(state.drawingOfMode('circular')).length,
      depthSeriesCount: circularAvailableDepthTrackCountForState(state)
    });
  };

  /** @param {unknown} renderer */
  const canAddCircularTrackRenderer = (renderer) => circularTrackAddBlockReason(renderer) === '';

  // Duplicate shares the availability of the other row controls (GX-01).
  /** @param {CircularTrackSlot | null | undefined} slot */
  const canDuplicateCircularTrackSlot = (slot) => {
    const drawing = state.drawingOfMode('circular');
    if (!slot || isManagedConservationSlot(slot) || state.sessionOperationAvailability?.()) return false;
    if (slot.enabled === false) return true;
    if (slot.renderer === 'annotations') return annotationSetIds(drawing).length > 0;
    if (slot.renderer === 'depth') return circularAvailableDepthTrackCountForState(state) > 0;
    return true;
  };

  /** @param {CircularTrackSlot | null | undefined} slot @param {string} mark */
  const circularAnnotationMarkSelected = (slot, mark) => {
    const selected = Array.isArray(slot?.params?.marks)
      ? slot.params.marks.map((value) => String(value).trim().toLowerCase()).filter(Boolean)
      : [];
    return selected.length === 0 || selected.includes(mark);
  };

  /** @param {CircularTrackSlot | null | undefined} slot @param {string} mark @param {unknown} checked */
  const setCircularAnnotationMarkSelected = (slot, mark, checked) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!slot || slot.renderer !== 'annotations' || !ANNOTATION_MARK_OPTIONS.includes(mark)) return;
    slot.params = cloneParams(slot.params);
    const current = Array.isArray(slot.params.marks) && slot.params.marks.length > 0
      ? Array.from(new Set(slot.params.marks.map((value) => String(value).trim().toLowerCase())))
      : ANNOTATION_MARK_OPTIONS.slice();
    const next = checked
      ? Array.from(new Set([...current, mark]))
      : current.filter((value) => value !== mark);
    if (next.length === 0 || next.length === ANNOTATION_MARK_OPTIONS.length) delete slot.params.marks;
    else slot.params.marks = next;
  };

  /** @param {CircularTrackSlot | null | undefined} slot @param {string} field @param {number} defaultValue */
  const circularAnnotationNumberValue = (slot, field, defaultValue) => {
    const raw = slot?.params?.[field];
    return raw === null || raw === undefined || raw === '' ? defaultValue : raw;
  };

  /** @param {CircularTrackSlot | null | undefined} slot @param {string} field @param {unknown} value @param {number} defaultValue */
  const setCircularAnnotationNumber = (slot, field, value, defaultValue) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!slot || slot.renderer !== 'annotations') return;
    slot.params = cloneParams(slot.params);
    if (value === null || value === undefined || value === '') {
      delete slot.params[field];
      return;
    }
    const numeric = Number(value);
    if (Number.isFinite(numeric) && numeric === defaultValue) delete slot.params[field];
    else slot.params[field] = numeric;
  };

  /** @param {CircularTrackSlot | null | undefined} slot */
  const circularAnnotationCoverAnchor = (slot) => slot?.params?.cover_anchor === true;

  /** @param {CircularTrackSlot | null | undefined} slot @param {unknown} checked */
  const setCircularAnnotationCoverAnchor = (slot, checked) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!slot || slot.renderer !== 'annotations') return;
    slot.params = cloneParams(slot.params);
    if (checked) slot.params.cover_anchor = true;
    else delete slot.params.cover_anchor;
  };

  /**
   * A managed Depth row for a series, added to `takenSlots`. Its ID is the
   * preferred one numbered with the first unused number, except that the first
   * series takes `depth` while its preferred `depth_1` is taken and `depth` is not.
   * @param {unknown} trackIndex @param {{ id?: unknown }[]} takenSlots @param {unknown} desiredCount
   */
  const makeDepthSlotForTrackIndex = (trackIndex, takenSlots, desiredCount) => {
    const normalizedIndex = Math.max(0, Number(trackIndex) || 0);
    const preferredId = Number(desiredCount) > 1 ? `depth_${normalizedIndex + 1}` : 'depth';
    const numberedId = nextNumberedId(preferredId, takenSlots);
    const id = numberedId !== preferredId && normalizedIndex === 0 && nextNumberedId('depth', takenSlots) === 'depth'
      ? 'depth'
      : numberedId;
    const slot = makeSlot({
      id,
      renderer: 'depth',
      params: {
        track_index: normalizedIndex
      }
    });
    takenSlots.push(slot);
    return slot;
  };

  const normalizeSlotsInPlace = () => {
    const drawing = state.drawingOfMode('circular');
    const normalized = normalizeCircularTrackSlots(
      drawing.adv.circular_track_slots,
      drawing.adv.nt,
      drawing.form.track_type
    );
    const axis = axisIndexForCurrentSlots(drawing, normalized);
    syncSlotsFromAxisIndex(normalized, axis, drawing.form.track_type);
    drawing.adv.circular_track_slots_axis_index = enforceSingleOnAxisSlot(
      normalized,
      axis,
      drawing.form.track_type
    );
    const suppressed = applyCircularSuppressControlsToSlots(normalized, drawing.form, drawing.adv.nt);
    const identityPreserving = suppressed.map((slot, index) => (
      replaceObjectContents(drawing.adv.circular_track_slots[index], slot)
    ));
    drawing.adv.circular_track_slots.splice(
      0,
      drawing.adv.circular_track_slots.length,
      ...identityPreserving
    );
  };

  /**
   * @param {DrawingState} drawing
   * @param {{ nextSlots: CircularTrackSlot[], newSlots?: CircularTrackSlot[], managedPredicate: (slot: CircularTrackSlot | null | undefined) => boolean, preferredInsertIndex?: number | null }} mutation
   */
  const commitManagedSlotMutation = (drawing, {
    nextSlots,
    newSlots = [],
    managedPredicate,
    preferredInsertIndex = null
  }) => {
    /** @type {CircularTrackSlot[]} */
    const currentSlots = Array.isArray(drawing.adv.circular_track_slots)
      ? drawing.adv.circular_track_slots
      : [];
    const currentAxis = axisIndexForCurrentSlots(drawing, currentSlots);
    const retainedManagedIds = new Set(
      (Array.isArray(nextSlots) ? nextSlots : [])
        .filter((slot) => managedPredicate(slot))
        .map((slot) => String(slot?.id || '').trim())
        .filter(Boolean)
    );
    const removedBeforeAxis = currentSlots
      .slice(0, currentAxis)
      .filter((slot) => (
        managedPredicate(slot) &&
        !retainedManagedIds.has(String(slot?.id || '').trim())
      )).length;
    const committed = Array.isArray(nextSlots) ? nextSlots.slice() : [];
    const axis = Math.max(0, currentAxis - removedBeforeAxis);
    const additions = (Array.isArray(newSlots) ? newSlots : []).map((slot) => {
      applyPlacementDefaults(slot, 'inside');
      return slot;
    });
    if (additions.length > 0) {
      const preferred = Number.isInteger(preferredInsertIndex)
        // Number.isInteger is true only for a number, so it is not null here.
        ? /** @type {number} */ (preferredInsertIndex)
        : axis;
      const insideFloor = committed.reduce((floor, slot, index) => {
        if (
          index < axis ||
          slot?.renderer === 'annotations' ||
          effectiveSlotPlacement(slot, drawing.form.track_type) !== 'overlay'
        ) {
          return floor;
        }
        return Math.max(floor, index + 1);
      }, axis);
      const insertIndex = Math.max(
        insideFloor,
        Math.min(preferred, committed.length)
      );
      committed.splice(insertIndex, 0, ...additions);
    }
    drawing.adv.circular_track_slots_axis_index = Math.min(axis, committed.length);
    drawing.adv.circular_track_slots.splice(
      0,
      drawing.adv.circular_track_slots.length,
      ...committed
    );
    normalizeSlotsInPlace();
  };

  // Reset is the preset reset for the current Track layout.
  const resetCircularTrackSlotsFromSimpleControls = () => (
    resetCircularTrackSlotsToPreset(state.drawingOfMode('circular').form.track_type)
  );

  // Managed Depth rows follow Depth sources (PD-OI-058). The saved stack is
  // reconciled too, so a later Use custom stack shows the row; an empty saved
  // stack is built by Reset when the stack is enabled.
  /**
   * What a Depth file removal does to the stack, from the Depth file rows it
   * would leave. It writes nothing.
   * @param {DrawingState} drawing @param {DepthFileSlots[]} nextDepth
   */
  const circularDepthSourceLoss = (drawing, nextDepth) => depthSourceLoss(
    Array.isArray(drawing.adv.circular_track_slots) ? drawing.adv.circular_track_slots : [],
    circularSourcedDepthTrackIndexesForState(state),
    circularSourcedDepthTrackIndexesForState({ files: { c_depth: nextDepth } }),
    (rows) => draftDrawsDepth('circular', { form: drawing.form, adv: { ...drawing.adv, circular_track_slots: rows } })
  );
  /** @param {DrawingState} drawing @param {number[]} previousSourced */
  const reconcileCircularDepthSlots = (drawing, previousSourced) => {
    const slots = Array.isArray(drawing.adv.circular_track_slots) ? drawing.adv.circular_track_slots : [];
    if (slots.length === 0 && !drawing.adv.circular_track_slots_enabled) return;
    // The reconciled rows are the drawing's own rows (R13: the service types them as DepthSlotRow).
    const { slots: nextSlots, additions } = /** @type {{ slots: CircularTrackSlot[], additions: number[] }} */ (
      reconcileManagedDepthSlots({ slots, previousSourced, sourced: circularSourcedDepthTrackIndexesForState(state) })
    );
    if (additions.length === 0 && nextSlots.length === slots.length) return;
    const seriesCount = circularAvailableDepthTrackCountForState(state);
    /** @type {{ id?: unknown }[]} */
    const takenSlots = [...nextSlots];
    // The shortcut copies of the new rows are rows (the service types them as Record<string, unknown>).
    const newSlots = /** @type {CircularTrackSlot[]} */ (applyCircularGeometryShortcuts(
      additions.map((trackIndex) => makeDepthSlotForTrackIndex(trackIndex, takenSlots, seriesCount)),
      circularGeometryShortcutsForState(drawing)
    ));
    let preferredInsertIndex = -1;
    nextSlots.forEach((slot, index) => {
      if (slot?.renderer === 'depth') preferredInsertIndex = index;
    });
    if (preferredInsertIndex < 0) {
      preferredInsertIndex = nextSlots.findIndex((slot) => slot?.renderer === 'ticks');
    }
    if (preferredInsertIndex < 0) {
      preferredInsertIndex = nextSlots.findIndex((slot) => slot?.renderer === 'features');
    }
    commitManagedSlotMutation(drawing, {
      nextSlots,
      newSlots,
      managedPredicate: isDefaultManagedDepthSlot,
      preferredInsertIndex: Math.max(0, preferredInsertIndex + 1)
    });
  };

  // The only entry for Depth source changes: run the change, then reconcile.
  // The Circular drawing's Show Depth goes off with its last Depth
  // source; the other mode's drawing is not touched (OV-82, OV-101).
  /** @param {() => unknown} mutate */
  const changeCircularDepthSources = (mutate) => {
    const drawing = state.drawingOfMode('circular');
    const previousSourced = circularSourcedDepthTrackIndexesForState(state);
    mutate();
    reconcileCircularDepthSlots(drawing, previousSourced);
    if (previousSourced.length > 0 && circularSourcedDepthTrackIndexesForState(state).length === 0) drawing.form.show_depth = false;
  };

  /**
   * The slot ids of the managed ring rows a change of the comparison sources
   * takes out of the stack, from the sources it would leave. It writes nothing.
   * The rows follow their sources only while the stack is in use
   * (`refreshCircularConservationSlots`).
   * @param {DrawingState} drawing @param {readonly ConservationSlotEntry[]} entries
   * @returns {Set<string>}
   */
  const circularConservationRowsLeaving = (drawing, entries) => (
    drawing.adv.circular_track_slots_enabled === true
      ? new Set([...managedRowsWithoutSource(drawing.adv.circular_track_slots, entries)].map((slot) => String(slot.id)))
      : new Set()
  );

  const syncCircularConservationSlots = () => {
    const drawing = state.drawingOfMode('circular');
    const entries = conservationEntriesForState(state, drawing);
    /** @type {CircularTrackSlot[]} */
    const slots = Array.isArray(drawing.adv.circular_track_slots)
      ? drawing.adv.circular_track_slots
      : [];
    const leaving = managedRowsWithoutSource(slots, entries);
    const existingIds = new Set(slots.map((slot) => String(slot?.id || '').trim()).filter(Boolean));
    const entryByKey = new Map(entries.map((entry) => [String(entry.sourceKey || ''), entry]));
    const entryOrderByKey = new Map(entries.map((entry, index) => [String(entry.sourceKey || ''), index]));

    /** @type {CircularTrackSlot[]} */
    const nextSlots = [];
    const presentKeys = new Set();
    slots.forEach((slot) => {
      if (leaving.has(slot)) return;
      if (!isManagedConservationSlot(slot)) {
        nextSlots.push(slot);
        return;
      }
      const key = managedConservationSlotKey(slot);
      refreshManagedConservationSlot(slot, entryByKey.get(key), entryOrderByKey.get(key) ?? presentKeys.size);
      presentKeys.add(key);
      nextSlots.push(slot);
    });

    /** @type {CircularTrackSlot[]} */
    const missingSlots = [];
    entries.forEach((entry, index) => {
      const key = String(entry.sourceKey || '');
      if (!key || presentKeys.has(key)) return;
      missingSlots.push(makeManagedConservationSlot(entry, index, existingIds));
      presentKeys.add(key);
    });

    let preferredInsertIndex = -1;
    nextSlots.forEach((slot, index) => {
      if (isManagedConservationSlot(slot)) preferredInsertIndex = index;
    });
    if (preferredInsertIndex < 0) {
      preferredInsertIndex = nextSlots.findIndex((slot) => slot?.renderer === 'features');
    }
    commitManagedSlotMutation(drawing, {
      nextSlots,
      newSlots: missingSlots,
      managedPredicate: isManagedConservationSlot,
      preferredInsertIndex: Math.max(0, preferredInsertIndex + 1)
    });
  };

  /** @param {unknown} preset */
  const resetCircularTrackSlotsToPreset = (preset) => {
    const drawing = state.drawingOfMode('circular');
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const normalizedPreset = normalizeCircularTrackPreset(preset);
    // One Depth row per loaded Depth series; Show Depth is not an input (TK-10).
    const depthTrackCount = circularAvailableDepthTrackCountForState(state);
    const templateSlots = applyCircularGeometryShortcuts(createDefaultCircularTrackSlots({
      nt: drawing.adv.nt,
      showDepth: depthTrackCount > 0,
      depthTrackCount,
      showGc: !drawing.form.suppress_gc,
      showSkew: !drawing.form.suppress_skew,
      showTicks: drawing.form.show_scale !== false,
      preset: normalizedPreset
    }), circularGeometryShortcutsForState(drawing));
    drawing.adv.circular_track_slots_axis_index = inferLegacyAxisIndexFromFeature(
      normalizeCircularTrackSlots(templateSlots, drawing.adv.nt, normalizedPreset),
      normalizedPreset
    );
    const normalized = applyCircularTrackOrderPlacements(
      templateSlots,
      drawing.adv.nt,
      normalizedPreset,
      drawing.adv.circular_track_slots_axis_index
    );
    drawing.form.track_type = normalizedPreset;
    drawing.adv.circular_track_slots.splice(0, drawing.adv.circular_track_slots.length, ...normalized);
    syncCircularConservationSlots();
  };

  /** @param {unknown} enabled */
  const setCircularTrackSlotsEnabled = (enabled) => {
    const drawing = state.drawingOfMode('circular');
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    drawing.adv.circular_track_slots_enabled = Boolean(enabled);
    if (!drawing.adv.circular_track_slots_enabled) return;
    if (!Array.isArray(drawing.adv.circular_track_slots) || drawing.adv.circular_track_slots.length === 0) {
      resetCircularTrackSlotsFromSimpleControls();
    } else {
      // The stack it turns on holds the managed comparison ring rows.
      syncCircularConservationSlots();
    }
  };

  /** @param {string} renderer @param {unknown} [placement] */
  const addCircularTrackSlot = (renderer, placement = null) => {
    const drawing = state.drawingOfMode('circular');
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!canAddCircularTrackRenderer(renderer)) return;
    normalizeSlotsInPlace();
    const slot = createCircularTrackSlotForRenderer(renderer, drawing.adv.circular_track_slots, drawing.adv.nt, placement);
    if (slot.renderer === 'annotations') {
      slot.params.set_id = String(drawing.annotationSets?.[0]?.id || '').trim();
    } else if (slot.renderer === 'depth') {
      const available = circularAvailableDepthTrackCountForState(state);
      const claimed = new Set(
        drawing.adv.circular_track_slots
          .filter((/** @type {CircularTrackSlot} */ candidate) => candidate?.enabled !== false && candidate?.renderer === 'depth')
          .map((/** @type {CircularTrackSlot} */ candidate) => normalizeTrackIndex(candidate?.params?.track_index))
          .filter((/** @type {number | null} */ trackIndex) => trackIndex !== null)
      );
      let trackIndex = 0;
      while (trackIndex < available && claimed.has(trackIndex)) trackIndex += 1;
      slot.params.track_index = trackIndex < available ? trackIndex : 0;
    }
    drawing.adv.circular_track_slots.push(slot);
    normalizeSlotsInPlace();
  };

  /** @param {unknown} index */
  const duplicateCircularTrackSlot = (index) => {
    const drawing = state.drawingOfMode('circular');
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    normalizeSlotsInPlace();
    const idx = Number(index);
    if (!Number.isInteger(idx) || idx < 0 || idx >= drawing.adv.circular_track_slots.length) return;
    if (!canDuplicateCircularTrackSlot(drawing.adv.circular_track_slots[idx])) return;
    const source = normalizeCircularTrackSlot(drawing.adv.circular_track_slots[idx], idx, drawing.adv.nt, drawing.form.track_type);
    const duplicate = createCircularTrackSlotForRenderer(source.renderer, drawing.adv.circular_track_slots, drawing.adv.nt);
    duplicate.id = nextCopyId(String(source.id || '') || duplicate.id, drawing.adv.circular_track_slots);
    duplicate.enabled = source.enabled;
    duplicate.width = source.width;
    duplicate.radius = source.radius;
    duplicate.inner_gap_px = source.inner_gap_px;
    duplicate.outer_gap_px = source.outer_gap_px;
    duplicate.side = source.side;
    duplicate.z = source.z;
    duplicate.params = cloneParams(source.params);
    drawing.adv.circular_track_slots.splice(idx + 1, 0, duplicate);
    const axis = axisIndexForCurrentSlots(drawing, drawing.adv.circular_track_slots);
    if (idx < axis) drawing.adv.circular_track_slots_axis_index = axis + 1;
    normalizeSlotsInPlace();
  };

  /**
   * @param {unknown} index
   * @param {(() => void) | Event} [turnOffReaching] Turns off the rows that start at or extend to it
   *   (`askTrackHide`); a template's trailing event never arrives (`featureSlotEdits` takes it).
   */
  const removeCircularTrackSlot = (index, turnOffReaching) => {
    const drawing = state.drawingOfMode('circular');
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const idx = Number(index);
    if (!Number.isInteger(idx) || idx < 0 || idx >= drawing.adv.circular_track_slots.length) return;
    const axis = axisIndexForCurrentSlots(drawing, drawing.adv.circular_track_slots);
    drawing.adv.circular_track_slots.splice(idx, 1);
    drawing.adv.circular_track_slots_axis_index = idx < axis ? Math.max(0, axis - 1) : Math.min(axis, drawing.adv.circular_track_slots.length);
    if (typeof turnOffReaching === 'function') turnOffReaching();
    normalizeSlotsInPlace();
  };

  /** @param {DrawingState} drawing @param {unknown} fromIndex @param {unknown} toIndex */
  const wouldCircularTrackSlotMoveCrossAxis = (drawing, fromIndex, toIndex) => {
    const from = Number(fromIndex);
    const to = Number(toIndex);
    const normalized = normalizedSlotsForCurrentState(drawing);
    if (
      !Number.isInteger(from) ||
      !Number.isInteger(to) ||
      from < 0 ||
      to < 0 ||
      from >= normalized.length ||
      to >= normalized.length ||
      from === to
    ) {
      return true;
    }

    const axis = axisIndexForCurrentSlots(drawing, normalized);
    const movedPlacement = effectiveSlotPlacement(normalized[from], drawing.form.track_type);
    if (movedPlacement === 'overlay') return true;
    return (from < axis) !== (to < axis);
  };

  /** @param {unknown} fromIndex @param {unknown} toIndex */
  const moveCircularTrackSlot = (fromIndex, toIndex) => {
    const drawing = state.drawingOfMode('circular');
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (wouldCircularTrackSlotMoveCrossAxis(drawing, fromIndex, toIndex)) return;
    normalizeSlotsInPlace();
    const from = Number(fromIndex);
    const to = Number(toIndex);
    if (
      !Number.isInteger(from) ||
      !Number.isInteger(to) ||
      from < 0 ||
      to < 0 ||
      from >= drawing.adv.circular_track_slots.length ||
      to >= drawing.adv.circular_track_slots.length ||
      from === to
    ) {
      return;
    }
    const [moved] = drawing.adv.circular_track_slots.splice(from, 1);
    drawing.adv.circular_track_slots.splice(to, 0, moved);
    normalizeSlotsInPlace();
  };

  /** @param {unknown} index @param {unknown} direction */
  const canMoveCircularTrackSlot = (index, direction) => {
    const drawing = state.drawingOfMode('circular');
    const idx = Number(index);
    const step = Number(direction);
    if (!Number.isInteger(idx) || !Number.isInteger(step) || step === 0) return false;
    const target = idx + Math.sign(step);
    return !wouldCircularTrackSlotMoveCrossAxis(drawing, idx, target);
  };

  /** @param {unknown} index */
  const canMoveCircularTrackSlotOutside = (index) => {
    const drawing = state.drawingOfMode('circular');
    const idx = Number(index);
    const normalized = normalizedSlotsForCurrentState(drawing);
    if (!Number.isInteger(idx) || idx < 0 || idx >= normalized.length) return false;
    if (effectiveSlotPlacement(normalized[idx], drawing.form.track_type) === 'overlay') return true;
    return idx >= axisIndexForCurrentSlots(drawing, normalized);
  };

  /** @param {unknown} index */
  const canMoveCircularTrackSlotInside = (index) => {
    const drawing = state.drawingOfMode('circular');
    const idx = Number(index);
    const normalized = normalizedSlotsForCurrentState(drawing);
    if (!Number.isInteger(idx) || idx < 0 || idx >= normalized.length) return false;
    if (effectiveSlotPlacement(normalized[idx], drawing.form.track_type) === 'overlay') return true;
    return idx < axisIndexForCurrentSlots(drawing, normalized);
  };

  /** @param {unknown} index */
  const canMoveCircularTrackSlotToAxis = (index) => {
    const drawing = state.drawingOfMode('circular');
    const idx = Number(index);
    const normalized = normalizedSlotsForCurrentState(drawing);
    if (!Number.isInteger(idx) || idx < 0 || idx >= normalized.length) return false;
    const slot = normalized[idx];
    return ['features', 'ticks', 'annotations'].includes(slot?.renderer) && effectiveSlotPlacement(slot, drawing.form.track_type) !== 'overlay';
  };

  /** @param {DrawingState} drawing @param {unknown} index @param {unknown} placement */
  const moveCircularTrackSlotToPlacement = (drawing, index, placement) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const idx = Number(index);
    if (!Number.isInteger(idx) || idx < 0 || idx >= drawing.adv.circular_track_slots.length) return;
    normalizeSlotsInPlace();
    if (idx >= drawing.adv.circular_track_slots.length) return;
    const targetPlacement = normalizePlacement(placement);

    const directOverlaySlot = drawing.adv.circular_track_slots[idx];
    if (directOverlaySlot?.renderer === 'annotations' && targetPlacement !== 'overlay') {
      directOverlaySlot.params = cloneParams(directOverlaySlot.params);
      OVERLAY_ANNOTATION_PARAMS.forEach((key) => { delete directOverlaySlot.params[key]; });
    }
    if (targetPlacement === 'overlay' && directOverlaySlot?.renderer === 'annotations') {
      syncSlotPlacementFromSide(directOverlaySlot, 'overlay');
      directOverlaySlot.params = cloneParams(directOverlaySlot.params);
      bindOnlyCircularAnnotationAnchor(drawing, directOverlaySlot);
      normalizeSlotsInPlace();
      return;
    }

    if (targetPlacement === 'overlay') {
      const movedSlot = drawing.adv.circular_track_slots[idx];
      if (!movedSlot) return;
      const movedPreviousPlacement = effectiveSlotPlacement(movedSlot, drawing.form.track_type);
      const existingAxisIndex = drawing.adv.circular_track_slots.findIndex((/** @type {CircularTrackSlot} */ slot, /** @type {number} */ slotIndex) => (
        slotIndex !== idx &&
        effectiveSlotPlacement(slot, drawing.form.track_type) === 'overlay'
      ));
      if (existingAxisIndex >= 0) {
        const existingAxisSlot = drawing.adv.circular_track_slots[existingAxisIndex];
        const demotedPlacement = movedPreviousPlacement === 'overlay' ? 'inside' : movedPreviousPlacement;
        syncSlotPlacementFromSide(existingAxisSlot, demotedPlacement);
        syncSlotPlacementFromSide(movedSlot, 'overlay');
        drawing.adv.circular_track_slots[existingAxisIndex] = movedSlot;
        drawing.adv.circular_track_slots[idx] = existingAxisSlot;
        drawing.adv.circular_track_slots_axis_index = existingAxisIndex;
        normalizeSlotsInPlace();
        return;
      }
    }

    let axis = axisIndexForCurrentSlots(drawing, drawing.adv.circular_track_slots);
    const [slot] = drawing.adv.circular_track_slots.splice(idx, 1);
    if (!slot) return;
    if (idx < axis) axis -= 1;
    syncSlotPlacementFromSide(slot, targetPlacement);
    if (targetPlacement === 'outside') {
      drawing.adv.circular_track_slots.splice(axis, 0, slot);
      drawing.adv.circular_track_slots_axis_index = axis + 1;
    } else if (targetPlacement === 'overlay') {
      drawing.adv.circular_track_slots.splice(axis, 0, slot);
      drawing.adv.circular_track_slots_axis_index = axis;
    } else {
      const onAxisIndex = drawing.adv.circular_track_slots.findIndex((/** @type {CircularTrackSlot} */ candidate) => (
        effectiveSlotPlacement(candidate, drawing.form.track_type) === 'overlay'
      ));
      const insertIndex = onAxisIndex >= 0 ? onAxisIndex + 1 : axis;
      drawing.adv.circular_track_slots.splice(insertIndex, 0, slot);
      drawing.adv.circular_track_slots_axis_index = onAxisIndex >= 0 ? onAxisIndex : axis;
    }
    normalizeSlotsInPlace();
  };

  /** @param {unknown} index */
  const moveCircularTrackSlotOutside = (index) => {
    const drawing = state.drawingOfMode('circular');
    if (!canMoveCircularTrackSlotOutside(index)) return;
    moveCircularTrackSlotToPlacement(drawing, index, 'outside');
  };

  /** @param {unknown} index */
  const moveCircularTrackSlotInside = (index) => {
    const drawing = state.drawingOfMode('circular');
    if (!canMoveCircularTrackSlotInside(index)) return;
    moveCircularTrackSlotToPlacement(drawing, index, 'inside');
  };

  /** @param {unknown} index */
  const moveCircularTrackSlotToAxis = (index) => {
    const drawing = state.drawingOfMode('circular');
    if (!canMoveCircularTrackSlotToAxis(index)) return;
    moveCircularTrackSlotToPlacement(drawing, index, 'overlay');
  };

  /** @param {CircularTrackSlot | null | undefined} slot @param {string | undefined} renderer */
  const updateCircularTrackSlotRenderer = (slot, renderer) => {
    const drawing = state.drawingOfMode('circular');
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    // A renderer that is still undefined is no member of the list, so the guard below returns.
    renderer = /** @type {string} */ (renderer || slot?.renderer);
    if (!slot || !SUPPORTED_RENDERERS.includes(renderer)) return;
    slot.params = renderer === slot.renderer ? cloneParams(slot.params) : paramsKeptOnRendererChange(slot.params, renderer);
    slot.renderer = renderer;
    slot.side = normalizeSlotSide(slot.side);
    if (renderer === 'ticks') {
      delete slot.params['axis'];
      slot.params.tick_label_layout = normalizeTickLabelLayout(
        slot.params.tick_label_layout ??
          tickLabelLayoutFromSides(slot.params.label_side, slot.params.tick_side)
      );
      delete slot.params.label_side;
      delete slot.params.tick_side;
      if (normalizeOptionalText(slot.params.preset) === null) delete slot.params.preset;
      else slot.params.preset = normalizeCircularTrackPreset(slot.params.preset);
    } else if (renderer === 'features') {
      slot.params.lane_direction = laneDirectionForSide(slot.side);
    } else if (renderer === 'dinucleotide_content' || renderer === 'dinucleotide_skew') {
      slot.params.nt = normalizeNt(slot.params.nt, circularDrawingDinucleotide());
    } else if (renderer === 'annotations') {
      slot.side = slot.side === 'overlay' ? 'overlay' : 'outside';
      slot.params.set_id = String(slot.params.set_id || drawing.annotationSets?.[0]?.id || '');
      slot.params.overflow = 'error';
      slot.params.show_labels = true;
      slot.params.layer = 'foreground';
    }
    normalizeSlotsInPlace();
  };

  /** @param {DrawingState} drawing @param {string} renderer */
  const activeCircularTrackSlotsForRenderer = (drawing, renderer) => {
    normalizeSlotsInPlace();
    return drawing.adv.circular_track_slots.filter((/** @type {CircularTrackSlot} */ slot) => (
      slot &&
      slot.enabled !== false &&
      String(slot.renderer || '').trim() === renderer &&
      slotUsesDrawingDinucleotide(slot, drawing.adv.nt)
    ));
  };

  /**
   * Hiding a GC track turns off the enabled stack rows that draw it, also with
   * the custom stack off. It asks first, naming the rows that start at or extend
   * to them (`askTrackHide`), with the stack on also when no row does; null when
   * nothing asks.
   * @param {DrawingState} drawing
   * @param {string} renderer
   * @param {(turnOffReaching: () => void) => unknown} apply
   * @param {() => boolean} unchanged
   * @returns {Promise<unknown> | null}
   */
  const askCircularSuppressOverride = (drawing, renderer, apply, unchanged) => {
    const hidden = activeCircularTrackSlotsForRenderer(drawing, renderer);
    // The Hide control of GC content or GC skew.
    const texts = ownEntry(SUPPRESS_TEXT_BY_RENDERER, renderer);
    if (hidden.length === 0 || !texts) return null;
    const stackInUse = Boolean(drawing.adv.circular_track_slots_enabled);
    return askTrackHide({ askChoice, runUndoable }, {
      rows: () => drawing.adv.circular_track_slots,
      mode: 'circular',
      hiding: new Set(hidden.map((/** @type {CircularTrackSlot} */ slot) => String(slot.id))),
      asksAlone: stackInUse,
      question: (numbers) => ({
        title: texts.title(),
        // "Your custom track stack has..." is not the drawing's state when the stack is off.
        ...(stackInUse ? { description: texts.description(hidden.length) } : {}),
        action: numbers.length ? texts.hideWith(trackRowsText(numbers)) : texts.hide(),
        keep: texts.keep()
      }),
      label: ownEntry(SUPPRESS_CONTROL_LABEL_BY_RENDERER, renderer) || msg('Hide track'),
      unchanged,
      apply
    });
  };

  /** @typedef {{ target?: { checked: boolean } | null } | null} CircularSuppressToggleEvent The checkbox change event of a Suppress control. */

  /**
   * @param {DrawingState} drawing
   * @param {string} key
   * @param {unknown} checked
   * @param {CircularSuppressToggleEvent} [event]
   */
  const setCircularSuppressControl = (drawing, key, checked, event = null) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const renderer = SUPPRESS_RENDERER_BY_KEY[key];
    if (!renderer) return;
    // Every renderer of SUPPRESS_RENDERER_BY_KEY has a form key.
    const formKey = /** @type {string} */ (ownEntry(SUPPRESS_FORM_KEY_BY_RENDERER, renderer));
    const nextChecked = Boolean(checked);
    const previousChecked = Boolean(drawing.form?.[formKey]);

    if (nextChecked === previousChecked) {
      if (event?.target) event.target.checked = previousChecked;
      return;
    }

    /** @param {() => void} [turnOffReaching] */
    const apply = (turnOffReaching = () => {}) => changeTrackVisibility(() => {
      drawing.form[formKey] = nextChecked;
      turnOffReaching();
      normalizeSlotsInPlace();
    });
    const asked = nextChecked ? askCircularSuppressOverride(drawing, renderer, apply,
      () => state.activeDrawing() === drawing && !state.sessionOperationAvailability?.() && Boolean(drawing.form[formKey]) === previousChecked) : null;
    if (asked) {
      // The checkbox shows the kept value until the choice; the chosen hide is one History step.
      if (event?.target) event.target.checked = previousChecked;
      return asked;
    }
    apply();
    if (event?.target) event.target.checked = nextChecked;
  };

  // The Dinucleotide setting and a row's dinucleotide are inputs of the Hide
  // GC projection (TK-09), so their edits run it here (R10).
  /** @param {DrawingState} drawing */
  const projectCircularSuppressControls = (drawing) => {
    const slots = drawing.adv.circular_track_slots;
    applyCircularSuppressControlsToSlots(slots, drawing.form, drawing.adv.nt)
      .forEach((slot, index) => replaceObjectContents(slots[index], slot));
  };

  /** @param {unknown} value */
  const setCircularDinucleotide = (value) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const drawing = state.drawingOfMode('circular');
    drawing.adv.nt = String(value ?? '');
    projectCircularSuppressControls(drawing);
  };

  /**
   * @param {CircularTrackSlot | null | undefined} slot
   * @param {unknown} value
   */
  const setCircularTrackSlotNt = (slot, value) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!slot) return;
    slot.params = cloneParams(slot.params);
    slot.params.nt = String(value ?? '').trim();
    projectCircularSuppressControls(state.drawingOfMode('circular'));
  };

  /** @param {unknown} checked @param {CircularSuppressToggleEvent} [event] */
  const setCircularGcSuppressed = (checked, event = null) => {
    const drawing = state.drawingOfMode('circular');
    return setCircularSuppressControl(drawing, 'gc_content', checked, event);
  };

  /** @param {unknown} checked @param {CircularSuppressToggleEvent} [event] */
  const setCircularSkewSuppressed = (checked, event = null) => {
    const drawing = state.drawingOfMode('circular');
    return setCircularSuppressControl(drawing, 'gc_skew', checked, event);
  };

  /**
   * @param {CircularTrackSlot | null | undefined} slot @param {unknown} enabled
   * @param {(() => void) | Event} [turnOffReaching] Turns off the rows that start at or extend to it
   *   (`askTrackHide`); a template's trailing event never arrives (`featureSlotEdits` takes it).
   */
  const setCircularTrackSlotEnabled = (slot, enabled, turnOffReaching) => {
    const drawing = state.drawingOfMode('circular');
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!slot || circularTrackSlotHiddenBySuppressForm(slot, drawing.form, drawing.adv.nt)) return;
    slot.enabled = Boolean(enabled);
    if (typeof turnOffReaching === 'function') turnOffReaching();
    if (slot.params && typeof slot.params === 'object' && !Array.isArray(slot.params)) {
      const token = ownEntry(SUPPRESS_KEY_BY_RENDERER, String(slot.renderer || '').trim());
      if (token && slot.params[GLOBAL_SUPPRESS_PARAM] === token) {
        slot.params = cloneParams(slot.params);
        delete slot.params[GLOBAL_SUPPRESS_PARAM];
      }
    }
    normalizeSlotsInPlace();
  };

  /** @param {CircularTrackSlot | null | undefined} slot */
  const circularTrackSlotHiddenBySuppress = (slot) => {
    const drawing = state.drawingOfMode('circular');
    return circularTrackSlotHiddenBySuppressForm(slot, drawing.form, drawing.adv.nt);
  };

  /** @param {CircularTrackSlot | null | undefined} slot */
  const circularTrackSlotEffectiveEnabled = (slot) => (
    Boolean(slot?.enabled !== false) && !circularTrackSlotHiddenBySuppress(slot)
  );

  /** @param {CircularTrackSlot | null | undefined} slot */
  const circularTrackSlotSuppressMessage = (slot) => {
    if (!circularTrackSlotHiddenBySuppress(slot)) return '';
    return ownEntry(SUPPRESS_TEXT_BY_RENDERER, String(slot?.renderer || '').trim())?.hiddenBy() ?? '';
  };

  // A placement select (the Annotations row's) moves the row across the Axis
  // boundary as the row menu does; writing its side alone made the boundary
  // follow it and flipped the other rows (OV-1263).
  /** @param {CircularTrackSlot | null | undefined} slot @param {unknown} placement */
  const updateCircularTrackSlotPlacement = (slot, placement) => {
    const drawing = state.drawingOfMode('circular');
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!slot || !SUPPORTED_RENDERERS.includes(slot.renderer) || normalizeOptionalText(placement) === null) return;
    const index = drawing.adv.circular_track_slots.indexOf(slot);
    if (index < 0 || effectiveSlotPlacement(slot, drawing.form.track_type) === normalizePlacement(placement)) return;
    return moveCircularTrackSlotToPlacement(drawing, index, placement);
  };

  /** @param {CircularTrackSlot | null | undefined} slot @param {unknown} laneDirection */
  const updateCircularTrackFeatureLane = (slot, laneDirection) => {
    const drawing = state.drawingOfMode('circular');
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!slot || slot.renderer !== 'features') return;
    slot.params = cloneParams(slot.params);
    const explicitLaneDirection = normalizeOptionalText(laneDirection);
    if (explicitLaneDirection === null) {
      delete slot.params.lane_direction;
      slot.side = null;
      normalizeSlotsInPlace();
      return;
    }
    const normalizedLaneDirection = normalizeLaneDirection(explicitLaneDirection);
    const index = drawing.adv.circular_track_slots.findIndex((/** @type {CircularTrackSlot} */ candidate) => candidate === slot);
    if (index >= 0) {
      moveCircularTrackSlotToPlacement(drawing, index, sideForLaneDirection(normalizedLaneDirection));
      return;
    }
    slot.params.lane_direction = normalizedLaneDirection;
    slot.side = sideForLaneDirection(slot.params.lane_direction);
    normalizeSlotsInPlace();
  };

  /** @param {string} renderer */
  const circularTrackRendererLabel = (renderer) => tm(ownEntry(RENDERER_LABELS, renderer) || renderer);
  /** @param {string} renderer */
  const supportsCircularTrackSlotPlacement = (renderer) => SUPPORTED_RENDERERS.includes(renderer);
  /** @param {CircularTrackSlot | null | undefined} slot */
  const isManagedCircularConservationSlot = (slot) => isManagedConservationSlot(slot);
  /** @param {CircularTrackSlot | null | undefined} slot */
  const circularTrackSlotDisplayLabel = (slot) => {
    if (isManagedConservationSlot(slot)) {
      return String(slot?.params?.label || slot?.params?.fileName || slot?.id || t('Conservation')).trim();
    }
    return circularTrackRendererLabel(/** @type {string} */ (slot?.renderer)); // no slot: no label, as the unkeyed label table gives none
  };
  /** @param {CircularTrackSlot | null | undefined} slot */
  const circularTrackSlotDisplayMeta = (slot) => {
    if (isManagedConservationSlot(slot)) {
      return String(slot?.params?.fileName || slot?.params?.series_key || '').trim();
    }
    return '';
  };
  /** @param {CircularTrackSlot | null | undefined} slot */
  const circularTrackSlotLegendLabelPlaceholder = (slot) => {
    const renderer = String(slot?.renderer || '').trim();
    if (renderer === 'dinucleotide_content' || renderer === 'dinucleotide_skew') {
      const nt = normalizeNt(slot?.params?.nt ?? slot?.params?.dinucleotide, circularDrawingDinucleotide());
      return renderer === 'dinucleotide_content' ? `${nt} content` : `${nt} skew`;
    }
    return t('Legend label');
  };
  /** @param {CircularTrackSlot | null | undefined} slot */
  const circularTrackSlotColor = (slot) => {
    if (isManagedConservationSlot(slot)) return String(slot?.params?.color || '').trim();
    return '';
  };

  /** @param {CircularTrackSlot | null | undefined} slot @param {string} key */
  const circularTrackSlotHasSkewColorOverride = (slot, key) => (
    slot?.renderer === 'dinucleotide_skew' &&
    normalizeOptionalText(slot?.params?.[key]) !== null
  );

  /** @param {CircularTrackSlot | null | undefined} slot @param {string} key */
  const circularTrackSlotSkewColorValue = (slot, key) => {
    const drawing = state.drawingOfMode('circular');
    return resolveTrackSlotSkewColorValue(/** @type {any} */ ({
      slot,
      key,
      currentColors: drawing.currentColors,
      paletteDefinitions: state.paletteDefinitions,
      selectedPalette: drawing.selectedPalette
    }));
  };

  /** @param {CircularTrackSlot | null | undefined} slot @param {string} key @param {unknown} value */
  const setCircularTrackSlotSkewColor = (slot, key, value) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!slot || slot.renderer !== 'dinucleotide_skew' || !['positive_color', 'negative_color'].includes(key)) return;
    slot.params = cloneParams(slot.params);
    const color = normalizeColorParam(value);
    if (color === null) delete slot.params[key];
    else slot.params[key] = color;
  };

  /** @param {CircularTrackSlot | null | undefined} slot @param {string} key */
  const clearCircularTrackSlotSkewColor = (slot, key) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!slot || slot.renderer !== 'dinucleotide_skew' || !['positive_color', 'negative_color'].includes(key)) return;
    slot.params = cloneParams(slot.params);
    delete slot.params[key];
  };

  const circularTrackSlots = () => {
    const drawing = state.drawingOfMode('circular');
    return (
      Array.isArray(drawing.adv.circular_track_slots)
        ? drawing.adv.circular_track_slots.map((/** @type {CircularTrackSlot} */ slot, /** @type {number} */ index) => ({ kind: STACK_ENTRY_SLOT, slot, index }))
        : []
    );
  };

  const circularTrackStackEntries = () => {
    const drawing = state.drawingOfMode('circular');
    /** @type {CircularTrackSlot[]} */
    const slots = Array.isArray(drawing.adv.circular_track_slots) ? drawing.adv.circular_track_slots : [];
    const axisIndex = axisIndexForCurrentSlots(drawing, slots);
    const entries = [];
    let axisRendered = false;
    slots.forEach((slot, index) => {
      const onAxis = effectiveSlotPlacement(slot, drawing.form.track_type) === 'overlay';
      if (index === axisIndex && !onAxis) {
        entries.push({ kind: STACK_ENTRY_AXIS, key: 'axis' });
        axisRendered = true;
      }
      if (onAxis && !axisRendered) {
        axisRendered = true;
      }
      entries.push({ kind: STACK_ENTRY_SLOT, slot, index, onAxis });
    });
    if (!axisRendered || axisIndex >= slots.length) {
      entries.push({ kind: STACK_ENTRY_AXIS, key: 'axis' });
    }
    return entries;
  };

  /** @param {CircularTrackSlot} slot */
  const circularTrackSlotCliSpec = (slot) => {
    const drawing = state.drawingOfMode('circular');
    normalizeSlotsInPlace();
    return buildCircularTrackSlotSpec(slot, drawing.adv.nt, drawing.form.track_type);
  };

  /** @param {CircularTrackSlot} slot @param {'width' | 'radius'} field @param {unknown} scalar */
  const updateCircularTrackSlotMeasure = (slot, field, scalar) => {
    const drawing = state.drawingOfMode('circular');
    if (!['width', 'radius'].includes(field) || !drawing.adv.circular_track_slots.includes(slot)) return;
    const numericLeaf = scalar && typeof scalar === 'object' && !Array.isArray(scalar)
      ? /** @type {{ value?: unknown }} */ (scalar).value : scalar; // a plain object past the guard
    if (typeof numericLeaf === 'number' && !Number.isFinite(numericLeaf)) {
      parseOptionalCircularScalar(scalar, `Circular track slot '${slot.id}' ${field}`);
    }
    if (slot[field] !== scalar) slot[field] = scalar;
  };

  /**
   * The gap fields' writer (on commit, so a half-typed `1.` is not rewritten).
   * @param {CircularTrackSlot} slot
   * @param {'inner_gap_px' | 'outer_gap_px'} field
   * @param {unknown} text
   */
  const updateCircularTrackSlotGap = (slot, field, text) => {
    if (!['inner_gap_px', 'outer_gap_px'].includes(field)) return;
    const normalized = normalizeGapText(typeof text === 'string' ? text.trim() : text, field);
    if (slot[field] !== normalized) slot[field] = /** @type {string | null} */ (normalized);
  };

  /** @param {CircularTrackSlot | null | undefined} slot @param {string} field */
  const circularSlotManualValue = (slot, field) => {
    if (!slot) return '';
    if (field === 'width') return slot.width;
    if (field === 'radius') return slot.radius;
    if (field === 'inner_gap_px') return slot.inner_gap_px;
    if (field === 'outer_gap_px') return slot.outer_gap_px;
    return '';
  };

  const selectedResultIndexValue = () => Number(state?.selectedResultIndex?.value ?? 0) || 0;

  /** @param {string | undefined} slotId */
  const circularSlotGeometryLookup = (slotId) => (/** @type {any} */ ({
    geometry: String(state?.trackSlotResolvedGeometry?.value?.mode || '') === 'circular'
      ? state.trackSlotResolvedGeometry.value
      : null,
    resultIndex: selectedResultIndexValue(),
    recordIndex: 0,
    slotId
  }));

  // Only a rendered row has resolved geometry. Any other row, and every row
  // before the first Generate, shows the Auto width and gaps it asks for as an
  // estimate (TK-15) and no radius: the layout places it (OV-1712).
  /** @param {CircularTrackSlot | null | undefined} slot @param {number} _slotIndex @param {string} field */
  const circularTrackSlotGeometryAutoText = (slot, _slotIndex, field) => {
    if (isManualSlotValue(circularSlotManualValue(slot, field))) return '';
    const lookup = circularSlotGeometryLookup(slot?.id);
    const resolved = circularTrackSlotEffectiveEnabled(slot) ? findTrackSlotGeometry(lookup) : null;
    if (!resolved) {
      const follows = field === 'width' || field === 'radius' ? overlayFollowsNote(slot) : '';
      if (follows) return follows;
      if (field === 'radius') return t('Auto');
      if (field === 'inner_gap_px' || field === 'outer_gap_px') return formatPxAuto(AUTO_SLOT_SIZES.gapPx, { estimate: true });
      if (field !== 'width') return '';
      const drawing = state.drawingOfMode('circular');
      const widthPx = autoWidthPx(String(slot?.renderer || '').trim(), getPreviewLengthParam(drawing, estimateRecordLengths()), drawing);
      return widthPx === null ? '' : formatPxAuto(widthPx, { estimate: true });
    }
    if (field === 'width') return formatPxAuto(resolved.widthPx);
    if (field === 'radius') {
      // `r` pins a ticks row's anchor, not the band centre the payload holds (GX-18).
      const tickAnchor = slot?.renderer === 'ticks'
        ? tickAnchorRadiusFactor(resolved, findTrackSlotGeometryRecord(lookup)?.axisRadiusPx, slot.params?.tick_label_layout)
        : null;
      return tickAnchor === null
        ? formatRadiusFactorAuto(resolved.radiusFactor)
        : formatRadiusFactorAuto(tickAnchor, { awayFromAxis: true });
    }
    if (field === 'inner_gap_px') return formatPxAuto(resolved.innerGapPx);
    if (field === 'outer_gap_px') return formatPxAuto(resolved.outerGapPx);
    return '';
  };

  /** @param {CircularTrackSlot | null | undefined} slot @param {string} field */
  const circularTrackSlotGeometryUnitSuffix = (slot, field) => {
    const text = String(circularSlotManualValue(slot, field) ?? '').trim();
    if (!text) return '';
    if (!['inner_gap_px', 'outer_gap_px'].includes(field)) return '';
    return /px$/i.test(text) ? '' : 'px';
  };

  /** @param {CircularTrackSlot | null | undefined} slot @param {string} field */
  const circularTrackSlotGeometryHasManual = (slot, field) => (
    isManualSlotValue(circularSlotManualValue(slot, field))
  );

  // What the preset draws where, in words; the radius ratios stay
  // with the slot fields' Auto values.
  const circularTrackPresetSummary = () => {
    const drawing = state.drawingOfMode('circular');
    const preset = normalizeCircularTrackPreset(drawing.form.track_type);
    const lane = laneDirectionForPreset(preset);
    const nt = circularDrawingDinucleotide();
    const sides = ownEntry(AUTO_SLOT_SIZES.presetNumericRowSides, preset);
    /** @type {string[]} */
    const inside = [];
    /** @type {string[]} */
    const outside = [];
    /** @param {string} name @param {string} rowId */
    const place = (name, rowId) => {
      const side = ownEntry(sides, rowId);
      if (side === 'inside') inside.push(name);
      else if (side === 'outside') outside.push(name);
    };
    if (Boolean(drawing.form.show_depth)) place(tm(RENDERER_LABELS.depth), 'depth');
    if (!Boolean(drawing.form.suppress_gc)) place(t('{nt} content', { nt }), 'gc_content');
    if (!Boolean(drawing.form.suppress_skew)) place(t('{nt} skew', { nt }), 'gc_skew');
    return {
      label: formatPresetName(preset),
      detail: [laneDirectionLabel(lane), inside.length ? t('{tracks} inside', { tracks: listText(inside) }) : '', outside.length ? t('{tracks} outside', { tracks: listText(outside) }) : '']
        .filter(Boolean).join(' · ')
    };
  };

  /** @param {CircularTrackSlot | null | undefined} slot */
  const circularTrackSlotUsesPresetGeometry = (slot) => {
    if (!slot || typeof slot !== 'object') return false;
    if (!SUPPORTED_RENDERERS.includes(slot.renderer)) return false;
    if (!slotHasManualGeometry(slot)) return true;
    if (normalizeOptionalPlacement(slot.side) === null) return true;
    if (slot.renderer === 'features' && normalizeOptionalText(slot.params?.lane_direction) === null) return true;
    if (slot.renderer === 'ticks') {
      return (
        normalizeOptionalText(slot.params?.tick_label_layout) === null
      );
    }
    return false;
  };

  return {
    circularTrackRenderers: UI_RENDERERS,
    circularTrackSlotEditorKey,
    circularTrackSlotOptionsOpen,
    setCircularTrackSlotOptionsOpen,
    circularTrackSlotOptionsSummary,
    circularDrawingDinucleotide,
    circularTrackRendererLabel,
    normalizeCircularTrackSlots: normalizeSlotsInPlace,
    syncCircularConservationSlots,
    changeCircularDepthSources,
    circularDepthSourceLoss,
    circularConservationRowsLeaving,
    ...featureSlotEdits(trackLayoutWithVisibility(changeTrackLayout, changeTrackVisibility), {
      resetCircularTrackSlotsFromSimpleControls,
      resetCircularTrackSlotsToPreset,
      setCircularTrackSlotsEnabled,
      addCircularTrackSlot,
      duplicateCircularTrackSlot,
      removeCircularTrackSlot,
      setCircularTrackSlotEnabled,
      moveCircularTrackSlot,
      moveCircularTrackSlotOutside,
      moveCircularTrackSlotInside,
      moveCircularTrackSlotToAxis,
      updateCircularTrackSlotRenderer,
      updateCircularTrackSlotPlacement,
      updateCircularTrackFeatureLane
    }),
    setCircularGcSuppressed,
    setCircularSkewSuppressed,
    setCircularDinucleotide,
    setCircularTrackSlotNt,
    canAddCircularTrackRenderer,
    circularTrackAddBlockReason,
    canDuplicateCircularTrackSlot,
    circularTrackSlotEffectiveEnabled,
    circularTrackSlotHiddenBySuppress,
    circularTrackSlotSuppressMessage,
    canMoveCircularTrackSlot,
    canMoveCircularTrackSlotOutside,
    canMoveCircularTrackSlotInside,
    canMoveCircularTrackSlotToAxis,
    updateCircularTrackSlotMeasure,
    updateCircularTrackSlotGap,
    circularTrackSlotIssue,
    circularTrackGlobalIssues,
    circularAnnotationAnchorOptions,
    circularAnnotationAnchorIsKnown,
    annotationTrackMarkOptions: ANNOTATION_MARK_OPTIONS,
    circularAnnotationMarkSelected,
    setCircularAnnotationMarkSelected,
    circularAnnotationLaneGapValue: (/** @type {CircularTrackSlot | null | undefined} */ slot) => circularAnnotationNumberValue(slot, 'lane_gap_px', 3),
    setCircularAnnotationLaneGap: (/** @type {CircularTrackSlot | null | undefined} */ slot, /** @type {unknown} */ value) => setCircularAnnotationNumber(slot, 'lane_gap_px', value, 3),
    circularAnnotationPaddingValue: (/** @type {CircularTrackSlot | null | undefined} */ slot) => circularAnnotationNumberValue(slot, 'padding_px', 2),
    setCircularAnnotationPadding: (/** @type {CircularTrackSlot | null | undefined} */ slot, /** @type {unknown} */ value) => setCircularAnnotationNumber(slot, 'padding_px', value, 2),
    circularAnnotationCoverAnchor,
    setCircularAnnotationCoverAnchor,
    setCircularOverlayAnnotationParam,
    supportsCircularTrackSlotPlacement,
    isManagedCircularConservationSlot,
    circularTrackSlots,
    circularTrackStackEntries,
    circularTrackSlotCliSpec,
    circularTrackSlotDisplayLabel,
    circularTrackSlotDisplayMeta,
    circularTrackSlotLegendLabelPlaceholder,
    circularTrackSlotColor,
    circularTrackSlotHasSkewColorOverride,
    circularTrackSlotSkewColorValue,
    circularTrackSlotGeometryAutoText,
    circularTrackSlotGeometryHasManual,
    circularTrackSlotGeometryUnitSuffix,
    setCircularTrackSlotSkewColor,
    clearCircularTrackSlotSkewColor,
    circularTrackPresetSummary,
    circularTrackSlotUsesPresetGeometry
  };
};
