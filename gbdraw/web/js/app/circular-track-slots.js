// @ts-check
/** @import { DrawingState } from '../state.js' */
/** @import { ChangeTrackLayout } from './track-slot-edits.js' */
import {
  CONSERVATION_SLOT_MANAGER,
  isManagedConservationSlot,
  normalizeFileList,
  orderedConservationSources,
  safeConservationSlotId
} from '../services/conservation-series.js';
import {
  isDefaultManagedDepthSlot,
  reconcileManagedDepthSlots,
  representativeDepthFiles
} from '../services/depth-track-state.js';
import { resolveTrackSlotSkewColorValue } from './track-slot-colors.js';
import {
  findTrackSlotGeometry,
  findTrackSlotGeometryRecord,
  formatPxAuto,
  formatRadiusFactorAuto,
  isManualSlotValue,
  normalizeOptionalText,
  tickAnchorRadiusFactor
} from '../services/track-slot-display.js';
import { featureSlotEdits } from './track-slot-edits.js';
import { parseOptionalCircularScalar, parseOptionalPixel, validateCustomTrackPlan } from '../services/track-slot-validation.js';
import { visibleFeatureUnderlaysForState } from '../utils/feature-rendering.js';
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
  tickLabelLayoutFromSides
} from '../services/circular-track-slot-model.js';

const UI_RENDERERS = SUPPORTED_RENDERERS.filter(
  (renderer) => renderer !== 'sequence_conservation'
);

const RENDERER_LABELS = {
  features: 'Features',
  ticks: 'Ticks',
  dinucleotide_content: 'Dinucleotide content',
  dinucleotide_skew: 'Dinucleotide skew',
  depth: 'Depth',
  sequence_conservation: 'Pairwise comparison',
  annotations: 'Annotations',
  spacer: 'Spacer'
};

const STACK_ENTRY_AXIS = 'axis';
const STACK_ENTRY_SLOT = 'slot';
const GLOBAL_SUPPRESS_PARAM = '_suppressed_by_global';
const SUPPRESS_RENDERER_BY_KEY = {
  gc_content: 'dinucleotide_content',
  gc_skew: 'dinucleotide_skew'
};
const SUPPRESS_KEY_BY_RENDERER = {
  dinucleotide_content: 'gc_content',
  dinucleotide_skew: 'gc_skew'
};
const SUPPRESS_FORM_KEY_BY_RENDERER = {
  dinucleotide_content: 'suppress_gc',
  dinucleotide_skew: 'suppress_skew'
};
const SUPPRESS_TRACK_LABEL_BY_RENDERER = {
  dinucleotide_content: 'GC content',
  dinucleotide_skew: 'GC skew'
};
const SUPPRESS_CONTROL_LABEL_BY_RENDERER = {
  dinucleotide_content: 'Hide GC Content',
  dinucleotide_skew: 'Hide GC Skew'
};
const PRESET_LABELS = {
  tuckin: 'Tuckin',
  middle: 'Middle',
  spreadout: 'Spreadout'
};
const PREVIEW_LENGTH_THRESHOLD_BP = 50000;
const PREVIEW_RADIUS_PX = 390;
const PREVIEW_TRACK_RATIO = 0.19;
const PREVIEW_TRACK_DICT = {
  short: {
    spreadout: { 1: 1.0, 2: 0.85, 3: 0.65, 4: 0.45 },
    middle: { 1: 1.0, 2: 0.75, 3: 0.55, 4: 0.35 },
    tuckin: { 1: 1.0, 2: 0.64, 3: 0.44, 4: 0.24 }
  },
  long: {
    spreadout: { 1: 1.0, 2: 0.80, 3: 0.60, 4: 0.40 },
    middle: { 1: 1.0, 2: 0.75, 3: 0.55, 4: 0.35 },
    tuckin: { 1: 1.0, 2: 0.70, 3: 0.50, 4: 0.30 }
  }
};
const PREVIEW_TRACK_RATIO_FACTORS = {
  short: [0.50, 1.0, 1.0],
  long: [0.25, 1.0, 1.0]
};
const ANNOTATION_MARK_OPTIONS = Object.freeze([
  'line',
  'bracket',
  'band',
  'highlight'
]);

const formatPresetName = (preset) => PRESET_LABELS[normalizeCircularTrackPreset(preset)] || PRESET_LABELS.tuckin;

const laneDirectionLabel = (laneDirection) => {
  const lane = normalizeLaneDirection(laneDirection);
  if (lane === 'outside') return 'feature stack outside axis';
  if (lane === 'split') return 'feature stack centered on axis';
  return 'feature stack inside axis';
};

const previewWidthPxForRenderer = (renderer, lengthParam) => {
  const base = PREVIEW_RADIUS_PX * PREVIEW_TRACK_RATIO;
  const factors = PREVIEW_TRACK_RATIO_FACTORS[lengthParam] || PREVIEW_TRACK_RATIO_FACTORS.long;
  if (renderer === 'features') return base * Number(factors[0]);
  if (renderer === 'sequence_conservation') return base * Number(factors[0]);
  if (renderer === 'depth') return base * Number(factors[1]) * 0.5;
  if (renderer === 'dinucleotide_skew') return base * Number(factors[2]);
  // An Auto ticks row draws marks of the default length max(6 px, 0.025 R)
  // (_default_tick_length_px, gbdraw/svg/circular_ticks.py), not 0 px (TK-15).
  if (renderer === 'ticks') return Math.max(6, 0.025 * PREVIEW_RADIUS_PX);
  return base * Number(factors[1]);
};

const previewSpacingPx = () => Math.max(1.0, 0.01 * PREVIEW_RADIUS_PX);

/** @param {DrawingState} drawing */
const previewFeatureLaneCount = (drawing) => {
  if (Boolean(drawing?.form?.separate_strands)) return 2;
  return 1;
};

/** @param {DrawingState} drawing */
const previewFeatureRadiusRatio = (preset, lengthParam, drawing) => {
  const normalized = normalizeCircularTrackPreset(preset);
  const laneWidth = previewWidthPxForRenderer('features', lengthParam);
  const laneCount = previewFeatureLaneCount(drawing);
  const spacing = previewSpacingPx();
  const bandWidth = (laneCount * laneWidth) + (Math.max(0, laneCount - 1) * spacing);
  if (normalized === 'tuckin') {
    return (PREVIEW_RADIUS_PX - spacing - (bandWidth / 2)) / PREVIEW_RADIUS_PX;
  }
  if (normalized === 'spreadout') {
    return (PREVIEW_RADIUS_PX + spacing + (bandWidth / 2)) / PREVIEW_RADIUS_PX;
  }
  return 1.0;
};

const getPreviewRecordEntries = (state) => {
  const recordsRef = state?.circularRecordList;
  const records = Array.isArray(recordsRef?.value)
    ? recordsRef.value
    : (Array.isArray(recordsRef) ? recordsRef : []);
  return records;
};

const getPreviewLengthParam = (state) => {
  const lengths = getPreviewRecordEntries(state)
    .map((entry) => Number(entry?.record_length ?? entry?.length ?? 0))
    .filter((value) => Number.isFinite(value) && value > 0);
  if (lengths.length === 0) return 'long';
  return Math.max(...lengths) < PREVIEW_LENGTH_THRESHOLD_BP ? 'short' : 'long';
};

/** @param {DrawingState} drawing */
const getBuiltinTrackId = (slot, renderer, drawing) => {
  const id = String(slot?.id || '').trim();
  const showDepth = Boolean(drawing?.form?.show_depth);
  const showGc = !Boolean(drawing?.form?.suppress_gc);
  const showSkew = !Boolean(drawing?.form?.suppress_skew);

  if (renderer === 'depth' && id === 'depth' && showDepth) return 2;
  if (renderer === 'dinucleotide_content' && id === 'gc_content' && showGc) {
    return showDepth ? 3 : 2;
  }
  if (renderer === 'dinucleotide_skew' && id === 'gc_skew' && showSkew) {
    if (showDepth) return showGc ? 4 : 3;
    return showGc ? 3 : 2;
  }
  return null;
};

/** @param {DrawingState} drawing */
const getPresetRadiusRatio = (slot, renderer, preset, lengthParam, drawing) => {
  if (renderer === 'features' && String(slot?.id || '').trim() === 'features') {
    return previewFeatureRadiusRatio(preset, lengthParam, drawing);
  }
  const trackId = getBuiltinTrackId(slot, renderer, drawing);
  if (trackId === null) return null;
  return PREVIEW_TRACK_DICT[lengthParam]?.[normalizeCircularTrackPreset(preset)]?.[trackId] ?? null;
};

const slotHasManualGeometry = (slot) => (
  normalizeOptionalText(slot?.width) !== null ||
  normalizeOptionalText(slot?.radius) !== null ||
  normalizeOptionalText(slot?.inner_gap_px) !== null ||
  normalizeOptionalText(slot?.outer_gap_px) !== null
);

const circularAvailableDepthTrackCountForState = (state) => {
  const files = representativeDepthFiles(state?.files?.c_depth);
  return files.some(Boolean) ? files.length : 0;
};

const circularSourcedDepthTrackIndexesForState = (state) => (
  representativeDepthFiles(state?.files?.c_depth)
    .flatMap((file, trackIndex) => (file ? [trackIndex] : []))
);

/** @param {DrawingState} drawing */
const circularDepthTrackCountForState = (state, drawing) => (
  Boolean(drawing?.form?.show_depth)
    ? circularAvailableDepthTrackCountForState(state)
    : 0
);

const applyPlacementDefaults = (slot, placement = 'inside') => {
  if (!slot) return;
  const requestedSide = normalizePlacement(placement);
  slot.side = requestedSide;
  slot.params = cloneParams(slot.params);
  if (slot.renderer === 'features') {
    slot.params.lane_direction = laneDirectionForSide(requestedSide);
  }
};

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
 * @param {Record<string, any>[]} [existingSlots]
 * @param {string} [nt]
 * @param {Record<string, any> | null} [placement]
 * @returns {Record<string, any>}
 */
export const createCircularTrackSlotForRenderer = (renderer, existingSlots = [], nt = 'GC', placement = null) => {
  const normalizedRenderer = SUPPORTED_RENDERERS.includes(renderer) ? renderer : 'dinucleotide_skew';
  const baseId = DEFAULT_SLOT_IDS[normalizedRenderer] || normalizedRenderer;
  const existingIds = new Set(
    (Array.isArray(existingSlots) ? existingSlots : [])
      .map((slot) => String(slot?.id || '').trim())
      .filter(Boolean)
  );
  let id = baseId;
  let suffix = 2;
  while (existingIds.has(id)) {
    id = `${baseId}_${suffix}`;
    suffix += 1;
  }

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

const appendOption = (options, key, value) => {
  const text = normalizeOptionalText(value);
  if (text === null) return;
  options.push(`${key}=${text}`);
};

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
    appendOption(options, 'layer', params.layer);
    if (params.cover_anchor === true) {
      appendOption(options, 'cover_anchor', 'true');
    }
  }
  appendOption(options, 'legend_label', params.legend_label);

  return `${normalized.id}:${normalized.renderer}${options.length ? `@${options.join(',')}` : ''}`;
};

export const hasEnabledCircularTrackRenderer = (slots, renderer) =>
  normalizeCircularTrackSlots(slots).some((slot) => slot.enabled && slot.renderer === renderer);

export const isCircularTrackRendererSuppressedByForm = (renderer, form = {}) => {
  const normalizedRenderer = String(renderer || '').trim();
  const formKey = SUPPRESS_FORM_KEY_BY_RENDERER[normalizedRenderer];
  return Boolean(formKey && form?.[formKey]);
};

export const circularTrackSlotHiddenBySuppressForm = (slot, form = {}) =>
  isCircularTrackRendererSuppressedByForm(slot?.renderer, form);

// Hide GC Content / Hide GC Skew own the visibility of their rows. A hidden
// renderer disables its enabled rows and marks them; once the form no longer
// hides it, the marked rows are enabled again. Rows the user disabled carry no
// mark and stay disabled. The projection depends only on the rows and the form.
export const applyCircularSuppressControlsToSlots = (slots, form = {}) => (
  (Array.isArray(slots) ? slots : []).map((slot) => {
    if (!slot || typeof slot !== 'object' || Array.isArray(slot)) return slot;
    const token = SUPPRESS_KEY_BY_RENDERER[String(slot.renderer || '').trim()];
    if (!token) return slot;
    const params = cloneParams(slot.params);
    if (circularTrackSlotHiddenBySuppressForm(slot, form)) {
      if (slot.enabled === false) return slot;
      params[GLOBAL_SUPPRESS_PARAM] = token;
      return { ...slot, enabled: false, params };
    }
    if (params[GLOBAL_SUPPRESS_PARAM] !== token) return slot;
    delete params[GLOBAL_SUPPRESS_PARAM];
    return { ...slot, enabled: true, params };
  })
);

/** @param {DrawingState} drawing */
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

/** @param {DrawingState} drawing */
const conservationEntriesForState = (state, drawing) => {
  if (drawing?.circularConservation?.enabled !== true) return [];
  return orderedConservationSources(
    conservationSourceFilesForState(state, drawing),
    drawing.circularConservation
  );
};

const positiveNumberOrNull = (value) => {
  if (value === null || value === undefined || value === '') return null;
  const parsed = Number(value);
  return Number.isFinite(parsed) && parsed > 0 ? parsed : null;
};

/** @param {DrawingState} drawing */
export const estimateCircularConservationLayoutWarning = (state, drawing) => {
  const currentMode = String(state?.mode?.value ?? state?.mode ?? '').trim().toLowerCase();
  if (currentMode && currentMode !== 'circular') return '';
  if (drawing?.circularConservation?.enabled !== true) return '';

  const entries = conservationEntriesForState(state, drawing);
  if (entries.length <= 0) return '';

  const preset = normalizeCircularTrackPreset(drawing?.form?.track_type);
  const lengthParam = getPreviewLengthParam(state);
  const axisRadius = PREVIEW_RADIUS_PX;
  const defaultSpacing = previewSpacingPx();
  const lane = laneDirectionForPreset(preset);
  const featureLaneCount = previewFeatureLaneCount(drawing);
  const featureWidth = previewWidthPxForRenderer('features', lengthParam);
  const featureBandWidth = (featureLaneCount * featureWidth) + (Math.max(0, featureLaneCount - 1) * defaultSpacing);
  let availableInsidePx = axisRadius - defaultSpacing;
  if (lane === 'inside') {
    const featureCenter = previewFeatureRadiusRatio(preset, lengthParam, drawing) * axisRadius;
    availableInsidePx = featureCenter - (featureBandWidth / 2) - defaultSpacing;
  } else if (lane === 'split') {
    availableInsidePx = axisRadius - (featureBandWidth / 2) - defaultSpacing;
  }
  availableInsidePx = Math.max(0, availableInsidePx);

  const ringWidth = positiveNumberOrNull(drawing?.circularConservation?.ring_width)
    ?? previewWidthPxForRenderer('sequence_conservation', lengthParam);
  const ringGap = positiveNumberOrNull(drawing?.circularConservation?.ring_gap) ?? defaultSpacing;
  const showGc = !Boolean(drawing?.form?.suppress_gc);
  const showSkew = !Boolean(drawing?.form?.suppress_skew);
  const showDepth = Boolean(drawing?.form?.show_depth);
  const numericAfterConservationPx =
    (showDepth ? previewWidthPxForRenderer('depth', lengthParam) + defaultSpacing : 0) +
    (showGc ? previewWidthPxForRenderer('dinucleotide_content', lengthParam) + defaultSpacing : 0) +
    (showSkew ? previewWidthPxForRenderer('dinucleotide_skew', lengthParam) + defaultSpacing : 0);
  const requestedStackPx =
    (entries.length * ringWidth) +
    (Math.max(0, entries.length - 1) * ringGap) +
    numericAfterConservationPx;
  const compressedStackPx =
    (entries.length * 4) +
    (Math.max(0, entries.length - 1) * 1) +
    (showDepth ? 10 + 1 : 0) +
    (showGc ? 10 + 1 : 0) +
    (showSkew ? 12 + 1 : 0);

  if (compressedStackPx > availableInsidePx) {
    return `There are ${entries.length} inside pairwise comparison ring(s) plus other circular tracks; even compressed rings may not fit. Reduce Ring Width/GAP, disable GC/skew/depth tracks, or move tracks outside before generating.`;
  }
  if (requestedStackPx > availableInsidePx * 0.9 || (entries.length >= 5 && (showGc || showSkew || showDepth))) {
    return `There are ${entries.length} inside pairwise comparison ring(s); gbdraw will auto-compress ring width/gap when needed. If generation still fails, reduce Ring Width/GAP or disable GC/skew/depth tracks.`;
  }
  return '';
};

const managedConservationSlotKey = (slot) => String(slot?.params?.series_key || '').trim();

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

const replaceObjectContents = (target, source) => {
  if (
    !target || typeof target !== 'object' || Array.isArray(target) ||
    !source || typeof source !== 'object' || Array.isArray(source)
  ) {
    return source;
  }
  Object.keys(target).forEach((key) => {
    if (!Object.prototype.hasOwnProperty.call(source, key)) delete target[key];
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
 */

// `changeTrackLayout` is the feature placement owner's transition, injected as
// a port (R10, Q3, R13): every stack edit that can change the feature slot
// runs through it.
/**
 * @param {CircularTrackSlotEditorOptions} options
 */
export const createCircularTrackSlotEditor = ({ state, changeTrackLayout = (apply) => apply() }) => {
  const editorKeys = new WeakMap();
  let nextEditorKey = 1;
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

  /** @param {DrawingState} drawing */
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
  const desiredCircularDepthTrackCount = (drawing) => circularDepthTrackCountForState(state, drawing);

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

  const circularTrackSlotIssue = (slot, index = null) => {
    const drawing = state.drawings.circular;
    const resolvedIndex = Number.isInteger(Number(index))
      ? Number(index)
      : drawing.adv.circular_track_slots.findIndex((candidate) => candidate === slot);
    if (resolvedIndex < 0) return '';
    return (circularTrackValidationPlan(drawing).rowIssues.get(resolvedIndex) || [])
      // The Width and Radius fields show their own invalid value (TK-12).
      .filter((issue) => !(issue.code === 'geometry_invalid' && ['width', 'radius'].includes(issue.field)))
      .map((issue) => issue.message)
      .join(' ');
  };

  const circularTrackGlobalIssues = () => {
    const drawing = state.drawings.circular;
    return (
      circularTrackValidationPlan(drawing).globalIssues.map((issue) => issue.message)
    );
  };

  const circularAnnotationAnchorOptions = (slot = null) => (
    state.drawings.circular.adv.circular_track_slots
      .filter((candidate) => (
        candidate &&
        candidate !== slot &&
        candidate.enabled !== false &&
        !['annotations', 'ticks', 'spacer'].includes(String(candidate.renderer || '').trim()) &&
        String(candidate.id || '').trim()
      ))
      .map((candidate) => ({
        id: String(candidate.id).trim(),
        label: `${String(candidate.id).trim()} · ${circularTrackRendererLabel(candidate.renderer)}`
      }))
  );

  const circularAnnotationAnchorIsKnown = (slot) => {
    const anchor = String(slot?.params?.anchor_slot || '').trim();
    return !anchor || circularAnnotationAnchorOptions(slot).some((option) => option.id === anchor);
  };

  const bindOnlyCircularAnnotationAnchor = (slot) => {
    if (!slot || slot.renderer !== 'annotations') return;
    slot.params = cloneParams(slot.params);
    const options = circularAnnotationAnchorOptions(slot);
    const current = String(slot.params.anchor_slot || '').trim();
    if (options.some((option) => option.id === current)) return;
    if (options.length === 1) slot.params.anchor_slot = options[0].id;
    else delete slot.params.anchor_slot;
  };

  const canAddCircularTrackRenderer = (renderer) => {
    const drawing = state.drawings.circular;
    const normalizedRenderer = String(renderer || '').trim();
    if (!UI_RENDERERS.includes(normalizedRenderer)) return false;
    if (normalizedRenderer === 'annotations') return annotationSetIds(drawing).length > 0;
    if (normalizedRenderer === 'depth') {
      return circularAvailableDepthTrackCountForState(state) > 0;
    }
    if (
      normalizedRenderer === 'features' &&
      visibleFeatureUnderlaysForState(state).length > 0
    ) {
      return !drawing.adv.circular_track_slots.some((slot) => (
        slot?.enabled !== false && slot?.renderer === 'features'
      ));
    }
    return true;
  };

  const canDuplicateCircularTrackSlot = (slot) => {
    const drawing = state.drawings.circular;
    if (!slot || isManagedConservationSlot(slot)) return false;
    if (slot.enabled === false) return true;
    if (
      slot.renderer === 'features' &&
      visibleFeatureUnderlaysForState(state).length > 0
    ) {
      return false;
    }
    if (slot.renderer === 'annotations') return annotationSetIds(drawing).length > 0;
    if (slot.renderer === 'depth') return circularAvailableDepthTrackCountForState(state) > 0;
    return true;
  };

  const circularAnnotationMarkSelected = (slot, mark) => {
    const selected = Array.isArray(slot?.params?.marks)
      ? slot.params.marks.map((value) => String(value).trim().toLowerCase()).filter(Boolean)
      : [];
    return selected.length === 0 || selected.includes(mark);
  };

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

  const circularAnnotationNumberValue = (slot, field, defaultValue) => {
    const raw = slot?.params?.[field];
    return raw === null || raw === undefined || raw === '' ? defaultValue : raw;
  };

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

  const circularAnnotationCoverAnchor = (slot) => slot?.params?.cover_anchor === true;

  const setCircularAnnotationCoverAnchor = (slot, checked) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!slot || slot.renderer !== 'annotations') return;
    slot.params = cloneParams(slot.params);
    if (checked) slot.params.cover_anchor = true;
    else delete slot.params.cover_anchor;
  };

  const makeDepthSlotForTrackIndex = (trackIndex, existingIds, desiredCount) => {
    const normalizedIndex = Math.max(0, Number(trackIndex) || 0);
    const preferredId = Number(desiredCount) > 1 ? `depth_${normalizedIndex + 1}` : 'depth';
    let id = preferredId;
    if (existingIds.has(id) && normalizedIndex === 0 && !existingIds.has('depth')) {
      id = 'depth';
    }
    let suffix = 2;
    while (existingIds.has(id)) {
      id = `${preferredId}_${suffix}`;
      suffix += 1;
    }
    existingIds.add(id);
    return makeSlot({
      id,
      renderer: 'depth',
      params: {
        track_index: normalizedIndex
      }
    });
  };

  const normalizeSlotsInPlace = () => {
    const drawing = state.drawings.circular;
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
    const suppressed = applyCircularSuppressControlsToSlots(normalized, drawing.form);
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
   * @param {{ nextSlots: any[], newSlots?: any[], managedPredicate: (slot: any) => boolean, preferredInsertIndex?: number | null }} mutation
   */
  const commitManagedSlotMutation = (drawing, {
    nextSlots,
    newSlots = [],
    managedPredicate,
    preferredInsertIndex = null
  }) => {
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
    resetCircularTrackSlotsToPreset(state.drawings.circular.form.track_type)
  );

  // Managed Depth rows follow Depth sources (PD-OI-058). The saved stack is
  // reconciled too, so a later Use custom stack shows the row; an empty saved
  // stack is built by Reset when the stack is enabled.
  /** @param {DrawingState} drawing */
  const reconcileCircularDepthSlots = (drawing, previousSourced) => {
    const slots = Array.isArray(drawing.adv.circular_track_slots) ? drawing.adv.circular_track_slots : [];
    if (slots.length === 0 && !drawing.adv.circular_track_slots_enabled) return;
    const { slots: nextSlots, additions } = reconcileManagedDepthSlots(/** @type {any} */ ({
      slots,
      previousSourced,
      sourced: circularSourcedDepthTrackIndexesForState(state),
      managedPredicate: isDefaultManagedDepthSlot
    }));
    if (additions.length === 0 && nextSlots.length === slots.length) return;
    const seriesCount = circularAvailableDepthTrackCountForState(state);
    const existingIds = new Set(nextSlots.map((slot) => String(slot?.id || '').trim()).filter(Boolean));
    const newSlots = applyCircularGeometryShortcuts(
      additions.map((trackIndex) => makeDepthSlotForTrackIndex(trackIndex, existingIds, seriesCount)),
      circularGeometryShortcutsForState(drawing)
    );
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
  const changeCircularDepthSources = (mutate) => {
    const drawing = state.drawings.circular;
    const previousSourced = circularSourcedDepthTrackIndexesForState(state);
    mutate();
    reconcileCircularDepthSlots(drawing, previousSourced);
    if (previousSourced.length > 0 && circularSourcedDepthTrackIndexesForState(state).length === 0) drawing.form.show_depth = false;
  };

  const syncCircularConservationSlots = () => {
    const drawing = state.drawings.circular;
    const entries = conservationEntriesForState(state, drawing);
    const slots = Array.isArray(drawing.adv.circular_track_slots)
      ? drawing.adv.circular_track_slots
      : [];
    const desiredKeys = new Set(entries.map((entry) => String(entry.sourceKey || '')));
    const existingIds = new Set(slots.map((slot) => String(slot?.id || '').trim()).filter(Boolean));
    const entryByKey = new Map(entries.map((entry) => [String(entry.sourceKey || ''), entry]));
    const entryOrderByKey = new Map(entries.map((entry, index) => [String(entry.sourceKey || ''), index]));

    const nextSlots = [];
    const presentKeys = new Set();
    slots.forEach((slot) => {
      if (!isManagedConservationSlot(slot)) {
        nextSlots.push(slot);
        return;
      }
      const key = managedConservationSlotKey(slot);
      if (!key || !desiredKeys.has(key) || presentKeys.has(key)) return;
      refreshManagedConservationSlot(slot, entryByKey.get(key), entryOrderByKey.get(key) ?? presentKeys.size);
      presentKeys.add(key);
      nextSlots.push(slot);
    });

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

  const resetCircularTrackSlotsToPreset = (preset) => {
    const drawing = state.drawings.circular;
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const normalizedPreset = normalizeCircularTrackPreset(preset);
    const templateSlots = applyCircularGeometryShortcuts(createDefaultCircularTrackSlots({
      nt: drawing.adv.nt,
      showDepth: Boolean(drawing.form.show_depth),
      depthTrackCount: desiredCircularDepthTrackCount(drawing),
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

  const setCircularTrackSlotsEnabled = (enabled) => {
    const drawing = state.drawings.circular;
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    drawing.adv.circular_track_slots_enabled = Boolean(enabled);
    if (
      drawing.adv.circular_track_slots_enabled &&
      (!Array.isArray(drawing.adv.circular_track_slots) || drawing.adv.circular_track_slots.length === 0)
    ) {
      resetCircularTrackSlotsFromSimpleControls();
    }
  };

  const addCircularTrackSlot = (renderer, placement = null) => {
    const drawing = state.drawings.circular;
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
          .filter((candidate) => candidate?.enabled !== false && candidate?.renderer === 'depth')
          .map((candidate) => normalizeTrackIndex(candidate?.params?.track_index))
          .filter((trackIndex) => trackIndex !== null)
      );
      let trackIndex = 0;
      while (trackIndex < available && claimed.has(trackIndex)) trackIndex += 1;
      slot.params.track_index = trackIndex < available ? trackIndex : 0;
    }
    drawing.adv.circular_track_slots.push(slot);
    normalizeSlotsInPlace();
  };

  const duplicateCircularTrackSlot = (index) => {
    const drawing = state.drawings.circular;
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    normalizeSlotsInPlace();
    const idx = Number(index);
    if (!Number.isInteger(idx) || idx < 0 || idx >= drawing.adv.circular_track_slots.length) return;
    if (!canDuplicateCircularTrackSlot(drawing.adv.circular_track_slots[idx])) return;
    const source = normalizeCircularTrackSlot(drawing.adv.circular_track_slots[idx], idx, drawing.adv.nt, drawing.form.track_type);
    const duplicate = createCircularTrackSlotForRenderer(source.renderer, drawing.adv.circular_track_slots, drawing.adv.nt);
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

  const removeCircularTrackSlot = (index) => {
    const drawing = state.drawings.circular;
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const idx = Number(index);
    if (!Number.isInteger(idx) || idx < 0 || idx >= drawing.adv.circular_track_slots.length) return;
    const axis = axisIndexForCurrentSlots(drawing, drawing.adv.circular_track_slots);
    drawing.adv.circular_track_slots.splice(idx, 1);
    drawing.adv.circular_track_slots_axis_index = idx < axis ? Math.max(0, axis - 1) : Math.min(axis, drawing.adv.circular_track_slots.length);
    normalizeSlotsInPlace();
  };

  /** @param {DrawingState} drawing */
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

  const moveCircularTrackSlot = (fromIndex, toIndex) => {
    const drawing = state.drawings.circular;
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

  const canMoveCircularTrackSlot = (index, direction) => {
    const drawing = state.drawings.circular;
    const idx = Number(index);
    const step = Number(direction);
    if (!Number.isInteger(idx) || !Number.isInteger(step) || step === 0) return false;
    const target = idx + Math.sign(step);
    return !wouldCircularTrackSlotMoveCrossAxis(drawing, idx, target);
  };

  const canMoveCircularTrackSlotOutside = (index) => {
    const drawing = state.drawings.circular;
    const idx = Number(index);
    const normalized = normalizedSlotsForCurrentState(drawing);
    if (!Number.isInteger(idx) || idx < 0 || idx >= normalized.length) return false;
    if (effectiveSlotPlacement(normalized[idx], drawing.form.track_type) === 'overlay') return true;
    return idx >= axisIndexForCurrentSlots(drawing, normalized);
  };

  const canMoveCircularTrackSlotInside = (index) => {
    const drawing = state.drawings.circular;
    const idx = Number(index);
    const normalized = normalizedSlotsForCurrentState(drawing);
    if (!Number.isInteger(idx) || idx < 0 || idx >= normalized.length) return false;
    if (effectiveSlotPlacement(normalized[idx], drawing.form.track_type) === 'overlay') return true;
    return idx < axisIndexForCurrentSlots(drawing, normalized);
  };

  const canMoveCircularTrackSlotToAxis = (index) => {
    const drawing = state.drawings.circular;
    const idx = Number(index);
    const normalized = normalizedSlotsForCurrentState(drawing);
    if (!Number.isInteger(idx) || idx < 0 || idx >= normalized.length) return false;
    const slot = normalized[idx];
    return ['features', 'ticks', 'annotations'].includes(slot?.renderer) && effectiveSlotPlacement(slot, drawing.form.track_type) !== 'overlay';
  };

  /** @param {DrawingState} drawing */
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
      delete directOverlaySlot.params.anchor_slot;
      delete directOverlaySlot.params.cover_anchor;
    }
    if (targetPlacement === 'overlay' && directOverlaySlot?.renderer === 'annotations') {
      syncSlotPlacementFromSide(directOverlaySlot, 'overlay');
      directOverlaySlot.params = cloneParams(directOverlaySlot.params);
      bindOnlyCircularAnnotationAnchor(directOverlaySlot);
      normalizeSlotsInPlace();
      return;
    }

    if (targetPlacement === 'overlay') {
      const movedSlot = drawing.adv.circular_track_slots[idx];
      if (!movedSlot) return;
      const movedPreviousPlacement = effectiveSlotPlacement(movedSlot, drawing.form.track_type);
      const existingAxisIndex = drawing.adv.circular_track_slots.findIndex((slot, slotIndex) => (
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
      const onAxisIndex = drawing.adv.circular_track_slots.findIndex((candidate) => (
        effectiveSlotPlacement(candidate, drawing.form.track_type) === 'overlay'
      ));
      const insertIndex = onAxisIndex >= 0 ? onAxisIndex + 1 : axis;
      drawing.adv.circular_track_slots.splice(insertIndex, 0, slot);
      drawing.adv.circular_track_slots_axis_index = onAxisIndex >= 0 ? onAxisIndex : axis;
    }
    normalizeSlotsInPlace();
  };

  const moveCircularTrackSlotOutside = (index) => {
    const drawing = state.drawings.circular;
    if (!canMoveCircularTrackSlotOutside(index)) return;
    moveCircularTrackSlotToPlacement(drawing, index, 'outside');
  };

  const moveCircularTrackSlotInside = (index) => {
    const drawing = state.drawings.circular;
    if (!canMoveCircularTrackSlotInside(index)) return;
    moveCircularTrackSlotToPlacement(drawing, index, 'inside');
  };

  const moveCircularTrackSlotToAxis = (index) => {
    const drawing = state.drawings.circular;
    if (!canMoveCircularTrackSlotToAxis(index)) return;
    moveCircularTrackSlotToPlacement(drawing, index, 'overlay');
  };

  const updateCircularTrackSlotRenderer = (slot, renderer) => {
    const drawing = state.drawings.circular;
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    renderer = renderer || slot?.renderer;
    if (!slot || !SUPPORTED_RENDERERS.includes(renderer)) return;
    slot.renderer = renderer;
    slot.params = cloneParams(slot.params);
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
      slot.params.nt = normalizeNt(slot.params.nt, normalizeNt(drawing.adv.nt));
    } else if (renderer === 'annotations') {
      slot.side = slot.side === 'overlay' ? 'overlay' : 'outside';
      slot.params.set_id = String(slot.params.set_id || drawing.annotationSets?.[0]?.id || '');
      slot.params.overflow = 'error';
      slot.params.show_labels = true;
      slot.params.layer = 'foreground';
    }
    normalizeSlotsInPlace();
  };

  /** @param {DrawingState} drawing */
  const activeCircularTrackSlotsForRenderer = (drawing, renderer) => {
    normalizeSlotsInPlace();
    return drawing.adv.circular_track_slots.filter((slot) => (
      slot &&
      slot.enabled !== false &&
      String(slot.renderer || '').trim() === renderer
    ));
  };

  /** @param {DrawingState} drawing */
  const confirmCircularSuppressOverride = (drawing, renderer) => {
    if (!drawing.adv.circular_track_slots_enabled) return true;
    const activeSlots = activeCircularTrackSlotsForRenderer(drawing, renderer);
    if (activeSlots.length === 0) return true;

    const trackLabel = SUPPRESS_TRACK_LABEL_BY_RENDERER[renderer] || 'selected';
    const plural = activeSlots.length === 1 ? '' : 's';
    const message = [
      `Custom Track Slots currently include ${activeSlots.length} enabled ${trackLabel} track${plural}.`,
      `Hiding ${trackLabel} will disable those custom track slot${plural} and override the custom track settings.`,
      '',
      'Continue?'
    ].join('\n');
    return globalThis.confirm ? globalThis.confirm(message) : true;
  };

  /** @typedef {{ target?: { checked: boolean } | null } | null} CircularSuppressToggleEvent The checkbox change event of a Suppress control. */

  /**
   * @param {DrawingState} drawing
   * @param {CircularSuppressToggleEvent} [event]
   */
  const setCircularSuppressControl = (drawing, key, checked, event = null) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const renderer = SUPPRESS_RENDERER_BY_KEY[key];
    if (!renderer) return;
    const formKey = SUPPRESS_FORM_KEY_BY_RENDERER[renderer];
    const nextChecked = Boolean(checked);
    const previousChecked = Boolean(drawing.form?.[formKey]);

    if (nextChecked === previousChecked) {
      if (event?.target) event.target.checked = previousChecked;
      return;
    }

    if (nextChecked && !confirmCircularSuppressOverride(drawing, renderer)) {
      if (event?.target) event.target.checked = previousChecked;
      return;
    }
    drawing.form[formKey] = nextChecked;
    normalizeSlotsInPlace();

    if (event?.target) event.target.checked = nextChecked;
  };

  /** @param {CircularSuppressToggleEvent} [event] */
  const setCircularGcSuppressed = (checked, event = null) => {
    const drawing = state.drawings.circular;
    setCircularSuppressControl(drawing, 'gc_content', checked, event);
  };

  /** @param {CircularSuppressToggleEvent} [event] */
  const setCircularSkewSuppressed = (checked, event = null) => {
    const drawing = state.drawings.circular;
    setCircularSuppressControl(drawing, 'gc_skew', checked, event);
  };

  const setCircularTrackSlotEnabled = (slot, enabled) => {
    const drawing = state.drawings.circular;
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!slot || circularTrackSlotHiddenBySuppressForm(slot, drawing.form)) return;
    slot.enabled = Boolean(enabled);
    if (slot.params && typeof slot.params === 'object' && !Array.isArray(slot.params)) {
      const token = SUPPRESS_KEY_BY_RENDERER[String(slot.renderer || '').trim()];
      if (token && slot.params[GLOBAL_SUPPRESS_PARAM] === token) {
        slot.params = cloneParams(slot.params);
        delete slot.params[GLOBAL_SUPPRESS_PARAM];
      }
    }
    normalizeSlotsInPlace();
  };

  const circularTrackSlotEffectiveEnabled = (slot) => (
    Boolean(slot?.enabled !== false) && !circularTrackSlotHiddenBySuppressForm(slot, state.drawings.circular.form)
  );

  const circularTrackSlotHiddenBySuppress = (slot) =>
    circularTrackSlotHiddenBySuppressForm(slot, state.drawings.circular.form);

  const circularTrackSlotSuppressMessage = (slot) => {
    if (!circularTrackSlotHiddenBySuppress(slot)) return '';
    const renderer = String(slot?.renderer || '').trim();
    const controlLabel = SUPPRESS_CONTROL_LABEL_BY_RENDERER[renderer] || 'Layout';
    return `Hidden by ${controlLabel}.`;
  };

  const updateCircularTrackSlotPlacement = (slot, placement) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!slot || !SUPPORTED_RENDERERS.includes(slot.renderer)) return;
    if (normalizeOptionalText(placement) === null) {
      slot.side = null;
      if (slot.renderer === 'features') {
        slot.params = cloneParams(slot.params);
        delete slot.params.lane_direction;
      }
      return;
    }
    applyPlacementDefaults(slot, placement);
    if (slot.renderer === 'annotations') {
      slot.params = cloneParams(slot.params);
      if (normalizePlacement(placement) === 'overlay') {
        bindOnlyCircularAnnotationAnchor(slot);
      } else {
        delete slot.params.anchor_slot;
        delete slot.params.cover_anchor;
      }
    }
    normalizeSlotsInPlace();
  };

  const updateCircularTrackFeatureLane = (slot, laneDirection) => {
    const drawing = state.drawings.circular;
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
    const index = drawing.adv.circular_track_slots.findIndex((candidate) => candidate === slot);
    if (index >= 0) {
      moveCircularTrackSlotToPlacement(drawing, index, sideForLaneDirection(normalizedLaneDirection));
      return;
    }
    slot.params.lane_direction = normalizedLaneDirection;
    slot.side = sideForLaneDirection(slot.params.lane_direction);
    normalizeSlotsInPlace();
  };

  const circularTrackRendererLabel = (renderer) => RENDERER_LABELS[renderer] || renderer;
  const supportsCircularTrackSlotPlacement = (renderer) => SUPPORTED_RENDERERS.includes(renderer);
  const isManagedCircularConservationSlot = (slot) => isManagedConservationSlot(slot);
  const circularTrackSlotDisplayLabel = (slot) => {
    if (isManagedConservationSlot(slot)) {
      return String(slot?.params?.label || slot?.params?.fileName || slot?.id || 'Conservation').trim();
    }
    return circularTrackRendererLabel(slot?.renderer);
  };
  const circularTrackSlotDisplayMeta = (slot) => {
    if (isManagedConservationSlot(slot)) {
      return String(slot?.params?.fileName || slot?.params?.series_key || '').trim();
    }
    return '';
  };
  const circularTrackSlotLegendLabelPlaceholder = (slot) => {
    const drawing = state.drawings.circular;
    const renderer = String(slot?.renderer || '').trim();
    if (renderer === 'dinucleotide_content' || renderer === 'dinucleotide_skew') {
      const nt = normalizeNt(slot?.params?.nt ?? slot?.params?.dinucleotide, normalizeNt(drawing.adv.nt));
      return renderer === 'dinucleotide_content' ? `${nt} content` : `${nt} skew`;
    }
    if (renderer === 'depth') return 'Depth';
    return 'Legend label';
  };
  const circularTrackSlotColor = (slot) => {
    if (isManagedConservationSlot(slot)) return String(slot?.params?.color || '').trim();
    return '';
  };

  const circularTrackSlotHasSkewColorOverride = (slot, key) => (
    slot?.renderer === 'dinucleotide_skew' &&
    normalizeOptionalText(slot?.params?.[key]) !== null
  );

  const circularTrackSlotSkewColorValue = (slot, key) => {
    const drawing = state.drawings.circular;
    return resolveTrackSlotSkewColorValue(/** @type {any} */ ({
      slot,
      key,
      currentColors: drawing.currentColors,
      paletteDefinitions: state.paletteDefinitions,
      selectedPalette: drawing.selectedPalette
    }));
  };

  const setCircularTrackSlotSkewColor = (slot, key, value) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!slot || slot.renderer !== 'dinucleotide_skew' || !['positive_color', 'negative_color'].includes(key)) return;
    slot.params = cloneParams(slot.params);
    const color = normalizeColorParam(value);
    if (color === null) delete slot.params[key];
    else slot.params[key] = color;
  };

  const clearCircularTrackSlotSkewColor = (slot, key) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!slot || slot.renderer !== 'dinucleotide_skew' || !['positive_color', 'negative_color'].includes(key)) return;
    slot.params = cloneParams(slot.params);
    delete slot.params[key];
  };

  const circularTrackSlots = () => {
    const drawing = state.drawings.circular;
    return (
      Array.isArray(drawing.adv.circular_track_slots)
        ? drawing.adv.circular_track_slots.map((slot, index) => ({ kind: STACK_ENTRY_SLOT, slot, index }))
        : []
    );
  };

  const circularTrackStackEntries = () => {
    const drawing = state.drawings.circular;
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

  const circularTrackSlotCliSpec = (slot) => {
    const drawing = state.drawings.circular;
    normalizeSlotsInPlace();
    return buildCircularTrackSlotSpec(slot, drawing.adv.nt, drawing.form.track_type);
  };

  const updateCircularTrackSlotMeasure = (slot, field, scalar) => {
    const drawing = state.drawings.circular;
    if (!['width', 'radius'].includes(field) || !drawing.adv.circular_track_slots.includes(slot)) return;
    const numericLeaf = scalar && typeof scalar === 'object' && !Array.isArray(scalar)
      ? scalar.value : scalar;
    if (typeof numericLeaf === 'number' && !Number.isFinite(numericLeaf)) {
      parseOptionalCircularScalar(scalar, `Circular track slot '${slot.id}' ${field}`);
    }
    if (slot[field] !== scalar) slot[field] = scalar;
  };

  const circularSlotManualValue = (slot, field) => {
    if (!slot) return '';
    if (field === 'width') return slot.width;
    if (field === 'radius') return slot.radius;
    if (field === 'inner_gap_px') return slot.inner_gap_px;
    if (field === 'outer_gap_px') return slot.outer_gap_px;
    return '';
  };

  const selectedResultIndexValue = () => Number(state?.selectedResultIndex?.value ?? 0) || 0;

  const circularSlotGeometryLookup = (slotId) => (/** @type {any} */ ({
    geometry: String(state?.trackSlotResolvedGeometry?.value?.mode || '') === 'circular'
      ? state.trackSlotResolvedGeometry.value
      : null,
    resultIndex: selectedResultIndexValue(),
    recordIndex: 0,
    slotId
  }));

  /** @param {DrawingState} drawing */
  const estimateCircularSlotGeometry = (drawing, slot, slotIndex) => {
    const preset = normalizeCircularTrackPreset(drawing.form.track_type);
    const lengthParam = getPreviewLengthParam(state);
    const renderer = String(slot?.renderer || '').trim();
    const widthPx = previewWidthPxForRenderer(renderer, lengthParam);
    const spacingPx = previewSpacingPx();
    let radiusFactor = getPresetRadiusRatio(slot, renderer, preset, lengthParam, drawing);
    if (radiusFactor === null) {
      const slots = Array.isArray(drawing.adv.circular_track_slots) ? drawing.adv.circular_track_slots : [];
      const axis = clampCircularTrackAxisIndex(drawing.adv.circular_track_slots_axis_index, slots.length)
        ?? inferLegacyAxisIndexFromFeature(slots, drawing.form.track_type);
      const placement = effectiveSlotPlacement(slot, preset);
      const distance = Math.max(1, Math.abs(Number(slotIndex) - Number(axis)) + 1);
      const step = (widthPx + spacingPx) / Math.max(1, PREVIEW_RADIUS_PX);
      if (placement === 'outside') radiusFactor = 1 + (distance * step);
      else if (placement === 'overlay') radiusFactor = 1.0;
      else radiusFactor = Math.max(0.05, 1 - (distance * step));
    }
    return {
      widthPx,
      radiusFactor,
      innerGapPx: spacingPx,
      outerGapPx: spacingPx
    };
  };

  // Only a rendered row has resolved geometry. Any other row, and every row
  // before the first Generate, shows an estimate marked as one (TK-15).
  const circularTrackSlotGeometryAutoText = (slot, slotIndex, field) => {
    if (isManualSlotValue(circularSlotManualValue(slot, field))) return '';
    const lookup = circularSlotGeometryLookup(slot?.id);
    const resolved = circularTrackSlotEffectiveEnabled(slot) ? findTrackSlotGeometry(lookup) : null;
    const estimate = !resolved;
    const geometry = resolved || estimateCircularSlotGeometry(state.drawings.circular, slot, slotIndex);
    if (field === 'width') return formatPxAuto(geometry.widthPx, { estimate });
    if (field === 'radius') {
      // `r` pins a ticks row's anchor, not the band centre the payload holds (GX-18).
      const tickAnchor = resolved && slot?.renderer === 'ticks'
        ? tickAnchorRadiusFactor(resolved, findTrackSlotGeometryRecord(lookup)?.axisRadiusPx, slot.params?.tick_label_layout)
        : null;
      return tickAnchor === null
        ? formatRadiusFactorAuto(geometry.radiusFactor, { estimate })
        : formatRadiusFactorAuto(tickAnchor, { awayFromAxis: true });
    }
    if (field === 'inner_gap_px') return formatPxAuto(geometry.innerGapPx, { estimate });
    if (field === 'outer_gap_px') return formatPxAuto(geometry.outerGapPx, { estimate });
    return '';
  };

  const circularTrackSlotGeometryUnitSuffix = (slot, field) => {
    const text = String(circularSlotManualValue(slot, field) ?? '').trim();
    if (!text) return '';
    if (!['inner_gap_px', 'outer_gap_px'].includes(field)) return '';
    return /px$/i.test(text) ? '' : 'px';
  };

  const circularTrackSlotGeometryHasManual = (slot, field) => (
    isManualSlotValue(circularSlotManualValue(slot, field))
  );

  const circularTrackPresetSummary = () => {
    const drawing = state.drawings.circular;
    const preset = normalizeCircularTrackPreset(drawing.form.track_type);
    const lengthParam = getPreviewLengthParam(state);
    const lane = laneDirectionForPreset(preset);
    const pieces = [
      laneDirectionLabel(lane),
      `feature r ${previewFeatureRadiusRatio(preset, lengthParam, drawing).toFixed(2)}x`
    ];
    const slotLike = (id, renderer) => ({ id, renderer });
    if (Boolean(drawing.form.show_depth)) {
      const depthRatio = getPresetRadiusRatio(slotLike('depth', 'depth'), 'depth', preset, lengthParam, drawing);
      if (depthRatio !== null) pieces.push(`depth r ${depthRatio.toFixed(2)}x`);
    }
    if (!Boolean(drawing.form.suppress_gc)) {
      const gcRatio = getPresetRadiusRatio(slotLike('gc_content', 'dinucleotide_content'), 'dinucleotide_content', preset, lengthParam, drawing);
      if (gcRatio !== null) pieces.push(`GC r ${gcRatio.toFixed(2)}x`);
    }
    if (!Boolean(drawing.form.suppress_skew)) {
      const skewRatio = getPresetRadiusRatio(slotLike('gc_skew', 'dinucleotide_skew'), 'dinucleotide_skew', preset, lengthParam, drawing);
      if (skewRatio !== null) pieces.push(`skew r ${skewRatio.toFixed(2)}x`);
    }
    return {
      label: `${formatPresetName(preset)} preset`,
      detail: `${lengthParam} defaults: ${pieces.join(' · ')}`,
    };
  };

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
    circularTrackRendererLabel,
    normalizeCircularTrackSlots: normalizeSlotsInPlace,
    syncCircularConservationSlots,
    changeCircularDepthSources,
    ...featureSlotEdits(changeTrackLayout, {
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
    canAddCircularTrackRenderer,
    canDuplicateCircularTrackSlot,
    circularTrackSlotEffectiveEnabled,
    circularTrackSlotHiddenBySuppress,
    circularTrackSlotSuppressMessage,
    canMoveCircularTrackSlot,
    canMoveCircularTrackSlotOutside,
    canMoveCircularTrackSlotInside,
    canMoveCircularTrackSlotToAxis,
    updateCircularTrackSlotMeasure,
    circularTrackSlotIssue,
    circularTrackGlobalIssues,
    circularAnnotationAnchorOptions,
    circularAnnotationAnchorIsKnown,
    annotationTrackMarkOptions: ANNOTATION_MARK_OPTIONS,
    circularAnnotationMarkSelected,
    setCircularAnnotationMarkSelected,
    circularAnnotationLaneGapValue: (slot) => circularAnnotationNumberValue(slot, 'lane_gap_px', 3),
    setCircularAnnotationLaneGap: (slot, value) => setCircularAnnotationNumber(slot, 'lane_gap_px', value, 3),
    circularAnnotationPaddingValue: (slot) => circularAnnotationNumberValue(slot, 'padding_px', 2),
    setCircularAnnotationPadding: (slot, value) => setCircularAnnotationNumber(slot, 'padding_px', value, 2),
    circularAnnotationCoverAnchor,
    setCircularAnnotationCoverAnchor,
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
