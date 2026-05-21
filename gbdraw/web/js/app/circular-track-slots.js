const SUPPORTED_RENDERERS = [
  'features',
  'ticks',
  'dinucleotide_content',
  'dinucleotide_skew',
  'depth',
  'spacer'
];

const RENDERER_LABELS = {
  features: 'Features',
  ticks: 'Ticks',
  dinucleotide_content: 'Dinucleotide content',
  dinucleotide_skew: 'Dinucleotide skew',
  depth: 'Depth',
  spacer: 'Spacer'
};

const DEFAULT_SLOT_IDS = {
  features: 'features',
  ticks: 'ticks',
  dinucleotide_content: 'gc_content',
  dinucleotide_skew: 'gc_skew',
  depth: 'depth',
  spacer: 'spacer'
};

const NUMERIC_RENDERERS = new Set(['dinucleotide_content', 'dinucleotide_skew', 'depth']);
const STACK_ENTRY_AXIS = 'axis';
const STACK_ENTRY_SLOT = 'slot';
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
const TICK_LABEL_LAYOUTS = [
  'label_out_tick_in',
  'label_in_tick_out',
  'tick_only',
  'label_only'
];
const DEFAULT_TICK_LABEL_LAYOUT = 'label_out_tick_in';

export const CIRCULAR_TRACK_PRESETS = ['tuckin', 'middle', 'spreadout'];

export const normalizeCircularTrackPreset = (value, fallback = 'tuckin') => {
  const text = String(value || fallback).trim().toLowerCase();
  return CIRCULAR_TRACK_PRESETS.includes(text) ? text : fallback;
};

const laneDirectionForPreset = (preset) => {
  const normalized = normalizeCircularTrackPreset(preset);
  if (normalized === 'middle') return 'split';
  if (normalized === 'spreadout') return 'outside';
  return 'inside';
};

const normalizeLaneDirection = (value, fallback = 'inside') => {
  const text = String(value || fallback).trim().toLowerCase();
  return ['inside', 'outside', 'split'].includes(text) ? text : fallback;
};

const sideForLaneDirection = (value) => {
  const lane = normalizeLaneDirection(value);
  return lane === 'split' ? 'overlay' : lane;
};

const laneDirectionForSide = (value) => {
  const side = normalizePlacement(value);
  if (side === 'outside') return 'outside';
  if (side === 'overlay') return 'split';
  return 'inside';
};

const normalizeNt = (value, fallback = 'GC') => {
  const text = String(value || '').trim().toUpperCase();
  return text || fallback;
};

const cleanToken = (value, fallback) => {
  const text = String(value || '').trim();
  return text || fallback;
};

const normalizeOptionalText = (value) => {
  const text = String(value ?? '').trim();
  return text.length > 0 ? text : null;
};

const normalizePlacement = (value, fallback = 'inside') => {
  const text = String(value || fallback).trim().toLowerCase();
  return ['inside', 'outside', 'overlay'].includes(text) ? text : fallback;
};

const normalizeOptionalPlacement = (value) => {
  if (value === null || value === undefined || value === '') return null;
  return normalizePlacement(value);
};

const normalizeTickLabelLayout = (value, fallback = DEFAULT_TICK_LABEL_LAYOUT) => {
  const text = String(value || fallback).trim().toLowerCase();
  return TICK_LABEL_LAYOUTS.includes(text) ? text : fallback;
};

const tickLabelLayoutFromSides = (labelSide, tickSide, fallback = DEFAULT_TICK_LABEL_LAYOUT) => {
  const label = String(labelSide || '').trim().toLowerCase();
  const tick = String(tickSide || '').trim().toLowerCase();
  if (label === 'outside' && tick === 'inside') return 'label_out_tick_in';
  if (label === 'inside' && tick === 'outside') return 'label_in_tick_out';
  if (label === 'none' && ['inside', 'outside', 'both'].includes(tick)) return 'tick_only';
  if (['inside', 'outside'].includes(label) && tick === 'none') return 'label_only';
  return fallback;
};

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
  if (renderer === 'depth') return base * Number(factors[1]) * 0.5;
  if (renderer === 'dinucleotide_skew') return base * Number(factors[2]);
  if (renderer === 'ticks') return 0;
  return base * Number(factors[1]);
};

const previewSpacingPx = () => Math.max(1.0, 0.01 * PREVIEW_RADIUS_PX);

const previewFeatureLaneCount = (state) => {
  if (Boolean(state?.form?.separate_strands)) return 2;
  return 1;
};

const previewFeatureRadiusRatio = (preset, lengthParam, state) => {
  const normalized = normalizeCircularTrackPreset(preset);
  const laneWidth = previewWidthPxForRenderer('features', lengthParam);
  const laneCount = previewFeatureLaneCount(state);
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

const getBuiltinTrackId = (slot, renderer, state) => {
  const id = String(slot?.id || '').trim();
  const showDepth = Boolean(state?.form?.show_depth);
  const showGc = !Boolean(state?.form?.suppress_gc);
  const showSkew = !Boolean(state?.form?.suppress_skew);

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

const getPresetRadiusRatio = (slot, renderer, preset, lengthParam, state) => {
  if (renderer === 'features' && String(slot?.id || '').trim() === 'features') {
    return previewFeatureRadiusRatio(preset, lengthParam, state);
  }
  const trackId = getBuiltinTrackId(slot, renderer, state);
  if (trackId === null) return null;
  return PREVIEW_TRACK_DICT[lengthParam]?.[normalizeCircularTrackPreset(preset)]?.[trackId] ?? null;
};

const slotHasManualGeometry = (slot) => (
  normalizeOptionalText(slot?.width) !== null ||
  normalizeOptionalText(slot?.radius) !== null ||
  normalizeOptionalText(slot?.spacing) !== null
);

const cloneParams = (params = {}) => {
  if (!params || typeof params !== 'object' || Array.isArray(params)) return {};
  return { ...params };
};

const normalizeSlotSide = (value) => normalizeOptionalPlacement(value);

const applyPlacementDefaults = (slot, placement = 'inside') => {
  if (!slot) return;
  const requestedSide = normalizePlacement(placement);
  slot.side = requestedSide;
  slot.params = cloneParams(slot.params);
  if (slot.renderer === 'features') {
    slot.params.lane_direction = laneDirectionForSide(requestedSide);
  }
};

const featureLaneForSlot = (slot, preset = 'tuckin') => {
  const params = cloneParams(slot?.params);
  const rawLane = normalizeOptionalText(params.lane_direction ?? params.lanes);
  return normalizeLaneDirection(rawLane, laneDirectionForPreset(preset));
};

const syncFeaturePlacement = (slot, preset = 'tuckin') => {
  if (!slot || slot.renderer !== 'features') return;
  const lane = featureLaneForSlot(slot, preset);
  slot.side = sideForLaneDirection(lane);
  slot.params = cloneParams(slot.params);
  slot.params.lane_direction = lane;
};

const syncTickParamsForPlacement = (slot, side) => {
  if (!slot || slot.renderer !== 'ticks') return;
  void side;
  slot.params = cloneParams(slot.params);
  slot.params.tick_label_layout = normalizeTickLabelLayout(slot.params.tick_label_layout);
  delete slot.params.label_side;
  delete slot.params.tick_side;
};

const syncSlotPlacementFromSide = (slot, side) => {
  if (!slot) return;
  const placement = normalizePlacement(side);
  slot.side = placement;
  if (slot.renderer === 'features') {
    slot.params = cloneParams(slot.params);
    slot.params.lane_direction = laneDirectionForSide(placement);
  } else if (slot.renderer === 'ticks') {
    syncTickParamsForPlacement(slot, placement);
  }
};

const effectiveSlotPlacement = (slot, preset = 'tuckin') => {
  if (!slot) return 'inside';
  if (slot.renderer === 'features') {
    return sideForLaneDirection(featureLaneForSlot(slot, preset));
  }
  return normalizePlacement(slot.side, 'inside');
};

export const inferLegacyAxisIndexFromFeature = (slots, preset = 'tuckin') => {
  const onAxisIndex = slots.findIndex((slot) => effectiveSlotPlacement(slot, preset) === 'overlay');
  if (onAxisIndex >= 0) return onAxisIndex;
  const featureIndex = slots.findIndex((slot) => slot?.enabled !== false && slot?.renderer === 'features');
  if (featureIndex < 0) {
    const firstInside = slots.findIndex((slot) => normalizePlacement(slot?.side, 'inside') !== 'outside');
    return firstInside < 0 ? slots.length : firstInside;
  }
  const featureSide = sideForLaneDirection(featureLaneForSlot(slots[featureIndex], preset));
  if (featureSide === 'outside') return featureIndex + 1;
  return featureIndex;
};

const axisIndexForSlots = inferLegacyAxisIndexFromFeature;

export const clampCircularTrackAxisIndex = (value, slotCount) => {
  const length = Math.max(0, Number(slotCount) || 0);
  const numeric = Number(value);
  if (!Number.isInteger(numeric)) return null;
  return Math.min(Math.max(numeric, 0), length);
};

const syncSlotsFromAxisIndex = (slots, axisIndex, preset = 'tuckin') => {
  const axis = clampCircularTrackAxisIndex(axisIndex, slots.length);
  const resolvedAxis = axis === null ? inferLegacyAxisIndexFromFeature(slots, preset) : axis;
  slots.forEach((slot, index) => {
    if (!slot) return;
    if (effectiveSlotPlacement(slot, preset) === 'overlay') {
      syncSlotPlacementFromSide(slot, 'overlay');
      return;
    }
    if (slot.renderer === 'features' && featureLaneForSlot(slot, preset) === 'split') {
      syncSlotPlacementFromSide(slot, 'overlay');
      return;
    }
    const side = index < resolvedAxis ? 'outside' : 'inside';
    syncSlotPlacementFromSide(slot, side);
  });
  return resolvedAxis;
};

const enforceSingleOnAxisSlot = (slots, axisIndex, preset = 'tuckin') => {
  const onAxisIndices = slots
    .map((slot, index) => (
      effectiveSlotPlacement(slot, preset) === 'overlay' ? index : null
    ))
    .filter((index) => Number.isInteger(index));
  if (onAxisIndices.length === 0) {
    return clampCircularTrackAxisIndex(axisIndex, slots.length) ?? inferLegacyAxisIndexFromFeature(slots, preset);
  }

  const clampedAxis = clampCircularTrackAxisIndex(axisIndex, slots.length);
  const keepIndex = onAxisIndices.includes(clampedAxis) ? clampedAxis : onAxisIndices[0];
  onAxisIndices.forEach((index) => {
    if (index === keepIndex) return;
    syncSlotPlacementFromSide(slots[index], index < keepIndex ? 'outside' : 'inside');
  });
  return keepIndex;
};

export const applyCircularTrackOrderPlacements = (slots, defaultNt = 'GC', preset = 'tuckin', axisIndex = null) => {
  const normalized = normalizeCircularTrackSlots(slots, defaultNt, preset);
  const resolvedAxis = syncSlotsFromAxisIndex(normalized, axisIndex, preset);
  enforceSingleOnAxisSlot(normalized, resolvedAxis, preset);
  return normalized;
};

const makeSlot = ({
  id,
  renderer,
  enabled = true,
  width = null,
  radius = null,
  spacing = null,
  side = null,
  z = 0,
  params = {}
}) => ({
  id,
  renderer,
  enabled,
  width,
  radius,
  spacing,
  side: normalizeSlotSide(side),
  z,
  params: cloneParams(params)
});

const paramsMatchExactly = (params, expected = {}) => {
  const actualEntries = Object.entries(cloneParams(params))
    .filter(([, value]) => normalizeOptionalText(value) !== null)
    .map(([key, value]) => [String(key), normalizeOptionalText(value)]);
  const expectedEntries = Object.entries(expected)
    .filter(([, value]) => normalizeOptionalText(value) !== null)
    .map(([key, value]) => [String(key), normalizeOptionalText(value)]);
  if (actualEntries.length !== expectedEntries.length) return false;
  const actualMap = Object.fromEntries(actualEntries);
  return expectedEntries.every(([key, value]) => actualMap[key] === value);
};

const hasBlankSlotGeometry = (source) =>
  normalizeOptionalText(source.width) === null &&
  normalizeOptionalText(source.radius) === null &&
  normalizeOptionalText(source.spacing) === null;

const hasDefaultSlotFlags = (source) =>
  source.enabled !== false &&
  Number(source.z || 0) === 0;

const isLegacyDefaultWebSlotShape = (source, renderer, defaultNt = 'GC', preset = 'tuckin') => {
  if (!source || typeof source !== 'object' || Array.isArray(source)) return false;
  if (!hasBlankSlotGeometry(source) || !hasDefaultSlotFlags(source)) return false;

  const normalizedPreset = normalizeCircularTrackPreset(preset);
  const normalizedId = String(source.id || '').trim();
  const side = normalizeOptionalPlacement(source.side);
  const params = cloneParams(source.params);

  if (renderer === 'features' && normalizedId === 'features') {
    const laneDirection = laneDirectionForPreset(normalizedPreset);
    return (
      side === sideForLaneDirection(laneDirection) &&
      paramsMatchExactly(params, { lane_direction: laneDirection })
    );
  }

  if (renderer === 'ticks' && normalizedId === 'ticks') {
    return (
      side === 'inside' &&
      paramsMatchExactly(params, {
        tick_label_layout: DEFAULT_TICK_LABEL_LAYOUT,
        preset: normalizedPreset
      })
    );
  }

  if (renderer === 'depth' && normalizedId === 'depth') {
    return side === 'inside' && paramsMatchExactly(params, {});
  }

  if (renderer === 'dinucleotide_content' && normalizedId === 'gc_content') {
    return (
      side === 'inside' &&
      paramsMatchExactly(params, { nt: normalizeNt(defaultNt) })
    );
  }

  if (renderer === 'dinucleotide_skew' && normalizedId === 'gc_skew') {
    return (
      side === 'inside' &&
      paramsMatchExactly(params, { nt: normalizeNt(defaultNt) })
    );
  }

  return false;
};

export const createDefaultCircularTrackSlots = ({
  nt = 'GC',
  showDepth = false,
  showGc = true,
  showSkew = true,
  preset = 'tuckin'
} = {}) => {
  void nt;
  void preset;
  const slots = [
    makeSlot({
      id: 'features',
      renderer: 'features'
    }),
    makeSlot({
      id: 'ticks',
      renderer: 'ticks'
    })
  ];
  if (showDepth) slots.push(makeSlot({ id: 'depth', renderer: 'depth' }));
  if (showGc) {
    slots.push(makeSlot({
      id: 'gc_content',
      renderer: 'dinucleotide_content'
    }));
  }
  if (showSkew) {
    slots.push(makeSlot({
      id: 'gc_skew',
      renderer: 'dinucleotide_skew'
    }));
  }
  return slots;
};

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
  }
  void nt;

  const slot = makeSlot({
    id,
    renderer: normalizedRenderer,
    side,
    params
  });
  if (side !== null) applyPlacementDefaults(slot, side);
  return slot;
};

export const normalizeCircularTrackSlot = (slot, index = 0, defaultNt = 'GC', preset = 'tuckin') => {
  const source = slot && typeof slot === 'object' && !Array.isArray(slot) ? slot : {};
  const renderer = SUPPORTED_RENDERERS.includes(source.renderer) ? source.renderer : 'dinucleotide_skew';
  const fallbackId = DEFAULT_SLOT_IDS[renderer] || `slot_${index + 1}`;
  const inheritsPresetDefaults = isLegacyDefaultWebSlotShape(source, renderer, defaultNt, preset);
  const params = inheritsPresetDefaults ? {} : cloneParams(source.params);
  [
    'side',
    'strict',
    'compress',
    'reserve',
    'r',
    'radius',
    'w',
    'width',
    'spacing'
  ].forEach((key) => {
    delete params[key];
  });

  let side = inheritsPresetDefaults ? null : normalizeSlotSide(source.side);
  const radius = source.radius ?? null;
  const spacing = source.spacing ?? null;

  if (renderer === 'dinucleotide_content' || renderer === 'dinucleotide_skew') {
    const nt = normalizeOptionalText(params.nt ?? params.dinucleotide);
    if (nt === null) {
      delete params.nt;
      delete params.dinucleotide;
    } else {
      params.nt = normalizeNt(nt);
      delete params.dinucleotide;
    }
  }
  if (renderer === 'ticks') {
    delete params['axis'];
    const legacyLayout = (
      normalizeOptionalText(params.tick_label_layout) === null
        ? tickLabelLayoutFromSides(params.label_side, params.tick_side)
        : params.tick_label_layout
    );
    params.tick_label_layout = normalizeTickLabelLayout(legacyLayout);
    delete params.label_side;
    delete params.tick_side;
    if (normalizeOptionalText(params.preset) === null) {
      delete params.preset;
    } else {
      params.preset = normalizeCircularTrackPreset(params.preset);
    }
  }
  if (renderer === 'features') {
    const rawLaneDirection = normalizeOptionalText(params.lane_direction ?? params.lanes);
    if (rawLaneDirection === null) {
      delete params.lane_direction;
      delete params.lanes;
      if (side !== null) params.lane_direction = laneDirectionForSide(side);
    } else {
      const normalizedLaneDirection = normalizeLaneDirection(rawLaneDirection);
      delete params.lanes;
      if (normalizedLaneDirection === 'inside') {
        params.lane_direction = 'inside';
        side = 'inside';
      } else {
        params.lane_direction = normalizedLaneDirection;
        side = sideForLaneDirection(params.lane_direction);
      }
    }
  }

  return makeSlot({
    id: cleanToken(source.id, fallbackId),
    renderer,
    enabled: source.enabled !== false,
    width: source.width ?? null,
    radius,
    spacing,
    side,
    z: Number.isFinite(Number(source.z)) ? Number(source.z) : 0,
    params
  });
};

export const normalizeCircularTrackSlots = (slots, defaultNt = 'GC', preset = 'tuckin') => {
  const base = Array.isArray(slots)
    ? slots
    : createDefaultCircularTrackSlots({ nt: defaultNt, preset });
  return base.map((slot, index) => normalizeCircularTrackSlot(slot, index, defaultNt, preset));
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
  const forceSplitLane = optionsOverride?.forceSplitLane === true;

  if (!normalized.enabled) options.push('enabled=false');
  appendOption(options, 'w', normalized.width);
  appendOption(options, 'r', normalized.radius);
  appendOption(options, 'spacing', normalized.spacing);
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
    if (
      normalizeOptionalText(params.lane_direction) !== null &&
      (
        (forceSplitLane && normalizeLaneDirection(params.lane_direction) === 'split') ||
        normalizeLaneDirection(params.lane_direction) !== laneDirectionForPreset(normalizedPreset)
      )
    ) {
      appendOption(options, 'lane_direction', params.lane_direction);
    }
  } else if (normalized.renderer === 'dinucleotide_content' || normalized.renderer === 'dinucleotide_skew') {
    const nt = normalizeOptionalText(params.nt);
    if (nt !== null && normalizeNt(nt) !== normalizeNt(defaultNt)) {
      options.push(`nt=${normalizeNt(nt)}`);
    }
  }
  appendOption(options, 'legend_label', params.legend_label);

  return `${normalized.id}:${normalized.renderer}${options.length ? `@${options.join(',')}` : ''}`;
};

export const hasEnabledCircularTrackRenderer = (slots, renderer) =>
  normalizeCircularTrackSlots(slots).some((slot) => slot.enabled && slot.renderer === renderer);

export const createCircularTrackSlotEditor = ({ state }) => {
  const axisIndexForCurrentSlots = (slots) => {
    const current = clampCircularTrackAxisIndex(state.adv.circular_track_slots_axis_index, slots.length);
    if (current !== null) {
      state.adv.circular_track_slots_axis_index = current;
      return current;
    }
    const inferred = inferLegacyAxisIndexFromFeature(slots, state.form.track_type);
    state.adv.circular_track_slots_axis_index = inferred;
    return inferred;
  };

  const normalizedSlotsForCurrentState = () => applyCircularTrackOrderPlacements(
    state.adv.circular_track_slots,
    state.adv.nt,
    state.form.track_type,
    state.adv.circular_track_slots_axis_index
  );

  const normalizeSlotsInPlace = () => {
    const normalized = normalizeCircularTrackSlots(
      state.adv.circular_track_slots,
      state.adv.nt,
      state.form.track_type
    );
    const axis = axisIndexForCurrentSlots(normalized);
    syncSlotsFromAxisIndex(normalized, axis, state.form.track_type);
    state.adv.circular_track_slots_axis_index = enforceSingleOnAxisSlot(
      normalized,
      axis,
      state.form.track_type
    );
    state.adv.circular_track_slots.splice(
      0,
      state.adv.circular_track_slots.length,
      ...normalized
    );
  };

  const resetCircularTrackSlotsFromSimpleControls = () => {
    const slots = createDefaultCircularTrackSlots({
      nt: state.adv.nt,
      showDepth: Boolean(state.form.show_depth),
      showGc: !state.form.suppress_gc,
      showSkew: !state.form.suppress_skew,
      preset: state.form.track_type
    });
    state.adv.circular_track_slots_axis_index = inferLegacyAxisIndexFromFeature(
      normalizeCircularTrackSlots(slots, state.adv.nt, state.form.track_type),
      state.form.track_type
    );
    const normalized = applyCircularTrackOrderPlacements(
      slots,
      state.adv.nt,
      state.form.track_type,
      state.adv.circular_track_slots_axis_index
    );
    state.adv.circular_track_slots.splice(0, state.adv.circular_track_slots.length, ...normalized);
  };

  const ensureCircularTrackDepthSlot = () => {
    normalizeSlotsInPlace();
    const existingDepthSlot = state.adv.circular_track_slots.find((slot) => slot?.renderer === 'depth');
    if (existingDepthSlot) {
      existingDepthSlot.enabled = true;
      return;
    }

    const currentSlots = JSON.stringify(normalizeCircularTrackSlots(
      state.adv.circular_track_slots,
      state.adv.nt,
      state.form.track_type
    ));
    const simpleSlotsWithoutDepth = JSON.stringify(normalizeCircularTrackSlots(
      createDefaultCircularTrackSlots({
        nt: state.adv.nt,
        showDepth: false,
        showGc: !state.form.suppress_gc,
        showSkew: !state.form.suppress_skew,
        preset: state.form.track_type
      }),
      state.adv.nt,
      state.form.track_type
    ));

    if (currentSlots === simpleSlotsWithoutDepth) {
      resetCircularTrackSlotsFromSimpleControls();
      return;
    }

    const slot = createCircularTrackSlotForRenderer('depth', state.adv.circular_track_slots, state.adv.nt);
    state.adv.circular_track_slots.push(slot);
    normalizeSlotsInPlace();
  };

  const resetCircularTrackSlotsToPreset = (preset) => {
    const normalizedPreset = normalizeCircularTrackPreset(preset);
    const templateSlots = createDefaultCircularTrackSlots({
      nt: state.adv.nt,
      showDepth: Boolean(state.form.show_depth),
      showGc: !state.form.suppress_gc,
      showSkew: !state.form.suppress_skew,
      preset: normalizedPreset
    });
    state.adv.circular_track_slots_axis_index = inferLegacyAxisIndexFromFeature(
      normalizeCircularTrackSlots(templateSlots, state.adv.nt, normalizedPreset),
      normalizedPreset
    );
    const normalized = applyCircularTrackOrderPlacements(
      templateSlots,
      state.adv.nt,
      normalizedPreset,
      state.adv.circular_track_slots_axis_index
    );
    state.form.track_type = normalizedPreset;
    state.adv.circular_track_slots.splice(0, state.adv.circular_track_slots.length, ...normalized);
  };

  const applyCircularTrackPreset = (preset) => resetCircularTrackSlotsToPreset(preset);

  const setCircularTrackSlotsEnabled = (enabled) => {
    state.adv.circular_track_slots_enabled = Boolean(enabled);
    if (state.adv.circular_track_slots_enabled) {
      resetCircularTrackSlotsFromSimpleControls();
    }
  };

  const addCircularTrackSlot = (renderer, placement = null) => {
    normalizeSlotsInPlace();
    const slot = createCircularTrackSlotForRenderer(renderer, state.adv.circular_track_slots, state.adv.nt, placement);
    state.adv.circular_track_slots.push(slot);
    normalizeSlotsInPlace();
  };

  const duplicateCircularTrackSlot = (index) => {
    normalizeSlotsInPlace();
    const idx = Number(index);
    if (!Number.isInteger(idx) || idx < 0 || idx >= state.adv.circular_track_slots.length) return;
    const source = normalizeCircularTrackSlot(state.adv.circular_track_slots[idx], idx, state.adv.nt, state.form.track_type);
    const duplicate = createCircularTrackSlotForRenderer(source.renderer, state.adv.circular_track_slots, state.adv.nt);
    duplicate.enabled = source.enabled;
    duplicate.width = source.width;
    duplicate.radius = source.radius;
    duplicate.spacing = source.spacing;
    duplicate.side = source.side;
    duplicate.z = source.z;
    duplicate.params = cloneParams(source.params);
    state.adv.circular_track_slots.splice(idx + 1, 0, duplicate);
    const axis = axisIndexForCurrentSlots(state.adv.circular_track_slots);
    if (idx < axis) state.adv.circular_track_slots_axis_index = axis + 1;
    normalizeSlotsInPlace();
  };

  const removeCircularTrackSlot = (index) => {
    const idx = Number(index);
    if (!Number.isInteger(idx) || idx < 0 || idx >= state.adv.circular_track_slots.length) return;
    const axis = axisIndexForCurrentSlots(state.adv.circular_track_slots);
    state.adv.circular_track_slots.splice(idx, 1);
    state.adv.circular_track_slots_axis_index = idx < axis ? Math.max(0, axis - 1) : Math.min(axis, state.adv.circular_track_slots.length);
    normalizeSlotsInPlace();
  };

  const wouldCircularTrackSlotMoveCrossAxis = (fromIndex, toIndex) => {
    const from = Number(fromIndex);
    const to = Number(toIndex);
    const normalized = normalizedSlotsForCurrentState();
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

    const axis = axisIndexForCurrentSlots(normalized);
    const movedPlacement = effectiveSlotPlacement(normalized[from], state.form.track_type);
    if (movedPlacement === 'overlay') return true;
    return (from < axis) !== (to < axis);
  };

  const moveCircularTrackSlot = (fromIndex, toIndex) => {
    if (wouldCircularTrackSlotMoveCrossAxis(fromIndex, toIndex)) return;
    normalizeSlotsInPlace();
    const from = Number(fromIndex);
    const to = Number(toIndex);
    if (
      !Number.isInteger(from) ||
      !Number.isInteger(to) ||
      from < 0 ||
      to < 0 ||
      from >= state.adv.circular_track_slots.length ||
      to >= state.adv.circular_track_slots.length ||
      from === to
    ) {
      return;
    }
    const [moved] = state.adv.circular_track_slots.splice(from, 1);
    state.adv.circular_track_slots.splice(to, 0, moved);
    normalizeSlotsInPlace();
  };

  const canMoveCircularTrackSlot = (index, direction) => {
    const idx = Number(index);
    const step = Number(direction);
    if (!Number.isInteger(idx) || !Number.isInteger(step) || step === 0) return false;
    const target = idx + Math.sign(step);
    return !wouldCircularTrackSlotMoveCrossAxis(idx, target);
  };

  const canMoveCircularTrackSlotOutside = (index) => {
    const idx = Number(index);
    const normalized = normalizedSlotsForCurrentState();
    if (!Number.isInteger(idx) || idx < 0 || idx >= normalized.length) return false;
    return idx >= axisIndexForCurrentSlots(normalized);
  };

  const canMoveCircularTrackSlotInside = (index) => {
    const idx = Number(index);
    const normalized = normalizedSlotsForCurrentState();
    if (!Number.isInteger(idx) || idx < 0 || idx >= normalized.length) return false;
    return idx < axisIndexForCurrentSlots(normalized);
  };

  const canMoveCircularTrackSlotToAxis = (index) => {
    const idx = Number(index);
    const normalized = normalizedSlotsForCurrentState();
    if (!Number.isInteger(idx) || idx < 0 || idx >= normalized.length) return false;
    const slot = normalized[idx];
    return ['features', 'ticks'].includes(slot?.renderer) && effectiveSlotPlacement(slot, state.form.track_type) !== 'overlay';
  };

  const moveCircularTrackSlotToPlacement = (index, placement) => {
    const idx = Number(index);
    if (!Number.isInteger(idx) || idx < 0 || idx >= state.adv.circular_track_slots.length) return;
    normalizeSlotsInPlace();
    if (idx >= state.adv.circular_track_slots.length) return;
    const targetPlacement = normalizePlacement(placement);

    if (targetPlacement === 'overlay') {
      const movedSlot = state.adv.circular_track_slots[idx];
      if (!movedSlot) return;
      const movedPreviousPlacement = effectiveSlotPlacement(movedSlot, state.form.track_type);
      const existingAxisIndex = state.adv.circular_track_slots.findIndex((slot, slotIndex) => (
        slotIndex !== idx &&
        effectiveSlotPlacement(slot, state.form.track_type) === 'overlay'
      ));
      if (existingAxisIndex >= 0) {
        const existingAxisSlot = state.adv.circular_track_slots[existingAxisIndex];
        const demotedPlacement = movedPreviousPlacement === 'overlay' ? 'inside' : movedPreviousPlacement;
        syncSlotPlacementFromSide(existingAxisSlot, demotedPlacement);
        syncSlotPlacementFromSide(movedSlot, 'overlay');
        state.adv.circular_track_slots[existingAxisIndex] = movedSlot;
        state.adv.circular_track_slots[idx] = existingAxisSlot;
        state.adv.circular_track_slots_axis_index = existingAxisIndex;
        normalizeSlotsInPlace();
        return;
      }
    }

    let axis = axisIndexForCurrentSlots(state.adv.circular_track_slots);
    const [slot] = state.adv.circular_track_slots.splice(idx, 1);
    if (!slot) return;
    if (idx < axis) axis -= 1;
    syncSlotPlacementFromSide(slot, targetPlacement);
    if (targetPlacement === 'outside') {
      state.adv.circular_track_slots.splice(axis, 0, slot);
      state.adv.circular_track_slots_axis_index = axis + 1;
    } else if (targetPlacement === 'overlay') {
      state.adv.circular_track_slots.splice(axis, 0, slot);
      state.adv.circular_track_slots_axis_index = axis;
    } else {
      state.adv.circular_track_slots.splice(axis, 0, slot);
      state.adv.circular_track_slots_axis_index = axis;
    }
    normalizeSlotsInPlace();
  };

  const moveCircularTrackSlotOutside = (index) => {
    if (!canMoveCircularTrackSlotOutside(index)) return;
    moveCircularTrackSlotToPlacement(index, 'outside');
  };

  const moveCircularTrackSlotInside = (index) => {
    if (!canMoveCircularTrackSlotInside(index)) return;
    moveCircularTrackSlotToPlacement(index, 'inside');
  };

  const moveCircularTrackSlotToAxis = (index) => {
    if (!canMoveCircularTrackSlotToAxis(index)) return;
    moveCircularTrackSlotToPlacement(index, 'overlay');
  };

  const updateCircularTrackSlotRenderer = (slot, renderer) => {
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
      slot.params.nt = normalizeNt(slot.params.nt, normalizeNt(state.adv.nt));
    }
    normalizeSlotsInPlace();
  };

  const updateCircularTrackSlotPlacement = (slot, placement) => {
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
    normalizeSlotsInPlace();
  };

  const updateCircularTrackFeatureLane = (slot, laneDirection) => {
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
    const index = state.adv.circular_track_slots.findIndex((candidate) => candidate === slot);
    if (index >= 0) {
      moveCircularTrackSlotToPlacement(index, sideForLaneDirection(normalizedLaneDirection));
      return;
    }
    slot.params.lane_direction = normalizedLaneDirection;
    slot.side = sideForLaneDirection(slot.params.lane_direction);
    normalizeSlotsInPlace();
  };

  const circularTrackRendererLabel = (renderer) => RENDERER_LABELS[renderer] || renderer;
  const supportsCircularTrackSlotPlacement = (renderer) => SUPPORTED_RENDERERS.includes(renderer);

  const circularTrackSlots = () => (
    Array.isArray(state.adv.circular_track_slots)
      ? state.adv.circular_track_slots.map((slot, index) => ({ kind: STACK_ENTRY_SLOT, slot, index }))
      : []
  );

  const circularTrackStackEntries = () => {
    const slots = Array.isArray(state.adv.circular_track_slots) ? state.adv.circular_track_slots : [];
    const axisIndex = axisIndexForCurrentSlots(slots);
    const entries = [];
    let axisRendered = false;
    slots.forEach((slot, index) => {
      const onAxis = effectiveSlotPlacement(slot, state.form.track_type) === 'overlay';
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
    normalizeSlotsInPlace();
    return buildCircularTrackSlotSpec(slot, state.adv.nt, state.form.track_type);
  };

  const circularTrackSlotAutoPlacementFromOrder = (slot) => {
    void slot;
    return null;
  };

  const circularTrackPresetSummary = () => {
    const preset = normalizeCircularTrackPreset(state.form.track_type);
    const lengthParam = getPreviewLengthParam(state);
    const lane = laneDirectionForPreset(preset);
    const pieces = [
      laneDirectionLabel(lane),
      `feature r ${previewFeatureRadiusRatio(preset, lengthParam, state).toFixed(2)}x`
    ];
    const slotLike = (id, renderer) => ({ id, renderer });
    if (Boolean(state.form.show_depth)) {
      const depthRatio = getPresetRadiusRatio(slotLike('depth', 'depth'), 'depth', preset, lengthParam, state);
      if (depthRatio !== null) pieces.push(`depth r ${depthRatio.toFixed(2)}x`);
    }
    if (!Boolean(state.form.suppress_gc)) {
      const gcRatio = getPresetRadiusRatio(slotLike('gc_content', 'dinucleotide_content'), 'dinucleotide_content', preset, lengthParam, state);
      if (gcRatio !== null) pieces.push(`GC r ${gcRatio.toFixed(2)}x`);
    }
    if (!Boolean(state.form.suppress_skew)) {
      const skewRatio = getPresetRadiusRatio(slotLike('gc_skew', 'dinucleotide_skew'), 'dinucleotide_skew', preset, lengthParam, state);
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
    circularTrackRenderers: SUPPORTED_RENDERERS,
    circularTrackRendererLabel,
    normalizeCircularTrackSlots: normalizeSlotsInPlace,
    resetCircularTrackSlotsFromSimpleControls,
    ensureCircularTrackDepthSlot,
    resetCircularTrackSlotsToPreset,
    applyCircularTrackPreset,
    setCircularTrackSlotsEnabled,
    addCircularTrackSlot,
    duplicateCircularTrackSlot,
    removeCircularTrackSlot,
    moveCircularTrackSlot,
    canMoveCircularTrackSlot,
    moveCircularTrackSlotOutside,
    moveCircularTrackSlotInside,
    moveCircularTrackSlotToAxis,
    canMoveCircularTrackSlotOutside,
    canMoveCircularTrackSlotInside,
    canMoveCircularTrackSlotToAxis,
    updateCircularTrackSlotRenderer,
    updateCircularTrackSlotPlacement,
    updateCircularTrackFeatureLane,
    supportsCircularTrackSlotPlacement,
    circularTrackSlots,
    circularTrackStackEntries,
    circularTrackSlotCliSpec,
    circularTrackPresetSummary,
    circularTrackSlotUsesPresetGeometry
  };
};
