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
const SECTION_DEFINITIONS = [
  {
    key: 'outer',
    label: 'Outer tracks',
    help: 'Tracks placed outside the genome axis.',
    addLabel: 'Add outer',
    emptyLabel: 'No outer tracks',
    renderers: SUPPORTED_RENDERERS
  },
  {
    key: 'axis',
    label: 'On-axis tracks',
    help: 'Tracks that overlap or attach directly to the genome axis.',
    addLabel: 'Add on-axis',
    emptyLabel: 'No on-axis tracks',
    renderers: SUPPORTED_RENDERERS
  },
  {
    key: 'inside',
    label: 'Inner tracks',
    help: 'Tracks placed inside the genome axis.',
    addLabel: 'Add inner',
    emptyLabel: 'No inner tracks',
    renderers: SUPPORTED_RENDERERS
  }
];
const SECTION_RANKS = {
  outer: 0,
  axis: 1,
  inside: 2
};
const PLACEMENT_LABELS = {
  inside: 'Inner',
  outside: 'Outer',
  overlay: 'On-axis'
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

const formatPresetName = (preset) => PRESET_LABELS[normalizeCircularTrackPreset(preset)] || PRESET_LABELS.tuckin;

const formatPx = (value) => {
  const number = Number(value);
  if (!Number.isFinite(number)) return String(value || '');
  const rounded = Math.round(number * 10) / 10;
  return `${Number.isInteger(rounded) ? String(rounded) : rounded.toFixed(1)}px`;
};

const formatRatioPx = (ratio) => {
  const number = Number(ratio);
  if (!Number.isFinite(number)) return 'auto';
  return `${number.toFixed(2)}x / ${formatPx(number * PREVIEW_RADIUS_PX)}`;
};

const laneDirectionLabel = (laneDirection) => {
  const lane = normalizeLaneDirection(laneDirection);
  if (lane === 'outside') return 'Outer lanes';
  if (lane === 'split') return 'Split lanes';
  return 'Inside lanes';
};

const rendererUsesNumericDefaults = (renderer) => NUMERIC_RENDERERS.has(renderer);

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
  if (renderer === 'features' && String(slot?.id || '').trim() === 'features') return 1.0;
  const trackId = getBuiltinTrackId(slot, renderer, state);
  if (trackId === null) return null;
  return PREVIEW_TRACK_DICT[lengthParam]?.[normalizeCircularTrackPreset(preset)]?.[trackId] ?? null;
};

const slotHasManualGeometry = (slot) => (
  normalizeOptionalText(slot?.width) !== null ||
  normalizeOptionalText(slot?.radius) !== null ||
  normalizeOptionalText(slot?.spacing) !== null
);

const getSlotRenderer = (slot) => {
  const renderer = String(slot?.renderer || '').trim();
  return SUPPORTED_RENDERERS.includes(renderer) ? renderer : 'dinucleotide_skew';
};

const getSectionKeyForSide = (side) => {
  if (side === 'outside') return 'outer';
  if (side === 'overlay') return 'axis';
  return 'inside';
};

const getSlotSectionKey = (slot, preset = 'tuckin') => {
  const renderer = getSlotRenderer(slot);
  const side = normalizeOptionalPlacement(slot?.side);
  if (side !== null) return getSectionKeyForSide(side);

  if (renderer === 'features') {
    const rawLane = normalizeOptionalText(slot?.params?.lane_direction ?? slot?.params?.lanes);
    const lane = rawLane === null ? laneDirectionForPreset(preset) : normalizeLaneDirection(rawLane);
    return getSectionKeyForSide(sideForLaneDirection(lane));
  }

  if (renderer === 'ticks') {
    const tickPreset = normalizeCircularTrackPreset(slot?.params?.preset, preset);
    return tickPreset === 'spreadout' ? 'axis' : 'inside';
  }

  return 'inside';
};

const getSectionDefinition = (sectionKey) => (
  SECTION_DEFINITIONS.find((definition) => definition.key === sectionKey) || SECTION_DEFINITIONS[2]
);

const normalizeRendererForSection = (sectionKey, renderer) => {
  const definition = getSectionDefinition(sectionKey);
  return definition.renderers.includes(renderer) ? renderer : definition.renderers[0];
};

const orderSlotsForSections = (slots, preset = 'tuckin') => {
  if (!Array.isArray(slots)) return [];
  return slots
    .map((slot, index) => ({ slot, index }))
    .sort((left, right) => {
      const leftRank = SECTION_RANKS[getSlotSectionKey(left.slot, preset)] ?? SECTION_RANKS.inside;
      const rightRank = SECTION_RANKS[getSlotSectionKey(right.slot, preset)] ?? SECTION_RANKS.inside;
      return leftRank === rightRank ? left.index - right.index : leftRank - rightRank;
    })
    .map((entry) => entry.slot);
};

const normalizeOptionalBool = (value) => {
  if (value === null || value === undefined || value === '') return null;
  if (typeof value === 'boolean') return value;
  const text = String(value).trim().toLowerCase();
  if (['1', 'true', 'yes', 'on'].includes(text)) return true;
  if (['0', 'false', 'no', 'off'].includes(text)) return false;
  return null;
};

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
    if (requestedSide === 'overlay' && slot.reserve === null) slot.reserve = true;
  }
  if ((!NUMERIC_RENDERERS.has(slot.renderer) || requestedSide !== 'inside') && slot.compress === true) {
    slot.compress = null;
  }
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
  strict = null,
  compress = null,
  reserve = null,
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
  strict,
  compress,
  reserve,
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
  Number(source.z || 0) === 0 &&
  normalizeOptionalBool(source.strict) === null;

const isLegacyDefaultWebSlotShape = (source, renderer, defaultNt = 'GC', preset = 'tuckin') => {
  if (!source || typeof source !== 'object' || Array.isArray(source)) return false;
  if (!hasBlankSlotGeometry(source) || !hasDefaultSlotFlags(source)) return false;

  const normalizedPreset = normalizeCircularTrackPreset(preset);
  const normalizedId = String(source.id || '').trim();
  const side = normalizeOptionalPlacement(source.side);
  const reserve = normalizeOptionalBool(source.reserve);
  const compress = normalizeOptionalBool(source.compress);
  const params = cloneParams(source.params);

  if (renderer === 'features' && normalizedId === 'features') {
    const laneDirection = laneDirectionForPreset(normalizedPreset);
    const expectedReserve = laneDirection === 'split' ? true : null;
    return (
      side === sideForLaneDirection(laneDirection) &&
      reserve === expectedReserve &&
      compress === null &&
      paramsMatchExactly(params, { lane_direction: laneDirection })
    );
  }

  if (renderer === 'ticks' && normalizedId === 'ticks') {
    return (
      side === 'inside' &&
      reserve === null &&
      compress === null &&
      paramsMatchExactly(params, {
        label_side: 'outside',
        tick_side: 'inside',
        preset: normalizedPreset
      })
    );
  }

  if (renderer === 'depth' && normalizedId === 'depth') {
    return side === 'inside' && reserve === null && compress === true && paramsMatchExactly(params, {});
  }

  if (renderer === 'dinucleotide_content' && normalizedId === 'gc_content') {
    return (
      side === 'inside' &&
      reserve === null &&
      compress === true &&
      paramsMatchExactly(params, { nt: normalizeNt(defaultNt) })
    );
  }

  if (renderer === 'dinucleotide_skew' && normalizedId === 'gc_skew') {
    return (
      side === 'inside' &&
      reserve === null &&
      compress === true &&
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

export const createCircularTrackSlotForRenderer = (renderer, existingSlots = [], nt = 'GC', placement = 'inside') => {
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
    params.label_side = 'outside';
    params.tick_side = 'inside';
  } else if (normalizedRenderer === 'features' && side !== null) {
    params.lane_direction = laneDirectionForSide(side);
  }
  void nt;

  const slot = makeSlot({
    id,
    renderer: normalizedRenderer,
    side,
    reserve: normalizedRenderer === 'features' && side === 'overlay' ? true : null,
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
  const strict = inheritsPresetDefaults ? null : normalizeOptionalBool(source.strict);
  let compress = inheritsPresetDefaults ? null : normalizeOptionalBool(source.compress);
  let reserve = inheritsPresetDefaults ? null : normalizeOptionalBool(source.reserve);
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
    if (!normalizeOptionalText(params.label_side) || params.label_side === 'legacy') {
      delete params.label_side;
    }
    if (!normalizeOptionalText(params.tick_side) || params.tick_side === 'legacy') {
      delete params.tick_side;
    }
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
    if (side === 'overlay' && reserve === null) reserve = true;
  }
  if (side !== null && !NUMERIC_RENDERERS.has(renderer) && compress === true) {
    compress = null;
  }
  if (!NUMERIC_RENDERERS.has(renderer) || side === 'outside' || side === 'overlay' || normalizeOptionalText(radius)) {
    if (compress === true) compress = null;
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
    strict,
    compress,
    reserve,
    params
  });
};

export const normalizeCircularTrackSlots = (slots, defaultNt = 'GC', preset = 'tuckin') => {
  const base = Array.isArray(slots)
    ? slots
    : createDefaultCircularTrackSlots({ nt: defaultNt, preset });
  return orderSlotsForSections(
    base.map((slot, index) => normalizeCircularTrackSlot(slot, index, defaultNt, preset)),
    preset
  );
};

const appendOption = (options, key, value) => {
  const text = normalizeOptionalText(value);
  if (text === null) return;
  options.push(`${key}=${text}`);
};

export const buildCircularTrackSlotSpec = (slot, defaultNt = 'GC', preset = 'tuckin') => {
  const normalized = normalizeCircularTrackSlot(slot, 0, defaultNt, preset);
  const options = [];
  const params = normalized.params || {};
  const normalizedPreset = normalizeCircularTrackPreset(preset);

  if (!normalized.enabled) options.push('enabled=false');
  appendOption(options, 'w', normalized.width);
  appendOption(options, 'r', normalized.radius);
  appendOption(options, 'spacing', normalized.spacing);
  appendOption(options, 'side', normalized.side);
  appendOption(options, 'strict', normalized.strict);
  appendOption(options, 'compress', normalized.compress);
  appendOption(options, 'reserve', normalized.reserve);
  if (Number.isFinite(Number(normalized.z)) && Number(normalized.z) !== 0) {
    options.push(`z=${Number(normalized.z)}`);
  }

  if (normalized.renderer === 'ticks') {
    appendOption(options, 'label_side', params.label_side);
    appendOption(options, 'tick_side', params.tick_side);
    if (normalizeOptionalText(params.preset) !== null && normalizeCircularTrackPreset(params.preset) !== normalizedPreset) {
      appendOption(options, 'preset', params.preset);
    }
  } else if (normalized.renderer === 'features') {
    if (
      normalizeOptionalText(params.lane_direction) !== null &&
      normalizeLaneDirection(params.lane_direction) !== laneDirectionForPreset(normalizedPreset)
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
  const normalizeSlotsInPlace = () => {
    const normalized = normalizeCircularTrackSlots(state.adv.circular_track_slots, state.adv.nt, state.form.track_type);
    state.adv.circular_track_slots.splice(
      0,
      state.adv.circular_track_slots.length,
      ...orderSlotsForSections(normalized, state.form.track_type)
    );
  };

  const findSlotIndexByIdentity = (slotOrId) => {
    const slots = state.adv.circular_track_slots;
    if (slotOrId && typeof slotOrId === 'object') {
      const objectIndex = slots.findIndex((slot) => slot === slotOrId);
      if (objectIndex >= 0) return objectIndex;
    }
    const slotId = String(
      slotOrId && typeof slotOrId === 'object'
        ? slotOrId.id
        : slotOrId
    ).trim();
    if (!slotId) return -1;
    return slots.findIndex((slot) => String(slot?.id || '').trim() === slotId);
  };

  const insertSlotInSection = (slot, sectionKey = getSlotSectionKey(slot, state.form.track_type)) => {
    const targetRank = SECTION_RANKS[getSectionDefinition(sectionKey).key] ?? SECTION_RANKS.inside;
    let insertIndex = 0;
    state.adv.circular_track_slots.forEach((existingSlot, index) => {
      const rank = SECTION_RANKS[getSlotSectionKey(existingSlot, state.form.track_type)] ?? SECTION_RANKS.inside;
      if (rank <= targetRank) insertIndex = index + 1;
    });
    state.adv.circular_track_slots.splice(insertIndex, 0, slot);
  };

  const moveSlotToCurrentSection = (slotOrId) => {
    const index = findSlotIndexByIdentity(slotOrId);
    if (index < 0) return;
    const [slot] = state.adv.circular_track_slots.splice(index, 1);
    insertSlotInSection(slot, getSlotSectionKey(slot, state.form.track_type));
  };

  const resetCircularTrackSlotsFromSimpleControls = () => {
    const slots = createDefaultCircularTrackSlots({
      nt: state.adv.nt,
      showDepth: Boolean(state.form.show_depth),
      showGc: !state.form.suppress_gc,
      showSkew: !state.form.suppress_skew,
      preset: state.form.track_type
    });
    state.adv.circular_track_slots.splice(0, state.adv.circular_track_slots.length, ...slots);
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

    const slot = createCircularTrackSlotForRenderer('depth', state.adv.circular_track_slots, state.adv.nt, 'inside');
    insertSlotInSection(slot);
  };

  const applyCircularTrackPreset = (preset) => {
    const normalizedPreset = normalizeCircularTrackPreset(preset);
    const current = JSON.stringify(normalizeCircularTrackSlots(state.adv.circular_track_slots, state.adv.nt, state.form.track_type));
    const replacement = createDefaultCircularTrackSlots({
      nt: state.adv.nt,
      showDepth: Boolean(state.form.show_depth),
      showGc: !state.form.suppress_gc,
      showSkew: !state.form.suppress_skew,
      preset: normalizedPreset
    });
    const next = JSON.stringify(normalizeCircularTrackSlots(replacement, state.adv.nt, normalizedPreset));
    if (
      current !== next &&
      Array.isArray(state.adv.circular_track_slots) &&
      state.adv.circular_track_slots.length > 0 &&
      typeof window !== 'undefined' &&
      !window.confirm('Replace the current custom circular track slots with this preset?')
    ) {
      return;
    }
    state.form.track_type = normalizedPreset;
    state.adv.circular_track_slots.splice(0, state.adv.circular_track_slots.length, ...replacement);
  };

  const addCircularTrackSlot = (renderer, placement = 'inside') => {
    normalizeSlotsInPlace();
    const slot = createCircularTrackSlotForRenderer(renderer, state.adv.circular_track_slots, state.adv.nt, placement);
    insertSlotInSection(slot);
  };

  const addCircularTrackSlotToSection = (sectionKey, renderer) => {
    normalizeSlotsInPlace();
    const normalizedSectionKey = getSectionDefinition(sectionKey).key;
    const normalizedRenderer = normalizeRendererForSection(normalizedSectionKey, renderer);
    const placement = normalizedSectionKey === 'outer'
      ? 'outside'
      : (normalizedSectionKey === 'axis' ? 'overlay' : 'inside');
    const slot = createCircularTrackSlotForRenderer(
      normalizedRenderer,
      state.adv.circular_track_slots,
      state.adv.nt,
      placement
    );
    insertSlotInSection(slot, normalizedSectionKey);
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
    duplicate.strict = source.strict;
    duplicate.compress = source.compress;
    duplicate.reserve = source.reserve;
    duplicate.params = cloneParams(source.params);
    state.adv.circular_track_slots.splice(idx + 1, 0, duplicate);
  };

  const removeCircularTrackSlot = (index) => {
    const idx = Number(index);
    if (!Number.isInteger(idx) || idx < 0 || idx >= state.adv.circular_track_slots.length) return;
    state.adv.circular_track_slots.splice(idx, 1);
  };

  const moveCircularTrackSlot = (fromIndex, toIndex) => {
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
  };

  const sameSectionSlotIndexes = (sectionKey) => state.adv.circular_track_slots
    .map((slot, index) => ({ slot, index }))
    .filter((entry) => getSlotSectionKey(entry.slot, state.form.track_type) === sectionKey)
    .map((entry) => entry.index);

  const canMoveCircularTrackSlotInSection = (index, direction) => {
    const idx = Number(index);
    const step = Number(direction);
    if (!Number.isInteger(idx) || !Number.isInteger(step) || step === 0) return false;
    const slot = state.adv.circular_track_slots[idx];
    if (!slot) return false;
    const indexes = sameSectionSlotIndexes(getSlotSectionKey(slot, state.form.track_type));
    const position = indexes.indexOf(idx);
    return position >= 0 && position + Math.sign(step) >= 0 && position + Math.sign(step) < indexes.length;
  };

  const moveCircularTrackSlotInSection = (index, direction) => {
    const idx = Number(index);
    const step = Math.sign(Number(direction));
    if (!canMoveCircularTrackSlotInSection(idx, step)) return;
    const slot = state.adv.circular_track_slots[idx];
    const indexes = sameSectionSlotIndexes(getSlotSectionKey(slot, state.form.track_type));
    const position = indexes.indexOf(idx);
    const target = indexes[position + step];
    const [moved] = state.adv.circular_track_slots.splice(idx, 1);
    const adjustedTarget = target > idx ? target - 1 : target;
    state.adv.circular_track_slots.splice(adjustedTarget, 0, moved);
  };

  const updateCircularTrackSlotRenderer = (slot, renderer) => {
    if (!slot || !SUPPORTED_RENDERERS.includes(renderer)) return;
    slot.renderer = renderer;
    slot.params = cloneParams(slot.params);
    slot.side = normalizeSlotSide(slot.side);
    if (renderer === 'ticks') {
      delete slot.params['axis'];
      if (normalizeOptionalText(slot.params.label_side) === null) delete slot.params.label_side;
      if (normalizeOptionalText(slot.params.tick_side) === null) delete slot.params.tick_side;
      if (normalizeOptionalText(slot.params.preset) === null) delete slot.params.preset;
      else slot.params.preset = normalizeCircularTrackPreset(slot.params.preset);
    } else if (renderer === 'features') {
      slot.params.lane_direction = laneDirectionForSide(slot.side);
      if (slot.side === 'overlay' && slot.reserve === null) slot.reserve = true;
    } else if (renderer === 'dinucleotide_content' || renderer === 'dinucleotide_skew') {
      slot.params.nt = normalizeNt(slot.params.nt, normalizeNt(state.adv.nt));
    }
    if (!NUMERIC_RENDERERS.has(renderer) && slot.compress === true) slot.compress = null;
    moveSlotToCurrentSection(slot);
  };

  const updateCircularTrackSlotPlacement = (slot, placement) => {
    if (!slot || !SUPPORTED_RENDERERS.includes(slot.renderer)) return;
    if (normalizeOptionalText(placement) === null) {
      slot.side = null;
      if (slot.compress === true) slot.compress = null;
      if (slot.renderer === 'features') {
        slot.params = cloneParams(slot.params);
        delete slot.params.lane_direction;
        slot.reserve = null;
      }
      moveSlotToCurrentSection(slot);
      return;
    }
    applyPlacementDefaults(slot, placement);
    moveSlotToCurrentSection(slot);
  };

  const updateCircularTrackFeatureLane = (slot, laneDirection) => {
    if (!slot || slot.renderer !== 'features') return;
    slot.params = cloneParams(slot.params);
    const explicitLaneDirection = normalizeOptionalText(laneDirection);
    if (explicitLaneDirection === null) {
      delete slot.params.lane_direction;
      slot.side = null;
      slot.reserve = null;
      moveSlotToCurrentSection(slot);
      return;
    }
    const normalizedLaneDirection = normalizeLaneDirection(explicitLaneDirection);
    if (normalizedLaneDirection === 'inside') {
      slot.params.lane_direction = 'inside';
      slot.side = 'inside';
      slot.reserve = null;
      moveSlotToCurrentSection(slot);
      return;
    }
    slot.params.lane_direction = normalizedLaneDirection;
    slot.side = sideForLaneDirection(slot.params.lane_direction);
    if (slot.side === 'overlay' && slot.reserve === null) slot.reserve = true;
    moveSlotToCurrentSection(slot);
  };

  const circularTrackRendererLabel = (renderer) => RENDERER_LABELS[renderer] || renderer;
  const circularTrackPlacementLabel = (placement) => (
    normalizeOptionalPlacement(placement) === null
      ? 'Auto'
      : (PLACEMENT_LABELS[normalizePlacement(placement)] || PLACEMENT_LABELS.inside)
  );
  const supportsCircularTrackSlotPlacement = (renderer) => SUPPORTED_RENDERERS.includes(renderer);

  const circularTrackRenderersForSection = (sectionKey) => [...getSectionDefinition(sectionKey).renderers];

  const circularTrackSlotSections = () => {
    const entries = Array.isArray(state.adv.circular_track_slots)
      ? state.adv.circular_track_slots.map((slot, index) => ({ slot, index }))
      : [];
    return SECTION_DEFINITIONS.map((definition) => ({
      ...definition,
      entries: entries.filter((entry) => getSlotSectionKey(entry.slot, state.form.track_type) === definition.key)
    }));
  };

  const circularTrackSlotSectionLabel = (slot) => getSectionDefinition(
    getSlotSectionKey(slot, state.form.track_type)
  ).label;

  const circularTrackSlotCliSpec = (slot) => buildCircularTrackSlotSpec(slot, state.adv.nt, state.form.track_type);

  const circularTrackPresetSummary = () => {
    const preset = normalizeCircularTrackPreset(state.form.track_type);
    const lengthParam = getPreviewLengthParam(state);
    const lane = laneDirectionForPreset(preset);
    const pieces = [`features ${laneDirectionLabel(lane)}`];
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

  const makeBadge = (key, label, value, source, title = '') => ({
    key,
    label,
    value,
    source,
    title
  });

  const circularTrackSlotEffectiveBadges = (slot) => {
    const normalized = normalizeCircularTrackSlot(slot, 0, state.adv.nt, state.form.track_type);
    const preset = normalizeCircularTrackPreset(state.form.track_type);
    const lengthParam = getPreviewLengthParam(state);
    const presetLabel = formatPresetName(preset);
    const badges = [];

    if (normalizeOptionalText(normalized.radius) !== null) {
      badges.push(makeBadge('radius', 'r', String(normalized.radius), 'manual', 'Manual radius override'));
    } else {
      const ratio = getPresetRadiusRatio(normalized, normalized.renderer, preset, lengthParam, state);
      const value = ratio === null
        ? (normalized.renderer === 'ticks' ? 'preset tick band' : 'auto pack')
        : formatRatioPx(ratio);
      badges.push(makeBadge('radius', 'r', value, 'preset', `Inherited from ${presetLabel}`));
    }

    if (normalizeOptionalText(normalized.width) !== null) {
      badges.push(makeBadge('width', 'w', String(normalized.width), 'manual', 'Manual width override'));
    } else {
      const widthPx = previewWidthPxForRenderer(normalized.renderer, lengthParam);
      const value = normalized.renderer === 'ticks' && widthPx <= 0
        ? 'preset ticks'
        : formatPx(widthPx);
      badges.push(makeBadge('width', 'w', value, 'preset', `Inherited from ${presetLabel}`));
    }

    if (normalizeOptionalText(normalized.spacing) !== null) {
      badges.push(makeBadge('spacing', 'gap', String(normalized.spacing), 'manual', 'Manual spacing override'));
    } else {
      badges.push(makeBadge('spacing', 'gap', formatPx(previewSpacingPx()), 'preset', `Inherited from ${presetLabel}`));
    }

    if (normalizeOptionalPlacement(slot?.side) !== null) {
      badges.push(makeBadge('side', 'side', circularTrackPlacementLabel(slot.side), 'manual', 'Manual placement override'));
    } else if (normalized.renderer === 'features') {
      badges.push(makeBadge('side', 'side', circularTrackPlacementLabel(sideForLaneDirection(laneDirectionForPreset(preset))), 'preset', `Inherited from ${presetLabel}`));
    } else if (normalized.renderer === 'ticks') {
      badges.push(makeBadge('side', 'side', 'preset ticks', 'preset', `Inherited from ${presetLabel}`));
    } else {
      badges.push(makeBadge('side', 'side', 'Inside', 'preset', `Inherited from ${presetLabel}`));
    }

    if (normalized.renderer === 'features') {
      const manualLane = normalizeOptionalText(slot?.params?.lane_direction);
      const lane = manualLane === null ? laneDirectionForPreset(preset) : normalizeLaneDirection(manualLane);
      badges.push(makeBadge(
        'lane',
        'lanes',
        laneDirectionLabel(lane),
        manualLane === null ? 'preset' : 'manual',
        manualLane === null ? `Inherited from ${presetLabel}` : 'Manual lane override'
      ));
    }

    if (normalized.renderer === 'ticks') {
      const manualLabelSide = normalizeOptionalText(slot?.params?.label_side);
      const manualTickSide = normalizeOptionalText(slot?.params?.tick_side);
      badges.push(makeBadge(
        'tick-label-side',
        'labels',
        manualLabelSide === null ? 'preset' : String(manualLabelSide),
        manualLabelSide === null ? 'preset' : 'manual',
        manualLabelSide === null ? `Inherited from ${presetLabel}` : 'Manual tick label placement'
      ));
      badges.push(makeBadge(
        'tick-side',
        'ticks',
        manualTickSide === null ? 'preset' : String(manualTickSide),
        manualTickSide === null ? 'preset' : 'manual',
        manualTickSide === null ? `Inherited from ${presetLabel}` : 'Manual tick placement'
      ));
    }

    if (rendererUsesNumericDefaults(normalized.renderer)) {
      const nt = normalizeOptionalText(slot?.params?.nt);
      if (nt !== null) {
        badges.push(makeBadge('nt', 'nt', normalizeNt(nt), 'manual', 'Manual dinucleotide override'));
      }
    }

    return badges;
  };

  const circularTrackSlotUsesPresetGeometry = (slot) => {
    if (!slot || typeof slot !== 'object') return false;
    if (!SUPPORTED_RENDERERS.includes(slot.renderer)) return false;
    if (!slotHasManualGeometry(slot)) return true;
    if (normalizeOptionalPlacement(slot.side) === null) return true;
    if (slot.renderer === 'features' && normalizeOptionalText(slot.params?.lane_direction) === null) return true;
    if (slot.renderer === 'ticks') {
      return (
        normalizeOptionalText(slot.params?.label_side) === null ||
        normalizeOptionalText(slot.params?.tick_side) === null
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
    applyCircularTrackPreset,
    addCircularTrackSlot,
    addCircularTrackSlotToSection,
    duplicateCircularTrackSlot,
    removeCircularTrackSlot,
    moveCircularTrackSlot,
    moveCircularTrackSlotInSection,
    canMoveCircularTrackSlotInSection,
    updateCircularTrackSlotRenderer,
    updateCircularTrackSlotPlacement,
    updateCircularTrackFeatureLane,
    circularTrackPlacementLabel,
    supportsCircularTrackSlotPlacement,
    circularTrackRenderersForSection,
    circularTrackSlotSections,
    circularTrackSlotSectionLabel,
    circularTrackSlotCliSpec,
    circularTrackPresetSummary,
    circularTrackSlotEffectiveBadges,
    circularTrackSlotUsesPresetGeometry
  };
};
