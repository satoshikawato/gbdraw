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
const PLACEMENT_LABELS = {
  inside: 'Inside',
  outside: 'Outside',
  overlay: 'Overlay'
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

export const findFeatureSlotIndex = (slots) => {
  if (!Array.isArray(slots)) return -1;
  return slots.findIndex((slot) => slot?.renderer === 'features');
};

const getSlotPlacement = (slot) => {
  if (!slot) return null;
  return normalizePlacement(slot.side, 'inside');
};

const getPlacementInsertionIndex = (slots, placement) => {
  const featureIndex = findFeatureSlotIndex(slots);
  if (featureIndex < 0) return Array.isArray(slots) ? slots.length : 0;

  if (placement === 'outside') {
    let insertIndex = featureIndex;
    for (let index = 0; index < featureIndex; index += 1) {
      if (getSlotPlacement(slots[index]) === 'outside') insertIndex = index + 1;
    }
    return insertIndex;
  }

  if (placement === 'inside') {
    let insertIndex = featureIndex + 1;
    for (let index = featureIndex + 1; index < slots.length; index += 1) {
      if (getSlotPlacement(slots[index]) === 'inside') insertIndex = index + 1;
    }
    return insertIndex;
  }

  return slots.length;
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
  const side = normalizePlacement(placement);
  slot.side = side;
  slot.params = cloneParams(slot.params);
  if (slot.renderer === 'features') {
    slot.params.lane_direction = laneDirectionForSide(side);
    if (side === 'overlay' && slot.reserve === null) slot.reserve = true;
  }
  if ((!NUMERIC_RENDERERS.has(slot.renderer) || side !== 'inside') && slot.compress === true) {
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
  } else if (normalizedRenderer === 'features') {
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
  applyPlacementDefaults(slot, side);
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
      params.lane_direction = normalizeLaneDirection(rawLaneDirection);
      side = sideForLaneDirection(params.lane_direction);
      delete params.lanes;
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
  return base.map((slot, index) => normalizeCircularTrackSlot(slot, index, defaultNt, preset));
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
    state.adv.circular_track_slots.splice(0, state.adv.circular_track_slots.length, ...normalized);
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

  const insertSlotRelativeToFeatures = (slot, placement) => {
    const normalizedPlacement = normalizePlacement(placement);
    if (normalizedPlacement === 'overlay') {
      state.adv.circular_track_slots.push(slot);
      return;
    }
    const insertIndex = getPlacementInsertionIndex(state.adv.circular_track_slots, normalizedPlacement);
    state.adv.circular_track_slots.splice(insertIndex, 0, slot);
  };

  const moveSlotRelativeToFeatures = (slotOrId, placement) => {
    const normalizedPlacement = normalizePlacement(placement);
    if (normalizedPlacement === 'overlay') return;
    const index = findSlotIndexByIdentity(slotOrId);
    if (index < 0) return;
    const [slot] = state.adv.circular_track_slots.splice(index, 1);
    insertSlotRelativeToFeatures(slot, normalizedPlacement);
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
    insertSlotRelativeToFeatures(slot, slot.side);
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
    moveSlotRelativeToFeatures(duplicate, duplicate.side);
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
    moveSlotRelativeToFeatures(slot, slot.side);
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
      return;
    }
    applyPlacementDefaults(slot, placement);
    moveSlotRelativeToFeatures(slot, slot.side);
  };

  const updateCircularTrackFeatureLane = (slot, laneDirection) => {
    if (!slot || slot.renderer !== 'features') return;
    slot.params = cloneParams(slot.params);
    const explicitLaneDirection = normalizeOptionalText(laneDirection);
    if (explicitLaneDirection === null) {
      delete slot.params.lane_direction;
      slot.side = null;
      slot.reserve = null;
      return;
    }
    slot.params.lane_direction = normalizeLaneDirection(explicitLaneDirection);
    slot.side = sideForLaneDirection(slot.params.lane_direction);
    if (slot.side === 'overlay' && slot.reserve === null) slot.reserve = true;
    moveSlotRelativeToFeatures(slot, slot.side);
  };

  const circularTrackRendererLabel = (renderer) => RENDERER_LABELS[renderer] || renderer;
  const circularTrackPlacementLabel = (placement) => (
    normalizeOptionalPlacement(placement) === null
      ? 'Preset'
      : (PLACEMENT_LABELS[normalizePlacement(placement)] || PLACEMENT_LABELS.inside)
  );
  const supportsCircularTrackSlotPlacement = (renderer) => SUPPORTED_RENDERERS.includes(renderer);

  const circularTrackSlotCliSpec = (slot) => buildCircularTrackSlotSpec(slot, state.adv.nt, state.form.track_type);

  return {
    circularTrackRenderers: SUPPORTED_RENDERERS,
    circularTrackRendererLabel,
    normalizeCircularTrackSlots: normalizeSlotsInPlace,
    resetCircularTrackSlotsFromSimpleControls,
    applyCircularTrackPreset,
    addCircularTrackSlot,
    duplicateCircularTrackSlot,
    removeCircularTrackSlot,
    moveCircularTrackSlot,
    updateCircularTrackSlotRenderer,
    updateCircularTrackSlotPlacement,
    updateCircularTrackFeatureLane,
    circularTrackPlacementLabel,
    supportsCircularTrackSlotPlacement,
    circularTrackSlotCliSpec
  };
};
