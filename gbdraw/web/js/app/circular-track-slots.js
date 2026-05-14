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

const PLACEMENT_RENDERERS = new Set(['dinucleotide_content', 'dinucleotide_skew', 'depth', 'spacer']);
const PLACEMENT_LABELS = {
  inside: 'Inside',
  outside: 'Outside',
  overlay: 'Overlay'
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

export const findFeatureSlotIndex = (slots) => {
  if (!Array.isArray(slots)) return -1;
  return slots.findIndex((slot) => slot?.renderer === 'features');
};

const getSlotPlacement = (slot) => {
  if (!slot || !PLACEMENT_RENDERERS.has(slot.renderer)) return null;
  return normalizePlacement(slot.params?.side, 'inside');
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

const applyPlacementDefaults = (params, placement = 'inside') => {
  const side = normalizePlacement(placement);
  params.side = side;
  if (side === 'inside') {
    params.strict = true;
    params.compress = true;
  } else if (side === 'outside') {
    params.strict = false;
    delete params.compress;
  } else {
    params.strict = false;
    delete params.compress;
  }
};

const makeSlot = ({
  id,
  renderer,
  enabled = true,
  width = null,
  radius = null,
  innerRadius = null,
  outerRadius = null,
  gapAfter = null,
  z = 0,
  params = {}
}) => ({
  id,
  renderer,
  enabled,
  width,
  radius,
  innerRadius,
  outerRadius,
  gapAfter,
  z,
  params: cloneParams(params)
});

export const createDefaultCircularTrackSlots = ({
  nt = 'GC',
  showDepth = false,
  showGc = true,
  showSkew = true
} = {}) => {
  const normalizedNt = normalizeNt(nt);
  const slots = [
    makeSlot({ id: 'features', renderer: 'features' }),
    makeSlot({
      id: 'ticks',
      renderer: 'ticks',
      params: { label_side: 'legacy', tick_side: 'legacy' }
    })
  ];
  if (showDepth) slots.push(makeSlot({ id: 'depth', renderer: 'depth' }));
  if (showGc) {
    slots.push(makeSlot({
      id: 'gc_content',
      renderer: 'dinucleotide_content',
      params: { nt: normalizedNt }
    }));
  }
  if (showSkew) {
    slots.push(makeSlot({
      id: 'gc_skew',
      renderer: 'dinucleotide_skew',
      params: { nt: normalizedNt }
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
  if (normalizedRenderer === 'ticks') {
    params.label_side = 'outside';
    params.tick_side = 'inside';
  } else if (normalizedRenderer === 'dinucleotide_content' || normalizedRenderer === 'dinucleotide_skew') {
    params.nt = normalizeNt(nt);
  }
  if (PLACEMENT_RENDERERS.has(normalizedRenderer)) {
    applyPlacementDefaults(params, placement);
  }

  return makeSlot({ id, renderer: normalizedRenderer, params });
};

export const normalizeCircularTrackSlot = (slot, index = 0, defaultNt = 'GC') => {
  const source = slot && typeof slot === 'object' && !Array.isArray(slot) ? slot : {};
  const renderer = SUPPORTED_RENDERERS.includes(source.renderer) ? source.renderer : 'dinucleotide_skew';
  const fallbackId = DEFAULT_SLOT_IDS[renderer] || `slot_${index + 1}`;
  const params = cloneParams(source.params);

  if ((renderer === 'dinucleotide_content' || renderer === 'dinucleotide_skew') && !normalizeOptionalText(params.nt)) {
    params.nt = normalizeNt(defaultNt);
  }
  if (renderer === 'ticks') {
    delete params['axis'];
    if (!normalizeOptionalText(params.label_side)) params.label_side = 'outside';
    if (!normalizeOptionalText(params.tick_side)) params.tick_side = 'inside';
  }
  if (PLACEMENT_RENDERERS.has(renderer)) {
    if (normalizeOptionalText(params.side)) params.side = normalizePlacement(params.side);
    const strict = normalizeOptionalBool(params.strict);
    const compress = normalizeOptionalBool(params.compress);
    const reserve = normalizeOptionalBool(params.reserve);
    if (strict === null) delete params.strict;
    else params.strict = strict;
    if (compress === null) delete params.compress;
    else params.compress = compress;
    if (reserve === null) delete params.reserve;
    else params.reserve = reserve;
  } else {
    delete params.side;
    delete params.strict;
    delete params.compress;
    delete params.reserve;
  }

  return makeSlot({
    id: cleanToken(source.id, fallbackId),
    renderer,
    enabled: source.enabled !== false,
    width: source.width ?? null,
    radius: source.radius ?? source.r ?? null,
    innerRadius: source.innerRadius ?? source.inner_radius ?? null,
    outerRadius: source.outerRadius ?? source.outer_radius ?? null,
    gapAfter: source.gapAfter ?? source.gap_after ?? null,
    z: Number.isFinite(Number(source.z)) ? Number(source.z) : 0,
    params
  });
};

export const normalizeCircularTrackSlots = (slots, defaultNt = 'GC') => {
  const base = Array.isArray(slots)
    ? slots
    : createDefaultCircularTrackSlots({ nt: defaultNt });
  return base.map((slot, index) => normalizeCircularTrackSlot(slot, index, defaultNt));
};

const appendOption = (options, key, value) => {
  const text = normalizeOptionalText(value);
  if (text === null) return;
  options.push(`${key}=${text}`);
};

export const buildCircularTrackSlotSpec = (slot, defaultNt = 'GC') => {
  const normalized = normalizeCircularTrackSlot(slot, 0, defaultNt);
  const options = [];
  const params = normalized.params || {};

  if (!normalized.enabled) options.push('enabled=false');
  appendOption(options, 'w', normalized.width);
  appendOption(options, 'r', normalized.radius);
  appendOption(options, 'ri', normalized.innerRadius);
  appendOption(options, 'ro', normalized.outerRadius);
  appendOption(options, 'gap', normalized.gapAfter);
  if (Number.isFinite(Number(normalized.z)) && Number(normalized.z) !== 0) {
    options.push(`z=${Number(normalized.z)}`);
  }

  if (normalized.renderer === 'ticks') {
    appendOption(options, 'label_side', params.label_side);
    appendOption(options, 'tick_side', params.tick_side);
  } else if (normalized.renderer === 'dinucleotide_content' || normalized.renderer === 'dinucleotide_skew') {
    options.push(`nt=${normalizeNt(params.nt, normalizeNt(defaultNt))}`);
  }
  if (PLACEMENT_RENDERERS.has(normalized.renderer)) {
    appendOption(options, 'side', params.side);
    appendOption(options, 'strict', params.strict);
    appendOption(options, 'compress', params.compress);
    appendOption(options, 'reserve', params.reserve);
  }
  appendOption(options, 'legend_label', params.legend_label);

  return `${normalized.id}:${normalized.renderer}${options.length ? `@${options.join(',')}` : ''}`;
};

export const hasEnabledCircularTrackRenderer = (slots, renderer) =>
  normalizeCircularTrackSlots(slots).some((slot) => slot.enabled && slot.renderer === renderer);

export const createCircularTrackSlotEditor = ({ state }) => {
  const normalizeSlotsInPlace = () => {
    const normalized = normalizeCircularTrackSlots(state.adv.circular_track_slots, state.adv.nt);
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
      showSkew: !state.form.suppress_skew
    });
    state.adv.circular_track_slots.splice(0, state.adv.circular_track_slots.length, ...slots);
  };

  const addCircularTrackSlot = (renderer, placement = 'inside') => {
    normalizeSlotsInPlace();
    const slot = createCircularTrackSlotForRenderer(renderer, state.adv.circular_track_slots, state.adv.nt, placement);
    if (PLACEMENT_RENDERERS.has(slot.renderer)) {
      insertSlotRelativeToFeatures(slot, slot.params.side);
    } else {
      state.adv.circular_track_slots.push(slot);
    }
  };

  const duplicateCircularTrackSlot = (index) => {
    normalizeSlotsInPlace();
    const idx = Number(index);
    if (!Number.isInteger(idx) || idx < 0 || idx >= state.adv.circular_track_slots.length) return;
    const source = normalizeCircularTrackSlot(state.adv.circular_track_slots[idx], idx, state.adv.nt);
    const duplicate = createCircularTrackSlotForRenderer(source.renderer, state.adv.circular_track_slots, state.adv.nt);
    duplicate.enabled = source.enabled;
    duplicate.width = source.width;
    duplicate.radius = source.radius;
    duplicate.innerRadius = source.innerRadius;
    duplicate.outerRadius = source.outerRadius;
    duplicate.gapAfter = source.gapAfter;
    duplicate.z = source.z;
    duplicate.params = cloneParams(source.params);
    state.adv.circular_track_slots.splice(idx + 1, 0, duplicate);
    if (PLACEMENT_RENDERERS.has(duplicate.renderer)) {
      moveSlotRelativeToFeatures(duplicate, duplicate.params.side);
    }
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
    let placement = null;
    if (renderer === 'ticks') {
      delete slot.params['axis'];
      slot.params.label_side = normalizeOptionalText(slot.params.label_side) || 'outside';
      slot.params.tick_side = normalizeOptionalText(slot.params.tick_side) || 'inside';
    } else if (renderer === 'dinucleotide_content' || renderer === 'dinucleotide_skew') {
      slot.params.nt = normalizeNt(slot.params.nt, normalizeNt(state.adv.nt));
    }
    if (PLACEMENT_RENDERERS.has(renderer) && !normalizeOptionalText(slot.params.side)) {
      applyPlacementDefaults(slot.params, 'inside');
    }
    if (PLACEMENT_RENDERERS.has(renderer)) {
      placement = normalizePlacement(slot.params.side);
    }
    if (!PLACEMENT_RENDERERS.has(renderer)) {
      delete slot.params.side;
      delete slot.params.strict;
      delete slot.params.compress;
      delete slot.params.reserve;
    }
    if (placement !== null) moveSlotRelativeToFeatures(slot, placement);
  };

  const updateCircularTrackSlotPlacement = (slot, placement) => {
    if (!slot || !PLACEMENT_RENDERERS.has(slot.renderer)) return;
    slot.params = cloneParams(slot.params);
    applyPlacementDefaults(slot.params, placement);
    moveSlotRelativeToFeatures(slot, slot.params.side);
  };

  const circularTrackRendererLabel = (renderer) => RENDERER_LABELS[renderer] || renderer;
  const circularTrackPlacementLabel = (placement) => PLACEMENT_LABELS[normalizePlacement(placement)] || PLACEMENT_LABELS.inside;
  const supportsCircularTrackSlotPlacement = (renderer) => PLACEMENT_RENDERERS.has(renderer);

  const circularTrackSlotCliSpec = (slot) => buildCircularTrackSlotSpec(slot, state.adv.nt);

  return {
    circularTrackRenderers: SUPPORTED_RENDERERS,
    circularTrackRendererLabel,
    normalizeCircularTrackSlots: normalizeSlotsInPlace,
    resetCircularTrackSlotsFromSimpleControls,
    addCircularTrackSlot,
    duplicateCircularTrackSlot,
    removeCircularTrackSlot,
    moveCircularTrackSlot,
    updateCircularTrackSlotRenderer,
    updateCircularTrackSlotPlacement,
    circularTrackPlacementLabel,
    supportsCircularTrackSlotPlacement,
    circularTrackSlotCliSpec
  };
};
