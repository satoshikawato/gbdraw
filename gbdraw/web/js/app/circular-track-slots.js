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

const cloneParams = (params = {}) => {
  if (!params || typeof params !== 'object' || Array.isArray(params)) return {};
  return { ...params };
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
      params: { axis: true, label_side: 'legacy', tick_side: 'legacy' }
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

export const createCircularTrackSlotForRenderer = (renderer, existingSlots = [], nt = 'GC') => {
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
    params.axis = true;
    params.label_side = 'outside';
    params.tick_side = 'inside';
  } else if (normalizedRenderer === 'dinucleotide_content' || normalizedRenderer === 'dinucleotide_skew') {
    params.nt = normalizeNt(nt);
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
    if (!normalizeOptionalText(params.label_side)) params.label_side = 'outside';
    if (!normalizeOptionalText(params.tick_side)) params.tick_side = 'inside';
    if (typeof params.axis !== 'boolean') params.axis = params.axis !== false;
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

const appendBooleanOption = (options, key, value) => {
  if (typeof value !== 'boolean') return;
  options.push(`${key}=${value ? 'true' : 'false'}`);
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
    appendBooleanOption(options, 'axis', params.axis);
    appendOption(options, 'label_side', params.label_side);
    appendOption(options, 'tick_side', params.tick_side);
  } else if (normalized.renderer === 'dinucleotide_content' || normalized.renderer === 'dinucleotide_skew') {
    options.push(`nt=${normalizeNt(params.nt, normalizeNt(defaultNt))}`);
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

  const resetCircularTrackSlotsFromSimpleControls = () => {
    const slots = createDefaultCircularTrackSlots({
      nt: state.adv.nt,
      showDepth: Boolean(state.form.show_depth),
      showGc: !state.form.suppress_gc,
      showSkew: !state.form.suppress_skew
    });
    state.adv.circular_track_slots.splice(0, state.adv.circular_track_slots.length, ...slots);
  };

  const addCircularTrackSlot = (renderer) => {
    normalizeSlotsInPlace();
    state.adv.circular_track_slots.push(
      createCircularTrackSlotForRenderer(renderer, state.adv.circular_track_slots, state.adv.nt)
    );
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
    if (renderer === 'ticks') {
      slot.params.axis = slot.params.axis !== false;
      slot.params.label_side = normalizeOptionalText(slot.params.label_side) || 'outside';
      slot.params.tick_side = normalizeOptionalText(slot.params.tick_side) || 'inside';
    } else if (renderer === 'dinucleotide_content' || renderer === 'dinucleotide_skew') {
      slot.params.nt = normalizeNt(slot.params.nt, normalizeNt(state.adv.nt));
    }
  };

  const circularTrackRendererLabel = (renderer) => RENDERER_LABELS[renderer] || renderer;

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
    circularTrackSlotCliSpec
  };
};
