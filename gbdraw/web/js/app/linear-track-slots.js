const SUPPORTED_RENDERERS = [
  'features',
  'dinucleotide_content',
  'dinucleotide_skew',
  'depth',
  'spacer'
];

const RENDERER_LABELS = {
  features: 'Features',
  dinucleotide_content: 'Dinucleotide content',
  dinucleotide_skew: 'Dinucleotide skew',
  depth: 'Depth',
  spacer: 'Spacer'
};

const RENDERER_ALIASES = {
  gc_content: 'dinucleotide_content',
  content: 'dinucleotide_content',
  gc_skew: 'dinucleotide_skew',
  skew: 'dinucleotide_skew'
};

const DEFAULT_SLOT_IDS = {
  features: 'features',
  dinucleotide_content: 'gc_content',
  dinucleotide_skew: 'gc_skew',
  depth: 'depth',
  spacer: 'spacer'
};

const STACK_ENTRY_AXIS = 'axis';
const STACK_ENTRY_SLOT = 'slot';

const cloneParams = (params = {}) => {
  if (!params || typeof params !== 'object' || Array.isArray(params)) return {};
  return { ...params };
};

const normalizeOptionalText = (value) => {
  const text = String(value ?? '').trim();
  return text.length > 0 ? text : null;
};

const normalizeRenderer = (value, fallback = 'features') => {
  const text = String(value || fallback).trim().toLowerCase();
  const renderer = RENDERER_ALIASES[text] || text;
  return SUPPORTED_RENDERERS.includes(renderer) ? renderer : fallback;
};

const normalizeSide = (value, fallback = 'below') => {
  const text = String(value || fallback).trim().toLowerCase();
  return ['above', 'below', 'overlay'].includes(text) ? text : fallback;
};

const sideForLinearTrackLayout = (trackLayout = 'middle') => {
  const normalized = String(trackLayout || 'middle').trim().toLowerCase();
  if (normalized === 'above' || normalized === 'spreadout') return 'above';
  if (normalized === 'below' || normalized === 'tuckin') return 'below';
  return 'overlay';
};

const normalizeNt = (value, fallback = 'GC') => {
  const text = String(value || '').trim().toUpperCase();
  return text || fallback;
};

const normalizeTrackIndex = (value) => {
  if (value === null || value === undefined || value === '') return null;
  const numeric = Number(value);
  if (!Number.isInteger(numeric) || numeric < 0) return null;
  return numeric;
};

const normalizePxText = (value) => {
  const text = normalizeOptionalText(value);
  if (text === null) return '';
  const withoutUnit = text.endsWith('px') ? text.slice(0, -2) : text;
  const numeric = Number(withoutUnit);
  return Number.isFinite(numeric) && numeric >= 0 ? String(text) : '';
};

const defaultSlot = (renderer, overrides = {}) => {
  const normalizedRenderer = normalizeRenderer(renderer);
  const params = cloneParams(overrides.params);
  return {
    id: String(overrides.id || DEFAULT_SLOT_IDS[normalizedRenderer] || normalizedRenderer),
    renderer: normalizedRenderer,
    enabled: overrides.enabled !== false,
    side: normalizeSide(overrides.side, normalizedRenderer === 'features' ? 'overlay' : 'below'),
    height: normalizePxText(overrides.height),
    spacing: normalizePxText(overrides.spacing),
    z: Number.isInteger(Number(overrides.z)) ? Number(overrides.z) : 0,
    params
  };
};

export const linearDepthTrackCountForState = (state) => {
  if (!Boolean(state?.form?.show_depth)) return 0;
  const seqs = Array.isArray(state?.linearSeqs) ? state.linearSeqs : [];
  const counts = seqs.map((seq) => {
    const value = seq?.depth;
    if (Array.isArray(value)) return value.filter(Boolean).length;
    return value ? 1 : 0;
  });
  return Math.max(1, ...counts);
};

export const createDefaultLinearTrackSlots = ({
  showDepth = false,
  depthTrackCount = 1,
  showGc = false,
  showSkew = false,
  nt = 'GC',
  trackLayout = 'middle'
} = {}) => {
  const slots = [
    defaultSlot('features', {
      id: 'features',
      side: sideForLinearTrackLayout(trackLayout)
    })
  ];
  if (showDepth) {
    const count = Math.max(1, Number(depthTrackCount) || 1);
    for (let index = 0; index < count; index += 1) {
      slots.push(defaultSlot('depth', {
        id: count === 1 ? 'depth' : `depth_${index + 1}`,
        side: 'below',
        params: { track_index: index }
      }));
    }
  }
  if (showGc) {
    slots.push(defaultSlot('dinucleotide_content', {
      id: 'gc_content',
      side: 'below',
      params: { nt: normalizeNt(nt) }
    }));
  }
  if (showSkew) {
    slots.push(defaultSlot('dinucleotide_skew', {
      id: 'gc_skew',
      side: 'below',
      params: { nt: normalizeNt(nt) }
    }));
  }
  return slots;
};

export const clampLinearTrackAxisIndex = (value, slotCount) => {
  if (value === null || value === undefined || value === '') return null;
  const numeric = Number(value);
  if (!Number.isInteger(numeric)) return null;
  return Math.max(0, Math.min(Number(slotCount) || 0, numeric));
};

export const normalizeLinearTrackSlots = (slots, nt = 'GC', trackLayout = 'middle') => {
  const source = Array.isArray(slots) && slots.length > 0
    ? slots
    : createDefaultLinearTrackSlots({ showGc: true, showSkew: true, nt, trackLayout });
  const usedIds = new Set();
  let hasFeature = false;
  return source
    .filter((slot) => slot && typeof slot === 'object' && !Array.isArray(slot))
    .map((slot, index) => {
      const renderer = normalizeRenderer(slot.renderer);
      const params = cloneParams(slot.params);
      if (renderer === 'features') {
        hasFeature = true;
      }
      if (renderer === 'depth') {
        const trackIndex = normalizeTrackIndex(params.track_index);
        params.track_index = trackIndex === null ? 0 : trackIndex;
      }
      if (renderer === 'dinucleotide_content' || renderer === 'dinucleotide_skew') {
        params.nt = normalizeNt(params.nt ?? params.dinucleotide, nt);
        delete params.dinucleotide;
      }
      let id = String(slot.id || DEFAULT_SLOT_IDS[renderer] || `slot_${index + 1}`).trim();
      if (!id) id = `slot_${index + 1}`;
      if (usedIds.has(id)) id = `${id}_${index + 1}`;
      usedIds.add(id);
      const side = renderer === 'features'
        ? normalizeSide(slot.side, sideForLinearTrackLayout(trackLayout))
        : normalizeSide(slot.side, 'below');
      return {
        id,
        renderer,
        enabled: slot.enabled !== false,
        side: renderer !== 'features' && side === 'overlay' ? 'below' : side,
        height: normalizePxText(slot.height),
        spacing: normalizePxText(slot.spacing),
        z: Number.isInteger(Number(slot.z)) ? Number(slot.z) : 0,
        params
      };
    })
    .filter((slot) => slot.renderer !== 'features' || (hasFeature && slot.id));
};

export const hasEnabledLinearTrackRenderer = (slots, renderer) => {
  const normalizedRenderer = normalizeRenderer(renderer);
  return normalizeLinearTrackSlots(slots).some(
    (slot) => slot.enabled !== false && slot.renderer === normalizedRenderer
  );
};

export const buildLinearTrackSlotSpec = (slot, { includeEnabled = false } = {}) => {
  const normalized = normalizeLinearTrackSlots([slot])[0];
  if (!normalized) return '';
  const parts = [];
  if (includeEnabled && normalized.enabled === false) parts.push('enabled=false');
  if (normalized.side) parts.push(`side=${normalized.side}`);
  if (normalizeOptionalText(normalized.height)) parts.push(`h=${normalized.height}`);
  if (normalizeOptionalText(normalized.spacing)) parts.push(`spacing=${normalized.spacing}`);
  if (Number(normalized.z) !== 0) parts.push(`z=${Number(normalized.z)}`);
  const params = cloneParams(normalized.params);
  if (normalized.renderer === 'depth') {
    const trackIndex = normalizeTrackIndex(params.track_index);
    parts.push(`track_index=${trackIndex === null ? 0 : trackIndex}`);
  }
  if (normalized.renderer === 'dinucleotide_content' || normalized.renderer === 'dinucleotide_skew') {
    parts.push(`nt=${normalizeNt(params.nt)}`);
  }
  if (normalizeOptionalText(params.legend_label)) {
    parts.push(`legend_label=${String(params.legend_label).trim()}`);
  }
  const suffix = parts.length > 0 ? `@${parts.join(',')}` : '';
  return `${normalized.id}:${normalized.renderer}${suffix}`;
};

export const applyLinearTrackOrderPlacements = (slots, axisIndex = null) => {
  const normalized = normalizeLinearTrackSlots(slots);
  const clampedAxis = clampLinearTrackAxisIndex(axisIndex, normalized.length);
  if (clampedAxis === null) return normalized;
  return normalized.map((slot, index) => ({
    ...slot,
    side: slot.side === 'overlay' ? 'overlay' : (index < clampedAxis ? 'above' : 'below')
  }));
};

export const createLinearTrackSlotEditor = ({ state }) => {
  const { computed } = window.Vue;
  const { adv, form } = state;

  const linearTrackRenderers = SUPPORTED_RENDERERS.slice();
  const linearTrackRendererLabel = (renderer) => RENDERER_LABELS[normalizeRenderer(renderer)] || String(renderer || '');

  const normalizeCurrentSlots = () => {
    const normalized = normalizeLinearTrackSlots(adv.linear_track_slots, adv.nt, form.linear_track_layout);
    adv.linear_track_slots.splice(0, adv.linear_track_slots.length, ...normalized);
    adv.linear_track_slots_axis_index = clampLinearTrackAxisIndex(
      adv.linear_track_slots_axis_index,
      adv.linear_track_slots.length
    );
  };

  const resetLinearTrackSlotsFromSimpleControls = () => {
    const slots = createDefaultLinearTrackSlots({
      showDepth: Boolean(form.show_depth),
      depthTrackCount: linearDepthTrackCountForState(state),
      showGc: Boolean(form.show_gc),
      showSkew: Boolean(form.show_skew),
      nt: adv.nt,
      trackLayout: form.linear_track_layout
    });
    adv.linear_track_slots.splice(0, adv.linear_track_slots.length, ...slots);
    adv.linear_track_slots_axis_index = slots.filter((slot) => slot.side === 'above').length;
  };

  const setLinearTrackSlotsEnabled = (enabled) => {
    adv.linear_track_slots_enabled = Boolean(enabled);
    if (adv.linear_track_slots_enabled) {
      resetLinearTrackSlotsFromSimpleControls();
    }
  };

  const addLinearTrackSlot = (renderer = 'spacer') => {
    normalizeCurrentSlots();
    const normalizedRenderer = normalizeRenderer(renderer, 'spacer');
    const baseId = DEFAULT_SLOT_IDS[normalizedRenderer] || 'slot';
    let nextId = baseId;
    let suffix = 2;
    const usedIds = new Set(adv.linear_track_slots.map((slot) => String(slot.id || '')));
    while (usedIds.has(nextId)) {
      nextId = `${baseId}_${suffix}`;
      suffix += 1;
    }
    adv.linear_track_slots.push(defaultSlot(normalizedRenderer, {
      id: nextId,
      side: normalizedRenderer === 'spacer' ? 'below' : undefined,
      height: normalizedRenderer === 'spacer' ? '12px' : ''
    }));
  };

  const duplicateLinearTrackSlot = (index) => {
    normalizeCurrentSlots();
    const idx = Number(index);
    const source = adv.linear_track_slots[idx];
    if (!source) return;
    const copy = {
      ...source,
      id: `${source.id || source.renderer}_copy`,
      params: cloneParams(source.params)
    };
    adv.linear_track_slots.splice(idx + 1, 0, copy);
    normalizeCurrentSlots();
  };

  const removeLinearTrackSlot = (index) => {
    normalizeCurrentSlots();
    const idx = Number(index);
    if (!Number.isInteger(idx) || idx < 0 || idx >= adv.linear_track_slots.length) return;
    adv.linear_track_slots.splice(idx, 1);
    normalizeCurrentSlots();
  };

  const canMoveLinearTrackSlot = (index, direction) => {
    const idx = Number(index);
    const target = idx + Math.sign(Number(direction));
    return Number.isInteger(idx) && target >= 0 && target < adv.linear_track_slots.length;
  };

  const moveLinearTrackSlot = (index, direction) => {
    normalizeCurrentSlots();
    const idx = Number(index);
    const target = idx + Math.sign(Number(direction));
    if (!canMoveLinearTrackSlot(idx, direction)) return;
    const [slot] = adv.linear_track_slots.splice(idx, 1);
    adv.linear_track_slots.splice(target, 0, slot);
    adv.linear_track_slots_axis_index = clampLinearTrackAxisIndex(
      adv.linear_track_slots_axis_index,
      adv.linear_track_slots.length
    );
  };

  const updateLinearTrackSlotRenderer = (slot) => {
    if (!slot) return;
    slot.renderer = normalizeRenderer(slot.renderer);
    slot.params = cloneParams(slot.params);
    if (slot.renderer === 'depth') {
      slot.params.track_index = normalizeTrackIndex(slot.params.track_index) ?? 0;
    }
    if (slot.renderer === 'dinucleotide_content' || slot.renderer === 'dinucleotide_skew') {
      slot.params.nt = normalizeNt(slot.params.nt, adv.nt);
    }
    if (slot.renderer !== 'features' && slot.side === 'overlay') slot.side = 'below';
    normalizeCurrentSlots();
  };

  const linearTrackStackEntries = computed(() => {
    const slots = Array.isArray(adv.linear_track_slots) ? adv.linear_track_slots : [];
    const axisIndex = clampLinearTrackAxisIndex(adv.linear_track_slots_axis_index, slots.length);
    const entries = [];
    slots.forEach((slot, index) => {
      if (axisIndex !== null && index === axisIndex) entries.push({ kind: STACK_ENTRY_AXIS });
      entries.push({ kind: STACK_ENTRY_SLOT, slot, index });
    });
    if (axisIndex === slots.length) entries.push({ kind: STACK_ENTRY_AXIS });
    if (axisIndex === null) entries.splice(1, 0, { kind: STACK_ENTRY_AXIS });
    return entries;
  });

  return {
    linearTrackRenderers,
    linearTrackRendererLabel,
    normalizeLinearTrackSlots: normalizeCurrentSlots,
    resetLinearTrackSlotsFromSimpleControls,
    setLinearTrackSlotsEnabled,
    addLinearTrackSlot,
    duplicateLinearTrackSlot,
    removeLinearTrackSlot,
    moveLinearTrackSlot,
    canMoveLinearTrackSlot,
    updateLinearTrackSlotRenderer,
    linearTrackSlots: () => {
      return Array.isArray(adv.linear_track_slots) ? adv.linear_track_slots : [];
    },
    linearTrackStackEntries: () => linearTrackStackEntries.value,
    linearTrackSlotCliSpec: (slot) => buildLinearTrackSlotSpec(slot),
    linearTrackSlotDisplayLabel: (slot) => linearTrackRendererLabel(slot?.renderer),
    linearTrackSlotDisplayMeta: (slot) => buildLinearTrackSlotSpec(slot)
  };
};

export { SUPPORTED_RENDERERS as LINEAR_TRACK_RENDERERS };
