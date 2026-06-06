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
const DEFAULT_LINEAR_SLOT_MANAGER = 'linear-default';

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

const normalizePlacement = (value, fallback = 'below') => normalizeSide(value, fallback);

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

export const effectiveLinearSlotPlacement = (slot) => {
  if (!slot) return 'below';
  const renderer = normalizeRenderer(slot.renderer);
  const side = normalizePlacement(slot.side, renderer === 'features' ? 'overlay' : 'below');
  if (renderer === 'features' && side === 'overlay') return 'overlay';
  return side === 'above' ? 'above' : 'below';
};

export const inferLinearTrackAxisIndexFromSlots = (slots) => {
  const normalizedSlots = Array.isArray(slots) ? slots : [];
  const featureIndex = normalizedSlots.findIndex((slot) => slot?.renderer === 'features');
  if (featureIndex >= 0) {
    const featurePlacement = effectiveLinearSlotPlacement(normalizedSlots[featureIndex]);
    if (featurePlacement === 'above') return featureIndex + 1;
    if (featurePlacement === 'overlay') return featureIndex;
    return featureIndex;
  }
  return normalizedSlots.filter((slot) => effectiveLinearSlotPlacement(slot) === 'above').length;
};

export const resolveLinearTrackAxisIndex = (slots, axisIndex = null) => {
  const normalizedSlots = Array.isArray(slots) ? slots : [];
  const clamped = clampLinearTrackAxisIndex(axisIndex, normalizedSlots.length);
  return clamped === null ? inferLinearTrackAxisIndexFromSlots(normalizedSlots) : clamped;
};

export const syncLinearSlotPlacementFromSide = (slot, placement) => {
  if (!slot) return;
  const renderer = normalizeRenderer(slot.renderer);
  const normalizedPlacement = normalizePlacement(placement, renderer === 'features' ? 'overlay' : 'below');
  slot.side = renderer === 'features' || normalizedPlacement !== 'overlay'
    ? normalizedPlacement
    : 'below';
};

export const syncLinearSlotsFromAxisIndex = (slots, axisIndex) => {
  const normalizedSlots = Array.isArray(slots) ? slots : [];
  const resolvedAxis = resolveLinearTrackAxisIndex(normalizedSlots, axisIndex);
  normalizedSlots.forEach((slot, index) => {
    if (!slot) return;
    const isOnAxisFeature = (
      index === resolvedAxis &&
      normalizeRenderer(slot.renderer) === 'features' &&
      effectiveLinearSlotPlacement(slot) === 'overlay'
    );
    syncLinearSlotPlacementFromSide(
      slot,
      isOnAxisFeature ? 'overlay' : (index < resolvedAxis ? 'above' : 'below')
    );
  });
  return resolvedAxis;
};

export const enforceSingleLinearOnAxisSlot = (slots, axisIndex) => {
  const normalizedSlots = Array.isArray(slots) ? slots : [];
  const onAxisIndices = normalizedSlots
    .map((slot, index) => (
      normalizeRenderer(slot?.renderer) === 'features' &&
      effectiveLinearSlotPlacement(slot) === 'overlay'
        ? index
        : null
    ))
    .filter((index) => Number.isInteger(index));
  if (onAxisIndices.length === 0) {
    return resolveLinearTrackAxisIndex(normalizedSlots, axisIndex);
  }

  const clampedAxis = clampLinearTrackAxisIndex(axisIndex, normalizedSlots.length);
  const keepIndex = onAxisIndices.includes(clampedAxis) ? clampedAxis : onAxisIndices[0];
  onAxisIndices.forEach((index) => {
    if (index === keepIndex) return;
    syncLinearSlotPlacementFromSide(normalizedSlots[index], index < keepIndex ? 'above' : 'below');
  });
  return keepIndex;
};

export const buildLinearTrackSlotSpec = (slot, { includeEnabled = false, includeSide = true } = {}) => {
  const normalized = normalizeLinearTrackSlots([slot])[0];
  if (!normalized) return '';
  const parts = [];
  if (includeEnabled && normalized.enabled === false) parts.push('enabled=false');
  if (normalized.side && (includeSide || normalized.side === 'overlay')) {
    parts.push(`side=${normalized.side}`);
  }
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

export const applyLinearTrackOrderPlacements = (slots, axisIndex = null, nt = 'GC', trackLayout = 'middle') => {
  const normalized = normalizeLinearTrackSlots(slots, nt, trackLayout);
  const resolvedAxis = syncLinearSlotsFromAxisIndex(normalized, axisIndex);
  enforceSingleLinearOnAxisSlot(normalized, resolvedAxis);
  return normalized;
};

const paramsMatchAllowedKeys = (params, allowedKeys) => {
  const allowed = new Set(allowedKeys);
  return Object.entries(cloneParams(params)).every(([key, value]) => (
    allowed.has(String(key)) || normalizeOptionalText(value) === null
  ));
};

const hasBlankLinearSlotGeometry = (slot) => (
  normalizeOptionalText(slot?.height) === null &&
  normalizeOptionalText(slot?.spacing) === null &&
  Number(slot?.z || 0) === 0
);

const isDefaultManagedLinearSlot = (slot, renderer = null) => {
  if (!slot || typeof slot !== 'object' || Array.isArray(slot)) return false;
  const normalizedRenderer = normalizeRenderer(slot.renderer);
  if (renderer !== null && normalizedRenderer !== normalizeRenderer(renderer)) return false;
  if (!hasBlankLinearSlotGeometry(slot)) return false;
  const id = String(slot.id || '').trim();
  const params = cloneParams(slot.params);
  if (normalizeOptionalText(params.legend_label) !== null) return false;
  if (params.managed === DEFAULT_LINEAR_SLOT_MANAGER) return true;

  if (normalizedRenderer === 'features') {
    return id === 'features' && paramsMatchAllowedKeys(params, []);
  }
  if (normalizedRenderer === 'dinucleotide_content') {
    return id === 'gc_content' && paramsMatchAllowedKeys(params, ['nt', 'dinucleotide']);
  }
  if (normalizedRenderer === 'dinucleotide_skew') {
    return id === 'gc_skew' && paramsMatchAllowedKeys(params, ['nt', 'dinucleotide']);
  }
  if (normalizedRenderer === 'depth') {
    return /^depth(?:_\d+)?$/.test(id) && paramsMatchAllowedKeys(params, ['track_index']);
  }
  return false;
};

export const createLinearTrackSlotEditor = ({ state }) => {
  const { adv, form } = state;

  const linearTrackRenderers = SUPPORTED_RENDERERS.slice();
  const linearTrackRendererLabel = (renderer) => RENDERER_LABELS[normalizeRenderer(renderer)] || String(renderer || '');

  const axisIndexForCurrentLinearSlots = (slots) => {
    const current = clampLinearTrackAxisIndex(adv.linear_track_slots_axis_index, slots.length);
    if (current !== null) {
      adv.linear_track_slots_axis_index = current;
      return current;
    }
    const inferred = inferLinearTrackAxisIndexFromSlots(slots);
    adv.linear_track_slots_axis_index = inferred;
    return inferred;
  };

  const normalizedSlotsForCurrentState = () => applyLinearTrackOrderPlacements(
    adv.linear_track_slots,
    adv.linear_track_slots_axis_index,
    adv.nt,
    form.linear_track_layout
  );

  const normalizeCurrentSlots = () => {
    const normalized = normalizeLinearTrackSlots(adv.linear_track_slots, adv.nt, form.linear_track_layout);
    const axis = axisIndexForCurrentLinearSlots(normalized);
    syncLinearSlotsFromAxisIndex(normalized, axis);
    adv.linear_track_slots_axis_index = enforceSingleLinearOnAxisSlot(
      normalized,
      axis
    );
    adv.linear_track_slots.splice(0, adv.linear_track_slots.length, ...normalized);
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
    const normalized = normalizeLinearTrackSlots(slots, adv.nt, form.linear_track_layout);
    adv.linear_track_slots_axis_index = inferLinearTrackAxisIndexFromSlots(normalized);
    syncLinearSlotsFromAxisIndex(normalized, adv.linear_track_slots_axis_index);
    adv.linear_track_slots_axis_index = enforceSingleLinearOnAxisSlot(
      normalized,
      adv.linear_track_slots_axis_index
    );
    adv.linear_track_slots.splice(0, adv.linear_track_slots.length, ...normalized);
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
    normalizeCurrentSlots();
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
    const axis = axisIndexForCurrentLinearSlots(adv.linear_track_slots);
    if (idx < axis) adv.linear_track_slots_axis_index = axis + 1;
    normalizeCurrentSlots();
  };

  const removeLinearTrackSlot = (index) => {
    normalizeCurrentSlots();
    const idx = Number(index);
    if (!Number.isInteger(idx) || idx < 0 || idx >= adv.linear_track_slots.length) return;
    const axis = axisIndexForCurrentLinearSlots(adv.linear_track_slots);
    adv.linear_track_slots.splice(idx, 1);
    adv.linear_track_slots_axis_index = idx < axis
      ? Math.max(0, axis - 1)
      : Math.min(axis, adv.linear_track_slots.length);
    normalizeCurrentSlots();
  };

  const wouldLinearTrackSlotMoveCrossAxis = (fromIndex, toIndex) => {
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

    const axis = axisIndexForCurrentLinearSlots(normalized);
    const movedPlacement = effectiveLinearSlotPlacement(normalized[from]);
    const targetPlacement = effectiveLinearSlotPlacement(normalized[to]);
    if (movedPlacement === 'overlay' || targetPlacement === 'overlay') return true;
    return (from < axis) !== (to < axis);
  };

  const canMoveLinearTrackSlot = (index, direction) => {
    const idx = Number(index);
    const step = Number(direction);
    if (!Number.isInteger(idx) || !Number.isInteger(step) || step === 0) return false;
    const target = idx + Math.sign(step);
    return !wouldLinearTrackSlotMoveCrossAxis(idx, target);
  };

  const moveLinearTrackSlot = (fromIndex, toIndex) => {
    if (wouldLinearTrackSlotMoveCrossAxis(fromIndex, toIndex)) return;
    normalizeCurrentSlots();
    const from = Number(fromIndex);
    const to = Number(toIndex);
    if (
      !Number.isInteger(from) ||
      !Number.isInteger(to) ||
      from < 0 ||
      to < 0 ||
      from >= adv.linear_track_slots.length ||
      to >= adv.linear_track_slots.length ||
      from === to
    ) {
      return;
    }
    const [slot] = adv.linear_track_slots.splice(from, 1);
    adv.linear_track_slots.splice(to, 0, slot);
    normalizeCurrentSlots();
  };

  const canMoveLinearTrackSlotAbove = (index) => {
    const idx = Number(index);
    const normalized = normalizedSlotsForCurrentState();
    if (!Number.isInteger(idx) || idx < 0 || idx >= normalized.length) return false;
    if (effectiveLinearSlotPlacement(normalized[idx]) === 'overlay') return true;
    return idx >= axisIndexForCurrentLinearSlots(normalized);
  };

  const canMoveLinearTrackSlotBelow = (index) => {
    const idx = Number(index);
    const normalized = normalizedSlotsForCurrentState();
    if (!Number.isInteger(idx) || idx < 0 || idx >= normalized.length) return false;
    if (effectiveLinearSlotPlacement(normalized[idx]) === 'overlay') return true;
    return idx < axisIndexForCurrentLinearSlots(normalized);
  };

  const canMoveLinearTrackSlotToAxis = (index) => {
    const idx = Number(index);
    const normalized = normalizedSlotsForCurrentState();
    if (!Number.isInteger(idx) || idx < 0 || idx >= normalized.length) return false;
    const slot = normalized[idx];
    return slot?.renderer === 'features' && effectiveLinearSlotPlacement(slot) !== 'overlay';
  };

  const moveLinearTrackSlotToPlacement = (index, placement) => {
    const idx = Number(index);
    if (!Number.isInteger(idx) || idx < 0 || idx >= adv.linear_track_slots.length) return;
    normalizeCurrentSlots();
    if (idx >= adv.linear_track_slots.length) return;
    const targetPlacement = normalizePlacement(placement);
    const movingSlot = adv.linear_track_slots[idx];
    if (targetPlacement === 'overlay' && normalizeRenderer(movingSlot?.renderer) !== 'features') return;

    if (targetPlacement === 'overlay') {
      const movedPreviousPlacement = effectiveLinearSlotPlacement(movingSlot);
      const existingAxisIndex = adv.linear_track_slots.findIndex((slot, slotIndex) => (
        slotIndex !== idx &&
        normalizeRenderer(slot?.renderer) === 'features' &&
        effectiveLinearSlotPlacement(slot) === 'overlay'
      ));
      if (existingAxisIndex >= 0) {
        const existingAxisSlot = adv.linear_track_slots[existingAxisIndex];
        const demotedPlacement = movedPreviousPlacement === 'overlay' ? 'below' : movedPreviousPlacement;
        syncLinearSlotPlacementFromSide(existingAxisSlot, demotedPlacement);
        syncLinearSlotPlacementFromSide(movingSlot, 'overlay');
        adv.linear_track_slots[existingAxisIndex] = movingSlot;
        adv.linear_track_slots[idx] = existingAxisSlot;
        adv.linear_track_slots_axis_index = existingAxisIndex;
        normalizeCurrentSlots();
        return;
      }
    }

    let axis = axisIndexForCurrentLinearSlots(adv.linear_track_slots);
    const onAxisIndex = adv.linear_track_slots.findIndex((slot) => (
      normalizeRenderer(slot?.renderer) === 'features' &&
      effectiveLinearSlotPlacement(slot) === 'overlay'
    ));
    const [slot] = adv.linear_track_slots.splice(idx, 1);
    if (!slot) return;
    if (idx < axis) axis -= 1;

    syncLinearSlotPlacementFromSide(slot, targetPlacement);
    if (targetPlacement === 'above') {
      adv.linear_track_slots.splice(axis, 0, slot);
      adv.linear_track_slots_axis_index = axis + 1;
    } else if (targetPlacement === 'overlay') {
      adv.linear_track_slots.splice(axis, 0, slot);
      adv.linear_track_slots_axis_index = axis;
    } else {
      const adjustedOnAxisIndex = onAxisIndex >= 0
        ? (idx === onAxisIndex ? -1 : (idx < onAxisIndex ? onAxisIndex - 1 : onAxisIndex))
        : -1;
      const insertIndex = adjustedOnAxisIndex >= 0 ? adjustedOnAxisIndex + 1 : axis;
      adv.linear_track_slots.splice(insertIndex, 0, slot);
      adv.linear_track_slots_axis_index = adjustedOnAxisIndex >= 0 ? adjustedOnAxisIndex : axis;
    }
    normalizeCurrentSlots();
  };

  const moveLinearTrackSlotAbove = (index) => {
    if (!canMoveLinearTrackSlotAbove(index)) return;
    moveLinearTrackSlotToPlacement(index, 'above');
  };

  const moveLinearTrackSlotBelow = (index) => {
    if (!canMoveLinearTrackSlotBelow(index)) return;
    moveLinearTrackSlotToPlacement(index, 'below');
  };

  const moveLinearTrackSlotToAxis = (index) => {
    if (!canMoveLinearTrackSlotToAxis(index)) return;
    moveLinearTrackSlotToPlacement(index, 'overlay');
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

  const updateLinearTrackSlotPlacement = (slot, placement) => {
    if (!slot) return;
    const index = adv.linear_track_slots.findIndex((candidate) => candidate === slot);
    if (index >= 0) {
      moveLinearTrackSlotToPlacement(index, placement);
      return;
    }
    syncLinearSlotPlacementFromSide(slot, placement);
    normalizeCurrentSlots();
  };

  const findDefaultFeatureIndex = () => {
    const slots = Array.isArray(adv.linear_track_slots) ? adv.linear_track_slots : [];
    const preferred = slots.findIndex((slot) => slot?.renderer === 'features' && String(slot?.id || '') === 'features');
    if (preferred >= 0) return preferred;
    return slots.findIndex((slot) => slot?.renderer === 'features');
  };

  const ensureFeatureSlot = () => {
    if (findDefaultFeatureIndex() >= 0) return;
    adv.linear_track_slots.unshift(defaultSlot('features', {
      id: 'features',
      side: sideForLinearTrackLayout(form.linear_track_layout)
    }));
    adv.linear_track_slots_axis_index = null;
  };

  const removeDefaultManagedSlots = (renderer) => {
    const next = adv.linear_track_slots.filter((slot) => !isDefaultManagedLinearSlot(slot, renderer));
    if (next.length !== adv.linear_track_slots.length) {
      adv.linear_track_slots.splice(0, adv.linear_track_slots.length, ...next);
      adv.linear_track_slots_axis_index = null;
    }
  };

  const ensureDefaultNumericSlot = (renderer, id, params = {}) => {
    const normalizedRenderer = normalizeRenderer(renderer);
    const entries = adv.linear_track_slots
      .map((slot, index) => ({ slot, index }))
      .filter((entry) => isDefaultManagedLinearSlot(entry.slot, normalizedRenderer));
    const primary = entries[0]?.slot || null;
    if (primary) {
      primary.id = id;
      primary.renderer = normalizedRenderer;
      primary.enabled = true;
      primary.params = cloneParams(params);
      entries.slice(1).reverse().forEach((entry) => {
        adv.linear_track_slots.splice(entry.index, 1);
      });
      return;
    }
    adv.linear_track_slots.push(defaultSlot(normalizedRenderer, {
      id,
      side: 'below',
      params
    }));
  };

  const ensureLinearTrackDepthSlots = () => {
    const desiredCount = linearDepthTrackCountForState(state);
    if (!Boolean(form.show_depth) || desiredCount <= 0) {
      removeDefaultManagedSlots('depth');
      normalizeCurrentSlots();
      return;
    }

    normalizeCurrentSlots();
    const slots = adv.linear_track_slots;
    const managedDepthEntries = slots
      .map((slot, index) => ({ slot, index }))
      .filter((entry) => isDefaultManagedLinearSlot(entry.slot, 'depth'));
    const claimed = new Set();
    const removeIndexes = [];
    managedDepthEntries.forEach((entry) => {
      const trackIndex = normalizeTrackIndex(entry.slot?.params?.track_index);
      if (trackIndex !== null && trackIndex < desiredCount && !claimed.has(trackIndex)) {
        entry.slot.enabled = true;
        entry.slot.params = { track_index: trackIndex };
        claimed.add(trackIndex);
      } else {
        removeIndexes.push(entry.index);
      }
    });
    removeIndexes.sort((a, b) => b - a).forEach((index) => {
      slots.splice(index, 1);
    });

    const existingIds = new Set(slots.map((slot) => String(slot?.id || '').trim()).filter(Boolean));
    for (let trackIndex = 0; trackIndex < desiredCount; trackIndex += 1) {
      if (claimed.has(trackIndex)) continue;
      const preferredId = desiredCount === 1 ? 'depth' : `depth_${trackIndex + 1}`;
      let id = preferredId;
      let suffix = 2;
      while (existingIds.has(id)) {
        id = `${preferredId}_${suffix}`;
        suffix += 1;
      }
      existingIds.add(id);
      slots.push(defaultSlot('depth', {
        id,
        side: 'below',
        params: { track_index: trackIndex }
      }));
      claimed.add(trackIndex);
    }
    normalizeCurrentSlots();
  };

  const syncLinearNumericSlotsFromSimpleControls = () => {
    normalizeCurrentSlots();
    if (Boolean(form.show_gc)) {
      ensureDefaultNumericSlot('dinucleotide_content', 'gc_content', { nt: normalizeNt(adv.nt) });
    } else {
      removeDefaultManagedSlots('dinucleotide_content');
    }
    if (Boolean(form.show_skew)) {
      ensureDefaultNumericSlot('dinucleotide_skew', 'gc_skew', { nt: normalizeNt(adv.nt) });
    } else {
      removeDefaultManagedSlots('dinucleotide_skew');
    }
    ensureLinearTrackDepthSlots();
    normalizeCurrentSlots();
  };

  const applyLinearTrackLayoutPreset = (trackLayout = form.linear_track_layout) => {
    ensureFeatureSlot();
    normalizeCurrentSlots();
    const featureIndex = findDefaultFeatureIndex();
    if (featureIndex < 0) return;
    moveLinearTrackSlotToPlacement(featureIndex, sideForLinearTrackLayout(trackLayout));
  };

  const reconcileLinearTrackSlotsFromSimpleControls = () => {
    ensureFeatureSlot();
    applyLinearTrackLayoutPreset(form.linear_track_layout);
    syncLinearNumericSlotsFromSimpleControls();
    normalizeCurrentSlots();
  };

  const linearTrackStackEntries = () => {
    const slots = Array.isArray(adv.linear_track_slots) ? adv.linear_track_slots : [];
    const axisIndex = axisIndexForCurrentLinearSlots(slots);
    const entries = [];
    let axisRendered = false;
    slots.forEach((slot, index) => {
      const onAxis = (
        index === axisIndex &&
        normalizeRenderer(slot?.renderer) === 'features' &&
        effectiveLinearSlotPlacement(slot) === 'overlay'
      );
      if (index === axisIndex && !onAxis) {
        entries.push({ kind: STACK_ENTRY_AXIS, key: 'axis' });
        axisRendered = true;
      }
      if (onAxis && !axisRendered) {
        axisRendered = true;
      }
      entries.push({ kind: STACK_ENTRY_SLOT, slot, index, onAxis });
    });
    if (!axisRendered || axisIndex >= slots.length) entries.push({ kind: STACK_ENTRY_AXIS, key: 'axis' });
    return entries;
  };

  const linearTrackSlotPlacementLabel = (slot) => {
    const placement = effectiveLinearSlotPlacement(slot);
    if (placement === 'overlay') return 'On Axis';
    return placement === 'above' ? 'Above Axis' : 'Below Axis';
  };

  const linearTrackSlotUsesPresetGeometry = (slot) => {
    if (!slot || typeof slot !== 'object') return false;
    const renderer = RENDERER_ALIASES[String(slot.renderer || '').trim().toLowerCase()] || String(slot.renderer || '').trim().toLowerCase();
    if (!SUPPORTED_RENDERERS.includes(renderer)) return false;
    return hasBlankLinearSlotGeometry(slot);
  };

  const parsePositivePxNumber = (value) => {
    const text = normalizeOptionalText(value);
    if (text === null) return null;
    const withoutUnit = text.endsWith('px') ? text.slice(0, -2) : text;
    const numeric = Number(withoutUnit);
    return Number.isFinite(numeric) && numeric > 0 ? numeric : null;
  };

  const ensureDepthTrackConfigForSlotIndex = (trackIndex) => {
    const idx = Math.max(0, Number(trackIndex) || 0);
    if (!Array.isArray(adv.depth_tracks)) adv.depth_tracks = [];
    while (adv.depth_tracks.length <= idx) {
      const nextIndex = adv.depth_tracks.length;
      adv.depth_tracks.push({
        label: nextIndex === 0 ? 'Depth' : `Depth ${nextIndex + 1}`,
        color: nextIndex === 0 ? String(adv.depth_color || '#4A90E2') : '',
        height: null,
        large_tick_interval: null,
        small_tick_interval: null,
        tick_font_size: null
      });
    }
    if (!adv.depth_tracks[idx] || typeof adv.depth_tracks[idx] !== 'object' || Array.isArray(adv.depth_tracks[idx])) {
      adv.depth_tracks[idx] = {
        label: idx === 0 ? 'Depth' : `Depth ${idx + 1}`,
        color: idx === 0 ? String(adv.depth_color || '#4A90E2') : '',
        height: null,
        large_tick_interval: null,
        small_tick_interval: null,
        tick_font_size: null
      };
    }
    return adv.depth_tracks[idx];
  };

  const linearDepthTrackIndexForSlot = (slot) => {
    if (!slot || normalizeRenderer(slot.renderer) !== 'depth') return null;
    const params = cloneParams(slot.params);
    return normalizeTrackIndex(params.track_index) ?? 0;
  };

  const heightTextFromDepthTrackConfig = (trackIndex) => {
    const config = ensureDepthTrackConfigForSlotIndex(trackIndex);
    const height = parsePositivePxNumber(config.height);
    return height === null ? '' : String(height);
  };

  const syncLinearDepthSlotHeightsFromDepthTracks = (trackIndex = null) => {
    const slots = Array.isArray(adv.linear_track_slots) ? adv.linear_track_slots : [];
    slots.forEach((slot) => {
      const slotTrackIndex = linearDepthTrackIndexForSlot(slot);
      if (slotTrackIndex === null) return;
      if (trackIndex !== null && Number(trackIndex) !== slotTrackIndex) return;
      slot.height = heightTextFromDepthTrackConfig(slotTrackIndex);
    });
  };

  const linearTrackSlotHeightValue = (slot) => {
    if (normalizeRenderer(slot?.renderer) !== 'depth') {
      return String(slot?.height || '');
    }
    const trackIndex = linearDepthTrackIndexForSlot(slot) ?? 0;
    const heightText = heightTextFromDepthTrackConfig(trackIndex);
    return heightText || String(slot?.height || '');
  };

  const setLinearTrackSlotHeight = (slot, value) => {
    if (!slot) return;
    const text = String(value ?? '').trim();
    slot.height = text;
    if (normalizeRenderer(slot.renderer) !== 'depth') return;
    const trackIndex = linearDepthTrackIndexForSlot(slot) ?? 0;
    const config = ensureDepthTrackConfigForSlotIndex(trackIndex);
    config.height = parsePositivePxNumber(text);
  };

  return {
    linearTrackRenderers,
    linearTrackRendererLabel,
    normalizeLinearTrackSlots: normalizeCurrentSlots,
    resetLinearTrackSlotsFromSimpleControls,
    reconcileLinearTrackSlotsFromSimpleControls,
    ensureLinearTrackDepthSlots,
    syncLinearNumericSlotsFromSimpleControls,
    applyLinearTrackLayoutPreset,
    setLinearTrackSlotsEnabled,
    addLinearTrackSlot,
    duplicateLinearTrackSlot,
    removeLinearTrackSlot,
    moveLinearTrackSlot,
    canMoveLinearTrackSlot,
    moveLinearTrackSlotAbove,
    moveLinearTrackSlotBelow,
    moveLinearTrackSlotToAxis,
    moveLinearTrackSlotToPlacement,
    canMoveLinearTrackSlotAbove,
    canMoveLinearTrackSlotBelow,
    canMoveLinearTrackSlotToAxis,
    updateLinearTrackSlotRenderer,
    updateLinearTrackSlotPlacement,
    linearTrackSlotHeightValue,
    setLinearTrackSlotHeight,
    syncLinearDepthSlotHeightsFromDepthTracks,
    linearTrackSlots: () => {
      return Array.isArray(adv.linear_track_slots) ? adv.linear_track_slots : [];
    },
    linearTrackStackEntries,
    linearTrackSlotCliSpec: (slot) => buildLinearTrackSlotSpec(slot),
    linearTrackSlotDisplayLabel: (slot) => linearTrackRendererLabel(slot?.renderer),
    linearTrackSlotDisplayMeta: (slot) => buildLinearTrackSlotSpec(slot),
    linearTrackSlotPlacementLabel,
    linearTrackSlotUsesPresetGeometry
  };
};

export { SUPPORTED_RENDERERS as LINEAR_TRACK_RENDERERS };
