// @ts-check
/** @import { DrawingState } from '../state.js' */
/** @import { ChangeTrackLayout } from './track-slot-edits.js' */
import {
  activeDepthTrackIndices,
  depthTrackMatrixWidth,
  reconcileManagedDepthSlots
} from '../services/depth-track-state.js';
import { resolveTrackSlotSkewColorValue } from './track-slot-colors.js';
import {
  findTrackSlotGeometry,
  formatPxAuto,
  isManualSlotValue,
  normalizeOptionalText
} from '../services/track-slot-display.js';
import { featureSlotEdits } from './track-slot-edits.js';
import { paramsKeptOnRendererChange, parseOptionalPixel, validateCustomTrackPlan } from '../services/track-slot-validation.js';
import { visibleFeatureUnderlaysForState } from '../utils/feature-rendering.js';

import {
  applyLinearTrackOrderPlacements,
  clampLinearTrackAxisIndex,
  cloneParams,
  createDefaultLinearTrackSlots,
  DEFAULT_SLOT_IDS,
  defaultSlot,
  effectiveLinearSlotPlacement,
  enforceSingleLinearOnAxisSlot,
  inferLinearTrackAxisIndexFromSlots,
  normalizeColorParam,
  normalizeLinearTrackSlots,
  normalizeNt,
  normalizePlacement,
  normalizeRenderer,
  normalizeTrackIndex,
  RENDERER_ALIASES,
  resolveLinearTrackAxisIndex,
  LINEAR_TRACK_RENDERERS as SUPPORTED_RENDERERS,
  syncLinearSlotPlacementFromSide,
  syncLinearSlotsFromAxisIndex
} from '../services/linear-track-slot-model.js';

const UI_RENDERERS = SUPPORTED_RENDERERS.slice();

const RENDERER_LABELS = {
  features: 'Features',
  dinucleotide_content: 'Dinucleotide content',
  dinucleotide_skew: 'Dinucleotide skew',
  depth: 'Depth',
  annotations: 'Annotations',
  spacer: 'Spacer'
};

const STACK_ENTRY_AXIS = 'axis';
const STACK_ENTRY_SLOT = 'slot';
const DEFAULT_LINEAR_SLOT_MANAGER = 'linear-default';
const ESTIMATED_LINEAR_GC_HEIGHT_PX = 20;
const ESTIMATED_LINEAR_DEPTH_HEIGHT_PX = 10;
const ESTIMATED_LINEAR_DEPTH_SPACING_PX = 8;
const ANNOTATION_MARK_OPTIONS = Object.freeze([
  'line',
  'bracket',
  'band',
  'highlight'
]);

export const linearAvailableDepthTrackCountForState = (state) => {
  const seqs = Array.isArray(state?.linearSeqs) ? state.linearSeqs : [];
  return depthTrackMatrixWidth(seqs.map((seq) => seq?.depth));
};

const linearSourcedDepthTrackIndexesForState = (state) => (
  activeDepthTrackIndices((Array.isArray(state?.linearSeqs) ? state.linearSeqs : []).map((seq) => seq?.depth))
);

export const hasEnabledLinearTrackRenderer = (slots, renderer) => {
  const normalizedRenderer = normalizeRenderer(renderer);
  return normalizeLinearTrackSlots(slots).some(
    (slot) => slot.enabled !== false && slot.renderer === normalizedRenderer
  );
};

export const linearTrackAxisIndexForEnabledSlots = (slots, axisIndex = null) => {
  const normalizedSlots = Array.isArray(slots) ? slots : [];
  const resolvedAxis = resolveLinearTrackAxisIndex(normalizedSlots, axisIndex);
  return normalizedSlots
    .slice(0, resolvedAxis)
    .filter((slot) => slot?.enabled !== false)
    .length;
};

export const buildLinearTrackSlotSpec = (slot, { includeEnabled = false, includeSide = true } = {}) => {
  const normalized = normalizeLinearTrackSlots([slot])[0];
  if (!normalized) return '';
  const parts = [];
  if (includeEnabled && normalized.enabled === false) parts.push('enabled=false');
  if (normalized.side && (includeSide || normalized.side === 'overlay')) {
    parts.push(`side=${normalized.side}`);
  }
  const height = parseOptionalPixel(normalized.height, `Linear track '${normalized.id}' height`, { allowZero: false });
  const spacing = parseOptionalPixel(normalized.spacing, `Linear track '${normalized.id}' spacing`, { allowZero: true });
  if (height !== null) parts.push(`h=${height}px`);
  if (spacing !== null) parts.push(`spacing=${spacing}px`);
  if (Number(normalized.z) !== 0) parts.push(`z=${Number(normalized.z)}`);
  const params = cloneParams(normalized.params);
  if (normalized.renderer === 'depth') {
    const trackIndex = normalizeTrackIndex(params.track_index);
    parts.push(`track_index=${trackIndex === null ? 0 : trackIndex}`);
  }
  if (normalized.renderer === 'dinucleotide_content' || normalized.renderer === 'dinucleotide_skew') {
    parts.push(`nt=${normalizeNt(params.nt)}`);
    if (normalized.renderer === 'dinucleotide_skew') {
      if (normalizeOptionalText(params.positive_color)) {
        parts.push(`positive_color=${String(params.positive_color).trim()}`);
      }
      if (normalizeOptionalText(params.negative_color)) {
        parts.push(`negative_color=${String(params.negative_color).trim()}`);
      }
    }
  }
  if (normalized.renderer === 'annotations') {
    if (Array.isArray(params.marks) && params.marks.length > 0) {
      parts.push(`marks=${params.marks.join('|')}`);
    }
    ['set_id', 'lane_gap_px', 'padding_px', 'overflow', 'anchor_slot', 'layer'].forEach((key) => {
      if (normalizeOptionalText(params[key])) parts.push(`${key}=${String(params[key]).trim()}`);
    });
    parts.push(`show_labels=${params.show_labels === false ? 'false' : 'true'}`);
    if (params.cover_anchor === true) parts.push('cover_anchor=true');
  }
  if (normalizeOptionalText(params.legend_label)) {
    parts.push(`legend_label=${String(params.legend_label).trim()}`);
  }
  const suffix = parts.length > 0 ? `@${parts.join(',')}` : '';
  return `${normalized.id}:${normalized.renderer}${suffix}`;
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

/** @param {string | null} [renderer] */
const isDefaultManagedLinearSlot = (slot, renderer = null) => {
  if (!slot || typeof slot !== 'object' || Array.isArray(slot)) return false;
  const normalizedRenderer = normalizeRenderer(slot.renderer);
  if (renderer !== null && normalizedRenderer !== normalizeRenderer(renderer)) return false;
  if (!hasBlankLinearSlotGeometry(slot)) return false;
  const id = String(slot.id || '').trim();
  const params = cloneParams(slot.params);
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
    return /^depth(?:_\d+)?$/.test(id) && paramsMatchAllowedKeys(params, ['track_index', 'legend_label']);
  }
  return false;
};

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

/**
 * @typedef {object} LinearTrackSlotEditorOptions
 * @property {Record<string, any>} state the Web state; its shape belongs to `state.js`
 * @property {ChangeTrackLayout} [changeTrackLayout] R10 port (default: apply directly)
 */

// `changeTrackLayout` is the feature placement owner's transition, injected as
// a port (R10, Q3, R13): every stack edit that can change the feature slot
// runs through it.
/**
 * @param {LinearTrackSlotEditorOptions} options
 */
export const createLinearTrackSlotEditor = ({ state, changeTrackLayout = (apply) => apply() }) => {
  const editorKeys = new WeakMap();
  let nextEditorKey = 1;
  const linearTrackSlotEditorKey = (slot) => {
    if (!slot || typeof slot !== 'object') return 'linear-slot-invalid';
    let key = editorKeys.get(slot);
    if (!key) {
      key = `linear-editor-slot-${nextEditorKey}`;
      nextEditorKey += 1;
      editorKeys.set(slot, key);
    }
    return key;
  };

  const linearTrackRenderers = UI_RENDERERS.slice();
  const linearTrackRendererLabel = (renderer) => RENDERER_LABELS[normalizeRenderer(renderer)] || String(renderer || '');
  /** @param {DrawingState} drawing */
  const annotationSetIds = (drawing) => (
    (Array.isArray(drawing.annotationSets) ? drawing.annotationSets : [])
      .map((set) => String(set?.id || '').trim())
      .filter(Boolean)
  );

  /** @param {DrawingState} drawing */
  const linearTrackValidationPlan = (drawing) => validateCustomTrackPlan({
    mode: 'linear',
    slots: drawing.adv.linear_track_slots,
    axisIndex: drawing.adv.linear_track_slots_axis_index,
    trackType: drawing.form.linear_track_layout,
    depthTrackCount: linearAvailableDepthTrackCountForState(state),
    depthSourcedTrackIndexes: linearSourcedDepthTrackIndexesForState(state),
    annotationSetIds: annotationSetIds(drawing),
    visibleFeatureUnderlays: visibleFeatureUnderlaysForState(state),
    conservationSeries: []
  });

  const linearTrackSlotIssue = (slot, index = null) => {
    const drawing = state.drawings.linear;
    const resolvedIndex = Number.isInteger(Number(index))
      ? Number(index)
      : drawing.adv.linear_track_slots.findIndex((candidate) => candidate === slot);
    if (resolvedIndex < 0) return '';
    return (linearTrackValidationPlan(drawing).rowIssues.get(resolvedIndex) || [])
      .map((issue) => issue.message)
      .join(' ');
  };

  const linearTrackGlobalIssues = () => {
    const drawing = state.drawings.linear;
    return (
      linearTrackValidationPlan(drawing).globalIssues.map((issue) => issue.message)
    );
  };

  const linearAnnotationAnchorOptions = (slot = null) => (
    state.drawings.linear.adv.linear_track_slots
      .filter((candidate) => (
        candidate &&
        candidate !== slot &&
        candidate.enabled !== false &&
        !['annotations', 'spacer'].includes(normalizeRenderer(candidate.renderer)) &&
        String(candidate.id || '').trim()
      ))
      .map((candidate) => ({
        id: String(candidate.id).trim(),
        label: `${String(candidate.id).trim()} · ${linearTrackRendererLabel(candidate.renderer)}`
      }))
  );

  const linearAnnotationAnchorIsKnown = (slot) => {
    const anchor = String(slot?.params?.anchor_slot || '').trim();
    return !anchor || linearAnnotationAnchorOptions(slot).some((option) => option.id === anchor);
  };

  const bindOnlyLinearAnnotationAnchor = (slot) => {
    if (!slot || normalizeRenderer(slot.renderer) !== 'annotations') return;
    slot.params = cloneParams(slot.params);
    const options = linearAnnotationAnchorOptions(slot);
    const current = String(slot.params.anchor_slot || '').trim();
    if (options.some((option) => option.id === current)) return;
    if (options.length === 1) slot.params.anchor_slot = options[0].id;
    else delete slot.params.anchor_slot;
  };

  const canAddLinearTrackRenderer = (renderer) => {
    const drawing = state.drawings.linear;
    const normalizedRenderer = normalizeRenderer(renderer, 'spacer');
    if (normalizedRenderer === 'annotations') return annotationSetIds(drawing).length > 0;
    if (normalizedRenderer === 'depth') {
      return linearAvailableDepthTrackCountForState(state) > 0;
    }
    if (normalizedRenderer === 'features') {
      return !drawing.adv.linear_track_slots.some((slot) => (
        slot?.enabled !== false && normalizeRenderer(slot?.renderer) === 'features'
      ));
    }
    return true;
  };

  // Duplicate shares the availability of the other row controls (GX-01).
  const canDuplicateLinearTrackSlot = (slot) => {
    const drawing = state.drawings.linear;
    if (!slot || state.sessionOperationAvailability?.()) return false;
    if (slot.enabled === false) return true;
    const renderer = normalizeRenderer(slot.renderer);
    if (renderer === 'features') return false;
    if (renderer === 'annotations') return annotationSetIds(drawing).length > 0;
    if (renderer === 'depth') return linearAvailableDepthTrackCountForState(state) > 0;
    return true;
  };

  const linearAnnotationMarkSelected = (slot, mark) => {
    const selected = Array.isArray(slot?.params?.marks)
      ? slot.params.marks.map((value) => String(value).trim().toLowerCase()).filter(Boolean)
      : [];
    return selected.length === 0 || selected.includes(mark);
  };

  const setLinearAnnotationMarkSelected = (slot, mark, checked) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!slot || normalizeRenderer(slot.renderer) !== 'annotations' || !ANNOTATION_MARK_OPTIONS.includes(mark)) return;
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

  const linearAnnotationNumberValue = (slot, field, defaultValue) => {
    const raw = slot?.params?.[field];
    return raw === null || raw === undefined || raw === '' ? defaultValue : raw;
  };

  const setLinearAnnotationNumber = (slot, field, value, defaultValue) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!slot || normalizeRenderer(slot.renderer) !== 'annotations') return;
    slot.params = cloneParams(slot.params);
    if (value === null || value === undefined || value === '') {
      delete slot.params[field];
      return;
    }
    const numeric = Number(value);
    if (Number.isFinite(numeric) && numeric === defaultValue) delete slot.params[field];
    else slot.params[field] = numeric;
  };

  const linearAnnotationCoverAnchor = (slot) => slot?.params?.cover_anchor === true;

  const setLinearAnnotationCoverAnchor = (slot, checked) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!slot || normalizeRenderer(slot.renderer) !== 'annotations') return;
    slot.params = cloneParams(slot.params);
    if (checked) slot.params.cover_anchor = true;
    else delete slot.params.cover_anchor;
  };

  /** @param {DrawingState} drawing */
  const axisIndexForCurrentLinearSlots = (drawing, slots) => {
    const current = clampLinearTrackAxisIndex(drawing.adv.linear_track_slots_axis_index, slots.length);
    if (current !== null) {
      drawing.adv.linear_track_slots_axis_index = current;
      return current;
    }
    const inferred = inferLinearTrackAxisIndexFromSlots(slots);
    drawing.adv.linear_track_slots_axis_index = inferred;
    return inferred;
  };

  /** @param {DrawingState} drawing */
  const normalizedSlotsForCurrentState = (drawing) => applyLinearTrackOrderPlacements(
    drawing.adv.linear_track_slots,
    drawing.adv.linear_track_slots_axis_index,
    drawing.adv.nt,
    drawing.form.linear_track_layout
  );

  const normalizeCurrentSlots = () => {
    const drawing = state.drawings.linear;
    const normalized = normalizeLinearTrackSlots(drawing.adv.linear_track_slots, drawing.adv.nt, drawing.form.linear_track_layout);
    const axis = axisIndexForCurrentLinearSlots(drawing, normalized);
    syncLinearSlotsFromAxisIndex(normalized, axis);
    drawing.adv.linear_track_slots_axis_index = enforceSingleLinearOnAxisSlot(
      normalized,
      axis
    );
    const identityPreserving = normalized.map((slot, index) => (
      replaceObjectContents(drawing.adv.linear_track_slots[index], slot)
    ));
    drawing.adv.linear_track_slots.splice(
      0,
      drawing.adv.linear_track_slots.length,
      ...identityPreserving
    );
  };

  const resetLinearTrackSlotsFromSimpleControls = () => {
    const drawing = state.drawings.linear;
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    // One Depth row per loaded Depth series; Show Depth is not an input (TK-10).
    const depthTrackCount = linearAvailableDepthTrackCountForState(state);
    const slots = createDefaultLinearTrackSlots({
      showDepth: depthTrackCount > 0,
      depthTrackCount,
      showGc: Boolean(drawing.form.show_gc),
      showSkew: Boolean(drawing.form.show_skew),
      nt: drawing.adv.nt,
      trackLayout: drawing.form.linear_track_layout
    });
    const normalized = normalizeLinearTrackSlots(slots, drawing.adv.nt, drawing.form.linear_track_layout);
    drawing.adv.linear_track_slots_axis_index = inferLinearTrackAxisIndexFromSlots(normalized);
    syncLinearSlotsFromAxisIndex(normalized, drawing.adv.linear_track_slots_axis_index);
    drawing.adv.linear_track_slots_axis_index = enforceSingleLinearOnAxisSlot(
      normalized,
      drawing.adv.linear_track_slots_axis_index
    );
    drawing.adv.linear_track_slots.splice(0, drawing.adv.linear_track_slots.length, ...normalized);
  };

  const setLinearTrackSlotsEnabled = (enabled) => {
    const drawing = state.drawings.linear;
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    drawing.adv.linear_track_slots_enabled = Boolean(enabled);
  };

  const addLinearTrackSlot = (renderer = 'spacer') => {
    const drawing = state.drawings.linear;
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const normalizedRenderer = normalizeRenderer(renderer, 'spacer');
    if (!canAddLinearTrackRenderer(normalizedRenderer)) return;
    normalizeCurrentSlots();
    const baseId = DEFAULT_SLOT_IDS[normalizedRenderer] || 'slot';
    let nextId = baseId;
    let suffix = 2;
    const usedIds = new Set(drawing.adv.linear_track_slots.map((slot) => String(slot.id || '')));
    while (usedIds.has(nextId)) {
      nextId = `${baseId}_${suffix}`;
      suffix += 1;
    }
    const slot = defaultSlot(normalizedRenderer, {
      id: nextId,
      side: normalizedRenderer === 'spacer' ? 'below' : undefined,
      height: normalizedRenderer === 'spacer' ? '12px' : '',
      params: normalizedRenderer === 'annotations'
        ? {
            set_id: annotationSetIds(drawing)[0] || '',
            overflow: 'error',
            show_labels: true,
            layer: 'foreground'
          }
        : {}
    });
    if (normalizedRenderer === 'depth') {
      const available = linearAvailableDepthTrackCountForState(state);
      const claimed = new Set(
        drawing.adv.linear_track_slots
          .filter((candidate) => candidate?.enabled !== false && normalizeRenderer(candidate?.renderer) === 'depth')
          .map((candidate) => normalizeTrackIndex(candidate?.params?.track_index))
          .filter((trackIndex) => trackIndex !== null)
      );
      let trackIndex = 0;
      while (trackIndex < available && claimed.has(trackIndex)) trackIndex += 1;
      slot.params.track_index = trackIndex < available ? trackIndex : 0;
    }
    drawing.adv.linear_track_slots.push(slot);
    normalizeCurrentSlots();
  };

  const duplicateLinearTrackSlot = (index) => {
    const drawing = state.drawings.linear;
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    normalizeCurrentSlots();
    const idx = Number(index);
    const source = drawing.adv.linear_track_slots[idx];
    if (!source) return;
    if (!canDuplicateLinearTrackSlot(source)) return;
    const copy = {
      ...source,
      id: `${source.id || source.renderer}_copy`,
      params: cloneParams(source.params)
    };
    drawing.adv.linear_track_slots.splice(idx + 1, 0, copy);
    const axis = axisIndexForCurrentLinearSlots(drawing, drawing.adv.linear_track_slots);
    if (idx < axis) drawing.adv.linear_track_slots_axis_index = axis + 1;
    normalizeCurrentSlots();
  };

  const setLinearTrackSlotEnabled = (slot, enabled) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!slot) return;
    slot.enabled = Boolean(enabled);
  };

  const removeLinearTrackSlot = (index) => {
    const drawing = state.drawings.linear;
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    normalizeCurrentSlots();
    const idx = Number(index);
    if (!Number.isInteger(idx) || idx < 0 || idx >= drawing.adv.linear_track_slots.length) return;
    const axis = axisIndexForCurrentLinearSlots(drawing, drawing.adv.linear_track_slots);
    drawing.adv.linear_track_slots.splice(idx, 1);
    drawing.adv.linear_track_slots_axis_index = idx < axis
      ? Math.max(0, axis - 1)
      : Math.min(axis, drawing.adv.linear_track_slots.length);
    normalizeCurrentSlots();
  };

  /** @param {DrawingState} drawing */
  const wouldLinearTrackSlotMoveCrossAxis = (drawing, fromIndex, toIndex) => {
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

    const axis = axisIndexForCurrentLinearSlots(drawing, normalized);
    const movedPlacement = effectiveLinearSlotPlacement(normalized[from]);
    const targetPlacement = effectiveLinearSlotPlacement(normalized[to]);
    if (movedPlacement === 'overlay' || targetPlacement === 'overlay') return true;
    return (from < axis) !== (to < axis);
  };

  const canMoveLinearTrackSlot = (index, direction) => {
    const drawing = state.drawings.linear;
    const idx = Number(index);
    const step = Number(direction);
    if (!Number.isInteger(idx) || !Number.isInteger(step) || step === 0) return false;
    const target = idx + Math.sign(step);
    return !wouldLinearTrackSlotMoveCrossAxis(drawing, idx, target);
  };

  const moveLinearTrackSlot = (fromIndex, toIndex) => {
    const drawing = state.drawings.linear;
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (wouldLinearTrackSlotMoveCrossAxis(drawing, fromIndex, toIndex)) return;
    normalizeCurrentSlots();
    const from = Number(fromIndex);
    const to = Number(toIndex);
    if (
      !Number.isInteger(from) ||
      !Number.isInteger(to) ||
      from < 0 ||
      to < 0 ||
      from >= drawing.adv.linear_track_slots.length ||
      to >= drawing.adv.linear_track_slots.length ||
      from === to
    ) {
      return;
    }
    const [slot] = drawing.adv.linear_track_slots.splice(from, 1);
    drawing.adv.linear_track_slots.splice(to, 0, slot);
    normalizeCurrentSlots();
  };

  const canMoveLinearTrackSlotAbove = (index) => {
    const drawing = state.drawings.linear;
    const idx = Number(index);
    const normalized = normalizedSlotsForCurrentState(drawing);
    if (!Number.isInteger(idx) || idx < 0 || idx >= normalized.length) return false;
    if (effectiveLinearSlotPlacement(normalized[idx]) === 'overlay') return true;
    return idx >= axisIndexForCurrentLinearSlots(drawing, normalized);
  };

  const canMoveLinearTrackSlotBelow = (index) => {
    const drawing = state.drawings.linear;
    const idx = Number(index);
    const normalized = normalizedSlotsForCurrentState(drawing);
    if (!Number.isInteger(idx) || idx < 0 || idx >= normalized.length) return false;
    if (effectiveLinearSlotPlacement(normalized[idx]) === 'overlay') return true;
    return idx < axisIndexForCurrentLinearSlots(drawing, normalized);
  };

  const canMoveLinearTrackSlotToAxis = (index) => {
    const drawing = state.drawings.linear;
    const idx = Number(index);
    const normalized = normalizedSlotsForCurrentState(drawing);
    if (!Number.isInteger(idx) || idx < 0 || idx >= normalized.length) return false;
    const slot = normalized[idx];
    return ['features', 'annotations'].includes(slot?.renderer) && effectiveLinearSlotPlacement(slot) !== 'overlay';
  };

  /** @param {DrawingState} drawing */
  const moveLinearTrackSlotToPlacement = (drawing, index, placement) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const idx = Number(index);
    if (!Number.isInteger(idx) || idx < 0 || idx >= drawing.adv.linear_track_slots.length) return;
    normalizeCurrentSlots();
    if (idx >= drawing.adv.linear_track_slots.length) return;
    const targetPlacement = normalizePlacement(placement);
    const movingSlot = drawing.adv.linear_track_slots[idx];
    const movingRenderer = normalizeRenderer(movingSlot?.renderer);
    if (targetPlacement === 'overlay' && !['features', 'annotations'].includes(movingRenderer)) return;
    if (movingRenderer === 'annotations' && targetPlacement !== 'overlay') {
      movingSlot.params = cloneParams(movingSlot.params);
      delete movingSlot.params.anchor_slot;
      delete movingSlot.params.cover_anchor;
    }
    if (targetPlacement === 'overlay' && movingRenderer === 'annotations') {
      syncLinearSlotPlacementFromSide(movingSlot, 'overlay');
      movingSlot.params = cloneParams(movingSlot.params);
      bindOnlyLinearAnnotationAnchor(movingSlot);
      normalizeCurrentSlots();
      return;
    }

    if (targetPlacement === 'overlay') {
      const movedPreviousPlacement = effectiveLinearSlotPlacement(movingSlot);
      const existingAxisIndex = drawing.adv.linear_track_slots.findIndex((slot, slotIndex) => (
        slotIndex !== idx &&
        normalizeRenderer(slot?.renderer) === 'features' &&
        effectiveLinearSlotPlacement(slot) === 'overlay'
      ));
      if (existingAxisIndex >= 0) {
        const existingAxisSlot = drawing.adv.linear_track_slots[existingAxisIndex];
        const demotedPlacement = movedPreviousPlacement === 'overlay' ? 'below' : movedPreviousPlacement;
        syncLinearSlotPlacementFromSide(existingAxisSlot, demotedPlacement);
        syncLinearSlotPlacementFromSide(movingSlot, 'overlay');
        drawing.adv.linear_track_slots[existingAxisIndex] = movingSlot;
        drawing.adv.linear_track_slots[idx] = existingAxisSlot;
        drawing.adv.linear_track_slots_axis_index = existingAxisIndex;
        normalizeCurrentSlots();
        return;
      }
    }

    let axis = axisIndexForCurrentLinearSlots(drawing, drawing.adv.linear_track_slots);
    const onAxisIndex = drawing.adv.linear_track_slots.findIndex((slot) => (
      normalizeRenderer(slot?.renderer) === 'features' &&
      effectiveLinearSlotPlacement(slot) === 'overlay'
    ));
    const [slot] = drawing.adv.linear_track_slots.splice(idx, 1);
    if (!slot) return;
    if (idx < axis) axis -= 1;

    syncLinearSlotPlacementFromSide(slot, targetPlacement);
    if (targetPlacement === 'above') {
      drawing.adv.linear_track_slots.splice(axis, 0, slot);
      drawing.adv.linear_track_slots_axis_index = axis + 1;
    } else if (targetPlacement === 'overlay') {
      drawing.adv.linear_track_slots.splice(axis, 0, slot);
      drawing.adv.linear_track_slots_axis_index = axis;
    } else {
      const adjustedOnAxisIndex = onAxisIndex >= 0
        ? (idx === onAxisIndex ? -1 : (idx < onAxisIndex ? onAxisIndex - 1 : onAxisIndex))
        : -1;
      const insertIndex = adjustedOnAxisIndex >= 0 ? adjustedOnAxisIndex + 1 : axis;
      drawing.adv.linear_track_slots.splice(insertIndex, 0, slot);
      drawing.adv.linear_track_slots_axis_index = adjustedOnAxisIndex >= 0 ? adjustedOnAxisIndex : axis;
    }
    normalizeCurrentSlots();
  };

  const moveLinearTrackSlotAbove = (index) => {
    const drawing = state.drawings.linear;
    if (!canMoveLinearTrackSlotAbove(index)) return;
    moveLinearTrackSlotToPlacement(drawing, index, 'above');
  };

  const moveLinearTrackSlotBelow = (index) => {
    const drawing = state.drawings.linear;
    if (!canMoveLinearTrackSlotBelow(index)) return;
    moveLinearTrackSlotToPlacement(drawing, index, 'below');
  };

  const moveLinearTrackSlotToAxis = (index) => {
    const drawing = state.drawings.linear;
    if (!canMoveLinearTrackSlotToAxis(index)) return;
    moveLinearTrackSlotToPlacement(drawing, index, 'overlay');
  };

  const updateLinearTrackSlotRenderer = (slot, renderer = slot?.renderer) => {
    const drawing = state.drawings.linear;
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!slot) return;
    const nextRenderer = normalizeRenderer(renderer);
    slot.params = nextRenderer === normalizeRenderer(slot.renderer)
      ? cloneParams(slot.params)
      : paramsKeptOnRendererChange(slot.params, nextRenderer);
    slot.renderer = nextRenderer;
    if (slot.renderer === 'depth') {
      slot.params.track_index = normalizeTrackIndex(slot.params.track_index) ?? 0;
    }
    if (slot.renderer === 'dinucleotide_content' || slot.renderer === 'dinucleotide_skew') {
      slot.params.nt = normalizeNt(slot.params.nt, drawing.adv.nt);
    }
    if (slot.renderer === 'annotations') {
      slot.side = slot.side === 'overlay' ? 'overlay' : 'above';
      slot.params.set_id = String(slot.params.set_id || drawing.annotationSets?.[0]?.id || '');
      slot.params.overflow = 'error';
      slot.params.show_labels = true;
      slot.params.layer = 'foreground';
    }
    if (slot.renderer !== 'features' && slot.renderer !== 'annotations' && slot.side === 'overlay') slot.side = 'below';
    normalizeCurrentSlots();
  };

  const updateLinearTrackSlotPlacement = (slot, placement) => {
    const drawing = state.drawings.linear;
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!slot) return;
    const index = drawing.adv.linear_track_slots.findIndex((candidate) => candidate === slot);
    if (index >= 0) {
      moveLinearTrackSlotToPlacement(drawing, index, placement);
      return;
    }
    syncLinearSlotPlacementFromSide(slot, placement);
    if (normalizeRenderer(slot.renderer) === 'annotations') {
      slot.params = cloneParams(slot.params);
      if (normalizePlacement(placement) === 'overlay') {
        bindOnlyLinearAnnotationAnchor(slot);
      } else {
        delete slot.params.anchor_slot;
        delete slot.params.cover_anchor;
      }
    }
    normalizeCurrentSlots();
  };

  // Managed Depth rows follow Depth sources (PD-OI-058), as in Circular. A
  // manual row left on a series without a source reports a row issue
  // (PD-OI-083); a logical series itself is kept (PD-OI-025).
  /** @param {DrawingState} drawing */
  const reconcileLinearDepthSlots = (drawing, previousSourced) => {
    const slots = Array.isArray(drawing.adv.linear_track_slots) ? drawing.adv.linear_track_slots : [];
    if (slots.length === 0 && !drawing.adv.linear_track_slots_enabled) return;
    const { slots: nextSlots, additions } = reconcileManagedDepthSlots(/** @type {any} */ ({
      slots,
      previousSourced,
      sourced: linearSourcedDepthTrackIndexesForState(state),
      managedPredicate: (slot) => isDefaultManagedLinearSlot(slot, 'depth')
    }));
    if (additions.length === 0 && nextSlots.length === slots.length) return;
    const previousAxisIndex = clampLinearTrackAxisIndex(drawing.adv.linear_track_slots_axis_index, slots.length);
    if (previousAxisIndex !== null) {
      const removedBeforeAxis = slots
        .slice(0, previousAxisIndex)
        .filter((slot) => !nextSlots.includes(slot)).length;
      drawing.adv.linear_track_slots_axis_index = previousAxisIndex - removedBeforeAxis;
    }
    const seriesCount = linearAvailableDepthTrackCountForState(state);
    const existingIds = new Set(nextSlots.map((slot) => String(slot?.id || '').trim()).filter(Boolean));
    const newSlots = additions.map((trackIndex) => {
      const preferredId = seriesCount <= 1 ? 'depth' : `depth_${trackIndex + 1}`;
      let id = preferredId;
      let suffix = 2;
      while (existingIds.has(id)) {
        id = `${preferredId}_${suffix}`;
        suffix += 1;
      }
      existingIds.add(id);
      return defaultSlot('depth', { id, side: 'below', params: { track_index: trackIndex } });
    });
    drawing.adv.linear_track_slots.splice(0, slots.length, ...nextSlots, ...newSlots);
    normalizeCurrentSlots();
  };

  // The only entry for Depth source changes: run the change, then reconcile.
  // The Linear drawing's Show Depth goes off with its last Depth
  // source; the other mode's drawing is not touched (OV-82, OV-101).
  const changeLinearDepthSources = (mutate) => {
    const drawing = state.drawings.linear;
    const previousSourced = linearSourcedDepthTrackIndexesForState(state);
    mutate();
    reconcileLinearDepthSlots(drawing, previousSourced);
    if (previousSourced.length > 0 && linearSourcedDepthTrackIndexesForState(state).length === 0) drawing.form.show_depth = false;
  };

  const linearTrackStackEntries = () => {
    const drawing = state.drawings.linear;
    const slots = Array.isArray(drawing.adv.linear_track_slots) ? drawing.adv.linear_track_slots : [];
    const axisIndex = axisIndexForCurrentLinearSlots(drawing, slots);
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

  const linearTrackSlotLegendLabelPlaceholder = (slot) => {
    const drawing = state.drawings.linear;
    const renderer = normalizeRenderer(slot?.renderer);
    if (renderer === 'dinucleotide_content' || renderer === 'dinucleotide_skew') {
      const nt = normalizeNt(slot?.params?.nt ?? slot?.params?.dinucleotide, drawing.adv.nt);
      return renderer === 'dinucleotide_content' ? `${nt} content` : `${nt} skew`;
    }
    if (renderer === 'depth') return 'Depth';
    return 'Legend label';
  };

  const linearTrackSlotUsesPresetGeometry = (slot) => {
    if (!slot || typeof slot !== 'object') return false;
    const renderer = RENDERER_ALIASES[String(slot.renderer || '').trim().toLowerCase()] || String(slot.renderer || '').trim().toLowerCase();
    if (!SUPPORTED_RENDERERS.includes(renderer)) return false;
    return hasBlankLinearSlotGeometry(slot);
  };

  const parsePositivePxNumber = (value) => {
    try {
      return parseOptionalPixel(value, 'Linear track height', { allowZero: false });
    } catch {
      return null;
    }
  };

  /** @param {DrawingState} drawing */
  const ensureDepthTrackConfigForSlotIndex = (drawing, trackIndex) => {
    const idx = Math.max(0, Number(trackIndex) || 0);
    if (!Array.isArray(drawing.adv.depth_tracks)) drawing.adv.depth_tracks = [];
    while (drawing.adv.depth_tracks.length <= idx) {
      const nextIndex = drawing.adv.depth_tracks.length;
      drawing.adv.depth_tracks.push({
        label: nextIndex === 0 ? 'Depth' : `Depth ${nextIndex + 1}`,
        color: nextIndex === 0 ? String(drawing.adv.depth_color || '#4A90E2') : '',
        height: null,
        large_tick_interval: null,
        small_tick_interval: null,
        tick_font_size: null
      });
    }
    if (!drawing.adv.depth_tracks[idx] || typeof drawing.adv.depth_tracks[idx] !== 'object' || Array.isArray(drawing.adv.depth_tracks[idx])) {
      drawing.adv.depth_tracks[idx] = {
        label: idx === 0 ? 'Depth' : `Depth ${idx + 1}`,
        color: idx === 0 ? String(drawing.adv.depth_color || '#4A90E2') : '',
        height: null,
        large_tick_interval: null,
        small_tick_interval: null,
        tick_font_size: null
      };
    }
    return drawing.adv.depth_tracks[idx];
  };

  const linearDepthTrackIndexForSlot = (slot) => {
    if (!slot || normalizeRenderer(slot.renderer) !== 'depth') return null;
    const params = cloneParams(slot.params);
    return normalizeTrackIndex(params.track_index) ?? 0;
  };

  /** @param {DrawingState} drawing */
  const heightTextFromDepthTrackConfig = (drawing, trackIndex) => {
    const config = ensureDepthTrackConfigForSlotIndex(drawing, trackIndex);
    const height = parsePositivePxNumber(config.height);
    return height === null ? '' : String(height);
  };

  const syncLinearDepthSlotHeightsFromDepthTracks = (trackIndex = null) => {
    const drawing = state.drawings.linear;
    const slots = Array.isArray(drawing.adv.linear_track_slots) ? drawing.adv.linear_track_slots : [];
    slots.forEach((slot) => {
      const slotTrackIndex = linearDepthTrackIndexForSlot(slot);
      if (slotTrackIndex === null) return;
      if (trackIndex !== null && Number(trackIndex) !== slotTrackIndex) return;
      slot.height = heightTextFromDepthTrackConfig(drawing, slotTrackIndex);
    });
  };

  const linearTrackSlotHeightValue = (slot) => {
    const drawing = state.drawings.linear;
    if (normalizeRenderer(slot?.renderer) !== 'depth') {
      return String(slot?.height || '');
    }
    const trackIndex = linearDepthTrackIndexForSlot(slot) ?? 0;
    if (isManualSlotValue(slot?.height)) return String(slot.height);
    return heightTextFromDepthTrackConfig(drawing, trackIndex);
  };

  const linearSlotManualValue = (slot, field) => {
    if (!slot) return '';
    if (field === 'height') return linearTrackSlotHeightValue(slot);
    if (field === 'spacing') return slot.spacing;
    return '';
  };

  const selectedResultIndexValue = () => Number(state?.selectedResultIndex?.value ?? 0) || 0;

  const resolvedLinearSlotGeometry = (slotId) => findTrackSlotGeometry(/** @type {any} */ ({
    geometry: String(state?.trackSlotResolvedGeometry?.value?.mode || '') === 'linear'
      ? state.trackSlotResolvedGeometry.value
      : null,
    resultIndex: selectedResultIndexValue(),
    recordIndex: 0,
    slotId
  }));

  /** @param {DrawingState} drawing */
  const estimateLinearSlotGeometry = (drawing, slot) => {
    const renderer = normalizeRenderer(slot?.renderer);
    const configuredHeight = parsePositivePxNumber(slot?.height);
    let heightPx = configuredHeight ?? ESTIMATED_LINEAR_GC_HEIGHT_PX;
    if (renderer === 'features') heightPx = 0;
    else if (renderer === 'spacer') heightPx = configuredHeight ?? 0;
    else if (renderer === 'depth') {
      heightPx = parsePositivePxNumber(drawing.adv.depth_height) ?? ESTIMATED_LINEAR_DEPTH_HEIGHT_PX;
    } else if (renderer === 'dinucleotide_content' || renderer === 'dinucleotide_skew') {
      heightPx = parsePositivePxNumber(drawing.adv.gc_height) ?? ESTIMATED_LINEAR_GC_HEIGHT_PX;
    }
    const spacingAfterPx = renderer === 'depth' ? ESTIMATED_LINEAR_DEPTH_SPACING_PX : 0;
    return {
      heightPx,
      spacingAfterPx,
      baseYOffsetPx: 0,
      finalYOffsetPx: 0
    };
  };

  // Only a rendered row has resolved geometry; a disabled row shows the estimate.
  /** @param {DrawingState} drawing */
  const linearTrackSlotDisplayGeometry = (drawing, slot) => {
    const resolved = slot?.enabled !== false ? resolvedLinearSlotGeometry(slot?.id) : null;
    return resolved || estimateLinearSlotGeometry(drawing, slot);
  };

  // The row index is part of the shared template call; Linear geometry is
  // found by slot ID alone.
  const linearTrackSlotGeometryAutoText = (slot, _slotIndex, field) => {
    const drawing = state.drawings.linear;
    if (isManualSlotValue(linearSlotManualValue(slot, field))) return '';
    const geometry = linearTrackSlotDisplayGeometry(drawing, slot);
    const value = field === 'height'
      ? geometry.heightPx
      : (field === 'spacing' ? geometry.spacingAfterPx : null);
    const autoText = formatPxAuto(value);
    if (normalizeRenderer(slot?.renderer) === 'features') {
      return autoText
        ? autoText.replace('(auto)', '(auto; varies by record)')
        : 'record-specific (auto)';
    }
    if (field === 'height' || field === 'spacing') return autoText;
    return '';
  };

  const linearTrackSlotGeometryUnitSuffix = (slot, field) => {
    const text = String(linearSlotManualValue(slot, field) ?? '').trim();
    if (!text || /px$/i.test(text)) return '';
    return 'px';
  };

  const linearTrackSlotGeometryHasManual = (slot, field) => (
    isManualSlotValue(linearSlotManualValue(slot, field))
  );

  const setLinearTrackSlotHeight = (slot, value) => {
    const drawing = state.drawings.linear;
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!slot) return;
    const text = String(value ?? '').trim();
    slot.height = text;
    if (normalizeRenderer(slot.renderer) !== 'depth') return;
    const trackIndex = linearDepthTrackIndexForSlot(slot) ?? 0;
    const config = ensureDepthTrackConfigForSlotIndex(drawing, trackIndex);
    try {
      config.height = parseOptionalPixel(text, 'Linear track height', { allowZero: false });
    } catch {
      // The slot retains the invalid text; keep the last valid shared config.
    }
  };

  const linearTrackSlotHasSkewColorOverride = (slot, key) => (
    normalizeRenderer(slot?.renderer) === 'dinucleotide_skew' &&
    normalizeOptionalText(slot?.params?.[key]) !== null
  );

  const linearTrackSlotSkewColorValue = (slot, key) => {
    const drawing = state.drawings.linear;
    return resolveTrackSlotSkewColorValue(/** @type {any} */ ({
      slot,
      key,
      currentColors: drawing.currentColors,
      paletteDefinitions: state.paletteDefinitions,
      selectedPalette: drawing.selectedPalette
    }));
  };

  const setLinearTrackSlotSkewColor = (slot, key, value) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!slot || normalizeRenderer(slot.renderer) !== 'dinucleotide_skew' || !['positive_color', 'negative_color'].includes(key)) return;
    slot.params = cloneParams(slot.params);
    const color = normalizeColorParam(value);
    if (color === null) delete slot.params[key];
    else slot.params[key] = color;
  };

  const clearLinearTrackSlotSkewColor = (slot, key) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!slot || normalizeRenderer(slot.renderer) !== 'dinucleotide_skew' || !['positive_color', 'negative_color'].includes(key)) return;
    slot.params = cloneParams(slot.params);
    delete slot.params[key];
  };

  return {
    linearTrackRenderers,
    linearTrackSlotEditorKey,
    linearTrackRendererLabel,
    normalizeLinearTrackSlots: normalizeCurrentSlots,
    changeLinearDepthSources,
    ...featureSlotEdits(changeTrackLayout, {
      resetLinearTrackSlotsFromSimpleControls,
      setLinearTrackSlotsEnabled,
      addLinearTrackSlot,
      duplicateLinearTrackSlot,
      removeLinearTrackSlot,
      setLinearTrackSlotEnabled,
      moveLinearTrackSlot,
      moveLinearTrackSlotAbove,
      moveLinearTrackSlotBelow,
      moveLinearTrackSlotToAxis,
      updateLinearTrackSlotRenderer,
      updateLinearTrackSlotPlacement
    }),
    canAddLinearTrackRenderer,
    canDuplicateLinearTrackSlot,
    canMoveLinearTrackSlot,
    canMoveLinearTrackSlotAbove,
    canMoveLinearTrackSlotBelow,
    canMoveLinearTrackSlotToAxis,
    linearTrackSlotIssue,
    linearTrackGlobalIssues,
    linearAnnotationAnchorOptions,
    linearAnnotationAnchorIsKnown,
    annotationTrackMarkOptions: ANNOTATION_MARK_OPTIONS,
    linearAnnotationMarkSelected,
    setLinearAnnotationMarkSelected,
    linearAnnotationLaneGapValue: (slot) => linearAnnotationNumberValue(slot, 'lane_gap_px', 3),
    setLinearAnnotationLaneGap: (slot, value) => setLinearAnnotationNumber(slot, 'lane_gap_px', value, 3),
    linearAnnotationPaddingValue: (slot) => linearAnnotationNumberValue(slot, 'padding_px', 2),
    setLinearAnnotationPadding: (slot, value) => setLinearAnnotationNumber(slot, 'padding_px', value, 2),
    linearAnnotationCoverAnchor,
    setLinearAnnotationCoverAnchor,
    linearTrackSlotHeightValue,
    linearTrackSlotGeometryAutoText,
    linearTrackSlotGeometryHasManual,
    linearTrackSlotGeometryUnitSuffix,
    setLinearTrackSlotHeight,
    linearTrackSlotHasSkewColorOverride,
    linearTrackSlotSkewColorValue,
    setLinearTrackSlotSkewColor,
    clearLinearTrackSlotSkewColor,
    syncLinearDepthSlotHeightsFromDepthTracks,
    linearTrackSlots: () => {
      const drawing = state.drawings.linear;
      return Array.isArray(drawing.adv.linear_track_slots) ? drawing.adv.linear_track_slots : [];
    },
    linearTrackStackEntries,
    linearTrackSlotCliSpec: (slot) => buildLinearTrackSlotSpec(slot),
    linearTrackSlotDisplayLabel: (slot) => linearTrackRendererLabel(slot?.renderer),
    linearTrackSlotDisplayMeta: (slot) => buildLinearTrackSlotSpec(slot),
    linearTrackSlotLegendLabelPlaceholder,
    linearTrackSlotPlacementLabel,
    linearTrackSlotUsesPresetGeometry
  };
};
