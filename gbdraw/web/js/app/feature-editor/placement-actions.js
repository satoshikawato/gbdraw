// @ts-check
/** @import { DrawingState } from '../../state.js' */
/** @import { FeaturePlacementDraftRow, FeaturePlacementTarget } from '../../services/feature-placement.js' */
/** @import { ChangeTrackLayout, LayoutControl } from '../track-slot-edits.js' */
import {
  canonicalFeaturePlacements,
  featureIdentityKeyOf,
  parseFeatureIdentityKey
} from '../../services/feature-placement.js';
import { resolveCircularTrackFeaturePlacement } from '../../services/circular-track-slot-model.js';
import { createDefaultLinearTrackSlots, effectiveLinearSlotPlacement } from '../../services/linear-track-slot-model.js';
import { validateCustomTrackPlan } from '../../services/track-slot-validation.js';

const SIDES = { circular: ['outward', 'inward'], linear: ['above', 'below'] };

// The placement targets of a draft's feature slot: the one availability
// predicate of the popup choices and of a layout change (R3, R10). It resolves
// the draft through the same slot owner used by request projection; Result
// geometry continues to describe the last successful Generate.
export const draftPlacementTargets = ({ mode, form, adv }) => {
  const trackType = mode === 'circular' ? form.track_type : form.linear_track_layout;
  const slot = adv[`${mode}_track_slots_enabled`]
    ? validateCustomTrackPlan({
      mode, slots: adv[`${mode}_track_slots`],
      axisIndex: adv[`${mode}_track_slots_axis_index`], trackType
    }).enabledSlots.find((entry) => entry.renderer === 'features')
    : (mode === 'circular' ? {} : createDefaultLinearTrackSlots({ trackLayout: trackType })[0]);
  if (!slot) return [];
  const bidirectional = mode === 'circular'
    ? resolveCircularTrackFeaturePlacement(slot, trackType).laneDirection === 'split'
    : effectiveLinearSlotPlacement(slot) === 'overlay' && !form.separate_strands;
  return [{ kind: 'main' }, ...(bidirectional
    ? SIDES[mode].map((side) => ({ kind: 'lane', side, level: 1 })) : [])];
};

// The draft inputs of draftPlacementTargets, saved and restored in place so
// the stack rows keep their identity (R10).
const LAYOUT_FORM_FIELDS = ['track_type', 'linear_track_layout', 'separate_strands'];
const copy = (value) => (Array.isArray(value) ? value.map(copy) : value && typeof value === 'object'
  ? Object.fromEntries(Object.entries(value).map(([key, entry]) => [key, copy(entry)])) : value);
/** @param {Record<string, any>} draft The draft's `form` and `adv`. */
export const saveTrackLayout = ({ form, adv }) => ({
  form: Object.fromEntries(LAYOUT_FORM_FIELDS.map((field) => [field, form[field]])),
  stacks: Object.keys(SIDES).map((mode) => {
    const slots = adv[`${mode}_track_slots`];
    return { mode, slots, enabled: adv[`${mode}_track_slots_enabled`], axis: adv[`${mode}_track_slots_axis_index`],
      rows: (Array.isArray(slots) ? slots : []).map((slot) => [slot, copy(slot)]) };
  })
});
/**
 * @param {Record<string, any>} draft The draft's `form` and `adv`.
 * @param {Record<string, any>} saved The result of `saveTrackLayout`.
 */
export const restoreTrackLayout = ({ form, adv }, saved) => {
  Object.assign(form, saved.form);
  for (const { mode, slots, enabled, axis, rows } of saved.stacks) {
    Object.assign(adv, { [`${mode}_track_slots_enabled`]: enabled, [`${mode}_track_slots_axis_index`]: axis,
      [`${mode}_track_slots`]: slots });
    for (const [slot, content] of rows) {
      Object.keys(slot).forEach((key) => delete slot[key]);
      Object.assign(slot, content);
    }
    if (Array.isArray(slots)) slots.splice(0, slots.length, ...rows.map(([slot]) => slot));
  }
};

/**
 * @typedef {object} FeaturePlacementActionsOptions
 * @property {Record<string, any>} state App state (state.js; not yet typed).
 * @property {(label: string, fn: () => any, options?: Record<string, any>) => any} runUndoable History's undoable step.
 * @property {() => ({ mode?: string } | null)} getCommittedRequest The committed canonical request.
 * @property {(feature: Record<string, any>) => boolean} isCurrentFeature Whether the feature belongs to the displayed Result.
 * @property {<T extends object>(value: T) => T} [reactive] Vue `reactive`
 * @property {() => Promise<any>} [nextTick] Vue `nextTick`
 */

/** @param {FeaturePlacementActionsOptions} options */
export const createFeaturePlacementActions = ({
  state, runUndoable, getCommittedRequest, isCurrentFeature,
  reactive = (value) => value, nextTick = () => Promise.resolve()
}) => {
  /** @param {DrawingState} drawing */
  const draft = (drawing) => ({ mode: state.mode.value, form: drawing.form, adv: drawing.adv });
  /** @param {DrawingState} drawing */
  const targetsFor = (drawing, features) => {
    const mode = state.mode.value;
    if (!features.length || getCommittedRequest()?.mode !== mode) return [];
    const items = state.featureCatalog.value?.items || [];
    if (!features.every((feature) => items.some((item) => item.recordKeys.includes(feature?.record_key)))) return [];
    return draftPlacementTargets(draft(drawing));
  };
  const choices = (features) => {
    const drawing = state.activeDrawing();
    const targets = targetsFor(drawing, features);
    return ['auto', 'main', ...SIDES[state.mode.value]].map((value) => {
      const enabled = features.length > 0 && features.every((feature) => {
        if (!featureIdentityKeyOf(feature) || !isCurrentFeature(feature)) return false;
        if (value === 'auto') return true;
        return targets.some((target) => value === 'main' ? target.kind === 'main' : target.side === value);
      });
      return { value, enabled, label: value === 'auto' ? 'Auto' : value === 'main' ? 'Main'
        : `${value[0].toUpperCase()}${value.slice(1)} lane 1`,
      reason: enabled ? '' : 'Unavailable in the current draft feature slot.' };
    });
  };
  const setPlacement = (features, value) => {
    const drawing = state.activeDrawing();
    const busy = state.sessionOperationAvailability?.();
    if (busy) return busy;
    const choice = choices(features).find((entry) => entry.value === value);
    if (!choice?.enabled) throw new Error(choice?.reason || 'Unknown feature placement.');
    return runUndoable(features.length === 1 ? 'Change feature placement' : 'Change selected feature placements', () => {
      for (const feature of features) {
        const key = featureIdentityKeyOf(feature);
        if (value === 'auto') delete drawing.featurePlacementOverrides[key];
        else {
          // `key` comes from featureIdentityKeyOf, so parseFeatureIdentityKey reads it back; canonicalFeaturePlacements validates the row.
          const row = /** @type {FeaturePlacementDraftRow} */ ({ ...parseFeatureIdentityKey(key), placement: /** @type {FeaturePlacementTarget} */ (value === 'main'
            ? { kind: 'main' } : { kind: 'lane', side: value, level: 1 }) });
          canonicalFeaturePlacements({ [key]: row }, state.mode.value);
          drawing.featurePlacementOverrides[key] = row;
        }
      }
    });
  };

  // Q3 (Owner, 2026-10-04) and R10: the one transition for an edit of a draft
  // feature-slot input. It applies the edit and, when a lane placement of the
  // drawing loses its lane, restores the inputs and asks. A drawing's
  // placements are of its own mode (PD-OI-086). Reset applies the
  // edit and removes those rows as one History step; Cancel keeps the control's
  // value and records none. The control's own History adapter records an edit
  // that loses nothing (R11). Restores (Undo/Redo, Session load, Import, Reset
  // Settings) install state as is; Generate names what they leave (R6).
  // The composition root passes `changeTrackLayout` to the track stack editors
  // as their one port (R13); each routes its feature-slot edits through it.
  const layoutChange = reactive({ open: false, count: 0, setting: '', value: '' });
  /** @type {{ control: LayoutControl | null, apply: () => any } | null} */
  let pendingLayout = null;
  /** @param {DrawingState} drawing */
  const laneSides = (drawing) => draftPlacementTargets({ mode: state.mode.value, form: drawing.form, adv: drawing.adv })
    .map((target) => target.side).filter(Boolean);
  /**
   * @param {DrawingState} drawing
   * @param {string[]} before
   */
  const lostLaneKeys = (drawing, before) => {
    const after = laneSides(drawing);
    return Object.keys(drawing.featurePlacementOverrides).filter((key) => {
      const { placement } = drawing.featurePlacementOverrides[key] || {};
      return before.includes(placement?.side) && !after.includes(placement?.side);
    });
  };
  const focusAfterRender = (find) => nextTick().then(() => find()?.focus?.());
  /**
   * The `ChangeTrackLayout` port of the track stack editors (app/track-slot-edits.js).
   * @type {ChangeTrackLayout}
   */
  const changeTrackLayout = (apply, control = null) => {
    const drawing = state.activeDrawing();
    const busy = state.sessionOperationAvailability?.();
    if (busy) return busy;
    if (!Object.values(drawing.featurePlacementOverrides).some((row) => row?.placement?.kind === 'lane')) return apply();
    const before = laneSides(drawing);
    const saved = saveTrackLayout(drawing);
    const result = apply();
    const lost = lostLaneKeys(drawing, before);
    if (!lost.length) return result;
    restoreTrackLayout(drawing, saved);
    const checkbox = control?.type === 'checkbox';
    const value = checkbox ? (control.checked ? 'On' : 'Off') : control?.selectedOptions?.[0]?.text || '';
    // The checkbox shows the kept value again; a bound select re-renders to it.
    if (checkbox) control.checked = !control.checked;
    pendingLayout = { control, apply };
    Object.assign(layoutChange, { open: true, count: lost.length, value,
      setting: control?.getAttribute?.('aria-label') || control?.labels?.[0]?.textContent.trim()
        || (control?.tagName === 'BUTTON' && control.textContent.trim()) || control?.title || 'This change' });
    void focusAfterRender(() => globalThis.document?.querySelector('[data-placement-layout-primary]'));
    return false;
  };
  const resolveLayoutChange = async (choice) => {
    const drawing = state.activeDrawing();
    const change = pendingLayout;
    pendingLayout = null;
    layoutChange.open = false;
    if (!change) return false;
    const result = choice === 'reset' && await runUndoable('Change setting and reset Feature placements', () => {
      const before = laneSides(drawing);
      change.apply();
      for (const key of lostLaneKeys(drawing, before)) delete drawing.featurePlacementOverrides[key];
    });
    void focusAfterRender(() => change.control);
    return result;
  };
  const changeLayoutSetting = (event, field) => {
    const drawing = state.activeDrawing();
    const control = event.target;
    const value = control.type === 'checkbox' ? control.checked : control.value;
    return changeTrackLayout(() => { drawing.form[field] = value; }, control);
  };

  return { choices, setPlacement, layoutChange, resolveLayoutChange, changeLayoutSetting, changeTrackLayout,
    valueFor: (feature) => {
      const drawing = state.activeDrawing();
      // The control lists the shown mode's drawing's placements (R2).
      const target = drawing.featurePlacementOverrides[featureIdentityKeyOf(feature)]?.placement;
      return target?.kind === 'main' ? 'main' : target?.side || 'auto';
    } };
};
