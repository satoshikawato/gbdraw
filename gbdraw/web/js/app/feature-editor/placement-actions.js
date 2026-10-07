// @ts-check
/** @import { FeaturePlacementTarget } from '../../services/feature-placement.js' */
/** @import { ChangeTrackLayout } from '../track-slot-edits.js' */
import {
  canonicalFeaturePlacements,
  featureIdentityKeyOf,
  parseFeatureIdentityKey
} from '../../services/feature-placement.js';
import { resolveCircularTrackFeaturePlacement } from '../circular-track-slots.js';
import { createDefaultLinearTrackSlots, effectiveLinearSlotPlacement } from '../linear-track-slots.js';
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
  const draft = () => ({ mode: state.mode.value, form: state.form, adv: state.adv });
  const targetsFor = (features) => {
    const mode = state.mode.value;
    if (!features.length || getCommittedRequest()?.mode !== mode) return [];
    const items = state.featureCatalog.value?.items || [];
    if (!features.every((feature) => items.some((item) => item.recordKeys.includes(feature?.record_key)))) return [];
    return draftPlacementTargets(draft());
  };
  const choices = (features) => {
    const targets = targetsFor(features);
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
    const busy = state.sessionOperationAvailability?.();
    if (busy) return busy;
    const choice = choices(features).find((entry) => entry.value === value);
    if (!choice?.enabled) throw new Error(choice?.reason || 'Unknown feature placement.');
    return runUndoable(features.length === 1 ? 'Change feature placement' : 'Change selected feature placements', () => {
      for (const feature of features) {
        const key = featureIdentityKeyOf(feature);
        if (value === 'auto') delete state.featurePlacementOverrides[key];
        else {
          const row = { ...parseFeatureIdentityKey(key), placement: /** @type {FeaturePlacementTarget} */ (value === 'main'
            ? { kind: 'main' } : { kind: 'lane', side: value, level: 1 }) };
          canonicalFeaturePlacements({ [key]: row });
          state.featurePlacementOverrides[key] = row;
        }
      }
    });
  };

  // Q3 (Owner, 2026-10-04) and R10: the one transition for an edit of a draft
  // feature-slot input. It applies the edit and, when a lane placement of
  // either mode loses its lane, restores the inputs and asks. Reset applies the
  // edit and removes those rows as one History step; Cancel keeps the control's
  // value and records none. The control's own History adapter records an edit
  // that loses nothing (R11). Restores (Undo/Redo, Session load, Import, Reset
  // Settings) install state as is; Generate names what they leave (R6).
  // The composition root passes `changeTrackLayout` to the track stack editors
  // as their one port (R13); each routes its feature-slot edits through it.
  const layoutChange = reactive({ open: false, count: 0, setting: '', value: '', scope: '' });
  let pendingLayout = null;
  const laneSides = () => Object.fromEntries(Object.keys(SIDES).map((mode) => [mode,
    draftPlacementTargets({ mode, form: state.form, adv: state.adv }).map((target) => target.side).filter(Boolean)]));
  const lostLaneKeys = (before) => {
    const after = laneSides();
    return Object.keys(state.featurePlacementOverrides).filter((key) => {
      const { scope, placement } = state.featurePlacementOverrides[key] || {};
      return before[scope]?.includes(placement?.side) && !after[scope].includes(placement?.side);
    });
  };
  const focusAfterRender = (find) => nextTick().then(() => find()?.focus?.());
  /**
   * The `ChangeTrackLayout` port of the track stack editors (app/track-slot-edits.js).
   * @type {ChangeTrackLayout}
   */
  const changeTrackLayout = (apply, control = null) => {
    const busy = state.sessionOperationAvailability?.();
    if (busy) return busy;
    if (!Object.values(state.featurePlacementOverrides).some((row) => row?.placement?.kind === 'lane')) return apply();
    const before = laneSides();
    const saved = saveTrackLayout(state);
    const result = apply();
    const lost = lostLaneKeys(before);
    if (!lost.length) return result;
    restoreTrackLayout(state, saved);
    const checkbox = control?.type === 'checkbox';
    const value = checkbox ? (control.checked ? 'On' : 'Off') : control?.selectedOptions?.[0]?.text || '';
    // The checkbox shows the kept value again; a bound select re-renders to it.
    if (checkbox) control.checked = !control.checked;
    const scopes = [...new Set(lost.map((key) => state.featurePlacementOverrides[key].scope))]
      .filter((scope) => scope !== state.mode.value);
    pendingLayout = { control, apply };
    Object.assign(layoutChange, { open: true, count: lost.length, value,
      scope: scopes.map((scope) => `${scope[0].toUpperCase()}${scope.slice(1)} `).join(''),
      setting: control?.getAttribute?.('aria-label') || control?.labels?.[0]?.textContent.trim()
        || (control?.tagName === 'BUTTON' && control.textContent.trim()) || control?.title || 'This change' });
    void focusAfterRender(() => globalThis.document?.querySelector('[data-placement-layout-primary]'));
    return false;
  };
  const resolveLayoutChange = async (choice) => {
    const change = pendingLayout;
    pendingLayout = null;
    layoutChange.open = false;
    if (!change) return false;
    const result = choice === 'reset' && await runUndoable('Change setting and reset Feature placements', () => {
      const before = laneSides();
      change.apply();
      for (const key of lostLaneKeys(before)) delete state.featurePlacementOverrides[key];
    });
    void focusAfterRender(() => change.control);
    return result;
  };
  const changeLayoutSetting = (event, field) => {
    const control = event.target;
    const value = control.type === 'checkbox' ? control.checked : control.value;
    return changeTrackLayout(() => { state.form[field] = value; }, control);
  };

  return { choices, setPlacement, layoutChange, resolveLayoutChange, changeLayoutSetting, changeTrackLayout,
    valueFor: (feature) => {
      // The control lists this mode's placements; a feature of the other
      // mode's Result reads as Auto until that mode is active (R2).
      const row = state.featurePlacementOverrides[featureIdentityKeyOf(feature)];
      const target = row?.scope === state.mode.value ? row.placement : null;
      return target?.kind === 'main' ? 'main' : target?.side || 'auto';
    } };
};
