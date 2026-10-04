import { canonicalFeaturePlacements, placementAppliesToMode } from '../../services/feature-placement.js';
import { resolveCircularTrackFeaturePlacement } from '../circular-track-slots.js';
import { createDefaultLinearTrackSlots, effectiveLinearSlotPlacement } from '../linear-track-slots.js';
import { validateCustomTrackPlan } from '../track-slot-validation.js';

const identity = (feature) => [feature?.record_key, feature?.biological_feature_id];
const keyFor = (feature) => JSON.stringify(identity(feature));
const SIDES = { circular: ['outward', 'inward'], linear: ['above', 'below'] };

// The placement targets of a draft's feature slot: the one availability
// predicate of the popup choices and of a layout change (R3, R10). It resolves
// the draft through the same slot owner used by request projection; Result
// geometry continues to describe the last successful Generate.
const draftPlacementTargets = ({ mode, form, adv }) => {
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

export const createFeaturePlacementActions = ({
  state, history, getCommittedRequest, isCurrentFeature,
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
        if (identity(feature).some((part) => !part) || !isCurrentFeature(feature)) return false;
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
    return history.runUndoable(features.length === 1 ? 'Change feature placement' : 'Change selected feature placements', () => {
      for (const feature of features) {
        const key = keyFor(feature);
        if (value === 'auto') delete state.featurePlacementOverrides[key];
        else {
          const [recordKey, biologicalFeatureId] = identity(feature);
          const row = { recordKey, biologicalFeatureId, placement: value === 'main'
            ? { kind: 'main' } : { kind: 'lane', side: value, level: 1 } };
          canonicalFeaturePlacements([row], state.mode.value);
          state.featurePlacementOverrides[key] = row;
        }
      }
    });
  };

  // Q3 (Owner, 2026-10-04): a layout edit that would leave this mode's lane
  // placements undrawable asks first. Reset applies the edit and removes those
  // rows as one History step; cancel keeps the setting and records no step. The
  // other mode's rows wait for their mode (R2). Session load, Undo/Redo, and
  // Reset Settings do not come here; Generate names what they leave (R6).
  const layoutChange = reactive({ open: false, count: 0, setting: '', value: '' });
  let pendingLayout = null;
  const lostLaneKeys = (next) => {
    const sides = (entry) => draftPlacementTargets(entry).map((target) => target.side).filter(Boolean);
    const [had, has] = [sides(draft()), sides(next(draft()))];
    return Object.keys(state.featurePlacementOverrides).filter((key) => {
      const side = state.featurePlacementOverrides[key]?.placement?.side;
      return had.includes(side) && !has.includes(side);
    });
  };
  const focusAfterRender = (find) => nextTick().then(() => find()?.focus?.());
  // R10: the control's own transition reconciles the rows before it commits.
  const changeLayout = (control, previous, next, apply) => {
    const busy = state.sessionOperationAvailability?.();
    if (busy) return busy;
    const count = lostLaneKeys(next).length;
    if (!count) return history.runUndoable('Change setting', apply);
    const checkbox = control.type === 'checkbox';
    const value = checkbox ? (control.checked ? 'On' : 'Off') : control.selectedOptions?.[0]?.text;
    control[checkbox ? 'checked' : 'value'] = previous;
    pendingLayout = { control, next, apply };
    Object.assign(layoutChange, { open: true, count, value,
      setting: control.getAttribute?.('aria-label') || control.labels?.[0]?.textContent.trim() || 'the setting' });
    void focusAfterRender(() => globalThis.document?.querySelector('[data-placement-layout-primary]'));
    return false;
  };
  const resolveLayoutChange = async (choice) => {
    const change = pendingLayout;
    pendingLayout = null;
    layoutChange.open = false;
    if (!change) return false;
    const result = choice === 'reset' && await history.runUndoable('Change setting and reset Feature placements', () => {
      const lost = lostLaneKeys(change.next);
      change.apply();
      for (const key of lost) delete state.featurePlacementOverrides[key];
    });
    void focusAfterRender(() => change.control);
    return result;
  };
  const changeLayoutSetting = (event, field) => {
    const control = event.target;
    const value = control.type === 'checkbox' ? control.checked : control.value;
    return changeLayout(control, state.form[field],
      (entry) => ({ ...entry, form: { ...entry.form, [field]: value } }),
      () => { state.form[field] = value; });
  };
  // A custom slot's side; `update` is the slot editor's transition.
  const changeFeatureSlotSide = (event, slot, update) => {
    const { value } = event.target;
    const mode = state.mode.value;
    const slots = `${mode}_track_slots`;
    const changed = mode === 'circular'
      ? { ...slot, side: null, params: { ...slot.params, lane_direction: value || undefined } }
      : { ...slot, side: value };
    return changeLayout(event.target, mode === 'circular' ? slot.params?.lane_direction ?? '' : slot.side,
      (entry) => ({ ...entry, adv: { ...entry.adv, [slots]: entry.adv[slots].map((item) => item === slot ? changed : item) } }),
      () => update(slot, value));
  };

  return { choices, setPlacement, layoutChange, resolveLayoutChange, changeLayoutSetting, changeFeatureSlotSide,
    valueFor: (feature) => {
      // The request projection skips another mode's lane; the popup reads it as Auto.
      const row = state.featurePlacementOverrides[keyFor(feature)];
      const target = placementAppliesToMode(row, state.mode.value) ? row?.placement : null;
      return target?.kind === 'main' ? 'main' : target?.side || 'auto';
    } };
};
