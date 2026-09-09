import { canonicalFeaturePlacements } from '../../services/feature-placement.js';
import { resolveCircularTrackFeaturePlacement } from '../circular-track-slots.js';
import { createDefaultLinearTrackSlots, effectiveLinearSlotPlacement } from '../linear-track-slots.js';
import { validateCustomTrackPlan } from '../track-slot-validation.js';

const identity = (feature) => [feature?.record_key, feature?.biological_feature_id];
const keyFor = (feature) => JSON.stringify(identity(feature));

export const createFeaturePlacementActions = ({ state, history, getCommittedRequest, isCurrentFeature }) => {
  const targetsFor = (features, sides) => {
    const mode = state.mode.value;
    if (!features.length || getCommittedRequest()?.mode !== mode) return [];
    const items = state.featureCatalog.value?.items || [];
    if (!features.every((feature) => items.some((item) => item.recordKeys.includes(feature?.record_key)))) return [];
    const { form, adv } = state;
    const trackType = mode === 'circular' ? form.track_type : form.linear_track_layout;
    // Resolve the draft through the same slot owner used by request projection.
    // Result geometry continues to describe the last successful Generate.
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
      ? sides.map((side) => ({ kind: 'lane', side, level: 1 })) : [])];
  };
  const choices = (features) => {
    const sides = state.mode.value === 'circular' ? ['outward', 'inward'] : ['above', 'below'];
    const targets = targetsFor(features, sides);
    return ['auto', 'main', ...sides].map((value) => {
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
  return { choices, setPlacement, valueFor: (feature) => {
    const target = state.featurePlacementOverrides[keyFor(feature)]?.placement;
    return target?.kind === 'main' ? 'main' : target?.side || 'auto';
  } };
};
