import { canonicalFeaturePlacements } from '../../services/feature-placement.js';

const identity = (feature) => [feature?.record_key, feature?.biological_feature_id];
const keyFor = (feature) => JSON.stringify(identity(feature));

export const createFeaturePlacementActions = ({ state, history, getCommittedRequest, isCurrentFeature }) => {
  const targetsFor = (feature) => {
    if (!getCommittedRequest()) return [];
    const items = state.featureCatalog.value?.items || [];
    const resultIndex = items.findIndex((item) => item.recordKeys.includes(feature?.record_key));
    const recordIndex = items[resultIndex]?.recordKeys.indexOf(feature?.record_key) ?? -1;
    const geometry = state.trackSlotResolvedGeometry.value;
    if (recordIndex < 0 || geometry?.mode !== state.mode.value) return [];
    return geometry.records?.find((record) => record.resultIndex === resultIndex
      && record.recordIndex === recordIndex)?.featurePlacementTargets || [];
  };
  const choices = (features) => {
    const sides = state.mode.value === 'circular' ? ['outward', 'inward'] : ['above', 'below'];
    return ['auto', 'main', ...sides].map((value) => {
      const enabled = features.length > 0 && features.every((feature) => {
        if (identity(feature).some((part) => !part) || !isCurrentFeature(feature)) return false;
        if (value === 'auto') return true;
        return targetsFor(feature).some((target) => value === 'main' ? target.kind === 'main' : target.side === value);
      });
      return { value, enabled, label: value === 'auto' ? 'Auto' : value === 'main' ? 'Main'
        : `${value[0].toUpperCase()}${value.slice(1)} lane 1`,
      reason: enabled ? '' : 'Unavailable in the resolved feature slot. Generate after changing the layout.' };
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
