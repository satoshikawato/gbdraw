// The identity by which a random-walk step names a feature
// (tests/web/helpers/live-random-walk.mjs). The walk records it and a replay
// finds the feature by it after the fixture has booted again, so it must not
// change between boots: a feature's `id` holds its record key, which for a
// Linear record is a UUID made at each boot (OV-2236). The record's position
// in the displayed Result (`record_idx`) and the biological feature ID (the
// key Python gives the feature in its record) do not change. The page imports
// this module (`/tests/web/helpers/live-random-walk-page.mjs`), so every
// reader of the identity uses this one derivation.

/**
 * A feature row of the displayed Result as the identity reads it.
 * @typedef {{ record_idx: unknown, biological_feature_id: unknown }} WalkFeatureRow
 */

/** @param {WalkFeatureRow} feature @returns {string} */
export const walkFeatureKey = (feature) => {
  const recordIndex = feature.record_idx;
  const featureId = String(feature.biological_feature_id ?? '');
  if (!Number.isInteger(recordIndex) || !featureId) {
    throw new Error(`a feature row without a record position or biological feature ID: ${JSON.stringify({ recordIndex, featureId })}`);
  }
  return `record ${recordIndex}, feature ${featureId}`;
};

/**
 * The feature of `features` that `key` names, or undefined.
 * @template {WalkFeatureRow} T
 * @param {readonly T[]} features
 * @param {string} key
 * @returns {T | undefined}
 */
export const findWalkFeature = (features, key) => features.find((feature) => walkFeatureKey(feature) === key);
