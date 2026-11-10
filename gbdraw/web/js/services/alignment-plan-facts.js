// @ts-check
// The committed Similarity alignment plan read against current facts. Pure: whether the plan's
// anchors still name one feature of the shown Result's catalog and of their Similarity Group (the
// check before Generate), and the placement that Reset alignment restores. The alignment owner
// (`app/similarity-alignment.js`) runs the helper, writes the plan and publishes the notices.
import { t } from '../utils/i18n.js';
import { cloneJsonData } from './json-clone.js';
import { canonicalRecordReverseComplement } from './record-display-model.js';
/**
 * @import { CanonicalAlignmentAnchorIdentity, CanonicalLinearRecordTranslation,
 *   CanonicalSimilarityAlignmentPlan } from '../contracts/request.generated.js'
 * @import { CanonicalRecordOrientation } from './record-display-model.js'
 * @import { AlignmentResetReceipt } from './session-active-config-contract.js'
 */

/** @typedef {CanonicalAlignmentAnchorIdentity} AnchorIdentity */

/**
 * A feature of the shown Result's catalog, as far as an anchor finds it.
 * @typedef {{ recordKey?: unknown, biologicalFeatureId?: unknown, sourceFeatureIndex?: unknown }} CatalogFeatureLike
 */

/**
 * Whether two anchors name the same feature: the same record and biological feature, and no
 * source index or stable SVG id that both carry and that differ.
 * @param {AnchorIdentity | null | undefined} left
 * @param {AnchorIdentity | null | undefined} right
 */
export const anchorsAgree = (left, right) => {
  if (!left || !right) return false;
  if (
    left.recordKey !== right.recordKey
    || left.biologicalFeatureId !== right.biologicalFeatureId
  ) return false;
  return /** @type {(keyof AnchorIdentity)[]} */ (['sourceFeatureIndex', 'stableFeatureSvgId']).every((field) => (
    left[field] === null
    || left[field] === undefined
    || right[field] === null
    || right[field] === undefined
    || left[field] === right[field]
  ));
};

/**
 * The one catalog feature that a plan anchor names, or null when none or
 * several match. Generate validation and Reset share this lookup.
 * @template {CatalogFeatureLike} T
 * @param {readonly T[]} features
 * @param {AnchorIdentity | null} anchor
 * @returns {T | null}
 */
const uniqueCatalogFeature = (features, anchor) => {
  if (!anchor) return null;
  const matches = features.filter((feature) => feature.recordKey === anchor.recordKey
    && feature.biologicalFeatureId === anchor.biologicalFeatureId
    && (anchor.sourceFeatureIndex === null || feature.sourceFeatureIndex === anchor.sourceFeatureIndex));
  return matches.length === 1 ? matches[0] : null;
};

/**
 * What the check before Generate makes of a committed plan.
 * - `current`: nothing to resolve. The group was regrouped by a comparison setting or left the
 *   catalog, and every anchor still resolves; the typed render validates every anchor.
 * - `stale-reference`: the reference is gone, or a target is gone from a regrouped group (target
 *   repair offers this group's members and keeps the other choices, so it needs the reference and
 *   every remaining target in the group).
 * - `resolve`: the helper checks the plan; `goneRecordKeys` are the aligned targets to repair.
 * @typedef {{ kind: 'current' }
 *   | { kind: 'stale-reference', reason: string }
 *   | { kind: 'resolve', goneRecordKeys: Set<string> }} PlanCatalogCheck
 */

/**
 * A plan names anchor features, not group membership (PD-OI-028, OV-341): it is stale only when
 * an anchor feature no longer resolves in the shown Result's feature catalog. A Result without a
 * catalog (an older Session) cannot be checked here; the typed render validates every anchor.
 * @param {object} options
 * @param {CanonicalSimilarityAlignmentPlan} options.plan
 * @param {readonly CatalogFeatureLike[]} options.catalogFeatures The shown Result's catalog features.
 * @param {readonly AnchorIdentity[]} options.groupAnchors The anchors of the plan's Similarity Group members now.
 * @returns {PlanCatalogCheck}
 */
export const classifyPlanAgainstCatalog = ({ plan, catalogFeatures, groupAnchors }) => {
  const resolves = (/** @type {AnchorIdentity | null} */ anchor) => (
    !catalogFeatures.length || uniqueCatalogFeature(catalogFeatures, anchor) !== null);
  if (!resolves(plan.reference)) {
    return { kind: 'stale-reference', reason: t('The alignment reference feature is no longer found in its record.') };
  }
  const goneRecordKeys = new Set(plan.records
    .filter(({ status, anchor }) => status === 'aligned' && !resolves(anchor))
    .map(({ recordKey }) => recordKey));
  const inGroup = (/** @type {AnchorIdentity | null} */ anchor) => (
    groupAnchors.filter((candidate) => anchorsAgree(candidate, anchor)).length === 1);
  const regrouped = !inGroup(plan.reference) || plan.records.some(({ recordKey, status, anchor }) => (
    status === 'aligned' && !goneRecordKeys.has(recordKey) && !inGroup(anchor)));
  if (!regrouped) return { kind: 'resolve', goneRecordKeys };
  return goneRecordKeys.size === 0 ? { kind: 'current' }
    : { kind: 'stale-reference', reason: t('An aligned target feature is no longer found in its record.') };
};

/**
 * The catalog features of the plan's anchored records, each with its anchor: the members a
 * Reset's helper request offers.
 * @template {CatalogFeatureLike} T
 * @param {CanonicalSimilarityAlignmentPlan} plan
 * @param {readonly T[]} catalogFeatures
 * @returns {(T & AnchorIdentity)[]}
 */
export const resetSourceMembers = (plan, catalogFeatures) => plan.records.flatMap(({ anchor }) => {
  if (!anchor) return [];
  const feature = uniqueCatalogFeature(catalogFeatures, anchor);
  if (!feature) throw new Error('Alignment Reset source anchor is unavailable or ambiguous.');
  return [{ ...feature, ...anchor }];
});

/** @typedef {{ recordKey: string, reverseComplement: boolean }} RecordOrientation */

/**
 * The placement Reset alignment restores, before the axes of records it turns are corrected: the
 * record translations with the reference record moved back by the receipt's shift, and each
 * record's direction (the one before the latest Align when the scope restores directions, else
 * the current one).
 * @param {object} options
 * @param {readonly CanonicalLinearRecordTranslation[]} options.translations The committed translations.
 * @param {readonly (CanonicalRecordOrientation & { recordKey: string })[]} options.records
 *   The committed request's records.
 * @param {AlignmentResetReceipt | null} options.receipt The validated receipt of the latest Align.
 * @param {string} options.scope `positions` or `positions-and-directions`.
 * @returns {{ translations: CanonicalLinearRecordTranslation[], orientations: RecordOrientation[],
 *   directionsChange: boolean }}
 */
export const restoredPlacement = ({ translations, records, receipt, scope }) => {
  const restored = /** @type {CanonicalLinearRecordTranslation[]} */ (cloneJsonData(translations));
  if (receipt?.referenceDeltaX) {
    const delta = receipt.referenceDeltaX;
    const reference = restored.find(({ recordKey }) => recordKey === delta.recordKey);
    if (!reference) throw new Error('Alignment reset reference placement is unavailable.');
    reference.x -= delta.deltaX;
  }
  const restore = new Map(scope === 'positions-and-directions'
    // The caller offers this scope only when the receipt is set.
    ? /** @type {AlignmentResetReceipt} */ (receipt).directions.map((delta) => [delta.recordKey, delta.before]) : []);
  /** @type {RecordOrientation[]} */
  const orientations = records.map((record) => ({ recordKey: record.recordKey,
    reverseComplement: restore.has(record.recordKey) ? Boolean(restore.get(record.recordKey))
      : canonicalRecordReverseComplement(record) }));
  const directionsChange = orientations.some((direction, index) => direction.reverseComplement
    !== canonicalRecordReverseComplement(records[index]));
  return { translations: restored, orientations, directionsChange };
};

/**
 * The helper's projection of the restored directions, as far as Reset reads it.
 * @typedef {object} RestoredAxisProjection
 * @property {Record<string, boolean>} geometryOrientations
 * @property {readonly { recordKey: string, beforeAxisY: number, variants: readonly { axisY: number }[] }[]} records
 */

/**
 * The translations Reset writes once the helper has projected the restored directions: each
 * record's axis stays where it was drawn before.
 * @param {object} options
 * @param {readonly CanonicalLinearRecordTranslation[]} options.translations From `restoredPlacement`.
 * @param {Record<string, boolean>} options.vector Each record's restored direction, by record key.
 * @param {RestoredAxisProjection | null} options.projection The helper's validated projection of `vector`.
 * @returns {CanonicalLinearRecordTranslation[]}
 */
export const restoredAxisTranslations = ({ translations, vector, projection }) => {
  if (!projection || Object.entries(vector).some(([recordKey, reverseComplement]) => (
    projection.geometryOrientations[recordKey] !== reverseComplement))) {
    throw new Error('Alignment Reset geometry is not validated.');
  }
  return translations.map((translation) => {
    // The helper's validation checked that the projection covers the request's records.
    const fact = /** @type {RestoredAxisProjection['records'][number]} */ (
      projection.records.find((record) => record.recordKey === translation.recordKey));
    return { ...translation, y: translation.y + (fact.beforeAxisY - fact.variants[Number(vector[translation.recordKey])].axisY) };
  });
};
