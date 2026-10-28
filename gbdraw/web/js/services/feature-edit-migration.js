// @ts-check
/** @import { FeaturePlacementDraft } from './feature-placement.js' */
import { canonicalFeaturePlacements, featureIdentityKey } from './feature-placement.js';
import { diagnosticError } from '../utils/error-normalization.js';
import { tn } from '../utils/i18n.js';

// Session 44 and older kept per-feature edits in four maps keyed by rendered
// SVG ID (`RENDERED_ID_FEATURE_EDIT_FIELDS`). Python alone converts them, once
// at Load (gbdraw/session_io.py); this module keeps the placement migration and
// the field list the Session authority reads.

// A migrated row's key: [scope, recordKey, biologicalFeatureId].
/**
 * @param {unknown} scope
 * @param {string} recordKey
 * @param {string} biologicalFeatureId
 * @returns {string}
 */
const scopedIdentityKey = (scope, recordKey, biologicalFeatureId) => (
  (scope === 'circular' || scope === 'linear') && featureIdentityKey(recordKey, biologicalFeatureId)
    ? JSON.stringify([scope, recordKey, biologicalFeatureId])
    : ''
);
/** @param {Record<string, any>} row */
const scopedIdentityKeyOf = (row) => scopedIdentityKey(row?.scope, row?.recordKey, row?.biologicalFeatureId);
export const RENDERED_ID_FEATURE_EDIT_FIELDS = Object.freeze([
  'featureVisibilityOverrides',
  'labelVisibilityOverrides',
  'labelTextFeatureOverrides',
  'labelTextFeatureOverrideSources'
]);

/**
 * @param {unknown} value
 * @returns {value is Record<string, any>}
 */
const isObject = (value) => value !== null && typeof value === 'object' && !Array.isArray(value);

/** @param {number} count */
export const ANNOTATION_TARGET_MIGRATION_NOTICE = (count) => tn(count,
  '{count} annotation from an older Session named a feature by hash=. It now names that feature by its source, so it stays on the feature when the crop or orientation changes.',
  '{count} annotations from an older Session named a feature by hash=. Each now names that feature by its source, so it stays on the feature when the crop or orientation changes.');

const LANE_MODES = new Map([['outward', 'circular'], ['inward', 'circular'], ['above', 'linear'], ['below', 'linear']]);

/**
 * The mode-scoped Feature placement drafts of a Session 41-44 draft keyed by
 * [recordKey, biologicalFeatureId]. Such a row reached every request with its
 * record key, so a Main row is kept for both modes and a lane row for the mode
 * of its side. A row whose key does not encode it is kept as is, and the draft
 * check rejects it.
 * @param {any} placements
 * @returns {any} The Feature placement draft, or the unvalidated input when it is not a map.
 */
export const migrateSessionFeaturePlacements = (placements) => {
  if (!isObject(placements)) return placements;
  // Entries, so a saved key such as `__proto__` stays a row of the draft.
  /** @type {Array<[string, any]>} */
  const migrated = [];
  Object.entries(placements).forEach(([key, row]) => {
    const target = row?.placement;
    const modes = key === JSON.stringify([row?.recordKey, row?.biologicalFeatureId])
      ? target?.kind === 'main' ? ['circular', 'linear'] : [LANE_MODES.get(target?.side)].filter(Boolean) : [];
    if (modes.length === 0) migrated.push([key, row]);
    modes.forEach((scope) => {
      const scoped = { scope, ...row };
      migrated.push([scopedIdentityKeyOf(scoped) || key, scoped]);
    });
  });
  return Object.fromEntries(migrated);
};

// A Session 41-44 draft's placement rows as `migrateSessionFeaturePlacements`
// names their mode (`[scope, recordKey, biologicalFeatureId]`); each is checked
// as a row of its mode. The Session 46 split moves each row into its mode's
// drawing.
/** @param {unknown} placements */
export const validateScopedFeaturePlacements = (placements) => {
  if (!isObject(placements)) throw diagnosticError('INPUT_INVALID', { field: 'schema', reason: 'FIELDS' });
  for (const [key, row] of Object.entries(/** @type {Record<string, any>} */ (placements))) {
    const { scope, ...rest } = isObject(row) ? row : {};
    if (!['circular', 'linear'].includes(scope) || key !== scopedIdentityKeyOf(row)) {
      throw diagnosticError('INPUT_INVALID', { field: 'schema', reason: 'FIELDS' });
    }
    canonicalFeaturePlacements(/** @type {FeaturePlacementDraft} */ ({
      [featureIdentityKey(rest.recordKey, rest.biologicalFeatureId)]: rest
    }), scope);
  }
};
