// @ts-check
import { CANONICAL_REQUEST_DESCRIPTION } from '../contracts/request-shape.generated.js';
import { readSpelling } from '../utils/contract-shape.js';
import { diagnosticError } from '../utils/error-normalization.js';
/** @import { CanonicalLosslessCollinearityParameters, CanonicalProteinComparisonSettings } from '../contracts/request.generated.js' */
/** @import { FieldReader } from '../utils/field-reading.js' */
/** @import { LiteralShape, ObjectShape } from '../utils/contract-shape.js' */

/**
 * The five comparison settings whose other spellings Python reads
 * (`read_spelling`); the Web reads them from the generated description, so it
 * reads exactly what Python reads (OV-776).
 * @typedef {'collinearityUnitMode' | 'collinearityAnchorMode' | 'collinearitySearchScope'
 *   | 'collinearityColorMode' | 'orthogroupMembershipMode'} ComparisonSettingName
 */
const COMPARISON_SETTINGS = /** @type {ObjectShape} */ (CANONICAL_REQUEST_DESCRIPTION.types.CanonicalProteinComparisonSettings).fields;
/**
 * Each of the five is a literal with spellings (tests/web/current-option-values.test.mjs).
 * @param {ComparisonSettingName} name
 */
const settingShape = (name) => /** @type {LiteralShape} */ (COMPARISON_SETTINGS[name]);

/**
 * The name of setting `name` that `value` spells, or null.
 * @template {ComparisonSettingName} K
 * @param {K} name
 * @param {unknown} value
 * @returns {CanonicalProteinComparisonSettings[K] | null}
 */
export const readComparisonSetting = (name, value) => (
  // `readSpelling` returns one of the setting's values, which its typedef lists.
  /** @type {CanonicalProteinComparisonSettings[K] | null} */ (readSpelling(value, settingShape(name)))
);

/**
 * The name of setting `name` that `value` spells; a value Python does not read
 * is refused, naming the setting by `label`.
 * @template {ComparisonSettingName} K
 * @param {K} name
 * @param {unknown} value
 * @param {string} label
 * @returns {CanonicalProteinComparisonSettings[K]}
 */
export const requireComparisonSetting = (name, value, label) => {
  const read = readComparisonSetting(name, value);
  if (read === null) throw new Error(`${label} must be one of: ${settingShape(name).values.join(', ')}.`);
  return read;
};

/**
 * Each reader also says whether the stored value was lost (BA-D25): a value
 * Python does not read is the setting's default (its blank spelling).
 * @template {ComparisonSettingName} K
 * @param {K} name
 * @returns {FieldReader<CanonicalProteinComparisonSettings[K]>}
 */
const comparisonSettingReader = (name) => {
  // Each of the five has a blank spelling, so an unset value reads as a name.
  const unset = /** @type {CanonicalProteinComparisonSettings[K]} */ (readComparisonSetting(name, null));
  return (value) => {
    const read = readComparisonSetting(name, value);
    return read === null ? { value: unset, lost: true } : { value: read, lost: false };
  };
};

export const readCollinearUnitMode = comparisonSettingReader('collinearityUnitMode');
export const readCollinearAnchorMode = comparisonSettingReader('collinearityAnchorMode');
export const readCollinearSearchScope = comparisonSettingReader('collinearitySearchScope');
export const readCollinearColorMode = comparisonSettingReader('collinearityColorMode');
export const readOrthogroupMembershipMode = comparisonSettingReader('orthogroupMembershipMode');

/** @typedef {CanonicalLosslessCollinearityParameters['mergeOrientation']} CollinearMergeOrientation */
// The merge orientation of the Collinear parameters, a literal with spellings
// ("either" reads as "strand", OV-2301).
const MERGE_ORIENTATION = /** @type {LiteralShape} */ (
  /** @type {ObjectShape} */ (CANONICAL_REQUEST_DESCRIPTION.types.CanonicalLosslessCollinearityParameters).fields.mergeOrientation
);
/**
 * The merge orientation `value` spells, or null.
 * @param {unknown} value
 * @returns {CollinearMergeOrientation | null}
 */
const readMergeOrientationSpelling = (value) => (
  // `readSpelling` returns one of the literal's values, which the typedef lists.
  /** @type {CollinearMergeOrientation | null} */ (readSpelling(value, MERGE_ORIENTATION))
);
const DEFAULT_MERGE_ORIENTATION = /** @type {CollinearMergeOrientation} */ (readMergeOrientationSpelling(null));
/** @type {FieldReader<CollinearMergeOrientation>} */
export const readCollinearMergeOrientation = (value) => {
  const read = readMergeOrientationSpelling(value);
  return read === null ? { value: DEFAULT_MERGE_ORIENTATION, lost: true } : { value: read, lost: false };
};
/**
 * The merge orientation `value` spells; a value Python does not read is refused.
 * @param {unknown} value
 * @returns {CollinearMergeOrientation}
 */
export const requireCollinearMergeOrientation = (value) => {
  const read = readMergeOrientationSpelling(value);
  if (read === null) throw diagnosticError('INPUT_INVALID', { field: 'collinear_merge_orientation', reason: 'COLLINEAR_MERGE_ORIENTATION' });
  return read;
};

/** @param {unknown} value */
export const collinearGroupScopeForEvidenceScope = (value) =>
  readCollinearSearchScope(value).value === 'all' ? 'global_collinear' : 'adjacent_local';

/** @param {unknown} value */
export const normalizeGroupMetadataScope = (value) => {
  const normalized = String(value || '').trim().toLowerCase().replace(/-/g, '_');
  if (normalized === 'global_collinear' || normalized === 'all_records') return 'global_collinear';
  if (normalized === 'adjacent_local' || normalized === 'local_collinear') return 'adjacent_local';
  if (normalized === 'record_local') return 'record_local';
  return 'global';
};

/** @param {unknown} value */
export const groupMetadataScopeLabel = (value) => {
  const scope = normalizeGroupMetadataScope(value);
  if (scope === 'global_collinear') return 'Collinearity-backed global evidence';
  if (scope === 'adjacent_local') return 'Local collinear group';
  if (scope === 'record_local') return 'Record-specific similarity group';
  return 'Cross-record similarity group';
};
