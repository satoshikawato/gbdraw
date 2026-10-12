// @ts-check
export const LINEAR_LABEL_VISIBILITY_MODES = Object.freeze(['auto', 'show', 'hide']);

/**
 * @param {unknown} value
 * @param {string} key
 */
const hasOwn = (value, key) => (
  Boolean(value)
  && typeof value === 'object'
  && !Array.isArray(value)
  && Object.prototype.hasOwnProperty.call(value, key)
);

/**
 * @param {unknown} value
 * @param {string} [label]
 * @returns {string}
 */
export const requireLinearLabelVisibilityMode = (value, label = 'Linear label visibility') => {
  const normalized = String(value ?? '').trim().toLowerCase();
  if (!LINEAR_LABEL_VISIBILITY_MODES.includes(normalized)) {
    throw new Error(`${label} must be one of: ${LINEAR_LABEL_VISIBILITY_MODES.join(', ')}.`);
  }
  return normalized;
};

/**
 * @param {unknown} mode
 * @param {{ hasSharedRow?: boolean }} [options]
 * @returns {boolean}
 */
export const resolveLinearLabelVisibility = (mode, { hasSharedRow = false } = {}) => {
  const selected = requireLinearLabelVisibilityMode(mode);
  if (selected === 'show') return true;
  if (selected === 'hide') return false;
  return !hasSharedRow;
};

/**
 * @param {unknown} mode
 * @param {{ hasSharedRow?: boolean }} [options]
 */
export const describeLinearLabelVisibility = (mode, options = {}) => {
  const selected = requireLinearLabelVisibilityMode(mode);
  if (selected !== 'auto') return selected === 'show' ? 'Show' : 'Hide';
  return `Auto · ${resolveLinearLabelVisibility(selected, options) ? 'Shown' : 'Hidden'}`;
};

/**
 * @param {Record<string, unknown>} [source]
 * @returns {Record<string, unknown>}
 */
export const migrateLegacyLinearLabelVisibility = (source = {}) => {
  if (!source || typeof source !== 'object' || Array.isArray(source)) return source;
  const migrated = { ...source };
  for (const [currentKey, legacyKey, label] of [
    ['linear_accession_visibility', 'linear_show_accession', 'Linear Accession visibility'],
    ['linear_length_visibility', 'linear_show_length', 'Linear Length / Coordinates visibility']
  ]) {
    if (hasOwn(migrated, currentKey)) {
      migrated[currentKey] = requireLinearLabelVisibilityMode(migrated[currentKey], label);
    } else if (hasOwn(migrated, legacyKey)) {
      if (typeof migrated[legacyKey] !== 'boolean') {
        throw new Error(`${label} legacy value must be a boolean.`);
      }
      migrated[currentKey] = migrated[legacyKey] ? 'show' : 'hide';
    } else {
      migrated[currentKey] = 'show';
    }
    delete migrated[legacyKey];
  }
  return migrated;
};
