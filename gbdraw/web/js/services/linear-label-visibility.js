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

/** Each current visibility field, its older boolean field, and its label. */
const LEGACY_VISIBILITY_FIELDS = Object.freeze([
  ['linear_accession_visibility', 'linear_show_accession', 'Linear Accession visibility'],
  ['linear_length_visibility', 'linear_show_length', 'Linear Length / Coordinates visibility']
]);

/**
 * The older boolean fields read as the current visibility (`show`, `hide`)
 * where the current field is not stored; with it, the older one is dropped.
 * A value that is no boolean keeps its older name (not read; Load lists it).
 * @param {Record<string, unknown>} source
 * @returns {Record<string, unknown>}
 */
export const renameLegacyLinearLabelVisibility = (source) => {
  if (!LEGACY_VISIBILITY_FIELDS.some(([, legacyKey]) => hasOwn(source, legacyKey))) return source;
  const migrated = { ...source };
  for (const [currentKey, legacyKey] of LEGACY_VISIBILITY_FIELDS) {
    if (!hasOwn(migrated, legacyKey) || typeof migrated[legacyKey] !== 'boolean') continue;
    if (!hasOwn(migrated, currentKey)) migrated[currentKey] = migrated[legacyKey] ? 'show' : 'hide';
    delete migrated[legacyKey];
  }
  return migrated;
};

/**
 * A Session 27-39 or older configuration's visibility: the older boolean
 * fields renamed, and Show where neither is stored. A stored current value is
 * kept as stored (the settings' readers read it, and Save requires one of the
 * modes). An older value that is no boolean refuses the file, except on Load
 * (`load`), which keeps it under its older name for the settings' readers to
 * list as not read, and shows the field.
 * @param {Record<string, unknown>} [source]
 * @param {{ load?: boolean }} [options]
 * @returns {Record<string, unknown>}
 */
export const migrateLegacyLinearLabelVisibility = (source = {}, { load = false } = {}) => {
  if (!source || typeof source !== 'object' || Array.isArray(source)) return source;
  const migrated = { ...renameLegacyLinearLabelVisibility(source) };
  for (const [currentKey, legacyKey, label] of LEGACY_VISIBILITY_FIELDS) {
    if (hasOwn(migrated, legacyKey) && !load) throw new Error(`${label} legacy value must be a boolean.`);
    if (!hasOwn(migrated, currentKey)) migrated[currentKey] = 'show';
  }
  return migrated;
};
