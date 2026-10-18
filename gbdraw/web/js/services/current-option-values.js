// @ts-check
import { DEFAULT_FEATURE_RENDERINGS, normalizeFeatureRendering, normalizedPersistedArrowGeometryState } from '../utils/feature-rendering.js';
import { renameLegacyLinearLabelVisibility, requireLinearLabelVisibilityMode } from './linear-label-visibility.js';
import { requireComparisonSetting } from './losat-normalization.js';

/**
 * @param {unknown} value
 * @param {string} fallback
 */
const normalizedValue = (value, fallback) => {
  const normalized = String(value ?? '').trim().toLowerCase();
  return normalized || fallback;
};

/**
 * @param {unknown} value
 * @returns {value is Record<string, unknown>}
 */
const isPlainObject = (value) => (
  Boolean(value) && typeof value === 'object' && !Array.isArray(value)
);

/**
 * @param {unknown} value
 * @param {string} key
 */
const hasOwn = (value, key) => (
  isPlainObject(value) && Object.prototype.hasOwnProperty.call(value, key)
);

/**
 * @template T
 * @param {T} source
 * @param {string} obsoleteName
 * @param {string} currentName
 * @returns {T}
 */
const migrateOwnField = (source, obsoleteName, currentName) => {
  if (!hasOwn(source, obsoleteName)) return source;
  /** @type {Record<string, unknown>} */
  // `hasOwn` above holds only for a plain object.
  const migrated = { .../** @type {Record<string, unknown>} */ (source) };
  if (!hasOwn(migrated, currentName)) {
    migrated[currentName] = migrated[obsoleteName];
  }
  delete migrated[obsoleteName];
  // The copy keeps the shape of `source`, with the obsolete name moved to the current one.
  return /** @type {T} */ (migrated);
};

/**
 * @template {string} T
 * @param {unknown} value
 * @param {string} fallback
 * @param {readonly T[]} supported
 * @param {string} label
 * @returns {T} One of the `supported` names.
 */
const requireCurrentValue = (value, fallback, supported, label) => {
  const normalized = normalizedValue(value, fallback);
  if (!(/** @type {readonly string[]} */ (supported)).includes(normalized)) {
    throw new Error(`${label} must be one of: ${supported.join(', ')}.`);
  }
  // `supported` holds it: checked above.
  return /** @type {T} */ (normalized);
};

// The values of each option a stored configuration names; Save refuses
// another value (`requireCurrent*`), and Load reads it as the default with a
// dialog row (`config-field-readers.js`). The older names each read as.
export const CIRCULAR_MULTI_RECORD_SIZE_MODES = /** @type {const} */ (['auto', 'linear', 'equal']);
export const LINEAR_TRACK_LAYOUTS = /** @type {const} */ (['above', 'middle', 'below']);
export const LINEAR_LABEL_PLACEMENTS = /** @type {const} */ (['auto', 'above_feature']);
export const PROTEIN_BLASTP_MODES = /** @type {const} */ (['pairwise', 'orthogroup', 'collinear']);
export const CIRCULAR_MULTI_RECORD_SIZE_MODE_ALIASES = Object.freeze({ sqrt: /** @type {const} */ ('auto') });
export const LINEAR_TRACK_LAYOUT_ALIASES = Object.freeze({ spreadout: /** @type {const} */ ('above'), tuckin: /** @type {const} */ ('below') });
export const LINEAR_LABEL_PLACEMENT_ALIASES = Object.freeze({ on_feature: /** @type {const} */ ('above_feature') });

/** @param {unknown} value */
export const requireCurrentCircularMultiRecordSizeMode = (value) => (
  requireCurrentValue(value, 'auto', CIRCULAR_MULTI_RECORD_SIZE_MODES, 'Circular multi-record size mode')
);

/** @param {unknown} value */
export const requireCurrentLinearTrackLayout = (value) => (
  requireCurrentValue(value, 'middle', LINEAR_TRACK_LAYOUTS, 'Linear track layout')
);

/** @param {unknown} value */
export const requireCurrentLinearLabelPlacement = (value) => (
  requireCurrentValue(value, 'auto', LINEAR_LABEL_PLACEMENTS, 'Linear label placement')
);

/** The Pairwise match styles (the Load reader and Generate read this list). */
export const PAIRWISE_MATCH_STYLES = /** @type {const} */ (['ribbon', 'curve']);

/**
 * @param {unknown} value
 * @param {unknown} fallback
 */
export const normalizeCurrentPairwiseMatchStyle = (value, fallback) => {
  const supported = /** @type {readonly string[]} */ (PAIRWISE_MATCH_STYLES);
  const normalizedFallback = String(fallback ?? '').trim().toLowerCase();
  if (!supported.includes(normalizedFallback)) {
    throw new Error('Pairwise match style fallback must be one of: ribbon, curve.');
  }
  const normalized = String(value ?? '').trim().toLowerCase();
  return supported.includes(normalized) ? normalized : normalizedFallback;
};

/** @param {unknown} value */
export const requireCurrentProteinBlastpMode = (value) => (
  requireCurrentValue(value, 'orthogroup', PROTEIN_BLASTP_MODES, 'Protein BLASTP mode')
);

/** @param {unknown} value */
export const requireCurrentCollinearSearchScope = (value) => (
  requireComparisonSetting('collinearitySearchScope', value, 'Collinear search scope')
);

// Older saved pipelines omit this field; fresh Web state explicitly supplies false.
/** @param {unknown} [value] */
export const requireCurrentCollinearInferOrthogroups = (value = true) => {
  if (typeof value !== 'boolean') {
    throw new Error('Collinear orthogroup inference must be a boolean.');
  }
  return value;
};

/** @param {unknown} value */
export const requireCurrentProteinBlastpCandidateLimit = (value) => {
  if (value === null || value === undefined || value === '') return null;
  const numeric = Number(value);
  if (!Number.isInteger(numeric) || numeric <= 0) {
    throw new Error('Protein BLASTP Max target seqs must be a positive integer or None.');
  }
  return numeric;
};

/**
 * @template {number | null} F
 * @param {unknown} value
 * @param {F} fallback
 * @param {number} minimum
 * @param {string} label
 * @returns {F | number}
 */
const requireCurrentInteger = (value, fallback, minimum, label) => {
  if (value === null || value === undefined || value === '') return fallback;
  const numeric = Number(value);
  if (!Number.isInteger(numeric) || numeric < minimum) {
    const domain = minimum === 0 ? 'a non-negative integer' : 'a positive integer';
    throw new Error(`${label} must be ${domain}.`);
  }
  return numeric;
};

/** @param {unknown} value */
export const requireCurrentProteinBlastpMaxHits = (value) => (
  requireCurrentInteger(value, 5, 1, 'Protein BLASTP Pairwise max hits')
);

/** @param {unknown} value */
export const requireCurrentOrthogroupMemberMaxHits = (value) => (
  requireCurrentInteger(value, null, 1, 'Protein BLASTP member hits per protein')
);

/** @param {unknown} value */
export const requireCurrentCollinearMinAnchors = (value) => (
  requireCurrentInteger(value, 1, 1, 'Collinear minimum anchors')
);

/** @param {unknown} value */
export const requireCurrentCollinearMaxUnitGap = (value) => (
  requireCurrentInteger(value, 0, 0, 'Collinear maximum unit gap')
);

/** @param {unknown} value */
export const requireCurrentCollinearMaxDiagonalDrift = (value) => (
  requireCurrentInteger(value, 0, 0, 'Collinear maximum diagonal drift')
);

/** @param {unknown} value */
export const requireCurrentCollinearMaxConflicts = (value) => (
  requireCurrentInteger(value, 1, 0, 'Collinear maximum conflicts in merge gap')
);

const DEFAULT_COLLINEAR_MAX_PARALOG_LINKS = 2;

/** @param {unknown} value */
export const requireCurrentCollinearMaxParalogLinks = (value) => (
  requireCurrentInteger(
    value, DEFAULT_COLLINEAR_MAX_PARALOG_LINKS, 1, 'Collinear maximum paralog links per orthogroup'
  )
);

/**
 * The paralog links per group a Generate draws with. Only Collinear mode reads
 * the stored setting; Similarity groups draws with the default (OV-1128).
 * @param {boolean} collinear Whether the Generate runs protein BLASTP in Collinear mode.
 * @param {unknown} value The stored setting.
 */
export const drawnCollinearMaxParalogLinks = (collinear, value) => (
  collinear ? requireCurrentCollinearMaxParalogLinks(value) : DEFAULT_COLLINEAR_MAX_PARALOG_LINKS
);

/** @param {unknown} value */
export const requireCurrentCollinearUnitMode = (value) => (
  requireComparisonSetting('collinearityUnitMode', value, 'Collinear unit mode')
);

/** @param {unknown} value */
export const requireCurrentCollinearAnchorMode = (value) => (
  requireComparisonSetting('collinearityAnchorMode', value, 'Collinear anchor mode')
);

/** @param {unknown} value */
export const requireCurrentCollinearMergeOrientation = (value) => (
  requireCurrentValue(
    value,
    'either',
    ['strand', 'order', 'either'],
    'Collinear merge orientation'
  )
);

/** @param {unknown} value */
export const requireCurrentCollinearColorMode = (value) => (
  requireComparisonSetting('collinearityColorMode', value, 'Collinear color mode')
);

/** @param {unknown} value */
export const requireCurrentOrthogroupMembershipMode = (value) => (
  requireComparisonSetting('orthogroupMembershipMode', value, 'Protein BLASTP orthogroup membership mode')
);

/** @param {unknown} value */
/**
 * @param {unknown} value
 * @param {string} fallback
 * @param {Readonly<Record<string, string>>} aliases
 */
const aliasedValue = (value, fallback, aliases) => {
  const normalized = normalizedValue(value, fallback);
  return Object.hasOwn(aliases, normalized) ? aliases[normalized] : normalized;
};

/** @param {unknown} value */
export const migratePersistedCircularMultiRecordSizeMode = (value) => requireCurrentCircularMultiRecordSizeMode(
  aliasedValue(value, 'auto', CIRCULAR_MULTI_RECORD_SIZE_MODE_ALIASES)
);

/** @param {unknown} value */
export const migratePersistedLinearTrackLayout = (value) => requireCurrentLinearTrackLayout(
  aliasedValue(value, 'middle', LINEAR_TRACK_LAYOUT_ALIASES)
);

/** @param {unknown} value */
export const migratePersistedLinearLabelPlacement = (value) => requireCurrentLinearLabelPlacement(
  aliasedValue(value, 'auto', LINEAR_LABEL_PLACEMENT_ALIASES)
);

/**
 * The sections of a Web config document that callers read after the migration. A document without a section, or
 * with a non-object one, holds `undefined` or that value there, so callers check a section with `isPlainObject`
 * before they write to it. The type states the shape callers rely on while `config.js` types its documents.
 * @typedef {{ form: Record<string, unknown>, adv: Record<string, unknown> }} PersistedWebSections
 */

/**
 * @template {{ adv?: unknown, losat?: unknown }} T
 * @param {T} config
 * @returns {T & PersistedWebSections}
 */
export const migratePersistedWebStateFieldNames = (config) => {
  if (!isPlainObject(config)) return /** @type {T & PersistedWebSections} */ (config);

  let adv = migrateOwnField(
    config.adv,
    'depth_tick_interval',
    'depth_large_tick_interval'
  );
  if (isPlainObject(adv)) adv = renameLegacyLinearLabelVisibility(adv);
  if (isPlainObject(adv) && Array.isArray(adv.depth_tracks)) {
    adv = {
      ...adv,
      depth_tracks: adv.depth_tracks.map((track) => (
        migrateOwnField(track, 'tick_interval', 'large_tick_interval')
      ))
    };
  }

  let losat = config.losat;
  if (isPlainObject(losat) && isPlainObject(losat.blastp)) {
    const blastp = migrateOwnField(
      losat.blastp,
      'collinearMaxGeneGap',
      'collinearMaxUnitGap'
    );
    if (blastp !== losat.blastp) {
      losat = { ...losat, blastp };
    }
  }

  return /** @type {T & PersistedWebSections} */ ({
    ...config,
    ...(adv === undefined ? {} : { adv }),
    ...(losat === undefined ? {} : { losat })
  });
};

/**
 * @template {{ adv?: { depth_tracks?: unknown } | null, losat?: { blastp?: unknown } | null } | null | undefined} T
 * @param {T} config
 * @returns {T}
 */
export const requireCurrentWebStateFieldNames = (config) => {
  const adv = isPlainObject(config) ? config.adv : null;
  for (const field of ['linear_show_accession', 'linear_show_length']) {
    if (hasOwn(adv, field)) {
      throw new Error(`Web state field adv.${field} is obsolete; use the selected visibility mode.`);
    }
  }
  if (hasOwn(adv, 'depth_tick_interval')) {
    throw new Error(
      'Web state field adv.depth_tick_interval is obsolete; use adv.depth_large_tick_interval.'
    );
  }
  if (Array.isArray(adv?.depth_tracks)) {
    adv.depth_tracks.forEach((track, index) => {
      if (hasOwn(track, 'tick_interval')) {
        throw new Error(
          `Web state field adv.depth_tracks[${index}].tick_interval is obsolete; ` +
          'use large_tick_interval.'
        );
      }
    });
  }
  const blastp = isPlainObject(config?.losat) ? config.losat.blastp : null;
  if (hasOwn(blastp, 'collinearMaxGeneGap')) {
    throw new Error(
      'Web state field losat.blastp.collinearMaxGeneGap is obsolete; ' +
      'use losat.blastp.collinearMaxUnitGap.'
    );
  }
  return config;
};

// The merge of stored settings installs text or null over a text setting; a
// value of another kind is not used, and Load lists it (`mergeStoredSettings`).
/** @param {unknown} value */
const installedOverText = (value) => typeof value === 'string' || value === null;

// The invariants of a configuration the app wrote (BA-OD96): no retired field
// name, no value that no current option names (the Linear label visibilities
// included), no default feature type's shape that names no shape, and no
// arrow ratio out of range. Only a value the merge installs is checked
// (BA-OD100). Save, History (`applyConfigData`) and the Load rollback hold
// them; Load reads a file's values with the reader table, which lists what it
// cannot restore (`config-field-readers.js`, E1a), and refuses none of them.
/**
 * @param {unknown} config
 */
export const requireCurrentConfigFields = (config) => {
  if (!isPlainObject(config)) return;
  requireCurrentWebStateFieldNames(config);
  const { form, adv } = config;
  if (isPlainObject(form) && hasOwn(form, 'linear_track_layout')) requireCurrentLinearTrackLayout(form.linear_track_layout);
  if (!isPlainObject(adv)) return;
  if (hasOwn(adv, 'label_placement')) requireCurrentLinearLabelPlacement(adv.label_placement);
  if (hasOwn(adv, 'multi_record_size_mode')) requireCurrentCircularMultiRecordSizeMode(adv.multi_record_size_mode);
  if (installedOverText(adv.linear_accession_visibility)) {
    requireLinearLabelVisibilityMode(adv.linear_accession_visibility, 'Linear Accession visibility');
  }
  if (installedOverText(adv.linear_length_visibility)) {
    requireLinearLabelVisibilityMode(adv.linear_length_visibility, 'Linear Length / Coordinates visibility');
  }
  // The arrow ratios and the default types' shapes `applyConfigData` installs.
  normalizedPersistedArrowGeometryState(adv);
  if (isPlainObject(adv.feature_shapes)) {
    for (const [featureType, shape] of Object.entries(adv.feature_shapes)) {
      if (Object.hasOwn(DEFAULT_FEATURE_RENDERINGS, featureType) && installedOverText(shape)) normalizeFeatureRendering(shape);
    }
  }
};
