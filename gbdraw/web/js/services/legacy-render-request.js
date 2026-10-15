// @ts-check
// A stored render request of schema 1 or 2 (Sessions 31-33) read as the
// schema-5 request that means the same (PR-D). The reader of stored requests
// then reads schema 5 and later only. Python's decoder owns the tables
// (`describe_legacy_request_tables` in gbdraw/session_request_codec.py,
// generated into LEGACY_REQUEST_TABLES); where the Web read a form more
// leniently than Python (OV-793, OV-795, OV-810) or with another meaning
// (OV-798), the lift keeps the Web's reading.
import { LEGACY_REQUEST_TABLES } from '../contracts/schema.generated.js';
import { diagnosticError } from '../utils/error-normalization.js';
import { loadRepair } from '../utils/load-repairs.js';
import { classifyOptionalPositiveNumber, DECIMAL_NUMBER_PATTERN } from '../utils/optional-positive-number.js';
import { migrateLegacyCircularTrackSlot, readLegacyCircularTrackSlotSpec } from './circular-track-slot-model.js';
import { parseOptionalPixel } from './track-slot-validation.js';

/** @import { LoadRepair, Repaired } from '../utils/load-repairs.js' */
/** @typedef {Record<string, unknown>} JsonObject */

const LIFTED_SCHEMA = 5;

/**
 * @param {unknown} value
 * @returns {value is JsonObject}
 */
const isObject = (value) => Boolean(value) && typeof value === 'object' && !Array.isArray(value);

/** @param {unknown} value */
const isUnset = (value) => value === null || value === undefined;

/** @param {unknown} value */
const foldedText = (value) => String(value ?? '').trim().toLowerCase();

/** @param {unknown} value A stored value as a repair row names it. */
const excerpt = (value) => (typeof value === 'string' ? value : JSON.stringify(value) ?? String(value));

/**
 * The value `table` maps a legacy text to (case and outer spaces ignored, as
 * the Web read it), or the value itself.
 * @param {unknown} value
 * @param {Readonly<Record<string, string>>} table
 */
const renamedToken = (value, table) => (
  typeof value === 'string' && Object.hasOwn(table, foldedText(value)) ? table[foldedText(value)] : value
);

/**
 * Explicit grouping (A2): a Linear request draws one figure; a Circular one a
 * grid when it has a layout or several records.
 * @param {JsonObject} request
 */
const inferredGrouping = (request) => {
  if (request.mode !== 'circular') return 'single';
  const layout = isObject(request.layout) ? request.layout : {};
  const records = Array.isArray(request.records) ? request.records : [];
  return Object.keys(layout).length > 0 || records.length > 1 ? 'grid' : 'single';
};

/**
 * The feature visibility table under its current names (A15, OV-792). A
 * current name that holds a value keeps it.
 * @param {JsonObject} options
 */
const liftFeatureVisibilityNames = (options) => {
  for (const [old, current] of Object.entries(LEGACY_REQUEST_TABLES.diagramOptionRenames)) {
    if (!Object.hasOwn(options, old)) continue;
    const stored = options[old];
    delete options[old];
    if (isUnset(options[current])) options[current] = stored;
  }
};

/**
 * The sparse defaults a schema 1-2 request omitted (A4, A5, A6), written out.
 * @param {JsonObject} options
 * @param {'circular' | 'linear'} mode
 */
const liftSparseDefaults = (options, mode) => {
  for (const [key, value] of Object.entries(LEGACY_REQUEST_TABLES.comparisonDefaults)) {
    if (isUnset(options[key])) options[key] = value;
  }
  if (isUnset(options.selectedFeaturesSet)) options.selectedFeaturesSet = [...LEGACY_REQUEST_TABLES.featureTypes];
  const shapes = isObject(options.featureShapes) ? options.featureShapes : {};
  for (const [type, shape] of Object.entries(LEGACY_REQUEST_TABLES.featureShapes)) {
    if (isUnset(shapes[type])) shapes[type] = shape;
  }
  options.featureShapes = shapes;
  if (!isUnset(options.config)) return;
  const explicit = options.configOverrides ?? {};
  if (!isObject(explicit)) return;
  const overrides = Object.fromEntries(Object.entries(explicit).filter(([, value]) => !isUnset(value)));
  /** @type {Readonly<Record<string, readonly string[]>>} */
  const flatPaths = LEGACY_REQUEST_TABLES.configOverridePaths;
  for (const [path, value] of Object.entries(LEGACY_REQUEST_TABLES.configOverrides[mode])) {
    // A flat name of the same leaf holds its value (OV-799).
    const named = [path, ...Object.keys(flatPaths).filter((name) => flatPaths[name].includes(path))];
    if (!named.some((name) => Object.hasOwn(overrides, name))) overrides[path] = value;
  }
  options.configOverrides = overrides;
};

/**
 * The retired names and values of configuration overrides (A7) and an
 * invalid comparison height (B17: a Linear request drew it before b1b48f6e0;
 * now unset with a row, and without one in a Circular request, which does
 * not use it: D-OD10).
 * @param {JsonObject} options
 * @param {'circular' | 'linear'} mode
 * @param {string} at the path of `options` in the file
 * @param {LoadRepair[]} repairs
 */
const liftConfigOverrides = (options, mode, at, repairs) => {
  const overrides = options.configOverrides;
  if (!isObject(overrides)) return;
  for (const [old, current] of Object.entries(LEGACY_REQUEST_TABLES.configOverrideRenames)) {
    if (!Object.hasOwn(overrides, old)) continue;
    const stored = overrides[old];
    delete overrides[old];
    if (isUnset(overrides[current])) overrides[current] = stored;
  }
  for (const key of ['comparison_height', ...LEGACY_REQUEST_TABLES.configOverridePaths.comparison_height]) {
    if (!Object.hasOwn(overrides, key) || classifyOptionalPositiveNumber(overrides[key]).status !== 'invalid') continue;
    if (mode === 'linear') {
      repairs.push(loadRepair({
        code: 'setting-replaced', path: `${at}.configOverrides.${key}`, effect: 'dropped', stored: excerpt(overrides[key])
      }));
    }
    delete overrides[key];
  }
  if (mode !== 'linear') return;
  /** @type {Readonly<Record<string, readonly string[]>>} */
  const canonicalPaths = LEGACY_REQUEST_TABLES.configOverridePaths;
  for (const [name, table] of Object.entries(LEGACY_REQUEST_TABLES.linearConfigOverrideValues)) {
    for (const key of [name, ...(canonicalPaths[name] ?? [])]) {
      if (Object.hasOwn(overrides, key)) overrides[key] = renamedToken(overrides[key], table);
    }
  }
};

/**
 * The Linear leaves of a stored full config the schema-5 reader requires in
 * their current form (A7: the two rewrites of Python's full-config
 * migration), and an invalid comparison height, which the config cannot
 * leave unset: its default is used (B17, F1), with a row in a Linear request
 * unless a valid override of the leaf was drawn instead (N3); a Circular
 * request does not use it and gets no row (D-OD10).
 * @param {JsonObject} options
 * @param {'circular' | 'linear'} mode
 * @param {string} at the path of `options` in the file
 * @param {LoadRepair[]} repairs
 */
const liftFullConfig = (options, mode, at, repairs) => {
  const config = options.config;
  if (!isObject(config)) return;
  const values = LEGACY_REQUEST_TABLES.linearConfigOverrideValues;
  const labels = isObject(config.labels) && isObject(config.labels.linear) ? config.labels.linear : null;
  const canvas = isObject(config.canvas) && isObject(config.canvas.linear) ? config.canvas.linear : null;
  if (mode === 'linear') {
    if (labels && Object.hasOwn(labels, 'placement')) labels.placement = renamedToken(labels.placement, values.label_placement);
    if (canvas && Object.hasOwn(canvas, 'track_layout')) {
      canvas.track_layout = renamedToken(canvas.track_layout, values.linear_track_layout);
    }
  }
  if (canvas && Object.hasOwn(canvas, 'comparison_height')
    && classifyOptionalPositiveNumber(canvas.comparison_height).status === 'invalid') {
    const overrides = isObject(options.configOverrides) ? options.configOverrides : {};
    const overridden = ['comparison_height', ...LEGACY_REQUEST_TABLES.configOverridePaths.comparison_height]
      .some((key) => Object.hasOwn(overrides, key) && classifyOptionalPositiveNumber(overrides[key]).status === 'valid');
    if (mode === 'linear' && !overridden) {
      repairs.push(loadRepair({
        code: 'setting-replaced', path: `${at}.config.canvas.linear.comparison_height`, effect: 'defaulted',
        stored: excerpt(canvas.comparison_height)
      }));
    }
    canvas.comparison_height = LEGACY_REQUEST_TABLES.comparisonHeightDefault;
  }
};

/**
 * A legacy pixel value: a number, text, or `{ value, unit: 'px' }`.
 * @param {unknown} value
 * @returns {number | null | undefined} undefined when it cannot be read.
 */
const legacyPixel = (value) => {
  try {
    if (isObject(value)) {
      return foldedText(value.unit) === 'px' ? parseOptionalPixel(value.value, 'spacing', { allowZero: true }) : undefined;
    }
    return parseOptionalPixel(value, 'spacing', { allowZero: true });
  } catch {
    return undefined;
  }
};

const PRIVATE_SPACING = LEGACY_REQUEST_TABLES.privateSpacingParam;
// A slot's legacy spacing, in both forms: its `spacing`, else the private
// spacing parameter Python's decoder read (N2, N3).
const LEGACY_SPACING_KEYS = ['spacing', PRIVATE_SPACING];
const GAP_NAMES = [['inner_gap_px', 'innerGapPx'], ['outer_gap_px', 'outerGapPx']];
const GAP_KEYS = new Set(GAP_NAMES.flat());

/**
 * A spacing in percent of the axis radius (`5%`), as Python reads it: the
 * factor 0.05 (`ScalarSpec.parse`).
 * @param {unknown} value
 */
const isPercent = (value) => typeof value === 'string' && value.trim().endsWith('%')
  && DECIMAL_NUMBER_PATTERN.test(value.trim().slice(0, -1).trim());

/**
 * How Python read a legacy slot spacing: in pixels (`4px`, `{ unit: 'px' }`),
 * as a multiple of the axis radius when it has no unit (`4`), a factor unit or
 * a percent sign, or not at all. The Web read a spacing without a unit as
 * pixels (OV-798).
 * @param {unknown} value
 * @returns {'unitless' | 'factor' | 'other'}
 */
const spacingUnit = (value) => {
  if (isObject(value)) return foldedText(value.unit) === 'factor' ? 'factor' : 'other';
  if (typeof value === 'number') return 'unitless';
  if (isPercent(value)) return 'factor';
  return typeof value === 'string' && DECIMAL_NUMBER_PATTERN.test(value.trim()) ? 'unitless' : 'other';
};

/**
 * A factor spacing of 0 (`{ value: 0, unit: 'factor' }`, `0%`): no gap, read
 * as 0 px with no row (N4, D-OD10). Any other factor fits no current
 * schema, which holds gaps in pixels only.
 * @param {unknown} value
 */
const isZeroFactor = (value) => spacingUnit(value) === 'factor' && (isObject(value)
  ? legacyPixel({ ...value, unit: 'px' }) === 0
  : Number(String(value).trim().slice(0, -1)) === 0);

/**
 * The row of a legacy spacing read with another meaning (OV-798): a
 * spacing without a unit is read as pixels, Python's factor form is
 * dropped. A spacing of 0 is no gap either way.
 * @param {unknown} spacing
 * @param {string} path
 * @param {string} stored
 * @returns {LoadRepair | null}
 */
const spacingRepair = (spacing, path, stored) => {
  const unit = spacingUnit(spacing);
  if (unit === 'factor') return loadRepair({ code: 'setting-replaced', path, effect: 'dropped', stored });
  const pixel = unit === 'unitless' ? legacyPixel(spacing) : undefined;
  return pixel ? loadRepair({ code: 'setting-replaced', path, effect: 'changed', stored, used: `${pixel} px` }) : null;
};

/**
 * A Circular object slot without the legacy spacing (A8): the spacing (or
 * the private spacing parameter Python's decoder read, N3) fills the gaps
 * the slot leaves unset, and its gaps are read as pixels as the Web read
 * them, text and `{ value, unit: 'px' }` included (S3). A spacing or gap
 * that cannot be read is left for the reader to refuse.
 * @param {JsonObject} slot
 * @param {string} at its path in the file
 * @param {LoadRepair[]} repairs
 * @returns {JsonObject}
 */
const liftCircularObjectSlot = (slot, at, repairs) => {
  let source = slot;
  let spacingPath = `${at}.spacing`;
  if (isObject(slot.params) && Object.hasOwn(slot.params, PRIVATE_SPACING) && isUnset(slot.spacing)) {
    const { [PRIVATE_SPACING]: spacing, ...params } = slot.params;
    source = { ...slot, spacing, params };
    spacingPath = `${at}.params.${PRIVATE_SPACING}`;
  }
  if (isZeroFactor(source.spacing)) source = { ...source, spacing: 0 };
  const gapUnset = isUnset(source.innerGapPx ?? source.inner_gap_px) || isUnset(source.outerGapPx ?? source.outer_gap_px);
  const repair = isUnset(source.spacing) || !gapUnset ? null : spacingRepair(source.spacing, spacingPath, excerpt(source.spacing));
  if (spacingUnit(source.spacing) === 'factor') source = { ...source, spacing: null };
  const migrated = /** @type {JsonObject} */ (migrateLegacyCircularTrackSlot(source));
  /** @type {JsonObject} */
  const lifted = {};
  for (const [key, value] of Object.entries(migrated)) {
    if (!GAP_KEYS.has(key)) lifted[key] = value;
  }
  for (const [legacy, current] of GAP_NAMES) {
    if (!Object.hasOwn(migrated, legacy) && !Object.hasOwn(migrated, current)) continue;
    const pixel = legacyPixel(migrated[current] ?? migrated[legacy]);
    if (pixel === undefined) return slot;
    lifted[current] = pixel;
  }
  if (repair) repairs.push(repair);
  return lifted;
};

/**
 * The gap a legacy spacing token fills: 0 for a factor of 0, none for any
 * other factor (dropped with a row), else the token as stored.
 * @param {string} value
 * @returns {string | null}
 */
const spacingGap = (value) => {
  if (isZeroFactor(value)) return '0';
  return spacingUnit(value) === 'factor' ? null : value;
};

/**
 * A Circular string slot under current tokens (A9), its legacy spacing read
 * as the object form reads it (N2), with the row of the spacing token the
 * parser consumed when it has no unit (OV-798) or is a factor (D-OD10).
 * @param {unknown} spec
 * @param {string} at its path in the file
 * @param {LoadRepair[]} repairs
 */
const liftCircularStringSlot = (spec, at, repairs) => {
  const { spec: lifted, spacing } = readLegacyCircularTrackSlotSpec(spec, LEGACY_SPACING_KEYS, spacingGap);
  const repair = spacing?.filled && !isZeroFactor(spacing.value) ? spacingRepair(spacing.value, at, spacing.stored) : null;
  if (repair) repairs.push(repair);
  return lifted;
};

/**
 * The Circular slots under current tokens (A8, A9).
 * @param {JsonObject} options
 * @param {string} optionsAt the path of `options` in the file
 * @param {LoadRepair[]} repairs
 */
const liftCircularSlots = (options, optionsAt, repairs) => {
  const tracks = options.tracks;
  if (!isObject(tracks) || !Array.isArray(tracks.circularTrackSlots)) return;
  const at = `${optionsAt}.tracks.circularTrackSlots`;
  tracks.circularTrackSlots = tracks.circularTrackSlots.map((/** @type {unknown} */ slot, /** @type {number} */ index) => (
    isObject(slot) ? liftCircularObjectSlot(slot, `${at}[${index}]`, repairs) : liftCircularStringSlot(slot, `${at}[${index}]`, repairs)
  ));
};

/**
 * A `kind=standard` collinearity payload as the lossless parameters Python
 * reads it as (A14, OV-794).
 * @param {JsonObject} comparison
 */
const liftStandardCollinearity = (comparison) => {
  const settings = comparison.settings;
  if (!isObject(settings)) return;
  const params = settings.collinearityParams;
  if (!isObject(params) || params.kind !== 'standard' || !isObject(params.parameters)) return;
  const standard = params.parameters;
  if (!LEGACY_REQUEST_TABLES.collinearityStandardParameters.every((key) => Object.hasOwn(standard, key))) return;
  /** @type {JsonObject} */
  const parameters = {};
  for (const [from, to] of Object.entries(LEGACY_REQUEST_TABLES.collinearityStandardFields)) parameters[to] = standard[from];
  parameters.mergeOrientation = LEGACY_REQUEST_TABLES.collinearityMergeOrientation;
  settings.collinearityParams = { kind: 'lossless', parameters };
  const paralogKey = LEGACY_REQUEST_TABLES.paralogLinksSetting;
  if (foldedText(comparison.mode) === 'collinear' && settings[paralogKey] === LEGACY_REQUEST_TABLES.paralogLinksDefault) {
    settings[paralogKey] = standard[LEGACY_REQUEST_TABLES.collinearityStandardParalogLinks];
  }
};

/**
 * Comparisons with explicit record indices (A12, A13) and lossless
 * collinearity parameters (A14). An index the Web read positionally becomes
 * that index: the n-th item of a kind joins records n and n + 1.
 * @param {JsonObject} request
 */
const liftComparisons = (request) => {
  if (!Array.isArray(request.comparisons)) return;
  /** @type {Map<unknown, number>} */
  const positions = new Map();
  request.comparisons = request.comparisons.map((/** @type {unknown} */ item) => {
    if (!isObject(item)) return item;
    const comparison = { ...item };
    if (comparison.kind === 'nucleotideBlast' || (comparison.kind === 'precomputedProteinComparison' && request.schema === 1)) {
      const position = positions.get(comparison.kind) ?? 0;
      positions.set(comparison.kind, position + 1);
      const index = (/** @type {unknown} */ value, /** @type {number} */ fallback) => (
        Number.isInteger(Number(value)) ? Number(value) : fallback);
      comparison.queryRecordIndex = index(comparison.queryRecordIndex, position);
      comparison.subjectRecordIndex = index(comparison.subjectRecordIndex, position + 1);
    }
    if (comparison.kind === 'generatedProteinComparison') {
      if (!Object.hasOwn(comparison, 'pairs')) comparison.pairs = [];
      liftStandardCollinearity(comparison);
    }
    return comparison;
  });
};

/**
 * A form no schema 1 or 2 writer stored, which Python's decoder refuses, as
 * the Web did: a field of a later schema (D-OD13 F3), or an output without
 * the prefix these schemas require (D-OD16 F2).
 * @param {JsonObject} request
 * @param {JsonObject} options
 */
const refuseUnwrittenForms = (request, options) => {
  if (Array.isArray(options.featureOverrides)) {
    throw diagnosticError('INPUT_INVALID', { field: 'featureOverrides', reason: 'SESSION_FORMAT' });
  }
  const layout = isObject(request.layout) ? request.layout : {};
  const records = Array.isArray(request.records) ? request.records : [];
  if ((isObject(options.output) && !Object.hasOwn(options.output, 'outputPrefix'))
    || Object.hasOwn(options, 'featurePlacements')
    || records.some((/** @type {unknown} */ record) => isObject(record) && Object.hasOwn(record, 'display'))
    || (request.mode === 'linear' && (Object.hasOwn(layout, 'recordTranslations') || Object.hasOwn(layout, 'similarityAlignment')))) {
    throw diagnosticError('INPUT_INVALID', { field: 'schema', reason: 'SESSION_FORMAT' });
  }
};

/**
 * A stored render request of schema 1 or 2 as the schema-5 request that
 * means the same; any other value unchanged. The repairs list what the
 * schema-5 request cannot carry or reads with another meaning (B17,
 * OV-798, D-OD10); a lossless lift has none. A value the reader refuses (not
 * an object, an unknown mode) is left to the reader; a form no writer
 * stored is refused here (D-OD13 F3, D-OD16 F2).
 * @template T
 * @param {T} value a stored `renderRequest` (parsed JSON)
 * @param {string} [path] its path in the file, which the rows name
 * @returns {Repaired<T | JsonObject>} `value` itself, or the lifted copy
 */
export const liftLegacyRenderRequest = (value, path = 'renderRequest') => {
  if (!isObject(value) || (value.schema !== 1 && value.schema !== 2)
    || (value.mode !== 'circular' && value.mode !== 'linear') || !isObject(value.diagramOptions)) {
    return { value, repairs: [] };
  }
  refuseUnwrittenForms(value, value.diagramOptions);
  /** @type {LoadRepair[]} */
  const repairs = [];
  const mode = value.mode;
  /** @type {JsonObject} */
  const request = JSON.parse(JSON.stringify(value));
  // The cast: `isObject(value.diagramOptions)` above; `request` is its copy.
  const options = /** @type {JsonObject} */ (request.diagramOptions);
  const optionsAt = `${path}.diagramOptions`;
  if (isObject(options.output)) delete options.output.outputPrefix;
  liftFeatureVisibilityNames(options);
  liftSparseDefaults(options, mode);
  liftConfigOverrides(options, mode, optionsAt, repairs);
  liftFullConfig(options, mode, optionsAt, repairs);
  if (mode === 'circular') {
    liftCircularSlots(options, optionsAt, repairs);
    if (isObject(request.layout) && Object.hasOwn(request.layout, 'multiRecordSizeMode')) {
      request.layout.multiRecordSizeMode = renamedToken(
        request.layout.multiRecordSizeMode, LEGACY_REQUEST_TABLES.multiRecordSizeModeRenames);
    }
  } else if (request.schema === 1) {
    // Schema 1 had no Linear layout; the Web ignored a stored one (OV-793).
    request.layout = {};
  }
  liftComparisons(request);
  if (request.schema === 1 && Array.isArray(request.records)) {
    request.records = request.records.map((/** @type {unknown} */ record, /** @type {number} */ index) => (
      isObject(record) && !Object.hasOwn(record, 'recordKey') ? { ...record, recordKey: `record-${index + 1}` } : record
    ));
  }
  request.grouping = inferredGrouping(request);
  request.schema = LIFTED_SCHEMA;
  return { value: request, repairs };
};
