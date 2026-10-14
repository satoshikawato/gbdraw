// @ts-check
// A stored render request of schema 1 or 2 (Sessions 31-33) read as the
// schema-5 request that means the same (PR-D). The reader of stored requests
// then reads schema 5 and later only. Python's decoder owns the tables
// (`describe_legacy_request_tables` in gbdraw/session_request_codec.py,
// generated into LEGACY_REQUEST_TABLES); where the Web read a form more
// leniently than Python (OV-793, OV-795, OV-810) or with another meaning
// (OV-798), the lift keeps the Web's reading.
import { LEGACY_REQUEST_TABLES } from '../contracts/schema.generated.js';
import { loadRepair } from '../utils/load-repairs.js';
import { classifyOptionalPositiveNumber, DECIMAL_NUMBER_PATTERN } from '../utils/optional-positive-number.js';
import { migrateLegacyCircularTrackSlot, readLegacyCircularTrackSlotSpec } from './circular-track-slot-model.js';
import { parseOptionalPixel } from './track-slot-validation.js';

/** @import { LoadRepair, Repaired } from '../utils/load-repairs.js' */
/** @typedef {Record<string, unknown>} JsonObject */

const LIFTED_SCHEMA = 5;
const LEGACY_SCHEMAS = new Set([1, 2]);
const OPTIONS_PATH = 'renderRequest.diagramOptions';

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
  if (isUnset(shapes.repeat_region)) shapes.repeat_region = 'rectangle';
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
 * invalid Linear comparison height (B17: drawn before b1b48f6e0, now unset
 * with a row).
 * @param {JsonObject} options
 * @param {'circular' | 'linear'} mode
 * @param {LoadRepair[]} repairs
 */
const liftConfigOverrides = (options, mode, repairs) => {
  const overrides = options.configOverrides;
  if (!isObject(overrides)) return;
  for (const [old, current] of Object.entries(LEGACY_REQUEST_TABLES.configOverrideRenames)) {
    if (!Object.hasOwn(overrides, old)) continue;
    const stored = overrides[old];
    delete overrides[old];
    if (isUnset(overrides[current])) overrides[current] = stored;
  }
  if (mode !== 'linear') return;
  /** @type {Readonly<Record<string, readonly string[]>>} */
  const canonicalPaths = LEGACY_REQUEST_TABLES.configOverridePaths;
  for (const [name, table] of Object.entries(LEGACY_REQUEST_TABLES.linearConfigOverrideValues)) {
    for (const key of [name, ...(canonicalPaths[name] ?? [])]) {
      if (Object.hasOwn(overrides, key)) overrides[key] = renamedToken(overrides[key], table);
    }
  }
  for (const key of ['comparison_height', ...LEGACY_REQUEST_TABLES.configOverridePaths.comparison_height]) {
    if (!Object.hasOwn(overrides, key) || classifyOptionalPositiveNumber(overrides[key]).status !== 'invalid') continue;
    repairs.push(loadRepair({
      code: 'setting-replaced', path: `${OPTIONS_PATH}.configOverrides.${key}`, effect: 'dropped', stored: excerpt(overrides[key])
    }));
    delete overrides[key];
  }
};

/**
 * The Linear leaves of a stored full config the schema-5 reader requires in
 * their current form (A7: the two rewrites of Python's full-config
 * migration), and an invalid comparison height, which the config cannot
 * leave unset: its default is used (B17, F1), with a row unless a valid
 * override of the leaf was drawn instead (N3).
 * @param {JsonObject} options
 * @param {LoadRepair[]} repairs
 */
const liftLinearFullConfig = (options, repairs) => {
  const config = options.config;
  if (!isObject(config)) return;
  const values = LEGACY_REQUEST_TABLES.linearConfigOverrideValues;
  const labels = isObject(config.labels) && isObject(config.labels.linear) ? config.labels.linear : null;
  if (labels && Object.hasOwn(labels, 'placement')) labels.placement = renamedToken(labels.placement, values.label_placement);
  const canvas = isObject(config.canvas) && isObject(config.canvas.linear) ? config.canvas.linear : null;
  if (!canvas) return;
  if (Object.hasOwn(canvas, 'track_layout')) canvas.track_layout = renamedToken(canvas.track_layout, values.linear_track_layout);
  if (Object.hasOwn(canvas, 'comparison_height') && classifyOptionalPositiveNumber(canvas.comparison_height).status === 'invalid') {
    const overrides = isObject(options.configOverrides) ? options.configOverrides : {};
    const overridden = ['comparison_height', ...LEGACY_REQUEST_TABLES.configOverridePaths.comparison_height]
      .some((key) => Object.hasOwn(overrides, key) && classifyOptionalPositiveNumber(overrides[key]).status === 'valid');
    if (!overridden) {
      repairs.push(loadRepair({
        code: 'setting-replaced', path: `${OPTIONS_PATH}.config.canvas.linear.comparison_height`, effect: 'defaulted',
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

const PRIVATE_SPACING = '__gbdraw_legacy_spacing';
// A slot's legacy spacing, in both forms: its `spacing`, else the private
// spacing parameter Python's decoder read (N2, N3).
const LEGACY_SPACING_KEYS = ['spacing', PRIVATE_SPACING];
const GAP_NAMES = [['inner_gap_px', 'innerGapPx'], ['outer_gap_px', 'outerGapPx']];
const GAP_KEYS = new Set(GAP_NAMES.flat());

/**
 * How Python read a legacy slot spacing: in pixels (`4px`, `{ unit: 'px' }`),
 * as a multiple of the axis radius when it has no unit (`4`) or a factor unit,
 * or not at all. The Web read a spacing without a unit as pixels (OV-798).
 * @param {unknown} value
 * @returns {'unitless' | 'factor' | 'other'}
 */
const spacingUnit = (value) => {
  if (isObject(value)) return foldedText(value.unit) === 'factor' ? 'factor' : 'other';
  if (typeof value === 'number') return 'unitless';
  return typeof value === 'string' && DECIMAL_NUMBER_PATTERN.test(value.trim()) ? 'unitless' : 'other';
};

/**
 * A factor spacing of 0: no gap, read as 0 px with no row (N4).
 * @param {unknown} value
 */
const isZeroFactor = (value) => isObject(value) && spacingUnit(value) === 'factor' && legacyPixel({ ...value, unit: 'px' }) === 0;

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
 * A Circular string slot under current tokens (A9), its legacy spacing read
 * as the object form reads it (N2), with the row of the spacing token the
 * parser consumed when it has no unit (OV-798).
 * @param {unknown} spec
 * @param {string} at its path in the file
 * @param {LoadRepair[]} repairs
 */
const liftCircularStringSlot = (spec, at, repairs) => {
  const { spec: lifted, spacing } = readLegacyCircularTrackSlotSpec(spec, LEGACY_SPACING_KEYS);
  const repair = spacing?.filled ? spacingRepair(spacing.value, at, spacing.stored) : null;
  if (repair) repairs.push(repair);
  return lifted;
};

/**
 * The Circular slots under current tokens (A8, A9).
 * @param {JsonObject} options
 * @param {LoadRepair[]} repairs
 */
const liftCircularSlots = (options, repairs) => {
  const tracks = options.tracks;
  if (!isObject(tracks) || !Array.isArray(tracks.circularTrackSlots)) return;
  const at = `${OPTIONS_PATH}.tracks.circularTrackSlots`;
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
 * A stored render request of schema 1 or 2 as the schema-5 request that
 * means the same; any other value unchanged. The repairs list what the
 * schema-5 request cannot carry or reads with another meaning (B17,
 * OV-798); a lossless lift has none. A value the
 * reader refuses (not an object, an unknown mode) is left to the reader.
 * @param {unknown} value a stored `renderRequest` (parsed JSON)
 * @returns {Repaired<unknown>}
 */
export const liftLegacyRenderRequest = (value) => {
  if (!isObject(value) || !LEGACY_SCHEMAS.has(/** @type {number} */ (value.schema))
    || (value.mode !== 'circular' && value.mode !== 'linear') || !isObject(value.diagramOptions)) {
    return { value, repairs: [] };
  }
  /** @type {LoadRepair[]} */
  const repairs = [];
  const mode = value.mode;
  const request = /** @type {JsonObject} */ (JSON.parse(JSON.stringify(value)));
  const options = /** @type {JsonObject} */ (request.diagramOptions);
  if (isObject(options.output)) delete options.output.outputPrefix;
  liftFeatureVisibilityNames(options);
  liftSparseDefaults(options, mode);
  liftConfigOverrides(options, mode, repairs);
  if (mode === 'circular') {
    liftCircularSlots(options, repairs);
    if (isObject(request.layout) && Object.hasOwn(request.layout, 'multiRecordSizeMode')
      && foldedText(request.layout.multiRecordSizeMode) === 'sqrt') {
      request.layout.multiRecordSizeMode = 'auto';
    }
  } else {
    liftLinearFullConfig(options, repairs);
    // Schema 1 had no Linear layout; the Web ignored a stored one (OV-793).
    if (request.schema === 1) request.layout = {};
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
