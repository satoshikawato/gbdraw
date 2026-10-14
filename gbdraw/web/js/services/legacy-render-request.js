// @ts-check
// A stored render request of schema 1 or 2 (Sessions 31-33) read as the
// schema-5 request that means the same (PR-D). The reader of stored requests
// then reads schema 5 and later only. Python's decoder owns the tables
// (`describe_legacy_request_tables` in gbdraw/session_request_codec.py,
// generated into LEGACY_REQUEST_TABLES); where the Web read a form more
// leniently than Python (OV-793, OV-795), the lift keeps the Web's reading.
import { LEGACY_REQUEST_TABLES } from '../contracts/schema.generated.js';
import { loadRepair } from '../utils/load-repairs.js';
import { classifyOptionalPositiveNumber } from '../utils/optional-positive-number.js';
import { migrateLegacyCircularTrackSlot, migrateLegacyCircularTrackSlotSpec } from './circular-track-slot-model.js';
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
  for (const [path, value] of Object.entries(LEGACY_REQUEST_TABLES.configOverrides[mode])) {
    if (!Object.hasOwn(overrides, path)) overrides[path] = value;
  }
  options.configOverrides = overrides;
};

/**
 * The retired names and values of configuration overrides (A7) and an
 * invalid Linear comparison height (B17: drawn before b1b48f6e0, now unset).
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
  const config = options.config;
  if (isObject(config)) {
    const placement = LEGACY_REQUEST_TABLES.linearConfigOverrideValues.label_placement;
    const trackLayout = LEGACY_REQUEST_TABLES.linearConfigOverrideValues.linear_track_layout;
    const labels = isObject(config.labels) && isObject(config.labels.linear) ? config.labels.linear : null;
    if (labels && Object.hasOwn(labels, 'placement')) labels.placement = renamedToken(labels.placement, placement);
    const canvas = isObject(config.canvas) && isObject(config.canvas.linear) ? config.canvas.linear : null;
    if (canvas && Object.hasOwn(canvas, 'track_layout')) canvas.track_layout = renamedToken(canvas.track_layout, trackLayout);
  }
  for (const key of ['comparison_height', ...LEGACY_REQUEST_TABLES.configOverridePaths.comparison_height]) {
    if (!Object.hasOwn(overrides, key) || classifyOptionalPositiveNumber(overrides[key]).status !== 'invalid') continue;
    repairs.push(loadRepair({
      code: 'setting-replaced', path: `${OPTIONS_PATH}.configOverrides.${key}`, effect: 'dropped',
      stored: String(overrides[key])
    }));
    delete overrides[key];
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

/**
 * A Circular object slot without the legacy spacing (A8): the spacing fills
 * the gaps the slot leaves unset. A spacing that cannot be read is left for
 * the reader to refuse.
 * @param {JsonObject} slot
 * @returns {JsonObject}
 */
const liftCircularObjectSlot = (slot) => {
  const migrated = /** @type {JsonObject} */ (migrateLegacyCircularTrackSlot(slot));
  /** @type {JsonObject} */
  const lifted = {};
  for (const [key, value] of Object.entries(migrated)) {
    if (key !== 'inner_gap_px' && key !== 'outer_gap_px') lifted[key] = value;
  }
  for (const [legacy, current] of /** @type {const} */ ([['inner_gap_px', 'innerGapPx'], ['outer_gap_px', 'outerGapPx']])) {
    if (!Object.hasOwn(migrated, legacy)) continue;
    const pixel = legacyPixel(migrated[legacy]);
    if (pixel === undefined) return slot;
    lifted[current] = pixel;
  }
  return lifted;
};

/**
 * The Circular slots under current tokens (A8, A9).
 * @param {JsonObject} options
 */
const liftCircularSlots = (options) => {
  const tracks = options.tracks;
  if (!isObject(tracks) || !Array.isArray(tracks.circularTrackSlots)) return;
  tracks.circularTrackSlots = tracks.circularTrackSlots.map((/** @type {unknown} */ slot) => (
    isObject(slot) ? liftCircularObjectSlot(slot) : migrateLegacyCircularTrackSlotSpec(slot)
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
 * schema-5 request cannot carry (B17); a lossless lift has none. A value the
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
    liftCircularSlots(options);
    if (isObject(request.layout) && Object.hasOwn(request.layout, 'multiRecordSizeMode')
      && foldedText(request.layout.multiRecordSizeMode) === 'sqrt') {
      request.layout.multiRecordSizeMode = 'auto';
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
