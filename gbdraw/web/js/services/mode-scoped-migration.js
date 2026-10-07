// @ts-check
// Session 46 keeps each diagram mode's settings and edits in its own slice
// (`modes.circular`, `modes.linear`, PD-OI-086). This state-free module owns
// the registry's slice shape (mode-scoped-settings.generated.js, written from
// gbdraw/web_support/mode_scoped_settings.py) and the split of a Session 27-44
// draft into the two slices. Its Python twin is `split_draft_into_modes` in
// gbdraw/session_io.py; tests/fixtures/sessions/mode-split-vectors.json pins
// both.
import { MODE_SCOPED_SETTINGS } from '../mode-scoped-settings.generated.js';
import { normalizePaletteColors } from '../utils/color-utils.js';
import { cloneJsonData } from './json-clone.js';
import { diagnosticError } from '../utils/error-normalization.js';
import { featureIdentityKey } from './feature-placement.js';

/** @typedef {'circular' | 'linear'} SliceMode */
/**
 * One registry row: a value that each mode's slice holds on its own.
 * @typedef {{ domain: string, path: string, modes: string, key: string | null, migrate: string }} ModeScopedRow
 */

/** @type {readonly SliceMode[]} */
export const SLICE_MODES = Object.freeze(['circular', 'linear']);
/** @type {readonly ModeScopedRow[]} */
export const MODE_SCOPED_ROWS = MODE_SCOPED_SETTINGS.rows;
/** @type {readonly string[]} */
export const LOSAT_EXECUTION_FIELDS = MODE_SCOPED_SETTINGS.losatExecutionFields;

const isObject = (value) => value !== null && typeof value === 'object' && !Array.isArray(value);
const fieldsError = () => diagnosticError('INPUT_INVALID', { field: 'schema', reason: 'FIELDS' });
/** @param {unknown} value */
const sliceMode = (value) => (value === 'circular' || value === 'linear' ? value : null);
/** @param {unknown} value */
const jsTruthy = (value) => Boolean(value);

// Each slice container ('' is the slice itself) and the keys it may hold.
/** @type {Readonly<Record<string, ReadonlySet<string>>>} */
export const MODE_SLICE_CONTAINERS = (() => {
  /** @type {Record<string, Set<string>>} */
  const containers = { '': new Set() };
  MODE_SCOPED_ROWS.forEach((row) => {
    const parts = row.domain.split('.');
    parts.forEach((part, depth) => {
      const parent = parts.slice(0, depth).join('.');
      (containers[parent] ||= new Set()).add(part);
      containers[parts.slice(0, depth + 1).join('.')] ||= new Set();
    });
    containers[row.domain].add(row.path);
  });
  return Object.freeze(containers);
})();

// A slice holds only registry fields: the writer and both validators reject
// any other key (plan 4.2).
/**
 * @param {unknown} value
 * @param {string} [container]
 */
export const validateModeSliceFields = (value, container = '') => {
  if (!isObject(value)) throw fieldsError();
  const allowed = MODE_SLICE_CONTAINERS[container];
  for (const [field, child] of Object.entries(/** @type {Record<string, unknown>} */ (value))) {
    if (!allowed.has(field)) throw fieldsError();
    const childContainer = container ? `${container}.${field}` : field;
    if (Object.hasOwn(MODE_SLICE_CONTAINERS, childContainer)) validateModeSliceFields(child, childContainer);
  }
};

// The modes a GUI-unmanaged config override leaf belongs to (`by-leaf`): the
// leaf of one mode's renderer prefixes belongs to that mode, any other leaf to
// both (gbdraw.api.options._MODE_CONFIG_OVERRIDE_PREFIXES).
/**
 * @param {string} path
 * @returns {readonly SliceMode[]}
 */
export const unmanagedConfigOverrideModes = (path) => {
  for (const mode of SLICE_MODES) {
    if (MODE_SCOPED_SETTINGS.unmanagedConfigOverrideModePrefixes[mode]
      .some((/** @type {string} */ prefix) => path === prefix || path.startsWith(`${prefix}.`))) {
      return [mode];
    }
  }
  return SLICE_MODES;
};

// depthFileSlotsFromValue: a list is the series slots; a value is one slot.
/** @param {unknown} value */
const depthSlots = (value) => (Array.isArray(value) ? [...value] : jsTruthy(value) ? [value] : []);

// Each mode's Depth series count in the Web file bindings, 0 when none is
// bound: the widest row when any slot holds a file.
/**
 * @param {unknown} bindings
 * @returns {Record<SliceMode, number>}
 */
export const sessionDepthSourceWidths = (bindings) => {
  const source = isObject(bindings) ? /** @type {Record<string, any>} */ (bindings) : {};
  const circular = source.c_depth;
  const circularRows = Array.isArray(circular) && circular.every(Array.isArray)
    ? circular.map(depthSlots)
    : [depthSlots(circular)];
  const linearRows = (Array.isArray(source.linearSeqs) ? source.linearSeqs : [])
    .filter(isObject).map((sequence) => depthSlots(sequence.depth));
  /** @param {unknown[][]} rows */
  const width = (rows) => (rows.some((row) => row.some(jsTruthy)) ? Math.max(...rows.map((row) => row.length)) : 0);
  return { circular: width(circularRows), linear: width(linearRows) };
};

// The mode whose values a flat draft holds for the mode-profile fields:
// `ui.mode`, else `renderRequest.mode`, else `modeProfiles.activeMode`.
/**
 * @param {Record<string, any>} draft
 * @param {unknown} modeProfiles
 * @returns {SliceMode | null}
 */
export const profileActiveMode = (draft, modeProfiles) => (
  sliceMode(draft?.ui?.mode)
  || sliceMode(draft?.renderRequest?.mode)
  || sliceMode(isObject(modeProfiles) ? /** @type {Record<string, any>} */ (modeProfiles).activeMode : null)
);

/**
 * @param {unknown} modeProfiles
 * @param {SliceMode} mode
 * @returns {Record<string, any>}
 */
const modeProfileValues = (modeProfiles, mode) => {
  const values = /** @type {Record<string, any>} */ (modeProfiles)?.profiles?.[mode]?.values;
  return isObject(values) ? values : {};
};

// The scoped draft rows of `mode` (a Session 41-44 migration's rows), keyed by
// identity pair and without `scope`; a row that names no mode is kept as it is.
/**
 * @param {unknown} rows
 * @param {SliceMode} mode
 */
export const unscopedDraftRows = (rows, mode) => {
  if (!isObject(rows)) return rows;
  /** @type {Array<[string, any]>} */
  const kept = [];
  Object.entries(/** @type {Record<string, any>} */ (rows)).forEach(([key, row]) => {
    if (isObject(row) && Object.hasOwn(row, 'scope')) {
      if (row.scope !== mode) return;
      const { scope: _scope, ...fields } = row;
      kept.push([featureIdentityKey(fields.recordKey, fields.biologicalFeatureId) || key, fields]);
    } else {
      kept.push([key, row]);
    }
  });
  return Object.fromEntries(kept);
};

const ANNOTATION_RECORD_BINDING_KEY = '_gbdraw_web_target_record_key';

// The mode whose record an annotation's target names, if it names one: a
// featureIdentity target's `scope`, or the mode its record binding starts with.
/** @param {unknown} annotation */
const annotationBindingMode = (annotation) => {
  if (!isObject(annotation)) return null;
  const { target, metadata } = /** @type {Record<string, any>} */ (annotation);
  if (isObject(target) && target.kind === 'featureIdentity') return sliceMode(target.scope);
  const binding = isObject(metadata) ? metadata[ANNOTATION_RECORD_BINDING_KEY] : null;
  if (typeof binding !== 'string') return null;
  const text = binding.trim();
  // A record key starts with its source key, a JSON array whose first item is the mode.
  for (let end = text.indexOf(']'); end >= 0; end = text.indexOf(']', end + 1)) {
    try {
      const sourceKey = JSON.parse(text.slice(0, end + 1));
      return Array.isArray(sourceKey) && sourceKey.length ? sliceMode(sourceKey[0]) : null;
    } catch {
      // Not yet the end of the first JSON value.
    }
  }
  return null;
};

/**
 * @param {unknown[]} sets
 * @param {SliceMode} mode
 */
const annotationSetsOfMode = (sets, mode) => sets.map((annotationSet) => {
  const annotations = isObject(annotationSet) ? /** @type {Record<string, any>} */ (annotationSet).annotations : null;
  if (!Array.isArray(annotations)) return cloneJsonData(annotationSet);
  const kept = annotations.flatMap((annotation) => {
    const bound = annotationBindingMode(annotation);
    if (bound !== null && bound !== mode) return [];
    const copy = cloneJsonData(annotation);
    if (isObject(copy?.target) && copy.target.kind === 'featureIdentity') delete copy.target.scope;
    return [copy];
  });
  return { ...cloneJsonData(annotationSet), annotations: kept };
});

/**
 * The value each slice takes for one saved registry value.
 * @param {ModeScopedRow} row
 * @param {any} value
 * @param {{ committed: SliceMode, widths: Record<SliceMode, number> }} context
 * @returns {Partial<Record<SliceMode, any>>}
 */
const splitSetting = (row, value, { committed, widths }) => {
  if (row.migrate === 'own') return { [/** @type {SliceMode} */ (row.modes)]: cloneJsonData(value) };
  if (row.migrate === 'result-mode') return { [committed]: unscopedDraftRows(cloneJsonData(value), committed) };
  if (row.migrate === 'layout') {
    const slots = isObject(value) ? value : {};
    return Object.fromEntries(SLICE_MODES.filter((mode) => Object.hasOwn(slots, mode))
      .map((mode) => [mode, cloneJsonData(slots[mode])]));
  }
  /** @type {Partial<Record<SliceMode, any>>} */
  const split = {};
  SLICE_MODES.forEach((mode) => {
    if (row.migrate === 'show-if-source') {
      split[mode] = jsTruthy(value) ? widths[mode] > 0 : value;
    } else if (row.migrate === 'depth' && Array.isArray(value)) {
      split[mode] = cloneJsonData(value.slice(0, Math.max(1, widths[mode])));
    } else if (row.migrate === 'by-scope' && Array.isArray(value)) {
      split[mode] = value.filter((draft) => isObject(draft) && draft.scope === mode).map((draft) => {
        const { scope: _scope, ...fields } = cloneJsonData(draft);
        return fields;
      });
    } else if (row.migrate === 'by-side') {
      split[mode] = unscopedDraftRows(cloneJsonData(value), mode);
    } else if (row.migrate === 'by-leaf' && isObject(value)) {
      split[mode] = Object.fromEntries(Object.entries(value)
        .filter(([path]) => unmanagedConfigOverrideModes(String(path)).includes(mode))
        .map(([path, leaf]) => [path, cloneJsonData(leaf)]));
    } else if (row.migrate === 'by-binding' && Array.isArray(value)) {
      split[mode] = annotationSetsOfMode(value, mode);
    } else {
      split[mode] = cloneJsonData(value);
    }
  });
  return split;
};

// A Session before the pairwise match style existed drew ribbons; the flat
// draft of such a Session holds that value for its active mode.
/** @type {Record<string, unknown>} */
const HISTORICAL_FLAT_PROFILE_VALUES = Object.freeze({ pairwise_match_style: 'ribbon' });
// A registry value that a Session 44 or older saved elsewhere.
/** @type {Record<string, [string, string]>} */
const FLAT_DRAFT_SOURCES = Object.freeze({ 'ui.selectedFeatureRecordIdx': ['features', 'selectedFeatureRecordIdx'] });
// The top-level homes that Session 46 moved into `modes`.
const FLAT_DRAFT_TOP_LEVEL_FIELDS = Object.freeze(['config', 'features']);
/** @type {ReadonlyArray<[string, ReadonlySet<string>]>} */
const RETIRED_NESTED_FIELDS = Object.freeze([
  ['editorState', new Set(['featureStrokes'])],
  ['editorState.legend', new Set(MODE_SCOPED_ROWS.filter((row) => row.domain === 'editorState.legend').map((row) => row.path))],
  ['ui', new Set(MODE_SCOPED_ROWS.filter((row) => row.domain === 'ui').map((row) => row.path))]
]);

/**
 * @param {unknown} source
 * @param {string} domain
 */
const containerAt = (source, domain) => domain.split('.').reduce(
  (/** @type {any} */ current, part) => (isObject(current) ? current[part] : undefined),
  source
);

// The draft `config` with `colorsAreOverrides` resolved and dropped: with the
// flag and colors, the stored colors override the palette's colors, as Load
// reads them (`applyConfigData`); otherwise the stored colors are complete.
/**
 * @param {Record<string, any>} config
 * @param {Record<string, any> | null | undefined} paletteColors
 */
const withResolvedOverrideColors = (config, paletteColors) => {
  const { colorsAreOverrides: _flag, ...resolved } = config;
  const colors = config.colors;
  if (jsTruthy(config.colorsAreOverrides) && isObject(colors) && Object.keys(colors).length > 0) {
    const overrides = Object.fromEntries(Object.entries(colors)
      .map(([key, value]) => [key, jsTruthy(value) ? String(value).trim() : '']));
    resolved.colors = normalizePaletteColors({ ...(paletteColors || {}), ...overrides });
  }
  return resolved;
};

/**
 * Splits the flat Web draft of a Session 27-44 into Session 46 mode slices.
 *
 * `draft` holds the Session's `config`, `features`, `editorState`, and `ui`
 * after the older normalizers (field names and placement rows, per-feature
 * edits, annotation targets); its other fields are kept. The result has
 * `modes` instead of the flat draft when the draft has a `config`: each
 * registry row fills the slices by its `migrate` token, a value a slice does
 * not take is absent there (that mode's default), and a field no row names is
 * dropped. App-level settings move to the top level: `config.losat`'s
 * execution settings to `ui.losatExecution`, `config.adv.rich_feature_popup`
 * to `ui.richFeaturePopup`, a boolean `config.paletteInstantPreviewEnabled` to
 * `ui` (over a saved `ui` value, as Load applies it last), and
 * `config.cliOptions` to `cliOptions`.
 *
 * @param {Record<string, any>} draft
 * @param {{
 *   committedMode: unknown,
 *   modeProfiles?: unknown,
 *   depthSources?: Partial<Record<SliceMode, number>> | null,
 *   paletteColors?: Record<string, any> | null
 * }} context `committedMode` is the saved Result's mode (`renderRequest.mode`,
 *   else `ui.mode`); `modeProfiles` defaults to `config.modeProfiles`;
 *   `depthSources` (each mode's Depth series count) to the counts of
 *   `webFiles.bindings`; `paletteColors` are the colors of the draft's
 *   palette, which override colors merge into.
 * @returns {Record<string, any>}
 */
export const splitDraftIntoModes = (draft, { committedMode, modeProfiles = undefined, depthSources = null, paletteColors = null }) => {
  const committed = sliceMode(committedMode);
  if (!committed) throw fieldsError();
  let configValue = draft.config;
  if (isObject(configValue) && Object.hasOwn(configValue, 'colorsAreOverrides')) {
    configValue = withResolvedOverrideColors(configValue, paletteColors);
    draft = { ...draft, config: configValue };
  }
  /** @type {Record<string, any>} */
  const config = isObject(configValue) ? configValue : {};
  const profiles = modeProfiles === undefined ? config.modeProfiles : modeProfiles;
  const active = profileActiveMode(draft, profiles) || committed;
  const widths = depthSources
    ? { circular: Number(depthSources.circular) || 0, linear: Number(depthSources.linear) || 0 }
    : sessionDepthSourceWidths(isObject(draft.webFiles) ? draft.webFiles.bindings : draft.files);
  const adv = isObject(config.adv) ? config.adv : {};

  /** @type {Record<SliceMode, Record<string, any>>} */
  const slices = { circular: {}, linear: {} };
  MODE_SCOPED_ROWS.forEach((row) => {
    const [sourceDomain, sourcePath] = FLAT_DRAFT_SOURCES[`${row.domain}.${row.path}`] || [row.domain, row.path];
    const container = containerAt(draft, sourceDomain);
    if (!isObject(container)) return;
    /** @type {Partial<Record<SliceMode, any>>} */
    let split;
    if (row.migrate === 'profile') {
      split = {};
      const flat = Object.hasOwn(container, sourcePath)
        ? container[sourcePath]
        : HISTORICAL_FLAT_PROFILE_VALUES[sourcePath];
      if (Object.hasOwn(container, sourcePath) || Object.hasOwn(HISTORICAL_FLAT_PROFILE_VALUES, sourcePath)) {
        split[active] = cloneJsonData(flat);
      }
      const other = active === 'circular' ? 'linear' : 'circular';
      const saved = modeProfileValues(profiles, other);
      if (Object.hasOwn(saved, sourcePath)) split[other] = cloneJsonData(saved[sourcePath]);
    } else if (Object.hasOwn(container, sourcePath)) {
      split = splitSetting(row, container[sourcePath], { committed, widths });
    } else {
      return;
    }
    /** @type {Array<[SliceMode, any]>} */ (Object.entries(split)).forEach(([mode, value]) => {
      let target = slices[mode];
      row.domain.split('.').forEach((part) => { target = (target[part] ||= {}); });
      target[row.path] = value;
    });
  });

  /** @type {Record<string, any>} */
  const result = Object.fromEntries(Object.entries(draft)
    .filter(([key]) => !FLAT_DRAFT_TOP_LEVEL_FIELDS.includes(key)));
  RETIRED_NESTED_FIELDS.forEach(([domain, fields]) => {
    const [head, rest] = domain.split(/\.(.*)/s);
    const container = containerAt(result, domain);
    if (!isObject(container) || !Object.keys(container).some((key) => fields.has(key))) return;
    const kept = Object.fromEntries(Object.entries(container).filter(([key]) => !fields.has(key)));
    if (rest) {
      const parent = { ...result[head] };
      if (Object.keys(kept).length) parent[rest] = kept;
      else delete parent[rest];
      result[head] = parent;
    } else {
      result[head] = kept;
    }
  });
  const losat = config.losat;
  const execution = Object.fromEntries(LOSAT_EXECUTION_FIELDS
    .filter((field) => isObject(losat) && Object.hasOwn(losat, field))
    .map((field) => [field, cloneJsonData(losat[field])]));
  // App-level settings leave the draft: LOSAT execution, the rich feature
  // popup (a missing value reads true), and Instant Preview.
  /** @type {Record<string, any>} */
  const appUi = Object.keys(execution).length ? { losatExecution: execution } : {};
  if (Object.hasOwn(adv, 'rich_feature_popup')) appUi.richFeaturePopup = cloneJsonData(adv.rich_feature_popup);
  // Load applies `ui`, then a boolean draft value, which wins.
  if (typeof config.paletteInstantPreviewEnabled === 'boolean') {
    appUi.paletteInstantPreviewEnabled = config.paletteInstantPreviewEnabled;
  }
  if (Object.keys(appUi).length) result.ui = { ...(result.ui || {}), ...appUi };
  // Session provenance, written only when the draft has it.
  if (Object.hasOwn(config, 'cliOptions')) result.cliOptions = cloneJsonData(config.cliOptions);
  if (isObject(configValue)) result.modes = { circular: slices.circular, linear: slices.linear };
  return result;
};
