// @ts-check
// Session 46 keeps each diagram mode's settings and edits in its own slice
// (`modes.circular`, `modes.linear`, PD-OI-086). This state-free module owns
// the registry's slice shape (mode-scoped-settings.generated.js, written from
// gbdraw/web_support/mode_scoped_settings.py) and the split of a Session 27-44
// draft into the two slices. Its Python twin is `split_draft_into_modes` in
// gbdraw/session_io.py; tests/fixtures/sessions/mode-split-vectors.json pins
// both. It also reads a Session 46 Legend's edits by row key (Session 47,
// `legendEditsByRowKey`; twin `legend_edits_by_row_key`, pinned by
// tests/fixtures/sessions/legend-row-key-vectors.json).
import { MODE_SCOPED_SETTINGS } from '../mode-scoped-settings.generated.js';
// `FLAT_DRAFT_TOP_LEVEL_FIELDS`: the top-level homes that Session 46 moved into `modes`.
import { FLAT_DRAFT_TOP_LEVEL_FIELDS, ROW_KEYED_LEGEND_SESSION_VERSION } from '../contracts/schema.generated.js';
import { legendEntryKey } from './legend-svg.js';
import { normalizeDefaultColor, normalizePaletteColors } from '../utils/color-utils.js';
import { normalizeLegendColorOverrides, normalizeStrokeOverrideMap } from './style-overrides.js';
/** @import { LegendEditsUnread } from './style-overrides.js' */
import { cloneJsonData } from './json-clone.js';
import { diagnosticError } from '../utils/error-normalization.js';
import { featureIdentityKey } from './feature-placement.js';
import { ANNOTATION_RECORD_BINDING_KEY } from './annotation-state.js';

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

/** @param {unknown} value */
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
 * The places of a slice that hold no registry field: a container that is no
 * object (`field` null), or a key it may not hold.
 * @param {unknown} value
 * @param {string} [container]
 * @returns {{ container: string, field: string | null }[]}
 */
export const modeSliceFieldViolations = (value, container = '') => {
  if (!isObject(value)) return [{ container, field: null }];
  const allowed = MODE_SLICE_CONTAINERS[container];
  return Object.entries(/** @type {Record<string, unknown>} */ (value)).flatMap(([field, child]) => {
    if (!allowed.has(field)) return [{ container, field }];
    const childContainer = container ? `${container}.${field}` : field;
    return Object.hasOwn(MODE_SLICE_CONTAINERS, childContainer) ? modeSliceFieldViolations(child, childContainer) : [];
  });
};
/** Save: a slice holds only registry fields. @param {unknown} value */
export const validateModeSliceFields = (value) => {
  if (modeSliceFieldViolations(value).length) throw fieldsError();
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
  // `annotations` is an array, so `annotationSet` is an object.
  return { ...cloneJsonData(/** @type {Record<string, unknown>} */ (annotationSet)), annotations: kept };
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

// The Stroke options disclosure (`showStroke`) that earlier Sessions saved on
// a Legend row is view state (OV-157): the split drops it from the Legend
// rows, as Load does (`normalizeSessionLegendEntries`).
const LEGEND_ROW_LISTS = ['entries', 'deletedEntries', 'dormantEntries'];
/** @param {Record<string, any>} draft */
const withoutStrokeDisclosure = (draft) => {
  const legend = isObject(draft.editorState) ? draft.editorState.legend : null;
  if (!isObject(legend)) return draft;
  const changed = { ...legend };
  for (const rows of LEGEND_ROW_LISTS) {
    if (!Array.isArray(legend[rows])) continue;
    changed[rows] = legend[rows].map((/** @type {unknown} */ entry) => {
      if (!isObject(entry)) return entry;
      const { showStroke: _showStroke, ...row } = /** @type {Record<string, any>} */ (entry);
      return row;
    });
  }
  return { ...draft, editorState: { ...draft.editorState, legend: changed } };
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
 * `config.cliOptions` to `cliOptions`. The Legend rows lose `showStroke`.
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
  draft = withoutStrokeDisclosure(draft);
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

/** @typedef {'colorOverrides' | 'strokeOverrides'} LegendEditStore */
/**
 * One Session 46 Legend edit that did not stay as it was (`legendEditsByRowKey`):
 * `moved` to the key of the row whose caption it was stored under; kept under
 * its key as `unmatched` (no row shows or holds the key) or as `unshown` (the
 * key of a renamed row that no row shows: Session 46 applied it to that row
 * through Python's key); or a `collision`: a later edit of the store reached
 * the same row (`rowKey`) and Session 46 applied it over this one, dropped
 * (`toKey: null`) unless a stroke field of it is kept (`toKey`: `rowKey`). `slice`: the mode slice, or `draft` (`sessionLegendEditsByRowKey`).
 * @typedef {{ store: LegendEditStore, fromKey: string, toKey: string, outcome: 'moved' | 'unmatched' | 'unshown', slice?: string }
 *   | { store: LegendEditStore, fromKey: string, toKey: string | null, outcome: 'collision', rowKey: string, slice?: string }} LegendEditRekey
 */

// The rows a Session 46 Legend's edits reach (Session 46 Load,
// `normalizeSessionLegendEntries` and `normalizeEditorStateData` at
// 9839de8e: an entry with a caption and a usable color, the first of each
// caption in its list), by the caption each row shows (`draftLegendRows` at
// 9839de8e): a listed row before a dormant one (the later dormant entry of a
// caption wins), and a deleted row as its Restore shows it. A dormant entry
// that is not renamed, or whose key a listed or deleted row holds, is no row.
// `renamedKeys`: the keys of the renamed rows.
/** @param {Record<string, unknown>} legend */
const session46LegendRows = (legend) => {
  /** @param {string} list @returns {{ caption: string, key: string, dormant: boolean }[]} */
  const loaded = (list) => {
    const entries = legend[list];
    /** @type {Set<string>} */
    const captions = new Set();
    return (Array.isArray(entries) ? entries : []).flatMap((/** @type {unknown} */ entry) => {
      if (!isObject(entry)) return [];
      const row = /** @type {{ caption?: unknown, originalCaption?: unknown, color?: unknown, dormant?: unknown }} */ (entry);
      const caption = String(row.caption || '').trim();
      if (!caption || !normalizeDefaultColor(row.color) || captions.has(caption)) return [];
      captions.add(caption);
      return [{ caption, key: legendEntryKey(row), dormant: row.dormant === true }];
    });
  };
  const entries = loaded('entries');
  const listed = entries.filter(({ dormant }) => !dormant);
  const deleted = loaded('deletedEntries');
  const listedKeys = new Set(listed.map(({ key }) => key));
  const deletedKeys = new Set(deleted.map(({ key }) => key));
  const dormant = [...loaded('dormantEntries'), ...entries.filter((entry) => entry.dormant)]
    .filter(({ caption, key }) => caption !== key && !listedKeys.has(key) && !deletedKeys.has(key));
  /** @type {Map<string, string>} */
  const rowOf = new Map([...dormant, ...listed].map(({ caption, key }) => [caption, key]));
  deleted.forEach(({ caption, key }) => { if (!rowOf.has(caption)) rowOf.set(caption, key); });
  const renamedKeys = new Set([...listed, ...dormant, ...deleted].filter(({ caption, key }) => caption !== key).map(({ key }) => key));
  return { rowOf, renamedKeys };
};

// A Session 46 Legend (a slice's, or a flat draft's, `editorState.legend`)
// read as Session 47, and what the reading did to each edit that did not stay
// as it was. Session 46 stored the Legend color and stroke edits under the
// caption a row shows, Session 47 under the row key (`legendEntryKey`). Each
// row key gets the value Session 46 applied to that row (`compileDirectEditorMutationPlan`
// at 9839de8e, app/candidate-render.js): the stores as Session 46 Load read
// them (style-overrides.js), each edit reaching the row whose caption it is
// stored under, else the row of Python's key of that name; the color edits in
// their order, then the strokes in the order of the color edits' keys and then
// their own; a later edit painted over an earlier one, a stroke field by field
// (`strokeColor`, `strokeWidth`). Other fields are kept as they are. The
// outcomes are data for Load's repair report, computed before Load commits;
// the twin is `legend_edits_by_row_key` (gbdraw/session_io.py).
/**
 * @template {Record<string, unknown>} Legend
 * @param {Legend} legend
 * @param {LegendEditsUnread | null} [unread] Written: the edits the stores' reading drops (Load lists them).
 * @returns {{ legend: Legend, outcomes: LegendEditRekey[] }}
 */
export const legendEditsByRowKey = (legend, unread = null) => {
  const { rowOf, renamedKeys } = session46LegendRows(legend);
  /** @param {string} key */
  const rowKeyOf = (key) => rowOf.get(key) || key;
  const colors = normalizeLegendColorOverrides(legend.colorOverrides, unread?.colors);
  const strokes = normalizeStrokeOverrideMap(legend.strokeOverrides, { requireOverride: true, losses: unread?.strokes });
  /** @type {Record<string, unknown>} */
  const read = { ...legend };
  /** @type {LegendEditRekey[]} */
  const outcomes = [];
  /**
   * @template T
   * @param {LegendEditStore} store
   * @param {Record<string, T>} edits
   * @param {string[]} order
   * @param {(earlier: T | undefined, later: T) => T} paint
   * @param {(edit: T) => string[]} fieldsOf The fields an edit sets; a row keeps each from its last edit that sets it.
   */
  const keyByRow = (store, edits, order, paint, fieldsOf) => {
    if (!isObject(legend[store])) return;
    const lastOf = new Map(order.map((key) => [rowKeyOf(key), key]));
    /** @param {string} key @param {string} field */
    const fieldOf = (key, field) => JSON.stringify([rowKeyOf(key), field]);
    const fieldLastOf = new Map(order.flatMap((key) => fieldsOf(edits[key]).map((field) => [fieldOf(key, field), key])));
    /** @type {Record<string, T>} */
    const keyed = {};
    order.forEach((key) => {
      const rowKey = rowKeyOf(key);
      keyed[rowKey] = paint(keyed[rowKey], edits[key]);
      if (lastOf.get(rowKey) !== key) {
        const kept = fieldsOf(edits[key]).some((field) => fieldLastOf.get(fieldOf(key, field)) === key);
        outcomes.push({ store, fromKey: key, toKey: kept ? rowKey : null, outcome: 'collision', rowKey });
      } else if (rowKey !== key) outcomes.push({ store, fromKey: key, toKey: rowKey, outcome: 'moved' });
      else if (!rowOf.has(key)) outcomes.push({ store, fromKey: key, toKey: key, outcome: renamedKeys.has(key) ? 'unshown' : 'unmatched' });
    });
    read[store] = keyed;
  };
  const colorKeys = Object.keys(colors);
  keyByRow('colorOverrides', colors, colorKeys, (_earlier, later) => later, () => ['color']);
  keyByRow('strokeOverrides', strokes, [...new Set([...colorKeys, ...Object.keys(strokes)])].filter((key) => Object.hasOwn(strokes, key)),
    (earlier, later) => ({ ...earlier, ...later }), (edit) => Object.keys(edit));
  return { legend: /** @type {Legend} */ (read), outcomes };
};

// A Session document of `version` with each Legend's edits keyed by row key:
// every mode slice's (Session 46), and the flat draft's of a Session 27-44;
// with the outcomes of `legendEditsByRowKey`, each naming its slice. A
// Session 47 or later is returned as it is.
/**
 * @template {Record<string, unknown>} Session
 * @param {Session} session
 * @param {number} version
 * @param {{ slice: string, unread: LegendEditsUnread }[] | null} [unread] Written: what each slice's reading drops.
 * @returns {{ session: Session, outcomes: LegendEditRekey[] }}
 */
export const sessionLegendEditsByRowKey = (session, version, unread = null) => {
  /** @type {LegendEditRekey[]} */
  const outcomes = [];
  if (version >= ROW_KEYED_LEGEND_SESSION_VERSION) return { session, outcomes };
  /** @template {Record<string, unknown>} Holder @param {Holder} holder @param {string} slice @returns {Holder} */
  const withRowKeyedLegend = (holder, slice) => {
    const editor = holder.editorState;
    if (!isObject(editor) || !isObject(/** @type {Record<string, unknown>} */ (editor).legend)) return holder;
    const editorState = /** @type {Record<string, unknown>} */ (editor);
    /** @type {LegendEditsUnread} */
    const sliceUnread = { colors: [], strokes: { leaves: [], dropped: 0 } };
    const read = legendEditsByRowKey(/** @type {Record<string, unknown>} */ (editorState.legend), unread ? sliceUnread : null);
    read.outcomes.forEach((outcome) => outcomes.push({ ...outcome, slice }));
    unread?.push({ slice, unread: sliceUnread });
    return { ...holder, editorState: { ...editorState, legend: read.legend } };
  };
  const read = withRowKeyedLegend(session, 'draft');
  if (!isObject(read.modes)) return { session: read, outcomes };
  return {
    session: {
      ...read,
      modes: Object.fromEntries(Object.entries(/** @type {Record<string, unknown>} */ (read.modes)).map(([mode, slice]) => [
        mode, isObject(slice) ? withRowKeyedLegend(/** @type {Record<string, unknown>} */ (slice), mode) : slice
      ]))
    },
    outcomes
  };
};
