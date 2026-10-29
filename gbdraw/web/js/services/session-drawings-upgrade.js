// @ts-check
// The Web's move of a Session 46-49 to the drawings of a Session 50 (the
// 46-49 fast path, DESIGN §4.2; Python twin `drawings_document_from_mode_scoped`
// in gbdraw/session_drawings.py, pinned by
// tests/fixtures/sessions/drawings-vectors.json). A state-free module: every
// function reads only the document it is given and two generated modules,
// `SESSION_DRAWINGS_LAYOUT` (the Python layout's field sets,
// tools/generate_web_contracts.py) and `MODE_SLICE_DEFAULTS` (the Web's
// packaged per-mode defaults, tools/generate_mode_slice_defaults.mjs).
//
// It owns the Session 46 Legend rekey (Session 47 keys a Legend's edits by
// row key): `legendEditsByRowKey` and `sessionLegendEditsByRowKey`, twins of
// `legend_edits_by_row_key` and `session_legend_edits_by_row_key`
// (gbdraw/session_io.py), pinned by
// tests/fixtures/sessions/legend-row-key-vectors.json.
import {
  FILE_IMPORT_ORIGIN_SESSION_VERSION,
  MODE_SCOPED_SESSION_VERSION,
  OTHER_MODE_RESULT_EDITOR_FIELDS,
  OTHER_MODE_RESULT_LEGEND_FIELDS,
  OTHER_MODE_RESULT_UI_FIELDS,
  ROW_KEYED_LEGEND_SESSION_VERSION,
  SESSION_DRAWINGS_LAYOUT
} from '../contracts/schema.generated.js';
import { MODE_SCOPED_SETTINGS } from '../mode-scoped-settings.generated.js';
import { MODE_SLICE_DEFAULTS } from '../mode-slice-defaults.generated.js';
import { MODE_SLICE_CONTAINERS } from './mode-scoped-migration.js';
import { legendEntryKey } from './drawing/paint.js';
import { withoutUnarrangedEntries } from './linear-record-layout.js';
import { normalizeDefaultColor } from '../utils/color-utils.js';
import { diagnosticError } from '../utils/error-normalization.js';
import { normalizeLegendColorOverrides, normalizeStrokeOverrideMap } from './style-overrides.js';
/** @import { LegendEditsUnread } from './style-overrides.js' */
/** @import { SessionChecked46 } from '../contracts/session.generated.js' */
/** @import { ResultSetPath } from './session-authority.js' */

/** @param {unknown} value */
const isObject = (value) => value !== null && typeof value === 'object' && !Array.isArray(value);

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

// ---------------------------------------------------------------------------
// The move of a Session 46-49 to drawings: `upgradeModeScopedSession`, the
// twin of `drawings_document_from_mode_scoped` (gbdraw/session_drawings.py),
// function by function as named below.

/** @typedef {'circular' | 'linear'} DrawingMode */
/** @typedef {Record<string, unknown>} Fields */
/**
 * What only the source version says (DESIGN §4.2): Session 48 card rows are
 * file rows when it has no Web draft; Sessions 46-47 keep protein Results in
 * the drawn frame; Sessions 46 and later saved records of Python's paint with
 * their Results (a Session 27-44 saved none, and Python's conversion adds
 * none: Load records them, `createSavedResultPlan`). No other reader compares
 * source versions. `resultsRecordPythonPaint` is the Web's own: Python's twin
 * (`ModeScopedSourceFacts`) reads no Results.
 * @typedef {{ cardRowsAreFileRows: boolean, proteinResultsInDrawnFrame: boolean,
 *   resultsRecordPythonPaint: boolean }} SourceFacts
 */
/**
 * The one reading of a Session file's source facts by its version, for the
 * Session 50 Load reads from it (a Session 46-49's here, a Session 27-44's
 * from Python's conversion, a Session 50 itself).
 * @param {number} version The file's Session version.
 * @param {boolean} hasDraft The file holds a Web draft (a Session 50 holds
 *   each drawing's own, which the drawing's projection reads: pass false).
 *   A Session 27-44's protein Results are in the search frame of the
 *   Session 50 Load reads: Python's conversion moves them (OV-1969).
 * @returns {SourceFacts}
 */
export const sessionSourceFactsOf = (version, hasDraft) => ({
  cardRowsAreFileRows: version >= FILE_IMPORT_ORIGIN_SESSION_VERSION && !hasDraft,
  proteinResultsInDrawnFrame: version >= MODE_SCOPED_SESSION_VERSION && version < FILE_IMPORT_ORIGIN_SESSION_VERSION,
  resultsRecordPythonPaint: version >= MODE_SCOPED_SESSION_VERSION
});
/**
 * One drawing of a Session 50 (`SessionDrawing`, gbdraw/session_contract.py).
 * @typedef {{
 *   id: string,
 *   name: string,
 *   mode: DrawingMode,
 *   renderRequest: Fields | null,
 *   results: unknown[],
 *   runMetadata?: Fields,
 *   cliInvocation?: Fields,
 *   cliOptions?: Fields,
 *   orthogroupState?: Fields,
 *   losatCache?: Fields,
 *   losatDerivedCache?: Fields,
 *   proteinIdentityManifest?: Fields,
 *   config?: Fields,
 *   features?: Fields,
 *   editorState?: Fields,
 *   ui: Fields,
 *   webFiles?: Fields
 * }} SessionDrawing50
 */
/**
 * A Session 50 (`GbdrawSession50`, gbdraw/session_contract.py; DESIGN §1).
 * @typedef {{
 *   format: 'gbdraw-session',
 *   version: number,
 *   createdAt: string,
 *   title?: string,
 *   resources: Fields,
 *   legacyArtifacts?: Fields,
 *   ui: Fields & { activeDrawingId: string },
 *   webFiles?: { resourceOriginalNames?: Fields },
 *   drawings: SessionDrawing50[]
 * }} Session50
 */

const LAYOUT = SESSION_DRAWINGS_LAYOUT;
const DRAWING_MODES = /** @type {readonly DrawingMode[]} */ (Object.keys(LAYOUT.drawingNames));
/** @param {readonly string[]} fields @param {...string} extra */
const fieldSet = (fields, ...extra) => new Set([...fields, ...extra]);
const SET_FIELDS = fieldSet(LAYOUT.setFields);
const SET_EDITOR_FIELDS = fieldSet(LAYOUT.setEditorFields);
const RESULT_EDITOR_FIELDS = fieldSet(OTHER_MODE_RESULT_EDITOR_FIELDS);
const RESULT_LEGEND_FIELDS = fieldSet(OTHER_MODE_RESULT_LEGEND_FIELDS);
const RESULT_UI_FIELDS = fieldSet(OTHER_MODE_RESULT_UI_FIELDS);
const PROJECT_FIELDS = fieldSet(LAYOUT.projectFields);
const PROJECT_UI_FIELDS = fieldSet(LAYOUT.projectUiFields);
const PROJECT_WEB_FILES_FIELDS = fieldSet(LAYOUT.projectWebFilesFields);
const DRAWING_FIELDS = fieldSet(LAYOUT.drawingFields);
const VIEW_FIELDS = fieldSet(LAYOUT.viewFields);
const DRAWING_WEB_FILE_KEYS = fieldSet(LAYOUT.drawingWebFileKeys);
const RESOURCE_REFERENCE_FIELDS = fieldSet(LAYOUT.resourceReferenceFields);
const DRAWING_UI_DROPPED = fieldSet(LAYOUT.projectUiFields, 'mode', 'activeDrawingId', ...Object.values(LAYOUT.inputTypeFields));
const LEGEND = new Set(['legend']);
const SLICE_DRAFT = new Set(['config', 'features']);
const ABSENT = Symbol('absent');
const layoutInvalid = () => diagnosticError('INPUT_INVALID', { field: 'schema', reason: 'FIELDS' });

/** @param {unknown} value @returns {Fields} */
const mapping = (value) => (isObject(value) ? /** @type {Fields} */ (value) : {});
/** @param {Fields} source @param {ReadonlySet<string>} fields @returns {Fields} */
const only = (source, fields) => Object.fromEntries(Object.entries(source).filter(([key]) => fields.has(key)));
/** @param {Fields} source @param {ReadonlySet<string>} fields @returns {Fields} */
const without = (source, fields) => Object.fromEntries(Object.entries(source).filter(([key]) => !fields.has(key)));
/** `_with_key`. @param {Fields} target @param {string} key @param {Fields} value @param {boolean} present */
const withKey = (target, key, value, present) => {
  if (present || Object.keys(value).length) target[key] = { ...value };
};
/** @param {unknown} mode @returns {mode is DrawingMode} */
const isDrawingMode = (mode) => typeof mode === 'string' && Object.hasOwn(LAYOUT.drawingNames, mode);
/** `_request_mode`. @param {Fields} view @returns {DrawingMode | null} */
const requestMode = (view) => {
  const mode = isObject(view.renderRequest) ? mapping(view.renderRequest).mode : null;
  return isDrawingMode(mode) ? mode : null;
};
/** `_view_mode`. @param {Fields} view @returns {DrawingMode} */
const viewMode = (view) => {
  const mode = requestMode(view) || mapping(view.ui).mode;
  if (!isDrawingMode(mode)) throw layoutInvalid();
  return mode;
};

/** `_shared`: a one-drawing view outside its committed Result set. @param {Fields} view @returns {Fields} */
const sharedOf = (view) => {
  const editor = mapping(view.editorState);
  const shared = without(view, new Set([...SET_FIELDS, 'otherModeResult']));
  shared.editorState = {
    ...without(editor, RESULT_EDITOR_FIELDS),
    legend: without(mapping(editor.legend), RESULT_LEGEND_FIELDS)
  };
  shared.ui = without(mapping(view.ui), RESULT_UI_FIELDS);
  return shared;
};

/** `_view`: the shared fields with one committed Result set. @param {Fields} shared @param {Fields} resultSet @returns {Fields} */
const viewOf = (shared, resultSet) => {
  const sharedEditor = mapping(shared.editorState);
  const setEditor = mapping(resultSet.editorState);
  return {
    ...without(shared, new Set(['editorState', 'ui'])),
    ...only(resultSet, SET_FIELDS),
    editorState: {
      ...without(sharedEditor, LEGEND),
      ...only(setEditor, SET_EDITOR_FIELDS),
      legend: { ...mapping(sharedEditor.legend), ...only(mapping(setEditor.legend), RESULT_LEGEND_FIELDS) }
    },
    ui: { ...mapping(shared.ui), ...only(mapping(resultSet.ui), RESULT_UI_FIELDS) }
  };
};

/**
 * Each mode's committed one-drawing view of a Session 46-49, in document
 * order (`read_session_drawings`, the 46-49 branch), with where the file
 * holds its Result set (`at`, the prefix of the set's paths).
 * @param {Fields} source
 * @returns {Map<DrawingMode, { view: Fields, at: ResultSetPath }>}
 */
const committedViews = (source) => {
  /** @type {Map<DrawingMode, { view: Fields, at: ResultSetPath }>} */
  const committed = new Map();
  if (source.renderRequest === null || source.renderRequest === undefined) {
    const shown = mapping(source.ui).mode;
    if (isDrawingMode(shown)) committed.set(shown, { view: source, at: '' });
    return committed;
  }
  const topMode = requestMode(source);
  if (!topMode) return committed;
  committed.set(topMode, { view: without(source, new Set(['otherModeResult'])), at: '' });
  if (isObject(source.otherModeResult)) {
    const view = viewOf(sharedOf(source), mapping(source.otherModeResult));
    const otherMode = requestMode(view);
    if (otherMode) committed.set(otherMode, { view, at: 'otherModeResult.' });
  }
  return committed;
};

/** `_undrawn_view`: the view of `mode` in a Session 46-49 with no Result of it. @param {Fields} source @param {DrawingMode} mode @returns {Fields} */
const undrawnView = (source, mode) => {
  const view = sharedOf(source);
  return {
    ...view,
    renderRequest: null,
    results: [],
    editorState: { ...mapping(view.editorState), featureCatalog: null },
    ui: { ...mapping(view.ui), mode }
  };
};

/** `canonical_resource_ids` (gbdraw/session_resources.py). @param {unknown} value @returns {Set<string>} */
const resourceIdsOf = (value) => {
  /** @type {Set<string>} */
  const found = new Set();
  const pending = [value];
  while (pending.length) {
    const item = pending.pop();
    if (Array.isArray(item)) pending.push(...item);
    else if (isObject(item)) {
      Object.entries(/** @type {Fields} */ (item)).forEach(([key, child]) => {
        if (RESOURCE_REFERENCE_FIELDS.has(key) && typeof child === 'string' && child.trim()) found.add(child);
        else pending.push(child);
      });
    }
  }
  return found;
};

/** `drawing_resource_ids` (gbdraw/session_resources.py). @param {Fields} drawing @returns {Set<string>} */
const drawingResourceIds = (drawing) => {
  const found = new Set(LAYOUT.drawingResourceFields.flatMap((field) => [...resourceIdsOf(drawing[field])]));
  const webFiles = mapping(drawing.webFiles);
  LAYOUT.resourceIdListFields.forEach((field) => {
    const ids = webFiles[field];
    if (Array.isArray(ids)) ids.forEach((id) => { if (typeof id === 'string') found.add(id); });
  });
  return found;
};

/** `_binds_input`, rule c: a binding of `mode`'s input slots names a file. @param {Fields} bindings @param {DrawingMode} mode */
const bindsInput = (bindings, mode) => LAYOUT.bindingModeSlots[mode].some((slot) => resourceIdsOf(bindings[slot]).size > 0);

/** @typedef {{ domain: string, path: string, modes: string, key: string | null }} RegistryRow */

/** `_canonical`: JSON with sorted keys. @param {unknown} value @returns {string} */
const canonical = (value) => JSON.stringify(value, (_key, item) => (
  isObject(item)
    ? Object.fromEntries(Object.entries(/** @type {Fields} */ (item)).sort(([a], [b]) => (a < b ? -1 : a > b ? 1 : 0)))
    : item
));

/** `_row_value`: a row's value in a mode slice, normalized; ABSENT when the slice leaves it out. @param {Fields} modeSlice @param {RegistryRow} row @returns {unknown} */
const rowValue = (modeSlice, row) => {
  /** @type {unknown} */
  let holder = modeSlice;
  for (const part of row.domain.split('.')) holder = isObject(holder) && Object.hasOwn(mapping(holder), part) ? mapping(holder)[part] : ABSENT;
  if (!isObject(holder) || !Object.hasOwn(mapping(holder), row.path)) return ABSENT;
  return withoutUnarrangedEntries(MODE_SLICE_DEFAULTS.normalize, row, mapping(holder)[row.path]);
};

/** `_default_value`: the packaged default of a row in `mode`; ABSENT when Save leaves it out. @param {DrawingMode} mode @param {RegistryRow} row @returns {unknown} */
const defaultValue = (mode, row) => {
  const held = MODE_SLICE_DEFAULTS.defaults[mode].find(({ domain, path }) => domain === row.domain && path === row.path);
  return held && Object.hasOwn(held, 'value') ? held.value : ABSENT;
};

/** `_keyed_entries`: a keyed row's entries, canonical, by entry key. @param {unknown} value @param {string} key @returns {Map<string, string>} */
const keyedEntries = (value, key) => {
  if (isObject(value)) return new Map(Object.entries(mapping(value)).map(([entryKey, entry]) => [entryKey, canonical(entry)]));
  const fields = key.startsWith('JSON[') ? key.slice('JSON['.length, -1).split(',') : [key];
  return new Map((Array.isArray(value) ? value : []).map((entry) => {
    const identity = isObject(entry)
      ? fields.map((field) => (Object.hasOwn(mapping(entry), field) ? mapping(entry)[field] : null))
      : [entry];
    return [canonical(identity), canonical(entry)];
  }));
};

/**
 * A registry row's value in `mode`'s slice when the slice holds it and it
 * differs from `mode`'s packaged default; ABSENT otherwise.
 * @param {Fields} modeSlice
 * @param {DrawingMode} mode
 * @param {RegistryRow} row
 * @returns {unknown}
 */
const changedRowValue = (modeSlice, mode, row) => {
  if (row.modes !== 'both' && row.modes !== mode) return ABSENT;
  const value = rowValue(modeSlice, row);
  if (value === ABSENT) return ABSENT;
  const fallback = defaultValue(mode, row);
  return fallback !== ABSENT && canonical(value) === canonical(fallback) ? ABSENT : value;
};

/**
 * Whether a drawing's slice (as Save writes it, `buildModeSliceData`) holds a
 * registry row that differs from its mode's packaged default: a changed
 * setting or an edit, which makes the drawing one Save writes (DESIGN §4.4).
 * @param {Fields} modeSlice
 * @param {DrawingMode} mode
 */
export const sliceDiffersFromDefault = (modeSlice, mode) => /** @type {readonly RegistryRow[]} */ (
  MODE_SCOPED_SETTINGS.rows
).some((row) => changedRowValue(modeSlice, mode, row) !== ABSENT);

/**
 * `_used_by_settings`, rule d: a registry row of `mode`'s slice that differs
 * from `mode`'s packaged default and, unless the row is `mode`'s alone, from
 * `other`'s value (entry by entry for a keyed row).
 * @param {Fields} slices
 * @param {DrawingMode} mode
 * @param {DrawingMode} other
 */
const usedBySettings = (slices, mode, other) => {
  const own = mapping(slices[mode]);
  const shown = mapping(slices[other]);
  return /** @type {readonly RegistryRow[]} */ (MODE_SCOPED_SETTINGS.rows).some((row) => {
    const value = changedRowValue(own, mode, row);
    if (value === ABSENT) return false;
    if (row.modes === mode) return true;
    const shownValue = rowValue(shown, row);
    if (shownValue === ABSENT) return true;
    if (row.key === null) return canonical(value) !== canonical(shownValue);
    const shownEntries = keyedEntries(shownValue, row.key);
    return [...keyedEntries(value, row.key)].some(([entryKey, entry]) => shownEntries.get(entryKey) !== entry);
  });
};

/**
 * `_migrated_drawing_modes` with its evidence: the modes of the migrated
 * drawings, the drawing that opens first. The mode shown at Save (`ui.mode`,
 * else the first committed mode) always gets a drawing; the other mode gets
 * one only when it is used: a Result, a non-null input binding of its mode
 * (rule c), or a registry row that differs from its default and from the
 * shown mode's (rule d). The shown drawing opens unless it has no stored
 * Result and the other has one: then the drawing with the Result opens
 * (the figure the released app showed), and the drawing that opens comes first.
 * @param {Fields} source
 * @param {Map<DrawingMode, { view: Fields }>} committed
 * @returns {DrawingMode[]}
 */
const migratedDrawingModes = (source, committed) => {
  const shownMode = mapping(source.ui).mode;
  const shown = isDrawingMode(shownMode) ? shownMode : [...committed.keys()][0];
  const other = DRAWING_MODES.find((mode) => mode !== shown) || shown;
  /** @param {DrawingMode} mode */
  const hasResult = (mode) => mapping(committed.get(mode)?.view).renderRequest != null;
  const slices = mapping(source.modes);
  const used = committed.has(other)
    || bindsInput(mapping(mapping(source.webFiles).bindings), other)
    || (Object.hasOwn(slices, other) && usedBySettings(slices, other, shown));
  if (!used) return [shown];
  return hasResult(other) && !hasResult(shown) ? [other, shown] : [shown, other];
};

/** `_project_of_view`: the Session 50 project fields of a view. @param {Fields} view @returns {Fields} */
const projectOf = (view) => {
  const project = { ...only(view, PROJECT_FIELDS), ui: only(mapping(view.ui), PROJECT_UI_FIELDS) };
  withKey(project, 'webFiles', only(mapping(view.webFiles), PROJECT_WEB_FILES_FIELDS), false);
  return project;
};

/**
 * `_drawing_of_view`: one Session 50 drawing from its one-drawing view; it
 * keeps its own mode's slice, input slots and request metadata. Its ID and
 * name are its mode's.
 * @param {Fields} view
 * @returns {SessionDrawing50}
 */
const drawingOf = (view) => {
  const mode = viewMode(view);
  const webFiles = mapping(view.webFiles);
  const unplaced = Object.keys(view).some((key) => !VIEW_FIELDS.has(key))
    || Object.keys(webFiles).some((key) => !DRAWING_WEB_FILE_KEYS.has(key));
  if (unplaced) throw layoutInvalid();
  const drawingSlice = mapping(mapping(view.modes)[mode]);
  /** @type {Fields} */
  const drawing = {
    id: mode,
    name: LAYOUT.drawingNames[mode],
    mode,
    ...only(view, DRAWING_FIELDS),
    ...only(drawingSlice, SLICE_DRAFT)
  };
  const editor = mapping(view.editorState);
  const sliceEditor = mapping(drawingSlice.editorState);
  if (Object.hasOwn(view, 'editorState') || Object.hasOwn(drawingSlice, 'editorState')) {
    const merged = { ...without(editor, LEGEND), ...without(sliceEditor, LEGEND) };
    withKey(merged, 'legend', { ...mapping(editor.legend), ...mapping(sliceEditor.legend) },
      Object.hasOwn(editor, 'legend') || Object.hasOwn(sliceEditor, 'legend'));
    drawing.editorState = merged;
  }
  const viewUi = mapping(view.ui);
  const ui = { ...without(viewUi, DRAWING_UI_DROPPED), ...mapping(drawingSlice.ui) };
  const inputTypeField = LAYOUT.inputTypeFields[mode];
  if (Object.hasOwn(viewUi, inputTypeField)) ui.inputType = viewUi[inputTypeField];
  drawing.ui = ui;
  const drawingWebFiles = only(webFiles, fieldSet(LAYOUT.drawingWebFilesFields[mode]));
  if (Object.hasOwn(webFiles, 'bindings')) {
    drawingWebFiles.bindings = only(
      mapping(webFiles.bindings),
      fieldSet(LAYOUT.bindingCardSlots, 'schema', ...LAYOUT.bindingModeSlots[mode])
    );
  }
  withKey(drawing, 'webFiles', drawingWebFiles, false);
  return /** @type {SessionDrawing50} */ (drawing);
};


/**
 * `_drawing_view` (the inverse of `drawingOf`): one drawing of a Session 50
 * as the one-drawing view a Session 46-49 gives, the project's fields with
 * the drawing's: its Result set and request metadata at the top level, its
 * draft as the slice of its mode (`modes.<mode>`), its input type under its
 * mode's `ui` key. Load projects a drawing from it as it projects the
 * committed set of a Session 46-49.
 * @param {Session50} project
 * @param {SessionDrawing50} drawing
 * @returns {Fields}
 */
export const drawingViewOf = (project, drawing) => {
  const { mode } = drawing;
  const sliceEditorFields = new Set([...MODE_SLICE_CONTAINERS.editorState].filter((field) => field !== 'legend'));
  const sliceLegendFields = MODE_SLICE_CONTAINERS['editorState.legend'];
  const sliceUiFields = MODE_SLICE_CONTAINERS.ui;
  const source = /** @type {Fields} */ (/** @type {unknown} */ (drawing));
  /** @type {Fields} */
  const view = { ...only(/** @type {Fields} */ (/** @type {unknown} */ (project)), PROJECT_FIELDS), ...only(source, DRAWING_FIELDS) };
  const editor = mapping(drawing.editorState);
  const legend = mapping(editor.legend);
  if (Object.hasOwn(drawing, 'editorState')) {
    const setEditor = without(editor, new Set([...sliceEditorFields, 'legend']));
    withKey(setEditor, 'legend', without(legend, sliceLegendFields), Object.hasOwn(editor, 'legend'));
    view.editorState = setEditor;
  }
  const drawingSlice = only(source, SLICE_DRAFT);
  const sliceEditor = only(editor, sliceEditorFields);
  withKey(sliceEditor, 'legend', only(legend, sliceLegendFields), false);
  withKey(drawingSlice, 'editorState', sliceEditor, false);
  const drawingUi = mapping(drawing.ui);
  withKey(drawingSlice, 'ui', only(drawingUi, sliceUiFields), false);
  view.modes = { [mode]: drawingSlice };
  /** @type {Fields} */
  const ui = {
    ...without(mapping(project.ui), new Set(['activeDrawingId'])),
    ...without(drawingUi, new Set([...sliceUiFields, 'inputType'])),
    mode
  };
  if (Object.hasOwn(drawingUi, 'inputType')) ui[LAYOUT.inputTypeFields[mode]] = drawingUi.inputType;
  view.ui = ui;
  withKey(view, 'webFiles', { ...mapping(project.webFiles), ...mapping(drawing.webFiles) }, false);
  return view;
};

/**
 * A Session 50 of one-drawing views (`_drawings_document`): each drawing from
 * its view (`drawingOf`) with its ID and name, in tab order; the views share
 * their project fields by construction (one source, or one Save), so the
 * first gives them. Save writes it; the move of a Session 46-49 builds it.
 * @param {readonly { id: string, name: string, view: Fields }[]} saved
 * @param {string} activeDrawingId
 * @returns {Session50}
 */
export const sessionDocumentOfViews = (saved, activeDrawingId) => {
  const project = projectOf(saved[0].view);
  project.version = LAYOUT.version;
  project.ui = { ...mapping(project.ui), activeDrawingId };
  project.drawings = saved.map(({ id, name, view }) => ({ ...drawingOf(view), id, name }));
  return /** @type {Session50} */ (/** @type {unknown} */ (project));
};

/**
 * A Session 46-49's Session 50, what only the source says, and where the
 * file holds each drawing's saved Result set (`setSources`, by drawing ID;
 * a drawing without a saved set has none), so Load finds a set's reading and
 * keeps its repair paths the file's own.
 * @typedef {{ document: Session50, sourceFacts: SourceFacts, setSources: ReadonlyMap<string, ResultSetPath> }} ModeScopedUpgrade
 */
/**
 * The Session 50 of an admitted Session 46-49, and what only the source says
 * (`drawings_document_from_mode_scoped`, gbdraw/session_drawings.py). Each
 * mode's committed set goes to its mode's drawing; the drawings, their
 * order and the one that opens are `migratedDrawingModes`'s. Each drawing keeps its own mode's slice, input
 * bindings, request metadata, the file cards' bindings and a copy of the
 * LOSAT artifacts; IDs and names are the modes'. A Session 46 gets its Legend
 * edits keyed by row key first (`sessionLegendEditsByRowKey`). Resources no
 * drawing names are left out.
 * @param {SessionChecked46} checked
 * @param {number} version The admitted source version.
 * @returns {ModeScopedUpgrade}
 */
export const upgradeModeScopedSession = (checked, version) => upgradeRowKeyedModeScopedSession(
  /** @type {SessionChecked46} */ (/** @type {unknown} */ (sessionLegendEditsByRowKey(/** @type {Fields} */ (checked), version).session)),
  version
);

/**
 * `upgradeModeScopedSession` of a Session 46-49 whose Legend edits are keyed
 * by row key already: Load reads them so before it plans (its preflight runs
 * `sessionLegendEditsByRowKey` and lists what that drops), and the reading
 * is not run twice.
 * @param {SessionChecked46} checked
 * @param {number} version The admitted source version.
 * @returns {ModeScopedUpgrade}
 */
export const upgradeRowKeyedModeScopedSession = (checked, version) => {
  if (!Number.isInteger(version) || version < MODE_SCOPED_SESSION_VERSION || version >= LAYOUT.version) {
    throw layoutInvalid();
  }
  const source = /** @type {Fields} */ (/** @type {unknown} */ (checked));
  const committed = committedViews(source);
  if (!committed.size) throw layoutInvalid();
  const order = migratedDrawingModes(source, committed);
  const views = order.map((mode) => committed.get(mode)?.view || undrawnView(source, mode));
  const project = sessionDocumentOfViews(views.map((view, index) => ({
    id: order[index], name: LAYOUT.drawingNames[order[index]], view
  })), order[0]);
  const migrated = project.drawings;
  const referenced = new Set(migrated.flatMap((drawing) => [...drawingResourceIds(drawing)]));
  project.resources = only(mapping(project.resources), new Set(Object.keys(mapping(project.resources)).filter((id) => referenced.has(id))));
  const hasDraft = Object.values(mapping(source.modes)).some((value) => isObject(mapping(value).config));
  /** @type {Map<string, ResultSetPath>} */
  const setSources = new Map();
  order.forEach((mode, index) => {
    const held = committed.get(mode);
    if (held && held.view.renderRequest != null) setSources.set(migrated[index].id, held.at);
  });
  return {
    setSources,
    document: /** @type {Session50} */ (/** @type {unknown} */ (project)),
    sourceFacts: sessionSourceFactsOf(version, hasDraft)
  };
};
