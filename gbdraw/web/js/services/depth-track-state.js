// @ts-check
import { resolveColorToHex } from '../utils/color-utils.js';
import { isUnsetValue, positiveNumberOrNull } from '../utils/field-reading.js';
import { normalizeOptionalText } from './track-slot-display.js';

/** @import { FileLike } from './session-resource-backing.js' */
/** @import { FieldLeaf, FieldReading } from '../utils/field-reading.js' */

/**
 * The Depth file slots of one record: one file, or a list indexed by series.
 * @typedef {FileLike | (FileLike | null)[] | null | undefined} DepthFileSlots
 */

/**
 * A Depth series config as a draft or a Session holds it: every field is unvalidated.
 * @typedef {{
 *   label?: unknown, color?: unknown, height?: unknown, large_tick_interval?: unknown,
 *   small_tick_interval?: unknown, tick_font_size?: unknown
 * }} RawDepthTrackConfig
 */

/**
 * A Depth series config after `readDepthTrackConfig`.
 * @typedef {{
 *   label: string, color: string, height: number | null, large_tick_interval: number | null,
 *   small_tick_interval: number | null, tick_font_size: number | null
 * }} DepthTrackConfig
 */

/**
 * A track slot row. Python owns its field set (R7); this module reads
 * `renderer`, `id`, `enabled`, `params`, and `depth_binding_error`.
 * @typedef {Record<string, any>} DepthSlotRow
 */

/**
 * The Depth track defaults of the calling stack. The values are a drawing's
 * stored single-series settings, unvalidated.
 * @typedef {object} DepthTrackDefaults
 * @property {unknown} [depthColor] The first series' color.
 * @property {unknown} [depthHeight]
 * @property {unknown} [largeTickInterval]
 * @property {unknown} [smallTickInterval]
 * @property {unknown} [tickFontSize]
 */

/**
 * A drawing's single-series Depth settings (`adv`), unvalidated.
 * @typedef {{ depth_color?: unknown, depth_height?: unknown, depth_large_tick_interval?: unknown,
 *   depth_small_tick_interval?: unknown, depth_tick_font_size?: unknown }} StoredDepthSettings
 */

/**
 * @typedef {object} ReconcileDepthTracksToFilesOptions
 * @property {DepthFileSlots} [files] Depth file slots of one record (an array or one file).
 * @property {RawDepthTrackConfig[]} [depthTracks] The Depth track configs before the change.
 * @property {number} [targetCount]
 * @property {DepthTrackDefaults} [defaults]
 */

/**
 * @typedef {object} DepthTrackSessionWidthOptions
 * @property {DepthFileSlots[]} [rows] Record-major Depth file rows.
 * @property {RawDepthTrackConfig[]} [depthTracks]
 * @property {DepthSlotRow[]} [slots]
 */

/**
 * @typedef {object} ReindexDepthSlotsOptions
 * @property {DepthSlotRow[]} [slots]
 * @property {number} [removedIndex] The logical series that was removed.
 * @property {number} [activeCount] The series count after the removal.
 * @property {(slot: DepthSlotRow) => boolean} [managedPredicate]
 */

/**
 * @typedef {object} DropInvalidManagedDepthSlotsOptions
 * @property {DepthSlotRow[]} [slots]
 * @property {number} [activeCount] The drawing's series count.
 * @property {Iterable<number>} [sourced] Series indexes with a source; without it every series has one.
 * @property {(slot: DepthSlotRow) => boolean} [managedPredicate]
 */

/**
 * @typedef {object} ReconcileManagedDepthSlotsOptions
 * @property {DepthSlotRow[]} [slots]
 * @property {Iterable<number>} [previousSourced] Series indexes with a source before.
 * @property {Iterable<number>} [sourced] Series indexes with a source after.
 * @property {(slot: DepthSlotRow) => boolean} [managedPredicate]
 */

/**
 * @typedef {object} SyncDepthSlotLabelsOptions
 * @property {DepthSlotRow[]} [slots]
 * @property {readonly RawDepthTrackConfig[]} [depthTracks]
 * @property {number | null} [activeCount]
 */

const DEPTH_TRACK_FALLBACK_COLORS = [
  '#4A90E2',
  '#E45756',
  '#2CA02C',
  '#F28E2B',
  '#9467BD',
  '#8C564B',
  '#17BECF',
  '#7F7F7F'
];


/**
 * The defaults of a Depth series an edit of a drawing gives it settings for
 * (`depthTrackConfigForEdit`): the series color, and the drawing's stored
 * height. Every writer of a drawing's series uses these; Generate leaves out
 * the height (`shapeDepthTracksForGenerate`).
 * @param {StoredDepthSettings} adv
 * @returns {DepthTrackDefaults}
 */
export const depthTrackEditDefaults = (adv) => ({
  depthColor: adv.depth_color,
  depthHeight: adv.depth_height,
  largeTickInterval: null,
  smallTickInterval: null,
  tickFontSize: null
});

/**
 * @param {unknown} [params]
 * @returns {Record<string, unknown>}
 */
const cloneParams = (params = {}) => {
  if (!params || typeof params !== 'object' || Array.isArray(params)) return {};
  return { ...params };
};

/**
 * The name the figure draws for a Depth series without a label: Depth for
 * the only series, Depth N for each of several (Python
 * `_default_depth_track_labels`, gbdraw/analysis/depth_tracks.py). Figure
 * text, so it stays English in every UI locale; the legend-name fields show it
 * as their placeholder.
 * @param {number} index The series' index.
 * @param {number} count The number of series.
 * @returns {string}
 */
export const defaultDepthSeriesLabel = (index, count) => (count === 1 ? 'Depth' : `Depth ${index + 1}`);

// The one rule for the color of a series without its own: the first series
// takes the drawing's Depth color, a later one the palette (OV-1510).
/**
 * @param {number} index
 * @param {DepthTrackDefaults} [defaults]
 * @returns {string}
 */
const fallbackColor = (index, defaults = {}) => {
  if (index === 0 && defaults.depthColor) return String(defaults.depthColor);
  return DEPTH_TRACK_FALLBACK_COLORS[index % DEPTH_TRACK_FALLBACK_COLORS.length];
};

/**
 * The color Depth series `index` of a drawing draws in while it has none of
 * its own: the editor shows it, an edit stores it, and the request sends it.
 * @param {number} index
 * @param {StoredDepthSettings} adv
 * @returns {string}
 */
export const depthSeriesDefaultColor = (index, adv) => fallbackColor(index, depthTrackEditDefaults(adv));

/**
 * The color an imported render request's Depth series takes: its own, else
 * the request's Depth color, in which Python draws every series without one
 * (`clone_depth_config` in gbdraw/analysis/depth_tracks.py). An import redraws
 * what the CLI drew; a series the Web creates takes `depthSeriesDefaultColor`.
 * @param {unknown} color The series' color in the request: text, or unset.
 * @param {string} depthColor The request's effective Depth color.
 * @returns {string}
 */
export const importedDepthSeriesColor = (color, depthColor) => (
  typeof color === 'string' && color ? color : depthColor
);

/** @param {FileLike | null | undefined} file */
const depthFileBaseName = (file) => {
  // split always returns at least one element, so pop is not undefined.
  const rawName = /** @type {string} */ (String(file?.name || '').split(/[\\/]/).pop()).trim();
  if (!rawName) return '';
  return rawName.replace(/\.[^.]+$/, '').trim() || rawName;
};

/**
 * @param {unknown} label
 * @param {FileLike | null | undefined} file
 */
const depthLabelMatchesFile = (label, file) => {
  const labelText = String(label || '').toLowerCase();
  const fileText = depthFileBaseName(file).toLowerCase();
  if (!labelText || !fileText) return false;
  if (labelText.includes(fileText) || fileText.includes(labelText)) return true;
  const accessionTokens = fileText.match(/[a-z]{2,5}\d{5,}/g) || [];
  return accessionTokens.some((token) => labelText.includes(token));
};

/**
 * @param {DepthFileSlots} value
 * @returns {(FileLike | null)[]}
 */
export const depthFileSlotsFromValue = (value) => {
  if (Array.isArray(value)) return value.slice();
  return value ? [value] : [];
};

/**
 * @param {unknown} value
 * @returns {value is DepthFileSlots[]}
 */
export const isRecordMajorDepthFileMatrix = (value) => (
  Array.isArray(value) && value.every((row) => Array.isArray(row))
);

/**
 * @param {unknown} value
 * @param {string} [fieldName]
 * @returns {number}
 */
export const parseDepthTrackIndexIdentity = (value, fieldName = 'track_index') => {
  /** @type {number | null} */
  let numeric = null;
  if (typeof value === 'number') {
    numeric = value;
  } else if (
    typeof value === 'string' &&
    /^[+-]?\d+$/.test(value.trim())
  ) {
    numeric = Number(value.trim());
  }
  if (numeric === null || !Number.isSafeInteger(numeric) || numeric < 0) {
    throw new Error(`${fieldName} must be a non-negative integer.`);
  }
  return numeric;
};

/**
 * @param {DepthFileSlots | DepthFileSlots[]} value
 * @param {unknown} [recordCount]
 * @returns {(FileLike | null)[][]}
 */
export const normalizeRecordMajorDepthFileRows = (value, recordCount = 1) => {
  const requestedCount = Math.max(1, Number(recordCount) || 1);
  const sourceIsMatrix = isRecordMajorDepthFileMatrix(value);
  const rows = sourceIsMatrix
    ? value.map((row) => depthFileSlotsFromValue(row))
    : Array.from(
        { length: requestedCount },
        () => depthFileSlotsFromValue(value)
      );
  const targetCount = Math.max(requestedCount, rows.length);
  while (rows.length < targetCount) rows.push([]);
  const width = depthTrackMatrixWidth(rows);
  return rows.map((row) => padDepthFileSlots(row, width));
};

/**
 * @param {DepthFileSlots | DepthFileSlots[]} value
 * @returns {(FileLike | null)[]}
 */
export const representativeDepthFiles = (value) => {
  const rows = isRecordMajorDepthFileMatrix(value)
    ? normalizeRecordMajorDepthFileRows(value, value.length)
    : [depthFileSlotsFromValue(value)];
  const width = depthTrackMatrixWidth(rows);
  return Array.from({ length: width }, (_, trackIndex) => (
    rows.map((row) => row[trackIndex]).find(Boolean) || null
  ));
};

/** @param {DepthFileSlots} slots */
export const compactDepthFileSlots = (slots) => {
  const next = depthFileSlotsFromValue(slots);
  while (next.length > 0 && !next[next.length - 1]) {
    next.pop();
  }
  return next;
};

/** @param {DepthFileSlots} value */
export const uploadedDepthFileCount = (value) => (
  depthFileSlotsFromValue(value).filter(Boolean).length
);

/** @param {DepthFileSlots[] | null | undefined} rows */
export const depthTrackMatrixWidth = (rows) => (
  (Array.isArray(rows) ? rows : []).reduce(
    (width, row) => Math.max(width, depthFileSlotsFromValue(row).length),
    0
  )
);

/** @param {DepthFileSlots[] | null | undefined} rows */
export const activeDepthTrackIndices = (rows) => {
  /** @type {Set<number>} */
  const active = new Set();
  (Array.isArray(rows) ? rows : []).forEach((row) => {
    depthFileSlotsFromValue(row).forEach((file, trackIndex) => {
      if (file) active.add(trackIndex);
    });
  });
  return Array.from(active).sort((left, right) => left - right);
};

/**
 * @param {DepthFileSlots[] | null | undefined} rows
 * @param {unknown} trackIndex
 */
export const depthTrackCoverageCount = (rows, trackIndex) => {
  const idx = Number(trackIndex);
  if (!Number.isInteger(idx) || idx < 0) return 0;
  return (Array.isArray(rows) ? rows : []).reduce(
    (count, row) => count + (depthFileSlotsFromValue(row)[idx] ? 1 : 0),
    0
  );
};

/**
 * @param {DepthFileSlots} value
 * @param {unknown} width
 */
export const padDepthFileSlots = (value, width) => {
  const next = depthFileSlotsFromValue(value);
  const targetWidth = Math.max(next.length, Math.max(0, Number(width) || 0));
  while (next.length < targetWidth) next.push(null);
  return next;
};

/**
 * @param {DepthFileSlots} value
 * @param {unknown} trackIndex
 * @param {number | null} [width]
 */
export const clearDepthTrackSourceAt = (value, trackIndex, width = null) => {
  const idx = Number(trackIndex);
  if (!Number.isInteger(idx) || idx < 0) return depthFileSlotsFromValue(value);
  const targetWidth = width === null || width === undefined
    ? idx + 1
    : Math.max(idx + 1, Number(width) || 0);
  const next = padDepthFileSlots(value, targetWidth);
  next[idx] = null;
  return next;
};

/**
 * @param {DepthFileSlots[] | null | undefined} rows
 * @param {unknown} trackIndex
 */
export const removeDepthTrackColumnAt = (rows, trackIndex) => {
  const idx = Number(trackIndex);
  const sourceRows = Array.isArray(rows) ? rows : [];
  if (!Number.isInteger(idx) || idx < 0) {
    return sourceRows.map((row) => depthFileSlotsFromValue(row));
  }
  return sourceRows.map((row) => {
    const next = depthFileSlotsFromValue(row);
    if (idx < next.length) next.splice(idx, 1);
    return next;
  });
};

// The one reader of a Depth series config: a stored series, a Session's, and
// the defaults of a series the app adds. Unset fields take `defaults`; a
// color name reads as its hex. `lost`: a stored value it could not use.
/**
 * @param {unknown} entry
 * @param {number} index
 * @param {DepthTrackDefaults} [defaults]
 * @returns {FieldReading<DepthTrackConfig>}
 */
export const readDepthTrackConfig = (entry, index, defaults = {}) => {
  const isEntry = Boolean(entry) && typeof entry === 'object' && !Array.isArray(entry);
  /** @type {RawDepthTrackConfig} */
  const source = isEntry ? /** @type {RawDepthTrackConfig} */ (entry) : {};
  const lost = !isEntry && !isUnsetValue(entry);
  /** @type {FieldLeaf[]} */
  const leaves = [];
  /** @param {string} field @param {unknown} value */
  const positive = (field, value) => {
    const reading = positiveNumberOrNull(value);
    if (reading.lost) leaves.push({ at: `.${field}`, stored: value, used: reading.value });
    return reading.value;
  };
  const value = {
    // A series nobody named keeps no label: its name follows the series count (OV-1508).
    label: String(source.label ?? ''),
    color: resolveColorToHex(String(source.color || fallbackColor(index, defaults))),
    height: positive('height', Object.hasOwn(source, 'height') ? source.height : defaults.depthHeight),
    large_tick_interval: positive('large_tick_interval', source.large_tick_interval ?? defaults.largeTickInterval),
    small_tick_interval: positive('small_tick_interval', source.small_tick_interval ?? defaults.smallTickInterval),
    tick_font_size: positive('tick_font_size', source.tick_font_size ?? defaults.tickFontSize)
  };
  return { value, lost: lost || leaves.length > 0, ...(leaves.length ? { leaves } : {}) };
};

/**
 * @param {unknown} entry
 * @param {number} index
 * @param {DepthTrackDefaults} [defaults]
 * @returns {DepthTrackConfig}
 */
export const normalizeDepthTrackConfig = (entry, index, defaults = {}) => (
  readDepthTrackConfig(entry, index, defaults).value
);

// A drawing's stored single-series Depth settings belong to its first series;
// every series takes the height.
/**
 * @param {StoredDepthSettings} adv
 * @param {number} index
 * @returns {DepthTrackDefaults}
 */
const storedSeriesDefaults = (adv, index) => ({
  depthColor: adv.depth_color,
  depthHeight: adv.depth_height,
  ...(index === 0 ? {
    largeTickInterval: adv.depth_large_tick_interval,
    smallTickInterval: adv.depth_small_tick_interval,
    tickFontSize: adv.depth_tick_font_size
  } : {})
});

/**
 * A drawing's stored Depth series (`adv.depth_tracks`), at least one.
 * @param {unknown} tracks
 * @param {StoredDepthSettings} [adv]
 * @returns {FieldReading<DepthTrackConfig[]>}
 */
export const readDepthTracks = (tracks, adv = {}) => {
  const readings = (Array.isArray(tracks) ? tracks : [])
    .map((entry, index) => readDepthTrackConfig(entry, index, storedSeriesDefaults(adv, index)));
  if (!readings.length) readings.push(readDepthTrackConfig(null, 0, storedSeriesDefaults(adv, 0)));
  const notList = !Array.isArray(tracks) && !isUnsetValue(tracks);
  const value = readings.map((reading) => reading.value);
  // An entry that is no series is replaced by a default series.
  const leaves = notList ? [] : readings.flatMap((reading, index) => (!reading.lost ? []
    : reading.leaves ? reading.leaves.map((leaf) => ({ ...leaf, at: `[${index}]${leaf.at}` }))
      : [{ at: `[${index}]`, stored: /** @type {unknown[]} */ (tracks)[index], used: reading.value }]));
  return { value, lost: notList || leaves.length > 0, ...(leaves.length ? { leaves } : {}) };
};

/**
 * @param {RawDepthTrackConfig[] | null | undefined} tracks
 * @param {unknown} count
 * @param {DepthTrackDefaults} [defaults]
 */
export const ensureDepthTrackConfigCount = (tracks, count, defaults = {}) => {
  const rawTracks = Array.isArray(tracks) ? tracks : [];
  const targetCount = Math.max(1, Number(count) || 1);
  const normalized = rawTracks.map((entry, index) => (
    normalizeDepthTrackConfig(entry, index, defaults)
  ));
  while (normalized.length < targetCount) {
    normalized.push(normalizeDepthTrackConfig(null, normalized.length, defaults));
  }
  return normalized;
};

/**
 * @param {RawDepthTrackConfig[] | null | undefined} tracks
 * @param {unknown} count
 * @param {DepthTrackDefaults} [defaults]
 * @returns {RawDepthTrackConfig[]}
 */
export const ensureDepthTrackConfigShape = (tracks, count, defaults = {}) => {
  const targetTracks = Array.isArray(tracks) ? tracks : [];
  const targetCount = Math.max(1, Number(count) || 1);
  for (let index = 0; index < targetCount; index += 1) {
    const current = targetTracks[index];
    if (!current || typeof current !== 'object' || Array.isArray(current)) {
      targetTracks[index] = normalizeDepthTrackConfig(null, index, defaults);
      continue;
    }
    const fallback = normalizeDepthTrackConfig(null, index, defaults);
    if (current.label === undefined || current.label === null) current.label = fallback.label;
    if (current.color === undefined || current.color === null) current.color = fallback.color;
    if (current.height === undefined) current.height = fallback.height;
    if (current.large_tick_interval === undefined) current.large_tick_interval = fallback.large_tick_interval;
    if (current.small_tick_interval === undefined) current.small_tick_interval = fallback.small_tick_interval;
    if (current.tick_font_size === undefined) current.tick_font_size = fallback.tick_font_size;
  }
  return targetTracks;
};

// Generate gives series 0..count-1 of a drawing the settings an edit gives
// them (`depthTrackEditDefaults`), except a height: a series without one
// keeps none (Auto). A Linear Depth row draws its series' stored height, so a
// height Generate stored would replace the row's own on the next Generate
// (OV-1510). Generate adds no series past `count`; a count of 0 shapes none.
/**
 * @param {RawDepthTrackConfig[]} tracks
 * @param {number} count
 * @param {StoredDepthSettings} adv
 * @returns {RawDepthTrackConfig[]}
 */
export const shapeDepthTracksForGenerate = (tracks, count, adv) => (
  count > 0
    ? ensureDepthTrackConfigShape(tracks, count, { ...depthTrackEditDefaults(adv), depthHeight: null })
    : tracks
);

// The settings of Depth series `index` for an edit (TK-03). An edit writes a
// series the drawing's mode has: `seriesCount` series, at least one. An index
// past them names no series, so there is nothing to write and the series list
// does not grow; the series up to `index` get their default fields.
/**
 * @param {RawDepthTrackConfig[]} tracks
 * @param {number} index
 * @param {number} seriesCount
 * @param {DepthTrackDefaults} [defaults]
 * @returns {RawDepthTrackConfig | null}
 */
export const depthTrackConfigForEdit = (tracks, index, seriesCount, defaults = {}) => {
  const idx = Math.max(0, Number(index) || 0);
  if (!Array.isArray(tracks) || idx >= Math.max(1, Number(seriesCount) || 0)) return null;
  ensureDepthTrackConfigShape(tracks, idx + 1, defaults);
  return tracks[idx];
};

// Each series keeps the config of its index, or a later config whose label
// names its file; a series whose config another series took takes the first
// config no series took, so no config is lost while there are as many series
// as configs. `sources`: the config index each series took.
/**
 * @param {ReconcileDepthTracksToFilesOptions} [options]
 * @returns {{ tracks: DepthTrackConfig[], sources: number[] }}
 */
export const reconcileDepthTracksToFiles = ({
  files,
  depthTracks,
  targetCount,
  defaults = {}
} = {}) => {
  const fileSlots = Array.isArray(files) ? files : depthFileSlotsFromValue(files);
  const rawTracks = Array.isArray(depthTracks) ? depthTracks : [];
  /** @type {Set<number>} */
  const usedTrackIndexes = new Set();
  const sources = Array.from({ length: Math.max(1, Number(targetCount) || 1) }, (_, index) => {
    const file = fileSlots[index] || null;
    let sourceIndex = index;
    if (
      file &&
      rawTracks[index] &&
      !depthLabelMatchesFile(rawTracks[index]?.label, file)
    ) {
      const matchingIndex = rawTracks.findIndex((entry, candidateIndex) => (
        candidateIndex > index &&
        !usedTrackIndexes.has(candidateIndex) &&
        depthLabelMatchesFile(entry?.label, file)
      ));
      if (matchingIndex >= 0) sourceIndex = matchingIndex;
    }
    if (usedTrackIndexes.has(sourceIndex)) {
      const unused = rawTracks.findIndex((_, candidateIndex) => !usedTrackIndexes.has(candidateIndex));
      if (unused >= 0) sourceIndex = unused;
    }
    usedTrackIndexes.add(sourceIndex);
    return sourceIndex;
  });
  return {
    tracks: sources.map((sourceIndex, index) => normalizeDepthTrackConfig(rawTracks[sourceIndex], index, defaults)),
    sources
  };
};

// The logical series a Depth row references: its track_index, else its
// depth_<n> ID, else 0 (the Python default). A row's position in the stack
// never selects a series. A row whose series was removed references none.
/**
 * @param {DepthSlotRow | null | undefined} slot
 * @returns {number | null}
 */
export const depthSlotTrackIndex = (slot) => {
  const rawTrackIndex = slot?.params?.track_index;
  if (rawTrackIndex !== null && rawTrackIndex !== undefined && rawTrackIndex !== '') {
    const parsed = Number(rawTrackIndex);
    if (Number.isInteger(parsed) && parsed >= 0) return parsed;
  }
  if (slot?.depth_binding_error) return null;
  const idMatch = String(slot?.id || '').trim().match(/^depth_(\d+)$/);
  if (idMatch) return Math.max(0, Number(idMatch[1]) - 1);
  return 0;
};

/** @param {DepthSlotRow[] | null | undefined} slots */
const referencedDepthTrackWidth = (slots) => {
  let width = 0;
  (Array.isArray(slots) ? slots : []).forEach((slot) => {
    if (!slot || String(slot.renderer || '') !== 'depth') return;
    const trackIndex = depthSlotTrackIndex(slot);
    if (trackIndex !== null) width = Math.max(width, trackIndex + 1);
  });
  return width;
};

/** @param {DepthTrackSessionWidthOptions} [options] */
export const depthTrackSessionWidth = ({ rows, depthTracks, slots } = {}) => Math.max(
  depthTrackMatrixWidth(rows),
  Array.isArray(depthTracks) ? depthTracks.length : 0,
  referencedDepthTrackWidth(slots)
);

/** @param {DepthSlotRow} slot */
const slotHasGeometry = (slot) => (
  normalizeOptionalText(slot?.width) !== null ||
  normalizeOptionalText(slot?.radius) !== null ||
  normalizeOptionalText(slot?.spacing) !== null ||
  normalizeOptionalText(slot?.inner_gap_px) !== null ||
  normalizeOptionalText(slot?.outer_gap_px) !== null ||
  Number(slot?.z || 0) !== 0
);

/**
 * @param {unknown} params
 * @param {string[]} allowedKeys
 */
const paramsOnlyAllowedKeys = (params, allowedKeys) => {
  const allowed = new Set(allowedKeys);
  return Object.entries(cloneParams(params)).every(([key, value]) => (
    allowed.has(String(key)) || normalizeOptionalText(value) === null
  ));
};

/** @param {DepthSlotRow | null | undefined} slot */
export const isDefaultManagedDepthSlot = (slot) => {
  if (!slot || typeof slot !== 'object' || Array.isArray(slot)) return false;
  if (String(slot.renderer || '') !== 'depth') return false;
  if (slotHasGeometry(slot)) return false;
  const id = String(slot.id || '').trim();
  if (!/^depth(?:_\d+)?$/.test(id)) return false;
  return paramsOnlyAllowedKeys(slot.params, ['track_index', 'legend_label', 'managed']);
};

/**
 * @typedef {DepthSlotRow & { params: Record<string, unknown> }} ClonedDepthSlot
 */

/**
 * @param {DepthSlotRow} slot
 * @returns {ClonedDepthSlot}
 */
const cloneSlot = (slot) => ({
  ...slot,
  params: cloneParams(slot?.params)
});

/**
 * @param {DepthSlotRow} slot
 * @param {number | null} [removedTrackIndex]
 */
const disableInvalidManualDepthSlot = (slot, removedTrackIndex = null) => {
  const next = cloneSlot(slot);
  next.enabled = false;
  if (/^depth(?:_\d+)?$/.test(String(next.id || '').trim())) {
    next.id = `disabled_${next.id}`;
  }
  const indexText = Number.isInteger(Number(removedTrackIndex))
    ? ` ${Number(removedTrackIndex)}`
    : '';
  next.depth_binding_error = `Depth logical track index${indexText} is no longer available.`;
  delete next.params.track_index;
  delete next.params.legend_label;
  return next;
};

/** @param {ReindexDepthSlotsOptions} [options] */
export const reindexDepthSlots = ({
  slots,
  removedIndex,
  activeCount,
  managedPredicate = isDefaultManagedDepthSlot
} = {}) => {
  const idx = Number(removedIndex);
  const count = Math.max(0, Number(activeCount) || 0);
  if (!Array.isArray(slots) || !Number.isInteger(idx) || idx < 0) {
    return Array.isArray(slots) ? slots.map((slot) => cloneSlot(slot)) : [];
  }
  /** @type {ClonedDepthSlot[]} */
  const nextSlots = [];
  slots.forEach((slot) => {
    if (!slot || String(slot.renderer || '') !== 'depth') {
      nextSlots.push(cloneSlot(slot));
      return;
    }
    const oldTrackIndex = depthSlotTrackIndex(slot);
    const isManaged = managedPredicate(slot);
    if (oldTrackIndex === null || oldTrackIndex === idx) {
      if (!isManaged) nextSlots.push(disableInvalidManualDepthSlot(slot, oldTrackIndex ?? idx));
      return;
    }
    const newTrackIndex = oldTrackIndex > idx ? oldTrackIndex - 1 : oldTrackIndex;
    if (newTrackIndex < 0 || newTrackIndex >= count) {
      if (!isManaged) nextSlots.push(disableInvalidManualDepthSlot(slot, oldTrackIndex));
      return;
    }
    const next = cloneSlot(slot);
    next.params.track_index = newTrackIndex;
    nextSlots.push(next);
  });
  return nextSlots;
};

/**
 * The ids of the Depth rows a Depth source change takes out of what the stack
 * can draw: a row on before and not on after, and a row of either state that
 * has gone. The change itself builds `after`, so this restates no rule of it.
 * @param {readonly DepthSlotRow[]} before
 * @param {readonly DepthSlotRow[]} after
 * @returns {Set<string>}
 */
export const depthRowsLeaving = (before, after) => {
  const ids = new Set(after.map((slot) => String(slot.id)));
  const on = new Set(after.filter((slot) => slot.enabled !== false).map((slot) => String(slot.id)));
  return new Set(before
    .filter((slot) => isDepthSlot(slot) && !on.has(String(slot.id)) && (slot.enabled !== false || !ids.has(String(slot.id))))
    .map((slot) => String(slot.id)));
};

/**
 * The rows the removal of Depth series `removedIndex` takes out of a stack,
 * `activeCount` series staying (`reindexDepthSlots`).
 * @param {readonly DepthSlotRow[]} slots
 * @param {number} removedIndex
 * @param {number} activeCount
 * @returns {Set<string>}
 */
export const depthSeriesRowsLeaving = (slots, removedIndex, activeCount) => depthRowsLeaving(
  slots, reindexDepthSlots({ slots: [...slots], removedIndex, activeCount })
);

/**
 * A Depth row that names none of the drawing's `count` series.
 * @param {DepthSlotRow} slot
 * @param {number} count
 */
const namesNoDepthSeries = (slot, count) => {
  if (!slot || String(slot.renderer || '') !== 'depth') return false;
  const trackIndex = depthSlotTrackIndex(slot);
  return trackIndex === null || trackIndex < 0 || trackIndex >= count;
};

/** @param {DepthSlotRow} slot */
const isDepthSlot = (slot) => Boolean(slot) && String(slot.renderer || '') === 'depth';

// What a Depth row keeps when its series has no Depth source (PD-OI-058,
// PD-OI-083, BA-OD94), the one rule of a Depth source change and of Load: a
// managed row goes with the source; a manual row stays as the user set it, and
// its row issue says that the series has no source.
/**
 * @param {DepthSlotRow} slot
 * @param {(slot: DepthSlotRow) => boolean} managedPredicate
 */
const keepsDepthRowWithoutSource = (slot, managedPredicate) => !managedPredicate(slot);

// A stack's Depth rows fit to the drawing's `activeCount` series: a row of a
// series without a source keeps or loses it by `keepsDepthRowWithoutSource`;
// a managed row that names no series is dropped and a manual one turned off.
/** @param {DropInvalidManagedDepthSlotsOptions} [options] */
export const dropInvalidManagedDepthSlots = ({
  slots,
  activeCount,
  sourced,
  managedPredicate = isDefaultManagedDepthSlot
} = {}) => {
  const count = Math.max(0, Number(activeCount) || 0);
  /** @type {Set<number | null> | null} */
  const withSource = sourced ? new Set(sourced) : null;
  if (!Array.isArray(slots)) return [];
  return /** @type {ClonedDepthSlot[]} */ (slots
    .map((slot) => {
      if (!isDepthSlot(slot)) return cloneSlot(slot);
      const trackIndex = depthSlotTrackIndex(slot);
      if (namesNoDepthSeries(slot, count)) {
        return managedPredicate(slot) ? null : disableInvalidManualDepthSlot(slot, trackIndex);
      }
      if (withSource && !withSource.has(trackIndex) && !keepsDepthRowWithoutSource(slot, managedPredicate)) return null;
      const next = cloneSlot(slot);
      next.params.track_index = trackIndex;
      return next;
    })
    .filter(Boolean)); // drops the null of a managed row
};

/**
 * Managed Depth row lifecycle (PD-OI-058), shared by both modes. Only the
 * transitions that change Depth sources call it, with the logical series
 * indexes that had a source before and after the change:
 * - a series that loses its last source loses its managed rows;
 * - a series that gains its first source gets one managed row, returned in
 *   `additions`, unless a row (enabled or disabled) already references it.
 * Manual rows are never changed. A manual row left on a series without a
 * source reports a row issue (PD-OI-083).
 * @param {ReconcileManagedDepthSlotsOptions} [options]
 */
export const reconcileManagedDepthSlots = ({
  slots,
  previousSourced = [],
  sourced = [],
  managedPredicate = isDefaultManagedDepthSlot
} = {}) => {
  // depthSlotTrackIndex returns null for a row bound to no series, so the lookups below may ask for null.
  /** @type {Set<number | null>} */
  const before = new Set(previousSourced);
  /** @type {Set<number | null>} */
  const after = new Set(sourced);
  const kept = (Array.isArray(slots) ? slots : []).filter((slot) => {
    if (!isDepthSlot(slot) || keepsDepthRowWithoutSource(slot, managedPredicate)) return true;
    const trackIndex = depthSlotTrackIndex(slot);
    return !before.has(trackIndex) || after.has(trackIndex);
  });
  const referenced = new Set(kept.filter(isDepthSlot).map((slot) => depthSlotTrackIndex(slot)));
  // `after` is built from `sourced`, which holds only numbers; null is only looked up in it.
  const additions = /** @type {number[]} */ (Array.from(after))
    .filter((trackIndex) => !before.has(trackIndex) && !referenced.has(trackIndex))
    .sort((left, right) => left - right);
  return { slots: kept, additions };
};

/** @param {SyncDepthSlotLabelsOptions} [options] */
export const syncDepthSlotLabels = ({ slots, depthTracks, activeCount = null } = {}) => {
  if (!Array.isArray(slots)) return;
  const count = activeCount === null || activeCount === undefined
    ? (Array.isArray(depthTracks) ? depthTracks.length : 0)
    : Math.max(0, Number(activeCount) || 0);
  slots.forEach((slot) => {
    if (!slot || String(slot.renderer || '') !== 'depth') return;
    const trackIndex = depthSlotTrackIndex(slot);
    if (trackIndex === null || trackIndex < 0 || trackIndex >= count) {
      return;
    }
    const label = String(depthTracks?.[trackIndex]?.label ?? '').trim();
    slot.params = cloneParams(slot.params);
    slot.params.track_index = trackIndex;
    if (label) {
      slot.params.legend_label = label;
    } else {
      delete slot.params.legend_label;
    }
  });
};

/**
 * The Legend captions the draft's Depth series would draw: each series label (the
 * default name for a series without one), and the legend label of each Depth row. The request carries no Depth source while
 * Show Depth is off, so Python cannot name these rows then (OV-81).
 * @param {{ depthTracks?: RawDepthTrackConfig[], slots?: DepthSlotRow[] }} [options]
 * @returns {string[]}
 */
export const depthSeriesLegendCaptions = ({ depthTracks, slots } = {}) => {
  const captions = [
    ...(Array.isArray(depthTracks) ? depthTracks : []).map((track, index, tracks) => (
      String(track?.label ?? '').trim() || defaultDepthSeriesLabel(index, tracks.length))),
    ...(Array.isArray(slots) ? slots : [])
      .filter((slot) => String(slot?.renderer || '') === 'depth')
      .map((slot) => slot?.params?.legend_label)
  ].map((caption) => String(caption ?? '').trim());
  return [...new Set(captions.filter(Boolean))];
};

/**
 * A drawing's Depth settings and track stacks.
 * @typedef {StoredDepthSettings & { depth_tracks?: unknown, circular_track_slots?: DepthSlotRow[],
 *   linear_track_slots?: DepthSlotRow[] }} DrawingDepthState
 */
/**
 * @typedef {object} DepthSourceFit
 * @property {DepthTrackConfig[]} depthTracks
 * @property {ClonedDepthSlot[]} circularSlots Before the stack's placements.
 * @property {ClonedDepthSlot[]} linearSlots Before the stack's placements.
 * @property {number} linearRowWidth The width each Linear record's Depth slots
 *   are padded to; 0 leaves them as saved.
 * @property {{ list: 'circular' | 'linear', id: string, enabled: boolean }[]} disabledRows
 *   The manual Depth rows that name no series of the drawing, which the fit
 *   turns off, and whether each was on (a row the fit marked before is not).
 */

/**
 * A loaded drawing's Depth series and Depth rows fit to the Depth sources
 * (plan 4.2 `depth`) by the rule the app keeps when a source goes (BA-OD94):
 * every stored series stays, with or without a Depth file (Circular adds a
 * series for each Circular Depth file past them, and a series takes the
 * stored config whose label names its file; Linear grows to the widest record
 * and pads the records that hold a file, OV-337); a Depth row keeps or loses
 * a series without a source by `keepsDepthRowWithoutSource`. Pure: Load lists
 * what it turns off before the commit (B45) and installs it after.
 * @param {{ mode: 'circular' | 'linear', adv: DrawingDepthState, circularFiles: DepthFileSlots | DepthFileSlots[],
 *   linearRows: (FileLike | null)[][] }} input
 * @returns {DepthSourceFit}
 */
export const fitDepthStateToSources = ({ mode, adv, circularFiles, linearRows }) => {
  const circularDepthFiles = representativeDepthFiles(circularFiles);
  const circularSlots = Array.isArray(adv.circular_track_slots) ? adv.circular_track_slots : [];
  const linearSlots = Array.isArray(adv.linear_track_slots) ? adv.linear_track_slots : [];
  const stored = readDepthTracks(adv.depth_tracks, adv).value;
  const linearDepthCount = depthTrackSessionWidth({ rows: linearRows, depthTracks: stored, slots: linearSlots });
  /** @type {DepthTrackConfig[]} */
  let depthTracks;
  if (mode === 'linear') {
    depthTracks = stored;
    while (depthTracks.length < linearDepthCount) {
      depthTracks.push(readDepthTrackConfig(null, depthTracks.length, storedSeriesDefaults(adv, depthTracks.length)).value);
    }
  } else {
    // Every Circular series takes the first series' stored settings as its
    // defaults, as this fit always read them.
    depthTracks = reconcileDepthTracksToFiles({
      files: circularDepthFiles,
      depthTracks: stored,
      targetCount: Math.max(stored.length, circularDepthFiles.some(Boolean) ? circularDepthFiles.length : 0),
      defaults: storedSeriesDefaults(adv, 0)
    }).tracks;
  }
  /**
   * @param {'circular' | 'linear'} list
   * @param {DepthSlotRow[]} slots
   * @param {number} seriesCount
   * @param {number[]} sourced
   */
  const fitStack = (list, slots, seriesCount, sourced) => {
    const fitted = dropInvalidManagedDepthSlots({ slots, activeCount: seriesCount, sourced });
    syncDepthSlotLabels({ slots: fitted, depthTracks, activeCount: seriesCount });
    const disabled = slots
      .filter((slot) => namesNoDepthSeries(slot, seriesCount) && !isDefaultManagedDepthSlot(slot) && !slot.depth_binding_error)
      .map((slot) => ({ list, id: String(slot.id || ''), enabled: slot.enabled !== false }));
    return { fitted, disabled };
  };
  // Each stack fits the series of its own mode's sources.
  const circular = fitStack('circular', circularSlots, mode === 'circular' ? depthTracks.length : Math.max(stored.length, circularDepthFiles.length),
    circularDepthFiles.flatMap((file, index) => (file ? [index] : [])));
  const linear = fitStack('linear', linearSlots, mode === 'linear' ? depthTracks.length : linearDepthCount, activeDepthTrackIndices(linearRows));
  return {
    depthTracks,
    circularSlots: circular.fitted,
    linearSlots: linear.fitted,
    linearRowWidth: mode === 'linear' && linearDepthCount > 0 && linearRows.some((row) => row.some(Boolean)) ? linearDepthCount : 0,
    disabledRows: [...circular.disabled, ...linear.disabled]
  };
};
