// @ts-check
import { resolveColorToHex } from '../utils/color-utils.js';
import {
  exactChoice, isNonNegativeInteger, isPositiveInteger, isUnsetValue, nonNegativeNumberOrNull, numberOr, onlyTrue,
  positiveNumberOrNull, textChoice, trueUnlessFalse, truthy
} from '../utils/field-reading.js';
import { loadRepair } from '../utils/load-repairs.js';
import { IMPORTED_COMPARISON_ACTIONS } from './imported-comparison-intent.js';
import { readCircularPlotTitlePosition } from './layout-preferences.js';
import { readCollinearAnchorMode, readCollinearSearchScope, readOrthogroupMembershipMode } from './losat-normalization.js';
import { WEB_UX_PROFILE } from '../web-ux-profile.js';

/** @import { FieldReader, FieldReading } from '../utils/field-reading.js' */
/** @import { LoadRepair } from '../utils/load-repairs.js' */

// The readers of the stored settings of a configuration (a drawing's draft, a
// legacy configuration file) and of a Session's app `ui`. `applyConfigData`
// and Load's apply read with these rows, History restores included; before
// Load applies a file it lists what they cannot restore (BA-D25, B36, B37,
// B43), so the apply makes no repair of its own.

/**
 * One stored field and its reader. A reader may read fields its table read
 * before it (`container`). `inPlace`: the list is replaced in its array.
 * `retired`: the field is no longer read and is deleted.
 * @typedef {object} FieldRow
 * @property {string} field
 * @property {(value: unknown, container: Record<string, unknown>) => FieldReading<unknown>} read
 * @property {boolean} [inPlace]
 * @property {boolean} [retired]
 */
/** @typedef {{ at: string, rows: readonly FieldRow[] }} FieldTable */

/** @param {unknown} value @returns {value is Record<string, unknown>} */
const isRecord = (value) => Boolean(value) && typeof value === 'object' && !Array.isArray(value);

/**
 * @param {string} field
 * @param {FieldReader<unknown>} read
 * @returns {FieldRow}
 */
const row = (field, read) => ({ field, read: (value) => read(value) });

/** @type {FieldRow['read']} */
const retiredField = (value) => ({ value: undefined, lost: !isUnsetValue(value) });

// A Session's Depth series as Load reads them; unset fields take the legacy
// single-series settings of `adv`.
const DEPTH_TRACK_FALLBACK_COLORS = [
  '#4A90E2', '#E45756', '#2CA02C', '#F28E2B', '#9467BD', '#8C564B', '#17BECF', '#7F7F7F'
];
/**
 * @typedef {{ label: string, color: string, height: number | null, large_tick_interval: number | null,
 *   small_tick_interval: number | null, tick_font_size: number | null }} SessionDepthTrack
 */
/**
 * @param {unknown} entry A saved Depth series, unvalidated.
 * @param {number} index
 * @param {Record<string, unknown>} [legacyAdv]
 * @returns {FieldReading<SessionDepthTrack>}
 */
export const readDepthTrack = (entry, index, legacyAdv = {}) => {
  const source = isRecord(entry) ? entry : {};
  let lost = !isRecord(entry) && !isUnsetValue(entry);
  /** @param {unknown} value */
  const positive = (value) => {
    const reading = positiveNumberOrNull(value);
    lost ||= reading.lost;
    return reading.value;
  };
  const fallbackColor = index === 0
    ? String(legacyAdv.depth_color || DEPTH_TRACK_FALLBACK_COLORS[0])
    : DEPTH_TRACK_FALLBACK_COLORS[index % DEPTH_TRACK_FALLBACK_COLORS.length];
  const value = {
    label: String(source.label ?? (index === 0 ? 'Depth' : `Depth ${index + 1}`)),
    color: resolveColorToHex(String(source.color || fallbackColor)),
    height: positive(Object.hasOwn(source, 'height') ? source.height : legacyAdv.depth_height),
    large_tick_interval: positive(source.large_tick_interval ?? (index === 0 ? legacyAdv.depth_large_tick_interval : null)),
    small_tick_interval: positive(source.small_tick_interval ?? (index === 0 ? legacyAdv.depth_small_tick_interval : null)),
    tick_font_size: positive(source.tick_font_size ?? (index === 0 ? legacyAdv.depth_tick_font_size : null))
  };
  return { value, lost };
};
/**
 * At least one series.
 * @param {unknown} tracks
 * @param {Record<string, unknown>} [legacyAdv]
 * @returns {FieldReading<SessionDepthTrack[]>}
 */
export const readDepthTracks = (tracks, legacyAdv = {}) => {
  const readings = (Array.isArray(tracks) ? tracks : []).map((entry, index) => readDepthTrack(entry, index, legacyAdv));
  if (!readings.length) readings.push(readDepthTrack(null, 0, legacyAdv));
  return {
    value: readings.map(({ value }) => value),
    lost: readings.some(({ lost }) => lost) || (!Array.isArray(tracks) && !isUnsetValue(tracks))
  };
};

/** @typedef {{ selector: string, row: number }} MultiRecordPosition */
const readGridRow = numberOr(isPositiveInteger, 1);
/**
 * Once per selector, ordered by row (the order within a row kept).
 * @type {FieldReader<MultiRecordPosition[]>}
 */
const readMultiRecordPositions = (value) => {
  let dropped = 0;
  let lost = !Array.isArray(value) && !isUnsetValue(value);
  /** @type {Set<string>} */
  const seen = new Set();
  /** @type {MultiRecordPosition[]} */
  const positions = [];
  (Array.isArray(value) ? value : []).forEach((entry) => {
    const selector = isRecord(entry) ? String(entry.selector ?? '').trim() : '';
    if (!isRecord(entry) || !selector || seen.has(selector)) {
      dropped += 1;
      return;
    }
    const gridRow = readGridRow(entry.row);
    lost ||= gridRow.lost;
    seen.add(selector);
    positions.push({ selector, row: gridRow.value });
  });
  const ordered = positions.map((entry, index) => ({ entry, index }))
    .sort((left, right) => (left.entry.row - right.entry.row) || (left.index - right.index))
    .map(({ entry }) => entry);
  return { value: ordered, lost: lost || dropped > 0, ...(dropped && !lost ? { dropped } : {}) };
};

/**
 * @param {unknown} value
 * @param {string} [fallback]
 * @returns {FieldReading<string>}
 */
const normalizeHexColor = (value, fallback = '#4e79a7') => {
  const resolved = String(resolveColorToHex(String(value || fallback).trim()) || fallback).trim();
  const shortMatch = resolved.match(/^#([0-9a-fA-F]{3})$/);
  if (shortMatch) {
    return { value: `#${shortMatch[1].split('').map((char) => char + char).join('').toLowerCase()}`, lost: false };
  }
  return /^#[0-9a-fA-F]{6}$/.test(resolved) ? { value: resolved.toLowerCase(), lost: false } : { value: fallback, lost: true };
};
const readGencode = numberOr(isPositiveInteger, 1);
/**
 * @typedef {{ sourceKey: string, fileName: string, sourceIndex: number, label: string, color: string,
 *   losat_gencode: number }} ConservationSeries
 */
/** @type {FieldReader<ConservationSeries[]>} */
export const readConservationSeries = (value) => {
  const entries = Array.isArray(value) ? value : [];
  let lost = !Array.isArray(value) && !isUnsetValue(value);
  const kept = entries.filter(isRecord);
  const series = kept.map((entry, index) => {
    const color = normalizeHexColor(entry.color);
    const gencode = readGencode(entry.losat_gencode);
    lost ||= color.lost || gencode.lost;
    return {
      sourceKey: String(entry.sourceKey || ''),
      fileName: String(entry.fileName || ''),
      sourceIndex: Number.isInteger(Number(entry.sourceIndex)) ? Number(entry.sourceIndex) : index,
      label: String(entry.label ?? entry.name ?? ''),
      color: color.value,
      losat_gencode: gencode.value
    };
  });
  const dropped = entries.length - kept.length;
  return { value: series, lost: lost || dropped > 0, ...(dropped && !lost ? { dropped } : {}) };
};

/**
 * @typedef {{ uid: string, row: number, canonicalCardinality?: 'exactly_one', canonicalRow?: number,
 *   canonicalColumn?: number }} LinearRecordRowEntry
 */
/** The rows that name a record and a row. @type {FieldReader<LinearRecordRowEntry[]>} */
export const readLinearRecordRows = (value) => {
  const entries = (Array.isArray(value) ? value : []).map((entry) => {
    const source = isRecord(entry) ? entry : {};
    return {
      uid: String(source.uid || ''), row: Number(source.row),
      ...(source.canonicalCardinality === 'exactly_one' ? { canonicalCardinality: /** @type {const} */ ('exactly_one') } : {}),
      ...(Number(source.canonicalRow) === Number(source.row)
        && typeof source.canonicalColumn === 'number' && Number.isInteger(source.canonicalColumn) && source.canonicalColumn > 0
        ? { canonicalRow: Number(source.canonicalRow), canonicalColumn: source.canonicalColumn }
        : {})
    };
  });
  const kept = entries.filter((entry) => entry.uid && Number.isInteger(entry.row) && entry.row > 0);
  const dropped = entries.length - kept.length;
  return { value: kept, lost: dropped > 0, ...(dropped ? { dropped } : {}) };
};

// The Web's lenient reading of a stored label rendering; Python refuses an
// unknown value (twin row `function:normalizeLabelRendering`).
const readLabelRenderingChoice = textChoice(['auto', 'embedded_only', 'external_only'], 'auto');
/** @type {FieldReader<string>} */
const normalizeLabelRendering = (value) => readLabelRenderingChoice(value);

/** @param {number} numeric */
const isRadiusRatio = (numeric) => numeric > 0 && numeric <= 1;
/** @param {number} numeric */
const isNonNegative = (numeric) => numeric >= 0;
const readNonNegativeNumberOrNull = nonNegativeNumberOrNull;

/** `adv`, in the order the rows depend on each other. @type {readonly FieldRow[]} */
const ADV_ROWS = Object.freeze([
  row('label_rendering', normalizeLabelRendering),
  row('circular_label_placement', textChoice(['radial', 'horizontal'], 'horizontal')),
  row('circular_label_spacing', positiveNumberOrNull),
  row('linear_label_spacing', positiveNumberOrNull),
  row('track_axis_gap', readNonNegativeNumberOrNull),
  row('multi_record_min_radius_ratio', numberOr(isRadiusRatio, 0.55)),
  row('multi_record_column_gap_ratio', numberOr(isNonNegative, 0.10, { coerceUnset: true })),
  row('multi_record_row_gap_ratio', numberOr(isNonNegative, 0.05, { coerceUnset: true })),
  row('multi_record_positions', readMultiRecordPositions),
  row('plot_title_font_size', positiveNumberOrNull),
  row('keep_full_definition_with_plot_title', onlyTrue),
  row('depth_color', (value) => ({ value: resolveColorToHex(String(value || '#4A90E2')), lost: false })),
  row('depth_normalize', onlyTrue),
  row('depth_show_axis', trueUnlessFalse),
  row('depth_show_ticks', trueUnlessFalse),
  row('depth_share_axis', onlyTrue),
  row('depth_height', positiveNumberOrNull),
  row('depth_width_circular', positiveNumberOrNull),
  row('depth_window_size', positiveNumberOrNull),
  row('depth_step_size', positiveNumberOrNull),
  row('depth_large_tick_interval', positiveNumberOrNull),
  row('depth_small_tick_interval', positiveNumberOrNull),
  row('depth_tick_font_size', positiveNumberOrNull),
  { field: 'depth_tracks', read: (value, adv) => readDepthTracks(value, adv), inPlace: true },
  row('gc_content_mode', textChoice(['percent', 'deviation'], 'deviation')),
  row('gc_content_show_axis', trueUnlessFalse),
  row('gc_content_show_ticks', trueUnlessFalse),
  row('gc_content_tick_interval', positiveNumberOrNull),
  row('gc_content_small_tick_interval', positiveNumberOrNull),
  row('gc_content_tick_font_size', positiveNumberOrNull),
  row('center_reserved_radius', readNonNegativeNumberOrNull),
  row('depth_min', readNonNegativeNumberOrNull),
  {
    field: 'depth_max',
    // Below the minimum: unset.
    read: (value, adv) => {
      const reading = readNonNegativeNumberOrNull(value);
      return typeof adv.depth_min === 'number' && reading.value !== null && adv.depth_min > reading.value
        ? { value: null, lost: true } : reading;
    }
  },
  row('gc_content_min_percent', numberOr(Number.isFinite, 0)),
  {
    field: 'gc_content_max_percent',
    // Below the minimum: the minimum.
    read: (value, adv) => {
      const reading = numberOr(Number.isFinite, 100)(value);
      return typeof adv.gc_content_min_percent === 'number' && adv.gc_content_min_percent > reading.value
        ? { value: adv.gc_content_min_percent, lost: !isUnsetValue(value) } : reading;
    }
  },
  row('linear_show_replicon', onlyTrue),
  row('pairwise_match_style', textChoice(['ribbon', 'curve'], 'ribbon'))
]);

/** `losat.blastp`: the protein search settings. @type {readonly FieldRow[]} */
const LOSAT_SEARCH_ROWS = Object.freeze([
  row('mode', textChoice(['pairwise', 'orthogroup', 'collinear'], 'orthogroup')),
  row('maxHits', numberOr(isPositiveInteger, 5)),
  row('candidateLimit', numberOr(isPositiveInteger, null)),
  row('orthogroupMembershipMode', readOrthogroupMembershipMode),
  row('orthogroupMemberMaxHits', numberOr(isPositiveInteger, null)),
  row('collinearMinAnchors', numberOr(isPositiveInteger, 1)),
  row('collinearMaxUnitGap', numberOr(isNonNegativeInteger, 0, { coerceUnset: true })),
  row('collinearMaxDiagonalDrift', numberOr(isNonNegativeInteger, 0, { coerceUnset: true })),
  row('collinearMaxConflictsInMergeGap', numberOr(isNonNegativeInteger, 1, { coerceUnset: true })),
  row('collinearMaxParalogLinksPerOrthogroup', numberOr(isPositiveInteger, 2)),
  row('collinearColorMode', textChoice(['average_identity', 'orientation', 'orientation_identity'], 'orientation', {
    aliases: { identity: 'average_identity' }, foldHyphens: true
  })),
  row('collinearUnitMode', textChoice(['auto', 'cds', 'locus'], 'auto')),
  row('collinearAnchorMode', readCollinearAnchorMode),
  row('collinearMergeOrientation', textChoice(['strand', 'order', 'either'], 'either')),
  row('collinearSearchScope', readCollinearSearchScope),
  { field: 'collinearBlockMergeGap', read: retiredField, retired: true },
  { field: 'collinearSingletonMergeGap', read: retiredField, retired: true },
  { field: 'orthogroupHitPolicy', read: retiredField, retired: true },
  // Read as `orthogroupMemberMaxHits` when that is unset.
  {
    field: 'orthogroupMaxHits',
    read: (value, blastp) => ({ value: undefined, lost: !isUnsetValue(value) && blastp.orthogroupMemberMaxHits !== undefined }),
    retired: true
  }
]);

/**
 * How LOSAT runs (an app setting; a Session 44 or older draft keeps it in
 * `config.losat`). Each field reads as the trimmed lowercase text.
 * @param {readonly string[]} words
 * @param {string} fallback
 * @param {(numeric: number) => boolean} [counts] accepts a count instead
 * @returns {FieldReader<string>}
 */
const executionWord = (words, fallback, counts = () => false) => (value) => {
  const raw = String(value ?? fallback).trim().toLowerCase();
  return words.includes(raw) || counts(Number(raw)) ? { value: raw, lost: false } : { value: fallback, lost: raw !== '' };
};
const readThreadBudget = executionWord(['safe', 'auto', 'available'], 'safe', isPositiveInteger);
/** @type {readonly FieldRow[]} */
export const LOSAT_EXECUTION_ROWS = Object.freeze([
  row('parallelWorkers', (value) => {
    const reading = executionWord([], '', isPositiveInteger)(value);
    return { value: reading.value || undefined, lost: reading.lost };
  }),
  row('executionMode', textChoice(['auto', 'serial', 'threaded'], 'auto')),
  row('threadsPerJob', executionWord(['auto'], 'auto', isPositiveInteger)),
  // `auto` is an older name of `safe` (`resolveLosatThreadPlan` reads both alike).
  row('totalThreadBudget', (value) => {
    const reading = readThreadBudget(value);
    return { ...reading, value: reading.value === 'auto' ? 'safe' : reading.value };
  })
]);

/** `circularConservation`. @type {readonly FieldRow[]} */
const CONSERVATION_ROWS = Object.freeze([
  row('enabled', onlyTrue),
  row('source', textChoice(['upload', 'losat'], 'losat')),
  row('losat_program', textChoice(['tblastx', 'blastn'], 'blastn')),
  row('subject_gencode', numberOr(isPositiveInteger, 1)),
  row('reference', textChoice(['auto', 'query', 'subject'], 'auto')),
  { field: 'series', read: (value) => readConservationSeries(value), inPlace: true },
  row('ring_width', positiveNumberOrNull),
  row('ring_gap', positiveNumberOrNull)
]);

/** Omitted: the fresh default in every Session version. @type {FieldReader<boolean>} */
export const readLinearRecordLayoutEnabled = (value) => (typeof value === 'boolean' ? { value, lost: false }
  : { value: WEB_UX_PROFILE.linear.arrangeInRowsByDefault, lost: !isUnsetValue(value) });
export const readLinearRecordGap = numberOr(isNonNegative, 24, { coerceUnset: true });
export const readImportedComparisonAction = exactChoice(Object.values(IMPORTED_COMPARISON_ACTIONS), null);
const readLosatProgramText = exactChoice(['blastn', 'tblastx', 'blastp'], 'blastn');
/** Read as text. @type {FieldReader<string>} */
export const readLosatProgram = (value) => readLosatProgramText(value ? String(value) : value);

/** @type {readonly FieldTable[]} */
const CONFIG_FIELD_TABLE = Object.freeze([
  { at: 'adv', rows: ADV_ROWS },
  { at: 'losat.blastp', rows: LOSAT_SEARCH_ROWS },
  { at: 'losat', rows: LOSAT_EXECUTION_ROWS },
  { at: 'circularConservation', rows: CONSERVATION_ROWS },
  {
    at: 'linearRecordLayout',
    rows: [
      row('enabled', readLinearRecordLayoutEnabled),
      row('recordGap', readLinearRecordGap),
      row('rows', readLinearRecordRows)
    ]
  },
  { at: 'importedComparisonResolution', rows: [row('action', readImportedComparisonAction)] },
  { at: '', rows: [row('losatProgram', readLosatProgram)] }
]);

export const readSessionUiMode = exactChoice(['circular', 'linear'], 'circular');
export const readFeaturePanelTab = exactChoice(['colors', 'labels'], 'colors');
/** A Session's app `ui`. @type {readonly FieldTable[]} */
const SESSION_UI_TABLE = Object.freeze([
  {
    at: '',
    rows: [
      row('mode', readSessionUiMode),
      row('featurePanelTab', readFeaturePanelTab),
      row('autoLabelReflow', truthy),
      row('paletteInstantPreviewEnabled', truthy),
      row('generatedMultiRecordCanvas', truthy),
      row('generatedCircularPlotTitlePosition', readCircularPlotTitlePosition)
    ]
  }
]);

/**
 * Reads the rows into `container` (an `adv`, `losat.blastp`, or
 * `circularConservation` the apply merged a stored configuration into).
 * @param {Record<string, unknown>} container
 * @param {'adv' | 'losat.blastp' | 'circularConservation'} at
 */
export const applyConfigFieldRows = (container, at) => {
  const table = CONFIG_FIELD_TABLE.find((entry) => entry.at === at);
  (table ? table.rows : []).forEach(({ field, read, inPlace, retired }) => {
    if (retired) {
      delete container[field];
      return;
    }
    const { value } = read(container[field], container);
    const current = container[field];
    if (inPlace && Array.isArray(current) && Array.isArray(value)) current.splice(0, current.length, ...value);
    else container[field] = value;
  });
};

/** @param {unknown} value */
const excerpt = (value) => (typeof value === 'string' ? value : JSON.stringify(value) ?? '');

/**
 * @param {unknown} value
 * @param {string} at Dotted; '' is the value itself.
 * @returns {Record<string, unknown> | null}
 */
const containerAt = (value, at) => {
  let container = value;
  for (const key of at.split('.').filter(Boolean)) container = isRecord(container) ? container[key] : undefined;
  return isRecord(container) ? container : null;
};

/**
 * What the rows of `tables` cannot restore of the values `stored` holds.
 * @param {readonly FieldTable[]} tables
 * @param {unknown} value
 * @param {string} path
 * @param {unknown} stored
 * @returns {LoadRepair[]}
 */
const fieldRepairs = (tables, value, path, stored) => tables.flatMap(({ at, rows }) => {
  const source = containerAt(value, at);
  const storedContainer = containerAt(stored, at);
  if (!isRecord(source) || !isRecord(storedContainer)) return [];
  /** @type {Record<string, unknown>} */
  const read = { ...source };
  return rows.flatMap(({ field, read: reader }) => {
    if (!Object.hasOwn(source, field)) return [];
    const reading = reader(source[field], read);
    read[field] = reading.value;
    if (!reading.lost || !Object.hasOwn(storedContainer, field)) return [];
    const fieldPath = [path, at, field].filter(Boolean).join('.');
    if (reading.dropped) {
      return [loadRepair({ code: 'setting-entries-dropped', path: fieldPath, effect: 'dropped', count: reading.dropped })];
    }
    return [loadRepair({
      code: 'setting-replaced', path: fieldPath, effect: 'defaulted', stored: excerpt(source[field]),
      ...(isUnsetValue(reading.value) ? {} : { used: excerpt(reading.value) })
    })];
  });
});

/**
 * The stored settings of a configuration that its readers cannot restore, as
 * Load lists them before `applyConfigData` reads the configuration with the
 * same rows. Only the values `stored` holds are listed (the file's draft; a
 * value the stored request supplies is read with the request).
 * @param {unknown} config
 * @param {string} path Where the configuration is in the file.
 * @param {unknown} [stored]
 */
export const configFieldRepairs = (config, path, stored = config) => fieldRepairs(CONFIG_FIELD_TABLE, config, path, stored);

/**
 * The same for a Session's app `ui`.
 * @param {unknown} ui
 * @param {string} path
 */
export const sessionUiFieldRepairs = (ui, path) => fieldRepairs(SESSION_UI_TABLE, ui, path, ui);

/**
 * The same for the LOSAT execution settings (`ui.losatExecution` of a Session
 * 46, `config.losat` of an older draft).
 * @param {unknown} execution
 * @param {string} path
 */
export const losatExecutionRepairs = (execution, path) => fieldRepairs([{ at: '', rows: LOSAT_EXECUTION_ROWS }], execution, path, execution);
