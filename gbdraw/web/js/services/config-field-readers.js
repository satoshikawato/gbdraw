// @ts-check
import { MODE_SCOPED_SETTINGS } from '../mode-scoped-settings.generated.js';
import { resolveColorToHex } from '../utils/color-utils.js';
import {
  changedLeaves, exactChoice, isNonNegativeInteger, isPositiveInteger, isUnsetValue, keepsMeaning, nonNegativeNumberOrNull, numberOr, onlyTrue,
  positiveNumberOrNull, textChoice, trueUnlessFalse, truthy
} from '../utils/field-reading.js';
import { loadRepair } from '../utils/load-repairs.js';
import { normalizeAnnotationSets } from './annotation-state.js';
import { PAIRWISE_MATCH_STYLES } from './current-option-values.js';
import { readDefinitionLineStyles } from './definition-line-style-state.js';
import { readDepthTracks } from './depth-track-state.js';
import { IMPORTED_COMPARISON_ACTIONS } from './imported-comparison-intent.js';
import { readCircularPlotTitlePosition } from './layout-preferences.js';
import { isPlainObject } from './losat-cache.js';
import { readCollinearAnchorMode, readCollinearSearchScope, readOrthogroupMembershipMode } from './losat-normalization.js';
import { normalizeOrthogroupDormantOverrides } from './orthogroup-feature-metadata.js';
import {
  createDefaultAdv, createDefaultCircularConservation, createDefaultForm, createDefaultLosat, createDefaultLosatpHitLimits,
  readLinearTypographyLinked
} from './session-active-config-contract.js';
import { WEB_UX_PROFILE } from '../web-ux-profile.js';
import { FEATURE_RENDERING_VALUES, normalizedPersistedArrowGeometryState } from '../utils/feature-rendering.js';
import { LINEAR_LABEL_VISIBILITY_MODES } from './linear-label-visibility.js';

/** @import { FieldLeaf, FieldReader, FieldReading } from '../utils/field-reading.js' */
/** @import { LoadRepair } from '../utils/load-repairs.js' */

// The readers of the stored settings of a configuration (a drawing's draft, a
// legacy configuration file) and of a Session's app `ui`, a drawing's `ui`
// and its feature edits. `applyConfigData` and Load's apply read with these
// readers, History restores included; before Load applies a file it reads the
// file with the same readers and lists what they cannot restore (BA-D25,
// B36, B37, B43), so the apply makes no repair of its own (G2).

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
/**
 * A stored value a reading did not restore: `at` is its place in the read
 * value (`adv.depth_tracks[0].height`); `dropped` counts the entries of a list
 * it dropped; `kept`: the stored value is not used, and the target keeps the
 * value it holds.
 * @typedef {{ at: string, stored: unknown, used: unknown, dropped?: number, kept?: boolean }} SettingLoss
 */

/**
 * @param {string} field
 * @param {FieldReader<unknown>} read
 * @returns {FieldRow}
 */
const row = (field, read) => ({ field, read: (value) => read(value) });

/** @type {FieldRow['read']} */
const retiredField = (value) => ({ value: undefined, lost: !isUnsetValue(value) });

/** An object, else empty. @type {FieldReader<Record<string, unknown>>} */
export const readPlainObject = (value) => (isPlainObject(value) ? { value, lost: false } : { value: {}, lost: !isUnsetValue(value) });
/** A list, else empty. @type {FieldReader<unknown[]>} */
export const readList = (value) => (Array.isArray(value) ? { value, lost: false } : { value: [], lost: !isUnsetValue(value) });
/**
 * A list `read` reads; any other value is not used, and the drawing keeps the
 * list it holds (`applyConfigData` installs a configuration's `rules` and
 * `whitelist` so).
 * @template T
 * @param {FieldReader<T[]>} read
 * @returns {FieldReader<T[]>}
 */
const listOrKept = (read) => (value) => (Array.isArray(value) ? read(value)
  : { value: [], lost: !isUnsetValue(value), kept: true });

/** @typedef {{ selector: string, row: number }} MultiRecordPosition */
const readGridRow = numberOr(isPositiveInteger, 1);
/**
 * Once per selector, ordered by row (the order within a row kept).
 * @type {FieldReader<MultiRecordPosition[]>}
 */
const readMultiRecordPositions = (value) => {
  let dropped = 0;
  /** @type {FieldLeaf[]} */
  const leaves = [];
  /** @type {Set<string>} */
  const seen = new Set();
  /** @type {MultiRecordPosition[]} */
  const positions = [];
  (Array.isArray(value) ? value : []).forEach((entry, index) => {
    const selector = isPlainObject(entry) ? String(entry.selector ?? '').trim() : '';
    if (!isPlainObject(entry) || !selector || seen.has(selector)) {
      dropped += 1;
      return;
    }
    const gridRow = readGridRow(entry.row);
    if (gridRow.lost) leaves.push({ at: `[${index}].row`, stored: entry.row, used: gridRow.value });
    seen.add(selector);
    positions.push({ selector, row: gridRow.value });
  });
  const ordered = positions.map((entry, index) => ({ entry, index }))
    .sort((left, right) => (left.entry.row - right.entry.row) || (left.index - right.index))
    .map(({ entry }) => entry);
  return listReading(ordered, value, dropped, leaves);
};

/**
 * The reading of a list: lost when the stored value is no list, when entries
 * were dropped, or when values inside it were replaced.
 * @template T
 * @param {T} read
 * @param {unknown} stored
 * @param {number} dropped
 * @param {readonly FieldLeaf[]} leaves
 * @returns {FieldReading<T>}
 */
const listReading = (read, stored, dropped, leaves) => {
  const notList = !Array.isArray(stored) && !isUnsetValue(stored);
  return {
    value: read,
    lost: notList || dropped > 0 || leaves.length > 0,
    ...(dropped && !notList ? { dropped } : {}),
    ...(leaves.length && !notList ? { leaves } : {})
  };
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
  /** @type {FieldLeaf[]} */
  const leaves = [];
  /** @type {ConservationSeries[]} */
  const series = [];
  entries.forEach((entry, index) => {
    if (!isPlainObject(entry)) return;
    const color = normalizeHexColor(entry.color);
    const gencode = readGencode(entry.losat_gencode);
    if (color.lost) leaves.push({ at: `[${index}].color`, stored: entry.color, used: color.value });
    if (gencode.lost) leaves.push({ at: `[${index}].losat_gencode`, stored: entry.losat_gencode, used: gencode.value });
    // A series without its file index takes its place in the list.
    const sourceIndex = numberOr(Number.isInteger, series.length, { coerceUnset: true })(entry.sourceIndex);
    if (sourceIndex.lost) leaves.push({ at: `[${index}].sourceIndex`, stored: entry.sourceIndex, used: sourceIndex.value });
    series.push({
      sourceKey: String(entry.sourceKey || ''),
      fileName: String(entry.fileName || ''),
      sourceIndex: sourceIndex.value,
      label: String(entry.label ?? entry.name ?? ''),
      color: color.value,
      losat_gencode: gencode.value
    });
  });
  return listReading(series, value, entries.length - series.length, leaves);
};

/**
 * @typedef {{ uid: string, row: number, canonicalCardinality?: 'exactly_one', canonicalRow?: number,
 *   canonicalColumn?: number }} LinearRecordRowEntry
 */
/**
 * The rows that name a record and a row. The canonical place of a row is kept
 * when it names the row's own row and a column; other canonical fields are
 * dropped.
 * @type {FieldReader<LinearRecordRowEntry[]>}
 */
export const readLinearRecordRows = (value) => {
  /** @type {FieldLeaf[]} */
  const leaves = [];
  let dropped = 0;
  /** @type {LinearRecordRowEntry[]} */
  const kept = [];
  (Array.isArray(value) ? value : []).forEach((entry, index) => {
    const source = isPlainObject(entry) ? entry : {};
    const uid = String(source.uid || '');
    const rowNumber = Number(source.row);
    if (!uid || !Number.isInteger(rowNumber) || rowNumber <= 0) {
      dropped += 1;
      return;
    }
    const canonical = Number(source.canonicalRow) === rowNumber
      && typeof source.canonicalColumn === 'number' && Number.isInteger(source.canonicalColumn) && source.canonicalColumn > 0;
    ['canonicalRow', 'canonicalColumn'].forEach((field) => {
      if (!canonical && !isUnsetValue(source[field])) leaves.push({ at: `[${index}].${field}`, stored: source[field], used: undefined });
    });
    if (!isUnsetValue(source.canonicalCardinality) && source.canonicalCardinality !== 'exactly_one') {
      leaves.push({ at: `[${index}].canonicalCardinality`, stored: source.canonicalCardinality, used: undefined });
    }
    kept.push({
      uid, row: rowNumber,
      ...(source.canonicalCardinality === 'exactly_one' ? { canonicalCardinality: /** @type {const} */ ('exactly_one') } : {}),
      ...(canonical ? { canonicalRow: rowNumber, canonicalColumn: /** @type {number} */ (source.canonicalColumn) } : {})
    });
  });
  return listReading(kept, value, dropped, leaves);
};

/** @typedef {{ feat: string, order: string }} QualifierPriorityRule */
/**
 * One rule per feature type, the later rule of a type winning; a rule without
 * a type or an order is dropped.
 * @type {FieldReader<QualifierPriorityRule[]>}
 */
export const readQualifierPriorityRules = (value) => {
  /** @type {QualifierPriorityRule[]} */
  const rules = [];
  let dropped = 0;
  (Array.isArray(value) ? value : []).forEach((rule) => {
    const feat = isPlainObject(rule) ? String(rule.feat ?? '').trim() : '';
    const order = isPlainObject(rule) ? String(rule.order ?? '').trim() : '';
    if (!feat || !order) {
      dropped += 1;
      return;
    }
    const existing = rules.find((entry) => entry.feat === feat);
    if (existing) {
      existing.order = order;
      dropped += 1;
    } else {
      rules.push({ feat, order });
    }
  });
  return listReading(rules, value, dropped, []);
};

/**
 * A list of rows read entry by entry: an entry that is no object is dropped
 * and counted; `readEntry` gives the row and the values of the entry it
 * replaced. A field the row does not name is not read.
 * @template T
 * @param {(entry: Record<string, unknown>) => { row: T, leaves: FieldLeaf[] }} readEntry
 * @returns {FieldReader<T[]>}
 */
const rowList = (readEntry) => (value) => {
  /** @type {T[]} */
  const rows = [];
  /** @type {FieldLeaf[]} */
  const leaves = [];
  let dropped = 0;
  (Array.isArray(value) ? value : []).forEach((entry, index) => {
    if (!isPlainObject(entry)) {
      dropped += 1;
      return;
    }
    const reading = readEntry(entry);
    rows.push(reading.row);
    reading.leaves.forEach((leaf) => leaves.push({ ...leaf, at: `[${index}]${leaf.at}` }));
  });
  return listReading(rows, value, dropped, leaves);
};
/** A rule's text field: an unset or falsy value is ''. @param {unknown} value */
const ruleText = (value) => String(value || '');
/**
 * @param {Record<string, unknown>} entry
 * @param {Record<string, unknown>} row
 * @param {readonly string[]} fields
 * @returns {FieldLeaf[]}
 */
const textFieldLeaves = (entry, row, fields) => fields.flatMap((field) => changedLeaves(entry[field], row[field], `.${field}`));

/** @typedef {{ feat: string, qual: string, val: string, color: string, cap: string, fromFile?: true }} StoredColorRule */
/**
 * The specific color rules of a configuration (`rules`). A color name reads as
 * its hex value; an unset color is black; `fromFile` is kept only when set.
 * A value that is no list is not used (`listOrKept`).
 * @type {FieldReader<StoredColorRule[]>}
 */
export const readColorRules = listOrKept(rowList((entry) => {
  const color = resolveColorToHex(String(entry.color || '#000000'));
  /** @type {StoredColorRule} */
  const row = {
    feat: ruleText(entry.feat), qual: ruleText(entry.qual), val: ruleText(entry.val), color, cap: ruleText(entry.cap),
    // The shape normalizeSpecificRule gives: `fromFile` only when true.
    ...(entry.fromFile ? { fromFile: true } : {})
  };
  const sameColor = typeof entry.color === 'string' && resolveColorToHex(entry.color) === color;
  return {
    row,
    leaves: [
      ...textFieldLeaves(entry, row, ['feat', 'qual', 'val', 'cap']),
      ...(!isUnsetValue(entry.color) && !sameColor && !keepsMeaning(entry.color, color)
        ? [{ at: '.color', stored: entry.color, used: color }] : []),
      ...(!isUnsetValue(entry.fromFile) && typeof entry.fromFile !== 'boolean' && !keepsMeaning(entry.fromFile, Boolean(entry.fromFile))
        ? [{ at: '.fromFile', stored: entry.fromFile, used: Boolean(entry.fromFile) }] : [])
    ]
  };
}));

/** @typedef {{ feat: string, qual: string, key: string }} StoredWhitelistRule */
/**
 * The label whitelist of a configuration (`whitelist`); a value that is no
 * list is not used (`listOrKept`).
 * @type {FieldReader<StoredWhitelistRule[]>}
 */
export const readWhitelistRules = listOrKept(rowList((entry) => {
  const row = { feat: ruleText(entry.feat), qual: ruleText(entry.qual), key: ruleText(entry.key) };
  return { row, leaves: textFieldLeaves(entry, row, ['feat', 'qual', 'key']) };
}));

/**
 * The annotation sets of a drawing (`normalizeAnnotationSets`), with each
 * stored value it replaces.
 * @type {FieldReader<ReturnType<typeof normalizeAnnotationSets>>}
 */
export const readAnnotationSets = (value) => {
  const sets = normalizeAnnotationSets(/** @type {Parameters<typeof normalizeAnnotationSets>[0]} */ (Array.isArray(value) ? value : []));
  return listReading(sets, value, 0, Array.isArray(value) ? changedLeaves(value, sets) : []);
};

/**
 * A map of text by key (similarity group names and descriptions, bulk label
 * text): a key without text, or a value that is no text (a number is its
 * text), is dropped.
 * @type {FieldReader<Record<string, string>>}
 */
export const readStringMap = (value) => {
  /** @type {Record<string, string>} */
  const map = {};
  let dropped = 0;
  Object.entries(isPlainObject(value) ? value : {}).forEach(([key, text]) => {
    const name = String(key || '').trim();
    // An entry without a key, or whose value is no text (an object, a
    // boolean; a number is its text), is dropped and counted.
    if (!name || (text !== null && !keepsMeaning(text, String(text)))) {
      dropped += 1;
      return;
    }
    map[name] = String(text ?? '');
  });
  const notObject = !isPlainObject(value) && !isUnsetValue(value);
  return {
    value: map,
    lost: notObject || dropped > 0,
    ...(dropped && !notObject ? { dropped } : {})
  };
};
/** The dormant group names (D-21): an entry with neither text is dropped. @type {FieldReader<Record<string, unknown>>} */
const readDormantOverrides = (value) => {
  const read = normalizeOrthogroupDormantOverrides(value);
  const dropped = isPlainObject(value) ? Object.keys(value).length - Object.keys(read).length : 0;
  const notObject = !isPlainObject(value) && !isUnsetValue(value);
  return { value: read, lost: notObject || dropped > 0, ...(dropped && !notObject ? { dropped } : {}) };
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

/**
 * The feature shapes after the merge (`adv.feature_shapes`): a shape is read
 * without case and outer spaces, a feature type without outer spaces. An
 * entry whose type is empty or whose shape names none is dropped (BA-OD100):
 * the merge copies the shape of a feature type the defaults lack as stored,
 * and a default type's shape of another kind is not used (a text one that
 * names no shape is refused before, `requireCurrentConfigFields`).
 * @type {FieldReader<Record<string, string>>}
 */
const readFeatureShapes = (value) => {
  /** @type {Record<string, string>} */
  const shapes = {};
  /** @type {FieldLeaf[]} */
  const leaves = [];
  Object.entries(isPlainObject(value) ? value : {}).forEach(([key, shape]) => {
    const featureType = key.trim();
    const name = typeof shape === 'string' ? shape.trim().toLowerCase() : '';
    if (featureType && FEATURE_RENDERING_VALUES.includes(name)) shapes[featureType] = name;
    // An entry of no feature type, or an unset shape, drew nothing: no row.
    else if (featureType && !isUnsetValue(shape)) leaves.push({ at: `.${key}`, stored: shape, used: undefined });
  });
  return { value: shapes, lost: leaves.length > 0, ...(leaves.length ? { leaves } : {}) };
};
// A stored visibility the merge installs names a mode (`requireCurrentConfigFields`).
const readLinearLabelVisibility = textChoice(LINEAR_LABEL_VISIBILITY_MODES, 'auto');

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
  row('feature_overlap_tolerance_bp', numberOr((numeric) => Number.isSafeInteger(numeric) && numeric >= 0, 0)),
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
  row('linear_accession_visibility', readLinearLabelVisibility),
  row('linear_length_visibility', readLinearLabelVisibility),
  row('feature_shapes', readFeatureShapes),
  row('pairwise_match_style', textChoice(PAIRWISE_MATCH_STYLES, 'ribbon')),
  row('linear_definition_line_styles', readDefinitionLineStyles)
]);

/**
 * `losat.blastp`: the protein search settings. `orthogroupMaxHits` had no
 * writer on `main` or in a release (OV-756): it is dropped like the other
 * retired keys.
 * @type {readonly FieldRow[]}
 */
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
  { field: 'orthogroupMaxHits', read: retiredField, retired: true }
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
/** The label filter, read without case. */
export const readFilterMode = textChoice(['None', 'Whitelist', 'Blacklist'], 'None', {
  aliases: { none: 'None', whitelist: 'Whitelist', blacklist: 'Blacklist' }
});
/** Read as text. @type {FieldReader<string>} */
export const readLosatProgram = (value) => readLosatProgramText(value ? String(value) : value);

/** The containers the merge reads first, then their rows (`installStoredSettings`). @type {readonly FieldTable[]} */
const CONTAINER_TABLES = Object.freeze([
  { at: 'form', rows: [] },
  { at: 'adv', rows: ADV_ROWS },
  { at: 'losat.blastp', rows: LOSAT_SEARCH_ROWS },
  { at: 'circularConservation', rows: CONSERVATION_ROWS }
]);
/** The fields read from the configuration as stored. @type {readonly FieldTable[]} */
const VALUE_TABLES = Object.freeze([
  { at: 'losat', rows: LOSAT_EXECUTION_ROWS },
  {
    at: 'linearRecordLayout',
    rows: [
      row('enabled', readLinearRecordLayoutEnabled),
      row('recordGap', readLinearRecordGap),
      row('rows', readLinearRecordRows)
    ]
  },
  { at: 'importedComparisonResolution', rows: [row('action', readImportedComparisonAction)] },
  {
    at: 'webEdits',
    rows: [
      row('orthogroupNameOverrides', readStringMap), row('orthogroupDescriptionOverrides', readStringMap),
      row('orthogroupDormantOverrides', readDormantOverrides)
    ]
  },
  {
    at: '',
    rows: [
      row('losatProgram', readLosatProgram),
      row('qualifierPriorityRules', readQualifierPriorityRules),
      // The older name is read only without the current one.
      {
        field: 'priorityRules',
        read: (value, container) => (Object.hasOwn(container, 'qualifierPriorityRules')
          ? { value, lost: false } : readQualifierPriorityRules(value))
      },
      row('annotationSets', readAnnotationSets),
      row('rules', readColorRules),
      row('whitelist', readWhitelistRules),
      row('filterMode', readFilterMode),
      row('colors', readPlainObject),
      row('featurePlacementOverrides', readPlainObject),
      row('unmanagedConfigOverrides', readPlainObject)
    ]
  }
]);

/** @param {unknown} value */
const storedText = (value) => typeof value === 'string' && value.trim() !== '';
/**
 * The plot title position the shown Result was drawn with: the stored one, or
 * else the older `circularPlotTitlePosition` of the same `ui`, which wins
 * without a loss over a stored value that is no text.
 * @param {unknown} value
 * @param {Record<string, unknown>} ui
 * @returns {FieldReading<string>}
 */
export const readGeneratedCircularPlotTitlePosition = (value, ui) => {
  if (storedText(value)) return readCircularPlotTitlePosition(value);
  const fallback = readCircularPlotTitlePosition(ui.circularPlotTitlePosition);
  return { value: fallback.value, lost: !isUnsetValue(value) && !storedText(ui.circularPlotTitlePosition) };
};
/**
 * A saved palette's colors by feature type, each read as a color; a value
 * that is no object is not read (`null`: the drawing's colors are used).
 * @type {FieldReader<Record<string, string> | null>}
 */
export const readPaletteColors = (value) => (isPlainObject(value)
  ? {
      value: Object.fromEntries(Object.entries(value).map(([key, color]) => [key, resolveColorToHex(String(color || '').trim())])),
      lost: false
    }
  : { value: null, lost: !isUnsetValue(value) });
/** Shown unless `false`. @type {FieldReader<boolean>} */
export const readRichFeaturePopup = trueUnlessFalse;
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
      { field: 'circularPlotTitlePosition', read: (value, ui) => (storedText(ui.generatedCircularPlotTitlePosition)
        ? { value, lost: false } : readCircularPlotTitlePosition(value)) },
      { field: 'generatedCircularPlotTitlePosition', read: readGeneratedCircularPlotTitlePosition },
      row('appliedPaletteColors', readPaletteColors)
    ]
  }
]);

/** @type {FieldReader<number>} */
const readPaddingSide = numberOr(Number.isFinite, 0);
/** A drawing's Canvas padding: a side that is no number is 0. @type {FieldReader<Record<'top' | 'right' | 'bottom' | 'left', number>>} */
export const readCanvasPadding = (value) => {
  const stored = isPlainObject(value) ? value : {};
  /** @type {FieldLeaf[]} */
  const leaves = [];
  const padding = { top: 0, right: 0, bottom: 0, left: 0 };
  /** @type {('top' | 'right' | 'bottom' | 'left')[]} */ (['top', 'right', 'bottom', 'left']).forEach((side) => {
    const reading = readPaddingSide(stored[side]);
    if (reading.lost) leaves.push({ at: `.${side}`, stored: stored[side], used: reading.value });
    padding[side] = reading.value;
  });
  const notObject = !isPlainObject(value) && !isUnsetValue(value);
  return { value: padding, lost: notObject || leaves.length > 0, ...(leaves.length ? { leaves } : {}) };
};
/** A drawing's `ui` and its feature edits (a Session 46 slice). @type {readonly FieldTable[]} */
const SLICE_TABLES = Object.freeze([
  {
    at: 'ui',
    rows: [
      row('canvasPadding', readCanvasPadding),
      row('linearTypographyLinked', readLinearTypographyLinked)
    ]
  },
  {
    at: 'features',
    rows: [
      row('featureColorOverrides', readPlainObject),
      row('featureVisibilityManualRules', readList),
      row('labelOverrideRows', readList),
      row('labelTextBulkOverrides', readStringMap)
    ]
  }
]);

/**
 * @param {unknown} value
 * @param {string} at Dotted; '' is the value itself.
 * @returns {Record<string, unknown> | null}
 */
const containerAt = (value, at) => {
  let container = value;
  for (const key of at.split('.').filter(Boolean)) container = isPlainObject(container) ? container[key] : undefined;
  return isPlainObject(container) ? container : null;
};

/**
 * Reads the rows into `container`; returns what they could not restore.
 * `stored` holds the retired fields (the merge does not copy them).
 * @param {Record<string, unknown>} container
 * @param {readonly FieldRow[]} rows
 * @param {string} at
 * @param {Record<string, unknown> | null} [stored]
 * @returns {SettingLoss[]}
 */
const readRows = (container, rows, at, stored = container) => rows.flatMap(({ field, read, inPlace, retired }) => {
  const fieldAt = [at, field].filter(Boolean).join('.');
  if (retired) {
    const storedValue = stored?.[field];
    delete container[field];
    return retiredField(storedValue, container).lost ? [{ at: fieldAt, stored: storedValue, used: undefined }] : [];
  }
  const before = container[field];
  const reading = read(before, container);
  if (inPlace && Array.isArray(before) && Array.isArray(reading.value)) before.splice(0, before.length, ...reading.value);
  else container[field] = reading.value;
  if (!reading.lost) return [];
  if (reading.kept) return [{ at: fieldAt, stored: before, used: undefined, kept: true }];
  if (reading.dropped || reading.leaves) {
    return [
      ...(reading.dropped ? [{ at: fieldAt, stored: before, used: undefined, dropped: reading.dropped }] : []),
      ...(reading.leaves ?? []).map((leaf) => ({ at: `${fieldAt}${leaf.at}`, stored: leaf.stored, used: leaf.used }))
    ];
  }
  return [{ at: fieldAt, stored: before, used: reading.value }];
});

/** @param {unknown} value A list, an object, or the `typeof` of anything else. */
const kindOf = (value) => (Array.isArray(value) ? 'list' : typeof value);
/** Stored maps whose keys the user names, merged with their own keys. @type {readonly string[]} */
const FREE_KEY_MAPS = Object.freeze(['adv.feature_shapes']);
/** @type {readonly string[]} */
const UNSAFE_KEYS = Object.freeze(['__proto__', 'constructor', 'prototype']);
/**
 * Merges a stored settings object into `target`: a key `target` lacks is
 * ignored, except in a map of `FREE_KEY_MAPS`, where it is copied; an object merges into an object; a list replaces a list; a value
 * of `target`'s kind, or `null` on either side, is copied. A value of another
 * kind (a list for an object included) is not used; it is returned as a loss
 * at `at`, and `target` keeps its value.
 * @param {Record<string, unknown>} target
 * @param {unknown} source
 * @param {string} [at]
 * @returns {SettingLoss[]}
 */
export const mergeStoredSettings = (target, source, at = '') => {
  if (!source || typeof source !== 'object') return [];
  const stored = /** @type {Record<string, unknown>} */ (source);
  return Object.keys(stored).flatMap((key) => {
    if (UNSAFE_KEYS.includes(key)) return [];
    if (!Object.hasOwn(target, key)) {
      // A shape of a feature type the defaults do not name is the user's own.
      if (FREE_KEY_MAPS.includes(at)) target[key] = copied(stored[key]);
      return [];
    }
    const targetValue = target[key];
    const sourceValue = stored[key];
    const keyAt = at ? `${at}.${key}` : key;
    if (isPlainObject(targetValue) && isPlainObject(sourceValue)) return mergeStoredSettings(targetValue, sourceValue, keyAt);
    if (Array.isArray(targetValue) && Array.isArray(sourceValue)) {
      targetValue.splice(0, targetValue.length, ...sourceValue);
      return [];
    }
    if (kindOf(targetValue) === kindOf(sourceValue) || targetValue === null || sourceValue === null) {
      target[key] = sourceValue;
      return [];
    }
    return [{ at: keyAt, stored: sourceValue, used: targetValue, kept: true }];
  });
};

/**
 * The settings containers of a drawing a configuration is installed into: the
 * drawing's own, or fresh defaults when Load reads a file before it applies it.
 * @typedef {{ form: Record<string, unknown>, adv: Record<string, unknown>,
 *   losat: { blastp: Record<string, unknown> }, circularConservation: Record<string, unknown> }} SettingsContainers
 */
/**
 * A drawing's settings as Load resets it (`resetSessionBaseline`).
 * @param {'circular' | 'linear'} mode
 * @returns {SettingsContainers}
 */
export const defaultSettingsContainers = (mode) => ({
  form: createDefaultForm(), adv: createDefaultAdv(mode), losat: createDefaultLosat(),
  circularConservation: createDefaultCircularConservation()
});
/** The `adv` fields copied as stored, whatever their kind, when the stored `adv` holds them: the two Linear font sizes. */
const RAW_ADV_FIELDS = Object.freeze(['scale_font_size', 'ruler_label_font_size']);
/**
 * Fields read as stored, whatever their kind (their merge losses are none):
 * the `RAW_ADV_FIELDS`; the arrow ratios, read as numbers below (a ratio
 * stored as its text reads as its number; one that is no number is refused
 * before, `requireCurrentConfigFields`); and the inference flag, copied with
 * the hit limits below.
 */
const RAW_FIELDS = Object.freeze([
  ...RAW_ADV_FIELDS.map((field) => `adv.${field}`), 'adv.arrow_head_length_ratio', 'adv.arrow_shaft_width_ratio',
  'losat.blastp.collinearInferOrthogroups'
]);
/** @param {unknown} value */
const copied = (value) => (value === undefined ? undefined : JSON.parse(JSON.stringify(value)));

/**
 * Installs a configuration's settings into `targets`: each container is
 * merged, then read by its rows. `applyConfigData` installs a drawing's
 * settings with it, and Load reads a file with it before it applies it, so
 * both lose the same values.
 * @param {SettingsContainers} targets
 * @param {Record<string, unknown>} config
 * @returns {SettingLoss[]}
 */
export const installStoredSettings = (targets, config) => {
  /** A stored container that is no object is not used (the drawing's stays). @param {Record<string, unknown>} target @param {string} at */
  const merged = (target, at) => (isPlainObject(config[at]) || isUnsetValue(config[at])
    ? mergeStoredSettings(target, config[at], at) : [{ at, stored: config[at], used: target, kept: true }]);
  const losses = [
    ...merged(targets.form, 'form'),
    ...merged(targets.adv, 'adv'),
    ...merged(targets.losat, 'losat'),
    ...merged(targets.circularConservation, 'circularConservation')
  ].filter(({ at }) => !RAW_FIELDS.includes(at) && !at.startsWith('losat.blastp.hitLimitsByMode.'));
  CONTAINER_TABLES.forEach(({ at, rows }) => {
    const container = at === 'losat.blastp' ? targets.losat.blastp
      : /** @type {Record<string, unknown>} */ (targets[/** @type {'form' | 'adv' | 'circularConservation'} */ (at)]);
    losses.push(...readRows(container, rows, at, containerAt(config, at)));
  });
  // A stored group the merge did not use (its row says the target keeps its
  // value) gives no value below either.
  /** @param {string} at */
  const used = (at) => !losses.some((loss) => loss.kept && loss.at === at);
  const storedAdv = used('adv') ? containerAt(config, 'adv') : null;
  Object.assign(targets.adv, Object.fromEntries(RAW_ADV_FIELDS
    .filter((field) => storedAdv && Object.hasOwn(storedAdv, field)).map((field) => [field, copied(storedAdv?.[field])])));
  // The arrow ratios the stored `adv` holds, else their defaults (an older
  // file draws the arrow it was drawn with); an `adv` not used keeps the target's.
  Object.assign(targets.adv, normalizedPersistedArrowGeometryState(used('adv') ? storedAdv : targets.adv));
  if (config.losat && used('losat') && used('losat.blastp')) {
    // The hit limits of each mode over their defaults, as stored, unless
    // they are not used; the inference flag as stored, else on for a collinear search.
    const blastp = containerAt(config, 'losat.blastp');
    const limits = blastp?.hitLimitsByMode;
    targets.losat.blastp.collinearInferOrthogroups = blastp?.collinearInferOrthogroups ?? (targets.losat.blastp.mode === 'collinear');
    if (used('losat.blastp.hitLimitsByMode')) {
      targets.losat.blastp.hitLimitsByMode = { ...createDefaultLosatpHitLimits(), ...(isPlainObject(limits) ? copied(limits) : {}) };
    }
  }
  return losses;
};

/**
 * Reads `tables` of `value` as stored; returns what they could not restore.
 * @param {readonly FieldTable[]} tables
 * @param {unknown} value
 * @returns {SettingLoss[]}
 */
const readValueTables = (tables, value) => tables.flatMap(({ at, rows }) => {
  const source = containerAt(value, at);
  if (!source) return [];
  /** @type {Record<string, unknown>} */
  const read = { ...source };
  return readRows(read, rows.filter(({ field }) => Object.hasOwn(source, field)), at, source);
});

/** Free text: Load names the place of such a value, not the value (12.3). */
const FREE_TEXT_KEYS = new Set([
  'label', 'name', 'fileName', 'sourceKey', 'uid', 'caption', 'prefix', 'species', 'strain', 'plot_title',
  'circular_record_label', 'circular_record_subtitle', 'labels', 'blacklistText', 'legendLabel', 'description', 'id'
]);
/** @param {string} at */
const isFreeText = (at) => FREE_TEXT_KEYS.has(at.replace(/\[\d+\]/g, '').split('.').pop() || '');
/** @param {unknown} value */
const excerpt = (value) => (typeof value === 'string' ? value : JSON.stringify(value) ?? '');

/**
 * Whether the drawing of `mode` draws the value at `at`: the mode-scoped
 * registry gives some settings to one mode only, and a value the drawing
 * never draws changes nothing when Load drops it.
 * @param {'circular' | 'linear' | null} mode
 * @param {string} root The registry domain of the read value (`config`, or '' for a slice).
 * @param {string} at
 */
const drawnByMode = (mode, root, at) => {
  if (!mode) return true;
  const [head, next] = at.replace(/\[\d+\]/g, '').split('.');
  const setting = MODE_SCOPED_SETTINGS.rows.find((entry) => entry.domain === [root, head].filter(Boolean).join('.') && entry.path === next)
    ?? MODE_SCOPED_SETTINGS.rows.find((entry) => entry.domain === root && entry.path === head);
  return !setting || setting.modes === 'both' || setting.modes === mode;
};

/**
 * @param {unknown} value
 * @param {string} at
 */
const holds = (value, at) => {
  let container = value;
  for (const key of at.replace(/\[(\d+)\]/g, '.$1').split('.').filter(Boolean)) {
    if (!container || typeof container !== 'object' || !Object.hasOwn(container, key)) return false;
    container = /** @type {Record<string, unknown>} */ (container)[key];
  }
  return true;
};

/**
 * The dialog rows of what a reading lost. A value the file does not hold
 * (`stored`), one the drawing never draws, and one read as its own text
 * (`"true"` as true) are no row. A row names the value used when it is one
 * value; free text is named by its place only.
 * @param {readonly SettingLoss[]} losses
 * @param {string} path Where the read value is in the file.
 * @param {unknown} stored
 * @param {{ mode?: 'circular' | 'linear' | null, root?: string, targets?: 'current' | 'default' }} [scope]
 *   The drawing's mode, the registry domain of the read value, and what a
 *   kept value is: the drawing's current one (a configuration file is
 *   installed over the drawing) or its default (Load resets a Session's drawing).
 * @returns {LoadRepair[]}
 */
const lossRepairs = (losses, path, stored, { mode = null, root = '', targets = 'default' } = {}) => losses
  .filter(({ at, stored: storedValue, used }) => holds(stored, at) && drawnByMode(mode, root, at)
    && !keepsMeaning(storedValue, used))
  .map(({ at, stored: storedValue, used, dropped, kept }) => {
    const fieldPath = [path, at].filter(Boolean).join('.');
    if (dropped) return loadRepair({ code: 'setting-entries-dropped', path: fieldPath, effect: 'dropped', count: dropped });
    // A kept value: the drawing's current one or its default; else nothing
    // used: the value is dropped; an unset value, a list or an object used:
    // the default is used. Free text is named by its place.
    const named = !isUnsetValue(used) && typeof used !== 'object';
    const keptEffect = kept ? (targets === 'current' ? 'kept' : 'defaulted') : null;
    if (isFreeText(at)) {
      return loadRepair({
        code: 'setting-replaced', path: fieldPath,
        effect: keptEffect ?? (used === undefined ? 'dropped' : named ? 'changed' : 'defaulted')
      });
    }
    return loadRepair({
      code: 'setting-replaced', path: fieldPath, effect: keptEffect ?? (used === undefined ? 'dropped' : 'defaulted'),
      stored: excerpt(storedValue), ...(named ? { used: excerpt(used) } : {})
    });
  });

/**
 * Reads the containers of `config` as `installStoredSettings` installs them
 * into `targets` (fresh defaults for a Session's drawing, a copy of the
 * drawing for a legacy configuration) and the other fields as
 * `applyConfigData` reads them; lists what they cannot restore. Only the
 * values `stored` holds are listed (the file's draft; a value the stored
 * request supplies is read with the request), and only those the drawing of
 * `mode` draws.
 * @param {SettingsContainers} targets Written.
 * @param {Record<string, unknown>} config
 * @param {string} path Where the configuration is in the file.
 * @param {{ stored?: unknown, mode?: 'circular' | 'linear' | null, current?: boolean }} [options]
 *   `current`: `targets` are the drawing's current settings (a configuration
 *   file), so a value not used keeps the drawing's.
 * @returns {LoadRepair[]}
 */
export const configFieldRepairs = (targets, config, path, { stored = config, mode = null, current = false } = {}) => lossRepairs([
  ...installStoredSettings(targets, config),
  ...readValueTables(VALUE_TABLES, config)
], path, stored, { mode, root: 'config', targets: current ? 'current' : 'default' });

/**
 * The same for a Session's app `ui`.
 * @param {unknown} ui
 * @param {string} path
 */
export const sessionUiFieldRepairs = (ui, path) => lossRepairs(readValueTables(SESSION_UI_TABLE, ui), path, ui);

/**
 * The same for a drawing's `ui` and feature edits (a Session 46 slice).
 * @param {unknown} slice
 * @param {string} path
 * @param {'circular' | 'linear'} mode
 */
export const sliceFieldRepairs = (slice, path, mode) => lossRepairs(readValueTables(SLICE_TABLES, slice), path, slice, { mode });

/**
 * The same for one value and its reader.
 * @param {FieldReader<unknown>} read
 * @param {unknown} value
 * @param {string} path
 */
export const valueRepairs = (read, value, path) => {
  const reading = read(value);
  return reading.lost ? lossRepairs([{ at: '', stored: value, used: reading.value }], path, value) : [];
};

/**
 * The same for the LOSAT execution settings (`ui.losatExecution` of a Session
 * 46, `config.losat` of an older draft).
 * @param {unknown} execution
 * @param {string} path
 */
export const losatExecutionRepairs = (execution, path) => lossRepairs(
  readValueTables([{ at: '', rows: LOSAT_EXECUTION_ROWS }], execution), path, execution
);
