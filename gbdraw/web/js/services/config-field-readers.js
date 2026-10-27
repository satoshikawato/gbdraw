// @ts-check
import { MODE_SCOPED_SETTINGS } from '../mode-scoped-settings.generated.js';
import { LABEL_STYLE_CHOICES, labelStyleForMode } from '../mode-profiles.js';
import { CANONICAL_REQUEST_DESCRIPTION } from '../contracts/request-shape.generated.js';
import { readCommittedPaint, resolveColorToHex } from '../utils/color-utils.js';
import { classifyOptionalNumber } from '../utils/optional-positive-number.js';
import {
  changedLeaves, exactChoice, isNonNegativeInteger, isPositiveInteger, isUnsetValue, keepsMeaning, nonNegativeNumberOrNull, numberOr, onlyTrue,
  positiveNumberOrNull, textChoice, trueUnlessFalse, truthy
} from '../utils/field-reading.js';
import { loadRepair } from '../utils/load-repairs.js';
import { DECIMAL_NUMBER_PATTERN } from '../utils/optional-positive-number.js';
import { normalizeAnnotationSets } from './annotation-state.js';
import {
  CIRCULAR_MULTI_RECORD_SIZE_MODE_ALIASES, CIRCULAR_MULTI_RECORD_SIZE_MODES, LINEAR_LABEL_PLACEMENT_ALIASES, LINEAR_LABEL_PLACEMENTS,
  LINEAR_TRACK_LAYOUT_ALIASES, LINEAR_TRACK_LAYOUTS, PAIRWISE_MATCH_STYLES, PROTEIN_BLASTP_MODES, migratePersistedWebStateFieldNames
} from './current-option-values.js';
import { readDefinitionLineStyles } from './definition-line-style-state.js';
import { readDepthTracks } from './depth-track-state.js';
import { IMPORTED_COMPARISON_ACTIONS } from './imported-comparison-intent.js';
import { readCircularPlotTitlePosition } from './layout-preferences.js';
import { isPlainObject } from './losat-cache.js';
import {
  readCollinearAnchorMode, readCollinearColorMode, readCollinearSearchScope, readCollinearUnitMode, readOrthogroupMembershipMode
} from './losat-normalization.js';
import { normalizeOrthogroupDormantOverrides } from './orthogroup-feature-metadata.js';
import {
  CIRCULAR_TRACK_SLOT_SCHEMA_VERSION, DEFAULT_FEATURE_OVERLAP_TOLERANCE_BP, DEFAULT_PNG_DPI, PNG_DPI_CHOICES, RETIRED_STORED_FIELDS,
  SESSION_INPUT_TYPES, STORED_CONFIG_DOMAINS, STORED_ROW_FIELDS, createDefaultAdv, createDefaultCircularConservation, createDefaultForm,
  createDefaultLosat, createDefaultLosatpHitLimits, readImportedTrackSlots, readLinearTypographyLinked
} from './session-active-config-contract.js';
import {
  applyCircularTrackOrderPlacements, clampCircularTrackAxisIndex, inferLegacyAxisIndexFromFeature, normalizeCircularTrackSlots
} from './circular-track-slot-model.js';
import {
  LINEAR_TRACK_SLOT_SCHEMA_VERSION, applyLinearTrackOrderPlacements, clampLinearTrackAxisIndex, normalizeLinearTrackSlots, resolveLinearTrackAxisIndex
} from './linear-track-slot-model.js';
import { normalizeLinearComparisonPlan } from './linear-comparisons.js';
import { MODE_SLICE_CONTAINERS } from './mode-scoped-migration.js';
import { readRecordDisplayDrafts } from './record-display-model.js';
import { isUnsafeObjectKey } from './safe-object-keys.js';
import { isRecordDrawKey } from './record-draw-selection.js';
import { isFeatureOverrideDraftRow, isFeaturePlacementDraftRow } from './feature-placement.js';
import { WEB_UX_PROFILE } from '../web-ux-profile.js';
import {
  DEFAULT_FEATURE_RENDERINGS, FEATURE_RENDERING_VALUES, readPersistedArrowHeadLengthRatio, readPersistedArrowShaftWidthRatio
} from '../utils/feature-rendering.js';
import { LINEAR_LABEL_VISIBILITY_MODES } from './linear-label-visibility.js';
import { featureColorEditColor } from './specific-color-rules.js';
import { isLabelOverrideRow } from './label-override-table.js';
import { isCardRow, markedRows } from './file-imports.js';

/** @import { FieldLeaf, FieldReader, FieldReading } from '../utils/field-reading.js' */
/** @import { CardMarkedRow, FileImports, PriorityRuleRow, WhitelistRuleRow } from './file-imports.js' */
/** @import { LoadRepair } from '../utils/load-repairs.js' */
/** @import { ObjectShape } from '../utils/contract-shape.js' */

// The readers of the stored settings of a configuration (a drawing's draft, a
// legacy configuration file) and of a Session's app `ui`, a drawing's `ui`
// and its feature edits. `applyConfigData` and Load's apply read with
// these readers, History restores included; before Load applies a file it
// reads the file with the same readers and lists what they cannot restore
// (BA-D25, B36, B37, B43), so the apply makes no repair of its own (G2). A
// value Save refuses is read here, not refused (E1a): the readers share their
// predicates with Save's checks (E-OD02).

/**
 * One stored field and its reader. A reader may read fields its table read
 * before it (`container`). `inPlace`: the list is replaced in its array.
 * `retired`: the field is no longer read and is deleted. `fileOnly`: only a
 * file's value is read; captured settings keep the value the app wrote (R11).
 * @typedef {object} FieldRow
 * @property {string} field
 * @property {(value: unknown, container: Record<string, unknown>) => FieldReading<unknown>} read
 * @property {boolean} [inPlace]
 * @property {boolean} [retired]
 * @property {boolean} [fileOnly]
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
    // A series without its file index (unset or blank) takes its place in the list.
    const sourceIndex = numberOr(Number.isInteger, series.length)(entry.sourceIndex);
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
 * A crop endpoint of a Linear card: a whole number from 1 (its text reads as
 * it), else unset (OV-568).
 * @type {FieldReader<number | null>}
 */
const readRegionEndpoint = (value) => {
  const classified = classifyOptionalNumber(value);
  if (classified.status === 'auto') return { value: null, lost: false };
  return classified.status === 'valid' && isPositiveInteger(classified.value)
    ? { value: classified.value, lost: false } : { value: null, lost: true };
};
/**
 * A Linear card's orientation: a boolean; any other value is not used, and
 * the card is drawn as stored (OV-568).
 * @type {FieldReader<boolean>}
 */
const readRegionReverse = (value) => (typeof value === 'boolean'
  ? { value, lost: false } : { value: false, lost: !isUnsetValue(value) });
/** A Linear card's text: text, else none (OV-568). @type {FieldReader<string>} */
const readCardText = (value) => (typeof value === 'string' ? { value, lost: false } : { value: '', lost: !isUnsetValue(value) });
/** A Linear card's genetic code, crop and texts. @type {readonly FieldRow[]} */
const LINEAR_SEQ_ROWS = Object.freeze([
  row('losat_gencode', readGencode),
  row('region_start', readRegionEndpoint),
  row('region_end', readRegionEndpoint),
  row('region_reverse', readRegionReverse),
  ...['definition', 'record_subtitle', 'file_definition', 'file_subtitle', 'inferred_definition', 'region_record_id']
    .map((field) => row(field, readCardText))
]);
/** The Linear card fields the Files read (`readLinearSeqFields`); a binding copies them as stored. */
export const LINEAR_SEQ_READ_FIELDS = Object.freeze(LINEAR_SEQ_ROWS.map(({ field }) => field));
/**
 * A Linear card's genetic code, crop and texts as the Files read them
 * (`files.linearSeqs[i]` of a Session 27-30, `webFiles.bindings.linearSeqs[i]`
 * after), with the rows of what they cannot restore at `path` (OV-568).
 * @typedef {{ losat_gencode: number, region_start: number | null, region_end: number | null, region_reverse: boolean,
 *   definition: string, record_subtitle: string, file_definition: string, file_subtitle: string,
 *   inferred_definition: string, region_record_id: string }} LinearSeqFields
 * @param {Record<string, unknown>} seq
 * @param {string} path
 * @returns {{ value: LinearSeqFields, repairs: LoadRepair[] }}
 */
export const readLinearSeqFields = (seq, path) => {
  /** @type {Record<string, unknown>} */
  const read = { ...seq };
  const losses = readRows(read, LINEAR_SEQ_ROWS, '', seq);
  return {
    value: /** @type {LinearSeqFields} */ (Object.fromEntries(LINEAR_SEQ_READ_FIELDS.map((field) => [field, read[field]]))),
    repairs: lossRepairs(losses, path, seq)
  };
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

/** @typedef {PriorityRuleRow} QualifierPriorityRule */
/**
 * One rule per feature type, the later rule of a type winning; a rule without
 * a type or an order is dropped. A rule the Priority File added keeps its
 * card and the user order it replaced (OV-1280); another mark is not read.
 * @type {FieldReader<QualifierPriorityRule[]>}
 */
export const readQualifierPriorityRules = (value) => {
  /** @type {QualifierPriorityRule[]} */
  const rules = [];
  /** @type {FieldLeaf[]} */
  const leaves = [];
  let dropped = 0;
  (Array.isArray(value) ? value : []).forEach((rule, index) => {
    const feat = isPlainObject(rule) ? String(rule.feat ?? '').trim() : '';
    const order = isPlainObject(rule) ? String(rule.order ?? '').trim() : '';
    if (!feat || !order) {
      dropped += 1;
      return;
    }
    const entry = /** @type {Record<string, unknown>} */ (rule);
    const imported = isCardRow(/** @type {CardMarkedRow} */ (entry), 'qualifier_priority');
    const handOrder = imported && typeof entry.handOrder === 'string' ? entry.handOrder : undefined;
    // A field no gbdraw writer wrote is not read (E-OD03).
    Object.keys(entry).filter((field) => field !== 'feat' && field !== 'order'
      && !(field === 'importedFrom' && imported) && !(field === 'handOrder' && handOrder !== undefined))
      .forEach((field) => leaves.push({ at: `[${index}].${field}`, stored: entry[field], used: undefined }));
    /** @type {QualifierPriorityRule} */
    const read = { ...(imported ? markedRows([{ feat, order }], 'qualifier_priority')[0] : { feat, order }),
      ...(handOrder !== undefined ? { handOrder } : {}) };
    const existing = rules.findIndex((other) => other.feat === feat);
    if (existing >= 0) {
      rules.splice(existing, 1, read);
      dropped += 1;
    } else {
      rules.push(read);
    }
  });
  return listReading(rules, value, dropped, leaves);
};

/**
 * A list of rows read entry by entry: an entry that is no object is dropped
 * and counted; `readEntry` gives the row and the values of the entry it
 * replaced. A field not in `fields` is not read (A17b).
 * @template T
 * @param {(entry: Record<string, unknown>) => { row: T, leaves: FieldLeaf[] }} readEntry
 * @param {readonly string[]} fields
 * @returns {FieldReader<T[]>}
 */
const rowList = (readEntry, fields) => (value) => {
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
    const unread = Object.keys(entry).filter((field) => !fields.includes(field))
      .map((field) => ({ at: `.${field}`, stored: entry[field], used: undefined }));
    [...reading.leaves, ...unread].forEach((leaf) => leaves.push({ ...leaf, at: `[${index}]${leaf.at}` }));
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
}, STORED_ROW_FIELDS.rules));

/** @typedef {WhitelistRuleRow} StoredWhitelistRule */
/**
 * The label whitelist of a configuration (`whitelist`); a value that is no
 * list is not used (`listOrKept`). A row the Whitelist File added keeps its
 * card (OV-1280).
 * @type {FieldReader<StoredWhitelistRule[]>}
 */
export const readWhitelistRules = listOrKept(rowList((entry) => {
  const imported = isCardRow(/** @type {CardMarkedRow} */ (entry), 'whitelist');
  const read = { feat: ruleText(entry.feat), qual: ruleText(entry.qual), key: ruleText(entry.key) };
  /** @type {StoredWhitelistRule} */
  const row = imported ? markedRows([read], 'whitelist')[0] : read;
  return { row, leaves: [...textFieldLeaves(entry, row, ['feat', 'qual', 'key']),
    ...(!isUnsetValue(entry.importedFrom) && !imported ? [{ at: '.importedFrom', stored: entry.importedFrom, used: undefined }] : [])] };
}, STORED_ROW_FIELDS.whitelist));

/** @param {unknown} value */
const isColorText = (value) => typeof value === 'string' && value.trim() !== '';
/**
 * What the Default colors and Blacklist files set (`fileImports`, OV-1280):
 * each Default color with the user color it replaced (`null`: the palette's),
 * and the Blacklist words with the text before. An entry of another shape is
 * not read.
 * @type {FieldReader<FileImports>}
 */
export const readFileImports = (value) => {
  if (!isPlainObject(value)) return { value: {}, lost: !isUnsetValue(value) };
  /** @type {FileImports} */
  const read = {};
  let lost = Object.keys(value).some((card) => card !== 'd_color' && card !== 'blacklist');
  if (isPlainObject(value.d_color)) {
    const colors = Object.entries(value.d_color).filter(([, entry]) => isPlainObject(entry) && isColorText(entry.color)
      && (entry.handColor === null || isColorText(entry.handColor)));
    lost ||= colors.length < Object.keys(value.d_color).length;
    if (colors.length) {
      read.d_color = Object.fromEntries(colors.map(([type, entry]) => [type, {
        color: String(/** @type {Record<string, unknown>} */ (entry).color),
        handColor: /** @type {string | null} */ (/** @type {Record<string, unknown>} */ (entry).handColor)
      }]));
    }
  } else lost ||= !isUnsetValue(value.d_color);
  const blacklist = value.blacklist;
  if (isPlainObject(blacklist) && Array.isArray(blacklist.words) && blacklist.words.every((word) => typeof word === 'string')
    && typeof blacklist.before === 'string') {
    read.blacklist = { words: [...blacklist.words], before: blacklist.before };
  } else lost ||= !isUnsetValue(blacklist);
  return { value: read, lost };
};

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
 * without case and outer spaces, a feature type without outer spaces. A
 * default feature type whose shape names none takes its default shape; any
 * other entry whose type is empty or whose shape names none is dropped
 * (BA-OD100, OV-757): the merge copies the shape of a feature type the
 * defaults lack as stored.
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
    if (featureType && FEATURE_RENDERING_VALUES.includes(name)) {
      shapes[featureType] = name;
      return;
    }
    // A default feature type keeps its default shape (OV-757).
    const fallback = Object.hasOwn(DEFAULT_FEATURE_RENDERINGS, featureType)
      ? /** @type {Readonly<Record<string, string>>} */ (DEFAULT_FEATURE_RENDERINGS)[featureType] : undefined;
    if (fallback) shapes[featureType] = fallback;
    // An entry of no feature type, or an unset shape, drew nothing: no row.
    if (featureType && !isUnsetValue(shape)) leaves.push({ at: `.${key}`, stored: shape, used: fallback });
  });
  return { value: shapes, lost: leaves.length > 0, ...(leaves.length ? { leaves } : {}) };
};
// A stored visibility the merge installs names a mode (`requireCurrentConfigFields`).
const readLinearLabelVisibility = textChoice(LINEAR_LABEL_VISIBILITY_MODES, 'auto');

/**
 * `adv` of a drawing of `mode`, in the order the rows depend on each other.
 * @param {'circular' | 'linear'} mode
 * @returns {readonly FieldRow[]}
 */
const advRows = (mode) => Object.freeze([
  row('label_rendering', normalizeLabelRendering),
  // A value of no style reads as the mode's packaged style (the generated mode profile), and is reported
  // (labels.<mode>.line, .protection).
  row('label_line', textChoice(LABEL_STYLE_CHOICES.line, labelStyleForMode(mode).label_line)),
  row('label_protection', textChoice(LABEL_STYLE_CHOICES.protection, labelStyleForMode(mode).label_protection)),
  row('label_placement', textChoice(LINEAR_LABEL_PLACEMENTS, 'auto', { aliases: LINEAR_LABEL_PLACEMENT_ALIASES })),
  row('multi_record_size_mode', textChoice(CIRCULAR_MULTI_RECORD_SIZE_MODES, 'auto', { aliases: CIRCULAR_MULTI_RECORD_SIZE_MODE_ALIASES })),
  row('circular_label_placement', textChoice(['radial', 'horizontal'], 'horizontal')),
  row('circular_label_spacing', positiveNumberOrNull),
  row('linear_label_spacing', positiveNumberOrNull),
  row('track_axis_gap', readNonNegativeNumberOrNull),
  row('multi_record_min_radius_ratio', numberOr(isRadiusRatio, 0.55)),
  row('multi_record_column_gap_ratio', numberOr(isNonNegative, 0.10)),
  row('multi_record_row_gap_ratio', numberOr(isNonNegative, 0.05)),
  row('multi_record_positions', readMultiRecordPositions),
  row('feature_overlap_tolerance_bp', numberOr((numeric) => Number.isSafeInteger(numeric) && numeric >= 0, DEFAULT_FEATURE_OVERLAP_TOLERANCE_BP)),
  row('plot_title_font_size', positiveNumberOrNull),
  row('scale_font_size', positiveNumberOrNull),
  row('ruler_label_font_size', positiveNumberOrNull),
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
  row('linear_definition_line_styles', readDefinitionLineStyles),
  // The track stacks (OV-758): each axis index is read before its stack, as a
  // stack read as the default unsets it; the stacks are read in the current
  // slot schema. A captured stack keeps the drafts its rows show as invalid
  // (TK-12), so only a file's stack is read.
  .../** @type {const} */ ([['circular', CIRCULAR_TRACK_SLOT_SCHEMA_VERSION], ['linear', LINEAR_TRACK_SLOT_SCHEMA_VERSION]]).flatMap(([mode, schema]) => [
    {
      field: `${mode}_track_slots_axis_index`,
      /** @type {FieldRow['read']} */
      read: (value, adv) => (readImportedTrackSlots(mode, adv[`${mode}_track_slots`]).reset
        ? { value: null, lost: !isUnsetValue(value) } : { value, lost: false }),
      fileOnly: true
    },
    { field: `${mode}_track_slots`, read: (/** @type {unknown} */ value) => readImportedTrackSlots(mode, value), inPlace: true, fileOnly: true },
    // A version up to the current one is migrated without a loss.
    {
      field: `${mode}_track_slots_schema_version`,
      read: (/** @type {unknown} */ value) => ({
        value: schema, lost: !isUnsetValue(value) && !(typeof value === 'number' && Number.isInteger(value) && value >= 1 && value <= schema)
      })
    }
  ])
]);

/**
 * `losat.blastp`: the protein search settings. `orthogroupMaxHits` had no
 * writer on `main` or in a release (OV-756): it is dropped like the other
 * retired keys.
 * @type {readonly FieldRow[]}
 */
const LOSAT_SEARCH_ROWS = Object.freeze([
  row('mode', textChoice(PROTEIN_BLASTP_MODES, 'orthogroup')),
  row('maxHits', numberOr(isPositiveInteger, 5)),
  row('candidateLimit', numberOr(isPositiveInteger, null)),
  row('orthogroupMembershipMode', readOrthogroupMembershipMode),
  row('orthogroupMemberMaxHits', numberOr(isPositiveInteger, null)),
  row('collinearMinAnchors', numberOr(isPositiveInteger, 1)),
  row('collinearMaxUnitGap', numberOr(isNonNegativeInteger, 0)),
  row('collinearMaxDiagonalDrift', numberOr(isNonNegativeInteger, 0)),
  row('collinearMaxConflictsInMergeGap', numberOr(isNonNegativeInteger, 1)),
  row('collinearMaxParalogLinksPerOrthogroup', numberOr(isPositiveInteger, 2)),
  row('collinearColorMode', readCollinearColorMode),
  row('collinearUnitMode', readCollinearUnitMode),
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
/** The Linear record gap's default, in px: Python's (`CanonicalLinearLayout.recordGapPx`). */
export const DEFAULT_LINEAR_RECORD_GAP_PX = /** @type {number} */ (
  /** @type {ObjectShape} */ (CANONICAL_REQUEST_DESCRIPTION.types.CanonicalLinearLayout).defaults?.recordGapPx);
/**
 * The Linear record gap a draft draws and saves with: a blank (the cleared
 * number input stores '') or an unset value is the default, as R7 reads a
 * blank (OV-828, BA-D30); any other value is the draft's own, which Generate
 * and Save convert as before.
 * @param {unknown} value
 * @returns {unknown}
 */
export const draftLinearRecordGap = (value) => (isUnsetValue(value) ? DEFAULT_LINEAR_RECORD_GAP_PX : value);
/** A blank or absent gap reads as the default, as `draftLinearRecordGap` gives it. @type {FieldReader<number>} */
export const readLinearRecordGap = numberOr(isNonNegative, DEFAULT_LINEAR_RECORD_GAP_PX);
/** The Linear alignment plan: an object, else unset (OV-578). @type {FieldReader<Record<string, unknown> | null>} */
export const readSimilarityAlignment = (value) => (isPlainObject(value) ? { value, lost: false } : { value: null, lost: !isUnsetValue(value) });
export const readImportedComparisonAction = exactChoice(Object.values(IMPORTED_COMPARISON_ACTIONS), null);
const readLosatProgramText = exactChoice(['blastn', 'tblastx', 'blastp'], 'blastn');
/** The label filter, read without case. */
export const readFilterMode = textChoice(['None', 'Whitelist', 'Blacklist'], 'None', {
  aliases: { none: 'None', whitelist: 'Whitelist', blacklist: 'Blacklist' }
});
/**
 * The Label blacklist, text: a number reads as its text; another value keeps
 * the text it always read as (`String`, nothing for an unset or false value),
 * with its row (E-OD26). @type {FieldReader<string>}
 */
export const readBlacklistText = (value) => (typeof value === 'string' || typeof value === 'number'
  ? { value: String(value), lost: false } : { value: String(value || ''), lost: !isUnsetValue(value) });
/** Read as text. @type {FieldReader<string>} */
export const readLosatProgram = (value) => readLosatProgramText(value ? String(value) : value);

/** Text, else unset (the drawing's palette stays). @type {FieldReader<string>} */
export const readPalette = (value) => (typeof value === 'string' ? { value: value.trim(), lost: false }
  : { value: '', lost: !isUnsetValue(value) });
/**
 * The records a drawing of `mode` turns OFF: each a record key of its mode,
 * once (a repeat names the same record: no loss).
 * @param {'circular' | 'linear'} mode
 * @returns {FieldReader<string[]>}
 */
export const recordsOffReader = (mode) => (value) => {
  const entries = Array.isArray(value) ? value : [];
  const keys = /** @type {string[]} */ (entries.filter((key) => isRecordDrawKey(mode, key)));
  // A key of the other mode's records (a Session 27-44 draft both read) turned nothing off here.
  const otherMode = mode === 'linear' ? 'circular' : 'linear';
  const unread = entries.filter((key) => !isRecordDrawKey(mode, key) && !isRecordDrawKey(otherMode, key)).length;
  return listReading([...new Set(keys)], value, unread, []);
};
/**
 * The rows of an identity-keyed draft map that `accepts`.
 * @param {(key: string, row: unknown) => boolean} accepts
 * @returns {FieldReader<Record<string, unknown>>}
 */
const draftRows = (accepts) => (value) => {
  /** @type {Record<string, unknown>} */
  const rows = {};
  let dropped = 0;
  Object.entries(isPlainObject(value) ? value : {}).forEach(([key, row]) => {
    if (accepts(key, row)) rows[key] = row;
    else dropped += 1;
  });
  const notObject = !isPlainObject(value) && !isUnsetValue(value);
  return { value: rows, lost: notObject || dropped > 0, ...(dropped && !notObject ? { dropped } : {}) };
};
/**
 * The Feature placement draft rows of a drawing of `mode` that
 * `canonicalFeaturePlacements` accepts (OV-764).
 * @param {'circular' | 'linear'} mode
 */
export const featurePlacementsReader = (mode) => draftRows((key, row) => isFeaturePlacementDraftRow(key, row, mode));
/** The feature edit draft rows `canonicalFeatureOverrides` accepts (A23'). */
export const readFeatureOverrides = draftRows(isFeatureOverrideDraftRow);
/**
 * A Linear comparison plan as `normalizeLinearComparisonPlan` reads it, with
 * what it does not use as stored (OV-764).
 * @type {FieldReader<ReturnType<typeof normalizeLinearComparisonPlan>>}
 */
export const readLinearComparisonPlan = (value) => {
  /** @type {FieldLeaf[]} */
  const leaves = [];
  const plan = normalizeLinearComparisonPlan(value, (leaf) => leaves.push(leaf));
  const notObject = !isPlainObject(value) && !isUnsetValue(value);
  return { value: plan, lost: notObject || leaves.length > 0, ...(leaves.length && !notObject ? { leaves } : {}) };
};
const readHitLimit = numberOr(isPositiveInteger, null);
/**
 * The hit limits of each protein search mode (`losat.blastp.hitLimitsByMode`)
 * over their defaults: a mode stored as null keeps its defaults; a mode or a
 * limit of another kind, or one no mode has, is not used (OV-764).
 * @type {FieldReader<Record<string, Record<string, number | null>>>}
 */
export const readHitLimitsByMode = (value) => {
  /** @type {Record<string, Record<string, number | null>>} */
  const limits = createDefaultLosatpHitLimits();
  /** @type {FieldLeaf[]} */
  const leaves = [];
  Object.entries(isPlainObject(value) ? value : {}).forEach(([mode, stored]) => {
    if (!Object.hasOwn(limits, mode)) {
      leaves.push({ at: `.${mode}`, stored, used: undefined });
      return;
    }
    if (stored === null) return;
    if (!isPlainObject(stored)) {
      leaves.push({ at: `.${mode}`, stored, used: limits[mode] });
      return;
    }
    // A stored mode replaces its defaults, as it always did.
    const defaults = limits[mode];
    limits[mode] = {};
    Object.entries(stored).forEach(([field, limit]) => {
      if (!Object.hasOwn(defaults, field)) {
        leaves.push({ at: `.${mode}.${field}`, stored: limit, used: undefined });
        return;
      }
      const reading = readHitLimit(limit);
      if (reading.lost) leaves.push({ at: `.${mode}.${field}`, stored: limit, used: reading.value });
      limits[mode][field] = reading.value;
    });
  });
  const notObject = !isPlainObject(value) && !isUnsetValue(value);
  return { value: limits, lost: notObject || leaves.length > 0, ...(leaves.length && !notObject ? { leaves } : {}) };
};

/**
 * The containers the merge reads first, then their rows (`installStoredSettings`), of a drawing of `mode`.
 * @param {'circular' | 'linear'} mode
 * @returns {readonly FieldTable[]}
 */
const containerTables = (mode) => Object.freeze([
  { at: 'form', rows: [
    row('linear_track_layout', textChoice(LINEAR_TRACK_LAYOUTS, 'middle', { aliases: LINEAR_TRACK_LAYOUT_ALIASES })),
    // A null lock is an omission, the historical OFF (`importSessionDocument`; E-OD24, no row); the merge keeps
    // the target's value for one of another kind, with its row.
    row('keep_definition_left_aligned', onlyTrue)
  ] },
  { at: 'adv', rows: advRows(mode) },
  { at: 'losat.blastp', rows: LOSAT_SEARCH_ROWS },
  { at: 'circularConservation', rows: CONSERVATION_ROWS }
]);
/** @type {Readonly<Record<'circular' | 'linear', readonly FieldTable[]>>} */
const CONTAINER_TABLES = Object.freeze({ circular: containerTables('circular'), linear: containerTables('linear') });
/** The containers `installStoredSettings` merges. */
const MERGED_CONTAINERS = Object.freeze(['form', 'adv', 'losat', 'circularConservation']);
/**
 * The fields of a drawing of `mode` read from the configuration as stored.
 * @type {Readonly<Record<'circular' | 'linear', readonly FieldTable[]>>}
 */
const DRAWING_TABLES = Object.freeze({
  circular: [{ at: '', rows: [row('recordsOff', recordsOffReader('circular')), row('featurePlacementOverrides', featurePlacementsReader('circular'))] }],
  linear: [{ at: '', rows: [row('recordsOff', recordsOffReader('linear')), row('featurePlacementOverrides', featurePlacementsReader('linear'))] }]
});
/** The fields read from the configuration as stored. @type {readonly FieldTable[]} */
const VALUE_TABLES = Object.freeze([
  { at: 'losat', rows: LOSAT_EXECUTION_ROWS },
  {
    at: 'linearRecordLayout',
    rows: [
      row('enabled', readLinearRecordLayoutEnabled),
      row('recordGap', readLinearRecordGap),
      row('rows', readLinearRecordRows),
      row('similarityAlignment', readSimilarityAlignment),
      row('recordTranslations', readList)
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
      row('blacklistText', readBlacklistText),
      row('fileImports', readFileImports),
      row('palette', readPalette),
      row('colors', readPlainObject),
      row('unmanagedConfigOverrides', readPlainObject),
      row('recordDisplayDrafts', (value) => readRecordDisplayDrafts(value)),
      row('linearComparisonPlan', readLinearComparisonPlan)
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
export const readPaletteColors = (value) => {
  if (!isPlainObject(value)) return { value: null, lost: !isUnsetValue(value) };
  // A color that names no color is dropped (OV-633).
  const entries = Object.entries(value).filter(([, color]) => readCommittedPaint(color) !== null);
  const dropped = Object.keys(value).length - entries.length;
  return {
    value: Object.fromEntries(entries.map(([key, color]) => [key, resolveColorToHex(String(color || '').trim())])),
    lost: dropped > 0, ...(dropped ? { dropped } : {})
  };
};
/** Shown unless `false`. @type {FieldReader<boolean>} */
export const readRichFeaturePopup = trueUnlessFalse;
export const readSessionUiMode = exactChoice(['circular', 'linear'], 'circular');
/** An input type, else GenBank (OV-594). */
export const readSessionInputType = exactChoice(SESSION_INPUT_TYPES, 'gb');
/**
 * The PNG DPI: one of the choices the DPI menu lists, kept as stored, else
 * the default (OV-595).
 * @type {FieldReader<unknown>}
 */
export const readDownloadDpi = (value) => {
  if (isUnsetValue(value)) return { value: DEFAULT_PNG_DPI, lost: false };
  const known = PNG_DPI_CHOICES.some((choice) => typeof value !== 'object' && choice.value === String(value).trim());
  return known ? { value: typeof value === 'number' ? value : String(value).trim(), lost: false } : { value: DEFAULT_PNG_DPI, lost: true };
};
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
      row('cInputType', readSessionInputType),
      row('lInputType', readSessionInputType),
      row('downloadDpi', readDownloadDpi),
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
/** Per-feature color edits: an edit whose paint names no color is dropped (RL-D). */
export const readFeatureColorEdits = draftRows((_key, edit) => readCommittedPaint(featureColorEditColor(edit)) !== null);
/** The label table rows the table serializes (RL-D). @type {FieldReader<unknown[]>} */
export const readLabelOverrideRows = (value) => {
  const entries = Array.isArray(value) ? value : [];
  const rows = entries.filter(isLabelOverrideRow);
  return listReading(rows, value, entries.length - rows.length, []);
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
      row('featureColorOverrides', readFeatureColorEdits),
      row('featureOverrides', readFeatureOverrides),
      row('featureVisibilityManualRules', readList),
      row('labelOverrideRows', readLabelOverrideRows),
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
  const container = valueAt(value, at);
  return isPlainObject(container) ? container : null;
};
/**
 * @param {unknown} value
 * @param {string} at Dotted; '' is the value itself.
 * @returns {unknown}
 */
const valueAt = (value, at) => {
  let container = value;
  for (const key of at.split('.').filter(Boolean)) container = isPlainObject(container) ? container[key] : undefined;
  return container;
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
  return readingLosses(reading, fieldAt, before);
});

/**
 * What a reading of the value stored at `at` lost.
 * @param {FieldReading<unknown>} reading
 * @param {string} at
 * @param {unknown} stored
 * @returns {SettingLoss[]}
 */
const readingLosses = (reading, at, stored) => {
  if (!reading.lost) return [];
  if (reading.kept) return [{ at, stored, used: undefined, kept: true }];
  if (reading.dropped || reading.leaves) {
    return [
      ...(reading.dropped ? [{ at, stored, used: undefined, dropped: reading.dropped }] : []),
      ...(reading.leaves ?? []).map((leaf) => ({ at: `${at}${leaf.at}`, stored: leaf.stored, used: leaf.used }))
    ];
  }
  return [{ at, stored, used: reading.value }];
};

/**
 * The comparison filters, whose writers store a number or its text
 * (`resolveComparisonThresholds` reads both): the e-value input is text, so any
 * text or number is its kind; the other filters are number inputs, so a number,
 * a text that reads as one, or a blank is theirs.
 * @type {readonly string[]}
 */
const TEXT_OR_NUMBER_FIELDS = Object.freeze(['adv.evalue']);
/** @type {readonly string[]} */
const NUMBER_OR_ITS_TEXT_FIELDS = Object.freeze(['adv.min_bitscore', 'adv.identity', 'adv.alignment_length']);
/**
 * A list, an object, or the `typeof` of anything else, but the kinds the
 * comparison filters hold count as one.
 * @param {unknown} value
 * @param {string} at
 */
const kindOf = (value, at) => {
  if (Array.isArray(value)) return 'list';
  if (TEXT_OR_NUMBER_FIELDS.includes(at) && (typeof value === 'string' || typeof value === 'number')) return 'text or number';
  if (NUMBER_OR_ITS_TEXT_FIELDS.includes(at) && (typeof value === 'number'
    || (typeof value === 'string' && (value.trim() === '' || DECIMAL_NUMBER_PATTERN.test(value.trim()))))) return 'number';
  return typeof value;
};
/** Stored maps whose keys the user names, merged with their own keys. @type {readonly string[]} */
const FREE_KEY_MAPS = Object.freeze(['adv.feature_shapes']);
/**
 * Stored keys the merge does not install and no loss of it: the retired and
 * moved fields, and the fields a row reads from the stored value (the LOSAT
 * execution settings of an older draft, a retired row).
 * @type {ReadonlySet<string>}
 */
const READ_ELSEWHERE = new Set([
  ...RETIRED_STORED_FIELDS,
  ...LOSAT_EXECUTION_ROWS.map(({ field }) => `losat.${field}`),
  ...LOSAT_SEARCH_ROWS.filter(({ retired }) => retired).map(({ field }) => `losat.blastp.${field}`)
]);
/**
 * Merges a stored settings object into `target`: a key `target` lacks is
 * not used, except in a map of `FREE_KEY_MAPS`, where it is copied; it is
 * returned as a loss unless a gbdraw writer wrote it and later retired or
 * moved it (`RETIRED_STORED_FIELDS`, E-OD03; TB19, TB32). An object merges
 * into an object, and `null` over an object merges nothing; a list replaces
 * a list; a value of `target`'s kind, or `null` on either side, is copied. An
 * empty text of another kind over a number (a cleared number input) is the
 * field's value in `defaults`, without a loss: a blank is the default (R7,
 * OV-827). Any other value of another kind (a list for an object included) is
 * not used; it is returned as a loss at `at`, and `target` keeps its value.
 * @param {Record<string, unknown>} target
 * @param {unknown} source
 * @param {string} at
 * @param {Record<string, unknown> | null} defaults The fresh settings of the same place (`defaultSettingsContainers`).
 * @returns {SettingLoss[]}
 */
const mergeStoredSettings = (target, source, at, defaults) => {
  if (!source || typeof source !== 'object') return [];
  const stored = /** @type {Record<string, unknown>} */ (source);
  return Object.keys(stored).flatMap((key) => {
    if (isUnsafeObjectKey(key)) return [];
    const keyAt = at ? `${at}.${key}` : key;
    if (!Object.hasOwn(target, key)) {
      // A shape of a feature type the defaults do not name is the user's own.
      if (FREE_KEY_MAPS.includes(at)) target[key] = copied(stored[key]);
      return FREE_KEY_MAPS.includes(at) || READ_ELSEWHERE.has(keyAt) ? [] : [{ at: keyAt, stored: stored[key], used: undefined }];
    }
    const targetValue = target[key];
    const sourceValue = stored[key];
    const fresh = defaults && Object.hasOwn(defaults, key) ? defaults[key] : undefined;
    // A null over a settings group holds nothing and merges nothing (OV-822, E-OD06).
    if (isPlainObject(targetValue) && (isPlainObject(sourceValue) || sourceValue === null)) {
      return mergeStoredSettings(targetValue, sourceValue, keyAt, isPlainObject(fresh) ? fresh : null);
    }
    if (Array.isArray(targetValue) && Array.isArray(sourceValue)) {
      targetValue.splice(0, targetValue.length, ...sourceValue);
      return [];
    }
    if (kindOf(targetValue, keyAt) === kindOf(sourceValue, keyAt) || targetValue === null || sourceValue === null) {
      target[key] = sourceValue;
      return [];
    }
    if (sourceValue === '' && typeof targetValue === 'number') {
      if (fresh !== undefined) target[key] = copied(fresh);
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
/**
 * Fields read from the stored value after the merge (their merge losses are
 * none): the arrow ratios (a ratio stored as its text reads as its number; one
 * out of range is its default, with a row), and the inference flag, read with
 * the hit limits below.
 */
const RAW_FIELDS = Object.freeze([
  'adv.arrow_head_length_ratio', 'adv.arrow_shaft_width_ratio', 'losat.blastp.collinearInferOrthogroups'
]);
/** @param {unknown} value */
const copied = (value) => (value === undefined ? undefined : JSON.parse(JSON.stringify(value)));

/** @typedef {{ arrow_head_length_ratio: number | null, arrow_shaft_width_ratio: number | null }} StoredArrowRatios */
/**
 * Installs a configuration's settings into `targets`: each container is
 * merged, then read by its rows. `applyConfigData` installs a drawing's
 * settings with it, and Load reads a file with it before it applies it, so
 * both lose the same values. The arrow ratios are read, not installed: the
 * caller that writes the drawing writes them (`applyConfigData`, R10).
 * @param {SettingsContainers} targets
 * @param {Record<string, unknown>} stored
 * @param {'circular' | 'linear'} mode The drawing's mode, whose defaults a blank reads as.
 * @param {{ captured?: boolean }} [options] `captured`: settings the app wrote (History Undo and Redo, the Load
 *   rollback), which the `fileOnly` rows leave as written (R11).
 * @returns {{ losses: SettingLoss[], arrowRatios: StoredArrowRatios | null }} `arrowRatios` is null when the stored
 *   `adv` is not used (the target keeps its ratios).
 */
export const installStoredSettings = (targets, stored, mode, { captured = false } = {}) => {
  // A retired field name reads as its current name at every version (A17g).
  const config = migratePersistedWebStateFieldNames(stored);
  const defaults = defaultSettingsContainers(mode);
  /**
   * A stored container that is no object is not used (the drawing's stays).
   * @param {'form' | 'adv' | 'losat' | 'circularConservation'} at
   */
  const merged = (at) => (isPlainObject(config[at]) || isUnsetValue(config[at])
    ? mergeStoredSettings(targets[at], config[at], at, defaults[at]) : [{ at, stored: config[at], used: targets[at], kept: true }]);
  const losses = [
    ...merged('form'),
    ...merged('adv'),
    ...merged('losat'),
    ...merged('circularConservation')
  ].filter(({ at }) => !RAW_FIELDS.includes(at) && !at.startsWith('losat.blastp.hitLimitsByMode.'));
  CONTAINER_TABLES[mode].forEach(({ at, rows }) => {
    const container = at === 'losat.blastp' ? targets.losat.blastp
      : /** @type {Record<string, unknown>} */ (targets[/** @type {'form' | 'adv' | 'circularConservation'} */ (at)]);
    losses.push(...readRows(container, captured ? rows.filter(({ fileOnly }) => !fileOnly) : rows, at, containerAt(config, at)));
  });
  // A stored group the merge did not use (its row says the target keeps its
  // value) gives no value below either.
  /** @param {string} at */
  const used = (at) => !losses.some((loss) => loss.kept && loss.at === at);
  // The arrow ratios the stored `adv` holds, else their defaults (an older
  // file draws the arrow it was drawn with); an `adv` not used keeps the target's.
  /** @type {StoredArrowRatios | null} */
  let arrowRatios = null;
  if (used('adv')) {
    const storedAdv = containerAt(config, 'adv');
    const head = readPersistedArrowHeadLengthRatio(storedAdv?.arrow_head_length_ratio);
    const shaft = readPersistedArrowShaftWidthRatio(storedAdv?.arrow_shaft_width_ratio);
    arrowRatios = { arrow_head_length_ratio: head.value, arrow_shaft_width_ratio: shaft.value };
    losses.push(...readingLosses(head, 'adv.arrow_head_length_ratio', storedAdv?.arrow_head_length_ratio),
      ...readingLosses(shaft, 'adv.arrow_shaft_width_ratio', storedAdv?.arrow_shaft_width_ratio));
  }
  if (config.losat && used('losat') && used('losat.blastp')) {
    // The inference flag as stored, else on for a collinear search; the hit
    // limits of each mode over their defaults, unless they are not used.
    const blastp = containerAt(config, 'losat.blastp');
    const inference = blastp?.collinearInferOrthogroups;
    targets.losat.blastp.collinearInferOrthogroups = typeof inference === 'boolean' ? inference : targets.losat.blastp.mode === 'collinear';
    if (!isUnsetValue(inference) && typeof inference !== 'boolean') {
      losses.push({ at: 'losat.blastp.collinearInferOrthogroups', stored: inference, used: targets.losat.blastp.collinearInferOrthogroups });
    }
    if (used('losat.blastp.hitLimitsByMode')) {
      const limits = readHitLimitsByMode(blastp?.hitLimitsByMode);
      targets.losat.blastp.hitLimitsByMode = limits.value;
      losses.push(...readingLosses(limits, 'losat.blastp.hitLimitsByMode', blastp?.hitLimitsByMode));
    }
  }
  return { losses, arrowRatios };
};

/**
 * Reads `tables` of `value` as stored; returns what they could not restore.
 * @param {readonly FieldTable[]} tables
 * @param {unknown} value
 * @returns {SettingLoss[]}
 */
const readValueTables = (tables, value) => tables.flatMap(({ at, rows }) => {
  const source = containerAt(value, at);
  if (!source) {
    // A container that is no object is not read (its fields take their
    // defaults); `installStoredSettings` lists the merged ones.
    const stored = at && !MERGED_CONTAINERS.includes(at) ? valueAt(value, at) : undefined;
    return isUnsetValue(stored) ? [] : [{ at, stored, used: {} }];
  }
  /** @type {Record<string, unknown>} */
  const read = { ...source };
  return readRows(read, rows.filter(({ field }) => Object.hasOwn(source, field)), at, source);
});

/** Free text: Load names the place of such a value, not the value (12.3). */
const FREE_TEXT_KEYS = new Set([
  'label', 'name', 'fileName', 'sourceKey', 'uid', 'caption', 'prefix', 'species', 'strain', 'plot_title',
  'circular_record_label', 'circular_record_subtitle', 'labels', 'blacklistText', 'legendLabel', 'description', 'id',
  // A track stack holds slot ids.
  'circular_track_slots', 'linear_track_slots'
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
 * @param {{ stored?: unknown, mode: 'circular' | 'linear', current?: boolean, resolveTrackPlacements?: boolean }} options
 *   `current`: `targets` are the drawing's current settings (a configuration
 *   file), so a value not used keeps the drawing's. `resolveTrackPlacements`:
 *   the apply resolves the track stacks (`resolvedTrackStacks`), so what that
 *   replaces is listed too.
 * @returns {LoadRepair[]}
 */
export const configFieldRepairs = (targets, config, path, {
  stored = config, mode, current = false, resolveTrackPlacements = true
}) => {
  const { losses } = installStoredSettings(targets, config, mode);
  if (resolveTrackPlacements) resolvedTrackStacks(targets.adv, targets.form, losses);
  return lossRepairs([
    ...losses,
    ...readValueTables(VALUE_TABLES, config),
    ...readValueTables(DRAWING_TABLES[mode], config),
    ...unknownDomainLosses(config)
  ], path, stored, { mode, root: 'config', targets: current ? 'current' : 'default' });
};

/**
 * The track stacks a drawing's settings resolve to: each stack with its axis
 * index read and its placements applied (`applyConfigData`'s
 * `resolveTrackPlacements`). Pure: Load reads the stacks it installs with it
 * before the commit, and `losses` collects what the slot readers replace
 * (OV-758; `at` under the settings).
 * @param {Record<string, unknown>} adv The drawing's settings, read.
 * @param {Record<string, unknown>} form
 * @param {SettingLoss[] | null} [losses]
 */
export const resolvedTrackStacks = (adv, form, losses = null) => {
  /** @type {FieldLeaf[]} */
  const circularLosses = [];
  /** @type {FieldLeaf[]} */
  const circularAxisLosses = [];
  /** @type {FieldLeaf[]} */
  const linearLosses = [];
  /** @type {FieldLeaf[]} */
  const linearAxisLosses = [];
  // The read settings hold text here (the reader rows); the casts name what the slot readers take.
  const nt = /** @type {string | undefined} */ (adv.nt);
  const trackType = /** @type {string | undefined} */ (form.track_type);
  const normalizedSlots = normalizeCircularTrackSlots(adv.circular_track_slots, nt, trackType, circularLosses);
  const importedAxis = clampCircularTrackAxisIndex(adv.circular_track_slots_axis_index, normalizedSlots.length, circularAxisLosses);
  const circularAxis = importedAxis === null
    ? inferLegacyAxisIndexFromFeature(normalizedSlots, trackType)
    : importedAxis;
  const linearStack = /** @type {Parameters<typeof normalizeLinearTrackSlots>[0]} */ (adv.linear_track_slots);
  const normalizedLinearSlots = normalizeLinearTrackSlots(linearStack, nt, form.linear_track_layout, linearLosses);
  const linearAxis = resolveLinearTrackAxisIndex(
    normalizedLinearSlots,
    clampLinearTrackAxisIndex(adv.linear_track_slots_axis_index, normalizedLinearSlots.length, linearAxisLosses)
  );
  losses?.push(...[
    ...circularLosses.map((leaf) => ({ ...leaf, at: `adv.circular_track_slots${leaf.at}` })),
    ...circularAxisLosses.map((leaf) => ({ ...leaf, at: 'adv.circular_track_slots_axis_index' })),
    ...linearLosses.map((leaf) => ({ ...leaf, at: `adv.linear_track_slots${leaf.at}` })),
    ...linearAxisLosses.map((leaf) => ({ ...leaf, at: 'adv.linear_track_slots_axis_index' }))
  ]);
  return {
    circularAxis,
    circularSlots: applyCircularTrackOrderPlacements(adv.circular_track_slots, nt, trackType, circularAxis),
    linearAxis,
    linearSlots: applyLinearTrackOrderPlacements(normalizedLinearSlots, linearAxis, nt, form.linear_track_layout)
  };
};

/**
 * The fields of a Session's stored draft its drawings do not install: they
 * install the registry fields of a container that is an object only
 * (`overlayModeSliceConfig`, the Session 27-44 split), so Load drops a field
 * no gbdraw writer wrote and a container of another kind, listed as
 * `configFieldRepairs` lists them in a configuration file (TB19, TB32).
 * @param {unknown} stored The stored draft (a slice's, or a Session 27-44's).
 * @param {string} path
 * @param {'circular' | 'linear'} mode
 * @param {{ slice?: boolean }} [options] `slice`: a Session 46 slice, whose
 *   `losat` holds no LOSAT execution setting (the app `ui` does; E-OD03).
 * @returns {LoadRepair[]}
 */
export const unreadFieldRepairs = (stored, path, mode, { slice = false } = {}) => {
  if (!isPlainObject(stored)) return [];
  const config = migratePersistedWebStateFieldNames(stored);
  const targets = defaultSettingsContainers(mode);
  const storedLosat = slice && isPlainObject(config.losat) ? config.losat : {};
  const losses = [
    ...MODE_SCOPED_SETTINGS.losatExecutionFields.filter((/** @type {string} */ field) => Object.hasOwn(storedLosat, field))
      .map((/** @type {string} */ field) => ({ at: `losat.${field}`, stored: storedLosat[field], used: undefined })),
    ...MERGED_CONTAINERS.flatMap((at) => {
      const stored = config[at];
      if (!isPlainObject(stored)) return isUnsetValue(stored) ? [] : [{ at, stored, used: {} }];
      const target = /** @type {Record<string, unknown>} */ (targets[/** @type {keyof SettingsContainers} */ (at)]);
      return mergeStoredSettings(target, stored, at, null).filter(({ used }) => used === undefined);
    }),
    ...unknownDomainLosses(config, slice ? SLICE_READ_DOMAINS : READ_DOMAINS),
    // A domain the reading replaces whole (one of another kind).
    ...readValueTables([...VALUE_TABLES, ...DRAWING_TABLES[mode]], config).filter(({ at, dropped }) => !dropped && !/[.[]/.test(at))
  ];
  return lossRepairs(losses, path, stored, { mode, root: 'config' });
};

const RETIRED_DOMAINS = RETIRED_STORED_FIELDS.filter((at) => !at.includes('.'));
/** The configuration domains Load reads. @type {ReadonlySet<string>} */
const READ_DOMAINS = new Set([...STORED_CONFIG_DOMAINS, ...RETIRED_DOMAINS]);
/**
 * The domains Load reads from a Session 46 slice: the registry's (a flat
 * draft's `modeProfiles` or `cliOptions` is no slice field).
 * @type {ReadonlySet<string>}
 */
const SLICE_READ_DOMAINS = new Set([...MODE_SLICE_CONTAINERS.config, ...RETIRED_DOMAINS]);
/**
 * A domain no gbdraw writer wrote is not read (TB19).
 * @param {Record<string, unknown>} config
 * @param {ReadonlySet<string>} [read]
 * @returns {SettingLoss[]}
 */
const unknownDomainLosses = (config, read = READ_DOMAINS) => Object.keys(config).filter((domain) => !read.has(domain))
  .map((domain) => ({ at: domain, stored: config[domain], used: undefined }));

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
 * The rows of what a normalizer noted it could not read (each loss at a place
 * the file holds, `at` under `path`; a key may hold a dot).
 * @param {readonly SettingLoss[]} losses
 * @param {string} path
 */
export const notedLossRepairs = (losses, path) => losses.flatMap((loss) => lossRepairs(
  [{ ...loss, at: '' }], [path, loss.at].filter(Boolean).join('.'), loss.stored
));

/**
 * The same for the LOSAT execution settings (`ui.losatExecution` of a Session
 * 46, `config.losat` of an older draft). `own`: the settings are a container
 * of their own (`ui.losatExecution`), so a container that is no object and a
 * key no execution setting names are not read either (A7e).
 * @param {unknown} execution
 * @param {string} path
 * @param {{ own?: boolean }} [options]
 */
export const losatExecutionRepairs = (execution, path, { own = false } = {}) => lossRepairs([
  ...(own && !isPlainObject(execution) && !isUnsetValue(execution) ? [{ at: '', stored: execution, used: undefined }] : []),
  ...(own && isPlainObject(execution) ? Object.keys(execution).filter((field) => !LOSAT_EXECUTION_ROWS.some((row) => row.field === field))
    .map((field) => ({ at: field, stored: execution[field], used: undefined })) : []),
  ...readValueTables([{ at: '', rows: LOSAT_EXECUTION_ROWS }], execution)
], path, execution);
