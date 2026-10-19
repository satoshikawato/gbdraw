// @ts-check
import { normalizeDefaultColor, normalizeSpecificRuleColor } from '../utils/color-utils.js';
import { diagnosticError } from '../utils/error-normalization.js';
import { normalizeTsvCell } from '../utils/tsv-cell.js';

/**
 * A row of a pre-normalization table (Sessions 31-39) that was repaired or dropped.
 * @typedef {{ row: number, repair: 'joined' | 'dropped' | 'invalid' }} TableRepair
 */
/** @typedef {{ feat: string, qual: string, val: string, color: string, cap: string, fromFile: true }} ParsedSpecificRule */
/**
 * What the serializers read of a row; a row that is not an object writes nothing.
 * @typedef {{ feat?: unknown, qual?: unknown, val?: unknown, color?: unknown, cap?: unknown, key?: unknown, order?: unknown }} EditorTableRow
 */

const SPECIFIC_RULE_COLUMNS = Object.freeze(['feature_type', 'qualifier', 'pattern', 'color']);

// The request fields of each editor table: the Web writes the `File` form, and
// the CLI and the Python API write the canonical `Table` form. Python reads the
// `Table` form first.
export const EDITOR_TABLE_OPTIONS = Object.freeze({
  visibility: ['featureVisibilityTable', 'featureVisibilityTableFile'],
  whitelist: ['labelWhitelistTable', 'labelWhitelistFile'],
  priority: ['qualifierPriorityTable', 'qualifierPriorityFile'],
  labelOverrides: ['labelOverrideTable', 'labelOverrideFile']
});
/**
 * The resource reference of one editor table in request `diagramOptions`.
 * @param {Record<string, any> | null | undefined} options
 * @param {keyof typeof EDITOR_TABLE_OPTIONS} table
 * @returns {{ resourceId: string, representation?: string } | null}
 */
export const editorTableRef = (options, table) => EDITOR_TABLE_OPTIONS[table]
  .map((key) => options?.[key]).find((ref) => ref?.resourceId) || null;

// Python reads these tables with a fixed column count and rejects any other row
// (gbdraw/io/table_text.py::read_literal_table), so the import does too (R4). The
// diagnostic is the one the Label override parser uses.
/**
 * @param {string[]} parts
 * @param {number} columnCount
 * @param {number} row
 */
const requireColumnCount = (parts, columnCount, row) => {
  if (parts.length !== columnCount) {
    throw diagnosticError('TABLE_INVALID', { row, columnCount, reason: 'FIELDS' });
  }
};

// Sessions 31–39 stored these tables as their Web writer wrote them, before cell
// values were normalized (normalizeTsvCell): a tab in a value made extra cells and
// a line break a short row. For them (`repairs`), a row is read as the current
// writer writes it: the extra cells join the last field with one space, and a row
// without its required fields (`complete`) is dropped; `repairs` lists each one.
/**
 * @param {string} line
 * @param {number} columnCount
 * @param {number} row
 * @param {TableRepair[] | null} repairs
 * @param {(cells: string[]) => unknown} complete
 * @returns {string[] | null}
 */
const tableCells = (line, columnCount, row, repairs, complete) => {
  const parts = line.split('\t');
  if (!repairs) {
    requireColumnCount(parts, columnCount, row);
    return parts.map((part) => part.trim());
  }
  const cells = [...parts.slice(0, columnCount - 1), parts.slice(columnCount - 1).join('\t')]
    .map(normalizeTsvCell);
  if (parts.length < columnCount || !complete(cells)) {
    repairs.push({ row, repair: 'dropped' });
    return null;
  }
  if (parts.length > columnCount) repairs.push({ row, repair: 'joined' });
  return cells;
};
/**
 * @param {string} key
 * @param {string} color
 */
const isColorHeader = (key, color) => key.toLowerCase() === 'feature_type' && color.toLowerCase() === 'color';

// A color outside the Default colors domain stops the import at its line, as
// Python stops (OV-272); a blank color cell keeps the palette color, as Python
// keeps the built-in one.
/**
 * @param {string} text
 * @param {{ legacyRows?: boolean }} [options]
 */
export const parseColorTable = (text, { legacyRows = false } = {}) => {
  /** @type {Record<string, string>} */
  const colors = {};
  let count = 0;
  /** @type {TableRepair[] | null} */
  const repairs = legacyRows ? [] : null;
  const lines = text.split(/\r?\n/);

  for (const [index, line] of lines.entries()) {
    if (!line.trim() || line.trim().startsWith('#') || line.trim().startsWith('[')) continue;
    const cells = tableCells(line, 2, index + 1, repairs, ([key, color]) => key && color);
    if (!cells) continue;
    const [key, color] = cells;
    if (isColorHeader(key, color) || !key || !color) continue;
    const normalized = normalizeDefaultColor(color);
    if (!normalized) {
      // A Sessions 31-39 table (`repairs`) lists the row instead of failing Load.
      if (repairs) { repairs.push({ row: index + 1, repair: 'invalid' }); continue; }
      throw diagnosticError('TABLE_INVALID', { row: index + 1, field: 'color', reason: 'COLOR' });
    }
    colors[key] = normalized;
    count++;
  }

  return { colors, count, repairs: repairs || [] };
};

/** @param {string | null | undefined} text */
export const parseSpecificRules = (text) => {
  /** @type {ParsedSpecificRule[]} */
  const rules = [];
  /** @type {ParsedSpecificRule[]} */
  const rulesWithCaptions = [];
  /** @type {Set<string>} */
  const identities = new Set();
  const lines = String(text ?? '').split(/\r?\n/);

  for (let index = 0; index < lines.length; index += 1) {
    const line = lines[index];
    const lineNo = index + 1;
    if (!line.trim() || line.trim().startsWith('#')) continue;
    const parts = line.split('\t');
    if (
      parts.length >= 2 &&
      parts[0].trim().toLowerCase() === 'feature_type' &&
      parts[1].trim().toLowerCase() === 'qualifier_key'
    ) continue;

    if (parts.length < 4 || parts.length > 5) {
      throw diagnosticError('TABLE_INVALID', { row: lineNo, reason: 'SPECIFIC_COLUMNS' });
    }

    const [feat, qual, val, colorRaw, captionRaw = ''] = parts.map((part) => part.trim());
    const required = [feat, qual, val, colorRaw];
    const missingIndex = required.findIndex((value) => !value);
    if (missingIndex >= 0) {
      throw diagnosticError('TABLE_INVALID', { row: lineNo, field: SPECIFIC_RULE_COLUMNS[missingIndex], reason: 'REQUIRED' });
    }
    const color = normalizeSpecificRuleColor(colorRaw);
    if (!color) {
      throw diagnosticError('TABLE_INVALID', { row: lineNo, field: 'color', reason: 'COLOR' });
    }

    /** @type {ParsedSpecificRule} */
    const rule = {
      feat,
      qual,
      val,
      color,
      cap: captionRaw,
      fromFile: true
    };
    const identity = JSON.stringify([rule.feat, rule.qual, rule.val, rule.color, rule.cap]);
    if (identities.has(identity)) continue;
    identities.add(identity);
    rules.push(rule);
    if (rule.cap) rulesWithCaptions.push(rule);
  }

  return { rules, rulesWithCaptions, count: rules.length };
};

/** @param {ReadonlyArray<EditorTableRow | null | undefined> | null | undefined} rules */
export const serializeSpecificRules = (rules) => {
  const rows = (Array.isArray(rules) ? rules : [])
    .map((rule) => [
      normalizeTsvCell(rule?.feat),
      normalizeTsvCell(rule?.qual),
      normalizeTsvCell(rule?.val),
      normalizeTsvCell(rule?.color),
      normalizeTsvCell(rule?.cap)
    ])
    .filter((fields) => fields.slice(0, 4).every(Boolean))
    .map((fields) => fields.join('\t'));

  return rows.length > 0 ? `${rows.join('\n')}\n` : '';
};

/** @param {ReadonlyArray<EditorTableRow | null | undefined> | null | undefined} rules */
export const serializeLabelWhitelistRules = (rules) => {
  const rows = (Array.isArray(rules) ? rules : [])
    .map((rule) => [rule?.feat, rule?.qual, rule?.key].map(normalizeTsvCell))
    .filter(([feat, qual]) => feat && qual)
    .map((fields) => fields.join('\t'));

  return rows.length > 0 ? `${rows.join('\n')}\n` : '';
};

/** @param {ReadonlyArray<EditorTableRow | null | undefined> | null | undefined} rules */
export const serializeQualifierPriorityRules = (rules) => {
  const rows = (Array.isArray(rules) ? rules : [])
    .map((rule) => [rule?.feat, rule?.order].map(normalizeTsvCell))
    .filter((fields) => fields.every(Boolean))
    .map((fields) => fields.join('\t'));

  return rows.length > 0 ? `${rows.join('\n')}\n` : '';
};

/**
 * @param {string} text
 * @param {{ legacyRows?: boolean }} [options]
 */
export const parsePriorityRules = (text, { legacyRows = false } = {}) => {
  /** @type {{ feat: string, order: string }[]} */
  const rules = [];
  /** @type {TableRepair[] | null} */
  const repairs = legacyRows ? [] : null;
  const lines = text.split(/\r?\n/);

  for (const [index, line] of lines.entries()) {
    if (!line.trim() || line.trim().startsWith('#')) continue;
    const cells = tableCells(line, 2, index + 1, repairs, ([feat, order]) => feat && order);
    if (!cells) continue;
    const [feat, order] = cells;
    if (feat.toLowerCase() === 'feature_type' && order.toLowerCase() === 'priorities') continue;
    rules.push({ feat, order });
  }

  return { rules, count: rules.length, repairs: repairs || [] };
};

/**
 * @param {string} text
 * @param {{ legacyRows?: boolean }} [options]
 */
export const parseWhitelistRules = (text, { legacyRows = false } = {}) => {
  /** @type {{ feat: string, qual: string, key: string }[]} */
  const rules = [];
  /** @type {TableRepair[] | null} */
  const repairs = legacyRows ? [] : null;
  const lines = text.split(/\r?\n/);

  for (const [index, line] of lines.entries()) {
    if (!line.trim() || line.trim().startsWith('#')) continue;
    const cells = tableCells(line, 3, index + 1, repairs, ([feat, qual]) => feat && qual);
    if (!cells) continue;
    const [feat, qual, key] = cells;
    // The CLI and Python store the table with this header.
    if (feat.toLowerCase() === 'feature_type' && qual.toLowerCase() === 'qualifier'
      && String(key).toLowerCase() === 'keyword') continue;
    rules.push({ feat, qual, key });
  }

  return { rules, count: rules.length, repairs: repairs || [] };
};

/** @param {string} text */
export const parseBlacklistWords = (text) => {
  const words = text
    .split(/[\r\n,]+/)
    .map((word) => word.trim())
    .filter((word) => word);

  return { words, count: words.length };
};

// OV-1280 (D-19): what a file card's file put into a drawing, so Remove on the
// card takes it back and leaves the drawing as if the file had never been
// loaded. A card is named by its `files` slot. A row a file adds to a list
// carries its card's mark (`cardRowMark`): a Specific rules row `fromFile`
// (on main since 2026-01-14; D-ORIGIN-1), a Whitelist or Qualifier priority row
// `importedFrom: <card>`. A priority row whose order the file replaced keeps
// the order it had as `handOrder`. The drawing's `fileImports` keeps what a
// file set in the Default colors (each type's color and the user color it
// replaced, `null` for the palette's) and in the Blacklist text (its words and
// the text before). An edit of a row makes it the user's (the mark goes), and
// so does a later change of a value the file set: Remove keeps both. A
// Whitelist or priority row without a mark (Sessions before 48, hand rows) is
// the user's.
/** @typedef {'d_color' | 't_color' | 'blacklist' | 'whitelist' | 'qualifier_priority'} FileImportCard */
/** @typedef {'t_color' | 'whitelist' | 'qualifier_priority'} RowFileCard A card whose file adds rows. */
/** @typedef {{ fromFile?: boolean, importedFrom?: string }} CardMarkedRow */
/** @typedef {{ color: string, handColor: string | null }} ImportedDefaultColor */
/** @typedef {{ words: string[], before: string }} ImportedBlacklist */
/**
 * What a drawing's Default colors and Blacklist files set (`fileImports`).
 * @typedef {{ d_color?: Record<string, ImportedDefaultColor>, blacklist?: ImportedBlacklist }} FileImports
 */
/** @typedef {{ feat: string, order: string, importedFrom?: 'qualifier_priority', handOrder?: string }} PriorityRuleRow */
/** @typedef {{ feat: string, qual: string, key: string, importedFrom?: 'whitelist' }} WhitelistRuleRow */

/**
 * The mark a row of each card's file carries.
 * @typedef {{ t_color: { fromFile: true }, whitelist: { importedFrom: 'whitelist' }, qualifier_priority: { importedFrom: 'qualifier_priority' } }} CardRowMarks
 */
/** @type {Readonly<CardRowMarks>} */
const CARD_ROW_MARKS = Object.freeze({
  t_color: Object.freeze({ fromFile: true }),
  whitelist: Object.freeze({ importedFrom: 'whitelist' }),
  qualifier_priority: Object.freeze({ importedFrom: 'qualifier_priority' })
});

/**
 * Whether `row` is a row `card`'s file added.
 * @param {Readonly<CardMarkedRow>} row
 * @param {RowFileCard} card
 */
export const isCardRow = (row, card) => card === 't_color'
  ? row.fromFile === true
  : row.importedFrom === CARD_ROW_MARKS[card].importedFrom;

/**
 * `rows` marked as rows of `card`'s file.
 * @template {object} T
 * @template {RowFileCard} C
 * @param {readonly T[]} rows
 * @param {C} card
 * @returns {(T & CardRowMarks[C])[]}
 */
export const markedRows = (rows, card) => rows.map((row) => ({ ...row, ...CARD_ROW_MARKS[card] }));

/**
 * Makes `row` the user's: it loses `card`'s mark (an edit of the row).
 * @param {CardMarkedRow} row Changed in place.
 * @param {RowFileCard} card
 */
export const unmarkCardRow = (row, card) => {
  Object.keys(CARD_ROW_MARKS[card]).forEach((key) => { delete row[/** @type {keyof CardMarkedRow} */ (key)]; });
};

/** @param {unknown} left @param {unknown} right */
const sameColor = (left, right) => String(left ?? '').toLowerCase() === String(right ?? '').toLowerCase();

/**
 * The Default colors with a file's colors, and what the file set.
 * @param {Readonly<Record<string, string>>} current
 * @param {Readonly<Record<string, string>>} colors The file's colors by feature type.
 * @param {(type: string) => string | null} userColorOf The user color a type has now (`null`: the palette's).
 * @returns {{ colors: Record<string, string>, imported: Record<string, ImportedDefaultColor> }}
 */
export const withImportedDefaultColors = (current, colors, userColorOf) => ({
  colors: { ...current, ...colors },
  imported: Object.fromEntries(Object.entries(colors)
    .map(([type, color]) => [type, { color, handColor: userColorOf(type) }]))
});

/**
 * The Default colors without what a file set: each type that still has the
 * file's color takes back its user color, else the palette's (none: no entry).
 * @param {Readonly<Record<string, string>>} current
 * @param {Readonly<Record<string, ImportedDefaultColor>>} imported
 * @param {Readonly<Record<string, string>>} paletteColors
 * @returns {Record<string, string>}
 */
export const withoutImportedDefaultColors = (current, imported, paletteColors) => {
  const colors = { ...current };
  Object.entries(imported).forEach(([type, { color, handColor }]) => {
    if (!sameColor(colors[type], color)) return;
    const restored = handColor ?? paletteColors[type];
    if (restored === undefined) delete colors[type];
    else colors[type] = restored;
  });
  return colors;
};

/**
 * The Blacklist text after a file's words: the words follow the text.
 * @param {string} before
 * @param {readonly string[]} words
 */
const blacklistWithWords = (before, words) => {
  const existing = before ? before.trim() : '';
  const separator = existing && !existing.endsWith(',') ? ', ' : '';
  return existing + separator + words.join(', ');
};

/**
 * @param {string} text
 * @param {readonly string[]} words The file's words (none: no change).
 * @returns {{ text: string, imported: ImportedBlacklist | null }}
 */
export const withImportedBlacklistWords = (text, words) => (words.length > 0
  ? { text: blacklistWithWords(text, words), imported: { words: [...words], before: text } }
  : { text, imported: null });

/**
 * The Blacklist text without a file's words: the text before when it is
 * unchanged since, else the text with each of the words (its last occurrence
 * and the separator before it) taken out.
 * @param {string} text
 * @param {Readonly<ImportedBlacklist>} imported
 */
export const withoutImportedBlacklistWords = (text, { words, before }) => {
  if (text === blacklistWithWords(before, words)) return before;
  const parts = text.split(/([,\n])/);
  [...words].reverse().forEach((word) => {
    let at = -1;
    for (let index = parts.length - 1; index >= 0 && at < 0; index -= 2) {
      if (parts[index].trim() === word) at = index;
    }
    if (at > 0) parts.splice(at - 1, 2);
    else if (at === 0 && parts.length > 1) {
      parts.splice(0, 2);
      parts[0] = parts[0].replace(/^\s+/, '');
    } else if (at === 0) parts[0] = '';
  });
  return parts.join('');
};

/**
 * Adds a Qualifier priority file's rules: a rule replaces the order of its
 * feature type's row, which keeps the user order it had (`handOrder`).
 * @param {PriorityRuleRow[]} rows The drawing's rows, changed in place.
 * @param {readonly { feat: string, order: string }[]} rules
 */
export const importPriorityRows = (rows, rules) => {
  rules.forEach(({ feat, order }) => {
    const row = rows.find((entry) => entry.feat === feat);
    if (!row) {
      rows.push(...markedRows([{ feat, order }], 'qualifier_priority'));
      return;
    }
    if (!isCardRow(row, 'qualifier_priority')) Object.assign(row, { handOrder: row.order }, CARD_ROW_MARKS.qualifier_priority);
    row.order = order;
  });
};

/**
 * Takes back a Qualifier priority file's rules: a row it added goes, a row
 * whose order it replaced takes back its user order.
 * @param {PriorityRuleRow[]} rows Changed in place.
 */
export const removeImportedPriorityRows = (rows) => {
  for (let index = rows.length - 1; index >= 0; index -= 1) {
    const { feat, handOrder } = rows[index];
    if (!isCardRow(rows[index], 'qualifier_priority')) continue;
    if (handOrder === undefined) rows.splice(index, 1);
    else rows.splice(index, 1, { feat, order: handOrder });
  }
};

/**
 * The rows of a list without those `card`'s file added.
 * @template {CardMarkedRow} T
 * @param {readonly T[]} rows
 * @param {RowFileCard} card
 * @returns {T[]}
 */
export const withoutImportedRows = (rows, card) => rows.filter((row) => !isCardRow(row, card));
