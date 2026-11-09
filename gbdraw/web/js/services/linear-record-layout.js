// @ts-check
/**
 * One record's place in the Linear layout, as the layout writes it.
 * @typedef {{
 *   uid: string,
 *   row: number,
 *   canonicalCardinality?: 'exactly_one',
 *   canonicalRow?: number,
 *   canonicalColumn?: number
 * }} LinearRecordRow
 */
/** @typedef {{ uid?: unknown } | null | undefined} LayoutSequence A record slot; the layout reads only its uid. */
/** @typedef {{ records: { sequence?: LayoutSequence | null }[] }} LinearSourceGroup One File's records, as the Record Layout groups them. */
/** @typedef {Partial<Omit<LinearRecordRow, 'uid' | 'row'>> & { uid?: unknown, row?: unknown }} LayoutEntry A saved layout row. */
/** @typedef {{ recordKey?: unknown, x?: unknown, y?: unknown } | null | undefined} StoredRecordTranslation A translation as the page holds it, before the request validates it. */

/**
 * The stored translations of a set of records, in the records' order: each
 * record's stored entries as they are, or one at (0, 0) for a record that has
 * none. An entry of a record outside the set is left out, so a request carries
 * the translations of the records it draws after a File was removed, replaced,
 * or added, or a record turned off (OV-2291); the alignment owner takes the
 * base of a plan the same way. Nothing is converted: the request's validation
 * still refuses a duplicate key or a value that is not a finite number.
 * @param {readonly StoredRecordTranslation[]} stored
 * @param {readonly string[]} recordKeys
 * @returns {StoredRecordTranslation[]}
 */
export const recordTranslationsFor = (stored, recordKeys) => recordKeys.flatMap((recordKey) => {
  const entries = stored.filter((entry) => String(entry?.recordKey) === recordKey);
  return entries.length ? entries : [{ recordKey, x: 0, y: 0 }];
});

/**
 * A row field that the Web fills from run-time identities, and that a reader
 * empties when its list arranges nothing: every entry holds exactly
 * `entryKeys`, and its `positionKey` is its 1-based position in the list
 * (`normalize` of gbdraw/data/mode-slice-defaults.json).
 * @typedef {{
 *   domain: string, path: string, field: string,
 *   emptyWhenUnarranged: { entryKeys: string[], positionKey: string },
 *   reason: string
 * }} UnarrangedFieldRule
 */

/**
 * `value` with the field of each rule of `row` emptied when its list arranges
 * nothing. The one Web evaluation of the rule: the defaults generator
 * (tools/generate_mode_slice_defaults.mjs) and the Session 46-48 fast path
 * (session-drawings-upgrade.js) read it; its Python twin is
 * `gbdraw.session_drawings._normalized`, held equal by the generator's `--check`.
 * @param {readonly UnarrangedFieldRule[]} rules
 * @param {{ domain: string, path: string }} row
 * @param {unknown} value
 * @returns {unknown}
 */
export const withoutUnarrangedEntries = (rules, row, value) => rules
  .filter((rule) => rule.domain === row.domain && rule.path === row.path)
  .reduce((current, rule) => {
    if (!current || typeof current !== 'object' || Array.isArray(current)) return current;
    const holder = /** @type {Record<string, unknown>} */ (current);
    const entries = holder[rule.field];
    const { entryKeys, positionKey } = rule.emptyWhenUnarranged;
    const keys = [...entryKeys].sort().join('\0');
    const unarranged = Array.isArray(entries) && entries.every((entry, index) => (
      entry !== null && typeof entry === 'object' && !Array.isArray(entry)
      && Object.keys(entry).sort().join('\0') === keys
      && /** @type {Record<string, unknown>} */ (entry)[positionKey] === index + 1
    ));
    return unarranged ? { ...holder, [rule.field]: [] } : holder;
  }, value);

/**
 * @param {unknown} value
 * @param {number} [fallback]
 */
const positiveRow = (value, fallback = 1) => {
  const row = Number(value);
  return Number.isInteger(row) && row > 0 ? row : fallback;
};

/**
 * @param {readonly LayoutSequence[] | null | undefined} sequences
 * @param {LayoutEntry[] | null | undefined} [entries]
 * @returns {LinearRecordRow[]}
 */
export const reconcileLinearRecordLayout = (sequences, entries = []) => {
  const previous = new Map(
    (Array.isArray(entries) ? entries : []).map((entry) => [String(entry?.uid || ''), entry])
  );
  return (Array.isArray(sequences) ? sequences : []).map((sequence, index) => {
    const saved = previous.get(String(sequence?.uid || ''));
    const row = positiveRow(saved?.row, index + 1);
    return {
      uid: String(sequence?.uid || ''),
      row,
      ...(saved?.canonicalCardinality === 'exactly_one'
        ? { canonicalCardinality: 'exactly_one' } : {}),
      ...(saved?.canonicalRow === row && Number.isInteger(saved?.canonicalColumn)
        // `Number.isInteger` above rejects undefined.
        && /** @type {number} */ (saved.canonicalColumn) > 0
        ? { canonicalRow: row, canonicalColumn: saved.canonicalColumn }
        : {})
    };
  });
};

/**
 * @param {readonly LayoutSequence[] | null | undefined} sequences
 * @param {LayoutEntry[] | null | undefined} [entries]
 * @param {{ enabled?: boolean }} [options]
 * @returns {{ uid: string, row: number }[]}
 */
export const resolveEffectiveLinearRecordRows = (
  sequences,
  entries = [],
  { enabled = true } = {}
) => (
  enabled
    ? reconcileLinearRecordLayout(sequences, entries)
    : (Array.isArray(sequences) ? sequences : []).map((sequence, index) => ({
        uid: String(sequence?.uid || ''),
        row: index + 1
      }))
);

/**
 * @param {readonly LayoutSequence[] | null | undefined} sequences
 * @param {LayoutEntry[] | null | undefined} [entries]
 * @param {{ enabled?: boolean }} [options]
 */
export const linearRecordLayoutHasSharedRow = (
  sequences,
  entries = [],
  { enabled = true } = {}
) => {
  /** @type {Set<number>} */
  const seenRows = new Set();
  return resolveEffectiveLinearRecordRows(sequences, entries, { enabled }).some(({ row }) => {
    if (seenRows.has(row)) return true;
    seenRows.add(row);
    return false;
  });
};

// File moves exchange whole row blocks. Each File must use its own
// consecutive rows (rendered rows are compacted, so empty row numbers between
// blocks do not count); one row per File is the block size 1 case.
/**
 * @param {{
 *   sourceGroups: LinearSourceGroup[] | null | undefined,
 *   entries: LayoutEntry[] | null | undefined,
 *   sourceIndex: number,
 *   direction: number
 * }} move
 * @returns {{ allowed: boolean, reason: string, rows: { uid: string, row: number }[] }}
 */
export const planLinearSourceRowMove = ({
  sourceGroups,
  entries,
  sourceIndex,
  direction
}) => {
  const groups = Array.isArray(sourceGroups) ? sourceGroups : [];
  const sequences = groups.flatMap((group) => (
    Array.isArray(group?.records) ? group.records.map((record) => record.sequence) : []
  ));
  const layout = reconcileLinearRecordLayout(sequences, entries);
  const index = sourceIndex;
  const offset = direction;

  const rowByUid = new Map(layout.map((entry) => [entry.uid, entry.row]));
  /** @param {LinearSourceGroup} group */
  const recordUids = (group) => group.records.map(({ sequence }) => String(sequence?.uid || ''));
  const sourceRows = groups.map((group) => (
    // Every uid of a group is a uid of `layout`, which reconciles these same sequences.
    [...new Set(/** @type {number[]} */ (recordUids(group).map((uid) => rowByUid.get(uid))))].sort((left, right) => left - right)
  ));
  const rowSlots = sourceRows.flat().sort((left, right) => left - right);
  const slotIndex = new Map(rowSlots.map((row, slot) => [row, slot]));
  // A group without rows has no slot; the subtraction is then NaN and the comparison false.
  const consecutiveBlocks = slotIndex.size === rowSlots.length && sourceRows.every((rows) => (
    /** @type {number} */ (slotIndex.get(/** @type {number} */ (rows.at(-1)))) - /** @type {number} */ (slotIndex.get(rows[0])) === rows.length - 1
  ));
  if (!consecutiveBlocks) {
    return { allowed: false, reason: 'custom-layout', rows: layout };
  }

  const target = index + offset;
  if (!Number.isInteger(index) || ![-1, 1].includes(offset)
      || index < 0 || index >= groups.length || target < 0 || target >= groups.length) {
    return { allowed: false, reason: 'boundary', rows: layout };
  }
  const order = groups.map((_, groupIndex) => groupIndex);
  [order[index], order[target]] = [order[target], order[index]];
  let nextSlot = 0;
  const rows = order.flatMap((groupIndex) => {
    const blockRows = sourceRows[groupIndex];
    const slots = rowSlots.slice(nextSlot, nextSlot + blockRows.length);
    nextSlot += blockRows.length;
    return recordUids(groups[groupIndex]).map((uid) => ({
      uid,
      row: slots[blockRows.indexOf(/** @type {number} */ (rowByUid.get(uid)))]
    }));
  });
  return { allowed: true, reason: '', rows };
};

/**
 * @param {readonly LayoutSequence[] | null | undefined} sequences
 * @param {LayoutEntry[] | null | undefined} entries
 * @returns {string[]}
 */
export const linearRecordPositionTokens = (sequences, entries) => {
  const rows = new Map(reconcileLinearRecordLayout(sequences, entries).map((entry) => [entry.uid, entry.row]));
  return (Array.isArray(sequences) ? sequences : [])
    .map((sequence, index) => ({ index, row: rows.get(String(sequence?.uid || '')) || index + 1 }))
    .sort((left, right) => left.row - right.row || left.index - right.index)
    .map((entry) => `#${entry.index + 1}@${entry.row}`);
};

/**
 * @param {{ uid?: unknown, row: number }[] | null | undefined} entries
 * @param {string} uid
 * @param {unknown} row
 */
export const setLinearRecordRow = (entries, uid, row) => {
  const target = (Array.isArray(entries) ? entries : []).find((entry) => entry.uid === uid);
  if (target) target.row = positiveRow(row, target.row);
};

/**
 * @param {{ uid?: unknown }[] | null | undefined} sequences swapped in place with the record's row neighbour
 * @param {LayoutEntry[] | null | undefined} entries
 * @param {string} uid
 * @param {number} direction
 * @returns {LinearRecordRow[]}
 */
export const moveLinearRecordInRow = (sequences, entries, uid, direction) => {
  const layout = reconcileLinearRecordLayout(sequences, entries);
  const target = layout.find((entry) => entry.uid === uid);
  if (!target) return layout;
  const rowUids = layout.filter((entry) => entry.row === target.row).map((entry) => entry.uid);
  const current = rowUids.indexOf(uid);
  const next = current + (direction < 0 ? -1 : 1);
  if (current < 0 || next < 0 || next >= rowUids.length) return layout;
  const sequenceList = Array.isArray(sequences) ? sequences : [];
  const leftIndex = sequenceList.findIndex((sequence) => sequence.uid === rowUids[current]);
  const rightIndex = sequenceList.findIndex((sequence) => sequence.uid === rowUids[next]);
  if (leftIndex >= 0 && rightIndex >= 0) {
    [sequenceList[leftIndex], sequenceList[rightIndex]] = [sequenceList[rightIndex], sequenceList[leftIndex]];
  }
  return reconcileLinearRecordLayout(sequenceList, layout);
};

/**
 * A record slot as the Record Layout cells read it: its uid and its windows (a card drawn in several windows has one
 * cell per window, PD-OI-124).
 * @typedef {{ uid?: unknown, windows: readonly { column: number | null }[] }} LayoutCard
 */
/**
 * One cell of the Record Layout: a card drawn whole or in one window (`window` null), or one window of a card drawn
 * in several (its index in the card's list). `row` is the card's row, shared by its cells; `position` is the cell's
 * place in its row (from 1), the order the request writes and the drawing numbers its windows by; `placed`, whether
 * the cell has a stored place, `column`.
 * @typedef {{ uid: string, cardIndex: number, window: number | null, row: number, position: number, placed: boolean, column: number | null }} LinearLayoutCell
 */

/**
 * The Record Layout's cells in drawing order: by row, then by their stored place in the row (a window's `column`, a
 * whole card's `canonicalColumn`), then cards and windows in list order. A cell without a stored place follows those
 * with one, so cards that never moved keep the list order of the rows they share.
 * @param {readonly LayoutCard[]} cards
 * @param {LayoutEntry[] | null | undefined} entries
 * @returns {LinearLayoutCell[]}
 */
export const linearLayoutCells = (cards, entries) => {
  const layout = reconcileLinearRecordLayout(cards, entries);
  const placed = cards.flatMap((card, cardIndex) => {
    const { uid, row, canonicalColumn } = layout[cardIndex];
    return card.windows.length > 1
      ? card.windows.map((window, index) => ({ uid, cardIndex, window: /** @type {number | null} */ (index), row, column: window.column }))
      : [{ uid, cardIndex, window: /** @type {number | null} */ (null), row, column: canonicalColumn ?? null }];
  }).sort((left, right) => left.row - right.row
    || (left.column ?? Infinity) - (right.column ?? Infinity)
    || left.cardIndex - right.cardIndex
    || (left.window ?? -1) - (right.window ?? -1));
  /** @type {Map<number, number>} */
  const filled = new Map();
  return placed.map(({ uid, cardIndex, window, row, column }) => {
    const position = (filled.get(row) || 0) + 1;
    filled.set(row, position);
    return { uid, cardIndex, window, row, position, placed: column !== null, column };
  });
};

/**
 * The place the request writes for each cell of an arranged layout (`grid_column` of a whole card, `column` of a
 * window), keyed by `layoutCellKey`, when the stored places cannot be written as they are: only some cells have one
 * (Python reads a place for every cell or none), or two cells of a row hold the same. Null otherwise: the cells keep
 * the places they hold, or none.
 * @param {readonly LinearLayoutCell[]} cells
 * @returns {Map<string, number> | null}
 */
export const linearLayoutPlaces = (cells) => {
  const placed = cells.filter((cell) => cell.placed);
  const taken = new Set(placed.map((cell) => `${cell.row}:${cell.column}`));
  return placed.length === 0 || (placed.length === cells.length && taken.size === cells.length) ? null
    : new Map(cells.map((cell) => [layoutCellKey(cell.uid, cell.window), cell.position]));
};
/** @param {string} uid @param {number | null} window */
export const layoutCellKey = (uid, window) => `${uid}\u0000${window ?? ''}`;

/**
 * Writes the places of `rowCells`, one row's cells in their new order: each window's `column` and a whole card's
 * `canonicalRow`/`canonicalColumn` in its layout entry.
 * @param {readonly { windows: { column: number | null }[] }[]} cards
 * @param {{ uid: string, row: number, canonicalRow?: number, canonicalColumn?: number }[]} entries reconciled entries
 * @param {readonly LinearLayoutCell[]} rowCells
 */
const placeRowCells = (cards, entries, rowCells) => rowCells.forEach((cell, index) => {
  if (cell.window !== null) {
    cards[cell.cardIndex].windows[cell.window].column = index + 1;
    return;
  }
  const entry = entries.find((candidate) => candidate.uid === cell.uid);
  if (entry) Object.assign(entry, { canonicalRow: cell.row, canonicalColumn: index + 1 });
});

/**
 * Moves one cell a place up or down within its row (Up and Down of a drawing that has a card drawn in several
 * windows; PD-OI-124): every cell of that row gets its place. A whole card's cell also swaps the card with the card of
 * the cell it passes in `cards`, so their "#n" swap as Up and Down of a whole card always did (PD-OI-018); a window's
 * cell changes only places. Returns the reconciled entries, or null when the cell is at the end of its row.
 * @param {{ uid?: unknown, windows: { column: number | null }[] }[]} cards
 * @param {LayoutEntry[] | null | undefined} entries
 * @param {{ uid: string, window: number | null }} target
 * @param {number} direction
 * @returns {LinearRecordRow[] | null}
 */
export const moveLinearLayoutCell = (cards, entries, target, direction) => {
  const layout = reconcileLinearRecordLayout(cards, entries);
  const cells = linearLayoutCells(cards, layout);
  const at = cells.find((cell) => cell.uid === target.uid && cell.window === target.window);
  if (!at) return null;
  const rowCells = cells.filter((cell) => cell.row === at.row);
  const from = rowCells.indexOf(at);
  const to = from + (direction < 0 ? -1 : 1);
  if (to < 0 || to >= rowCells.length) return null;
  const passed = rowCells[to];
  [rowCells[from], rowCells[to]] = [rowCells[to], rowCells[from]];
  placeRowCells(cards, layout, rowCells);
  if (at.window !== null || passed.uid === at.uid) return layout;
  [cards[at.cardIndex], cards[passed.cardIndex]] = [cards[passed.cardIndex], cards[at.cardIndex]];
  return reconcileLinearRecordLayout(cards, layout);
};

/**
 * Puts the windows of card `uid`, a card drawn in several windows, at the end of its row in their order, after a Row
 * change moved the record there (PD-OI-124): the row's other cells keep their order and get their places. Returns the
 * reconciled entries, or null for a card drawn whole or in one window (its place goes with its row, as before).
 * @param {readonly { uid?: unknown, windows: { column: number | null }[] }[]} cards
 * @param {LayoutEntry[] | null | undefined} entries
 * @param {string} uid
 * @returns {LinearRecordRow[] | null}
 */
export const placeRecordAtRowEnd = (cards, entries, uid) => {
  const card = cards.find((candidate) => candidate.uid === uid);
  if (!card || card.windows.length < 2) return null;
  const layout = reconcileLinearRecordLayout(cards, entries);
  const cells = linearLayoutCells(cards, layout);
  const row = cells.find((cell) => cell.uid === uid)?.row;
  const rowCells = cells.filter((cell) => cell.row === row);
  placeRowCells(cards, layout, [...rowCells.filter((cell) => cell.uid !== uid), ...rowCells.filter((cell) => cell.uid === uid)]);
  return layout;
};

/**
 * The place a new window of card `uid` takes in its row: right after the card's last cell, the cells after it moving
 * one place on (PD-OI-124 "+ Add window"). Writes the places of that row's other cells and returns the new window's
 * column, `lastColumn`, the place of the card's last cell (the column of its window when the card held one), and the
 * reconciled entries.
 * @param {readonly { uid?: unknown, windows: { column: number | null }[] }[]} cards
 * @param {LayoutEntry[] | null | undefined} entries
 * @param {string} uid
 * @returns {{ column: number, lastColumn: number, entries: LinearRecordRow[] } | null}
 */
export const placeNewLayoutWindow = (cards, entries, uid) => {
  const layout = reconcileLinearRecordLayout(cards, entries);
  const cells = linearLayoutCells(cards, layout);
  const last = cells.filter((cell) => cell.uid === uid).at(-1);
  if (!last) return null;
  const rowCells = cells.filter((cell) => cell.row === last.row);
  placeRowCells(cards, layout, rowCells);
  // The cells after the new window move one place on.
  rowCells.slice(last.position).forEach((cell) => {
    if (cell.window !== null) cards[cell.cardIndex].windows[cell.window].column = cell.position + 1;
    else {
      const entry = layout.find((candidate) => candidate.uid === cell.uid);
      if (entry) entry.canonicalColumn = cell.position + 1;
    }
  });
  return { column: last.position + 1, lastColumn: last.position, entries: layout };
};
