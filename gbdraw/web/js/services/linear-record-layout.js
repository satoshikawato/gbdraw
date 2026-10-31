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
