// @ts-check
const positiveRow = (value, fallback = 1) => {
  const row = Number(value);
  return Number.isInteger(row) && row > 0 ? row : fallback;
};

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
        && saved.canonicalColumn > 0
        ? { canonicalRow: row, canonicalColumn: saved.canonicalColumn }
        : {})
    };
  });
};

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

export const linearRecordLayoutHasSharedRow = (
  sequences,
  entries = [],
  { enabled = true } = {}
) => {
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
  const recordUids = (group) => group.records.map(({ sequence }) => String(sequence?.uid || ''));
  const sourceRows = groups.map((group) => (
    [...new Set(recordUids(group).map((uid) => rowByUid.get(uid)))].sort((left, right) => left - right)
  ));
  const rowSlots = sourceRows.flat().sort((left, right) => left - right);
  const slotIndex = new Map(rowSlots.map((row, slot) => [row, slot]));
  const consecutiveBlocks = slotIndex.size === rowSlots.length && sourceRows.every((rows) => (
    slotIndex.get(rows.at(-1)) - slotIndex.get(rows[0]) === rows.length - 1
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
      row: slots[blockRows.indexOf(rowByUid.get(uid))]
    }));
  });
  return { allowed: true, reason: '', rows };
};

export const linearRecordPositionTokens = (sequences, entries) => {
  const rows = new Map(reconcileLinearRecordLayout(sequences, entries).map((entry) => [entry.uid, entry.row]));
  return (Array.isArray(sequences) ? sequences : [])
    .map((sequence, index) => ({ index, row: rows.get(String(sequence?.uid || '')) || index + 1 }))
    .sort((left, right) => left.row - right.row || left.index - right.index)
    .map((entry) => `#${entry.index + 1}@${entry.row}`);
};

export const setLinearRecordRow = (entries, uid, row) => {
  const target = (Array.isArray(entries) ? entries : []).find((entry) => entry.uid === uid);
  if (target) target.row = positiveRow(row, target.row);
};

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
