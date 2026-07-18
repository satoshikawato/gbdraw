const positiveRow = (value, fallback = 1) => {
  const row = Number(value);
  return Number.isInteger(row) && row > 0 ? row : fallback;
};

export const reconcileLinearRecordLayout = (sequences, entries = []) => {
  const previous = new Map(
    (Array.isArray(entries) ? entries : []).map((entry) => [String(entry?.uid || ''), entry])
  );
  return (Array.isArray(sequences) ? sequences : []).map((sequence, index) => ({
    uid: String(sequence?.uid || ''),
    row: positiveRow(previous.get(String(sequence?.uid || ''))?.row, index + 1)
  }));
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
