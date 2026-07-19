let comparisonCounter = 0;

const keyFor = (queryUid, subjectUid) => `${queryUid}->${subjectUid}`;

export const hasLinearComparisonIntent = ({
  layoutEnabled = false,
  comparisons = [],
  sequences = [],
  blastSource = 'upload'
} = {}) => {
  const records = Array.isArray(sequences) ? sequences : [];
  if (layoutEnabled) return Array.isArray(comparisons) && comparisons.length > 0;
  if (String(blastSource || '') === 'losat') return records.length > 1;
  return records.slice(0, -1).some((sequence) => Boolean(sequence?.blast));
};

export const reconcileLinearComparisons = (sequences, comparisons = []) => {
  const valid = new Set((Array.isArray(sequences) ? sequences : []).map((sequence) => String(sequence?.uid || '')));
  const seen = new Set();
  return (Array.isArray(comparisons) ? comparisons : []).filter((comparison) => {
    const queryUid = String(comparison?.queryUid || '');
    const subjectUid = String(comparison?.subjectUid || '');
    const key = keyFor(queryUid, subjectUid);
    if (!valid.has(queryUid) || !valid.has(subjectUid) || queryUid === subjectUid || seen.has(key)) return false;
    seen.add(key);
    return true;
  });
};

export const createLinearComparison = (queryUid, subjectUid, source = 'upload') => ({
  id: `linear-comparison-${Date.now()}-${++comparisonCounter}`,
  queryUid: String(queryUid || ''),
  subjectUid: String(subjectUid || ''),
  source,
  file: null
});

export const addLinearComparison = (comparisons, queryUid, subjectUid, source = 'upload') => {
  if (!queryUid || !subjectUid || queryUid === subjectUid) return false;
  const key = keyFor(queryUid, subjectUid);
  if (comparisons.some((item) => keyFor(item.queryUid, item.subjectUid) === key)) return false;
  comparisons.push(createLinearComparison(queryUid, subjectUid, source));
  return true;
};

export const adjacentRowPairs = (sequences, layout, allPairs = false) => {
  const rows = new Map((Array.isArray(layout) ? layout : []).map((entry) => [entry.uid, Number(entry.row)]));
  const grouped = new Map();
  (Array.isArray(sequences) ? sequences : []).forEach((sequence, index) => {
    const row = rows.get(sequence.uid) || index + 1;
    if (!grouped.has(row)) grouped.set(row, []);
    grouped.get(row).push(sequence.uid);
  });
  const orderedRows = [...grouped.keys()].sort((a, b) => a - b);
  const pairs = [];
  for (let index = 0; index < orderedRows.length - 1; index += 1) {
    const upper = grouped.get(orderedRows[index]);
    const lower = grouped.get(orderedRows[index + 1]);
    if (allPairs) {
      upper.forEach((queryUid) => lower.forEach((subjectUid) => pairs.push([queryUid, subjectUid])));
    } else {
      for (let pairIndex = 0; pairIndex < Math.min(upper.length, lower.length); pairIndex += 1) {
        pairs.push([upper[pairIndex], lower[pairIndex]]);
      }
    }
  }
  return pairs;
};

export const validateLinearComparisons = (sequences, layout, comparisons) => {
  const rows = new Map((Array.isArray(layout) ? layout : []).map((entry) => [entry.uid, Number(entry.row)]));
  const validUids = new Set((Array.isArray(sequences) ? sequences : []).map((sequence) => sequence.uid));
  for (const comparison of comparisons) {
    if (!validUids.has(comparison.queryUid) || !validUids.has(comparison.subjectUid)) return 'A comparison references a removed record.';
    if (comparison.queryUid === comparison.subjectUid) return 'A comparison cannot connect a record to itself.';
    const queryRow = rows.get(comparison.queryUid);
    const subjectRow = rows.get(comparison.subjectUid);
    if (queryRow === subjectRow || Math.abs(queryRow - subjectRow) !== 1) return 'Comparisons must connect different adjacent rows.';
    if (comparison.source === 'upload' && !comparison.file) return 'Every uploaded comparison requires a BLAST TSV file.';
  }
  return '';
};
