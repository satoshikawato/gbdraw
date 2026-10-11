// @ts-check
const interactivePairName = (name) => {
  const normalized = String(name || '').trim();
  const match = normalized.match(/^(.*)\.interactive\.svg$/i);
  return match ? `${match[1]}.svg` : null;
};

export const normalizeLogicalResults = (results) => {
  if (!Array.isArray(results)) return [];
  const normalized = [];
  const indexByLogicalName = new Map();
  const sourceWasPlain = new Set();

  for (const result of results) {
    if (!result || typeof result !== 'object' || Array.isArray(result)) continue;
    const sourceName = String(result.name || '').trim();
    if (!sourceName) continue;
    const pairedPlainName = interactivePairName(sourceName);
    const logicalName = pairedPlainName || sourceName;
    const existingIndex = indexByLogicalName.get(logicalName);

    if (existingIndex === undefined) {
      indexByLogicalName.set(logicalName, normalized.length);
      normalized.push({
        ...result,
        name: pairedPlainName ? sourceName : logicalName
      });
      if (!pairedPlainName) sourceWasPlain.add(logicalName);
      continue;
    }

    if (!pairedPlainName && !sourceWasPlain.has(logicalName)) {
      normalized[existingIndex] = { ...result, name: logicalName };
      sourceWasPlain.add(logicalName);
    }
  }

  return normalized;
};

/**
 * Keep the names of the Results a rerender draws again (OV-136). The engine
 * names each Result after the request's output prefix (`out.svg`), while a
 * loaded Session's Result keeps the name its writer saved (the CLI saves the
 * output stem; a Gallery Session, its Gallery id). Every row of the Worker
 * reply that names a Result also carries its index, so all of them take the
 * kept name before admission checks each row against its Result. The caller
 * owns the reply. A reply with another number of Results keeps the engine's names.
 * @param {{ results?: unknown, metadata?: unknown }} response One Worker reply.
 * @param {readonly string[] | null} names The names of the Results drawn again.
 */
export const keepResultNames = (response, names) => {
  const results = response.results;
  if (
    !Array.isArray(names) || !Array.isArray(results) || names.length !== results.length
    || names.some((name) => typeof name !== 'string' || !name.trim())
  ) return;
  /** @type {Record<string, any>} */
  const metadata = response.metadata && typeof response.metadata === 'object' ? response.metadata : {};
  results.forEach((result, index) => {
    if (result && typeof result === 'object') result.name = names[index];
  });
  [
    metadata.featureCatalog?.items,
    metadata.legendRows,
    metadata.trackSlotGeometry?.records,
    metadata.annotationWarnings,
    metadata.comparisonWarnings
  ].forEach((rows) => {
    if (!Array.isArray(rows)) return;
    rows.forEach((row) => {
      if (row && typeof row === 'object' && Object.hasOwn(row, 'resultName')
        && Number.isInteger(row.resultIndex) && row.resultIndex >= 0 && row.resultIndex < names.length) {
        row.resultName = names[row.resultIndex];
      }
    });
  });
};

