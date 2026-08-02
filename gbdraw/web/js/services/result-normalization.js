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

