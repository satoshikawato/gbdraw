export const DEFINITION_LINE_STYLE_KINDS = Object.freeze([
  'name',
  'subtitle',
  'replicon',
  'accession',
  'length'
]);

const normalizeDefinitionLineFontSize = (value) => {
  if (value === null || value === undefined) return null;
  const raw = String(value).trim();
  if (!raw || ['auto', 'none', 'null', 'default'].includes(raw.toLowerCase())) return null;
  const numeric = Number(raw);
  return Number.isFinite(numeric) && numeric > 0 ? numeric : null;
};

const normalizeDefinitionLineFontWeight = (value) => {
  if (value === null || value === undefined) return null;
  const normalized = String(value).trim().toLowerCase();
  if (!normalized || ['auto', 'none', 'null', 'default', 'normal'].includes(normalized)) return null;
  if (['bold', 'lighter', 'bolder'].includes(normalized)) return normalized;
  return /^(100|200|300|400|500|600|700|800|900)$/.test(normalized) ? normalized : null;
};

const normalizeDefinitionLineFill = (value) => {
  if (value === null || value === undefined) return null;
  const normalized = String(value).trim();
  return normalized || null;
};

export const createDefaultDefinitionLineStyle = () => ({
  font_size: null,
  font_weight: null,
  fill: null
});

export const createDefaultLinearDefinitionLineStyles = () => Object.fromEntries(
  DEFINITION_LINE_STYLE_KINDS.map((kind) => [kind, createDefaultDefinitionLineStyle()])
);

export const normalizeDefinitionLineStyleState = (source) => {
  const normalized = createDefaultLinearDefinitionLineStyles();
  if (!source || typeof source !== 'object' || Array.isArray(source)) return normalized;
  DEFINITION_LINE_STYLE_KINDS.forEach((kind) => {
    const entry = source[kind];
    if (!entry || typeof entry !== 'object' || Array.isArray(entry)) return;
    normalized[kind] = {
      font_size: normalizeDefinitionLineFontSize(entry.font_size),
      font_weight: normalizeDefinitionLineFontWeight(entry.font_weight),
      fill: normalizeDefinitionLineFill(entry.fill)
    };
  });
  return normalized;
};
