// @ts-check
/** @import { FieldLeaf, FieldReading } from '../utils/field-reading.js' */
export const DEFINITION_LINE_STYLE_KINDS = Object.freeze([
  'name',
  'subtitle',
  'replicon',
  'accession',
  'length'
]);

/** @typedef {{ font_size: number | null, font_weight: string | null, fill: string | null }} DefinitionLineStyle */
/** @typedef {typeof DEFINITION_LINE_STYLE_KINDS[number]} DefinitionLineStyleKind */

// Each reader returns the value used and whether reading it lost the stored
// one; unset, `auto` and the other names of the default lose nothing.
/** @param {unknown} value @returns {FieldReading<number | null>} */
const normalizeDefinitionLineFontSize = (value) => {
  if (value === null || value === undefined) return { value: null, lost: false };
  const raw = String(value).trim();
  if (!raw || ['auto', 'none', 'null', 'default'].includes(raw.toLowerCase())) return { value: null, lost: false };
  const numeric = Number(raw);
  return Number.isFinite(numeric) && numeric > 0 ? { value: numeric, lost: false } : { value: null, lost: true };
};

/** @param {unknown} value @returns {FieldReading<string | null>} */
const normalizeDefinitionLineFontWeight = (value) => {
  if (value === null || value === undefined) return { value: null, lost: false };
  const normalized = String(value).trim().toLowerCase();
  if (!normalized || ['auto', 'none', 'null', 'default', 'normal'].includes(normalized)) return { value: null, lost: false };
  return ['bold', 'lighter', 'bolder'].includes(normalized) || /^(100|200|300|400|500|600|700|800|900)$/.test(normalized)
    ? { value: normalized, lost: false } : { value: null, lost: true };
};

/** @param {unknown} value @returns {FieldReading<string | null>} */
const normalizeDefinitionLineFill = (value) => {
  if (value === null || value === undefined) return { value: null, lost: false };
  return { value: String(value).trim() || null, lost: false };
};

const STYLE_READERS = Object.freeze({
  font_size: normalizeDefinitionLineFontSize, font_weight: normalizeDefinitionLineFontWeight, fill: normalizeDefinitionLineFill
});

/** @returns {DefinitionLineStyle} */
export const createDefaultDefinitionLineStyle = () => ({
  font_size: null,
  font_weight: null,
  fill: null
});

/** @returns {Record<DefinitionLineStyleKind, DefinitionLineStyle>} */
export const createDefaultLinearDefinitionLineStyles = () => Object.fromEntries(
  DEFINITION_LINE_STYLE_KINDS.map((kind) => [kind, createDefaultDefinitionLineStyle()])
);

/**
 * The Linear definition line styles a drawing stores; a kind that is not an
 * object takes the default style.
 * @param {unknown} source
 * @returns {FieldReading<Record<DefinitionLineStyleKind, DefinitionLineStyle>>}
 */
export const readDefinitionLineStyles = (source) => {
  const value = createDefaultLinearDefinitionLineStyles();
  if (!source || typeof source !== 'object' || Array.isArray(source)) {
    return { value, lost: source !== null && source !== undefined };
  }
  const stored = /** @type {Record<string, unknown>} */ (source);
  /** @type {FieldLeaf[]} */
  const leaves = [];
  DEFINITION_LINE_STYLE_KINDS.forEach((kind) => {
    const entry = stored[kind];
    if (entry === null || entry === undefined || entry === '') return;
    if (typeof entry !== 'object' || Array.isArray(entry)) {
      leaves.push({ at: `.${kind}`, stored: entry, used: undefined });
      return;
    }
    const fields = /** @type {Record<string, unknown>} */ (entry);
    const style = createDefaultDefinitionLineStyle();
    /** @type {(keyof DefinitionLineStyle)[]} */ (['font_size', 'font_weight', 'fill']).forEach((field) => {
      const reading = STYLE_READERS[field](fields[field]);
      if (reading.lost) leaves.push({ at: `.${kind}.${field}`, stored: fields[field], used: reading.value });
      Object.assign(style, { [field]: reading.value });
    });
    value[kind] = style;
  });
  return { value, lost: leaves.length > 0, ...(leaves.length ? { leaves } : {}) };
};

/**
 * @param {unknown} source
 * @returns {Record<DefinitionLineStyleKind, DefinitionLineStyle>}
 */
export const normalizeDefinitionLineStyleState = (source) => readDefinitionLineStyles(source).value;
