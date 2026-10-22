// @ts-check
import { parseOptionalCircularScalar } from './track-slot-validation.js';

/** @import { CircularMeasureObject } from './track-slot-validation.js' */

/**
 * The editor's view of a Width or Radius draft. An invalid draft's sentence
 * belongs to the validator (`circularMeasureInvalidParts`).
 * @typedef {object} CircularMeasureView
 * @property {string} valueText The numeric text the field shows.
 * @property {string | null} selectedUnit 'px' or 'factor'; null for Auto.
 * @property {boolean} isAuto
 * @property {boolean} invalid
 */

/**
 * @param {unknown} value
 * @returns {value is CircularMeasureObject} A non-array object: a typed `{ value, unit }` draft (its fields are checked on read).
 */
const isTypedMeasure = (value) => value !== null && typeof value === 'object' && !Array.isArray(value);

/** @type {(unit: unknown) => asserts unit is 'px' | 'factor'} */
const requireEditorUnit = (unit) => {
  if (unit !== 'px' && unit !== 'factor') throw new Error('Circular measure editor unit must be px or factor.');
};

/**
 * Read a detached view. The slot scalar remains the only persisted draft.
 * @param {unknown} scalar A slot Width or Radius draft.
 * @returns {CircularMeasureView}
 */
export const readCircularMeasure = (scalar) => {
  const typed = isTypedMeasure(scalar);
  let valueText = String((typed ? scalar.value : scalar) ?? '');
  /** @type {string | null} */
  let selectedUnit = typed ? String(scalar.unit ?? '').trim().toLowerCase() : 'factor';
  let invalid = false;
  let isAuto = false;
  try {
    const pair = parseOptionalCircularScalar(scalar);
    isAuto = pair === null;
    selectedUnit = pair?.unit ?? null;
    if (pair === null) valueText = '';
    else if (!typed) {
      const text = valueText.trim();
      valueText = /%$/.test(text) ? String(pair.value)
        : pair.unit === 'px' ? text.slice(0, -2).trim() : text;
    }
  } catch {
    invalid = true;
  }
  return { valueText, selectedUnit, isAuto, invalid };
};

/**
 * Numeric edit and approved suffix adapter, returning one existing draft scalar.
 * @param {unknown} valueText The edited text; anything else throws.
 * @param {unknown} selectedUnit The unit selector; anything but px or factor throws.
 */
export const writeCircularMeasureValue = (valueText, selectedUnit) => {
  requireEditorUnit(selectedUnit);
  if (typeof valueText !== 'string') throw new Error('Circular measure editor value must be text.');
  const text = valueText.trim();
  if (!text) return null;
  // Plain text belongs to the selector. Only a complete valid suffix overrides it.
  if (/(?:px|%)$/i.test(text)) {
    const view = readCircularMeasure(text);
    if (!view.invalid) return { value: view.valueText, unit: view.selectedUnit };
  }
  return { value: text, unit: selectedUnit };
};

/**
 * Reinterpret the same numeric lexeme; do not convert using a radius or Result.
 * @param {unknown} scalar A slot Width or Radius draft.
 * @param {unknown} nextUnit The unit selector; anything but px or factor throws.
 */
export const changeCircularMeasureUnit = (scalar, nextUnit) => {
  requireEditorUnit(nextUnit);
  const view = readCircularMeasure(scalar);
  if (view.isAuto) return null;
  if (view.selectedUnit === nextUnit) return scalar;
  return isTypedMeasure(scalar)
    ? { ...scalar, unit: nextUnit }
    : { value: typeof scalar === 'string' ? view.valueText : scalar, unit: nextUnit };
};

/**
 * Request-to-draft projection shares the editor view; invalid objects stay invalid.
 * @param {unknown} scalar A slot Width or Radius draft.
 */
export const projectCircularMeasureDraft = (scalar) => {
  if (scalar === null || scalar === undefined) return null;
  if (!isTypedMeasure(scalar)) return scalar;
  const view = readCircularMeasure(scalar);
  if (view.invalid) return scalar;
  return view.selectedUnit === 'px' ? `${view.valueText}px` : view.valueText;
};
