// @ts-check
import { parseOptionalCircularScalar } from '../track-slot-validation.js';

const isTypedMeasure = (value) => value !== null && typeof value === 'object' && !Array.isArray(value);

const requireEditorUnit = (unit) => {
  if (unit !== 'px' && unit !== 'factor') throw new Error('Circular measure editor unit must be px or factor.');
};

/** Read a detached view. The slot scalar remains the only persisted draft. */
export const readCircularMeasure = (scalar) => {
  const typed = isTypedMeasure(scalar);
  let valueText = String((typed ? scalar.value : scalar) ?? '');
  let selectedUnit = typed ? String(scalar.unit ?? '').trim().toLowerCase() : 'factor';
  let error = null;
  let isAuto = false;
  try {
    const pair = parseOptionalCircularScalar(scalar);
    isAuto = pair === null;
    selectedUnit = pair?.unit ?? null;
    if (isAuto) valueText = '';
    else if (!typed) {
      const text = valueText.trim();
      valueText = /%$/.test(text) ? String(pair.value)
        : pair.unit === 'px' ? text.slice(0, -2).trim() : text;
    }
  } catch (cause) {
    error = cause.message;
  }
  return { valueText, selectedUnit, isAuto, error };
};

/** Numeric edit and approved suffix adapter, returning one existing draft scalar. */
export const writeCircularMeasureValue = (valueText, selectedUnit) => {
  requireEditorUnit(selectedUnit);
  if (typeof valueText !== 'string') throw new Error('Circular measure editor value must be text.');
  const text = valueText.trim();
  if (!text) return null;
  // Plain text belongs to the selector. Only a complete valid suffix overrides it.
  if (/(?:px|%)$/i.test(text)) {
    const view = readCircularMeasure(text);
    if (!view.error) return { value: view.valueText, unit: view.selectedUnit };
  }
  return { value: text, unit: selectedUnit };
};

/** Reinterpret the same numeric lexeme; do not convert using a radius or Result. */
export const changeCircularMeasureUnit = (scalar, nextUnit) => {
  requireEditorUnit(nextUnit);
  const view = readCircularMeasure(scalar);
  if (view.isAuto) return null;
  if (view.selectedUnit === nextUnit) return scalar;
  return isTypedMeasure(scalar)
    ? { ...scalar, unit: nextUnit }
    : { value: typeof scalar === 'string' ? view.valueText : scalar, unit: nextUnit };
};

/** Request-to-draft projection shares the editor view; invalid objects stay invalid. */
export const projectCircularMeasureDraft = (scalar) => {
  if (scalar === null || scalar === undefined) return null;
  if (!isTypedMeasure(scalar)) return scalar;
  const view = readCircularMeasure(scalar);
  if (view.error) return scalar;
  return view.selectedUnit === 'px' ? `${view.valueText}px` : view.valueText;
};
