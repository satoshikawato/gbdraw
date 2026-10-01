import { diagnosticError } from '../services/error-normalization.js';

export const DECIMAL_NUMBER_PATTERN = /^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$/;

/**
 * JSON-representable classification of an optional number: blank or Auto is
 * 'auto', a finite number (any sign) is 'valid', anything else is 'invalid'.
 */
export const classifyOptionalNumber = (value) => {
  if (value === null || value === undefined) return { status: 'auto', value: null };
  if (typeof value === 'string') {
    const normalized = value.trim();
    if (!normalized || normalized.toLowerCase() === 'auto') return { status: 'auto', value: null };
    const numeric = DECIMAL_NUMBER_PATTERN.test(normalized) ? Number(normalized) : NaN;
    return Number.isFinite(numeric) ? { status: 'valid', value: numeric } : { status: 'invalid', raw: value };
  }
  return typeof value === 'number' && Number.isFinite(value)
    ? { status: 'valid', value }
    : { status: 'invalid', raw: value };
};

export const classifyOptionalPositiveNumber = (value) => {
  const classified = classifyOptionalNumber(value);
  return classified.status === 'valid' && !(classified.value > 0)
    ? { status: 'invalid', raw: value }
    : classified;
};

/**
 * R7 request projection: blank -> null, number -> the same number, anything
 * else -> typed INPUT_INVALID. It never converts or range-checks a value; the
 * Python typed layer owns the domain. `context` names the field or configPath.
 */
export const projectOptionalNumber = (value, context = {}) => {
  const classified = classifyOptionalNumber(value);
  if (classified.status === 'invalid') throw diagnosticError('INPUT_INVALID', { ...context, reason: 'FINITE' });
  return classified.value;
};
