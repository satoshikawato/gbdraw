const DECIMAL_NUMBER_PATTERN = /^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$/;

export const classifyOptionalPositiveNumber = (value) => {
  if (value === null || value === undefined) {
    return { status: 'auto', value: null };
  }

  if (typeof value === 'string') {
    const normalized = value.trim();
    if (!normalized || normalized.toLowerCase() === 'auto') {
      return { status: 'auto', value: null };
    }
    if (!DECIMAL_NUMBER_PATTERN.test(normalized)) {
      return { status: 'invalid', raw: value };
    }
    const numeric = Number(normalized);
    return Number.isFinite(numeric) && numeric > 0
      ? { status: 'valid', value: numeric }
      : { status: 'invalid', raw: value };
  }

  if (typeof value === 'number' && Number.isFinite(value) && value > 0) {
    return { status: 'valid', value };
  }

  return { status: 'invalid', raw: value };
};

