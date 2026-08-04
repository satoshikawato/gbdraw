const SVG_NUMBER = '[+-]?(?:\\d+(?:\\.\\d*)?|\\.\\d+)(?:[eE][+-]?\\d+)?';
const LEADING_TRANSLATE = new RegExp(
  `^(\\s*)translate\\(\\s*(${SVG_NUMBER})(?:\\s*,\\s*|\\s+)(${SVG_NUMBER})?\\s*\\)`
);
const ANY_TRANSLATE = new RegExp(
  `translate\\(\\s*(${SVG_NUMBER})(?:\\s*,\\s*|\\s+)(${SVG_NUMBER})?\\s*\\)`
);

const finiteNumber = (value, fallback = 0) => {
  const number = Number(value);
  return Number.isFinite(number) ? number : fallback;
};

export const readLeadingTranslate = (transformStr) => {
  const source = String(transformStr || '');
  const match = source.match(LEADING_TRANSLATE);
  if (!match) {
    return { found: false, x: 0, y: 0, prefix: '', tail: source };
  }
  return {
    found: true,
    x: finiteNumber(match[2]),
    y: finiteNumber(match[3]),
    prefix: match[1] || '',
    tail: source.slice(match[0].length)
  };
};

export const replaceLeadingTranslate = (transformStr, x, y) => {
  const source = String(transformStr || '');
  const leading = readLeadingTranslate(source);
  const translation = `translate(${finiteNumber(x)},${finiteNumber(y)})`;
  if (!leading.found) {
    return `${translation}${source.trim() ? ` ${source}` : ''}`;
  }
  return `${leading.prefix}${translation}${leading.tail}`;
};

export const prependTranslate = (transformStr, x, y) => {
  const source = String(transformStr || '');
  const translation = `translate(${finiteNumber(x)},${finiteNumber(y)})`;
  return `${translation}${source.trim() ? ` ${source}` : ''}`;
};

export const parseTransform = (transformStr) => {
  const source = String(transformStr || '');
  const leading = readLeadingTranslate(source);
  if (leading.found) return { x: leading.x, y: leading.y };
  const match = source.match(ANY_TRANSLATE);
  return match
    ? { x: finiteNumber(match[1]), y: finiteNumber(match[2]) }
    : { x: 0, y: 0 };
};
