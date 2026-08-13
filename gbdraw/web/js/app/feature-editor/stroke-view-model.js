const HEX_COLOR_PATTERN = /^#(?:[0-9a-f]{3}|[0-9a-f]{4}|[0-9a-f]{6}|[0-9a-f]{8})$/i;

const text = (value) => String(value ?? '').trim();
const hasOwn = (value, key) => Object.prototype.hasOwnProperty.call(value || {}, key);

export const normalizeFeatureStrokeColor = (value) => {
  const normalized = text(value).toLowerCase();
  return HEX_COLOR_PATTERN.test(normalized) ? normalized : null;
};

export const normalizeFeatureStrokeWidth = (value) => {
  if (value === null || value === undefined || value === '') return null;
  const normalized = Number(value);
  return Number.isFinite(normalized) && normalized >= 0 ? normalized : null;
};

export const normalizeFeatureStrokeValue = (value) => {
  const kind = value && typeof value === 'object'
    ? text(value.kind).toLowerCase()
    : '';
  if (value === null || kind === 'inherit') return Object.freeze({ kind: 'inherit' });

  const rawColor = value && typeof value === 'object'
    ? (value.strokeColor ?? value.color)
    : value;
  const rawWidth = value && typeof value === 'object'
    ? (value.strokeWidth ?? value.width)
    : undefined;
  const strokeColor = rawColor === undefined ? null : normalizeFeatureStrokeColor(rawColor);
  const strokeWidth = normalizeFeatureStrokeWidth(rawWidth);
  if ((rawColor !== undefined && !strokeColor) || (rawWidth !== undefined && strokeWidth === null)) {
    return null;
  }
  if (!strokeColor && strokeWidth === null) return null;
  return Object.freeze({
    kind: 'stroke',
    ...(strokeColor ? { strokeColor } : {}),
    ...(strokeWidth !== null ? { strokeWidth } : {})
  });
};

export const normalizeStrokeOverride = (value) => {
  if (!value || typeof value !== 'object' || Array.isArray(value)) return Object.freeze({});
  const strokeColor = hasOwn(value, 'strokeColor')
    ? normalizeFeatureStrokeColor(value.strokeColor)
    : null;
  const strokeWidth = hasOwn(value, 'strokeWidth')
    ? normalizeFeatureStrokeWidth(value.strokeWidth)
    : null;
  return Object.freeze({
    ...(strokeColor ? { strokeColor } : {}),
    ...(strokeWidth !== null ? { strokeWidth } : {})
  });
};

const ownedColor = (override) => {
  const normalized = normalizeStrokeOverride(override);
  return hasOwn(normalized, 'strokeColor') ? normalized.strokeColor : null;
};

const ownedWidth = (override) => {
  const normalized = normalizeStrokeOverride(override);
  return hasOwn(normalized, 'strokeWidth') ? normalized.strokeWidth : null;
};

const propertyOrigin = ({ exact, caption, fallback, property }) => {
  const read = property === 'strokeColor' ? ownedColor : ownedWidth;
  const exactValue = read(exact);
  if (exactValue !== null) return { value: exactValue, origin: 'local' };
  const captionValue = read(caption);
  if (captionValue !== null) return { value: captionValue, origin: 'legend-caption' };
  return { value: fallback, origin: 'svg-default' };
};

const originLabel = (origin, caption) => ({
  local: 'Feature override',
  'legend-caption': caption
    ? `Inherited from legend: ${caption}`
    : 'Inherited from legend caption',
  'svg-default': 'Inherited from SVG default'
}[origin]);

/**
 * Resolve the Feature stroke without turning inheritance into a disabled mode.
 * Color and width fall through independently; the primary origin describes the
 * color shown by the picker.
 */
export const resolveFeatureStrokeViewModel = ({
  exactOverride = null,
  captionOverride = null,
  legendCaption = '',
  svgDefaultStroke = null,
  svgDefaultColor = '#000000',
  svgDefaultWidth = 0.5
} = {}) => {
  const fallbackColor = normalizeFeatureStrokeColor(
    svgDefaultStroke?.color ?? svgDefaultStroke?.strokeColor ?? svgDefaultColor
  ) || '#000000';
  const fallbackWidth = normalizeFeatureStrokeWidth(
    svgDefaultStroke?.width ?? svgDefaultStroke?.strokeWidth ?? svgDefaultWidth
  ) ?? 0.5;
  const color = propertyOrigin({
    exact: exactOverride,
    caption: captionOverride,
    fallback: fallbackColor,
    property: 'strokeColor'
  });
  const width = propertyOrigin({
    exact: exactOverride,
    caption: captionOverride,
    fallback: fallbackWidth,
    property: 'strokeWidth'
  });
  const normalizedExact = normalizeStrokeOverride(exactOverride);

  return Object.freeze({
    effectiveColor: color.value,
    effectiveWidth: width.value,
    explicitValue: hasOwn(normalizedExact, 'strokeColor')
      ? normalizedExact.strokeColor
      : null,
    explicitWidth: hasOwn(normalizedExact, 'strokeWidth')
      ? normalizedExact.strokeWidth
      : null,
    origin: color.origin,
    widthOrigin: width.origin,
    originLabel: originLabel(color.origin, text(legendCaption)),
    canReset: Object.keys(normalizedExact).length > 0,
    allowNone: false,
    pickerEnabled: true
  });
};

export const getFeatureStrokeViewModel = resolveFeatureStrokeViewModel;
