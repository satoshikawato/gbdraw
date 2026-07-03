const LENGTH_THRESHOLD_BP = 50000;

const DEFAULTS = Object.freeze({
  circular: {
    radiusPx: 390,
    trackRatio: 0.19,
    trackRatioFactors: {
      short: [0.50, 1.0, 1.0],
      long: [0.25, 1.0, 1.0]
    },
    trackDict: {
      short: {
        spreadout: { 1: 1.0, 2: 0.85, 3: 0.65, 4: 0.45 },
        middle: { 1: 1.0, 2: 0.75, 3: 0.55, 4: 0.35 },
        tuckin: { 1: 1.0, 2: 0.64, 3: 0.44, 4: 0.24 }
      },
      long: {
        spreadout: { 1: 1.0, 2: 0.80, 3: 0.60, 4: 0.40 },
        middle: { 1: 1.0, 2: 0.75, 3: 0.55, 4: 0.35 },
        tuckin: { 1: 1.0, 2: 0.70, 3: 0.50, 4: 0.30 }
      }
    },
    definitionFontSizePx: 18,
    plotTitleFontSizePx: 32,
    tickLabelFontSizePx: 14,
    axisStrokeWidthPx: { short: 3, long: 1 },
    labelFontSizePx: { short: 14, long: 8 }
  },
  linear: {
    defaultCdsHeightPx: { short: 80, long: 20 },
    definitionFontSizePx: { short: 24, long: 10 },
    labelFontSizePx: { short: 24, long: 5 },
    axisStrokeWidthPx: { short: 5, long: 2 },
    scaleFontSizePx: { short: 24, long: 16 },
    rulerLabelFontSizePx: { short: 20, long: 12 },
    defaultGcHeightPx: 20,
    depthHeightPx: 10,
    comparisonHeightPx: 60
  },
  feature: {
    blockStrokeWidthPx: { short: 2, long: 0 },
    lineStrokeWidthPx: { short: 5, long: 1 }
  },
  legend: {
    boxSizePx: { short: 24, long: 20 },
    fontSizePx: { short: 20, long: 16 }
  },
  labels: {
    spacingPx: 3,
    radiusOffset: 1
  },
  scale: {
    strokeWidthPx: 3
  }
});

const SLIDING_WINDOW_PRESETS = Object.freeze([
  { maxExclusive: 1_000_000, window: 1000, step: 100 },
  { maxExclusive: 10_000_000, window: 10_000, step: 1000 },
  { maxExclusive: Infinity, window: 100_000, step: 10_000 }
]);

const LINEAR_TICK_THRESHOLDS = Object.freeze([
  [2000, 100],
  [20_000, 1000],
  [50_000, 5000],
  [150_000, 10_000],
  [250_000, 50_000],
  [1_000_000, 100_000],
  [2_000_000, 200_000],
  [5_000_000, 500_000],
  [Infinity, 1_000_000]
]);

const normalizeOptionalText = (value) => {
  const text = String(value ?? '').trim();
  return text.length > 0 ? text : null;
};

const parsePositiveNumber = (value) => {
  const text = normalizeOptionalText(value);
  if (text === null) return null;
  const numeric = Number(String(text).replace(/px$/i, ''));
  return Number.isFinite(numeric) && numeric > 0 ? numeric : null;
};

const roundDisplayNumber = (value, digits = 1) => {
  const numeric = Number(value);
  if (!Number.isFinite(numeric)) return '';
  const fixed = numeric.toFixed(digits);
  return fixed.replace(/\.0+$/, '').replace(/(\.\d*?)0+$/, '$1');
};

const formatBp = (value) => {
  const numeric = Number(value);
  if (!Number.isFinite(numeric)) return '';
  if (Math.abs(numeric) >= 1_000_000 && numeric % 1_000_000 === 0) {
    return `${roundDisplayNumber(numeric / 1_000_000, 1)} Mbp`;
  }
  if (Math.abs(numeric) >= 1000 && numeric % 1000 === 0) {
    return `${roundDisplayNumber(numeric / 1000, 1)} kbp`;
  }
  return `${roundDisplayNumber(numeric, 0)} bp`;
};

const formatPx = (value) => {
  const text = roundDisplayNumber(value, 1);
  return text ? `${text} px` : '';
};

const formatFactor = (value) => {
  const text = roundDisplayNumber(value, 2);
  return text ? `${text} R` : '';
};

const withAuto = (text) => (String(text || '').trim() ? `${text} (auto)` : '');

const formatShortLong = (values, formatter, activeLength = null) => {
  if (!values || typeof values !== 'object') return '';
  if (activeLength === 'short' || activeLength === 'long') {
    return withAuto(formatter(values[activeLength]));
  }
  const shortText = formatter(values.short);
  const longText = formatter(values.long);
  if (!shortText || !longText) return '';
  if (shortText === longText) return withAuto(shortText);
  const unitMatch = shortText.match(/^(.+?)\s+([A-Za-z%]+)$/);
  const longUnitMatch = longText.match(/^(.+?)\s+([A-Za-z%]+)$/);
  if (unitMatch && longUnitMatch && unitMatch[2] === longUnitMatch[2]) {
    return `${unitMatch[1]}/${longUnitMatch[1]} ${unitMatch[2]} (s/l auto)`;
  }
  return `${shortText}/${longText} (s/l auto)`;
};

const normalizeLengthParam = (value) => {
  const text = String(value || '').trim().toLowerCase();
  return text === 'short' || text === 'long' ? text : null;
};

const circularRecordLengths = (state) => {
  const recordsRef = state?.circularRecordList;
  const records = Array.isArray(recordsRef?.value)
    ? recordsRef.value
    : (Array.isArray(recordsRef) ? recordsRef : []);
  return records
    .map((entry) => Number(entry?.record_length ?? entry?.length ?? 0))
    .filter((value) => Number.isFinite(value) && value > 0);
};

const activeCircularLengthParam = (state) => {
  const lengths = circularRecordLengths(state);
  if (lengths.length === 0) return null;
  return Math.max(...lengths) < LENGTH_THRESHOLD_BP ? 'short' : 'long';
};

const activeLengthParam = (state) => {
  const mode = String(state?.mode?.value || state?.mode || 'circular');
  if (mode === 'circular') return activeCircularLengthParam(state);
  return null;
};

const maxKnownRecordLength = (state) => {
  const mode = String(state?.mode?.value || state?.mode || 'circular');
  if (mode === 'circular') {
    const lengths = circularRecordLengths(state);
    return lengths.length ? Math.max(...lengths) : null;
  }
  return null;
};

const slidingWindowForLength = (length) => {
  const numeric = Number(length);
  if (!Number.isFinite(numeric) || numeric <= 0) return null;
  return SLIDING_WINDOW_PRESETS.find((entry) => numeric < entry.maxExclusive) || SLIDING_WINDOW_PRESETS[2];
};

const autoLinearTickInterval = (length) => {
  const numeric = Number(length);
  if (!Number.isFinite(numeric) || numeric <= 0) return null;
  const match = LINEAR_TICK_THRESHOLDS.find(([threshold]) => numeric < threshold);
  return match ? match[1] : 1_000_000;
};

const autoCircularTickInterval = (length) => {
  const numeric = Number(length);
  if (!Number.isFinite(numeric) || numeric <= 0) return null;
  if (numeric <= 30_000) return 1000;
  if (numeric <= 50_000) return 5000;
  if (numeric <= 150_000) return 10_000;
  if (numeric <= 1_000_000) return 50_000;
  if (numeric <= 10_000_000) return 500_000;
  return 1_000_000;
};

const normalizeTrackPreset = (value) => {
  const text = String(value || '').trim().toLowerCase();
  return ['tuckin', 'middle', 'spreadout'].includes(text) ? text : 'tuckin';
};

const circularBaseTrackWidthPx = (lengthParam) => {
  const base = DEFAULTS.circular.radiusPx * DEFAULTS.circular.trackRatio;
  const factors = DEFAULTS.circular.trackRatioFactors[normalizeLengthParam(lengthParam) || 'long'];
  return { base, factors };
};

const circularWidthForRenderer = (renderer, lengthParam) => {
  const { base, factors } = circularBaseTrackWidthPx(lengthParam);
  if (renderer === 'features' || renderer === 'sequence_conservation') return base * Number(factors[0]);
  if (renderer === 'depth') return base * Number(factors[1]) * 0.5;
  return base * Number(factors[1]);
};

const previewSpacingPx = () => Math.max(1.0, 0.01 * DEFAULTS.circular.radiusPx);

const circularLengthValues = (state, resolver) => {
  const active = activeCircularLengthParam(state);
  if (active) return resolver(active);
  return {
    short: resolver('short'),
    long: resolver('long')
  };
};

const circularWidthText = (state, renderer) => {
  const value = circularLengthValues(state, (lengthParam) => circularWidthForRenderer(renderer, lengthParam));
  return typeof value === 'number'
    ? withAuto(formatPx(value))
    : formatShortLong(value, formatPx, null);
};

const circularFeatureRadius = (state, lengthParam) => {
  const preset = normalizeTrackPreset(state?.form?.track_type);
  const laneWidth = circularWidthForRenderer('features', lengthParam);
  const laneCount = state?.form?.separate_strands ? 2 : 1;
  const spacing = previewSpacingPx();
  const bandWidth = (laneCount * laneWidth) + (Math.max(0, laneCount - 1) * spacing);
  if (preset === 'tuckin') {
    return (DEFAULTS.circular.radiusPx - spacing - (bandWidth / 2)) / DEFAULTS.circular.radiusPx;
  }
  if (preset === 'spreadout') {
    return (DEFAULTS.circular.radiusPx + spacing + (bandWidth / 2)) / DEFAULTS.circular.radiusPx;
  }
  return 1.0;
};

const circularBuiltinTrackId = (renderer, state) => {
  const showDepth = Boolean(state?.form?.show_depth);
  const showGc = !Boolean(state?.form?.suppress_gc);
  if (renderer === 'depth') return showDepth ? 2 : 2;
  if (renderer === 'dinucleotide_content') return showDepth ? 3 : 2;
  if (renderer === 'dinucleotide_skew') {
    if (showDepth) return showGc ? 4 : 3;
    return showGc ? 3 : 2;
  }
  return null;
};

const circularRadiusForRenderer = (state, renderer, lengthParam) => {
  if (renderer === 'features') return circularFeatureRadius(state, lengthParam);
  const trackId = circularBuiltinTrackId(renderer, state);
  if (trackId === null) return null;
  const preset = normalizeTrackPreset(state?.form?.track_type);
  return DEFAULTS.circular.trackDict[normalizeLengthParam(lengthParam) || 'long']?.[preset]?.[trackId] ?? null;
};

const circularRadiusText = (state, renderer) => {
  const value = circularLengthValues(state, (lengthParam) => circularRadiusForRenderer(state, renderer, lengthParam));
  return typeof value === 'number'
    ? withAuto(formatFactor(value))
    : formatShortLong(value, formatFactor, null);
};

const depthTickFontSize = (state, trackConfig = null) => {
  const mode = String(state?.mode?.value || state?.mode || 'circular');
  if (mode === 'linear') {
    const trackHeight = parsePositiveNumber(trackConfig?.height)
      ?? parsePositiveNumber(state?.adv?.depth_height)
      ?? DEFAULTS.linear.depthHeightPx;
    return Math.max(5, Math.min(8, trackHeight * 0.7));
  }
  const trackWidth = parsePositiveNumber(state?.adv?.depth_width_circular)
    ?? circularWidthForRenderer('depth', activeCircularLengthParam(state) || 'long');
  return Math.max(5, Math.min(8, trackWidth * 0.22));
};

const gcTickFontSize = (state) => {
  const mode = String(state?.mode?.value || state?.mode || 'circular');
  if (mode === 'linear') {
    const trackHeight = parsePositiveNumber(state?.adv?.gc_height) ?? DEFAULTS.linear.defaultGcHeightPx;
    return Math.max(5, Math.min(8, trackHeight * 0.7));
  }
  const trackWidth = parsePositiveNumber(state?.adv?.gc_content_width_circular)
    ?? circularWidthForRenderer('dinucleotide_content', activeCircularLengthParam(state) || 'long');
  return Math.max(5, Math.min(8, trackWidth * 0.22));
};

const lengthDependentPx = (state, values) => formatShortLong(values, formatPx, activeLengthParam(state));

const dinucleotideWindowText = (state, field) => {
  const manualWindow = parsePositiveNumber(state?.adv?.window_size);
  const manualStep = parsePositiveNumber(state?.adv?.step_size);
  if (field === 'window' && manualWindow !== null) return withAuto(formatBp(manualWindow));
  if (field === 'step' && manualStep !== null) return withAuto(formatBp(manualStep));

  const length = maxKnownRecordLength(state);
  const preset = slidingWindowForLength(length);
  if (preset) return withAuto(formatBp(preset[field]));
  const values = SLIDING_WINDOW_PRESETS.map((entry) => formatBp(entry[field]));
  return `${values.join('/')} (auto)`;
};

const depthWindowText = (state, field) => {
  const sourceField = field === 'window' ? 'window_size' : 'step_size';
  const manualSource = parsePositiveNumber(state?.adv?.[sourceField]);
  if (manualSource !== null) {
    const divisor = field === 'window' ? 10 : 10;
    const minValue = field === 'window' ? 100 : 1;
    return withAuto(formatBp(Math.max(minValue, manualSource / divisor)));
  }

  const length = maxKnownRecordLength(state);
  const preset = slidingWindowForLength(length);
  if (preset) {
    const source = field === 'window' ? preset.window : preset.step;
    const minValue = field === 'window' ? 100 : 1;
    return withAuto(formatBp(Math.max(minValue, source / 10)));
  }
  return field === 'window'
    ? '100/1k/10k bp (auto)'
    : '10/100/1k bp (auto)';
};

const linearAxisGapText = (state) => {
  const layout = String(state?.form?.linear_track_layout || 'middle').toLowerCase();
  if (layout === 'middle') return 'ignored';
  const manualHeight = parsePositiveNumber(state?.adv?.feature_height);
  const multiplier = state?.form?.separate_strands ? 0.25 : 0.30;
  if (manualHeight !== null) return withAuto(formatPx(manualHeight * multiplier));
  return formatShortLong(
    {
      short: DEFAULTS.linear.defaultCdsHeightPx.short * multiplier,
      long: DEFAULTS.linear.defaultCdsHeightPx.long * multiplier
    },
    formatPx,
    null
  );
};

const linearLabelPlacementText = (state) => {
  const layout = String(state?.form?.linear_track_layout || 'middle').toLowerCase();
  return layout === 'below' || layout === 'tuckin' ? 'below' : 'above';
};

const scaleIntervalText = (state) => {
  const length = maxKnownRecordLength(state);
  const mode = String(state?.mode?.value || state?.mode || 'circular');
  if (mode === 'circular') {
    const interval = autoCircularTickInterval(length);
    return interval ? withAuto(formatBp(interval)) : 'by record length (auto)';
  }
  const interval = autoLinearTickInterval(length);
  return interval ? withAuto(formatBp(interval)) : 'by record length (auto)';
};

const plotTitleText = (state) => {
  if (String(state?.mode?.value || state?.mode || 'circular') !== 'circular') return '';
  if (String(state?.adv?.plot_title_position || 'none') === 'none') return '';
  const pieces = [state?.form?.species, state?.form?.strain]
    .map((value) => String(value || '').trim())
    .filter(Boolean);
  return pieces.length ? `${pieces.join(' ')} (auto)` : 'species + strain (auto)';
};

const autoTextByKey = (state, key, context = null) => {
  switch (key) {
    case 'depthWindow': return depthWindowText(state, 'window');
    case 'depthStep': return depthWindowText(state, 'step');
    case 'depthMin': return '0 (auto)';
    case 'depthMax': return 'from depth data (auto)';
    case 'circularDepthWidth': return circularWidthText(state, 'depth');
    case 'linearDepthHeight': return withAuto(formatPx(DEFAULTS.linear.depthHeightPx));
    case 'linearDepthTrackHeight': return withAuto(formatPx(parsePositiveNumber(state?.adv?.depth_height) ?? DEFAULTS.linear.depthHeightPx));
    case 'depthLargeTick': return 'endpoints only';
    case 'depthTickFont': return withAuto(formatPx(depthTickFontSize(state, context)));
    case 'conservationReference': return state?.circularConservation?.source === 'losat' ? 'subject (auto)' : 'detect side (auto)';
    case 'conservationRingWidth': return circularWidthText(state, 'sequence_conservation');
    case 'conservationRingGap': return withAuto(formatPx(previewSpacingPx()));
    case 'plotTitle': return plotTitleText(state);
    case 'definitionFontSize':
      return String(state?.mode?.value || state?.mode || 'circular') === 'circular'
        ? withAuto(formatPx(DEFAULTS.circular.definitionFontSizePx))
        : lengthDependentPx(state, DEFAULTS.linear.definitionFontSizePx);
    case 'legendBoxSize': return lengthDependentPx(state, DEFAULTS.legend.boxSizePx);
    case 'legendFontSize': return lengthDependentPx(state, DEFAULTS.legend.fontSizePx);
    case 'linearFeatureHeight': return lengthDependentPx(state, DEFAULTS.linear.defaultCdsHeightPx);
    case 'linearAxisGap': return linearAxisGapText(state);
    case 'circularFeatureWidth': return circularWidthText(state, 'features');
    case 'blockStrokeWidth': return lengthDependentPx(state, DEFAULTS.feature.blockStrokeWidthPx);
    case 'lineStrokeWidth': return lengthDependentPx(state, DEFAULTS.feature.lineStrokeWidthPx);
    case 'labelFontSize':
      return String(state?.mode?.value || state?.mode || 'circular') === 'circular'
        ? lengthDependentPx(state, DEFAULTS.circular.labelFontSizePx)
        : lengthDependentPx(state, DEFAULTS.linear.labelFontSizePx);
    case 'labelRendering': return 'per label';
    case 'linearLabelPlacement': return linearLabelPlacementText(state);
    case 'circularLabelSpacing': return withAuto(formatPx(DEFAULTS.labels.spacingPx));
    case 'linearLabelSpacing': return withAuto(formatPx(DEFAULTS.labels.spacingPx));
    case 'labelRadiusOffset': return withAuto(roundDisplayNumber(DEFAULTS.labels.radiusOffset, 1));
    case 'axisStrokeWidth':
      return String(state?.mode?.value || state?.mode || 'circular') === 'circular'
        ? lengthDependentPx(state, DEFAULTS.circular.axisStrokeWidthPx)
        : lengthDependentPx(state, DEFAULTS.linear.axisStrokeWidthPx);
    case 'scaleStrokeWidth': return withAuto(formatPx(DEFAULTS.scale.strokeWidthPx));
    case 'scaleFontSize':
      return String(state?.form?.scale_style || 'bar') === 'ruler'
        ? lengthDependentPx(state, DEFAULTS.linear.rulerLabelFontSizePx)
        : lengthDependentPx(state, DEFAULTS.linear.scaleFontSizePx);
    case 'scaleInterval': return scaleIntervalText(state);
    case 'tickLabelFontSize': return withAuto(formatPx(DEFAULTS.circular.tickLabelFontSizePx));
    case 'dinucleotideWindow': return dinucleotideWindowText(state, 'window');
    case 'dinucleotideStep': return dinucleotideWindowText(state, 'step');
    case 'gcTickFont': return withAuto(formatPx(gcTickFontSize(state)));
    case 'linearGcHeight': return withAuto(formatPx(DEFAULTS.linear.defaultGcHeightPx));
    case 'circularGcWidth': return circularWidthText(state, 'dinucleotide_content');
    case 'circularGcRadius': return circularRadiusText(state, 'dinucleotide_content');
    case 'circularSkewWidth': return circularWidthText(state, 'dinucleotide_skew');
    case 'circularSkewRadius': return circularRadiusText(state, 'dinucleotide_skew');
    case 'pairwiseMatchHeight': return withAuto(formatPx(DEFAULTS.linear.comparisonHeightPx));
    default: return '';
  }
};

export const createAutoValueDisplay = (state) => {
  const autoValueText = (key, context = null) => autoTextByKey(state, key, context);
  const autoValueVisible = (value, key, context = null) => (
    normalizeOptionalText(value) === null && autoValueText(key, context) !== ''
  );

  return {
    autoValueText,
    autoValueVisible
  };
};
