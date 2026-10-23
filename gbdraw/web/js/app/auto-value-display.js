// @ts-check
/** @import { DrawingState } from '../state.js' */
import {
  linearLabelSizeLengths,
  normalizeOptionalText,
  sharedSizeClass,
  sizeClassForRecordLengths
} from '../services/track-slot-display.js';
import { requireCurrentLinearTrackLayout } from '../services/current-option-values.js';
import { t } from '../utils/i18n.js';
import { MODE_PROFILE_DATA } from '../mode-profiles.generated.js';

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

// The length tiers of the Auto sliding window and scale ticks are Python's
// (gbdraw/auto_sizes.py), carried by the generated mode-profile data.
const { slidingWindow, circularTicks, linearTicks } = MODE_PROFILE_DATA.autoSizes;
// One object per tier: `sharedLengthValue` compares the values of the drawn records by identity.
const SLIDING_WINDOW_PRESETS = Object.freeze([
  ...slidingWindow.tiers.map(([maxExclusive, window, step]) => Object.freeze({ maxExclusive, window, step })),
  Object.freeze({ maxExclusive: Infinity, window: slidingWindow.above[0], step: slidingWindow.above[1] })
]);

/** @param {unknown} value */
const parsePositiveNumber = (value) => {
  const text = normalizeOptionalText(value);
  if (text === null) return null;
  const numeric = Number(String(text).replace(/px$/i, ''));
  return Number.isFinite(numeric) && numeric > 0 ? numeric : null;
};

/**
 * @param {unknown} value
 * @param {number} [digits]
 */
const roundDisplayNumber = (value, digits = 1) => {
  const numeric = Number(value);
  if (!Number.isFinite(numeric)) return '';
  const fixed = numeric.toFixed(digits);
  return fixed.replace(/\.0+$/, '').replace(/(\.\d*?)0+$/, '$1');
};

/** @param {unknown} value */
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

/** @param {unknown} value */
const formatPx = (value) => {
  const text = roundDisplayNumber(value, 1);
  return text ? `${text} px` : '';
};

/** @param {unknown} value */
const formatFactor = (value) => {
  const text = roundDisplayNumber(value, 2);
  return text ? `${text} R` : '';
};

/** @param {unknown} text */
const withAuto = (text) => (String(text || '').trim() ? t('{value} (auto)', { value: String(text) }) : '');

// An Auto value that depends on genome size shows the value for the
// current genome; before its length is known, it says so instead of listing
// the short- and long-genome values side by side.
/**
 * @typedef {'short' | 'long'} SizeClass
 */

/**
 * One value per genome size class.
 * @template T
 * @typedef {{ short: T, long: T }} ByLength
 */

/**
 * @param {ByLength<number | null> | null | undefined} values
 * @param {(value: number | null) => string} formatter
 * @param {string | null} [activeLength]
 */
const formatShortLong = (values, formatter, activeLength = null) => {
  if (!values || typeof values !== 'object') return '';
  if (activeLength === 'short' || activeLength === 'long') {
    return withAuto(formatter(values[activeLength]));
  }
  const shortText = formatter(values.short);
  const longText = formatter(values.long);
  if (!shortText || !longText) return '';
  if (shortText === longText) return withAuto(shortText);
  return t('Auto (by genome size)');
};

/**
 * @param {unknown} value
 * @returns {SizeClass | null}
 */
const normalizeLengthParam = (value) => {
  const text = String(value || '').trim().toLowerCase();
  return text === 'short' || text === 'long' ? text : null;
};

/**
 * What the record-length readers take from the Web state (`state.js`): the
 * shown mode, a Vue ref.
 * @typedef {object} RecordLengthState
 * @property {{ value?: string }} [mode]
 */

/**
 * The track whose Auto text is asked for: only a Linear Depth track's height is read.
 * @typedef {{ height?: unknown }} AutoValueContext
 */

/**
 * The records the next Generate draws, by length. The notes do not choose
 * records: the app registers these readers in `createAutoValueDisplay`, and
 * each returns an empty list while a drawn length is not known.
 * @typedef {object} DrawnLengthSources
 * @property {() => number[][]} circularDiagramLengths the records of each Circular diagram
 * @property {() => number[]} linearRecordLengths the records of the Linear diagram, in drawing order
 */
/** @type {WeakMap<RecordLengthState, DrawnLengthSources>} */
const drawnLengthSources = new WeakMap();

/** @param {RecordLengthState} state */
const shownMode = (state) => (
  String(state?.mode?.value || state?.mode || 'circular') === 'circular' ? 'circular' : 'linear'
);

/** @param {RecordLengthState} state */
const circularDiagramLengths = (state) => drawnLengthSources.get(state)?.circularDiagramLengths() || [];

/** @param {RecordLengthState} state */
const linearRecordLengths = (state) => drawnLengthSources.get(state)?.linearRecordLengths() || [];

/**
 * The record lengths of each diagram the shown mode's next Generate draws, or
 * [] while unknown.
 * @param {RecordLengthState} state
 */
const drawnDiagramLengths = (state) => {
  if (shownMode(state) === 'circular') return circularDiagramLengths(state);
  const lengths = linearRecordLengths(state);
  return lengths.length ? [lengths] : [];
};

/** @param {RecordLengthState} state @param {DrawingState} drawing */
const activeCircularLengthParam = (state, drawing) => (
  sharedSizeClass('circular', circularDiagramLengths(state), drawing)
);

/**
 * The size class of the shown mode's next Generate, or null while it is
 * unknown or separate diagrams take different ones.
 * @param {RecordLengthState} state
 * @param {DrawingState} drawing the drawing of the shown mode
 */
const activeLengthParam = (state, drawing) => (
  sharedSizeClass(shownMode(state), drawnDiagramLengths(state), drawing)
);

/**
 * The value a sliding window or a scale interval takes in the next Generate:
 * each Circular record follows its own length, and a Linear diagram its
 * longest record. Null while no length is known or the records differ.
 * @template T
 * @param {RecordLengthState} state
 * @param {(length: number) => T} valueForLength
 * @returns {T | null}
 */
const sharedLengthValue = (state, valueForLength) => {
  const diagrams = drawnDiagramLengths(state);
  const lengths = shownMode(state) === 'circular' ? diagrams.flat() : diagrams.map((records) => Math.max(...records));
  const values = new Set(lengths.map(valueForLength));
  return values.size === 1 ? [...values][0] : null;
};

/**
 * The row of `tiers` whose bound the length is under (`atMost` also accepts the bound).
 * @template {number[]} Row
 * @param {unknown} length
 * @param {readonly Row[]} tiers
 * @param {boolean} [inclusive]
 * @returns {Row | null}
 */
const tierRow = (length, tiers, inclusive = false) => {
  const numeric = Number(length);
  if (!Number.isFinite(numeric) || numeric <= 0) return null;
  return tiers.find(([bound]) => (inclusive ? numeric <= bound : numeric < bound)) || null;
};

/** @param {unknown} length */
const slidingWindowForLength = (length) => {
  const numeric = Number(length);
  if (!Number.isFinite(numeric) || numeric <= 0) return null;
  return SLIDING_WINDOW_PRESETS.find((entry) => numeric < entry.maxExclusive) || null;
};

/** @param {unknown} length */
const autoLinearTickInterval = (length) => {
  const row = tierRow(length, linearTicks.tiers);
  return row ? row[1] : (Number(length) > 0 ? linearTicks.above : null);
};

/** @param {unknown} length */
const autoCircularTickInterval = (length) => {
  const row = tierRow(length, circularTicks.tiers, true);
  return row ? row[1] : (Number(length) > 0 ? circularTicks.above[0] : null);
};

/**
 * @param {unknown} value
 * @returns {'tuckin' | 'middle' | 'spreadout'}
 */
const normalizeTrackPreset = (value) => {
  const text = String(value || '').trim().toLowerCase();
  // `includes` checks the three names, which `text` equals when it is true.
  return ['tuckin', 'middle', 'spreadout'].includes(text) ? /** @type {'tuckin' | 'middle' | 'spreadout'} */ (text) : 'tuckin';
};

/** @param {SizeClass} lengthParam */
const circularBaseTrackWidthPx = (lengthParam) => {
  const base = DEFAULTS.circular.radiusPx * DEFAULTS.circular.trackRatio;
  const factors = DEFAULTS.circular.trackRatioFactors[normalizeLengthParam(lengthParam) || 'long'];
  return { base, factors };
};

/**
 * @param {string} renderer
 * @param {SizeClass} lengthParam
 */
const circularWidthForRenderer = (renderer, lengthParam) => {
  const { base, factors } = circularBaseTrackWidthPx(lengthParam);
  if (renderer === 'features' || renderer === 'sequence_conservation') return base * Number(factors[0]);
  if (renderer === 'depth') return base * Number(factors[1]) * 0.5;
  return base * Number(factors[1]);
};

const previewSpacingPx = () => Math.max(1.0, 0.01 * DEFAULTS.circular.radiusPx);

/**
 * @template T
 * @param {RecordLengthState} state
 * @param {DrawingState} drawing
 * @param {(lengthParam: SizeClass) => T} resolver
 * @returns {T | ByLength<T>}
 */
const circularLengthValues = (state, drawing, resolver) => {
  const active = activeCircularLengthParam(state, drawing);
  if (active) return resolver(active);
  return {
    short: resolver('short'),
    long: resolver('long')
  };
};

/**
 * @param {RecordLengthState} state
 * @param {DrawingState} drawing
 * @param {string} renderer
 */
const circularWidthText = (state, drawing, renderer) => {
  const value = circularLengthValues(state, drawing, (lengthParam) => circularWidthForRenderer(renderer, lengthParam));
  return typeof value === 'number'
    ? withAuto(formatPx(value))
    : formatShortLong(value, formatPx, null);
};

/**
 * @param {DrawingState} drawing
 * @param {SizeClass} lengthParam
 */
const circularFeatureRadius = (drawing, lengthParam) => {
  const preset = normalizeTrackPreset(drawing?.form?.track_type);
  const laneWidth = circularWidthForRenderer('features', lengthParam);
  const laneCount = drawing?.form?.separate_strands ? 2 : 1;
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

/**
 * @param {string} renderer
 * @param {DrawingState} drawing
 * @returns {2 | 3 | 4 | null}
 */
const circularBuiltinTrackId = (renderer, drawing) => {
  const showDepth = Boolean(drawing?.form?.show_depth);
  const showGc = !Boolean(drawing?.form?.suppress_gc);
  if (renderer === 'depth') return showDepth ? 2 : 2;
  if (renderer === 'dinucleotide_content') return showDepth ? 3 : 2;
  if (renderer === 'dinucleotide_skew') {
    if (showDepth) return showGc ? 4 : 3;
    return showGc ? 3 : 2;
  }
  return null;
};

/**
 * @param {DrawingState} drawing
 * @param {string} renderer
 * @param {SizeClass} lengthParam
 */
const circularRadiusForRenderer = (drawing, renderer, lengthParam) => {
  if (renderer === 'features') return circularFeatureRadius(drawing, lengthParam);
  const trackId = circularBuiltinTrackId(renderer, drawing);
  if (trackId === null) return null;
  const preset = normalizeTrackPreset(drawing?.form?.track_type);
  return DEFAULTS.circular.trackDict[normalizeLengthParam(lengthParam) || 'long']?.[preset]?.[trackId] ?? null;
};

/**
 * @param {RecordLengthState} state
 * @param {DrawingState} drawing
 * @param {string} renderer
 */
const circularRadiusText = (state, drawing, renderer) => {
  const value = circularLengthValues(state, drawing, (lengthParam) => circularRadiusForRenderer(drawing, renderer, lengthParam));
  return typeof value === 'number'
    ? withAuto(formatFactor(value))
    : formatShortLong(value, formatFactor, null);
};

/**
 * @param {RecordLengthState} state
 * @param {DrawingState} drawing
 * @param {AutoValueContext | null} [trackConfig]
 */
const depthTickFontSize = (state, drawing, trackConfig = null) => {
  const mode = String(state?.mode?.value || state?.mode || 'circular');
  if (mode === 'linear') {
    const trackHeight = parsePositiveNumber(trackConfig?.height)
      ?? parsePositiveNumber(drawing?.adv?.depth_height)
      ?? DEFAULTS.linear.depthHeightPx;
    return Math.max(5, Math.min(8, trackHeight * 0.7));
  }
  const trackWidth = parsePositiveNumber(drawing?.adv?.depth_width_circular)
    ?? circularWidthForRenderer('depth', activeCircularLengthParam(state, drawing) || 'long');
  return Math.max(5, Math.min(8, trackWidth * 0.22));
};

/**
 * @param {RecordLengthState} state
 * @param {DrawingState} drawing
 */
const gcTickFontSize = (state, drawing) => {
  const mode = String(state?.mode?.value || state?.mode || 'circular');
  if (mode === 'linear') {
    const trackHeight = parsePositiveNumber(drawing?.adv?.gc_height) ?? DEFAULTS.linear.defaultGcHeightPx;
    return Math.max(5, Math.min(8, trackHeight * 0.7));
  }
  const trackWidth = parsePositiveNumber(drawing?.adv?.gc_content_width_circular)
    ?? circularWidthForRenderer('dinucleotide_content', activeCircularLengthParam(state, drawing) || 'long');
  return Math.max(5, Math.min(8, trackWidth * 0.22));
};

/**
 * @param {RecordLengthState} state
 * @param {DrawingState} drawing
 * @param {ByLength<number>} values
 */
const lengthDependentPx = (state, drawing, values) => formatShortLong(values, formatPx, activeLengthParam(state, drawing));

/**
 * @param {RecordLengthState} state
 * @param {DrawingState} drawing
 * @param {'window' | 'step'} field
 */
const dinucleotideWindowText = (state, drawing, field) => {
  const manualWindow = parsePositiveNumber(drawing?.adv?.window_size);
  const manualStep = parsePositiveNumber(drawing?.adv?.step_size);
  if (field === 'window' && manualWindow !== null) return withAuto(formatBp(manualWindow));
  if (field === 'step' && manualStep !== null) return withAuto(formatBp(manualStep));

  const preset = sharedLengthValue(state, slidingWindowForLength);
  if (preset) return withAuto(formatBp(preset[field]));
  return t('Auto (by genome size)');
};

/**
 * @param {RecordLengthState} state
 * @param {DrawingState} drawing
 * @param {'window' | 'step'} field
 */
const depthWindowText = (state, drawing, field) => {
  const sourceField = field === 'window' ? 'window_size' : 'step_size';
  const manualSource = parsePositiveNumber(drawing?.adv?.[sourceField]);
  if (manualSource !== null) {
    const divisor = field === 'window' ? 10 : 10;
    const minValue = field === 'window' ? 100 : 1;
    return withAuto(formatBp(Math.max(minValue, manualSource / divisor)));
  }

  const preset = sharedLengthValue(state, slidingWindowForLength);
  if (preset) {
    const source = field === 'window' ? preset.window : preset.step;
    const minValue = field === 'window' ? 100 : 1;
    return withAuto(formatBp(Math.max(minValue, source / 10)));
  }
  return t('Auto (by genome size)');
};

/** @param {RecordLengthState} state @param {DrawingState} drawing */
const linearAxisGapText = (state, drawing) => {
  const layout = String(drawing?.form?.linear_track_layout || 'middle').toLowerCase();
  if (layout === 'middle') return 'ignored';
  const manualHeight = parsePositiveNumber(drawing?.adv?.feature_height);
  const multiplier = drawing?.form?.separate_strands ? 0.25 : 0.30;
  if (manualHeight !== null) return withAuto(formatPx(manualHeight * multiplier));
  return formatShortLong(
    {
      short: DEFAULTS.linear.defaultCdsHeightPx.short * multiplier,
      long: DEFAULTS.linear.defaultCdsHeightPx.long * multiplier
    },
    formatPx,
    activeLengthParam(state, drawing)
  );
};

/** @param {DrawingState} drawing */
const linearLabelPlacementText = (drawing) => {
  const layout = requireCurrentLinearTrackLayout(drawing?.form?.linear_track_layout);
  return layout === 'below' ? 'below' : 'above';
};

/** @param {RecordLengthState} state */
const scaleIntervalText = (state) => {
  const interval = sharedLengthValue(
    state,
    shownMode(state) === 'circular' ? autoCircularTickInterval : autoLinearTickInterval
  );
  return interval ? withAuto(formatBp(interval)) : t('by record length (auto)');
};

// Linear labels take the size class of the records in the label scope, which
// can differ from the class of the diagram (OV-493).
/** @param {RecordLengthState} state @param {DrawingState} drawing */
const linearLabelFontSizeText = (state, drawing) => formatShortLong(
  DEFAULTS.linear.labelFontSizePx,
  formatPx,
  sizeClassForRecordLengths(
    'linear',
    linearLabelSizeLengths(linearRecordLengths(state), drawing.form.show_labels_linear),
    drawing
  )
);

/**
 * @param {RecordLengthState} state
 * @param {DrawingState} drawing
 */
const plotTitleText = (state, drawing) => {
  if (String(state?.mode?.value || state?.mode || 'circular') !== 'circular') return '';
  if (String(drawing?.adv?.plot_title_position || 'none') === 'none') return '';
  const pieces = [drawing?.form?.species, drawing?.form?.strain]
    .map((value) => String(value || '').trim())
    .filter(Boolean);
  return pieces.length ? t('{value} (auto)', { value: pieces.join(' ') }) : t('species + strain (auto)');
};

/**
 * @param {RecordLengthState} state
 * @param {DrawingState} drawing
 * @param {string} key
 * @param {AutoValueContext | null} [context]
 */
const autoTextByKey = (state, drawing, key, context = null) => {
  switch (key) {
    case 'depthWindow': return depthWindowText(state, drawing, 'window');
    case 'depthStep': return depthWindowText(state, drawing, 'step');
    case 'depthMin': return t('0 (auto)');
    case 'depthMax': return t('from depth data (auto)');
    case 'circularDepthWidth': return circularWidthText(state, drawing, 'depth');
    case 'linearDepthHeight': return withAuto(formatPx(DEFAULTS.linear.depthHeightPx));
    case 'linearDepthTrackHeight': return withAuto(formatPx(parsePositiveNumber(drawing?.adv?.depth_height) ?? DEFAULTS.linear.depthHeightPx));
    case 'depthLargeTick': return 'endpoints only';
    case 'depthTickFont': return withAuto(formatPx(depthTickFontSize(state, drawing, context)));
    case 'conservationReference': return drawing?.circularConservation?.source === 'losat' ? t('subject (auto)') : t('detect side (auto)');
    case 'conservationRingWidth': return circularWidthText(state, drawing, 'sequence_conservation');
    case 'conservationRingGap': return withAuto(formatPx(previewSpacingPx()));
    case 'plotTitle': return plotTitleText(state, drawing);
    case 'definitionFontSize':
      return String(state?.mode?.value || state?.mode || 'circular') === 'circular'
        ? withAuto(formatPx(DEFAULTS.circular.definitionFontSizePx))
        : lengthDependentPx(state, drawing, DEFAULTS.linear.definitionFontSizePx);
    case 'legendBoxSize': return lengthDependentPx(state, drawing, DEFAULTS.legend.boxSizePx);
    case 'legendFontSize': return lengthDependentPx(state, drawing, DEFAULTS.legend.fontSizePx);
    case 'linearFeatureHeight': return lengthDependentPx(state, drawing, DEFAULTS.linear.defaultCdsHeightPx);
    case 'linearAxisGap': return linearAxisGapText(state, drawing);
    case 'circularFeatureWidth': return circularWidthText(state, drawing, 'features');
    case 'blockStrokeWidth': return lengthDependentPx(state, drawing, DEFAULTS.feature.blockStrokeWidthPx);
    case 'lineStrokeWidth': return lengthDependentPx(state, drawing, DEFAULTS.feature.lineStrokeWidthPx);
    case 'labelFontSize':
      return String(state?.mode?.value || state?.mode || 'circular') === 'circular'
        ? lengthDependentPx(state, drawing, DEFAULTS.circular.labelFontSizePx)
        : linearLabelFontSizeText(state, drawing);
    case 'labelRendering': return t('per label');
    case 'linearLabelPlacement': return linearLabelPlacementText(drawing);
    case 'circularLabelSpacing': return withAuto(formatPx(DEFAULTS.labels.spacingPx));
    case 'linearLabelSpacing': return withAuto(formatPx(DEFAULTS.labels.spacingPx));
    case 'labelRadiusOffset': return withAuto(roundDisplayNumber(DEFAULTS.labels.radiusOffset, 1));
    case 'axisStrokeWidth':
      return String(state?.mode?.value || state?.mode || 'circular') === 'circular'
        ? lengthDependentPx(state, drawing, DEFAULTS.circular.axisStrokeWidthPx)
        : lengthDependentPx(state, drawing, DEFAULTS.linear.axisStrokeWidthPx);
    case 'scaleStrokeWidth': return withAuto(formatPx(DEFAULTS.scale.strokeWidthPx));
    case 'scaleFontSize': return lengthDependentPx(state, drawing, DEFAULTS.linear.scaleFontSizePx);
    case 'rulerLabelFontSize': return lengthDependentPx(state, drawing, DEFAULTS.linear.rulerLabelFontSizePx);
    case 'scaleInterval': return scaleIntervalText(state);
    case 'tickLabelFontSize': return withAuto(formatPx(DEFAULTS.circular.tickLabelFontSizePx));
    case 'dinucleotideWindow': return dinucleotideWindowText(state, drawing, 'window');
    case 'dinucleotideStep': return dinucleotideWindowText(state, drawing, 'step');
    case 'gcTickFont': return withAuto(formatPx(gcTickFontSize(state, drawing)));
    case 'linearGcHeight': return withAuto(formatPx(DEFAULTS.linear.defaultGcHeightPx));
    case 'circularGcWidth': return circularWidthText(state, drawing, 'dinucleotide_content');
    case 'circularGcRadius': return circularRadiusText(state, drawing, 'dinucleotide_content');
    case 'circularSkewWidth': return circularWidthText(state, drawing, 'dinucleotide_skew');
    case 'circularSkewRadius': return circularRadiusText(state, drawing, 'dinucleotide_skew');
    case 'pairwiseMatchHeight': return withAuto(formatPx(DEFAULTS.linear.comparisonHeightPx));
    default: return '';
  }
};

/**
 * @param {RecordLengthState & { activeDrawing: () => DrawingState }} state the Web state of `state.js`
 * @param {Partial<DrawnLengthSources>} [sources] the records the next Generate draws
 */
export const createAutoValueDisplay = (state, sources = {}) => {
  drawnLengthSources.set(state, { circularDiagramLengths: () => [], linearRecordLengths: () => [], ...sources });
  /**
   * @param {string} key
   * @param {AutoValueContext | null} [context]
   */
  const autoValueText = (key, context = null) => {
    const drawing = state.activeDrawing();
    return autoTextByKey(state, drawing, key, context);
  };
  /**
   * @param {unknown} value
   * @param {string} key
   * @param {AutoValueContext | null} [context]
   */
  const autoValueVisible = (value, key, context = null) => (
    normalizeOptionalText(value) === null && autoValueText(key, context) !== ''
  );

  return {
    autoValueText,
    autoValueVisible
  };
};
