// @ts-check
/** @import { DrawingState } from '../state.js' */
/** @import { TrackSlotRow } from '../services/track-slot-validation.js' */
import {
  linearLabelSizeLengths,
  normalizeOptionalText,
  sharedSizeClass,
  sizeClassForRecordLengths
} from '../services/track-slot-display.js';
import { requireCurrentLinearTrackLayout } from '../services/current-option-values.js';
import { depthSlotTrackIndex } from '../services/depth-track-state.js';
import { draftDrawsDepth } from '../services/track-legend-rows.js';
import { parseOptionalCircularScalar } from '../services/track-slot-validation.js';
import { t } from '../utils/i18n.js';
import { MODE_PROFILE_DATA } from '../mode-profiles.generated.js';

const DEFAULTS = Object.freeze({
  circular: {
    radiusPx: 390,
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

// The length tiers and divisor of the Auto sliding window, the Depth window rule,
// the scale tick tiers, the Auto track widths and the scalar-axis tick font rule
// are Python's (gbdraw/auto_sizes.py), carried by the generated mode-profile data.
const { slidingWindow, depthWindow, circularTicks, linearTicks, circularSlots, scalarAxisTickFont } = MODE_PROFILE_DATA.autoSizes;
const SLIDING_WINDOW_TIERS = [...slidingWindow.tiers, [Infinity, ...slidingWindow.above]];

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

/** A length in bp as the Auto texts show it (`10 kbp`). @param {unknown} value */
export const formatBp = (value) => {
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
 * The Depth series whose Auto text is asked for (a row of the Depth panel): its index, and its config, of which only
 * a Linear series' height is read.
 * @typedef {{ index: number, config: { height?: unknown } }} AutoValueContext
 */

/**
 * The records the next Generate draws, by length. The notes do not choose
 * records: the app registers these readers in `createAutoValueDisplay`, and
 * each returns an empty list while a drawn length is not known.
 * @typedef {object} DrawnLengthSources
 * @property {() => number[][]} circularDiagramLengths the records of each Circular diagram
 * @property {() => number[]} linearRecordLengths the records of the Linear diagram, in drawing order
 * @property {() => number | null | undefined} linearWrappedAxisLength while Wrap Rows is on, the bases on the
 *   longest drawn axis of the shown Result that drew this draft (Python's `drawnAxisLength`), null while it is not
 *   known; undefined while no record is cut, so the axis is the longest record
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

/**
 * The lengths sizes and ticks follow: a wrapped Linear drawing's drawn axis (sliding windows keep the whole
 * records, `drawnDiagramLengths`). Python owns the axis; [] while it is not known.
 * @param {RecordLengthState} state
 */
const drawnAxisDiagramLengths = (state) => {
  const axis = shownMode(state) === 'linear' ? drawnLengthSources.get(state)?.linearWrappedAxisLength() : undefined;
  if (axis === undefined) return drawnDiagramLengths(state);
  return axis === null ? [] : [[axis]];
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
  sharedSizeClass(shownMode(state), drawnAxisDiagramLengths(state), drawing)
);

/**
 * The value a sliding window or a scale interval takes in the next Generate:
 * each Circular record follows its own length, and a Linear diagram its
 * longest record. Null while no length is known or the records differ.
 * @template T
 * @param {RecordLengthState} state
 * @param {(length: number) => T} valueForLength
 * @param {number[][]} [diagrams] The lengths it follows: the whole records, or the drawn axis for ticks.
 * @returns {T | null}
 */
const sharedLengthValue = (state, valueForLength, diagrams = drawnDiagramLengths(state)) => {
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

/**
 * Python's `sliding_window`: the tier's window and step, the window at most the length
 * divided by the divisor and the step at most the window divided by it, each at least 1 bp.
 * @param {unknown} length
 * @returns {Readonly<{ window: number, step: number }> | null}
 */
const slidingWindowForLength = (length) => {
  const numeric = Number(length);
  if (!Number.isFinite(numeric) || numeric <= 0) return null;
  const tier = SLIDING_WINDOW_TIERS.find(([below]) => numeric < below) ?? SLIDING_WINDOW_TIERS[SLIDING_WINDOW_TIERS.length - 1];
  const window = Math.max(1, Math.min(tier[1], Math.floor(numeric / slidingWindow.divisor)));
  const step = Math.max(1, Math.min(tier[2], Math.floor(window / slidingWindow.divisor)));
  return Object.freeze({ window, step });
};

/**
 * Python's `depth_sliding_window` for one field: the GC value divided by the divisor, in whole base pairs,
 * the window at least the minimum window and the step at least 1 bp.
 * @param {number} source The GC window or step.
 * @param {'window' | 'step'} field
 */
const depthWindowValue = (source, field) =>
  Math.max(field === 'window' ? depthWindow.minWindow : 1, Math.floor(source / depthWindow.divisor));

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

/** @typedef {keyof typeof circularSlots.widthPx.long} CircularAutoWidthRenderer */

/**
 * The Auto width of a Circular track (Python's `circular_track_width_px`).
 * @param {CircularAutoWidthRenderer} renderer
 * @param {SizeClass} lengthParam
 */
const circularWidthForRenderer = (renderer, lengthParam) => (
  circularSlots.widthPx[normalizeLengthParam(lengthParam) || 'long'][renderer]
);

/**
 * The Auto tick font of a GC content or Depth axis on a track this thick (Python's `scalar_axis_tick_font_size`).
 * @param {'circular' | 'linear'} mode
 * @param {number} thicknessPx
 */
const scalarAxisTickFontSize = (mode, thicknessPx) => {
  const [low, high] = scalarAxisTickFont.rangePx;
  return Math.max(low, Math.min(high, thicknessPx * scalarAxisTickFont.thicknessFraction[mode]));
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
 * @param {CircularAutoWidthRenderer} renderer
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
 * The first enabled Circular Depth row of the custom stack that draws Depth series `trackIndex` (by the series it
 * references), or the first that draws a series when `trackIndex` is null; null when the stack is off or has none.
 * @param {DrawingState} drawing
 * @param {number | null} trackIndex
 */
const customStackDepthSlot = ({ adv }, trackIndex) => {
  if (!adv.circular_track_slots_enabled || !Array.isArray(adv.circular_track_slots)) return null;
  /** @type {TrackSlotRow[]} */
  const slots = adv.circular_track_slots;
  return slots.find((slot) => (
    Boolean(slot)
    && slot?.enabled !== false
    && slot?.renderer === 'depth'
    && (trackIndex === null ? depthSlotTrackIndex(slot) !== null : depthSlotTrackIndex(slot) === trackIndex)
  )) ?? null;
};

/**
 * The width (px) a Circular Depth track of series `trackIndex` is drawn at, as Python sizes its tick font from it (the
 * slot's `requested_width_px`): with Use custom stack on, the width of the series' first Depth row (a ×R width of the
 * packaged radius); otherwise the Depth Width shortcut; else the Auto width.
 * @param {RecordLengthState} state
 * @param {DrawingState} drawing
 * @param {number} trackIndex
 */
const circularDepthWidthPx = (state, drawing, trackIndex) => {
  const autoWidth = circularWidthForRenderer('depth', activeCircularLengthParam(state, drawing) || 'long');
  if (!drawing.adv.circular_track_slots_enabled) return parsePositiveNumber(drawing.adv.depth_width_circular) ?? autoWidth;
  const slot = customStackDepthSlot(drawing, trackIndex);
  let width = null;
  try {
    width = slot ? parseOptionalCircularScalar(slot.width) : null;
  } catch {
    // A width the row cannot use draws nothing; the note keeps the Auto width.
  }
  if (width === null) return autoWidth;
  return width.unit === 'px' ? width.value : width.value * DEFAULTS.circular.radiusPx;
};

/**
 * The tick font Circular Depth series `trackIndex` draws its numbers in: its own Tick Font Size, else the drawing's
 * Depth tick font, else the Auto size of its track's width (Python's `depth_track_config` and
 * `depth_axis_tick_font_size_px`).
 * @param {RecordLengthState} state
 * @param {DrawingState} drawing
 * @param {number} trackIndex
 */
const circularDepthSeriesTickFontPx = (state, drawing, trackIndex) => {
  const series = Array.isArray(drawing.adv.depth_tracks) ? drawing.adv.depth_tracks[trackIndex] : null;
  return parsePositiveNumber(series?.tick_font_size)
    ?? drawingDepthTick(drawing, 'tick_font_size')
    ?? scalarAxisTickFontSize('circular', circularDepthWidthPx(state, drawing, trackIndex));
};

/**
 * The drawing's own Depth tick (`adv.depth_<field>`, which a configuration or a Session from the CLI or an older
 * writer sets), or null. A Depth series without its own tick is drawn with it (Python's `DepthConfigurator` reads a
 * series' null tick as the drawing's), so it is the Auto value of each series' field (OV-3496).
 * @param {DrawingState} drawing
 * @param {'large_tick_interval' | 'small_tick_interval' | 'tick_font_size'} field
 */
const drawingDepthTick = (drawing, field) => parsePositiveNumber(drawing.adv[`depth_${field}`]);

/**
 * The Auto text of a Depth series' tick interval field: the drawing's own tick, else `fallback`.
 * @param {DrawingState} drawing
 * @param {'large_tick_interval' | 'small_tick_interval'} field
 * @param {string} fallback
 */
const depthSeriesTickText = (drawing, field, fallback) => {
  const tick = drawingDepthTick(drawing, field);
  return tick === null ? fallback : withAuto(roundDisplayNumber(tick, 2));
};

/**
 * The font a Depth series without its own Tick Font Size draws its numbers in: the drawing's Depth tick font, else the
 * Auto size of its track (Linear: the series' height; Circular: its track's width).
 * @param {RecordLengthState} state
 * @param {DrawingState} drawing
 * @param {AutoValueContext | null} [trackConfig]
 */
const depthTickFontSize = (state, drawing, trackConfig = null) => {
  const drawingFont = drawingDepthTick(drawing, 'tick_font_size');
  if (drawingFont !== null) return drawingFont;
  if (shownMode(state) === 'linear') {
    const trackHeight = parsePositiveNumber(trackConfig?.config.height)
      ?? parsePositiveNumber(drawing?.adv?.depth_height)
      ?? DEFAULTS.linear.depthHeightPx;
    return scalarAxisTickFontSize('linear', trackHeight);
  }
  // The Auto value of the series' own Tick Font Size: the size of its track's width.
  return scalarAxisTickFontSize('circular', circularDepthWidthPx(state, drawing, trackConfig?.index ?? 0));
};

/**
 * The font a Percent GC axis takes from the depth axis drawn with ticks beside it, as Python's
 * `CircularScalarAxes.gc_tick_font_size` draws it: the font of the first Depth track the stack draws (pinned by
 * `gcAxisTickFonts` of tests/fixtures/auto_size_vectors.json); null when the GC axis keeps the Auto size of its own
 * track.
 * @param {RecordLengthState} state
 * @param {DrawingState} drawing
 */
const gcTickFontOfDepth = (state, drawing) => {
  const { adv } = drawing;
  if (adv.gc_content_mode !== 'percent' || !draftDrawsDepth('circular', drawing)) return null;
  if (!(adv.depth_show_axis && adv.depth_show_ticks)) return null;
  const firstSlot = customStackDepthSlot(drawing, null);
  return circularDepthSeriesTickFontPx(state, drawing, (firstSlot && depthSlotTrackIndex(firstSlot)) ?? 0);
};

/**
 * @param {RecordLengthState} state
 * @param {DrawingState} drawing
 */
const gcTickFontSize = (state, drawing) => {
  if (shownMode(state) === 'linear') {
    const trackHeight = parsePositiveNumber(drawing?.adv?.gc_height) ?? DEFAULTS.linear.defaultGcHeightPx;
    return scalarAxisTickFontSize('linear', trackHeight);
  }
  const trackWidth = parsePositiveNumber(drawing?.adv?.gc_content_width_circular)
    ?? circularWidthForRenderer('dinucleotide_content', activeCircularLengthParam(state, drawing) || 'long');
  return gcTickFontOfDepth(state, drawing) ?? scalarAxisTickFontSize('circular', trackWidth);
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

  // Each field by itself: the records may share a window and not a step (Python's record_auto_values).
  const value = sharedLengthValue(state, (length) => slidingWindowForLength(length)?.[field] ?? null);
  return value !== null ? withAuto(formatBp(value)) : t('Auto (by genome size)');
};

/**
 * @param {RecordLengthState} state
 * @param {DrawingState} drawing
 * @param {'window' | 'step'} field
 */
const depthWindowText = (state, drawing, field) => {
  const sourceField = field === 'window' ? 'window_size' : 'step_size';
  const manualSource = parsePositiveNumber(drawing?.adv?.[sourceField]);
  if (manualSource !== null) return withAuto(formatBp(depthWindowValue(manualSource, field)));

  // Records of different GC windows can share a Depth window (at least 100 bp; Python's record_auto_values).
  const value = sharedLengthValue(state, (length) => {
    const preset = slidingWindowForLength(length);
    return preset ? depthWindowValue(preset[field], field) : null;
  });
  return value !== null ? withAuto(formatBp(value)) : t('Auto (by genome size)');
};

/** @param {RecordLengthState} state @param {DrawingState} drawing */
const linearAxisGapText = (state, drawing) => {
  const layout = String(drawing?.form?.linear_track_layout || 'middle').toLowerCase();
  if (layout === 'middle') return t('ignored');
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
  return layout === 'below' ? t('below') : t('above');
};

/** @param {RecordLengthState} state */
const scaleIntervalText = (state) => {
  const interval = sharedLengthValue(
    state,
    shownMode(state) === 'circular' ? autoCircularTickInterval : autoLinearTickInterval,
    drawnAxisDiagramLengths(state)
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
    drawnLengthSources.get(state)?.linearWrappedAxisLength() === undefined
      ? linearLabelSizeLengths(linearRecordLengths(state), drawing.form.show_labels_linear)
      // Wrapped: the labels follow the drawn axis (Python's `linear_label_font_size`).
      : drawnAxisDiagramLengths(state).flat(),
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
    case 'depthLargeTick': return depthSeriesTickText(drawing, 'large_tick_interval', t('endpoints only'));
    case 'depthSmallTick': return depthSeriesTickText(drawing, 'small_tick_interval', t('None'));
    case 'depthTickFont': return withAuto(formatPx(depthTickFontSize(state, drawing, context)));
    case 'conservationReference': return drawing?.circularConservation?.source === 'losat' ? t('subject') : t('detect side');
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
  drawnLengthSources.set(state, {
    circularDiagramLengths: () => [], linearRecordLengths: () => [], linearWrappedAxisLength: () => undefined, ...sources
  });
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
