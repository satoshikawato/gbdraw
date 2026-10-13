// @ts-check
import { LEGEND_POSITIONS } from '../contracts/schema.generated.js';
import { WEB_UX_PROFILE } from '../web-ux-profile.js';

/**
 * A legend side and a title side, as a diagram of one mode shows them.
 * @typedef {{ legend: string, plotTitlePosition: string }} ActiveLayoutPreference
 */
/**
 * The circular multi-record sides; null follows the single-record side.
 * @typedef {{ legend: string | null, plotTitlePosition: string | null }} MultiLayoutPreference
 */
/**
 * The layout of both modes (`drawing.layoutPreferences`).
 * @typedef {{
 *   circular: { single: ActiveLayoutPreference, multi: MultiLayoutPreference },
 *   linear: ActiveLayoutPreference
 * }} LayoutPreferences
 */
/**
 * Stored layout preferences before normalization: a branch or a side may be absent.
 * @typedef {{
 *   circular?: { single?: StoredLayoutSides | null, multi?: StoredLayoutSides | null } | null,
 *   linear?: StoredLayoutSides | null
 * }} StoredLayoutPreferences
 */
/** @typedef {{ legend?: unknown, plotTitlePosition?: unknown }} StoredLayoutSides */

/** @param {unknown} value */
export const normalizeCircularPlotTitlePosition = (value) => {
  const normalized = String(value || '').trim().toLowerCase();
  return ['none', 'top', 'bottom'].includes(normalized) ? normalized : 'none';
};

/** @param {unknown} value */
export const normalizeLinearPlotTitlePosition = (value) => {
  const normalized = String(value || '').trim().toLowerCase();
  return ['center', 'top', 'bottom'].includes(normalized) ? normalized : 'bottom';
};

/**
 * The legend position of a mode: `value` without case and outer spaces when it
 * is one of the mode's positions (Python names them), otherwise null. A value
 * outside the set reads like an absent one, so each caller keeps the fallback
 * it has for an unset position (OV-652).
 * @param {unknown} value
 * @param {string} mode 'linear'; any other mode reads the Circular positions.
 * @returns {string | null}
 */
export const legendPositionOfMode = (value, mode) => {
  const position = typeof value === 'string' ? value.trim().toLowerCase() : '';
  return LEGEND_POSITIONS[mode === 'linear' ? 'linear' : 'circular'].includes(position) ? position : null;
};

/**
 * @param {unknown} value
 * @returns {value is string}
 */
const storedText = (value) => typeof value === 'string' && value.trim() !== '';

/**
 * The legend position of a mode, or the mode's default legend position for a
 * value outside the mode's positions (OV-652, OV-671).
 * @param {unknown} value
 * @param {string} mode 'linear'; any other mode reads the Circular positions.
 * @returns {string}
 */
export const legendPositionOrDefault = (value, mode) => (
  legendPositionOfMode(value, mode) ?? WEB_UX_PROFILE[mode === 'linear' ? 'linear' : 'circular'].legend
);

/**
 * The side a Result saved without composition metadata (Session 27 to 30) is
 * laid out with: `none` when its Session stored no side (no value, or one the
 * old reader skipped: `null`, `false`, blank text), otherwise the stored
 * position of the Result's mode or the mode's default (OV-671, BA-D26).
 * @param {unknown} value The Session's `config.form.legend`.
 * @param {string} mode The Result's mode.
 * @returns {string}
 */
export const drawnLegendSideOfMode = (value, mode) => (
  !value || (typeof value === 'string' && value.trim() === '') ? 'none' : legendPositionOrDefault(value, mode)
);

/** @returns {LayoutPreferences} */
export const createDefaultLayoutPreferences = () => ({
  circular: {
    single: {
      legend: WEB_UX_PROFILE.circular.legend,
      plotTitlePosition: WEB_UX_PROFILE.circular.plotTitlePosition
    },
    multi: {
      legend: null,
      plotTitlePosition: null
    }
  },
  linear: {
    legend: WEB_UX_PROFILE.linear.legend,
    plotTitlePosition: WEB_UX_PROFILE.linear.plotTitlePosition
  }
});

/**
 * @param {StoredLayoutPreferences | null | undefined} preferences
 * @param {boolean} [useMultiRecord]
 * @returns {ActiveLayoutPreference}
 */
export const resolveCircularLayoutPreference = (preferences, useMultiRecord = false) => {
  const defaults = createDefaultLayoutPreferences();
  const circular = preferences?.circular || defaults.circular;
  const single = {
    legend: legendPositionOrDefault(circular.single?.legend, 'circular'),
    plotTitlePosition: normalizeCircularPlotTitlePosition(
      circular.single?.plotTitlePosition
    )
  };
  if (!useMultiRecord) return single;
  return {
    legend: legendPositionOfMode(circular.multi?.legend, 'circular') ?? single.legend,
    plotTitlePosition: storedText(circular.multi?.plotTitlePosition)
      ? normalizeCircularPlotTitlePosition(circular.multi.plotTitlePosition)
      : single.plotTitlePosition
  };
};

/**
 * @param {StoredLayoutPreferences | null | undefined} preferences
 * @param {string} mode
 * @param {boolean} [useMultiRecord]
 * @returns {ActiveLayoutPreference}
 */
export const resolveActiveLayoutPreference = (
  preferences,
  mode,
  useMultiRecord = false
) => {
  if (mode === 'linear') {
    return {
      legend: legendPositionOrDefault(preferences?.linear?.legend, 'linear'),
      plotTitlePosition: normalizeLinearPlotTitlePosition(
        preferences?.linear?.plotTitlePosition
      )
    };
  }
  return resolveCircularLayoutPreference(preferences, useMultiRecord);
};

/**
 * @param {LayoutPreferences} preferences
 * @param {string} mode
 * @param {boolean} useMultiRecord
 * @param {StoredLayoutSides} patch
 */
export const updateActiveLayoutPreference = (
  preferences,
  mode,
  useMultiRecord,
  patch
) => {
  if (mode === 'linear') {
    if (Object.prototype.hasOwnProperty.call(patch, 'legend')) {
      preferences.linear.legend = legendPositionOrDefault(patch.legend, 'linear');
    }
    if (Object.prototype.hasOwnProperty.call(patch, 'plotTitlePosition')) {
      preferences.linear.plotTitlePosition = normalizeLinearPlotTitlePosition(
        patch.plotTitlePosition
      );
    }
    return;
  }

  const key = useMultiRecord ? 'multi' : 'single';
  if (Object.prototype.hasOwnProperty.call(patch, 'legend')) {
    preferences.circular[key].legend = legendPositionOrDefault(patch.legend, 'circular');
  }
  if (Object.prototype.hasOwnProperty.call(patch, 'plotTitlePosition')) {
    preferences.circular[key].plotTitlePosition = normalizeCircularPlotTitlePosition(
      patch.plotTitlePosition
    );
  }
};

/**
 * @param {unknown} source Saved layout preferences, validated here.
 * @returns {LayoutPreferences}
 */
export const normalizeLayoutPreferences = (source) => {
  const defaults = createDefaultLayoutPreferences();
  if (!source || typeof source !== 'object' || Array.isArray(source)) {
    return defaults;
  }
  // The guard above makes `source` a non-array object, the stored shape (`StoredLayoutPreferences`).
  const stored = /** @type {StoredLayoutPreferences} */ (source);
  const rawMulti = stored.circular?.multi;
  return {
    circular: {
      single: resolveCircularLayoutPreference(stored, false),
      multi: {
        legend: legendPositionOfMode(rawMulti?.legend, 'circular'),
        plotTitlePosition: storedText(rawMulti?.plotTitlePosition)
          ? normalizeCircularPlotTitlePosition(rawMulti.plotTitlePosition)
          : null
      }
    },
    linear: resolveActiveLayoutPreference(stored, 'linear')
  };
};

/**
 * @param {LayoutPreferences} target
 * @param {unknown} source
 */
export const replaceLayoutPreferences = (target, source) => {
  const normalized = normalizeLayoutPreferences(source);
  Object.assign(target.circular.single, normalized.circular.single);
  Object.assign(target.circular.multi, normalized.circular.multi);
  Object.assign(target.linear, normalized.linear);
};

/**
 * @param {Record<string, unknown> | null | undefined} ui
 * @param {{ mode?: string, multiRecord?: boolean, activeLegend?: string, activePlotTitlePosition?: string }} [context]
 * @returns {LayoutPreferences}
 */
export const migrateLegacyLayoutPreferences = (
  ui,
  {
    mode = 'circular',
    multiRecord = false,
    activeLegend = mode === 'linear'
      ? WEB_UX_PROFILE.linear.legend
      : WEB_UX_PROFILE.circular.legend,
    activePlotTitlePosition = mode === 'linear'
      ? WEB_UX_PROFILE.linear.plotTitlePosition
      : WEB_UX_PROFILE.circular.plotTitlePosition
  } = {}
) => {
  if (ui?.layoutPreferences) {
    return normalizeLayoutPreferences(ui.layoutPreferences);
  }

  // Each legend slot takes its first source that holds a position of its mode.
  const activeCircularLegend = legendPositionOfMode(activeLegend, 'circular');
  const circularLegend = legendPositionOfMode(ui?.circularLegendPosition, 'circular')
    ?? legendPositionOfMode(ui?.legend, 'circular')
    ?? activeCircularLegend
    ?? 'left';
  const circularPlotTitle = storedText(ui?.circularPlotTitlePosition)
    ? normalizeCircularPlotTitlePosition(ui.circularPlotTitlePosition)
    : normalizeCircularPlotTitlePosition(activePlotTitlePosition);
  const single = {
    legend: legendPositionOfMode(ui?.circularSingleRecordLegendPosition, 'circular')
      ?? (multiRecord ? null : activeCircularLegend)
      ?? circularLegend,
    plotTitlePosition: storedText(ui?.circularSingleRecordPlotTitlePosition)
      ? normalizeCircularPlotTitlePosition(ui.circularSingleRecordPlotTitlePosition)
      : multiRecord
        ? circularPlotTitle
        : normalizeCircularPlotTitlePosition(activePlotTitlePosition)
  };
  const multi = {
    legend: legendPositionOfMode(ui?.circularMultiRecordLegendPosition, 'circular')
      ?? (multiRecord ? activeCircularLegend : null)
      ?? circularLegend,
    plotTitlePosition: storedText(ui?.circularMultiRecordPlotTitlePosition)
      ? normalizeCircularPlotTitlePosition(ui.circularMultiRecordPlotTitlePosition)
      : multiRecord
        ? normalizeCircularPlotTitlePosition(activePlotTitlePosition)
        : circularPlotTitle
  };
  return {
    circular: { single, multi },
    linear: {
      legend: legendPositionOfMode(ui?.linearLegendPosition, 'linear')
        ?? (mode === 'linear' ? legendPositionOfMode(activeLegend, 'linear') : null)
        ?? 'bottom',
      plotTitlePosition: storedText(ui?.linearPlotTitlePosition)
        ? normalizeLinearPlotTitlePosition(ui.linearPlotTitlePosition)
        : mode === 'linear'
          ? normalizeLinearPlotTitlePosition(activePlotTitlePosition)
          : 'bottom'
    }
  };
};

// The `ui` fields of the layout before Session 44 kept `ui.layoutPreferences`.
export const LEGACY_LAYOUT_PREFERENCE_FIELDS = Object.freeze([
  'legend',
  'circularLegendPosition',
  'linearLegendPosition',
  'circularPlotTitlePosition',
  'linearPlotTitlePosition',
  'circularSingleRecordLegendPosition',
  'circularSingleRecordPlotTitlePosition',
  'circularMultiRecordLegendPosition',
  'circularMultiRecordPlotTitlePosition'
]);

/**
 * @param {unknown} value
 * @returns {value is Record<string, unknown>}
 */
const isRecord = (value) => Boolean(value) && typeof value === 'object' && !Array.isArray(value);

// The legacy legend fields and the mode whose positions each one holds. The
// migration reads `legend` as a Circular position in every mode (and as the
// Linear one in Linear, whose positions are Circular's).
/** @type {Readonly<Record<string, 'circular' | 'linear'>>} */
const LEGACY_LEGEND_FIELD_MODES = Object.freeze({
  legend: 'circular',
  circularLegendPosition: 'circular',
  circularSingleRecordLegendPosition: 'circular',
  circularMultiRecordLegendPosition: 'circular',
  linearLegendPosition: 'linear'
});

// Partial objects remain authoritative for session compatibility; normalization
// supplies the current defaults for omitted branches. A legend field outside
// its mode's positions is not a stored layout (OV-652).
/** @param {Record<string, unknown> | null | undefined} ui */
const hasStoredLayoutPreferences = (ui) => (
  isRecord(ui?.layoutPreferences) ||
  LEGACY_LAYOUT_PREFERENCE_FIELDS.some((field) => (
    Object.hasOwn(LEGACY_LEGEND_FIELD_MODES, field)
      ? legendPositionOfMode(ui?.[field], LEGACY_LEGEND_FIELD_MODES[field]) !== null
      : storedText(ui?.[field])
  ))
);

/**
 * The layout preferences a Session 27-44 restores (Session Load and Gallery
 * publication). A saved layout owner (current or legacy `ui` fields) wins.
 * Without one, a canonical Session takes the layout projected from its
 * committed request (`projected`). Legacy fields migrate with the committed
 * values (canonical) or `active` (other payloads) as their fallback.
 * @param {Record<string, unknown> | null | undefined} ui
 * @param {{ mode: string, multiRecord: boolean, projected?: StoredLayoutPreferences | null,
 *   active: ActiveLayoutPreference }} context
 * @returns {LayoutPreferences}
 */
export const restoredLayoutPreferences = (ui, { mode, multiRecord, projected = null, active }) => {
  if (projected && !hasStoredLayoutPreferences(ui)) return normalizeLayoutPreferences(projected);
  const fallback = projected ? resolveActiveLayoutPreference(projected, mode, multiRecord) : active;
  const migrationUi = (
    !isRecord(ui?.layoutPreferences) &&
    mode === 'linear' &&
    legendPositionOfMode(ui?.linearLegendPosition, 'linear') === null &&
    storedText(ui?.legend)
  )
    ? { ...ui, linearLegendPosition: ui.legend }
    : ui;
  return migrateLegacyLayoutPreferences(migrationUi, {
    mode,
    multiRecord,
    activeLegend: fallback.legend,
    activePlotTitlePosition: fallback.plotTitlePosition
  });
};
