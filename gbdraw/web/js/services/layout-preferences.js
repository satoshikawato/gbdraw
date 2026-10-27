// @ts-check
import { LEGEND_POSITIONS } from '../contracts/schema.generated.js';
import { WEB_UX_PROFILE } from '../web-ux-profile.js';
import { loadRepair } from '../utils/load-repairs.js';
import { changedLeaves, textChoice } from '../utils/field-reading.js';

/** @import { LoadRepair, Repaired } from '../utils/load-repairs.js' */

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

export const readCircularPlotTitlePosition = textChoice(['none', 'top', 'bottom'], 'none');
/** @param {unknown} value */
export const normalizeCircularPlotTitlePosition = (value) => readCircularPlotTitlePosition(value).value;

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
 * A stored value Load reads as a legend position: a falsy value (absent,
 * `null`, `false`, 0, empty text) and blank text are unset (no repair); case
 * and outer spaces read as the position.
 * @param {unknown} value
 */
const storedLegendValue = (value) => Boolean(value) && !(typeof value === 'string' && value.trim() === '');

/**
 * The legend position of a mode that Load reads from a file: the mode's
 * default for a value outside the mode's positions, listed in the Load dialog
 * at `path` (OV-652, BA-D11, BA-D25).
 * @param {unknown} value
 * @param {string} mode 'linear'; any other mode reads the Circular positions.
 * @param {string} path Where the value is in the file.
 * @param {'legend-position-defaulted' | 'legend-side-defaulted'} [code] The side a saved Result is drawn with
 *   is `legend-side-defaulted`.
 * @returns {Repaired<string>}
 */
export const readLegendPositionOfMode = (value, mode, path, code = 'legend-position-defaulted') => {
  const position = legendPositionOrDefault(value, mode);
  return {
    value: position,
    repairs: storedLegendValue(value) && legendPositionOfMode(value, mode) === null
      ? [loadRepair({ code, path, effect: 'defaulted', stored: String(value), used: position })]
      : []
  };
};

/**
 * The side a Result saved without composition metadata (Session 27 to 30) is
 * laid out with: `none` when its Session stored no side (no value, a falsy
 * value the old `||` chain read as `none`, or blank text, which failed Load
 * before), otherwise the stored position of the Result's mode or the mode's
 * default (OV-671, BA-D26).
 * @param {unknown} value The Session's `config.form.legend`.
 * @param {string} mode The Result's mode.
 * @returns {string}
 */
export const drawnLegendSideOfMode = (value, mode) => (
  !value || (typeof value === 'string' && value.trim() === '') ? 'none' : legendPositionOrDefault(value, mode)
);

/**
 * `drawnLegendSideOfMode` as Load reads it from a file: a stored side outside
 * the Result's mode is listed in the Load dialog at `path` (BA-D26).
 * @param {unknown} value
 * @param {string} mode
 * @param {string} path
 * @returns {Repaired<string>}
 */
export const readDrawnLegendSideOfMode = (value, mode, path) => {
  const side = drawnLegendSideOfMode(value, mode);
  return {
    value: side,
    repairs: storedLegendValue(value) && legendPositionOfMode(value, mode) === null
      ? [loadRepair({ code: 'legend-side-defaulted', path, effect: 'defaulted', stored: String(value), used: side })]
      : []
  };
};

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
 * Whether a circular multi-record side is stored (a null side follows the
 * single-record side).
 * @param {LayoutPreferences} preferences
 */
export const hasStoredCircularMultiRecordLayout = (preferences) => (
  storedText(preferences.circular.multi.legend) || storedText(preferences.circular.multi.plotTitlePosition)
);

/**
 * The multi-record sides the Multi-Record Canvas toggle stores when it turns
 * the canvas on: the single-record sides, with a left Legend and no title
 * moved below the grid. Null when a multi-record side is stored, so the
 * toggle stores them once.
 * @param {LayoutPreferences} preferences
 * @returns {ActiveLayoutPreference | null}
 */
export const circularMultiRecordDefaults = (preferences) => {
  if (hasStoredCircularMultiRecordLayout(preferences)) return null;
  const single = resolveCircularLayoutPreference(preferences, false);
  return {
    legend: single.legend === 'left' ? 'bottom' : single.legend,
    plotTitlePosition: single.plotTitlePosition === 'none' ? 'bottom' : single.plotTitlePosition
  };
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
 * The layout values stored in a Session's `ui` that Load does not use as
 * stored. `ui.layoutPreferences`, when stored, is the layout: each value it
 * reads otherwise (a legend outside its mode's positions, OV-652; an unknown
 * plot title position) is listed with the value used. Otherwise the legacy
 * legend fields are, and a legend outside its mode's positions is not used.
 * A Session 46 mode slice stores its mode's slot only (`slotOf`); a slot that
 * is no object is not read, and the committed request's layout, when it has
 * one (`projected`), is used instead without a loss.
 * @param {unknown} ui
 * @param {string} path Where `ui` is in the file.
 * @param {{ slotOf?: 'circular' | 'linear' | null, projected?: boolean }} [options]
 * @returns {LoadRepair[]}
 */
export const storedLegendPositionRepairs = (ui, path, { slotOf = null, projected = false } = {}) => {
  if (!isRecord(ui)) return [];
  const layout = ui.layoutPreferences;
  if (isRecord(layout)) {
    const stored = slotOf ? { [slotOf]: layout } : layout;
    const read = normalizeLayoutPreferences(stored);
    return changedLeaves(slotOf ? layout : stored, slotOf ? read[slotOf] : read).map(({ at, stored: value, used }) => {
      const leafPath = `${path}.layoutPreferences${at}`;
      return at.endsWith('.legend') && (used === null || typeof used === 'string')
        ? loadRepair({
            code: 'legend-position-defaulted', path: leafPath, effect: used === null ? 'dropped' : 'defaulted',
            stored: String(value), ...(used === null ? {} : { used })
          })
        : loadRepair({
            code: 'setting-replaced', path: leafPath, effect: used === undefined ? 'dropped' : 'defaulted',
            stored: typeof value === 'string' ? value : JSON.stringify(value) ?? '',
            ...(typeof used === 'string' ? { used } : {})
          });
    });
  }
  if (slotOf) {
    return layout === undefined || layout === null || projected ? []
      : [loadRepair({ code: 'setting-replaced', path: `${path}.layoutPreferences`, effect: 'defaulted', stored: JSON.stringify(layout) ?? '' })];
  }
  return Object.entries(LEGACY_LEGEND_FIELD_MODES)
    .filter(([field, mode]) => storedLegendValue(ui[field]) && legendPositionOfMode(ui[field], mode) === null)
    .map(([field]) => loadRepair({ code: 'legend-position-defaulted', path: `${path}.${field}`, effect: 'dropped', stored: String(ui[field]) }));
};

/**
 * A Session 27-44's `ui` as Load reads its layout, at the file's own path,
 * before Python converts the file (D-15-129): what it does not use as stored
 * (`storedLegendPositionRepairs`), and a stored `ui.layoutPreferences` as
 * Load uses it, so Python receives the read value.
 * @param {Record<string, unknown>} ui
 * @param {string} path Where `ui` is in the file.
 * @returns {{ value: Record<string, unknown>, repairs: LoadRepair[] }}
 */
export const readStoredLayoutUi = (ui, path) => ({
  value: isRecord(ui.layoutPreferences) ? { ...ui, layoutPreferences: normalizeLayoutPreferences(ui.layoutPreferences) } : ui,
  repairs: storedLegendPositionRepairs(ui, path)
});

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
