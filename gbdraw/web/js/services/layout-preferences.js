// @ts-check
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
 * @param {unknown} value
 * @param {string} fallback
 */
const normalizeLegendPosition = (value, fallback) => {
  const normalized = String(value ?? '').trim().toLowerCase();
  return normalized || fallback;
};

/**
 * @param {unknown} value
 * @returns {value is string}
 */
const storedText = (value) => typeof value === 'string' && value.trim() !== '';

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
    legend: normalizeLegendPosition(
      circular.single?.legend,
      WEB_UX_PROFILE.circular.legend
    ),
    plotTitlePosition: normalizeCircularPlotTitlePosition(
      circular.single?.plotTitlePosition
    )
  };
  if (!useMultiRecord) return single;
  return {
    legend: storedText(circular.multi?.legend)
      ? normalizeLegendPosition(circular.multi.legend, single.legend)
      : single.legend,
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
      legend: normalizeLegendPosition(
        preferences?.linear?.legend,
        WEB_UX_PROFILE.linear.legend
      ),
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
      preferences.linear.legend = normalizeLegendPosition(
        patch.legend,
        WEB_UX_PROFILE.linear.legend
      );
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
    preferences.circular[key].legend = normalizeLegendPosition(
      patch.legend,
      WEB_UX_PROFILE.circular.legend
    );
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
  const single = resolveCircularLayoutPreference(/** @type {StoredLayoutPreferences} */ (source), false);
  const rawMulti = /** @type {StoredLayoutPreferences} */ (source).circular?.multi;
  return {
    circular: {
      single,
      multi: {
        legend: storedText(rawMulti?.legend)
          ? normalizeLegendPosition(rawMulti.legend, single.legend)
          : null,
        plotTitlePosition: storedText(rawMulti?.plotTitlePosition)
          ? normalizeCircularPlotTitlePosition(rawMulti.plotTitlePosition)
          : null
      }
    },
    linear: {
      legend: normalizeLegendPosition(
        /** @type {StoredLayoutPreferences} */ (source).linear?.legend,
        WEB_UX_PROFILE.linear.legend
      ),
      plotTitlePosition: normalizeLinearPlotTitlePosition(
        /** @type {StoredLayoutPreferences} */ (source).linear?.plotTitlePosition
      )
    }
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

  const circularLegend = storedText(ui?.circularLegendPosition)
    ? normalizeLegendPosition(ui.circularLegendPosition, activeLegend)
    : storedText(ui?.legend)
      ? normalizeLegendPosition(ui.legend, activeLegend)
      : normalizeLegendPosition(activeLegend, 'left');
  const circularPlotTitle = storedText(ui?.circularPlotTitlePosition)
    ? normalizeCircularPlotTitlePosition(ui.circularPlotTitlePosition)
    : normalizeCircularPlotTitlePosition(activePlotTitlePosition);
  const single = {
    legend: storedText(ui?.circularSingleRecordLegendPosition)
      ? normalizeLegendPosition(ui.circularSingleRecordLegendPosition, circularLegend)
      : multiRecord
        ? circularLegend
        : normalizeLegendPosition(activeLegend, circularLegend),
    plotTitlePosition: storedText(ui?.circularSingleRecordPlotTitlePosition)
      ? normalizeCircularPlotTitlePosition(ui.circularSingleRecordPlotTitlePosition)
      : multiRecord
        ? circularPlotTitle
        : normalizeCircularPlotTitlePosition(activePlotTitlePosition)
  };
  const multi = {
    legend: storedText(ui?.circularMultiRecordLegendPosition)
      ? normalizeLegendPosition(ui.circularMultiRecordLegendPosition, circularLegend)
      : multiRecord
        ? normalizeLegendPosition(activeLegend, circularLegend)
        : circularLegend,
    plotTitlePosition: storedText(ui?.circularMultiRecordPlotTitlePosition)
      ? normalizeCircularPlotTitlePosition(ui.circularMultiRecordPlotTitlePosition)
      : multiRecord
        ? normalizeCircularPlotTitlePosition(activePlotTitlePosition)
        : circularPlotTitle
  };
  return {
    circular: { single, multi },
    linear: {
      legend: storedText(ui?.linearLegendPosition)
        ? normalizeLegendPosition(ui.linearLegendPosition, 'bottom')
        : mode === 'linear'
          ? normalizeLegendPosition(activeLegend, 'bottom')
          : 'bottom',
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

// Partial objects remain authoritative for session compatibility; normalization
// supplies the current defaults for omitted branches.
/** @param {Record<string, unknown> | null | undefined} ui */
const hasStoredLayoutPreferences = (ui) => (
  isRecord(ui?.layoutPreferences) ||
  LEGACY_LAYOUT_PREFERENCE_FIELDS.some((field) => storedText(ui?.[field]))
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
    !storedText(ui?.linearLegendPosition) &&
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
