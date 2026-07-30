import {
  normalizeCircularPlotTitlePosition,
  normalizeLinearPlotTitlePosition
} from './plot-title-position.js';
import { WEB_UX_PROFILE } from '../web-ux-profile.js';

const normalizeLegendPosition = (value, fallback) => {
  const normalized = String(value ?? '').trim().toLowerCase();
  return normalized || fallback;
};

const storedText = (value) => typeof value === 'string' && value.trim() !== '';

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

export const normalizeLayoutPreferences = (source) => {
  const defaults = createDefaultLayoutPreferences();
  if (!source || typeof source !== 'object' || Array.isArray(source)) {
    return defaults;
  }
  const single = resolveCircularLayoutPreference(source, false);
  const rawMulti = source.circular?.multi;
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
        source.linear?.legend,
        WEB_UX_PROFILE.linear.legend
      ),
      plotTitlePosition: normalizeLinearPlotTitlePosition(
        source.linear?.plotTitlePosition
      )
    }
  };
};

export const replaceLayoutPreferences = (target, source) => {
  const normalized = normalizeLayoutPreferences(source);
  Object.assign(target.circular.single, normalized.circular.single);
  Object.assign(target.circular.multi, normalized.circular.multi);
  Object.assign(target.linear, normalized.linear);
};

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
