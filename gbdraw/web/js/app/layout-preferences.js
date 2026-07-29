import {
  normalizeCircularPlotTitlePosition,
  normalizeLinearPlotTitlePosition
} from './plot-title-position.js';

const normalizeLegendPosition = (value, fallback) => {
  const normalized = String(value ?? '').trim().toLowerCase();
  return normalized || fallback;
};

const storedText = (value) => typeof value === 'string' && value.trim() !== '';

export const createDefaultLayoutPreferences = () => ({
  circular: {
    single: {
      legend: 'left',
      plotTitlePosition: 'none'
    },
    multi: {
      legend: null,
      plotTitlePosition: null
    }
  },
  linear: {
    legend: 'bottom',
    plotTitlePosition: 'bottom'
  }
});

export const resolveCircularLayoutPreference = (preferences, useMultiRecord = false) => {
  const defaults = createDefaultLayoutPreferences();
  const circular = preferences?.circular || defaults.circular;
  const single = {
    legend: normalizeLegendPosition(circular.single?.legend, 'left'),
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
      legend: normalizeLegendPosition(preferences?.linear?.legend, 'bottom'),
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
      preferences.linear.legend = normalizeLegendPosition(patch.legend, 'bottom');
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
    preferences.circular[key].legend = normalizeLegendPosition(patch.legend, 'left');
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
      legend: normalizeLegendPosition(source.linear?.legend, 'bottom'),
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
    activeLegend = mode === 'linear' ? 'bottom' : 'left',
    activePlotTitlePosition = mode === 'linear' ? 'bottom' : 'none'
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
