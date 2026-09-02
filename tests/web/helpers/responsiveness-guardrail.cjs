const median = (values) => {
  if (!Array.isArray(values) || values.length === 0) {
    throw new TypeError('median requires at least one numeric value.');
  }
  if (values.some((value) => !Number.isFinite(value))) {
    throw new TypeError('median accepts only finite numeric values.');
  }
  const sorted = [...values].sort((left, right) => left - right);
  const middle = Math.floor(sorted.length / 2);
  return sorted.length % 2
    ? sorted[middle]
    : (sorted[middle - 1] + sorted[middle]) / 2;
};

const aggregateResponsivenessSamples = (samples) => {
  if (!Array.isArray(samples) || samples.length === 0) {
    throw new TypeError('At least one responsiveness sample is required.');
  }
  return Object.freeze({
    repetitions: samples.length,
    medianLongTaskTotalMs: median(samples.map(({ longTasks }) => longTasks.totalDurationMs)),
    medianLongestLongTaskMs: median(
      samples.map(({ longTasks }) => longTasks.longestDurationMs)
    ),
    medianRafP95Ms: median(samples.map(({ raf }) => raf.p95Ms)),
    maximumHeartbeatGapMs: Math.max(
      ...samples.map(({ heartbeat }) => heartbeat.maximumGapMs)
    )
  });
};

const evaluateResponsivenessGuardrail = (aggregate, thresholds) => Object.freeze([
  aggregate.medianLongTaskTotalMs > thresholds.maximumMedianLongTaskTotalMs
    ? `median Long Task total ${aggregate.medianLongTaskTotalMs} ms exceeded `
      + `${thresholds.maximumMedianLongTaskTotalMs} ms`
    : null,
  aggregate.medianLongestLongTaskMs > thresholds.maximumMedianLongestLongTaskMs
    ? `median longest Long Task ${aggregate.medianLongestLongTaskMs} ms exceeded `
      + `${thresholds.maximumMedianLongestLongTaskMs} ms`
    : null
].filter(Boolean));

module.exports = {
  aggregateResponsivenessSamples,
  evaluateResponsivenessGuardrail,
  median
};
