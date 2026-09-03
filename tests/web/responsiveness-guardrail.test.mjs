import assert from 'node:assert/strict';
import { createRequire } from 'node:module';

const require = createRequire(import.meta.url);
const {
  aggregateResponsivenessSamples,
  evaluateResponsivenessGuardrail,
  median
} = require('./helpers/responsiveness-guardrail.cjs');

assert.equal(median([9, 1, 5]), 5);
assert.equal(median([8, 2, 4, 6]), 5);
assert.throws(() => median([]), /at least one/);
assert.throws(() => median([1, Number.NaN]), /finite/);

const aggregate = aggregateResponsivenessSamples([
  {
    longTasks: { totalDurationMs: 20, longestDurationMs: 12 },
    raf: { p95Ms: 22 },
    heartbeat: { maximumGapMs: 55 }
  },
  {
    longTasks: { totalDurationMs: 10, longestDurationMs: 8 },
    raf: { p95Ms: 18 },
    heartbeat: { maximumGapMs: 65 }
  },
  {
    longTasks: { totalDurationMs: 30, longestDurationMs: 16 },
    raf: { p95Ms: 20 },
    heartbeat: { maximumGapMs: 60 }
  }
]);
assert.deepEqual(aggregate, {
  repetitions: 3,
  medianLongTaskTotalMs: 20,
  medianLongestLongTaskMs: 12,
  medianRafP95Ms: 20,
  maximumHeartbeatGapMs: 65
});
assert.deepEqual(evaluateResponsivenessGuardrail(aggregate, {
  maximumMedianLongestLongTaskMs: 15
}), []);
assert.deepEqual(evaluateResponsivenessGuardrail(aggregate, {
  maximumMedianLongestLongTaskMs: 11
}), [
  'median longest Long Task 12 ms exceeded 11 ms'
]);

const slowRunnerAggregate = aggregateResponsivenessSamples([
  {
    longTasks: { totalDurationMs: 2755, longestDurationMs: 196 },
    raf: { p95Ms: 116.7 },
    heartbeat: { maximumGapMs: 200.8 }
  },
  {
    longTasks: { totalDurationMs: 2598, longestDurationMs: 104 },
    raf: { p95Ms: 100 },
    heartbeat: { maximumGapMs: 112.5 }
  },
  {
    longTasks: { totalDurationMs: 2767, longestDurationMs: 115 },
    raf: { p95Ms: 116.6 },
    heartbeat: { maximumGapMs: 122.4 }
  }
]);
assert.deepEqual(evaluateResponsivenessGuardrail(slowRunnerAggregate, {
  maximumMedianLongestLongTaskMs: 250
}), []);

const historicalRegressionAggregate = aggregateResponsivenessSamples([
  {
    longTasks: { totalDurationMs: 2977, longestDurationMs: 905 },
    raf: { p95Ms: 0 },
    heartbeat: { maximumGapMs: 0 }
  },
  {
    longTasks: { totalDurationMs: 2505, longestDurationMs: 964 },
    raf: { p95Ms: 0 },
    heartbeat: { maximumGapMs: 0 }
  },
  {
    longTasks: { totalDurationMs: 3361, longestDurationMs: 959 },
    raf: { p95Ms: 0 },
    heartbeat: { maximumGapMs: 0 }
  }
]);
assert.deepEqual(evaluateResponsivenessGuardrail(historicalRegressionAggregate, {
  maximumMedianLongestLongTaskMs: 250
}), [
  'median longest Long Task 959 ms exceeded 250 ms'
]);

console.log(JSON.stringify({
  status: 'responsiveness guardrail aggregate tests passed',
  assertions: 10
}));
