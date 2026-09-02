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
  maximumMedianLongTaskTotalMs: 25,
  maximumMedianLongestLongTaskMs: 15
}), []);
assert.deepEqual(evaluateResponsivenessGuardrail(aggregate, {
  maximumMedianLongTaskTotalMs: 19,
  maximumMedianLongestLongTaskMs: 11
}), [
  'median Long Task total 20 ms exceeded 19 ms',
  'median longest Long Task 12 ms exceeded 11 ms'
]);

console.log(JSON.stringify({
  status: 'responsiveness guardrail aggregate tests passed',
  assertions: 8
}));
