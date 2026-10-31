// The minimal-reproduction search of the random walk
// (tests/web/helpers/live-random-walk-shrink.mjs), over a stubbed predicate
// instead of browser replays.
import assert from 'node:assert/strict';
import test from 'node:test';
import { pairPhaseBudget, shrinkWalk } from './helpers/live-random-walk-shrink.mjs';

const range = (from, to) => Array.from({ length: to - from + 1 }, (_, offset) => from + offset);
// A predicate that reproduces when the candidate holds every step of `needed`.
const needs = (needed, calls = []) => async (indices) => {
  calls.push(indices);
  return needed.every((index) => indices.includes(index));
};

test('the failing step alone is tried first', async () => {
  const calls = [];
  const result = await shrinkWalk({ prefix: range(1, 9), failing: 10, budget: 40, reproduces: needs([10], calls) });
  assert.deepEqual(result, { indices: [10], replays: 1, finished: true, confirmed: true });
  assert.deepEqual(calls, [[10]]);
});

test('then each earlier step with it, nearest first', async () => {
  const calls = [];
  const result = await shrinkWalk({ prefix: range(1, 12), failing: 13, budget: 40, reproduces: needs([2, 13], calls) });
  assert.deepEqual(result.indices, [2, 13]);
  assert.equal(result.finished, true);
  assert.deepEqual(calls, [[13], ...range(2, 12).reverse().map((index) => [index, 13])]);
});

test('then ddmin over the prefix finds a three-step interaction', async () => {
  const result = await shrinkWalk({ prefix: range(1, 19), failing: 20, budget: 200, reproduces: needs([3, 11, 20]) });
  assert.deepEqual(result.indices, [3, 11, 20]);
  assert.equal(result.finished, true);
  assert.equal(result.confirmed, true);
});

test('a candidate is replayed at most once', async () => {
  const calls = [];
  await shrinkWalk({ prefix: range(1, 15), failing: 16, budget: 500, reproduces: needs([4, 9, 14, 16], calls) });
  const keys = calls.map((indices) => indices.join(','));
  assert.equal(new Set(keys).size, keys.length);
  for (const indices of calls) {
    assert.equal(indices.at(-1), 16, 'every candidate ends with the failing step');
    assert.deepEqual(indices, [...indices].sort((left, right) => left - right), 'candidates keep walk order');
  }
});

test('the budget stops the search with the best candidate so far', async () => {
  const unbounded = await shrinkWalk({ prefix: range(1, 29), failing: 30, budget: 1000, reproduces: needs([5, 17, 30]) });
  assert.ok(unbounded.replays > 30);
  // 20 replays for the step alone and the pairs (half the budget), then 20 for ddmin.
  const bounded = await shrinkWalk({ prefix: range(1, 29), failing: 30, budget: 40, reproduces: needs([5, 17, 30]) });
  assert.equal(bounded.replays, 40);
  assert.equal(bounded.finished, false);
  assert.ok(bounded.indices.includes(5) && bounded.indices.includes(17) && bounded.indices.at(-1) === 30);
  assert.ok(bounded.indices.length < 30, 'ddmin had cut the prefix before the budget ran out');
});

test('a budget of 0 replays nothing and keeps the whole walk', async () => {
  const calls = [];
  const result = await shrinkWalk({ prefix: [1, 2, 3], failing: 4, budget: 0, reproduces: needs([4], calls) });
  assert.deepEqual(result, { indices: [1, 2, 3, 4], replays: 0, finished: false, confirmed: false });
  assert.deepEqual(calls, []);
});

test('when no replay reproduces, the whole walk stays the reproduction, unconfirmed', async () => {
  const result = await shrinkWalk({ prefix: range(1, 7), failing: 8, budget: 100, reproduces: async () => false });
  assert.deepEqual(result.indices, range(1, 8));
  assert.equal(result.finished, true);
  assert.equal(result.confirmed, false);
});

test('the first step of a walk has no prefix to cut', async () => {
  const result = await shrinkWalk({ prefix: [], failing: 1, budget: 40, reproduces: async () => false });
  assert.deepEqual(result, { indices: [1], replays: 1, finished: true, confirmed: false });
});

test('the step alone and the pairs use at most half the budget; ddmin gets the rest', async () => {
  assert.deepEqual([0, 1, 2, 9, 10, 40].map(pairPhaseBudget), [0, 1, 1, 5, 5, 20]);
  const calls = [];
  // The one step needed lies beyond the pairs the cap allows, so ddmin finds it.
  const result = await shrinkWalk({ prefix: range(1, 29), failing: 30, budget: 10, reproduces: needs([3, 30], calls) });
  assert.deepEqual(calls.slice(0, 5), [[30], [29, 30], [28, 30], [27, 30], [26, 30]]);
  // ddmin starts with the first half of the prefix: subsets before complements.
  assert.deepEqual(calls[5], [...range(1, 14), 30]);
  assert.equal(result.replays, 10);
  assert.ok(result.confirmed && result.indices.includes(3) && result.indices.length <= 3, `${result.indices}`);
});

test('a candidate tried before costs no replay', async () => {
  // At granularity 2 the complements of the two halves are the halves again.
  const calls = [];
  const result = await shrinkWalk({ prefix: range(1, 7), failing: 8, budget: 100, reproduces: needs([2, 5, 8], calls) });
  assert.deepEqual(result.indices, [2, 5, 8]);
  assert.equal(result.replays, calls.length);
  assert.equal(result.replays, 19);
  // With a budget of exactly those replays the search still finishes.
  const exact = await shrinkWalk({ prefix: range(1, 7), failing: 8, budget: 19, reproduces: needs([2, 5, 8]) });
  assert.deepEqual(exact, { indices: [2, 5, 8], replays: 19, finished: true, confirmed: true });
});

test('ddmin follows Zeller: subsets, then complements with granularity n - 1, then 2n', async () => {
  const calls = [];
  const result = await shrinkWalk({ prefix: range(1, 8), failing: 9, budget: 100, reproduces: needs([1, 3, 5, 7, 9], calls) });
  assert.deepEqual(result, { indices: [1, 3, 5, 7, 9], replays: 36, finished: true, confirmed: true });
  // After the step alone and the 8 pairs (none reproduces).
  assert.deepEqual(calls.slice(9).map((indices) => indices.join(',')), [
    // n = 2: subsets (their complements are the same candidates)
    '1,2,3,4,9', '5,6,7,8,9',
    // n = 4: subsets, then complements
    '1,2,9', '3,4,9', '5,6,9', '7,8,9', '3,4,5,6,7,8,9', '1,2,5,6,7,8,9', '1,2,3,4,7,8,9', '1,2,3,4,5,6,9',
    // n = 8: the complement without 2 reproduces, so n = 7 on 7 steps
    '2,3,4,5,6,7,8,9', '1,3,4,5,6,7,8,9',
    // n = 7: the complement without 4 reproduces, so n = 6
    '1,4,5,6,7,8,9', '1,3,5,6,7,8,9',
    // n = 6: without 6, so n = 5
    '3,5,6,7,8,9', '1,5,6,7,8,9', '1,3,6,7,8,9', '1,3,5,7,8,9',
    // n = 5: without 8, so n = 4
    '3,5,7,8,9', '1,5,7,8,9', '1,3,7,8,9', '1,3,5,8,9', '1,3,5,7,9',
    // n = 4: no complement reproduces, and n is the size: done
    '3,5,7,9', '1,5,7,9', '1,3,7,9', '1,3,5,9'
  ]);
});
