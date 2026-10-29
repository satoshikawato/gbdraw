// The search for a minimal reproduction of a failing random-walk step
// (tests/web/helpers/live-random-walk.mjs). A pure module: the walk passes a
// `reproduces` predicate that replays a candidate in the browser, and the node
// test tests/web/live-random-walk-shrink.test.mjs drives it with a stub.
//
// A candidate is a list of step indices in walk order that ends with the
// failing step. The search tries, in order: the failing step alone; each
// earlier step with it, nearest first, until those replays reach half the
// budget (`pairPhaseBudget`), so ddmin always gets the other half; then ddmin
// (Zeller and Hildebrandt 2002) over the steps before it, with the failing
// step always appended. The whole walk is taken to reproduce (it failed in the
// walk itself) and is not replayed. Each predicate call on a new candidate is
// a replay; a candidate tried before costs none. The search stops at `budget`
// replays and returns the shortest reproducing candidate found so far.

/**
 * @typedef {object} ShrinkResult
 * @property {number[]} indices the shortest candidate that reproduced, or the whole walk
 * @property {number} replays the predicate calls made
 * @property {boolean} finished false when the search stopped at the budget
 * @property {boolean} confirmed true when `indices` reproduced in a replay (the whole walk is never replayed)
 */

// The parts of `items` cut into `count` runs as even as possible.
/** @param {number[]} items @param {number} count */
const split = (items, count) => {
  const parts = [];
  let start = 0;
  for (let part = 0; part < count; part += 1) {
    const end = start + Math.floor((items.length - start) / (count - part));
    parts.push(items.slice(start, end));
    start = end;
  }
  return parts;
};

const OUT_OF_BUDGET = Symbol('out of budget');

/** The replays the failing step alone and the pairs may use: half the budget, rounded up. @param {number} budget */
export const pairPhaseBudget = (budget) => Math.ceil(budget / 2);

/**
 * @param {{ prefix: number[], failing: number, budget: number, reproduces: (indices: number[]) => Promise<boolean> }} options
 *   `prefix`: the indices of the steps before the failing one, in walk order.
 * @returns {Promise<ShrinkResult>}
 */
export const shrinkWalk = async ({ prefix, failing, budget, reproduces }) => {
  const tried = new Map();
  let replays = 0;
  let best = [...prefix, failing];
  let confirmed = false;
  /** @param {number[]} subset a subset of `prefix`, in walk order */
  const test = async (subset) => {
    const key = subset.join(',');
    if (tried.has(key)) return tried.get(key);
    if (replays >= budget) throw OUT_OF_BUDGET;
    replays += 1;
    const result = Boolean(await reproduces([...subset, failing]));
    tried.set(key, result);
    if (result && (!confirmed || subset.length + 1 < best.length)) {
      best = [...subset, failing];
      confirmed = true;
    }
    return result;
  };
  const done = (finished) => ({ indices: best, replays, finished, confirmed });
  try {
    if (await test([])) return done(true);
    for (const index of [...prefix].reverse()) {
      if (replays >= pairPhaseBudget(budget)) break;
      if (await test([index])) return done(true);
    }
    // ddmin over the prefix: the whole prefix reproduces (the walk failed).
    let current = prefix;
    let granularity = 2;
    while (current.length >= 2) {
      const parts = split(current, granularity);
      let reduced = false;
      for (const part of parts) {
        if (await test(part)) {
          current = part;
          granularity = 2;
          reduced = true;
          break;
        }
      }
      if (!reduced) {
        for (const part of parts) {
          const complement = current.filter((index) => !part.includes(index));
          if (await test(complement)) {
            current = complement;
            granularity = Math.max(granularity - 1, 2);
            reduced = true;
            break;
          }
        }
      }
      if (!reduced) {
        if (granularity >= current.length) break;
        granularity = Math.min(granularity * 2, current.length);
      }
    }
    if (!confirmed) best = [...current, failing];
    return done(true);
  } catch (error) {
    if (error === OUT_OF_BUDGET) return done(false);
    throw error;
  }
};
