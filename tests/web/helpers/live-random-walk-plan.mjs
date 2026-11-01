// The kind schedule of the live-vs-Generate random walk (gbdraw/web/CLAUDE.md
// R3; tests/web/helpers/live-random-walk.mjs runs it). A pure module, so the
// node test tests/web/live-random-walk-plan.test.mjs checks the coverage
// without a browser.
//
// A plan is a list of walks, each on one fixture, each a list of planned steps
// `{ kind, popup }`. `planWalk` builds it greedily from a seed until
// every coverage target of its level is covered:
//
// - level `pr`: every fixture; every feasible ordered pair of kinds (A, B), A = B included,
//   with A before B in one walk and no barrier step (Reset Settings, Session
//   Load: they change every drawing) strictly between them; every kind in each
//   popup state; and the steps the caller requires (fixture, kind, popup).
// - level `promotion`: the `pr` targets, every feasible ordered pair (A, B) as
//   adjacent steps of one walk, and every fixture x kind x popup.
//
// At each position the planner takes the step that covers the most uncovered
// targets, then the one that leaves the most uncovered pairs to start, then a
// seeded draw. The step arguments (which feature, color, index) are drawn
// later by the walk from its own generator, so runs with different seeds
// still differ while the coverage is fixed.

// The step kinds.
export const WALK_KINDS = Object.freeze([
  'fill', 'stroke', 'featureVisibility', 'labelText', 'labelVisibility',
  'visibilityRuleAdd', 'visibilityRuleEdit', 'visibilityRuleRemove', 'colorRuleAdd', 'colorRuleEdit', 'colorRuleRemove',
  'legendRename', 'legendColor', 'legendSort', 'legendSortDefault', 'legendDelete', 'legendRestore', 'legendAdd',
  'trackData', 'palette', 'undo', 'undoRedo', 'resultSwitch', 'reset', 'load', 'drawingSwitch'
]);

// The kinds that change every drawing by design (PD-OI-090 Reset, Save and
// Load): a pair across one of them is not an interaction.
export const BARRIER_KINDS = Object.freeze(['reset', 'load']);

// The fixtures: the shown drawing with its Result, and the other drawing's
// Result (Linear always draws the one-record file). The name is the fixture's
// identity.
export const FIXTURES = Object.freeze([
  { name: 'Circular, one Result', mode: 'circular', results: 'single' },
  { name: 'Linear, one Result', mode: 'linear', results: 'single' },
  { name: 'Circular, two-Result batch', mode: 'circular', results: 'batch' }
]);

/** @typedef {'pr' | 'promotion'} PlanLevel */
/** @typedef {{ kind: string, popup: boolean }} PlannedStep */
/** @typedef {{ fixtureIndex: number, steps: PlannedStep[] }} PlannedWalk */
/** @typedef {{ walks: PlannedWalk[] }} WalkPlan */
/** A step a caller requires: the fixture by name. */
/** @typedef {{ fixture: string, kind: string, popup: boolean }} RequiredStep */

export const DEFAULT_WALK_LENGTH = Object.freeze({ pr: 40, promotion: 50 });

// The pull-request walk (tests/web/live-generate-random-walk.playwright.spec.js):
// a fixed seed, and the case-matrix cases of
// live-generate-parity.playwright.spec.js that the walk replaced (R3), each a
// required step in the fixture and popup state the case
// pinned. Dropping one from the plan fails the walk spec and the node test:
// restore the case, or keep a plan that draws it.
export const PR_SEED = 20261018;
/** @type {readonly (RequiredStep & { edit: string })[]} */
export const REPLACED_MATRIX_CASES = Object.freeze([
  { edit: 'Specific color rule add', kind: 'colorRuleAdd', fixture: 'Linear, one Result', popup: true }
]);
// The most walks a plan may hold; a plan that reaches it leaves targets uncovered.
const MAX_WALKS = 60;
// The greedy plans tried per seed; the plan with the fewest missing targets,
// then walks, then steps wins. A short plan is worth more tries.
const ATTEMPTS = Object.freeze({ pr: 16, promotion: 1 });

// Only the batch has two Results to switch between. Reset clears its batch
// settings and its drawing switch would leave it, so neither runs there.
/** @param {string} kind @param {{ results: string }} fixture */
export const allowedIn = (kind, fixture) => {
  if (kind === 'resultSwitch') return fixture.results === 'batch';
  if (kind === 'reset' || kind === 'drawingSwitch') return fixture.results !== 'batch';
  return true;
};

// mulberry32: a 32-bit seeded generator; the same seed gives the same draws.
/** @param {number} seed */
export const createRandom = (seed) => {
  let state = seed >>> 0;
  const next = () => {
    state = (state + 0x6D2B79F5) >>> 0;
    let value = state;
    value = Math.imul(value ^ (value >>> 15), value | 1);
    value ^= value + Math.imul(value ^ (value >>> 7), value | 61);
    return ((value ^ (value >>> 14)) >>> 0) / 4294967296;
  };
  /** @param {number} count */
  const int = (count) => Math.floor(next() * count);
  /** @template T @param {readonly T[]} items @returns {T} */
  const pick = (items) => items[int(items.length)];
  return { next, int, pick };
};

const POPUP_STATES = Object.freeze([false, true]);
// The pair keys, built once: the planner reads them in its inner loop.
const PAIR_KEYS = new Map(WALK_KINDS.map((first) => [first, new Map(WALK_KINDS.map((second) => [second, `${first}>${second}`]))]));
const pairKey = (first, second) => PAIR_KEYS.get(first).get(second);
const flagKey = (kind, popup) => `${kind}|${popup ? 'popup' : 'none'}`;
const fixtureFlagKey = (fixtureIndex, kind, popup) => `${FIXTURES[fixtureIndex].name}|${flagKey(kind, popup)}`;
const isBarrier = (kind) => BARRIER_KINDS.includes(kind);
const KINDS_IN = FIXTURES.map((fixture) => WALK_KINDS.filter((kind) => allowedIn(kind, fixture)));
const kindsIn = (fixtureIndex) => KINDS_IN[fixtureIndex];
const pairFeasibleIn = (fixtureIndex, first, second) => allowedIn(first, FIXTURES[fixtureIndex]) && allowedIn(second, FIXTURES[fixtureIndex]);

/**
 * The coverage targets of a level, in groups: `fixtures` (a walk on each
 * fixture), `sequence` (ordered pairs as a
 * subsequence within a barrier-free stretch), `flags` (kind x popup),
 * `required`, and for `promotion` also `adjacent` (ordered pairs as adjacent
 * steps) and `fixtureFlags` (fixture x kind x popup).
 * @param {PlanLevel} level
 * @param {readonly RequiredStep[]} [required]
 * @returns {Record<string, Set<string>>}
 */
export const coverageTargets = (level, required = []) => {
  if (level !== 'pr' && level !== 'promotion') throw new Error(`unknown plan level ${level}`);
  const pairs = new Set();
  for (const first of WALK_KINDS) {
    for (const second of WALK_KINDS) {
      if (FIXTURES.some((_, fixtureIndex) => pairFeasibleIn(fixtureIndex, first, second))) pairs.add(pairKey(first, second));
    }
  }
  const flags = new Set(WALK_KINDS.flatMap((kind) => POPUP_STATES.map((popup) => flagKey(kind, popup))));
  const requiredKeys = new Set(required.map(({ fixture, kind, popup }) => {
    const fixtureIndex = FIXTURES.findIndex(({ name }) => name === fixture);
    if (fixtureIndex < 0) throw new Error(`no fixture ${fixture}`);
    if (!WALK_KINDS.includes(kind) || !allowedIn(kind, FIXTURES[fixtureIndex])) throw new Error(`${kind} cannot run in ${fixture}`);
    return fixtureFlagKey(fixtureIndex, kind, popup);
  }));
  const fixtures = new Set(FIXTURES.map(({ name }) => name));
  const targets = { fixtures, sequence: pairs, flags, required: requiredKeys };
  if (level === 'pr') return targets;
  const fixtureFlags = new Set(FIXTURES.flatMap((_, fixtureIndex) => kindsIn(fixtureIndex)
    .flatMap((kind) => POPUP_STATES.map((popup) => fixtureFlagKey(fixtureIndex, kind, popup)))));
  return { ...targets, adjacent: new Set(pairs), fixtureFlags };
};

// The targets a step at the end of a walk covers, by group. `segment` holds
// the kinds since the last barrier (the barrier included); `previous` is the
// kind before the step, or null at the start of a walk.
const coveredBy = (fixtureIndex, segment, previous, { kind, popup }) => ({
  fixtures: [FIXTURES[fixtureIndex].name],
  sequence: [...segment].map((first) => pairKey(first, kind)),
  flags: [flagKey(kind, popup)],
  required: [fixtureFlagKey(fixtureIndex, kind, popup)],
  adjacent: previous === null ? [] : [pairKey(previous, kind)],
  fixtureFlags: [fixtureFlagKey(fixtureIndex, kind, popup)]
});
// The segment after a step: a barrier starts a new one.
const nextSegment = (segment, kind) => (isBarrier(kind) ? new Set([kind]) : new Set([...segment, kind]));

/**
 * Which targets of a level a plan covers and which it misses.
 * @param {WalkPlan} plan
 * @param {PlanLevel} level
 * @param {{ required?: readonly RequiredStep[] }} [options]
 * @returns {{ covered: Record<string, number>, missing: string[] }}
 */
export const coverageOf = (plan, level, { required = [] } = {}) => {
  const targets = coverageTargets(level, required);
  const hit = Object.fromEntries(Object.keys(targets).map((group) => [group, new Set()]));
  for (const { fixtureIndex, steps } of plan.walks) {
    let segment = new Set();
    let previous = null;
    for (const step of steps) {
      if (!allowedIn(step.kind, FIXTURES[fixtureIndex])) throw new Error(`${step.kind} cannot run in ${FIXTURES[fixtureIndex].name}`);
      for (const [group, keys] of Object.entries(coveredBy(fixtureIndex, segment, previous, step))) {
        if (hit[group]) for (const key of keys) hit[group].add(key);
      }
      segment = nextSegment(segment, step.kind);
      previous = step.kind;
    }
  }
  const covered = {};
  const missing = [];
  for (const [group, keys] of Object.entries(targets)) {
    covered[group] = [...keys].filter((key) => hit[group].has(key)).length;
    for (const key of keys) if (!hit[group].has(key)) missing.push(`${group} ${key}`);
  }
  return { covered, missing };
};

/**
 * The walk plan of a seed: walks of at most `walkLength` steps until every
 * target of `level` is covered (or MAX_WALKS is reached; `coverageOf` then
 * names what is missing).
 * @param {number} seed
 * @param {{ level: PlanLevel, walkLength?: number, required?: readonly RequiredStep[] }} options
 * @returns {WalkPlan}
 */
export const planWalk = (seed, { level, walkLength = DEFAULT_WALK_LENGTH[level], required = [] }) => {
  if (!Number.isInteger(walkLength) || walkLength < 2) throw new Error(`a walk needs at least 2 steps; got ${walkLength}`);
  let best = null;
  let bestRank = null;
  for (let attempt = 0; attempt < ATTEMPTS[level]; attempt += 1) {
    const plan = greedyPlan(createRandom((seed ^ 0x5EED5EED) + Math.imul(attempt, 0x9E3779B1)), { level, walkLength, required });
    const rank = [coverageOf(plan, level, { required }).missing.length, plan.walks.length, stepCount(plan)];
    const differs = bestRank ? rank.findIndex((value, index) => value !== bestRank[index]) : -1;
    if (!bestRank || (differs >= 0 && rank[differs] < bestRank[differs])) {
      best = plan;
      bestRank = rank;
    }
  }
  return best;
};

/** The number of steps of a plan. @param {WalkPlan} plan */
export const stepCount = (plan) => plan.walks.reduce((sum, { steps }) => sum + steps.length, 0);

// One greedy plan from a generator.
/**
 * @param {ReturnType<typeof createRandom>} random
 * @param {{ level: PlanLevel, walkLength: number, required: readonly RequiredStep[] }} options
 * @returns {WalkPlan}
 */
const greedyPlan = (random, { level, walkLength, required }) => {
  const uncovered = coverageTargets(level, required);
  const remaining = () => Object.values(uncovered).reduce((sum, keys) => sum + keys.size, 0);
  const gainOf = (keysByGroup) => Object.entries(keysByGroup)
    .reduce((sum, [group, keys]) => sum + (uncovered[group] ? keys.filter((key) => uncovered[group].has(key)).length : 0), 0);
  // The pairs a kind could still start in a fixture: they break ties toward
  // the kinds whose later steps cover more.
  const openPairsFrom = (fixtureIndex, kind) => kindsIn(fixtureIndex).reduce((sum, second) => {
    const key = pairKey(kind, second);
    return sum + (uncovered.sequence.has(key) ? 1 : 0) + (uncovered.adjacent?.has(key) ? 1 : 0);
  }, 0);
  // The uncovered targets a fixture can cover.
  const openIn = (fixtureIndex) => {
    const kinds = kindsIn(fixtureIndex);
    let count = 0;
    for (const group of ['sequence', 'adjacent']) {
      for (const key of uncovered[group] || []) {
        const [first, second] = key.split('>');
        if (kinds.includes(first) && kinds.includes(second)) count += 1;
      }
    }
    for (const key of uncovered.flags) if (kinds.includes(key.split('|')[0])) count += 1;
    if (uncovered.fixtures.has(FIXTURES[fixtureIndex].name)) count += 1;
    for (const group of ['required', 'fixtureFlags']) {
      for (const key of uncovered[group] || []) if (key.startsWith(`${FIXTURES[fixtureIndex].name}|`)) count += 1;
    }
    return count;
  };

  /** @type {PlannedWalk[]} */
  const walks = [];
  while (remaining() > 0 && walks.length < MAX_WALKS) {
    const open = FIXTURES.map((_, fixtureIndex) => openIn(fixtureIndex));
    const most = Math.max(...open);
    if (most === 0) break;
    const fixtureIndex = random.pick(open.flatMap((count, index) => (count === most ? [index] : [])));
    const kinds = kindsIn(fixtureIndex);
    /** @type {PlannedStep[]} */
    const steps = [];
    let segment = new Set();
    let previous = null;
    while (steps.length < walkLength && remaining() > 0) {
      const bestOf = (candidates) => {
        let best = [];
        let bestScore = [-1, -1];
        for (const kind of candidates) {
          const potential = openPairsFrom(fixtureIndex, kind);
          // The pairs a kind covers do not depend on its popup state.
          const pairGain = gainOf({
            sequence: [...segment].map((first) => pairKey(first, kind)),
            adjacent: previous === null ? [] : [pairKey(previous, kind)]
          });
          for (const popup of POPUP_STATES) {
            const step = { kind, popup };
            const gain = pairGain + gainOf({
              flags: [flagKey(kind, popup)],
              required: [fixtureFlagKey(fixtureIndex, kind, popup)],
              fixtureFlags: [fixtureFlagKey(fixtureIndex, kind, popup)]
            });
            if (gain > bestScore[0] || (gain === bestScore[0] && potential > bestScore[1])) {
              best = [step];
              bestScore = [gain, potential];
            } else if (gain === bestScore[0] && potential === bestScore[1]) {
              best.push(step);
            }
          }
        }
        return { best, gain: bestScore[0] };
      };
      // A barrier ends the pairs of the steps before it, so it goes at the
      // start or the end of a walk, or where no other step covers anything.
      const atEdge = steps.length === 0 || steps.length === walkLength - 1;
      let { best, gain } = bestOf(atEdge ? kinds : kinds.filter((kind) => !isBarrier(kind)));
      if (!atEdge && gain === 0) ({ best } = bestOf(kinds));
      const step = random.pick(best);
      for (const [group, keys] of Object.entries(coveredBy(fixtureIndex, segment, previous, step))) {
        if (uncovered[group]) for (const key of keys) uncovered[group].delete(key);
      }
      steps.push(step);
      segment = nextSegment(segment, step.kind);
      previous = step.kind;
    }
    walks.push({ fixtureIndex, steps });
  }
  return { walks };
};
