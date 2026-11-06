// @ts-check
const AUTO_PARALLEL_RUN_LIMIT = 4;

/**
 * @param {unknown} value
 * @returns {number | null}
 */
const positiveInteger = (value) => {
  const parsed = Number(value);
  return Number.isInteger(parsed) && parsed >= 1 ? parsed : null;
};

export const getLosatHardwareThreads = () =>
  positiveInteger(globalThis.navigator?.hardwareConcurrency) || 4;

/**
 * How LOSAT runs: one app-level setting (`state.losatExecution`, Session `ui.losatExecution`) for both drawings, as
 * saved: 'auto', 'serial' or 'threaded' for the mode, 'safe', 'available' or a count for the total thread budget, and
 * 'auto' (or blank) or a count for the threads per job and the pair workers.
 * @typedef {{ executionMode?: unknown, totalThreadBudget?: unknown, threadsPerJob?: unknown, parallelWorkers?: unknown }} LosatExecutionSettings
 */

/** @param {unknown} value @param {string} fallback The saved word of a setting, in lower case. */
const settingWord = (value, fallback) => String(value || fallback).trim().toLowerCase();

// What a run asks LOSAT for, read from the saved settings: the count or mode the user chose, or undefined where the
// runner plans it (`resolveLosatThreadPlan`).
/** @param {LosatExecutionSettings} execution */
export const requestedLosatParallelWorkers = (execution) => positiveInteger(settingWord(execution.parallelWorkers, 'auto')) ?? undefined;

/** @param {LosatExecutionSettings} execution */
export const requestedLosatExecutionMode = (execution) => {
  const mode = settingWord(execution.executionMode, 'auto');
  return ['auto', 'serial', 'threaded'].includes(mode) ? mode : 'auto';
};

/** blastp runs the threads per job the user chose; another program runs one. @param {LosatExecutionSettings} execution @param {unknown} program */
export const requestedLosatThreadsPerJob = (execution, program) => (
  program === 'blastp' ? positiveInteger(settingWord(execution.threadsPerJob, 'auto')) ?? undefined : 1
);

/** 'available' is every hardware thread and a count is at most that; 'safe' is the runner's plan. @param {LosatExecutionSettings} execution */
export const requestedLosatTotalThreadBudget = (execution) => {
  const budget = settingWord(execution.totalThreadBudget, 'safe');
  if (budget === 'available') return getLosatHardwareThreads();
  const count = positiveInteger(budget);
  return count === null ? undefined : Math.min(count, getLosatHardwareThreads());
};

/**
 * The requested settings. The three request values are the saved Settings
 * values: 'safe' or 'available' or a count for the total budget, and a count or
 * blank (Auto) for the other two.
 * @typedef {object} LosatThreadPlanRequest
 * @property {number} [jobCount]
 * @property {number} [hardwareThreads]
 * @property {string | number | null} [totalThreadBudget]
 * @property {string | number | null} [threadsPerJob]
 * @property {string | number | null} [parallelWorkers]
 */

/**
 * The resolved plan. `totalBudget` resolves the request's `totalThreadBudget`.
 * @typedef {object} LosatThreadPlan
 * @property {number} totalBudget
 * @property {number} threadsPerJob
 * @property {number} pairWorkers
 * @property {number} autoPairWorkers
 * @property {number} maxPairWorkers
 */

// Resolve requested settings without replacing the saved Auto or manual values.
// The controls use the estimated source-job count; execution uses pending jobs.
/**
 * @param {LosatThreadPlanRequest} [request]
 * @returns {LosatThreadPlan}
 */
export const resolveLosatThreadPlan = ({
  jobCount = 0,
  hardwareThreads = getLosatHardwareThreads(),
  totalThreadBudget,
  threadsPerJob,
  parallelWorkers
} = {}) => {
  const hardware = positiveInteger(hardwareThreads) || 4;
  const budgetMode = String(totalThreadBudget ?? 'safe').trim().toLowerCase();
  const requestedBudget = budgetMode === 'available' ? hardware : positiveInteger(totalThreadBudget);
  const totalBudget = Math.min(hardware, requestedBudget || Math.max(1, Math.floor(hardware / 2)));
  const jobs = positiveInteger(jobCount) || 0;
  const requestedThreads = positiveInteger(threadsPerJob);
  const fixedThreads = requestedThreads === null ? null : Math.min(requestedThreads, totalBudget);
  const maxPairWorkers = Math.min(jobs, Math.floor(totalBudget / (fixedThreads || 1)));
  const autoPairWorkers = ['safe', 'auto', ''].includes(budgetMode)
    ? Math.min(AUTO_PARALLEL_RUN_LIMIT, maxPairWorkers)
    : maxPairWorkers;
  const pairWorkers = Math.min(maxPairWorkers, positiveInteger(parallelWorkers) || autoPairWorkers);
  const effectiveThreads = fixedThreads || (pairWorkers > 0 ? Math.floor(totalBudget / pairWorkers) : 1);
  return { totalBudget, threadsPerJob: effectiveThreads, pairWorkers, autoPairWorkers, maxPairWorkers };
};
