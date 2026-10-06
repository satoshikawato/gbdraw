// @ts-check
const AUTO_PARALLEL_RUN_LIMIT = 4;

const positiveInteger = (value) => {
  const parsed = Number(value);
  return Number.isInteger(parsed) && parsed >= 1 ? parsed : null;
};

export const getLosatHardwareThreads = () =>
  positiveInteger(globalThis.navigator?.hardwareConcurrency) || 4;

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
