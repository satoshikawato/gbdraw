// @ts-check
import { normalizeCollinearSearchScope } from '../services/losat-normalization.js';
import { buildLosatJobSpecs } from './linear-comparisons.js';
import { losatRecordGencode, planLosatSourceJobs } from './linear-sources.js';
import { getLosatHardwareThreads, resolveLosatThreadPlan } from '../services/losat-thread-plan.js';
import { losatThreadingPrecondition } from '../services/losat.js';

const { computed, ref, watch, onMounted } = window.Vue;

const parsePositiveInteger = (value) => {
  const parsed = Number(value);
  return Number.isInteger(parsed) && parsed >= 1 ? parsed : null;
};

const createPositiveIntegerOptions = (maxValue) =>
  Array.from({ length: Math.max(0, Math.floor(maxValue)) }, (_unused, index) => {
    const value = index + 1;
    return {
      value: String(value),
      label: `${value}`
    };
  });

const appendRequestedIntegerOption = (options, requestedValue, effectiveValue) => {
  const requested = parsePositiveInteger(requestedValue);
  if (requested === null) return options;
  const requestedOption = {
    value: String(requested),
    label: requested === effectiveValue ? `${requested}` : `${requested} (${effectiveValue} effective)`
  };
  return options.some((option) => option.value === requestedOption.value)
    ? options
    : [...options, requestedOption];
};

/**
 * The saved LOSAT Settings values. A thread value is 'safe', 'available', a
 * count, or blank (Auto), as `resolveLosatThreadPlan` reads it.
 * @typedef {object} LosatSettingsValues
 * @property {string} [executionMode]
 * @property {string | number | null} [threadsPerJob]
 * @property {string | number | null} [parallelWorkers]
 * @property {string | number | null} [totalThreadBudget]
 * @property {{ mode?: string, collinearInferOrthogroups?: boolean, collinearSearchScope?: string }} [blastp]
 */

/**
 * The refs and reactive objects of `state.js` this module reads. A ref is
 * `{ value }` because Vue comes from `window.Vue`, which is `any`.
 * @typedef {object} LosatSettingsState
 * @property {Record<string, any>[]} linearSeqs
 * @property {{ value: any }} linearComparisonResolution
 * @property {LosatSettingsValues} losat
 * @property {{ value: string }} losatProgram
 * @property {{ value?: { state?: string } }} [losatThreadingStatus]
 */

/**
 * @param {{ state: LosatSettingsState }} options
 */
export const createLosatSettings = ({ state }) => {
  const {
    linearSeqs,
    linearComparisonResolution,
    losat,
    losatProgram
  } = state;

  // `state.js` always provides the computed plan.
  const readResolution = () => linearComparisonResolution.value || {};

  const losatHardwareThreads = ref(getLosatHardwareThreads());
  onMounted(() => {
    losatHardwareThreads.value = getLosatHardwareThreads();
  });

  const losatThreadsPerJobFixed = computed(() => losatProgram.value !== 'blastp' || losat.executionMode === 'serial');

  // Generate plans its jobs with the same two functions (N-10). Only the
  // translation tables vary per record, so they are the only arguments that
  // can split a source batch.
  const losatEstimatedJobCount = computed(() => {
    const resolution = readResolution();
    if (resolution.valid === false || !resolution.hasLosatIntent) return 0;
    const program = losatProgram.value;
    try {
      const specs = buildLosatJobSpecs({
        resolution,
        recordCount: linearSeqs.length,
        recordUids: linearSeqs.map((sequence) => sequence?.uid),
        program,
        blastpMode: String(losat.blastp?.mode || 'orthogroup'),
        collinearInferOrthogroups: losat.blastp?.collinearInferOrthogroups !== false,
        collinearSearchScope: normalizeCollinearSearchScope(losat.blastp?.collinearSearchScope)
      });
      return planLosatSourceJobs({
        sequences: linearSeqs,
        specs,
        buildArgs: (query, subject) => (program === 'tblastx'
          ? [losatRecordGencode(linearSeqs[query]), losatRecordGencode(linearSeqs[subject])]
          : [])
      }).jobs.length;
    } catch {
      return 0;
    }
  });

  const losatThreadPlan = computed(() => resolveLosatThreadPlan({
    hardwareThreads: losatHardwareThreads.value,
    jobCount: losatEstimatedJobCount.value,
    totalThreadBudget: losat.totalThreadBudget,
    threadsPerJob: losatThreadsPerJobFixed.value ? 1 : losat.threadsPerJob,
    parallelWorkers: losat.parallelWorkers
  }));
  const losatSafeThreadBudget = computed(() =>
    resolveLosatThreadPlan({ hardwareThreads: losatHardwareThreads.value }).totalBudget
  );
  const losatTotalThreadBudget = computed(() => losatThreadPlan.value.totalBudget);
  const losatTotalThreadBudgetOptions = computed(() =>
    createPositiveIntegerOptions(losatHardwareThreads.value)
  );
  const losatEffectiveThreadsPerJob = computed(() => losatThreadPlan.value.threadsPerJob);

  const losatThreadOptions = computed(() => {
    if (losatThreadsPerJobFixed.value) {
      return appendRequestedIntegerOption([{ value: '1', label: 'Fixed (1)' }], losat.threadsPerJob, 1);
    }
    const maxThreads = Math.max(1, losatTotalThreadBudget.value);
    return appendRequestedIntegerOption(
      createPositiveIntegerOptions(maxThreads),
      losat.threadsPerJob,
      losatEffectiveThreadsPerJob.value
    );
  });

  const losatMaxPairWorkers = computed(() => losatThreadPlan.value.maxPairWorkers);
  const losatAutoPairWorkers = computed(() => losatThreadPlan.value.autoPairWorkers);
  const losatPairWorkerOptions = computed(() => losatEstimatedJobCount.value === 0 ? [] : appendRequestedIntegerOption(
    createPositiveIntegerOptions(losatMaxPairWorkers.value).map((option) => ({
      ...option, label: `${option.value} ${option.value === '1' ? 'run' : 'runs'}`
    })),
    losat.parallelWorkers,
    losatThreadPlan.value.pairWorkers
  ));

  const losatEffectiveExecutionMode = computed(() => {
    const raw = String(losat.executionMode || 'auto').trim().toLowerCase();
    if (raw === 'serial') return 'serial';
    if (raw === 'threaded') return 'threaded';
    const threadingState = String(state?.losatThreadingStatus?.value?.state || '').trim().toLowerCase();
    if (threadingState && threadingState !== 'available' && threadingState !== 'running') return 'serial';
    if (losatMaxPairWorkers.value > 1 || losatEffectiveThreadsPerJob.value > 1) return 'threaded if useful';
    return 'serial';
  });

  // Threaded stays strict (PD-OI-018); the option says when this browser
  // environment cannot run it.
  const losatThreadedOptionLabel = computed(() => (
    losatThreadingPrecondition().state === 'available' ? 'Threaded' : 'Threaded (unavailable here)'
  ));

  const losatThreadingPlanSummary = computed(() =>
    'By default, LOSAT can use up to half the number of cores available.'
  );

  const hasValidLosatIntent = () => {
    const resolution = readResolution();
    return resolution.valid === true && resolution.hasLosatIntent === true;
  };

  watch(
    losatTotalThreadBudgetOptions,
    (options) => {
      if (!hasValidLosatIntent()) return;
      const raw = String(losat.totalThreadBudget || 'safe').trim().toLowerCase();
      if (['safe', 'available'].includes(raw)) return;
      const values = options.map((option) => option.value);
      if (values.includes(raw)) return;
      const parsed = parsePositiveInteger(raw);
      losat.totalThreadBudget = parsed !== null && parsed >= losatHardwareThreads.value
        ? 'available'
        : 'safe';
    },
    { immediate: true }
  );

  return {
    losatHardwareThreads,
    losatSafeThreadBudget,
    losatTotalThreadBudget,
    losatTotalThreadBudgetOptions,
    losatThreadsPerJobFixed,
    losatThreadOptions,
    losatEffectiveThreadsPerJob,
    losatEffectiveExecutionMode,
    losatThreadedOptionLabel,
    losatEstimatedJobCount,
    losatMaxPairWorkers,
    losatAutoPairWorkers,
    losatPairWorkerOptions,
    losatThreadingPlanSummary
  };
};
