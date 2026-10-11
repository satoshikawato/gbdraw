// @ts-check
import { normalizeCollinearSearchScope } from '../services/losat-normalization.js';
import { buildLosatJobSpecs } from '../services/linear-comparisons.js';
import { losatRecordGencode, losatSearchSequences, planLosatSourceJobs } from '../services/linear-sources.js';
import { drawnLinearSequences } from '../services/record-draw-selection.js';
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
 * How LOSAT runs (`state.losatExecution`, one app-level setting for both
 * drawings). A thread value is 'safe', 'available', a count, or blank (Auto),
 * as `resolveLosatThreadPlan` reads it.
 * @typedef {object} LosatExecutionValues
 * @property {string} [executionMode]
 * @property {string | number | null} [threadsPerJob]
 * @property {string | number | null} [parallelWorkers]
 * @property {string | number | null} [totalThreadBudget]
 */

/**
 * A drawing's LOSAT search settings this module reads.
 * @typedef {object} LosatSettingsValues
 * @property {{ mode?: string, collinearInferOrthogroups?: boolean, collinearSearchScope?: string }} [blastp]
 */

/**
 * The refs and reactive objects of `state.js` this module reads. A ref is
 * `{ value }` because Vue comes from `window.Vue`, which is `any`.
 * @typedef {object} LosatSettingsState
 * @property {Record<string, any>[]} linearSeqs
 * @property {() => LosatSettingsDrawing} activeDrawing
 * @property {LosatExecutionValues} losatExecution
 * @property {{ value?: { state?: string } }} [losatThreadingStatus]
 */

/**
 * The members of a drawing (`DrawingState` of state.js) this module reads.
 * @typedef {object} LosatSettingsDrawing
 * @property {{ value: any }} linearComparisonResolution
 * @property {LosatSettingsValues} losat
 * @property {{ value: string }} losatProgram
 * @property {readonly string[]} [recordsOff]
 */

/**
 * @param {{ state: LosatSettingsState }} options
 */
export const createLosatSettings = ({ state }) => {
  const {
    linearSeqs,
    losatExecution
  } = state;

  // `state.js` always provides the computed plan.
  /** @param {LosatSettingsDrawing} drawing */
  const readResolution = (drawing) => drawing.linearComparisonResolution.value || {};

  const losatHardwareThreads = ref(getLosatHardwareThreads());
  onMounted(() => {
    losatHardwareThreads.value = getLosatHardwareThreads();
  });

  const losatThreadsPerJobFixed = computed(() => {
    return state.activeDrawing().losatProgram.value !== 'blastp' || losatExecution.executionMode === 'serial';
  });

  // Generate plans its jobs with the same two functions (N-10). Only the
  // translation tables vary per record, so they are the only arguments that
  // can split a source batch. Pairs join drawn records, and OFF records stay
  // in their source database (record selection D-07).
  const losatEstimatedJobCount = computed(() => {
    const drawing = state.activeDrawing();
    const resolution = readResolution(drawing);
    if (resolution.valid === false || !resolution.hasLosatIntent) return 0;
    const program = drawing.losatProgram.value;
    try {
      const drawn = drawnLinearSequences(linearSeqs, drawing.recordsOff);
      const search = losatSearchSequences(drawn, linearSeqs);
      const specs = buildLosatJobSpecs({
        resolution,
        recordCount: drawn.length,
        recordUids: drawn.map((sequence) => sequence?.uid),
        program,
        blastpMode: String(drawing.losat.blastp?.mode || 'orthogroup'),
        collinearInferOrthogroups: drawing.losat.blastp?.collinearInferOrthogroups !== false,
        collinearSearchScope: normalizeCollinearSearchScope(drawing.losat.blastp?.collinearSearchScope)
      });
      return planLosatSourceJobs({
        ...search,
        specs,
        buildArgs: (query, subject) => (program === 'tblastx'
          ? [losatRecordGencode(search.sequences[query]), losatRecordGencode(search.sequences[subject])]
          : [])
      }).jobs.length;
    } catch {
      return 0;
    }
  });

  const losatThreadPlan = computed(() => resolveLosatThreadPlan({
    hardwareThreads: losatHardwareThreads.value,
    jobCount: losatEstimatedJobCount.value,
    totalThreadBudget: losatExecution.totalThreadBudget,
    threadsPerJob: losatThreadsPerJobFixed.value ? 1 : losatExecution.threadsPerJob,
    parallelWorkers: losatExecution.parallelWorkers
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
      return appendRequestedIntegerOption([{ value: '1', label: 'Fixed (1)' }], losatExecution.threadsPerJob, 1);
    }
    const maxThreads = Math.max(1, losatTotalThreadBudget.value);
    return appendRequestedIntegerOption(
      createPositiveIntegerOptions(maxThreads),
      losatExecution.threadsPerJob,
      losatEffectiveThreadsPerJob.value
    );
  });

  const losatMaxPairWorkers = computed(() => losatThreadPlan.value.maxPairWorkers);
  const losatAutoPairWorkers = computed(() => losatThreadPlan.value.autoPairWorkers);
  const losatPairWorkerOptions = computed(() => losatEstimatedJobCount.value === 0 ? [] : appendRequestedIntegerOption(
    createPositiveIntegerOptions(losatMaxPairWorkers.value).map((option) => ({
      ...option, label: `${option.value} ${option.value === '1' ? 'run' : 'runs'}`
    })),
    losatExecution.parallelWorkers,
    losatThreadPlan.value.pairWorkers
  ));

  const losatEffectiveExecutionMode = computed(() => {
    const raw = String(losatExecution.executionMode || 'auto').trim().toLowerCase();
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

  /** @param {LosatSettingsDrawing} drawing */
  const hasValidLosatIntent = (drawing) => {
    const resolution = readResolution(drawing);
    return resolution.valid === true && resolution.hasLosatIntent === true;
  };

  watch(
    losatTotalThreadBudgetOptions,
    (options) => {
      if (!hasValidLosatIntent(state.activeDrawing())) return;
      const raw = String(losatExecution.totalThreadBudget || 'safe').trim().toLowerCase();
      if (['safe', 'available'].includes(raw)) return;
      const values = options.map((option) => option.value);
      if (values.includes(raw)) return;
      const parsed = parsePositiveInteger(raw);
      losatExecution.totalThreadBudget = parsed !== null && parsed >= losatHardwareThreads.value
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
