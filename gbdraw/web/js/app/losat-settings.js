import { normalizeCollinearSearchScope } from './losat-normalization.js';
import { groupLinearSourceRecords } from './linear-sources.js';
import { getLosatHardwareThreads, resolveLosatThreadPlan } from '../services/losat-thread-plan.js';

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

export const createLosatSettings = ({ state }) => {
  const {
    linearSeqs,
    linearComparisonResolution,
    losat,
    losatProgram
  } = state;

  const losatHardwareThreads = ref(getLosatHardwareThreads());
  onMounted(() => {
    losatHardwareThreads.value = getLosatHardwareThreads();
  });

  const losatThreadsPerJobFixed = computed(() => losatProgram.value !== 'blastp' || losat.executionMode === 'serial');

  const losatEstimatedJobCount = computed(() => {
    const resolution = linearComparisonResolution?.value || linearComparisonResolution || {};
    if (resolution.valid === false || !resolution.hasLosatIntent) return 0;
    const sources = new Map();
    groupLinearSourceRecords(linearSeqs).forEach((group) => {
      group.records.forEach(({ index }) => sources.set(index, group.uid));
    });
    const jobs = new Set();
    const addPair = (query, subject) => jobs.add(JSON.stringify([
      sources.get(query), sources.get(subject),
      ...(losatProgram.value === 'tblastx'
        ? [linearSeqs[query]?.losat_gencode, linearSeqs[subject]?.losat_gencode] : [])
    ]));
    const edges = (resolution.edges || []).filter((edge) => edge.source === 'losat');
    edges.forEach((edge) => addPair(edge.queryIndex, edge.subjectIndex));
    const blastpMode = String(losat.blastp?.mode || 'orthogroup').trim().toLowerCase();
    if (losatProgram.value === 'blastp' && ['orthogroup', 'collinear'].includes(blastpMode)
      && resolution.mode === 'adjacent' && resolution.defaultSource === 'losat') {
      linearSeqs.forEach((_, index) => addPair(index, index));
      if (blastpMode === 'orthogroup'
        || normalizeCollinearSearchScope(losat.blastp?.collinearSearchScope) === 'all') {
        const representatives = [...new Set(sources.values())].map((uid) =>
          [...sources].find(([, sourceUid]) => sourceUid === uid)[0]);
        representatives.forEach((query) => representatives.forEach((subject) => addPair(query, subject)));
      } else {
        edges.forEach((edge) => addPair(edge.subjectIndex, edge.queryIndex));
      }
    }
    return jobs.size;
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

  const losatThreadingPlanSummary = computed(() =>
    'By default, LOSAT can use up to half the number of cores available.'
  );

  const hasValidLosatIntent = () => {
    const resolution = linearComparisonResolution?.value || linearComparisonResolution || {};
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
    losatEstimatedJobCount,
    losatMaxPairWorkers,
    losatAutoPairWorkers,
    losatPairWorkerOptions,
    losatThreadingPlanSummary
  };
};
