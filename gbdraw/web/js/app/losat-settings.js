import { normalizeCollinearSearchScope } from './losat-normalization.js';

const { computed, ref, watch, onMounted } = window.Vue;

const DEFAULT_LOSAT_PAIR_WORKER_AUTO_LIMIT = 4;

const getBrowserHardwareThreads = () =>
  Math.max(1, Number(globalThis.navigator?.hardwareConcurrency || 4) || 4);

const getSafeLosatThreadBudget = (hardwareThreads) =>
  Math.max(1, Math.floor(hardwareThreads / 2));

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

const appendRequestedIntegerOption = (options, requestedValue) => {
  const requested = parsePositiveInteger(requestedValue);
  if (requested === null) return options;
  const requestedOption = {
    value: String(requested),
    label: `${requested}`
  };
  return options.some((option) => option.value === requestedOption.value)
    ? options
    : [...options, requestedOption];
};

export const createLosatSettings = ({ state }) => {
  const {
    linearSeqs,
    losat,
    losatProgram
  } = state;

  const losatHardwareThreads = ref(getBrowserHardwareThreads());
  onMounted(() => {
    losatHardwareThreads.value = getBrowserHardwareThreads();
  });

  const losatThreadsPerJobFixed = computed(() => losatProgram.value !== 'blastp');

  const losatEstimatedJobCount = computed(() => {
    const recordCount = Math.max(0, Array.isArray(linearSeqs) ? linearSeqs.length : 0);
    if (recordCount < 2) return 1;
    if (losatProgram.value !== 'blastp') return recordCount - 1;

    const blastpMode = String(losat.blastp?.mode || 'orthogroup').trim().toLowerCase();
    if (blastpMode === 'orthogroup') return recordCount * recordCount;
    if (blastpMode === 'collinear') {
      const scope = normalizeCollinearSearchScope(losat.blastp?.collinearSearchScope);
      const pairCount = scope === 'all'
        ? Math.floor((recordCount * (recordCount - 1)) / 2)
        : recordCount - 1;
      return Math.max(1, recordCount + pairCount * 2);
    }
    return recordCount - 1;
  });

  const losatSafeThreadBudget = computed(() =>
    getSafeLosatThreadBudget(losatHardwareThreads.value)
  );

  const losatTotalThreadBudget = computed(() => {
    const raw = String(losat.totalThreadBudget || 'safe').trim().toLowerCase();
    if (raw === 'available') return losatHardwareThreads.value;
    const parsed = parsePositiveInteger(raw);
    if (parsed !== null) return Math.min(parsed, losatHardwareThreads.value);
    return losatSafeThreadBudget.value;
  });

  const losatTotalThreadBudgetOptions = computed(() => {
    return createPositiveIntegerOptions(losatHardwareThreads.value);
  });

  const getLosatAutoThreadsPerJob = () => {
    if (losatThreadsPerJobFixed.value) return 1;
    const hardwareBudget = losatSafeThreadBudget.value;
    if (losatEstimatedJobCount.value !== 1) return Math.max(1, Math.min(2, hardwareBudget));
    return Math.max(1, hardwareBudget);
  };

  const losatEffectiveThreadsPerJob = computed(() => {
    if (losatThreadsPerJobFixed.value) return 1;
    const raw = String(losat.threadsPerJob || 'auto').trim().toLowerCase();
    const requested = raw === 'auto'
      ? getLosatAutoThreadsPerJob()
      : parsePositiveInteger(raw) || getLosatAutoThreadsPerJob();
    return Math.max(1, Math.min(requested, Math.max(1, losatTotalThreadBudget.value)));
  });

  const losatThreadOptions = computed(() => {
    if (losatThreadsPerJobFixed.value) {
      return [{ value: '1', label: 'Fixed (1)' }];
    }
    const maxThreads = Math.max(1, losatTotalThreadBudget.value);
    return appendRequestedIntegerOption(
      createPositiveIntegerOptions(maxThreads),
      losat.threadsPerJob
    );
  });

  const losatMaxPairWorkers = computed(() => {
    const perJobSlots = losatEffectiveThreadsPerJob.value;
    const budgetLimited = Math.max(1, Math.floor(losatTotalThreadBudget.value / perJobSlots));
    return Math.max(1, Math.min(losatEstimatedJobCount.value, budgetLimited));
  });

  const losatAutoPairWorkers = computed(() => {
    const budgetMode = String(losat.totalThreadBudget || 'safe').trim().toLowerCase();
    if (!['safe', 'auto'].includes(budgetMode)) return losatMaxPairWorkers.value;
    const hardwareLimit = Math.max(1, losatHardwareThreads.value);
    const defaultConcurrency = Math.min(
      losatEstimatedJobCount.value,
      hardwareLimit,
      DEFAULT_LOSAT_PAIR_WORKER_AUTO_LIMIT
    );
    return Math.max(1, Math.min(defaultConcurrency, losatMaxPairWorkers.value));
  });

  const losatPairWorkerOptions = computed(() => {
    return Array.from({ length: losatMaxPairWorkers.value }, (_, index) => {
      const value = index + 1;
      return {
        value: String(value),
        label: `${value} ${value === 1 ? 'run' : 'runs'}`
      };
    });
  });

  const losatThreadingPlanSummary = computed(() =>
    'By default, LOSAT can use up to half the number of cores available.'
  );

  watch(
    losatTotalThreadBudgetOptions,
    (options) => {
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

  watch(
    losatPairWorkerOptions,
    (options) => {
      if (losat.parallelWorkers === undefined || losat.parallelWorkers === null) return;
      const values = options.map((option) => option.value);
      if (values.includes(String(losat.parallelWorkers))) return;
      losat.parallelWorkers = options[options.length - 1]?.value || undefined;
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
    losatAutoPairWorkers,
    losatPairWorkerOptions,
    losatThreadingPlanSummary
  };
};
