const { computed, ref, watch, onMounted } = window.Vue;

const DEFAULT_LOSAT_PAIR_WORKER_AUTO_LIMIT = 4;
const LOSAT_THREAD_OPTION_BASES = [1, 2, 4, 8, 16, 24, 32, 48, 64, 96, 128];

const getBrowserHardwareThreads = () =>
  Math.max(1, Number(globalThis.navigator?.hardwareConcurrency || 4) || 4);

const getSafeLosatThreadBudget = (hardwareThreads) =>
  Math.max(1, Math.floor(hardwareThreads / 2));

const parsePositiveInteger = (value) => {
  const parsed = Number(value);
  return Number.isInteger(parsed) && parsed >= 1 ? parsed : null;
};

const uniqueSortedPositiveIntegers = (values) =>
  Array.from(new Set(values.filter((value) => Number.isInteger(value) && value >= 1)))
    .sort((left, right) => left - right);

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

  const losatEstimatedJobCount = computed(() => {
    const recordCount = Math.max(0, Array.isArray(linearSeqs) ? linearSeqs.length : 0);
    if (recordCount < 2) return 1;
    if (losatProgram.value !== 'blastp') return recordCount - 1;

    const blastpMode = String(losat.blastp?.mode || 'orthogroup').trim().toLowerCase();
    if (blastpMode === 'orthogroup') return recordCount * (recordCount - 1);
    if (blastpMode === 'collinear') {
      const scope = String(losat.blastp?.collinearSearchScope || 'adjacent').trim().toLowerCase();
      const pairCount = scope === 'all'
        ? Math.floor((recordCount * (recordCount - 1)) / 2)
        : recordCount - 1;
      const directionMultiplier = String(losat.blastp?.collinearAnchorMode || '').trim().toLowerCase() === 'rbh'
        ? 2
        : 1;
      return Math.max(1, pairCount * directionMultiplier);
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
    const hardware = losatHardwareThreads.value;
    const reservedValues = new Set([losatSafeThreadBudget.value, hardware]);
    const optionValues = uniqueSortedPositiveIntegers([
      ...LOSAT_THREAD_OPTION_BASES.filter((value) => value <= hardware && !reservedValues.has(value))
    ]);
    return optionValues.map((value) => ({
      value: String(value),
      label: `${value}`
    }));
  });

  const getLosatAutoThreadsPerJob = () => {
    const hardwareBudget = losatSafeThreadBudget.value;
    if (losatEstimatedJobCount.value !== 1) return Math.max(1, Math.min(2, hardwareBudget));
    return Math.max(1, Math.min(hardwareBudget, Math.max(1, losatHardwareThreads.value - 1)));
  };

  const losatEffectiveThreadsPerJob = computed(() => {
    const raw = String(losat.threadsPerJob || 'auto').trim().toLowerCase();
    const requested = raw === 'auto'
      ? getLosatAutoThreadsPerJob()
      : parsePositiveInteger(raw) || getLosatAutoThreadsPerJob();
    return Math.max(1, Math.min(requested, Math.max(1, losatTotalThreadBudget.value - 1)));
  });

  const losatThreadOptions = computed(() => {
    const maxThreads = Math.max(1, losatTotalThreadBudget.value - 1);
    const optionValues = uniqueSortedPositiveIntegers([
      maxThreads,
      ...LOSAT_THREAD_OPTION_BASES.filter((value) => value <= maxThreads)
    ]);
    return optionValues.map((value) => ({
      value: String(value),
      label: `${value}`
    }));
  });

  const losatMaxPairWorkers = computed(() => {
    const perJobSlots = losatEffectiveThreadsPerJob.value + 1;
    const budgetLimited = Math.max(1, Math.floor(losatTotalThreadBudget.value / perJobSlots));
    return Math.max(1, Math.min(losatEstimatedJobCount.value, budgetLimited));
  });

  const losatAutoPairWorkers = computed(() => {
    const budgetMode = String(losat.totalThreadBudget || 'safe').trim().toLowerCase();
    if (!['safe', 'auto'].includes(budgetMode)) return losatMaxPairWorkers.value;
    const hardwareLimit = Math.max(1, losatHardwareThreads.value - 1);
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
        label: `${value} ${value === 1 ? 'pair' : 'pairs'}`
      };
    });
  });

  const losatThreadingPlanSummary = computed(() => {
    const selectedWorkers = parsePositiveInteger(losat.parallelWorkers) || losatAutoPairWorkers.value;
    const effectiveWorkers = Math.min(selectedWorkers, losatMaxPairWorkers.value);
    return `${losatHardwareThreads.value} cores reported; LOSAT can use up to ${losatTotalThreadBudget.value} cores. ` +
      `Selected: ${effectiveWorkers} parallel ${effectiveWorkers === 1 ? 'pair' : 'pairs'}, ` +
      `${losatEffectiveThreadsPerJob.value} LOSAT thread${losatEffectiveThreadsPerJob.value === 1 ? '' : 's'} per pair.`;
  });

  watch(
    losatThreadOptions,
    (options) => {
      if (String(losat.threadsPerJob || 'auto').trim().toLowerCase() === 'auto') return;
      const values = options.map((option) => option.value);
      if (values.includes(String(losat.threadsPerJob))) return;
      losat.threadsPerJob = options[options.length - 1]?.value || 'auto';
    },
    { immediate: true }
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
    losatThreadOptions,
    losatEffectiveThreadsPerJob,
    losatAutoPairWorkers,
    losatPairWorkerOptions,
    losatThreadingPlanSummary
  };
};
