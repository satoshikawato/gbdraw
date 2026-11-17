// @ts-check
import { readCollinearSearchScope } from '../services/losat-normalization.js';
import { requireCurrentCollinearInferOrthogroups } from '../services/current-option-values.js';
import { buildLosatJobSpecs } from '../services/linear-comparisons.js';
import { drawsSeveralUnits, linearUnitPlanPayload, planLinearUnitSearch } from '../services/linear-units.js';
import { losatSearchSequences, nucleotideLosatArgs, planLosatSourceJobs, proteinLosatArgs } from '../services/linear-sources.js';
import { losatRecordGencode } from '../services/drawing/inputs.js';
import { drawnLinearSequences } from '../services/record-draw-selection.js';
import { getLosatHardwareThreads, resolveLosatThreadPlan } from '../services/losat-thread-plan.js';
import { losatThreadingPrecondition } from '../services/losat.js';
import { t, tn } from '../utils/i18n.js';
/** @import { LinearComparisonResolution } from '../services/linear-comparisons.js' */
/** @import { DrawingState, LinearSeq, WebState } from '../state.js' */
/** @import { LinearUnitPlanPayload, LinearUnitSearch, LinearCardPlacesOf } from '../services/linear-units.js' */

const { computed, ref, watch, onMounted } = window.Vue;

/** @typedef {{ value: string, label: string }} IntegerOption */

/** @param {unknown} value */
const parsePositiveInteger = (value) => {
  const parsed = Number(value);
  return Number.isInteger(parsed) && parsed >= 1 ? parsed : null;
};

/**
 * @param {number} maxValue
 * @returns {IntegerOption[]}
 */
const createPositiveIntegerOptions = (maxValue) =>
  Array.from({ length: Math.max(0, Math.floor(maxValue)) }, (_unused, index) => {
    const value = index + 1;
    return {
      value: String(value),
      label: `${value}`
    };
  });

/**
 * @param {IntegerOption[]} options
 * @param {unknown} requestedValue
 * @param {number} effectiveValue
 */
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
 * The members of `state.js` this module reads. `losatExecution` is how LOSAT
 * runs, one app-level setting for both drawings; a thread value in it is
 * 'safe', 'available', a count, or blank (Auto), as `resolveLosatThreadPlan`
 * reads it.
 * @typedef {Pick<WebState, 'activeDrawing' | 'losatExecution' | 'losatThreadingStatus'>} LosatSettingsState
 */

/**
 * The members of a drawing this module reads: its comparison plan as resolved, its LOSAT settings and program, and
 * its record rows with the records it leaves out.
 * @typedef {Pick<DrawingState, 'linearComparisonResolution' | 'losat' | 'losatProgram' | 'recordsOff' | 'linearSeqs'>}
 *   LosatSettingsDrawing
 */

/**
 * How the estimate plans drawn units: the Worker helper that plans them, and the places a drawing's request writes for
 * its drawn cards (`linearDrawingUnitPlaces` of services/config.js).
 * @typedef {object} LosatUnitPlanPorts
 * @property {(payload: LinearUnitPlanPayload) => Promise<{ result?: unknown }>} planLinearUnits
 * @property {() => LinearCardPlacesOf} shownPlacesOf The places the drawn cards of the shown drawing write in its
 *   request.
 */

/**
 * The drawn units of the shown drawing's plan when a drawn card draws several windows, from Python as Generate plans
 * them (`planLinearUnitSearch`, OV-2741); otherwise the cards are the units. `unitPlanKey` is the payload Python
 * answers ('' when no plan is needed), so a plan of another payload is not read.
 * @param {LosatSettingsState} state
 * @param {(drawing: LosatSettingsDrawing) => LinearComparisonResolution} readResolution
 * @param {LosatUnitPlanPorts} ports
 */
const watchLosatUnitPlan = (state, readResolution, { planLinearUnits, shownPlacesOf }) => {
  const unitPlanKey = computed(() => {
    const drawing = state.activeDrawing();
    const resolution = readResolution(drawing);
    const drawn = drawnLinearSequences(drawing.linearSeqs, drawing.recordsOff);
    if (resolution.valid === false || !resolution.hasLosatIntent || !drawsSeveralUnits(drawn)) return '';
    try {
      return JSON.stringify(linearUnitPlanPayload(drawn, resolution, shownPlacesOf()));
    } catch {
      return 'refused'; // A window Generate refuses while it is typed: no plan, as when the helper fails.
    }
  });
  /** @type {{ value: { key: string, search: LinearUnitSearch | null } | null }} */
  const unitPlan = ref(null);
  watch(() => unitPlanKey.value, /** @param {string} key */ async (key) => {
    if (!key) return;
    const drawing = state.activeDrawing();
    let search = null;
    try {
      search = await planLinearUnitSearch(drawnLinearSequences(drawing.linearSeqs, drawing.recordsOff),
        readResolution(drawing), shownPlacesOf(), planLinearUnits);
    } catch {
      search = null;
    }
    if (unitPlanKey.value === key) unitPlan.value = { key, search };
  }, { immediate: true });
  return { unitPlanKey, unitPlan };
};

/**
 * The arguments the job estimate plans a drawing's searches with: Generate's builders with no E-value scale. A scale
 * is one per genome pair, and a source batch is one, so it never splits a batch; the translation tables of TLOSATX
 * records do.
 * @param {LosatSettingsDrawing} drawing
 * @param {readonly (LinearSeq | undefined)[]} sequences The searched records, by record index.
 * @returns {(query: number, subject: number) => string[]}
 */
const estimateLosatArgs = (drawing, sequences) => (query, subject) => {
  const program = drawing.losatProgram.value;
  if (program !== 'blastn' && program !== 'tblastx') {
    return proteinLosatArgs({ pairwise: drawing.losat.blastp?.mode === 'pairwise', candidateLimit: null, databaseSize: null });
  }
  // The typed tables as Generate passes them; Generate refuses one LOSAT does not read (OV-3326).
  return nucleotideLosatArgs({
    program, task: drawing.losat.blastn.task, searchSpace: null,
    gencodes: () => [losatRecordGencode(sequences[query]), losatRecordGencode(sequences[subject])]
  });
};

/**
 * @param {{ state: LosatSettingsState } & LosatUnitPlanPorts} options
 */
export const createLosatSettings = ({ state, ...unitPlanPorts }) => {
  const { losatExecution } = state;

  // `state.js` always provides the computed plan.
  /** @param {LosatSettingsDrawing} drawing */
  const readResolution = (drawing) => drawing.linearComparisonResolution.value;

  const losatHardwareThreads = ref(getLosatHardwareThreads());
  onMounted(() => {
    losatHardwareThreads.value = getLosatHardwareThreads();
  });

  const losatThreadsPerJobFixed = computed(() => {
    return state.activeDrawing().losatProgram.value !== 'blastp' || losatExecution.executionMode === 'serial';
  });

  const { unitPlanKey, unitPlan } = watchLosatUnitPlan(state, readResolution, unitPlanPorts);
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
      const key = unitPlanKey.value;
      const planned = unitPlan.value;
      const units = key && planned && planned.key === key ? planned.search : null;
      // A plan Python has not answered yet, or could not make, counts no job.
      if (key && !units) return 0;
      const cards = drawnLinearSequences(drawing.linearSeqs, drawing.recordsOff);
      const drawn = units ? units.sequences : cards;
      const search = losatSearchSequences(drawn, drawing.linearSeqs, cards);
      const specs = buildLosatJobSpecs({
        resolution: units ? units.resolution : resolution,
        recordCount: drawn.length,
        recordUids: drawn.map((sequence) => sequence?.uid),
        program,
        blastpMode: String(drawing.losat.blastp?.mode || 'orthogroup'),
        collinearInferOrthogroups: requireCurrentCollinearInferOrthogroups(drawing.losat.blastp?.collinearInferOrthogroups),
        collinearSearchScope: readCollinearSearchScope(drawing.losat.blastp?.collinearSearchScope).value
      });
      return planLosatSourceJobs({
        ...search,
        specs,
        buildArgs: estimateLosatArgs(drawing, search.sequences)
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
      ...option, label: tn(Number(option.value), '{count} parallel run', '{count} parallel runs')
    })),
    losatExecution.parallelWorkers,
    losatThreadPlan.value.pairWorkers
  ));

  // The shown name of the mode Auto resolves to (UI text, compared nowhere): threaded wherever this browser can
  // run it (D-51), serial when it cannot or the last run fell back.
  const losatEffectiveExecutionMode = computed(() => {
    const raw = String(losatExecution.executionMode || 'auto').trim().toLowerCase();
    if (raw === 'serial') return t('serial');
    if (raw === 'threaded') return t('threaded');
    const threadingState = String(state?.losatThreadingStatus?.value?.state || '').trim().toLowerCase();
    if (threadingState && threadingState !== 'available' && threadingState !== 'running') return t('serial');
    return losatThreadingPrecondition().state === 'available' ? t('threaded') : t('serial');
  });

  // Threaded stays strict (PD-OI-018); the option says when this browser
  // environment cannot run it.
  const losatThreadedOptionLabel = computed(() => (
    losatThreadingPrecondition().state === 'available' ? t('Threaded') : t('Threaded (unavailable here)')
  ));

  /** @param {LosatSettingsDrawing} drawing */
  const hasValidLosatIntent = (drawing) => {
    const resolution = readResolution(drawing);
    return resolution.valid === true && resolution.hasLosatIntent === true;
  };

  watch(
    losatTotalThreadBudgetOptions,
    /** @param {IntegerOption[]} options */
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
    losatPairWorkerOptions
  };
};
