import assert from 'node:assert/strict';
import { mkdtemp, readFile, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';

globalThis.window = {
  Vue: {
    ref: (value) => ({ value }),
    computed: (getter) => ({
      get value() {
        return getter();
      }
    }),
    watch: (source, callback, options = {}) => {
      if (!options.immediate) return;
      const value = typeof source === 'function' ? source() : source.value;
      callback(value);
    },
    onMounted: (callback) => {
      callback();
    }
  }
};

Object.defineProperty(globalThis, 'navigator', {
  value: { hardwareConcurrency: 8 },
  configurable: true
});

const sourceRoot = new URL('../../gbdraw/web/js/app/', import.meta.url);
const tempDir = await mkdtemp(join(tmpdir(), 'gbdraw-losat-settings-'));
const tempModulePath = join(tempDir, 'losat-settings.mjs');
const tempNormalizationPath = join(tempDir, 'losat-normalization.mjs');
const source = await readFile(new URL('losat-settings.js', sourceRoot), 'utf8');
await writeFile(
  tempModulePath,
  source.replace("./losat-normalization.js", "./losat-normalization.mjs")
    .replace("./linear-sources.js", new URL("linear-sources.js", sourceRoot).href)
    .replace("../services/losat-thread-plan.js", new URL("../services/losat-thread-plan.js", sourceRoot).href)
);
await writeFile(
  tempNormalizationPath,
  await readFile(new URL('losat-normalization.js', sourceRoot), 'utf8')
);

const { createLosatSettings } = await import(pathToFileURL(tempModulePath));

const resolved = ({ mode = 'adjacent', defaultSource = 'losat', sources = [], valid = true } = {}) => ({
  value: {
    mode,
    defaultSource,
    valid,
    hasLosatIntent: sources.includes('losat'),
    edges: sources.map((source, ordinal) => ({ source, ordinal, queryIndex: ordinal, subjectIndex: ordinal + 1 }))
  }
});

const state = {
  linearSeqs: Array.from({ length: 5 }, (_, index) => ({ uid: `record-${index}` })),
  linearComparisonResolution: resolved({ sources: ['losat', 'losat', 'losat', 'losat'] }),
  losat: {
    totalThreadBudget: 'safe',
    threadsPerJob: '32',
    parallelWorkers: undefined,
    blastp: {
      mode: 'orthogroup',
      collinearSearchScope: 'adjacent'
    }
  },
  losatProgram: { value: 'blastp' }
};

const settings = createLosatSettings({ state });

assert.equal(state.losat.threadsPerJob, '32');
assert.equal(settings.losatEffectiveThreadsPerJob.value, 4);
assert(settings.losatThreadOptions.value.some((option) => option.value === '32'));

const loadingOrderState = {
  linearSeqs: Array.from({ length: 5 }, (_, index) => ({ uid: `record-${index}` })),
  linearComparisonResolution: resolved({ sources: ['losat', 'losat', 'losat', 'losat'] }),
  losat: {
    totalThreadBudget: 'safe',
    threadsPerJob: '32',
    parallelWorkers: undefined,
    blastp: {
      mode: 'orthogroup',
      collinearSearchScope: 'adjacent'
    }
  },
  losatProgram: { value: 'blastn' }
};

createLosatSettings({ state: loadingOrderState });
assert.equal(loadingOrderState.losat.threadsPerJob, '32');

const noneState = {
  linearSeqs: Array.from({ length: 3 }, (_, index) => ({ uid: `record-${index}` })),
  linearComparisonResolution: resolved({ mode: 'none', sources: [] }),
  losat: {
    totalThreadBudget: '17',
    threadsPerJob: '32',
    parallelWorkers: '9',
    blastp: { mode: 'pairwise', collinearSearchScope: 'adjacent' }
  },
  losatProgram: { value: 'blastn' }
};
const noneSettings = createLosatSettings({ state: noneState });
assert.equal(noneSettings.losatEstimatedJobCount.value, 0);
assert.equal(noneSettings.losatMaxPairWorkers.value, 0);
assert.equal(noneSettings.losatAutoPairWorkers.value, 0);
assert.deepEqual(noneSettings.losatPairWorkerOptions.value, []);
assert.equal(noneState.losat.totalThreadBudget, '17');
assert.equal(noneState.losat.threadsPerJob, '32');
assert.equal(noneState.losat.parallelWorkers, '9');

const mixedState = {
  linearSeqs: Array.from({ length: 4 }, (_, index) => ({ uid: `record-${index}` })),
  linearComparisonResolution: resolved({
    mode: 'selected',
    defaultSource: 'upload',
    sources: ['losat', 'upload', 'losat']
  }),
  losat: {
    totalThreadBudget: 'safe',
    threadsPerJob: 'auto',
    parallelWorkers: undefined,
    blastp: { mode: 'pairwise', collinearSearchScope: 'adjacent' }
  },
  losatProgram: { value: 'blastn' }
};
const mixedSettings = createLosatSettings({ state: mixedState });
assert.equal(mixedSettings.losatEstimatedJobCount.value, 2);
mixedState.losatProgram.value = 'blastp';
assert.equal(mixedSettings.losatEstimatedJobCount.value, 2);

const expansionState = {
  linearSeqs: Array.from({ length: 5 }, (_, index) => ({ uid: `record-${index}` })),
  linearComparisonResolution: resolved({ sources: ['losat', 'losat', 'losat', 'losat'] }),
  losat: {
    totalThreadBudget: 'safe',
    threadsPerJob: 'auto',
    parallelWorkers: undefined,
    blastp: { mode: 'orthogroup', collinearSearchScope: 'all' }
  },
  losatProgram: { value: 'blastp' }
};
const expansionSettings = createLosatSettings({ state: expansionState });
assert.equal(
  expansionSettings.losatEstimatedJobCount.value,
  25,
  'Similarity groups must plan all five self and ordered cross-record searches'
);
expansionState.losat.blastp.mode = 'collinear';
assert.equal(
  expansionSettings.losatEstimatedJobCount.value,
  25,
  'fresh Collinear all-vs-all must plan all five self and ordered cross-record searches'
);
expansionState.losat.blastp.collinearSearchScope = 'adjacent';
assert.equal(
  expansionSettings.losatEstimatedJobCount.value,
  13,
  'an explicit adjacent Collinear scope must retain the narrower job plan'
);
expansionState.losat.blastp.collinearSearchScope = 'all';
assert.equal(expansionSettings.losatEstimatedJobCount.value, 25);

expansionState.linearComparisonResolution.value = {
  ...expansionState.linearComparisonResolution.value,
  valid: false
};
assert.equal(expansionSettings.losatEstimatedJobCount.value, 0);

console.log('losat-settings tests passed');

const sourceFiles = [{ name: 'upper.gb' }, { name: 'lower.gb' }];
expansionState.linearComparisonResolution.value.valid = true;
expansionState.linearSeqs.forEach((sequence, index) => { sequence.gb = sourceFiles[index < 2 ? 0 : 1]; });
expansionState.losat.blastp.mode = 'orthogroup';
assert.equal(expansionSettings.losatEstimatedJobCount.value, 4);
expansionState.losat.blastp.mode = 'pairwise';
expansionState.linearComparisonResolution = { value: {
  mode: 'adjacent', defaultSource: 'losat', valid: true, hasLosatIntent: true,
  edges: [0, 1].flatMap((queryIndex) => [2, 3, 4].map((subjectIndex) => ({
    source: 'losat', queryIndex, subjectIndex
  })))
} };
// The computed owner retains its original resolution ref.
const batchedSettings = createLosatSettings({ state: expansionState });
assert.equal(batchedSettings.losatEstimatedJobCount.value, 1);

// Four source jobs must use the selected total, rather than stay at two threads per job.
Object.defineProperty(globalThis, 'navigator', {
  value: { hardwareConcurrency: 64 }, configurable: true
});
const budgetState = {
  linearSeqs: [{ uid: 'a' }, { uid: 'b' }],
  linearComparisonResolution: resolved({ sources: ['losat'] }),
  losat: { totalThreadBudget: '32', threadsPerJob: 'auto', parallelWorkers: undefined,
    blastp: { mode: 'orthogroup', collinearSearchScope: 'all' } },
  losatProgram: { value: 'blastp' }
};
const budgetSettings = createLosatSettings({ state: budgetState });
for (const [budget, runs, threads] of [['32', 4, 8], ['16', 4, 4], ['2', 2, 1], ['safe', 4, 8], ['available', 4, 16]]) {
  budgetState.losat.totalThreadBudget = budget;
  assert.equal(budgetSettings.losatAutoPairWorkers.value, runs);
  assert.equal(budgetSettings.losatEffectiveThreadsPerJob.value, threads);
  assert(runs * threads <= budgetSettings.losatTotalThreadBudget.value);
  assert.equal(budgetState.losat.threadsPerJob, 'auto');
  assert.equal(budgetState.losat.parallelWorkers, undefined);
}
budgetState.losat.totalThreadBudget = '32';
budgetState.losat.parallelWorkers = '2';
assert.equal(budgetSettings.losatEffectiveThreadsPerJob.value, 16);
assert.equal(budgetSettings.losatMaxPairWorkers.value, 4);
budgetState.losat.threadsPerJob = '8';
budgetState.losat.totalThreadBudget = '4';
assert.equal(budgetSettings.losatEffectiveThreadsPerJob.value, 4);
assert(budgetSettings.losatThreadOptions.value.some(option => option.value === '8' && option.label.includes('4 effective')));
assert(budgetSettings.losatPairWorkerOptions.value.some(option => option.value === '2' && option.label.includes('1 effective')));
assert.equal(budgetState.losat.threadsPerJob, '8');
assert.equal(budgetState.losat.parallelWorkers, '2');
budgetState.losat.totalThreadBudget = '32';
assert.equal(budgetSettings.losatEffectiveThreadsPerJob.value, 8);
assert.equal(budgetSettings.losatMaxPairWorkers.value, 4);
for (const mode of ['blastn', 'tblastx']) {
  budgetState.losatProgram.value = mode;
  assert.equal(budgetSettings.losatEffectiveThreadsPerJob.value, 1);
  assert.equal(budgetState.losat.threadsPerJob, '8');
}
console.log('losat total-budget allocation tests passed');
