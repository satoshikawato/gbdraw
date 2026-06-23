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
);
await writeFile(
  tempNormalizationPath,
  await readFile(new URL('losat-normalization.js', sourceRoot), 'utf8')
);

const { createLosatSettings } = await import(pathToFileURL(tempModulePath));

const state = {
  linearSeqs: [{}, {}, {}, {}, {}],
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
  linearSeqs: [{}, {}, {}, {}, {}],
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

console.log('losat-settings tests passed');
