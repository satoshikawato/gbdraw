// Runs the app's Generate adopt port once with a Linear canonical Session over
// Circular rings that are off with a ring FASTA kept and no rows, and prints
// the rings before and after as JSON. The module
// mocks need `--experimental-test-module-mocks`, so a test runs this file in a
// child process (A2 review-I2, Medium 1).
import { File } from 'node:buffer';
import { mock } from 'node:test';
import { installFakeSvgDom } from '../fake-svg-dom.mjs';
import { modeDrawing } from './drawing-state.mjs';

const at = (path) => new URL(path, import.meta.url).href;
globalThis.window = {
  Vue: {
    ref: (value) => ({ value }), shallowRef: (value) => ({ value }), reactive: (value) => value,
    isRef: (value) => Boolean(value && typeof value === 'object' && 'value' in value),
    computed: (getter) => ({ get value() { return getter(); } }), nextTick: async () => {},
    watch: () => () => {}, onMounted: () => {}, onUnmounted: () => {}, provide: () => {}
  },
  DOMPurify: { sanitize: (value) => value }
};
globalThis.document = { createElement: () => ({ click: () => {} }) };
installFakeSvgDom();
globalThis.File = File;

const runAnalysisUrl = at('../../../gbdraw/web/js/app/run-analysis.js');
const configUrl = at('../../../gbdraw/web/js/services/config.js');
const real = await import(`${runAnalysisUrl}?real`);
const realConfig = await import(`${configUrl}?real`);
/** @type {{ adoptCanonicalRenderArtifacts: (canonical: unknown, options: unknown) => void }} */
let ports = /** @type {never} */ (null);
// A Linear adopt installs only the Linear comparisons and the committed Session.
mock.module(configUrl, { exports: { ...realConfig, adoptCanonicalRenderArtifacts: () => false } });
mock.module(runAnalysisUrl, {
  exports: { ...real, createRunAnalysis: (options) => { ports = options; return real.createRunAnalysis(options); } }
});
const { state } = await import(at('../../../gbdraw/web/js/state.js'));
const { createAppSetup } = await import(at('../../../gbdraw/web/js/app/app-setup.js'));
createAppSetup();
const rings = modeDrawing(state, 'circular').circularConservation;
state.files.c_conservation_fastas = [new File(['>r\nACGT\n'], 'ring.fa')];
rings.source = 'losat'; rings.labels = ''; rings.enabled = false; rings.series.splice(0);
const before = { enabled: rings.enabled, series: rings.series.length };
ports.adoptCanonicalRenderArtifacts({ renderRequest: { mode: 'linear', records: [] }, resources: {} }, { adoptOwnedRequest: false });
console.log(JSON.stringify({ before, after: { enabled: rings.enabled, series: rings.series.length } }));
