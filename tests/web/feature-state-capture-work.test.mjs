// J5 (perf 0.14.x, counts not timings): a capture of the feature and
// orthogroup state (a History checkpoint, twice per Reset Settings; a Session
// save; an import rollback) reads the catalog features and groups from their
// raw arrays, so its reactive reads do not grow with the 10k-40k catalog rows,
// and the copy holds the same values.
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { test } from 'node:test';
import vm from 'node:vm';

const context = vm.createContext({ console });
vm.runInContext(readFileSync(new URL('../../gbdraw/web/vendor/vue/vue.global.js', import.meta.url), 'utf8'), context);
globalThis.window = { Vue: context.Vue, setTimeout, clearTimeout };
const { effect, toRaw } = context.Vue;
const { state } = await import('../../gbdraw/web/js/state.js');
const { buildFeatureStateData, buildOrthogroupStateData } = await import('../../gbdraw/web/js/services/config.js');

const feature = (index) => ({
  svg_id: `f${index}`, type: 'CDS', start: index, end: index + 9, strand: 1, nucleotide_sequence: 'ATG',
  qualifiers: { locus_tag: [`L${index}`], product: [`p${index}`] }, selector: { hash: `h${index}`, qualifiers: {} }
});
const capture = (count) => {
  state.extractedFeatures.value = Array.from({ length: count }, (_, index) => feature(index));
  state.biologicalFeatures.value = Array.from({ length: count }, (_, index) => feature(index));
  state.orthogroups.value = Array.from({ length: count }, (_, index) => ({ id: `og${index}`, members: [`f${index}`] }));
  const drawing = state.activeDrawing();
  /** @type {Record<string, any>} */
  const copy = {};
  // The reactive values the capture reads: the dependencies of an effect that runs it.
  const runner = effect(() => {
    copy.features = buildFeatureStateData(drawing);
    copy.orthogroups = buildOrthogroupStateData(drawing);
  });
  let reads = 0;
  for (let link = runner.effect.deps; link; link = link.nextDep) reads += 1;
  runner.effect.stop();
  return { reads, copy };
};

test('the reactive reads of a feature state capture do not grow with the catalog', () => {
  const small = capture(20);
  const large = capture(200);
  assert.equal(large.reads, small.reads, `${small.reads} reads for 20 rows, ${large.reads} for 200`);
});

test('the capture holds the catalog values without their sequences', () => {
  const { copy } = capture(3);
  const { nucleotide_sequence: _sequence, ...expected } = feature(1);
  assert.equal(JSON.stringify(copy.features.extractedFeatures[1]), JSON.stringify(expected));
  assert.equal(JSON.stringify(copy.features.biologicalFeatures[2].qualifiers), JSON.stringify(feature(2).qualifiers));
  assert.equal(JSON.stringify(copy.orthogroups.groups), JSON.stringify(toRaw(state.orthogroups.value)));
});
