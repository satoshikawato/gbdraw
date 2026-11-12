// Circular similarity rings search with one E-value scale per genome (PD-OI-018 revision 6, OV-3319): the Web's
// Worker plans each ring's -searchsp with the Python owner the command line's ring search reads, from the reference
// staged as the Generate's request names it and the lengths the ring reader answered, and the ring's LOSAT search
// and its raw cache entry pass it. Real Worker, real Pyodide, real LOSAT Wasm.
const { test, expect } = require('@playwright/test');
const { writeFileSync } = require('node:fs');
const { gzipSync } = require('node:zlib');
const { evaluateWithRetainedPromise } = require('./helpers/app-lifecycle.cjs');
const { openWithGenBank } = require('./helpers/audit-browser.cjs');

/** A reproducible sequence of `length` bases. @param {number} seed @param {number} length */
const bases = (seed, length) => {
  let state = seed;
  return Array.from({ length }, () => {
    state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
    return 'acgt'[state >>> 30];
  }).join('');
};

/** A GenBank record of `sequence`. @param {string} id @param {string} sequence */
const genbankRecord = (id, sequence) => [
  `LOCUS       ${id.padEnd(16)}${String(sequence.length).padStart(11)} bp    DNA     circular UNK 01-JAN-1980`,
  `DEFINITION  ${id} synthetic record.`,
  `ACCESSION   ${id}`,
  `VERSION     ${id}`,
  'FEATURES             Location/Qualifiers',
  `     source          1..${sequence.length}`,
  'ORIGIN',
  ...Array.from({ length: Math.ceil(sequence.length / 60) }, (_, line) => `${String(line * 60 + 1).padStart(9)} ${
    sequence.slice(line * 60, line * 60 + 60).match(/.{1,10}/g).join(' ')}`),
  '//'
].join('\n');

// A gzip reference is staged for the plan as Generate's request names it (its gzip bytes, PD-OI-121 revision 3); Python
// expands it, so the search space is that of the expanded records.
for (const gzip of [false, true]) test(`a ring of two files searches with the search space of the genome pair Python plans${
  gzip ? ', from a gzip reference' : ''}`, async ({ page }, testInfo) => {
  test.setTimeout(300_000);
  // The reference: one file of two records (360 + 240 bp); the ring: a chromosome file and a plasmid file.
  const chromosome = bases(977, 360);
  const plasmid = bases(431, 240);
  const reference = testInfo.outputPath(gzip ? 'reference.gbk.gz' : 'reference.gbk');
  const text = `${genbankRecord('REFCHR', chromosome)}\n${genbankRecord('REFPLS', plasmid)}\n`;
  writeFileSync(reference, gzip ? gzipSync(text) : text);
  await openWithGenBank(page, reference);
  await page.evaluate(() => {
    window.__GBDRAW_RING_ARGS__ = [];
    window.__GBDRAW_LOSAT_EXECUTOR__ = async (jobs, options) => {
      window.__GBDRAW_RING_ARGS__.push(jobs.map((job) => job.extraArgs));
      const { runLosatPairsParallel } = await import('/gbdraw/web/js/services/losat.js');
      return runLosatPairsParallel(jobs, options);
    };
  });
  await page.evaluate(async ({ ringChromosome, ringPlasmid }) => {
    const app = window.__GBDRAW_APP__;
    Object.assign(app.circularConservation, { enabled: true, source: 'losat', losat_program: 'blastn' });
    Object.assign(app.adv, { min_bitscore: 0, evalue: 10, identity: 0, alignment_length: 0 });
    // The serial engine reads -searchsp through LOSAT's Web ABI, and needs no cross-origin isolation.
    app.losatExecution.executionMode = 'serial';
    const file = (name, text) => new File([text], name, { type: 'text/plain', lastModified: 1 });
    app.addCircularConservationComparisonFile({ target: { files: [file('ring-chromosome.fasta', ringChromosome)], value: '' } });
    app.addRingRowFiles(0, { target: { files: [file('ring-plasmid.fasta', ringPlasmid)], value: '' } });
    await window.Vue.nextTick();
  }, {
    ringChromosome: `>RINGCHR\n${chromosome.slice(0, 300)}\n`,
    ringPlasmid: `>RINGPLS\n${plasmid.slice(0, 120)}\n`
  });
  const result = await evaluateWithRetainedPromise(page, () => window.__GBDRAW_APP__.runAnalysis());
  expect(result, JSON.stringify(await page.evaluate(() => window.__GBDRAW_APP__.errorLog))).toEqual({ status: 'ok' });

  // The reference genome is both records (600 bp), the ring's genome both of its files (420 bp).
  const space = String((360 + 240) * (300 + 120));
  const [[searched]] = await page.evaluate(() => window.__GBDRAW_RING_ARGS__);
  expect(searched).toEqual(['--task', 'megablast', '--searchsp', space]);
  const stored = await page.evaluate(async () => {
    const { state } = await import('/gbdraw/web/js/state.js');
    return [...state.losatCache.value.values()].filter((entry) => entry.flow === 'circular-conservation')
      .map((entry) => entry.args);
  });
  expect(stored).toEqual([['--task', 'megablast', '--searchsp', space]]);
});
