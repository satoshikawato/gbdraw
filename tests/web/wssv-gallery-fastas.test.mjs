import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';

import { buildRestoredMatchSequenceSources } from '../../gbdraw/web/js/app/match-sequences.js';
import { projectCanonicalSessionRequest } from '../../gbdraw/web/js/services/session-request.js';


const sessionPath = new URL(
  '../../gbdraw/web/gallery/sessions/WSSV_genome_comparison.gbdraw-session.json',
  import.meta.url
);
const session = JSON.parse(await readFile(sessionPath, 'utf8'));
const projection = projectCanonicalSessionRequest({
  renderRequest: session.renderRequest,
  resources: session.resources,
  webFiles: session.webFiles,
  legacyFiles: session.files,
  storedConfig: session.config,
  fileBindings: session.cliInvocation?.fileBindings
});

const expectedComparisons = [
  ['NC_003225.3', 309286],
  ['AF440570.1', 307287],
  ['NC_075105.1', 305119],
  ['AF369029.2', 292967],
  ['AP027278.1', 299976],
  ['AP027279.1', 293923],
  ['AP027284.1', 298496],
  ['AP027286.1', 289353],
  ['AP027288.1', 288252],
  ['KT995471.1', 284148],
  ['KY827813.1', 281054],
  ['MF768985.1', 285973],
  ['SRR14509867', 282445],
  ['SRR12919258', 289862],
  ['SRR17256726', 296593],
  ['Shantou2020', 284424],
  ['SRR8144089', 291020],
  ['SRR8144084', 289337],
  ['SRR22022264', 285171],
  ['ERR5659803', 283826]
];

assert.equal(projection.mode, 'circular');
assert.equal(projection.inputType, 'gb');
assert.equal(projection.files.c_conservation_blasts_source, 'losat-cache');
assert.equal(projection.files.c_conservation_fastas.length, 20);
assert.deepEqual(
  projection.files.c_conservation_fastas.map((file) => file.name),
  [
    'CN01.fasta',
    'WSSV-TW.fasta',
    'WSSV-CN.fasta',
    'WSSV-TH.fasta',
    'JP01A.fasta',
    'JP01B.fasta',
    'Pc2020.fasta',
    'E1.fasta',
    '0722-1.fasta',
    'CN03.fasta',
    'CN04.fasta',
    'WSSV-AU.fasta',
    'EU129.fa',
    'GCF7.fa',
    'MES-753.fa',
    'Shantou2019.fa',
    'POMZ1.fa',
    'POMZ4.fa',
    'MG18PR-0187-N40S.fa',
    'Angostura2013.fa'
  ]
);

const withText = (file) => ({
  ...file,
  text: async () => Buffer.from(file.data, 'base64').toString('utf8')
});
const restoredSources = await buildRestoredMatchSequenceSources({
  mode: projection.mode,
  cInputType: projection.inputType,
  files: {
    ...projection.files,
    c_gb: withText(projection.files.c_gb),
    c_conservation_fastas: projection.files.c_conservation_fastas.map(withText)
  },
  circularConservation: session.config.circularConservation
});
const comparisons = restoredSources.filter(
  (source) => source.origin === 'homology-comparison'
);

assert.equal(
  restoredSources.filter((source) => source.origin === 'circular-reference').length,
  1
);
assert.deepEqual(
  comparisons.map((source) => [source.recordId, source.sequence.length]),
  expectedComparisons
);
assert.deepEqual(
  comparisons.map((source) => source.sourceIndex),
  Array.from({ length: 20 }, (_, index) => index)
);
