import assert from 'node:assert/strict';
import { createHash } from 'node:crypto';
import { readFile } from 'node:fs/promises';

import { buildRestoredMatchSequenceSources } from '../../gbdraw/web/js/app/match-sequences.js';
import { readFileText } from '../../gbdraw/web/js/services/file-content-cache.js';
import {
  adoptCurrentSessionResources,
  sessionResourceSource
} from '../../gbdraw/web/js/services/session-resource-backing.js';
import { projectCanonicalSessionRequest } from '../../gbdraw/web/js/services/session-request.js';

const sessionPath = new URL(
  '../../gbdraw/web/gallery/sessions/WSSV_genome_comparison.gbdraw-session.json',
  import.meta.url
);
const session = JSON.parse(await readFile(sessionPath, 'utf8'));
const resourceMetrics = [];
globalThis.__GBDRAW_TEST_HOOKS__ = {
  onStructuralMetric(metric) {
    resourceMetrics.push({ ...metric });
  }
};
const sessionResourceTable = adoptCurrentSessionResources(session.resources);
const projection = projectCanonicalSessionRequest({
  renderRequest: session.renderRequest,
  resources: session.resources,
  webFiles: session.webFiles,
  legacyFiles: session.files,
  storedConfig: session.config,
  fileBindings: session.cliInvocation?.fileBindings,
  sessionResourceTable
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

const lazyViews = [projection.files.c_gb, ...projection.files.c_conservation_fastas];
for (const file of lazyViews) {
  assert.equal(Object.isFrozen(file), true);
  for (const field of ['text', 'arrayBuffer', 'data', 'resourceId']) {
    assert.equal(Object.hasOwn(file, field), false, `${file.name} exposed ${field}`);
  }
}
assert.equal(resourceMetrics.length, 0, 'resource projection must remain metadata-only');

const restoredSources = await buildRestoredMatchSequenceSources({
  mode: projection.mode,
  cInputType: projection.inputType,
  files: projection.files,
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

const expectedResourceIds = session.webFiles.conservationLosatFastaSources;
assert.deepEqual(
  projection.files.c_conservation_fastas.map((file) => (
    sessionResourceSource(file)?.resourceId
  )),
  expectedResourceIds
);
assert.equal(
  sessionResourceSource(projection.files.c_gb)?.resourceId,
  session.renderRequest.records[0].source.resourceId
);
assert.deepEqual(
  projection.files.c_conservation_fastas.map((file) => sessionResourceSource(file)?.descriptor),
  expectedResourceIds.map((resourceId) => session.resources[resourceId])
);

const fastaTexts = await Promise.all(
  projection.files.c_conservation_fastas.map((file) => readFileText(file))
);
const displayedCacheEntries = session.losatCache.entries.filter((entry) => (
  entry.flow === 'circular-conservation' && entry.display === true
));
assert.equal(displayedCacheEntries.length, 20);
assert.deepEqual(
  fastaTexts.map((text) => createHash('sha256').update(text).digest('hex')),
  displayedCacheEntries.map((entry) => entry.queryCanonicalHash)
);
assert.equal(
  resourceMetrics.filter(({ name }) => name === 'resourceTextReadCount').length,
  21
);
assert.equal(
  resourceMetrics.filter(({ name }) => name === 'base64DecodeCount').length,
  21
);
