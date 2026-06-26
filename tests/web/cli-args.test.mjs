import assert from 'node:assert/strict';
import { readFile, writeFile, mkdtemp } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';

const repoRoot = process.cwd();
const sourcePath = join(repoRoot, 'gbdraw', 'web', 'js', 'app', 'cli-args.js');
const tempDir = await mkdtemp(join(tmpdir(), 'gbdraw-cli-args-'));
const modulePath = join(tempDir, 'cli-args.mjs');
await writeFile(modulePath, await readFile(sourcePath, 'utf8'), 'utf8');

const {
  DEFAULT_CIRCULAR_CONSERVATION_BLAST_FILTERS,
  buildDefinitionLineStyleAssignments,
  buildBlastFilterArgs,
  buildFeatureShapeAssignments,
  buildRecordSelectorArgs,
  buildReverseComplementArgs,
  createDefaultLinearDefinitionLineStyles,
  isCliDefaultFeatureList,
  normalizeDefinitionLineStyleState
} = await import(pathToFileURL(modulePath));

assert.equal(
  isCliDefaultFeatureList(['CDS', 'rRNA', 'tRNA', 'tmRNA', 'ncRNA', 'misc_RNA', 'repeat_region']),
  true
);
assert.equal(
  isCliDefaultFeatureList(['CDS', 'rRNA', 'tRNA', 'tmRNA', 'ncRNA', 'repeat_region']),
  false
);

assert.deepEqual(
  buildFeatureShapeAssignments(
    ['CDS', 'rRNA', 'tRNA', 'tmRNA', 'ncRNA', 'repeat_region'],
    {
      CDS: 'arrow',
      rRNA: 'arrow',
      tRNA: 'arrow',
      tmRNA: 'arrow',
      ncRNA: 'arrow',
      repeat_region: 'rectangle'
    }
  ),
  []
);

assert.deepEqual(
  buildFeatureShapeAssignments(
    ['CDS', 'repeat_region', 'custom_feature'],
    {
      CDS: 'rectangle',
      repeat_region: 'arrow',
      custom_feature: 'arrow'
    }
  ),
  ['CDS=rectangle', 'repeat_region=arrow', 'custom_feature=arrow']
);

assert.deepEqual(createDefaultLinearDefinitionLineStyles(), {
  name: { font_size: null, font_weight: null, fill: null },
  subtitle: { font_size: null, font_weight: null, fill: null },
  replicon: { font_size: null, font_weight: null, fill: null },
  accession: { font_size: null, font_weight: null, fill: null },
  length: { font_size: null, font_weight: null, fill: null }
});

assert.deepEqual(
  normalizeDefinitionLineStyleState({
    name: { font_size: '12', font_weight: 'Bold', fill: '#111111' },
    accession: { font_size: '', font_weight: 'Default', fill: '' },
    replicon: { font_weight: 'normal' },
    length: { font_size: '9', font_weight: '700', fill: 'rgb(1,2,3)' }
  }),
  {
    name: { font_size: 12, font_weight: 'bold', fill: '#111111' },
    subtitle: { font_size: null, font_weight: null, fill: null },
    replicon: { font_size: null, font_weight: null, fill: null },
    accession: { font_size: null, font_weight: null, fill: null },
    length: { font_size: 9, font_weight: '700', fill: 'rgb(1,2,3)' }
  }
);

assert.deepEqual(
  buildDefinitionLineStyleAssignments({
    name: { font_size: 12, font_weight: 'bold', fill: '#111111' },
    accession: { font_weight: 'normal' },
    length: { fill: 'rgb(1,2,3)' }
  }),
  ['name:size=12,weight=bold,color=#111111', 'length:color=rgb(1,2,3)']
);

assert.deepEqual(
  buildBlastFilterArgs({
    bitscore: 50,
    evalue: '0.01',
    identity: 0,
    alignment_length: 0
  }),
  []
);

assert.deepEqual(
  buildBlastFilterArgs({
    bitscore: 50,
    evalue: '1e-2',
    identity: 0,
    alignment_length: 0
  }, DEFAULT_CIRCULAR_CONSERVATION_BLAST_FILTERS),
  ['--evalue', '1e-2', '--identity', 0]
);

assert.deepEqual(buildRecordSelectorArgs(['', '', '']), []);
assert.deepEqual(buildRecordSelectorArgs(['', 'NC_2', '']), ['--record_id', '', '--record_id', 'NC_2']);

assert.deepEqual(buildReverseComplementArgs([false, false]), []);
assert.deepEqual(
  buildReverseComplementArgs([false, true, false]),
  ['--reverse_complement', '0', '--reverse_complement', '1']
);
