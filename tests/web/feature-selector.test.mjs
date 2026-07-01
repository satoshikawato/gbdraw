import assert from 'node:assert/strict';
import { readFile, writeFile, mkdtemp } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';

const repoRoot = process.cwd();
const sourcePath = join(repoRoot, 'gbdraw', 'web', 'js', 'app', 'feature-selector.js');
const tempDir = await mkdtemp(join(tmpdir(), 'gbdraw-feature-selector-'));
await writeFile(join(tempDir, 'package.json'), '{"type":"module"}\n', 'utf8');
await writeFile(
  join(tempDir, 'feature-selector.js'),
  await readFile(sourcePath, 'utf8'),
  'utf8'
);

const {
  buildFeatureSelectorUniquenessIndex,
  buildSelectorSafetyUniquenessIndex,
  exactRegexValue,
  selectFeatureSelector
} = await import(pathToFileURL(join(tempDir, 'feature-selector.js')));

const makeFeature = (overrides = {}) => ({
  svg_id: 'h1',
  record_id: 'rec1',
  type: 'CDS',
  selector: {
    hash: 'h1',
    record_location: 'rec1:10..30:+',
    qualifiers: {
      protein_id: ['P1'],
      locus_tag: ['L1']
    }
  },
  ...overrides
});

{
  const scope = [
    makeFeature().selector,
  ].map((selector) => ({ record_id: 'rec1', feature_type: 'CDS', selector }));
  const index = buildSelectorSafetyUniquenessIndex(scope);
  const selected = selectFeatureSelector(
    makeFeature(),
    index,
    { priority: ['protein_id', 'locus_tag'], requireSelector: true, requireSafetyScope: true }
  );
  assert.deepEqual(selected, { qualifier: 'protein_id', value: 'P1', isFallbackHash: false });
}

{
  const feature = makeFeature();
  const scope = [
    { record_id: 'rec1', feature_type: 'CDS', selector: feature.selector },
    {
      record_id: 'rec1',
      feature_type: 'CDS',
      selector: {
        hash: 'h2',
        record_location: 'rec1:40..60:+',
        qualifiers: { protein_id: ['P1'], locus_tag: ['L2'] }
      }
    }
  ];
  const selected = selectFeatureSelector(
    feature,
    buildSelectorSafetyUniquenessIndex(scope),
    { priority: ['protein_id', 'locus_tag'], requireSelector: true, requireSafetyScope: true }
  );
  assert.deepEqual(selected, { qualifier: 'locus_tag', value: 'L1', isFallbackHash: false });
}

{
  const feature = makeFeature({
    selector: {
      hash: 'h1',
      record_location: 'rec1:10..30:+',
      qualifiers: { protein_id: ['P1'] }
    }
  });
  const scope = [
    { record_id: 'rec1', feature_type: 'CDS', selector: feature.selector },
    {
      record_id: 'rec1',
      feature_type: 'CDS',
      selector: {
        hash: 'h2',
        record_location: 'rec1:40..60:+',
        qualifiers: { protein_id: ['P1'] }
      }
    }
  ];
  const selected = selectFeatureSelector(
    feature,
    buildSelectorSafetyUniquenessIndex(scope),
    { priority: ['protein_id', 'locus_tag'], requireSelector: true, requireSafetyScope: true }
  );
  assert.deepEqual(selected, {
    qualifier: 'record_location',
    value: 'rec1:10..30:+',
    isFallbackHash: false
  });
}

{
  const selected = selectFeatureSelector(
    makeFeature(),
    buildSelectorSafetyUniquenessIndex(null),
    { priority: ['protein_id'], requireSelector: true, requireSafetyScope: true }
  );
  assert.deepEqual(selected, { qualifier: 'hash', value: 'h1', isFallbackHash: true });
}

{
  const feature = makeFeature({
    selector: {
      hash: 'h1',
      qualifiers: {
        product: ['specific product'],
        gene: ['geneA'],
        note: ['note A'],
        id: ['qualifier-id'],
        name: ['qualifier-name']
      }
    }
  });
  const scope = [{ record_id: 'rec1', feature_type: 'CDS', selector: feature.selector }];
  const selected = selectFeatureSelector(
    feature,
    buildSelectorSafetyUniquenessIndex(scope),
    { priority: ['protein_id', 'locus_tag', 'gene_id', 'old_locus_tag'], requireSelector: true, requireSafetyScope: true }
  );
  assert.deepEqual(selected, { qualifier: 'hash', value: 'h1', isFallbackHash: true });
}

{
  const feature = makeFeature({
    selector: {
      hash: 'h1',
      record_location: 'rec1:10..30:+',
      qualifiers: { protein_id: ['P1'], gene: ['geneA'] }
    }
  });
  const index = buildFeatureSelectorUniquenessIndex([feature], { preferSelector: false });
  assert.equal(
    selectFeatureSelector(feature, index, { priority: ['gene'], preferSelector: false }).qualifier,
    'gene'
  );
  assert.equal(
    selectFeatureSelector(feature, index, { priority: ['protein_id'], preferSelector: false }).qualifier,
    'protein_id'
  );
}

assert.equal(exactRegexValue('YP_009725295.1'), '^YP_009725295\\.1$');

console.log('feature selector tests passed');
