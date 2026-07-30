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
  SPECIFIC_COLOR_QUALIFIER_PRESETS,
  buildFeatureSelectorUniquenessIndex,
  buildSelectorSafetyUniquenessIndex,
  collectSpecificColorQualifierSuggestions,
  exactRegexValue,
  resolveFeatureLabelSelector,
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
  const feature = {
    svg_id: 'ff51a6081_record_2',
    stable_svg_id: 'ff51a6081',
    record_id: 'LC921558.1',
    type: 'mRNA',
    selector: {
      hash: 'ff51a6081',
      record_location: 'LC921558.1:775..13422:+',
      qualifiers: { gene: ['penF'] }
    }
  };
  const scope = [
    { record_id: 'LC921558.1', feature_type: 'mRNA', selector: feature.selector },
    {
      record_id: 'LC921558.1',
      feature_type: 'mRNA',
      selector: {
        hash: 'f8468d457',
        record_location: 'LC921558.1:775..13422:+',
        qualifiers: { gene: ['penF'] }
      }
    }
  ];
  const selected = selectFeatureSelector(
    feature,
    buildSelectorSafetyUniquenessIndex(scope),
    { priority: ['protein_id', 'locus_tag'], requireSelector: true, requireSafetyScope: true }
  );
  assert.deepEqual(selected, { qualifier: 'hash', value: 'ff51a6081', isFallbackHash: true });
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

{
  const feature = makeFeature({
    selector: {
      hash: 'h1',
      qualifiers: {
        product: ['wsv360-like protein'],
        gene: ['wsv360']
      }
    }
  });
  assert.deepEqual(resolveFeatureLabelSelector(feature, 'wsv360-like protein'), {
    qualifier: 'product',
    value: 'wsv360-like protein',
    pattern: '^wsv360-like protein$'
  });
}

{
  const feature = makeFeature({
    qualifiers: {
      product: ['  ATPase (A)+ [x]?  ']
    },
    selector: {
      hash: 'h1',
      qualifiers: {
        gene: ['atpA']
      }
    }
  });
  assert.equal(
    resolveFeatureLabelSelector(feature, 'ATPase (A)+ [x]?')?.pattern,
    '^  ATPase \\(A\\)\\+ \\[x\\]\\?  $'
  );
  assert.equal(resolveFeatureLabelSelector(feature, 'manually edited label'), null);
}

{
  const suggestions = collectSpecificColorQualifierSuggestions(
    [
      { qualifiers: { custom_tag: ['one'], color: ['#ff0000'] } },
      { selector: { qualifiers: { VendorKey: ['two'] } } }
    ],
    [{ qual: 'rule_only' }, { qual: ' custom_tag ' }]
  );
  assert.equal(SPECIFIC_COLOR_QUALIFIER_PRESETS.includes('color'), true);
  assert.deepEqual(suggestions.slice(0, SPECIFIC_COLOR_QUALIFIER_PRESETS.length), SPECIFIC_COLOR_QUALIFIER_PRESETS);
  assert.equal(suggestions.includes('custom_tag'), true);
  assert.equal(suggestions.includes('VendorKey'), true);
  assert.equal(suggestions.includes('rule_only'), true);
  assert.equal(suggestions.filter((value) => value === 'custom_tag').length, 1);
}

console.log('feature selector tests passed');
