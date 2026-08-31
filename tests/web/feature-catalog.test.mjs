import assert from 'node:assert/strict';
import test from 'node:test';

import {
  admitFeatureCatalog,
  biologicalFeatureKey,
  featureStateFromCatalog,
  stableFeatureOverrideKey,
  validateFeatureCatalog
} from '../../gbdraw/web/js/services/feature-catalog.js';

const results = [{ name: 'diagram.svg', content: '<svg />' }];
const fullNote = `${'x'.repeat(49)}😀tail`;
const compactNote = `${'x'.repeat(49)}😀`;
const catalog = {
  schema: 3,
  items: [{
    resultIndex: 0,
    resultName: 'diagram.svg',
    recordKeys: ['record-instance-a'],
    features: [{
      svgId: 'f0001',
      recordKey: 'record-instance-a',
      biologicalFeatureId: 'source-feature-a',
      fillColor: '#abcdef'
    }],
    biologicalFeatures: [{
      recordKey: 'record-instance-a',
      biologicalFeatureId: 'source-feature-a',
      stableFeatureId: 'stable-source-a',
      record_idx: 0,
      sourceFeatureIndex: 4,
      record_id: 'duplicated-accession',
      type: 'CDS',
      start: 1,
      end: 5,
      strand: -1,
      sequenceSourceIndex: 0,
      gene: 'curated-gene',
      amino_acid_sequence: 'MPEPTIDE',
      translationFromAminoAcidSequence: true,
      qualifiers: {
        protein_id: ['public-protein'],
        locus_tag: ['LOCUS_001'],
        old_locus_tag: ['OLD_LOCUS_001'],
        gene: ['derived-gene'],
        product: ['enzyme'],
        note: [fullNote]
      }
    }],
    orthogroups: [{
      id: 'group-1',
      members: [{
        recordKey: 'record-instance-a',
        biologicalFeatureId: 'source-feature-a',
        representative: true
      }]
    }],
    annotations: [],
    comparisonMatches: [],
    sequenceSources: [{
      key: 'linear:record:0',
      origin: 'linear-record',
      recordIndex: 0,
      sequence: 'AACCGG'
    }]
  }]
};

test('schema-3 catalog validates and expands stable biological identities', () => {
  const validated = validateFeatureCatalog(catalog, results);
  const state = featureStateFromCatalog(validated);
  const expectedKey = biologicalFeatureKey('record-instance-a', 'source-feature-a');

  assert.equal(state.extractedFeatures.length, 1);
  assert.equal(state.biologicalFeatures.length, 1);
  assert.equal(state.extractedFeatures[0].id, expectedKey);
  assert.equal(stableFeatureOverrideKey(state.extractedFeatures[0]), expectedKey);
  assert.equal(state.extractedFeatures[0].svg_id, 'f0001');
  assert.equal(state.extractedFeatures[0].stable_svg_id, 'stable-source-a');
  assert.equal(state.biologicalFeatures[0].svg_id, 'stable-source-a');
  assert.equal(state.biologicalFeatures[0].sourceFeatureIndex, 4);
  assert.equal(state.biologicalFeatures[0].feature_index, 4);
  assert.equal(state.biologicalFeatures[0].nucleotide_sequence, 'CGGT');
  assert.equal(state.biologicalFeatures[0].protein_id, 'public-protein');
  assert.equal(state.biologicalFeatures[0].locus_tag, 'LOCUS_001');
  assert.equal(state.biologicalFeatures[0].old_locus_tag, 'OLD_LOCUS_001');
  assert.equal(state.biologicalFeatures[0].gene, 'curated-gene');
  assert.equal(state.biologicalFeatures[0].product, 'enzyme');
  assert.equal(state.biologicalFeatures[0].note, compactNote);
  assert.deepEqual(
    state.biologicalFeatures[0].qualifiers.translation,
    ['MPEPTIDE']
  );
  assert.equal(
    state.biologicalFeatures[0].translationFromAminoAcidSequence,
    undefined
  );
  assert.deepEqual(state.biologicalFeatures[0].location_parts, [{
    start: 1,
    end: 5,
    strand: -1,
    display: '2..5'
  }]);
  assert.deepEqual(state.featureRecordIds, ['duplicated-accession']);
  assert.equal(state.orthogroups[0].members[0].renderedFeatureSvgId, 'f0001');
  assert.equal(state.orthogroups[0].members[0].proteinId, 'public-protein');
  assert.equal(state.orthogroups[0].members[0].featureIndex, 4);
});

test('duplicate-location source indexes survive catalog validation and expansion', () => {
  const duplicated = structuredClone(catalog);
  const second = {
    ...structuredClone(duplicated.items[0].biologicalFeatures[0]),
    biologicalFeatureId: 'source-feature-b',
    sourceFeatureIndex: 7
  };
  duplicated.items[0].biologicalFeatures.push(second);
  duplicated.items[0].features.push({
    svgId: 'f0002',
    recordKey: 'record-instance-a',
    biologicalFeatureId: 'source-feature-b',
    fillColor: '#fedcba'
  });
  duplicated.items[0].orthogroups[0].members.push({
    recordKey: 'record-instance-a',
    biologicalFeatureId: 'source-feature-b'
  });
  duplicated.items[0].comparisonMatches.push({
    id: 'duplicate-location-match',
    queryFeatureReferences: [{
      recordKey: 'record-instance-a',
      biologicalFeatureId: 'source-feature-a'
    }, {
      recordKey: 'record-instance-a',
      biologicalFeatureId: 'source-feature-b'
    }]
  });

  const state = featureStateFromCatalog(
    validateFeatureCatalog(duplicated, results)
  );

  assert.deepEqual(
    state.biologicalFeatures.map((feature) => feature.sourceFeatureIndex),
    [4, 7]
  );
  assert.deepEqual(
    state.extractedFeatures.map((feature) => feature.feature_index),
    [4, 7]
  );
  assert.deepEqual(
    state.orthogroups[0].members.map((member) => member.featureIndex),
    [4, 7]
  );

  const missingSourceIndex = structuredClone(duplicated);
  delete missingSourceIndex.items[0].biologicalFeatures[1].sourceFeatureIndex;
  assert.throws(
    () => validateFeatureCatalog(missingSourceIndex, results),
    /Reload the page and Generate again/
  );
});

test('catalog rejects invalid or conflicting source feature indexes', () => {
  for (const invalidValue of [-1, 0.5, '4', true]) {
    const invalid = structuredClone(catalog);
    invalid.items[0].biologicalFeatures[0].sourceFeatureIndex = invalidValue;
    assert.throws(
      () => validateFeatureCatalog(invalid, results),
      /Reload the page and Generate again/
    );
  }
  const conflicting = structuredClone(catalog);
  conflicting.items[0].biologicalFeatures[0].feature_index = 7;
  assert.throws(
    () => validateFeatureCatalog(conflicting, results),
    /Reload the page and Generate again/
  );
});

test('orthogroup members do not select one of multiple rendered references', () => {
  const repeatedRendering = structuredClone(catalog);
  repeatedRendering.items[0].features.push({
    svgId: 'f0001-alternate',
    recordKey: 'record-instance-a',
    biologicalFeatureId: 'source-feature-a',
    fillColor: '#fedcba'
  });

  const state = featureStateFromCatalog(
    validateFeatureCatalog(repeatedRendering, results)
  );
  const member = state.orthogroups[0].members[0];

  assert.deepEqual(
    state.extractedFeatures.map((feature) => feature.svg_id),
    ['f0001', 'f0001-alternate']
  );
  assert.equal(member.recordKey, 'record-instance-a');
  assert.equal(member.biologicalFeatureId, 'source-feature-a');
  assert.equal(member.stableFeatureSvgId, 'stable-source-a');
  assert.equal(member.renderedFeatureSvgId, '');
});

test('catalog presentation scope partitions local collinear groups without changing core scope', () => {
  const localCatalog = structuredClone(catalog);
  Object.assign(localCatalog.items[0].orthogroups[0], {
    scope: 'cross_record',
    presentationScope: 'adjacent_local',
    collinearGroupScope: 'adjacent_local',
    groupKind: 'collinear_gene_group'
  });

  const state = featureStateFromCatalog(
    validateFeatureCatalog(localCatalog, results)
  );

  assert.deepEqual(state.orthogroups, []);
  assert.equal(state.collinearGroups.length, 1);
  assert.equal(state.collinearGroups[0].scope, 'cross_record');
  assert.equal(state.collinearGroups[0].presentationScope, 'adjacent_local');
  assert.equal(state.collinearGroups[0].members[0].proteinId, 'public-protein');
});

test('catalog result topology must match the logical Results exactly', () => {
  assert.throws(
    () => validateFeatureCatalog(catalog, [{ ...results[0], name: 'other.svg' }]),
    /Reload the page and Generate again/
  );
  assert.throws(
    () => validateFeatureCatalog({ ...catalog, items: [] }, results),
    /Reload the page and Generate again/
  );
});

test('catalog rejects empty, duplicate, and normalized-colliding rendered IDs', () => {
  for (const svgId of ['', '   ']) {
    const malformed = structuredClone(catalog);
    malformed.items[0].features[0].svgId = svgId;
    assert.throws(
      () => validateFeatureCatalog(malformed, results),
      /Reload the page and Generate again/
    );
  }

  const duplicate = structuredClone(catalog);
  duplicate.items[0].features.push({
    ...duplicate.items[0].features[0]
  });
  assert.throws(
    () => validateFeatureCatalog(duplicate, results),
    /Reload the page and Generate again/
  );

  const normalizedCollision = structuredClone(catalog);
  normalizedCollision.items[0].features[0].svgId = 'f0001__part1';
  normalizedCollision.items[0].features.push({
    ...normalizedCollision.items[0].features[0],
    svgId: 'f0001'
  });
  assert.throws(
    () => validateFeatureCatalog(normalizedCollision, results),
    /Reload the page and Generate again/
  );
});

test('catalog rejects dangling rendered and group references', () => {
  const danglingRendered = structuredClone(catalog);
  danglingRendered.items[0].features[0].biologicalFeatureId = 'missing';
  assert.throws(
    () => validateFeatureCatalog(danglingRendered, results),
    /Reload the page and Generate again/
  );

  const danglingMember = structuredClone(catalog);
  danglingMember.items[0].orthogroups[0].members[0].recordKey = 'other';
  assert.throws(
    () => validateFeatureCatalog(danglingMember, results),
    /Reload the page and Generate again/
  );

  const danglingMatchGroup = structuredClone(catalog);
  danglingMatchGroup.items[0].comparisonMatches.push({
    id: 'match-with-missing-group',
    orthogroup_ids: ['missing-group']
  });
  assert.throws(
    () => validateFeatureCatalog(danglingMatchGroup, results),
    /Reload the page and Generate again/
  );
});

test('catalog validates ordered plural comparison endpoint references', () => {
  const plural = structuredClone(catalog);
  plural.items[0].features.push({
    svgId: 'f0002',
    recordKey: 'record-instance-a',
    biologicalFeatureId: 'source-feature-b'
  });
  plural.items[0].biologicalFeatures.push({
    recordKey: 'record-instance-a',
    biologicalFeatureId: 'source-feature-b',
    type: 'CDS',
    start: 5,
    end: 6,
    strand: 1
  });
  plural.items[0].comparisonMatches.push({
    id: 'plural-match',
    queryFeatureReferences: [{
      recordKey: 'record-instance-a',
      biologicalFeatureId: 'source-feature-a'
    }, {
      recordKey: 'record-instance-a',
      biologicalFeatureId: 'source-feature-b'
    }]
  });
  validateFeatureCatalog(plural, results);

  const duplicated = structuredClone(plural);
  duplicated.items[0].comparisonMatches[0].queryFeatureReferences[1] = {
    ...duplicated.items[0].comparisonMatches[0].queryFeatureReferences[0]
  };
  assert.throws(
    () => validateFeatureCatalog(duplicated, results),
    /Reload the page and Generate again/
  );

  const conflictingSingular = structuredClone(plural);
  Object.assign(conflictingSingular.items[0].comparisonMatches[0], {
    queryRecordKey: 'record-instance-a',
    queryBiologicalFeatureId: 'source-feature-a'
  });
  assert.throws(
    () => validateFeatureCatalog(conflictingSingular, results),
    /Reload the page and Generate again/
  );
});

test('catalog rejects duplicate match, member, and partial presentation metadata', () => {
  const duplicateMember = structuredClone(catalog);
  duplicateMember.items[0].orthogroups[0].members.push({
    ...duplicateMember.items[0].orthogroups[0].members[0]
  });
  assert.throws(
    () => validateFeatureCatalog(duplicateMember, results),
    /Reload the page and Generate again/
  );

  const duplicateMatch = structuredClone(catalog);
  duplicateMatch.items[0].comparisonMatches.push(
    { id: 'duplicate' },
    { id: 'duplicate' }
  );
  assert.throws(
    () => validateFeatureCatalog(duplicateMatch, results),
    /Reload the page and Generate again/
  );

  const partialPresentation = structuredClone(catalog);
  partialPresentation.items[0].orthogroups[0].presentationScope = 'adjacent_local';
  partialPresentation.items[0].orthogroups[0].groupKind = 'collinear_gene_group';
  assert.throws(
    () => validateFeatureCatalog(partialPresentation, results),
    /Reload the page and Generate again/
  );
});

test('catalog rejects invalid sequence-source references', () => {
  for (const invalidIndex of [-1, 1, 0.5, '0']) {
    const invalid = structuredClone(catalog);
    invalid.items[0].biologicalFeatures[0].sequenceSourceIndex = invalidIndex;
    assert.throws(
      () => validateFeatureCatalog(invalid, results),
      /Reload the page and Generate again/
    );
  }
  const invalidCoordinates = structuredClone(catalog);
  invalidCoordinates.items[0].biologicalFeatures[0].end = 99;
  assert.throws(
    () => validateFeatureCatalog(invalidCoordinates, results),
    /Reload the page and Generate again/
  );
  const wrongRecord = structuredClone(catalog);
  wrongRecord.items[0].sequenceSources[0].recordIndex = 1;
  assert.throws(
    () => validateFeatureCatalog(wrongRecord, results),
    /Reload the page and Generate again/
  );
  for (const conflict of [
    'DIFFERENT', '', ' ', null, [], [null], 0, false
  ]) {
    const typedConflict = structuredClone(catalog);
    typedConflict.items[0].biologicalFeatures[0]
      .qualifiers.translation = conflict;
    assert.throws(
      () => validateFeatureCatalog(typedConflict, results),
      /Reload the page and Generate again/
    );
  }
  const missingAminoAcid = structuredClone(catalog);
  delete missingAminoAcid.items[0].biologicalFeatures[0].amino_acid_sequence;
  assert.throws(
    () => validateFeatureCatalog(missingAminoAcid, results),
    /Reload the page and Generate again/
  );
  for (const invalidAminoAcid of [0, false, {}, []]) {
    const typedAminoAcid = structuredClone(catalog);
    typedAminoAcid.items[0].biologicalFeatures[0]
      .amino_acid_sequence = invalidAminoAcid;
    assert.throws(
      () => validateFeatureCatalog(typedAminoAcid, results),
      /Reload the page and Generate again/
    );
  }
  const camelAminoAcid = structuredClone(catalog);
  const camelFeature = camelAminoAcid.items[0].biologicalFeatures[0];
  camelFeature.aminoAcidSequence = camelFeature.amino_acid_sequence;
  delete camelFeature.amino_acid_sequence;
  validateFeatureCatalog(camelAminoAcid, results);
  for (const shadowingValue of [null, '']) {
    const shadowedAminoAcid = structuredClone(camelAminoAcid);
    shadowedAminoAcid.items[0].biologicalFeatures[0]
      .amino_acid_sequence = shadowingValue;
    assert.throws(
      () => validateFeatureCatalog(shadowedAminoAcid, results),
      /Reload the page and Generate again/
    );
  }
  for (const invalidSequence of [123, {}, 'AAC CGG']) {
    const typedSequence = structuredClone(catalog);
    typedSequence.items[0].sequenceSources[0].sequence = invalidSequence;
    assert.throws(
      () => validateFeatureCatalog(typedSequence, results),
      /Reload the page and Generate again/
    );
  }
  const invalidUnreferencedSource = structuredClone(catalog);
  invalidUnreferencedSource.items[0].sequenceSources.push({
    origin: 'linear-record',
    recordIndex: 0,
    sequence: 'AAC CGG'
  });
  assert.throws(
    () => validateFeatureCatalog(invalidUnreferencedSource, results),
    /Reload the page and Generate again/
  );
  const coexistingSequence = structuredClone(catalog);
  coexistingSequence.items[0].biologicalFeatures[0]
    .nucleotide_sequence = 'CGGT';
  assert.throws(
    () => validateFeatureCatalog(coexistingSequence, results),
    /Reload the page and Generate again/
  );
});

test('inline qualifier differences remain authoritative', () => {
  const inline = structuredClone(catalog);
  const biological = inline.items[0].biologicalFeatures[0];
  delete biological.translationFromAminoAcidSequence;
  biological.qualifiers.translation = ['DIFFERENT'];
  const state = featureStateFromCatalog(validateFeatureCatalog(inline, results));

  assert.equal(state.biologicalFeatures[0].gene, 'curated-gene');
  assert.deepEqual(
    state.biologicalFeatures[0].qualifiers.translation,
    ['DIFFERENT']
  );
});

test('whole-record source validation is bounded by source count', () => {
  const sourceSequence = `ATG${'A'.repeat(1_000_000)}`;
  const featureCount = 512;
  const largeCatalog = {
    schema: 3,
    items: [{
      resultIndex: 0,
      resultName: 'diagram.svg',
      recordKeys: ['record-key'],
      features: [],
      biologicalFeatures: Array.from(
        { length: featureCount },
        (_, index) => ({
          recordKey: 'record-key',
          biologicalFeatureId: `feature-${index}`,
          type: 'CDS',
          start: 0,
          end: 3,
          strand: '+',
          sequenceSourceIndex: 0
        })
      ),
      orthogroups: [],
      annotations: [],
      comparisonMatches: [],
      sequenceSources: [{
        origin: 'linear-record',
        recordIndex: 0,
        sequence: sourceSequence
      }]
    }]
  };
  const nativeTest = RegExp.prototype.test;
  let sourceScans = 0;
  RegExp.prototype.test = function (value) {
    if (this.source === '\\s' && value === sourceSequence) sourceScans += 1;
    return nativeTest.call(this, value);
  };
  try {
    const validated = validateFeatureCatalog(largeCatalog, results);
    assert.equal(sourceScans, 1);
    const state = featureStateFromCatalog(validated);
    assert.equal(state.biologicalFeatures.length, featureCount);
    assert.equal(state.biologicalFeatures[0].nucleotide_sequence, 'ATG');
    assert.equal(sourceScans, 1);
  } finally {
    RegExp.prototype.test = nativeTest;
  }
});

test('admitted projections do not reread later catalog mutations', () => {
  const initialSequence = `ATG${'A'.repeat(1_000_000)}`;
  const updatedSequence = `CCC${'A'.repeat(1_000_000)}`;
  const largeCatalog = {
    schema: 3,
    items: [{
      resultIndex: 0,
      resultName: 'diagram.svg',
      recordKeys: ['record-key'],
      features: [],
      biologicalFeatures: [{
        recordKey: 'record-key',
        biologicalFeatureId: 'feature',
        type: 'CDS',
        start: 0,
        end: 3,
        strand: '+',
        sequenceSourceIndex: 0
      }],
      orthogroups: [],
      annotations: [],
      comparisonMatches: [],
      sequenceSources: [{
        origin: 'linear-record',
        recordIndex: 0,
        sequence: initialSequence
      }]
    }]
  };
  const sourceSequences = new Set([initialSequence, updatedSequence]);
  const nativeTest = RegExp.prototype.test;
  let sourceScans = 0;
  RegExp.prototype.test = function (value) {
    if (this.source === '\\s' && sourceSequences.has(value)) sourceScans += 1;
    return nativeTest.call(this, value);
  };
  try {
    const validated = validateFeatureCatalog(largeCatalog, results);
    assert.equal(
      featureStateFromCatalog(validated)
        .biologicalFeatures[0].nucleotide_sequence,
      'ATG'
    );
    validated.items[0].sequenceSources[0].sequence = updatedSequence;
    assert.equal(
      featureStateFromCatalog(validated)
        .biologicalFeatures[0].nucleotide_sequence,
      'ATG'
    );
  } finally {
    RegExp.prototype.test = nativeTest;
  }
  assert.equal(sourceScans, 1);
});

test('compound identity isolates duplicate accessions and positional indices', () => {
  const otherRecord = biologicalFeatureKey('record-instance-b', 'source-feature-a');
  assert.notEqual(
    biologicalFeatureKey('record-instance-a', 'source-feature-a'),
    otherRecord
  );
});

test('schema-3 compact defaults restore stable IDs and record indexes', () => {
  const compactCatalog = structuredClone(catalog);
  const biological = compactCatalog.items[0].biologicalFeatures[0];
  delete biological.stableFeatureId;
  delete biological.record_idx;

  const state = featureStateFromCatalog(
    validateFeatureCatalog(compactCatalog, results)
  );
  assert.equal(state.biologicalFeatures[0].stable_feature_id, 'source-feature-a');
  assert.equal(state.biologicalFeatures[0].record_idx, 0);
  assert.equal(state.extractedFeatures[0].stable_svg_id, 'source-feature-a');
  assert.equal(state.orthogroups[0].members[0].recordIndex, 0);
});

test('linear catalog state uses global record indexes and matching display labels', () => {
  const multiRecordCatalog = structuredClone(catalog);
  multiRecordCatalog.items[0].recordKeys.push('record-instance-b');
  multiRecordCatalog.items[0].features.push({
    svgId: 'f0002',
    recordKey: 'record-instance-b',
    biologicalFeatureId: 'source-feature-b'
  });
  multiRecordCatalog.items[0].biologicalFeatures.push({
    recordKey: 'record-instance-b',
    biologicalFeatureId: 'source-feature-b',
    record_id: 'second-accession',
    type: 'CDS',
    start: 60,
    end: 90,
    strand: 1
  });

  const state = featureStateFromCatalog(
    validateFeatureCatalog(multiRecordCatalog, results, { mode: 'linear' }),
    { mode: 'linear' }
  );

  assert.deepEqual(state.featureRecordIds, [
    'File 1: duplicated-accession',
    'File 2: second-accession'
  ]);
  assert.deepEqual(
    state.extractedFeatures.map((feature) => ({
      recordIdx: feature.record_idx,
      fileIdx: feature.fileIdx,
      displayRecordId: feature.displayRecordId
    })),
    [{
      recordIdx: 0,
      fileIdx: 0,
      displayRecordId: 'File 1: duplicated-accession'
    }, {
      recordIdx: 1,
      fileIdx: 1,
      displayRecordId: 'File 2: second-accession'
    }]
  );
});

test('circular batch catalog state rebases item-local record indexes globally', () => {
  const batchResults = [
    { name: 'first.svg', content: '<svg />' },
    { name: 'second.svg', content: '<svg />' }
  ];
  const firstItem = structuredClone(catalog.items[0]);
  firstItem.resultName = 'first.svg';
  delete firstItem.biologicalFeatures[0].record_idx;
  const secondItem = {
    resultIndex: 1,
    resultName: 'second.svg',
    recordKeys: ['record-instance-b'],
    features: [{
      svgId: 'f0002',
      recordKey: 'record-instance-b',
      biologicalFeatureId: 'source-feature-b'
    }],
    biologicalFeatures: [{
      recordKey: 'record-instance-b',
      biologicalFeatureId: 'source-feature-b',
      record_id: 'second-accession',
      type: 'CDS',
      start: 60,
      end: 90,
      strand: 1
    }],
    orthogroups: [],
    annotations: [],
    comparisonMatches: []
  };
  const batchCatalog = { schema: 3, items: [firstItem, secondItem] };

  const state = featureStateFromCatalog(
    validateFeatureCatalog(batchCatalog, batchResults, { mode: 'circular' }),
    { mode: 'circular' }
  );

  assert.deepEqual(state.featureRecordIds, [
    'duplicated-accession',
    'second-accession'
  ]);
  assert.deepEqual(
    state.extractedFeatures.map((feature) => feature.record_idx),
    [0, 1]
  );
  assert.deepEqual(
    state.extractedFeatures.map((feature) => feature.displayRecordId),
    ['duplicated-accession', 'second-accession']
  );
});

test('one catalog admission produces identities, indexes, scalar footprint, and zero secondary traversal', () => {
  const metrics = [];
  const events = [];
  globalThis.__GBDRAW_TEST_HOOKS__ = {
    onStructuralMetric: (metric) => metrics.push(metric),
    onSessionLifecycleEvent: (event) => events.push(event)
  };
  try {
    const source = structuredClone(catalog);
    const beforeAdmission = structuredClone(source);
    const admission = admitFeatureCatalog(source, results, {
      adopt: true,
      mode: 'linear'
    });
    assert.equal(admission.catalog, source);
    assert.deepEqual(source, beforeAdmission);
    assert.equal(admission.featureState.extractedFeatures.length, 1);
    assert.equal(admission.featureOrthogroupIndex.size > 0, true);
    assert.equal(admission.renderedIdentitiesByResult.length, 1);
    assert.equal(
      admission.renderedIdentitiesByResult[0].byRenderedId.get('f0001').stableId,
      'stable-source-a'
    );
    assert.deepEqual(admission.scalarMetrics, {
      resultCount: 1,
      itemCount: 1,
      recordCount: 1,
      renderedFeatureCount: 1,
      biologicalFeatureCount: 1,
      orthogroupRecordCount: 1,
      comparisonMatchCount: 0,
      annotationCount: 0,
      sequenceSourceCount: 1,
      sequenceCharacters: 6
    });
    assert.equal(
      metrics.filter((metric) => metric.name === 'featureCatalogAdmissionCount')
        .reduce((total, metric) => total + metric.value, 0),
      1
    );
    assert.equal(
      metrics.filter((metric) => metric.name === 'featureCatalogSecondaryTraversalCount')
        .reduce((total, metric) => total + metric.value, 0),
      0
    );
    assert.deepEqual(
      events.filter((event) => event.name.startsWith('catalog.')).map((event) => event.name),
      ['catalog.admission-started', 'catalog.admission-completed']
    );
  } finally {
    delete globalThis.__GBDRAW_TEST_HOOKS__;
  }
});
