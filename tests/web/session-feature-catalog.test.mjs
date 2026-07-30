import assert from 'node:assert/strict';
import test from 'node:test';

import {
  compactSessionFeatureCatalog,
  expandSessionFeatureCatalog
} from '../../gbdraw/web/js/services/session-feature-catalog.js';

const richFeatureState = () => {
  const longNote = `${'😀'.repeat(30)}${'x'.repeat(30)}`;
  const qualifiers = {
    protein_id: ['protein-1'],
    locus_tag: ['locus-1'],
    product: ['example protein'],
    translation: ['M'],
    note: [longNote]
  };
  const biological = {
    id: 'f0',
    svg_id: 'f_stable',
    stable_svg_id: 'f_stable',
    stable_feature_id: 'f_stable',
    record_id: 'record',
    record_idx: 0,
    feature_index: 0,
    organism: 'Example',
    type: 'CDS',
    start: 1,
    end: 3,
    strand: '+',
    protein_id: 'protein-1',
    source_protein_id: 'protein-1',
    locus_tag: 'locus-1',
    gene_id: '',
    old_locus_tag: '',
    gene: '',
    product: 'example protein',
    note: Array.from(longNote).slice(0, 50).join(''),
    qualifiers,
    selector: {
      qualifiers: structuredClone(qualifiers),
      hash: 'f_stable',
      location: '1..3',
      record_location: 'record:1..3:+'
    },
    location_parts: [{ start: 1, end: 3, strand: 1 }],
    nucleotide_sequence: 'ATG',
    amino_acid_sequence: 'M',
    sequence_warnings: []
  };
  const extracted = structuredClone(biological);
  delete extracted.feature_index;
  extracted.rendered_feature_svg_id = 'f_stable_record_1';
  return {
    extractedFeatures: [extracted],
    biologicalFeatures: [biological],
    selectedFeatureRecordIdx: 0
  };
};

test('current session feature catalog is compact and lossless', () => {
  const features = richFeatureState();
  const compact = compactSessionFeatureCatalog({
    format: 'gbdraw-session',
    version: 39,
    features
  });

  assert.equal(compact.features.extractedFeatures, undefined);
  assert.deepEqual(compact.features.featureCatalog, {
    schema: 1,
    encoding: 'biological-authority-v1',
    profile: 'rich-v1',
    extracted: [[0, 'f0', 'f_stable_record_1']]
  });
  assert.equal(compact.features.biologicalFeatures[0].stable_svg_id, undefined);
  assert.deepEqual(expandSessionFeatureCatalog(compact).features, features);
});

test('sanitized browser feature catalog retains its sequence-free shape', () => {
  const features = richFeatureState();
  for (const collection of [features.extractedFeatures, features.biologicalFeatures]) {
    for (const feature of collection) {
      delete feature.nucleotide_sequence;
      delete feature.amino_acid_sequence;
    }
  }

  const compact = compactSessionFeatureCatalog({ version: 39, features });

  assert.equal(compact.features.featureCatalog.profile, 'sanitized-v1');
  assert.deepEqual(expandSessionFeatureCatalog(compact).features, features);
});

test('compact feature catalog rejects an out-of-range biological reference', () => {
  const compact = compactSessionFeatureCatalog({
    version: 39,
    features: richFeatureState()
  });
  compact.features.featureCatalog.extracted[0][0] = 99;

  assert.throws(
    () => expandSessionFeatureCatalog(compact),
    /Invalid compact extracted-feature reference/
  );
});

test('compact feature catalog bounds expansion and requires string qualifiers', () => {
  const compact = compactSessionFeatureCatalog({
    version: 39,
    features: richFeatureState()
  });
  const duplicate = structuredClone(compact);
  duplicate.features.featureCatalog.extracted.push(
    structuredClone(duplicate.features.featureCatalog.extracted[0])
  );
  assert.throws(
    () => expandSessionFeatureCatalog(duplicate),
    /Invalid compact extracted-feature reference/
  );

  const malformedQualifier = structuredClone(compact);
  malformedQualifier.features.biologicalFeatures[0].qualifiers.protein_id = [null];
  assert.throws(
    () => expandSessionFeatureCatalog(malformedQualifier),
    /feature qualifiers must contain string arrays/
  );
});
