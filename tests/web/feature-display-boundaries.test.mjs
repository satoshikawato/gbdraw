import assert from 'node:assert/strict';
import test from 'node:test';
import { mkdtemp, readFile, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';

const tempDir = await mkdtemp(join(tmpdir(), 'gbdraw-feature-display-boundaries-'));
await writeFile(join(tempDir, 'package.json'), '{"type":"module"}\n', 'utf8');
await writeFile(
  join(tempDir, 'feature-utils.js'),
  await readFile('gbdraw/web/js/app/feature-utils.js', 'utf8'),
  'utf8'
);
await writeFile(
  join(tempDir, 'feature-sequence-fasta.js'),
  await readFile('gbdraw/web/js/app/feature-sequence-fasta.js', 'utf8'),
  'utf8'
);

const { buildFeatureSequenceFastas } = await import(
  pathToFileURL(join(tempDir, 'feature-sequence-fasta.js'))
);
const { getFeatureCaption } = await import(
  pathToFileURL(join(tempDir, 'feature-utils.js'))
);

const runtimeHandle = `h_${'a'.repeat(26)}`;
const featureAnalysisId = `f_${'b'.repeat(64)}`;
const unsupportedHistoricalTransportId = `record@instance|alias~f_${'c'.repeat(64)}`;
const legacyId = 'p_r_old_0_9_1_deadbeefdead';
const generatedId = 'gbd_r0001_cds000001';

test('protein FASTA skips internal IDs and uses biological aliases', () => {
  const { aminoAcidFasta } = buildFeatureSequenceFastas({
    record_id: 'record-a',
    start: 0,
    end: 9,
    strand: '+',
    source_protein_id: '',
    protein_id: runtimeHandle,
    locus_tag: 'LOCUS_SAFE',
    amino_acid_sequence: 'MPEPTIDE'
  });
  assert.match(aminoAcidFasta, /^>LOCUS_SAFE\b/);
  assert.doesNotMatch(aminoAcidFasta, /h_[a-z2-7]{26}/);

  for (const proteinId of [
    runtimeHandle,
    featureAnalysisId,
    unsupportedHistoricalTransportId,
    legacyId,
    generatedId
  ]) {
    const fasta = buildFeatureSequenceFastas({
      protein_id: proteinId,
      locus_tag: 'BIOLOGICAL_ALIAS',
      amino_acid_sequence: 'MK'
    }).aminoAcidFasta;
    assert.match(fasta, /^>BIOLOGICAL_ALIAS\b/);
    assert.doesNotMatch(fasta, new RegExp(proteinId.replace(/[.*+?^${}()|[\]\\]/g, '\\$&')));
  }
});

test('feature captions ignore internal labels and use safe fallbacks', () => {
  for (const label of [
    runtimeHandle,
    featureAnalysisId,
    unsupportedHistoricalTransportId,
    legacyId,
    generatedId
  ]) {
    assert.equal(getFeatureCaption({
      label,
      product: 'safe product',
      type: 'CDS',
      start: 0,
      end: 9
    }), 'safe product');
  }

  assert.equal(getFeatureCaption({
    label: generatedId,
    type: 'CDS',
    start: 0,
    end: 9
  }), 'CDS at 0..9');
  assert.equal(getFeatureCaption({
    label: featureAnalysisId,
    displayLabel: runtimeHandle,
    locus_tag: 'LOCUS_CAPTION_SAFE'
  }), 'LOCUS_CAPTION_SAFE');
  assert.equal(getFeatureCaption({
    product: runtimeHandle,
    gene: unsupportedHistoricalTransportId,
    locus_tag: 'LOCUS_AFTER_INTERNAL_PRODUCT'
  }), 'LOCUS_AFTER_INTERNAL_PRODUCT');
});

test('FASTA descriptions skip internal metadata and retain safe later candidates', () => {
  const safeFallback = buildFeatureSequenceFastas({
    source_protein_id: 'WP_SAFE_DESCRIPTION.1',
    product: runtimeHandle,
    protein: featureAnalysisId,
    gene: unsupportedHistoricalTransportId,
    locus_tag: 'LOCUS_DESCRIPTION_SAFE',
    amino_acid_sequence: 'MK'
  }).aminoAcidFasta;
  assert.match(
    safeFallback,
    /^>WP_SAFE_DESCRIPTION\.1 LOCUS_DESCRIPTION_SAFE\n/
  );
  assert.doesNotMatch(safeFallback, /(?:h_[a-z2-7]{26}|f_[0-9a-f]{64}|~f_)/);

  const internalOnly = buildFeatureSequenceFastas({
    source_protein_id: 'WP_SAFE_HEADER.1',
    product: runtimeHandle,
    protein: featureAnalysisId,
    gene: unsupportedHistoricalTransportId,
    locus_tag: legacyId,
    name: generatedId,
    amino_acid_sequence: 'MK'
  }).aminoAcidFasta;
  assert.equal(internalOnly, '>WP_SAFE_HEADER.1\nMK');
});
