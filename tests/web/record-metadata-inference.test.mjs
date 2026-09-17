import test from 'node:test';
import assert from 'node:assert/strict';

import {
  formatInferredOrganismStrain,
  formatInferredSubtitle,
  extractGenBankMetadata,
  parseSequenceRecordText
} from '../../gbdraw/web/js/app/record-discovery.js';

test('formatInferredOrganismStrain formats binomial names with italics and appends strain', () => {
  assert.equal(
    formatInferredOrganismStrain({
      organism: 'Escherichia coli O157:H7',
      strain: 'Sakai'
    }),
    '<i>Escherichia coli</i> O157:H7 Sakai'
  );

  assert.equal(
    formatInferredOrganismStrain({
      organism: 'Streptomyces lividus',
      strain: 'CBS 844.73'
    }),
    '<i>Streptomyces lividus</i> CBS 844.73'
  );

  assert.equal(
    formatInferredOrganismStrain({
      organism: 'Homo sapiens'
    }),
    '<i>Homo sapiens</i>'
  );

  assert.equal(
    formatInferredOrganismStrain({
      organism: 'Candidatus Tyloplasma litorale',
      strain: 'Fukuoka2020'
    }),
    'Candidatus <i>Tyloplasma litorale</i> Fukuoka2020'
  );

  assert.equal(
    formatInferredOrganismStrain({
      organism: 'Escherichia coli str. K-12 substr. MG1655',
      strain: 'K-12'
    }),
    '<i>Escherichia coli</i> str. K-12 substr. MG1655'
  );
});

test('formatInferredSubtitle detects complete genome, plasmids, chromosomes, and cluster titles', () => {
  assert.equal(
    formatInferredSubtitle({
      definition: 'Escherichia coli O157:H7 str. Sakai DNA, complete genome.'
    }),
    'Complete genome'
  );

  assert.equal(
    formatInferredSubtitle({
      definition: 'Homo sapiens mitochondrion, complete genome.'
    }),
    'Mitochondrion, complete genome'
  );

  assert.equal(
    formatInferredSubtitle({
      plasmid: 'pOSAK1'
    }),
    'Plasmid pOSAK1'
  );

  assert.equal(
    formatInferredSubtitle({
      chromosome: '1'
    }),
    'Chromosome 1'
  );

  assert.equal(
    formatInferredSubtitle({
      definition: 'Streptomyces lividus lividomycin biosynthesis gene cluster.',
      organism: 'Streptomyces lividus'
    }),
    'Lividomycin biosynthesis gene cluster'
  );
});

test('extractGenBankMetadata parses organism, strain, plasmid, and definition from text chunk', () => {
  const sampleChunk = `LOCUS       NC_002695            5498578 bp    DNA     circular CON 12-FEB-2021
DEFINITION  Escherichia coli O157:H7 str. Sakai DNA, complete genome.
ACCESSION   NC_002695
VERSION     NC_002695.2
SOURCE      Escherichia coli O157:H7 str. Sakai
  ORGANISM  Escherichia coli O157:H7 str. Sakai
            Bacteria; Pseudomonadati; Pseudomonadota.
FEATURES             Location/Qualifiers
     source          1..5498578
                     /organism="Escherichia coli O157:H7 str. Sakai"
                     /mol_type="genomic DNA"
                     /strain="Sakai"
                     /sub_strain="RIMD 0509952"
                     /serovar="O157:H7"
`;

  const meta = extractGenBankMetadata(sampleChunk);
  assert.equal(meta.organism, 'Escherichia coli O157:H7 str. Sakai');
  assert.equal(meta.strain, 'Sakai');
  assert.equal(meta.inferredDefinition, '<i>Escherichia coli</i> O157:H7 str. Sakai');
  assert.equal(meta.inferredSubtitle, 'Complete genome');
});

test('parseSequenceRecordText attaches inferred metadata to all records in multi-record file', () => {
  const multiRecordGenBank = `LOCUS       chr1                 1000 bp    DNA     circular BCT 01-JAN-2020
DEFINITION  Escherichia coli K-12 chromosome, complete genome.
ACCESSION   NC_000001
VERSION     NC_000001.1
FEATURES             Location/Qualifiers
     source          1..1000
                     /organism="Escherichia coli"
                     /strain="K-12"
//
LOCUS       plasmid1              500 bp    DNA     circular BCT 01-JAN-2020
DEFINITION  Escherichia coli K-12 plasmid pTEST, complete sequence.
ACCESSION   NC_000002
VERSION     NC_000002.1
FEATURES             Location/Qualifiers
     source          1..500
                     /organism="Escherichia coli"
                     /strain="K-12"
                     /plasmid="pTEST"
//
`;

  const records = parseSequenceRecordText(multiRecordGenBank, 'genbank');
  assert.equal(records.length, 2);
  assert.equal(records[0].inferredDefinition, '<i>Escherichia coli</i> K-12');
  assert.equal(records[0].inferredSubtitle, 'Complete genome');
  assert.equal(records[1].inferredDefinition, '<i>Escherichia coli</i> K-12');
  assert.equal(records[1].inferredSubtitle, 'Plasmid pTEST');
});

