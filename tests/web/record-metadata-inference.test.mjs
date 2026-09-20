import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';

import {
  formatInferredOrganismStrain,
  extractGenBankMetadata,
  parseSequenceRecordText
} from '../../gbdraw/web/js/app/record-discovery.js';

// gbdraw/core/record_metadata.py owns this inference for the Worker and the
// renderer; record-discovery.js reimplements it for the no-Worker upload path.
// tests/test_record_metadata.py asserts the same table, so the two cannot drift.
const CASES = JSON.parse(
  readFileSync(new URL('../fixtures/record_metadata_inference_cases.json', import.meta.url), 'utf8')
);

test('formatInferredOrganismStrain matches the shared inference table', () => {
  for (const testCase of CASES.definition) {
    assert.equal(
      formatInferredOrganismStrain({ organism: testCase.organism, strain: testCase.strain }),
      testCase.expected,
      testCase.name
    );
  }
});

test('extractGenBankMetadata parses organism and strain without automatic subtitles', () => {
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
  assert.equal(Object.hasOwn(meta, 'inferredSubtitle'), false);
  assert.equal(Object.hasOwn(meta, 'inferredSubtitleFromReplicon'), false);
});

test('extractGenBankMetadata reads /isolate ahead of /strain, like infer_record_source_metadata', () => {
  const chunk = `LOCUS       TEST                 1000 bp    DNA     circular BCT 01-JAN-2020
DEFINITION  Bacillus subtilis genomic DNA.
FEATURES             Location/Qualifiers
     source          1..1000
                     /organism="Bacillus subtilis"
                     /strain="ignored"
                     /isolate="168"
`;
  assert.equal(extractGenBankMetadata(chunk).strain, '168');
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

  assert.equal(records[1].inferredDefinition, '<i>Escherichia coli</i> K-12');
  for (const record of records) {
    assert.equal(Object.hasOwn(record, 'inferredSubtitle'), false);
    assert.equal(Object.hasOwn(record, 'inferredSubtitleFromReplicon'), false);
  }
});
