// The Web record count follows gbdraw.io.input_text.flat_file_format: the first
// non-blank line names the format, and an EMBL record starts at "ID   ".
import assert from 'node:assert/strict';
import test from 'node:test';

import { countGenBankRecords } from '../../gbdraw/web/js/services/genbank-header.js';

test('record counts read EMBL by its first non-blank line', () => {
  const record = 'LOCUS       GOOD1                   4 bp    DNA     linear   UNK 01-JAN-1980\nORIGIN\n        1 acgt\n//\n';
  const embl = 'ID   EMBL1; SV 1; linear; DNA; STD; PRO; 4 BP.\nXX\nSQ   Sequence 4 BP;\n     acgt 4\n//\n';
  assert.equal(countGenBankRecords(`\n${embl}${embl.replace('EMBL1', 'EMBL2')}`), 2);
  assert.equal(countGenBankRecords(`${record}ID   NOT_A_RECORD\n`), 1);
  assert.equal(countGenBankRecords(` ${embl}`), 0);
});
