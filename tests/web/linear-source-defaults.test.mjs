import assert from 'node:assert/strict';
import test from 'node:test';

globalThis.window = {
  Vue: {
    ref: (value) => ({ value }),
    reactive: (obj) => obj,
    computed: (getter) => ({
      get value() {
        return getter();
      }
    }),
    watch: () => {},
    onMounted: () => {}
  }
};

const {
  groupLinearSourceRecords,
  getLinearSourceDefaultDefinition,
  setLinearSourceDefaultDefinition,
  getLinearSourceDefaultSubtitle,
  setLinearSourceDefaultSubtitle,
  resolveLinearRecordEffectiveDefinition,
  resolveLinearRecordEffectiveSubtitle
} = await import('../../gbdraw/web/js/app/linear-sources.js');
const { createLinearSeq } = await import('../../gbdraw/web/js/state.js');

test('createLinearSeq handles file_definition and file_subtitle', () => {
  const seq = createLinearSeq({
    file_definition: '<i>Streptomyces coelicolor</i>',
    file_subtitle: 'A3(2)'
  });
  assert.equal(seq.file_definition, '<i>Streptomyces coelicolor</i>');
  assert.equal(seq.file_subtitle, 'A3(2)');
  assert.equal(seq.definition, '');
  assert.equal(seq.record_subtitle, '');

  const defSeq = createLinearSeq();
  assert.equal(defSeq.file_definition, '');
  assert.equal(defSeq.file_subtitle, '');
});

test('linear sources default getters and setters operate on grouped records', () => {
  const sharedFile = { name: 'multi_record.gb' };
  const seq1 = createLinearSeq({ uid: 'seq-1', gb: sharedFile });
  const seq2 = createLinearSeq({ uid: 'seq-2', gb: sharedFile });
  const seq3 = createLinearSeq({ uid: 'seq-3', gb: { name: 'other.gb' } });

  const groups = groupLinearSourceRecords([seq1, seq2, seq3]);
  assert.equal(groups.length, 2);
  const multiGroup = groups[0];
  const singleGroup = groups[1];

  assert.equal(getLinearSourceDefaultDefinition(multiGroup), '');
  assert.equal(getLinearSourceDefaultSubtitle(multiGroup), '');

  setLinearSourceDefaultDefinition(multiGroup, '<i>Escherichia coli</i> O157:H7');
  setLinearSourceDefaultSubtitle(multiGroup, 'Complete genome');

  assert.equal(getLinearSourceDefaultDefinition(multiGroup), '<i>Escherichia coli</i> O157:H7');
  assert.equal(getLinearSourceDefaultSubtitle(multiGroup), 'Complete genome');
  assert.equal(seq1.file_definition, '<i>Escherichia coli</i> O157:H7');
  assert.equal(seq2.file_definition, '<i>Escherichia coli</i> O157:H7');
  assert.equal(seq1.file_subtitle, 'Complete genome');
  assert.equal(seq2.file_subtitle, 'Complete genome');

  // other group remains unaffected
  assert.equal(getLinearSourceDefaultDefinition(singleGroup), '');
  assert.equal(seq3.file_definition, '');
});

test('resolveLinearRecordEffectiveDefinition falls back to file default when record definition is empty', () => {
  const sharedFile = { name: 'multi_record.gb' };
  const seq1 = createLinearSeq({
    uid: 'seq-1',
    gb: sharedFile,
    file_definition: '<i>Escherichia coli</i>'
  });
  const seq2 = createLinearSeq({
    uid: 'seq-2',
    gb: sharedFile,
    definition: '<i>Custom E. coli</i>',
    file_definition: '<i>Escherichia coli</i>'
  });
  const groups = groupLinearSourceRecords([seq1, seq2]);

  assert.equal(resolveLinearRecordEffectiveDefinition(seq1, groups[0]), '<i>Escherichia coli</i>');
  assert.equal(resolveLinearRecordEffectiveDefinition(seq2, groups[0]), '<i>Custom E. coli</i>');

  // Test without passing source group explicitly
  assert.equal(resolveLinearRecordEffectiveDefinition(seq1), '<i>Escherichia coli</i>');
  assert.equal(resolveLinearRecordEffectiveDefinition(seq2), '<i>Custom E. coli</i>');
});

test('resolveLinearRecordEffectiveSubtitle falls back to file default when record subtitle is empty', () => {
  const sharedFile = { name: 'multi_record.gb' };
  const seq1 = createLinearSeq({
    uid: 'seq-1',
    gb: sharedFile,
    file_subtitle: 'Default subtitle'
  });
  const seq2 = createLinearSeq({
    uid: 'seq-2',
    gb: sharedFile,
    record_subtitle: 'Overridden subtitle',
    file_subtitle: 'Default subtitle'
  });
  const groups = groupLinearSourceRecords([seq1, seq2]);

  assert.equal(resolveLinearRecordEffectiveSubtitle(seq1, groups[0]), 'Default subtitle');
  assert.equal(resolveLinearRecordEffectiveSubtitle(seq2, groups[0]), 'Overridden subtitle');

  assert.equal(resolveLinearRecordEffectiveSubtitle(seq1), 'Default subtitle');
  assert.equal(resolveLinearRecordEffectiveSubtitle(seq2), 'Overridden subtitle');
});
