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
  inferredDefinitionForRecord,
  linearFileCardRecord,
  resolveLinearRecordEffectiveDefinition,
  resolveLinearRecordEffectiveSubtitle
} = await import('../../gbdraw/web/js/services/linear-sources.js');
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

// IN-06 (D-12, PD-OI-067): a record draws its own Definition, else the File
// default the user entered, else the definition inferred from that record.
test('each record falls back to its own inferred definition after the File default', () => {
  const file = { name: 'two_organisms.gb' };
  const recordA = createLinearSeq({ uid: 'a', gb: file, inferred_definition: '<i>Escherichia coli</i> K-12' });
  const recordB = createLinearSeq({ uid: 'b', gb: file, inferred_definition: '<i>Bacillus subtilis</i> 168' });
  assert.equal(resolveLinearRecordEffectiveDefinition(recordA), '<i>Escherichia coli</i> K-12');
  assert.equal(resolveLinearRecordEffectiveDefinition(recordB), '<i>Bacillus subtilis</i> 168');
  const [group] = groupLinearSourceRecords([recordA, recordB]);
  setLinearSourceDefaultDefinition(group, 'Shared organism');
  assert.equal(resolveLinearRecordEffectiveDefinition(recordA), 'Shared organism');
  assert.equal(resolveLinearRecordEffectiveDefinition(recordB), 'Shared organism');
  recordB.definition = 'Own value';
  assert.equal(resolveLinearRecordEffectiveDefinition(recordB), 'Own value');
  setLinearSourceDefaultDefinition(group, '');
  recordB.definition = '';
  assert.equal(resolveLinearRecordEffectiveDefinition(recordB), '<i>Bacillus subtilis</i> 168');
});

test('the inferred definition belongs to the record a row selects', () => {
  const records = [
    { selector: '#1', recordId: 'REC_A.1', inferredDefinition: 'A' },
    { selector: '#2', recordId: 'REC_B.1', inferredDefinition: 'B' }
  ];
  assert.equal(inferredDefinitionForRecord(records, 'REC_B.1'), 'B');
  assert.equal(inferredDefinitionForRecord(records, '#1'), 'A');
  assert.equal(inferredDefinitionForRecord(records, ''), '', 'an automatic selector over two records names none');
  assert.equal(inferredDefinitionForRecord(records.slice(1), ''), 'B');
  assert.equal(inferredDefinitionForRecord(records, 'missing'), '');
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

test('sessions without per-record values fall back to comparing label and fileDefinition', async () => {
  const { projectCanonicalSessionRequest, CANONICAL_REQUEST_SCHEMA } = await import('../../gbdraw/web/js/services/session-request.js');
  const canonicalRequest = {
    renderRequest: {
      schema: CANONICAL_REQUEST_SCHEMA,
      mode: 'linear',
      grouping: 'single',
      records: [
        {
          recordKey: 'rec-1',
          source: { kind: 'genbank', resourceId: 'res-1' },
          presentation: { label: '<i>Escherichia coli</i>', subtitle: 'Chromosome 1' },
          display: { isCircular: null, startCoordinate: null }
        },
        {
          recordKey: 'rec-2',
          source: { kind: 'genbank', resourceId: 'res-1' },
          presentation: { label: '<i>Custom Override</i>', subtitle: 'Overridden Subtitle' },
          display: { isCircular: null, startCoordinate: null }
        }
      ],
      diagramOptions: { featurePlacements: [], featureOverrides: [] },
      output: { prefix: 'test' }
    },
    webFiles: {
      linearRecordMetadata: [
        { recordKey: 'rec-1', fileDefinition: '<i>Escherichia coli</i>', fileSubtitle: 'Chromosome 1' },
        { recordKey: 'rec-2', fileDefinition: '<i>Escherichia coli</i>', fileSubtitle: 'Chromosome 1' }
      ]
    },
    resources: {
      'res-1': { kind: 'file', name: 'records.gb' }
    }
  };

  const projected = projectCanonicalSessionRequest(canonicalRequest);
  assert.equal(projected.files.linearSeqs.length, 2);

  // rec-1 inherited file default: definition and record_subtitle must be empty string
  assert.equal(projected.files.linearSeqs[0].file_definition, '<i>Escherichia coli</i>');
  assert.equal(projected.files.linearSeqs[0].file_subtitle, 'Chromosome 1');
  assert.equal(projected.files.linearSeqs[0].definition, '');
  assert.equal(projected.files.linearSeqs[0].record_subtitle, '');

  // rec-2 had an explicit override: definition and record_subtitle must be preserved
  assert.equal(projected.files.linearSeqs[1].file_definition, '<i>Escherichia coli</i>');
  assert.equal(projected.files.linearSeqs[1].file_subtitle, 'Chromosome 1');
  assert.equal(projected.files.linearSeqs[1].definition, '<i>Custom Override</i>');
  assert.equal(projected.files.linearSeqs[1].record_subtitle, 'Overridden Subtitle');
});


test('an explicit override equal to the file default survives a save and reload', async () => {
  const { projectCanonicalSessionRequest, CANONICAL_REQUEST_SCHEMA } = await import('../../gbdraw/web/js/services/session-request.js');
  const shared = '<i>Escherichia coli</i>';
  const canonicalRequest = {
    renderRequest: {
      schema: CANONICAL_REQUEST_SCHEMA,
      mode: 'linear',
      grouping: 'single',
      records: [
        {
          recordKey: 'rec-1',
          source: { kind: 'genbank', resourceId: 'res-1' },
          presentation: { label: shared, subtitle: 'Chromosome 1' },
          display: { isCircular: null, startCoordinate: null }
        },
        {
          recordKey: 'rec-2',
          source: { kind: 'genbank', resourceId: 'res-1' },
          presentation: { label: shared, subtitle: 'Chromosome 1' },
          display: { isCircular: null, startCoordinate: null }
        }
      ],
      diagramOptions: { featurePlacements: [], featureOverrides: [] },
      output: { prefix: 'test' }
    },
    webFiles: {
      linearRecordMetadata: [
        // Inherited: the record carries no value of its own.
        {
          recordKey: 'rec-1',
          fileDefinition: shared,
          fileSubtitle: 'Chromosome 1',
          recordDefinition: '',
          recordSubtitle: ''
        },
        // Pinned: the user typed the same text the file default happens to have.
        {
          recordKey: 'rec-2',
          fileDefinition: shared,
          fileSubtitle: 'Chromosome 1',
          recordDefinition: shared,
          recordSubtitle: 'Chromosome 1'
        }
      ]
    },
    resources: { 'res-1': { kind: 'file', name: 'records.gb' } }
  };

  const projected = projectCanonicalSessionRequest(canonicalRequest);
  assert.equal(projected.files.linearSeqs[0].definition, '');
  assert.equal(projected.files.linearSeqs[0].record_subtitle, '');
  // Both records render the same text, but only rec-2 keeps following the file
  // default when the user later edits it.
  assert.equal(projected.files.linearSeqs[1].definition, shared);
  assert.equal(projected.files.linearSeqs[1].record_subtitle, 'Chromosome 1');
});

// Which name fields a File card shows is one decision of the File's records
// and stored File values (web.linear.file-default-definition).
test('a one-record File shows its record names on the card and no empty File names block', () => {
  const record = createLinearSeq({ uid: 'only', gb: { name: 'one.gb' }, inferred_definition: 'Inferred' });
  const [source] = groupLinearSourceRecords([record]);
  const cardRecord = linearFileCardRecord(source);
  assert.equal(cardRecord, source.records[0], 'the record entry itself, found by identity');
  assert.equal(cardRecord.sequence, record);
});

test('a one-record File that holds a File-level value keeps the File names block alone', () => {
  for (const stored of [{ file_definition: 'File organism' }, { file_subtitle: 'File subtitle' }]) {
    const record = createLinearSeq({ uid: 'only', gb: { name: 'one.gb' }, definition: 'Own', ...stored });
    const [source] = groupLinearSourceRecords([record]);
    // One names block: the record's names stay in Record options.
    assert.equal(linearFileCardRecord(source), null, JSON.stringify(stored));
    // Clearing the stored value puts the record's names on the card instead.
    setLinearSourceDefaultDefinition(source, '');
    setLinearSourceDefaultSubtitle(source, '');
    assert.equal(linearFileCardRecord(source), source.records[0]);
  }
});

test('a File with two or more records keeps the File names block and per-record names', () => {
  const file = { name: 'two.gb' };
  const records = [createLinearSeq({ uid: 'a', gb: file }), createLinearSeq({ uid: 'b', gb: file })];
  const [source] = groupLinearSourceRecords(records);
  assert.equal(linearFileCardRecord(source), null);
  setLinearSourceDefaultDefinition(source, 'Shared');
  assert.equal(linearFileCardRecord(source), null);
});

test('a File card without an inspected or uploaded file follows the rows it has', () => {
  // No file and no discovery state is read: one row is the one-record card.
  const [empty] = groupLinearSourceRecords([createLinearSeq({ uid: 'empty' })]);
  assert.equal(linearFileCardRecord(empty), empty.records[0]);
});

test('the File card template takes its name fields from the one decision and one component', async () => {
  const { readFileSync } = await import('node:fs');
  const html = readFileSync(new URL('../../gbdraw/web/index.html', import.meta.url), 'utf8');
  const lines = html.split('\n');
  const block = lines.find((line) => line.includes(':data-linear-source-defaults="source.uid"'));
  assert.match(block, /'linear-file-names-unused': linearFileCardRecord\(source\) \}/);
  const uses = lines.filter((line) => line.includes('<linear-record-names'));
  assert.equal(uses.length, 2, 'on the File card and in Record options');
  assert.match(uses[0], /v-if="linearFileCardRecord\(source\)"/);
  assert.match(uses[1], /v-if="!linearFileCardRecord\(source\)"/);
  // The record's fields exist once: in the component both places render.
  for (const label of ["t('Definition for sequence {n}', { n: idx + 1 })", "t('Subtitle / title for sequence {n}', { n: idx + 1 })"]) {
    assert.equal(lines.filter((line) => line.includes(label)).length, 1, label);
  }
  const template = html.slice(html.indexOf('id="linear-record-names-template"'));
  assert.ok(template.indexOf("t('Definition for sequence {n}', { n: idx + 1 })") > 0, 'the fields live in the component template');
  assert.match(template.slice(0, template.indexOf('</script>')), /v-model="seq\.definition"[\s\S]*v-model="seq\.record_subtitle"/);
});

test('a blank record Definition placeholder previews the inherited definition Generate draws', async () => {
  const { readFileSync } = await import('node:fs');
  const html = readFileSync(new URL('../../gbdraw/web/index.html', import.meta.url), 'utf8');
  const input = html.split('\n').find((line) => line.includes("t('Definition for sequence {n}', { n: idx + 1 })"));
  assert.ok(input, 'record Definition input not found');
  // One rule: the placeholder uses the function that resolves the drawn text,
  // so File default, then the inferred definition, then the example text.
  assert.match(
    input,
    /:placeholder="resolveLinearRecordEffectiveDefinition\(seq, source\) \|\| t\('e\.g\., \{example\}', \{ example: 'Escherichia coli O157:H7 Sakai' \}\)"/
  );

  const inferredOnly = createLinearSeq({ inferred_definition: '<i>Vibrio</i> sp.' });
  const source = { records: [{ sequence: inferredOnly }], sequence: inferredOnly };
  assert.equal(resolveLinearRecordEffectiveDefinition(inferredOnly, source), '<i>Vibrio</i> sp.');
  inferredOnly.file_definition = 'File default';
  assert.equal(resolveLinearRecordEffectiveDefinition(inferredOnly, source), 'File default');
  assert.equal(resolveLinearRecordEffectiveDefinition(createLinearSeq(), { records: [] }), '');
});
