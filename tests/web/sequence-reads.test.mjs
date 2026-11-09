// The popups' sequence actions (D5): the feature popup opens without reading anything and reads its sequences when
// the Sequence tab is shown; the match popup's Copy and FASTA read the spans when pressed; a group's member FASTA
// reads the members when pressed. Each shows "Reading the sequence…" while Python reads, keeps a failure with its
// Retry, and writes nothing (popup, clipboard, file) for a read the door discarded.
import assert from 'node:assert/strict';
import { test } from 'node:test';
import { createSequenceReads } from '../../gbdraw/web/js/app/sequence-reads.js';
import { createCatalogSequences, featureSequenceKey } from '../../gbdraw/web/js/services/catalog-sequences.js';
import { buildFeaturePopup } from '../../gbdraw/web/js/services/feature-popup.js';
import { matchSpanUnavailableReason } from '../../gbdraw/web/js/services/pairwise-match-popup.js';
import { MATCH_SPAN_REASONS } from '../../gbdraw/web/js/contracts/schema.generated.js';

const feat = {
  id: 'rec-1\u0000kinase', recordKey: 'rec-1', biologicalFeatureId: 'kinase', record_id: 'REC', type: 'CDS',
  start: 0, end: 6, strand: 1, qualifiers: { product: ['kinase'] }
};
const popupOf = () => buildFeaturePopup(feat, {
  svgId: 'svg-kinase', legendName: 'CDS', color: '#000000', strokeColor: '#000000', strokeWidth: 0.5,
  featureVisibility: 'auto', identityEditable: true
});
const sequencesOf = (read) => ({ ...read, nucleotideSequence: 'ATGAAA', aminoAcidSequence: 'MK', sequenceWarnings: ['partial'] });

const harness = ({ readFeatures, readSpans, isShown } = {}) => {
  const copied = [];
  const downloaded = [];
  const reads = createSequenceReads({
    isShown: isShown || (() => true),
    readFeatures: readFeatures || (async (asked) => new Map(asked.map((read) => [featureSequenceKey(read), sequencesOf(read)]))),
    readSpans: readSpans || (async (spans) => spans.map((span) => ({
      status: 'ok', reason: null, recordId: 'REC', orientation: '+', length: 4, sequence: 'ACGT', sourceLength: 8
    }))),
    copyPending: async (pending) => {
      const text = await pending;
      copied.push(text);
      return text;
    },
    downloadText: (filename, text) => downloaded.push({ filename, text })
  });
  return { reads, copied, downloaded };
};

test('the feature popup opens unread; the Sequence tab reads once and shows the FASTA and warnings', async () => {
  let calls = 0;
  const { reads } = harness({ readFeatures: async (asked) => {
    calls += 1;
    return new Map(asked.map((read) => [featureSequenceKey(read), sequencesOf(read)]));
  } });
  const popup = popupOf();
  assert.equal(popup.sequence.status, 'idle');
  assert.equal(calls, 0, 'opening the popup reads nothing');
  const shown = reads.showFeatureSequence(popup);
  assert.equal(popup.sequence.status, 'reading');
  await shown;
  assert.equal(popup.sequence.status, 'ready');
  assert.match(popup.sequence.nucleotideFasta, /\nATGAAA$/);
  assert.match(popup.sequence.aminoAcidFasta, /\nMK$/);
  assert.deepEqual(popup.sequence.warnings, ['partial']);
  await reads.showFeatureSequence(popup);
  assert.equal(calls, 1, 'a shown tab is not read again');
});

test('a failed tab read shows the error, and Retry reads again', async () => {
  let fail = true;
  const { reads } = harness({ readFeatures: async (asked) => {
    if (fail) throw Object.assign(new Error('x'), { code: 'SEQUENCE_REQUIRES_GENERATE', stage: 'helper', context: {} });
    return new Map(asked.map((read) => [featureSequenceKey(read), sequencesOf(read)]));
  } });
  const popup = popupOf();
  await reads.showFeatureSequence(popup);
  assert.equal(popup.sequence.status, 'failed');
  assert.match(popup.sequence.error.summary, /read after the next Generate/);
  fail = false;
  await reads.showFeatureSequence(popup);
  assert.equal(popup.sequence.status, 'ready');
});

test('a tab read the door discarded is read again for the Result now shown', async () => {
  let calls = 0;
  const { reads } = harness({ readFeatures: async (asked) => {
    calls += 1;
    return calls === 1 ? null : new Map(asked.map((read) => [featureSequenceKey(read), sequencesOf(read)]));
  } });
  const popup = popupOf();
  await reads.showFeatureSequence(popup);
  assert.equal(calls, 2);
  assert.equal(popup.sequence.status, 'ready');
});

// D5 review r1 (LOW-7): the tab reads again for the Result now shown only while its popup is shown; a closed popup
// (or one showing another feature) stops the loop and keeps nothing it would show as reading.
test('a discarded tab read stops reading when its popup is no longer shown', async () => {
  let calls = 0;
  let shown = true;
  const popup = popupOf();
  const { reads } = harness({
    readFeatures: async () => {
      calls += 1;
      shown = false;
      return null;
    },
    isShown: (asked) => shown && asked === popup
  });
  await reads.showFeatureSequence(popup);
  assert.equal(calls, 1);
  assert.equal(popup.sequence.status, 'idle', 'the tab reads again when it is shown again');
});

const entryOf = (role, read) => ({
  span: { role, recordId: 'REC', start: 1, end: 4, displayRole: role === 'query' ? 'Query' : 'Subject' },
  orientation: '+', sequenceLength: 4, read, available: Boolean(read), unavailableReason: read ? '' : 'never supplied',
  matchId: 'm1', filename: `m1_${role}.fna`
});
const bundleOf = (entries) => ({ entries, combinedFilename: 'm1_both.fna', read: { status: '', error: null, retry: null } });

test('span Copy reads the span when pressed and copies its FASTA; FASTA downloads it', async () => {
  const asked = [];
  const { reads, copied, downloaded } = harness({ readSpans: async (spans) => {
    asked.push(spans);
    return spans.map(() => ({ status: 'ok', reason: null, recordId: 'REC', orientation: '+', length: 4, sequence: 'ACGT', sourceLength: 8 }));
  } });
  const bundle = bundleOf([entryOf('query', { kind: 'record', recordKey: 'rec-1', recordId: '', start: 1, end: 4 }), entryOf('subject', null)]);
  assert.equal(await reads.spanAction(bundle, [bundle.entries[0]], 'copy'), true);
  assert.deepEqual(asked, [[{ kind: 'record', recordKey: 'rec-1', recordId: '', start: 1, end: 4 }]]);
  assert.equal(copied[0], '>m1_query|record=REC|coords=1..4|strand=+\nACGT\n');
  assert.equal(bundle.read.status, '');
  const both = bundleOf([
    entryOf('query', { kind: 'record', recordKey: 'rec-1', recordId: '', start: 1, end: 4 }),
    entryOf('subject', { kind: 'record', recordKey: 'rec-2', recordId: '', start: 1, end: 4 })
  ]);
  assert.equal(reads.spansReadable(both.entries), true);
  await reads.spanAction(both, both.entries, 'download');
  assert.deepEqual(downloaded.map(({ filename }) => filename), ['m1_both.fna'], 'both spans download as one file');
  assert.equal(downloaded[0].text, '>m1_query|record=REC|coords=1..4|strand=+\nACGT\n>m1_subject|record=REC|coords=1..4|strand=+\nACGT\n');
});

// D5 review r1 (MEDIUM-3): Copy both spans and Download both spans write both spans or nothing, as before D5, when
// the popup opens with a span it cannot read (a BLAST source without its comparison FASTA) and when Python cannot
// read one when pressed: the unreadable span's row shows its reason.
test('both spans are offered only when every span can be read, and a span unreadable when pressed writes nothing', async () => {
  const replies = [];
  const { reads, copied, downloaded } = harness({ readSpans: async (spans) => replies.shift()(spans) });
  const partial = bundleOf([entryOf('query', { kind: 'record', recordKey: 'rec-1', recordId: '', start: 1, end: 4 }), entryOf('subject', null)]);
  assert.equal(reads.spansReadable(partial.entries), false);
  assert.equal(await reads.spanAction(partial, partial.entries, 'download'), false);
  assert.deepEqual(downloaded, [], 'no file holds the query span alone');
  assert.equal(partial.entries[1].unavailableReason, 'never supplied');

  const both = bundleOf([
    entryOf('query', { kind: 'record', recordKey: 'rec-1', recordId: '', start: 1, end: 4 }),
    entryOf('subject', { kind: 'record', recordKey: 'rec-2', recordId: '', start: 1, end: 9 })
  ]);
  replies.push((spans) => spans.map((span) => (span.recordKey === 'rec-2'
    ? { status: 'unavailable', reason: 'OUT_OF_RANGE', recordId: 'REC', orientation: '+', length: 0, sequence: '', sourceLength: 8 }
    : { status: 'ok', reason: null, recordId: 'REC', orientation: '+', length: 4, sequence: 'ACGT', sourceLength: 8 })));
  assert.equal(await reads.spanAction(both, both.entries, 'copy'), false);
  assert.deepEqual(copied, [], 'the clipboard gets no part of the pair');
  assert.deepEqual(both.entries.map((entry) => [entry.available, entry.unavailableReason]), [
    [true, ''], [false, 'Coordinates exceed the sequence length (8 bp).']
  ]);
  assert.equal(both.read.status, '', 'the reason shows in the span\'s row');
  assert.equal(reads.spansReadable(both.entries), false);
});

test('a span Python cannot read shows its reason in its row and copies nothing', async () => {
  const { reads, copied } = harness({ readSpans: async (spans) => spans.map(() => ({
    status: 'unavailable', reason: 'OUT_OF_RANGE', recordId: 'REC', orientation: '+', length: 0, sequence: '', sourceLength: 8
  })) });
  const bundle = bundleOf([entryOf('query', { kind: 'record', recordKey: 'rec-1', recordId: '', start: 1, end: 4 })]);
  assert.equal(await reads.spanAction(bundle, bundle.entries, 'copy'), false);
  assert.equal(bundle.entries[0].available, false);
  assert.equal(bundle.entries[0].unavailableReason, 'Coordinates exceed the sequence length (8 bp).');
  assert.equal(bundle.read.status, '', 'the reason shows in the span\'s row');
  assert.deepEqual(copied, []);
});

// Each reason Python gives shows the text the match popup showed for that case before D5, with the record and length.
test('each reason a span cannot be read shows the match popup\'s text from before D5', () => {
  const entry = entryOf('query', { kind: 'reference', recordId: 'NC_1', recordKey: null, start: 1, end: 4 });
  entry.span.recordId = 'NC_1';
  const textOf = (reason) => matchSpanUnavailableReason(entry, {
    status: 'unavailable', reason, recordId: '', orientation: '+', length: 0, sequence: '', sourceLength: 12345
  });
  assert.deepEqual(Object.fromEntries(MATCH_SPAN_REASONS.map((reason) => [reason, textOf(reason)])), {
    NEVER_SUPPLIED: 'Comparison sequence was not supplied for this BLAST source.',
    SOURCE_UNAVAILABLE: 'The sequence source for this match is unavailable.',
    SOURCE_KEY_CONFLICT: 'The sequence source key conflicts with the match identity.',
    RECORD_UNKNOWN: 'No sequence source matched record “NC_1”.',
    RECORD_AMBIGUOUS: 'Record ID “NC_1” is ambiguous in the available sequence sources.',
    RECORD_ALIAS_AMBIGUOUS: 'Record alias “NC_1” is ambiguous in the available sequence sources.',
    COORDINATE_NOT_WHOLE: 'Coordinates must be whole 1-based values.',
    COORDINATE_BELOW_ONE: 'Coordinates must be greater than or equal to 1.',
    OUT_OF_RANGE: `Coordinates exceed the sequence length (${(12345).toLocaleString()} bp).`
  });
});

test('a discarded span read writes nothing; a failed one keeps its error and Retry', async () => {
  let reply = null;
  const { reads, copied } = harness({ readSpans: async (spans) => {
    if (reply === 'fail') throw Object.assign(new Error('x'), { code: 'RESOURCE_INVALID', stage: 'helper' });
    return reply;
  } });
  const bundle = bundleOf([entryOf('query', { kind: 'record', recordKey: 'rec-1', recordId: '', start: 1, end: 4 })]);
  assert.equal(await reads.spanAction(bundle, bundle.entries, 'copy'), false);
  assert.equal(bundle.read.status, '');
  assert.deepEqual(copied, []);
  reply = 'fail';
  await reads.spanAction(bundle, bundle.entries, 'copy');
  assert.equal(bundle.read.status, 'failed');
  assert.equal(bundle.read.error.operation, 'readMatchSequences');
  reply = [{ status: 'ok', reason: null, recordId: 'REC', orientation: '+', length: 4, sequence: 'ACGT', sourceLength: 8 }];
  assert.equal(await bundle.read.retry(), true);
  assert.equal(bundle.read.status, '');
  assert.equal(copied.length, 1);
});

// D5 review r1 (MEDIUM-4): a refused clipboard write is told apart from a failed read: the text was read, so the
// popup says "Copy failed" with Retry, as the Similarity group buttons said before D5, not an unknown failure.
test('a clipboard write the browser refuses says Copy failed with Retry; a later action supersedes an earlier one', async () => {
  let refuse = true;
  const reads = createSequenceReads({
    isShown: () => true,
    readFeatures: async (asked) => new Map(asked.map((read) => [featureSequenceKey(read), sequencesOf(read)])),
    readSpans: async () => null,
    copyPending: async (pending) => {
      const text = await pending;
      if (refuse) throw new Error('denied');
      return text;
    },
    downloadText: () => {}
  });
  const state = { status: '', error: null, retry: null };
  const row = { read: { recordKey: 'rec-1', biologicalFeatureId: 'kinase', bases: 6 }, fastaFeature: feat };
  assert.equal(await reads.memberAction(state, [row], 'nt', 'copy', 'x.fna'), false);
  assert.equal(state.status, 'copy-failed');
  assert.equal(state.error, null);
  refuse = false;
  assert.equal(await state.retry(), true);
  assert.equal(state.status, '');
  let release;
  const slow = createSequenceReads({
    isShown: () => true,
    readFeatures: () => new Promise((_, reject) => { release = () => reject(new Error('worker down')); }),
    readSpans: async () => null, copyPending: async (pending) => pending, downloadText: () => {}
  });
  const first = slow.memberAction(state, [row], 'nt', 'copy', 'x.fna');
  const firstRelease = release;
  const second = slow.memberAction(state, [row], 'nt', 'copy', 'x.fna');
  firstRelease();
  assert.equal(await first, false);
  assert.equal(state.status, 'reading', 'the earlier action writes nothing once a later one runs');
  release();
  assert.equal(await second, false);
  assert.equal(state.status, 'failed');
});

test('member FASTA reads the members when pressed and copies the kind asked', async () => {
  const asked = [];
  const { reads, copied, downloaded } = harness({ readFeatures: async (reads_) => {
    asked.push(reads_);
    return new Map(reads_.map((read) => [featureSequenceKey(read), sequencesOf(read)]));
  } });
  const read = { recordKey: 'rec-1', biologicalFeatureId: 'kinase', bases: 6 };
  const rows = [{ read, fastaFeature: feat }, { read: null, fastaFeature: { id: 'orphan' } }];
  const state = { status: '', error: null, retry: null };
  await reads.memberAction(state, rows, 'aa', 'copy', 'g_aa.faa');
  assert.deepEqual(asked, [[read]]);
  assert.match(copied[0], /\nMK\n$/);
  await reads.memberAction(state, rows, 'nt', 'download', 'g_nt.fna');
  assert.equal(downloaded[0].filename, 'g_nt.fna');
  assert.match(downloaded[0].text, /\nATGAAA\n$/);
});

test('a row without an identity in an older Session shows the reason to Generate, not "No sequence"', async () => {
  // A Session 27-39 that kept only its Results has no request to read from; its rows may lack the identity a read
  // names. The door, not the row, says why nothing can be read.
  const door = createCatalogSequences({ getCommittedSession: () => null, helpers: /** @type {never} */ ({}) });
  const { reads, copied } = harness({ readFeatures: door.readFeatures, readSpans: door.readSpans });
  const popup = buildFeaturePopup({ ...feat, biologicalFeatureId: undefined }, {
    svgId: 'svg-kinase', legendName: 'CDS', color: '#000000', strokeColor: '#000000', strokeWidth: 0.5,
    featureVisibility: 'auto', identityEditable: true
  });
  await reads.showFeatureSequence(popup);
  assert.equal(popup.sequence.status, 'failed');
  assert.equal(popup.sequence.error.code, 'SEQUENCE_REQUIRES_GENERATE');
  // A Session 27-39 that kept its request draws Results from SVG metadata without catalog identities: the
  // Sequence tab gives the same reason, not "No sequence available".
  let asked = 0;
  const withRequest = harness({ readFeatures: async (features) => { asked += features.length; return new Map(); } });
  const legacy = buildFeaturePopup({ ...feat, recordKey: undefined, biologicalFeatureId: undefined }, {
    svgId: 'svg-kinase', legendName: 'CDS', color: '#000000', strokeColor: '#000000', strokeWidth: 0.5,
    featureVisibility: 'auto', identityEditable: true
  });
  await withRequest.reads.showFeatureSequence(legacy);
  assert.equal(legacy.sequence.status, 'failed');
  assert.equal(legacy.sequence.error.code, 'SEQUENCE_REQUIRES_GENERATE');
  assert.equal(asked, 0, 'Python is not asked for a feature it cannot find');
  const state = { status: '', error: null, retry: null };
  await reads.memberAction(state, [{ read: null, fastaFeature: feat }], 'nt', 'copy', 'members.fasta');
  assert.equal(state.status, 'failed');
  assert.equal(state.error.code, 'SEQUENCE_REQUIRES_GENERATE');
  assert.deepEqual(copied, []);
});
