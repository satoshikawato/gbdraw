// Copy and Download of a Similarity group member's or a whole group's FASTA read the sequences from Python through
// the catalog-sequence door (D5): the members carry none, a member is available when its catalog row has an
// identity, one read asks every member once, and the read runs through the one owner of a popup's member read state
// (app/sequence-reads.js), which the Similarity group popup of a match uses too: "Reading the sequence…" while
// Python reads, "No sequence available", or the failure with Details, Copy diagnostics and Retry in place. A stale
// reply writes nothing, and a read the popup no longer shows writes nothing to the state it shows.
import assert from 'node:assert/strict';
import test from 'node:test';
import { withDrawings } from './helpers/drawing-state.mjs';
import { loadPageVue } from './helpers/vue-runtime.mjs';

const Vue = loadPageVue();
globalThis.window = { Vue, setTimeout, clearTimeout };
globalThis.isSecureContext = true;
let copied = null;
const downloads = [];

const { createOrthogroupEditor } = await import('../../gbdraw/web/js/app/orthogroups.js');
const { createSequenceReads } = await import('../../gbdraw/web/js/app/sequence-reads.js');
const { featureSequenceKey } = await import('../../gbdraw/web/js/services/catalog-sequences.js');
const { diagnosticError } = await import('../../gbdraw/web/js/utils/error-normalization.js');
const { ref, reactive } = Vue;

const row = (id, start, end, extra = {}) => ({
  recordKey: 'rec-1', biologicalFeatureId: id, record_id: 'REC', type: 'CDS', start, end, strand: 1,
  locus_tag: id.toUpperCase(), ...extra
});
// A schema-5 row that still carries its sequences: the members never take them.
const alpha = row('alpha', 0, 12, { nucleotide_sequence: 'CCCCCC', amino_acid_sequence: 'WWW' });
const beta = row('beta', 30, 39);
const member = (feature) => ({ recordKey: feature.recordKey, biologicalFeatureId: feature.biologicalFeatureId });
const orphan = { recordKey: 'rec-1', biologicalFeatureId: 'not-in-catalog' };
const group = { id: 'og_1', name: 'LivB', members: [member(alpha), member(beta), orphan] };

const sequencesOf = (feature, nucleotideSequence, aminoAcidSequence) => [
  featureSequenceKey(member(feature)),
  { ...member(feature), nucleotideSequence, aminoAcidSequence, sequenceWarnings: [] }
];

/** An editor whose door answers with `answer(features)`; each call is recorded. */
const setup = (answer) => {
  const reads = [];
  const state = {
    mode: 'linear', orthogroups: ref([group]),
    orthogroupNameOverrides: reactive({}), orthogroupDescriptionOverrides: reactive({}),
    selectedOrthogroupId: ref('og_1'), orthogroupSearch: ref(''), orthogroupSortMode: ref('id'),
    clickedFeature: ref(null), showRightDrawer: ref(false), rightDrawerTab: ref('features'),
    svgContainer: ref(null), linearSeqs: reactive([]), extractedFeatures: ref([]),
    biologicalFeatures: ref([alpha, beta])
  };
  const sequenceReads = createSequenceReads({
    readFeatures: (features) => {
      reads.push(features);
      return answer(features);
    },
    readSpans: async () => null,
    copyPending: async (pending) => {
      copied = await pending;
      return copied;
    },
    downloadText: (filename) => downloads.push(filename)
  });
  const editor = createOrthogroupEditor({ state: withDrawings(state), memberAction: sequenceReads.memberAction });
  return { editor, state, reads };
};
const label = (editor, kind, original, memberOrNull = null, action = 'copy') =>
  editor.orthogroupCopyFeedbackLabel(group, kind, original, memberOrNull, 'drawer', action);
const readState = (editor, context = 'drawer') => editor.orthogroupMemberRead(context);

test('members carry no sequences, and a member with a catalog identity is available', () => {
  const { editor } = setup(async () => new Map());
  const members = editor.getEnrichedOrthogroupMembers(group);
  for (const enriched of members) {
    for (const key of ['nucleotideSequence', 'aminoAcidSequence', 'nucleotide_sequence', 'amino_acid_sequence',
      'sequenceWarnings']) {
      assert.equal(Object.hasOwn(enriched, key), false, key);
    }
  }
  // Both kinds follow one rule: whether a member holds a sequence of a kind is known only by reading it.
  assert.equal(editor.getOrthogroupSequenceCount(group), 2);
  assert.deepEqual(members.map((enriched) => editor.hasOrthogroupMemberSequence(enriched)), [true, true, false]);
  const detail = editor.selectedOrthogroupDetail.value;
  assert.equal(detail.hasSequence, true);
  assert.deepEqual(detail.members.map((enriched) => detail.hasMemberSequence(enriched)), [true, true, false]);
});

test('a member Copy reads that member, shows the read in place, and copies the reply', async () => {
  let release;
  const { editor, reads } = setup((features) => new Promise((resolve) => {
    release = () => resolve(new Map([sequencesOf(alpha, 'ATGAAATGA', 'MK')]));
  }));
  copied = null;
  const [first] = editor.getEnrichedOrthogroupMembers(group);
  const copying = editor.copyOrthogroupMemberSequence(first, 'nt', group);
  assert.deepEqual(reads, [[{ recordKey: 'rec-1', biologicalFeatureId: 'alpha', bases: 12 }]]);
  assert.equal(readState(editor).status, 'reading');
  assert.equal(readState(editor, 'feature-popup').status, '');
  release();
  assert.equal(await copying, true);
  assert.equal(copied, '>REC:1-12 ALPHA\nATGAAATGA\n');
  assert.equal(readState(editor).status, '');
  assert.equal(label(editor, 'nt', 'Copy nt', first), 'Copied!');
});

test('a group FASTA reads every member once and skips members without a sequence', async () => {
  const { editor, reads } = setup(async () => new Map([sequencesOf(alpha, 'ATG', 'MK'), sequencesOf(beta, 'TTT', '')]));
  copied = null;
  assert.equal(await editor.copyOrthogroupSequences(group, 'aa'), true);
  assert.equal(reads.length, 1);
  assert.deepEqual(reads[0].map((feature) => feature.biologicalFeatureId), ['alpha', 'beta']);
  assert.equal(copied, '>ALPHA ALPHA\nMK\n');
  assert.equal(label(editor, 'aa', 'Copy aa'), 'Copied!');
});

test('a read with no sequence copies nothing and says so in place', async () => {
  const { editor } = setup(async () => new Map([sequencesOf(alpha, '', ''), sequencesOf(beta, '', '')]));
  copied = null;
  assert.equal(await editor.copyOrthogroupSequences(group, 'nt'), false);
  assert.equal(copied, null);
  assert.equal(readState(editor).status, 'empty');
  assert.equal(label(editor, 'nt', 'Copy nt'), 'Copy nt');
});

test('a stale reply writes nothing', async () => {
  const { editor } = setup(async () => null);
  copied = null;
  assert.equal(await editor.copyOrthogroupSequences(group, 'nt'), false);
  assert.equal(copied, null);
  assert.equal(label(editor, 'nt', 'Copy nt'), 'Copy nt');
  assert.equal(readState(editor).status, '');
  downloads.length = 0;
  assert.equal(await editor.downloadOrthogroupSequences(group, 'nt'), false);
  assert.deepEqual(downloads, []);
});

test('a failed read shows its error with Retry in place, and Retry reads again', async () => {
  let fail = true;
  const { editor, reads } = setup(async () => {
    if (fail) throw diagnosticError('SEQUENCE_REQUIRES_GENERATE', {}, { stage: 'helper', operation: 'readFeatureSequences' });
    return new Map([sequencesOf(alpha, 'ATG', 'M')]);
  });
  copied = null;
  assert.equal(await editor.copyOrthogroupSequences(group, 'nt', 'feature-popup'), false);
  assert.equal(copied, null);
  const shown = readState(editor, 'feature-popup');
  assert.equal(shown.status, 'failed');
  assert.equal(shown.error?.code, 'SEQUENCE_REQUIRES_GENERATE');
  assert.equal(typeof shown.retry, 'function');
  assert.equal(readState(editor).status, '');
  fail = false;
  assert.equal(await shown.retry(), true);
  assert.equal(reads.length, 2);
  assert.equal(copied, '>REC:1-12 ALPHA\nATG\n');
  assert.equal(readState(editor, 'feature-popup').status, '');
});

test('Download reads through the port, names the file, and shows a failure in place', async () => {
  const { editor } = setup(async () => new Map([sequencesOf(beta, 'ACGT', 'MP')]));
  downloads.length = 0;
  const [, second] = editor.getEnrichedOrthogroupMembers(group);
  assert.equal(await editor.downloadOrthogroupMemberSequence(second, 'aa', group), true);
  assert.deepEqual(downloads, ['og_1_BETA_aa.faa']);
  assert.equal(label(editor, 'aa', 'Download aa', second, 'download'), 'Download aa');
  const empty = setup(async () => new Map([sequencesOf(alpha, '', ''), sequencesOf(beta, '', '')]));
  assert.equal(await empty.editor.downloadOrthogroupSequences(group, 'nt'), false);
  assert.equal(readState(empty.editor).status, 'empty');
  const failing = setup(async () => { throw new Error('worker down'); });
  assert.equal(await failing.editor.downloadOrthogroupSequences(group, 'nt'), false);
  assert.equal(readState(failing.editor).status, 'failed');
  assert.ok(failing.editor.orthogroupMemberRead('drawer').error);
  assert.equal(label(failing.editor, 'nt', 'Download nt', null, 'download'), 'Download nt');
});

test('a read the popup no longer shows writes nothing to the state it shows', async () => {
  let reject;
  const { editor, state } = setup(() => new Promise((_, rejectRead) => { reject = rejectRead; }));
  const reading = editor.copyOrthogroupSequences(group, 'nt', 'feature-popup');
  assert.equal(readState(editor, 'feature-popup').status, 'reading');
  // Another feature is clicked: the popup shows the new feature, with no read.
  state.clickedFeature.value = { feat: {} };
  assert.equal(readState(editor, 'feature-popup').status, '');
  reject(new Error('worker down'));
  assert.equal(await reading, false);
  assert.equal(readState(editor, 'feature-popup').status, '');
});

test('the editor reports no read failure to the page error log', async () => {
  const { readFileSync } = await import('node:fs');
  const source = readFileSync(new URL('../../gbdraw/web/js/app/orthogroups.js', import.meta.url), 'utf8');
  assert.doesNotMatch(source, /errorLog/);
  assert.doesNotMatch(source, /copyPendingTextToClipboard|downloadTextFile/);
});

test('the orthogroup panel and the match popup take a member FASTA header from one owner', async () => {
  const { readFileSync } = await import('node:fs');
  const { buildFeatureSequenceFastas, memberFastaFeature } = await import('../../gbdraw/web/js/services/feature-sequence-fasta.js');
  // The catalog row is read first, the member fills what it lacks, and the header's protein ID is the one the
  // popups show (an internal ID never).
  const header = memberFastaFeature(
    { recordId: 'REC', start: 0, end: 12, strand: -1, sourceProteinId: 'WP_1.1', product: 'kinase' },
    { recordKey: 'rec-1', biologicalFeatureId: 'alpha', sourceProteinId: 'h_abcdefghijklmnopqrstuvwxyz', gene: 'abc' },
    '1-12(-)'
  );
  assert.equal(header.record_id, 'REC');
  assert.equal(header.strand, '-');
  assert.equal(header.protein_id, 'WP_1.1');
  assert.match(buildFeatureSequenceFastas(header, { nucleotideSequence: 'ATG', aminoAcidSequence: 'M' }).aminoAcidFasta,
    /^>WP_1\.1 /);
  assert.equal(header.product, 'kinase');
  assert.equal(header.gene, 'abc');
  for (const file of ['gbdraw/web/js/app/orthogroups.js', 'gbdraw/web/js/services/pairwise-match-popup.js']) {
    const source = readFileSync(new URL(`../../${file}`, import.meta.url), 'utf8');
    assert.doesNotMatch(source, /const (buildMemberFeaturePayload|normalizeMemberStrand) =/, file);
  }
});
