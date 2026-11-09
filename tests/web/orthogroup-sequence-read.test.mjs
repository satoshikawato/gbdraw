// Copy and Download of a Similarity group member's or a whole group's FASTA read the sequences from Python through
// the catalog-sequence door (D5): the members carry none, a member is available when its catalog row has an
// identity, one read asks every member once, a stale reply writes nothing, an empty read says so, and a failure
// says so and reaches the page's error log.
import assert from 'node:assert/strict';
import test from 'node:test';
import { withDrawings } from './helpers/drawing-state.mjs';
import { loadPageVue } from './helpers/vue-runtime.mjs';

const Vue = loadPageVue();
globalThis.window = { Vue, setTimeout, clearTimeout };
globalThis.isSecureContext = true;
let copied = null;
Object.defineProperty(globalThis, 'navigator', {
  configurable: true,
  value: { clipboard: { writeText: async (value) => { copied = String(value); } } }
});
const downloads = [];
globalThis.document = {
  createElement: () => ({ addEventListener() {}, click() { downloads.push(this.download); } }),
  body: { appendChild() {} }
};
globalThis.URL.createObjectURL = () => 'blob:orthogroup-sequence-read';
globalThis.URL.revokeObjectURL = () => {};

const { createOrthogroupEditor } = await import('../../gbdraw/web/js/app/orthogroups.js');
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
    biologicalFeatures: ref([alpha, beta]), errorLog: ref(null)
  };
  const editor = createOrthogroupEditor({
    state: withDrawings(state),
    readFeatureSequences: (features) => {
      reads.push(features);
      return answer(features);
    }
  });
  return { editor, state, reads };
};
const label = (editor, kind, original, memberOrNull = null, action = 'copy') =>
  editor.orthogroupCopyFeedbackLabel(group, kind, original, memberOrNull, 'drawer', action);

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

test('a member Copy reads that member through the port and copies the reply', async () => {
  let release;
  const { editor, reads } = setup((features) => new Promise((resolve) => {
    release = () => resolve(new Map([sequencesOf(alpha, 'ATGAAATGA', 'MK')]));
  }));
  copied = null;
  const [first] = editor.getEnrichedOrthogroupMembers(group);
  const copying = editor.copyOrthogroupMemberSequence(first, 'nt', group);
  assert.deepEqual(reads, [[{ recordKey: 'rec-1', biologicalFeatureId: 'alpha', bases: 12 }]]);
  assert.equal(label(editor, 'nt', 'Copy nt', first), 'Reading the sequence…');
  release();
  assert.equal(await copying, true);
  assert.equal(copied, '>REC:1-12 ALPHA\nATGAAATGA\n');
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

test('a read with no sequence copies nothing and says so', async () => {
  const { editor } = setup(async () => new Map([sequencesOf(alpha, '', ''), sequencesOf(beta, '', '')]));
  copied = null;
  assert.equal(await editor.copyOrthogroupSequences(group, 'nt'), false);
  assert.equal(copied, null);
  assert.equal(label(editor, 'nt', 'Copy nt'), 'No sequence available');
});

test('a stale reply writes nothing', async () => {
  const { editor, state } = setup(async () => null);
  copied = null;
  assert.equal(await editor.copyOrthogroupSequences(group, 'nt'), false);
  assert.equal(copied, null);
  assert.equal(label(editor, 'nt', 'Copy nt'), 'Copy nt');
  assert.equal(state.errorLog.value, null);
  downloads.length = 0;
  assert.equal(await editor.downloadOrthogroupSequences(group, 'nt'), false);
  assert.deepEqual(downloads, []);
});

test('a failed read says Copy failed and reports the error', async () => {
  const { editor, state } = setup(async () => {
    throw diagnosticError('SEQUENCE_REQUIRES_GENERATE', {}, { stage: 'helper', operation: 'readFeatureSequences' });
  });
  copied = null;
  assert.equal(await editor.copyOrthogroupSequences(group, 'nt'), false);
  assert.equal(copied, null);
  assert.equal(label(editor, 'nt', 'Copy nt'), 'Copy failed');
  assert.equal(state.errorLog.value?.code, 'SEQUENCE_REQUIRES_GENERATE');
});

test('Download reads through the port, names the file, and gives the same feedback', async () => {
  const { editor, state } = setup(async () => new Map([sequencesOf(beta, 'ACGT', 'MP')]));
  downloads.length = 0;
  const [, second] = editor.getEnrichedOrthogroupMembers(group);
  assert.equal(await editor.downloadOrthogroupMemberSequence(second, 'aa', group), true);
  assert.deepEqual(downloads, ['og_1_BETA_aa.faa']);
  assert.equal(label(editor, 'aa', 'Download aa', second, 'download'), 'Download aa');
  const empty = setup(async () => new Map([sequencesOf(alpha, '', ''), sequencesOf(beta, '', '')]));
  assert.equal(await empty.editor.downloadOrthogroupSequences(group, 'nt'), false);
  assert.equal(label(empty.editor, 'nt', 'Download nt', null, 'download'), 'No sequence available');
  const failing = setup(async () => { throw new Error('worker down'); });
  assert.equal(await failing.editor.downloadOrthogroupSequences(group, 'nt'), false);
  assert.equal(label(failing.editor, 'nt', 'Download nt', null, 'download'), 'Download failed');
  assert.ok(failing.state.errorLog.value);
  assert.equal(state.errorLog.value, null);
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
