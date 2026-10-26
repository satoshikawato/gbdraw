// UX-1b U1: the views of the similarity popups (UX-B01 collinear block,
// UX-B02 similarity-group ribbon) built by app/pairwise-match-popup.js, and
// OV-420 (a pairwise match with no similarity group showed a "Similarity
// group" section holding the query product).
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import test from 'node:test';

import {
  buildMatchPopupPayload,
  buildPairwiseMatchHoverRows,
  groupCountText
} from '../../gbdraw/web/js/app/pairwise-match-popup.js';
import { SOURCE_LOCALE, setUiLocale } from '../../gbdraw/web/js/utils/i18n.js';

const elementFrom = (attributes) => ({
  style: { fill: '' },
  getAttribute: (name) => attributes[name] ?? ''
});

// Features resolve from an element only with the record and source-feature
// indexes the renderer writes (see pairwise-payload-performance.test.mjs).
const featurePair = (queryRecordId, subjectRecordId) => new Map([
  ['fq1', {
    fileIdx: 0, sourceFeatureIndex: 3, svg_id: 'fq1', stable_feature_id: 'fq1', record_id: queryRecordId,
    start: 99, end: 400, strand: '+', protein_id: 'QUERY.1', product: 'query enzyme'
  }],
  ['fs1', {
    fileIdx: 1, sourceFeatureIndex: 7, svg_id: 'fs1', stable_feature_id: 'fs1', record_id: subjectRecordId,
    start: 199, end: 500, strand: '+', protein_id: 'SUBJECT.1', product: 'subject enzyme'
  }]
]);
const featureAttributes = {
  'data-query-record-index': '0',
  'data-subject-record-index': '1',
  'data-query-feature-index': '3',
  'data-subject-feature-index': '7',
  'data-query-feature-svg-id': 'fq1',
  'data-subject-feature-svg-id': 'fs1'
};


test('OV-420: a pairwise match with no similarity group shows no similarity group', () => {
  const payload = buildMatchPopupPayload(elementFrom({
    'data-gbdraw-match-id': 'comparison1_match1',
    'data-match-kind': 'pairwise',
    'data-query-record-id': 'record_a',
    'data-subject-record-id': 'record_b',
    'data-qstart': '100',
    'data-qend': '400',
    'data-sstart': '200',
    'data-send': '500',
    'data-identity': '91.2',
    ...featureAttributes
  }), { featureLookup: featurePair('record_a', 'record_b') });
  assert.deepEqual(
    payload.sections.map((section) => section.title).filter((title) => title === 'Similarity group'),
    []
  );
  assert.equal(payload.group ?? null, null);
});

test('UX-B02: the group view names the group, counts members, and keeps computed facts technical', () => {
  const members = [
    { recordIndex: 0, recordId: 'record_a', start: 99, end: 400, strand: 1, proteinId: 'QUERY.1', product: 'query enzyme', role: 'anchor', confidence: 'high' },
    { recordIndex: 1, recordId: 'record_b', start: 199, end: 500, strand: 1, proteinId: 'SUBJECT.1', product: 'subject enzyme', role: 'local_paralog', confidence: 'medium' }
  ];
  const orthogroups = [{
    id: 'og_7',
    displayName: 'enzyme',
    description: 'No informative product/gene/note consensus was found.',
    member_count: 2,
    record_coverage_count: 2,
    rbhOrthogroupIds: ['og_7'],
    orthologPaths: [{}, {}],
    members
  }];
  const element = elementFrom({
    'data-gbdraw-match-id': 'comparison1_match2',
    'data-match-kind': 'orthogroup',
    'data-orthogroup-id': 'og_7'
  });
  const featureLookup = featurePair('record_a', 'record_b');
  const payload = buildMatchPopupPayload(element, { featureLookup, orthogroups });
  assert.equal(payload.title, 'og_7:enzyme', 'the hover summary title is unchanged');
  assert.deepEqual(payload.sections, []);
  const { group } = payload;
  assert.equal(group.title, 'enzyme');
  assert.equal(group.badge, 'og_7');
  assert.equal(group.countText, '2 genes in 2 records');
  assert.equal(group.description, '', 'a computed description is not the lead');
  assert.deepEqual(group.technicalRows.map((row) => row.label), [
    'Similarity group ID', 'Description', 'Scope', 'RBH seeds', 'Group paths'
  ]);
  assert.deepEqual(group.memberRows.map((row) => row.roleNote), ['', 'local paralog · medium confidence']);

  const described = buildMatchPopupPayload(element, {
    featureLookup,
    orthogroups,
    orthogroupDescriptionOverrides: { og_7: 'Aminotransferases' }
  }).group;
  assert.equal(described.description, 'Aminotransferases', 'a description the user wrote leads');
  assert.equal(described.technicalRows.some((row) => row.label === 'Description'), false);
});

test('UX-B02: the count line names genes and records once', () => {
  assert.equal(groupCountText(5, 5), '5 genes in 5 records');
  assert.equal(groupCountText('1', '1'), '1 gene in 1 record');
  assert.equal(groupCountText(3, ''), '3 genes');
  assert.equal(groupCountText('', ''), '');
});

test('UX-B01: a collinear block leads with one summary and lists gene pairs by group', () => {
  const groups = [
    { id: 'og_2', displayName: 'rsmA', members: [
      { recordIndex: 0, stableFeatureSvgId: 'fq1', recordId: 'record_a', start: 99, end: 400 },
      { recordIndex: 1, stableFeatureSvgId: 'fs1', recordId: 'record_b', start: 199, end: 500 }
    ] }
  ];
  const payload = buildMatchPopupPayload(elementFrom({
    'data-gbdraw-match-id': 'comparison1_match3',
    'data-match-kind': 'collinear',
    'data-collinearity-block-id': 'block_0386',
    'data-collinearity-block-kind': 'cluster',
    'data-collinearity-orientation': 'plus',
    'data-collinearity-color-mode': 'orientation_identity',
    'data-collinearity-block-score': '8448',
    'data-collinearity-anchor-index': '1',
    'data-collinearity-anchor-count': '15',
    'data-alignment-length': '5831',
    'data-identity': '76.0591457726',
    'data-orthogroup-id': 'og_2',
    'data-query-record-id': 'AP027133.1',
    'data-subject-record-id': 'AP027132.1',
    'data-qstart': '10124',
    'data-qend': '28017',
    'data-sstart': '16057',
    'data-send': '34176',
    ...featureAttributes
  }), { featureLookup: featurePair('AP027133.1', 'AP027132.1'), orthogroups: groups });
  const { block } = payload;
  assert.deepEqual(payload.sections, []);
  assert.equal(block.query.text, 'AP027133.1 10,124..28,017 (17.9 kb)');
  assert.equal(block.subject.text, 'AP027132.1 16,057..34,176 (18.1 kb)');
  assert.deepEqual(block.facts, ['Same orientation', '76% identity', '1 gene pair']);
  assert.equal(
    block.summaryText,
    'Collinear block: AP027133.1 10,124..28,017 (17.9 kb) vs AP027132.1 16,057..34,176 (18.1 kb); Same orientation; 76% identity; 1 gene pair'
  );
  assert.deepEqual(block.genePairs.map((pair) => [pair.query?.label, pair.subject?.label, pair.group?.id]), [
    ['QUERY.1', 'SUBJECT.1', 'og_2']
  ]);
  assert.deepEqual(Object.fromEntries(block.detailRows.map((row) => [row.label, row.value])), {
    'Block ID': 'block_0386',
    Kind: 'cluster',
    Orientation: 'plus',
    'Color mode': 'orientation_identity',
    'Average identity': '76.06%',
    'Aligned length': '5,831 aa',
    'Block score': '8448',
    Anchor: '1 / 15'
  });
  assert.match(block.note, /^Block envelopes may include/);
  assert.equal(block.pairsCopyText.split('\n').length, 2);
  assert.equal(payload.blockOrthogroups[0].title, 'rsmA');
});

// OV-322: the orientation of a collinear block is shown as its own text unless it is
// `plus` or `minus`. A hand-edited Session SVG can name an inherited member of the
// text table (`constructor`, `__proto__`), which is not an orientation.
test('OV-322: a collinear block shows an orientation other than plus or minus as its own text', () => {
  const factsFor = (orientation) => buildMatchPopupPayload(elementFrom({
    'data-gbdraw-match-id': 'comparison1_match3',
    'data-match-kind': 'collinear',
    'data-collinearity-block-id': 'block_0386',
    'data-collinearity-orientation': orientation,
    'data-identity': '76.0591457726',
    'data-query-record-id': 'AP027133.1',
    'data-subject-record-id': 'AP027132.1',
    'data-qstart': '10124',
    'data-qend': '28017',
    'data-sstart': '16057',
    'data-send': '34176'
  })).block.facts;
  assert.deepEqual(factsFor('plus')[0], 'Same orientation');
  assert.deepEqual(factsFor('minus')[0], 'Inverted');
  assert.deepEqual(factsFor('unknown')[0], 'unknown');
  assert.deepEqual(factsFor('constructor')[0], 'constructor');
  assert.deepEqual(factsFor('__proto__')[0], '__proto__');
});

// OV-526 (same pattern as OV-322): the name and description overrides of a similarity group are
// read by its ID, which a hand-edited Session SVG can set to an inherited member of a
// plain object (`constructor`). Only an own entry is an override.
test('OV-526: a similarity group ID that names an inherited member has no override', () => {
  const payload = buildMatchPopupPayload(elementFrom({
    'data-gbdraw-match-id': 'og_match',
    'data-match-kind': 'orthogroup',
    'data-orthogroup-id': 'constructor'
  }), { orthogroupNameOverrides: {}, orthogroupDescriptionOverrides: {} });
  assert.equal(payload.group.displayName, '');
  assert.equal(payload.group.description, '');
  const named = buildMatchPopupPayload(elementFrom({
    'data-gbdraw-match-id': 'og_match',
    'data-match-kind': 'orthogroup',
    'data-orthogroup-id': 'og_1'
  }), { orthogroupNameOverrides: { og_1: 'rsmA' }, orthogroupDescriptionOverrides: new Map([['og_1', 'A note']]) });
  assert.equal(named.group.displayName, 'rsmA');
  assert.equal(named.group.description, 'A note');
});

test('a popup built before a language switch keeps its hover rows: sections and rows are found by identity', async () => {
  const element = elementFrom({
    'data-gbdraw-match-id': 'comparison1_match1',
    'data-match-kind': 'pairwise',
    'data-query-record-id': 'record_a',
    'data-subject-record-id': 'record_b',
    'data-qstart': '100',
    'data-qend': '400',
    'data-sstart': '200',
    'data-send': '500',
    'data-identity': '91.2',
    ...featureAttributes
  });
  const options = { featureLookup: featurePair('record_a', 'record_b') };
  const catalog = Object.freeze({
    Summary: 'S-ja', Alignment: 'A-ja', Identity: 'I-ja', 'Query interval': 'QI-ja',
    'Subject interval': 'SI-ja', Kind: 'K-ja', Query: 'Q-ja', Subject: 'Sub-ja'
  });
  const english = buildMatchPopupPayload(element, options);
  try {
    await setUiLocale('ja', async () => ({ default: catalog }));
    const japanese = buildMatchPopupPayload(element, options);
    assert.deepEqual(japanese.sections.slice(0, 2).map((entry) => [entry.id, entry.title]), [['Summary', 'S-ja'], ['Alignment', 'A-ja']]);
    const expected = [['Kind', 'K-ja', 'pairwise'], ['Identity', 'I-ja', '91.2'], ['Query', 'Q-ja', '100..400'], ['Subject', 'Sub-ja', '200..500']];
    const rowsOf = (payload) => buildPairwiseMatchHoverRows(payload).map((row) => [row.id, row.label, row.value]);
    assert.deepEqual(rowsOf(english), expected, 'an English payload shown after the switch');
    assert.deepEqual(rowsOf(japanese), expected);
  } finally {
    await setUiLocale(SOURCE_LOCALE);
  }
  assert.deepEqual(
    buildPairwiseMatchHoverRows(english).map((row) => [row.label, row.value]),
    [['Kind', 'pairwise'], ['Identity', '91.2'], ['Query', '100..400'], ['Subject', '200..500']]
  );
});

test('the collinear block summary shows "vs" in the muted style before the subject record', () => {
  const html = readFileSync(new URL('../../gbdraw/web/index.html', import.meta.url), 'utf8');
  const summary = html.slice(html.indexOf('data-collinear-summary'));
  const subjectLine = summary.split('\n').find((line) => line.includes('block.subject.recordId'));
  assert.match(subjectLine, /<span class="match-muted">\{\{ t\('vs'\) \}\}\{\{ ' ' \}\}<\/span><span class="font-semibold">\{\{ clickedPairwiseMatch\.block\.subject\.recordId \}\}<\/span>/);
});

test('OV-1390: member roles and confidence levels show in the UI language; an unknown value shows as written', async () => {
  const members = [
    { recordIndex: 0, recordId: 'record_a', start: 99, end: 400, strand: 1, proteinId: 'QUERY.1', role: 'anchor', confidence: 'high' },
    { recordIndex: 1, recordId: 'record_b', start: 199, end: 500, strand: 1, proteinId: 'SUBJECT.1', role: 'local_paralog', confidence: 'medium' },
    { recordIndex: 1, recordId: 'record_b', start: 599, end: 900, strand: 1, proteinId: 'SUBJECT.2', role: 'future_role', confidence: 'low' }
  ];
  const orthogroups = [{ id: 'og_7', displayName: 'enzyme', member_count: 3, record_coverage_count: 2, members }];
  const groupElement = elementFrom({
    'data-gbdraw-match-id': 'comparison1_match2',
    'data-match-kind': 'orthogroup',
    'data-orthogroup-id': 'og_7'
  });
  const matchElement = elementFrom({
    'data-gbdraw-match-id': 'comparison1_match1',
    'data-match-kind': 'pairwise',
    'data-qstart': '100',
    'data-qend': '400',
    'data-sstart': '200',
    'data-send': '500',
    'data-query-orthogroup-role': 'low_confidence',
    'data-subject-orthogroup-role': 'coortholog',
    'data-query-orthogroup-confidence': 'low',
    'data-subject-orthogroup-confidence': 'future_level',
    ...featureAttributes
  });
  const options = { featureLookup: featurePair('record_a', 'record_b'), orthogroups };
  const shown = () => {
    const alignment = buildMatchPopupPayload(matchElement, options).sections.find((entry) => entry.id === 'Alignment');
    const value = (id) => alignment.rows.find((row) => row.id === id)?.value;
    return {
      notes: buildMatchPopupPayload(groupElement, options).group.memberRows.map((row) => row.roleNote),
      rows: ['Query member role', 'Subject member role', 'Query confidence', 'Subject confidence'].map(value)
    };
  };
  assert.deepEqual(shown(), {
    notes: ['', 'local paralog · medium confidence', 'future_role · low confidence'],
    rows: ['low confidence', 'coortholog', 'low', 'future_level']
  });
  try {
    await setUiLocale('ja');
    assert.deepEqual(shown(), {
      notes: ['', 'レコード内パラログ · 信頼度：中', 'future_role · 信頼度：低'],
      // ja keeps the English word coortholog: no established Japanese term (catalog decision, no coined terms).
      rows: ['低信頼度', 'coortholog', '低', 'future_level']
    });
  } finally {
    await setUiLocale(SOURCE_LOCALE);
  }
});
