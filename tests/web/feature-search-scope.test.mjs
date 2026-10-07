import assert from 'node:assert/strict';
import { test } from 'node:test';
import * as featureUtils from '../../gbdraw/web/js/services/feature-utils.js';
import { runFeatureSearch } from '../../gbdraw/web/js/app/feature-search/search-core.js';
import { buildMatchPopupPayload } from '../../gbdraw/web/js/app/pairwise-match-popup.js';
import { STANDALONE_INTERACTIVE_SCRIPT } from '../../gbdraw/web/js/services/standalone-interactivity-assets.js';
import {
  FEATURE_CATALOG_SCHEMA,
  featureStateFromCatalog,
  validateFeatureCatalog
} from '../../gbdraw/web/js/services/feature-catalog.js';

// Feature Search and the Interactive SVG runtime are two implementations of one
// contract (the runtime cannot import modules). These cases bind them together.
const features = [
  {
    svg_id: 'f-cytb', type: 'CDS', record_id: 'REC', gene: 'CYTB', product: 'cytochrome b',
    start: 700, end: 1000, strand: '-',
    location_parts: [
      { start: 900, end: 1000, strand: '-', display: '901..1000' },
      { start: 700, end: 800, strand: '-', display: '701..800' }
    ],
    qualifiers: { gene: ['CYTB'], product: ['cytochrome b'], translation: ['MKKK'] },
    nucleotide_sequence: 'ATGAAAAAAAAATAA', amino_acid_sequence: 'MKKK'
  },
  {
    svg_id: 'f-trnf', type: 'tRNA', record_id: 'REC', product: 'tRNA-Phe',
    start: 1010, end: 1031, strand: '+',
    location_parts: [{ start: 1010, end: 1031, strand: '+', display: '1011..1031' }],
    qualifiers: { product: ['tRNA-Phe'] },
    nucleotide_sequence: 'GTTTATGTAGCTTACCTCCAA'
  },
  {
    svg_id: 'f-dnaa', type: 'CDS', record_id: 'REC', gene: 'dnaA', product: 'replication initiator',
    start: 1100, end: 1109, strand: '+',
    location_parts: [{ start: 1100, end: 1109, strand: '+', display: '1101..1109' }],
    qualifiers: { gene: ['dnaA'], product: ['replication initiator'] },
    nucleotide_sequence: 'ATGCCCTGA', amino_acid_sequence: 'MP'
  },
  {
    svg_id: 'f-gyra', type: 'CDS', record_id: 'REC', gene: 'gyrA', product: 'DNA gyrase subunit A',
    start: 1200, end: 1212, strand: '+',
    location_parts: [{ start: 1200, end: 1212, strand: '+', display: '1201..1212' }],
    qualifiers: { gene: ['gyrA'], product: ['DNA gyrase subunit A'], translation: ['MSDN'] },
    nucleotide_sequence: 'ATGTCCGACAAC', amino_acid_sequence: 'MSDN'
  },
  {
    svg_id: 'f-hns', type: 'CDS', record_id: 'REC', gene: 'hns', product: 'DNA-binding protein H-NS',
    start: 1300, end: 1327, strand: '+',
    location_parts: [{ start: 1300, end: 1327, strand: '+', display: '1301..1327' }],
    qualifiers: { gene: ['hns'], product: ['DNA-binding protein H-NS'], translation: ['MSGYRAKLE'] },
    nucleotide_sequence: 'ATGTCCGGCTACCGCGCCAAGCTGGAG'
  },
  {
    svg_id: 'f-origin', type: 'CDS', record_id: 'REC', gene: 'oriX', product: 'origin spanning protein',
    start: 0, end: 4000, strand: '+',
    location_parts: [
      { start: 3900, end: 4000, strand: '+', display: '3901..4000' },
      { start: 0, end: 200, strand: '+', display: '1..200' }
    ],
    qualifiers: { gene: ['oriX'], product: ['origin spanning protein'] }
  },
  {
    svg_id: 'f-simple', type: 'CDS', record_id: 'REC', gene: 'dupA', product: 'duplicate protein',
    start: 300, end: 600, strand: '+',
    location_parts: [{ start: 300, end: 600, strand: '+', display: '301..600' }],
    qualifiers: { gene: ['dupA'], product: ['duplicate protein'] }
  }
];
const renderedFeatureIds = new Set(features.map((feature) => feature.svg_id));
const appSearch = (query, field, popupMode = 'rich', options = {}) => runFeatureSearch({
  features,
  renderedFeatureIds,
  query,
  field,
  popupMode,
  ...options
});

const locationCases = [
  [features.find((feature) => feature.svg_id === 'f-origin'), '3901..4000, 1..200 (+)', '300 bp'],
  [features.find((feature) => feature.svg_id === 'f-cytb'), '901..1000, 701..800 (-)', '200 bp'],
  [features.find((feature) => feature.svg_id === 'f-simple'), '301..600 (+)', '300 bp'],
  [{ start: 300, end: 600, strand: '+' }, '301..600 (+)', '300 bp'],
  [{ start: 10, end: 20, strand: '-', location_parts: [{ start: 10, end: 20, strand: '-' }] }, '11..20 (-)', '10 bp'],
  [{ strand: '+' }, '', '']
];

const embeddedSlice = (startMarker, endMarker) => {
  const start = STANDALONE_INTERACTIVE_SCRIPT.indexOf(startMarker);
  const end = STANDALONE_INTERACTIVE_SCRIPT.indexOf(endMarker, start);
  assert.ok(start >= 0 && end > start, `missing embedded source ${startMarker}`);
  return STANDALONE_INTERACTIVE_SCRIPT.slice(start, end);
};
const embeddedFunction = (name) => {
  const start = STANDALONE_INTERACTIVE_SCRIPT.indexOf(`function ${name}(`);
  assert.notEqual(start, -1, `missing embedded function ${name}`);
  let depth = 0;
  for (let index = STANDALONE_INTERACTIVE_SCRIPT.indexOf('{', start); index < STANDALONE_INTERACTIVE_SCRIPT.length; index += 1) {
    if (STANDALONE_INTERACTIVE_SCRIPT[index] === '{') depth += 1;
    if (STANDALONE_INTERACTIVE_SCRIPT[index] === '}') depth -= 1;
    if (depth === 0) return STANDALONE_INTERACTIVE_SCRIPT.slice(start, index + 1);
  }
  return assert.fail(`unterminated embedded function ${name}`);
};
const createEmbeddedSearch = (popupMode) => new Function('features', 'popupMode', `
  var orthogroupsById = new Map();
  function consistentTextIdentity() { return { valid: false, value: '' }; }
  function biologicalFeatureForMember() { return null; }
  ${embeddedFunction('normalizeArray')}
  ${embeddedFunction('featureLocationParts')}
  ${embeddedFunction('locationText')}
  ${embeddedFunction('featureLengthText')}
  ${embeddedSlice('var searchFieldIds = {', 'var featuresById = new Map();')}
  ${embeddedSlice('function normalizeSearchText(', 'preparedSearchIndex = buildPreparedSearchIndex();')}
  var index = buildPreparedSearchIndex();
  return {
    locationText: locationText,
    featureLengthText: featureLengthText,
    search: function (query, field, qualifierKey, useRegex) {
      var normalizedField = normalizeSearchField(field);
      var matcher = compileSearchMatcher(query, useRegex);
      if (!matcher.active || matcher.error) return [];
      var candidateIds = index.featureOrder;
      var key = normalizeSearchText(qualifierKey).trim();
      if (normalizedField === 'qualifier-value' && key) {
        candidateIds = index.qualifierFeatureIdsByKey.get(key) || [];
      }
      return candidateIds.filter(function (svgId) {
        return preparedFeatureSearchMatches(
          index.byId.get(svgId), matcher, normalizedField, qualifierKey
        ).length > 0;
      });
    }
  };
`)(features, popupMode);

test('All matches names and annotations, not sequence content or /translation (D-01)', () => {
  for (const [query, expected] of [['CYTB', ['f-cytb']], ['dnaA', ['f-dnaa']], ['gyrA', ['f-gyra']]]) {
    assert.deepEqual(appSearch(query, 'all').matches, expected, query);
  }
  assert.deepEqual(appSearch('polymerase', 'all').matches, []);
  assert.deepEqual(appSearch('H-NS', 'all').matches, ['f-hns']);
  // The dedicated sequence and qualifier fields keep sequence search, IUPAC included.
  assert.ok(appSearch('CYTB', 'nucleotide').matches.includes('f-trnf'));
  assert.deepEqual(appSearch('ATGTCCGGN', 'nucleotide').matches, ['f-hns']);
  assert.deepEqual(appSearch('GYRA', 'amino-acid').matches, ['f-hns']);
  assert.deepEqual(appSearch('GYRA', 'qualifier-value', 'rich', { qualifierKey: 'translation' }).matches, ['f-hns']);
});

test('Location search uses the 1-based INSDC location only (D-16)', () => {
  assert.deepEqual(appSearch('300', 'location').matches, []);
  assert.deepEqual(appSearch('301..600', 'location').matches, ['f-simple']);
  assert.deepEqual(appSearch('3901..4000', 'location').matches, ['f-origin']);
  assert.deepEqual(appSearch('1..200', 'location').matches, ['f-origin']);
  const details = Object.values(appSearch('0', 'location').matchDetails).flat();
  assert.ok(details.length > 0);
  assert.deepEqual([...new Set(details.map((detail) => detail.label))], ['Location']);
});

test('feature location and length come from the location parts', () => {
  for (const [feature, location, length] of locationCases) {
    assert.equal(featureUtils.formatFeatureLocation(feature), location);
    assert.equal(featureUtils.formatFeatureLength(feature), length);
  }
});

test('Interactive SVG search and location text match the app', () => {
  const queries = [
    'CYTB', 'dnaA', 'gyrA', 'GYRA', 'CCTC', 'ATGAAA', 'MSGYRA', '300', '301', '0',
    '3901..4000', '1..200', '(+)', 'H-NS', 'tRNA', 'translation', 'cytochrome'
  ];
  const fields = [
    'all', 'label', 'type', 'record-id', 'location', 'strand', 'qualifier-key',
    'qualifier-value', 'nucleotide', 'amino-acid'
  ];
  for (const popupMode of ['rich', 'simple']) {
    const embedded = createEmbeddedSearch(popupMode);
    for (const query of queries) {
      for (const field of fields) {
        assert.deepEqual(
          embedded.search(query, field, '', false),
          appSearch(query, field, popupMode).matches,
          `${popupMode} ${field} ${query}`
        );
      }
    }
    assert.deepEqual(
      embedded.search('^cyt', 'all', '', true),
      appSearch('^cyt', 'all', popupMode, { useRegex: true }).matches
    );
    for (const [feature, location, length] of locationCases) {
      assert.equal(embedded.locationText(feature), location);
      assert.equal(embedded.featureLengthText(feature), length);
    }
  }
});

// The popup title and the Details Protein ID of the Interactive SVG follow the
// app popup (getFeatureCaption and resolveFeatureProteinId): the app reads the
// admitted catalog, the runtime the catalog entry with its rendered ID.
test('Interactive SVG popup title and Protein ID match the app popup', () => {
  const internal = `f_${'0'.repeat(64)}`;
  const qualifierSets = [
    { gene: ['TRNF'], product: ['tRNA-Phe'], note: ['NAR: 1455'] },
    { gene: ['COX1'], product: ['cytochrome c oxidase subunit I'], protein_id: ['YP_003024028.1'] },
    { gene: ['TRNF'] },
    { locus_tag: ['LOC_0001'], protein_id: [internal] },
    { note: [`${'n'.repeat(49)}😀tail`] },
    {},
    { product: [internal], gene: ['dnaA'] }
  ];
  const ids = qualifierSets.map((_, index) => `f${String(index).padStart(4, '0')}`);
  const catalog = {
    schema: FEATURE_CATALOG_SCHEMA,
    items: [{
      resultIndex: 0,
      resultName: 'diagram.svg',
      recordKeys: ['record-1'],
      features: ids.map((svgId, index) => ({
        svgId,
        recordKey: 'record-1',
        biologicalFeatureId: `b${index}`,
        fillColor: '#abcdef',
        drawnSelector: null
      })),
      biologicalFeatures: qualifierSets.map((qualifiers, index) => ({
        recordKey: 'record-1',
        biologicalFeatureId: `b${index}`,
        record_idx: 0,
        sourceFeatureIndex: index,
        record_id: 'REC',
        type: index === 0 ? 'tRNA' : 'CDS',
        start: 10 * index,
        end: 10 * index + 5,
        strand: 1,
        anchorProfile: { precision: 'exact', operator: 'single', partOrder: 'biological', strand: '+' },
        qualifiers
      })),
      orthogroups: [],
      annotations: [],
      comparisonMatches: []
    }]
  };
  const results = [{ name: 'diagram.svg', content: '<svg />' }];
  const appFeatures = featureStateFromCatalog(validateFeatureCatalog(structuredClone(catalog), results))
    .extractedFeatures;
  const embedded = new Function(`
    ${embeddedFunction('normalizeArray')}
    ${embeddedFunction('getFeatureQualifiers')}
    ${embeddedFunction('isInternalProteinDisplayId')}
    ${embeddedFunction('firstNonInternalDisplayText')}
    ${embeddedFunction('qualifierDisplayValue')}
    ${embeddedFunction('featureCaption')}
    ${embeddedFunction('featureProteinId')}
    return { featureCaption: featureCaption, featureProteinId: featureProteinId };
  `)();
  const runtimeFeatures = catalog.items[0].biologicalFeatures.map((feature, index) => ({
    ...feature,
    svg_id: ids[index]
  }));
  const titles = runtimeFeatures.map((feature) => embedded.featureCaption(feature));
  assert.deepEqual(titles, appFeatures.map((feature) => featureUtils.getFeatureCaption(feature)));
  assert.deepEqual(titles.slice(0, 5), [
    'tRNA-Phe', 'cytochrome c oxidase subunit I', 'TRNF', 'LOC_0001', `${'n'.repeat(49)}😀`
  ]);
  const proteinIds = runtimeFeatures.map((feature) => embedded.featureProteinId(feature, null));
  assert.deepEqual(proteinIds, appFeatures.map((feature) => featureUtils.resolveFeatureProteinId(feature, null)));
  // Only a protein ID names the Protein ID row; the gene or locus tag of a tRNA does not.
  assert.deepEqual(proteinIds, ['', 'YP_003024028.1', '', '', '', '', '']);
  // An edited label is the title in both rules.
  const edited = { ...runtimeFeatures[0], display_label: 'Phe' };
  assert.equal(embedded.featureCaption(edited), 'Phe');
  assert.equal(featureUtils.getFeatureCaption({ ...appFeatures[0], display_label: 'Phe' }), 'Phe');
  // A similarity-group member's protein ID still counts.
  const member = { proteinId: 'WP_000001.1' };
  assert.equal(embedded.featureProteinId(runtimeFeatures[0], member), 'WP_000001.1');
  assert.equal(featureUtils.resolveFeatureProteinId(appFeatures[0], member), 'WP_000001.1');
});

test('match popup feature sections show the shared split location', () => {
  const origin = features.find((feature) => feature.svg_id === 'f-origin');
  const rendered = { ...origin, fileIdx: 0, sourceFeatureIndex: 5, stable_feature_id: 'stable-origin' };
  const attributes = {
    'data-gbdraw-pairwise-match-id': 'match-origin',
    'data-match-kind': 'pairwise',
    'data-query-record-id': 'REC',
    'data-subject-record-id': 'OTHER',
    'data-query-record-index': '0',
    'data-query-feature-index': '5',
    'data-query-stable-feature-svg-id': 'stable-origin',
    'data-query-feature-svg-id': 'f-origin',
    'data-qstart': '3950',
    'data-qend': '100',
    'data-sstart': '1',
    'data-send': '150'
  };
  const payload = buildMatchPopupPayload(
    { style: { fill: '' }, getAttribute: (name) => attributes[name] || '' },
    {
      featureLookup: new Map([['f-origin', rendered]]),
      sourceFeatures: [{ ...rendered, svg_id: 'stable-origin' }]
    }
  );
  const query = payload.sections.find((section) => section.title === 'Query');
  assert.equal(query.rows.find((row) => row.label === 'Location').value, '3901..4000, 1..200 (+)');
  assert.deepEqual(query.featureRows.map((row) => row.location), ['3901..4000, 1..200 (+)']);
  assert.equal(createEmbeddedSearch('rich').locationText(rendered), '3901..4000, 1..200 (+)');
});
