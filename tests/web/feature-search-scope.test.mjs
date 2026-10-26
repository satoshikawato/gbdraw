import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { test } from 'node:test';
import * as featureUtils from '../../gbdraw/web/js/services/feature-utils.js';
import {
  buildFeatureSearchIndex,
  formatSearchMatchLine,
  runFeatureSearch
} from '../../gbdraw/web/js/app/feature-search/search-core.js';
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
  [{ strand: '+' }, '', ''],
  // OV-968: a null or blank start or end is missing, not 0.
  [{ strand: '+', location_parts: [{ start: null, end: 20, strand: '+' }] }, '', ''],
  [{ strand: '+', location_parts: [{ start: '  ', end: 20, strand: '+' }] }, '', ''],
  [{ strand: '+', location_parts: [{ start: 5, end: null, strand: '+' }] }, '', ''],
  [{ start: null, end: 20, strand: '+' }, '', ''],
  // PT-OD23: hand-edited bounds read by the writers' rule (PT-OD14).
  [{ strand: '+', location_parts: [{ start: '0x10', end: 20, strand: '+' }] }, '', ''],
  [{ strand: '+', location_parts: [{ start: 1.5, end: 20, strand: '+' }] }, '2..20 (+)', '19 bp'],
  [{ start: 0, end: 2 ** 53 + 1, strand: '+' }, '', '']
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
  ${embeddedFunction('locationPartBound')}
  ${embeddedFunction('featureLocationParts')}
  ${embeddedFunction('rangeText')}
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

// PT-OD23: the app, the Interactive SVG viewer, and both interactive-SVG writers
// read a start or end by one rule (../fixtures/location-parts-vectors.json).
test('the app and the Interactive SVG read a start or end by the writers\' rule', () => {
  const vectors = JSON.parse(readFileSync(new URL('../fixtures/location-parts-vectors.json', import.meta.url), 'utf8'));
  const embeddedBound = new Function(`${embeddedFunction('locationPartBound')} return locationPartBound;`)();
  /** @type {[unknown, number | null][]} */
  const values = [
    [true, null], ['0x10', null], ['1_0', null], ['\u0661', null], [1.5, 1], [2 ** 53 + 1, null],
    [' 20 ', 20], ['-0.5', 0], [undefined, null]
  ];
  for (const { input, expected } of vectors.featureCoordinateCases) {
    for (const field of ['start', 'end']) values.push([input[field], expected[field]]);
  }
  for (const [value, expected] of values) {
    const label = String(value);
    assert.equal(featureUtils.locationPartNumber(value), expected, label);
    assert.deepEqual(featureUtils.locationPartBound({ start: value }, 'start'), expected ?? NaN, label);
    assert.deepEqual(embeddedBound({ start: value }, 'start'), expected ?? NaN, label);
  }
});

// OV-1010 / OV-1017: the runtime shows a member's location as the app popup
// does: a null or blank bound is missing (not 0), a bound that does not read
// shows its stored text, and one side shows alone when the other is empty.
// The app text is `locationRangeText` plus the popup's strand rule
// (pairwise-match-popup.js memberLocationText).
test('Interactive SVG member location text equals the app member text', () => {
  const viewerText = new Function(`
    ${embeddedFunction('locationPartBound')}
    ${embeddedFunction('rawBoundText')}
    ${embeddedFunction('strandText')}
    ${embeddedFunction('memberLocationText')}
    return memberLocationText;
  `)();
  const appStrand = (strand) => {
    const text = String(strand).trim();
    if (text === '-1') return '-';
    if (text === '1') return '+';
    return text;
  };
  const appText = (member) => {
    const range = featureUtils.locationRangeText(member);
    const strand = appStrand(member.strand);
    return range && strand ? `${range} (${strand})` : range;
  };
  assert.equal(viewerText({ start: 9, end: 18, strand: '-1' }), '10..18 (-)');
  const bounds = [9, 1.5, null, '', '  ', '0x10', '20', 1e400, 2 ** 53];
  for (const start of bounds) {
    for (const end of [38, undefined]) {
      for (const strand of ['-1', '1', '+', '-', '']) {
        const member = end === undefined ? { start, strand } : { start, end, strand };
        assert.equal(viewerText(member), appText(member), JSON.stringify([String(start), end, strand]));
      }
    }
  }
});

// OV-1015: the Similarity group member rows (feature popup table, Similarity
// drawer) and the member popup show one location text, read with the shared
// reader. The drawer groups digits (en-US); the table does not.
test('Similarity group member rows read their location with the shared reader', () => {
  const { locationRangeText } = featureUtils;
  for (const [start, plain, grouped] of [
    [10, '11..20', '11..20'],
    [null, '20', '20'],
    [1.5, '2..20', '2..20'],
    ['0x10', '0x10..20', '0x10..20'],
    [100, '101..20', '101..20'],
    [1233, '1234..20', '1,234..20']
  ]) {
    assert.equal(locationRangeText({ start, end: 20 }), plain, JSON.stringify(start));
    assert.equal(locationRangeText({ start, end: 20 }, 'en-US'), grouped, JSON.stringify(start));
  }
  assert.equal(locationRangeText({ start: 9999, end: 20000 }, 'en-US'), '10,000..20,000');
  assert.equal(locationRangeText({ start: null, end: undefined }), '');
  assert.equal(locationRangeText(null), '');
  const html = readFileSync(new URL('../../gbdraw/web/index.html', import.meta.url), 'utf8');
  assert.doesNotMatch(html, /Number\(member\.(start|end)\)/);
  assert.equal(html.match(/\{\{ orthogroupMemberRangeText\(member\) \}\}/g)?.length, 1);
  assert.equal(html.match(/\{\{ orthogroupMemberRangeText\(member, 'en-US'\) \}\}/g)?.length, 1);
});

// OV-211: the active match names its record, in the app and the Interactive SVG.
test('the match line names the record of the active match', () => {
  const embeddedLine = new Function(`
    ${embeddedFunction('collapseWhitespace')}
    ${embeddedFunction('searchMatchSnippet')}
    ${embeddedFunction('formatSearchMatchDetail')}
    ${embeddedFunction('formatSearchMatchLine')}
    return formatSearchMatchLine;
  `)();
  const gene = { label: 'Qualifier gene', value: 'dnaA', match: 'dnaA' };
  const product = { label: 'Product', value: 'chromosomal replication initiator protein DnaA', match: 'DnaA' };
  for (const [feature, details, expected] of [
    [{ record_id: 'NC_022349.1' }, [gene], 'Record NC_022349.1 · Matched Qualifier gene: dnaA'],
    [{ record_id: 'NC_022349.1' }, [product, gene], 'Record NC_022349.1 · Matched Product: ...ation initiator protein DnaA'],
    [{ record_id: '' }, [gene], 'Matched Qualifier gene: dnaA'],
    [null, [gene], 'Matched Qualifier gene: dnaA'],
    [{ record_id: 'NC_022349.1' }, [], '']
  ]) {
    assert.equal(formatSearchMatchLine(feature, details), expected);
    assert.equal(embeddedLine(feature, details), expected);
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
    ${embeddedFunction('locationPartBound')}
    ${embeddedFunction('featureLocationParts')}
    ${embeddedFunction('rangeText')}
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
  // An unnamed feature falls back to its 1-based range (start is 0-based), in both rules.
  assert.equal(titles[5], 'CDS at 51..55');
  const split = { ...runtimeFeatures[5], location_parts: [{ start: 50, end: 52 }, { start: 53, end: 55, display: '54..55' }] };
  assert.equal(embedded.featureCaption(split), 'CDS at 51..52, 54..55');
  assert.equal(featureUtils.getFeatureCaption({ ...appFeatures[5], location_parts: split.location_parts }), 'CDS at 51..52, 54..55');
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
    'data-gbdraw-match-id': 'match-origin',
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

test('match popup span rows use the input-file coordinates of the Interval rows (CI-05)', () => {
  // Record 0 is source 3000..8000 (crop); record 1 is source 1000..12000
  // reverse-complemented. data-qstart/qend/sstart/send are record-local.
  const groups = {
    0: { 'data-gbdraw-record-source-start': '3000', 'data-gbdraw-record-source-end': '8000', 'data-gbdraw-record-source-step': '1' },
    1: { 'data-gbdraw-record-source-start': '1000', 'data-gbdraw-record-source-end': '12000', 'data-gbdraw-record-source-step': '-1' }
  };
  const svg = { querySelectorAll: () => Object.entries(groups).map(([index, group]) => ({
    getAttribute: (name) => (name === 'data-gbdraw-record-index' ? index : group[name] ?? null)
  })) };
  const attributes = {
    'data-gbdraw-match-id': 'm1', 'data-match-kind': 'pairwise', 'data-source-index': '0',
    'data-query-record-id': 'Q', 'data-subject-record-id': 'S',
    'data-query-record-index': '0', 'data-subject-record-index': '1',
    'data-qstart': '1', 'data-qend': '4', 'data-sstart': '4', 'data-send': '1'
  };
  const sequences = { 'linear:record:0': 'ACGTACGT', 'linear:record:1': 'TTGGCCAA' };
  const payload = buildMatchPopupPayload(
    { ownerSVGElement: svg, style: { fill: '' }, getAttribute: (name) => attributes[name] || '' },
    { resolveSequenceSource: (key) => ({ source: { recordId: key, sequence: sequences[key] } }) }
  );
  const rows = Object.fromEntries(payload.sections.flatMap((section) => section.rows).map((row) => [row.label, row.value]));
  const spans = payload.sequenceBundle.entries.map(({ span }) => `${span.start}..${span.end}`);
  assert.deepEqual(spans, [rows['Query interval'], rows['Subject interval']]);
  assert.deepEqual(spans, ['3000..3003', '11997..12000']);
  // The bases and the FASTA headers are unchanged.
  assert.match(payload.sequenceBundle.entries[0].fasta, /coords=3000\.\.3003\|strand=\+\nACGT\n$/);
  // A record shown as in its file (rotate only) keeps its coordinates.
  const plain = buildMatchPopupPayload(
    { ownerSVGElement: { querySelectorAll: () => [] }, style: { fill: '' }, getAttribute: (name) => attributes[name] || '' },
    { resolveSequenceSource: (key) => ({ source: { recordId: key, sequence: sequences[key] } }) }
  );
  assert.deepEqual(plain.sequenceBundle.entries.map(({ span }) => `${span.start}..${span.end}`), ['1..4', '4..1']);
});

// A search index reads each similarity-group member once, not once per feature
// and field (OV-914); a feature still takes the protein ID of its one member.
test('the search index reads each group member a bounded number of times (OV-914)', () => {
  const build = (count, { ambiguous = false } = {}) => {
    let memberReads = 0;
    const member = (index) => new Proxy({
      record_index: 0, feature_index: index, rendered_svg_id: `f${index}`, source_protein_id: `PX${index}.1`
    }, {
      get(target, property, receiver) {
        memberReads += 1;
        return Reflect.get(target, property, receiver);
      }
    });
    const members = Array.from({ length: count }, (_, index) => member(index));
    if (ambiguous) members.push(member(0));
    const index = buildFeatureSearchIndex({
      features: Array.from({ length: count }, (_, index) => ({
        svg_id: `f${index}`, rendered_svg_id: `f${index}`, record_index: 0, feature_index: index,
        type: 'CDS', record_id: 'REC', start: index * 10, end: index * 10 + 9, strand: '+', orthogroup_id: 'OG1'
      })),
      orthogroups: [{ id: 'OG1', name: 'Group one', members }]
    });
    const proteinIds = index.featureOrder.map((svgId) => (
      index.byId.get(svgId).itemsByField.get('orthogroup').find((item) => item.label === 'Protein ID')?.value || ''
    ));
    return { memberReads, proteinIds };
  };
  const small = build(10);
  const large = build(40);
  assert.deepEqual(large.proteinIds, Array.from({ length: 40 }, (_, index) => `PX${index}.1`));
  // Linear in the group size: 4x the members cost about 4x the reads, not 16x.
  assert.ok(large.memberReads < 6 * small.memberReads, `${small.memberReads} reads for 10, ${large.memberReads} for 40`);
  // A feature that matches two members takes no member's protein ID.
  const ambiguous = build(3, { ambiguous: true });
  assert.deepEqual(ambiguous.proteinIds, ['', 'PX1.1', 'PX2.1']);
});
