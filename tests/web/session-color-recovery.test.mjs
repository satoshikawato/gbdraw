import assert from 'node:assert/strict';
import { createHash } from 'node:crypto';
import { readFile } from 'node:fs/promises';
import test from 'node:test';

import {
  parseSpecificRules,
  serializeSpecificRules
} from '../../gbdraw/web/js/app/file-imports.js';
import { resolveOrderedSpecificRule } from '../../gbdraw/web/js/app/feature-editor/fill-view-model.js';
import { recoverV40ColorAuthority } from '../../gbdraw/web/js/services/session-color-recovery.js';

const MAIN_GALLERY_SHA256 = '5e566c305535bb4caccf50723158462471ecf82ce5f71ecc2fda57f91e80f937';
const GALLERY_SESSION = 'gbdraw/web/gallery/sessions/BGC0000708-BGC0000713.gbdraw-session.json';
const TOBACCO_GALLERY_SESSION =
  'gbdraw/web/gallery/sessions/tobacco-chloroplast.gbdraw-session.json';

const resource = (kind, name, value) => ({
  kind,
  name,
  type: 'text/tab-separated-values',
  size: Buffer.byteLength(value),
  lastModified: 0,
  encoding: 'base64',
  data: Buffer.from(value).toString('base64')
});

const rule = {
  feat: 'CDS',
  qual: 'gene_kind',
  val: 'transport',
  color: '#577edb',
  cap: 'Transport-related genes'
};

const minimalV40Session = () => ({
  format: 'gbdraw-session',
  version: 40,
  renderRequest: {
    schema: 5,
    mode: 'linear',
    records: [{ recordKey: 'record-1' }],
    diagramOptions: {
      colors: {
        colorTable: { resourceId: 'canonical-specific', representation: 'canonicalTsv' },
        colorTableFile: null,
        defaultColors: { resourceId: 'canonical-defaults', representation: 'canonicalTsv' },
        defaultColorsFile: null,
        defaultColorsPalette: 'orange'
      }
    },
    comparisons: []
  },
  resources: {
    'canonical-specific': resource(
      'canonical-tsv',
      'canonical-specific.tsv',
      `feature_type\tqualifier_key\tvalue\tcolor\tcaption\n${serializeSpecificRules([rule], { schema: 5 })}`
    ),
    'canonical-defaults': resource(
      'canonical-tsv',
      'canonical-defaults.tsv',
      'feature_type\tcolor\nCDS\t#54bcf8\n'
    ),
    'saved-default-file': resource(
      'colors-default-colors-file',
      'default-colors.tsv',
      'CDS\t#dddddd\n'
    )
  },
  config: {
    palette: 'orange',
    colors: { CDS: '#dddddd' },
    rules: [structuredClone(rule)]
  },
  ui: {
    appliedPaletteName: 'orange',
    appliedPaletteColors: { CDS: '#dddddd' }
  },
  results: [{
    name: 'diagram.svg',
    content: [
      '<svg>',
      '<path fill="#577edb" id="transport-feature"/>',
      '<path fill="#54bcf8" id="default-feature"/>',
      '</svg>'
    ].join('')
  }],
  features: {
    featureColorOverrides: {
      file0_f1: { color: '#577edb', caption: 'Transport-related genes' }
    }
  },
  editorState: {
    legend: {
      entries: [{ caption: 'Transport-related genes', color: '#577edb' }, {
        caption: 'other proteins', color: '#dddddd'
      }]
    },
    featureCatalog: {
      schema: 3,
      items: [{
        resultIndex: 0,
        resultName: 'diagram.svg',
        recordKeys: ['record-1'],
        biologicalFeatures: [{
          recordKey: 'record-1',
          biologicalFeatureId: 'transport-bio',
          stableFeatureId: 'transport-feature',
          sourceFeatureIndex: 1,
          type: 'CDS',
          qualifiers: { gene_kind: ['transport'] }
        }, {
          recordKey: 'record-1',
          biologicalFeatureId: 'default-bio',
          stableFeatureId: 'default-feature',
          sourceFeatureIndex: 2,
          type: 'CDS',
          qualifiers: {}
        }],
        features: [{
          recordKey: 'record-1',
          biologicalFeatureId: 'transport-bio',
          svgId: 'transport-feature',
          fillColor: '#577edb'
        }, {
          recordKey: 'record-1',
          biologicalFeatureId: 'default-bio',
          svgId: 'default-feature',
          fillColor: '#54bcf8'
        }]
      }]
    }
  }
});

test('main-backed v40 Gallery recovery keeps the consistent specific table and repairs only defaults', async () => {
  const raw = await readFile(GALLERY_SESSION);
  assert.equal(createHash('sha256').update(raw).digest('hex'), MAIN_GALLERY_SHA256);
  const session = JSON.parse(raw.toString('utf8'));
  const before = structuredClone(session);
  const resources = session.resources;
  const records = session.renderRequest.records;
  const comparisons = session.renderRequest.comparisons;

  const recovered = recoverV40ColorAuthority(session);

  assert.equal(recovered.recovered, true);
  assert.equal(recovered.specificStatus, 'consistent');
  assert.equal(recovered.defaultStatus, 'recovered');
  assert.equal(recovered.session.renderRequest.schema, 6);
  assert.deepEqual(
    recovered.session.renderRequest.diagramOptions.colors.colorTable,
    before.renderRequest.diagramOptions.colors.colorTable
  );
  assert.equal(recovered.session.renderRequest.diagramOptions.colors.colorTableFile, null);
  assert.equal(recovered.session.renderRequest.diagramOptions.colors.defaultColors, null);
  assert.deepEqual(
    recovered.session.renderRequest.diagramOptions.colors.defaultColorsFile,
    { resourceId: 'colors-default-colors-file', representation: 'file' }
  );
  assert.equal(recovered.session.resources, resources);
  assert.equal(recovered.session.renderRequest.records, records);
  assert.equal(recovered.session.renderRequest.comparisons, comparisons);
  assert.match(
    recovered.session.results[0].content,
    /fill="#dddddd" id="fd4ab5f0c_record_1"/
  );
  assert.deepEqual(session, before, 'recovery must not mutate the v40 envelope or resources');
});

test('v40 Gallery empty editor rules leave the canonical specific table authoritative', async () => {
  const session = JSON.parse(await readFile(TOBACCO_GALLERY_SESSION, 'utf8'));
  const before = structuredClone(session);
  const canonicalSpecificRef = session.renderRequest.diagramOptions.colors.colorTable;

  assert.deepEqual(session.config.rules, []);
  const recovered = recoverV40ColorAuthority(session);

  assert.equal(recovered.specificStatus, 'consistent');
  assert.deepEqual(
    recovered.session.renderRequest.diagramOptions.colors.colorTable,
    canonicalSpecificRef
  );
  assert.equal(recovered.session.renderRequest.diagramOptions.colors.colorTableFile, null);
  assert.deepEqual(session, before, 'compatibility inspection must not mutate Gallery input');

  const missingRules = structuredClone(session);
  delete missingRules.config.rules;
  const recoveredMissingRules = recoverV40ColorAuthority(missingRules);
  assert.equal(recoveredMissingRules.specificStatus, 'consistent');
  assert.deepEqual(
    recoveredMissingRules.session.renderRequest.diagramOptions.colors.colorTable,
    canonicalSpecificRef
  );
});

test('minimal v40 recovery has deterministic copy-on-write behavior', () => {
  const session = minimalV40Session();
  const before = structuredClone(session);
  const resources = session.resources;
  const catalog = session.editorState.featureCatalog;
  const results = session.results;

  const recovered = recoverV40ColorAuthority(session);

  assert.equal(recovered.session.renderRequest.schema, 6);
  assert.equal(recovered.session.resources, resources);
  assert.notEqual(recovered.session.editorState.featureCatalog, catalog);
  assert.notEqual(recovered.session.results, results);
  assert.equal(
    recovered.session.editorState.featureCatalog.items[0].features[1].fillColor,
    '#dddddd'
  );
  assert.match(
    recovered.session.results[0].content,
    /fill="#dddddd" id="default-feature"/
  );
  assert.deepEqual(
    recovered.session.renderRequest.diagramOptions.colors.defaultColorsFile,
    { resourceId: 'saved-default-file', representation: 'file' }
  );
  assert.deepEqual(session, before);
});

test('v40 recovery rejects ambiguous support, conflicting authorities, and unexplained artifacts', () => {
  const ambiguous = minimalV40Session();
  ambiguous.resources['second-saved-default-file'] = resource(
    'colors-default-colors-file',
    'second-default-colors.tsv',
    'CDS\t#dddddd\n'
  );
  assert.throws(
    () => recoverV40ColorAuthority(ambiguous),
    /ambiguous saved resource support/
  );

  const conflictingAuthorities = minimalV40Session();
  conflictingAuthorities.ui.appliedPaletteColors.CDS = '#eeeeee';
  assert.throws(
    () => recoverV40ColorAuthority(conflictingAuthorities),
    /conflicting applied color authorities/
  );

  const conflictingCaption = minimalV40Session();
  conflictingCaption.editorState.legend.entries[0].color = '#000000';
  assert.throws(
    () => recoverV40ColorAuthority(conflictingCaption),
    /conflicting saved legend state/
  );

  const conflictingResult = minimalV40Session();
  conflictingResult.results[0].content = conflictingResult.results[0].content.replace(
    'fill="#577edb"',
    'fill="#000000"'
  );
  assert.throws(
    () => recoverV40ColorAuthority(conflictingResult),
    /conflicting catalogue and Result Feature fills/
  );

  const duplicateCanonicalPair = minimalV40Session();
  const secondItem = structuredClone(duplicateCanonicalPair.editorState.featureCatalog.items[0]);
  secondItem.resultIndex = 1;
  secondItem.resultName = 'second.svg';
  secondItem.features.forEach((feature) => {
    feature.svgId = `second-${feature.svgId}`;
  });
  duplicateCanonicalPair.editorState.featureCatalog.items.push(secondItem);
  duplicateCanonicalPair.results.push({
    name: 'second.svg',
    content: duplicateCanonicalPair.results[0].content
      .replaceAll('id="', 'id="second-')
  });
  assert.throws(
    () => recoverV40ColorAuthority(duplicateCanonicalPair),
    /ambiguous biological Feature identity/
  );
});

test('a derived specific-table mismatch needs one supporting resource and complete rule evidence', () => {
  const session = minimalV40Session();
  session.resources['canonical-specific'] = resource(
    'canonical-tsv',
    'canonical-specific.tsv',
    serializeSpecificRules([{ ...rule, color: '#111111' }], { schema: 5 })
  );
  session.resources['saved-specific-file'] = resource(
    'colors-color-table-file',
    'specific-colors.tsv',
    serializeSpecificRules([rule], { schema: 5 })
  );

  const recovered = recoverV40ColorAuthority(session);

  assert.equal(recovered.specificStatus, 'recovered');
  assert.deepEqual(
    recovered.session.renderRequest.diagramOptions.colors.colorTableFile,
    { resourceId: 'saved-specific-file', representation: 'file' }
  );
});

test('schema-5 biological instance_hash remains regex while the reserved pseudo-selector is rejected', () => {
  const session = minimalV40Session();
  const biologicalRule = {
    feat: 'CDS',
    qual: 'instance_hash',
    val: '^bio-[0-9]+$',
    color: '#577edb',
    cap: 'Biological instance qualifier'
  };
  session.config.rules = [biologicalRule];
  session.resources['canonical-specific'] = resource(
    'canonical-tsv',
    'canonical-specific.tsv',
    serializeSpecificRules([biologicalRule], { schema: 5 })
  );
  session.resources['canonical-defaults'] = resource(
    'canonical-tsv',
    'canonical-defaults.tsv',
    'CDS\t#dddddd\n'
  );
  session.editorState.featureCatalog.items[0].biologicalFeatures[0].qualifiers = {
    instance_hash: ['bio-7']
  };

  const recovered = recoverV40ColorAuthority(session);
  const projectedRules = parseSpecificRules(
    Buffer.from(recovered.session.resources['canonical-specific'].data, 'base64').toString('utf8'),
    { schema: recovered.session.renderRequest.schema }
  ).rules;
  assert.equal(projectedRules[0].qual, 'instance_hash');
  assert.equal(projectedRules[0].match, undefined);
  assert.equal(
    resolveOrderedSpecificRule(
      session.editorState.featureCatalog.items[0].biologicalFeatures[0],
      projectedRules
    )?.ruleIndex,
    0
  );

  const reserved = minimalV40Session();
  const reservedRule = {
    feat: 'CDS',
    qual: '__gbdraw_instance_hash__',
    val: '^fi1_',
    color: '#577edb',
    cap: 'Unsafe legacy interpretation'
  };
  reserved.config.rules = [reservedRule];
  reserved.resources['canonical-specific'] = resource(
    'canonical-tsv',
    'canonical-specific.tsv',
    'CDS\t__gbdraw_instance_hash__\t^fi1_\t#577edb\tUnsafe legacy interpretation\n'
  );
  assert.throws(
    () => recoverV40ColorAuthority(reserved),
    /cannot safely promote a schema-5 reserved instance selector/
  );

  const semantic = minimalV40Session();
  const semanticRule = {
    feat: '*',
    qual: '__gbdraw_semantic_scope__',
    val: '^fs1:',
    color: '#577edb',
    cap: 'Unsafe semantic interpretation'
  };
  semantic.config.rules = [semanticRule];
  semantic.resources['canonical-specific'] = resource(
    'canonical-tsv',
    'canonical-specific.tsv',
    '*\t__gbdraw_semantic_scope__\t^fs1:\t#577edb\tUnsafe semantic interpretation\n'
  );
  assert.throws(
    () => recoverV40ColorAuthority(semantic),
    /cannot safely promote a schema-5 reserved semantic selector/
  );
});
