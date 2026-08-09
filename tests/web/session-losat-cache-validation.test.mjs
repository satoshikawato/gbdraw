import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import test from 'node:test';
import { gunzipSync } from 'node:zlib';

globalThis.window = {
  Vue: {
    ref: (value) => ({ value }),
    reactive: (value) => value,
    computed: (getter) => ({ get value() { return getter(); } }),
    nextTick: async () => {}
  },
  DOMPurify: { sanitize: (value) => value }
};
globalThis.document = {};
let restoredPrimaryTextReads = 0;
let restoredPrimaryArrayBufferReads = 0;
globalThis.File = class File extends Blob {
  constructor(parts, name, options = {}) {
    super(parts, options);
    this.name = String(name || 'file');
    this.lastModified = options.lastModified ?? Date.now();
  }

  async text() {
    restoredPrimaryTextReads += 1;
    return super.text();
  }

  async arrayBuffer() {
    restoredPrimaryArrayBufferReads += 1;
    return super.arrayBuffer();
  }
};
const alerts = [];
globalThis.alert = (message) => alerts.push(String(message));

const {
  adoptCanonicalRenderArtifacts,
  importSession,
  serializeActiveRenderFiles,
  validateSessionLosatArtifacts
} = await import('../../gbdraw/web/js/services/config.js');
const { state } = await import('../../gbdraw/web/js/state.js');

const rawEntry = (key) => ({
  schema: 2,
  kind: 'raw-losat',
  key,
  program: 'blastn',
  identityKind: 'nucleotide',
  text: ''
});

const proteinRawEntry = (key) => ({
  schema: 4,
  kind: 'raw-losat',
  idEncoding: 'runtime-handle-v1',
  key,
  program: 'blastp',
  identityKind: 'protein',
  queryProteinSetHash: 'sha256:query',
  subjectProteinSetHash: 'sha256:subject',
  queryRuntimeBindingHash: 'sha256:query-binding',
  subjectRuntimeBindingHash: 'sha256:subject-binding',
  queryRecordInstanceKey: 'query',
  subjectRecordInstanceKey: 'subject',
  text: ''
});

const derivedEntry = (key) => ({
  schema: 3,
  kind: 'derived-losatp-payload',
  idEncoding: 'runtime-handle-v1',
  key,
  payload: {}
});

const session = (rawEntries, derivedEntries) => ({
  losatCache: { entries: rawEntries },
  losatDerivedCache: { entries: derivedEntries },
  proteinIdentityManifest: {
    schema: 2,
    proteinSets: {},
    recordAnalyses: {},
    recordInstances: {}
  }
});

const canonicalResource = (kind, name, text) => ({
  kind,
  name,
  type: 'text/plain',
  size: new TextEncoder().encode(text).byteLength,
  lastModified: 0,
  encoding: 'base64',
  data: btoa(text)
});

test('canonical protein comparisons rehydrate typed files and pipeline state', async () => {
  alerts.length = 0;
  const genbank = `LOCUS       TEST                        4 bp    DNA     linear   UNK 01-JAN-1980
FEATURES             Location/Qualifiers
ORIGIN
        1 atgc
//
`;
  const proteinTable = `query\tsubject\tidentity\talignment_length\tmismatches\tgap_opens\tqstart\tqend\tsstart\tsend\tevalue\tbitscore
protein-a\tprotein-b\t95\t20\t1\t0\t10\t30\t50\t70\t1e-20\t120
`;
  const orthogroupsJson = JSON.stringify({
    schema: 1,
    valueKind: 'orthogroupResult',
    value: { orthogroups: {} }
  });
  const collinearityJson = JSON.stringify({
    schema: 1,
    valueKind: 'blocks',
    value: []
  });
  const generatedProteinComparison = {
    kind: 'generatedProteinComparison',
    mode: 'none',
    pairs: [],
    settings: {
      collinearityParams: {
        kind: 'lossless',
        parameters: {
          minAnchors: 3,
          maxUnitGap: 4,
          maxDiagonalDrift: 5,
          maxConflicts: 6,
          mergeOrientation: 'either'
        }
      },
      collinearityUnitMode: 'cds',
      collinearityAnchorMode: 'rbh',
      collinearitySearchScope: 'all',
      collinearityColorMode: 'average_identity',
      losatpBin: '/custom/losat',
      ncbiBlastpBin: '/custom/blastp',
      losatpThreads: 2,
      proteinBlastpMaxHits: 9,
      proteinBlastpCandidateLimit: 17,
      orthogroupMembershipMode: 'anchor_core_v1',
      orthogroupMemberMaxHits: 7,
      collinearMaxParalogLinksPerOrthogroup: 8,
      alignOrthogroupFeature: 'feature-anchor'
    }
  };
  const renderRequest = {
    schema: 5,
    mode: 'linear',
    grouping: 'single',
    records: [0, 1, 2].map((index) => ({
      recordKey: `record-${index + 1}`,
      source: { kind: 'genbank', resourceId: `record-${index + 1}` },
      selector: null,
      region: null,
      presentation: {
        label: null,
        subtitle: null,
        reverseComplement: false,
        gridRow: null,
        gridColumn: null
      }
    })),
    diagramOptions: {
      configOverrides: {},
      colors: {
        colorTable: null,
        colorTableFile: null,
        defaultColors: null,
        defaultColorsPalette: 'default',
        defaultColorsFile: null
      },
      selectedFeaturesSet: ['CDS'],
      featureShapes: { CDS: 'arrow' },
      plotTitle: '',
      evalue: 1e-5,
      bitscore: 50,
      identity: 70,
      alignmentLength: 0,
      tracks: {
        circularTrackSlots: null,
        circularTrackAxisIndex: null,
        linearTrackSlots: null,
        linearTrackAxisIndex: null,
        centerReservedRadius: null
      },
      output: { legend: 'bottom', plotTitlePosition: 'bottom' },
      pairwiseMatchStyle: 'ribbon'
    },
    layout: {},
    comparisons: [
      {
        kind: 'precomputedProteinComparison',
        resourceId: 'resolved-protein',
        encoding: 'canonicalTsv',
        queryRecordIndex: 0,
        subjectRecordIndex: 2
      },
      {
        kind: 'orthogroupResult',
        resourceId: 'orthogroups',
        encoding: 'canonicalJson'
      },
      {
        kind: 'collinearityResult',
        resourceId: 'collinearity',
        encoding: 'canonicalJson',
        valueKind: 'blocks'
      },
      generatedProteinComparison
    ],
    output: {
      prefix: 'typed-round-trip',
      formats: ['interactive_svg'],
      overwrite: false,
      interactiveMetadataPolicy: 'auto'
    }
  };
  const resources = {
    'record-1': canonicalResource('genbank', 'record-1.gb', genbank),
    'record-2': canonicalResource('genbank', 'record-2.gb', genbank),
    'record-3': canonicalResource('genbank', 'record-3.gb', genbank),
    'resolved-protein': canonicalResource(
      'canonical-tsv',
      'resolved-protein.tsv',
      proteinTable
    ),
    orthogroups: canonicalResource(
      'orthogroup-result',
      'orthogroups.json',
      orthogroupsJson
    ),
    collinearity: canonicalResource(
      'collinearity-result',
      'collinearity.json',
      collinearityJson
    )
  };
  const payload = {
    format: 'gbdraw-session',
    version: 39,
    renderRequest,
    resources,
    config: {
      losatProgram: 'blastn',
      losat: {
        outfmt: '6',
        parallelWorkers: '3',
        executionMode: 'threaded',
        totalThreadBudget: '12',
        threadsPerJob: '4',
        blastn: { task: 'dc-megablast' },
        blastp: { mode: 'collinear', maxHits: 1 }
      }
    },
    losatCache: { entries: [] },
    losatDerivedCache: { entries: [] },
    proteinIdentityManifest: {
      schema: 2,
      proteinSets: {},
      recordAnalyses: {},
      recordInstances: {}
    }
  };
  const file = new Blob([JSON.stringify(payload)], { type: 'application/json' });
  const event = { target: { files: [file], value: 'selected' } };

  const result = await importSession(event);

  assert.equal(result.status, 'ok');
  assert.deepEqual(
    state.files.linearCanonicalComparisons.map((comparison) => comparison.kind),
    [
      'precomputedProteinComparison',
      'orthogroupResult',
      'collinearityResult',
      'generatedProteinComparison'
    ]
  );
  assert.equal(state.linearComparisonPlan.edges.length, 0);
  assert.equal(
    state.files.linearCanonicalComparisons[0].file.name,
    'resolved-protein.tsv'
  );
  assert.deepEqual(
    state.files.linearCanonicalComparisons[3],
    generatedProteinComparison
  );
  assert.equal(
    state.files.linearCanonicalComparisons[1].file.name,
    'orthogroups.json'
  );
  assert.equal(
    state.files.linearCanonicalComparisons[2].valueKind,
    'blocks'
  );
  assert.equal(state.losatProgram.value, 'blastp');
  assert.equal(state.losat.executionMode, 'threaded');
  assert.equal(state.losat.parallelWorkers, '3');
  assert.equal(state.losat.totalThreadBudget, '12');
  assert.equal(state.losat.blastn.task, 'dc-megablast');
  assert.equal(state.losat.blastp.mode, 'collinear');
  assert.equal(state.losat.blastp.maxHits, 9);
  assert.equal(state.losat.blastp.collinearSearchScope, 'all');
  assert.equal(state.selectedOrthogroupAlignmentFeature.value, 'feature-anchor');

  state.files.linearCanonicalComparisons = [];
  adoptCanonicalRenderArtifacts({ renderRequest, resources });
  assert.deepEqual(
    state.files.linearCanonicalComparisons.map((comparison) => comparison.kind),
    [
      'precomputedProteinComparison',
      'orthogroupResult',
      'collinearityResult',
      'generatedProteinComparison'
    ]
  );

  const serialized = await serializeActiveRenderFiles(state.mode.value, state);
  assert.equal(serialized.linearCanonicalComparisons[0].file.data, btoa(proteinTable));
  assert.equal(
    serialized.linearCanonicalComparisons[1].file.data,
    btoa(orthogroupsJson)
  );
  assert.equal(
    serialized.linearCanonicalComparisons[2].file.data,
    btoa(collinearityJson)
  );
  assert.equal(serialized.linearCanonicalComparisons[2].valueKind, 'blocks');
  assert.deepEqual(
    serialized.linearCanonicalComparisons[3],
    generatedProteinComparison
  );
});

test('current session restores match sequences from the catalog without rereading primary files', async () => {
  alerts.length = 0;
  restoredPrimaryTextReads = 0;
  restoredPrimaryArrayBufferReads = 0;
  const genbank = `LOCUS       CATALOG                     8 bp    DNA     linear   UNK 01-JAN-1980
ACCESSION   CATALOG
VERSION     CATALOG.1
FEATURES             Location/Qualifiers
ORIGIN
        1 aaccggtt
//
`;
  const renderRequest = {
    schema: 5,
    mode: 'linear',
    grouping: 'single',
    records: [{
      recordKey: 'catalog-record',
      source: { kind: 'genbank', resourceId: 'catalog-record' },
      selector: null,
      region: null,
      presentation: {
        label: null,
        subtitle: null,
        reverseComplement: false,
        gridRow: null,
        gridColumn: null
      }
    }],
    diagramOptions: {
      configOverrides: {},
      colors: {
        colorTable: null,
        colorTableFile: null,
        defaultColors: null,
        defaultColorsPalette: 'default',
        defaultColorsFile: null
      },
      selectedFeaturesSet: ['CDS'],
      featureShapes: { CDS: 'arrow' },
      plotTitle: '',
      evalue: 1e-5,
      bitscore: 50,
      identity: 70,
      alignmentLength: 0,
      tracks: {
        circularTrackSlots: null,
        circularTrackAxisIndex: null,
        linearTrackSlots: null,
        linearTrackAxisIndex: null,
        centerReservedRadius: null
      },
      output: { legend: 'bottom', plotTitlePosition: 'bottom' },
      pairwiseMatchStyle: 'ribbon'
    },
    layout: {},
    comparisons: [],
    output: {
      prefix: 'catalog-session',
      formats: ['svg'],
      overwrite: false,
      interactiveMetadataPolicy: 'omit'
    }
  };
  const sequenceSource = {
    key: 'linear:record:0',
    recordId: 'CATALOG.1',
    aliases: ['CATALOG', 'CATALOG.1'],
    sequence: 'AACCGGTT',
    origin: 'linear-record',
    recordIndex: 0
  };
  const payload = {
    format: 'gbdraw-session',
    version: 40,
    renderRequest,
    resources: {
      'catalog-record': canonicalResource(
        'genbank',
        'catalog-primary.gbk',
        genbank
      )
    },
    webFiles: {
      bindings: {
        schema: 1,
        linearSeqs: [{
          uid: 'catalog-record',
          gb: {
            resourceId: 'catalog-record',
            name: 'catalog-primary.gbk',
            type: 'application/genbank',
            lastModified: 7
          },
          gff: null,
          fasta: null,
          depth: null,
          losat_gencode: 1,
          definition: '',
          record_subtitle: '',
          region_record_id: '',
          region_start: null,
          region_end: null,
          region_reverse: false
        }]
      }
    },
    config: {
      form: { legend: 'bottom' },
      adv: {
        circular_track_slots_enabled: false,
        circular_track_slots: [],
        circular_track_slots_axis_index: null,
        circular_track_slots_schema_version: 4,
        linear_track_slots_enabled: false,
        linear_track_slots: [],
        linear_track_slots_axis_index: null,
        linear_track_slots_schema_version: 2
      }
    },
    ui: { mode: 'linear', selectedResultIndex: 0 },
    results: [{
      name: 'catalog-session',
      content: '<svg xmlns="http://www.w3.org/2000/svg"></svg>'
    }],
    features: {},
    editorState: {
      featureCatalog: {
        schema: 3,
        items: [{
          resultIndex: 0,
          resultName: 'catalog-session',
          recordKeys: ['catalog-record'],
          features: [],
          biologicalFeatures: [],
          orthogroups: [],
          annotations: [],
          comparisonMatches: [],
          sequenceSources: [sequenceSource]
        }]
      }
    },
    orthogroupState: {},
    losatCache: { entries: [] },
    losatDerivedCache: { entries: [] },
    proteinIdentityManifest: {
      schema: 2,
      proteinSets: {},
      recordAnalyses: {},
      recordInstances: {}
    }
  };
  const file = new Blob([JSON.stringify(payload)], { type: 'application/json' });
  const event = { target: { files: [file], value: 'selected' } };

  const result = await importSession(event);

  assert.equal(result.status, 'ok');
  assert.equal(state.linearSeqs[0].gb.name, 'catalog-primary.gbk');
  assert.equal(restoredPrimaryTextReads, 0);
  assert.equal(restoredPrimaryArrayBufferReads, 0);
  assert.deepEqual(state.matchSequenceRegistry.values(), [{
    ...sequenceSource,
    aliases: ['CATALOG.1', 'CATALOG'],
    sourceIndex: null
  }]);

  const secondGenbank = genbank
    .replaceAll('CATALOG', 'SECOND')
    .replace('aaccggtt', 'ttggccaa');
  const incompletePayload = structuredClone(payload);
  incompletePayload.renderRequest.records.push({
    ...structuredClone(renderRequest.records[0]),
    recordKey: 'second-record',
    source: { kind: 'genbank', resourceId: 'second-record' }
  });
  incompletePayload.resources['second-record'] = canonicalResource(
    'genbank',
    'second-primary.gbk',
    secondGenbank
  );
  incompletePayload.webFiles.bindings.linearSeqs.push({
    ...structuredClone(incompletePayload.webFiles.bindings.linearSeqs[0]),
    uid: 'second-record',
    gb: {
      resourceId: 'second-record',
      name: 'second-primary.gbk',
      type: 'application/genbank',
      lastModified: 8
    }
  });
  incompletePayload.editorState.featureCatalog.items[0].recordKeys.push(
    'second-record'
  );
  restoredPrimaryTextReads = 0;
  const incompleteResult = await importSession({
    target: {
      files: [new Blob([JSON.stringify(incompletePayload)], { type: 'application/json' })],
      value: 'selected'
    }
  });

  assert.equal(incompleteResult.status, 'ok');
  assert.equal(restoredPrimaryTextReads, 2);
  assert.deepEqual(
    state.matchSequenceRegistry.values().map((source) => source.recordId),
    ['CATALOG.1', 'SECOND.1']
  );
});

test('fresh Circular LOSAT resources become session-owned live files', async () => {
  const genbank = `LOCUS       TEST                        4 bp    DNA     linear   UNK 01-JAN-1980
FEATURES             Location/Qualifiers
ORIGIN
        1 atgc
//
`;
  const fasta = '>comparison\natgc\n';
  const resources = {
    record: canonicalResource('genbank', 'record.gb', genbank),
    'blast-a': canonicalResource('conservation-blast-file', 'comparison-a.tsv', ''),
    'blast-b': canonicalResource(
      'conservation-blast-file',
      'comparison-b.tsv',
      'record\tcomparison\t99\t4\t0\t0\t1\t4\t1\t4\t1e-20\t50\n'
    ),
    fasta: canonicalResource('conservation-fasta-file', 'comparison.fna', fasta)
  };
  const renderRequest = {
    schema: 5,
    mode: 'circular',
    grouping: 'single',
    records: [{
      recordKey: 'record',
      source: { kind: 'genbank', resourceId: 'record' },
      selector: null,
      region: null,
      presentation: {
        label: null,
        subtitle: null,
        reverseComplement: false,
        gridRow: null,
        gridColumn: null
      }
    }],
    diagramOptions: {
      configOverrides: {},
      colors: {
        colorTable: null,
        colorTableFile: null,
        defaultColors: null,
        defaultColorsPalette: 'default',
        defaultColorsFile: null
      },
      selectedFeaturesSet: ['CDS'],
      featureShapes: { CDS: 'arrow' },
      plotTitle: '',
      evalue: 1e-5,
      bitscore: 50,
      identity: 70,
      alignmentLength: 0,
      tracks: {
        circularTrackSlots: null,
        circularTrackAxisIndex: null,
        linearTrackSlots: null,
        linearTrackAxisIndex: null,
        centerReservedRadius: null
      },
      output: { legend: 'left', plotTitlePosition: 'none' },
      conservationBlastFiles: [
        { resourceId: 'blast-a', representation: 'file' },
        { resourceId: 'blast-b', representation: 'file' }
      ],
      conservationFastaFiles: [
        null,
        { resourceId: 'fasta', representation: 'file' }
      ],
      conservationReference: 'subject',
      conservationLabels: ['Missing FASTA', 'Saved comparison'],
      conservationColors: ['#e15759', '#4e79a7'],
      conservationRingWidth: 14,
      conservationRingGap: 3
    },
    layout: {},
    comparisons: [],
    output: {
      prefix: 'circular-losat',
      formats: ['interactive_svg'],
      overwrite: false,
      interactiveMetadataPolicy: 'auto'
    }
  };
  state.circularConservation.series.splice(
    0,
    state.circularConservation.series.length,
    {
      sourceIndex: 4,
      label: 'Old label',
      color: '#000000',
      losat_gencode: 11
    },
    {
      sourceIndex: 9,
      label: 'Second old label',
      color: '#111111',
      losat_gencode: 4
    }
  );

  adoptCanonicalRenderArtifacts({
    renderRequest,
    resources,
    webFiles: {
      conservationBlastSource: 'losat-cache',
      conservationLosatFastaSources: [null, 'fasta']
    }
  });

  assert.equal(state.files.c_conservation_blasts_source, 'losat-cache');
  assert.equal(state.files.c_conservation_blasts.length, 2);
  assert.equal(state.files.c_conservation_blasts[0].size, 0);
  assert.equal(state.files.c_conservation_fastas.length, 2);
  assert.equal(state.files.c_conservation_fastas[0], null);
  assert.equal(state.files.c_conservation_fastas[1].name, 'comparison.fna');
  assert.equal(state.circularConservation.source, 'losat');
  assert.deepEqual(
    {
      sourceIndex: state.circularConservation.series[0].sourceIndex,
      label: state.circularConservation.series[0].label,
      color: state.circularConservation.series[0].color,
      losatGencode: state.circularConservation.series[0].losat_gencode
    },
    {
      sourceIndex: 0,
      label: 'Missing FASTA',
      color: '#e15759',
      losatGencode: 11
    }
  );
  assert.equal(state.circularConservation.series[1].sourceIndex, 1);
  assert.equal(state.circularConservation.series[1].losat_gencode, 4);

  state.mode.value = 'circular';
  const serialized = await serializeActiveRenderFiles(state.mode.value, state);
  assert.equal(serialized.c_conservation_blasts[0].data, '');
  assert.equal(serialized.c_conservation_fastas[0], null);
  assert.equal(serialized.c_conservation_fastas[1].data, btoa(fasta));
});

test('legacy standalone config import migrates retired values without a writer envelope', async () => {
  alerts.length = 0;
  const file = new Blob([JSON.stringify({
    form: {
      prefix: 'legacy-standalone',
      linear_track_layout: 'tuckin'
    },
    adv: {
      label_placement: 'on_feature',
      multi_record_size_mode: 'sqrt'
    }
  })], { type: 'application/json' });
  const event = { target: { files: [file], value: 'selected' } };

  const result = await importSession(event);

  assert.equal(result.status, 'legacy');
  assert.equal(state.form.prefix, 'legacy-standalone');
  assert.equal(state.form.linear_track_layout, 'below');
  assert.equal(state.adv.label_placement, 'above_feature');
  assert.equal(state.adv.multi_record_size_mode, 'auto');
  assert.deepEqual(alerts, [
    'Legacy configuration loaded. Save as a session to use the current format.'
  ]);
});

const frozenV39Session = () => JSON.parse(gunzipSync(readFileSync(new URL(
  '../fixtures/sessions/BGC0000708-BGC0000713.v39.gbdraw-session.json.gz',
  import.meta.url
))));

test('visible protein cache entries from the frozen v39 session recover stable edge identities', async () => {
  alerts.length = 0;
  const sessionData = frozenV39Session();
  const sourcePairs = sessionData.losatCache.entries
    .filter((entry) => entry.display !== false)
    .map((entry) => [entry.queryRecordInstanceKey, entry.subjectRecordInstanceKey]);
  const file = new Blob([JSON.stringify(sessionData)], { type: 'application/json' });
  const result = await importSession({ target: { files: [file], value: 'selected' } });

  assert.equal(result.status, 'ok');
  assert.deepEqual(
    state.losatCacheInfo.value.map((entry) => [
      entry.edgeKey,
      entry.queryUid,
      entry.subjectUid,
      entry.queryIndex,
      entry.subjectIndex,
      entry.ordinal
    ]),
    sourcePairs.map(([queryUid, subjectUid], ordinal) => [
      `${queryUid}->${subjectUid}`,
      queryUid,
      subjectUid,
      ordinal,
      ordinal + 1,
      ordinal
    ])
  );
  const firstRawEntry = state.losatCache.value.get(state.losatCacheInfo.value[0].key);
  assert.deepEqual(
    [firstRawEntry.queryRecordInstanceKey, firstRawEntry.subjectRecordInstanceKey],
    sourcePairs[0]
  );
});

test('a conflicting restored edge key is not exposed under the wrong endpoints', async () => {
  alerts.length = 0;
  const sessionData = frozenV39Session();
  const visibleEntries = sessionData.losatCache.entries.filter(
    (entry) => entry.display !== false
  );
  const conflictingEdgeKey = [
    visibleEntries[1].queryRecordInstanceKey,
    visibleEntries[1].subjectRecordInstanceKey
  ].join('->');
  visibleEntries[0].edgeKey = conflictingEdgeKey;
  const file = new Blob([JSON.stringify(sessionData)], { type: 'application/json' });
  const result = await importSession({ target: { files: [file], value: 'selected' } });

  assert.equal(result.status, 'ok');
  assert.equal(state.losatCacheInfo.value[0].edgeKey, undefined);
  assert.equal(
    state.losatCacheInfo.value.filter((entry) => entry.edgeKey === conflictingEdgeKey).length,
    1
  );
});

for (const version of [39, 40]) {
  test(`version-${version} import rejects duplicate raw LOSAT cache keys`, () => {
    assert.throws(
      () => validateSessionLosatArtifacts(
        session([rawEntry('raw-key'), rawEntry('raw-key')], []),
        version
      ),
      /Duplicate LOSAT cache key/
    );
  });

  test(`version-${version} import rejects duplicate derived LOSATP cache keys`, () => {
    assert.throws(
      () => validateSessionLosatArtifacts(
          session(
            [],
            [derivedEntry('derived-key'), derivedEntry('derived-key')]
        ),
        version
      ),
      /Duplicate derived LOSATP cache key/
    );
  });

  for (const invalidKey of [undefined, '', null, 1]) {
    test(`version-${version} import rejects a nucleotide raw cache key of ${String(invalidKey)}`, () => {
      assert.throws(
        () => validateSessionLosatArtifacts(
          session([rawEntry(invalidKey)], []),
          version
        ),
        /LOSAT cache entry at losatCache\.entries\[0\] requires a key/
      );
    });
  }

  test(`version-${version} import rejects an empty protein raw cache key`, () => {
    assert.throws(
      () => validateSessionLosatArtifacts(
        session([proteinRawEntry('')], []),
        version
      ),
      /LOSAT cache entry at losatCache\.entries\[0\] requires a key/
    );
  });

  for (const field of ['losatCache', 'losatDerivedCache']) {
    test(`version-${version} import rejects a non-object ${field}`, () => {
      const malformed = session([], []);
      malformed[field] = [];
      assert.throws(
        () => validateSessionLosatArtifacts(malformed, version),
        new RegExp(`Session ${field} must be an object when present`)
      );
    });

    test(`version-${version} import rejects a non-array ${field}.entries`, () => {
      const malformed = session([], []);
      malformed[field] = { entries: null };
      assert.throws(
        () => validateSessionLosatArtifacts(malformed, version),
        new RegExp(`Session ${field}\\.entries must be an array`)
      );
    });
  }
  test(`version-${version} import rejects branch-internal raw schema 3`, () => {
    const malformed = session([{ ...proteinRawEntry('raw-key'), schema: 3 }], []);
    assert.throws(
      () => validateSessionLosatArtifacts(malformed, version),
      new RegExp(`Session version ${version} contains a non-current raw LOSAT entry`)
    );
  });

  test(`version-${version} import rejects branch-internal derived schema 2`, () => {
    const malformed = session([], [{ ...derivedEntry('derived-key'), schema: 2 }]);
    assert.throws(
      () => validateSessionLosatArtifacts(malformed, version),
      new RegExp(`Session version ${version} contains an invalid derived LOSATP entry`)
    );
  });
}

for (const unsupportedVersion of [34, 35, 36, 37, 38]) {
  test(`session import rejects branch-internal version ${unsupportedVersion}`, async () => {
    alerts.length = 0;
    const file = new Blob([JSON.stringify({
      format: 'gbdraw-session',
      version: unsupportedVersion
    })], { type: 'application/json' });
    const event = { target: { files: [file], value: 'selected' } };
    const result = await importSession(event);
    assert.equal(result.status, 'error');
    assert.match(
      String(result.error?.message || ''),
      new RegExp(`Unsupported session version: ${unsupportedVersion}`)
    );
    assert.deepEqual(alerts, [
      `Failed to load session: Unsupported session version: ${unsupportedVersion}.`
    ]);
    assert.equal(event.target.value, '');
  });
}
